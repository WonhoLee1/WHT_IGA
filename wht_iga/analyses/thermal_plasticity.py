# -*- coding: utf-8 -*-
"""Thermo-plastic precursor for E03 on the S28 path (Track T2, spike S30).

Kratos 10.4 has no thermo-plastic law with thermal expansion (S28). The thermal
strain is therefore imposed as ``InitialState`` strain on the generic J2 laws
(``SmallStrainIsotropicPlasticity3DVonMisesVonMises`` /
``FiniteStrainIsotropicPlasticity3DVonMisesVonMises``) or on the linear elastic
law (thermoelastic glass), and the temperature-dependent ``YOUNG_MODULUS`` and
``YIELD_STRESS`` are written into per-group ``Properties`` at every level.
Both are re-read by the laws on every call (S28 for the yield stress; the
elastic matrix is rebuilt from the Properties in ``CalculateElasticMatrix``),
so the stress is the total form ``sigma = D(T) : (eps - eps_th(T) - eps_p)``,
which is also Abaqus' form for temperature-dependent elasticity.

Pieces
------
* :class:`TempTable` - piecewise-linear property vs temperature, constant
  beyond the end points (Abaqus table behaviour).
* :class:`ThermoMaterial` - E(T), nu, secant expansion alpha(T) about the
  zero temperature, yield(T) or ``None`` (thermoelastic).
  ``thermal_strain(T, T_init) = alpha(T)(T - T_zero) - alpha(T_init)(T_init - T_zero)``
  (Abaqus ``*EXPANSION, ZERO=`` with an initial-temperature field: the
  initial state is strain free).
* :class:`CellThermalGroups` - groups elements (= cells; one element per cell
  from ``attach_cell_elements``) by (material, rounded temperature history
  over ALL thermal states). Cells with the same key have the same thermal
  strain and properties at every level, so one sub model part + one
  ``Properties`` per group is built ONCE and reused for the whole history; a
  level costs one ``SetInitialStateProcess3D`` run and two ``SetValue`` per
  group (``Geometry.SetValue`` is not bound in Python, S28).
* :func:`run_thermal_steps` - thermal-step driver on top of
  :func:`wht_iga.analyses.nonlinear.run_load_stepping`: level ``s`` in
  ``[k, k+1]`` is the linear interpolation between thermal states ``k`` and
  ``k+1``; one strategy for the whole history.

Limits (S30): the temperature is constant per element (the process writes the
same state to every integration point of an element); hardening is perfect
plasticity (``HARDENING_CURVE = 3``) - the generic laws harden on plastic
dissipation, not on PEEQ, so E03's linear / multilinear hardening is NOT
represented; groups and Properties must be built before the strategy is
initialized (the elements clone ``CONSTITUTIVE_LAW`` from their Properties in
``Initialize``).
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend.geometry import cell_gauss_points

__all__ = [
    "TempTable",
    "ThermoMaterial",
    "CellThermalGroups",
    "default_law_factory",
    "thermal_levels",
    "run_thermal_steps",
    "e03_materials_initial_yield",
    "attach_region_elements",
    "fix_isostatic_321",
]


# --------------------------------------------------------------------------- data

@dataclass(frozen=True)
class TempTable:
    """Property vs temperature, linear between points, constant outside (Abaqus)."""

    temperatures: tuple[float, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        t = np.asarray(self.temperatures, dtype=float)
        v = np.asarray(self.values, dtype=float)
        if t.ndim != 1 or t.size == 0 or t.shape != v.shape:
            raise ValueError("TempTable needs equally long, non-empty 1D temperatures/values")
        if np.any(np.diff(t) <= 0.0):
            raise ValueError("TempTable temperatures must be strictly increasing")
        object.__setattr__(self, "temperatures", tuple(float(x) for x in t))
        object.__setattr__(self, "values", tuple(float(x) for x in v))

    def __call__(self, temperature):
        return np.interp(temperature, self.temperatures, self.values)


def _eval(prop: float | TempTable, temperature):
    if isinstance(prop, TempTable):
        return prop(temperature)
    if np.ndim(temperature):
        return np.full(np.shape(temperature), float(prop))
    return float(prop)


@dataclass(frozen=True)
class ThermoMaterial:
    """Temperature-dependent isotropic material (E03 data forms).

    ``expansion`` is the secant (mean) coefficient about ``t_zero`` (Abaqus
    ``ZERO`` temperature, default 0 as in Abaqus). ``yield_stress=None`` means
    thermoelastic (linear elastic law, e.g. the S-8061 glass).
    """

    name: str
    young_modulus: float | TempTable
    poisson_ratio: float
    expansion: float | TempTable
    t_zero: float = 0.0
    yield_stress: float | TempTable | None = None

    @property
    def plastic(self) -> bool:
        return self.yield_stress is not None

    def young(self, temperature):
        return _eval(self.young_modulus, temperature)

    def yield_at(self, temperature):
        if self.yield_stress is None:
            raise ValueError(f"material {self.name!r} is thermoelastic (no yield stress)")
        return _eval(self.yield_stress, temperature)

    def thermal_strain(self, temperature, t_init):
        """alpha(T)(T - T_zero) - alpha(T_init)(T_init - T_zero); zero at T = T_init."""
        a = _eval(self.expansion, temperature)
        a0 = _eval(self.expansion, t_init)
        return a * (np.asarray(temperature, dtype=float) - self.t_zero) \
            - a0 * (np.asarray(t_init, dtype=float) - self.t_zero)


def default_law_factory(mods: dict[str, Any], finite_strain: bool = False
                        ) -> Callable[[ThermoMaterial], Any]:
    """Law prototype per material: linear elastic if thermoelastic, else
    generic von Mises J2 (small or finite strain). ``mods`` from ``load_kratos``."""
    sma = mods["StructuralMechanicsApplication"]
    cla = mods["ConstitutiveLawsApplication"]

    def factory(mat: ThermoMaterial):
        if not mat.plastic:
            return sma.LinearElastic3DLaw()
        if finite_strain:
            return cla.FiniteStrainIsotropicPlasticity3DVonMisesVonMises()
        return cla.SmallStrainIsotropicPlasticity3DVonMisesVonMises()

    return factory


# --------------------------------------------------------------------------- groups

class CellThermalGroups:
    """Cell groups with equal (material, temperature history) -> one sub model
    part + one ``Properties`` each, reused for every thermal level.

    Args:
        KM, mods: Kratos core module and the ``load_kratos`` dict.
        mp: root model part holding ``elements``.
        elements: Kratos elements (one per cell), not yet initialized.
        materials: material list; ``material_index`` (scalar or per element)
            indexes it.
        states: thermal states, ``states[0]`` = stress-free initial
            temperature; each entry a scalar or a per-element array.
        law_factory: ``ThermoMaterial -> law prototype`` (default:
            :func:`default_law_factory` small strain).
        quantum: temperature rounding for grouping; the group temperature is
            the member mean, so the per-cell temperature error is < quantum
            (<= quantum at bin-edge ties).
        fracture_energy: ``FRACTURE_ENERGY`` for the plastic laws (required by
            Kratos, irrelevant for perfect plasticity).
    """

    def __init__(self, KM: Any, mods: dict[str, Any], mp: Any, elements: Sequence[Any],
                 materials: Sequence[ThermoMaterial], material_index,
                 states: Sequence[Any],
                 law_factory: Callable[[ThermoMaterial], Any] | None = None,
                 quantum: float = 1e-6, name_prefix: str = "ThermalGroup",
                 fracture_energy: float = 1e20) -> None:
        if len(states) < 1:
            raise ValueError("need at least the initial thermal state")
        if quantum <= 0.0:
            raise ValueError("quantum must be positive")
        t0 = time.perf_counter()
        self._KM = KM
        self._cl = mods["ConstitutiveLawsApplication"]
        self.mp = mp
        self.materials = list(materials)
        n_el = len(elements)
        self.n_states = len(states)
        hist = np.empty((self.n_states, n_el))
        for s, st in enumerate(states):
            hist[s] = np.broadcast_to(np.asarray(st, dtype=float), (n_el,))
        mat_idx = np.broadcast_to(np.asarray(material_index, dtype=np.int64), (n_el,))
        if mat_idx.size and (mat_idx.min() < 0 or mat_idx.max() >= len(self.materials)):
            raise ValueError("material_index out of range")

        keys = np.column_stack([mat_idx, np.round(hist.T / quantum).astype(np.int64)])
        _, inverse = np.unique(keys, axis=0, return_inverse=True)
        inverse = np.asarray(inverse).reshape(-1)
        self.n_groups = int(inverse.max()) + 1 if n_el else 0
        counts = np.bincount(inverse, minlength=self.n_groups)
        self.group_T = np.vstack([np.bincount(inverse, weights=hist[s], minlength=self.n_groups)
                                  / counts for s in range(self.n_states)])
        self.group_material = np.zeros(self.n_groups, dtype=np.int64)
        self.group_material[inverse] = mat_idx
        self.group_of_element = {int(el.Id): int(g) for el, g in zip(elements, inverse)}

        factory = law_factory if law_factory is not None else default_law_factory(mods)
        laws = [factory(m) for m in self.materials]
        pid0 = max([int(p.Id) for p in mp.Properties] + [0]) + 1
        order = np.argsort(inverse, kind="stable")
        bounds = np.concatenate([[0], np.cumsum(counts)])
        self.sub_model_parts: list[Any] = []
        self.properties: list[Any] = []
        for g in range(self.n_groups):
            mat = self.materials[self.group_material[g]]
            props = mp.CreateNewProperties(pid0 + g)
            props.SetValue(KM.POISSON_RATIO, float(mat.poisson_ratio))
            props.SetValue(KM.CONSTITUTIVE_LAW, laws[self.group_material[g]])
            if mat.plastic:
                props.SetValue(KM.FRACTURE_ENERGY, float(fracture_energy))
                props.SetValue(self._cl.HARDENING_CURVE, 3)  # perfect plasticity
            members = [elements[e] for e in order[bounds[g]:bounds[g + 1]]]
            for el in members:
                el.Properties = props
            smp = mp.CreateSubModelPart(f"{name_prefix}_{g}")
            smp.AddElements([int(el.Id) for el in members])
            self.sub_model_parts.append(smp)
            self.properties.append(props)
        self._zero6 = KM.Vector(6)
        for i in range(6):
            self._zero6[i] = 0.0
        self._eye3 = KM.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                self._eye3[i, j] = 1.0 if i == j else 0.0
        self.build_seconds = time.perf_counter() - t0
        self.apply_seconds = 0.0
        self.n_apply = 0
        # set the level-0 properties so the laws can be initialized consistently
        self._set_properties(self.group_values(0.0))

    @property
    def n_steps(self) -> int:
        return self.n_states - 1

    def temperatures(self, level: float) -> np.ndarray:
        """Group temperatures at ``level`` (linear between thermal states)."""
        s = float(level)
        if self.n_steps == 0:
            if abs(s) > 1e-12:
                raise ValueError("only one thermal state: level must be 0")
            return self.group_T[0].copy()
        if s < -1e-12 or s > self.n_steps + 1e-12:
            raise ValueError(f"level {s} outside [0, {self.n_steps}]")
        k = min(max(int(math.floor(s)), 0), self.n_steps - 1)
        w = min(max(s - k, 0.0), 1.0)
        return (1.0 - w) * self.group_T[k] + w * self.group_T[k + 1]

    def group_values(self, level: float) -> dict[str, np.ndarray]:
        """Per-group ``T``, ``eps_th``, ``E``, ``sigma_y`` (NaN if elastic) at ``level``."""
        T = self.temperatures(level)
        eps = np.empty(self.n_groups)
        E = np.empty(self.n_groups)
        sy = np.full(self.n_groups, np.nan)
        for m, mat in enumerate(self.materials):
            sel = self.group_material == m
            if not np.any(sel):
                continue
            eps[sel] = mat.thermal_strain(T[sel], self.group_T[0][sel])
            E[sel] = mat.young(T[sel])
            if mat.plastic:
                sy[sel] = mat.yield_at(T[sel])
        return {"T": T, "eps_th": eps, "E": E, "sigma_y": sy}

    def _set_properties(self, vals: dict[str, np.ndarray]) -> None:
        KM = self._KM
        for g, props in enumerate(self.properties):
            props.SetValue(KM.YOUNG_MODULUS, float(vals["E"][g]))
            if not math.isnan(vals["sigma_y"][g]):
                props.SetValue(KM.YIELD_STRESS, float(vals["sigma_y"][g]))

    def apply(self, level: float) -> dict[str, np.ndarray]:
        """Write thermal strain (InitialState) and E/yield for ``level``.

        Absolute and idempotent: safe to call again for a cutback retry. The
        strategy must already be initialized (the process writes into the
        elements' constitutive laws)."""
        t0 = time.perf_counter()
        KM = self._KM
        vals = self.group_values(level)
        self._set_properties(vals)
        strain = KM.Vector(6)
        for g, smp in enumerate(self.sub_model_parts):
            e = float(vals["eps_th"][g])
            for i in range(6):
                strain[i] = e if i < 3 else 0.0
            KM.SetInitialStateProcess3D(smp, strain, self._zero6, self._eye3).ExecuteInitializeSolutionStep()
        self.apply_seconds += time.perf_counter() - t0
        self.n_apply += 1
        return vals


# --------------------------------------------------------------------------- driver

class _ConvergedHook:
    """Strategy proxy: calls ``callback()`` after each finalized (= converged,
    see ``nonlinear.solve_increment``) increment."""

    def __init__(self, strategy: Any, callback: Callable[[], None]) -> None:
        self._s = strategy
        self._cb = callback

    def __getattr__(self, name: str):
        return getattr(self._s, name)

    def FinalizeSolutionStep(self) -> None:  # noqa: N802 (Kratos API name)
        self._s.FinalizeSolutionStep()
        self._cb()


def thermal_levels(n_steps: int, increments: int | Sequence[int]) -> list[float]:
    """Absolute levels ``k + j / n_k`` (j = 1..n_k) for steps k = 0..n_steps-1;
    the last level of step k is exactly ``k + 1``."""
    incs = [int(increments)] * n_steps if isinstance(increments, (int, np.integer)) \
        else [int(n) for n in increments]
    if len(incs) != n_steps or any(n < 1 for n in incs):
        raise ValueError(f"need {n_steps} positive increment counts, got {incs}")
    levels: list[float] = []
    for k, n in enumerate(incs):
        levels.extend(float(k + 1) if j == n else k + j / n for j in range(1, n + 1))
    return levels


def run_thermal_steps(mp: Any, KM: Any, make_strategy: Callable[[], Any],
                      groups: CellThermalGroups,
                      increments: int | Sequence[int] = 10,
                      on_converged: Callable[[float], None] | None = None,
                      cutback: float = 0.5, min_frac: float = 1e-3) -> dict:
    """Run the thermal history ``groups.states[0] -> ... -> states[-1]``.

    ``on_converged(level)`` is called after every converged increment
    (including cutback sub-increments), with the element state committed.
    Returns the ``run_load_stepping`` dict plus ``levels`` (targets) and
    ``apply_seconds`` (time spent writing thermal state)."""
    if groups.n_steps < 1:
        raise ValueError("need at least two thermal states")
    levels = thermal_levels(groups.n_steps, increments)
    current = {"level": 0.0}

    def apply(level: float) -> None:
        current["level"] = float(level)
        groups.apply(level)

    def make():
        strat = make_strategy()
        if on_converged is None:
            return strat
        return _ConvergedHook(strat, lambda: on_converged(current["level"]))

    out = NL.run_load_stepping(mp, KM, make, levels, apply, cutback=cutback,
                               min_frac=min_frac, start_level=0.0)
    out["levels"] = levels
    out["apply_seconds"] = groups.apply_seconds
    return out


# --------------------------------------------------------------------------- E03 data

def _read_tsv(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter="\t", ndmin=2)


def e03_materials_initial_yield(data_dir: str | Path) -> dict[str, ThermoMaterial]:
    """E03 (glass-to-metal seal) materials from the Coreform Abaqus tables.

    Source: ``research-cf-IGA/user-guide/example-problems/glass-to-metal-seal``
    (Tables 9-16 and ``data/*.tsv``); all expansions are secant about 20 C.
    APPROXIMATION (name says it): plasticity is perfect at the INITIAL yield
    (304L: eps_p = 0 row per temperature; Alloy 52: 234.4 MPa) - the Alloy 52
    linear hardening (1924 MPa) and the 304L multilinear hardening are not
    represented (generic Kratos laws harden on plastic dissipation, S30).
    """
    d = Path(data_dir)
    glass = _read_tsv(d / "abaqus_s8061_elastic.tsv")        # E, nu, T
    ss_el = _read_tsv(d / "abaqus_ss304l_elastic.tsv")       # E, nu, T
    ss_ex = _read_tsv(d / "abaqus_ss304l_expansion.tsv")     # alpha, T
    ss_pl = _read_tsv(d / "abaqus_ss304l_plastic.tsv")       # sigma_y, eps_p, T
    a52_pl = _read_tsv(d / "abaqus_alloy52_plastic.tsv")     # sigma_y, eps_p
    for name, arr, nu in (("S-8061", glass, 0.22), ("304L", ss_el, 0.3)):
        if not np.allclose(arr[:, 1], nu):
            raise ValueError(f"{name}: Poisson ratio not constant {nu}")
    y0 = ss_pl[ss_pl[:, 1] == 0.0]
    return {
        "S8061_glass": ThermoMaterial(
            "S8061_glass", TempTable(tuple(glass[:, 2]), tuple(glass[:, 0])), 0.22,
            9.30e-6, t_zero=20.0, yield_stress=None),
        "Alloy52_pins": ThermoMaterial(
            "Alloy52_pins", 165500.0, 0.29, 9.95e-6, t_zero=20.0,
            yield_stress=float(a52_pl[0, 0])),
        "SS304L_header": ThermoMaterial(
            "SS304L_header", TempTable(tuple(ss_el[:, 2]), tuple(ss_el[:, 0])), 0.3,
            TempTable(tuple(ss_ex[:, 1]), tuple(ss_ex[:, 0])), t_zero=20.0,
            yield_stress=TempTable(tuple(y0[:, 2]), tuple(y0[:, 0]))),
    }


# --------------------------------------------------------------------------- multi-body helpers

def attach_region_elements(
    KM: Any,
    model_part: Any,
    volume: Any,
    box: Any,
    regions: Sequence[tuple[Mapping[tuple[int, int, int], np.ndarray], Iterable[tuple[int, int, int]]]],
    element_name: str = "SmallDisplacementElement3D4N",
    start_id: int = 1,
    n_gauss: int | None = None,
) -> tuple[list[Any], np.ndarray, set[tuple[int, int, int]]]:
    """Create one QP element per ``(cell, material_region)`` on a shared spline.

    Each entry in ``regions`` is ``(qp_table, out_cells)`` produced by
    :func:`wht_iga.solver.trimmed.build_trimmed_table` for one bonded solid or
    material domain. Non-OUT cells use ``qp_table[(i, j, k)]`` when cut or
    ``cell_gauss_points`` when strictly interior. A temporary placeholder
    ``Properties`` is attached and later replaced by :class:`CellThermalGroups`.

    Returns:
        ``(elements, material_index_per_element, active_cells_union)``.
    """
    from wht_iga.kratos_backend.elements import check_qp_element_name

    check_qp_element_name(element_name)
    nx, ny, nz = (int(v) for v in box.n_cells)
    pid = max([int(p.Id) for p in model_part.Properties] + [0]) + 1
    tmp_props = model_part.CreateNewProperties(pid)
    elements: list[Any] = []
    mat_indices: list[int] = []
    active_union: set[tuple[int, int, int]] = set()
    eid = int(start_id)

    for m_idx, (qp_table, out_cells) in enumerate(regions):
        out_set = {tuple(int(v) for v in c) for c in out_cells}
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    cell = (i, j, k)
                    if cell in out_set:
                        continue
                    if cell in qp_table:
                        pts = np.asarray(qp_table[cell], dtype=float).tolist()
                    else:
                        pts = cell_gauss_points(box, i, j, k, n_gauss=n_gauss)
                    if not pts:
                        continue
                    active_union.add(cell)
                    qps = KM.GeometriesVector()
                    volume.CreateQuadraturePointGeometries(qps, 2, pts)
                    if len(qps) != 1:
                        raise RuntimeError(
                            f"region {m_idx} cell {cell}: expected 1 QP geometry, got {len(qps)}"
                        )
                    elements.append(
                        model_part.CreateNewElement(element_name, eid, qps[0], tmp_props)
                    )
                    mat_indices.append(m_idx)
                    eid += 1

    return elements, np.asarray(mat_indices, dtype=np.int64), active_union


def fix_isostatic_321(
    KM: Any,
    model_part: Any,
    box: Any,
    inactive_mask: np.ndarray | None = None,
    axis: int = 0,
    side: int = 0,
) -> tuple[int, int, int]:
    """Apply a minimal 6-DOF 3-2-1 kinematic restraint on active face control points.

    Selects three non-collinear active control points ``(p0, p1, p2)`` on the
    boundary face ``(axis, side)`` where ``p1 - p0`` aligns with the first
    in-plane axis ``t1`` and ``p2 - p0`` spans the second in-plane axis ``t2``:
    ``p0`` fixes ``(u_axis, u_t1, u_t2)``, ``p1`` fixes ``(u_axis, u_t2)``, and
    ``p2`` fixes ``(u_axis)``. This removes all 6 rigid-body modes while
    allowing stress-free uniform thermal expansion.
    """
    comps = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    t1, t2 = (d for d in range(3) if d != axis)
    face_ids = [int(I) for I in box.boundary_control_points(axis, side)]
    if inactive_mask is not None:
        face_ids = [I for I in face_ids if not bool(inactive_mask[I])]
    if len(face_ids) < 3:
        raise ValueError(f"need at least 3 active control points on face ({axis}, {side}), got {len(face_ids)}")

    coords = box.control_points()
    center = coords[face_ids].mean(axis=0)
    p0 = min(face_ids, key=lambda I: float(np.sum((coords[I] - center) ** 2)))
    p1_cands = [I for I in face_ids if I != p0 and abs(coords[I, t2] - coords[p0, t2]) < 1e-9]
    if not p1_cands:
        raise ValueError(f"no collinear partner along axis {t1} on face ({axis}, {side})")
    p1 = max(p1_cands, key=lambda I: abs(float(coords[I, t1] - coords[p0, t1])))
    p2_cands = [I for I in face_ids if abs(coords[I, t2] - coords[p0, t2]) > 1e-6]
    if not p2_cands:
        raise ValueError(f"no out-of-line partner along axis {t2} on face ({axis}, {side})")
    v1 = coords[p1] - coords[p0]
    p2 = max(p2_cands, key=lambda I: abs(float(np.cross(v1, coords[I] - coords[p0])[axis])))

    for c in (axis, t1, t2):
        model_part.GetNode(p0 + 1).Fix(comps[c])
    for c in (axis, t2):
        model_part.GetNode(p1 + 1).Fix(comps[c])
    model_part.GetNode(p2 + 1).Fix(comps[axis])
    return p0, p1, p2
