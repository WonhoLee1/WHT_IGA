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
* :class:`HardeningTable` - isotropic hardening sigma_y(PEEQ, T) in the
  Abaqus ``*PLASTIC`` form (S42).
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

Hardening (S42). Without a :class:`HardeningTable` the plastic laws are
perfectly plastic (``HARDENING_CURVE = 3``). With one, the group Properties
use ``HARDENING_CURVE = 6`` (CurveDefinedByPoints) with
``PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE`` = the PEEQ grid plus one far
point that continues the last stress flat (Abaqus holds it constant; Kratos
would soften linearly after the last point). The generic laws harden on the
normalised plastic dissipation ``kappa = W_p * l_c / FRACTURE_ENERGY`` with
``W_p = sum sigma : d eps_p`` and ``l_c`` = largest distance of the element's
control points from their centroid (reference configuration). On every
segment Kratos uses ``S = a sqrt(1 - b kappa)``, i.e. ``S^2`` linear in
``W_p``, which is exactly ``S`` linear in PEEQ with ``W_p = int S dPEEQ``;
``l_c`` and ``FRACTURE_ENERGY`` cancel (no mesh-size dependence). Two
discrete defects remain and are removed by a PEEQ sync: Kratos integrates
``W_p`` with the end-of-iteration stress (excess about ``H dPEEQ^2 / 2`` per
increment, accumulating), and a temperature-dependent curve would be read at
the plastic WORK of the history instead of its PEEQ. So
:meth:`CellThermalGroups.commit` accumulates PEEQ from the converged
``PLASTIC_STRAIN_VECTOR`` increments (``sqrt(2/3 |d eps_p|^2)``) and
:meth:`CellThermalGroups.apply` writes ``kappa = l_c W_T(PEEQ) / G_f`` at
every integration point before each increment (``SetValuesOnIntegrationPoints
(PLASTIC_DISSIPATION)``), so the threshold is sigma_y(PEEQ, T) at the start
of every increment (Abaqus semantics). Residual error: Kratos' return
tolerance ``F <= 1e-4 sigma_y`` plus the within-increment excess, about
``H^2 dPEEQ^2 / sigma_y`` in stress.

Limits (S30, S42): the temperature is constant per element (the process
writes the same state to every integration point of an element); hardening is
isotropic only; groups and Properties must be built before the strategy is
initialized (the elements clone ``CONSTITUTIVE_LAW`` from their Properties in
``Initialize``); the PEEQ sync is a Python loop over the elements of
hardening groups (two element calls per element and increment).
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend.geometry import cell_gauss_points

__all__ = [
    "TempTable",
    "HardeningTable",
    "ThermoMaterial",
    "CellThermalGroups",
    "characteristic_length",
    "default_law_factory",
    "thermal_levels",
    "run_thermal_steps",
    "e03_materials_initial_yield",
    "e03_materials",
    "attach_region_elements",
    "fix_isostatic_321",
]

# Length of the flat segment appended after the last PEEQ point of a hardening
# table (Kratos softens after its last point; Abaqus holds the stress).
FLAT_EXTENSION = 10.0


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


def _curve_work(grid: np.ndarray, stress: np.ndarray, eps_p) -> np.ndarray:
    """``int_0^eps_p S de`` for S piecewise linear on ``grid``, flat beyond the end."""
    e = np.asarray(eps_p, dtype=float)
    if np.any(e < 0.0):
        raise ValueError("equivalent plastic strain must be >= 0")
    seg = np.concatenate([[0.0], np.cumsum(0.5 * (stress[1:] + stress[:-1]) * np.diff(grid))])
    slope = np.append(np.diff(stress) / np.diff(grid), 0.0)
    j = np.clip(np.searchsorted(grid, e, side="right") - 1, 0, grid.size - 1)
    de = e - grid[j]
    return seg[j] + stress[j] * de + 0.5 * slope[j] * de * de


@dataclass(frozen=True)
class HardeningTable:
    """Isotropic hardening sigma_y(PEEQ, T) - Abaqus ``*PLASTIC`` with temperature.

    ``stresses[k][j]`` is the yield stress at ``temperatures[k]`` and
    equivalent plastic strain ``plastic_strains[j]`` (common grid, first point
    0). Linear in PEEQ between points and constant after the last one, linear
    in T between rows and constant outside (Abaqus, equal strain grids).
    Curves must not soften (Kratos' point curve would, but S42 verified only
    hardening).
    """

    plastic_strains: tuple[float, ...]
    stresses: tuple[tuple[float, ...], ...]
    temperatures: tuple[float, ...] = (0.0,)

    def __post_init__(self) -> None:
        e = np.asarray(self.plastic_strains, dtype=float)
        s = np.asarray(self.stresses, dtype=float)
        t = np.asarray(self.temperatures, dtype=float)
        if e.ndim != 1 or e.size < 2 or e[0] != 0.0 or np.any(np.diff(e) <= 0.0):
            raise ValueError("plastic_strains: >= 2 strictly increasing values starting at 0")
        if t.ndim != 1 or t.size < 1 or np.any(np.diff(t) <= 0.0):
            raise ValueError("temperatures must be strictly increasing")
        if s.shape != (t.size, e.size):
            raise ValueError(f"stresses must be {t.size} rows x {e.size} values, got {s.shape}")
        if np.any(s <= 0.0) or np.any(np.diff(s, axis=1) < 0.0):
            raise ValueError("stresses must be positive and non-decreasing in PEEQ")
        object.__setattr__(self, "plastic_strains", tuple(float(x) for x in e))
        object.__setattr__(self, "stresses", tuple(tuple(float(x) for x in row) for row in s))
        object.__setattr__(self, "temperatures", tuple(float(x) for x in t))

    def curve(self, temperature: float) -> np.ndarray:
        """Yield stress at every grid point for one temperature."""
        s = np.asarray(self.stresses)
        return np.array([np.interp(float(temperature), self.temperatures, s[:, j])
                         for j in range(s.shape[1])])

    def initial_yield(self, temperature):
        return np.interp(temperature, self.temperatures, np.asarray(self.stresses)[:, 0])

    def stress(self, eps_p, temperature: float):
        """sigma_y(PEEQ, T); constant after the last PEEQ point."""
        return np.interp(eps_p, self.plastic_strains, self.curve(temperature))

    def work(self, eps_p, temperature: float):
        """Plastic work ``int_0^PEEQ sigma_y(e, T) de`` per unit volume at fixed T."""
        return _curve_work(np.asarray(self.plastic_strains), self.curve(temperature), eps_p)


@dataclass(frozen=True)
class ThermoMaterial:
    """Temperature-dependent isotropic material (E03 data forms).

    ``expansion`` is the secant (mean) coefficient about ``t_zero`` (Abaqus
    ``ZERO`` temperature, default 0 as in Abaqus). Plasticity: either
    ``yield_stress`` (perfect plasticity) or ``hardening`` (isotropic
    hardening table, S42); neither means thermoelastic (linear elastic law,
    e.g. the S-8061 glass).
    """

    name: str
    young_modulus: float | TempTable
    poisson_ratio: float
    expansion: float | TempTable
    t_zero: float = 0.0
    yield_stress: float | TempTable | None = None
    hardening: HardeningTable | None = None

    def __post_init__(self) -> None:
        if self.yield_stress is not None and self.hardening is not None:
            raise ValueError(f"material {self.name!r}: give yield_stress OR hardening, not both")

    @property
    def plastic(self) -> bool:
        return self.yield_stress is not None or self.hardening is not None

    def young(self, temperature):
        return _eval(self.young_modulus, temperature)

    def yield_at(self, temperature):
        """Initial yield stress (PEEQ = 0) at ``temperature``."""
        if self.hardening is not None:
            return self.hardening.initial_yield(temperature)
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


def characteristic_length(element: Any) -> float:
    """Kratos ``CalculateCharacteristicLengthOnReferenceConfiguration`` for an
    element: largest distance of its geometry points (for IGA QP elements: the
    control points with support on the cell) from their centroid, reference
    coordinates. p = 1 unit cell: sqrt(3)/2."""
    pts = np.array([[p.X0, p.Y0, p.Z0] for p in element.GetGeometry()], dtype=float)
    return float(np.sqrt(np.max(np.sum((pts - pts.mean(axis=0)) ** 2, axis=1))))


def _kvec(KM: Any, values) -> Any:
    v = KM.Vector(len(values))
    for i, x in enumerate(values):
        v[i] = float(x)
    return v


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
        quantum: temperature rounding for grouping. A rounding bin is one
            quantum wide and the group temperature is the member MEAN (not
            the bin centre), so the per-cell temperature error is below one
            quantum, not quantum/2.
        fracture_energy: ``FRACTURE_ENERGY`` for the plastic laws (required by
            Kratos; irrelevant for perfect plasticity, cancels for hardening
            tables but must exceed ``l_c`` x the work of the extended curve).

    Hardening groups (material with a :class:`HardeningTable`) additionally
    track PEEQ per integration point: :meth:`commit` after every converged
    increment, :meth:`apply` writes the matching plastic dissipation.
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
        self.fracture_energy = float(fracture_energy)
        self.sub_model_parts: list[Any] = []
        self.properties: list[Any] = []
        self._members: list[list[Any]] = []
        for g in range(self.n_groups):
            mat = self.materials[self.group_material[g]]
            props = mp.CreateNewProperties(pid0 + g)
            props.SetValue(KM.POISSON_RATIO, float(mat.poisson_ratio))
            props.SetValue(KM.CONSTITUTIVE_LAW, laws[self.group_material[g]])
            if mat.plastic:
                props.SetValue(KM.FRACTURE_ENERGY, self.fracture_energy)
                if mat.hardening is None:
                    props.SetValue(self._cl.HARDENING_CURVE, 3)  # perfect plasticity
                else:
                    props.SetValue(self._cl.HARDENING_CURVE, 6)  # CurveDefinedByPoints
                    grid = mat.hardening.plastic_strains
                    props.SetValue(self._cl.PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE,
                                   _kvec(KM, grid + (grid[-1] + FLAT_EXTENSION,)))
            members = [elements[e] for e in order[bounds[g]:bounds[g + 1]]]
            for el in members:
                el.Properties = props
            smp = mp.CreateSubModelPart(f"{name_prefix}_{g}")
            smp.AddElements([int(el.Id) for el in members])
            self.sub_model_parts.append(smp)
            self.properties.append(props)
            self._members.append(members)
        self.hardening_groups = [g for g in range(self.n_groups)
                                 if self.materials[self.group_material[g]].hardening is not None]
        self._hard: dict[int, dict[str, Any]] | None = None  # PEEQ state, built on first use
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

    def _curve(self, g: int, temperature: float) -> np.ndarray:
        return self.materials[self.group_material[g]].hardening.curve(temperature)

    def _set_properties(self, vals: dict[str, np.ndarray]) -> None:
        KM = self._KM
        for g, props in enumerate(self.properties):
            props.SetValue(KM.YOUNG_MODULUS, float(vals["E"][g]))
            if not math.isnan(vals["sigma_y"][g]):
                props.SetValue(KM.YIELD_STRESS, float(vals["sigma_y"][g]))
        for g in self.hardening_groups:
            c = self._curve(g, vals["T"][g])  # c[0] == sigma_y (Kratos Check)
            self.properties[g].SetValue(self._cl.EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE,
                                        _kvec(KM, np.append(c, c[-1])))

    def _plastic_strain(self, el: Any) -> np.ndarray:
        return np.array([list(v) for v in el.CalculateOnIntegrationPoints(
            self._cl.PLASTIC_STRAIN_VECTOR, self.mp.ProcessInfo)], dtype=float).reshape(-1, 6)

    def _hardening_state(self) -> dict[int, dict[str, Any]]:
        if self._hard is None:  # laws must exist (strategy initialized)
            self._hard = {}
            for g in self.hardening_groups:
                els = self._members[g]
                eps_p = [self._plastic_strain(el) for el in els]
                self._hard[g] = {"lc": [characteristic_length(el) for el in els],
                                 "eps_p": eps_p, "peeq": [np.zeros(len(e)) for e in eps_p]}
        return self._hard

    def commit(self) -> None:
        """Accumulate PEEQ of hardening groups after a CONVERGED increment:
        ``dPEEQ = sqrt(2/3 |d eps_p|^2)`` (Voigt, engineering shear)."""
        for g, st in self._hardening_state().items():
            for i, el in enumerate(self._members[g]):
                ep = self._plastic_strain(el)
                d = ep - st["eps_p"][i]
                st["peeq"][i] = st["peeq"][i] + np.sqrt(
                    (2.0 / 3.0) * (np.sum(d[:, :3] ** 2, axis=1) + 0.5 * np.sum(d[:, 3:] ** 2, axis=1)))
                st["eps_p"][i] = ep

    def equivalent_plastic_strain(self) -> dict[int, np.ndarray]:
        """Tracked PEEQ per integration point, ``{element id: array}`` (hardening groups)."""
        return {int(el.Id): st["peeq"][i].copy()
                for g, st in self._hardening_state().items()
                for i, el in enumerate(self._members[g])}

    def _write_dissipation(self, vals: dict[str, np.ndarray]) -> None:
        """kappa = l_c W_T(PEEQ) / G_f at every integration point (threshold = sigma_y(PEEQ, T))."""
        KM = self._KM
        for g, st in self._hardening_state().items():
            grid = np.asarray(self.materials[self.group_material[g]].hardening.plastic_strains)
            c = self._curve(g, vals["T"][g])
            for i, el in enumerate(self._members[g]):
                peeq = st["peeq"][i]
                if peeq.size and peeq.max() >= grid[-1] + FLAT_EXTENSION:
                    raise RuntimeError(f"PEEQ {peeq.max():.4g} beyond the flat extension of group {g}")
                kappa = st["lc"][i] / self.fracture_energy * _curve_work(grid, c, peeq)
                el.SetValuesOnIntegrationPoints(KM.PLASTIC_DISSIPATION, [float(k) for k in kappa],
                                                self.mp.ProcessInfo)

    def apply(self, level: float) -> dict[str, np.ndarray]:
        """Write thermal strain (InitialState), E/yield/hardening curve and, for
        hardening groups, the plastic dissipation of the tracked PEEQ.

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
        self._write_dissipation(vals)
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
                      cutback: float = 0.5, min_frac: float = 1e-3,
                      mechanical_load: Callable[[float], None] | None = None) -> dict:
    """Run the thermal history ``groups.states[0] -> ... -> states[-1]``.

    ``mechanical_load(level)`` (optional) sets absolute mechanical loads /
    prescribed displacements for ``level`` after the thermal state is written.
    After every converged increment (including cutback sub-increments) the
    hardening PEEQ is committed (:meth:`CellThermalGroups.commit`) and then
    ``on_converged(level)`` is called, with the element state committed.
    Returns the ``run_load_stepping`` dict plus ``levels`` (targets) and
    ``apply_seconds`` (time spent writing thermal / hardening state)."""
    if groups.n_steps < 1:
        raise ValueError("need at least two thermal states")
    levels = thermal_levels(groups.n_steps, increments)
    current = {"level": 0.0}

    def apply(level: float) -> None:
        current["level"] = float(level)
        groups.apply(level)
        if mechanical_load is not None:
            mechanical_load(level)

    def converged() -> None:
        groups.commit()
        if on_converged is not None:
            on_converged(current["level"])

    def make():
        return _ConvergedHook(make_strategy(), converged)

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
    (304L: eps_p = 0 row per temperature; Alloy 52: 234.4 MPa). The tables
    with hardening are :func:`e03_materials` (S42).
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


def e03_materials(data_dir: str | Path) -> dict[str, ThermoMaterial]:
    """E03 materials with the Coreform isotropic hardening tables (S42).

    Alloy 52: ``abaqus_alloy52_plastic.tsv`` (234.4 MPa at PEEQ 0, 580.72 at
    0.18: linear 1924 MPa, constant afterwards), temperature independent.
    304L: ``abaqus_ss304l_plastic.tsv``, 12 PEEQ points (0 .. 0.023, common to
    all rows) at 7 temperatures (-55 .. 500 C), linear in T between rows.
    Elastic / expansion data as :func:`e03_materials_initial_yield`.
    """
    base = e03_materials_initial_yield(data_dir)
    d = Path(data_dir)
    ss_pl = _read_tsv(d / "abaqus_ss304l_plastic.tsv")       # sigma_y, eps_p, T
    a52_pl = _read_tsv(d / "abaqus_alloy52_plastic.tsv")     # sigma_y, eps_p
    temps = np.unique(ss_pl[:, 2])
    rows = [ss_pl[ss_pl[:, 2] == t] for t in temps]
    grid = rows[0][:, 1]
    if any(r.shape != rows[0].shape or not np.array_equal(r[:, 1], grid) for r in rows):
        raise ValueError("304L: PEEQ grids differ between temperatures")
    ss_h = HardeningTable(tuple(grid), tuple(tuple(r[:, 0]) for r in rows), tuple(temps))
    a52_h = HardeningTable(tuple(a52_pl[:, 1]), (tuple(a52_pl[:, 0]),))
    return {
        "S8061_glass": base["S8061_glass"],
        "Alloy52_pins": replace(base["Alloy52_pins"], yield_stress=None, hardening=a52_h),
        "SS304L_header": replace(base["SS304L_header"], yield_stress=None, hardening=ss_h),
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
