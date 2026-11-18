# -*- coding: utf-8 -*-
"""S42: isotropic hardening on the S28/S30 InitialState path (``thermal_plasticity``).

Mapping (module doc, S42): ``HARDENING_CURVE = 6`` (CurveDefinedByPoints) with
``PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE`` = PEEQ grid (+ a flat far
point). Kratos hardens on ``kappa = W_p l_c / G_f``; on each segment
``S = a sqrt(1 - b kappa)``, i.e. S^2 linear in W_p = int S dPEEQ, which is
exactly S linear in PEEQ, so ``l_c`` and ``G_f`` cancel. ``CellThermalGroups``
tracks PEEQ from the converged plastic strain increments and writes the
matching ``PLASTIC_DISSIPATION`` before every increment (Abaqus semantics
sigma_y(PEEQ, T) for temperature-dependent tables).

A-priori stress bound used by every hardening test (S42):

    |sigma - sigma_exact| <= 1e-4 S + 2 H^2 dPEEQ^2 / S_min  [+ 2 (h2 - h1)+ dPEEQ]

* 1e-4 S: Kratos' return tolerance (``F <= 1e-4 threshold``; trials in that
  band are also treated as elastic).
* Kratos integrates ``W_p = sum q dlambda`` over its return iterations with
  the end-of-iteration equivalent stress q (decreasing from q1), so within
  one increment ``W_p`` exceeds ``int S dPEEQ`` by at most
  ``(q1 - S_n) dPEEQ``, and ``q1 - S_n = A3 dlambda_1 <= H (q0 / S_n) dPEEQ``
  (q0 trial, S_n threshold at the start, A3 = H q0 / S_n Kratos' hardening
  term). dS/dW = H / S turns it into stress. Factor 2: q0 / S_n <= 1.2 here
  and Kratos' increment vs the analytic one.
* Only when an increment crosses a CONVEX kink (slope rises h1 -> h2):
  + 2 (h2 - h1) dPEEQ, because Kratos accepts a negative F after an
  over-return (test (i)).
* With the PEEQ sync the excess does not accumulate: the bound of the LAST
  plastic increment applies (monotone loading: per level; cyclic: the
  maximum over the history). dPEEQ, H (largest slope touched), S from the
  analytic solution.

Models: bars of cells [0, n L] x [0, L]^2, one element per cell, u_x = 0 at
x = 0 and prescribed (or 0) at x = n L. Lateral FREE (symmetry planes y = 0,
z = 0) gives uniaxial STRESS, lateral FIXED uniaxial STRAIN (S30).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses import nonlinear as NL
from wht_iga.analyses import thermal_plasticity as TP
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_DATA = (Path(__file__).resolve().parents[2] / "research-cf-IGA" / "user-guide" / "example-problems"
         / "glass-to-metal-seal" / "data")


# --------------------------------------------------------------------------- helpers

def _bar(KM, nx, lateral_fixed, L=1.0, degree=1, element_name="SmallDisplacementElement3D4N"):
    model = KM.Model()
    mp = model.CreateModelPart("S42")
    for v in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(v)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L * nx, L, L), n_cells=(nx, 1, 1), degree=degree)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = attach_cell_elements(KM, mp, vol, box, mp.CreateNewProperties(1),
                               element_name=element_name)
    mp.CloneTimeStep(0.0)
    right = []
    for n in mp.Nodes:
        if abs(n.X0) < 1e-12 or abs(n.X0 - L * nx) < 1e-12:
            n.Fix(KM.DISPLACEMENT_X)
        if abs(n.X0 - L * nx) < 1e-12:
            right.append(n)
        if lateral_fixed or abs(n.Y0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Y)
        if lateral_fixed or abs(n.Z0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Z)
    return model, mp, els, right


def _ip(KM, mp, el, var):
    return np.array([list(v) for v in el.CalculateOnIntegrationPoints(var, mp.ProcessInfo)])


def _ux_at(KM, mp, x):
    vals = [n.GetSolutionStepValue(KM.DISPLACEMENT_X) for n in mp.Nodes if abs(n.X0 - x) < 1e-12]
    assert vals and np.ptp(vals) < 1e-14
    return vals[0]


def _run(KM, mods, mp, els, materials, mat_idx, states, increments, *, load=None,
         var=None, top=1.0, finite=False):
    """CellThermalGroups + run_thermal_steps; record per converged level."""
    var = KM.CAUCHY_STRESS_VECTOR if var is None else var
    groups = TP.CellThermalGroups(KM, mods, mp, els, materials, mat_idx, states,
                                  law_factory=TP.default_law_factory(mods, finite_strain=finite))
    xs = sorted({round(n.X0, 12) for n in mp.Nodes})
    rec: list[dict] = []

    def obs(level):
        pq = groups.equivalent_plastic_strain()
        rec.append({
            "level": level,
            "stress": [_ip(KM, mp, el, var) for el in els],
            "peeq": [pq.get(int(el.Id)) for el in els],
            "u": [_ux_at(KM, mp, x) for x in xs],
            "uy_top": [n.GetSolutionStepValue(KM.DISPLACEMENT_Y) for n in mp.Nodes
                       if abs(n.Y0 - top) < 1e-12],
        })

    out = TP.run_thermal_steps(mp, KM, lambda: NL.make_newton_strategy(KM, mp), groups,
                               increments, on_converged=obs, mechanical_load=load)
    assert out["n_solves"] == len(out["levels"]), "unexpected cutback (levels would not align)"
    return groups, rec


def _pull(KM, nodes, disp):
    """mechanical_load: u_x = disp(level) on ``nodes``."""
    def load(level):
        for n in nodes:
            n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(disp(level)))
    return load


def _bound(dpeeq, H, s_min, s_max, dh_up=0.0):
    """A-priori stress bound (module docstring); ``dh_up``: slope increase at a
    kink crossed by the increment (over-return term)."""
    return 1e-4 * s_max + 2.0 * H ** 2 * dpeeq ** 2 / s_min + 2.0 * dh_up * dpeeq


def _peeq_uniaxial_stress(eps, E, grid, S):
    """Monotone uniaxial stress: PEEQ with E (eps - PEEQ) = S(PEEQ), S piecewise
    linear on ``grid``, flat after the last point (f decreasing -> one root)."""
    if E * eps <= S[0]:
        return 0.0
    for j in range(len(grid)):
        last = j + 1 == len(grid)
        h = 0.0 if last else (S[j + 1] - S[j]) / (grid[j + 1] - grid[j])
        e = (E * eps - S[j] + h * grid[j]) / (E + h)
        if grid[j] - 1e-15 <= e <= (np.inf if last else grid[j + 1] + 1e-15):
            return e
    raise RuntimeError("no segment")


def _slopes_touched(grid, S, a, b):
    """(largest slope, largest slope INCREASE between consecutive segments) on
    [a, b]; slope 0 after the last point."""
    h = np.append(np.diff(S) / np.diff(grid), 0.0)
    lo = max(np.searchsorted(grid, a, side="right") - 1, 0)
    hi = max(np.searchsorted(grid, b, side="right") - 1, 0)
    hs = h[lo:hi + 1]
    return float(hs.max()), float(max([0.0] + list(np.diff(hs))))


def _uniaxial_stress_thermal_reference(T_path, E, alpha, S0, H):
    """Incremental 1D return (backward Euler in PEEQ, = Abaqus), clamped cell,
    uniaxial stress: trial sigma = -E (alpha T + p); yield |sigma| = S0(T) +
    H(T) PEEQ with linear hardening at each T. Returns per level dict."""
    p = eb = 0.0
    out = []
    for T in T_path:
        trial = -E * (alpha * T + p)
        f = abs(trial) - (S0(T) + H(T) * eb)
        d = max(f, 0.0) / (E + H(T))
        c = float(np.sign(trial))
        p += c * d
        eb += d
        out.append({"T": T, "sigma": trial - c * E * d, "p": p, "peeq": eb, "dpeeq": d,
                    "H": H(T), "S0": S0(T)})
    return out


def _check_uniaxial_stress_levels(rec, ref, tols, E, nu, L=1.0, alpha=0.0):
    for r, q, tol in zip(rec, ref, tols):
        for st in r["stress"]:
            assert st[:, 0] == pytest.approx(q["sigma"], abs=tol), r["level"]
            assert np.abs(st[:, 1:]).max() < tol, r["level"]
        for pq in r["peeq"]:
            assert pq == pytest.approx(q["peeq"], abs=tol / E), r["level"]
        eps_y = alpha * q["T"] - nu * q["sigma"] / E - q["p"] / 2.0
        assert r["uy_top"] == pytest.approx([L * eps_y] * len(r["uy_top"]), abs=L * tol / E)


# --------------------------------------------------------------------------- (a)

_SQ3 = np.sqrt(3.0) / 2.0


@pytest.mark.parametrize("L, degree, nx, lc", [
    (1.0, 1, 1, [_SQ3]),
    (0.1, 1, 1, [0.1 * _SQ3]),
    (1.0, 2, 3, [np.sqrt(0.5 + 25.0 / 36.0), np.sqrt(1.5), np.sqrt(0.5 + 25.0 / 36.0)]),
], ids=["p1_unit_cell", "p1_cell_0.1", "p2_three_cells"])
def test_a_linear_hardening_uniaxial_stress_mesh_independent(L, degree, nx, lc):
    """Uniaxial stress, E = 100, nu = 0.3, sigma_y = 1 + 10 PEEQ (table (0, 1),
    (1, 11)), eps = u/l: 0 -> 0.052 in 20 increments.

    1D: PEEQ = (E eps - S0) / (E + H) after eps = S0/E = 0.01, sigma = S0 + H
    PEEQ (tangent E H/(E + H) = 9.09), eps_yy = -nu sigma/E - PEEQ/2; at the
    end PEEQ = 4.2/110 = 0.0381818, sigma = 1.381818.
    l_c (Kratos formula, largest control-point distance from the centroid):
    p = 1 cell of size L: L sqrt(3)/2; p = 2, 3 cells on [0, 3]: Greville x
    points (0, .5, 1.5), (.5, 1.5, 2.5), (1.5, 2.5, 3), y, z (0, .5, 1) ->
    sqrt(1/2 + 25/36), sqrt(3/2). l_c spans 0.087 .. 1.22: the response must
    not depend on it (and a wrong l_c in the PEEQ sync would shift S)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, right = _bar(KM, nx, lateral_fixed=False, L=L, degree=degree)
    assert [TP.characteristic_length(el) for el in els] == pytest.approx(lc, abs=1e-12)
    E, nu, S0, H, eps_max = 100.0, 0.3, 1.0, 10.0, 0.052
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable((0.0, 1.0), ((S0, S0 + H),)))
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [0.0, 0.0], 20,
                   load=_pull(KM, right, lambda s: s * eps_max * L * nx), top=L)
    ref, tols, prev = [], [], 0.0
    for r in rec:
        e = r["level"] * eps_max
        eb = _peeq_uniaxial_stress(e, E, [0.0, 1.0], [S0, S0 + H])
        ref.append({"T": 0.0, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        tols.append(_bound(eb - prev, H, S0, E * (e - eb)))
        prev = eb
    assert ref[-1]["peeq"] == pytest.approx(4.2 / 110.0, abs=1e-14)
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 42.0 / 110.0, abs=1e-12)
    assert max(tols) < 1.5e-3
    _check_uniaxial_stress_levels(rec, ref, tols, E, nu, L=L)


# --------------------------------------------------------------------------- (b)

def test_b_linear_hardening_uniaxial_strain():
    """Uniaxial strain (lateral fixed), E = 100, nu = 0.3, sigma_y = 1 + 10
    PEEQ, eps: 0 -> 0.041 in 20 increments.

    lambda = E nu / ((1 + nu)(1 - 2 nu)) = 57.692, G = E / (2 (1 + nu)) =
    38.462. Isochoric plastic strain (p, -p/2, -p/2), PEEQ = p:
        sigma_x = (lambda + 2G) eps - 2G p,  sigma_y = sigma_z = lambda eps + G p,
        q = sigma_x - sigma_y = 2G eps - 3G p = S0 + H p
        -> p = (2G eps - S0) / (3G + H) after eps = S0 / (2G) = 0.013,
    tangent d sigma_x / d eps = K + 4 G H / (3 (3G + H)) = 83.333 + 4.090.
    Errors from the bound r: sigma_x <= 2/3 r, sigma_y <= r/3, PEEQ <= r/(3G + H)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, right = _bar(KM, 1, lateral_fixed=True)
    E, nu, S0, H, eps_max = 100.0, 0.3, 1.0, 10.0, 0.041
    lam, G = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable((0.0, 1.0), ((S0, S0 + H),)))
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [0.0, 0.0], 20,
                   load=_pull(KM, right, lambda s: s * eps_max))
    prev = 0.0
    hist = []
    for r in rec:
        e = r["level"] * eps_max
        p = max(0.0, (2 * G * e - S0) / (3 * G + H))
        sx, sy = (lam + 2 * G) * e - 2 * G * p, lam * e + G * p
        tol = _bound(p - prev, H, S0, sx - sy)
        prev = p
        st = r["stress"][0]
        assert st[:, 0] == pytest.approx(sx, abs=tol), r["level"]
        assert st[:, 1] == pytest.approx(sy, abs=tol) and st[:, 2] == pytest.approx(sy, abs=tol)
        assert np.abs(st[:, 3:]).max() < tol
        assert r["peeq"][0] == pytest.approx(p, abs=tol / (3 * G + H))
        hist.append((e, sx, float(st[:, 0].mean()), tol))
    (e1, s1, k1, t1), (e2, s2, k2, t2) = hist[-2], hist[-1]
    tangent = E / 3 / (1 - 2 * nu) + 4 * G * H / (3 * (3 * G + H))
    assert (s2 - s1) / (e2 - e1) == pytest.approx(tangent, rel=1e-12)   # closed form = 1D answer
    assert (k2 - k1) / (e2 - e1) == pytest.approx(tangent, abs=(t1 + t2) / (e2 - e1))
    assert np.all(rec[-1]["peeq"][0] > 0.01)


# --------------------------------------------------------------------------- (c)

_GRID_C = (0.0, 0.005, 0.02, 0.04)
_S_C = (1.0, 1.4, 1.7, 1.8)


def test_c_multilinear_curve_uniaxial_stress():
    """Uniaxial stress, E = 100, nu = 0.3, sigma_y points (PEEQ, S) = (0, 1.0),
    (0.005, 1.4), (0.02, 1.7), (0.04, 1.8) (slopes 80, 20, 5) and constant
    after 0.04 (Abaqus); eps: 0 -> 0.07 in 70 increments.

    Per segment j: PEEQ = (E eps - S_j + H_j e_j) / (E + H_j); segment ends
    at eps = e_j + S_j / E: 0.01, 0.019, 0.037, 0.058; at eps = 0.07 PEEQ =
    0.052 in the flat extension, sigma = 1.8. Per-level bound with the
    largest slope touched by the increment (80 near yield)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, right = _bar(KM, 1, lateral_fixed=False)
    E, nu, eps_max = 100.0, 0.3, 0.07
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable(_GRID_C, (_S_C,)))
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [0.0, 0.0], 70,
                   load=_pull(KM, right, lambda s: s * eps_max))
    ref, tols, prev = [], [], 0.0
    for r in rec:
        e = r["level"] * eps_max
        eb = _peeq_uniaxial_stress(e, E, _GRID_C, _S_C)
        ref.append({"T": 0.0, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        h_max, dh_up = _slopes_touched(_GRID_C, _S_C, prev, eb)
        assert dh_up == 0.0  # concave curve: no over-return term
        tols.append(_bound(eb - prev, h_max, _S_C[0], E * (e - eb)))
        prev = eb
    ends = [g + s / E for g, s in zip(_GRID_C, _S_C)]
    assert ends == pytest.approx([0.01, 0.019, 0.037, 0.058], abs=1e-15)
    assert ref[-1]["peeq"] == pytest.approx(0.052, abs=1e-14) and ref[-1]["sigma"] == pytest.approx(1.8)
    # every segment and the flat extension is visited
    hit = {int(np.searchsorted(_GRID_C, q["peeq"], side="right")) for q in ref if q["peeq"] > 0}
    assert hit == {1, 2, 3, 4}
    _check_uniaxial_stress_levels(rec, ref, tols, E, nu)


# --------------------------------------------------------------------------- (d)

def _heat_hold_cool(KM, mods, mat, element_name="SmallDisplacementElement3D4N", var=None,
                    finite=False):
    model, mp, els, _right = _bar(KM, 1, lateral_fixed=False, element_name=element_name)
    return _run(KM, mods, mp, els, [mat], 0, [0.0, 40.0, 40.0, 0.0], [16, 2, 16],
                var=var, finite=finite)


_T_HHC = [2.5 * k for k in range(1, 17)] + [40.0, 40.0] + [40.0 - 2.5 * k for k in range(1, 17)]


def test_d_thermal_heat_hold_cool_linear_hardening():
    """One clamped cell, lateral free (uniaxial stress), E = 100, nu = 0.3,
    alpha = 1e-3, sigma_y = 1 + 10 PEEQ, T: 0 -> 40 (16 inc) -> 40 (hold, 2)
    -> 0 (16 inc); eps_xx = 0, sigma = -E (alpha T + p).

    Heat: elastic to T = 10, then compressive flow p = -PEEQ,
        PEEQ = (E alpha T - S0) / (E + H): T = 40 -> 3/110, sigma = -1.272727.
    Hold: no change. Cool: elastic sigma = E (PEEQ_h - alpha T) until it
    reaches +S(PEEQ_h) = 1.272727 at T = 14.545; then reverse flow
    (isotropic: PEEQ keeps growing while p decreases),
        dP = (E (PEEQ_h - alpha T) - S0 - H PEEQ_h) / (E + H):
    T = 0 -> dP = 1.454545/110 = 0.0132231, PEEQ = 0.0404959,
    residual sigma = 1 + 10 PEEQ = +1.404959 (perfect plasticity: +1.0).
    Lateral u_y(top) = alpha T - nu sigma/E - p/2."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, S0, H = 100.0, 1.0, 10.0
    mat = TP.ThermoMaterial("m", E, 0.3, 1e-3, hardening=TP.HardeningTable((0.0, 1.0), ((S0, S0 + H),)))
    _g, rec = _heat_hold_cool(KM, mods, mat)
    ref = _uniaxial_stress_thermal_reference(_T_HHC, E, 1e-3, lambda T: S0, lambda T: H)
    T_rec = [40.0 * min(r["level"], 1.0) if r["level"] <= 2.0 else 40.0 * (3.0 - r["level"])
             for r in rec]
    assert T_rec == pytest.approx(_T_HHC, abs=1e-12)
    assert ref[15]["sigma"] == pytest.approx(-(1.0 + 30.0 / 110.0), abs=1e-12)
    assert ref[15]["peeq"] == pytest.approx(3.0 / 110.0, abs=1e-14)
    assert ref[-1]["peeq"] == pytest.approx(3.0 / 110.0 + 1.6 / 121.0, abs=1e-14)
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 10.0 * (3.0 / 110.0 + 1.6 / 121.0), abs=1e-12)
    assert ref[-1]["p"] == pytest.approx(-3.0 / 110.0 + 1.6 / 121.0, abs=1e-14)
    tol = max(_bound(q["dpeeq"], H, S0, 1.5) for q in ref)
    assert tol < 1.2e-3
    _check_uniaxial_stress_levels(rec, ref, [tol] * len(ref), E, 0.3, alpha=1e-3)


# --------------------------------------------------------------------------- (e)

_HT_E = TP.HardeningTable((0.0, 0.1), ((1.0, 2.0), (0.6, 1.1)), temperatures=(0.0, 40.0))


def test_e_temperature_dependent_hardening_table():
    """(d) with a temperature-dependent table: rows (PEEQ 0, 0.1) = (1.0, 2.0)
    at 0 C and (0.6, 1.1) at 40 C, i.e. S0(T) = 1 - 0.01 T, H(T) = 10 -
    0.125 T; the curve is rewritten every increment and read at the tracked
    PEEQ (Abaqus sigma_y(PEEQ, T)).

    Heat: yield at 0.1 T = S0(T) -> T = 9.09; PEEQ(T) = (0.1 T - S0(T)) /
    (E + H(T)) (monotone, so every backward-Euler level is exact):
    T = 40 -> PEEQ_h = 3.4/105 = 0.0323810, sigma = -(0.6 + 5 PEEQ_h) = -0.761905.
    Cool: elastic sigma = E PEEQ_h - 0.1 T until sigma = S0(T) + H(T) PEEQ_h:
    T* = 1.914286 / 0.0859524 = 22.27 C (22.5 C still elastic by 0.0196);
    then dP = (E (PEEQ_h - alpha T) - S0(T) - H(T) PEEQ_h) / (E + H(T)):
    T = 0 -> dP = 1.914286 / 110 = 0.0174026, sigma = 1 + 10 (PEEQ_h + dP)
    = +1.497835. A plastic-WORK internal variable (Kratos without the PEEQ
    sync) reads the 40 C curve at a larger strain (work done on hotter, lower
    curves... and the reverse on cooling) and misses these values."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E = 100.0
    mat = TP.ThermoMaterial("m", E, 0.3, 1e-3, hardening=_HT_E)
    _g, rec = _heat_hold_cool(KM, mods, mat)
    S0 = lambda T: 1.0 - 0.01 * T  # noqa: E731 (hand-written, not the table)
    H = lambda T: 10.0 - 0.125 * T  # noqa: E731
    ref = _uniaxial_stress_thermal_reference(_T_HHC, E, 1e-3, S0, H)
    peeq_h = 3.4 / 105.0
    assert ref[15]["peeq"] == pytest.approx(peeq_h, abs=1e-14)
    assert ref[15]["sigma"] == pytest.approx(-(0.6 + 5.0 * peeq_h), abs=1e-12)
    i_rev = next(i for i, q in enumerate(ref) if i > 17 and q["dpeeq"] > 0)
    assert ref[i_rev]["T"] == 20.0 and ref[i_rev - 1]["T"] == 22.5
    dP = (E * peeq_h - 1.0 - 10.0 * peeq_h) / 110.0
    assert ref[-1]["peeq"] == pytest.approx(peeq_h + dP, abs=1e-14)
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 10.0 * (peeq_h + dP), abs=1e-12)
    tol = max(_bound(q["dpeeq"], q["H"], 0.6, 1.6) for q in ref)
    assert tol < 2e-3
    _check_uniaxial_stress_levels(rec, ref, [tol] * len(ref), E, 0.3, alpha=1e-3)


# --------------------------------------------------------------------------- (f)

def _series_hardening_reference(levels_data):
    """Series bar of unit cells, uniaxial strain, nu = 0, linear hardening per
    cell (S0 = inf: elastic). s_i = sigma + E (th - p/2) (S30); plastic cell
    (c = +-1): c s_i' = S0 + H (PEEQ + D) with p' = p + c D gives
        D = (c (sigma + E th - E p/2) - S0 - H PEEQ) / (E/2 + H),
        eps_i = sigma/E + th + p + (sigma + E th - E p/2 - c (S0 + H PEEQ)) / (E/2 + H),
    elastic eps_i = sigma/E + th + p; sum eps_i = 0 is linear in sigma;
    active set as in S30. Per level ``(E, th, S0, H)``."""
    n = len(levels_data[0][0])
    p, eb = np.zeros(n), np.zeros(n)
    out = []
    for E, th, S0, H in levels_data:
        E, th, S0, H = (np.asarray(a, dtype=float) for a in (E, th, S0, H))
        act: dict[int, float] = {}
        for _ in range(50):
            k = E / 2.0 + H
            a = sum(1.0 / E[i] + (1.0 / k[i] if i in act else 0.0) for i in range(n))
            b = sum(th[i] + p[i] + ((E[i] * th[i] - E[i] * p[i] / 2.0
                                     - act[i] * (S0[i] + H[i] * eb[i])) / k[i] if i in act else 0.0)
                    for i in range(n))
            sig = -b / a
            s = sig + E * (th - p / 2.0)
            D = {i: (c * s[i] - S0[i] - H[i] * eb[i]) / k[i] for i, c in act.items()}
            changed = False
            for i in range(n):
                if i not in act and np.isfinite(S0[i]) and abs(s[i]) > S0[i] + H[i] * eb[i]:
                    act[i] = float(np.sign(s[i]))
                    changed = True
            for i in list(act):
                if D.get(i, 0.0) < 0.0:
                    del act[i]
                    changed = True
            if not changed:
                break
        else:
            raise RuntimeError("reference active set did not settle")
        dpeeq = np.zeros(n)
        for i, c in act.items():
            p[i] += c * D[i]
            eb[i] += D[i]
            dpeeq[i] = D[i]
        out.append({"sigma": sig, "p": p.copy(), "peeq": eb.copy(), "dpeeq": dpeeq,
                    "eps": sig / E + th + p, "E": E, "th": th, "S0": S0, "H": H})
    return out


def test_f_three_cells_two_hardening_groups_residual():
    """Lateral fixed, nu = 0, E = 100; cells A1, A2 metal (alpha 2e-3, table
    rows (PEEQ 0, 0.1) = (1.0, 2.0) at 0 C, (0.4, 0.9) at 20 C: S0 = 1 -
    0.03 T, H = 10 - 0.25 T) and B glass-like (linear elastic law, alpha
    0.5e-3). T: 0 -> (20, 10, 20) (16 inc) -> hold (2) -> 0 (16 inc):
    3 groups (A at 20 C, A at 10 C, B), per-group curve, per-cell PEEQ.

    Reference: the active-set series solution above (independent of Kratos
    and of the module). Path: A1 yields from 10 C (s > 0), A2 yields
    (s < 0) at the end of heating; hold keeps the state; A1 yields in reverse
    on the last cooling increment (PEEQ grows while p shrinks); residual
    sigma_x = -0.4329. Errors from the bound r of each cell: sigma_x <= r,
    sigma_yy <= 2 r, PEEQ <= 4 r/E, u <= 5 r/E per cell."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    ht = TP.HardeningTable((0.0, 0.1), ((1.0, 2.0), (0.4, 0.9)), temperatures=(0.0, 20.0))
    metal = TP.ThermoMaterial("A", 100.0, 0.0, 2e-3, hardening=ht)
    glass = TP.ThermoMaterial("B", 100.0, 0.0, 0.5e-3)
    states = [0.0, np.array([20.0, 10.0, 20.0]), np.array([20.0, 10.0, 20.0]), 0.0]
    model, mp, els, _right = _bar(KM, 3, lateral_fixed=True)
    groups, rec = _run(KM, mods, mp, els, [metal, glass], [0, 0, 1], states, [16, 2, 16])
    assert groups.n_groups == 3 and len(groups.hardening_groups) == 2

    def cell_T(level):
        k = min(int(np.floor(level)), 2)
        S = [np.broadcast_to(np.asarray(s, dtype=float), (3,)) for s in states]
        return (1.0 - (level - k)) * S[k] + (level - k) * S[k + 1]

    data = []
    for r in rec:
        T = cell_T(r["level"])
        S0 = np.append(1.0 - 0.03 * T[:2], np.inf)
        H = np.append(10.0 - 0.25 * T[:2], 0.0)
        data.append((np.full(3, 100.0), np.array([2e-3, 2e-3, 0.5e-3]) * T, S0, H))
    ref = _series_hardening_reference(data)
    i_heat, i_end = 15, len(ref) - 1
    assert np.all(ref[i_heat]["peeq"][:2] > 0.0) and ref[i_heat]["p"][1] < 0.0 < ref[i_heat]["p"][0]
    assert ref[i_end]["dpeeq"][0] > 0.0 and ref[i_end]["p"][0] < ref[i_heat]["p"][0]
    assert ref[i_end]["sigma"] == pytest.approx(-0.43292, abs=1e-5)
    s_max = max(q["S0"][i] + q["H"][i] * q["peeq"][i] for q in ref for i in range(2))
    assert s_max < 1.5
    tol = max(_bound(q["dpeeq"][i], q["H"][i], 0.4, 1.5) for q in ref for i in range(2))
    assert tol < 2e-3
    for r, q in zip(rec, ref):
        for i, st in enumerate(r["stress"]):
            assert st[:, 0] == pytest.approx(q["sigma"], abs=tol), (r["level"], i)
            syy = q["E"][i] * (-q["th"][i] + q["p"][i] / 2.0)
            assert st[:, 1] == pytest.approx(syy, abs=2 * tol) and st[:, 2] == pytest.approx(syy, abs=2 * tol)
        for i in range(2):
            assert r["peeq"][i] == pytest.approx(q["peeq"][i], abs=4 * tol / 100.0), (r["level"], i)
        assert r["u"] == pytest.approx(np.concatenate([[0.0], np.cumsum(q["eps"])]), abs=15 * tol / 100.0)
    for k in (i_heat + 1, i_heat + 2):  # hold
        for st_h, st_0 in zip(rec[k]["stress"], rec[i_heat]["stress"]):
            assert st_h[:, 0] == pytest.approx(st_0[:, 0], abs=2 * tol)


# --------------------------------------------------------------------------- (g)

def test_g_total_lagrangian_finite_strain_hardening():
    """(d) with TotalLagrangianElement3D4N + FiniteStrainIsotropicPlasticity3D
    VonMisesVonMises (Hencky strain, same integrator), scaled: E x 1e4,
    alpha / 1e4, PEEQ grid / 1e4 (H x 1e4) -> same stresses, strains / 1e4
    (eps_max = 4e-6). Tolerance: the hardening bound (scale invariant) plus
    10 eps_max sigma_max for the PK2 / Hencky vs small-strain difference
    (S30 (e))."""
    s = 1e4
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, S0, H = 100.0, 1.0, 10.0
    mat = TP.ThermoMaterial("m", E * s, 0.3, 1e-3 / s,
                            hardening=TP.HardeningTable((0.0, 1.0 / s), ((S0, S0 + H),)))
    _g, rec = _heat_hold_cool(KM, mods, mat, element_name="TotalLagrangianElement3D4N",
                              var=KM.PK2_STRESS_VECTOR, finite=True)
    ref = _uniaxial_stress_thermal_reference(_T_HHC, E, 1e-3, lambda T: S0, lambda T: H)
    eps_max = max(max(abs(1e-3 * q["T"]), abs(q["p"])) for q in ref) / s
    assert eps_max < 5e-6
    tol = max(_bound(q["dpeeq"], H, S0, 1.5) for q in ref) + 10.0 * eps_max * 1.5
    for r, q in zip(rec, ref):
        st = r["stress"][0]
        assert st[:, 0] == pytest.approx(q["sigma"], abs=tol), r["level"]
        assert r["peeq"][0] == pytest.approx(q["peeq"] / s, abs=tol / (E * s)), r["level"]
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 10.0 * (3.0 / 110.0 + 1.6 / 121.0), abs=1e-12)


# --------------------------------------------------------------------------- (h) E03 data

def test_h_e03_hardening_tables():
    """``e03_materials`` against the Coreform TSV values (``data/*.tsv``):
    Alloy 52 (234.4, 0), (580.72, 0.18): slope (580.72 - 234.4)/0.18 = 1924,
    constant after 0.18; 304L rows at 7 temperatures -55..500 C on one 12-point grid, linear
    in T (60 C = mean of the 20 and 100 C rows), constant outside. Plastic
    work by hand trapezoids. Initial yield = e03_materials_initial_yield."""
    mats = TP.e03_materials(_DATA)
    base = TP.e03_materials_initial_yield(_DATA)
    a52 = mats["Alloy52_pins"].hardening
    assert a52.stress([0.0, 0.09, 0.18, 0.3], 20.0) == pytest.approx(
        [234.4, 234.4 + 1924.0 * 0.09, 580.72, 580.72], abs=1e-9)
    assert a52.stress(0.05, -55.0) == a52.stress(0.05, 450.0)   # T independent
    assert a52.work([0.09, 0.18, 0.28], 20.0) == pytest.approx(
        [234.4 * 0.09 + 0.5 * 1924.0 * 0.09 ** 2, 0.18 * (234.4 + 580.72) / 2.0,
         0.18 * (234.4 + 580.72) / 2.0 + 0.1 * 580.72], rel=1e-12)
    ss = mats["SS304L_header"].hardening
    assert ss.temperatures == (-55.0, 20.0, 100.0, 200.0, 300.0, 400.0, 500.0)
    assert len(ss.plastic_strains) == 12 and ss.plastic_strains[-1] == 0.023
    assert ss.stress(0.01, 20.0) == pytest.approx(227.0) and ss.stress(0.01, 500.0) == pytest.approx(128.5)
    assert ss.stress(0.01, 60.0) == pytest.approx((227.0 + 189.5) / 2.0, abs=1e-12)
    assert ss.stress(0.03, 20.0) == pytest.approx(265.2) and ss.stress(0.0, -100.0) == pytest.approx(194.4)
    assert ss.work(0.0001, 20.0) == pytest.approx(1e-4 * (162.0 + 176.2) / 2.0, rel=1e-12)
    for T in (-55.0, 20.0, 100.0, 250.0, 450.0, 500.0):
        for name in ("Alloy52_pins", "SS304L_header"):
            assert mats[name].yield_at(T) == pytest.approx(base[name].yield_at(T), abs=1e-12)
    assert mats["S8061_glass"] == base["S8061_glass"] and not mats["S8061_glass"].plastic
    with pytest.raises(ValueError):
        TP.ThermoMaterial("x", 1.0, 0.3, 0.0, yield_stress=1.0, hardening=a52)
    with pytest.raises(ValueError):
        TP.HardeningTable((0.0, 0.1), ((2.0, 1.0),))   # softening not verified


def test_i_e03_304l_curve_uniaxial_stress():
    """E03 304L curve at 20 C (12 points, first slope 14.2/1e-4 = 142000 MPa,
    flat after PEEQ 0.023), E = 193000 (E(20 C) = E(27 C)), nu = 0.3,
    uniaxial stress, eps: 0 -> 0.002 in 40 increments, -> 0.03 in 56. Answer
    per level from the piecewise solution; per-level bound with the largest
    slope touched (0.4 MPa near yield, steep first segment). At eps = 0.03:
    PEEQ = 0.03 - 265.2 / 193000 > 0.023 -> sigma = 265.2 MPa.

    The curve has one convex kink (slope 3800 -> 5700 at PEEQ 0.001), crossed
    by the increment ending at level 1.0. Kratos' first return iteration uses
    the slope of the start segment, overshoots dlambda by about
    (h2 - h1) db / (3G + h1) and accepts the negative F (``F <= 1e-4 S``):
    sigma ends BELOW the curve by up to (h2 - h1) dPEEQ (x (3G + h2)/(3G + h1)
    ~ 1); hence the term 2 (h2 - h1)+ dPEEQ of ``_bound`` (S42: the bound
    without it failed exactly there, -0.0398 vs 0.0197 MPa)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    ss = TP.e03_materials(_DATA)["SS304L_header"]
    mat = TP.ThermoMaterial("ss20", 193000.0, 0.3, 0.0, hardening=TP.HardeningTable(
        ss.hardening.plastic_strains, (tuple(ss.hardening.curve(20.0)),)))
    grid = [0.0, 0.0001, 0.00025, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.022, 0.023]
    S = [162.0, 176.2, 181.3, 185.7, 187.6, 193.3, 208.0, 227.0, 243.0, 257.6, 262.7, 265.2]
    assert mat.hardening.curve(20.0) == pytest.approx(S, abs=1e-12)
    disp = lambda s: 0.002 * s if s <= 1.0 else 0.002 + 0.028 * (s - 1.0)  # noqa: E731
    model, mp, els, right = _bar(KM, 1, lateral_fixed=False)
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [20.0, 20.0, 20.0], [40, 56],
                   load=_pull(KM, right, disp))
    E = 193000.0
    ref, tols, kinks, prev = [], [], [], 0.0
    for r in rec:
        e = disp(r["level"])
        eb = _peeq_uniaxial_stress(e, E, grid, S)
        ref.append({"T": 0.0, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        h_max, dh_up = _slopes_touched(np.array(grid), np.array(S), prev, eb)
        tols.append(_bound(eb - prev, h_max, S[0], E * (e - eb), dh_up))
        kinks.append(dh_up > 0.0)
        prev = eb
    assert ref[-1]["sigma"] == pytest.approx(265.2, abs=1e-9) and ref[-1]["peeq"] > 0.023
    assert max(tols) < 0.5
    # the only convex kink (3800 -> 5700 MPa at PEEQ 0.001) is crossed once
    assert sum(kinks) == 1 and rec[kinks.index(True)]["level"] == pytest.approx(1.0)
    _check_uniaxial_stress_levels(rec, ref, tols, E, 0.3)
