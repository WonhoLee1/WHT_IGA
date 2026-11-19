# -*- coding: utf-8 -*-
"""S44: isotropic hardening on the S28/S30 InitialState path (``thermal_plasticity``).

Mapping (module doc, S44): ``HARDENING_CURVE = 6`` (CurveDefinedByPoints) with
``PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE`` = PEEQ grid (+ a flat far
point). Kratos hardens on ``kappa = W_p l_c / G_f``; on segment j it uses
``S = a sqrt(1 - b kappa)``, i.e. S^2 linear in W_p = int S dPEEQ, which is
exactly S linear in PEEQ (S^2 = S_j^2 + 2 h_j (W - W_j)), so ``l_c`` and
``G_f`` cancel; the form needs a^2 = S_j^2 - 2 h_j W_j > 0 (test (l)).
``CellThermalGroups`` tracks PEEQ from the converged plastic strain
increments and writes the matching ``PLASTIC_DISSIPATION`` before every
increment (Abaqus semantics sigma_y(PEEQ, T) for temperature-dependent tables).

A-priori bound on the yield error of one Kratos increment, delta = q - S(PEEQ)
(S44, ``_bound``):

    |delta| <= b = 1e-4 S + 2 rho H^2 dP^2 / S_min  [+ 2 (h2 - h1)+ dP],
    rho = 1 + (H + 3G) dP / S_min

* 1e-4 S: Kratos' return tolerance (``F <= 1e-4 threshold``; trials in that
  band are also treated as elastic).
* Kratos integrates ``W_p = sum q dlambda`` over its return iterations with
  the end-of-iteration equivalent stress q (decreasing from q1), so within
  one increment ``W_p`` exceeds ``int S dPEEQ`` by at most ``(q1 - S_n) dP``
  with ``q1 - S_n = A3 dlambda_1 <= H (q0 / S_n) dP`` (A3 = H q0 / S_n is
  Kratos' hardening term, q0 the trial). The radial return gives
  q0 = q_{n+1} + 3G dP <= S_n + (H + 3G) dP, hence q0 / S_n <= rho.
  dS/dW = H / S turns the excess into stress. Factor 2: the analytic dP
  stands for Kratos' (they differ by O(delta / E)) and higher-order terms.
* Only when an increment crosses a CONVEX kink (slope rises h1 -> h2): the
  first return iteration uses h1, overshoots, and Kratos accepts the negative
  F = -(h2 - h1)(dlambda_1 - d) (``F <= 1e-4 S``): + 2 (h2 - h1)+ dP (test (i)).
* With the PEEQ sync the work excess does not accumulate. dP, H (largest
  slope touched) and S come from the analytic solution.

From delta to the stress error (one cell, total strain prescribed, uniaxial
stress): a flow phase is memoryless (sigma_err = E delta / (E + h)), elastic
levels freeze the plastic strain error, and a reversal re-reads the earlier
PEEQ error; :func:`_chain_tolerances` carries this (tests (d), (e), (g), (k)).

Models: bars of cells [0, n L] x [0, L]^2, one element per cell, u_x = 0 at
x = 0 and prescribed (or 0) at x = n L. Lateral FREE (symmetry planes y = 0,
z = 0) gives uniaxial STRESS, lateral FIXED uniaxial STRAIN (S30). (j) is a
2 x 2 x 2 patch in simple shear.
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
_SQRT3 = float(np.sqrt(3.0))


# --------------------------------------------------------------------------- helpers

def _model(KM, name="S44"):
    model = KM.Model()
    mp = model.CreateModelPart(name)
    for v in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(v)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def _bar(KM, nx, lateral_fixed, L=1.0, degree=1, element_name="SmallDisplacementElement3D4N"):
    model, mp = _model(KM)
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


def _bound(dpeeq, H, s_min, s_max, G, dh_up=0.0):
    """A-priori yield-error bound b of one increment (module docstring);
    ``dh_up``: slope increase at a convex kink crossed by the increment."""
    rho = 1.0 + (H + 3.0 * G) * dpeeq / s_min
    return 1e-4 * s_max + 2.0 * rho * H ** 2 * dpeeq ** 2 / s_min + 2.0 * dh_up * dpeeq


def _chain_tolerances(ref, bounds):
    """Per-level (stress, plastic-strain) tolerances for ONE cell in uniaxial
    stress with prescribed total strain (strain-driven or clamped thermal)
    from the per-level yield-error bounds ``bounds`` (b = 1e-4 |sigma| on
    elastic levels).

    In a flow phase of sign c that started with errors (dPEEQ_0, dp_0),
    sigma = -E (eps_th + p) [or E (eps - p)] and c sigma = S(PEEQ) + delta
    give, with X = dPEEQ_0 - c dp_0,
        dPEEQ = (E X - delta) / (E + h),  c dsigma = h dPEEQ + delta = E (h X + delta) / (E + h).
    First phase (X = 0): |dsigma| <= b, |dPEEQ| <= b / E, memoryless.
    A new phase (reversal) starts with |X| <= |dPEEQ_0| + |dp_0| <= 2 pt =:
    carry (pt bounds |dp| and |dPEEQ| at its start), so |dsigma| <= b + h
    carry, |dPEEQ| <= carry + b/E and |dp| <= |dp_0| + |dPEEQ_0| + |dPEEQ|
    <= carry + |dPEEQ|. Elastic levels freeze dp: |dsigma| <= b + E pt.
    ``ref`` entries need E, H (plastic levels), p (signed) and dpeeq."""
    tols, pts = [], []
    pt = carry = phase = p_prev = 0.0
    for q, b in zip(ref, bounds):
        if q["dpeeq"] > 0.0:
            c = float(np.sign(q["p"] - p_prev))
            if c != phase:
                carry, phase = 2.0 * pt, c
            dP = carry + b / q["E"]
            pt = carry + dP
            tols.append(b + q["H"] * carry)
        else:
            tols.append(b + q["E"] * pt)
        pts.append(pt)
        p_prev = q["p"]
    return tols, pts


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
        out.append({"T": T, "th": alpha * T, "E": E, "sigma": trial - c * E * d, "p": p,
                    "peeq": eb, "dpeeq": d, "H": H(T), "S0": S0(T)})
    return out


def _check_uniaxial_stress_levels(rec, ref, tols, nu, L=1.0, peeq_tols=None):
    """sigma_xx to tol, other components below tol, PEEQ to the plastic-strain
    tolerance (tol / E unless given), u_y(top) = L eps_yy with eps_yy = th -
    nu sigma / E - p / 2 to L (nu tol / E + pt / 2)."""
    assert len(rec) == len(ref) == len(tols)
    for k, (r, q, tol) in enumerate(zip(rec, ref, tols)):
        E = q["E"]
        pt = tol / E if peeq_tols is None else peeq_tols[k]
        for st in r["stress"]:
            assert st[:, 0] == pytest.approx(q["sigma"], abs=tol), r["level"]
            assert np.abs(st[:, 1:]).max() < tol, r["level"]
        for pq in r["peeq"]:
            assert pq == pytest.approx(q["peeq"], abs=pt), r["level"]
        eps_y = q["th"] - nu * q["sigma"] / E - q["p"] / 2.0
        assert r["uy_top"] == pytest.approx([L * eps_y] * len(r["uy_top"]),
                                            abs=L * (nu * tol / E + 0.5 * pt)), r["level"]


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
    end PEEQ = 4.2/110 = 0.0381818, sigma = 1.381818. Monotone, strain
    prescribed: per-level tolerance b (``_chain_tolerances`` first phase,
    memoryless).
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
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable((0.0, 1.0), ((S0, S0 + H),)))
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [0.0, 0.0], 20,
                   load=_pull(KM, right, lambda s: s * eps_max * L * nx), top=L)
    ref, tols, prev = [], [], 0.0
    for r in rec:
        e = r["level"] * eps_max
        eb = _peeq_uniaxial_stress(e, E, [0.0, 1.0], [S0, S0 + H])
        ref.append({"T": 0.0, "th": 0.0, "E": E, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        tols.append(_bound(eb - prev, H, S0, E * (e - eb), G))
        prev = eb
    assert ref[-1]["peeq"] == pytest.approx(4.2 / 110.0, abs=1e-14)
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 42.0 / 110.0, abs=1e-12)
    assert max(tols) < 2e-3   # the bound is meaningful: 0.15 % of sigma
    _check_uniaxial_stress_levels(rec, ref, tols, nu, L=L)


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
    All strains prescribed, so a yield error delta (|delta| <= b) gives
    dp = -delta / (3G + H), dq = 3G delta / (3G + H): sigma_x <= 2/3 b,
    sigma_y <= b/3, PEEQ <= b/(3G + H) (asserted: sigma to b)."""
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
        tol = _bound(p - prev, H, S0, sx - sy, G)
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
    after 0.04 (Abaqus); eps: 0 -> 0.07 in 200 increments.

    Per segment j: PEEQ = (E eps - S_j + H_j e_j) / (E + H_j); segment ends
    at eps = e_j + S_j / E: 0.01, 0.019, 0.037, 0.058; at eps = 0.07 PEEQ =
    0.052 in the flat extension, sigma = 1.8. Per-level bound with the
    largest slope touched by the increment (80 near yield)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, right = _bar(KM, 1, lateral_fixed=False)
    E, nu, eps_max = 100.0, 0.3, 0.07
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable(_GRID_C, (_S_C,)))
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [0.0, 0.0], 200,
                   load=_pull(KM, right, lambda s: s * eps_max))
    ref, tols, prev = [], [], 0.0
    for r in rec:
        e = r["level"] * eps_max
        eb = _peeq_uniaxial_stress(e, E, _GRID_C, _S_C)
        ref.append({"T": 0.0, "th": 0.0, "E": E, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        h_max, dh_up = _slopes_touched(_GRID_C, _S_C, prev, eb)
        assert dh_up == 0.0  # concave curve: no over-return term
        tols.append(_bound(eb - prev, h_max, _S_C[0], E * (e - eb), G))
        prev = eb
    ends = [g + s / E for g, s in zip(_GRID_C, _S_C)]
    assert ends == pytest.approx([0.01, 0.019, 0.037, 0.058], abs=1e-15)
    assert ref[-1]["peeq"] == pytest.approx(0.052, abs=1e-14) and ref[-1]["sigma"] == pytest.approx(1.8)
    # every segment and the flat extension is visited
    hit = {int(np.searchsorted(_GRID_C, q["peeq"], side="right")) for q in ref if q["peeq"] > 0}
    assert hit == {1, 2, 3, 4}
    assert max(tols) < 2e-3
    _check_uniaxial_stress_levels(rec, ref, tols, nu)


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
    T = 0 -> dP = 1.454545/110 = 1.6/121 = 0.0132231, PEEQ = 0.0404959,
    residual sigma = 1 + 10 PEEQ = +1.404959 (perfect plasticity: +1.0).
    Each phase is monotone, so every backward-Euler level is exact.
    Lateral u_y(top) = alpha T - nu sigma/E - p/2. Tolerances: per-level b
    carried through hold, cooling and the reversal (``_chain_tolerances``)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, nu, S0, H = 100.0, 0.3, 1.0, 10.0
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 1e-3, hardening=TP.HardeningTable((0.0, 1.0), ((S0, S0 + H),)))
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
    tols, pts = _chain_tolerances(
        ref, [_bound(q["dpeeq"], H, S0, abs(q["sigma"]), G) for q in ref])
    assert max(tols) < 2e-3
    _check_uniaxial_stress_levels(rec, ref, tols, nu, peeq_tols=pts)


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
    = +1.497835. Without the PEEQ sync Kratos would read each curve at the
    PEEQ whose work ON THAT CURVE equals the work done along the path; the
    work was done on hotter, lower curves, so the 0 C curve would be read at a
    smaller PEEQ (lower threshold, earlier reverse yield)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, nu = 100.0, 0.3
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 1e-3, hardening=_HT_E)
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
    tols, pts = _chain_tolerances(
        ref, [_bound(q["dpeeq"], q["H"], q["S0"], abs(q["sigma"]), G) for q in ref])
    assert max(tols) < 2e-3
    _check_uniaxial_stress_levels(rec, ref, tols, nu, peeq_tols=pts)


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
        sign = np.zeros(n)
        for i, c in act.items():
            p[i] += c * D[i]
            eb[i] += D[i]
            dpeeq[i] = D[i]
            sign[i] = c
        out.append({"sigma": sig, "p": p.copy(), "peeq": eb.copy(), "dpeeq": dpeeq, "sign": sign,
                    "eps": sig / E + th + p, "E": E, "th": th, "S0": S0, "H": H})
    return out


def test_f_three_cells_two_hardening_groups_residual():
    """Lateral fixed, nu = 0, E = 100; cells A1, A2 metal (alpha 2e-3, table
    rows (PEEQ 0, 0.1) = (1.0, 2.0) at 0 C, (0.4, 0.9) at 20 C: S0 = 1 -
    0.03 T, H = 10 - 0.25 T) and B glass-like (linear elastic law, alpha
    0.5e-3). T: 0 -> (20, 10, 20) (16 inc) -> hold (2) -> 0 (16 inc):
    3 groups (A at 20 C, A at 10 C, B), per-group curve, per-cell PEEQ.

    Reference: the active-set series solution above (independent of Kratos
    and of the module). Hand check, end of heating: A1 (20 C: S0 0.4, H 5,
    th 0.04, s > 0) and A2 (10 C: S0 0.7, H 7.5, th 0.02, s < 0) both on the
    yield surface, B elastic (th 0.01):
        sigma = -3.6 + 55 p1 = -2.7 + 57.5 p2,  0.03 sigma + 0.07 + p1 + p2 = 0
        -> sigma = -0.182411 / 0.0655731 = -2.781796, p1 = 0.014876, p2 = -0.001423.
    Hold keeps the state. Cooling is elastic until A1 yields in reverse on
    the last increment (PEEQ grows while p shrinks); at 0 C (S0 1, H 10, th 0):
        -sigma + 50 (p1 - D) = 1 + 10 (p1 + D),  0.03 sigma + p1 - D + p2 = 0
        -> residual sigma = -0.0202029 / 0.0466667 = -0.432920, D = 4.66e-4
    (A2 stays elastic: |s2| = 0.362 < 1.014).

    Tolerances from the per-increment bound r = max b (series, n = 3 cells,
    E/2 + h is the plastic modulus of s): before the reversal |dsigma| <= r
    (weighted mean of the cell errors), |dPEEQ_i| <= 2 r / (E/2 + h) <= 4 r/E,
    |dsigma_yy| = E |dp| / 2 <= 2 r, |du| <= 5 r/E per cell. After A1's
    reversal its earlier PEEQ error re-enters: with h = 10 at 0 C and
    |dPEEQ_1,0| <= 2 r / 55, |dsigma| <= (0.03 + 0.0242 + 0.0167) r / 0.0467
    = 1.52 r, |dPEEQ_1| <= 10.3 r/E, |dsigma_yy,1| <= 5.2 r, |du| <= 16.8 r/E;
    asserted with 2 r, 11 r/E, 6 r, 20 r/E."""
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
    assert ref[i_heat]["sigma"] == pytest.approx(
        -(0.07 + 3.6 / 55.0 + 2.7 / 57.5) / (0.03 + 1.0 / 55.0 + 1.0 / 57.5), abs=1e-12)
    assert ref[i_heat]["p"][:2] == pytest.approx([0.0148764316, -0.0014225437], abs=1e-9)
    assert ref[i_end]["dpeeq"][0] > 0.0 and ref[i_end]["p"][0] < ref[i_heat]["p"][0]
    assert ref[i_end]["sigma"] == pytest.approx(-0.0202029335 / (0.03 + 1.0 / 60.0), abs=1e-8)
    # reversal: first level where a metal cell flows against its earlier sign
    seen = np.zeros(2)
    i_rev = None
    for k, q in enumerate(ref):
        for i in range(2):
            if q["sign"][i] != 0.0:
                if seen[i] != 0.0 and q["sign"][i] != seen[i] and i_rev is None:
                    i_rev = k
                seen[i] = q["sign"][i]
    assert i_rev == i_end and ref[i_end]["sign"][1] == 0.0   # only A1 reverses, last increment
    s_max = max(q["S0"][i] + q["H"][i] * q["peeq"][i] for q in ref for i in range(2))
    assert s_max < 1.5
    r_b = max(_bound(q["dpeeq"][i], q["H"][i], 0.4, 1.5, 50.0) for q in ref for i in range(2))
    assert r_b < 2e-3
    for k, (r, q) in enumerate(zip(rec, ref)):
        fs, fyy, fp, fu = (1.0, 2.0, 4.0, 15.0) if k < i_rev else (2.0, 6.0, 11.0, 20.0)
        for i, st in enumerate(r["stress"]):
            assert st[:, 0] == pytest.approx(q["sigma"], abs=fs * r_b), (r["level"], i)
            syy = q["E"][i] * (-q["th"][i] + q["p"][i] / 2.0)
            assert st[:, 1] == pytest.approx(syy, abs=fyy * r_b), (r["level"], i)
            assert st[:, 2] == pytest.approx(syy, abs=fyy * r_b), (r["level"], i)
        for i in range(2):
            assert r["peeq"][i] == pytest.approx(q["peeq"][i], abs=fp * r_b / 100.0), (r["level"], i)
        assert r["u"] == pytest.approx(np.concatenate([[0.0], np.cumsum(q["eps"])]), abs=fu * r_b / 100.0)
    for k in (i_heat + 1, i_heat + 2):  # hold
        for st_h, st_0 in zip(rec[k]["stress"], rec[i_heat]["stress"]):
            assert st_h[:, 0] == pytest.approx(st_0[:, 0], abs=2 * r_b)


# --------------------------------------------------------------------------- (g)

def test_g_total_lagrangian_finite_strain_hardening():
    """(d) with TotalLagrangianElement3D4N + FiniteStrainIsotropicPlasticity3D
    VonMisesVonMises (Hencky strain, same integrator), scaled: E x 1e4,
    alpha / 1e4, PEEQ grid / 1e4 (H x 1e4) -> same stresses, strains / 1e4
    (eps_max = 4e-6). Tolerance: the (d) chain tolerances (b is scale
    invariant: H dP and G dP are) plus 10 eps_max sigma_max for the PK2 /
    Hencky vs small-strain difference (S30 (e))."""
    s = 1e4
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, nu, S0, H = 100.0, 0.3, 1.0, 10.0
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E * s, nu, 1e-3 / s,
                            hardening=TP.HardeningTable((0.0, 1.0 / s), ((S0, S0 + H),)))
    _g, rec = _heat_hold_cool(KM, mods, mat, element_name="TotalLagrangianElement3D4N",
                              var=KM.PK2_STRESS_VECTOR, finite=True)
    ref = _uniaxial_stress_thermal_reference(_T_HHC, E, 1e-3, lambda T: S0, lambda T: H)
    eps_max = max(max(abs(1e-3 * q["T"]), abs(q["p"])) for q in ref) / s
    assert eps_max < 5e-6
    geo = 10.0 * eps_max * 1.5
    tols, pts = _chain_tolerances(
        ref, [_bound(q["dpeeq"], H, S0, abs(q["sigma"]), G) for q in ref])
    assert len(rec) == len(ref)
    for r, q, tol, pt in zip(rec, ref, tols, pts):
        st = r["stress"][0]
        assert st[:, 0] == pytest.approx(q["sigma"], abs=tol + geo), r["level"]
        assert r["peeq"][0] == pytest.approx(q["peeq"] / s, abs=(pt + geo / E) / s), r["level"]
    assert ref[-1]["sigma"] == pytest.approx(1.0 + 10.0 * (3.0 / 110.0 + 1.6 / 121.0), abs=1e-12)


# --------------------------------------------------------------------------- (h) E03 data

def test_h_e03_hardening_tables():
    """``e03_materials`` against the Coreform TSV values (``data/*.tsv``):
    Alloy 52 (234.4, 0), (580.72, 0.18): slope (580.72 - 234.4)/0.18 = 1924,
    constant after 0.18; 304L rows at 7 temperatures -55..500 C on one
    12-point grid, linear in T (60 C = mean of the 20 and 100 C rows),
    constant outside. Plastic work by hand trapezoids. Initial yield =
    e03_materials_initial_yield."""
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
    slope touched (steep first segment near yield). At eps = 0.03:
    PEEQ = 0.03 - 265.2 / 193000 > 0.023 -> sigma = 265.2 MPa.

    The curve has one convex kink (slope 3800 -> 5700 at PEEQ 0.001), crossed
    by the increment ending at level 1.0: at eps = 0.00195 PEEQ = (376.35 -
    185.7 + 1.9) / 196800 = 0.000978 (slope 3800), at eps = 0.002 PEEQ =
    (386 - 187.6 + 5.7) / 198700 = 0.001027 (slope 5700). Kratos' first
    return iteration uses the slope of the start segment, overshoots dlambda
    and accepts the negative F (``F <= 1e-4 S``): sigma ends BELOW the curve
    by up to (h2 - h1) dPEEQ; hence the term 2 (h2 - h1)+ dPEEQ of
    ``_bound``."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    ss = TP.e03_materials(_DATA)["SS304L_header"]
    E, nu = 193000.0, 0.3
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("ss20", E, nu, 0.0, hardening=TP.HardeningTable(
        ss.hardening.plastic_strains, (tuple(ss.hardening.curve(20.0)),)))
    grid = [0.0, 0.0001, 0.00025, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.022, 0.023]
    S = [162.0, 176.2, 181.3, 185.7, 187.6, 193.3, 208.0, 227.0, 243.0, 257.6, 262.7, 265.2]
    assert mat.hardening.curve(20.0) == pytest.approx(S, abs=1e-12)
    disp = lambda s: 0.002 * s if s <= 1.0 else 0.002 + 0.028 * (s - 1.0)  # noqa: E731
    model, mp, els, right = _bar(KM, 1, lateral_fixed=False)
    _g, rec = _run(KM, mods, mp, els, [mat], 0, [20.0, 20.0, 20.0], [40, 56],
                   load=_pull(KM, right, disp))
    ref, tols, kinks, prev = [], [], [], 0.0
    for r in rec:
        e = disp(r["level"])
        eb = _peeq_uniaxial_stress(e, E, grid, S)
        ref.append({"T": 20.0, "th": 0.0, "E": E, "sigma": E * (e - eb), "p": eb, "peeq": eb})
        h_max, dh_up = _slopes_touched(np.array(grid), np.array(S), prev, eb)
        tols.append(_bound(eb - prev, h_max, S[0], E * (e - eb), G, dh_up))
        kinks.append(dh_up > 0.0)
        prev = eb
    assert ref[-1]["sigma"] == pytest.approx(265.2, abs=1e-9) and ref[-1]["peeq"] > 0.023
    assert max(tols) < 0.5   # 0.3 % of the initial yield
    # the only convex kink (3800 -> 5700 MPa at PEEQ 0.001) is crossed once
    assert sum(kinks) == 1 and rec[kinks.index(True)]["level"] == pytest.approx(1.0)
    _check_uniaxial_stress_levels(rec, ref, tols, nu)


# --------------------------------------------------------------------------- (j) shear

def test_j_simple_shear_patch_multilinear():
    """Simple shear on a 2 x 2 x 2 patch of p = 1 cells on [0, 2]^3: every
    boundary control point has u = (gamma y, 0, 0), the centre one is free;
    the homogeneous field is an exact discrete equilibrium (constant stress,
    int grad N_I dV = 0 for the interior function), so every integration
    point sees eps_xy = gamma / 2 and nothing else. E = 100, nu = 0.3
    (G = 38.4615), curve of (c), gamma: 0 -> 0.12 in 120 increments.

    J2 simple shear: engineering plastic shear gamma_p, PEEQ = gamma_p /
    sqrt(3) (sqrt(2/3 eps_p : eps_p) with eps_p,xy = gamma_p / 2),
    tau = G (gamma - sqrt(3) PEEQ), sqrt(3) tau = S(PEEQ):
        PEEQ = (sqrt(3) G gamma - S_j + h_j e_j) / (3G + h_j)  on segment j,
    i.e. the uniaxial formula with E -> 3G, eps -> gamma / sqrt(3); yield at
    gamma = 1 / (sqrt(3) G) = 0.0150; at gamma = 0.12 PEEQ = (7.99398 - 1.8)
    / 115.385 = 0.05368 > 0.04 (flat), tau = 1.8 / sqrt(3) = 1.03923.
    The PEEQ sync must count the shear plastic strain (Voigt, engineering
    gamma) with weight 1/2 in sqrt(2/3 (sum eps_ii^2 + 1/2 sum gamma_ij^2));
    uniaxial tests (a)-(i) never exercise that term. Errors from the yield
    error b (all strains prescribed): q_err = 3G delta / (3G + h) <= b, so
    tau to b / sqrt(3), PEEQ to b / (3G); other stress components (exactly 0
    for the radial return of a pure-shear trial) below b; the free centre
    node at (gamma, 0, 0) to 1e-6 gamma (Newton displacement criterion 1e-8
    relative)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp = _model(KM, "S44shear")
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 2.0, 2.0), n_cells=(2, 2, 2), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = attach_cell_elements(KM, mp, vol, box, mp.CreateNewProperties(1))
    mp.CloneTimeStep(0.0)
    boundary, centre = [], []
    for n in mp.Nodes:
        if any(abs(c - 1.0) > 1e-12 for c in (n.X0, n.Y0, n.Z0)):
            for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
                n.Fix(v)
            boundary.append(n)
        else:
            centre.append(n)
    assert len(els) == 8 and len(boundary) == 26 and len(centre) == 1
    E, nu, g_max = 100.0, 0.3, 0.12
    G = E / (2.0 * (1.0 + nu))
    mat = TP.ThermoMaterial("m", E, nu, 0.0, hardening=TP.HardeningTable(_GRID_C, (_S_C,)))
    groups = TP.CellThermalGroups(KM, mods, mp, els, [mat], 0, [0.0, 0.0])

    def load(level):
        for n in boundary:
            n.SetSolutionStepValue(KM.DISPLACEMENT_X, g_max * level * n.Y0)

    rec: list[dict] = []

    def obs(level):
        pq = groups.equivalent_plastic_strain()
        rec.append({"level": level,
                    "stress": np.vstack([_ip(KM, mp, el, KM.CAUCHY_STRESS_VECTOR) for el in els]),
                    "peeq": np.concatenate([pq[int(el.Id)] for el in els]),
                    "uc": np.array(centre[0].GetSolutionStepValue(KM.DISPLACEMENT))})

    out = TP.run_thermal_steps(mp, KM, lambda: NL.make_newton_strategy(KM, mp), groups, 120,
                               on_converged=obs, mechanical_load=load)
    assert out["n_solves"] == len(out["levels"]) == len(rec) == 120
    prev = 0.0
    for r in rec:
        gam = g_max * r["level"]
        P = _peeq_uniaxial_stress(gam / _SQRT3, 3.0 * G, _GRID_C, _S_C)
        tau = G * (gam - _SQRT3 * P)
        if P > 0.0:
            assert _SQRT3 * tau == pytest.approx(np.interp(P, _GRID_C, _S_C), abs=1e-12)
        h_max, _dh = _slopes_touched(_GRID_C, _S_C, prev, P)
        b = _bound(P - prev, h_max, _S_C[0], _SQRT3 * tau, G)
        prev = P
        st = r["stress"]
        assert st[:, 3] == pytest.approx(tau, abs=b / _SQRT3), r["level"]
        assert np.abs(st[:, [0, 1, 2, 4, 5]]).max() < b, r["level"]
        assert r["peeq"] == pytest.approx(P, abs=b / (3.0 * G)), r["level"]
        assert r["uc"] == pytest.approx([gam, 0.0, 0.0], abs=1e-6 * g_max), r["level"]
    assert P == pytest.approx((_SQRT3 * G * g_max - 1.8) / (3.0 * G), abs=1e-14) and P > 0.04
    assert tau == pytest.approx(1.8 / _SQRT3, abs=1e-12)


# --------------------------------------------------------------------------- (k) E03 history

_E03_STATES = [450.0, 20.0, -55.0, 200.0]
_E03_INCS = [43, 15, 26]


def _e03_curves(name):
    """E, nu, secant alpha, hardening curve at T straight from the Coreform
    data (TSV tables; Alloy 52 constants from doc Table 11: E 165500, nu
    0.29, alpha 9.95e-6), independent of ``thermal_plasticity``."""
    ld = lambda n: np.loadtxt(_DATA / n, delimiter="\t", ndmin=2)  # noqa: E731
    if name == "SS304L_header":
        el, ex, pl = (ld(f"abaqus_ss304l_{k}.tsv") for k in ("elastic", "expansion", "plastic"))
        temps = np.unique(pl[:, 2])
        grid = pl[pl[:, 2] == temps[0]][:, 1]
        rows = np.array([pl[pl[:, 2] == t][:, 0] for t in temps])
        return (lambda T: float(np.interp(T, el[:, 2], el[:, 0])), 0.3,
                lambda T: float(np.interp(T, ex[:, 1], ex[:, 0])), grid,
                lambda T: np.array([np.interp(T, temps, rows[:, j]) for j in range(grid.size)]))
    a52 = ld("abaqus_alloy52_plastic.tsv")
    return (lambda T: 165500.0, 0.29, lambda T: 9.95e-6, a52[:, 1], lambda T: a52[:, 0].copy())


def _e03_reference(name, T_path, t_init=450.0, t_zero=20.0):
    """Clamped cell, uniaxial stress, backward Euler at the level temperature
    (Abaqus): eps_th = alpha(T)(T - 20) - alpha(450)(450 - 20), trial
    sigma = -E(T) (eps_th + p), threshold S(PEEQ_n, T); flow solves
    E dP + S(PEEQ_n + dP, T) = |trial| by bisection (monotone)."""
    E_of, nu, a_of, grid, S_of = _e03_curves(name)
    p = pb = 0.0
    out = []
    for T in T_path:
        E, S = E_of(T), S_of(T)
        th = a_of(T) * (T - t_zero) - a_of(t_init) * (t_init - t_zero)
        trial = -E * (th + p)
        f = lambda d: abs(trial) - E * d - np.interp(pb + d, grid, S)  # noqa: E731
        d = 0.0
        if f(0.0) > 0.0:
            lo, hi = 0.0, abs(trial) / E
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(mid) > 0.0 else (lo, mid)
            d = 0.5 * (lo + hi)
        c = float(np.sign(trial))
        h_max, dh_up = _slopes_touched(grid, S, pb, pb + d)
        out.append({"T": T, "E": E, "nu": nu, "th": th, "sigma": trial - c * E * d,
                    "p": p + c * d, "peeq": pb + d, "dpeeq": d, "H": h_max, "dh_up": dh_up,
                    "S0": float(S[0])})
        p += c * d
        pb += d
    return out


def _e03_T(level):
    k = min(int(np.floor(level)), len(_E03_STATES) - 2)
    w = level - k
    return (1.0 - w) * _E03_STATES[k] + w * _E03_STATES[k + 1]


@pytest.mark.parametrize("name", ["SS304L_header", "Alloy52_pins"])
def test_k_e03_materials_thermal_history(name):
    """One clamped cell (uniaxial stress) of an E03 metal with the full
    ``e03_materials`` data (E(T), secant alpha(T) about 20 C, sigma_y(PEEQ, T))
    through the E03 history 450 (stress free) -> 20 -> -55 -> 200 C in 43 / 15
    / 26 increments; reference ``_e03_reference`` from the raw tables.

    Each flow phase is monotone, so its levels solve E(T)(m(T) - PEEQ) =
    S(PEEQ, T) exactly (m = -eps_th, increment independent).
    304L: m(20) = 1.83112744186e-5 * 430 = 7.873848e-3, segment slope 3800:
        PEEQ = (193000 m - 208 + 19) / 196800 = 6.761446e-3, sigma = 214.6935;
    -55 C: m = 7.873848e-3 + 1.72e-5 * 75 = 9.163848e-3, slope 4560:
        PEEQ = (193000 m - 249.6 + 22.8) / 197560 = 7.804326e-3, sigma = 262.3877;
    200 C: E = 183290.91, m = 7.873848e-3 - 1.77697666667e-5 * 180 =
    4.675290e-3, reverse flow D from E (PEEQ_-55 - D - m) = 168.5 + 2680
    (PEEQ_-55 + D - 0.01): D = 2.209530e-3, PEEQ = 1.0013856e-2,
    sigma = -168.5371.
    Alloy 52 (165500, linear 1924): yield from 307.7 C, PEEQ = (E m - 234.4)
    / 167424 with m = 9.95e-6 (450 - T): 20 C 2.829294e-3, sigma 239.8436;
    -55 C 3.566968e-3, sigma 241.2628; 200 C elastic, sigma = 165500
    (2.4875e-3 - 3.566968e-3) = -178.6520 (> -241.26: no reverse flow).
    Tolerances: per-level b with the slopes touched on the curve at T, the
    convex-kink term, chained through elastic stretches and the reversal
    (``_chain_tolerances``)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mat = TP.e03_materials(_DATA)[name]
    model, mp, els, _right = _bar(KM, 1, lateral_fixed=False)
    _g, rec = _run(KM, mods, mp, els, [mat], 0, _E03_STATES, _E03_INCS)
    ref = _e03_reference(name, [_e03_T(r["level"]) for r in rec])
    i1, i2, i3 = 42, 57, 83
    assert [rec[i]["level"] for i in (i1, i2, i3)] == pytest.approx([1.0, 2.0, 3.0], abs=1e-12)
    if name == "SS304L_header":
        m20 = 1.83112744186e-5 * 430.0
        P1 = (193000.0 * m20 - 189.0) / 196800.0
        m55 = m20 + 1.72e-5 * 75.0
        P2 = (193000.0 * m55 - 226.8) / 197560.0
        E3, m3, h3 = 187000.0 - 4000.0 * 51.0 / 55.0, m20 - 1.77697666667e-5 * 180.0, 2680.0
        D = (E3 * (P2 - m3) - 168.5 - h3 * (P2 - 0.01)) / (E3 + h3)
        anchors = [(P1, 208.0 + 3800.0 * (P1 - 0.005)), (P2, 249.6 + 4560.0 * (P2 - 0.005)),
                   (P2 + D, -(168.5 + h3 * (P2 + D - 0.01)))]
        assert [a[1] for a in anchors] == pytest.approx([214.6935, 262.3877, -168.5371], abs=1e-4)
        assert ref[i3]["p"] == pytest.approx(P2 - D, abs=1e-12)   # reversed flow in step 3
    else:
        P1 = (165500.0 * 9.95e-6 * 430.0 - 234.4) / 167424.0
        P2 = (165500.0 * 9.95e-6 * 505.0 - 234.4) / 167424.0
        anchors = [(P1, 234.4 + 1924.0 * P1), (P2, 234.4 + 1924.0 * P2),
                   (P2, 165500.0 * (9.95e-6 * 250.0 - P2))]
        assert [a[1] for a in anchors] == pytest.approx([239.8436, 241.2628, -178.6520], abs=1e-4)
        assert all(q["dpeeq"] == 0.0 for q in ref[i2 + 1:])          # step 3 elastic
    for i, (P, sig) in zip((i1, i2, i3), anchors):
        assert ref[i]["peeq"] == pytest.approx(P, rel=1e-9)
        assert ref[i]["sigma"] == pytest.approx(sig, rel=1e-9)
    bounds = [_bound(q["dpeeq"], q["H"], q["S0"], abs(q["sigma"]), q["E"] / (2.0 * (1.0 + q["nu"])),
                     q["dh_up"]) for q in ref]
    tols, pts = _chain_tolerances(ref, bounds)
    # meaningful bound: below 1 % of the 20 C stress (largest at first yield on
    # the steep first 304L segments, 34000-142000 MPa)
    assert max(tols) < 0.01 * anchors[0][1]
    _check_uniaxial_stress_levels(rec, ref, tols, ref[0]["nu"], peeq_tols=pts)


# --------------------------------------------------------------------------- (l) guard

def test_l_point_curve_representability_guard():
    """Kratos' segment form needs a_j^2 = S_j^2 - 2 h_j W_j > 0
    (``point_curve_margin`` = a_j^2 / S_j^2). By hand:
    (0, 1), (1, 2), (1.1, 10): W_1 = 1.5, h = 80 -> 1 - 240/4 = -59 (S44:
    Kratos' return loop never converges on it);
    curve (c): W_1 = 0.006, W_2 = 0.02925 -> 1 - 0.24/1.96 = 0.877551,
    1 - 0.2925/2.89 = 0.898789.
    Rows (1, 1.2, 1.8) at 0 C and (1, 9, 17) at 10 C on (0, 1, 2) are each
    valid (1 - 1.32/1.44 = 0.083333, 1 - 80/81 = 0.012346) but the 5 C
    curve (1, 5.1, 9.4) is not (1 - 2 * 4.3 * 3.05 / 26.01 = -0.008458): a^2
    is quadratic in T, so the curve is checked again whenever it is written.
    E03: every 304L row and the 1 C sweep -55..500 C stay valid; the
    smallest margin (independent trapezoid sum) is 0.0872, 500 C, PEEQ
    0.022 -> 0.023."""
    assert TP.point_curve_margin((0.0, 1.0, 1.1), (1.0, 2.0, 10.0)) == pytest.approx([1.0, -59.0], abs=1e-12)
    assert TP.point_curve_margin(_GRID_C, _S_C) == pytest.approx(
        [1.0, 1.0 - 0.24 / 1.96, 1.0 - 0.2925 / 2.89], abs=1e-12)
    with pytest.raises(ValueError, match="segment 1"):
        TP.HardeningTable((0.0, 1.0, 1.1), ((1.0, 2.0, 10.0),))
    ht = TP.HardeningTable((0.0, 1.0, 2.0), ((1.0, 1.2, 1.8), (1.0, 9.0, 17.0)), temperatures=(0.0, 10.0))
    assert TP.point_curve_margin(ht.plastic_strains, ht.curve(0.0))[1] == pytest.approx(1 - 1.32 / 1.44)
    assert TP.point_curve_margin(ht.plastic_strains, ht.curve(10.0))[1] == pytest.approx(1 - 80 / 81)
    assert TP.point_curve_margin(ht.plastic_strains, ht.curve(5.0))[1] == pytest.approx(
        1 - 2 * 4.3 * 3.05 / 26.01, abs=1e-12)
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _right = _bar(KM, 1, lateral_fixed=True)
    mat = TP.ThermoMaterial("m", 100.0, 0.3, 0.0, hardening=ht)
    groups = TP.CellThermalGroups(KM, mods, mp, els, [mat], 0, [0.0, 5.0])
    NL.make_newton_strategy(KM, mp).Initialize()
    with pytest.raises(ValueError, match="T = 5"):
        groups.apply(1.0)
    groups.apply(0.0)   # 0 C curve is fine
    # E03 data
    pl = np.loadtxt(_DATA / "abaqus_ss304l_plastic.tsv", delimiter="\t")
    temps = np.unique(pl[:, 2])
    grid = pl[pl[:, 2] == temps[0]][:, 1]
    rows = np.array([pl[pl[:, 2] == t][:, 0] for t in temps])
    worst = (np.inf, None, None)
    for T in np.arange(-55.0, 500.5, 1.0):
        S = np.array([np.interp(T, temps, rows[:, j]) for j in range(grid.size)])
        W = 0.0
        for j in range(grid.size - 1):
            h = (S[j + 1] - S[j]) / (grid[j + 1] - grid[j])
            m = 1.0 - 2.0 * h * W / S[j] ** 2
            worst = min(worst, (m, T, j))
            W += 0.5 * (S[j] + S[j + 1]) * (grid[j + 1] - grid[j])
        assert TP.point_curve_margin(grid, S) == pytest.approx(
            [1.0 - 2.0 * ((S[j + 1] - S[j]) / (grid[j + 1] - grid[j]))
             * np.sum(0.5 * (S[1:j + 1] + S[:j]) * np.diff(grid[:j + 1])) / S[j] ** 2
             for j in range(grid.size - 1)], abs=1e-12)
    assert worst[0] == pytest.approx(0.0872, abs=1e-4) and worst[1:] == (500.0, 10)
    mats = TP.e03_materials(_DATA)   # rows checked on construction (raise if not representable)
    assert mats["SS304L_header"].hardening.temperatures == tuple(temps)
