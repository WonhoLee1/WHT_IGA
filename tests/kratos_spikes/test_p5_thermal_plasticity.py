# -*- coding: utf-8 -*-
"""P5 / S30: thermo-plastic precursor (``analyses.thermal_plasticity``) for E03.

Thermal strain is InitialState strain per cell group, E(T) and sigma_y(T) are
refreshed per group at every level, and the thermal-step driver runs on
``nonlinear.run_load_stepping`` (one strategy for the whole history).

Model: clamped bars of unit cells [0, n] x [0, 1]^2, p = 1, one element per
cell, u_x = 0 at x = 0 and x = n, nu = 0, perfect plasticity unless stated.
Two lateral settings, chosen so the FE field is exactly 1D:

* lateral FREE (symmetry planes y = 0, z = 0 only) -> uniaxial STRESS. Only
  used with one cell / equal thermal strain in all cells: neighbouring cells
  with different lateral (thermal or plastic) strain cannot share a
  continuous u_y without interface shear, so a multi-cell bar with a
  cell-wise temperature is NOT 1D when its sides are free.
* lateral FIXED (u_y = u_z = 0 everywhere) -> uniaxial STRAIN, exactly 1D for
  any cell-wise data. With nu = 0 the axial and lateral stresses decouple:
  sigma_x = E (eps - th - p), sigma_y = sigma_z = E (-th + p/2) with the
  isochoric plastic strain (p, -p/2, -p/2); von Mises is
  s = sigma_x - sigma_y = sigma + E (th - p/2), |s| <= sigma_y.
  Equilibrium: sigma_x = sigma in every cell; compatibility: sum eps_i = 0.
  Elastic cell: eps = sigma/E + th + p_old. Plastic cell (s = c sigma_y,
  c = +-1): p = 2 (sigma/E + th - c sigma_y/E), eps = 3 (sigma/E + th) - 2 c sigma_y/E.
  Compatibility is then linear in sigma; :func:`_uniaxial_strain_reference`
  solves it per level with an active-set loop (independent of Kratos).

Kratos treats a trial state with 0 < F <= 1e-4 sigma_y as elastic
(``generic_small_strain_isotropic_plasticity.cpp``); every path is checked to
stay out of that band (:func:`_assert_clear_of_yield_band`) so the references
are unambiguous. Stresses are compared to 1e-9 (exact 1D field, linear
return), unless stated.
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

_DATA = (Path(__file__).resolve().parents[2] / "research-cf-IGA" / "user-guide" / "example-problems"
         / "glass-to-metal-seal" / "data")

_TOL_S = 1e-9      # stress, O(1) values, exact 1D field
_TOL_U = 1e-11     # displacement, O(1e-2) values


# --------------------------------------------------------------------------- helpers

def _bar(KM, nx: int, lateral_fixed: bool, element_name: str = "SmallDisplacementElement3D4N",
         ny: int = 1, nz: int = 1):
    model = KM.Model()
    mp = model.CreateModelPart("T2")
    for v in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(v)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(float(nx), 1.0, 1.0), n_cells=(nx, ny, nz), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = attach_cell_elements(KM, mp, vol, box, mp.CreateNewProperties(1),
                               element_name=element_name)
    mp.CloneTimeStep(0.0)
    for n in mp.Nodes:
        if abs(n.X0) < 1e-12 or abs(n.X0 - nx) < 1e-12:
            n.Fix(KM.DISPLACEMENT_X)
        if lateral_fixed or abs(n.Y0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Y)
        if lateral_fixed or abs(n.Z0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Z)
    return model, mp, els, vol


def _ip_stress(KM, mp, el, var):
    """(n_ip, 6) stress at all integration points of one element."""
    return np.array([list(v) for v in el.CalculateOnIntegrationPoints(var, mp.ProcessInfo)])


def _ux_at(KM, mp, x):
    vals = [n.GetSolutionStepValue(KM.DISPLACEMENT_X) for n in mp.Nodes if abs(n.X0 - x) < 1e-12]
    assert vals and np.ptp(vals) < 1e-14
    return vals[0]


def _record(KM, mp, els, var, nx):
    """Observer storing (level, per-element IP stresses, interface u_x, x=0 reaction)."""
    rec: list[dict] = []

    def obs(level):
        rec.append({
            "level": level,
            "stress": [_ip_stress(KM, mp, el, var) for el in els],
            "u": [_ux_at(KM, mp, float(x)) for x in range(nx + 1)],
            "rx0": sum(n.GetSolutionStepValue(KM.REACTION_X) for n in mp.Nodes if abs(n.X0) < 1e-12),
            "uy_top": [n.GetSolutionStepValue(KM.DISPLACEMENT_Y) for n in mp.Nodes
                       if abs(n.Y0 - 1.0) < 1e-12],
        })
    return rec, obs


def _uniaxial_strain_reference(levels_data):
    """Series bar of unit cells, uniaxial strain, nu = 0, perfect J2 (see module doc).

    ``levels_data``: per level ``(E, th, sy)`` arrays (sy = inf: elastic).
    Returns per level dict(sigma, eps, p, s, s_trial)."""
    n = len(levels_data[0][0])
    p = np.zeros(n)
    out = []
    for E, th, sy in levels_data:
        E, th, sy = (np.asarray(a, dtype=float) for a in (E, th, sy))
        s_trial = None
        act: dict[int, float] = {}
        for _ in range(50):
            A = sum(3.0 / E[i] if i in act else 1.0 / E[i] for i in range(n))
            B = sum(3.0 * th[i] - 2.0 * act[i] * sy[i] / E[i] if i in act else th[i] + p[i]
                    for i in range(n))
            sig = -B / A
            s = sig + E * (th - p / 2.0)
            if s_trial is None:
                s_trial = s.copy()
            p_new = p.copy()
            for i, c in act.items():
                p_new[i] = 2.0 * (sig / E[i] + th[i] - c * sy[i] / E[i])
            changed = False
            for i in range(n):
                if i not in act and abs(s[i]) > sy[i]:
                    act[i] = float(np.sign(s[i]))
                    changed = True
            for i, c in list(act.items()):
                if c * (p_new[i] - p[i]) < 0.0:
                    del act[i]
                    changed = True
            if not changed:
                break
        else:
            raise RuntimeError("reference active set did not settle")
        eps = np.array([3.0 * (sig / E[i] + th[i]) - 2.0 * act[i] * sy[i] / E[i] if i in act
                        else sig / E[i] + th[i] + p[i] for i in range(n)])
        p = p_new
        out.append({"sigma": sig, "eps": eps, "p": p.copy(), "s": sig + E * (th - p / 2.0),
                    "s_trial": s_trial, "sy": sy, "E": E, "th": th})
    return out


def _assert_clear_of_yield_band(ref, rel=1e-3):
    """No trial state in (1e-9, rel] * sigma_y beyond yield (Kratos elastic band 1e-4)."""
    for r in ref:
        fin = np.isfinite(r["sy"])  # elastic cells (sy = inf) cannot yield
        over = (np.abs(r["s_trial"][fin]) - r["sy"][fin]) / r["sy"][fin]
        bad = (over > 1e-9) & (over <= rel)
        assert not np.any(bad), f"trial within Kratos yield band: {over}"


def _check_uniaxial_strain(rec, ref, tol_s=_TOL_S, tol_u=_TOL_U):
    assert len(rec) == len(ref)
    for r, q in zip(rec, ref):
        for i, st in enumerate(r["stress"]):
            assert st[:, 0] == pytest.approx(q["sigma"], abs=tol_s), (r["level"], i)
            sig_y = q["E"][i] * (-q["th"][i] + q["p"][i] / 2.0)
            assert st[:, 1] == pytest.approx(sig_y, abs=tol_s), (r["level"], i)
            assert st[:, 2] == pytest.approx(sig_y, abs=tol_s), (r["level"], i)
        assert r["u"] == pytest.approx(np.concatenate([[0.0], np.cumsum(q["eps"])]), abs=tol_u)
        assert r["rx0"] == pytest.approx(-q["sigma"], abs=tol_s * 10)  # unit section


def _groups(KM, mods, mp, els, materials, mat_idx, states, finite=False):
    return TP.CellThermalGroups(KM, mods, mp, els, materials, mat_idx, states,
                                law_factory=TP.default_law_factory(mods, finite_strain=finite))


def _run(KM, mp, groups, increments, rec_obs):
    out = TP.run_thermal_steps(mp, KM, lambda: NL.make_newton_strategy(KM, mp), groups,
                               increments, on_converged=rec_obs)
    assert out["n_solves"] == len(out["levels"]), "unexpected cutback (levels would not align)"
    return out


def _cell_T(states, level):
    """Independent per-cell temperature at a level (linear between states)."""
    S = [np.asarray(s, dtype=float) for s in states]
    k = min(int(np.floor(level)), len(S) - 2)
    w = level - k
    return (1.0 - w) * S[k] + w * S[k + 1]


# --------------------------------------------------------------------------- (a)

def test_a_cellwise_temperature_uniaxial_strain_bar():
    """Six cells, T: 0 -> (0, 2, 2, 4, 4, 14), E = 100, alpha = 1e-3,
    sigma_y(T) = 1 - 0.025 T (table (0, 1), (20, 0.5)), lateral fixed, 5 increments.

    Elastic (levels 0.2-0.6): sigma = -E mean(th) (compatibility), e.g. level
    0.6: sigma = -100 * 0.0026 = -0.26; s_i = sigma + E th_i, the hottest cell
    has the largest s and the lowest sigma_y. Final level, cell 6 alone at
    s = +sigma_y(14) = 0.65, cells 1-5 elastic:
        sigma (n_e + 3) / E = -sum_e th_e - 3 th_6 + 2 sigma_y/E
        sigma = (-1.2 - 4.2 + 1.3) / 8 = -0.5125,
        p_6 = 2 (sigma/E + th_6 - sigma_y/E) = 0.00475,
    elastic cells |s| = (0.5125, 0.3125, 0.3125, 0.1125, 0.1125) < sigma_y(T_i).
    Cells with equal temperature history share one group (4 groups for 6 cells),
    and the sub model parts are reused for all levels."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _vol = _bar(KM, 6, lateral_fixed=True)
    Tf = np.array([0.0, 2.0, 2.0, 4.0, 4.0, 14.0])
    mat = TP.ThermoMaterial("m", 100.0, 0.0, 1e-3, t_zero=0.0,
                            yield_stress=TP.TempTable((0.0, 20.0), (1.0, 0.5)))
    states = [0.0, Tf]
    groups = _groups(KM, mods, mp, els, [mat], 0, states)
    assert groups.n_groups == 4
    gid = [groups.group_of_element[el.Id] for el in els]
    assert gid[1] == gid[2] and gid[3] == gid[4] and len(set(gid)) == 4
    n_smp = mp.NumberOfSubModelParts()

    rec, obs = _record(KM, mp, els, KM.CAUCHY_STRESS_VECTOR, 6)
    _run(KM, mp, groups, 5, obs)
    assert mp.NumberOfSubModelParts() == n_smp == 4

    ref = _uniaxial_strain_reference([
        (np.full(6, 100.0), 1e-3 * _cell_T(states, r["level"]),
         np.interp(_cell_T(states, r["level"]), [0.0, 20.0], [1.0, 0.5])) for r in rec])
    _assert_clear_of_yield_band(ref)
    assert ref[2]["p"].max() == 0.0 and ref[-1]["p"][5] > 0.0  # both regimes covered
    assert ref[2]["sigma"] == pytest.approx(-0.26, abs=1e-12)
    assert ref[-1]["sigma"] == pytest.approx(-0.5125, abs=1e-12)
    assert ref[-1]["p"] == pytest.approx([0, 0, 0, 0, 0, 0.00475], abs=1e-12)
    _check_uniaxial_strain(rec, ref)


# --------------------------------------------------------------------------- (b)

_SY_B = TP.TempTable((0.0, 10.0, 30.0, 60.0), (1.0, 0.8, 0.5, 0.3))
# hand-derived (docstring of test_b), sigma_xx at T = 5, 10, ..., 40, then 35, ..., 0
_SIG_B = [-0.5, -0.8, -0.725, -0.65, -0.575, -0.5, -(0.5 - 0.2 * 5 / 30), -(0.5 - 0.2 * 10 / 30),
          0.1 * 5 - (0.5 - 0.2 * 10 / 30), 0.5, 0.575, 0.65, 0.725, 0.8, 0.9, 1.0]


def test_b_yield_stress_table_heat_and_cool():
    """One clamped cell, lateral free (uniaxial stress), E = 100, alpha = 1e-3,
    sigma_y table (0, 1.0), (10, 0.8), (30, 0.5), (60, 0.3); T: 0 -> 40 -> 0,
    8 increments each (5 C).

    eps_xx = 0, so sigma = -E (alpha T + p). Heating: elastic sigma = -0.1 T
    until 0.1 T = sigma_y(T) (T = 8.33); afterwards sigma = -sigma_y(T)
    (flow stays monotone: dp/dT = -alpha + sigma_y'/E < 0):
    T = 5..40 -> -0.5, -0.8, -0.725, -0.65, -0.575, -0.5, -0.4667, -0.4333.
    At 40 C p = -(0.04 - 0.4333/E) = -0.035667. Cooling: elastic
    sigma = 0.1 (40 - T) - 0.4333 (T = 35: +0.0667) until it reaches
    +sigma_y(T) at T = 30.71; then sigma = +sigma_y(T) down to 0 C:
    T = 30..0 -> 0.5, 0.575, 0.65, 0.725, 0.8, 0.9, 1.0.
    Lateral: u_y(top) = eps_yy = alpha T - p/2."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _vol = _bar(KM, 1, lateral_fixed=False)
    mat = TP.ThermoMaterial("m", 100.0, 0.0, 1e-3, t_zero=0.0, yield_stress=_SY_B)
    groups = _groups(KM, mods, mp, els, [mat], 0, [0.0, 40.0, 0.0])
    rec, obs = _record(KM, mp, els, KM.CAUCHY_STRESS_VECTOR, 1)
    _run(KM, mp, groups, 8, obs)

    T_path = [5.0 * k for k in range(1, 9)] + [40.0 - 5.0 * k for k in range(1, 9)]
    assert [_cell_T([0.0, 40.0, 0.0], r["level"]).item() for r in rec] == pytest.approx(T_path)
    # independent incremental 1D return map (checks the hand list too)
    p, sig_ref, p_ref = 0.0, [], []
    for T in T_path:
        trial = -100.0 * (1e-3 * T + p)
        sy = float(np.interp(T, _SY_B.temperatures, _SY_B.values))
        assert not 0.0 < abs(trial) - sy <= 1e-3 * sy
        s = trial if abs(trial) <= sy else np.sign(trial) * sy
        p = -1e-3 * T - s / 100.0
        sig_ref.append(s)
        p_ref.append(p)
    assert sig_ref == pytest.approx(_SIG_B, abs=1e-12)
    for r, s, pp, T in zip(rec, sig_ref, p_ref, T_path):
        st = r["stress"][0]
        assert st[:, 0] == pytest.approx(s, abs=_TOL_S), T
        assert np.abs(st[:, 1:3]).max() < _TOL_S
        assert r["uy_top"] == pytest.approx([1e-3 * T - pp / 2.0] * len(r["uy_top"]), abs=_TOL_U)


# --------------------------------------------------------------------------- (c)

def _cycle_materials(scale: float = 1.0):
    """Metal A (plastic) + glass-like B (thermoelastic), both with E(T).
    ``scale`` multiplies E and divides alpha (same stresses, strains / scale)."""
    a = TP.ThermoMaterial("A", TP.TempTable((0.0, 20.0), (100.0 * scale, 60.0 * scale)), 0.0,
                          2e-3 / scale, t_zero=0.0,
                          yield_stress=TP.TempTable((0.0, 20.0), (1.0, 0.3)))
    b = TP.ThermoMaterial("B", TP.TempTable((0.0, 20.0), (100.0 * scale, 80.0 * scale)), 0.0,
                          0.5e-3 / scale, t_zero=0.0, yield_stress=None)
    return [a, b]


_CYCLE_STATES = [0.0, 20.0, 20.0, 0.0]
_CYCLE_INCS = [4, 2, 4]


def _cycle_reference(rec, scale: float = 1.0):
    data = []
    for r in rec:
        T = _cell_T(_CYCLE_STATES, r["level"]).item()
        E = np.array([np.interp(T, [0, 20], [100.0, 60.0]), np.interp(T, [0, 20], [100.0, 80.0])])
        data.append((E * scale, np.array([2e-3, 0.5e-3]) * T / scale,
                     np.array([np.interp(T, [0, 20], [1.0, 0.3]), np.inf])))
    return _uniaxial_strain_reference(data)


def _run_cycle(element_name, var, finite, scale):
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    var = KM.CAUCHY_STRESS_VECTOR if var is None else var
    model, mp, els, _vol = _bar(KM, 2, lateral_fixed=True, element_name=element_name)
    groups = _groups(KM, mods, mp, els, _cycle_materials(scale), [0, 1], _CYCLE_STATES,
                     finite=finite)
    assert groups.n_groups == 2
    rec, obs = _record(KM, mp, els, var, 2)
    _run(KM, mp, groups, _CYCLE_INCS, obs)
    ref = _cycle_reference(rec, scale)
    _assert_clear_of_yield_band(ref)
    return rec, ref


def test_c_three_step_cycle_two_materials_residual_stress():
    """Cells A (metal: alpha 2e-3, E 100 -> 60, sigma_y 1.0 -> 0.3 over 0-20 C)
    and B (glass-like, linear elastic law: alpha 0.5e-3, E 100 -> 80), lateral
    fixed, T: 0 -> 20 (heat, 4 inc) -> 20 (hold, 2 inc) -> 0 (cool, 4 inc).

    Heat: A yields between 10 and 15 C (s_A = +sigma_y); at 20 C with
    E_A = 60, E_B = 80, th = (0.04, 0.01):
        eps_A = 3 (sigma/60 + 0.04) - 2 * 0.3/60 = sigma/20 + 0.11,
        eps_B = sigma/80 + 0.01,  eps_A + eps_B = 0
        -> sigma = -0.12 / 0.0625 = -1.92, p_A = 0.006, u(1) = 0.014.
    Hold at the high-temperature properties: no change (no drift, no double
    thermal strain). Cool (A elastic, p_A frozen; E back to 100, th = 0):
        sigma (1/100 + 1/100) = -p_A -> residual sigma = -0.3,
        sigma_y,A = E_A p_A / 2 = +0.3, sigma_y,B = 0, u(1) = eps_A = 0.003,
        s_A = -0.6 inside sigma_y(0) = 1.
    With E frozen at 100 the heated state would differ (E(T) enters p_A)."""
    rec, ref = _run_cycle("SmallDisplacementElement3D4N", None, False, 1.0)
    i_heat, i_hold, i_end = 3, 5, len(ref) - 1
    assert ref[i_heat]["sigma"] == pytest.approx(-1.92, abs=1e-12)
    assert ref[i_heat]["p"] == pytest.approx([0.006, 0.0], abs=1e-12)
    assert ref[i_end]["sigma"] == pytest.approx(-0.3, abs=1e-12)
    assert ref[i_end]["eps"][0] == pytest.approx(0.003, abs=1e-12)
    _check_uniaxial_strain(rec, ref)
    # hold step: state identical to the end of heating
    for k in (4, i_hold):
        for st_h, st_0 in zip(rec[k]["stress"], rec[i_heat]["stress"]):
            assert st_h == pytest.approx(st_0, abs=_TOL_S)
    # residual stress state, written out
    sA, sB = rec[i_end]["stress"]
    assert sA[:, :3] == pytest.approx(np.tile([-0.3, 0.3, 0.3], (len(sA), 1)), abs=_TOL_S)
    assert sB[:, :3] == pytest.approx(np.tile([-0.3, 0.0, 0.0], (len(sB), 1)), abs=_TOL_S)


# --------------------------------------------------------------------------- (d)

_E_D_PTS = ((0.0, 50.0, 100.0), (100.0, 50.0, 40.0))
_E_D = TP.TempTable(*_E_D_PTS)


@pytest.mark.parametrize("plastic_law", [False, True], ids=["linear_elastic_law", "j2_law_elastic"])
def test_d_temperature_dependent_young_modulus(plastic_law):
    """E(T) table (0, 100), (50, 50), (100, 40), alpha = 1e-3, sigma_y = 1e6
    (J2 law kept elastic) or the linear elastic law (glass path).

    (1) one clamped cell, lateral free, T: 0 -> 50 -> 100 (2 inc each):
        sigma = -E(T) alpha T = -1.875, -2.5, -3.375, -4.0 at T = 25..100
        (frozen E = 100 would give -2.5, -5, -7.5, -10).
    (2) two cells, lateral fixed, T: 0 -> (50, 100) (2 inc): per-group E,
        sigma = -(th_A + th_B) / (1/E_A + 1/E_B) = -0.15/0.045 = -3.3333 and
        sigma_y,i = -E_i th_i = (-2.5, -4.0) at the end."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mat = TP.ThermoMaterial("m", _E_D, 0.0, 1e-3, t_zero=0.0,
                            yield_stress=1e6 if plastic_law else None)

    model, mp, els, _vol = _bar(KM, 1, lateral_fixed=False)
    groups = _groups(KM, mods, mp, els, [mat], 0, [0.0, 50.0, 100.0])
    rec, obs = _record(KM, mp, els, KM.CAUCHY_STRESS_VECTOR, 1)
    _run(KM, mp, groups, 2, obs)
    got = [r["stress"][0][:, 0] for r in rec]
    for g, ref in zip(got, (-1.875, -2.5, -3.375, -4.0)):
        assert g == pytest.approx(ref, abs=_TOL_S)

    model2, mp2, els2, _vol2 = _bar(KM, 2, lateral_fixed=True)
    states = [0.0, np.array([50.0, 100.0])]
    groups2 = _groups(KM, mods, mp2, els2, [mat], 0, states)
    rec2, obs2 = _record(KM, mp2, els2, KM.CAUCHY_STRESS_VECTOR, 2)
    _run(KM, mp2, groups2, 2, obs2)
    ref2 = _uniaxial_strain_reference([
        (np.interp(_cell_T(states, r["level"]), *_E_D_PTS), 1e-3 * _cell_T(states, r["level"]),
         np.full(2, np.inf))
        for r in rec2])
    assert ref2[-1]["sigma"] == pytest.approx(-0.15 / 0.045, abs=1e-12)
    _check_uniaxial_strain(rec2, ref2)
    sA, sB = rec2[-1]["stress"]
    assert sA[:, 1] == pytest.approx(-2.5, abs=_TOL_S) and sB[:, 1] == pytest.approx(-4.0, abs=_TOL_S)


# --------------------------------------------------------------------------- (e)

def test_e_total_lagrangian_finite_strain_matches_small_strain():
    """The (c) cycle with TotalLagrangianElement3D4N, the finite-strain J2 law
    (FiniteStrainIsotropicPlasticity3DVonMisesVonMises, Hencky strain, yield on
    Kirchhoff stress) for A and the linear elastic law (Green-Lagrange) for B.
    E x 1e4 and alpha / 1e4: the small-strain stresses are unchanged and
    |eps| <= eps_max = 4e-6 (th_A at 20 C).

    PK2 vs small strain differs only by geometric terms: 1D TL equilibrium
    holds for P = tau / (1 + eps) (not for S = tau / (1 + eps)^2), and
    ln(1 + eps) - eps, GL - eps = O(eps^2). So
    |S - sigma| <= 10 eps_max sigma_max (generous constant; the leading term
    is 2 eps sigma); displacements to that stress tolerance / E_min. No yield
    band term: ``_assert_clear_of_yield_band`` keeps every trial > 1e-3
    sigma_y away from the 1e-4 band and the perfectly plastic J2 return is
    exact in one iteration."""
    scale = 1e4
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    rec, ref = _run_cycle("TotalLagrangianElement3D4N", KM.PK2_STRESS_VECTOR, True, scale)
    eps_max = max(max(np.abs(q["eps"]).max(), np.abs(q["th"]).max(), np.abs(q["p"]).max())
                  for q in ref)
    sig_max = max(abs(q["sigma"]) for q in ref)
    assert eps_max < 5e-6
    tol_s = 10.0 * eps_max * sig_max
    _check_uniaxial_strain(rec, ref, tol_s=tol_s, tol_u=2 * tol_s / (60.0 * scale))
    assert ref[3]["sigma"] == pytest.approx(-1.92, abs=1e-12)  # same plastic path as (c)
    assert ref[-1]["sigma"] == pytest.approx(-0.3, abs=1e-12)


# --------------------------------------------------------------------------- data + scale

def test_f_e03_secant_expansion_matches_recommended_thermal_strain():
    """E03 data form (Abaqus *EXPANSION, ZERO = 20 C): the secant tables give
    the Coreform 'recommended thermal strain relative to 20 C' CSV (an
    independent column-wise derivation of the same supplier data) at all 12
    temperatures, for all three materials; table values outside the range are
    held constant (304L E at -55 C = E(27 C) = 193000); the 304L -55 C initial
    yield is the 20 C value x 1.2 (doc statement); with the 450 C stress-free
    initial state, eps_th(20 C) = -eps_rel20(450 C)."""
    mats = TP.e03_materials_initial_yield(_DATA)
    csv = np.genfromtxt(_DATA / "thermal_strain_recommended.csv", delimiter=",", names=True)
    cols = {"S8061_glass": "s8061_strain_rel_20C", "Alloy52_pins": "alloy52_strain_rel_20C",
            "SS304L_header": "ss304_strain_rel_20C"}
    T = csv["temperature_C"]
    for name, col in cols.items():
        got = mats[name].thermal_strain(T, 20.0)
        assert got == pytest.approx(csv[col], abs=1e-9), name   # CSV printed to 1e-9
        assert mats[name].thermal_strain(450.0, 450.0) == pytest.approx(0.0, abs=1e-18)
        i450 = int(np.flatnonzero(T == 450.0)[0])
        assert mats[name].thermal_strain(20.0, 450.0) == pytest.approx(-csv[col][i450], abs=1e-9)
    ss = mats["SS304L_header"]
    assert ss.young(-55.0) == 193000.0 and ss.young(600.0) == 155000.0
    assert ss.yield_at(-55.0) == pytest.approx(1.2 * ss.yield_at(20.0), abs=1e-9)
    assert mats["S8061_glass"].young(25.0) == 69000.0 and not mats["S8061_glass"].plastic


def test_g_grouping_many_cells_sets_cellwise_initial_strain():
    """2000 cells (20 x 10 x 10), linear elastic, alpha = 1e-3, T_zero = 0,
    random cell temperatures on a 0.05 C grid in [0, 100] plus jitter in
    +-4e-7 C (< quantum/2 = 5e-7: every cell stays in the rounding bin of its
    grid value), two states (T_init = 20 uniform).

    Reference group count: number of distinct rounded keys (Python set, not
    numpy.unique). After ``apply(1.0)`` every element's INITIAL_STRAIN_VECTOR
    (read back from its constitutive law) equals alpha (T_cell - 20) to
    alpha * 8e-7 = 8e-10 (group T = member mean, members within 8e-7 C) in
    the three normal components, 0 in shear."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _vol = _bar(KM, 20, lateral_fixed=True, ny=10, nz=10)
    assert len(els) == 2000
    rng = np.random.default_rng(30)
    T_round = 0.05 * rng.integers(0, 2001, size=2000)
    T_cell = T_round + rng.uniform(-4e-7, 4e-7, size=2000)
    mat = TP.ThermoMaterial("glass", 100.0, 0.0, 1e-3, t_zero=0.0, yield_stress=None)
    groups = _groups(KM, mods, mp, els, [mat], 0, [20.0, T_cell])
    assert groups.n_groups == len({round(t / 1e-6) for t in T_cell}) == len(set(T_round))
    strat = NL.make_newton_strategy(KM, mp)
    strat.Initialize()
    groups.apply(1.0)
    for el, t in zip(els, T_cell):
        eps0 = np.array([list(v) for v in el.CalculateOnIntegrationPoints(
            KM.INITIAL_STRAIN_VECTOR, mp.ProcessInfo)])
        assert eps0[:, :3] == pytest.approx(1e-3 * (t - 20.0), abs=1e-9)
        assert np.abs(eps0[:, 3:]).max() == 0.0


def test_h_intermediate_state_separates_groups():
    """Two cells share initial (T=0) and final (T=20) states and material, but
    differ at intermediate state 1 (T=(10, 0)), lateral fixed, E=100, alpha=1e-3,
    sigma_y=0.3, 2 increments per step.

    Grouping must inspect the full temperature history (not only hist.T[:, -1:]):
    two groups are built, and at state 1 (level 1.0) both cells yield in
    opposite directions (sigma=-0.5, p=(+0.004, -0.004), sigma_yy=(-0.8, -0.2),
    u(1)=0.009). At the shared final state T=20 (level 2.0), retained plastic
    strains keep sigma_yy=(-1.8, -2.2) and u(1)=0.004 distinct."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _vol = _bar(KM, 2, lateral_fixed=True)
    mat = TP.ThermoMaterial("m", 100.0, 0.0, 1e-3, t_zero=0.0, yield_stress=0.3)
    states = [0.0, np.array([10.0, 0.0]), 20.0]
    groups = _groups(KM, mods, mp, els, [mat], 0, states)
    assert groups.n_groups == 2
    assert groups.group_of_element[els[0].Id] != groups.group_of_element[els[1].Id]

    rec, obs = _record(KM, mp, els, KM.CAUCHY_STRESS_VECTOR, 2)
    _run(KM, mp, groups, 2, obs)
    ref = _uniaxial_strain_reference([
        (np.full(2, 100.0), 1e-3 * _cell_T(states, r["level"]), np.full(2, 0.3)) for r in rec])
    _assert_clear_of_yield_band(ref)
    _check_uniaxial_strain(rec, ref)
    assert rec[1]["level"] == pytest.approx(1.0)
    assert rec[1]["stress"][0][:, 1] == pytest.approx(-0.8, abs=_TOL_S)
    assert rec[1]["stress"][1][:, 1] == pytest.approx(-0.2, abs=_TOL_S)
    assert rec[1]["u"][1] == pytest.approx(0.009, abs=_TOL_U)
    assert rec[1]["rx0"] == pytest.approx(0.5, abs=_TOL_S * 10)
    assert rec[-1]["stress"][0][:, 1] == pytest.approx(-1.8, abs=_TOL_S)
    assert rec[-1]["stress"][1][:, 1] == pytest.approx(-2.2, abs=_TOL_S)
    assert rec[-1]["u"][1] == pytest.approx(0.004, abs=_TOL_U)


def test_i_group_temperature_error_below_one_quantum():
    """Documented grouping bound: quantum = 1, one cell at 99.51 C and 20 at
    100.49 C round to the same bin (100) -> one group at the member MEAN
    (99.51 + 20 * 100.49) / 21 = 100.443333; the lone cell is off by
    0.933333 quantum: more than quantum/2 (bin centre would give <= 0.49),
    less than one quantum (bin width)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    model, mp, els, _vol = _bar(KM, 21, lateral_fixed=True)
    T = np.array([99.51] + [100.49] * 20)
    mat = TP.ThermoMaterial("glass", 100.0, 0.0, 1e-3, t_zero=0.0, yield_stress=None)
    groups = TP.CellThermalGroups(KM, mods, mp, els, [mat], 0, [20.0, T], quantum=1.0)
    assert groups.n_groups == 1
    assert groups.group_T[1, 0] == pytest.approx((99.51 + 20 * 100.49) / 21.0, abs=1e-12)
    err = np.abs(T - groups.group_T[1, 0]).max()
    assert err == pytest.approx(0.98 * 20.0 / 21.0, abs=1e-12)
    assert 0.5 < err < 1.0
