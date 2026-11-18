# -*- coding: utf-8 -*-
"""S24-hybrid spike: Kratos MechanicalExplicitStrategy with Python HRZ overwrite.

Verifies the S24-correction hybrid path: after ``Initialize()``, master
``NODAL_MASS`` is overwritten with HRZ(T^T M_c T) values and the Kratos
explicit strategy is used unchanged (no Python time loop).

Gates (evidence-calibrated, probe 2026-11-16):
- HRZ reduced masses strictly positive for p=1/2 x eta=1/1e-2/1e-4/1e-6.
- Overwritten masses persist across Solve (Kratos does not refresh them).
- Wave accuracy: untrimmed tip L2 < 2% (S24 Q2 Python loop: 0.85%);
  trimmed+T tip L2 < 2% (Python loop: 1.16%; Kratos MPC path: 1.2%).
- No blow-up over the full run.

Review 2026-11-17 (dev_log/track2_review_20261116_opus.md, H2): the trimmed
case used to probe the eliminated slave control point (extrapolated to
x = 10, not the physical tip x = 9.75) against u = s*c/E*t, which is only
valid until the reflection returns at 2*L_trim/c = 1.95 s < t_end. That gave
4.36 % and a 5 % gate. It now probes the field at the physical tip and uses
the fixed-free reference with the reflection.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.kratos_spikes.test_s24_explicit_dynamics import (
    _assemble_trimmed_consistent_mass,
)
from wht_iga.dynamics import explicit as _EX
from wht_iga.dynamics.explicit import build_x_extension as _bxe
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference import assemble_elastic as ref
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_L, _E, _NU, _RHO = 10.0, 100.0, 0.0, 1.0
_C = float(np.sqrt(_E / _RHO))
_SIGMA0 = 1.0


def _explicit_model(KM, SM, name: str):
    model = KM.Model()
    mp = model.CreateModelPart(name)
    for v in (KM.DISPLACEMENT, KM.VELOCITY, KM.ACCELERATION, KM.REACTION,
              SM.MIDDLE_VELOCITY, KM.NODAL_MASS, KM.FORCE_RESIDUAL,
              KM.RESIDUAL_VECTOR, SM.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(v)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.COMPUTE_LUMPED_MASS_MATRIX, True)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_ALPHA, 0.0)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_BETA, 0.0)
    return model, mp


def _build_wave_model(KM, SM, box, eta_last, use_ties, name):
    mp_holder = {}
    model, mp = _explicit_model(KM, SM, name)
    vol, _ = make_nurbs_volume(KM, mp, box)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, _E)
    props.SetValue(KM.POISSON_RATIO, _NU)
    props.SetValue(KM.DENSITY, _RHO)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    nx = box.n_cells[0]
    if eta_last == 1.0 and not use_ties:
        attach_cell_elements(KM, mp, vol, box, props)
    else:
        from wht_iga.kratos_backend.geometry import cell_gauss_points
        eid = 1
        for k in range(box.n_cells[2]):
            for j in range(box.n_cells[1]):
                for i in range(nx):
                    pts = cell_gauss_points(box, i, j, k)
                    if eta_last < 1.0 and i == nx - 1:
                        u0 = float(box.breakpoints(0)[i])
                        pts = [[u0 + eta_last * (x - u0), y, z, eta_last * w]
                               for x, y, z, w in pts]
                    qps = KM.GeometriesVector()
                    vol.CreateQuadraturePointGeometries(qps, 2, pts)
                    mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                    eid += 1
    nu_cp = box.n_cp[0]
    assert mp.NumberOfNodes() == box.n_ctrl
    for nd in mp.Nodes:
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
    for I in range(box.n_ctrl):
        if I % nu_cp == 0:
            mp.GetNode(I + 1).Fix(KM.DISPLACEMENT_X)
    mp_holder["mp"] = mp
    mp_holder["model"] = model
    return mp_holder


def _apply_tip_traction(KM, SM, mp, box, nx, eta, use_ties, degree):
    nu_cp = box.n_cp[0]
    if eta == 1.0 and not use_ties:
        f_ext = ref.assemble_face_traction(box, axis=0, side=1,
                                           traction=[_SIGMA0, 0.0, 0.0])
        tip_f = {I: float(f_ext[3 * I]) for I in range(box.n_ctrl)
                 if I % nu_cp == nu_cp - 1}
    else:
        xi_tip = np.array([-1.0 + 2.0 * eta])
        fu_t, Nu_t, _, _, _ = box.axis_tables(0, xi_tip)
        w_x_tip = np.zeros(box.n_cp[0])
        w_x_tip[fu_t[nx - 1]:fu_t[nx - 1] + degree + 1] = Nu_t[nx - 1, 0, :]
        f_unit = ref.assemble_face_traction(box, axis=0, side=1,
                                            traction=[1.0, 0.0, 0.0])
        yz_w = f_unit[3 * (box.n_cp[0] - 1)::3 * box.n_cp[0]]
        tip_f = {}
        for kc in range(box.n_cp[2]):
            for jc in range(box.n_cp[1]):
                yz = jc + box.n_cp[1] * kc
                for ic in range(box.n_cp[0]):
                    I = ic + box.n_cp[0] * yz
                    tip_f[I] = _SIGMA0 * float(w_x_tip[ic]) * float(yz_w[yz])
    props2 = mp.GetProperties()[7]
    lcid = 200000
    for I, f in tip_f.items():
        if abs(f) < 1e-300:
            continue
        mp.CreateNewCondition("PointLoadCondition3D1N", lcid, [I + 1], props2)
        mp.GetNode(I + 1).SetSolutionStepValue(SM.POINT_LOAD_X, float(f))
        lcid += 1


def _run_kratos_wave(KM, SM, mp, box, nx, dt, t_end):
    nu_cp = box.n_cp[0]
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    tip_ids = [I + 1 for I in range(box.n_ctrl) if I % nu_cp == nu_cp - 1]
    n_steps = int(round(t_end / dt))
    hist = []
    ok = True
    t_abs = 0.0
    for _s in range(n_steps):
        t_abs += dt
        mp.CloneTimeStep(t_abs)
        mp.ProcessInfo[KM.STEP] = _s + 1
        strat.Solve()
        u = float(np.mean([mp.GetNode(i).GetSolutionStepValue(KM.DISPLACEMENT_X)
                           for i in tip_ids]))
        hist.append(u)
        if not np.isfinite(u) or abs(u) > 5.0:
            ok = False
            break
    hist = np.array(hist)
    t_hist = (np.arange(len(hist)) + 1) * dt
    exact = (_SIGMA0 * _C / _E) * t_hist
    rel = float(np.linalg.norm(hist - exact) / np.linalg.norm(exact))
    return ok, rel, len(hist)


def test_hybrid_hrz_masses_positive():
    """HRZ(T^T M_c T) diagonal strictly positive for p=1/2 x eta sweep."""
    from wht_iga.dynamics.explicit import build_x_extension, hrz_lump_diagonal

    for degree, nx in ((1, 4), (2, 4)):
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(float(nx), 1.0, 1.0),
                        n_cells=(nx, 1, 1), degree=degree)
        info = build_x_extension(box)
        T = info["T"]
        for eta in (1.0, 1e-2, 1e-4, 1e-6):
            M_c = _assemble_trimmed_consistent_mass(box, rho=1.0, eta_last_x=eta)
            m = hrz_lump_diagonal(T.T @ M_c @ T, total_mass=float(nx - 1 + eta))
            assert np.all(m > 0.0), f"non-positive HRZ mass p={degree} eta={eta}"


def test_hybrid_kratos_wave_untrimmed():
    """Kratos explicit on untrimmed bar (Kratos lumped mass, no overwrite): tip L2 < 2%."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    nx, degree, dt, t_end = 40, 2, 0.005, 2.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(nx, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, 1.0, False, "HYB_A")
    mp = mh["mp"]
    _apply_tip_traction(KM, SM, mp, box, nx, 1.0, False, degree)
    ok, rel, n = _run_kratos_wave(KM, SM, mp, box, nx, dt, t_end)
    assert ok and n == int(round(t_end / dt))
    assert rel < 0.02, f"untrimmed hybrid tip L2 too large: {rel}"


@pytest.mark.parametrize("eta", [1e-4, 1e-6])
def test_hybrid_kratos_wave_trimmed_ties(eta: float):
    """Kratos explicit + T ties + HRZ overwrite on trimmed bar: stable, < 2%.

    Overwritten master masses persist across Solve (Kratos does not refresh
    NODAL_MASS every step); trimmed tip stays stable at dt = 0.005 (CFL 0.2).
    Probe: displacement field at the physical tip x = (nx - 1 + eta) h.
    Reference: fixed-free bar under a step end load, u_tip = v0 t until
    2 L_trim / c, then v0 (2 L_trim / c) - v0 (t - 2 L_trim / c).
    """
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    nx, degree, dt, t_end = 40, 2, 0.005, 2.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(nx, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, eta, True, f"HYB_T{eta:.0e}")
    mp = mh["mp"]
    nu_cp = box.n_cp[0]
    info = _bxe(box)
    T = info["T"]
    src = [int(s) for s in info["src"]]
    wmap = dict(zip(src, [float(v) for v in info["masters"][0][1]]))
    cid = 100000
    for kc in range(box.n_cp[2]):
        for jc in range(box.n_cp[1]):
            s_id = (nu_cp - 1) + nu_cp * (jc + box.n_cp[1] * kc) + 1
            masters = [mp.GetNode(s_ + nu_cp * (jc + box.n_cp[1] * kc) + 1).GetDof(KM.DISPLACEMENT_X)
                       for s_ in src]
            add_linear_tie(KM, mp, cid,
                           [mp.GetNode(s_id).GetDof(KM.DISPLACEMENT_X)],
                           masters, np.array([[wmap[s_] for s_ in src]]))
            cid += 1
    M_c = _assemble_trimmed_consistent_mass(box, rho=_RHO, eta_last_x=eta)
    M_red = T.T @ M_c @ T
    L_trim = (nx - 1 + eta) * (_L / nx)
    hrz = _EX.hrz_lump_diagonal(M_red, total_mass=_RHO * L_trim)
    lex = [I for I in range(box.n_ctrl) if I % nu_cp != nu_cp - 1]
    assert len(lex) == T.shape[1] // 3
    ov = {lex[Im] + 1: float(hrz[3 * Im]) for Im in range(len(lex))
          if lex[Im] % nu_cp != 0}
    assert all(v > 0.0 for v in ov.values())
    _apply_tip_traction(KM, SM, mp, box, nx, eta, True, degree)
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    _EX.overwrite_nodal_mass(KM, mp, ov)
    before = _EX.read_nodal_mass(KM, mp, sorted(ov.keys()))
    mp.CloneTimeStep(dt)
    strat.Solve()
    after = _EX.read_nodal_mass(KM, mp, sorted(ov.keys()))
    assert bool(np.array_equal(before, after)), "Kratos refreshed NODAL_MASS"
    # main loop with the SAME strategy object (already stepped once above)
    # physical-tip field weights: x basis at the trimmed end, every (y, z) row
    fu_t, Nu_t, _, _, _ = box.axis_tables(0, np.array([-1.0 + 2.0 * eta]))
    w_x_tip = np.zeros(nu_cp)
    w_x_tip[fu_t[nx - 1]:fu_t[nx - 1] + degree + 1] = Nu_t[nx - 1, 0, :]
    rows = [nu_cp * (jc + box.n_cp[1] * kc)
            for kc in range(box.n_cp[2]) for jc in range(box.n_cp[1])]
    nodes = [mp.GetNode(I + 1) for I in range(box.n_ctrl)]
    n_steps = int(round(t_end / dt)) - 1
    hist = []
    ok = True
    t_abs = dt
    for _s in range(n_steps):
        t_abs += dt
        mp.CloneTimeStep(t_abs)
        mp.ProcessInfo[KM.STEP] = _s + 2
        strat.Solve()
        ux = np.array([nd.GetSolutionStepValue(KM.DISPLACEMENT_X) for nd in nodes])
        u = float(np.mean([w_x_tip @ ux[r:r + nu_cp] for r in rows]))
        hist.append(u)
        if not np.isfinite(u) or abs(u) > 5.0:
            ok = False
            break
    hist = np.array(hist)
    t_hist = (np.arange(len(hist)) + 2) * dt
    v0 = _SIGMA0 * _C / _E
    t_ref = 2.0 * L_trim / _C  # reflection back at the free end (1.95 s)
    exact = np.where(t_hist <= t_ref, v0 * t_hist, v0 * t_ref - v0 * (t_hist - t_ref))
    rel = float(np.linalg.norm(hist - exact) / np.linalg.norm(exact))
    assert ok and len(hist) == n_steps
    assert rel < 0.02, f"trimmed hybrid tip L2 too large (eta={eta}): {rel}"
