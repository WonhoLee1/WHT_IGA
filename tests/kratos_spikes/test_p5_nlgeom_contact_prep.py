# -*- coding: utf-8 -*-
"""Track 3 spike A: NLGEOM + Neo-Hookean + J2 via analyses.nonlinear (prep).

Gates (evidence: probe 2026-11-16, all machine precision):
- NH uniaxial (UL + TL, 1x1x1 p=1, lam=1.2): Cauchy sigma_11 vs analytic < 1e-12.
- J2 perfect plasticity (SmallDisplacement, 2% strain): plateau at yield < 1e-9.
- J2 linear hardening: tangent slope E*H/(E+H) within 2%.
- Load stepping (3-step ramp): all converged, final stress = analytic.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga import analyses  # noqa: F401  (package import check)
from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_E_NH, _NU_NH, _LAM1 = 3.0, 0.3, 1.2
_E_J2, _NU_J2, _YS = 100.0, 0.0, 0.5


def _build_uniaxial(KM, SM, box, law, element_name, lam1, props_setter):
    model = KM.Model()
    mp = model.CreateModelPart("P5")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props_setter(props)
    props.SetValue(KM.CONSTITUTIVE_LAW, law)
    els = attach_cell_elements(KM, mp, vol, box, props, element_name=element_name)
    mp.CloneTimeStep(0.0)
    mp.CloneTimeStep(1.0)
    for node in mp.Nodes:
        if abs(node.X0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        if abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, lam1 - 1.0)
        if abs(node.Y0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        if abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
    return model, mp, els


def _cauchy_xx(KM, mp, els):
    return float(np.mean([el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR,
                                                          mp.ProcessInfo)[0][0]
                          for el in els]))


def _solve(KM, mp):
    return NL.solve_static(mp, NL.make_newton_strategy(KM, mp))


@pytest.mark.parametrize("element_name", ["UpdatedLagrangianElement3D4N",
                                          "TotalLagrangianElement3D4N"])
def test_p5_neohookean_uniaxial_matches_analytic(element_name: str):
    """NH uniaxial via Alpha driver matches S19 analytic (< 1e-12)."""
    from tests.kratos_spikes.test_s19_constitutive_models import (
        _analytic_uniaxial_neohookean,
    )
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    _, sig_ref, _ = _analytic_uniaxial_neohookean(_E_NH, _NU_NH, _LAM1)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(1, 1, 1), degree=1)
    _, mp, els = _build_uniaxial(
        KM, mods["StructuralMechanicsApplication"], box,
        NL.hyperelastic_neohookean_law(CL), element_name, _LAM1,
        lambda p: (p.SetValue(KM.YOUNG_MODULUS, _E_NH),
                   p.SetValue(KM.POISSON_RATIO, _NU_NH)))
    assert _solve(KM, mp)
    sig = _cauchy_xx(KM, mp, els)
    assert abs(sig - sig_ref) / abs(sig_ref) < 1e-12, (sig, sig_ref)


def test_p5_j2_perfect_plasticity_plateau():
    """J2 perfect plasticity caps uniaxial stress at yield (< 1e-9)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(1, 1, 1), degree=1)
    _, mp, els = _build_uniaxial(
        KM, mods["StructuralMechanicsApplication"], box, NL.j2_law(CL),
        "SmallDisplacementElement3D4N", 1.02,
        lambda p: (p.SetValue(KM.YOUNG_MODULUS, _E_J2),
                   p.SetValue(KM.POISSON_RATIO, _NU_J2),
                   NL.set_j2_properties(p, KM, _YS, 0.0)))
    assert _solve(KM, mp)
    sig = _cauchy_xx(KM, mp, els)
    assert abs(sig - _YS) / _YS < 1e-9, (sig, _YS)


def test_p5_j2_linear_hardening_tangent():
    """J2 linear hardening: tangent slope E*H/(E+H) within 2% (two levels)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    E, H, ys = 100.0, 100.0, 0.5
    sigs = []
    for lam1 in (1.015, 1.025):
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                        n_cells=(1, 1, 1), degree=1)
        _, mp, els = _build_uniaxial(
            KM, SM, box, NL.j2_law(CL), "SmallDisplacementElement3D4N", lam1,
            lambda p: (p.SetValue(KM.YOUNG_MODULUS, E),
                       p.SetValue(KM.POISSON_RATIO, 0.0),
                       NL.set_j2_properties(p, KM, ys, H)))
        assert _solve(KM, mp)
        sigs.append(_cauchy_xx(KM, mp, els))
    tangent = (sigs[1] - sigs[0]) / (1.025 - 1.015)
    tangent_ref = E * H / (E + H)
    assert abs(tangent - tangent_ref) / tangent_ref < 0.02, (tangent, tangent_ref)


def test_p5_tied_surface_layer_build():
    """tied_surface_layer builds coincident, outward-oriented, tied triangles."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    from wht_iga.analyses import contact as _C
    from wht_iga.kratos_backend.elements import attach_cell_elements
    from wht_iga.kratos_backend.geometry import make_nurbs_volume

    model = KM.Model()
    mp = model.CreateModelPart("P5C")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(1, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 100.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    attach_cell_elements(KM, mp, vol, box, props)
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    cprop = mp.GetProperties()[2]
    lay = _C.tied_surface_layer(KM, SM, mp, box, axis=2, side=1,
                                node_id_start=1000, cond_id_start=100,
                                tie_id_start=5000, props=cprop)
    assert len(lay["nodes"]) == 4
    assert mp.NumberOfConditions() == 2
    # coincidence with spline face nodes
    for nd in lay["nodes"]:
        assert any(abs(n.X0 - nd.X0) < 1e-12 and abs(n.Y0 - nd.Y0) < 1e-12
                   and abs(n.Z0 - 1.0) < 1e-12
                   for n in mp.Nodes if n.Id < 100)
    # outward orientation (+z): signed area z-component positive
    import numpy as np
    for cond in lay["sub_model_part"].Conditions:
        p = [np.array([nd.X0, nd.Y0, nd.Z0]) for nd in cond.GetNodes()]
        n = np.cross(p[1] - p[0], p[2] - p[0])
        assert n[2] > 0.0, n


def test_p5_load_stepping_ramp():
    """run_load_stepping: 2-step NH ramp converges, final = analytic."""
    from tests.kratos_spikes.test_s19_constitutive_models import (
        _analytic_uniaxial_neohookean,
    )
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    E, nu, lam_full = 3.0, 0.3, 1.2
    _, sig_ref, _ = _analytic_uniaxial_neohookean(E, nu, lam_full)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(1, 1, 1), degree=1)
    _, mp, els = _build_uniaxial(
        KM, SM, box, NL.hyperelastic_neohookean_law(CL),
        "TotalLagrangianElement3D4N", 1.0,
        lambda p: (p.SetValue(KM.YOUNG_MODULUS, E),
                   p.SetValue(KM.POISSON_RATIO, nu)))
    tip = [n for n in mp.Nodes if abs(n.X0 - 1.0) < 1e-12]
    out = NL.run_load_stepping(
        mp, KM, lambda: NL.make_newton_strategy(KM, mp), [0.5, 1.0],
        lambda f: [n.SetSolutionStepValue(KM.DISPLACEMENT_X, (lam_full - 1.0) * f)
                   for n in tip])
    assert out["steps"] == 2
    assert all(h["converged"] for h in out["history"])
    sig = _cauchy_xx(KM, mp, els)
    assert abs(sig - sig_ref) / abs(sig_ref) < 1e-9, (sig, sig_ref)


def test_p5_bimaterial_thermal_clamped():
    """E03 reduced: 2-material clamped bar, uniform dT -> analytic stress (< 1e-9)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    from wht_iga.kratos_backend.elements import attach_cell_elements  # noqa (parity)
    from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume

    E_A, a_A = 2e5, 1e-5
    E_B, a_B = 1e5, 0.5e-5
    dT = -100.0
    sig_exact = -(a_A + a_B) * dT / (1.0 / E_A + 1.0 / E_B)

    model = KM.Model()
    mp = model.CreateModelPart("P5E")
    for v in (KM.DISPLACEMENT, KM.REACTION, KM.TEMPERATURE):
        mp.AddNodalSolutionStepVariable(v)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0),
                    n_cells=(2, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = []
    for i, (Ee, aa) in enumerate(((E_A, a_A), (E_B, a_B))):
        pr = mp.GetProperties()[i + 1]
        pr.SetValue(KM.YOUNG_MODULUS, Ee)
        pr.SetValue(KM.POISSON_RATIO, 0.0)
        pr.SetValue(KM.THERMAL_EXPANSION_COEFFICIENT, aa)
        pr.SetValue(KM.REFERENCE_TEMPERATURE, 0.0)
        pr.SetValue(KM.CONSTITUTIVE_LAW, CL.ThermalElasticIsotropic3D())
        qps = KM.GeometriesVector()
        vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, 0, 0))
        els.append(mp.CreateNewElement("SmallDisplacementElement3D4N", i + 1, qps[0], pr))
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    for nd in mp.Nodes:
        nd.SetSolutionStepValue(KM.TEMPERATURE, dT)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
        if abs(nd.X0) < 1e-12 or abs(nd.X0 - 2.0) < 1e-12:
            nd.Fix(KM.DISPLACEMENT_X)
    mp.CloneTimeStep(0.0)
    mp.CloneTimeStep(1.0)
    assert NL.solve_static(mp, NL.make_newton_strategy(KM, mp))
    sig = np.mean([el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR,
                                                   mp.ProcessInfo)[0][0]
                   for el in els])
    assert abs(sig - sig_exact) / abs(sig_exact) < 1e-9, (sig, sig_exact)


def test_p5_hertz_ssn_coarse():
    """Hertz SSN (coarse NX=12): converges, F_N within 35% of analytic quarter."""
    from wht_iga.kratos_backend.elements import (
        equation_permutation,
        solve_linear,
        system_matrix_csr,
    )
    from wht_iga.kratos_backend.geometry import make_nurbs_volume

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    E, nu = 1000.0, 0.3
    R_sph, delta, k_pen, NX, NZ = 0.5, 0.05, 1e4, 12, 2
    F_ref = (4.0 / 3.0) * (E / (1.0 - nu * nu)) * np.sqrt(R_sph) * delta ** 1.5 / 4.0

    model = KM.Model()
    mp = model.CreateModelPart("P5H")
    for v in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(v)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, -1.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(NX, NX, 2 * NZ), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    from wht_iga.kratos_backend.elements import attach_cell_elements
    attach_cell_elements(KM, mp, vol, box, props,
                         element_name="SmallDisplacementElement3D4N")
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    z0 = np.array([nd.Z0 for nd in mp.Nodes])
    x0 = np.array([nd.X0 for nd in mp.Nodes])
    y0 = np.array([nd.Y0 for nd in mp.Nodes])
    fixed = np.zeros(box.n_dofs, dtype=bool)
    for I in range(box.n_ctrl):
        if abs(z0[I] + 1.0) < 1e-12:
            fixed[3 * I:3 * I + 3] = True
        if abs(x0[I]) < 1e-12:
            fixed[3 * I] = True
        if abs(y0[I]) < 1e-12:
            fixed[3 * I + 1] = True
    free = ~fixed
    strat = solve_linear(KM, mp)
    Kfull = system_matrix_csr(strat).tocsr()
    perm = equation_permutation(KM, mp, box)
    K = Kfull[np.ix_(np.argsort(perm), np.argsort(perm))].tocsr()
    Kf = K[np.ix_(free, free)].tocsc()
    zc = 1.0 + R_sph - delta
    top = [I for I in range(box.n_ctrl) if abs(z0[I] - 1.0) < 1e-12]

    def _zs(x, y):
        r2 = x * x + y * y
        return 1e30 if r2 >= R_sph * R_sph else zc - float(np.sqrt(R_sph * R_sph - r2))

    zsq = np.array([_zs(x0[I], y0[I]) for I in top])
    cand = [I for I, z in zip(top, zsq) if np.isfinite(z) and z < 1e29]
    zcands = np.array([_zs(x0[I], y0[I]) for I in cand])
    free_index = -np.ones(box.n_dofs, dtype=int)
    free_index[np.flatnonzero(free)] = np.arange(int(free.sum()))
    zfree = free_index[np.array([3 * I + 2 for I in cand])]
    assert np.all(zfree >= 0)
    from scipy.sparse.linalg import spsolve
    u = np.zeros(int(free.sum()))
    prev_active, rnorm = None, 1.0
    for _ in range(100):
        uf = np.zeros(box.n_dofs)
        uf[free] = u
        gap = z0[np.array(cand)] + uf[np.array([3 * I + 2 for I in cand])] - zcands
        active = gap > 0
        f_pen = np.zeros(int(free.sum()))
        np.add.at(f_pen, zfree[np.flatnonzero(active)],
                  k_pen * np.maximum(gap[active], 0.0))
        rvec = Kf @ u + f_pen
        rnorm = float(np.linalg.norm(rvec) / (1.0 + np.linalg.norm(f_pen)))
        done = prev_active is not None and np.array_equal(active, prev_active)
        prev_active = active.copy()
        if done and rnorm < 1e-9:
            break
        if not np.any(active):
            continue
        K_aug = Kf.tolil()
        for r in zfree[np.flatnonzero(active)]:
            K_aug[r, r] += k_pen
        u = u + spsolve(K_aug.tocsc(), -rvec)
    uf = np.zeros(box.n_dofs)
    uf[free] = u
    pen = np.maximum(z0[np.array(cand)] + uf[np.array([3 * I + 2 for I in cand])] - zcands, 0.0)
    fn = float(k_pen * pen.sum())
    assert int(active.sum()) > 0
    assert rnorm < 1e-6, rnorm
    assert abs(fn - F_ref) / F_ref < 0.35, (fn, F_ref)


def test_p5_friction_plateau_mu():
    """Coulomb return map: stick ramp then F_T/F_N plateaus at mu (probe parity)."""
    from wht_iga.dynamics import explicit as EX

    k_t, mu, f_n = 1e3, 0.3, 200.0
    ref = None
    ratios = []
    slipped_any = False
    for step in range(1, 41):
        g_t = np.array([0.002 * step, 0.0])
        ft, ref, slipped = EX.coulomb_friction_update(g_t, k_t, mu, f_n, ref)
        slipped_any |= bool(slipped)
        ratios.append(float(np.linalg.norm(ft)) / f_n)
    assert ratios[-1] == pytest.approx(mu, rel=1e-12)
    assert all(r <= mu * (1 + 1e-12) for r in ratios)
    assert slipped_any
    assert ratios[9] == pytest.approx(0.1, rel=1e-12)  # stick leg parity with probe
