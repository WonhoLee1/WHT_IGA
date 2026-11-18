# -*- coding: utf-8 -*-
"""S28: thermal strain as Kratos initial strain on J2 plasticity (IGA QP elements).

Kratos 10.4 has no thermo-plastic constitutive law (thermal laws are
``ThermalElasticIsotropic3D`` and ``SmallStrainThermalIsotropicDamage*``
only). The generic plasticity laws subtract the ``InitialState`` strain in
both ``CalculateMaterialResponse`` and ``FinalizeMaterialResponse``, so the
thermal strain alpha*dT can be imposed as initial strain via
``SetInitialStateProcess3D`` (whole model part or one sub model part per cell
group: ``Geometry.SetValue`` is not bound in Python).

All references are analytic (E = 100, nu = 0, sigma_y = 1, perfect plasticity).
"""

from __future__ import annotations

import pytest

from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_E, _SY, _ALPHA = 100.0, 1.0, 1e-3


def _mods():
    m = load_kratos()
    return m["KratosMultiphysics"], m["ConstitutiveLawsApplication"]


def _bar(KM, CL, law, element_name, nx=1, lateral_fixed=False):
    """Clamped bar [0, nx] x [0,1]^2 (x fixed at both ends), symmetry planes y=0, z=0."""
    model = KM.Model()
    mp = model.CreateModelPart("S28")
    for v in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(v)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(float(nx), 1.0, 1.0), n_cells=(nx, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    pr = mp.GetProperties()[1]
    for var, val in ((KM.YOUNG_MODULUS, _E), (KM.POISSON_RATIO, 0.0),
                     (KM.YIELD_STRESS, _SY), (KM.FRACTURE_ENERGY, 1e10)):
        pr.SetValue(var, val)
    pr.SetValue(CL.HARDENING_CURVE, 3)  # perfect plasticity
    pr.SetValue(KM.CONSTITUTIVE_LAW, law)
    els = attach_cell_elements(KM, mp, vol, box, pr, element_name=element_name)
    mp.CloneTimeStep(0.0)
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    for n in mp.Nodes:
        if abs(n.X0) < 1e-12 or abs(n.X0 - nx) < 1e-12:
            n.Fix(KM.DISPLACEMENT_X)
        if lateral_fixed or abs(n.Y0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Y)
        if lateral_fixed or abs(n.Z0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Z)
    return mp, els


def _set_thermal_strain(KM, model_part, eps_th):
    strain, stress, F = KM.Vector(6), KM.Vector(6), KM.Matrix(3, 3)
    for i in range(6):
        strain[i] = eps_th if i < 3 else 0.0
        stress[i] = 0.0
    for i in range(3):
        for j in range(3):
            F[i, j] = 1.0 if i == j else 0.0
    KM.SetInitialStateProcess3D(model_part, strain, stress, F).ExecuteInitializeSolutionStep()


def _stress(KM, mp, el, var=None):
    var = KM.CAUCHY_STRESS_VECTOR if var is None else var
    return list(el.CalculateOnIntegrationPoints(var, mp.ProcessInfo)[0])


def _run_path(KM, mp, els, levels, apply_temperature, var=None):
    """run_load_stepping over temperature levels; stress after each converged level."""
    seen, applied = [], []

    def load(dT):
        if applied:
            seen.append(_stress(KM, mp, els[0], var)[0])
        applied.append(dT)
        apply_temperature(dT)

    out = NL.run_load_stepping(mp, KM, lambda: NL.make_newton_strategy(KM, mp), levels, load)
    seen.append(_stress(KM, mp, els[0], var)[0])
    assert out["steps"] == len(levels) and out["n_solves"] == len(levels)
    return seen


def test_heat_yield_cool_residual_stress():
    """Clamped bar, uniform dT path 5 -> 15 -> 0 (E*alpha*dT = 0.5, 1.5, 0):
    sigma_xx = -0.5 (elastic), -1.0 (yield), then +0.5 residual after cooling
    (plastic strain -(alpha*15 - sigma_y/E) = -0.005 left behind)."""
    KM, CL = _mods()
    mp, els = _bar(KM, CL, CL.SmallStrainIsotropicPlasticity3DVonMisesVonMises(),
                   "SmallDisplacementElement3D4N")
    got = _run_path(KM, mp, els, [5.0, 15.0, 0.0],
                    lambda dT: _set_thermal_strain(KM, mp, _ALPHA * dT))
    for g, ref in zip(got, (-0.5, -1.0, 0.5)):
        assert g == pytest.approx(ref, abs=1e-9), got


@pytest.mark.parametrize("dT, ref_a, ref_b, ux_ref", [
    (10.0, (-0.5, -1.0, -1.0), (-0.5, 0.0, 0.0), 0.005),   # elastic
    (40.0, (-2.0, -3.0, -3.0), (-2.0, -1.0, -1.0), 0.04),  # both cells at yield
])
def test_per_cell_thermal_strain_two_materials(dT, ref_a, ref_b, ux_ref):
    """Two cells, alpha_A = 1e-3, alpha_B = 0, lateral faces fixed, clamped ends.

    Elastic: eps_A = alpha*dT/2, sigma_xx = -E*alpha*dT/2, sigma_yy,A = -E*alpha*dT.
    dT = 40: |sigma_xx - sigma_yy| = sigma_y in both cells with isochoric plastic
    strain (2a, -a, -a); solving equilibrium + compatibility gives eps_A = 0.04,
    sigma_A = (-2, -3, -3), sigma_B = (-2, -1, -1)."""
    KM, CL = _mods()
    mp, els = _bar(KM, CL, CL.SmallStrainIsotropicPlasticity3DVonMisesVonMises(),
                   "SmallDisplacementElement3D4N", nx=2, lateral_fixed=True)
    cell_a = mp.CreateSubModelPart("CellA")
    cell_a.AddElements([els[0].Id])
    cell_b = mp.CreateSubModelPart("CellB")
    cell_b.AddElements([els[1].Id])

    def apply(t):
        _set_thermal_strain(KM, cell_a, 1e-3 * t)
        _set_thermal_strain(KM, cell_b, 0.0)

    out = NL.run_load_stepping(mp, KM, lambda: NL.make_newton_strategy(KM, mp), [dT], apply)
    assert out["steps"] == 1
    for el, ref in ((els[0], ref_a), (els[1], ref_b)):
        assert _stress(KM, mp, el)[:3] == pytest.approx(list(ref), abs=1e-9)
    for n in mp.Nodes:
        if abs(n.X0 - 1.0) < 1e-12:
            assert n.GetSolutionStepValue(KM.DISPLACEMENT_X) == pytest.approx(ux_ref, abs=1e-12)


def test_yield_stress_change_between_steps_is_honoured():
    """T-dependent yield proxy: lowering YIELD_STRESS 1.0 -> 0.6 between two
    steps at the same dT = 15 gives sigma_xx = -0.6 (the law re-reads it)."""
    KM, CL = _mods()
    mp, els = _bar(KM, CL, CL.SmallStrainIsotropicPlasticity3DVonMisesVonMises(),
                   "SmallDisplacementElement3D4N")
    calls = []

    def apply(dT):
        calls.append(dT)
        if len(calls) == 2:
            mp.GetProperties()[1].SetValue(KM.YIELD_STRESS, 0.6)
        _set_thermal_strain(KM, mp, _ALPHA * dT)

    got = _run_path(KM, mp, els, [15.0, 15.0], apply)
    assert got == pytest.approx([-1.0, -0.6], abs=1e-9)


def test_finite_strain_total_lagrangian_path():
    """NLGEOM route for E03: TL element + FiniteStrainIsotropicPlasticity3DVonMisesVonMises.
    With the thermal strain added to the Green-Lagrange strain, the PK2 path
    5 -> 15 -> 0 equals the small-strain answers (-0.5, -1.0, +0.5)."""
    KM, CL = _mods()
    mp, els = _bar(KM, CL, CL.FiniteStrainIsotropicPlasticity3DVonMisesVonMises(),
                   "TotalLagrangianElement3D4N")
    got = _run_path(KM, mp, els, [5.0, 15.0, 0.0],
                    lambda dT: _set_thermal_strain(KM, mp, _ALPHA * dT), KM.PK2_STRESS_VECTOR)
    assert got == pytest.approx([-0.5, -1.0, 0.5], abs=1e-9)
