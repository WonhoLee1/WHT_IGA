"""S19: Hyperelastic / ConstitutiveModelsApplication laws on cell-wise QP geometries.

Tests:
1. Python binding availability in ``ConstitutiveModelsApplication`` (MooneyRivlin / Yeoh
   exist in C++ headers only; NeoHookean / Ogden / SVK are bound).
2. ``ConstitutiveModelsApplication.LargeStrain3DLaw`` on ``NurbsVolume`` +
   ``CreateQuadraturePointGeometries`` fails in ``Constitutive3DLaw::CalculateDomainVariables``
   because ``QuadraturePointGeometry`` does not override ``Geometry::Volume()``.
3. ``ConstitutiveLawsApplication.HyperElastic3DLaw`` (compressible Neo-Hookean)
   succeeds on both ``TotalLagrangianElement3D4N`` and ``UpdatedLagrangianElement3D4N``
   (1-cell p=1, multi-cell p=2, and compression/tension stretch sweep), matching
   the analytical uniaxial solution to < 1e-12 (target 1e-8).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize import brentq

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _load_apps():
    return load_kratos(
        (
            "IgaApplication",
            "StructuralMechanicsApplication",
            "LinearSolversApplication",
            "ConstitutiveLawsApplication",
            "ConstitutiveModelsApplication",
        )
    )


def _analytic_uniaxial_neohookean(E: float, nu: float, lam1: float) -> tuple[float, float, float]:
    """Return (lam_t, cauchy_11, pk2_11) for HyperElasticIsotropicNeoHookean3D.

    Kratos ``HyperElastic3DLaw`` uses:
        S = lambda * ln(J) * C^{-1} + mu * (I - C^{-1})
        sigma = (lambda * ln(J) * I + mu * (B - I)) / J
    where J = lam1 * lam_t^2 and lateral stress-free condition sigma_22 = 0 gives:
        lambda * ln(lam1 * lam_t^2) + mu * (lam_t^2 - 1) = 0.
    """
    mu = E / (2.0 * (1.0 + nu))
    lame_lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    lam_t = brentq(
        lambda lt: lame_lam * math.log(lam1 * lt * lt) + mu * (lt * lt - 1.0),
        0.1,
        5.0,
        xtol=1e-15,
        rtol=1e-15,
    )
    J = lam1 * lam_t * lam_t
    cauchy_11 = mu / J * (lam1 * lam1 - lam_t * lam_t)
    pk2_11 = mu * (1.0 - (lam_t / lam1) ** 2)
    return lam_t, cauchy_11, pk2_11


def _solve_uniaxial(
    KM,
    law,
    element_name: str,
    degree: int,
    n_cells: tuple[int, int, int],
    lam1: float,
    props_setter,
):
    model = KM.Model()
    mp = model.CreateModelPart("S19")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.SetBufferSize(2)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=degree)
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

    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    conv = KM.DisplacementCriteria(1e-13, 1e-13)
    conv.SetEchoLevel(0)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, conv, builder, 25, True, False, True)
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    strategy.Solve()
    return model, mp, els


def test_cm_python_bindings_inventory():
    mods = _load_apps()
    CM = mods["ConstitutiveModelsApplication"]
    assert not hasattr(CM, "MooneyRivlinModel")
    assert not hasattr(CM, "IsochoricMooneyRivlinModel")
    assert not hasattr(CM, "YeohModel")
    for name in (
        "LargeStrain3DLaw",
        "NeoHookeanModel",
        "IncompressibleNeoHookeanModel",
        "OgdenModel",
        "SaintVenantKirchhoffModel",
    ):
        assert hasattr(CM, name), name


@pytest.mark.parametrize(
    "model_name,use_nh_props",
    [
        ("NeoHookeanModel", True),
        ("IncompressibleNeoHookeanModel", True),
        ("OgdenModel", True),
        ("SaintVenantKirchhoffModel", False),
    ],
)
@pytest.mark.parametrize("element_name", ["TotalLagrangianElement3D4N", "UpdatedLagrangianElement3D4N"])
def test_cm_large_strain_fails_on_qp_geometry_volume(model_name, use_nh_props, element_name):
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    CM = mods["ConstitutiveModelsApplication"]
    law = CM.LargeStrain3DLaw(getattr(CM, model_name)())

    def _props(p):
        if use_nh_props:
            p.SetValue(KM.C10, 0.5)
            p.SetValue(KM.BULK_MODULUS, 10.0)
        else:
            p.SetValue(KM.YOUNG_MODULUS, 3.0)
            p.SetValue(KM.POISSON_RATIO, 0.3)

    with pytest.raises(RuntimeError, match="Calling base class 'Volume' method instead of derived class one"):
        _solve_uniaxial(KM, law, element_name, degree=1, n_cells=(1, 1, 1), lam1=1.2, props_setter=_props)


# TotalLagrangian was left out because of intermittent DETJ0=-0.5/NaN stresses. Root cause
# (REVIEW_P0_P1 §1 #8): the NurbsVolume was freed when ``_solve_uniaxial`` returned, and TL's
# reference Jacobian re-evaluates the parent geometry. ``make_nurbs_volume`` now registers the
# volume on the model part; repeated runs guard the regression.
@pytest.mark.parametrize("repeat", range(4))
@pytest.mark.parametrize("element_name", ["UpdatedLagrangianElement3D4N", "TotalLagrangianElement3D4N"])
@pytest.mark.parametrize("degree,n_cells", [(1, (1, 1, 1)), (2, (2, 2, 2))])
def test_cl_hyperelastic_uniaxial_matches_analytic(element_name, degree, n_cells, repeat):
    import gc

    gc.collect()
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    E, nu, lam1 = 3.0, 0.3, 1.2

    _model, mp, els = _solve_uniaxial(
        KM,
        CL.HyperElastic3DLaw(),
        element_name,
        degree=degree,
        n_cells=n_cells,
        lam1=lam1,
        props_setter=lambda p: (p.SetValue(KM.YOUNG_MODULUS, E), p.SetValue(KM.POISSON_RATIO, nu)),
    )

    lam_t_ref, sig11_ref, pk2_11_ref = _analytic_uniaxial_neohookean(E, nu, lam1)
    sig_all = np.vstack(
        [np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)) for el in els]
    )
    pk2_all = np.vstack(
        [np.array(el.CalculateOnIntegrationPoints(KM.PK2_STRESS_VECTOR, mp.ProcessInfo)) for el in els]
    )
    uy_top = np.array([n.GetSolutionStepValue(KM.DISPLACEMENT_Y) for n in mp.Nodes if abs(n.Y0 - 1.0) < 1e-12])

    sig_err = float(np.max(np.abs(sig_all[:, 0] - sig11_ref)))
    lat_err = float(np.max(np.abs(sig_all[:, 1:])))
    pk2_err = float(np.max(np.abs(pk2_all[:, 0] - pk2_11_ref)))
    lam_t_err = float(np.max(np.abs((1.0 + uy_top) - lam_t_ref)))

    print(
        f"S19 {element_name} p={degree} cells={n_cells}: "
        f"sig_err={sig_err:.3e}, pk2_err={pk2_err:.3e}, lat_err={lat_err:.3e}, lam_t_err={lam_t_err:.3e}"
    )
    assert sig_err < 1e-11
    assert pk2_err < 1e-11
    assert lat_err < 1e-11
    assert lam_t_err < 1e-11


@pytest.mark.parametrize("lam1", [0.85, 1.05, 1.2, 1.35, 1.5])
def test_cl_hyperelastic_stretch_curve(lam1):
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    E, nu = 3.0, 0.3

    _model, mp, els = _solve_uniaxial(
        KM,
        CL.HyperElastic3DLaw(),
        "UpdatedLagrangianElement3D4N",
        degree=1,
        n_cells=(1, 1, 1),
        lam1=lam1,
        props_setter=lambda p: (p.SetValue(KM.YOUNG_MODULUS, E), p.SetValue(KM.POISSON_RATIO, nu)),
    )
    _, sig11_ref, _ = _analytic_uniaxial_neohookean(E, nu, lam1)
    sig_all = np.array(els[0].CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo))
    err = float(np.max(np.abs(sig_all[:, 0] - sig11_ref)))
    print(f"S19 stretch lam1={lam1}: sig11={sig_all[0, 0]:.12f}, ref={sig11_ref:.12f}, err={err:.3e}")
    assert err < 1e-11
