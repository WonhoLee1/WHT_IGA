"""S15: temperature-dependent properties via PiecewiseLinearTable.

Reference: analytic table interpolation E(T) = 1 + T/100 and the stiffness
scaling K(T=100)/K(T=0) = 2 exactly (linear elasticity is homogeneous in E).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _build(E_table: bool, T: float):
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S15")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.TEMPERATURE, T)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    cl = SM.LinearElastic3DLaw()
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    if E_table:
        tab = KM.PiecewiseLinearTable()
        tab.AddRow(0.0, 1.0)
        tab.AddRow(100.0, 2.0)
        props.SetValue(KM.YOUNG_MODULUS, 1.0)
        props.SetTable(KM.YOUNG_MODULUS, KM.TEMPERATURE, tab)
    else:
        props.SetValue(KM.YOUNG_MODULUS, 1.0)
    return KM, mp, box, vol, props


def test_table_set_get_interpolates():
    KM, mp, _, _, props = _build(False, 0.0)
    tab = KM.PiecewiseLinearTable()
    tab.AddRow(0.0, 1.0)
    tab.AddRow(100.0, 2.0)
    props.SetTable(KM.YOUNG_MODULUS, KM.TEMPERATURE, tab)
    assert props.HasTable(KM.YOUNG_MODULUS, KM.TEMPERATURE)
    assert abs(props.GetTable(KM.YOUNG_MODULUS, KM.TEMPERATURE).GetValue(20.0) - 1.2) < 1e-14
    assert abs(props.GetTable(KM.YOUNG_MODULUS, KM.TEMPERATURE).GetValue(50.0) - 1.5) < 1e-14


def test_table_valued_E_does_not_reach_stiffness():
    """Documents the limitation: K(T=100, table E=1->2) == K(E=1).

    Neither LinearElastic3DLaw nor ThermalElasticIsotropic3D consumes the
    YOUNG_MODULUS table in this build; both assemble with the base value.
    """
    import KratosMultiphysics.scipy_conversion_tools  # noqa: F401

    Ks = {}
    for tag, T, Ebase in (("table100", 100.0, 1.0), ("base", 0.0, 1.0)):
        KM, mp, box, vol, props = _build(True, T)
        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, T)
        attach_cell_elements(KM, mp, vol, box, props)
        strat = KM.ResidualBasedLinearStrategy(
            mp, KM.ResidualBasedIncrementalUpdateStaticScheme(),
            KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
            False, False, False, False,
        )
        strat.SetEchoLevel(0)
        strat.Initialize()
        mp.CloneTimeStep(1.0)
        strat.Solve()
        A = system_matrix_csr(strat).toarray()
        perm = equation_permutation(KM, mp, box)
        Ks[tag] = A[np.ix_(np.argsort(perm), np.argsort(perm))]
    np.testing.assert_allclose(Ks["table100"], Ks["base"], rtol=1e-12)


def test_properties_update_failover_scales_stiffness():
    """Failover from the implementation request: per-interval Properties update."""
    import KratosMultiphysics.scipy_conversion_tools  # noqa: F401

    Ks = {}
    for E in (1.0, 2.0):
        mods = load_kratos()
        KM = mods["KratosMultiphysics"]
        SM = mods["StructuralMechanicsApplication"]
        model = KM.Model()
        mp = model.CreateModelPart("S15F")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
        vol, _ = make_nurbs_volume(KM, mp, box)
        props = mp.GetProperties()[1]
        props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
        props.SetValue(KM.POISSON_RATIO, 0.0)
        props.SetValue(KM.YOUNG_MODULUS, E)
        attach_cell_elements(KM, mp, vol, box, props)
        strat = KM.ResidualBasedLinearStrategy(
            mp, KM.ResidualBasedIncrementalUpdateStaticScheme(),
            KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
            False, False, False, False,
        )
        strat.SetEchoLevel(0)
        strat.Initialize()
        mp.CloneTimeStep(1.0)
        strat.Solve()
        A = system_matrix_csr(strat).toarray()
        perm = equation_permutation(KM, mp, box)
        Ks[E] = A[np.ix_(np.argsort(perm), np.argsort(perm))]
    np.testing.assert_allclose(Ks[2.0], 2.0 * Ks[1.0], rtol=1e-12)
