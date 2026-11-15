"""S10: nodal/spring elements for the BC-B / BC-C virtual-node path.

Question (v2 §5.8): do SM ``NodalConcentratedElement3D1N`` and
``SpringDamperElement3D`` work on slave (virtual) nodes?

Reference: analytic spring extensions (u = F/k). Both solves are linear
static (NewtonRaphson, SkylineLU) with tol 1e-9 on |u - F/k|.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.kratos


def _solve(KM, mp) -> None:  # type: ignore[no-untyped-def]
    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()


def test_nodal_concentrated_stiffness() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S10a")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    node = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    node.Fix(KM.DISPLACEMENT_Y)
    node.Fix(KM.DISPLACEMENT_Z)

    props = mp.CreateNewProperties(1)
    elem = mp.CreateNewElement("NodalConcentratedElement3D1N", 1, [1], props)
    elem.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, [2.5, 0.0, 0.0])
    cond = mp.CreateNewCondition("PointLoadCondition3D1N", 1, [1], props)
    node.SetSolutionStepValue(SM.POINT_LOAD_X, 5.0)
    _solve(KM, mp)
    ux = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
    print(f"\nS10a NodalConcentrated: u_x = {ux:.12f} (expect 5.0/2.5 = 2.0)")
    assert abs(ux - 2.0) < 1e-9


def test_spring_damper_between_nodes() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S10b")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
    mp.AddNodalSolutionStepVariable(KM.ROTATION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    node_a = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
    node_b = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, mp)
    KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, mp)
    KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, mp)
    for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z,
                KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z):
        node_a.Fix(var)

    props = mp.CreateNewProperties(1)
    spring = mp.CreateNewElement("SpringDamperElement3D", 1, [1, 2], props)
    spring.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, [3.0, 1.0e6, 1.0e6])
    spring.SetValue(SM.NODAL_ROTATIONAL_STIFFNESS, [1.0e6, 1.0e6, 1.0e6])
    mp.CreateNewCondition("PointLoadCondition3D1N", 1, [2], props)
    node_b.SetSolutionStepValue(SM.POINT_LOAD_X, 4.0)
    _solve(KM, mp)
    ux = node_b.GetSolutionStepValue(KM.DISPLACEMENT_X)
    uy = node_b.GetSolutionStepValue(KM.DISPLACEMENT_Y)
    print(f"\nS10b SpringDamper: u_x = {ux:.12f} (expect 4.0/3.0 = 1.333...), u_y = {uy:.3e}")
    assert abs(ux - 4.0 / 3.0) < 1e-9
    assert abs(uy) < 1e-9
