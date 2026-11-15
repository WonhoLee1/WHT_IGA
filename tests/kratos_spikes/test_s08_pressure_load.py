"""S08: pressure / traction / follower loads on boundary QP geometry.

Question (v2 §5.8): does IgaApp ``LoadCondition`` with PRESSURE /
SURFACE_LOAD / PRESSURE_FOLLOWER_LOAD on surface QP geometry integrate the
correct total force, and what is the sign convention?

Reference: analytic totals on the unit-cube top face (area 1.0):
|F| = |p| * A. Bottom-face reaction sum must balance it (tol 1e-9).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _top_face_triangles() -> list[np.ndarray]:
    return [
        np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]),
        np.array([[0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]]),
    ]


def _build(KM, SM, name: str, E: float, nu: float):  # type: ignore[no-untyped-def]
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    volume, _ = make_nurbs_volume(KM, mp, box)
    vol_props = mp.CreateNewProperties(1)
    vol_props.SetValue(KM.YOUNG_MODULUS, E)
    vol_props.SetValue(KM.POISSON_RATIO, nu)
    vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    qp_vol = KM.GeometriesVector()
    volume.CreateQuadraturePointGeometries(qp_vol, 2, cell_gauss_points(box, 0, 0, 0))
    mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qp_vol[0], vol_props)
    for node in mp.Nodes:
        if abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)
    # NOTE: volume + surface wrappers must stay alive for the whole analysis
    # (lifetime rule, see wht_iga/kratos_backend/surface.py). Returned so the
    # caller retains them.
    qp_geoms, wrappers, _, _ = make_surface_quadrature(KM, mp, volume, _top_face_triangles(), rule_order=2)
    return mp, box, qp_geoms, volume, wrappers


def _solve(KM, mp, move_mesh: bool = False) -> None:  # type: ignore[no-untyped-def]
    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(
        mp, scheme, crit, builder, 15, True, False, move_mesh
    )
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()


def _bottom_reaction_z(KM, mp) -> float:  # type: ignore[no-untyped-def]
    return sum(
        node.GetSolutionStepValue(KM.REACTION_Z) for node in mp.Nodes if abs(node.Z0) < 1e-12
    )


def test_pressure_total_and_sign() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]
    mp, _box, qp_geoms, volume_keepalive, wrappers_keepalive = _build(KM, SM, "S08p", 1000.0, 0.3)
    cond_props = mp.CreateNewProperties(2)
    for i, qp in enumerate(qp_geoms):
        cond = mp.CreateNewCondition("LoadCondition", 100 + i, qp, cond_props)
        cond.SetValue(KM.PRESSURE, 2.0)
    _solve(KM, mp)
    fz = _bottom_reaction_z(KM, mp)
    print(f"\nS08 PRESSURE=2.0 bottom reaction sum z = {fz:.12f} (expect -2.0: +PRESSURE pushes along outward normal)")
    np.testing.assert_allclose(fz, -2.0, rtol=0.0, atol=1e-9)


def test_surface_load_total() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]
    mp, _box, qp_geoms, volume_keepalive, wrappers_keepalive = _build(KM, SM, "S08t", 1000.0, 0.3)
    cond_props = mp.CreateNewProperties(2)
    traction = KM.Vector(3)
    traction[0], traction[1], traction[2] = 0.0, 0.0, -3.0
    for i, qp in enumerate(qp_geoms):
        cond = mp.CreateNewCondition("LoadCondition", 100 + i, qp, cond_props)
        cond.SetValue(SM.SURFACE_LOAD, traction)
    _solve(KM, mp)
    fz = _bottom_reaction_z(KM, mp)
    print(f"\nS08 SURFACE_LOAD=(0,0,-3) bottom reaction sum z = {fz:.12f} (expect +3.0)")
    np.testing.assert_allclose(fz, 3.0, rtol=0.0, atol=1e-9)


def test_pressure_follower_matches_dead_load() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]
    mp, _box, qp_geoms, volume_keepalive, wrappers_keepalive = _build(KM, SM, "S08f", 1.0e6, 0.3)
    cond_props = mp.CreateNewProperties(2)
    for i, qp in enumerate(qp_geoms):
        cond = mp.CreateNewCondition("LoadCondition", 100 + i, qp, cond_props)
        cond.SetValue(IGA.PRESSURE_FOLLOWER_LOAD, 2.0)
    _solve(KM, mp, move_mesh=True)
    fz = _bottom_reaction_z(KM, mp)
    print(f"\nS08 FOLLOWER=2.0 bottom reaction sum z = {fz:.12f} (expect -2.0, same sign as dead load)")
    assert fz < 0.0
    assert abs(fz + 2.0) / 2.0 < 1e-5
    np.testing.assert_allclose(fz, -2.0, rtol=0.0, atol=2.0e-5)
