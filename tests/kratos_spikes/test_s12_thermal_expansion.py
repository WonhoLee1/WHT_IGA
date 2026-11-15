"""S12: thermal expansion (ConstitutiveLaws ThermalElasticIsotropic3D).

Question (v2 §5.8): does ``ThermalElasticIsotropic3D`` compute thermal strain
from nodal TEMPERATURE on immersed volume QP geometry?

Reference: unit cube, uniform ΔT=100 (T=120, T_ref=20), α=1e-5, E=1000, nu=0.3.
Free expansion with compatible minimal constraints: corner (1,1,1) must displace
α·ΔT·(1,1,1) = 1e-3·(1,1,1) and all fixed-DOF reactions must vanish (stress-free).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.thermal import ensure_thermal_variables, make_thermal_elastic_properties
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def test_free_thermal_expansion() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, CLA = (
        mods["KratosMultiphysics"],
        mods["StructuralMechanicsApplication"],
        mods["ConstitutiveLawsApplication"],
    )
    model = KM.Model()
    mp = model.CreateModelPart("S12")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    ensure_thermal_variables(mp, KM)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    alpha, T_ref, T_hot = 1.0e-5, 20.0, 120.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    volume, _ = make_nurbs_volume(KM, mp, box)
    props = make_thermal_elastic_properties(mp, KM, SM, CLA, E=1000.0, nu=0.3, alpha=alpha, T_ref=T_ref)
    qv = KM.GeometriesVector()
    volume.CreateQuadraturePointGeometries(qv, 2, cell_gauss_points(box, 0, 0, 0))
    mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qv[0], props)

    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, T_hot)

    def _node_at(x: float, y: float, z: float):  # type: ignore[no-untyped-def]
        return next(n for n in mp.Nodes if abs(n.X0 - x) < 1e-12 and abs(n.Y0 - y) < 1e-12 and abs(n.Z0 - z) < 1e-12)

    # Compatible minimal constraints (all consistent with uniform expansion).
    fixed: list[tuple] = []
    for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
        _node_at(0.0, 0.0, 0.0).Fix(var)
        fixed.append((_node_at(0.0, 0.0, 0.0), var))
    for var in (KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
        _node_at(1.0, 0.0, 0.0).Fix(var)
        fixed.append((_node_at(1.0, 0.0, 0.0), var))
    _node_at(0.0, 1.0, 0.0).Fix(KM.DISPLACEMENT_Z)
    fixed.append((_node_at(0.0, 1.0, 0.0), KM.DISPLACEMENT_Z))

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()

    corner = _node_at(1.0, 1.0, 1.0)
    u = np.array([corner.GetSolutionStepValue(v) for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)])
    expected = alpha * (T_hot - T_ref) * np.ones(3)
    print(f"\nS12 corner u = {u} (expect {expected})")
    np.testing.assert_allclose(u, expected, rtol=1e-9, atol=0.0)
    rmap = {KM.DISPLACEMENT_X: KM.REACTION_X, KM.DISPLACEMENT_Y: KM.REACTION_Y, KM.DISPLACEMENT_Z: KM.REACTION_Z}
    reactions = [node.GetSolutionStepValue(rmap[var]) for node, var in fixed]
    print(f"S12 max |reaction| = {max(abs(r) for r in reactions):.3e} (expect 0; E*a*dT = 1.0)")
    assert max(abs(r) for r in reactions) < 1e-6
    _ = (volume, SM)
