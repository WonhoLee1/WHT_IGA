"""S14: contact on independent surface layers (ContactStructuralMechanicsApplication).

Question (v2 §5.8): does ContactSM penalty/ALM mortar work on independent
free-DOF surface layers (spring-coupled)?

Scope (time-boxed spike): (1) ALM/mesh-tying mortar conditions instantiate on
plain Triangle3D3 layers; (2) spring-tie force transmission between two stacked
immersed blocks is exact (equilibrium + interface continuity); full frictional
contact solve is deferred (partial).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
    "ContactStructuralMechanicsApplication",
)


def _load():  # type: ignore[no-untyped-def]
    from wht_iga.kratos_backend import load_kratos

    return load_kratos(_APPS)


def test_mortar_conditions_instantiate() -> None:
    """Registration evidence: both mortar families exist and construct on a
    Triangle3D3 pair geometry (full pairing needs Phase-5 setup)."""
    mods = _load()
    KM = mods["KratosMultiphysics"]
    model = KM.Model()
    mp = model.CreateModelPart("S14reg")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    a = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
    b = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
    c = mp.CreateNewNode(3, 0.0, 1.0, 0.0)
    tri = KM.Triangle3D3(a, b, c)
    props = mp.CreateNewProperties(1)
    out = []
    for name in ("ALMNVFrictionlessMortarContactCondition3D3N", "MeshTyingMortarCondition3D3N"):
        try:
            cond = mp.CreateNewCondition(name, len(out) + 1, tri, props)
            out.append(f"{name}: created ({type(cond).__name__})")
        except Exception as exc:  # noqa: BLE001 - the message IS the evidence
            out.append(f"{name}: {type(exc).__name__}: {exc}")
    print("\nS14 registration:\n  " + "\n  ".join(out))
    # Both families must construct (S14.md: registration evidence for both).
    assert len(out) == 2
    assert all(": created" in line for line in out), out


def _block(KM, SM, mp, box, z0: float, first_id: int, E: float, prop_id: int):  # type: ignore[no-untyped-def]
    volume, _ = make_nurbs_volume(KM, mp, box, first_node_id=first_id)
    props = mp.CreateNewProperties(prop_id)
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, 0.3)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    nx, ny, nz = box.n_cells
    eid0 = mp.NumberOfElements()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                qv = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qv, 2, cell_gauss_points(box, i, j, k))
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid0 + 1, qv[0], props)
                eid0 += 1
    return volume


def test_spring_tie_transmits_force() -> None:
    """Two stacked p=1 cubes tied node-to-node by stiff springs (BC-C surrogate):
    pressure 2.0 on the top face must balance bottom reactions exactly and the
    interface jump must be O(F/k)."""
    mods = _load()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S14tie")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.ROTATION)
    mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    lo_box = BoxSpline(lo=(0, 0, 0), hi=(1, 1, 1), n_cells=(1, 1, 1), degree=1)
    hi_box = BoxSpline(lo=(0, 0, 1), hi=(1, 1, 2), n_cells=(1, 1, 1), degree=1)
    vol_lo = _block(KM, SM, mp, lo_box, 0.0, 1, 1000.0, 1)
    vol_hi = _block(KM, SM, mp, hi_box, 1.0, 100, 1000.0, 2)
    KM.VariableUtils().AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X, mp)
    KM.VariableUtils().AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y, mp)
    KM.VariableUtils().AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z, mp)
    for node in mp.Nodes:
        if abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)
    for node in mp.Nodes:
        if node.Id < 50_000:
            node.Fix(KM.ROTATION_X)
            node.Fix(KM.ROTATION_Y)
            node.Fix(KM.ROTATION_Z)

    # Coincident interface pairs -> stiff springs (normal + tangential).
    stiff = mp.CreateNewProperties(3)
    pairs = 0
    lo_face = [n for n in mp.Nodes if abs(n.Z0 - 1.0) < 1e-12 and n.Id < 100]
    hi_face = [n for n in mp.Nodes if abs(n.Z0 - 1.0) < 1e-12 and n.Id >= 100]
    for nl in lo_face:
        partner = next(n for n in hi_face if abs(n.X0 - nl.X0) < 1e-12 and abs(n.Y0 - nl.Y0) < 1e-12)
        el = mp.CreateNewElement("SpringDamperElement3D", mp.NumberOfElements() + 1, [nl.Id, partner.Id], stiff)
        el.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, [1.0e10, 1.0e10, 1.0e10])
        el.SetValue(SM.NODAL_ROTATIONAL_STIFFNESS, [1.0e10, 1.0e10, 1.0e10])
        pairs += 1
    assert pairs == 4

    # Pressure on the very top face.
    tris = [
        np.array([[0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [1.0, 1.0, 2.0]]),
        np.array([[0.0, 0.0, 2.0], [1.0, 1.0, 2.0], [0.0, 1.0, 2.0]]),
    ]
    qp_geoms, wrappers, _, _ = make_surface_quadrature(KM, mp, vol_hi, tris, rule_order=2, first_node_id=50_000)
    cp = mp.CreateNewProperties(4)
    for i, qp in enumerate(qp_geoms):
        cond = mp.CreateNewCondition("LoadCondition", 500 + i, qp, cp)
        cond.SetValue(KM.PRESSURE, 2.0)

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()

    fz = sum(node.GetSolutionStepValue(KM.REACTION_Z) for node in mp.Nodes if abs(node.Z0) < 1e-12)
    jumps = []
    for nl in lo_face:
        partner = next(n for n in hi_face if abs(n.X0 - nl.X0) < 1e-12 and abs(n.Y0 - nl.Y0) < 1e-12)
        jumps.append(
            max(
                abs(nl.GetSolutionStepValue(v) - partner.GetSolutionStepValue(v))
                for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
            )
        )
    print(f"\nS14 spring-tie: bottom reaction z = {fz:.10f} (expect -2.0), max interface jump = {max(jumps):.3e}")
    np.testing.assert_allclose(fz, -2.0, rtol=0.0, atol=1e-6)
    assert max(jumps) < 1e-6
    _ = (wrappers, lo_box, hi_box)
