"""S22: ContactStructuralMechanicsApplication processes on spline-coupled surface triangles.

Question: Do mesh_tying_process (mortar) and penalty_contact_process work on
surface-layer Triangle3D3 nodes coupled to splines via SpringDamperElement3D
or LinearMasterSlaveConstraint?

Findings tested directly:
1. Both mesh_tying_process and penalty_contact_process fail in ExecuteInitialize()
   when the main ModelPart contains IGA QuadraturePointGeometry elements, because
   InterfacePreprocessCondition::CheckOnTheFace calls Geometry::GenerateFaces(),
   which is not implemented on QuadraturePointGeometry.
2. Running them on an element-free surface ModelPart also fails in
   InterfacePreprocessCondition::CheckAndCreateProperties ("It was not possible
   to add a property") because no parent solid element exists.
3. As a working surface-layer baseline (surface Triangle3D3 nodes tied to spline
   control points via LinearMasterSlaveConstraint + interface SpringDamperElement3D),
   a stacked two-block tie matches a monolithic 2-cell block under uniaxial
   compression (displacement and reaction force verification).
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
    "ContactStructuralMechanicsApplication",
)


def _setup_two_blocks_with_surface_layers():
    """Create two stacked p=1 blocks (z in [0,1] and [1,2]) with surface Triangle3D3 layers at z=1."""
    mods = load_kratos(_APPS)
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    model = KM.Model()
    mp = model.CreateModelPart("Structure")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.ROTATION)
    mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
    mp.AddNodalSolutionStepVariable(KM.NORMAL)
    mp.AddNodalSolutionStepVariable(KM.NODAL_H)
    mp.AddNodalSolutionStepVariable(KM.NODAL_AREA)
    mp.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    lo_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    hi_box = BoxSpline(lo=(0.0, 0.0, 1.0), hi=(1.0, 1.0, 2.0), n_cells=(1, 1, 1), degree=1)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1000.0)
    props.SetValue(KM.POISSON_RATIO, 0.3)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    vol_lo, _ = make_nurbs_volume(KM, mp, lo_box, first_node_id=1)
    attach_cell_elements(KM, mp, vol_lo, lo_box, props, start_id=1)
    vol_hi, _ = make_nurbs_volume(KM, mp, hi_box, first_node_id=100)
    attach_cell_elements(KM, mp, vol_hi, hi_box, props, start_id=100)

    smp_lo = mp.CreateSubModelPart("InterfaceLo")
    smp_hi = mp.CreateSubModelPart("InterfaceHi")
    pts = [(0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    s_lo = [smp_lo.CreateNewNode(1000 + i, x, y, z) for i, (x, y, z) in enumerate(pts)]
    s_hi = [smp_hi.CreateNewNode(2000 + i, x, y, z) for i, (x, y, z) in enumerate(pts)]

    cprop = mp.GetProperties()[2]
    smp_lo.CreateNewCondition("SurfaceCondition3D3N", 1000, [s_lo[0].Id, s_lo[1].Id, s_lo[2].Id], cprop)
    smp_lo.CreateNewCondition("SurfaceCondition3D3N", 1001, [s_lo[0].Id, s_lo[2].Id, s_lo[3].Id], cprop)
    smp_hi.CreateNewCondition("SurfaceCondition3D3N", 2000, [s_hi[0].Id, s_hi[2].Id, s_hi[1].Id], cprop)
    smp_hi.CreateNewCondition("SurfaceCondition3D3N", 2001, [s_hi[0].Id, s_hi[3].Id, s_hi[2].Id], cprop)

    for v, r in zip(
        (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z, KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z),
        (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z, KM.REACTION_MOMENT_X, KM.REACTION_MOMENT_Y, KM.REACTION_MOMENT_Z),
    ):
        KM.VariableUtils().AddDof(v, r, mp)

    return mods, model, mp, s_lo, s_hi, props, vol_lo, vol_hi


def test_s22_mesh_tying_process_fails_on_qp_geometry() -> None:
    """mesh_tying_process fails on IGA QuadraturePointGeometry elements (GenerateFaces unimplemented)."""
    mods, model, mp, *_ = _setup_two_blocks_with_surface_layers()
    import KratosMultiphysics.ContactStructuralMechanicsApplication.mesh_tying_process as mtp

    KM = mods["KratosMultiphysics"]
    mp.SetBufferSize(2)
    cfg = {
        "Parameters": {
            "model_part_name": "Structure",
            "mesh_tying_model_part": {"0": ["InterfaceLo", "InterfaceHi"]},
            "assume_master_slave": {"0": ["InterfaceLo"]},
            "mesh_tying_property_ids": {"0": 3},
            "consider_static_condensation": True,
        }
    }
    proc = mtp.Factory(KM.Parameters(json.dumps(cfg)), model)
    with pytest.raises(RuntimeError, match="GenerateFaces") as exc_info:
        proc.ExecuteInitialize()
    print("S22 mesh_tying_process expected failure:", str(exc_info.value).splitlines()[0])


def test_s22_penalty_contact_process_fails_on_qp_geometry() -> None:
    """penalty_contact_process fails on IGA QuadraturePointGeometry elements (GenerateFaces unimplemented)."""
    mods, model, mp, *_ = _setup_two_blocks_with_surface_layers()
    import KratosMultiphysics.ContactStructuralMechanicsApplication.penalty_contact_process as pcp

    KM = mods["KratosMultiphysics"]
    mp.SetBufferSize(2)
    cfg = {
        "Parameters": {
            "model_part_name": "Structure",
            "contact_model_part": {"0": ["InterfaceLo", "InterfaceHi"]},
            "assume_master_slave": {"0": ["InterfaceLo"]},
            "contact_property_ids": {"0": 3},
        }
    }
    proc = pcp.Factory(KM.Parameters(json.dumps(cfg)), model)
    with pytest.raises(RuntimeError, match="GenerateFaces") as exc_info:
        proc.ExecuteInitialize()
    print("S22 penalty_contact_process expected failure:", str(exc_info.value).splitlines()[0])


def test_s22_mesh_tying_process_fails_without_parent_elements() -> None:
    """Running mesh_tying_process on an element-free surface ModelPart fails in CheckAndCreateProperties."""
    mods = load_kratos(_APPS)
    KM = mods["KratosMultiphysics"]
    import KratosMultiphysics.ContactStructuralMechanicsApplication.mesh_tying_process as mtp

    model = KM.Model()
    mp = model.CreateModelPart("ContactSurf")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.NORMAL)
    mp.AddNodalSolutionStepVariable(KM.NODAL_H)
    mp.AddNodalSolutionStepVariable(KM.NODAL_AREA)
    mp.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    smp_lo = mp.CreateSubModelPart("InterfaceLo")
    smp_hi = mp.CreateSubModelPart("InterfaceHi")
    pts = [(0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    n1 = [smp_lo.CreateNewNode(1000 + i, x, y, z) for i, (x, y, z) in enumerate(pts)]
    n2 = [smp_hi.CreateNewNode(2000 + i, x, y, z) for i, (x, y, z) in enumerate(pts)]
    cprop = mp.GetProperties()[2]
    smp_lo.CreateNewCondition("SurfaceCondition3D3N", 1000, [n1[0].Id, n1[1].Id, n1[2].Id], cprop)
    smp_hi.CreateNewCondition("SurfaceCondition3D3N", 2000, [n2[0].Id, n2[2].Id, n2[1].Id], cprop)
    mp.SetBufferSize(2)

    cfg = {
        "Parameters": {
            "model_part_name": "ContactSurf",
            "mesh_tying_model_part": {"0": ["InterfaceLo", "InterfaceHi"]},
            "assume_master_slave": {"0": ["InterfaceLo"]},
            "consider_static_condensation": True,
            "scale_factor_parameters": {"manual_scale_factor": True, "scale_factor": 1.0e3},
        }
    }
    proc = mtp.Factory(KM.Parameters(json.dumps(cfg)), model)
    with pytest.raises(RuntimeError, match="not possible to add a property") as exc_info:
        proc.ExecuteInitialize()
    print("S22 element-free ContactSurf failure:", str(exc_info.value).splitlines()[0])


def test_s22_surface_layer_spring_tie_vs_monolithic_block() -> None:
    """Fallback verification: surface Triangle3D3 nodes tied to spline via LinearMasterSlaveConstraint
    and coupled across interface with stiff SpringDamperElement3D vs monolithic 2-cell block."""
    mods, model, mp, s_lo, s_hi, props, vol_lo, vol_hi = _setup_two_blocks_with_surface_layers()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    lo_face = [n for n in mp.Nodes if n.Id < 100 and abs(n.Z0 - 1.0) < 1e-12]
    hi_face = [n for n in mp.Nodes if 100 <= n.Id < 1000 and abs(n.Z0 - 1.0) < 1e-12]
    cid = 1
    for sn in s_lo:
        mn = next(n for n in lo_face if abs(n.X0 - sn.X0) < 1e-12 and abs(n.Y0 - sn.Y0) < 1e-12)
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            add_linear_tie(KM, mp, cid, [sn.GetDof(var)], [mn.GetDof(var)], np.array([[1.0]]))
            cid += 1
    for sn in s_hi:
        mn = next(n for n in hi_face if abs(n.X0 - sn.X0) < 1e-12 and abs(n.Y0 - sn.Y0) < 1e-12)
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            add_linear_tie(KM, mp, cid, [sn.GetDof(var)], [mn.GetDof(var)], np.array([[1.0]]))
            cid += 1

    stiff = mp.CreateNewProperties(3)
    for nl, nh in zip(s_lo, s_hi):
        el = mp.CreateNewElement("SpringDamperElement3D", mp.NumberOfElements() + 1, [nl.Id, nh.Id], stiff)
        el.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, [1.0e10, 1.0e10, 1.0e10])
        el.SetValue(SM.NODAL_ROTATIONAL_STIFFNESS, [1.0e10, 1.0e10, 1.0e10])

    for node in mp.Nodes:
        if node.Id < 100 and abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)
        elif 100 <= node.Id < 1000 and abs(node.Z0 - 2.0) < 1e-12:
            node.SetSolutionStepValue(SM.POINT_LOAD_Z, -0.5)
            mp.CreateNewCondition("PointLoadCondition3D1N", mp.NumberOfConditions() + 1, [node.Id], props)
        node.Fix(KM.ROTATION_X)
        node.Fix(KM.ROTATION_Y)
        node.Fix(KM.ROTATION_Z)

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-10, 1e-12)
    strat = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strat.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strat.Solve()
    assert strat.IsConverged()

    # Monolithic 2-cell reference block
    model_m = KM.Model()
    mp_m = model_m.CreateModelPart("Mono")
    mp_m.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp_m.AddNodalSolutionStepVariable(KM.REACTION)
    mp_m.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp_m.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mono_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 2.0), n_cells=(1, 1, 2), degree=1)
    props_m = mp_m.GetProperties()[1]
    props_m.SetValue(KM.YOUNG_MODULUS, 1000.0)
    props_m.SetValue(KM.POISSON_RATIO, 0.3)
    props_m.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    vol_m, _ = make_nurbs_volume(KM, mp_m, mono_box)
    attach_cell_elements(KM, mp_m, vol_m, mono_box, props_m)
    for node in mp_m.Nodes:
        if abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)
        elif abs(node.Z0 - 2.0) < 1e-12:
            node.SetSolutionStepValue(SM.POINT_LOAD_Z, -0.5)
            mp_m.CreateNewCondition("PointLoadCondition3D1N", mp_m.NumberOfConditions() + 1, [node.Id], props_m)
    mp_m.SetBufferSize(2)
    strat_m = KM.ResidualBasedNewtonRaphsonStrategy(
        mp_m,
        KM.ResidualBasedIncrementalUpdateStaticScheme(),
        KM.ResidualCriteria(1e-10, 1e-12),
        KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
        15,
        True,
        False,
        False,
    )
    strat_m.SetEchoLevel(0)
    mp_m.CloneTimeStep(1.0)
    strat_m.Solve()

    uz_top_tied = np.mean([n.GetSolutionStepValue(KM.DISPLACEMENT_Z) for n in mp.Nodes if 100 <= n.Id < 1000 and abs(n.Z0 - 2.0) < 1e-12])
    uz_top_mono = np.mean([n.GetSolutionStepValue(KM.DISPLACEMENT_Z) for n in mp_m.Nodes if abs(n.Z0 - 2.0) < 1e-12])
    jump_z = max(abs(nl.GetSolutionStepValue(KM.DISPLACEMENT_Z) - nh.GetSolutionStepValue(KM.DISPLACEMENT_Z)) for nl, nh in zip(s_lo, s_hi))
    rz_sum = sum(n.GetSolutionStepValue(KM.REACTION_Z) for n in mp.Nodes if n.Id < 100 and abs(n.Z0) < 1e-12)
    diff_top = abs(uz_top_tied - uz_top_mono)
    print(f"S22 surface spring tie: uz_top_tied={uz_top_tied:.9f}, uz_top_mono={uz_top_mono:.9f}, diff={diff_top:.3e}, jump_z={jump_z:.3e}, sum_Rz={rz_sum:.9f}")
    assert abs(rz_sum - 2.0) < 1e-6
    assert jump_z < 1e-9
    assert diff_top < 1e-9
