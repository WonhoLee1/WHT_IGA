"""S23: Non-matching interface MPC weights via MappingApplication (nearest_element)
and AssignMasterSlaveConstraintsToNeighboursUtility.

Question: Can AssignMasterSlaveConstraintsToNeighboursUtility or MappingApplication's
nearest_element mapper build non-matching interface MPC weights that pass
constant and linear displacement patch tests (target: 1e-12)?

Findings tested directly:
1. MappingApplication ``nearest_element`` mapper produces sparse weights W that
   satisfy partition of unity (err = 0.0) and linear coordinate reproduction
   W @ X_master = X_slave (err ~ 1e-16 < 1e-12).
2. Two non-matching blocks (2x2x1 vs 3x3x1) tied by these MPC weights pass the
   constant displacement patch test (zero stress) to machine precision (~ 1e-16).
3. Kinematic linear displacement transfer (master face prescribed to u=(x,2y,3z),
   slave block driven through W) passes the linear patch test to ~ 1e-15.
4. However, when *both* elastic blocks have free interface DOFs under non-zero
   stress (u=(x,2y,3z), sigma_zz != 0), pure nodal collocation MPC (W) violates
   variational traction consistency (int N_m dA != W^T int N_s dA), leaving a
   measured ~0.0860 displacement error at the master interior node.
5. AssignMasterSlaveConstraintsToNeighboursUtility creates 48 constraints on the
   16 slave nodes x 3 DOFs, and we measure its constant/linear reproduction error.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements, solve_linear
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
    "MappingApplication",
)


def _create_tri_surface_on_box_face(KM, mp, smp, box: BoxSpline, node_map: dict[int, object], z_val: float, first_cid: int):
    """Create Triangle3D3 SurfaceCondition3D3N elements on the z=z_val face using existing control-point nodes."""
    nu, nv, nw = box.n_cp
    X = box.control_points()
    k_match = None
    for k in range(nw):
        I = k * nu * nv
        if abs(X[I, 2] - z_val) < 1e-12:
            k_match = k
            break
    assert k_match is not None
    cprop = mp.GetProperties()[2]
    cid = first_cid
    for j in range(nv - 1):
        for i in range(nu - 1):
            n00 = node_map[i + nu * (j + nv * k_match)]
            n10 = node_map[(i + 1) + nu * (j + nv * k_match)]
            n11 = node_map[(i + 1) + nu * ((j + 1) + nv * k_match)]
            n01 = node_map[i + nu * ((j + 1) + nv * k_match)]
            for n in (n00, n10, n11, n01):
                if not smp.HasNode(n.Id):
                    smp.AddNode(n, 0)
            smp.CreateNewCondition("SurfaceCondition3D3N", cid, [n00.Id, n10.Id, n11.Id], cprop)
            smp.CreateNewCondition("SurfaceCondition3D3N", cid + 1, [n00.Id, n11.Id, n01.Id], cprop)
            cid += 2


def _build_nonmatching_model(fix_master_interface: bool = False):
    """Build two non-matching p=1 blocks (2x2x1 vs 3x3x1) tied at z=1 via nearest_element mapper."""
    mods = load_kratos(_APPS)
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    import KratosMultiphysics.scipy_conversion_tools as sct

    model = KM.Model()
    mp = model.CreateModelPart("S23")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    lo_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 1), degree=1)
    hi_box = BoxSpline(lo=(0.0, 0.0, 1.0), hi=(1.0, 1.0, 2.0), n_cells=(3, 3, 1), degree=1)

    vol_lo, _ = make_nurbs_volume(KM, mp, lo_box, first_node_id=1)
    attach_cell_elements(KM, mp, vol_lo, lo_box, props, start_id=1)
    lo_nodes = {n.Id - 1: n for n in mp.Nodes if n.Id <= lo_box.n_ctrl}

    vol_hi, _ = make_nurbs_volume(KM, mp, hi_box, first_node_id=100)
    attach_cell_elements(KM, mp, vol_hi, hi_box, props, start_id=100)
    hi_nodes = {n.Id - 100: n for n in mp.Nodes if n.Id >= 100}

    smp_master = mp.CreateSubModelPart("MasterFace")
    smp_slave = mp.CreateSubModelPart("SlaveFace")
    _create_tri_surface_on_box_face(KM, mp, smp_master, lo_box, lo_nodes, 1.0, first_cid=1000)
    _create_tri_surface_on_box_face(KM, mp, smp_slave, hi_box, hi_nodes, 1.0, first_cid=2000)

    mapper = KM.MapperFactory.CreateMapper(
        smp_master,
        smp_slave,
        KM.Parameters(json.dumps({"mapper_type": "nearest_element"})),
    )
    W = sct.to_csr(mapper.GetMappingMatrix())

    master_ordered = list(smp_master.Nodes)
    slave_ordered = list(smp_slave.Nodes)
    xm = np.array([[n.X0, n.Y0, n.Z0] for n in master_ordered])
    xs = np.array([[n.X0, n.Y0, n.Z0] for n in slave_ordered])
    pou_err = float(np.max(np.abs(W @ np.ones(len(master_ordered)) - 1.0)))
    coord_err = float(np.max(np.abs(W @ xm - xs)))

    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    cid = 1
    for r, snode in enumerate(slave_ordered):
        row = W.getrow(r)
        keep = np.abs(row.data) > 1e-14
        m_nodes = [master_ordered[c] for c in row.indices[keep]]
        vals = row.data[keep]
        for var in vars3:
            add_linear_tie(
                KM,
                mp,
                cid,
                [snode.GetDof(var)],
                [mn.GetDof(var) for mn in m_nodes],
                vals.reshape(1, -1),
                np.zeros(1),
            )
            cid += 1

    return KM, mp, {n.Id for n in slave_ordered}, {n.Id for n in master_ordered}, pou_err, coord_err, (vol_lo, vol_hi)


def _solve_patch(field_fn, fix_master_interface: bool = False) -> tuple[float, float, float]:
    KM, mp, slave_ids, master_ids, pou_err, coord_err, _vols = _build_nonmatching_model(fix_master_interface)
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for node in mp.Nodes:
        if node.Id in slave_ids:
            continue
        on_outer = (
            abs(node.X0) < 1e-12
            or abs(node.X0 - 1.0) < 1e-12
            or abs(node.Y0) < 1e-12
            or abs(node.Y0 - 1.0) < 1e-12
            or abs(node.Z0) < 1e-12
            or abs(node.Z0 - 2.0) < 1e-12
        )
        if on_outer or (fix_master_interface and node.Id in master_ids):
            u_exact = field_fn(node.X0, node.Y0, node.Z0)
            for var, val in zip(vars3, u_exact):
                node.Fix(var)
                node.SetSolutionStepValue(var, float(val))

    solve_linear(KM, mp)
    max_err = 0.0
    for node in mp.Nodes:
        u_num = np.array(node.GetSolutionStepValue(KM.DISPLACEMENT))
        u_exact = np.array(field_fn(node.X0, node.Y0, node.Z0))
        max_err = max(max_err, float(np.max(np.abs(u_num - u_exact))))
    return pou_err, coord_err, max_err


def test_s23_nearest_element_weights_completeness() -> None:
    """nearest_element mapping matrix W satisfies partition of unity and linear completeness (< 1e-12)."""
    _, _, _, _, pou_err, coord_err, _ = _build_nonmatching_model()
    print(f"S23 nearest_element completeness: pou_err={pou_err:.3e}, coord_err={coord_err:.3e}")
    assert pou_err < 1e-12
    assert coord_err < 1e-12


def test_s23_nearest_element_constant_patch() -> None:
    """Constant displacement patch u=(0.1, -0.2, 0.3) across non-matching 2x2 / 3x3 interface."""
    pou_err, coord_err, max_err = _solve_patch(lambda x, y, z: (0.1, -0.2, 0.3), fix_master_interface=False)
    print(f"S23 constant patch: pou_err={pou_err:.3e}, coord_err={coord_err:.3e}, max_u_err={max_err:.3e}")
    assert max_err < 1e-12


def test_s23_nearest_element_kinematic_linear_patch() -> None:
    """Kinematic linear patch u=(x, 2y, 3z) with master interface prescribed and slave driven by W."""
    pou_err, coord_err, max_err = _solve_patch(lambda x, y, z: (x, 2.0 * y, 3.0 * z), fix_master_interface=True)
    print(f"S23 kinematic linear patch: pou_err={pou_err:.3e}, coord_err={coord_err:.3e}, max_u_err={max_err:.3e}")
    assert max_err < 1e-12


def test_s23_nearest_element_coupled_linear_patch_collocation_gap() -> None:
    """Free-interface coupled linear patch under non-zero stress exposes collocation non-variational gap."""
    _, _, max_err = _solve_patch(lambda x, y, z: (x, 2.0 * y, 3.0 * z), fix_master_interface=False)
    print(f"S23 coupled free-interface linear patch (non-variational collocation gap): max_u_err={max_err:.6f}")
    # Documents the collocation traction mismatch (S23.md: 0.085987); band guards the recorded value.
    assert 0.08 < max_err < 0.09


def test_s23_assign_neighbours_utility_weights() -> None:
    """Measure AssignMasterSlaveConstraintsToNeighboursUtility on the non-matching interface."""
    mods = load_kratos(_APPS)
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    model = KM.Model()
    mp = model.CreateModelPart("S23Util")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    lo_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 1), degree=1)
    hi_box = BoxSpline(lo=(0.0, 0.0, 1.0), hi=(1.0, 1.0, 2.0), n_cells=(3, 3, 1), degree=1)
    vol_lo, _ = make_nurbs_volume(KM, mp, lo_box, first_node_id=1)
    attach_cell_elements(KM, mp, vol_lo, lo_box, props, start_id=1)
    vol_hi, _ = make_nurbs_volume(KM, mp, hi_box, first_node_id=100)
    attach_cell_elements(KM, mp, vol_hi, hi_box, props, start_id=100)

    smp_m = mp.CreateSubModelPart("MasterFace")
    smp_s = mp.CreateSubModelPart("SlaveFace")
    for n in mp.Nodes:
        if abs(n.Z0 - 1.0) < 1e-12:
            if n.Id < 100:
                smp_m.AddNode(n, 0)
            else:
                smp_s.AddNode(n, 0)

    util = KM.AssignMasterSlaveConstraintsToNeighboursUtility(smp_m.Nodes)
    util.AssignMasterSlaveConstraintsToNodes(
        smp_s.Nodes,
        0.75,
        mp,
        [KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z],
        1e-9,
    )
    n_cons = mp.NumberOfMasterSlaveConstraints()

    slave_ids = {n.Id for n in smp_s.Nodes}
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for node in mp.Nodes:
        if node.Id in slave_ids:
            continue
        on_outer = (
            abs(node.X0) < 1e-12
            or abs(node.X0 - 1.0) < 1e-12
            or abs(node.Y0) < 1e-12
            or abs(node.Y0 - 1.0) < 1e-12
            or abs(node.Z0) < 1e-12
            or abs(node.Z0 - 2.0) < 1e-12
            or abs(node.Z0 - 1.0) < 1e-12
        )
        if on_outer:
            for var, val in zip(vars3, (node.X0, 2.0 * node.Y0, 3.0 * node.Z0)):
                node.Fix(var)
                node.SetSolutionStepValue(var, float(val))

    solve_linear(KM, mp)
    lin_err = max(
        float(np.max(np.abs(np.array(n.GetSolutionStepValue(KM.DISPLACEMENT)) - np.array([n.X0, 2.0 * n.Y0, 3.0 * n.Z0]))))
        for n in smp_s.Nodes
    )
    print(f"S23 AssignMasterSlaveConstraintsToNeighboursUtility: n_cons={n_cons}, kinematic_linear_err={lin_err:.3e}")
    assert n_cons == 48
    # Decision evidence (S23.md: 7.942e-2): the utility misses the 1e-12 linear patch, so the
    # mapper route is chosen. Band guards the recorded value.
    assert 0.075 < lin_err < 0.085
