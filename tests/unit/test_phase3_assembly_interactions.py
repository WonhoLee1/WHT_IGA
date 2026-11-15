"""Unit tests (Kratos-free) for Track P3-C assembly, interactions, and Native FE helpers."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.grid.hierarchy import THBHierarchy
from wht_iga.kratos_backend.assembly_builder import (
    AssemblySpec,
    IdOffsetTable,
    InstanceSpec,
    PartSpec,
    compose_thb_operator,
    transform_stress_voigt,
    validate_assembly_spec,
    validate_rotation_matrix,
)
from wht_iga.kratos_backend.interactions import (
    LinearRelation,
    RBE2Spec,
    RBE3Spec,
    SurfaceMesh,
    TieSpec,
    build_rbe2,
    build_rbe3,
    build_tie,
    compile_constraint_graph,
    distribute_rbe3_loads,
    validate_constraint_graph,
)
from wht_iga.kratos_backend.native_fe import (
    check_element_jacobian,
    ensure_positive_jacobian,
    to_kratos_connectivity,
)
from wht_iga.reference.box import BoxSpline


def _make_grid_tri_surface(
    nx: int,
    ny: int,
    z_val: float = 0.0,
    diag_flip: bool = False,
    start_node_id: int = 1,
) -> SurfaceMesh:
    """Create a triangulated unit square surface [0, 1] x [0, 1] x {z_val}."""
    xs = np.linspace(0.0, 1.0, nx + 1)
    ys = np.linspace(0.0, 1.0, ny + 1)
    coords = np.array(
        [[x, y, z_val] for y in ys for x in xs],
        dtype=np.float64,
    )
    node_ids = np.arange(start_node_id, start_node_id + len(coords), dtype=np.int64)
    tris: list[tuple[int, int, int]] = []
    for j in range(ny):
        for i in range(nx):
            n00 = i + (nx + 1) * j
            n10 = (i + 1) + (nx + 1) * j
            n01 = i + (nx + 1) * (j + 1)
            n11 = (i + 1) + (nx + 1) * (j + 1)
            if not diag_flip:
                tris.append((n00, n10, n11))
                tris.append((n00, n11, n01))
            else:
                tris.append((n00, n10, n01))
                tris.append((n10, n11, n01))
    return SurfaceMesh(
        name=f"grid_{nx}x{ny}",
        kind="fe",
        coords=coords,
        triangles=np.asarray(tris, dtype=np.int64),
        node_ids=node_ids,
        dof_node_ids=node_ids,
    )


def test_id_offset_table_monotonic_and_utf8_roundtrip(tmp_path: Path) -> None:
    """ID allocator enforces non-overlapping ranges and preserves UTF-8 JSON roundtrip."""
    tbl = IdOffsetTable(job_name="트랙_P3_조립")
    n1 = tbl.allocate_nodes("Parts.beam_a", 12, provenance="control_point")
    n2 = tbl.allocate_nodes("Parts.beam_b", 8, provenance="fe_node")
    e1 = tbl.allocate_elements("Parts.beam_a", 4, kind="SmallDisplacementElement3D4N")
    c1 = tbl.allocate_constraints("Couplings.tie1", 9)
    p1 = tbl.allocate_properties("Parts.beam_a", 2, provenance="bimetallic")

    assert (n1.start_id, n1.end_id) == (1, 12)
    assert (n2.start_id, n2.end_id) == (13, 20)
    assert (e1.start_id, e1.end_id) == (1, 4)
    assert (c1.start_id, c1.end_id) == (1, 9)
    assert (p1.start_id, p1.end_id) == (1, 2)
    tbl.metadata["description"] = "한글 메타데이터 검증"

    out_file = tbl.to_json(tmp_path / "job.resolved.json")
    raw = json.loads(out_file.read_text(encoding="utf-8"))
    assert raw["job_name"] == "트랙_P3_조립"
    reloaded = IdOffsetTable.from_json(out_file)
    assert reloaded.next_id("nodes") == 21
    assert reloaded.metadata["description"] == "한글 메타데이터 검증"

    # Inject intentional overlap and verify validator rejects it.
    tbl.allocate_nodes("Parts.bad", 5)
    from wht_iga.kratos_backend.assembly_builder import IdRangeEntry

    tbl.entries.append(IdRangeEntry("nodes", "Parts.overlap", 10, 15, 6, "fe_node"))
    with pytest.raises(ValueError, match="Overlapping"):
        tbl.validate_non_overlapping()


def test_assembly_spec_validation_and_transforms() -> None:
    """Validate part ownership rules, duplicate rejection, and Voigt tensor rotation."""
    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 1.0, 1.0), (2, 1, 1), degree=2)
    part = PartSpec(name="beam", kind="iga", box=box)

    ang = np.deg2rad(35.0)
    R = np.array(
        [[np.cos(ang), -np.sin(ang), 0.0], [np.sin(ang), np.cos(ang), 0.0], [0.0, 0.0, 1.0]]
    )
    inst1 = InstanceSpec("beam", "left", dependent=True, translation=(1.0, 2.0, 3.0), rotation_matrix=R)
    inst2 = InstanceSpec("beam", "right", dependent=True, translation=(5.0, 0.0, 0.0))
    validate_assembly_spec(AssemblySpec("ok", {"beam": part}, [inst1, inst2]))

    # Voigt rotation matches direct 3x3 R S R^T
    s_loc = np.array([120.0, -40.0, 15.0, 30.0, -12.0, 25.0])
    s_rot = transform_stress_voigt(s_loc, R)
    S_mat = np.array([[120.0, 30.0, 25.0], [30.0, -40.0, -12.0], [25.0, -12.0, 15.0]])
    S_exp = R @ S_mat @ R.T
    s_exp = np.array([S_exp[0, 0], S_exp[1, 1], S_exp[2, 2], S_exp[0, 1], S_exp[1, 2], S_exp[0, 2]])
    assert np.allclose(s_rot, s_exp, atol=1e-12)

    # Improper reflection rejected
    R_reflect = np.diag([-1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match="det="):
        validate_rotation_matrix(R_reflect)

    # Mixed dependent + independent on same part rejected
    bad_mix = AssemblySpec(
        "mix",
        {"beam": part},
        [inst1, InstanceSpec("beam", "indep", dependent=False)],
    )
    with pytest.raises(ValueError, match="mixes dependent and independent"):
        validate_assembly_spec(bad_mix)

    # Dependent instance cannot override grid
    bad_override = AssemblySpec(
        "ovr",
        {"beam": part},
        [InstanceSpec("beam", "bad", dependent=True, box_override=box)],
    )
    with pytest.raises(ValueError, match="override"):
        validate_assembly_spec(bad_override)


def test_compose_thb_operator_partition_of_unity() -> None:
    """Rectangular C_THB and E compose accurately with partition of unity == 1.0."""
    thb = THBHierarchy.from_uniform((0.0, 4.0), (0.0, 1.0), (0.0, 1.0), (4, 1, 1), degree=2)
    thb.mark_box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), target_level=1)
    thb.enforce_2_to_1_balance()
    C_thb = thb.build_c_thb_matrix()
    n_thb = C_thb.shape[1]

    # Non-trivial column-preserving extension mock test on exact identity subspace
    E_id = sp.eye(n_thb, format="csr")
    T = compose_thb_operator(E_id, C_thb)
    assert T is not None
    assert T.shape == C_thb.shape
    row_sums = np.asarray(T.sum(axis=1)).ravel()
    assert np.allclose(row_sums, 1.0, atol=1e-13)


def test_constraint_graph_validator_and_chain_flattening() -> None:
    """Enforce S04 c==0, double-slave, cycle rejection, and exact chain flattening."""
    # 1. S04 nonzero constant rejected
    with pytest.raises(ValueError, match="S04 violation"):
        compile_constraint_graph(
            [LinearRelation((2, "DISPLACEMENT_X"), {(1, "DISPLACEMENT_X"): 1.0}, constant=0.5)]
        )

    # 2. Double slave rejected
    with pytest.raises(ValueError, match="Double-slave"):
        compile_constraint_graph(
            [
                LinearRelation((3, "DISPLACEMENT_X"), {(1, "DISPLACEMENT_X"): 1.0}),
                LinearRelation((3, "DISPLACEMENT_X"), {(2, "DISPLACEMENT_X"): 1.0}),
            ]
        )

    # 3. Cycle rejected even when flattening enabled
    with pytest.raises(ValueError, match="Cyclic"):
        compile_constraint_graph(
            [
                LinearRelation((1, "DISPLACEMENT_X"), {(2, "DISPLACEMENT_X"): 0.5}),
                LinearRelation((2, "DISPLACEMENT_X"), {(3, "DISPLACEMENT_X"): 1.0}),
                LinearRelation((3, "DISPLACEMENT_X"), {(1, "DISPLACEMENT_X"): 0.25}),
            ],
            chain_policy="flatten",
        )

    # 4. Chained relation s -> m -> r
    chained = [
        LinearRelation(
            (10, "DISPLACEMENT_X"),
            {(20, "DISPLACEMENT_X"): 0.4, (21, "DISPLACEMENT_X"): 0.6},
        ),
        LinearRelation(
            (20, "DISPLACEMENT_X"),
            {(30, "DISPLACEMENT_X"): 0.25, (31, "DISPLACEMENT_X"): 0.75},
        ),
        LinearRelation(
            (21, "DISPLACEMENT_X"),
            {(31, "DISPLACEMENT_X"): 0.5, (32, "DISPLACEMENT_X"): 0.5},
        ),
    ]
    with pytest.raises(ValueError, match="chained"):
        compile_constraint_graph(chained, chain_policy="reject")

    flat = compile_constraint_graph(chained, chain_policy="flatten")
    rep = validate_constraint_graph(flat, allow_chains=False)
    assert rep.valid
    row10 = next(r for r in flat if r.slave_dof == (10, "DISPLACEMENT_X"))
    # Exact AB composition: 0.4*0.25=0.10, 0.4*0.75 + 0.6*0.5=0.60, 0.6*0.5=0.30
    assert abs(row10.terms[(30, "DISPLACEMENT_X")] - 0.10) < 1e-14
    assert abs(row10.terms[(31, "DISPLACEMENT_X")] - 0.60) < 1e-14
    assert abs(row10.terms[(32, "DISPLACEMENT_X")] - 0.30) < 1e-14


def test_tie_matching_mortar_and_penalty_patch_accuracy() -> None:
    """Verify matching, clipped L2 mortar virtual-work patch, and FE-IGA penalty tie."""
    # Non-matching unit squares with different counts and opposite diagonals
    master = _make_grid_tri_surface(3, 3, z_val=0.0, diag_flip=False, start_node_id=1)
    slave = _make_grid_tri_surface(5, 4, z_val=0.0, diag_flip=True, start_node_id=100)

    # Arbitrary linear displacement field u(x, y) = 2.0 + 0.5 x - 0.3 y
    def u_lin(pts: np.ndarray) -> np.ndarray:
        return 2.0 + 0.5 * pts[:, 0] - 0.3 * pts[:, 1]

    u_m = u_lin(master.coords)
    u_s_exact = u_lin(slave.coords)

    # 1. mpc_matching fallback barycentric projection
    res_mpc = build_tie(TieSpec("bary", master=master, slave=slave, method="mpc_matching"))
    assert res_mpc.max_gap < 1e-13
    u_s_mpc = np.zeros(len(slave.coords))
    s_lookup = {int(nid): idx for idx, nid in enumerate(slave.node_ids)}
    m_lookup = {int(nid): idx for idx, nid in enumerate(master.node_ids)}
    for rel in res_mpc.relations:
        if rel.slave_dof[1] == "DISPLACEMENT_X":
            si = s_lookup[rel.slave_dof[0]]
            u_s_mpc[si] = sum(w * u_m[m_lookup[k[0]]] for k, w in rel.terms.items())
    assert np.max(np.abs(u_s_mpc - u_s_exact)) < 1e-13

    # 2. Clipped L2 mortar projection: exact partition of unity + linear patch + virtual work!
    res_mortar = build_tie(TieSpec("mortar", master=master, slave=slave, method="mortar_l2", gauss_order=3))
    assert res_mortar.resolved_method == "mortar_l2"
    u_s_mortar = np.zeros(len(slave.coords))
    P_mat = np.zeros((len(slave.coords), len(master.coords)))
    for rel in res_mortar.relations:
        if rel.slave_dof[1] == "DISPLACEMENT_X":
            si = s_lookup[rel.slave_dof[0]]
            for k, w in rel.terms.items():
                P_mat[si, m_lookup[k[0]]] = w
            u_s_mortar[si] = float(P_mat[si] @ u_m)
    assert np.max(np.abs(P_mat.sum(axis=1) - 1.0)) < 1e-12
    assert np.max(np.abs(u_s_mortar - u_s_exact)) < 1e-12

    # 3. Penalty spring tie between FE surface and IGA surface on z=0
    box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2, 2, 1), degree=2)
    iga_nids = np.arange(500, 500 + box.n_ctrl, dtype=np.int64)
    iga_face = SurfaceMesh.from_iga_face("iga_z0", box, axis=2, side="lo", dof_node_ids=iga_nids)
    res_pen = build_tie(
        TieSpec("fe_iga", master=iga_face, slave=slave, method="penalty_spring", E_max=1e5)
    )
    assert res_pen.resolved_method == "penalty_spring"
    assert res_pen.penalty_matrix is not None
    K = res_pen.penalty_matrix
    assert np.max(np.abs((K - K.T).toarray())) < 1e-8 * np.max(np.abs(K.data))
    # Linear displacement on BOTH sides must produce zero interface penalty energy
    all_nids = res_pen.penalty_node_ids
    cp_coords = box.control_points()
    u_all = np.zeros(3 * len(all_nids))
    for pos, nid in enumerate(all_nids):
        if nid in s_lookup:
            p = slave.coords[s_lookup[nid]]
        else:
            p = cp_coords[int(nid) - 500]
        u_all[3 * pos : 3 * pos + 3] = [u_lin(p[None, :])[0], 0.5 * p[2] - 0.1 * p[0], 0.3 * p[1]]
    res_vec = K @ u_all
    pen_energy = float(u_all @ res_vec)
    scale = float(K.diagonal().max())
    assert np.max(np.abs(res_vec)) / scale < 1e-13
    assert abs(pen_energy) / scale < 1e-12


def test_rbe2_and_rbe3_kinematics_and_static_equivalence() -> None:
    """Verify RBE2 rigid motion and RBE3 exact F/M virtual-work distribution without stiffening."""
    xs = np.linspace(-0.05, 0.05, 4)
    ys = np.linspace(-0.05, 0.05, 4)
    coords = np.array([[x, y, 1.0] for y in ys for x in xs], dtype=np.float64)
    nids = np.arange(1, len(coords) + 1, dtype=np.int64)
    ref_coord = (0.0, 0.0, 1.0)

    # Arbitrary rigid motion
    u0 = np.array([0.02, -0.01, 0.03])
    omega = np.array([0.004, -0.003, 0.002])

    # RBE2
    r2 = build_rbe2(RBE2Spec("r2", ref_node_id=999, ref_coord=ref_coord, slave_node_ids=nids, slave_coords=coords))
    q_ref = np.concatenate([u0, omega])
    u_r2 = (r2.operator @ q_ref).reshape(-1, 3)
    u_exact = np.array([u0 + np.cross(omega, p - np.asarray(ref_coord)) for p in coords])
    assert np.max(np.abs(u_r2 - u_exact)) < 1e-15

    # RBE3 (both uniform and inverse_distance modes)
    for wmode in ("uniform", "inverse_distance"):
        r3 = build_rbe3(
            RBE3Spec(
                f"r3_{wmode}",
                ref_node_id=999,
                ref_coord=ref_coord,
                slave_coords=coords,
                slave_node_ids=nids,
                weight_mode=wmode,
            )
        )
        q_rec = r3.operator @ u_exact.ravel()
        assert np.max(np.abs(q_rec - q_ref)) < 1e-13

        # Virtual-work force/moment equilibrium check: f = W^T [F; M]
        F_app = (12.0, -7.5, 40.0)
        M_app = (5.0, -2.5, 1.8)
        nodal_f = distribute_rbe3_loads(r3, force=F_app, moment=M_app)
        F_sum = nodal_f.sum(axis=0)
        M_sum = sum(
            np.cross(p - np.asarray(ref_coord), fi) for p, fi in zip(coords, nodal_f)
        )
        assert np.max(np.abs(F_sum - np.asarray(F_app))) < 1e-12
        assert np.max(np.abs(M_sum - np.asarray(M_app))) < 1e-12


def test_native_fe_jacobian_check_and_auto_flip_all_elements() -> None:
    """Test Jacobian sign checks and inverted orientation auto-flip across all 5 element types."""
    hex8_pts = np.array(
        [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ],
        dtype=np.float64,
    )
    hex20_pts = np.vstack(
        [
            hex8_pts,
            [
                [0.5, 0, 0], [1, 0.5, 0], [0.5, 1, 0], [0, 0.5, 0],
                [0.5, 0, 1], [1, 0.5, 1], [0.5, 1, 1], [0, 0.5, 1],
                [0, 0, 0.5], [1, 0, 0.5], [1, 1, 0.5], [0, 1, 0.5],
            ],
        ]
    )
    tet4_pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    tet10_pts = np.vstack(
        [tet4_pts, [[0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0, 0.5], [0, 0.5, 0.5]]]
    )
    wedge6_pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]],
        dtype=np.float64,
    )

    cases = [
        ("hex8", hex8_pts, np.array([4, 5, 6, 7, 0, 1, 2, 3])),
        (
            "hex20",
            hex20_pts,
            np.array([4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 16, 17, 18, 19]),
        ),
        ("tet4", tet4_pts, np.array([0, 2, 1, 3])),
        ("tet10", tet10_pts, np.array([0, 2, 1, 3, 6, 5, 4, 7, 9, 8])),
        ("wedge6", wedge6_pts, np.array([3, 4, 5, 0, 1, 2])),
    ]

    for etype, coords, inv_perm in cases:
        valid_conn = np.arange(len(coords), dtype=np.int64)
        min_d, _ = check_element_jacobian(etype, coords[valid_conn])
        assert min_d > 0.0

        inv_conn = valid_conn[inv_perm]
        inv_min_d, _ = check_element_jacobian(etype, coords[inv_conn])
        assert inv_min_d < 0.0

        with pytest.raises(ValueError, match="Inverted"):
            ensure_positive_jacobian(etype, inv_conn, coords, auto_flip=False)

        fixed_conn, was_flipped = ensure_positive_jacobian(etype, inv_conn, coords, auto_flip=True)
        assert was_flipped
        fixed_min_d, _ = check_element_jacobian(etype, coords[fixed_conn])
        assert fixed_min_d > 0.0

        # Check midside midpoint geometric consistency on flipped quadratic elements
        if etype == "tet10":
            c = coords[fixed_conn]
            assert np.allclose(c[4], 0.5 * (c[0] + c[1]))
            assert np.allclose(c[5], 0.5 * (c[1] + c[2]))
            assert np.allclose(c[6], 0.5 * (c[2] + c[0]))
            assert np.allclose(c[8], 0.5 * (c[1] + c[3]))
            assert np.allclose(c[9], 0.5 * (c[2] + c[3]))
        elif etype == "hex20":
            c = coords[fixed_conn]
            assert np.allclose(c[8], 0.5 * (c[0] + c[1]))
            assert np.allclose(c[12], 0.5 * (c[4] + c[5]))
            assert np.allclose(c[16], 0.5 * (c[0] + c[4]))
            kr = to_kratos_connectivity("hex20", fixed_conn)
            assert kr[12] == fixed_conn[16] and kr[16] == fixed_conn[12]
