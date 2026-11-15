"""Kratos integration tests for Track P3-C (THB, multi-cell, dependent instances,
mortar L2 Tie, RBE2, RBE3 V03-moment, Native FE 5 element types, and IGA-FE hybrid).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from wht_iga.grid.bspline import basis_1d
from wht_iga.grid.hierarchy import THBHierarchy
from wht_iga.kratos_backend.assembly_builder import (
    AssemblyBuilder,
    AssemblySpec,
    CellDomainSpec,
    InstanceSpec,
    PartSpec,
    transform_stress_voigt,
)
from wht_iga.kratos_backend.elements import solve_linear
from wht_iga.kratos_backend.interactions import (
    RBE2Spec,
    RBE3Spec,
    SurfaceMesh,
    TieSpec,
    apply_rbe2,
    apply_rbe3,
    apply_tie,
)
from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe
from wht_iga.meshio_ext.model import ElemBlock, MeshModel
from wht_iga.reference.box import BoxSpline


_DATA_DIR = Path(__file__).resolve().parents[1] / "verification" / "data"


def _integrate_iga_face_traction(
    box: BoxSpline,
    axis: int,
    side: str,
    traction: tuple[float, float, float],
    transform_vec: np.ndarray | None = None,
) -> dict[int, np.ndarray]:
    """Compute consistent nodal traction loads sum_q N_I(x_q) * t * w_q on a box face."""
    in_axes = [a for a in (0, 1, 2) if a != axis]
    a0, a1 = in_axes[0], in_axes[1]
    bps0 = box.breakpoints(a0)
    bps1 = box.breakpoints(a1)
    p = int(box.degree)
    gx, gw = np.polynomial.legendre.leggauss(p + 1)
    nu, nv, nw = (int(v) for v in box.n_cp)
    ku, kv, kw = box.knots
    c_fixed = float(box.lo[axis] if side == "lo" else box.hi[axis])
    t_vec = np.asarray(traction, dtype=np.float64)
    if transform_vec is not None:
        t_vec = transform_vec @ t_vec

    loads: dict[int, np.ndarray] = {}
    knots_list = [ku, kv, kw]
    s_fix, n_fix_arr = basis_1d(knots_list[axis], p, c_fixed, 0)
    i_fix0 = s_fix - p
    n_fix = n_fix_arr[0]

    for j_cell in range(len(bps1) - 1):
        v0, v1 = float(bps1[j_cell]), float(bps1[j_cell + 1])
        pts1 = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gx
        wts1 = 0.5 * (v1 - v0) * gw
        for i_cell in range(len(bps0) - 1):
            u0, u1 = float(bps0[i_cell]), float(bps0[i_cell + 1])
            pts0 = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gx
            wts0 = 0.5 * (u1 - u0) * gw
            for q1, w1 in zip(pts1, wts1):
                s1, n1_arr = basis_1d(knots_list[a1], p, float(q1), 0)
                i1_0 = s1 - p
                n1 = n1_arr[0]
                for q0, w0 in zip(pts0, wts0):
                    s0, n0_arr = basis_1d(knots_list[a0], p, float(q0), 0)
                    i0_0 = s0 - p
                    n0 = n0_arr[0]
                    w_area = float(w0 * w1)
                    for lf in range(p + 1):
                        if abs(float(n_fix[lf])) < 1e-16:
                            continue
                        idx_fix = i_fix0 + lf
                        for l1 in range(p + 1):
                            idx1 = i1_0 + l1
                            for l0 in range(p + 1):
                                idx0 = i0_0 + l0
                                weight = float(n_fix[lf] * n1[l1] * n0[l0]) * w_area
                                coords = [0, 0, 0]
                                coords[axis] = idx_fix
                                coords[a0] = idx0
                                coords[a1] = idx1
                                flat = coords[0] + nu * (coords[1] + nv * coords[2])
                                loads[flat] = loads.get(flat, np.zeros(3)) + weight * t_vec
    return loads


def _apply_iga_loads(
    KM: object,
    SM: object,
    model_part: object,
    node_ids: list[int],
    loads_by_flat: dict[int, np.ndarray],
    props: object,
    start_cond_id: int = 1,
) -> int:
    """Attach consistent nodal forces as PointLoadCondition3D1N."""
    cid = int(start_cond_id)
    for flat, fvec in loads_by_flat.items():
        if np.linalg.norm(fvec) <= 1e-18:
            continue
        nid = int(node_ids[flat])
        nd = model_part.GetNode(nid)
        cur = np.array(nd.GetSolutionStepValue(SM.POINT_LOAD, 0), dtype=np.float64)
        nd.SetSolutionStepValue(SM.POINT_LOAD, 0, [float(x) for x in (cur + fvec)])
        model_part.CreateNewCondition("PointLoadCondition3D1N", cid, [nid], props)
        cid += 1
    return cid


def _make_structured_hex_mesh(
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    ncells: tuple[int, int, int],
    etype: str = "hex8",
    part_id: int = 1,
) -> MeshModel:
    """Build structured MeshModel of hex8 or hex20 (or split into tet4/wedge6/tet10)."""
    nx, ny, nz = ncells
    x0, y0, z0 = lo
    x1, y1, z1 = hi

    if etype == "hex8":
        xs = np.linspace(x0, x1, nx + 1)
        ys = np.linspace(y0, y1, ny + 1)
        zs = np.linspace(z0, z1, nz + 1)
        nodes = np.array(
            [[x, y, z] for z in zs for y in ys for x in xs], dtype=np.float64
        )
        def _nid(i: int, j: int, k: int) -> int:
            return i + (nx + 1) * (j + (ny + 1) * k)
        conn = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    conn.append([
                        _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                        _nid(i, j, k + 1), _nid(i + 1, j, k + 1), _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                    ])
        return MeshModel(nodes=nodes, blocks=[ElemBlock(etype, np.asarray(conn, dtype=np.int64), part_id)])

    if etype == "hex20":
        hx, hy, hz = 2 * nx, 2 * ny, 2 * nz
        xs = np.linspace(x0, x1, hx + 1)
        ys = np.linspace(y0, y1, hy + 1)
        zs = np.linspace(z0, z1, hz + 1)
        coord_to_id: dict[tuple[int, int, int], int] = {}
        nodes_list: list[list[float]] = []
        def _pt(i: int, j: int, k: int) -> int:
            key = (i, j, k)
            if key not in coord_to_id:
                coord_to_id[key] = len(nodes_list)
                nodes_list.append([float(xs[i]), float(ys[j]), float(zs[k])])
            return coord_to_id[key]

        conn = []
        for kz in range(nz):
            for jy in range(ny):
                for ix in range(nx):
                    cx, cy, cz = 2 * ix, 2 * jy, 2 * kz
                    row = [
                        _pt(cx, cy, cz), _pt(cx + 2, cy, cz), _pt(cx + 2, cy + 2, cz), _pt(cx, cy + 2, cz),
                        _pt(cx, cy, cz + 2), _pt(cx + 2, cy, cz + 2), _pt(cx + 2, cy + 2, cz + 2), _pt(cx, cy + 2, cz + 2),
                        _pt(cx + 1, cy, cz), _pt(cx + 2, cy + 1, cz), _pt(cx + 1, cy + 2, cz), _pt(cx, cy + 1, cz),
                        _pt(cx + 1, cy, cz + 2), _pt(cx + 2, cy + 1, cz + 2), _pt(cx + 1, cy + 2, cz + 2), _pt(cx, cy + 1, cz + 2),
                        _pt(cx, cy, cz + 1), _pt(cx + 2, cy, cz + 1), _pt(cx + 2, cy + 2, cz + 1), _pt(cx, cy + 2, cz + 1),
                    ]
                    conn.append(row)
        return MeshModel(
            nodes=np.asarray(nodes_list, dtype=np.float64),
            blocks=[ElemBlock("hex20", np.asarray(conn, dtype=np.int64), part_id)],
        )

    base = _make_structured_hex_mesh(lo, hi, ncells, "hex8", part_id)
    hconn = base.blocks[0].conn
    if etype == "wedge6":
        wconn = []
        for c in hconn:
            wconn.append([c[0], c[1], c[2], c[4], c[5], c[6]])
            wconn.append([c[0], c[2], c[3], c[4], c[6], c[7]])
        return MeshModel(nodes=base.nodes, blocks=[ElemBlock("wedge6", np.asarray(wconn, dtype=np.int64), part_id)])

    if etype in ("tet4", "tet10"):
        # Kuhn 6-tet triangulation of structured hexahedral cells
        tets4 = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    def _n(ii: int, jj: int, kk: int) -> int:
                        return (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk))
                    v000 = _n(0, 0, 0)
                    v100 = _n(1, 0, 0)
                    v010 = _n(0, 1, 0)
                    v110 = _n(1, 1, 0)
                    v001 = _n(0, 0, 1)
                    v101 = _n(1, 0, 1)
                    v011 = _n(0, 1, 1)
                    v111 = _n(1, 1, 1)
                    raw_tets = [
                        [v000, v100, v110, v111],
                        [v000, v110, v010, v111],
                        [v000, v010, v011, v111],
                        [v000, v011, v001, v111],
                        [v000, v001, v101, v111],
                        [v000, v101, v100, v111],
                    ]
                    pts = base.nodes
                    for rt in raw_tets:
                        J = np.column_stack([pts[rt[1]] - pts[rt[0]], pts[rt[2]] - pts[rt[0]], pts[rt[3]] - pts[rt[0]]])
                        if np.linalg.det(J) < 0:
                            rt = [rt[0], rt[2], rt[1], rt[3]]
                        tets4.append(rt)
        if etype == "tet4":
            return MeshModel(nodes=base.nodes, blocks=[ElemBlock("tet4", np.asarray(tets4, dtype=np.int64), part_id)])

        # Subdivide into tet10 by inserting unique edge midpoints
        nodes_list = base.nodes.tolist()
        edge_mid: dict[tuple[int, int], int] = {}
        def _mid(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edge_mid:
                edge_mid[key] = len(nodes_list)
                nodes_list.append((0.5 * (base.nodes[a] + base.nodes[b])).tolist())
            return edge_mid[key]
        tets10 = []
        for a, b, c, d in tets4:
            tets10.append([a, b, c, d, _mid(a, b), _mid(b, c), _mid(c, a), _mid(a, d), _mid(b, d), _mid(c, d)])
        return MeshModel(nodes=np.asarray(nodes_list, dtype=np.float64), blocks=[ElemBlock("tet10", np.asarray(tets10, dtype=np.int64), part_id)])

    raise ValueError(f"Unsupported etype {etype}")


@pytest.mark.kratos
def test_p3_1_thb_hierarchical_solve_and_cantilever(kratos) -> None:
    """T-P3.1 Path A: Verify THB algebraic reduction and Kratos MPC solve match & converge."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]

    # V03 cantilever geometry with 2-level local THB refinement near root z in [0, 0.4]
    thb = THBHierarchy.from_uniform(
        (-0.05, 0.05), (-0.05, 0.05), (0.0, 1.0), (1, 1, 4), degree=2
    )
    thb.mark_box((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4), target_level=1)
    thb.enforce_2_to_1_balance()

    E_mod, nu = 2.0e11, 0.0
    part = PartSpec(name="beam", kind="thb", thb=thb, young_modulus=E_mod, poisson_ratio=nu)
    spec = AssemblySpec(
        job_name="thb_cantilever",
        parts={"beam": part},
        instances=[InstanceSpec("beam", "i1")],
    )
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(spec)
    rt = built.get_instance("beam_i1")
    box_fine = rt.box
    assert box_fine is not None and rt.T_basis is not None

    # 1. Exact affine patch test under THB reduction: u(x,y,z) = [0.002*z, -0.001*x, 0.003*y]
    # Note: pure shear/rigid displacement -> zero internal stress imbalance! Actually apply exact linear uniaxial strain:
    # u = [0, 0, eps * z] -> nu=0 => sigma_zz = E*eps => traction [0, 0, E*eps] on z=1, u_z=0 on z=0
    eps = 1e-4
    pz = int(box_fine.degree)
    nu_cp, nv_cp, nw_cp = box_fine.n_cp
    cp_coords = box_fine.control_points()
    for idx, nid in enumerate(rt.node_ids):
        nd = built.structure.GetNode(nid)
        x, y, z = cp_coords[idx]
        # Fix z=0 layer (exact zero u)
        k_layer = idx // (nu_cp * nv_cp)
        if k_layer == 0:
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)

    f_axial = _integrate_iga_face_traction(box_fine, axis=2, side="hi", traction=(0.0, 0.0, E_mod * eps))
    _apply_iga_loads(KM, SM, built.structure, rt.node_ids, f_axial, rt.properties_map[0], start_cond_id=1)

    res_red = builder.assemble_and_solve_thb_reduced(built, "beam_i1")
    u_hier = res_red["u_hier"].reshape(-1, 3)
    u_z_exact = eps * cp_coords[:, 2]
    assert np.max(np.abs(u_hier[:, 2] - u_z_exact)) < 1e-11
    assert np.max(np.abs(u_hier[:, :2])) < 1e-12

    # 2. Cantilever tip shear V03: P = 5.0 N in -y over A = 0.01 m^2 => t_y = -500.0 Pa
    v03 = json.loads((_DATA_DIR / "v03.json").read_text(encoding="utf-8"))
    eb_tip = float(v03["reference_values"]["euler_bernoulli_tip_uy_m"]["traction"])  # -1.0e-6 m

    built2 = builder.build(spec)
    rt2 = built2.get_instance("beam_i1")
    for idx, nid in enumerate(rt2.node_ids):
        k_layer = idx // (nu_cp * nv_cp)
        if k_layer == 0:
            nd = built2.structure.GetNode(nid)
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
    f_shear = _integrate_iga_face_traction(box_fine, axis=2, side="hi", traction=(0.0, -500.0, 0.0))
    _apply_iga_loads(KM, SM, built2.structure, rt2.node_ids, f_shear, rt2.properties_map[0], start_cond_id=1)
    res_bend = builder.assemble_and_solve_thb_reduced(built2, "beam_i1")
    u_red_bend = res_bend["u_hier"].reshape(-1, 3)

    # Now build a 3rd identical assembly solved NATIVELY inside Kratos via LinearMasterSlaveConstraints
    built3 = builder.build(spec)
    rt3 = built3.get_instance("beam_i1")
    builder.attach_thb_master_constraints(built3, "beam_i1")
    # Fix z=0 THB masters (those with non-zero support at z=0) & apply load on fine z=1 CPs?
    # Note in Kratos BlockBuilder, loads on slave nodes are transposed onto masters automatically,
    # whereas Fix() must be on THB master DOFs at z=0!
    T_csr = rt3.T_basis.tocsr()
    z0_rows = [idx for idx in range(len(rt3.node_ids)) if idx // (nu_cp * nv_cp) == 0]
    z0_masters = set(int(c) for c in T_csr[z0_rows, :].indices)
    for midx in z0_masters:
        nd = built3.structure.GetNode(rt3.thb_master_node_ids[midx])
        nd.Fix(KM.DISPLACEMENT_X)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
    _apply_iga_loads(KM, SM, built3.structure, rt3.node_ids, f_shear, rt3.properties_map[0], start_cond_id=1)
    solve_linear(KM, built3.structure, use_pardiso=False, echo=0)
    u_mpc = np.array(
        [list(built3.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
         for nid in rt3.node_ids],
        dtype=np.float64,
    )

    # Algebraic reduction vs Kratos MasterSlaveConstraint solve must agree closely!
    assert np.max(np.abs(u_mpc - u_red_bend)) / abs(eb_tip) < 1e-8
    # Evaluate centerline tip displacement at (0, 0, 1.0)
    N_tip = SurfaceMesh.from_iga_face("tip", box_fine, axis=2, side="hi", dof_node_ids=rt2.node_ids).eval_basis_fn(
        np.array([[0.0, 0.0, 1.0]])
    )
    tip_uy = float((N_tip @ u_red_bend[:, 1])[0])
    rel_err = abs(tip_uy - eb_tip) / abs(eb_tip)
    assert rel_err < 0.025


@pytest.mark.kratos
def test_p3_1_multicell_v04_bimetallic_and_v05_partitioned(kratos) -> None:
    """T-P3.1 Path B: V04 Cu/W bimetallic strip (incl. straddling cells) and V05 partitioned beam."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]

    v04 = json.loads((_DATA_DIR / "v04.json").read_text(encoding="utf-8"))
    v04_ref_merge = float(v04["reference_values"]["tie_comparison_traction_case"]["p2"]["merge_tip_uy_mm"])
    v04_ref_eb = float(v04["reference_values"]["transformed_section_tip_uy_mm"]["traction"])

    # V04 with 3 cells in y so the middle cell STRADDLES the y=0 Cu/W interface!
    box_v04 = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 3, 10), degree=2)
    dom_cu = CellDomainSpec(
        cell_id=1,
        name="Cu",
        young_modulus=1.0e5,
        poisson_ratio=0.0,
        bounds=((-5.0, 5.0), (-5.0, 0.0), (0.0, 100.0)),
    )
    dom_w = CellDomainSpec(
        cell_id=2,
        name="W",
        young_modulus=4.0e5,
        poisson_ratio=0.0,
        bounds=((-5.0, 5.0), (0.0, 5.0), (0.0, 100.0)),
    )
    part_v04 = PartSpec(
        name="bimetallic",
        kind="iga",
        box=box_v04,
        cell_domains=[dom_cu, dom_w],
    )
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built_v04 = builder.build(
        AssemblySpec("v04_shared", {"bimetallic": part_v04}, [InstanceSpec("bimetallic", "m1")])
    )
    rt = built_v04.get_instance("bimetallic_m1")
    # Single shared CP field (no duplicate DOFs!) across both material domains
    assert len(rt.node_ids) == box_v04.n_ctrl
    assert len(rt.properties_map) == 2

    nu_cp, nv_cp, _ = box_v04.n_cp
    for idx, nid in enumerate(rt.node_ids):
        if idx // (nu_cp * nv_cp) == 0:
            nd = built_v04.structure.GetNode(nid)
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
    f_v04 = _integrate_iga_face_traction(box_v04, axis=2, side="hi", traction=(0.0, -0.01, 0.0))
    _apply_iga_loads(KM, SM, built_v04.structure, rt.node_ids, f_v04, rt.properties_map[1], start_cond_id=1)
    solve_linear(KM, built_v04.structure, use_pardiso=False, echo=0)

    u_v04 = np.array(
        [list(built_v04.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
         for nid in rt.node_ids],
        dtype=np.float64,
    )
    N_tip = SurfaceMesh.from_iga_face("v04_tip", box_v04, axis=2, side="hi", dof_node_ids=rt.node_ids).eval_basis_fn(
        np.array([[0.0, 0.0, 100.0]])
    )
    tip_uy_v04 = float((N_tip @ u_v04[:, 1])[0])
    assert abs(tip_uy_v04 - v04_ref_merge) / abs(v04_ref_merge) < 0.01
    assert abs(tip_uy_v04 - v04_ref_eb) / abs(v04_ref_eb) < 0.015

    # V05 multiplanar 7-cell partition vs monolithic single-domain reference
    v05 = json.loads((_DATA_DIR / "v05.json").read_text(encoding="utf-8"))
    box_v05 = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 7), degree=2)
    z_bps = box_v05.breakpoints(2)
    domains_v05 = [
        CellDomainSpec(
            cell_id=c + 1,
            name=f"zone_{c}",
            young_modulus=1.0e11,
            poisson_ratio=0.3,
            bounds=((-5.0, 5.0), (-5.0, 5.0), (float(z_bps[c]), float(z_bps[c + 1]))),
        )
        for c in range(7)
    ]
    part_mono = PartSpec("v05_mono", "iga", box=box_v05, young_modulus=1.0e11, poisson_ratio=0.3)
    part_part = PartSpec("v05_part", "iga", box=box_v05, cell_domains=domains_v05)

    tips: list[float] = []
    for pname, part_obj in [("mono", part_mono), ("part", part_part)]:
        b = builder.build(
            AssemblySpec(f"v05_{pname}", {part_obj.name: part_obj}, [InstanceSpec(part_obj.name, "i1")])
        )
        r = b.get_instance(f"{part_obj.name}_i1")
        n_u, n_v, _ = box_v05.n_cp
        for idx, nid in enumerate(r.node_ids):
            if idx // (n_u * n_v) == 0:
                nd = b.structure.GetNode(nid)
                nd.Fix(KM.DISPLACEMENT_X)
                nd.Fix(KM.DISPLACEMENT_Y)
                nd.Fix(KM.DISPLACEMENT_Z)
        f_v05 = _integrate_iga_face_traction(box_v05, axis=2, side="hi", traction=(0.0, -0.01, 0.0))
        _apply_iga_loads(
            KM, SM, b.structure, r.node_ids, f_v05, next(iter(r.properties_map.values())), start_cond_id=1
        )
        solve_linear(KM, b.structure, use_pardiso=False, echo=0)
        u_arr = np.array(
            [list(b.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
             for nid in r.node_ids],
            dtype=np.float64,
        )
        N_t = SurfaceMesh.from_iga_face("t", box_v05, 2, "hi", r.node_ids).eval_basis_fn(np.array([[0.0, 0.0, 100.0]]))
        tips.append(float((N_t @ u_arr[:, 1])[0]))

    rel_v05_diff = abs(tips[1] - tips[0]) / abs(tips[0])
    assert rel_v05_diff < 1e-12
    assert v05["id"] == "V05"


@pytest.mark.kratos
def test_p3_2_dependent_instances_rotation_translation_invariance(kratos, tmp_path: Path) -> None:
    """T-P3.2: Three rotated & translated dependent instances match within 1e-10."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]

    box = BoxSpline((-0.05, -0.05, 0.0), (0.05, 0.05, 1.0), (2, 2, 6), degree=2)
    part = PartSpec("cant", "iga", box=box, young_modulus=2.0e11, poisson_ratio=0.3)

    a1, a2 = np.deg2rad(42.0), np.deg2rad(-25.0)
    Rz = np.array([[np.cos(a1), -np.sin(a1), 0.0], [np.sin(a1), np.cos(a1), 0.0], [0.0, 0.0, 1.0]])
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(a2), -np.sin(a2)], [0.0, np.sin(a2), np.cos(a2)]])
    R1 = Rz @ Rx
    b1 = np.deg2rad(-60.0)
    R2 = np.array([[np.cos(b1), 0.0, np.sin(b1)], [0.0, 1.0, 0.0], [-np.sin(b1), 0.0, np.cos(b1)]])

    instances = [
        InstanceSpec("cant", "ref", dependent=True, translation=(0.0, 0.0, 0.0)),
        InstanceSpec("cant", "rot1", dependent=True, translation=(10.0, -4.0, 2.5), rotation_matrix=R1),
        InstanceSpec("cant", "rot2", dependent=True, translation=(-6.0, 7.0, -3.0), rotation_matrix=R2),
    ]
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(AssemblySpec("multi_inst_job", {"cant": part}, instances))
    resolved_path = built.save_resolved_json(tmp_path / "multi_inst_job.resolved.json")
    assert resolved_path.exists()

    nu_cp, nv_cp, _ = box.n_cp
    f_local = _integrate_iga_face_traction(box, axis=2, side="hi", traction=(100.0, -500.0, 200.0))
    cond_id = 1
    for rt in built.instances.values():
        for idx, nid in enumerate(rt.node_ids):
            if idx // (nu_cp * nv_cp) == 0:
                nd = built.structure.GetNode(nid)
                nd.Fix(KM.DISPLACEMENT_X)
                nd.Fix(KM.DISPLACEMENT_Y)
                nd.Fix(KM.DISPLACEMENT_Z)
        R_inst = rt.spec.rotation()
        f_rot = {k: R_inst @ v for k, v in f_local.items()}
        cond_id = _apply_iga_loads(
            KM, SM, built.structure, rt.node_ids, f_rot, rt.properties_map[0], start_cond_id=cond_id
        )

    solve_linear(KM, built.structure, use_pardiso=False, echo=0)

    local_disps: dict[str, np.ndarray] = {}
    for key, rt in built.instances.items():
        u_glob = np.array(
            [list(built.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
             for nid in rt.node_ids],
            dtype=np.float64,
        )
        local_disps[key] = rt.spec.inverse_transform_vector(u_glob)

    ref_u = local_disps["cant_ref"]
    ref_max = float(np.max(np.abs(ref_u)))
    assert ref_max > 1e-8
    for key in ("cant_rot1", "cant_rot2"):
        rel_diff = float(np.max(np.abs(local_disps[key] - ref_u))) / ref_max
        assert rel_diff < 1e-10

    # Voigt stress tensor transformation verification for rotated instance
    stress_local = np.array([10.0, -3.0, 25.0, 4.0, -2.0, 7.0])
    stress_glob = transform_stress_voigt(stress_local, R1)
    stress_back = transform_stress_voigt(stress_glob, R1.T)
    assert np.max(np.abs(stress_back - stress_local)) < 1e-13


@pytest.mark.kratos
def test_p3_2_nonmatching_mortar_rbe2_and_rbe3_v03_moment(kratos) -> None:
    """T-P3.2: Free-interface mortar L2 stress patch, RBE2 kinematics, and RBE3 V03-moment."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]
    builder = AssemblyBuilder(KM=KM, SM=SM)

    # 1. Free-interface non-matching FE-FE mortar_l2 tension patch (with Poisson nu=0.3!)
    m_bot = _make_structured_hex_mesh((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2, 2, 1), "hex8")
    m_top = _make_structured_hex_mesh((0.0, 0.0, 1.0), (1.0, 1.0, 2.0), (3, 3, 1), "hex8")
    E_p, nu_p = 1.0e5, 0.3
    spec_p = AssemblySpec(
        "mortar_patch",
        {
            "bot": PartSpec("bot", "native_fe", mesh=m_bot, young_modulus=E_p, poisson_ratio=nu_p),
            "top": PartSpec("top", "native_fe", mesh=m_top, young_modulus=E_p, poisson_ratio=nu_p),
        },
        [InstanceSpec("bot", "b"), InstanceSpec("top", "t")],
    )
    built_p = builder.build(spec_p)
    rt_b = built_p.get_instance("bot_b")
    rt_t = built_p.get_instance("top_t")

    # Minimal 3-2-1 BCs on bottom z=0 leaving lateral Poisson contraction completely free
    for nid in rt_b.node_ids:
        nd = built_p.structure.GetNode(nid)
        if abs(nd.Z0) < 1e-12:
            nd.Fix(KM.DISPLACEMENT_Z)
            if abs(nd.X0) < 1e-12 and abs(nd.Y0) < 1e-12:
                nd.Fix(KM.DISPLACEMENT_X)
                nd.Fix(KM.DISPLACEMENT_Y)
            elif abs(nd.X0 - 1.0) < 1e-12 and abs(nd.Y0) < 1e-12:
                nd.Fix(KM.DISPLACEMENT_Y)

    # Apply uniform normal stress sigma_zz = 100 on z=2 quad faces via SurfaceLoadCondition3D4N
    props_t = rt_t.properties_map[0]
    cid = 1
    for row in m_top.blocks[0].conn:
        if np.all(np.abs(m_top.nodes[row[[4, 5, 6, 7]], 2] - 2.0) < 1e-10):
            knids = [rt_t.fe_handle.mesh_to_kratos_node[int(row[k])] for k in (4, 5, 6, 7)]
            cond = built_p.structure.CreateNewCondition("SurfaceLoadCondition3D4N", cid, knids, props_t)
            cond.SetValue(SM.SURFACE_LOAD, [0.0, 0.0, 100.0])
            cid += 1

    surf_m = SurfaceMesh.from_fe("bot_z1", m_bot, rt_b.fe_handle.mesh_to_kratos_node, lambda c: abs(c[2] - 1.0) < 1e-10)
    surf_s = SurfaceMesh.from_fe("top_z1", m_top, rt_t.fe_handle.mesh_to_kratos_node, lambda c: abs(c[2] - 1.0) < 1e-10)
    tie_res = apply_tie(
        KM, built_p.structure, TieSpec("mortar_z1", master=surf_m, slave=surf_s, method="mortar_l2", gauss_order=3),
        id_table=built_p.id_table, couplings_smp=built_p.couplings_smp,
    )
    assert tie_res.resolved_method == "mortar_l2"
    solve_linear(KM, built_p.structure, use_pardiso=False, echo=0)

    eps_zz = 100.0 / E_p
    eps_lat = -nu_p * eps_zz
    max_patch_err = 0.0
    for nid in (*rt_b.node_ids, *rt_t.node_ids):
        nd = built_p.structure.GetNode(nid)
        u = np.array(nd.GetSolutionStepValue(KM.DISPLACEMENT, 0))
        u_ex = np.array([eps_lat * nd.X0, eps_lat * nd.Y0, eps_zz * nd.Z0])
        max_patch_err = max(max_patch_err, float(np.max(np.abs(u - u_ex))))
    assert max_patch_err < 1e-10

    # 2. RBE2 rigid motion check in Kratos
    m_r2 = _make_structured_hex_mesh((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2, 2, 2), "hex8")
    built_r2 = builder.build(
        AssemblySpec("rbe2_test", {"b": PartSpec("b", "native_fe", mesh=m_r2, young_modulus=1e5, poisson_ratio=0.3)}, [InstanceSpec("b", "1")])
    )
    rt_r2 = built_r2.get_instance("b_1")
    z1_verts = [v for v in range(len(m_r2.nodes)) if abs(m_r2.nodes[v, 2] - 1.0) < 1e-10]
    z1_nids = [rt_r2.fe_handle.mesh_to_kratos_node[v] for v in z1_verts]
    ref_nid = built_r2.id_table.next_id("nodes")
    apply_rbe2(
        KM,
        built_r2.structure,
        RBE2Spec("r2", ref_node_id=ref_nid, ref_coord=(0.5, 0.5, 1.2), slave_node_ids=z1_nids, slave_coords=m_r2.nodes[z1_verts]),
        id_table=built_r2.id_table,
        couplings_smp=built_r2.couplings_smp,
    )
    # Minimal support on non-RBE2 nodes? Or fix all DOFs of ref node AND prescribe rigid rotation on RP while rest is free -> Wait: 3D continuum is free on z<1, so z=1 rigid motion moves z=0 rigidly with zero strain! Wait: we must fix at least z=1 via RP -> since RP has 6 fixed DOFs, the 3D mesh anchored to z=1 via RBE2 has ZERO rigid body modes!
    rp = built_r2.structure.GetNode(ref_nid)
    u0 = np.array([0.05, -0.02, 0.01])
    th0 = np.array([0.002, -0.003, 0.004])
    for dv, val in zip((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z), u0):
        rp.Fix(dv); rp.SetSolutionStepValue(dv, 0, float(val))
    for rv, val in zip((KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z), th0):
        rp.Fix(rv); rp.SetSolutionStepValue(rv, 0, float(val))
    solve_linear(KM, built_r2.structure, use_pardiso=False, echo=0)
    # Entire free block must follow exact rigid motion u0 + th0 x (x - x_ref)
    for nid in rt_r2.node_ids:
        nd = built_r2.structure.GetNode(nid)
        u_sol = np.array(nd.GetSolutionStepValue(KM.DISPLACEMENT, 0))
        u_ex = u0 + np.cross(th0, np.array([nd.X0 - 0.5, nd.Y0 - 0.5, nd.Z0 - 1.2]))
        assert np.max(np.abs(u_sol - u_ex)) < 1e-12

    # 3. RBE3 V03-moment benchmark (quadratic hex20 beam with nu=0.3 to prove zero lateral stiffening!)
    v03 = json.loads((_DATA_DIR / "v03.json").read_text(encoding="utf-8"))
    eb_moment_uy = float(v03["reference_values"]["euler_bernoulli_tip_uy_m"]["moment"])  # -1.5e-6 m
    m_v03 = _make_structured_hex_mesh((-0.05, -0.05, 0.0), (0.05, 0.05, 1.0), (2, 2, 8), "hex20")
    built_r3 = builder.build(
        AssemblySpec("v03_rbe3", {"beam": PartSpec("beam", "native_fe", mesh=m_v03, young_modulus=2.0e11, poisson_ratio=0.3)}, [InstanceSpec("beam", "1")])
    )
    rt_r3 = built_r3.get_instance("beam_1")
    for nid in rt_r3.node_ids:
        nd = built_r3.structure.GetNode(nid)
        if abs(nd.Z0) < 1e-12:
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)

    z1_v03 = [v for v in range(len(m_v03.nodes)) if abs(m_v03.nodes[v, 2] - 1.0) < 1e-10]
    z1_nids_v03 = [rt_r3.fe_handle.mesh_to_kratos_node[v] for v in z1_v03]
    rp3_id = built_r3.id_table.next_id("nodes")
    apply_rbe3(
        KM,
        SM,
        built_r3.structure,
        RBE3Spec("r3_tip", ref_node_id=rp3_id, ref_coord=(0.0, 0.0, 1.0), slave_coords=m_v03.nodes[z1_v03], slave_node_ids=z1_nids_v03),
        moment=(-5.0, 0.0, 0.0),  # Mx=-5.0 bends beam toward -y (right-hand rule about x puts y<0 in compression -> wait check sign!)
        mode="mpc",
        id_table=built_r3.id_table,
        couplings_smp=built_r3.couplings_smp,
    )
    solve_linear(KM, built_r3.structure, use_pardiso=False, echo=0)
    rp3 = built_r3.structure.GetNode(rp3_id)
    uy_r3 = float(rp3.GetSolutionStepValue(KM.DISPLACEMENT_Y, 0))
    # Compare magnitude & check anticlastic lateral expansion at tip corner (x=0.05, y=0.05, z=1)
    assert abs(abs(uy_r3) - abs(eb_moment_uy)) / abs(eb_moment_uy) < 0.02
    corner_ux = max(
        abs(float(built_r3.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT_X, 0)))
        for nid in z1_nids_v03
    )
    assert corner_ux > 1e-9  # Non-stiffening RBE3 preserves Poisson anticlastic deformation!


@pytest.mark.kratos
def test_p3_3_native_fe_all_types_and_hybrid_iga_fe_v04(kratos) -> None:
    """T-P3.3: Patch tests for hex8/tet4/tet10/hex20/wedge6 (w/ auto-flip) & IGA-FE V04 Hybrid."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]
    builder = AssemblyBuilder(KM=KM, SM=SM)

    # 1. All 5 Native FE element types with intentionally inverted elements auto-flipped
    E_fe, nu_fe = 1.0e5, 0.25
    eps = 1.0e-3
    inv_perms = {
        "hex8": [4, 5, 6, 7, 0, 1, 2, 3],
        "hex20": [4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 16, 17, 18, 19],
        "tet4": [0, 2, 1, 3],
        "tet10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
        "wedge6": [3, 4, 5, 0, 1, 2],
    }
    for etype in ("hex8", "tet4", "tet10", "hex20", "wedge6"):
        mesh = _make_structured_hex_mesh((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2, 2, 2), etype)
        # Invert element 0 cleanly to prove auto_flip works in Kratos
        mesh.blocks[0].conn[0] = mesh.blocks[0].conn[0][inv_perms[etype]]
        built = builder.build(
            AssemblySpec(f"fe_{etype}", {"p": PartSpec("p", "native_fe", mesh=mesh, young_modulus=E_fe, poisson_ratio=nu_fe)}, [InstanceSpec("p", "1")])
        )
        rt = built.get_instance("p_1")
        assert rt.fe_handle.flipped_count >= 1
        # Prescribe exact linear patch displacement u = [eps*x, -nu*eps*y, -nu*eps*z] on all outer faces
        for nid in rt.node_ids:
            nd = built.structure.GetNode(nid)
            on_bnd = any(abs(c) < 1e-10 or abs(c - 1.0) < 1e-10 for c in (nd.X0, nd.Y0, nd.Z0))
            if on_bnd:
                u_bc = (eps * nd.X0, -nu_fe * eps * nd.Y0, -nu_fe * eps * nd.Z0)
                for dv, val in zip((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z), u_bc):
                    nd.Fix(dv); nd.SetSolutionStepValue(dv, 0, float(val))
        solve_linear(KM, built.structure, use_pardiso=False, echo=0)
        err = 0.0
        for nid in rt.node_ids:
            nd = built.structure.GetNode(nid)
            u = np.array(nd.GetSolutionStepValue(KM.DISPLACEMENT, 0))
            u_ex = np.array([eps * nd.X0, -nu_fe * eps * nd.Y0, -nu_fe * eps * nd.Z0])
            err = max(err, float(np.max(np.abs(u - u_ex))))
        assert err < 1e-11, f"{etype} linear patch error {err}"

    # 2. IGA-FE Hybrid V04 Cu/W bimetallic cantilever (non-matching FE Cu + IGA W!)
    v04 = json.loads((_DATA_DIR / "v04.json").read_text(encoding="utf-8"))
    v04_ref = float(v04["reference_values"]["tie_comparison_traction_case"]["p2"]["merge_tip_uy_mm"])

    mesh_cu = _make_structured_hex_mesh((-5.0, -5.0, 0.0), (5.0, 0.0, 100.0), (2, 2, 20), "hex20")
    box_w = BoxSpline((-5.0, 0.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 16), degree=2)

    spec_hyb = AssemblySpec(
        "v04_hybrid",
        {
            "cu_fe": PartSpec("cu_fe", "native_fe", mesh=mesh_cu, young_modulus=1.0e5, poisson_ratio=0.0),
            "w_iga": PartSpec("w_iga", "iga", box=box_w, young_modulus=4.0e5, poisson_ratio=0.0),
        },
        [InstanceSpec("cu_fe", "1"), InstanceSpec("w_iga", "1")],
    )
    built_h = builder.build(spec_hyb)
    rt_fe = built_h.get_instance("cu_fe_1")
    rt_iga = built_h.get_instance("w_iga_1")

    # Clamp z=0 on both parts
    for nid in rt_fe.node_ids:
        nd = built_h.structure.GetNode(nid)
        if abs(nd.Z0) < 1e-10:
            nd.Fix(KM.DISPLACEMENT_X); nd.Fix(KM.DISPLACEMENT_Y); nd.Fix(KM.DISPLACEMENT_Z)
    nu_w, nv_w, _ = box_w.n_cp
    for idx, nid in enumerate(rt_iga.node_ids):
        if idx // (nu_w * nv_w) == 0:
            nd = built_h.structure.GetNode(nid)
            nd.Fix(KM.DISPLACEMENT_X); nd.Fix(KM.DISPLACEMENT_Y); nd.Fix(KM.DISPLACEMENT_Z)

    # Apply uniform shear traction ty = -0.01 MPa on z=100 of BOTH Cu (FE) and W (IGA)
    from wht_iga.meshio_ext.model import outer_faces
    # For FE hex20: distribute ty=-0.01 across the z=100 face using corner/quad SurfaceLoadCondition3D8N?
    # Note: outer_faces returns 4 corner indices for hex20? Check SurfaceLoadCondition3D4N on 4 corners or 2 tris:
    # Even simpler & exact for uniform traction on [-5,5]x[-5,0]: integrate bilinear/tri FE face or consistent hex8 split,
    # Actually let's check: MeshModel outer_faces gives quad4 corners; for exact quadratic hex20 face traction,
    # create SurfaceLoadCondition3D8N using [0,1,2,3,8,9,10,11] OR refine to hex20 end-face!
    # Wait: let's compute exact Serendipity quad8 surface nodal forces via 3x3 Gauss quadrature on each quad!
    gx3, gw3 = np.polynomial.legendre.leggauss(3)
    fe_z100_loads: dict[int, np.ndarray] = {}
    hconn_20 = mesh_cu.blocks[0].conn
    # Identify top z=100 quad8 face of each hex20: local nodes [4,5,6,7, 12,13,14,15] in MeshModel hex20!
    for row in hconn_20:
        if abs(mesh_cu.nodes[row[4], 2] - 100.0) < 1e-10:
            q_idx = [int(row[k]) for k in (4, 5, 6, 7, 12, 13, 14, 15)]
            q_pts = mesh_cu.nodes[q_idx]
            for eta, we in zip(gx3, gw3):
                for xi, wx in zip(gx3, gw3):
                    # Standard 8-node serendipity shape functions in (xi, eta)
                    N8 = np.array([
                        0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1),
                        0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1),
                        0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1),
                        0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1),
                        0.5 * (1 - xi * xi) * (1 - eta),
                        0.5 * (1 + xi) * (1 - eta * eta),
                        0.5 * (1 - xi * xi) * (1 + eta),
                        0.5 * (1 - xi) * (1 - eta * eta),
                    ])
                    dx = 0.5 * float(np.linalg.norm(q_pts[1] - q_pts[0]))
                    dy = 0.5 * float(np.linalg.norm(q_pts[3] - q_pts[0]))
                    w_tot = float(wx * we * dx * dy)
                    for a, mi in enumerate(q_idx):
                        knid = rt_fe.fe_handle.mesh_to_kratos_node[mi]
                        fe_z100_loads[knid] = fe_z100_loads.get(knid, np.zeros(3)) + N8[a] * w_tot * np.array([0.0, -0.01, 0.0])

    cid_h = 1
    props_cu = rt_fe.properties_map[0]
    for knid, fv in fe_z100_loads.items():
        if np.linalg.norm(fv) > 1e-18:
            nd = built_h.structure.GetNode(knid)
            nd.SetSolutionStepValue(SM.POINT_LOAD, 0, [float(x) for x in fv])
            built_h.structure.CreateNewCondition("PointLoadCondition3D1N", cid_h, [knid], props_cu)
            cid_h += 1

    f_w_z100 = _integrate_iga_face_traction(box_w, axis=2, side="hi", traction=(0.0, -0.01, 0.0))
    _apply_iga_loads(KM, SM, built_h.structure, rt_iga.node_ids, f_w_z100, rt_iga.properties_map[0], start_cond_id=cid_h)

    # Interface Tie on y=0: non-matching FE (Cu) vs IGA (W)
    surf_fe_y0 = SurfaceMesh.from_fe(
        "cu_y0", mesh_cu, rt_fe.fe_handle.mesh_to_kratos_node, lambda c: abs(c[1]) < 1e-10
    )
    surf_iga_y0 = SurfaceMesh.from_iga_face(
        "w_y0", box_w, axis=1, side="lo", dof_node_ids=rt_iga.node_ids, n_div=(4, 24)
    )
    tie_hyb = TieSpec("hyb_y0", master=surf_iga_y0, slave=surf_fe_y0, method="auto", E_max=4.0e5, penalty_beta=1e3)
    res_h = solve_hybrid_iga_fe(KM, SM, built_h, ties=[tie_hyb])
    assert res_h.tie_results[0].resolved_method == "penalty_spring"

    # Evaluate tip deflection on IGA W side at (0, 0, 100) and interface jump at (0, 0, 50)
    u_w_sol = np.array(
        [list(built_h.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
         for nid in rt_iga.node_ids],
        dtype=np.float64,
    )
    N_tip_w = SurfaceMesh.from_iga_face("w_tip", box_w, 2, "hi", rt_iga.node_ids).eval_basis_fn(
        np.array([[0.0, 0.0, 100.0]])
    )
    tip_uy_hyb = float((N_tip_w @ u_w_sol[:, 1])[0])
    assert abs(tip_uy_hyb - v04_ref) / abs(v04_ref) < 0.02

    # Check interface displacement jump at (0, 0, 100) between FE node at (0,0,100) and IGA W
    fe_tip_matches = [
        nid for nid in rt_fe.node_ids
        if np.linalg.norm([
            built_h.structure.GetNode(nid).X0 - 0.0,
            built_h.structure.GetNode(nid).Y0 - 0.0,
            built_h.structure.GetNode(nid).Z0 - 100.0,
        ]) < 1e-8
    ]
    assert len(fe_tip_matches) == 1
    fe_tip_uy = float(built_h.structure.GetNode(fe_tip_matches[0]).GetSolutionStepValue(KM.DISPLACEMENT_Y, 0))
    assert abs(fe_tip_uy - tip_uy_hyb) / abs(tip_uy_hyb) < 1e-3
