"""Phase 3 Track P3-B verification tests:

1. Lossless mesh read/write roundtrips across .inp, .bdf, .fem, .rad, .k, .msh
   (nodes, blocks, parts, node_sets, elem_sets, face_sets, SelectionSet export).
2. Body-fit structured all-hex -> NURBS volume converter (single & multi-block
   Cartesian and curvilinear patches, p=1..3, exact reproduction, plus rejection
   of T-junction, irregular valence, pinch-out, incomplete L-block, and mixed
   meshes with UNSTRUCTURED_HEX_UNSUPPORTED).
3. Unified OpenCASCADE STEP/IGES/BREP loader & B-Rep healing pipeline
   (surface/mode-0 IGES -> BRepBuilderAPI_Sewing + MakeSolid + ShapeFix_Solid
   with volume & area conservation).
4. Gmsh Surface Loop / Volume safeguard on shell-only IGES CAD files.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from wht_iga.geometry.cad import brep_area, brep_volume, tessellate
from wht_iga.geometry.occ_io import CADGeometry, HealingReport, load_cad_file
from wht_iga.meshgen.gmsh_runner import GmshOptions, build_geo_script, generate_mesh_with_gmsh
from wht_iga.meshio_ext import (
    BodyFitConversionResult,
    BodyFitSplinePatch,
    BodyFitValidationError,
    ElemBlock,
    MeshModel,
    ValidationIssue,
    mesh_to_bodyfit_patches,
    read_mesh,
    write_mesh,
)
from wht_iga.selection.selector import CadAdapter, MeshAdapter
from wht_iga.selection.sets import SelectionSet, SetRegistry


def _make_rich_mixed_solid_mesh() -> MeshModel:
    """Two-part mixed solid mesh (hex8 + wedge6 + tet4) with arbitrary coords and sets."""
    nodes = np.array(
        [
            [0.0, 0.0, 0.0],               # 0
            [1.12345678901234, 0.0, 0.0],  # 1 (precision check)
            [1.0, 1.05, 0.0],              # 2
            [0.0, 1.0, 0.0],               # 3
            [0.0, 0.0, 1.0],               # 4
            [1.0, 0.0, 1.0],               # 5
            [1.0, 1.0, 1.0],               # 6
            [0.0, 1.0, 1.0],               # 7
            [2.25, 0.0, 0.0],              # 8 (hex8 #2 / wedge6)
            [2.25, 1.0, 0.0],              # 9
            [2.25, 0.0, 1.0],              # 10
            [2.25, 1.0, 1.0],              # 11
            [1.5, 0.5, 2.345],             # 12 (tet4 apex)
        ],
        dtype=np.float64,
    )
    # Two blocks sorted canonically by (elem_type, part_id) so block index ordering matches readers.
    hex_conn = np.array(
        [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 8, 9, 2, 5, 10, 11, 6],
        ],
        dtype=np.int64,
    )
    tet_conn = np.array([[4, 5, 6, 12]], dtype=np.int64)
    wedge_conn = np.array([[8, 9, 2, 10, 11, 6]], dtype=np.int64)

    return MeshModel(
        nodes=nodes,
        blocks=[
            ElemBlock("hex8", hex_conn, part_id=1),
            ElemBlock("tet4", tet_conn, part_id=2),
            ElemBlock("wedge6", wedge_conn, part_id=2),
        ],
        parts={1: "STEEL_CORE", 2: "ALUM_CAP"},
        node_sets={
            "BASE_NODES": [0, 1, 2, 3, 8, 9],
            "TOP_NODES": [4, 5, 6, 7, 10, 11, 12],
            "EMPTY_NSET": [],
        },
        elem_sets={
            "CORE_ELEMS": [0, 1],
            "CAP_ELEMS": [2, 3],
            "ALL_ELEMS": [0, 1, 2, 3],
            "EMPTY_ELSET": [],
        },
        face_sets={
            "BOTTOM_Z": [(0, 0), (1, 0), (3, 0)],
            "SHARED_INTERFACE_LEFT_OWNER": [(0, 3)],  # x+ face of elem 0 (shared with elem 1!)
            "SHARED_INTERFACE_RIGHT_OWNER": [(1, 2)],  # x- face of elem 1 (opposite winding!)
            "TET_SLOPE": [(2, 2)],
        },
        unit_system="mm",
    )


@pytest.mark.parametrize("ext", [".inp", ".bdf", ".fem", ".rad", ".k", ".msh"])
def test_mesh_roundtrip_all_formats_lossless(tmp_path: Path, ext: str) -> None:
    """Write and re-read a multi-part mixed solid mesh across all 6 file extensions."""
    model = _make_rich_mixed_solid_mesh()
    out_path = tmp_path / f"roundtrip_mesh{ext}"

    written = write_mesh(model, out_path)
    assert written == out_path and out_path.is_file()

    back, issues = read_mesh(out_path, unit_system="mm")
    assert not any(i["severity"] == "error" for i in issues), f"Unexpected errors for {ext}: {issues}"
    assert back.unit_system == "mm"
    np.testing.assert_allclose(back.nodes, model.nodes, rtol=0.0, atol=1e-12)

    assert len(back.blocks) == len(model.blocks)
    for b_orig, b_back in zip(model.blocks, back.blocks):
        assert b_back.elem_type == b_orig.elem_type
        assert b_back.part_id == b_orig.part_id
        np.testing.assert_array_equal(b_back.conn, b_orig.conn)

    assert back.parts == model.parts
    assert back.node_sets == model.node_sets
    assert back.elem_sets == model.elem_sets
    assert back.face_sets == model.face_sets


@pytest.mark.parametrize("ext", [".inp", ".bdf", ".rad", ".k", ".msh"])
def test_second_order_hex20_tet10_roundtrip(tmp_path: Path, ext: str) -> None:
    """Ensure 20-node hex and 10-node tet connectivities + sets roundtrip cleanly."""
    # Reference unit hex (8 corners + 12 edge midpoints) + 10-node tet attached on top.
    c8 = np.array(
        [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ],
        dtype=np.float64,
    )
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    mids = np.array([0.5 * (c8[a] + c8[b]) for a, b in edges], dtype=np.float64)
    # Extra tet10 vertices sharing face (4, 5, 7) with apex at (0, 0, 2)
    apex = np.array([[0.0, 0.0, 2.0]], dtype=np.float64)
    tet_mids = np.array(
        [
            0.5 * (c8[4] + c8[5]),
            0.5 * (c8[5] + c8[7]),
            0.5 * (c8[7] + c8[4]),
            0.5 * (c8[4] + apex[0]),
            0.5 * (c8[5] + apex[0]),
            0.5 * (c8[7] + apex[0]),
        ],
        dtype=np.float64,
    )
    nodes = np.vstack([c8, mids, apex, tet_mids])
    model = MeshModel(
        nodes=nodes,
        blocks=[
            ElemBlock("hex20", np.arange(20, dtype=np.int64)[None, :], part_id=1),
            ElemBlock(
                "tet10",
                np.array([[4, 5, 7, 20, 21, 22, 23, 24, 25, 26]], dtype=np.int64),
                part_id=1,
            ),
        ],
        parts={1: "SECOND_ORDER"},
        node_sets={"CORNERS": list(range(8))},
        elem_sets={"QUADRATIC": [0, 1]},
        face_sets={"BOT": [(0, 0)]},
        unit_system="m",
    )

    p = tmp_path / f"quad{ext}"
    write_mesh(model, p)
    back, issues = read_mesh(p, unit_system="m")
    assert not any(i["severity"] == "error" for i in issues)
    np.testing.assert_allclose(back.nodes, model.nodes, atol=1e-12)
    assert [b.elem_type for b in back.blocks] == ["hex20", "tet10"]
    np.testing.assert_array_equal(back.blocks[0].conn, model.blocks[0].conn)
    np.testing.assert_array_equal(back.blocks[1].conn, model.blocks[1].conn)
    assert back.node_sets == model.node_sets
    assert back.elem_sets == model.elem_sets
    assert back.face_sets == model.face_sets


def test_write_mesh_merges_selection_registry(tmp_path: Path) -> None:
    """Verify SelectionSet / SetRegistry node & element sets export cleanly."""
    model = _make_rich_mixed_solid_mesh()
    reg = SetRegistry()
    reg.define_set("SEL_CORNER_NODES", SelectionSet.from_ids("node", [0, 4, 8]))
    reg.define_set("SEL_ENDS", SelectionSet.from_ids("element", [0, 3]))

    p = tmp_path / "with_registry.inp"
    write_mesh(model, p, extra_sets=reg)
    back, _ = read_mesh(p, "mm")
    assert back.node_sets["SEL_CORNER_NODES"] == [0, 4, 8]
    assert back.elem_sets["SEL_ENDS"] == [0, 3]


def test_mesh_dispatch_rejects_unknown_extension(tmp_path: Path) -> None:
    model = _make_rich_mixed_solid_mesh()
    with pytest.raises(ValueError, match="Unsupported mesh format"):
        write_mesh(model, tmp_path / "bad.xyz")
    with pytest.raises(ValueError, match="Unsupported mesh format"):
        read_mesh(tmp_path / "bad.xyz", unit_system="mm")


def _build_structured_hex_grid(
    nx: int,
    ny: int,
    nz: int,
    part_id: int = 1,
    part_name: str = "BLOCK",
    curvilinear: bool = False,
    offset_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> MeshModel:
    """Construct a structured Cartesian or polar-curved hex8 MeshModel."""
    rs = np.linspace(1.0, 2.5, nx + 1)
    ts = np.linspace(0.0, np.pi / 2.0, ny + 1)
    xs = np.linspace(0.0, float(nx), nx + 1) + offset_xyz[0]
    ys = np.linspace(0.0, float(ny), ny + 1) + offset_xyz[1]
    zs = np.linspace(0.0, float(nz), nz + 1) + offset_xyz[2]

    nodes = np.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3), dtype=np.float64)

    def nid(i: int, j: int, k: int) -> int:
        return i + (nx + 1) * (j + (ny + 1) * k)

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                if curvilinear:
                    r, th = rs[i], ts[j]
                    nodes[nid(i, j, k)] = [r * np.cos(th), r * np.sin(th), zs[k]]
                else:
                    nodes[nid(i, j, k)] = [xs[i], ys[j], zs[k]]

    conn = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn.append(
                    [
                        nid(i, j, k),
                        nid(i + 1, j, k),
                        nid(i + 1, j + 1, k),
                        nid(i, j + 1, k),
                        nid(i, j, k + 1),
                        nid(i + 1, j, k + 1),
                        nid(i + 1, j + 1, k + 1),
                        nid(i, j + 1, k + 1),
                    ]
                )
    return MeshModel(
        nodes=nodes,
        blocks=[ElemBlock("hex8", np.asarray(conn, dtype=np.int64), part_id=part_id)],
        parts={part_id: part_name},
        unit_system="mm",
    )


@pytest.mark.parametrize("curvilinear", [False, True])
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_bodyfit_single_block_structured_exact(curvilinear: bool, degree: int) -> None:
    """Convert structured Cartesian and curved hex8 blocks into NURBS patches (p=1..3)."""
    model = _build_structured_hex_grid(3, 4, 2, curvilinear=curvilinear)
    # Shuffle element order to verify topological flood-fill is independent of element ordering.
    rng_perm = np.array([23, 5, 0, 12, 7, 19, 1, 14, 8, 20, 2, 16, 10, 21, 3, 18, 9, 22, 4, 15, 11, 17, 6, 13])
    shuffled = MeshModel(
        nodes=model.nodes,
        blocks=[ElemBlock("hex8", model.blocks[0].conn[rng_perm], part_id=1)],
        parts=model.parts,
        unit_system="mm",
    )

    res = mesh_to_bodyfit_patches(shuffled, degree=degree, continuity="c0")
    assert isinstance(res, BodyFitConversionResult)
    assert res.ok and res.issues == []
    assert len(res.patches) == 1

    patch = res.patches[0]
    assert isinstance(patch, BodyFitSplinePatch)
    assert patch.n_cells == (3, 4, 2)
    assert patch.degree == (degree, degree, degree)

    # Clamped knot lengths + Kratos reduced knot conventions match knots.to_kratos.
    for d in range(3):
        ku = patch.knots[d]
        kk = patch.kratos_knots()[d]
        assert len(ku) == patch.shape[d] + degree + 1
        assert len(kk) == patch.shape[d] + degree - 1

    # Check Fortran flat ordering I = a + nu * (b + nv * c).
    flat = patch.control_points_flat()
    nu, nv, nw = patch.shape
    np.testing.assert_allclose(flat[1 + nu * (2 + nv * 1)], patch.control_points[1, 2, 1])

    # Evaluate spline at all mesh vertices (u=i/nx, v=j/ny, w=k/nz).
    nx, ny, nz = patch.n_cells
    us, vs, ws = np.meshgrid(
        np.linspace(0.0, 1.0, nx + 1),
        np.linspace(0.0, 1.0, ny + 1),
        np.linspace(0.0, 1.0, nz + 1),
        indexing="ij",
    )
    eval_nodes = patch.evaluate(us.ravel(), vs.ravel(), ws.ravel())
    ref_nodes = shuffled.nodes[patch.node_grid.ravel()]
    np.testing.assert_allclose(eval_nodes, ref_nodes, rtol=0.0, atol=1e-12)

    # Check arbitrary interior points against exact p=1 trilinear interpolation inside a cell.
    res_p1 = mesh_to_bodyfit_patches(shuffled, degree=1)
    sample_u = np.array([0.17, 0.42, 0.83])
    sample_v = np.array([0.21, 0.58, 0.91])
    sample_w = np.array([0.13, 0.50, 0.77])
    np.testing.assert_allclose(
        patch.evaluate(sample_u, sample_v, sample_w),
        res_p1.patches[0].evaluate(sample_u, sample_v, sample_w),
        rtol=0.0,
        atol=1e-12,
    )


def test_bodyfit_multiblock_conforming_patches() -> None:
    """Two abutting structured parts produce two NURBS patches with shared interface."""
    m1 = _build_structured_hex_grid(2, 2, 2, part_id=1, part_name="LEFT", offset_xyz=(0.0, 0.0, 0.0))
    m2 = _build_structured_hex_grid(2, 2, 2, part_id=2, part_name="RIGHT", offset_xyz=(2.0, 0.0, 0.0))

    # Merge coincident interface vertices along x = 2.0 into shared global node IDs.
    nodes = list(m1.nodes)
    tol = 1e-12
    m2_map: dict[int, int] = {}
    for idx2, pt in enumerate(m2.nodes):
        found = None
        for idx1, ref in enumerate(nodes):
            if np.linalg.norm(pt - ref) <= tol:
                found = idx1
                break
        if found is None:
            found = len(nodes)
            nodes.append(pt)
        m2_map[idx2] = found

    conn2 = np.vectorize(m2_map.__getitem__)(m2.blocks[0].conn)
    combined = MeshModel(
        nodes=np.asarray(nodes, dtype=np.float64),
        blocks=[
            ElemBlock("hex8", m1.blocks[0].conn, part_id=1),
            ElemBlock("hex8", conn2, part_id=2),
        ],
        parts={1: "LEFT", 2: "RIGHT"},
        unit_system="mm",
    )

    res = mesh_to_bodyfit_patches(combined, degree=2, continuity="smooth")
    assert res.ok and len(res.patches) == 2
    assert len(res.interfaces) == 1
    assert len(res.interfaces[0].shared_nodes) == 9  # (2+1)*(2+1) shared nodes on x=2 plane

    # Both patches agree on their shared boundary plane u_left=1.0 vs u_right=0.0.
    v_s = np.linspace(0.0, 1.0, 5)
    w_s = np.linspace(0.0, 1.0, 5)
    vg, wg = np.meshgrid(v_s, w_s, indexing="ij")
    pts_left = res.patches[0].evaluate(np.ones(vg.size), vg.ravel(), wg.ravel())
    pts_right = res.patches[1].evaluate(np.zeros(vg.size), vg.ravel(), wg.ravel())
    # Verify exact set of physical coordinates spans x=2, y in [0,2], z in [0,2].
    assert np.allclose(pts_left[:, 0], 2.0, atol=1e-12)
    assert np.allclose(pts_right[:, 0], 2.0, atol=1e-12)
    np.testing.assert_allclose(pts_left, pts_right, rtol=0.0, atol=1e-12)


def test_bodyfit_rejects_unstructured_variants() -> None:
    """Detect T-junctions, irregular valence, pinch-outs, L-shape, and mixed types."""
    # 1) Mixed hex8 + tet4
    mixed = _make_rich_mixed_solid_mesh()
    with pytest.raises(BodyFitValidationError) as ei:
        mesh_to_bodyfit_patches(mixed, strict=True)
    assert ei.value.issues[0].code == "UNSTRUCTURED_HEX_UNSUPPORTED"
    assert "immersed IGA" in ei.value.issues[0].message

    # 2) T-junction: 1 coarse hex [0,1]x[0,2]x[0,1] next to 2 stacked fine hexes on [1,2]
    t_nodes = np.array(
        [
            # Coarse hex [0,1] x [0,2] x [0,1] -> 0..7
            [0, 0, 0], [1, 0, 0], [1, 2, 0], [0, 2, 0],
            [0, 0, 1], [1, 0, 1], [1, 2, 1], [0, 2, 1],
            # Mid-plane nodes at x=1, y=1 -> 8, 9 (hanging on edges (1,2) and (5,6))
            [1, 1, 0], [1, 1, 1],
            # Fine nodes at x=2 -> 10..15
            [2, 0, 0], [2, 1, 0], [2, 2, 0],
            [2, 0, 1], [2, 1, 1], [2, 2, 1],
        ],
        dtype=np.float64,
    )
    t_mesh = MeshModel(
        nodes=t_nodes,
        blocks=[
            ElemBlock("hex8", np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int64), part_id=1),
            ElemBlock(
                "hex8",
                np.array(
                    [
                        [1, 10, 11, 8, 5, 13, 14, 9],
                        [8, 11, 12, 2, 9, 14, 15, 6],
                    ],
                    dtype=np.int64,
                ),
                part_id=2,
            ),
        ],
        parts={1: "COARSE", 2: "FINE"},
        unit_system="mm",
    )
    r_tj = mesh_to_bodyfit_patches(t_mesh, strict=False)
    assert not r_tj.ok and r_tj.issues[0].code == "UNSTRUCTURED_HEX_UNSUPPORTED"

    # 3) Irregular interior vertex / edge valence = 5 (pentagonal core around central axis)
    angles = np.linspace(0.0, 2.0 * np.pi, 5, endpoint=False)
    ring_inner = np.column_stack([np.cos(angles), np.sin(angles)])
    ring_outer = 2.0 * ring_inner
    # 5 sectors around center: center (0,0), mid-ray ring_inner[i], mid-outer (ring_inner[i]+ring_inner[i+1]), etc.
    pts2d = [(0.0, 0.0)]
    for i in range(5):
        pts2d.append((float(ring_inner[i, 0]), float(ring_inner[i, 1])))
    for i in range(5):
        mid = ring_inner[i] + ring_inner[(i + 1) % 5]
        pts2d.append((float(mid[0]), float(mid[1])))
    pts2d_arr = np.array(pts2d, dtype=np.float64)
    # Extrude 2 layers in z (z=0, z=1, z=2) so center nodes at z=1 are TRUE interior nodes with valence 10 != 8!
    n_2d = len(pts2d_arr)
    pts3d = np.vstack(
        [
            np.column_stack([pts2d_arr, np.zeros(n_2d)]),
            np.column_stack([pts2d_arr, np.ones(n_2d)]),
            np.column_stack([pts2d_arr, np.full(n_2d, 2.0)]),
        ]
    )
    pent_conn = []
    for layer in range(2):
        off0 = layer * n_2d
        off1 = (layer + 1) * n_2d
        for i in range(5):
            a = 1 + i
            b = 6 + i
            c = 1 + ((i + 1) % 5)
            pent_conn.append(
                [
                    off0 + 0, off0 + a, off0 + b, off0 + c,
                    off1 + 0, off1 + a, off1 + b, off1 + c,
                ]
            )
    pent_mesh = MeshModel(
        nodes=pts3d,
        blocks=[ElemBlock("hex8", np.array(pent_conn, dtype=np.int64), part_id=1)],
        unit_system="mm",
    )
    r_pent = mesh_to_bodyfit_patches(pent_mesh, strict=False)
    assert not r_pent.ok and r_pent.issues[0].code == "UNSTRUCTURED_HEX_UNSUPPORTED"

    # 4) Pinch-out (collapsed edge where two vertices coincide in connectivity)
    pinch_nodes = _build_structured_hex_grid(1, 1, 1).nodes
    pinch_mesh = MeshModel(
        nodes=pinch_nodes,
        blocks=[ElemBlock("hex8", np.array([[0, 1, 2, 3, 4, 5, 6, 6]], dtype=np.int64), part_id=1)],
        unit_system="mm",
    )
    r_pinch = mesh_to_bodyfit_patches(pinch_mesh, strict=False)
    assert not r_pinch.ok and r_pinch.issues[0].code == "UNSTRUCTURED_HEX_UNSUPPORTED"

    # 5) Incomplete L-shaped 3-cell single block (cannot fill 2x2x1 lattice)
    m_2x2 = _build_structured_hex_grid(2, 2, 1)
    l_mesh = MeshModel(
        nodes=m_2x2.nodes,
        blocks=[ElemBlock("hex8", m_2x2.blocks[0].conn[:3], part_id=1)],
        unit_system="mm",
    )
    r_l = mesh_to_bodyfit_patches(l_mesh, strict=False)
    assert not r_l.ok and r_l.issues[0].code == "UNSTRUCTURED_HEX_UNSUPPORTED"


def test_occ_io_iges_mode0_sewing_and_make_solid(tmp_path: Path) -> None:
    """Export real OCP curved filleted box to IGES (mode 0 faces), heal to watertight solid."""
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.IGESControl import IGESControl_Controller, IGESControl_Writer
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape()
    fillet = BRepFilletAPI_MakeFillet(box)
    exp = TopExp_Explorer(box, TopAbs_EDGE)
    fillet.Add(2.0, TopoDS.Edge(exp.Current()))
    fillet.Build()
    fbox = fillet.Shape()

    ref_vol = brep_volume(fbox)
    ref_area = brep_area(fbox)
    assert ref_vol > 5900.0

    # Write IGES in mode 0 (trimmed surfaces/faces only, no TopAbs_SOLID).
    iges_path = tmp_path / "filleted_box.igs"
    IGESControl_Controller.Init_s()
    iw = IGESControl_Writer("MM", 0)
    iw.AddShape(fbox)
    iw.ComputeModel()
    assert iw.Write(str(iges_path)) == 1

    raw_geom = load_cad_file(iges_path, heal=False, unit_system="mm")
    assert isinstance(raw_geom, CADGeometry)
    assert isinstance(raw_geom.healing_report, HealingReport)
    assert raw_geom.file_format == "iges"
    assert raw_geom.healing_report.solids_after == 0
    assert raw_geom.healing_report.free_edges_after > 0
    assert not raw_geom.healing_report.watertight
    assert not raw_geom.has_solid

    healed_geom = load_cad_file(iges_path, heal=True, sewing_tol=1e-6, unit_system="mm")
    rep = healed_geom.healing_report
    assert healed_geom.has_solid and len(healed_geom.solids) == 1
    assert rep.solids_before == 0
    assert rep.solids_after == 1
    assert rep.solids_constructed == 1
    assert rep.free_edges_after == 0
    assert rep.watertight
    assert abs(healed_geom.volume - ref_vol) / ref_vol < 1e-9
    assert abs(healed_geom.area - ref_area) / ref_area < 1e-9

    # Verify compatibility with CadAdapter and geometry.cad.tessellate
    adapter = CadAdapter(healed_geom.shape)
    assert len(adapter._ensure_solids()) == 1
    assert len(adapter._ensure_faces()) == len(healed_geom.faces)
    tess = tessellate(healed_geom.shape, deflection=1.0)
    assert len(tess.tris) > 0

    # Also test STEP loading via unified load_cad_file
    step_path = tmp_path / "filleted_box.step"
    sw = STEPControl_Writer()
    sw.Transfer(fbox, STEPControl_AsIs)
    assert sw.Write(str(step_path)) == 1
    step_geom = load_cad_file(step_path, heal=True, unit_system="mm")
    assert step_geom.file_format == "step"
    assert step_geom.healing_report.watertight
    assert abs(step_geom.volume - ref_vol) / ref_vol < 1e-9


def test_gmsh_surface_loop_safeguard_meshes_iges_mode0_3d(tmp_path: Path) -> None:
    """Verify gmsh_runner Surface Loop safeguard generates 3D tetrahedra on shell-only IGES."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.IGESControl import IGESControl_Controller, IGESControl_Writer

    script = build_geo_script(
        cad_path=tmp_path / "shell_only.igs",
        out_msh_path=tmp_path / "out.msh",
        options=GmshOptions(element_type="tet4", mesh_size_max=10.0),
    )
    assert 'SetFactory("OpenCASCADE");' in script
    assert "Surface Loop(_wht_sl)" in script
    assert "Volume(1) = {_wht_sl};" in script

    gmsh_exe = shutil.which("gmsh")
    assert gmsh_exe is not None, "External gmsh binary expected in test environment"

    box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape()
    iges_path = tmp_path / "shell_only.igs"
    IGESControl_Controller.Init_s()
    iw = IGESControl_Writer("MM", 0)
    iw.AddShape(box)
    iw.ComputeModel()
    assert iw.Write(str(iges_path)) == 1

    out_msh = tmp_path / "iges_tet4.msh"
    res = generate_mesh_with_gmsh(
        cad_path=iges_path,
        out_msh_path=out_msh,
        options=GmshOptions(element_type="tet4", mesh_size_max=8.0, format="msh2"),
        physical_groups={"volume": {"SOLID": [1]}},
        geo_path=tmp_path / "iges_tet4.geo",
    )
    assert res.ok, f"gmsh failed on IGES mode 0: {res.stderr}"
    model, read_issues = read_mesh(out_msh, unit_system="mm")
    assert not any(i["severity"] == "error" for i in read_issues)
    assert model.n_elem > 0
    assert any(b.elem_type == "tet4" for b in model.blocks)
    assert "SOLID" in model.parts.values()
