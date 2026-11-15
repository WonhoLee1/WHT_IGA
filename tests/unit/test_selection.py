"""Unit tests for wht_iga.selection (v2.1 §8.5).

Verifies:
- SelectionSet boolean algebra (|, &, -, ^, ~) and SetRegistry
- OCP CAD geometry (box + cylindrical hole + STEP/IGES roundtrip):
  by_geom_type, within_cylinder, by_normal, flood_fill, measure() vs analytic
- 2x2x2 hex8 mesh: free_surface_faces (24), free_surface_elements (8),
  free_surface_nodes (26, excluding center node), flood_fill, CAD<->FEM bridge
- inspect_contact(): 0.05 mm penetration and 0.1 mm clearance cases
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from wht_iga.selection import (
    CadAdapter,
    MeshAdapter,
    SelectionSet,
    SetRegistry,
    inspect_contact,
    measure_cad,
    measure_mesh,
)


def _make_2x2x2_hex8_mesh(length: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    """Create a structured 2x2x2 hex8 mesh on [0, length]^3 (27 nodes, 8 elements)."""
    xs = np.linspace(0.0, length, 3)
    coords = np.array([(x, y, z) for z in xs for y in xs for x in xs], dtype=np.float64)

    def nid(ix: int, iy: int, iz: int) -> int:
        return iz * 9 + iy * 3 + ix

    elems = []
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                elems.append(
                    [
                        nid(ix, iy, iz),
                        nid(ix + 1, iy, iz),
                        nid(ix + 1, iy + 1, iz),
                        nid(ix, iy + 1, iz),
                        nid(ix, iy, iz + 1),
                        nid(ix + 1, iy, iz + 1),
                        nid(ix + 1, iy + 1, iz + 1),
                        nid(ix, iy + 1, iz + 1),
                    ]
                )
    return coords, np.asarray(elems, dtype=np.int64)


def _make_ocp_box_with_hole(
    lx: float = 40.0,
    ly: float = 20.0,
    lz: float = 10.0,
    radius: float = 5.0,
):
    """Build an OCP box [0,lx]x[0,ly]x[0,lz] with a centered Z-through cylindrical hole."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    box = BRepPrimAPI_MakeBox(lx, ly, lz).Shape()
    ax = gp_Ax2(gp_Pnt(0.5 * lx, 0.5 * ly, -1.0), gp_Dir(0.0, 0.0, 1.0))
    cyl = BRepPrimAPI_MakeCylinder(ax, radius, lz + 2.0).Shape()
    cut = BRepAlgoAPI_Cut(box, cyl)
    cut.Build()
    assert cut.IsDone()
    return cut.Shape()


def _make_surface_box_tris(
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    z0: float,
    z1: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Create a closed outward-oriented triangle mesh of an axis-aligned box."""
    nodes = np.array(
        [
            [x0, y0, z0],  # 0
            [x1, y0, z0],  # 1
            [x1, y1, z0],  # 2
            [x0, y1, z0],  # 3
            [x0, y0, z1],  # 4
            [x1, y0, z1],  # 5
            [x1, y1, z1],  # 6
            [x0, y1, z1],  # 7
        ],
        dtype=np.float64,
    )
    tris = np.array(
        [
            # bottom (z0, normal -z)
            [0, 2, 1],
            [0, 3, 2],
            # top (z1, normal +z)
            [4, 5, 6],
            [4, 6, 7],
            # front (y0, normal -y)
            [0, 1, 5],
            [0, 5, 4],
            # back (y1, normal +y)
            [2, 3, 7],
            [2, 7, 6],
            # left (x0, normal -x)
            [0, 4, 7],
            [0, 7, 3],
            # right (x1, normal +x)
            [1, 2, 6],
            [1, 6, 5],
        ],
        dtype=np.int64,
    )
    return nodes, tris


def test_selection_set_boolean_algebra_and_registry() -> None:
    """Verify SelectionSet |, &, -, ^, ~ operators and SetRegistry derivations."""
    uni = range(6)
    a = SelectionSet.from_ids("node", [0, 1, 2, 2], name="a", universe=uni)
    b = SelectionSet.from_ids("node", [2, 3, 4], name="b", universe=uni)

    assert a.ids == (0, 1, 2)
    assert (a | b).ids == (0, 1, 2, 3, 4)
    assert (a & b).ids == (2,)
    assert (a - b).ids == (0, 1)
    assert (a ^ b).ids == (0, 1, 3, 4)
    assert (~a).ids == (3, 4, 5)
    assert a.union(b) == (a | b)
    assert a.intersection(b) == (a & b)
    assert a.difference(b) == (a - b)
    assert a.symmetric_difference(b) == (a ^ b)
    assert a.complement() == (~a)

    # Kind mismatch must raise ValueError
    e = SelectionSet.from_ids("element", [0, 1])
    with pytest.raises(ValueError, match="different kinds"):
        _ = a | e

    reg = SetRegistry()
    reg.define_set("a", a)
    reg.define_set("b", b)
    reg.derive_union("u", "a", "b")
    reg.derive_intersection("i", "a", "b")
    reg.derive_difference("d", "a", "b")
    reg.derive_symmetric_difference("sd", "a", "b")
    reg.derive_complement("ca", "a")

    assert reg["u"].ids == (0, 1, 2, 3, 4)
    assert reg["i"].ids == (2,)
    assert reg["d"].ids == (0, 1)
    assert reg["sd"].ids == (0, 1, 3, 4)
    assert reg["ca"].ids == (3, 4, 5)
    assert set(reg.by_kind("node").keys()) == {"a", "b", "u", "i", "d", "sd", "ca"}


def test_ocp_cad_selection_roundtrip_and_measure(tmp_path: Path) -> None:
    """Verify OCP box+hole selection, STEP/IGES roundtrip, and analytic measure()."""
    from OCP.IGESControl import IGESControl_Writer
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

    lx, ly, lz, r = 40.0, 20.0, 10.0, 5.0
    shape = _make_ocp_box_with_hole(lx, ly, lz, r)

    # STEP roundtrip
    step_path = tmp_path / "plate_hole.stp"
    sw = STEPControl_Writer()
    sw.Transfer(shape, STEPControl_AsIs)
    assert int(sw.Write(str(step_path))) == 1

    # IGES roundtrip (BRep mode=1 preserves solids and analytical surface types)
    from OCP.IGESControl import IGESControl_Controller

    IGESControl_Controller.Init_s()
    iges_path = tmp_path / "plate_hole.igs"
    iw = IGESControl_Writer("MM", 1)
    iw.AddShape(shape)
    iw.ComputeModel()
    assert iw.Write(str(iges_path))

    for cad_path in (step_path, iges_path):
        cad = CadAdapter.from_file(cad_path)
        cyl_faces = cad.faces().by_geom_type("cylinder", radius=r)
        plane_faces = cad.faces().by_geom_type("plane")
        assert len(cyl_faces) >= 1
        assert len(plane_faces) == 6

        # Within cylinder + normal perpendicular to Z selects the hole wall
        hole_q = cad.faces().within_cylinder((20.0, 10.0, -1.0), (20.0, 10.0, 11.0), r=5.2).by_normal(
            perp=(0.0, 0.0, 1.0), tol_deg=5.0
        )
        assert hole_q.ids == cyl_faces.ids

        # Top planar face
        top_face = cad.faces().by_normal((0.0, 0.0, 1.0), tol_deg=5.0)
        assert len(top_face) == 1

        # Edge geometry types: 2 circle edges (R=5), 12 box line edges (+ cylinder seam line edge)
        circles = cad.edges().by_geom_type("circle", radius=r)
        lines = cad.edges().by_geom_type("line")
        assert len(circles) >= 2
        assert len(lines) >= 12

        # Flood fill with 30 deg on 90-deg box stays on top face; 95 deg expands to side planes
        assert len(top_face.flood_fill(angle_deg=30.0)) == 1
        assert len(top_face.flood_fill(angle_deg=95.0)) >= 5

        # Analytical volume, area, centroid, and inertia tensor
        m_all = cad.solids().measure()
        v_box = lx * ly * lz
        v_cyl = math.pi * (r ** 2) * lz
        v_exact = v_box - v_cyl
        a_exact = 2.0 * (lx * ly + ly * lz + lx * lz) - 2.0 * math.pi * (r ** 2) + 2.0 * math.pi * r * lz
        ixx_exact = (v_box / 12.0) * (ly ** 2 + lz ** 2) - v_cyl * ((r ** 2) / 4.0 + (lz ** 2) / 12.0)
        iyy_exact = (v_box / 12.0) * (lx ** 2 + lz ** 2) - v_cyl * ((r ** 2) / 4.0 + (lz ** 2) / 12.0)
        izz_exact = (v_box / 12.0) * (lx ** 2 + ly ** 2) - v_cyl * ((r ** 2) / 2.0)

        assert m_all.volume() == pytest.approx(v_exact, rel=1e-6)
        assert m_all.area() == pytest.approx(a_exact, rel=1e-6)
        assert np.allclose(m_all.centroid(), [20.0, 10.0, 5.0], atol=1e-6)
        inertia = m_all.inertia_tensor()
        assert np.allclose(np.diag(inertia), [ixx_exact, iyy_exact, izz_exact], rtol=1e-6)
        assert np.allclose(inertia - np.diag(np.diag(inertia)), 0.0, atol=1e-5)

        bbox = m_all.bbox()
        assert np.allclose(bbox.min, [0.0, 0.0, 0.0], atol=1e-6)
        assert np.allclose(bbox.max, [lx, ly, lz], atol=1e-6)
        assert np.allclose(bbox.extents, [lx, ly, lz], atol=1e-6)
        assert bbox.diagonal == pytest.approx(math.sqrt(lx ** 2 + ly ** 2 + lz ** 2), rel=1e-6)


def test_ocp_flood_fill_across_tangent_fillet() -> None:
    """Verify flood_fill(angle_deg=30) propagates across tangent fillet faces."""
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    box = BRepPrimAPI_MakeBox(20.0, 20.0, 20.0).Shape()
    fillet = BRepFilletAPI_MakeFillet(box)
    exp = TopExp_Explorer(box, TopAbs_EDGE)
    # Fillet the first edge (radius 2.0)
    edge0 = TopoDS.Edge(exp.Current())
    fillet.Add(2.0, edge0)
    fillet.Build()
    assert fillet.IsDone()

    cad = CadAdapter(fillet.Shape())
    cyls = cad.faces().by_geom_type("cylinder", radius=2.0)
    assert len(cyls) == 1
    # Starting from the tangent fillet with 30 deg threshold should reach the 2 adjacent planes
    filled = cyls.flood_fill(angle_deg=30.0)
    assert len(filled) == 3


def test_mesh_2x2x2_hex8_free_surface_and_cad_bridge() -> None:
    """Verify 2x2x2 hex8 free_surface_faces (24), nodes (26), measure, and CAD<->FEM bridge."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    coords, elems = _make_2x2x2_hex8_mesh(length=2.0)
    mesh = MeshAdapter({"nodes": coords, "hex8": elems})

    fs_faces = mesh.free_surface_faces()
    fs_elems = mesh.free_surface_elements()
    fs_nodes = mesh.free_surface_nodes()

    assert len(fs_faces) == 24
    assert len(fs_elems) == 8
    assert len(fs_nodes) == 26
    # Center node at (1, 1, 1) is node index 13 and must be excluded from free_surface_nodes
    center_node_idx = 1 * 9 + 1 * 3 + 1
    assert center_node_idx not in fs_nodes

    # Top faces (z = 2.0, normal +z) -> 4 faces; flood_fill(30 deg) stays on those 4 faces
    top_faces = fs_faces.by_normal((0.0, 0.0, 1.0), tol_deg=5.0)
    assert len(top_faces) == 4
    assert len(top_faces.flood_fill(angle_deg=30.0)) == 4
    assert len(top_faces.flood_fill(angle_deg=95.0)) == 24

    # Measure full 2x2x2 hex8 mesh: volume=8, centroid=(1,1,1), I = (1/6)*8*4*I3 = (16/3)*I3
    m_mesh = mesh.elements().measure()
    assert m_mesh.volume() == pytest.approx(8.0, rel=1e-12)
    assert np.allclose(m_mesh.centroid(), [1.0, 1.0, 1.0], atol=1e-12)
    assert np.allclose(m_mesh.inertia_tensor(), (16.0 / 3.0) * np.eye(3), atol=1e-10)

    # CAD <-> FEM bridge on matching [0, 2]^3 box
    cad = CadAdapter(BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape())
    cad.bind_mesh(mesh, tol=1e-6)
    cad_top = cad.faces().by_normal((0.0, 0.0, 1.0), tol_deg=5.0)
    top_nset = cad_top.to_node_set("top_nodes")
    top_sset = cad_top.to_surface_set("top_surf")
    top_eset = cad_top.to_element_set("top_elems")

    assert len(top_nset) == 9
    assert len(top_sset) == 4
    assert len(top_eset) == 4

    # Set algebra between derived sets
    mesh.define_set("top_surf", top_sset)
    outer_minus_top = fs_faces - mesh.sets["top_surf"]
    assert len(outer_minus_top) == 20


def test_inspect_contact_penetration_and_clearance() -> None:
    """Verify inspect_contact() on 0.05 mm penetration and 0.1 mm clearance cases."""
    m_nodes, m_tris = _make_surface_box_tris(0.0, 10.0, 0.0, 10.0, -10.0, 0.0)

    # Case 1: 0.05 mm penetration (slave bottom at z = -0.05, facing master top at z = 0.0)
    s_pen_nodes, s_pen_tris = _make_surface_box_tris(2.0, 8.0, 2.0, 8.0, -0.05, 9.95)
    # Use only the bottom face of slave (z = -0.05, tris 0..1) and top face of master (z = 0.0, tris 2..3)
    rep_pen = inspect_contact((m_nodes, m_tris[2:4]), (s_pen_nodes, s_pen_tris[0:2]))
    assert rep_pen.max_penetration == pytest.approx(0.05, abs=1e-9)
    assert rep_pen.min_clearance == pytest.approx(0.0, abs=1e-9)
    assert rep_pen.normal_opposition_deg == pytest.approx(0.0, abs=1e-9)
    assert rep_pen.n_penetrating_points == rep_pen.n_sample_points
    assert rep_pen.contact_area_ratio == pytest.approx(1.0, abs=1e-9)
    assert rep_pen.suggested_search_radius > 0.05

    # Case 2: 0.10 mm clearance (slave bottom at z = +0.10)
    s_gap_nodes, s_gap_tris = _make_surface_box_tris(2.0, 8.0, 2.0, 8.0, 0.10, 10.10)
    rep_gap = inspect_contact((m_nodes, m_tris[2:4]), (s_gap_nodes, s_gap_tris[0:2]))
    assert rep_gap.max_penetration == pytest.approx(0.0, abs=1e-9)
    assert rep_gap.min_clearance == pytest.approx(0.10, abs=1e-9)
    assert rep_gap.max_gap == pytest.approx(0.10, abs=1e-9)
    assert rep_gap.mean_gap == pytest.approx(0.10, abs=1e-9)
    assert rep_gap.normal_opposition_deg == pytest.approx(0.0, abs=1e-9)
    assert rep_gap.n_penetrating_points == 0
    assert rep_gap.suggested_search_radius > 0.10

    # Also verify full closed solids (with back faces filtered out by opposing normals)
    rep_full_pen = inspect_contact((m_nodes, m_tris), (s_pen_nodes, s_pen_tris[0:2]))
    assert rep_full_pen.max_penetration == pytest.approx(0.05, abs=1e-9)

    # Verify measure_cad distance_to between two separated CAD boxes
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Vec
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Trsf

    b1 = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(0.0, 0.0, 12.5))
    b2 = BRepBuilderAPI_Transform(b1, trsf, True).Shape()
    assert measure_cad(b1).distance_to(measure_cad(b2)) == pytest.approx(2.5, abs=1e-9)
