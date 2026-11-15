"""S1 CAD: OCP read/heal/tessellate/validate (v2 §7.1).

Pipeline: ``STEPControl_Reader`` → ``ShapeFix_Shape`` →
``BRepBuilderAPI_Sewing`` (tol = 1e-6 × bbox diagonal) → per-shape
``BRepMesh_IncrementalMesh`` (linear deflection = tessellation_size,
angular deviation 0.3). Face triangulations come from ``BRep_Tool``;
shared-edge vertices are welded with a tolerance and ``face_id`` kept.
Default tessellation_size = min(h)/4 (caller passes grid spacing).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .issues import IssueLog


@dataclass
class Tessellation:
    """Welded triangle soup with per-triangle face tags."""

    nodes: np.ndarray  # (N, 3) float64
    tris: np.ndarray  # (M, 3) int64 into nodes
    face_ids: np.ndarray  # (M,) int64
    deflection: float = 0.0


@dataclass
class ValidationReport:
    """Geometric error record for one tessellated shape."""

    tess_area: float = 0.0
    brep_area: float = 0.0
    area_err: float = 0.0  # |tess - brep| / brep
    tess_volume: float = 0.0  # signed volume of closed mesh
    brep_volume: float = 0.0
    max_dev: float = 0.0  # sampled BRepExtrema vertex-to-shape distance
    watertight: bool = False
    issues: IssueLog = field(default_factory=IssueLog)


def read_step(path: str):  # type: ignore[no-untyped-def]
    """Read a STEP file into a TopoDS_Shape."""
    from OCP.STEPControl import STEPControl_Reader

    reader = STEPControl_Reader()
    if int(reader.ReadFile(path)) != 1:
        raise IOError(f"STEP read failed: {path}")
    reader.TransferRoots()
    return reader.OneShape()


def write_step(shape, path: str) -> None:  # type: ignore[no-untyped-def]
    """Write a TopoDS_Shape to STEP."""
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs

    writer = STEPControl_Writer()
    writer.Transfer(shape, STEPControl_AsIs)
    if writer.Write(path) != 1:
        raise IOError(f"STEP write failed: {path}")


def shape_bbox(shape) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[no-untyped-def]
    """Axis-aligned bbox of a shape as (lo, hi)."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box, True)
    lo = np.array([box.GetXMin(), box.GetYMin(), box.GetZMin()])
    hi = np.array([box.GetXMax(), box.GetYMax(), box.GetZMax()])
    return lo, hi


def bbox_diag(shape) -> float:  # type: ignore[no-untyped-def]
    """Bbox diagonal length (sewing-tolerance reference)."""
    lo, hi = shape_bbox(shape)
    return float(np.linalg.norm(hi - lo))


def heal(shape, tol: float | None = None):  # type: ignore[no-untyped-def]
    """ShapeFix + Sewing; default tol = 1e-6 × bbox diagonal (v2 §7.1)."""
    from OCP.ShapeFix import ShapeFix_Shape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing

    if tol is None:
        tol = 1e-6 * bbox_diag(shape)
    fixer = ShapeFix_Shape(shape)
    fixer.SetPrecision(tol)
    fixer.SetMinTolerance(1e-9)
    fixer.SetMaxTolerance(0.1 * tol if tol > 0 else 1e-4)
    fixer.Perform()
    sewing = BRepBuilderAPI_Sewing(tol)
    sewing.Add(fixer.Shape())
    sewing.Perform()
    return sewing.SewedShape()


def brep_area(shape) -> float:  # type: ignore[no-untyped-def]
    """Exact BRep surface area via BRepGProp."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props)
    return float(props.Mass())


def brep_volume(shape) -> float:  # type: ignore[no-untyped-def]
    """Exact BRep volume via BRepGProp."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return float(props.Mass())


def tessellate(
    shape,  # type: ignore[no-untyped-def]
    deflection: float,
    angular: float = 0.3,
    weld_tol: float | None = None,
) -> Tessellation:
    """Mesh shape and extract welded per-face triangulations.

    One ``BRepMesh_IncrementalMesh`` call per shape (per cell); face
    triangles are read with ``BRep_Tool.Triangulation_s`` + location
    transform, then welded by coordinate quantization so shared-edge
    vertices coincide exactly. ``face_id`` is the 0-based face order.
    """
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.BRep import BRep_Tool
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopoDS import TopoDS
    from OCP.TopLoc import TopLoc_Location

    BRepMesh_IncrementalMesh(shape, deflection, False, angular, True)
    lo, hi = shape_bbox(shape)
    scale = float(np.linalg.norm(hi - lo)) or 1.0
    if weld_tol is None:
        weld_tol = 1e-7 * scale
    key_scale = 1.0 / weld_tol

    node_keys: dict[tuple[int, int, int], int] = {}
    nodes: list[list[float]] = []
    tris: list[tuple[int, int, int]] = []
    fids: list[int] = []

    def _nid(p: tuple[float, float, float]) -> int:
        key = (round(p[0] * key_scale), round(p[1] * key_scale), round(p[2] * key_scale))
        idx = node_keys.get(key)
        if idx is None:
            idx = len(nodes)
            node_keys[key] = idx
            nodes.append([p[0], p[1], p[2]])
        return idx

    exp = TopExp_Explorer(shape, TopAbs_FACE)
    face_id = 0
    n_faces = 0
    from OCP.TopAbs import TopAbs_Orientation

    reversed_flag = TopAbs_Orientation.TopAbs_REVERSED
    while exp.More():
        face = TopoDS.Face(exp.Current())
        flip = face.Orientation() == reversed_flag
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is not None:
            trsf = loc.Transformation()
            from OCP.gp import gp_Pnt

            pts: list[tuple[float, float, float]] = []
            for i in range(1, tri.NbNodes() + 1):
                nd = tri.Node(i)
                p = gp_Pnt(nd.X(), nd.Y(), nd.Z())
                p.Transform(trsf)  # identity-safe: keeps location handling exact
                pts.append((p.X(), p.Y(), p.Z()))
            for i in range(1, tri.NbTriangles() + 1):
                t = tri.Triangle(i)
                a, b, c = t.Value(1) - 1, t.Value(2) - 1, t.Value(3) - 1
                if flip:  # REVERSED face: swap winding → outward normal
                    b, c = c, b
                na, nb, nc = _nid(pts[a]), _nid(pts[b]), _nid(pts[c])
                if len({na, nb, nc}) < 3:
                    continue  # pole-degenerate triangle: carries no area
                tris.append((na, nb, nc))
                fids.append(face_id)
            face_id += 1
        n_faces += 1
        exp.Next()
    if n_faces == 0:
        raise ValueError("tessellate: shape has no faces")
    return Tessellation(
        nodes=np.asarray(nodes),
        tris=np.asarray(tris, dtype=np.int64),
        face_ids=np.asarray(fids, dtype=np.int64),
        deflection=deflection,
    )


def tess_area(tess: Tessellation) -> float:
    """Sum of triangle areas."""
    p = tess.nodes[tess.tris]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


def signed_volume(tess: Tessellation) -> float:
    """Signed mesh volume Σ v0·(v1×v2)/6 (shape function of orientation)."""
    p = tess.nodes[tess.tris]
    return float((p[:, 0] * np.cross(p[:, 1], p[:, 2])).sum(axis=1).sum() / 6.0)


def boundary_edge_count(tess: Tessellation) -> int:
    """Edges used by exactly one triangle (0 = watertight; skips self-loops)."""
    t = tess.tris
    edges = np.sort(
        np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1
    )
    edges = edges[edges[:, 0] != edges[:, 1]]
    if len(edges) == 0:
        return 0
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return int((counts == 1).sum())


def sampled_max_deviation(shape, tess: Tessellation, n_samples: int = 200) -> float:  # type: ignore[no-untyped-def]
    """Max vertex/edge-midpoint-to-BRep distance over a strided sample.

    Mesh nodes sit exactly on the surface, so vertices alone cannot see
    chordal error — edge midpoints are included for that reason.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.gp import gp_Pnt

    pts = [tess.nodes]
    e = tess.nodes[tess.tris]
    pts.append((e[:, 0] + e[:, 1]) / 2)
    pts.append((e[:, 1] + e[:, 2]) / 2)
    pts.append((e[:, 2] + e[:, 0]) / 2)
    allp = np.vstack(pts)
    idx = np.linspace(0, len(allp) - 1, min(n_samples, len(allp))).astype(int)
    worst = 0.0
    for i in idx:
        x, y, z = allp[i]
        v = BRepBuilderAPI_MakeVertex(gp_Pnt(float(x), float(y), float(z))).Vertex()
        worst = max(worst, float(BRepExtrema_DistShapeShape(v, shape).Value()))
    return worst


def validate_shape(
    shape,  # type: ignore[no-untyped-def]
    tess: Tessellation,
    n_samples: int = 200,
) -> ValidationReport:
    """Fill geometric-error record; error issues block ``validate()``."""
    log = IssueLog()
    ta = tess_area(tess)
    ba = brep_area(shape)
    area_err = abs(ta - ba) / ba if ba > 0 else float("inf")
    sv = signed_volume(tess)
    bv = brep_volume(shape)
    n_bnd = boundary_edge_count(tess)
    watertight = n_bnd == 0
    if not watertight:
        log.add("error", "NON_WATERTIGHT", f"{n_bnd} boundary edges", "tessellation")
    if not sv > 0:
        log.add("error", "BAD_VOLUME_SIGN", f"signed volume {sv:.6g} <= 0", "tessellation")
    diag = bbox_diag(shape)
    vol_tol = max(1e-6, 10 * tess.deflection / diag) if tess.deflection > 0 else 1e-6
    if abs(sv - bv) / abs(bv) > vol_tol if bv != 0 else True:
        log.add(
            "warning", "VOLUME_MISMATCH",
            f"tess {sv:.9g} vs brep {bv:.9g}", "tessellation",
        )
    dev = sampled_max_deviation(shape, tess, n_samples)
    return ValidationReport(
        tess_area=ta, brep_area=ba, area_err=area_err,
        tess_volume=sv, brep_volume=bv, max_dev=dev,
        watertight=watertight, issues=log,
    )


def fuse_cells(shapes: list) -> object:  # type: ignore[no-untyped-def]
    """General fuse of cell solids via ``BOPAlgo_Builder`` (v2 §7.1).

    OCP (cadquery-ocp) exposes shape lists as
    ``OCP.collections.List_TopoDS_Shape``; the base ``BOPAlgo_Builder``
    with arguments only performs the general fuse.
    """
    from OCP.BOPAlgo import BOPAlgo_Builder
    from OCP.collections import List_TopoDS_Shape

    args = List_TopoDS_Shape()
    for s in shapes:
        args.Append(s)
    builder = BOPAlgo_Builder()
    builder.SetArguments(args)
    builder.Perform()
    if builder.HasErrors():
        raise RuntimeError("BOPAlgo_Builder general fuse failed")
    return builder.Shape()


def cut_shape(base, tool):  # type: ignore[no-untyped-def]
    """Boolean cut (base minus tool) via BRepAlgoAPI."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    op = BRepAlgoAPI_Cut(base, tool)
    op.Build()
    if not op.IsDone():
        raise RuntimeError("BRepAlgoAPI_Cut failed")
    return op.Shape()
