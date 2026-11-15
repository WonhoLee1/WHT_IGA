"""Unified CAD and Mesh entity selector, free-surface extractor, and CAD-FEM bridge (v2.1 §8.5).

Provides :class:`CadAdapter` (OpenCASCADE ``TopoDS_Shape`` / STEP / IGES),
:class:`MeshAdapter` (node/element connectivity or ``MeshModel``), and
:class:`EntityQuery` with fluent spatial/topological/geometric filters,
boolean algebra, free-surface extraction, and CAD-to-FEM set mapping.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from .contact_check import _closest_point_on_triangle
from .measure import MeasureResult, measure_cad, measure_mesh
from .sets import SelectionSet, SetKind, SetRegistry

# Canonical local face vertex indices matching wht_iga.meshio_ext.model
_SOLID_FACE_TABLES: dict[str, tuple[tuple[int, ...], ...]] = {
    "hex8": (
        (0, 3, 2, 1),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
    ),
    "hex20": (
        (0, 3, 2, 1),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
    ),
    "tet4": (
        (0, 2, 1),
        (0, 1, 3),
        (1, 2, 3),
        (2, 0, 3),
    ),
    "tet10": (
        (0, 2, 1),
        (0, 1, 3),
        (1, 2, 3),
        (2, 0, 3),
    ),
    "wedge6": (
        (0, 2, 1),
        (3, 4, 5),
        (0, 1, 4, 3),
        (1, 2, 5, 4),
        (2, 0, 3, 5),
    ),
}


@dataclass(frozen=True)
class _EntityRecord:
    """Cached geometric metadata for a single CAD or mesh entity."""

    entity_id: int
    centroid: np.ndarray
    normal: np.ndarray | None = None
    area: float = 0.0
    length: float = 0.0
    volume: float = 0.0
    geom_type: str = ""
    radius: float | None = None
    bbox_lo: np.ndarray | None = None
    bbox_hi: np.ndarray | None = None


def _dedup_ocp_shapes(shape: Any, topabs_type: Any) -> list[Any]:
    """Extract unique OCP sub-shapes of ``topabs_type`` preserving traversal order."""
    from OCP.TopExp import TopExp_Explorer

    buckets: dict[int, list[Any]] = {}
    unique: list[Any] = []
    exp = TopExp_Explorer(shape, topabs_type)
    while exp.More():
        cur = exp.Current()
        h = hash(cur)
        bucket = buckets.setdefault(h, [])
        if not any(cur.IsSame(existing) for existing in bucket):
            bucket.append(cur)
            unique.append(cur)
        exp.Next()
    return unique


def _ocp_shape_index_lookup(shapes: Sequence[Any]) -> tuple[dict[int, list[tuple[int, Any]]], list[Any]]:
    """Build a hash-bucketed lookup ``shape -> index`` for ``IsSame`` matching."""
    buckets: dict[int, list[tuple[int, Any]]] = {}
    for idx, shp in enumerate(shapes):
        h = hash(shp)
        buckets.setdefault(h, []).append((idx, shp))
    return buckets, list(shapes)


def _find_ocp_shape_index(buckets: dict[int, list[tuple[int, Any]]], target: Any) -> int | None:
    h = hash(target)
    for idx, shp in buckets.get(h, ()):
        if target.IsSame(shp):
            return idx
    return None


def _surface_type_and_radius(face: Any) -> tuple[str, float | None]:
    """Inspect an OCP ``TopoDS_Face`` via ``BRepAdaptor_Surface``."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import (
        GeomAbs_BSplineSurface,
        GeomAbs_BezierSurface,
        GeomAbs_Cone,
        GeomAbs_Cylinder,
        GeomAbs_OffsetSurface,
        GeomAbs_Plane,
        GeomAbs_Sphere,
        GeomAbs_SurfaceOfExtrusion,
        GeomAbs_SurfaceOfRevolution,
        GeomAbs_Torus,
    )

    surf = BRepAdaptor_Surface(face)
    stype = surf.GetType()
    if stype == GeomAbs_Plane:
        return "plane", None
    if stype == GeomAbs_Cylinder:
        return "cylinder", float(surf.Cylinder().Radius())
    if stype == GeomAbs_Sphere:
        return "sphere", float(surf.Sphere().Radius())
    if stype == GeomAbs_Cone:
        return "cone", float(surf.Cone().RefRadius())
    if stype == GeomAbs_Torus:
        return "torus", float(surf.Torus().MinorRadius())
    if stype == GeomAbs_BSplineSurface:
        return "bspline", None
    if stype == GeomAbs_BezierSurface:
        return "bezier", None
    if stype == GeomAbs_SurfaceOfRevolution:
        return "revolution", None
    if stype == GeomAbs_SurfaceOfExtrusion:
        return "extrusion", None
    if stype == GeomAbs_OffsetSurface:
        return "offset", None
    return "other", None


def _curve_type_and_radius(edge: Any) -> tuple[str, float | None]:
    """Inspect an OCP ``TopoDS_Edge`` via ``BRepAdaptor_Curve``."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GeomAbs import (
        GeomAbs_BSplineCurve,
        GeomAbs_BezierCurve,
        GeomAbs_Circle,
        GeomAbs_Ellipse,
        GeomAbs_Hyperbola,
        GeomAbs_Line,
        GeomAbs_OffsetCurve,
        GeomAbs_Parabola,
    )

    crv = BRepAdaptor_Curve(edge)
    ctype = crv.GetType()
    if ctype == GeomAbs_Line:
        return "line", None
    if ctype == GeomAbs_Circle:
        return "circle", float(crv.Circle().Radius())
    if ctype == GeomAbs_Ellipse:
        return "ellipse", None
    if ctype == GeomAbs_Hyperbola:
        return "hyperbola", None
    if ctype == GeomAbs_Parabola:
        return "parabola", None
    if ctype == GeomAbs_BezierCurve:
        return "bezier", None
    if ctype == GeomAbs_BSplineCurve:
        return "bspline", None
    if ctype == GeomAbs_OffsetCurve:
        return "offset", None
    return "other", None


def _face_normal_at_uv(face: Any, u: float, v: float) -> np.ndarray | None:
    """Compute outward unit normal on ``face`` at parameter ``(u, v)``."""
    from OCP.BRepGProp import BRepGProp_Face
    from OCP.gp import gp_Pnt, gp_Vec

    gface = BRepGProp_Face(face)
    pnt = gp_Pnt()
    vec = gp_Vec()
    gface.Normal(u, v, pnt, vec)
    arr = np.array([vec.X(), vec.Y(), vec.Z()], dtype=np.float64)
    nrm = float(np.linalg.norm(arr))
    if nrm <= 1e-15:
        return None
    return arr / nrm


def _face_edge_normal(face: Any, edge: Any) -> np.ndarray | None:
    """Evaluate outward unit normal of ``face`` at the midpoint of ``edge``."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve2d

    try:
        c2d = BRepAdaptor_Curve2d(edge, face)
        u_mid = 0.5 * (c2d.FirstParameter() + c2d.LastParameter())
        uv = c2d.Value(u_mid)
        return _face_normal_at_uv(face, uv.X(), uv.Y())
    except Exception:
        return None


class CadAdapter:
    """Adapter wrapping an OpenCASCADE ``TopoDS_Shape`` for unified querying.

    Parameters
    ----------
    shape:
        OpenCASCADE ``TopoDS_Shape`` instance.
    """

    def __init__(self, shape: Any) -> None:
        self.shape = shape
        self.sets = SetRegistry()
        self._bound_mesh: MeshAdapter | None = None
        self._bound_tol: float | None = None

        self._faces: list[Any] | None = None
        self._edges: list[Any] | None = None
        self._vertices: list[Any] | None = None
        self._solids: list[Any] | None = None

        self._face_records: list[_EntityRecord] | None = None
        self._edge_records: list[_EntityRecord] | None = None
        self._vertex_records: list[_EntityRecord] | None = None
        self._solid_records: list[_EntityRecord] | None = None
        self._face_adj: dict[int, list[tuple[int, float]]] | None = None
        self._tess_cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_step(cls, path: str | Path) -> CadAdapter:
        """Load a STEP (``.stp`` / ``.step``) file into a :class:`CadAdapter`."""
        from OCP.STEPControl import STEPControl_Reader

        reader = STEPControl_Reader()
        if int(reader.ReadFile(str(path))) != 1:
            raise IOError(f"STEP read failed: {path}")
        reader.TransferRoots()
        return cls(reader.OneShape())

    @classmethod
    def from_iges(cls, path: str | Path) -> CadAdapter:
        """Load an IGES (``.igs`` / ``.iges``) file into a :class:`CadAdapter`."""
        from OCP.IGESControl import IGESControl_Reader

        reader = IGESControl_Reader()
        if int(reader.ReadFile(str(path))) != 1:
            raise IOError(f"IGES read failed: {path}")
        reader.TransferRoots()
        return cls(reader.OneShape())

    @classmethod
    def from_file(cls, path: str | Path) -> CadAdapter:
        """Load a STEP or IGES file based on file extension."""
        p = Path(path)
        ext = p.suffix.lower()
        if ext in (".stp", ".step"):
            return cls.from_step(p)
        if ext in (".igs", ".iges"):
            return cls.from_iges(p)
        raise ValueError(f"Unsupported CAD extension {ext!r}; expected .stp/.step/.igs/.iges")

    def bind_mesh(self, mesh: MeshAdapter, tol: float | None = None) -> None:
        """Bind a :class:`MeshAdapter` for CAD-to-FEM set bridging."""
        self._bound_mesh = mesh
        self._bound_tol = tol

    def define_set(self, name: str, selection_or_query: SelectionSet | EntityQuery) -> SelectionSet:
        """Register a named :class:`SelectionSet` in this adapter's registry."""
        return self.sets.define_set(name, selection_or_query)

    def _ensure_faces(self) -> list[Any]:
        if self._faces is None:
            from OCP.TopAbs import TopAbs_FACE
            from OCP.TopoDS import TopoDS

            self._faces = [TopoDS.Face(s) for s in _dedup_ocp_shapes(self.shape, TopAbs_FACE)]
        return self._faces

    def _ensure_edges(self) -> list[Any]:
        if self._edges is None:
            from OCP.TopAbs import TopAbs_EDGE
            from OCP.TopoDS import TopoDS

            self._edges = [TopoDS.Edge(s) for s in _dedup_ocp_shapes(self.shape, TopAbs_EDGE)]
        return self._edges

    def _ensure_vertices(self) -> list[Any]:
        if self._vertices is None:
            from OCP.TopAbs import TopAbs_VERTEX
            from OCP.TopoDS import TopoDS

            self._vertices = [TopoDS.Vertex(s) for s in _dedup_ocp_shapes(self.shape, TopAbs_VERTEX)]
        return self._vertices

    def _ensure_solids(self) -> list[Any]:
        if self._solids is None:
            from OCP.TopAbs import TopAbs_SOLID
            from OCP.TopoDS import TopoDS

            self._solids = [TopoDS.Solid(s) for s in _dedup_ocp_shapes(self.shape, TopAbs_SOLID)]
        return self._solids

    def _get_records(self, kind: SetKind) -> list[_EntityRecord]:
        from OCP.BRep import BRep_Tool
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.BRepBndLib import BRepBndLib
        from OCP.BRepGProp import BRepGProp
        from OCP.Bnd import Bnd_Box
        from OCP.GProp import GProp_GProps

        if kind == "cad_face":
            if self._face_records is None:
                recs: list[_EntityRecord] = []
                for idx, face in enumerate(self._ensure_faces()):
                    gtype, rad = _surface_type_and_radius(face)
                    props = GProp_GProps()
                    BRepGProp.SurfaceProperties_s(face, props)
                    area = float(props.Mass())
                    com = props.CentreOfMass()
                    centroid = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
                    surf = BRepAdaptor_Surface(face)
                    u_mid = 0.5 * (surf.FirstUParameter() + surf.LastUParameter())
                    v_mid = 0.5 * (surf.FirstVParameter() + surf.LastVParameter())
                    normal = _face_normal_at_uv(face, u_mid, v_mid)
                    box = Bnd_Box()
                    BRepBndLib.Add_s(face, box, True)
                    lo = np.array([box.GetXMin(), box.GetYMin(), box.GetZMin()], dtype=np.float64)
                    hi = np.array([box.GetXMax(), box.GetYMax(), box.GetZMax()], dtype=np.float64)
                    recs.append(
                        _EntityRecord(
                            entity_id=idx,
                            centroid=centroid,
                            normal=normal,
                            area=area,
                            geom_type=gtype,
                            radius=rad,
                            bbox_lo=lo,
                            bbox_hi=hi,
                        )
                    )
                self._face_records = recs
            return self._face_records

        if kind == "cad_edge":
            if self._edge_records is None:
                recs = []
                for idx, edge in enumerate(self._ensure_edges()):
                    gtype, rad = _curve_type_and_radius(edge)
                    props = GProp_GProps()
                    BRepGProp.LinearProperties_s(edge, props)
                    length = float(props.Mass())
                    com = props.CentreOfMass()
                    centroid = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
                    box = Bnd_Box()
                    BRepBndLib.Add_s(edge, box, True)
                    lo = np.array([box.GetXMin(), box.GetYMin(), box.GetZMin()], dtype=np.float64)
                    hi = np.array([box.GetXMax(), box.GetYMax(), box.GetZMax()], dtype=np.float64)
                    recs.append(
                        _EntityRecord(
                            entity_id=idx,
                            centroid=centroid,
                            length=length,
                            geom_type=gtype,
                            radius=rad,
                            bbox_lo=lo,
                            bbox_hi=hi,
                        )
                    )
                self._edge_records = recs
            return self._edge_records

        if kind == "cad_vertex":
            if self._vertex_records is None:
                recs = []
                for idx, vtx in enumerate(self._ensure_vertices()):
                    p = BRep_Tool.Pnt_s(vtx)
                    pt = np.array([p.X(), p.Y(), p.Z()], dtype=np.float64)
                    recs.append(_EntityRecord(entity_id=idx, centroid=pt, bbox_lo=pt, bbox_hi=pt))
                self._vertex_records = recs
            return self._vertex_records

        if kind == "cad_solid":
            if self._solid_records is None:
                recs = []
                for idx, sld in enumerate(self._ensure_solids()):
                    props = GProp_GProps()
                    BRepGProp.VolumeProperties_s(sld, props)
                    vol = float(props.Mass())
                    com = props.CentreOfMass()
                    centroid = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
                    box = Bnd_Box()
                    BRepBndLib.Add_s(sld, box, True)
                    lo = np.array([box.GetXMin(), box.GetYMin(), box.GetZMin()], dtype=np.float64)
                    hi = np.array([box.GetXMax(), box.GetYMax(), box.GetZMax()], dtype=np.float64)
                    recs.append(
                        _EntityRecord(
                            entity_id=idx,
                            centroid=centroid,
                            volume=vol,
                            bbox_lo=lo,
                            bbox_hi=hi,
                        )
                    )
                self._solid_records = recs
            return self._solid_records

        raise ValueError(f"Unsupported CAD entity kind: {kind!r}")

    def _get_shapes(self, kind: SetKind, ids: Sequence[int]) -> list[Any]:
        table = {
            "cad_face": self._ensure_faces,
            "cad_edge": self._ensure_edges,
            "cad_vertex": self._ensure_vertices,
            "cad_solid": self._ensure_solids,
        }[kind]()
        return [table[i] for i in ids]

    def _get_face_adjacency(self) -> dict[int, list[tuple[int, float]]]:
        """Return ``{face_id: [(neighbor_face_id, dihedral_angle_deg), ...]}``."""
        if self._face_adj is not None:
            return self._face_adj

        from OCP.TopAbs import TopAbs_EDGE
        from OCP.TopExp import TopExp_Explorer
        from OCP.TopoDS import TopoDS

        faces = self._ensure_faces()
        edges = self._ensure_edges()
        edge_buckets, _ = _ocp_shape_index_lookup(edges)
        recs = self._get_records("cad_face")

        edge_to_faces: dict[int, list[int]] = {}
        for fid, face in enumerate(faces):
            exp = TopExp_Explorer(face, TopAbs_EDGE)
            seen_e: set[int] = set()
            while exp.More():
                e = TopoDS.Edge(exp.Current())
                eid = _find_ocp_shape_index(edge_buckets, e)
                if eid is not None and eid not in seen_e:
                    seen_e.add(eid)
                    edge_to_faces.setdefault(eid, []).append(fid)
                exp.Next()

        adj: dict[int, list[tuple[int, float]]] = {i: [] for i in range(len(faces))}
        for eid, fids in edge_to_faces.items():
            if len(fids) < 2:
                continue
            edge = edges[eid]
            for i in range(len(fids)):
                for j in range(i + 1, len(fids)):
                    f1, f2 = fids[i], fids[j]
                    n1 = _face_edge_normal(faces[f1], edge)
                    if n1 is None:
                        n1 = recs[f1].normal
                    n2 = _face_edge_normal(faces[f2], edge)
                    if n2 is None:
                        n2 = recs[f2].normal
                    if n1 is None or n2 is None:
                        ang = 180.0
                    else:
                        dot = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
                        ang = float(np.degrees(np.arccos(dot)))
                    adj[f1].append((f2, ang))
                    adj[f2].append((f1, ang))

        self._face_adj = adj
        return adj

    def tessellate_faces(
        self,
        face_ids: Sequence[int] | None = None,
        deflection: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Tessellate selected CAD faces into ``(nodes, tris)``."""
        from OCP.BRep import BRep_Tool
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        from OCP.TopAbs import TopAbs_Orientation
        from OCP.TopLoc import TopLoc_Location
        from OCP.gp import gp_Pnt

        if deflection is None:
            box = Bnd_Box()
            BRepBndLib.Add_s(self.shape, box, True)
            diag = float(
                np.linalg.norm(
                    [
                        box.GetXMax() - box.GetXMin(),
                        box.GetYMax() - box.GetYMin(),
                        box.GetZMax() - box.GetZMin(),
                    ]
                )
            )
            deflection = max(1e-3 * diag, 1e-4)

        key = round(float(deflection), 12)
        if key not in self._tess_cache:
            BRepMesh_IncrementalMesh(self.shape, float(deflection), False, 0.3, True)
            nodes_list: list[tuple[float, float, float]] = []
            tris_list: list[tuple[int, int, int]] = []
            tri_fids: list[int] = []
            reversed_flag = TopAbs_Orientation.TopAbs_REVERSED

            for fid, face in enumerate(self._ensure_faces()):
                flip = face.Orientation() == reversed_flag
                loc = TopLoc_Location()
                tri = BRep_Tool.Triangulation_s(face, loc)
                if tri is None:
                    continue
                offset = len(nodes_list)
                trsf = loc.Transformation()
                for i in range(1, tri.NbNodes() + 1):
                    nd = tri.Node(i)
                    gp = gp_Pnt(nd.X(), nd.Y(), nd.Z())
                    gp.Transform(trsf)
                    nodes_list.append((gp.X(), gp.Y(), gp.Z()))
                for i in range(1, tri.NbTriangles() + 1):
                    t = tri.Triangle(i)
                    a, b, c = t.Value(1) - 1, t.Value(2) - 1, t.Value(3) - 1
                    if flip:
                        b, c = c, b
                    tris_list.append((offset + a, offset + b, offset + c))
                    tri_fids.append(fid)

            self._tess_cache[key] = (
                np.asarray(nodes_list, dtype=np.float64).reshape(-1, 3),
                np.asarray(tris_list, dtype=np.int64).reshape(-1, 3),
                np.asarray(tri_fids, dtype=np.int64),
            )

        nodes_arr, tris_arr, fids_arr = self._tess_cache[key]
        if face_ids is None:
            return nodes_arr, tris_arr
        mask = np.isin(fids_arr, np.asarray(face_ids, dtype=np.int64))
        return nodes_arr, tris_arr[mask]

    def faces(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all CAD faces."""
        n = len(self._ensure_faces())
        ids = tuple(range(n))
        return EntityQuery(adapter=self, kind="cad_face", ids=ids, universe=frozenset(ids))

    def edges(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all CAD edges."""
        n = len(self._ensure_edges())
        ids = tuple(range(n))
        return EntityQuery(adapter=self, kind="cad_edge", ids=ids, universe=frozenset(ids))

    def vertices(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all CAD vertices."""
        n = len(self._ensure_vertices())
        ids = tuple(range(n))
        return EntityQuery(adapter=self, kind="cad_vertex", ids=ids, universe=frozenset(ids))

    def solids(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all CAD solids."""
        n = len(self._ensure_solids())
        ids = tuple(range(n))
        return EntityQuery(adapter=self, kind="cad_solid", ids=ids, universe=frozenset(ids))


@dataclass(frozen=True)
class _SurfaceFaceEntry:
    """Internal representation of an exterior or explicit surface face."""

    face_nodes: tuple[int, ...]  # 3 (tri) or 4 (quad) corner node indices
    parent_elem: int  # -1 if explicit surface tri without parent element
    normal: np.ndarray
    area: float
    centroid: np.ndarray


class MeshAdapter:
    """Adapter wrapping a finite element mesh (``MeshModel`` or dict/arrays).

    Parameters
    ----------
    mesh:
        A ``MeshModel`` instance, a ``dict`` containing ``"nodes"`` and
        connectivity (``"elements"`` + ``"elem_type"``, ``"blocks"``, or
        canonical type keys such as ``"hex8"``, ``"tet4"``, ``"wedge6"``),
        or an ``(N, 3)`` node coordinate array when ``elements``/``tris``
        are passed as keyword arguments.
    elements:
        Optional ``(M, K)`` element connectivity array when ``mesh`` is a node array.
    elem_type:
        Optional canonical element type (``"hex8"``, ``"tet4"``, ``"tet10"``,
        ``"hex20"``, ``"wedge6"``).
    tris:
        Optional explicit surface triangle connectivity ``(T, 3)``.
    """

    def __init__(
        self,
        mesh: Any,
        *,
        elements: np.ndarray | None = None,
        elem_type: str | None = None,
        tris: np.ndarray | None = None,
    ) -> None:
        self.sets = SetRegistry()
        self._blocks: list[tuple[str, np.ndarray]] = []
        self._explicit_tris: np.ndarray | None = None

        if isinstance(mesh, dict):
            self._nodes = np.asarray(mesh["nodes"], dtype=np.float64).reshape(-1, 3)
            if "blocks" in mesh:
                for blk in mesh["blocks"]:
                    et = str(blk["elem_type"]).lower() if isinstance(blk, dict) else str(blk.elem_type).lower()
                    cn = np.asarray(blk["conn"] if isinstance(blk, dict) else blk.conn, dtype=np.int64)
                    self._blocks.append((et, cn))
            elif "elements" in mesh:
                et = str(mesh.get("elem_type", elem_type or "")).lower()
                if not et:
                    raise ValueError("dict mesh with 'elements' must specify 'elem_type'")
                self._blocks.append((et, np.asarray(mesh["elements"], dtype=np.int64)))
            else:
                for et_key in _SOLID_FACE_TABLES:
                    if et_key in mesh:
                        self._blocks.append((et_key, np.asarray(mesh[et_key], dtype=np.int64)))
            if "tris" in mesh:
                self._explicit_tris = np.asarray(mesh["tris"], dtype=np.int64).reshape(-1, 3)
        elif hasattr(mesh, "nodes") and hasattr(mesh, "blocks"):
            # MeshModel or compatible dataclass
            self._nodes = np.asarray(mesh.nodes, dtype=np.float64).reshape(-1, 3)
            for blk in mesh.blocks:
                self._blocks.append((str(blk.elem_type).lower(), np.asarray(blk.conn, dtype=np.int64)))
            for sname, nids in getattr(mesh, "node_sets", {}).items():
                self.sets.define_set(
                    sname,
                    SelectionSet.from_ids("node", nids, name=sname, universe=range(len(self._nodes))),
                )
            total_e = sum(len(c) for _, c in self._blocks)
            for sname, eids in getattr(mesh, "elem_sets", {}).items():
                self.sets.define_set(
                    sname,
                    SelectionSet.from_ids("element", eids, name=sname, universe=range(total_e)),
                )
        else:
            self._nodes = np.asarray(mesh, dtype=np.float64).reshape(-1, 3)
            if elements is not None:
                et = (elem_type or "").lower()
                if not et:
                    raise ValueError("elem_type is required when passing elements array")
                self._blocks.append((et, np.asarray(elements, dtype=np.int64)))
            if tris is not None:
                self._explicit_tris = np.asarray(tris, dtype=np.int64).reshape(-1, 3)

        # Global element table (for uniform 0-based element indexing)
        self._elem_ranges: list[tuple[str, np.ndarray, int, int]] = []
        offset = 0
        for et, cn in self._blocks:
            n_e = len(cn)
            self._elem_ranges.append((et, cn, offset, offset + n_e))
            offset += n_e
        self._n_elements = offset

        self._node_records: list[_EntityRecord] | None = None
        self._elem_records: list[_EntityRecord] | None = None
        self._free_surfaces: list[_SurfaceFaceEntry] | None = None
        self._surf_records: list[_EntityRecord] | None = None
        self._surf_adj: dict[int, list[tuple[int, float]]] | None = None

    @property
    def node_coords(self) -> np.ndarray:
        """Return ``(N, 3)`` node coordinates array."""
        return self._nodes

    def define_set(self, name: str, selection_or_query: SelectionSet | EntityQuery) -> SelectionSet:
        """Register a named :class:`SelectionSet` in this adapter's registry."""
        return self.sets.define_set(name, selection_or_query)

    def _ensure_free_surfaces(self) -> list[_SurfaceFaceEntry]:
        if self._free_surfaces is not None:
            return self._free_surfaces

        if not self._elem_ranges and self._explicit_tris is not None:
            entries: list[_SurfaceFaceEntry] = []
            for t in self._explicit_tris:
                pts = self._nodes[t]
                cross = np.cross(pts[1] - pts[0], pts[2] - pts[0])
                nrm = float(np.linalg.norm(cross))
                area = 0.5 * nrm
                unit_n = cross / nrm if nrm > 1e-15 else np.array([0.0, 0.0, 1.0])
                entries.append(
                    _SurfaceFaceEntry(
                        face_nodes=(int(t[0]), int(t[1]), int(t[2])),
                        parent_elem=-1,
                        normal=unit_n,
                        area=area,
                        centroid=pts.mean(axis=0),
                    )
                )
            self._free_surfaces = entries
            return entries

        face_counts: dict[tuple[int, ...], int] = {}
        face_meta: dict[tuple[int, ...], tuple[tuple[int, ...], int]] = {}

        for et, conn, start_eid, _ in self._elem_ranges:
            if et not in _SOLID_FACE_TABLES:
                raise ValueError(f"Unsupported solid element type {et!r} for free surface extraction")
            table = _SOLID_FACE_TABLES[et]
            for local_e, row in enumerate(conn):
                geid = start_eid + local_e
                for local_face in table:
                    fnodes = tuple(int(row[k]) for k in local_face)
                    key = tuple(sorted(fnodes))
                    face_counts[key] = face_counts.get(key, 0) + 1
                    if key not in face_meta:
                        face_meta[key] = (fnodes, geid)

        entries = []
        for key, count in face_counts.items():
            if count != 1:
                continue
            fnodes, geid = face_meta[key]
            pts = self._nodes[list(fnodes)]
            if len(fnodes) == 3:
                cross = np.cross(pts[1] - pts[0], pts[2] - pts[0])
                nrm = float(np.linalg.norm(cross))
                area = 0.5 * nrm
                unit_n = cross / nrm if nrm > 1e-15 else np.array([0.0, 0.0, 1.0])
            else:
                # Quad face: split into (0, 1, 2) and (0, 2, 3)
                c1 = np.cross(pts[1] - pts[0], pts[2] - pts[0])
                c2 = np.cross(pts[2] - pts[0], pts[3] - pts[0])
                c_sum = c1 + c2
                nrm = float(np.linalg.norm(c_sum))
                area = 0.5 * (float(np.linalg.norm(c1)) + float(np.linalg.norm(c2)))
                unit_n = c_sum / nrm if nrm > 1e-15 else np.array([0.0, 0.0, 1.0])
            entries.append(
                _SurfaceFaceEntry(
                    face_nodes=fnodes,
                    parent_elem=geid,
                    normal=unit_n,
                    area=area,
                    centroid=pts.mean(axis=0),
                )
            )

        self._free_surfaces = entries
        return entries

    def _get_records(self, kind: SetKind) -> list[_EntityRecord]:
        if kind == "node":
            if self._node_records is None:
                self._node_records = [
                    _EntityRecord(entity_id=i, centroid=self._nodes[i], bbox_lo=self._nodes[i], bbox_hi=self._nodes[i])
                    for i in range(len(self._nodes))
                ]
            return self._node_records

        if kind == "element":
            if self._elem_records is None:
                recs: list[_EntityRecord] = []
                for _, conn, start_eid, _ in self._elem_ranges:
                    for local_e, row in enumerate(conn):
                        pts = self._nodes[row]
                        recs.append(
                            _EntityRecord(
                                entity_id=start_eid + local_e,
                                centroid=pts.mean(axis=0),
                                bbox_lo=pts.min(axis=0),
                                bbox_hi=pts.max(axis=0),
                            )
                        )
                self._elem_records = recs
            return self._elem_records

        if kind == "surface_tri":
            if self._surf_records is None:
                entries = self._ensure_free_surfaces()
                recs = []
                for idx, ent in enumerate(entries):
                    pts = self._nodes[list(ent.face_nodes)]
                    recs.append(
                        _EntityRecord(
                            entity_id=idx,
                            centroid=ent.centroid,
                            normal=ent.normal,
                            area=ent.area,
                            bbox_lo=pts.min(axis=0),
                            bbox_hi=pts.max(axis=0),
                        )
                    )
                self._surf_records = recs
            return self._surf_records

        raise ValueError(f"Unsupported Mesh entity kind: {kind!r}")

    def _get_surface_adjacency(self) -> dict[int, list[tuple[int, float]]]:
        if self._surf_adj is not None:
            return self._surf_adj
        entries = self._ensure_free_surfaces()
        edge_to_faces: dict[tuple[int, int], list[int]] = {}
        for fid, ent in enumerate(entries):
            fn = ent.face_nodes
            m = len(fn)
            for i in range(m):
                a, b = fn[i], fn[(i + 1) % m]
                ekey = (a, b) if a < b else (b, a)
                edge_to_faces.setdefault(ekey, []).append(fid)

        adj: dict[int, list[tuple[int, float]]] = {i: [] for i in range(len(entries))}
        for fids in edge_to_faces.values():
            if len(fids) < 2:
                continue
            for i in range(len(fids)):
                for j in range(i + 1, len(fids)):
                    f1, f2 = fids[i], fids[j]
                    dot = float(np.clip(np.dot(entries[f1].normal, entries[f2].normal), -1.0, 1.0))
                    ang = float(np.degrees(np.arccos(dot)))
                    adj[f1].append((f2, ang))
                    adj[f2].append((f1, ang))
        self._surf_adj = adj
        return adj

    def to_surface_tris(self, face_ids: Sequence[int] | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(nodes, tris)`` for selected surface face IDs (splitting quads)."""
        entries = self._ensure_free_surfaces()
        ids = range(len(entries)) if face_ids is None else face_ids
        tris: list[tuple[int, int, int]] = []
        for fid in ids:
            fn = entries[fid].face_nodes
            if len(fn) == 3:
                tris.append((fn[0], fn[1], fn[2]))
            elif len(fn) == 4:
                tris.append((fn[0], fn[1], fn[2]))
                tris.append((fn[0], fn[2], fn[3]))
        return self._nodes, np.asarray(tris, dtype=np.int64).reshape(-1, 3)

    def nodes(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all mesh nodes."""
        ids = tuple(range(len(self._nodes)))
        return EntityQuery(adapter=self, kind="node", ids=ids, universe=frozenset(ids))

    def elements(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over all mesh elements."""
        ids = tuple(range(self._n_elements))
        return EntityQuery(adapter=self, kind="element", ids=ids, universe=frozenset(ids))

    def surface_tris(self) -> EntityQuery:
        """Return an :class:`EntityQuery` over exterior/surface faces."""
        entries = self._ensure_free_surfaces()
        ids = tuple(range(len(entries)))
        return EntityQuery(adapter=self, kind="surface_tri", ids=ids, universe=frozenset(ids))

    def free_surface_faces(self) -> EntityQuery:
        """Extract all exterior (non-shared) element faces from the solid mesh."""
        return self.surface_tris()

    def free_surface_elements(self) -> EntityQuery:
        """Return unique solid elements that border the exterior free surface."""
        entries = self._ensure_free_surfaces()
        eids = tuple(sorted({e.parent_elem for e in entries if e.parent_elem >= 0}))
        return EntityQuery(
            adapter=self,
            kind="element",
            ids=eids,
            universe=frozenset(range(self._n_elements)),
        )

    def free_surface_nodes(self) -> EntityQuery:
        """Return unique nodes lying on the exterior free surface."""
        entries = self._ensure_free_surfaces()
        nids = tuple(sorted({nid for e in entries for nid in e.face_nodes}))
        return EntityQuery(
            adapter=self,
            kind="node",
            ids=nids,
            universe=frozenset(range(len(self._nodes))),
        )


class EntityQuery:
    """Fluent, domain-unified query over CAD or Mesh entities.

    Parameters
    ----------
    adapter:
        Owning :class:`CadAdapter` or :class:`MeshAdapter`.
    kind:
        Entity kind (``"cad_face"``, ``"cad_edge"``, ``"cad_vertex"``,
        ``"cad_solid"``, ``"node"``, ``"element"``, ``"surface_tri"``).
    ids:
        Currently selected 0-based entity IDs.
    universe:
        Full set of valid entity IDs for complement (``~``) operations.
    """

    def __init__(
        self,
        adapter: CadAdapter | MeshAdapter,
        kind: SetKind,
        ids: Iterable[int],
        universe: frozenset[int],
    ) -> None:
        self.adapter = adapter
        self.kind = kind
        self.ids: tuple[int, ...] = tuple(sorted({int(i) for i in ids}))
        self.universe = universe

    @property
    def count(self) -> int:
        """Number of selected entities."""
        return len(self.ids)

    def __len__(self) -> int:
        return len(self.ids)

    def __iter__(self):
        return iter(self.ids)

    def __contains__(self, item: object) -> bool:
        return isinstance(item, (int, np.integer)) and int(item) in set(self.ids)

    def __bool__(self) -> bool:
        return len(self.ids) > 0

    def __getitem__(self, idx: int | slice) -> int | tuple[int, ...]:
        return self.ids[idx]

    def _spawn(self, new_ids: Iterable[int]) -> EntityQuery:
        return EntityQuery(
            adapter=self.adapter,
            kind=self.kind,
            ids=new_ids,
            universe=self.universe,
        )

    def _records(self) -> list[_EntityRecord]:
        return self.adapter._get_records(self.kind)

    # ------------------------------------------------------------------
    # Geometric & spatial filters
    # ------------------------------------------------------------------

    def by_geom_type(
        self,
        geom_type: str,
        radius: float | tuple[float, float] | None = None,
        *,
        tol: float = 1e-6,
    ) -> EntityQuery:
        """Filter entities by underlying OCP surface/curve type and optional radius.

        Parameters
        ----------
        geom_type:
            Target geometry type name (e.g. ``"plane"``, ``"cylinder"``,
            ``"sphere"``, ``"cone"``, ``"torus"``, ``"bspline"``,
            ``"circle"``, ``"line"``).
        radius:
            Optional target radius (scalar ``float`` matched within ``tol``)
            or inclusive ``(r_min, r_max)`` range.
        tol:
            Absolute tolerance when ``radius`` is a scalar.
        """
        target = geom_type.strip().lower()
        recs = self._records()
        out: list[int] = []
        for eid in self.ids:
            r = recs[eid]
            if r.geom_type != target:
                continue
            if radius is not None:
                if r.radius is None:
                    continue
                if isinstance(radius, (tuple, list)):
                    if not (float(radius[0]) - tol <= r.radius <= float(radius[1]) + tol):
                        continue
                else:
                    if abs(r.radius - float(radius)) > max(tol, 1e-6 * abs(float(radius))):
                        continue
            out.append(eid)
        return self._spawn(out)

    def within_box(
        self,
        lo: Sequence[float] | np.ndarray,
        hi: Sequence[float] | np.ndarray,
        *,
        tol: float = 1e-9,
    ) -> EntityQuery:
        """Filter entities whose centroid lies inside axis-aligned box ``[lo, hi]``."""
        lo_a = np.asarray(lo, dtype=np.float64).reshape(3) - tol
        hi_a = np.asarray(hi, dtype=np.float64).reshape(3) + tol
        recs = self._records()
        out = [
            eid
            for eid in self.ids
            if np.all(recs[eid].centroid >= lo_a) and np.all(recs[eid].centroid <= hi_a)
        ]
        return self._spawn(out)

    def within_sphere(
        self,
        center: Sequence[float] | np.ndarray,
        r: float,
        *,
        tol: float = 1e-9,
    ) -> EntityQuery:
        """Filter entities whose centroid lies within distance ``r`` of ``center``."""
        c = np.asarray(center, dtype=np.float64).reshape(3)
        r_lim = float(r) + tol
        recs = self._records()
        out = [eid for eid in self.ids if float(np.linalg.norm(recs[eid].centroid - c)) <= r_lim]
        return self._spawn(out)

    def within_cylinder(
        self,
        p0: Sequence[float] | np.ndarray,
        p1: Sequence[float] | np.ndarray,
        r: float,
        *,
        tol: float = 1e-9,
    ) -> EntityQuery:
        """Filter entities whose centroid lies inside finite cylinder ``(p0, p1, r)``."""
        a = np.asarray(p0, dtype=np.float64).reshape(3)
        b = np.asarray(p1, dtype=np.float64).reshape(3)
        axis = b - a
        axis_len = float(np.linalg.norm(axis))
        if axis_len <= 1e-15:
            raise ValueError("within_cylinder requires distinct endpoints p0 and p1")
        u = axis / axis_len
        r_lim = float(r) + tol
        recs = self._records()
        out: list[int] = []
        for eid in self.ids:
            v = recs[eid].centroid - a
            proj = float(np.dot(v, u))
            if proj < -tol or proj > axis_len + tol:
                continue
            radial = float(np.linalg.norm(v - proj * u))
            if radial <= r_lim:
                out.append(eid)
        return self._spawn(out)

    def on_plane(
        self,
        origin: Sequence[float] | np.ndarray,
        normal: Sequence[float] | np.ndarray,
        *,
        tol: float = 1e-6,
    ) -> EntityQuery:
        """Filter entities whose centroid lies on plane ``(origin, normal)`` within ``tol``."""
        o = np.asarray(origin, dtype=np.float64).reshape(3)
        n = np.asarray(normal, dtype=np.float64).reshape(3)
        n_norm = float(np.linalg.norm(n))
        if n_norm <= 1e-15:
            raise ValueError("on_plane requires a non-zero normal vector")
        n_hat = n / n_norm
        recs = self._records()
        out = [eid for eid in self.ids if abs(float(np.dot(recs[eid].centroid - o, n_hat))) <= tol]
        return self._spawn(out)

    def between_planes(
        self,
        origin: Sequence[float] | np.ndarray,
        normal: Sequence[float] | np.ndarray,
        d_min: float,
        d_max: float,
        *,
        tol: float = 1e-9,
    ) -> EntityQuery:
        """Filter entities whose signed distance to ``(origin, normal)`` is in ``[d_min, d_max]``."""
        o = np.asarray(origin, dtype=np.float64).reshape(3)
        n = np.asarray(normal, dtype=np.float64).reshape(3)
        n_norm = float(np.linalg.norm(n))
        if n_norm <= 1e-15:
            raise ValueError("between_planes requires a non-zero normal vector")
        n_hat = n / n_norm
        lo_d = float(d_min) - tol
        hi_d = float(d_max) + tol
        recs = self._records()
        out = [
            eid
            for eid in self.ids
            if lo_d <= float(np.dot(recs[eid].centroid - o, n_hat)) <= hi_d
        ]
        return self._spawn(out)

    def by_normal(
        self,
        direction: Sequence[float] | np.ndarray | None = None,
        *,
        tol_deg: float = 10.0,
        perp: Sequence[float] | np.ndarray | None = None,
    ) -> EntityQuery:
        """Filter faces/surfaces by normal alignment or perpendicularity.

        Parameters
        ----------
        direction:
            Target unit normal direction; selects entities within ``tol_deg``.
        tol_deg:
            Angular tolerance in degrees.
        perp:
            If provided, selects entities whose normal is perpendicular to ``perp``
            within ``tol_deg`` (i.e. ``|dot(n, perp_hat)| <= sin(tol_deg)``).
        """
        if direction is None and perp is None:
            raise ValueError("by_normal requires either direction or perp")
        recs = self._records()
        out: list[int] = []

        if direction is not None:
            d = np.asarray(direction, dtype=np.float64).reshape(3)
            d_norm = float(np.linalg.norm(d))
            if d_norm <= 1e-15:
                raise ValueError("by_normal direction must be non-zero")
            d_hat = d / d_norm
            cos_thr = float(np.cos(np.radians(tol_deg)))
            for eid in self.ids:
                n = recs[eid].normal
                if n is not None and float(np.dot(n, d_hat)) >= cos_thr:
                    out.append(eid)
        else:
            p = np.asarray(perp, dtype=np.float64).reshape(3)
            p_norm = float(np.linalg.norm(p))
            if p_norm <= 1e-15:
                raise ValueError("by_normal perp must be non-zero")
            p_hat = p / p_norm
            sin_thr = float(np.sin(np.radians(tol_deg)))
            for eid in self.ids:
                n = recs[eid].normal
                if n is not None and abs(float(np.dot(n, p_hat))) <= sin_thr:
                    out.append(eid)
        return self._spawn(out)

    def by_area(self, min_area: float = 0.0, max_area: float = float("inf")) -> EntityQuery:
        """Filter faces/surfaces by area in ``[min_area, max_area]``."""
        recs = self._records()
        out = [eid for eid in self.ids if min_area <= recs[eid].area <= max_area]
        return self._spawn(out)

    def by_length(self, min_length: float = 0.0, max_length: float = float("inf")) -> EntityQuery:
        """Filter edges by length in ``[min_length, max_length]``."""
        recs = self._records()
        out = [eid for eid in self.ids if min_length <= recs[eid].length <= max_length]
        return self._spawn(out)

    def nearest_to(self, pt: Sequence[float] | np.ndarray, k: int = 1) -> EntityQuery:
        """Return the ``k`` nearest entities to point ``pt``."""
        if not self.ids or k <= 0:
            return self._spawn(())
        p = np.asarray(pt, dtype=np.float64).reshape(3)
        recs = self._records()
        dists = [(float(np.linalg.norm(recs[eid].centroid - p)), eid) for eid in self.ids]
        dists.sort(key=lambda item: item[0])
        return self._spawn(eid for _, eid in dists[: int(k)])

    def flood_fill(self, angle_deg: float = 30.0) -> EntityQuery:
        """Expand selection across adjacent faces with dihedral angle ``<= angle_deg``."""
        if self.kind not in ("cad_face", "surface_tri"):
            raise ValueError(f"flood_fill requires cad_face or surface_tri, got {self.kind!r}")
        if not self.ids:
            return self._spawn(())

        if isinstance(self.adapter, CadAdapter):
            adj = self.adapter._get_face_adjacency()
        else:
            adj = self.adapter._get_surface_adjacency()

        visited = set(self.ids)
        queue: deque[int] = deque(self.ids)
        while queue:
            cur = queue.popleft()
            for nbr, ang in adj.get(cur, ()):
                if nbr not in visited and ang <= angle_deg:
                    visited.add(nbr)
                    queue.append(nbr)
        return self._spawn(visited)

    # ------------------------------------------------------------------
    # Boolean operators
    # ------------------------------------------------------------------

    def _other_ids(self, other: EntityQuery | SelectionSet) -> frozenset[int]:
        if isinstance(other, EntityQuery):
            if other.kind != self.kind:
                raise ValueError(f"Kind mismatch: {self.kind!r} vs {other.kind!r}")
            return frozenset(other.ids)
        if isinstance(other, SelectionSet):
            if other.kind != self.kind:
                raise ValueError(f"Kind mismatch: {self.kind!r} vs {other.kind!r}")
            return other.id_set
        raise TypeError(f"Expected EntityQuery or SelectionSet, got {type(other).__name__}")

    def __or__(self, other: EntityQuery | SelectionSet) -> EntityQuery:
        return self._spawn(frozenset(self.ids) | self._other_ids(other))

    def __and__(self, other: EntityQuery | SelectionSet) -> EntityQuery:
        return self._spawn(frozenset(self.ids) & self._other_ids(other))

    def __sub__(self, other: EntityQuery | SelectionSet) -> EntityQuery:
        return self._spawn(frozenset(self.ids) - self._other_ids(other))

    def __xor__(self, other: EntityQuery | SelectionSet) -> EntityQuery:
        return self._spawn(frozenset(self.ids) ^ self._other_ids(other))

    def __invert__(self) -> EntityQuery:
        return self._spawn(self.universe - frozenset(self.ids))

    # ------------------------------------------------------------------
    # Conversion & CAD <-> FEM bridge
    # ------------------------------------------------------------------

    def to_set(self, name: str = "") -> SelectionSet:
        """Convert this query into a :class:`SelectionSet` of the same ``kind``."""
        return SelectionSet.from_ids(self.kind, self.ids, name=name, universe=self.universe)

    def _resolve_mesh_and_tol(
        self,
        mesh: MeshAdapter | None,
        tol: float | None,
    ) -> tuple[MeshAdapter, float]:
        if isinstance(self.adapter, MeshAdapter):
            m = self.adapter
        else:
            m = mesh if mesh is not None else self.adapter._bound_mesh
            if m is None:
                raise ValueError("CAD-to-FEM mapping requires a MeshAdapter (pass mesh=... or call bind_mesh)")
        if tol is None:
            if isinstance(self.adapter, CadAdapter) and self.adapter._bound_tol is not None:
                tol = self.adapter._bound_tol
            else:
                pts = m.node_coords
                diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0))) if len(pts) > 0 else 1.0
                tol = max(1e-3 * diag, 1e-4)
        return m, float(tol)

    def _map_cad_to_node_ids(self, mesh: MeshAdapter, tol: float) -> list[int]:
        """Find mesh node IDs within ``tol`` of the selected CAD entities."""
        assert isinstance(self.adapter, CadAdapter)
        if not self.ids:
            return []

        node_pts = mesh.node_coords
        if len(node_pts) == 0:
            return []

        if self.kind == "cad_vertex":
            recs = self._records()
            vpts = np.asarray([recs[i].centroid for i in self.ids], dtype=np.float64)
            tree = cKDTree(vpts)
            dists, _ = tree.query(node_pts)
            return [int(i) for i in np.where(dists <= tol)[0]]

        if self.kind == "cad_face":
            # Coarse filter by face tessellation + exact point-to-triangle distance
            t_nodes, t_tris = self.adapter.tessellate_faces(self.ids)
            if len(t_tris) == 0:
                return []
            centroids = t_nodes[t_tris].mean(axis=1)
            tree = cKDTree(centroids)
            k = min(8, len(t_tris))
            _, cand = tree.query(node_pts, k=k)
            if k == 1:
                cand = cand[:, None]
            matched: list[int] = []
            for nid, p in enumerate(node_pts):
                for ti in cand[nid]:
                    t = t_tris[ti]
                    q = _closest_point_on_triangle(p, t_nodes[t[0]], t_nodes[t[1]], t_nodes[t[2]])
                    if float(np.linalg.norm(p - q)) <= tol:
                        matched.append(nid)
                        break
            return matched

        # cad_edge or cad_solid: BRepExtrema_DistShapeShape per candidate vertex
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.gp import gp_Pnt

        shapes = self.adapter._get_shapes(self.kind, self.ids)
        matched_ids: list[int] = []
        for nid, p in enumerate(node_pts):
            vtx = BRepBuilderAPI_MakeVertex(gp_Pnt(float(p[0]), float(p[1]), float(p[2]))).Vertex()
            for shp in shapes:
                dcalc = BRepExtrema_DistShapeShape(vtx, shp)
                dcalc.Perform()
                if dcalc.IsDone() and float(dcalc.Value()) <= tol:
                    matched_ids.append(nid)
                    break
        return matched_ids

    def to_node_set(
        self,
        name: str = "",
        *,
        mesh: MeshAdapter | None = None,
        tol: float | None = None,
    ) -> SelectionSet:
        """Map selected CAD or Mesh entities to a ``"node"`` :class:`SelectionSet`."""
        m, resolved_tol = self._resolve_mesh_and_tol(mesh, tol)
        node_uni = frozenset(range(len(m.node_coords)))

        if isinstance(self.adapter, MeshAdapter):
            if self.kind == "node":
                return SelectionSet.from_ids("node", self.ids, name=name, universe=node_uni)
            if self.kind == "surface_tri":
                entries = m._ensure_free_surfaces()
                nids = {nid for fid in self.ids for nid in entries[fid].face_nodes}
                return SelectionSet.from_ids("node", nids, name=name, universe=node_uni)
            if self.kind == "element":
                selected = set(self.ids)
                nids: set[int] = set()
                for _, conn, start_eid, _ in m._elem_ranges:
                    for local_e, row in enumerate(conn):
                        if (start_eid + local_e) in selected:
                            nids.update(int(x) for x in row)
                return SelectionSet.from_ids("node", nids, name=name, universe=node_uni)

        nids_cad = self._map_cad_to_node_ids(m, resolved_tol)
        return SelectionSet.from_ids("node", nids_cad, name=name, universe=node_uni)

    def to_element_set(
        self,
        name: str = "",
        *,
        mesh: MeshAdapter | None = None,
        tol: float | None = None,
    ) -> SelectionSet:
        """Map selected CAD or Mesh entities to an ``"element"`` :class:`SelectionSet`."""
        m, resolved_tol = self._resolve_mesh_and_tol(mesh, tol)
        elem_uni = frozenset(range(m._n_elements))

        if isinstance(self.adapter, MeshAdapter):
            if self.kind == "element":
                return SelectionSet.from_ids("element", self.ids, name=name, universe=elem_uni)
            if self.kind == "surface_tri":
                entries = m._ensure_free_surfaces()
                eids = {entries[fid].parent_elem for fid in self.ids if entries[fid].parent_elem >= 0}
                return SelectionSet.from_ids("element", eids, name=name, universe=elem_uni)
            if self.kind == "node":
                nset = set(self.ids)
                eids = set()
                for _, conn, start_eid, _ in m._elem_ranges:
                    for local_e, row in enumerate(conn):
                        if any(int(x) in nset for x in row):
                            eids.add(start_eid + local_e)
                return SelectionSet.from_ids("element", eids, name=name, universe=elem_uni)

        node_set = self.to_node_set(mesh=m, tol=resolved_tol).id_set
        eids_cad: set[int] = set()
        for _, conn, start_eid, _ in m._elem_ranges:
            for local_e, row in enumerate(conn):
                if any(int(x) in node_set for x in row):
                    eids_cad.add(start_eid + local_e)
        return SelectionSet.from_ids("element", eids_cad, name=name, universe=elem_uni)

    def to_surface_set(
        self,
        name: str = "",
        *,
        mesh: MeshAdapter | None = None,
        tol: float | None = None,
    ) -> SelectionSet:
        """Map selected entities to a ``"surface_tri"`` (or ``"cad_face"``) :class:`SelectionSet`."""
        if isinstance(self.adapter, CadAdapter) and mesh is None and self.adapter._bound_mesh is None:
            return SelectionSet.from_ids(self.kind, self.ids, name=name, universe=self.universe)

        m, resolved_tol = self._resolve_mesh_and_tol(mesh, tol)
        entries = m._ensure_free_surfaces()
        surf_uni = frozenset(range(len(entries)))

        if isinstance(self.adapter, MeshAdapter):
            if self.kind == "surface_tri":
                return SelectionSet.from_ids("surface_tri", self.ids, name=name, universe=surf_uni)
            if self.kind == "node":
                nset = set(self.ids)
                fids = [fid for fid, ent in enumerate(entries) if all(n in nset for n in ent.face_nodes)]
                return SelectionSet.from_ids("surface_tri", fids, name=name, universe=surf_uni)
            if self.kind == "element":
                eset = set(self.ids)
                fids = [fid for fid, ent in enumerate(entries) if ent.parent_elem in eset]
                return SelectionSet.from_ids("surface_tri", fids, name=name, universe=surf_uni)

        node_set = self.to_node_set(mesh=m, tol=resolved_tol).id_set
        fids_cad = [fid for fid, ent in enumerate(entries) if all(n in node_set for n in ent.face_nodes)]
        return SelectionSet.from_ids("surface_tri", fids_cad, name=name, universe=surf_uni)

    # ------------------------------------------------------------------
    # Measurement & Contact integration
    # ------------------------------------------------------------------

    def to_surface_tris(self, deflection: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(nodes, tris)`` representation of the selected surface entities."""
        if isinstance(self.adapter, CadAdapter):
            if self.kind != "cad_face":
                raise ValueError(f"to_surface_tris requires cad_face, got {self.kind!r}")
            return self.adapter.tessellate_faces(self.ids, deflection=deflection)
        if self.kind != "surface_tri":
            raise ValueError(f"to_surface_tris requires surface_tri, got {self.kind!r}")
        return self.adapter.to_surface_tris(self.ids)

    def measure(self) -> MeasureResult:
        """Compute geometric measurements (:class:`MeasureResult`) for selected entities."""
        if not self.ids:
            raise ValueError("Cannot measure an empty EntityQuery")
        if isinstance(self.adapter, CadAdapter):
            shapes = self.adapter._get_shapes(self.kind, self.ids)
            return measure_cad(shapes)

        # MeshAdapter
        if self.kind == "node":
            return measure_mesh(self.adapter.node_coords, node_indices=np.asarray(self.ids, dtype=np.int64))
        if self.kind == "surface_tri":
            nodes, tris = self.adapter.to_surface_tris(self.ids)
            return measure_mesh(nodes, tris=tris)
        if self.kind == "element":
            if len(self.adapter._elem_ranges) == 1:
                et, conn, _, _ = self.adapter._elem_ranges[0]
                return measure_mesh(
                    self.adapter.node_coords,
                    elements=conn[list(self.ids)],
                    elem_type=et,
                )
            # Multi-block fallback: measure surface or first matching block
            et, conn, start_eid, end_eid = self.adapter._elem_ranges[0]
            local = [i - start_eid for i in self.ids if start_eid <= i < end_eid]
            return measure_mesh(self.adapter.node_coords, elements=conn[local], elem_type=et)
        raise ValueError(f"Unsupported kind for measure(): {self.kind!r}")
