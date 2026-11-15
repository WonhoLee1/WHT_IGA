"""Geometric measurement utilities for CAD shapes and finite element meshes (v2.1 §8.5).

Provides :class:`BBoxResult`, :class:`MeasureResult`, :func:`measure_cad`,
and :func:`measure_mesh` for computing bounding boxes, lengths, areas,
volumes, centroids, inertia tensors, and minimum distances.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class BBoxResult:
    """Axis-aligned bounding box summary.

    Parameters
    ----------
    min:
        Minimum ``(x, y, z)`` coordinates as shape ``(3,)`` array.
    max:
        Maximum ``(x, y, z)`` coordinates as shape ``(3,)`` array.
    center:
        Bounding box center ``0.5 * (min + max)`` as shape ``(3,)`` array.
    extents:
        Side lengths ``max - min`` as shape ``(3,)`` array.
    diagonal:
        Euclidean diagonal length ``||max - min||``.
    """

    min: np.ndarray
    max: np.ndarray
    center: np.ndarray
    extents: np.ndarray
    diagonal: float

    @classmethod
    def from_bounds(cls, lo: np.ndarray | Sequence[float], hi: np.ndarray | Sequence[float]) -> BBoxResult:
        """Construct a :class:`BBoxResult` from ``lo`` and ``hi`` coordinate vectors."""
        lo_arr = np.asarray(lo, dtype=np.float64).reshape(3)
        hi_arr = np.asarray(hi, dtype=np.float64).reshape(3)
        extents = hi_arr - lo_arr
        return cls(
            min=lo_arr,
            max=hi_arr,
            center=0.5 * (lo_arr + hi_arr),
            extents=extents,
            diagonal=float(np.linalg.norm(extents)),
        )


class MeasureResult:
    """Unified geometric measurement container for CAD and Mesh selections.

    Parameters
    ----------
    bbox_result:
        Precomputed :class:`BBoxResult`.
    length_val:
        Total curve/edge length (``0.0`` if not applicable).
    area_val:
        Total surface area (``0.0`` if not applicable).
    volume_val:
        Total solid volume (``0.0`` if not applicable).
    centroid_val:
        Geometric centroid ``(3,)`` array.
    inertia_val:
        ``(3, 3)`` inertia tensor about the centroid (unit density).
    sample_points:
        Optional ``(N, 3)`` point cloud used as fallback for distance queries.
    ocp_shape:
        Optional OpenCASCADE ``TopoDS_Shape`` for exact BRep distance queries.
    """

    def __init__(
        self,
        *,
        bbox_result: BBoxResult,
        length_val: float = 0.0,
        area_val: float = 0.0,
        volume_val: float = 0.0,
        centroid_val: np.ndarray | Sequence[float] = (0.0, 0.0, 0.0),
        inertia_val: np.ndarray | None = None,
        sample_points: np.ndarray | None = None,
        ocp_shape: Any = None,
    ) -> None:
        self._bbox = bbox_result
        self._length = float(length_val)
        self._area = float(area_val)
        self._volume = float(volume_val)
        self._centroid = np.asarray(centroid_val, dtype=np.float64).reshape(3)
        self._inertia = (
            np.asarray(inertia_val, dtype=np.float64).reshape(3, 3)
            if inertia_val is not None
            else np.zeros((3, 3), dtype=np.float64)
        )
        self._sample_points = (
            np.asarray(sample_points, dtype=np.float64).reshape(-1, 3)
            if sample_points is not None and len(sample_points) > 0
            else self._centroid.reshape(1, 3)
        )
        self._ocp_shape = ocp_shape

    def bbox(self) -> BBoxResult:
        """Return the axis-aligned bounding box :class:`BBoxResult`."""
        return self._bbox

    def length(self) -> float:
        """Return total edge/curve length."""
        return self._length

    def area(self) -> float:
        """Return total surface area."""
        return self._area

    def volume(self) -> float:
        """Return total solid volume."""
        return self._volume

    def centroid(self) -> np.ndarray:
        """Return geometric centroid as a ``(3,)`` float64 array."""
        return self._centroid.copy()

    def inertia_tensor(self) -> np.ndarray:
        """Return ``(3, 3)`` inertia tensor about the centroid (unit density)."""
        return self._inertia.copy()

    def distance_to(self, other: MeasureResult | Any) -> float:
        """Compute minimum Euclidean distance to ``other``.

        Parameters
        ----------
        other:
            Another :class:`MeasureResult`, an :class:`EntityQuery` (with ``.measure()``),
            an OCP ``TopoDS_Shape``, or an array-like of 3D coordinates ``(N, 3)`` or ``(3,)``.
        """
        other_mr = _coerce_to_measure_result(other)
        if self._ocp_shape is not None and other_mr._ocp_shape is not None:
            from OCP.BRepExtrema import BRepExtrema_DistShapeShape

            dist_calc = BRepExtrema_DistShapeShape(self._ocp_shape, other_mr._ocp_shape)
            dist_calc.Perform()
            if dist_calc.IsDone():
                return float(dist_calc.Value())

        tree = cKDTree(self._sample_points)
        dists, _ = tree.query(other_mr._sample_points)
        return float(np.min(dists))


def _coerce_to_measure_result(obj: Any) -> MeasureResult:
    if isinstance(obj, MeasureResult):
        return obj
    if hasattr(obj, "measure") and callable(obj.measure):
        return obj.measure()
    if isinstance(obj, (np.ndarray, list, tuple)):
        pts = np.asarray(obj, dtype=np.float64).reshape(-1, 3)
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
        return MeasureResult(
            bbox_result=BBoxResult.from_bounds(lo, hi),
            centroid_val=pts.mean(axis=0),
            sample_points=pts,
        )
    # Assume TopoDS_Shape
    return measure_cad(obj)


def _gprop_inertia_to_numpy(props: Any) -> np.ndarray:
    mat = props.MatrixOfInertia()
    return np.array(
        [
            [mat.Value(1, 1), mat.Value(1, 2), mat.Value(1, 3)],
            [mat.Value(2, 1), mat.Value(2, 2), mat.Value(2, 3)],
            [mat.Value(3, 1), mat.Value(3, 2), mat.Value(3, 3)],
        ],
        dtype=np.float64,
    )


def _build_compound(shapes: Sequence[Any]) -> Any:
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    if len(shapes) == 1:
        return shapes[0]
    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for s in shapes:
        builder.Add(comp, s)
    return comp


def _sample_shape_points(shape: Any, deflection: float | None = None) -> np.ndarray:
    from OCP.BRep import BRep_Tool
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopAbs import TopAbs_FACE, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.TopLoc import TopLoc_Location
    from OCP.gp import gp_Pnt

    pts: list[tuple[float, float, float]] = []
    # Vertices first (always available)
    vexp = TopExp_Explorer(shape, TopAbs_VERTEX)
    while vexp.More():
        p = BRep_Tool.Pnt_s(TopoDS.Vertex(vexp.Current()))
        pts.append((p.X(), p.Y(), p.Z()))
        vexp.Next()

    # Face triangulation nodes if faces exist
    fexp = TopExp_Explorer(shape, TopAbs_FACE)
    if fexp.More():
        if deflection is None:
            from OCP.Bnd import Bnd_Box
            from OCP.BRepBndLib import BRepBndLib

            box = Bnd_Box()
            BRepBndLib.Add_s(shape, box, True)
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
        BRepMesh_IncrementalMesh(shape, deflection, False, 0.3, True)
        fexp2 = TopExp_Explorer(shape, TopAbs_FACE)
        while fexp2.More():
            face = TopoDS.Face(fexp2.Current())
            loc = TopLoc_Location()
            tri = BRep_Tool.Triangulation_s(face, loc)
            if tri is not None:
                trsf = loc.Transformation()
                for i in range(1, tri.NbNodes() + 1):
                    nd = tri.Node(i)
                    gp = gp_Pnt(nd.X(), nd.Y(), nd.Z())
                    gp.Transform(trsf)
                    pts.append((gp.X(), gp.Y(), gp.Z()))
            fexp2.Next()

    if not pts:
        return np.zeros((0, 3), dtype=np.float64)
    return np.asarray(pts, dtype=np.float64)


def measure_cad(shapes: Any) -> MeasureResult:
    """Measure geometric properties of one or more OpenCASCADE ``TopoDS_Shape`` objects.

    Parameters
    ----------
    shapes:
        A single ``TopoDS_Shape`` or a non-empty sequence of ``TopoDS_Shape`` instances.
    """
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    if isinstance(shapes, (list, tuple)):
        if len(shapes) == 0:
            raise ValueError("measure_cad requires at least one TopoDS_Shape")
        compound = _build_compound(shapes)
    else:
        compound = shapes

    box = Bnd_Box()
    BRepBndLib.Add_s(compound, box, True)
    lo = np.array([box.GetXMin(), box.GetYMin(), box.GetZMin()], dtype=np.float64)
    hi = np.array([box.GetXMax(), box.GetYMax(), box.GetZMax()], dtype=np.float64)
    bbox_res = BBoxResult.from_bounds(lo, hi)

    has_solid = TopExp_Explorer(compound, TopAbs_SOLID).More()
    has_face = TopExp_Explorer(compound, TopAbs_FACE).More()
    has_edge = TopExp_Explorer(compound, TopAbs_EDGE).More()

    length_val = 0.0
    area_val = 0.0
    volume_val = 0.0
    centroid_val = bbox_res.center.copy()
    inertia_val = np.zeros((3, 3), dtype=np.float64)

    if has_edge:
        lprops = GProp_GProps()
        BRepGProp.LinearProperties_s(compound, lprops)
        length_val = float(lprops.Mass())
        if not has_face and not has_solid and length_val > 0.0:
            com = lprops.CentreOfMass()
            centroid_val = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
            inertia_val = _gprop_inertia_to_numpy(lprops)

    if has_face:
        sprops = GProp_GProps()
        BRepGProp.SurfaceProperties_s(compound, sprops)
        area_val = float(sprops.Mass())
        if not has_solid and area_val > 0.0:
            com = sprops.CentreOfMass()
            centroid_val = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
            inertia_val = _gprop_inertia_to_numpy(sprops)

    if has_solid:
        vprops = GProp_GProps()
        BRepGProp.VolumeProperties_s(compound, vprops)
        volume_val = float(vprops.Mass())
        if volume_val > 0.0:
            com = vprops.CentreOfMass()
            centroid_val = np.array([com.X(), com.Y(), com.Z()], dtype=np.float64)
            inertia_val = _gprop_inertia_to_numpy(vprops)

    pts = _sample_shape_points(compound)
    return MeasureResult(
        bbox_result=bbox_res,
        length_val=length_val,
        area_val=area_val,
        volume_val=volume_val,
        centroid_val=centroid_val,
        inertia_val=inertia_val,
        sample_points=pts,
        ocp_shape=compound,
    )


def _hex8_to_tets(conn: np.ndarray) -> np.ndarray:
    """Decompose ``(M, 8)`` hex8 connectivity into ``(6*M, 4)`` tetrahedra."""
    # Standard 6-tet decomposition along diagonal (0, 6)
    tet_idx = np.array(
        [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ],
        dtype=np.int64,
    )
    return conn[:, tet_idx].reshape(-1, 4)


def _wedge6_to_tets(conn: np.ndarray) -> np.ndarray:
    """Decompose ``(M, 6)`` wedge6 connectivity into ``(3*M, 4)`` tetrahedra."""
    tet_idx = np.array(
        [
            [0, 1, 2, 3],
            [1, 2, 4, 3],
            [2, 4, 5, 3],
        ],
        dtype=np.int64,
    )
    return conn[:, tet_idx].reshape(-1, 4)


def _tet_volume_centroid_inertia(nodes: np.ndarray, tets: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute exact volume, centroid, and inertia tensor for a tetrahedral mesh."""
    pts = nodes[tets]  # (T, 4, 3)
    v0 = pts[:, 0, :]
    e1 = pts[:, 1, :] - v0
    e2 = pts[:, 2, :] - v0
    e3 = pts[:, 3, :] - v0
    vols = np.abs(np.einsum("ij,ij->i", e1, np.cross(e2, e3))) / 6.0
    total_vol = float(vols.sum())
    if total_vol <= 0.0:
        c = pts.reshape(-1, 3).mean(axis=0)
        return 0.0, c, np.zeros((3, 3), dtype=np.float64)

    tet_centroids = pts.mean(axis=1)  # (T, 3)
    global_c = (vols[:, None] * tet_centroids).sum(axis=0) / total_vol

    # Covariance about each tet centroid: C_tet = (V / 20) * sum_{k=0..3} p_k p_k^T
    rel = pts - tet_centroids[:, None, :]  # (T, 4, 3)
    cov_local = np.einsum("t,tki,tkj->ij", vols / 20.0, rel, rel)
    # Parallel axis contribution from tet centroid to global centroid
    d = tet_centroids - global_c[None, :]  # (T, 3)
    cov_shift = np.einsum("t,ti,tj->ij", vols, d, d)
    cov_total = cov_local + cov_shift
    inertia = np.trace(cov_total) * np.eye(3, dtype=np.float64) - cov_total
    return total_vol, global_c, inertia


def _tri_area_centroid_inertia(nodes: np.ndarray, tris: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute exact area, centroid, and surface inertia tensor for a triangle mesh."""
    pts = nodes[tris]  # (T, 3, 3)
    cross = np.cross(pts[:, 1, :] - pts[:, 0, :], pts[:, 2, :] - pts[:, 0, :])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    total_area = float(areas.sum())
    if total_area <= 0.0:
        c = pts.reshape(-1, 3).mean(axis=0)
        return 0.0, c, np.zeros((3, 3), dtype=np.float64)

    tri_centroids = pts.mean(axis=1)
    global_c = (areas[:, None] * tri_centroids).sum(axis=0) / total_area

    rel = pts - tri_centroids[:, None, :]
    cov_local = np.einsum("t,tki,tkj->ij", areas / 12.0, rel, rel)
    d = tri_centroids - global_c[None, :]
    cov_shift = np.einsum("t,ti,tj->ij", areas, d, d)
    cov_total = cov_local + cov_shift
    inertia = np.trace(cov_total) * np.eye(3, dtype=np.float64) - cov_total
    return total_area, global_c, inertia


def measure_mesh(
    nodes: np.ndarray,
    *,
    elements: np.ndarray | None = None,
    elem_type: str | None = None,
    tris: np.ndarray | None = None,
    edges: np.ndarray | None = None,
    node_indices: np.ndarray | None = None,
) -> MeasureResult:
    """Measure geometric properties of a mesh or mesh entity subset.

    Parameters
    ----------
    nodes:
        Node coordinate array of shape ``(N, 3)``.
    elements:
        Optional element connectivity array ``(M, K)`` for solid elements.
    elem_type:
        Canonical element type (``"hex8"``, ``"tet4"``, ``"tet10"``, ``"hex20"``,
        ``"wedge6"``). Required when ``elements`` is provided.
    tris:
        Optional surface triangle connectivity ``(T, 3)``.
    edges:
        Optional 1D edge connectivity ``(E, 2)``.
    node_indices:
        Optional 1D array of active node indices when measuring a node subset.
    """
    nodes_arr = np.asarray(nodes, dtype=np.float64).reshape(-1, 3)

    if elements is not None and len(elements) > 0:
        conn = np.asarray(elements, dtype=np.int64)
        used_nodes = np.unique(conn.ravel())
        pts = nodes_arr[used_nodes]
        bbox_res = BBoxResult.from_bounds(pts.min(axis=0), pts.max(axis=0))

        etype = (elem_type or "").lower()
        if etype in ("tet4", "tet10"):
            tets = conn[:, :4]
        elif etype in ("hex8", "hex20"):
            tets = _hex8_to_tets(conn[:, :8])
        elif etype == "wedge6":
            tets = _wedge6_to_tets(conn[:, :6])
        else:
            raise ValueError(f"Unsupported solid elem_type for measure_mesh: {elem_type!r}")

        vol, c, inertia = _tet_volume_centroid_inertia(nodes_arr, tets)
        area = 0.0
        if tris is not None and len(tris) > 0:
            area, _, _ = _tri_area_centroid_inertia(nodes_arr, np.asarray(tris, dtype=np.int64))
        return MeasureResult(
            bbox_result=bbox_res,
            area_val=area,
            volume_val=vol,
            centroid_val=c,
            inertia_val=inertia,
            sample_points=pts,
        )

    if tris is not None and len(tris) > 0:
        tris_arr = np.asarray(tris, dtype=np.int64).reshape(-1, 3)
        used_nodes = np.unique(tris_arr.ravel())
        pts = nodes_arr[used_nodes]
        bbox_res = BBoxResult.from_bounds(pts.min(axis=0), pts.max(axis=0))
        area, c, inertia = _tri_area_centroid_inertia(nodes_arr, tris_arr)
        return MeasureResult(
            bbox_result=bbox_res,
            area_val=area,
            centroid_val=c,
            inertia_val=inertia,
            sample_points=pts,
        )

    if edges is not None and len(edges) > 0:
        edges_arr = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
        used_nodes = np.unique(edges_arr.ravel())
        pts = nodes_arr[used_nodes]
        bbox_res = BBoxResult.from_bounds(pts.min(axis=0), pts.max(axis=0))
        seg_vecs = nodes_arr[edges_arr[:, 1]] - nodes_arr[edges_arr[:, 0]]
        seg_lens = np.linalg.norm(seg_vecs, axis=1)
        total_len = float(seg_lens.sum())
        mids = 0.5 * (nodes_arr[edges_arr[:, 0]] + nodes_arr[edges_arr[:, 1]])
        c = (seg_lens[:, None] * mids).sum(axis=0) / total_len if total_len > 0 else pts.mean(axis=0)
        return MeasureResult(
            bbox_result=bbox_res,
            length_val=total_len,
            centroid_val=c,
            sample_points=pts,
        )

    idx = np.asarray(node_indices, dtype=np.int64) if node_indices is not None else np.arange(len(nodes_arr))
    pts = nodes_arr[idx]
    if len(pts) == 0:
        raise ValueError("measure_mesh: empty selection")
    bbox_res = BBoxResult.from_bounds(pts.min(axis=0), pts.max(axis=0))
    return MeasureResult(
        bbox_result=bbox_res,
        centroid_val=pts.mean(axis=0),
        sample_points=pts,
    )
