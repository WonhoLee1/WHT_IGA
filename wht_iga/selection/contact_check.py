"""Pre-solve contact gap, penetration, and normal alignment inspection (v2.1 §8.5).

Uses ``scipy.spatial.cKDTree`` and exact 3D point-to-triangle projection to
compute :class:`ContactGapReport` between a master surface and a slave surface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class ContactGapReport:
    """Summary of geometric contact state between master and slave surfaces.

    Parameters
    ----------
    max_penetration:
        Maximum penetration depth (``>= 0.0``). Positive when slave nodes/points
        lie behind the master outward normal.
    min_clearance:
        Minimum positive clearance (``>= 0.0``). Equals ``0.0`` when any
        penetration or exact contact occurs.
    max_gap:
        Maximum signed gap across sampled slave points.
    mean_gap:
        Mean signed gap across sampled slave points.
    min_signed_gap:
        Minimum signed gap (negative indicates penetration).
    normal_opposition_deg:
        Mean angle in degrees between ``master_normal`` and ``-slave_normal``
        (``0.0`` indicates perfectly opposing normals).
    contact_area_ratio:
        Fraction of slave surface area within the search radius / active contact zone.
    suggested_search_radius:
        Recommended contact search radius for solver setup.
    n_penetrating_points:
        Number of sampled slave points with negative signed gap.
    n_sample_points:
        Total number of sampled slave points evaluated.
    """

    max_penetration: float
    min_clearance: float
    max_gap: float
    mean_gap: float
    min_signed_gap: float
    normal_opposition_deg: float
    contact_area_ratio: float
    suggested_search_radius: float
    n_penetrating_points: int
    n_sample_points: int


def _closest_point_on_triangle(
    p: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
) -> np.ndarray:
    """Return the exact closest point on triangle ``(a, b, c)`` to point ``p`` (Ericson)."""
    ab = b - a
    ac = c - a
    ap = p - a

    d1 = float(np.dot(ab, ap))
    d2 = float(np.dot(ac, ap))
    if d1 <= 0.0 and d2 <= 0.0:
        return a

    bp = p - b
    d3 = float(np.dot(ab, bp))
    d4 = float(np.dot(ac, bp))
    if d3 >= 0.0 and d4 <= d3:
        return b

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3)
        return a + v * ab

    cp = p - c
    d5 = float(np.dot(ab, cp))
    d6 = float(np.dot(ac, cp))
    if d6 >= 0.0 and d5 <= d6:
        return c

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6)
        return a + w * ac

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b)

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    return a + ab * v + ac * w


def _tri_normals_and_areas(nodes: np.ndarray, tris: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute unit normals ``(T, 3)`` and areas ``(T,)`` for triangles."""
    v0 = nodes[tris[:, 0]]
    v1 = nodes[tris[:, 1]]
    v2 = nodes[tris[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)
    norms = np.linalg.norm(cross, axis=1)
    areas = 0.5 * norms
    safe = np.where(norms > 1e-15, norms, 1.0)
    unit_normals = cross / safe[:, None]
    return unit_normals, areas


def _coerce_surface_mesh(
    surface: Any,
    tess_deflection: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a surface specification into ``(nodes, tris)`` arrays."""
    if isinstance(surface, tuple) and len(surface) == 2:
        nodes = np.asarray(surface[0], dtype=np.float64).reshape(-1, 3)
        tris = np.asarray(surface[1], dtype=np.int64).reshape(-1, 3)
        return nodes, tris

    if hasattr(surface, "to_surface_tris") and callable(surface.to_surface_tris):
        return surface.to_surface_tris(deflection=tess_deflection)

    if hasattr(surface, "nodes") and hasattr(surface, "tris"):
        nodes = np.asarray(surface.nodes, dtype=np.float64).reshape(-1, 3)
        tris = np.asarray(surface.tris, dtype=np.int64).reshape(-1, 3)
        return nodes, tris

    # Fallback: assume OCP TopoDS_Shape
    from OCP.Bnd import Bnd_Box
    from OCP.BRep import BRep_Tool
    from OCP.BRepBndLib import BRepBndLib
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopAbs import TopAbs_FACE, TopAbs_Orientation
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.TopLoc import TopLoc_Location
    from OCP.gp import gp_Pnt

    if tess_deflection is None:
        box = Bnd_Box()
        BRepBndLib.Add_s(surface, box, True)
        diag = float(
            np.linalg.norm(
                [
                    box.GetXMax() - box.GetXMin(),
                    box.GetYMax() - box.GetYMin(),
                    box.GetZMax() - box.GetZMin(),
                ]
            )
        )
        tess_deflection = max(1e-3 * diag, 1e-4)

    BRepMesh_IncrementalMesh(surface, tess_deflection, False, 0.3, True)
    nodes_list: list[tuple[float, float, float]] = []
    tris_list: list[tuple[int, int, int]] = []
    reversed_flag = TopAbs_Orientation.TopAbs_REVERSED

    exp = TopExp_Explorer(surface, TopAbs_FACE)
    while exp.More():
        face = TopoDS.Face(exp.Current())
        flip = face.Orientation() == reversed_flag
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is not None:
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
        exp.Next()

    if not tris_list:
        raise ValueError("Surface shape produced no triangles")
    return np.asarray(nodes_list, dtype=np.float64), np.asarray(tris_list, dtype=np.int64)


def inspect_contact(
    master: Any,
    slave: Any,
    *,
    search_radius: float | None = None,
    tess_deflection: float | None = None,
    k_candidates: int = 16,
) -> ContactGapReport:
    """Inspect geometric contact gap, penetration, and normal opposition.

    Parameters
    ----------
    master:
        Master surface as ``(nodes, tris)``, an :class:`EntityQuery`, a
        ``Tessellation``, or an OCP ``TopoDS_Shape``.
    slave:
        Slave surface in any of the same formats as ``master``.
    search_radius:
        Optional contact search radius threshold for computing
        ``contact_area_ratio``. If ``None``, a radius is suggested automatically.
    tess_deflection:
        Optional linear deflection used when tessellating CAD shapes.
    k_candidates:
        Number of nearest master triangles queried per slave sample point.
    """
    m_nodes, m_tris = _coerce_surface_mesh(master, tess_deflection=tess_deflection)
    s_nodes, s_tris = _coerce_surface_mesh(slave, tess_deflection=tess_deflection)

    if len(m_tris) == 0 or len(s_tris) == 0:
        raise ValueError("inspect_contact requires non-empty master and slave triangulations")

    m_normals, _ = _tri_normals_and_areas(m_nodes, m_tris)
    s_normals, s_areas = _tri_normals_and_areas(s_nodes, s_tris)

    m_centroids = m_nodes[m_tris].mean(axis=1)
    s_centroids = s_nodes[s_tris].mean(axis=1)

    # Evaluate at slave triangle centroids (area-weighted) and unique slave vertices
    s_vert_ids = np.unique(s_tris.ravel())
    # Compute vertex normals by area-weighted averaging of incident triangle normals
    s_vnormals = np.zeros((len(s_nodes), 3), dtype=np.float64)
    for idx_corner in range(3):
        np.add.at(s_vnormals, s_tris[:, idx_corner], s_normals * s_areas[:, None])
    vn_lens = np.linalg.norm(s_vnormals, axis=1)
    vn_safe = np.where(vn_lens > 1e-15, vn_lens, 1.0)
    s_vnormals = s_vnormals / vn_safe[:, None]

    sample_pts = np.vstack([s_centroids, s_nodes[s_vert_ids]])
    sample_normals = np.vstack([s_normals, s_vnormals[s_vert_ids]])

    tree = cKDTree(m_centroids)
    k = min(int(k_candidates), len(m_tris))
    _, cand_indices = tree.query(sample_pts, k=k)
    if k == 1:
        cand_indices = cand_indices[:, None]

    n_samples = len(sample_pts)
    signed_gaps = np.empty(n_samples, dtype=np.float64)
    opp_angles_deg = np.empty(n_samples, dtype=np.float64)

    for i in range(n_samples):
        p = sample_pts[i]
        ns = sample_normals[i]
        cands = cand_indices[i]

        # Prefer opposing-normal master triangles when available
        dots = m_normals[cands] @ ns
        opposing_mask = dots < 0.0
        pool = cands[opposing_mask] if np.any(opposing_mask) else cands

        best_dist2 = float("inf")
        best_q = p
        best_nm = m_normals[pool[0]]
        for tri_idx in pool:
            t = m_tris[tri_idx]
            q = _closest_point_on_triangle(p, m_nodes[t[0]], m_nodes[t[1]], m_nodes[t[2]])
            d2 = float(np.dot(p - q, p - q))
            if d2 < best_dist2:
                best_dist2 = d2
                best_q = q
                best_nm = m_normals[tri_idx]

        signed_gaps[i] = float(np.dot(p - best_q, best_nm))
        cos_opp = float(np.clip(np.dot(best_nm, -ns), -1.0, 1.0))
        opp_angles_deg[i] = float(np.degrees(np.arccos(cos_opp)))

    min_signed = float(np.min(signed_gaps))
    max_pen = float(max(0.0, -min_signed))
    min_clr = float(max(0.0, min_signed))
    max_gap = float(np.max(signed_gaps))
    mean_gap = float(np.mean(signed_gaps))
    n_pen = int(np.sum(signed_gaps < -1e-12))
    mean_opp_deg = float(np.mean(opp_angles_deg))

    # Characteristic master edge length for suggested search radius
    e01 = np.linalg.norm(m_nodes[m_tris[:, 1]] - m_nodes[m_tris[:, 0]], axis=1)
    char_len = float(np.median(e01)) if len(e01) > 0 else 1.0
    suggested = float(max(2.0 * max(max_pen, min_clr), 0.25 * char_len, 1e-6))

    radius_tol = float(search_radius) if search_radius is not None else suggested
    centroid_gaps = signed_gaps[: len(s_tris)]
    total_s_area = float(s_areas.sum())
    if total_s_area > 0.0:
        active_mask = centroid_gaps <= radius_tol
        area_ratio = float(s_areas[active_mask].sum() / total_s_area)
    else:
        area_ratio = 0.0

    return ContactGapReport(
        max_penetration=max_pen,
        min_clearance=min_clr,
        max_gap=max_gap,
        mean_gap=mean_gap,
        min_signed_gap=min_signed,
        normal_opposition_deg=mean_opp_deg,
        contact_area_ratio=area_ratio,
        suggested_search_radius=suggested,
        n_penetrating_points=n_pen,
        n_sample_points=n_samples,
    )
