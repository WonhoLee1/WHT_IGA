# -*- coding: utf-8 -*-
"""3D curved triangulated surface contact projector, NLGEOM/J2 body operator,
and multi-body SSN contact coupling.

Extends :mod:`wht_iga.analyses.contact` and :mod:`wht_iga.analyses.contact_nlgeom`
without modifying them:

* :class:`DeformableTriangleSurfaceMaster` — drop-in replacement for
  :class:`~wht_iga.analyses.contact.DeformableBoxFaceMaster` that projects slave
  Gauss points onto an arbitrary deformed triangulated B-Rep surface
  ``(V_def = tri_nodes + B_verts @ U_m, tris)`` via Numba Voronoi closest-point
  projection with AABB bounding-box pruning, crease-preserving smooth corner
  normals, and exact sparse CSR basis interpolation
  ``B_m = b0 * B_v0 + b1 * B_v1 + b2 * B_v2``.
* :func:`filter_triangles` / :func:`subdivide_long_edges` — surface mesh helpers
  for selecting contact zones and refining straight cylinder generators without
  moving curved boundaries.
* :class:`NonlinearBodyOperator` — multi-step Kratos ``TotalLagrangianElement3D4N``
  bulk operator supporting single or multi-region trimmed solids,
  ``HyperElastic3DLaw``, and ``FiniteStrainIsotropicPlasticity3DVonMisesVonMises``
  2-point / piecewise-linear isotropic hardening (``HARDENING_CURVE = 6``) with
  pure trial evaluations (``eval``) and step-end commit (``commit_step``).
* :func:`solve_multibody_ssn_contact` and
  :func:`solve_multibody_nlgeom_contact_step` — coupled N-body, P-pair frictional
  SSN contact solver and outer NLGEOM/plasticity Newton loop.
"""

from __future__ import annotations

import gc
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
from numba import njit

from wht_iga.analyses.contact import (
    _solve_linear_step,
    bspline_basis_matrix,
    expand_vector_basis,
)
from wht_iga.analyses.thermal_plasticity import (
    FLAT_EXTENSION,
    HardeningTable,
    _check_point_curve,
    _curve_work,
    characteristic_length,
)
from wht_iga.kratos_backend.elements import check_qp_element_name, equation_permutation
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume


def release_solver_memory() -> None:
    """Release cached MKL Pardiso factorizations and collect Python/Kratos objects."""
    try:
        from pypardiso.scipy_aliases import pypardiso_solver

        pypardiso_solver.free_memory(everything=True)
    except Exception:
        pass
    gc.collect()


# --------------------------------------------------------------------------------------
# Surface triangle mesh helpers
# --------------------------------------------------------------------------------------


def filter_triangles(
    nodes: np.ndarray,
    tris: np.ndarray,
    keep_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a triangle subset and compact its vertex array.

    Parameters
    ----------
    nodes:
        Vertex coordinates ``(N, 3)``.
    tris:
        Triangle connectivity ``(K, 3)``.
    keep_mask:
        Boolean mask of length ``K`` selecting triangles to keep.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Compacted ``(sub_nodes, sub_tris)`` with contiguous 0-based vertex indices.
    """
    pts = np.asarray(nodes, dtype=np.float64)
    tri_arr = np.asarray(tris, dtype=np.int64)
    mask = np.asarray(keep_mask, dtype=bool)
    sub_tris = tri_arr[mask]
    if sub_tris.shape[0] == 0:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64)
    used, inv = np.unique(sub_tris.ravel(), return_inverse=True)
    return np.ascontiguousarray(pts[used], dtype=np.float64), np.ascontiguousarray(
        inv.reshape(-1, 3), dtype=np.int64
    )


def subdivide_long_edges(
    nodes: np.ndarray,
    tris: np.ndarray,
    max_edge: float,
    max_passes: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Conformingly subdivide triangles whose edges exceed ``max_edge`` (preserves winding).

    Uses edge-shared Red-Green (1->2 / 1->3 / 1->4) subdivision so every split
    edge is bisected identically by both incident triangles. Watertight manifold
    meshes remain strictly conforming and watertight (zero T-junction cracks),
    preserving continuous vertex normals and ``exclude_boundary=True`` topology.
    """
    pts = np.asarray(nodes, dtype=np.float64)
    tri_arr = np.asarray(tris, dtype=np.int64)
    if tri_arr.shape[0] == 0 or float(max_edge) <= 0.0:
        return pts.copy(), tri_arr.copy()

    limit_sq = float(max_edge) ** 2
    for _ in range(int(max_passes)):
        p = pts[tri_arr]  # (K, 3, 3)
        e0 = np.sum((p[:, 1] - p[:, 0]) ** 2, axis=1)
        e1 = np.sum((p[:, 2] - p[:, 1]) ** 2, axis=1)
        e2 = np.sum((p[:, 0] - p[:, 2]) ** 2, axis=1)

        m0 = e0 > limit_sq
        m1 = e1 > limit_sq
        m2 = e2 > limit_sq
        if not (np.any(m0) or np.any(m1) or np.any(m2)):
            break

        # Gather every undirected edge exceeding max_edge across all triangles
        parts = []
        if np.any(m0):
            parts.append(tri_arr[m0][:, [0, 1]])
        if np.any(m1):
            parts.append(tri_arr[m1][:, [1, 2]])
        if np.any(m2):
            parts.append(tri_arr[m2][:, [2, 0]])
        raw_edges = np.vstack(parts)
        canon_edges = np.sort(raw_edges, axis=1)
        uniq_edges = np.unique(canon_edges, axis=0)

        n0 = pts.shape[0]
        mids = 0.5 * (pts[uniq_edges[:, 0]] + pts[uniq_edges[:, 1]])
        pts = np.vstack([pts, mids])

        edge_mid: dict[tuple[int, int], int] = {
            (int(u), int(v)): n0 + idx for idx, (u, v) in enumerate(uniq_edges)
        }

        out_tris: list[tuple[int, int, int]] = []
        for row in tri_arr:
            v0, v1, v2 = int(row[0]), int(row[1]), int(row[2])
            k01 = (v0, v1) if v0 < v1 else (v1, v0)
            k12 = (v1, v2) if v1 < v2 else (v2, v1)
            k20 = (v2, v0) if v2 < v0 else (v0, v2)
            s01 = edge_mid.get(k01, -1)
            s12 = edge_mid.get(k12, -1)
            s20 = edge_mid.get(k20, -1)

            n_split = (s01 >= 0) + (s12 >= 0) + (s20 >= 0)
            if n_split == 0:
                out_tris.append((v0, v1, v2))
            elif n_split == 3:
                out_tris.append((v0, s01, s20))
                out_tris.append((s01, v1, s12))
                out_tris.append((s20, s12, v2))
                out_tris.append((s01, s12, s20))
            elif n_split == 1:
                if s01 >= 0:
                    a, b, c, m = v0, v1, v2, s01
                elif s12 >= 0:
                    a, b, c, m = v1, v2, v0, s12
                else:
                    a, b, c, m = v2, v0, v1, s20
                out_tris.append((a, m, c))
                out_tris.append((m, b, c))
            else:  # n_split == 2
                # Rotate so edges (a, b) and (b, c) are split, (c, a) is unsplit
                if s20 < 0:
                    a, b, c, mab, mbc = v0, v1, v2, s01, s12
                elif s01 < 0:
                    a, b, c, mab, mbc = v1, v2, v0, s12, s20
                else:
                    a, b, c, mab, mbc = v2, v0, v1, s20, s01
                out_tris.append((mab, b, mbc))
                # Split quad (a, mab, mbc, c) along its shorter diagonal
                d_a_mbc = float(np.sum((pts[a] - pts[mbc]) ** 2))
                d_mab_c = float(np.sum((pts[mab] - pts[c]) ** 2))
                if d_a_mbc <= d_mab_c:
                    out_tris.append((a, mab, mbc))
                    out_tris.append((a, mbc, c))
                else:
                    out_tris.append((a, mab, c))
                    out_tris.append((mab, mbc, c))

        tri_arr = np.asarray(out_tris, dtype=np.int64)

    return np.ascontiguousarray(pts, dtype=np.float64), np.ascontiguousarray(
        tri_arr, dtype=np.int64
    )


# --------------------------------------------------------------------------------------
# Numba closest-point & smooth corner normal kernels
# --------------------------------------------------------------------------------------


@njit(cache=True)
def _compute_corner_normals_njit(
    V: np.ndarray,
    tris: np.ndarray,
    v_tri_indptr: np.ndarray,
    v_tri_ids: np.ndarray,
    v_corner_idx: np.ndarray,
    outward_flip: bool,
    smooth_normals: bool,
    crease_cos: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute unit face normals and crease-preserving angle-weighted corner normals."""
    K = tris.shape[0]
    fn = np.empty((K, 3), dtype=np.float64)
    angles = np.empty((K, 3), dtype=np.float64)
    N_corner = np.empty((K, 3, 3), dtype=np.float64)
    sign = -1.0 if outward_flip else 1.0

    for k in range(K):
        i0, i1, i2 = tris[k, 0], tris[k, 1], tris[k, 2]
        p0x, p0y, p0z = V[i0, 0], V[i0, 1], V[i0, 2]
        p1x, p1y, p1z = V[i1, 0], V[i1, 1], V[i1, 2]
        p2x, p2y, p2z = V[i2, 0], V[i2, 1], V[i2, 2]

        e01x, e01y, e01z = p1x - p0x, p1y - p0y, p1z - p0z
        e02x, e02y, e02z = p2x - p0x, p2y - p0y, p2z - p0z
        e12x, e12y, e12z = p2x - p1x, p2y - p1y, p2z - p1z

        cx = sign * (e01y * e02z - e01z * e02y)
        cy = sign * (e01z * e02x - e01x * e02z)
        cz = sign * (e01x * e02y - e01y * e02x)
        clen = np.sqrt(cx * cx + cy * cy + cz * cz)
        if clen > 1e-30:
            fn[k, 0] = cx / clen
            fn[k, 1] = cy / clen
            fn[k, 2] = cz / clen
        else:
            fn[k, 0] = 0.0
            fn[k, 1] = 0.0
            fn[k, 2] = 1.0

        l01 = np.sqrt(e01x * e01x + e01y * e01y + e01z * e01z)
        l02 = np.sqrt(e02x * e02x + e02y * e02y + e02z * e02z)
        l12 = np.sqrt(e12x * e12x + e12y * e12y + e12z * e12z)

        if l01 > 1e-30 and l02 > 1e-30:
            d0 = (e01x * e02x + e01y * e02y + e01z * e02z) / (l01 * l02)
            angles[k, 0] = np.arccos(min(1.0, max(-1.0, d0)))
        else:
            angles[k, 0] = 0.0

        if l01 > 1e-30 and l12 > 1e-30:
            d1 = -(e01x * e12x + e01y * e12y + e01z * e12z) / (l01 * l12)
            angles[k, 1] = np.arccos(min(1.0, max(-1.0, d1)))
        else:
            angles[k, 1] = 0.0

        if l02 > 1e-30 and l12 > 1e-30:
            d2 = (e02x * e12x + e02y * e12y + e02z * e12z) / (l02 * l12)
            angles[k, 2] = np.arccos(min(1.0, max(-1.0, d2)))
        else:
            angles[k, 2] = 0.0

    for k in range(K):
        fnx, fny, fnz = fn[k, 0], fn[k, 1], fn[k, 2]
        if not smooth_normals:
            for c in range(3):
                N_corner[k, c, 0] = fnx
                N_corner[k, c, 1] = fny
                N_corner[k, c, 2] = fnz
            continue

        for c in range(3):
            vid = tris[k, c]
            s0 = 0.0
            s1 = 0.0
            s2 = 0.0
            for idx in range(v_tri_indptr[vid], v_tri_indptr[vid + 1]):
                k2 = v_tri_ids[idx]
                c2 = v_corner_idx[idx]
                dot_f = fnx * fn[k2, 0] + fny * fn[k2, 1] + fnz * fn[k2, 2]
                if dot_f >= crease_cos:
                    w = angles[k2, c2]
                    s0 += w * fn[k2, 0]
                    s1 += w * fn[k2, 1]
                    s2 += w * fn[k2, 2]
            slen = np.sqrt(s0 * s0 + s1 * s1 + s2 * s2)
            if slen > 1e-30:
                N_corner[k, c, 0] = s0 / slen
                N_corner[k, c, 1] = s1 / slen
                N_corner[k, c, 2] = s2 / slen
            else:
                N_corner[k, c, 0] = fnx
                N_corner[k, c, 1] = fny
                N_corner[k, c, 2] = fnz

    return fn, N_corner


@njit(cache=True)
def _closest_point_triangle_njit(
    px: float,
    py: float,
    pz: float,
    ax: float,
    ay: float,
    az: float,
    bx: float,
    by: float,
    bz: float,
    cx: float,
    cy: float,
    cz: float,
) -> tuple[float, float, float, float, float, float, float]:
    """Exact Voronoi closest point ``Q`` and barycentric coordinates ``(b0, b1, b2)``."""
    abx, aby, abz = bx - ax, by - ay, bz - az
    acx, acy, acz = cx - ax, cy - ay, cz - az
    apx, apy, apz = px - ax, py - ay, pz - az

    d1 = abx * apx + aby * apy + abz * apz
    d2 = acx * apx + acy * apy + acz * apz
    if d1 <= 0.0 and d2 <= 0.0:
        dx, dy, dz = px - ax, py - ay, pz - az
        return ax, ay, az, 1.0, 0.0, 0.0, dx * dx + dy * dy + dz * dz

    bpx, bpy, bpz = px - bx, py - by, pz - bz
    d3 = abx * bpx + aby * bpy + abz * bpz
    d4 = acx * bpx + acy * bpy + acz * bpz
    if d3 >= 0.0 and d4 <= d3:
        dx, dy, dz = px - bx, py - by, pz - bz
        return bx, by, bz, 0.0, 1.0, 0.0, dx * dx + dy * dy + dz * dz

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        denom = d1 - d3
        v = d1 / denom if denom > 1e-30 else 0.0
        qx = ax + v * abx
        qy = ay + v * aby
        qz = az + v * abz
        dx, dy, dz = px - qx, py - qy, pz - qz
        return qx, qy, qz, 1.0 - v, v, 0.0, dx * dx + dy * dy + dz * dz

    cpx, cpy, cpz = px - cx, py - cy, pz - cz
    d5 = abx * cpx + aby * cpy + abz * cpz
    d6 = acx * cpx + acy * cpy + acz * cpz
    if d6 >= 0.0 and d5 <= d6:
        dx, dy, dz = px - cx, py - cy, pz - cz
        return cx, cy, cz, 0.0, 0.0, 1.0, dx * dx + dy * dy + dz * dz

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        denom = d2 - d6
        w = d2 / denom if denom > 1e-30 else 0.0
        qx = ax + w * acx
        qy = ay + w * acy
        qz = az + w * acz
        dx, dy, dz = px - qx, py - qy, pz - qz
        return qx, qy, qz, 1.0 - w, 0.0, w, dx * dx + dy * dy + dz * dz

    va = d3 * d6 - d5 * d4
    d4_d3 = d4 - d3
    d5_d6 = d5 - d6
    if va <= 0.0 and d4_d3 >= 0.0 and d5_d6 >= 0.0:
        denom = d4_d3 + d5_d6
        w = d4_d3 / denom if denom > 1e-30 else 0.0
        qx = bx + w * (cx - bx)
        qy = by + w * (cy - by)
        qz = bz + w * (cz - bz)
        dx, dy, dz = px - qx, py - qy, pz - qz
        return qx, qy, qz, 0.0, 1.0 - w, w, dx * dx + dy * dy + dz * dz

    s = va + vb + vc
    if s > 1e-30:
        inv_s = 1.0 / s
        v = vb * inv_s
        w = vc * inv_s
        u = 1.0 - v - w
        qx = ax + v * abx + w * acx
        qy = ay + v * aby + w * acy
        qz = az + v * abz + w * acz
        dx, dy, dz = px - qx, py - qy, pz - qz
        return qx, qy, qz, u, v, w, dx * dx + dy * dy + dz * dz

    # Degenerate zero-area triangle fallback: pick nearest vertex
    d_a = apx * apx + apy * apy + apz * apz
    d_b = bpx * bpx + bpy * bpy + bpz * bpz
    d_c = cpx * cpx + cpy * cpy + cpz * cpz
    if d_a <= d_b and d_a <= d_c:
        return ax, ay, az, 1.0, 0.0, 0.0, d_a
    if d_b <= d_c:
        return bx, by, bz, 0.0, 1.0, 0.0, d_b
    return cx, cy, cz, 0.0, 0.0, 1.0, d_c


@njit(cache=True)
def _project_surface_triangles_njit(
    x_s: np.ndarray,
    V: np.ndarray,
    tris: np.ndarray,
    fn: np.ndarray,
    N_corner: np.ndarray,
    hint_tri: np.ndarray,
    tri_edge_is_boundary: np.ndarray,
    vert_is_boundary: np.ndarray,
    exclude_boundary: bool,
    boundary_cos_tol: float,
    max_penetration: float,
    max_distance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Exact AABB-pruned closest-point projection onto ``(V, tris)``."""
    K = tris.shape[0]
    M = x_s.shape[0]

    tri_min = np.empty((K, 3), dtype=np.float64)
    tri_max = np.empty((K, 3), dtype=np.float64)
    tri_cx = np.empty(K, dtype=np.float64)
    tri_cy = np.empty(K, dtype=np.float64)
    tri_cz = np.empty(K, dtype=np.float64)

    for k in range(K):
        i0, i1, i2 = tris[k, 0], tris[k, 1], tris[k, 2]
        for d in range(3):
            v0, v1, v2 = V[i0, d], V[i1, d], V[i2, d]
            lo = v0
            if v1 < lo:
                lo = v1
            if v2 < lo:
                lo = v2
            hi = v0
            if v1 > hi:
                hi = v1
            if v2 > hi:
                hi = v2
            tri_min[k, d] = lo
            tri_max[k, d] = hi
        tri_cx[k] = (V[i0, 0] + V[i1, 0] + V[i2, 0]) / 3.0
        tri_cy[k] = (V[i0, 1] + V[i1, 1] + V[i2, 1]) / 3.0
        tri_cz[k] = (V[i0, 2] + V[i1, 2] + V[i2, 2]) / 3.0

    tri_idx = np.empty(M, dtype=np.int64)
    bary = np.empty((M, 3), dtype=np.float64)
    x_m = np.empty((M, 3), dtype=np.float64)
    n_m = np.empty((M, 3), dtype=np.float64)
    g_n = np.empty(M, dtype=np.float64)

    for q in range(M):
        px, py, pz = x_s[q, 0], x_s[q, 1], x_s[q, 2]

        seed = hint_tri[q]
        if seed < 0 or seed >= K:
            best_c2 = 1e300
            seed = 0
            for k in range(K):
                dx = px - tri_cx[k]
                dy = py - tri_cy[k]
                dz = pz - tri_cz[k]
                c2 = dx * dx + dy * dy + dz * dz
                if c2 < best_c2:
                    best_c2 = c2
                    seed = k

        i0, i1, i2 = tris[seed, 0], tris[seed, 1], tris[seed, 2]
        best_qx, best_qy, best_qz, best_b0, best_b1, best_b2, best_d2 = _closest_point_triangle_njit(
            px, py, pz,
            V[i0, 0], V[i0, 1], V[i0, 2],
            V[i1, 0], V[i1, 1], V[i1, 2],
            V[i2, 0], V[i2, 1], V[i2, 2],
        )
        best_k = seed
        best_align = abs((px - best_qx) * fn[seed, 0] + (py - best_qy) * fn[seed, 1] + (pz - best_qz) * fn[seed, 2])

        for k in range(K):
            if k == seed:
                continue
            dx_box = 0.0
            if px < tri_min[k, 0]:
                dx_box = tri_min[k, 0] - px
            elif px > tri_max[k, 0]:
                dx_box = px - tri_max[k, 0]
            dy_box = 0.0
            if py < tri_min[k, 1]:
                dy_box = tri_min[k, 1] - py
            elif py > tri_max[k, 1]:
                dy_box = py - tri_max[k, 1]
            dz_box = 0.0
            if pz < tri_min[k, 2]:
                dz_box = tri_min[k, 2] - pz
            elif pz > tri_max[k, 2]:
                dz_box = pz - tri_max[k, 2]

            if dx_box * dx_box + dy_box * dy_box + dz_box * dz_box > best_d2 * (1.0 + 1e-12) + 1e-24:
                continue

            i0, i1, i2 = tris[k, 0], tris[k, 1], tris[k, 2]
            qx, qy, qz, b0, b1, b2, d2 = _closest_point_triangle_njit(
                px, py, pz,
                V[i0, 0], V[i0, 1], V[i0, 2],
                V[i1, 0], V[i1, 1], V[i1, 2],
                V[i2, 0], V[i2, 1], V[i2, 2],
            )
            tie_tol = 1e-22 * max(1.0, best_d2)
            if d2 < best_d2 - tie_tol:
                best_d2 = d2
                best_k = k
                best_qx, best_qy, best_qz = qx, qy, qz
                best_b0, best_b1, best_b2 = b0, b1, b2
                best_align = abs((px - qx) * fn[k, 0] + (py - qy) * fn[k, 1] + (pz - qz) * fn[k, 2])
            elif d2 <= best_d2 + tie_tol:
                align = abs((px - qx) * fn[k, 0] + (py - qy) * fn[k, 1] + (pz - qz) * fn[k, 2])
                if align > best_align:
                    best_d2 = d2
                    best_k = k
                    best_qx, best_qy, best_qz = qx, qy, qz
                    best_b0, best_b1, best_b2 = b0, b1, b2
                    best_align = align

        nx = best_b0 * N_corner[best_k, 0, 0] + best_b1 * N_corner[best_k, 1, 0] + best_b2 * N_corner[best_k, 2, 0]
        ny = best_b0 * N_corner[best_k, 0, 1] + best_b1 * N_corner[best_k, 1, 1] + best_b2 * N_corner[best_k, 2, 1]
        nz = best_b0 * N_corner[best_k, 0, 2] + best_b1 * N_corner[best_k, 1, 2] + best_b2 * N_corner[best_k, 2, 2]
        nlen = np.sqrt(nx * nx + ny * ny + nz * nz)
        if nlen > 1e-30:
            nx /= nlen
            ny /= nlen
            nz /= nlen
        else:
            nx, ny, nz = fn[best_k, 0], fn[best_k, 1], fn[best_k, 2]

        rx, ry, rz = px - best_qx, py - best_qy, pz - best_qz
        gn_val = rx * nx + ry * ny + rz * nz
        dist = np.sqrt(max(best_d2, 0.0))

        # Optional off-patch / out-of-range active candidate filtering
        if max_distance > 0.0 and dist > max_distance:
            gn_val = max(abs(gn_val), dist)
        elif max_penetration > 0.0 and gn_val < -max_penetration:
            gn_val = max(abs(gn_val), dist)
        elif exclude_boundary and dist > 1e-12:
            on_bnd = False
            if best_b0 <= 1e-10 and tri_edge_is_boundary[best_k, 1]:
                on_bnd = True
            elif best_b1 <= 1e-10 and tri_edge_is_boundary[best_k, 2]:
                on_bnd = True
            elif best_b2 <= 1e-10 and tri_edge_is_boundary[best_k, 0]:
                on_bnd = True
            elif best_b0 >= 1.0 - 1e-10 and vert_is_boundary[tris[best_k, 0]]:
                on_bnd = True
            elif best_b1 >= 1.0 - 1e-10 and vert_is_boundary[tris[best_k, 1]]:
                on_bnd = True
            elif best_b2 >= 1.0 - 1e-10 and vert_is_boundary[tris[best_k, 2]]:
                on_bnd = True
            if on_bnd and abs(gn_val) < boundary_cos_tol * dist:
                gn_val = max(abs(gn_val), dist)

        tri_idx[q] = best_k
        bary[q, 0] = best_b0
        bary[q, 1] = best_b1
        bary[q, 2] = best_b2
        x_m[q, 0] = best_qx
        x_m[q, 1] = best_qy
        x_m[q, 2] = best_qz
        n_m[q, 0] = nx
        n_m[q, 1] = ny
        n_m[q, 2] = nz
        g_n[q] = gn_val

    return tri_idx, bary, x_m, n_m, g_n


# --------------------------------------------------------------------------------------
# DeformableTriangleSurfaceMaster
# --------------------------------------------------------------------------------------


class DeformableTriangleSurfaceMaster:
    """Closest-point projector and basis evaluator on a deformed 3D triangle surface.

    Compatible with :class:`~wht_iga.analyses.contact.DeformableBoxFaceMaster`:
    exposes ``n_ctrl`` and ``evaluate(x_s, U_m) -> (g_n, n_m, x_m, ab_m, B_m)``,
    plugging directly into :func:`~wht_iga.analyses.contact.solve_deformable_ssn_contact`,
    :func:`~wht_iga.analyses.contact_nlgeom.solve_nlgeom_contact_step`, and
    :func:`solve_multibody_nlgeom_contact_step`.

    Normal & gap semantics:
    * Each triangle ``(v0, v1, v2)`` has outward unit normal
      ``sign * normalize((V1 - V0) x (V2 - V0))`` where ``sign = -1`` when
      ``outward_flip=True`` and ``+1`` otherwise (matching OCC ``tessellate``).
    * When ``smooth_normals=True`` (default), corner normals are angle-weighted
      averages across smooth neighbors (``dot(n_1, n_2) >= crease_cos``), giving
      ``C^0``-continuous normals across curved facets while preserving sharp
      creases (> 60 deg by default).
    * Signed gap is ``g_n = dot(x_s - x_m, n_m)`` (negative on penetration).
    * When ``exclude_boundary=True`` (or ``max_penetration`` / ``max_distance`` is
      set), slave points overhanging an open boundary edge (``|g_n| < boundary_cos_tol * dist``)
      or exceeding distance limits are safely masked out by setting ``g_n > 0``.
    * ``ab_m`` returns triangle barycentric coordinates ``(b1, b2)`` (with
      ``b0 = 1 - b1 - b2``), and ``B_m`` is the sparse CSR row interpolation
      ``b0 * B_verts[v0] + b1 * B_verts[v1] + b2 * B_verts[v2]``.

    Parameters
    ----------
    box:
        ``BoxSpline`` discretizing the master body.
    tri_nodes:
        Reference vertex coordinates ``(N_v, 3)``.
    tris:
        Triangle connectivity ``(K, 3)`` indexing ``tri_nodes``.
    outward_flip:
        Invert triangle winding normal if ``True``.
    smooth_normals:
        Interpolate crease-preserving angle-weighted vertex normals for smooth
        curved-surface SSN Newton convergence.
    crease_cos:
        Cosine threshold for smoothing across adjacent triangles (default ``0.5``
        for 60 deg).
    exclude_boundary:
        If ``True``, deactivate projections overhanging an open patch boundary edge.
    boundary_cos_tol:
        Minimum ``|dot(x_s - x_m, n_m)| / ||x_s - x_m||`` required at boundary
        edges/vertices when ``exclude_boundary=True``.
    max_penetration:
        Optional positive cap; projections with ``g_n < -max_penetration`` are
        marked inactive (prevents grabbing the back wall of thin shells).
    max_distance:
        Optional maximum projection distance ``||x_s - x_m||`` for candidate
        activation.
    """

    def __init__(
        self,
        box: Any,
        tri_nodes: np.ndarray,
        tris: np.ndarray,
        outward_flip: bool = False,
        *,
        smooth_normals: bool = True,
        crease_cos: float = 0.5,
        exclude_boundary: bool = False,
        boundary_cos_tol: float = 0.5,
        max_penetration: float | None = None,
        max_distance: float | None = None,
    ) -> None:
        pts = np.ascontiguousarray(tri_nodes, dtype=np.float64)
        tri_arr = np.ascontiguousarray(tris, dtype=np.int64)
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError(f"tri_nodes must have shape (N, 3), got {pts.shape}")
        if tri_arr.ndim != 2 or tri_arr.shape[1] != 3 or tri_arr.shape[0] == 0:
            raise ValueError(f"tris must have non-empty shape (K, 3), got {tri_arr.shape}")
        if tri_arr.min() < 0 or tri_arr.max() >= pts.shape[0]:
            raise ValueError("tris vertex index out of range for tri_nodes")

        self.box = box
        self.n_ctrl = int(box.n_ctrl)
        self.tri_nodes = pts
        self.tris = tri_arr
        self.n_verts = int(pts.shape[0])
        self.n_tris = int(tri_arr.shape[0])
        self.outward_flip = bool(outward_flip)
        self.smooth_normals = bool(smooth_normals)
        self.crease_cos = float(crease_cos)
        self.exclude_boundary = bool(exclude_boundary)
        self.boundary_cos_tol = float(boundary_cos_tol)
        self.max_penetration = -1.0 if max_penetration is None else float(max_penetration)
        self.max_distance = -1.0 if max_distance is None else float(max_distance)

        self.B_verts = bspline_basis_matrix(box, self.tri_nodes).tocsr()
        self._build_topology()
        self._hint_tri: np.ndarray | None = None
        self.last_tri_idx: np.ndarray | None = None

    def _build_topology(self) -> None:
        """Build vertex-to-triangle CSR adjacency and open-boundary edge flags."""
        K = self.n_tris
        N = self.n_verts
        tris = self.tris

        # Vertex-to-triangle CSR
        counts = np.bincount(tris.ravel(), minlength=N).astype(np.int64)
        indptr = np.empty(N + 1, dtype=np.int64)
        indptr[0] = 0
        np.cumsum(counts, out=indptr[1:])
        v_tri_ids = np.empty(3 * K, dtype=np.int64)
        v_corner_idx = np.empty(3 * K, dtype=np.int64)
        cursor = indptr[:-1].copy()
        for k in range(K):
            for c in range(3):
                v = int(tris[k, c])
                pos = int(cursor[v])
                v_tri_ids[pos] = k
                v_corner_idx[pos] = c
                cursor[v] = pos + 1

        self._v_tri_indptr = indptr
        self._v_tri_ids = v_tri_ids
        self._v_corner_idx = v_corner_idx

        # Boundary edges: edge 0=(v0,v1), edge 1=(v1,v2), edge 2=(v2,v0)
        e_pairs = np.vstack(
            [
                np.sort(tris[:, [0, 1]], axis=1),
                np.sort(tris[:, [1, 2]], axis=1),
                np.sort(tris[:, [2, 0]], axis=1),
            ]
        )
        _, inv, edge_counts = np.unique(
            e_pairs, axis=0, return_inverse=True, return_counts=True
        )
        is_bnd_flat = edge_counts[inv] == 1
        self._tri_edge_is_boundary = np.ascontiguousarray(
            np.column_stack(
                [is_bnd_flat[:K], is_bnd_flat[K : 2 * K], is_bnd_flat[2 * K :]]
            ),
            dtype=np.bool_,
        )
        vert_bnd = np.zeros(N, dtype=np.bool_)
        bnd_edges = e_pairs[is_bnd_flat]
        if bnd_edges.size:
            vert_bnd[bnd_edges.ravel()] = True
        self._vert_is_boundary = vert_bnd

    def evaluate(
        self,
        x_s: np.ndarray,
        U_m: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Any]:
        """Project deformed slave points ``x_s (M, 3)`` onto deformed master triangles.

        Returns ``(g_n (M,), n_m (M, 3), x_m (M, 3), ab_m (M, 2), B_m (M, n_ctrl))``.
        """
        import scipy.sparse as sp

        xs = np.ascontiguousarray(x_s, dtype=np.float64)
        if xs.ndim == 1:
            xs = xs.reshape(1, 3)
        M = int(xs.shape[0])
        if M == 0:
            return (
                np.empty(0, dtype=np.float64),
                np.empty((0, 3), dtype=np.float64),
                np.empty((0, 3), dtype=np.float64),
                np.empty((0, 2), dtype=np.float64),
                sp.csr_matrix((0, self.n_ctrl), dtype=np.float64),
            )

        um = np.asarray(U_m, dtype=np.float64).reshape(self.n_ctrl, 3)
        V_def = np.ascontiguousarray(self.tri_nodes + self.B_verts @ um, dtype=np.float64)

        fn, N_corner = _compute_corner_normals_njit(
            V_def,
            self.tris,
            self._v_tri_indptr,
            self._v_tri_ids,
            self._v_corner_idx,
            self.outward_flip,
            self.smooth_normals,
            self.crease_cos,
        )

        if self._hint_tri is None or self._hint_tri.shape[0] != M:
            hint = np.full(M, -1, dtype=np.int64)
        else:
            hint = self._hint_tri

        tri_idx, bary, x_m, n_m, g_n = _project_surface_triangles_njit(
            xs,
            V_def,
            self.tris,
            fn,
            N_corner,
            hint,
            self._tri_edge_is_boundary,
            self._vert_is_boundary,
            self.exclude_boundary,
            self.boundary_cos_tol,
            self.max_penetration,
            self.max_distance,
        )
        self._hint_tri = tri_idx
        self.last_tri_idx = tri_idx

        cols = self.tris[tri_idx].ravel()
        vals = bary.ravel()
        indptr = np.arange(0, 3 * M + 1, 3, dtype=np.int64)
        P = sp.csr_matrix((vals, cols, indptr), shape=(M, self.n_verts), dtype=np.float64)
        B_m = (P @ self.B_verts).tocsr()
        B_m.sort_indices()
        ab_m = np.ascontiguousarray(bary[:, 1:3], dtype=np.float64)
        return g_n, n_m, x_m, ab_m, B_m


# --------------------------------------------------------------------------------------
# Multi-region NLGEOM + J2 isotropic hardening Kratos body operator
# --------------------------------------------------------------------------------------


class NonlinearBodyOperator:
    """Multi-step NLGEOM + J2 isotropic hardening operator on one ``BoxSpline``.

    Supports single-material trimmed solids and multi-region partitioned solids
    (e.g. E05 pin shaft + end caps) under ``TotalLagrangianElement3D4N``.

    ``regions`` is a sequence of dicts, each containing:
    * ``trim_soup``: ``(nodes, tris)`` or ``None`` (untrimmed body)
    * ``E``: Young's modulus
    * ``nu``: Poisson's ratio
    * ``hardening``: optional :class:`~wht_iga.analyses.thermal_plasticity.HardeningTable`
      (uses ``FiniteStrainIsotropicPlasticity3DVonMisesVonMises`` with
      ``HARDENING_CURVE = 6`` when provided)
    * ``law_factory``: optional Kratos law constructor when ``hardening is None``
      (defaults to ``HyperElastic3DLaw``)
    """

    def __init__(
        self,
        KM: Any,
        mods: Mapping[str, Any],
        box: Any,
        regions: Sequence[Mapping[str, Any]],
        *,
        element: str = "TotalLagrangianElement3D4N",
        zero_fn: Callable[[np.ndarray], tuple[int, ...]] | None = None,
        presc_fn: Callable[[np.ndarray], tuple[int, ...]] | None = None,
        q: int = 2,
        ghost_alpha: float = 0.0,
        stiffness_reg: np.ndarray | None = None,
        surface_bc_B: Any | None = None,
        surface_bc_weights: np.ndarray | None = None,
        surface_bc_penalty: float = 0.0,
        surface_bc_target: tuple[float, float, float] | np.ndarray | None = None,
        fracture_energy: float = 1e20,
    ) -> None:
        from wht_iga.solver.active import inactive_ctrl_mask
        from wht_iga.solver.trimmed import build_trimmed_table

        check_qp_element_name(element)
        if not regions:
            raise ValueError("regions must be non-empty")

        self.KM = KM
        self._SM = mods["StructuralMechanicsApplication"]
        self._CL = mods["ConstitutiveLawsApplication"]
        self.box = box
        self.n_ctrl = int(box.n_ctrl)
        self.n_dofs = int(box.n_dofs)
        self.fracture_energy = float(fracture_energy)

        nx, ny, nz = (int(v) for v in box.n_cells)
        n_gauss = int(q) + 1
        all_cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
        out_union: set[tuple[int, int, int]] = set(all_cells)

        region_tables: list[tuple[dict, set[tuple[int, int, int]]]] = []
        n_fallback_total = 0
        n_cut_total = 0
        for reg in regions:
            soup = reg.get("trim_soup")
            if soup is not None:
                qp_table, out_cells, _issues, _rec, diag = build_trimmed_table(
                    soup[0], soup[1], box, q=int(q), cache_dir="off"
                )
                out_set = {tuple(int(v) for v in c) for c in out_cells}
                n_fallback_total += int(diag["n_fallback"])
                n_cut_total += int(diag["n_cut"])
            else:
                qp_table, out_set = {}, set()
            out_union &= out_set
            region_tables.append((qp_table, out_set))

        self.diag = {"n_cut": n_cut_total, "n_fallback": n_fallback_total}
        inact = inactive_ctrl_mask(box, sorted(out_union))

        self.model = KM.Model()
        mp = self.model.CreateModelPart("NonlinearBody")
        self.mp = mp
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.SetBufferSize(2)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        vol, _ = make_nurbs_volume(KM, mp, box)
        self._hardening_groups: list[dict[str, Any]] = []
        self._elements: list[Any] = []
        eid = 1

        for r_idx, (reg, (qp_table, out_set)) in enumerate(zip(regions, region_tables)):
            E = float(reg["E"])
            nu = float(reg["nu"])
            hardening: HardeningTable | None = reg.get("hardening")
            props = mp.CreateNewProperties(r_idx + 1)
            props.SetValue(KM.YOUNG_MODULUS, E)
            props.SetValue(KM.POISSON_RATIO, nu)

            if hardening is not None:
                grid = hardening.plastic_strains
                curve = hardening.curve(0.0)
                _check_point_curve(grid, curve, f"region {r_idx}")
                props.SetValue(KM.YIELD_STRESS, float(curve[0]))
                props.SetValue(KM.FRACTURE_ENERGY, self.fracture_energy)
                props.SetValue(self._CL.HARDENING_CURVE, 6)
                v_eps = KM.Vector(len(grid) + 1)
                for j, val in enumerate(grid):
                    v_eps[j] = float(val)
                v_eps[len(grid)] = float(grid[-1] + FLAT_EXTENSION)
                v_sig = KM.Vector(len(curve) + 1)
                for j, val in enumerate(curve):
                    v_sig[j] = float(val)
                v_sig[len(curve)] = float(curve[-1])
                props.SetValue(self._CL.PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE, v_eps)
                props.SetValue(self._CL.EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE, v_sig)
                props.SetValue(
                    KM.CONSTITUTIVE_LAW,
                    self._CL.FiniteStrainIsotropicPlasticity3DVonMisesVonMises(),
                )
            else:
                factory = reg.get("law_factory") or self._CL.HyperElastic3DLaw
                props.SetValue(KM.CONSTITUTIVE_LAW, factory())

            reg_elems = []
            for i, j, k in all_cells:
                cell = (i, j, k)
                if cell in out_set:
                    continue
                if cell in qp_table:
                    pts = np.asarray(qp_table[cell], dtype=float).tolist()
                else:
                    pts = cell_gauss_points(box, i, j, k, n_gauss=n_gauss)
                if not pts:
                    continue
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                el = mp.CreateNewElement(element, eid, qps[0], props)
                reg_elems.append(el)
                self._elements.append(el)
                eid += 1

            if hardening is not None:
                self._hardening_groups.append(
                    {
                        "elements": reg_elems,
                        "grid": np.asarray(hardening.plastic_strains, dtype=np.float64),
                        "curve": np.asarray(hardening.curve(0.0), dtype=np.float64),
                        "lc": [characteristic_length(el) for el in reg_elems],
                        "eps_p": None,
                        "peeq": None,
                    }
                )

        if float(ghost_alpha) > 0.0 and out_union:
            E0 = float(regions[0]["E"])
            nu0 = float(regions[0]["nu"])
            pg = mp.CreateNewProperties(len(regions) + 1)
            pg.SetValue(KM.YOUNG_MODULUS, E0 * float(ghost_alpha))
            pg.SetValue(KM.POISSON_RATIO, nu0)
            pg.SetValue(KM.CONSTITUTIVE_LAW, self._SM.LinearElastic3DLaw())
            for i, j, k in sorted(out_union):
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(
                    qps, 2, cell_gauss_points(box, i, j, k, n_gauss=n_gauss)
                )
                mp.CreateNewElement(element, eid, qps[0], pg)
                eid += 1

        for v, r in (
            (KM.DISPLACEMENT_X, KM.REACTION_X),
            (KM.DISPLACEMENT_Y, KM.REACTION_Y),
            (KM.DISPLACEMENT_Z, KM.REACTION_Z),
        ):
            KM.VariableUtils().AddDof(v, r, mp)

        self._scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        self._builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
        self._strat = KM.ResidualBasedLinearStrategy(
            mp, self._scheme, self._builder, False, False, False, False
        )
        self._strat.SetEchoLevel(0)
        self._strat.Initialize()

        # Assign equation IDs via step 0, then leave step 0 finalized so begin_step(t) starts cleanly
        mp.CloneTimeStep(0.0)
        self._strat.InitializeSolutionStep()
        self._A = self._strat.GetSystemMatrix()
        self._b = self._strat.GetSystemVector()
        self._Dx = KM.Vector(self._b.Size(), 0.0)
        self._space = KM.UblasSparseSpace()
        perm = equation_permutation(KM, mp, box)
        n_eq = int(self._A.Size1())
        if n_eq == self.n_dofs:
            self._P_act = None
            self._act_dof = None
            self._eq_ids = perm
            used_ctrl_set = set(range(box.n_ctrl))
        else:
            import scipy.sparse as sp

            used_ctrl = sorted(
                {
                    n.Id - 1
                    for el in mp.Elements
                    for n in el.GetNodes()
                    if 0 <= n.Id - 1 < box.n_ctrl
                }
            )
            used_ctrl_set = set(used_ctrl)
            act_dof = np.array([3 * i + d for i in used_ctrl for d in range(3)], dtype=np.int64)
            eq_ids = perm[act_dof]
            self._act_dof = act_dof
            self._eq_ids = eq_ids
            self._P_act = sp.csr_matrix(
                (np.ones(act_dof.size, dtype=np.float64), (act_dof, eq_ids)),
                shape=(self.n_dofs, n_eq),
            )
        self._strat.FinalizeSolutionStep()

        self._ctrl_nodes = sorted(
            (n for n in mp.Nodes if 0 <= n.Id - 1 < box.n_ctrl), key=lambda n: n.Id
        )
        self._X0 = np.array([[n.X0, n.Y0, n.Z0] for n in self._ctrl_nodes], dtype=np.float64)
        self._step_open = False
        self._step_clock = 0.0

        Xcp = box.control_points()
        zero = np.zeros(box.n_dofs, dtype=bool)
        presc = np.zeros(box.n_dofs, dtype=bool)
        for I in range(box.n_ctrl):
            if inact[I] or I not in used_ctrl_set:
                zero[3 * I : 3 * I + 3] = True
        if zero_fn is not None:
            for I in range(box.n_ctrl):
                for d in zero_fn(Xcp[I]):
                    zero[3 * I + d] = True
        if presc_fn is not None:
            for I in range(box.n_ctrl):
                if inact[I]:
                    continue
                for d in presc_fn(Xcp[I]):
                    if not zero[3 * I + d]:
                        presc[3 * I + d] = True
        self.zero_mask = zero
        self.presc_mask = presc
        self.free_mask = ~(zero | presc)

        if stiffness_reg is not None:
            reg_arr = np.asarray(stiffness_reg, dtype=np.float64).reshape(self.n_dofs)
            self._stiffness_reg = reg_arr
        else:
            self._stiffness_reg = None

        if surface_bc_B is not None and surface_bc_weights is not None and float(surface_bc_penalty) > 0.0:
            import scipy.sparse as sp

            w_bc = np.asarray(surface_bc_weights, dtype=np.float64).ravel() * float(surface_bc_penalty)
            gram_scalar = (surface_bc_B.T @ sp.diags(w_bc, format="csr") @ surface_bc_B).tocsr()
            self._K_sbc = sp.kron(gram_scalar, sp.eye(3, format="csr"), format="csr")
            self._u_sbc_target = np.zeros(self.n_dofs, dtype=np.float64)
            if surface_bc_target is not None:
                self.set_surface_bc_target(surface_bc_target)
        else:
            self._K_sbc = None
            self._u_sbc_target = None

    def set_surface_bc_target(
        self, target: tuple[float, float, float] | np.ndarray
    ) -> None:
        """Set uniform 3D displacement target ``(ux, uy, uz)`` or full vector on surface BC."""
        if self._K_sbc is None:
            raise RuntimeError("set_surface_bc_target called without surface_bc_B configured")
        arr = np.asarray(target, dtype=np.float64)
        if arr.size == 3:
            tgt = np.tile(arr.reshape(3), self.n_ctrl)
        elif arr.size == self.n_dofs:
            tgt = arr.reshape(self.n_dofs).copy()
        else:
            raise ValueError(f"Expected 3-vector or ({self.n_dofs},) array, got shape {arr.shape}")
        tgt[self.zero_mask] = 0.0
        self._u_sbc_target = tgt

    def surface_bc_reaction(self, u_full: np.ndarray) -> np.ndarray:
        """Return total 3D reaction force ``-K_sbc @ (u_full - u_target)`` on surface BC."""
        if self._K_sbc is None:
            return np.zeros(3, dtype=np.float64)
        u_arr = np.asarray(u_full, dtype=np.float64).reshape(self.n_dofs)
        r_vec = -(self._K_sbc @ (u_arr - self._u_sbc_target))
        return r_vec.reshape(-1, 3).sum(axis=0)

    def begin_step(self, time_val: float | None = None, *, sync_peeq: bool = True) -> None:
        """Open a new load step on Kratos ModelPart (must precede ``eval``)."""
        if self._step_open:
            raise RuntimeError("begin_step called while previous step is still open")
        t = self._step_clock + 1.0 if time_val is None else float(time_val)
        if t <= self._step_clock:
            t = self._step_clock + 1.0
        self._step_clock = t
        self.mp.CloneTimeStep(t)
        self._strat.InitializeSolutionStep()
        if sync_peeq and self._hardening_groups:
            self._write_dissipation()
        self._step_open = True

    def eval(self, u_full: np.ndarray) -> tuple[Any, np.ndarray]:
        """Evaluate pure trial ``(K_T csr, f_int)`` at ``u_full`` (does not commit history)."""
        import KratosMultiphysics.scipy_conversion_tools as sct

        if not self._step_open:
            self.begin_step(sync_peeq=True)

        KM = self.KM
        U = np.asarray(u_full, dtype=np.float64).reshape(-1, 3)
        for n, x0, uu in zip(self._ctrl_nodes, self._X0, U):
            n.SetSolutionStepValue(
                KM.DISPLACEMENT, KM.Array3([float(uu[0]), float(uu[1]), float(uu[2])])
            )
            n.X, n.Y, n.Z = x0[0] + uu[0], x0[1] + uu[1], x0[2] + uu[2]

        self._scheme.InitializeNonLinIteration(self.mp, self._A, self._Dx, self._b)
        self._space.SetToZeroMatrix(self._A)
        self._space.SetToZeroVector(self._b)
        self._builder.Build(self._scheme, self.mp, self._A, self._b)
        self._scheme.FinalizeNonLinIteration(self.mp, self._A, self._Dx, self._b)

        K_eq = sct.to_csr(self._A)
        b_arr = -np.array(self._b, copy=True)
        if self._P_act is None:
            K = K_eq[self._eq_ids, :][:, self._eq_ids].tocsr()
            f_int = b_arr[self._eq_ids]
        else:
            K = (self._P_act @ K_eq @ self._P_act.T).tocsr()
            f_int = np.zeros(self.n_dofs, dtype=np.float64)
            f_int[self._act_dof] = b_arr[self._eq_ids]
        if self._stiffness_reg is not None:
            import scipy.sparse as sp

            K = (K + sp.diags(self._stiffness_reg, format="csr")).tocsr()
            f_int = f_int + self._stiffness_reg * np.asarray(u_full, dtype=np.float64).reshape(-1)
        if self._K_sbc is not None:
            u_arr = np.asarray(u_full, dtype=np.float64).reshape(self.n_dofs)
            K = (K + self._K_sbc).tocsr()
            f_int = f_int + self._K_sbc @ (u_arr - self._u_sbc_target)
        return K, f_int

    def _read_plastic_strain(self, el: Any) -> np.ndarray:
        return np.array(
            [
                list(v)
                for v in el.CalculateOnIntegrationPoints(
                    self._CL.PLASTIC_STRAIN_VECTOR, self.mp.ProcessInfo
                )
            ],
            dtype=np.float64,
        ).reshape(-1, 6)

    def _ensure_peeq_state(self, grp: dict[str, Any]) -> None:
        if grp["eps_p"] is None:
            eps_p = [self._read_plastic_strain(el) for el in grp["elements"]]
            grp["eps_p"] = eps_p
            grp["peeq"] = [np.zeros(len(e), dtype=np.float64) for e in eps_p]

    def _write_dissipation(self) -> None:
        """Synchronize ``PLASTIC_DISSIPATION`` from accumulated PEEQ before a step."""
        KM = self.KM
        for grp in self._hardening_groups:
            self._ensure_peeq_state(grp)
            grid = grp["grid"]
            curve = grp["curve"]
            for i, el in enumerate(grp["elements"]):
                peeq = grp["peeq"][i]
                kappa = grp["lc"][i] / self.fracture_energy * _curve_work(grid, curve, peeq)
                el.SetValuesOnIntegrationPoints(
                    KM.PLASTIC_DISSIPATION, [float(k) for k in kappa], self.mp.ProcessInfo
                )

    def commit_step(self) -> float:
        """Finalize converged step in Kratos and accumulate PEEQ. Returns max PEEQ."""
        if not self._step_open:
            raise RuntimeError("commit_step called with no open step")
        self._strat.FinalizeSolutionStep()
        self._step_open = False

        max_peeq = 0.0
        for grp in self._hardening_groups:
            self._ensure_peeq_state(grp)
            for i, el in enumerate(grp["elements"]):
                ep = self._read_plastic_strain(el)
                d = ep - grp["eps_p"][i]
                grp["peeq"][i] = grp["peeq"][i] + np.sqrt(
                    (2.0 / 3.0)
                    * (np.sum(d[:, :3] ** 2, axis=1) + 0.5 * np.sum(d[:, 3:] ** 2, axis=1))
                )
                grp["eps_p"][i] = ep
                if grp["peeq"][i].size:
                    max_peeq = max(max_peeq, float(grp["peeq"][i].max()))
        return max_peeq

    def abort_step(self, u_restore: np.ndarray | None = None) -> None:
        """Abort an uncommitted step, restoring nodal coordinates and committed dissipation."""
        if not self._step_open:
            return
        if u_restore is not None:
            KM = self.KM
            U = np.asarray(u_restore, dtype=np.float64).reshape(-1, 3)
            for n, x0, uu in zip(self._ctrl_nodes, self._X0, U):
                n.SetSolutionStepValue(
                    KM.DISPLACEMENT, KM.Array3([float(uu[0]), float(uu[1]), float(uu[2])])
                )
                n.X, n.Y, n.Z = x0[0] + uu[0], x0[1] + uu[1], x0[2] + uu[2]
        if self._hardening_groups:
            self._write_dissipation()
        self._step_open = False

    def max_peeq(self) -> float:
        """Return maximum accumulated PEEQ across all hardening elements."""
        val = 0.0
        for grp in self._hardening_groups:
            if grp["peeq"] is None:
                continue
            for arr in grp["peeq"]:
                if arr.size:
                    val = max(val, float(arr.max()))
        return val

    def max_von_mises(self) -> float:
        """Return maximum Cauchy von Mises stress across all body elements."""
        vm_max = 0.0
        for el in self._elements:
            st = np.array(
                [
                    list(v)
                    for v in el.CalculateOnIntegrationPoints(
                        self.KM.CAUCHY_STRESS_VECTOR, self.mp.ProcessInfo
                    )
                ],
                dtype=np.float64,
            ).reshape(-1, 6)
            if st.size == 0:
                continue
            vm = np.sqrt(
                0.5
                * (
                    (st[:, 0] - st[:, 1]) ** 2
                    + (st[:, 1] - st[:, 2]) ** 2
                    + (st[:, 2] - st[:, 0]) ** 2
                    + 6.0 * (st[:, 3] ** 2 + st[:, 4] ** 2 + st[:, 5] ** 2)
                )
            )
            vm_max = max(vm_max, float(vm.max()))
        return vm_max

    def close(self) -> None:
        """Release Kratos objects to minimize peak RSS across sequential tests."""
        self._elements.clear()
        self._hardening_groups.clear()
        self._ctrl_nodes.clear()
        self._A = None
        self._b = None
        self._Dx = None
        self._strat = None
        self._builder = None
        self._scheme = None
        self.mp = None
        self.model = None
        gc.collect()


# --------------------------------------------------------------------------------------
# Multi-body frictional SSN contact + outer NLGEOM/J2 Newton loop
# --------------------------------------------------------------------------------------


def solve_multibody_ssn_contact(
    K_free_list: Sequence[Any],
    f_ext_list: Sequence[np.ndarray],
    free_masks: Sequence[np.ndarray],
    u_fixed_list: Sequence[np.ndarray],
    pairs: Sequence[Mapping[str, Any]],
    eps_n: float,
    *,
    u_init_list: Sequence[np.ndarray] | None = None,
    mu: float = 0.0,
    eps_t: float | None = None,
    slip_states: Sequence[dict[str, np.ndarray] | None] | None = None,
    max_iter: int = 60,
    rtol: float = 1e-8,
    use_pardiso: bool = True,
) -> dict[str, Any]:
    """Solve coupled N-body, P-pair surface-Gauss SSN contact equilibrium.

    Each pair dict specifies ``name``, ``m_idx``, ``s_idx``, ``master_face``
    (e.g. :class:`DeformableTriangleSurfaceMaster` or
    :class:`~wht_iga.analyses.contact.DeformableBoxFaceMaster`), ``X_s``, ``w_q``,
    and ``B_s``. Convected tangential slip uses the step-start ``slip_states``.
    """
    import scipy.sparse as sp

    if float(eps_n) <= 0.0:
        raise ValueError(f"eps_n must be > 0, got {eps_n}")
    mu = float(mu)
    if mu < 0.0:
        raise ValueError(f"mu must be >= 0, got {mu}")
    if mu > 0.0:
        eps_t = float(eps_n if eps_t is None else eps_t)
        if eps_t <= 0.0:
            raise ValueError(f"eps_t must be > 0 when mu > 0, got {eps_t}")

    n_bodies = len(K_free_list)
    free_arrs = [np.asarray(m, dtype=bool) for m in free_masks]
    n_free = [int(m.sum()) for m in free_arrs]
    offsets = np.cumsum([0] + n_free)
    n_tot = int(offsets[-1])

    K_diag = sp.block_diag([K.tocsr() for K in K_free_list], format="csr")
    f_ext_vec = np.concatenate(
        [np.asarray(f, dtype=np.float64).reshape(nf) for f, nf in zip(f_ext_list, n_free)]
    )

    u_bases = []
    for b_idx, (m, u_fix) in enumerate(zip(free_arrs, u_fixed_list)):
        u_b = np.zeros(m.size, dtype=np.float64)
        u_fix_arr = np.asarray(u_fix, dtype=np.float64).reshape(m.size)
        u_b[~m] = u_fix_arr[~m]
        u_bases.append(u_b)

    if u_init_list is None:
        u = np.zeros(n_tot, dtype=np.float64)
    else:
        u = np.concatenate(
            [np.asarray(ui, dtype=np.float64).reshape(nf) for ui, nf in zip(u_init_list, n_free)]
        ).copy()

    pair_prep = []
    for p_idx, p in enumerate(pairs):
        si = int(p["s_idx"])
        mi = int(p["m_idx"])
        X_s = np.asarray(p["X_s"], dtype=np.float64)
        w_q = np.asarray(p["w_q"], dtype=np.float64)
        M = int(X_s.shape[0])
        B_s_vec = expand_vector_basis(p["B_s"], ndim=3)
        B_s_free = B_s_vec[:, free_arrs[si]].tocsr()
        x_s_base = X_s + (B_s_vec[:, ~free_arrs[si]] @ u_bases[si][~free_arrs[si]]).reshape(M, 3)

        st = None if slip_states is None else slip_states[p_idx]
        n_ctrl_m = int(p["master_face"].n_ctrl)
        if st is None:
            x_s_prev = X_s.copy()
            U_m_prev = np.zeros((n_ctrl_m, 3), dtype=np.float64)
            t_t_prev = np.zeros((M, 3), dtype=np.float64)
        else:
            x_s_prev = np.asarray(st["x_s_prev"], dtype=np.float64).reshape(M, 3)
            U_m_prev = np.asarray(st["U_m_prev"], dtype=np.float64).reshape(n_ctrl_m, 3)
            t_t_prev = np.asarray(st["t_t_prev"], dtype=np.float64).reshape(M, 3)

        pair_prep.append(
            {
                "name": str(p["name"]),
                "mi": mi,
                "si": si,
                "master_face": p["master_face"],
                "M": M,
                "w_q": w_q,
                "B_s_vec": B_s_vec,
                "B_s_free": B_s_free,
                "x_s_base": x_s_base,
                "x_s_prev": x_s_prev,
                "U_m_prev": U_m_prev,
                "t_t_prev": t_t_prev,
                "n_ctrl_m": n_ctrl_m,
            }
        )

    I3 = np.eye(3, dtype=np.float64)

    def _unpack_full(b_idx: int, u_vec: np.ndarray) -> np.ndarray:
        full = u_bases[b_idx].copy()
        full[free_arrs[b_idx]] = u_vec[offsets[b_idx] : offsets[b_idx + 1]]
        return full

    def _evaluate(u_vec: np.ndarray, need_tangent: bool = True):
        f_c_blocks = [np.zeros(nf, dtype=np.float64) for nf in n_free]
        K_blocks: list[list[Any]] = [[None] * n_bodies for _ in range(n_bodies)]
        pair_out: dict[str, dict[str, Any]] = {}
        active_sigs = []
        slip_sigs = []

        for pp in pair_prep:
            mi, si, M = pp["mi"], pp["si"], pp["M"]
            us_f = u_vec[offsets[si] : offsets[si + 1]]
            U_m_curr = _unpack_full(mi, u_vec).reshape(pp["n_ctrl_m"], 3)
            x_s = pp["x_s_base"] + (pp["B_s_free"] @ us_f).reshape(M, 3)

            g_n, n_m, _x_m, _ab_m, B_m = pp["master_face"].evaluate(x_s, U_m_curr)
            active = g_n < 0.0
            act_idx = np.flatnonzero(active)
            slip_mask = np.zeros(M, dtype=bool)
            t_res = np.zeros((M, 3), dtype=np.float64)
            t_t_res = np.zeros((M, 3), dtype=np.float64)

            if act_idx.size > 0:
                w_act = pp["w_q"][act_idx]
                g_act = g_n[act_idx]
                n_act = n_m[act_idx]
                nn = np.einsum("qi,qj->qij", n_act, n_act)
                P_t = I3[None, :, :] - nn
                t_n = (eps_n * g_act)[:, None] * n_act
                k_loc = (w_act * eps_n)[:, None, None] * nn if need_tangent else None

                B_m_act_scalar = B_m[act_idx, :]
                if mu > 0.0:
                    du_m_act = B_m_act_scalar @ (U_m_curr - pp["U_m_prev"])
                    du_rel = (x_s[act_idx] - pp["x_s_prev"][act_idx]) - du_m_act
                    du_t = np.einsum("qij,qj->qi", P_t, du_rel)
                    t_prev_proj = np.einsum("qij,qj->qi", P_t, pp["t_t_prev"][act_idx])
                    t_trial = t_prev_proj + float(eps_t) * du_t
                    tn = np.linalg.norm(t_trial, axis=1)
                    cap = -mu * eps_n * g_act
                    slip = tn > cap
                    slip_mask[act_idx] = slip

                    t_t_act = t_trial.copy()
                    if need_tangent:
                        k_t = (w_act * float(eps_t))[:, None, None] * P_t
                    if np.any(slip):
                        s_idx = np.flatnonzero(slip)
                        d_hat = t_trial[s_idx] / tn[s_idx, None]
                        t_t_act[s_idx] = cap[s_idx, None] * d_hat
                        if need_tangent:
                            dd = np.einsum("qi,qj->qij", d_hat, d_hat)
                            dn = np.einsum("qi,qj->qij", d_hat, n_act[s_idx])
                            k_t[s_idx] = w_act[s_idx, None, None] * (
                                (float(eps_t) * (cap[s_idx] / tn[s_idx]))[:, None, None]
                                * (P_t[s_idx] - dd)
                                - (mu * eps_n) * dn
                            )
                    t_t_res[act_idx] = t_t_act
                    t_res[act_idx] = w_act[:, None] * (t_n + t_t_act)
                    if need_tangent:
                        k_loc = k_loc + k_t
                else:
                    t_res[act_idx] = w_act[:, None] * t_n

                rows = (3 * act_idx[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
                B_s_act = pp["B_s_free"][rows, :]
                B_m_act = expand_vector_basis(B_m_act_scalar, ndim=3)[:, free_arrs[mi]].tocsr()
                t_act_flat = t_res[act_idx].ravel()
                f_c_blocks[mi] -= B_m_act.T @ t_act_flat
                f_c_blocks[si] += B_s_act.T @ t_act_flat

                if need_tangent:
                    n_a = act_idx.size
                    K_blk = sp.bsr_matrix(
                        (k_loc, np.arange(n_a, dtype=np.int64), np.arange(n_a + 1, dtype=np.int64)),
                        shape=(3 * n_a, 3 * n_a),
                    ).tocsr()
                    K_mm = (B_m_act.T @ K_blk @ B_m_act).tocsr()
                    K_ms = -(B_m_act.T @ K_blk @ B_s_act).tocsr()
                    K_sm = -(B_s_act.T @ K_blk @ B_m_act).tocsr()
                    K_ss = (B_s_act.T @ K_blk @ B_s_act).tocsr()
                    for r_i, c_i, blk in ((mi, mi, K_mm), (mi, si, K_ms), (si, mi, K_sm), (si, si, K_ss)):
                        K_blocks[r_i][c_i] = blk if K_blocks[r_i][c_i] is None else (K_blocks[r_i][c_i] + blk)

            active_sigs.append(active)
            slip_sigs.append(slip_mask)
            pair_out[pp["name"]] = {
                "n_active": int(act_idx.size),
                "n_slip": int(slip_mask.sum()),
                "F_total": t_res.sum(axis=0),
                "F_t": (pp["w_q"][:, None] * t_t_res).sum(axis=0),
                "normal_pressure_integral": float(np.sum(-pp["w_q"][act_idx] * eps_n * g_n[act_idx]))
                if act_idx.size
                else 0.0,
                "x_s": x_s,
                "U_m_curr": U_m_curr,
                "t_t_res": t_t_res,
                "B_m": B_m,
                "act_idx": act_idx,
                "t_res": t_res,
            }

        f_c = np.concatenate(f_c_blocks)
        res = K_diag @ u_vec + f_c - f_ext_vec
        denom = max(1.0, float(np.linalg.norm(f_ext_vec)), float(np.linalg.norm(f_c)))
        rnorm = float(np.linalg.norm(res) / denom)

        K_tan = None
        if need_tangent:
            for r_i in range(n_bodies):
                for c_i in range(n_bodies):
                    if K_blocks[r_i][c_i] is None:
                        K_blocks[r_i][c_i] = sp.csr_matrix(
                            (n_free[r_i], n_free[c_i]), dtype=np.float64
                        )
            K_tan = (K_diag + sp.bmat(K_blocks, format="csr")).tocsr()

        return res, rnorm, K_tan, pair_out, active_sigs, slip_sigs

    prev_active: list[np.ndarray] | None = None
    prev_slip: list[np.ndarray] | None = None
    prev2_active: list[np.ndarray] | None = None
    prev2_slip: list[np.ndarray] | None = None
    converged = False
    n_solves = 0
    rnorm_init = 0.0

    for it_ssn in range(1, int(max_iter) + 1):
        res, rnorm, K_tan, pair_out, active_sigs, slip_sigs = _evaluate(u, need_tangent=True)
        if it_ssn == 1:
            rnorm_init = float(rnorm)
        same_act = prev_active is not None and all(
            np.array_equal(a, b) for a, b in zip(active_sigs, prev_active)
        )
        same_slp = prev_slip is not None and all(
            np.array_equal(a, b) for a, b in zip(slip_sigs, prev_slip)
        )
        cyc2_act = (
            it_ssn >= 25
            and prev2_active is not None
            and all(np.array_equal(a, b) for a, b in zip(active_sigs, prev2_active))
        )
        cyc2_slp = (
            it_ssn >= 25
            and prev2_slip is not None
            and all(np.array_equal(a, b) for a, b in zip(slip_sigs, prev2_slip))
        )
        if rnorm <= rtol and ((same_act and same_slp) or (cyc2_act and cyc2_slp)):
            converged = True
            break

        prev2_active = prev_active
        prev2_slip = prev_slip
        prev_active = [a.copy() for a in active_sigs]
        prev_slip = [s.copy() for s in slip_sigs]

        du = _solve_linear_step(K_tan, -res, use_pardiso=use_pardiso)
        n_solves += 1
        u_next = u + du

        if mu > 0.0 and float(np.linalg.norm(du)) > 0.0:
            res_norm0 = float(np.linalg.norm(res))
            alpha = 1.0
            best_alpha = 1.0
            best_norm = np.inf
            accepted = False
            while alpha >= 1.0 / 64.0:
                res_try, _, _, _, _, _ = _evaluate(u + alpha * du, need_tangent=False)
                norm_try = float(np.linalg.norm(res_try))
                if norm_try < best_norm:
                    best_norm = norm_try
                    best_alpha = alpha
                if norm_try < (1.0 - 1e-4 * alpha) * res_norm0:
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                alpha = best_alpha if (it_ssn >= 20 or best_norm < 2.0 * res_norm0) else 0.5
            u_next = u + alpha * du

        u = u_next

    if not converged:
        raise RuntimeError(
            f"Multi-body SSN contact did not converge in {max_iter} iterations (rnorm={rnorm:.3e})"
        )

    u_free_list = [u[offsets[i] : offsets[i + 1]].copy() for i in range(n_bodies)]
    u_full_list = [_unpack_full(i, u) for i in range(n_bodies)]
    f_c_full_list = [np.zeros(m.size, dtype=np.float64) for m in free_arrs]
    next_slip_states = []
    pair_stats = {}

    for pp in pair_prep:
        po = pair_out[pp["name"]]
        act_idx = po["act_idx"]
        mi, si = pp["mi"], pp["si"]
        if act_idx.size > 0:
            t_act_flat = po["t_res"][act_idx].ravel()
            B_m_vec_act = expand_vector_basis(po["B_m"][act_idx, :], ndim=3)
            f_c_full_list[mi] -= B_m_vec_act.T @ t_act_flat
            rows_full = (3 * act_idx[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
            f_c_full_list[si] += pp["B_s_vec"][rows_full, :].T @ t_act_flat
        next_slip_states.append(
            {
                "x_s_prev": po["x_s"].copy(),
                "U_m_prev": po["U_m_curr"].copy(),
                "t_t_prev": po["t_t_res"].copy(),
            }
        )
        pair_stats[pp["name"]] = {
            "n_active": po["n_active"],
            "n_slip": po["n_slip"],
            "F_total": po["F_total"],
            "F_t": po["F_t"],
            "normal_pressure_integral": po["normal_pressure_integral"],
        }

    return {
        "u_free_list": u_free_list,
        "u_full_list": u_full_list,
        "f_c_full_list": f_c_full_list,
        "pair_stats": pair_stats,
        "slip_states": next_slip_states,
        "iterations": int(n_solves),
        "rnorm": float(rnorm),
        "rnorm_init": float(rnorm_init),
    }


def solve_multibody_nlgeom_contact_step(
    bodies: Sequence[NonlinearBodyOperator],
    u_full_list: Sequence[np.ndarray],
    pairs: Sequence[Mapping[str, Any]],
    eps_n: float,
    *,
    step_time: float | None = None,
    mu: float = 0.0,
    eps_t: float | None = None,
    slip_states: Sequence[dict[str, np.ndarray] | None] | None = None,
    outer_tol: float = 1e-6,
    max_outer: int = 30,
    rtol: float = 1e-8,
    ssn_max_iter: int = 150,
    use_pardiso: bool = True,
    sync_peeq: bool = True,
) -> dict[str, Any]:
    """Run one coupled NLGEOM / J2 plasticity load step across ``bodies`` and ``pairs``.

    Each body must already have this step's prescribed DOFs placed in
    ``u_full_list[i][body.presc_mask]`` (and zero DOFs at 0). Calls
    ``begin_step`` on every body, iterates outer Newton wrapping
    :func:`solve_multibody_ssn_contact` until ``outer_rel <= outer_tol``, and
    calls ``commit_step`` on convergence.
    """
    u_orig = [np.asarray(u, dtype=np.float64).copy() for u in u_full_list]
    u_curr = [u.copy() for u in u_orig]
    free_masks = [b.free_mask for b in bodies]
    n_free_list = [int(fm.sum()) for fm in free_masks]
    offsets_f = np.cumsum([0] + n_free_list)
    for b in bodies:
        b.begin_step(step_time, sync_peeq=sync_peeq)

    rel = np.inf
    r0 = np.inf
    res: dict[str, Any] | None = None
    prev_x: np.ndarray | None = None
    prev_du: np.ndarray | None = None
    aa_dX: list[np.ndarray] = []
    aa_dG: list[np.ndarray] = []
    aa_m = 4
    omega = 1.0
    try:
        for n_outer in range(1, int(max_outer) + 1):
            K_free_list = []
            f_ext_list = []
            u_free_curr = []
            for b, u_f in zip(bodies, u_curr):
                K_b, f_int_b = b.eval(u_f)
                fm = b.free_mask
                Kf_b = K_b[fm][:, fm].tocsr()
                u_fb = u_f[fm].copy()
                K_free_list.append(Kf_b)
                f_ext_list.append(Kf_b @ u_fb - f_int_b[fm])
                u_free_curr.append(u_fb)

            res = solve_multibody_ssn_contact(
                K_free_list,
                f_ext_list,
                free_masks,
                u_curr,
                pairs,
                eps_n,
                u_init_list=u_free_curr,
                mu=mu,
                eps_t=eps_t,
                slip_states=slip_states,
                rtol=rtol,
                max_iter=int(ssn_max_iter),
                use_pardiso=use_pardiso,
            )

            x_vec = np.concatenate(u_free_curr)
            du_parts = [
                res["u_free_list"][i] - u_free_curr[i] for i in range(len(bodies))
            ]
            du_vec = np.concatenate(du_parts)
            du_norm = float(np.linalg.norm(du_vec))
            u_norm = float(np.linalg.norm(np.concatenate(u_curr)))
            rel = du_norm / max(u_norm, 1e-300)
            r0 = float(res.get("rnorm_init", np.inf))

            if rel <= outer_tol or (n_outer >= 2 and r0 <= outer_tol and rel <= 10.0 * outer_tol):
                for i, b in enumerate(bodies):
                    u_curr[i][b.free_mask] = res["u_free_list"][i]
                break

            if prev_x is not None and prev_du is not None:
                aa_dX.append(x_vec - prev_x)
                aa_dG.append(du_vec - prev_du)
                if len(aa_dX) > aa_m:
                    aa_dX.pop(0)
                    aa_dG.pop(0)
                if n_outer >= 6:
                    dot_prev = float(np.dot(du_vec, prev_du))
                    if dot_prev < 0.0:
                        diff = du_vec - prev_du
                        d2 = float(np.dot(diff, diff))
                        if d2 > 1e-30:
                            omega = float(np.clip(-omega * float(np.dot(prev_du, diff)) / d2, 0.2, 1.0))
                    elif omega < 1.0:
                        omega = min(1.0, 1.25 * omega)

            prev_x = x_vec
            prev_du = du_vec

            x_next: np.ndarray | None = None
            if n_outer >= 4 and rel < 5e-2 and len(aa_dX) >= 2:
                dG = np.column_stack(aa_dG)
                dX = np.column_stack(aa_dX)
                G_mat = dG.T @ dG
                reg = 1e-12 * max(float(np.trace(G_mat)), 1e-30)
                try:
                    gamma = np.linalg.solve(G_mat + reg * np.eye(G_mat.shape[0]), dG.T @ du_vec)
                    cand = (x_vec - dX @ gamma) + (du_vec - dG @ gamma)
                    step_aa = cand - x_vec
                    if np.all(np.isfinite(step_aa)) and float(np.linalg.norm(step_aa)) <= 3.0 * max(du_norm, 1e-30):
                        x_next = cand
                except np.linalg.LinAlgError:
                    x_next = None

            if x_next is not None:
                for i, b in enumerate(bodies):
                    u_curr[i][b.free_mask] = x_next[offsets_f[i] : offsets_f[i + 1]]
            else:
                for i, b in enumerate(bodies):
                    u_curr[i][b.free_mask] = u_free_curr[i] + omega * du_parts[i]
        else:
            raise RuntimeError(
                f"Multi-body NLGEOM outer Newton did not converge in {max_outer} iterations "
                f"(rel={rel:.2e}, r0={r0:.2e})"
            )
    except Exception:
        for b, u_0 in zip(bodies, u_orig):
            b.abort_step(u_0)
        raise

    f_int_list = []
    peeq_max_list = []
    for b, u_f in zip(bodies, u_curr):
        _, f_int_b = b.eval(u_f)
        f_int_list.append(f_int_b)
        peeq_max_list.append(b.commit_step())

    return {
        **res,
        "u_full_list": u_curr,
        "f_int_list": f_int_list,
        "n_outer": int(n_outer),
        "outer_rel": float(rel),
        "outer_rnorm": float(r0),
        "peeq_max_list": peeq_max_list,
    }
