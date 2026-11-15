"""Framework-wide Numba spatial acceleration for triangulated surfaces.

Provides ``TriSoupAccel`` with two compiled data structures built once per
``(nodes, tris)`` mesh:

1. **Uniform 3D grid CSR index** (``tri_offsets``, ``tri_ids``): narrows any
   axis-aligned bounding box query from ``O(M)`` triangles to ``O(k)`` local
   overlapping candidates using an ``@njit`` prefix-sum CSR table.
2. **Flat Binary BVH with Dipole + Quadrupole Winding Moments** (Barill et al.
   2018): accelerates ``winding_inside`` and far-field solid-angle queries from
   ``O(M)`` per point to ``O(log M)``, falling back to exact Van Oosterom–
   Strackee evaluation for near-boundary or near-threshold points.
"""

from __future__ import annotations

from collections import OrderedDict
import threading

import numpy as np
from numba import njit, prange  # type: ignore[import-untyped]


@njit(cache=True)
def _mesh_fingerprint_njit(nodes: np.ndarray, tris: np.ndarray) -> int:
    """Fast deterministic uint64-style int64 hash of node and tri contents."""
    h = np.uint64(14695981039346656037)
    p = np.uint64(1099511628211)
    nr = nodes.shape[0]
    nc = nodes.shape[1]
    step_n = 1 if nr <= 4096 else max(1, nr // 4096)
    for i in range(0, nr, step_n):
        for j in range(nc):
            bits = np.float64(nodes[i, j]).view(np.uint64)
            h = (h ^ bits) * p
    tr = tris.shape[0]
    tc = tris.shape[1]
    step_t = 1 if tr <= 4096 else max(1, tr // 4096)
    for i in range(0, tr, step_t):
        for j in range(tc):
            h = (h ^ np.uint64(tris[i, j])) * p
    h = (h ^ np.uint64(nr)) * p
    h = (h ^ np.uint64(tr)) * p
    return np.int64(h)


@njit(cache=True, parallel=True)
def _pack_tris_njit(nodes: np.ndarray, tris: np.ndarray):
    """Pack triangle vertices, AABBs, signed area vectors, and centroids."""
    m = tris.shape[0]
    soup = np.empty((m, 3, 3), dtype=np.float64)
    tmin = np.empty((m, 3), dtype=np.float64)
    tmax = np.empty((m, 3), dtype=np.float64)
    area_vec = np.empty((m, 3), dtype=np.float64)
    centroid = np.empty((m, 3), dtype=np.float64)
    area_abs = np.empty(m, dtype=np.float64)

    for i in prange(m):
        i0 = tris[i, 0]
        i1 = tris[i, 1]
        i2 = tris[i, 2]
        for d in range(3):
            v0 = nodes[i0, d]
            v1 = nodes[i1, d]
            v2 = nodes[i2, d]
            soup[i, 0, d] = v0
            soup[i, 1, d] = v1
            soup[i, 2, d] = v2
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
            tmin[i, d] = lo
            tmax[i, d] = hi
            centroid[i, d] = (v0 + v1 + v2) / 3.0

        ax = soup[i, 1, 0] - soup[i, 0, 0]
        ay = soup[i, 1, 1] - soup[i, 0, 1]
        az = soup[i, 1, 2] - soup[i, 0, 2]
        bx = soup[i, 2, 0] - soup[i, 0, 0]
        by = soup[i, 2, 1] - soup[i, 0, 1]
        bz = soup[i, 2, 2] - soup[i, 0, 2]
        nx = 0.5 * (ay * bz - az * by)
        ny = 0.5 * (az * bx - ax * bz)
        nz = 0.5 * (ax * by - ay * bx)
        area_vec[i, 0] = nx
        area_vec[i, 1] = ny
        area_vec[i, 2] = nz
        area_abs[i] = np.sqrt(nx * nx + ny * ny + nz * nz)

    return soup, tmin, tmax, area_vec, centroid, area_abs


@njit(cache=True, nogil=True)
def _tri_bin_bounds(tmin_val: float, tmax_val: float, glo: float, inv_h: float, g_max: int) -> tuple[int, int]:
    """Compute inclusive integer grid bin range [b0, b1] for a 1D interval."""
    b0 = int(np.floor((tmin_val - glo) * inv_h))
    b1 = int(np.floor((tmax_val - glo) * inv_h))
    if b0 < 0:
        b0 = 0
    elif b0 > g_max:
        b0 = g_max
    if b1 < 0:
        b1 = 0
    elif b1 > g_max:
        b1 = g_max
    return b0, b1


@njit(cache=True)
def _build_grid_csr_njit(tmin: np.ndarray, tmax: np.ndarray, glo: np.ndarray, inv_h: np.ndarray,
                         gx: int, gy: int, gz: int):
    """Build 3D uniform-grid CSR index (offsets, tri_ids) over triangle AABBs."""
    m = tmin.shape[0]
    nvox = gx * gy * gz
    counts = np.zeros(nvox, dtype=np.int64)
    gx1, gy1, gz1 = gx - 1, gy - 1, gz - 1
    gl0, gl1, gl2 = glo[0], glo[1], glo[2]
    ih0, ih1, ih2 = inv_h[0], inv_h[1], inv_h[2]

    for t in range(m):
        x0, x1 = _tri_bin_bounds(tmin[t, 0], tmax[t, 0], gl0, ih0, gx1)
        y0, y1 = _tri_bin_bounds(tmin[t, 1], tmax[t, 1], gl1, ih1, gy1)
        z0, z1 = _tri_bin_bounds(tmin[t, 2], tmax[t, 2], gl2, ih2, gz1)
        for kz in range(z0, z1 + 1):
            base_z = kz * gy
            for jy in range(y0, y1 + 1):
                base_y = (base_z + jy) * gx
                for ix in range(x0, x1 + 1):
                    counts[base_y + ix] += 1

    offsets = np.empty(nvox + 1, dtype=np.int64)
    offsets[0] = 0
    for v in range(nvox):
        offsets[v + 1] = offsets[v] + counts[v]

    tri_ids = np.empty(offsets[nvox], dtype=np.int64)
    cursor = offsets[:-1].copy()

    for t in range(m):
        x0, x1 = _tri_bin_bounds(tmin[t, 0], tmax[t, 0], gl0, ih0, gx1)
        y0, y1 = _tri_bin_bounds(tmin[t, 1], tmax[t, 1], gl1, ih1, gy1)
        z0, z1 = _tri_bin_bounds(tmin[t, 2], tmax[t, 2], gl2, ih2, gz1)
        for kz in range(z0, z1 + 1):
            base_z = kz * gy
            for jy in range(y0, y1 + 1):
                base_y = (base_z + jy) * gx
                for ix in range(x0, x1 + 1):
                    v = base_y + ix
                    tri_ids[cursor[v]] = t
                    cursor[v] += 1

    return offsets, tri_ids


@njit(cache=True, nogil=True)
def query_box_candidates_njit(blo: np.ndarray, bhi: np.ndarray,
                              tmin: np.ndarray, tmax: np.ndarray,
                              glo: np.ndarray, inv_h: np.ndarray,
                              gx: int, gy: int, gz: int,
                              offsets: np.ndarray, tri_ids: np.ndarray,
                              stamp_buf: np.ndarray, stamp_val: int,
                              out_buf: np.ndarray) -> int:
    """Collect deduplicated candidate triangle IDs whose AABB overlaps [blo, bhi].

    Caller supplies pre-allocated ``stamp_buf`` (len M, int32/int64) and ``out_buf``
    and increments ``stamp_val`` on each call (never zeroes ``stamp_buf``).
    Returns number of candidates written to ``out_buf``. If ``out_buf`` fills,
    caller should grow and retry (handled in safe wrapper).
    """
    x0, x1 = _tri_bin_bounds(blo[0], bhi[0], glo[0], inv_h[0], gx - 1)
    y0, y1 = _tri_bin_bounds(blo[1], bhi[1], glo[1], inv_h[1], gy - 1)
    z0, z1 = _tri_bin_bounds(blo[2], bhi[2], glo[2], inv_h[2], gz - 1)

    ax0, ay0, az0 = blo[0], blo[1], blo[2]
    ax1, ay1, az1 = bhi[0], bhi[1], bhi[2]
    cap = out_buf.shape[0]
    n = 0
    for kz in range(z0, z1 + 1):
        base_z = kz * gy
        for jy in range(y0, y1 + 1):
            base_y = (base_z + jy) * gx
            for ix in range(x0, x1 + 1):
                v = base_y + ix
                s = offsets[v]
                e = offsets[v + 1]
                for p in range(s, e):
                    tid = tri_ids[p]
                    if stamp_buf[tid] == stamp_val:
                        continue
                    stamp_buf[tid] = stamp_val
                    if (tmax[tid, 0] >= ax0 and tmin[tid, 0] <= ax1
                            and tmax[tid, 1] >= ay0 and tmin[tid, 1] <= ay1
                            and tmax[tid, 2] >= az0 and tmin[tid, 2] <= az1):
                        if n < cap:
                            out_buf[n] = tid
                        n += 1
    return n


@njit(cache=True)
def _partition_median_njit(order: np.ndarray, vals: np.ndarray, start: int, end: int) -> int:
    """Quickselect-style median partition of ``order[start:end]`` by ``vals[order[i]]``."""
    mid = (start + end) // 2
    lo = start
    hi = end - 1
    while lo < hi:
        piv_idx = order[(lo + hi) // 2]
        piv = vals[piv_idx]
        i = lo
        j = hi
        while i <= j:
            while vals[order[i]] < piv:
                i += 1
            while vals[order[j]] > piv:
                j -= 1
            if i <= j:
                tmp = order[i]
                order[i] = order[j]
                order[j] = tmp
                i += 1
                j -= 1
        if j < mid:
            lo = i
        if i > mid:
            hi = j
    return mid


@njit(cache=True)
def _build_bvh_njit(soup: np.ndarray, tmin: np.ndarray, tmax: np.ndarray,
                    area_vec: np.ndarray, centroid: np.ndarray, area_abs: np.ndarray,
                    leaf_size: int = 8):
    """Build top-down binary BVH with cluster dipole + quadrupole winding moments."""
    m = soup.shape[0]
    order = np.arange(m, dtype=np.int64)
    max_nodes = max(1, 2 * m)

    bmin = np.empty((max_nodes, 3), dtype=np.float64)
    bmax = np.empty((max_nodes, 3), dtype=np.float64)
    left = np.full(max_nodes, -1, dtype=np.int64)
    right = np.full(max_nodes, -1, dtype=np.int64)
    tri_start = np.empty(max_nodes, dtype=np.int64)
    tri_count = np.empty(max_nodes, dtype=np.int64)
    center = np.empty((max_nodes, 3), dtype=np.float64)
    rad = np.empty(max_nodes, dtype=np.float64)
    dipole = np.empty((max_nodes, 3), dtype=np.float64)
    quad = np.empty((max_nodes, 9), dtype=np.float64)

    # Explicit build stack of (node_idx, start, end)
    st_node = np.empty(max_nodes, dtype=np.int64)
    st_s = np.empty(max_nodes, dtype=np.int64)
    st_e = np.empty(max_nodes, dtype=np.int64)

    n_nodes = 1
    sp = 0
    st_node[0] = 0
    st_s[0] = 0
    st_e[0] = m
    sp = 1

    # Reusable axis-coord buffer
    cx_buf = centroid[:, 0]
    cy_buf = centroid[:, 1]
    cz_buf = centroid[:, 2]

    while sp > 0:
        sp -= 1
        nid = st_node[sp]
        s = st_s[sp]
        e = st_e[sp]
        cnt = e - s

        # Bounding box and multipole moments
        x0 = 1e300
        y0 = 1e300
        z0 = 1e300
        x1 = -1e300
        y1 = -1e300
        z1 = -1e300
        wsum = 0.0
        px = 0.0
        py = 0.0
        pz = 0.0
        d0 = 0.0
        d1 = 0.0
        d2 = 0.0

        for i in range(s, e):
            tid = order[i]
            if tmin[tid, 0] < x0:
                x0 = tmin[tid, 0]
            if tmin[tid, 1] < y0:
                y0 = tmin[tid, 1]
            if tmin[tid, 2] < z0:
                z0 = tmin[tid, 2]
            if tmax[tid, 0] > x1:
                x1 = tmax[tid, 0]
            if tmax[tid, 1] > y1:
                y1 = tmax[tid, 1]
            if tmax[tid, 2] > z1:
                z1 = tmax[tid, 2]
            aa = area_abs[tid]
            wsum += aa
            px += aa * centroid[tid, 0]
            py += aa * centroid[tid, 1]
            pz += aa * centroid[tid, 2]
            d0 += area_vec[tid, 0]
            d1 += area_vec[tid, 1]
            d2 += area_vec[tid, 2]

        bmin[nid, 0] = x0
        bmin[nid, 1] = y0
        bmin[nid, 2] = z0
        bmax[nid, 0] = x1
        bmax[nid, 1] = y1
        bmax[nid, 2] = z1

        if wsum > 0.0:
            px /= wsum
            py /= wsum
            pz /= wsum
        else:
            px = 0.5 * (x0 + x1)
            py = 0.5 * (y0 + y1)
            pz = 0.5 * (z0 + z1)

        center[nid, 0] = px
        center[nid, 1] = py
        center[nid, 2] = pz
        dipole[nid, 0] = d0
        dipole[nid, 1] = d1
        dipole[nid, 2] = d2

        # Compute bounding radius from weighted centroid & quadrupole matrix A2[j, k]
        max_r2 = 0.0
        q00 = 0.0
        q01 = 0.0
        q02 = 0.0
        q10 = 0.0
        q11 = 0.0
        q12 = 0.0
        q20 = 0.0
        q21 = 0.0
        q22 = 0.0
        for i in range(s, e):
            tid = order[i]
            for k in range(3):
                vx = soup[tid, k, 0] - px
                vy = soup[tid, k, 1] - py
                vz = soup[tid, k, 2] - pz
                r2 = vx * vx + vy * vy + vz * vz
                if r2 > max_r2:
                    max_r2 = r2
            dcx = centroid[tid, 0] - px
            dcy = centroid[tid, 1] - py
            dcz = centroid[tid, 2] - pz
            ax = area_vec[tid, 0]
            ay = area_vec[tid, 1]
            az = area_vec[tid, 2]
            q00 += ax * dcx
            q01 += ax * dcy
            q02 += ax * dcz
            q10 += ay * dcx
            q11 += ay * dcy
            q12 += ay * dcz
            q20 += az * dcx
            q21 += az * dcy
            q22 += az * dcz

        rad[nid] = np.sqrt(max_r2)
        quad[nid, 0] = q00
        quad[nid, 1] = q01
        quad[nid, 2] = q02
        quad[nid, 3] = q10
        quad[nid, 4] = q11
        quad[nid, 5] = q12
        quad[nid, 6] = q20
        quad[nid, 7] = q21
        quad[nid, 8] = q22
        tri_start[nid] = s
        tri_count[nid] = cnt

        if cnt <= leaf_size:
            continue

        # Split along longest centroid extent (handles coincident bounding vertices too)
        dx = x1 - x0
        dy = y1 - y0
        dz = z1 - z0
        ax = 0
        if dy > dx and dy >= dz:
            ax = 1
        elif dz > dx and dz > dy:
            ax = 2

        if ax == 0:
            mid = _partition_median_njit(order, cx_buf, s, e)
        elif ax == 1:
            mid = _partition_median_njit(order, cy_buf, s, e)
        else:
            mid = _partition_median_njit(order, cz_buf, s, e)

        if mid <= s or mid >= e:
            mid = (s + e) // 2

        l_id = n_nodes
        r_id = n_nodes + 1
        n_nodes += 2
        left[nid] = l_id
        right[nid] = r_id

        st_node[sp] = r_id
        st_s[sp] = mid
        st_e[sp] = e
        sp += 1
        st_node[sp] = l_id
        st_s[sp] = s
        st_e[sp] = mid
        sp += 1

    # Reorder soup in BVH order for sequential leaf memory access
    bvh_soup = np.empty((m, 3, 3), dtype=np.float64)
    for i in range(m):
        tid = order[i]
        for k in range(3):
            for d in range(3):
                bvh_soup[i, k, d] = soup[tid, k, d]

    return (bvh_soup, left[:n_nodes], right[:n_nodes],
            tri_start[:n_nodes], tri_count[:n_nodes],
            center[:n_nodes], rad[:n_nodes],
            dipole[:n_nodes], quad[:n_nodes])


@njit(cache=True, inline="always", nogil=True)
def _winding_bvh_one(px: float, py: float, pz: float,
                     bvh_soup: np.ndarray, left: np.ndarray, right: np.ndarray,
                     tri_start: np.ndarray, tri_count: np.ndarray,
                     center: np.ndarray, rad: np.ndarray,
                     dipole: np.ndarray, quad: np.ndarray,
                     beta: float) -> float:
    """Fast BVH winding number (Barill et al. 2018 dipole + quadrupole) for 1 point."""
    inv_4pi = 1.0 / (4.0 * np.pi)
    beta2 = beta * beta
    stack = np.empty(128, dtype=np.int64)
    stack[0] = 0
    sp = 1
    acc = 0.0

    while sp > 0:
        sp -= 1
        nid = stack[sp]
        dx = center[nid, 0] - px
        dy = center[nid, 1] - py
        dz = center[nid, 2] - pz
        r2 = dx * dx + dy * dy + dz * dz
        rr = rad[nid]
        if rr == 0.0 or r2 > beta2 * rr * rr:
            # Far-field multipole expansion
            if r2 > 0.0:
                inv_r = 1.0 / np.sqrt(r2)
                inv_r3 = inv_r / r2
                inv_r5 = inv_r3 / r2
                # Dipole: A1 . r0 / r^3
                term1 = (dipole[nid, 0] * dx + dipole[nid, 1] * dy + dipole[nid, 2] * dz) * inv_r3
                # Quadrupole: tr(A2) / r^3 - 3 * r0^T A2 r0 / r^5
                tr_q = quad[nid, 0] + quad[nid, 4] + quad[nid, 8]
                rqr = (
                    dx * (quad[nid, 0] * dx + quad[nid, 1] * dy + quad[nid, 2] * dz)
                    + dy * (quad[nid, 3] * dx + quad[nid, 4] * dy + quad[nid, 5] * dz)
                    + dz * (quad[nid, 6] * dx + quad[nid, 7] * dy + quad[nid, 8] * dz)
                )
                acc += (term1 + tr_q * inv_r3 - 3.0 * rqr * inv_r5) * inv_4pi
            continue

        lid = left[nid]
        if lid < 0:
            # Leaf: exact Van Oosterom-Strackee solid angle
            s = tri_start[nid]
            e = s + tri_count[nid]
            loc_s = 0.0
            for m in range(s, e):
                ax = bvh_soup[m, 0, 0] - px
                ay = bvh_soup[m, 0, 1] - py
                az = bvh_soup[m, 0, 2] - pz
                bx = bvh_soup[m, 1, 0] - px
                by = bvh_soup[m, 1, 1] - py
                bz = bvh_soup[m, 1, 2] - pz
                cx = bvh_soup[m, 2, 0] - px
                cy = bvh_soup[m, 2, 1] - py
                cz = bvh_soup[m, 2, 2] - pz
                la = np.sqrt(ax * ax + ay * ay + az * az)
                lb = np.sqrt(bx * bx + by * by + bz * bz)
                lc = np.sqrt(cx * cx + cy * cy + cz * cz)
                num = (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz)
                       + az * (bx * cy - by * cx))
                den = (la * lb * lc + (ax * bx + ay * by + az * bz) * lc
                       + (bx * cx + by * cy + bz * cz) * la
                       + (cx * ax + cy * ay + cz * az) * lb)
                loc_s += 2.0 * np.arctan2(num, den)
            acc += loc_s * inv_4pi
        else:
            rid = right[nid]
            stack[sp] = rid
            stack[sp + 1] = lid
            sp += 2

    return acc


@njit(cache=True, parallel=True)
def _winding_bvh_batch_parallel_njit(
    pts: np.ndarray,
    bvh_soup: np.ndarray, left: np.ndarray, right: np.ndarray,
    tri_start: np.ndarray, tri_count: np.ndarray,
    center: np.ndarray, rad: np.ndarray,
    dipole: np.ndarray, quad: np.ndarray,
    beta: float, out: np.ndarray,
) -> None:
    for k in prange(pts.shape[0]):
        out[k] = _winding_bvh_one(
            pts[k, 0], pts[k, 1], pts[k, 2],
            bvh_soup, left, right, tri_start, tri_count,
            center, rad, dipole, quad, beta,
        )


@njit(cache=True, nogil=True)
def _winding_bvh_batch_serial_njit(
    pts: np.ndarray,
    bvh_soup: np.ndarray, left: np.ndarray, right: np.ndarray,
    tri_start: np.ndarray, tri_count: np.ndarray,
    center: np.ndarray, rad: np.ndarray,
    dipole: np.ndarray, quad: np.ndarray,
    beta: float, out: np.ndarray,
) -> None:
    for k in range(pts.shape[0]):
        out[k] = _winding_bvh_one(
            pts[k, 0], pts[k, 1], pts[k, 2],
            bvh_soup, left, right, tri_start, tri_count,
            center, rad, dipole, quad, beta,
        )


@njit(cache=True, nogil=True)
def _exact_solid_angles_subset_njit(
    pts: np.ndarray, bvh_soup: np.ndarray, out: np.ndarray, indices: np.ndarray
) -> None:
    """Exact Van Oosterom-Strackee solid angle ONLY for uncertain ``indices``."""
    inv_4pi = 1.0 / (4.0 * np.pi)
    m_tot = bvh_soup.shape[0]
    for idx in range(indices.shape[0]):
        k = indices[idx]
        px, py, pz = pts[k, 0], pts[k, 1], pts[k, 2]
        s = 0.0
        for m in range(m_tot):
            ax = bvh_soup[m, 0, 0] - px
            ay = bvh_soup[m, 0, 1] - py
            az = bvh_soup[m, 0, 2] - pz
            bx = bvh_soup[m, 1, 0] - px
            by = bvh_soup[m, 1, 1] - py
            bz = bvh_soup[m, 1, 2] - pz
            cx = bvh_soup[m, 2, 0] - px
            cy = bvh_soup[m, 2, 1] - py
            cz = bvh_soup[m, 2, 2] - pz
            la = np.sqrt(ax * ax + ay * ay + az * az)
            lb = np.sqrt(bx * bx + by * by + bz * bz)
            lc = np.sqrt(cx * cx + cy * cy + cz * cz)
            num = (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz)
                   + az * (bx * cy - by * cx))
            den = (la * lb * lc + (ax * bx + ay * by + az * bz) * lc
                   + (bx * cx + by * cy + bz * cz) * la
                   + (cx * ax + cy * ay + cz * az) * lb)
            s += 2.0 * np.arctan2(num, den)
        out[k] = s * inv_4pi


class TriSoupAccel:
    """Compiled Numba spatial grid + BVH winding accelerator for a triangle mesh."""

    def __init__(self, nodes: np.ndarray, tris: np.ndarray, *, leaf_size: int = 8) -> None:
        self.nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        self.tris = np.ascontiguousarray(tris, dtype=np.int64)
        self.n_tris = int(self.tris.shape[0])
        self.fingerprint = int(_mesh_fingerprint_njit(self.nodes, self.tris))

        if self.n_tris == 0:
            self.soup = np.zeros((0, 3, 3), dtype=np.float64)
            self.tmin = np.zeros((0, 3), dtype=np.float64)
            self.tmax = np.zeros((0, 3), dtype=np.float64)
            self.area_vec = np.zeros((0, 3), dtype=np.float64)
            self.centroid = np.zeros((0, 3), dtype=np.float64)
            self.area_abs = np.zeros(0, dtype=np.float64)
            self.glo = np.zeros(3, dtype=np.float64)
            self.inv_h = np.ones(3, dtype=np.float64)
            self.gdims = (1, 1, 1)
            self.offsets = np.zeros(2, dtype=np.int64)
            self.tri_ids = np.zeros(0, dtype=np.int64)
            self.bvh_soup = self.soup
            self.bvh_left = np.full(1, -1, dtype=np.int64)
            self.bvh_right = np.full(1, -1, dtype=np.int64)
            self.bvh_start = np.zeros(1, dtype=np.int64)
            self.bvh_count = np.zeros(1, dtype=np.int64)
            self.bvh_center = np.zeros((1, 3), dtype=np.float64)
            self.bvh_rad = np.zeros(1, dtype=np.float64)
            self.bvh_dipole = np.zeros((1, 3), dtype=np.float64)
            self.bvh_quad = np.zeros((1, 9), dtype=np.float64)
        else:
            (self.soup, self.tmin, self.tmax, self.area_vec,
             self.centroid, self.area_abs) = _pack_tris_njit(self.nodes, self.tris)
            glo = self.tmin.min(axis=0)
            ghi = self.tmax.max(axis=0)
            span = np.maximum(ghi - glo, 1e-12)
            # Balance voxel count with triangle count: ~(4 * M)**(1/3), clamped [4, 96]
            target = max(4, min(96, int(np.ceil((4.0 * self.n_tris) ** (1.0 / 3.0)))))
            max_span = float(span.max())
            dims = np.maximum(1, np.minimum(96, np.round(target * (span / max_span)).astype(np.int64)))
            self.gdims = (int(dims[0]), int(dims[1]), int(dims[2]))
            pad = 1e-9 * span
            self.glo = glo - pad
            gspan = span + 2.0 * pad
            self.inv_h = dims.astype(np.float64) / gspan
            self.offsets, self.tri_ids = _build_grid_csr_njit(
                self.tmin, self.tmax, self.glo, self.inv_h, *self.gdims
            )
            (self.bvh_soup, self.bvh_left, self.bvh_right,
             self.bvh_start, self.bvh_count,
             self.bvh_center, self.bvh_rad,
             self.bvh_dipole, self.bvh_quad) = _build_bvh_njit(
                self.soup, self.tmin, self.tmax,
                self.area_vec, self.centroid, self.area_abs,
                leaf_size=int(leaf_size),
            )

        # Thread-local scratch buffers so concurrent workers share TriSoupAccel safely
        self._thread_local = threading.local()

    def _get_query_state(self) -> list:
        """Return per-thread [stamp_buf, stamp_val, out_buf] scratch buffers."""
        state = getattr(self._thread_local, "state", None)
        if state is None:
            state = [
                np.zeros(self.n_tris, dtype=np.int32),
                0,
                np.empty(min(max(256, self.n_tris), 16384), dtype=np.int64),
            ]
            self._thread_local.state = state
        return state

    def query_box_ids(self, blo: np.ndarray, bhi: np.ndarray) -> np.ndarray:
        """Return deduplicated triangle indices whose AABB overlaps ``[blo, bhi]``."""
        if self.n_tris == 0:
            return np.zeros(0, dtype=np.int64)
        state = self._get_query_state()
        stamp_buf, stamp_val, out_buf = state
        stamp_val += 1
        if stamp_val >= 2_000_000_000:
            stamp_buf.fill(0)
            stamp_val = 1
        blo_a = np.asarray(blo, dtype=np.float64)
        bhi_a = np.asarray(bhi, dtype=np.float64)
        n = query_box_candidates_njit(
            blo_a, bhi_a, self.tmin, self.tmax,
            self.glo, self.inv_h, *self.gdims,
            self.offsets, self.tri_ids,
            stamp_buf, stamp_val, out_buf,
        )
        if n > out_buf.shape[0]:
            out_buf = np.empty(n, dtype=np.int64)
            state[2] = out_buf
            stamp_val += 1
            n = query_box_candidates_njit(
                blo_a, bhi_a, self.tmin, self.tmax,
                self.glo, self.inv_h, *self.gdims,
                self.offsets, self.tri_ids,
                stamp_buf, stamp_val, out_buf,
            )
        state[1] = stamp_val
        return out_buf[:n].copy()

    def has_box_aabb_overlap(self, blo: np.ndarray, bhi: np.ndarray) -> bool:
        """Fast check whether any triangle AABB overlaps ``[blo, bhi]``."""
        if self.n_tris == 0:
            return False
        blo_a = np.asarray(blo, dtype=np.float64)
        bhi_a = np.asarray(bhi, dtype=np.float64)
        return _has_box_aabb_overlap_njit(
            blo_a, bhi_a, self.tmin, self.tmax,
            self.glo, self.inv_h, *self.gdims,
            self.offsets, self.tri_ids,
        )

    def has_boxes_aabb_overlap(self, blos: np.ndarray, h: np.ndarray) -> np.ndarray:
        """Vectorized/parallel check for N boxes ``[blos[m], blos[m] + h]``."""
        blos_a = np.ascontiguousarray(blos, dtype=np.float64).reshape(-1, 3)
        h_a = np.ascontiguousarray(h, dtype=np.float64)
        out = np.zeros(blos_a.shape[0], dtype=np.bool_)
        if self.n_tris == 0 or blos_a.shape[0] == 0:
            return out
        _has_boxes_aabb_overlap_batch_njit(
            blos_a, h_a, self.tmin, self.tmax,
            self.glo, self.inv_h, *self.gdims,
            self.offsets, self.tri_ids, out,
        )
        return out

    def winding_values(self, pts: np.ndarray, *, beta: float = 2.6,
                       refine_thresh: float | None = None,
                       refine_band: float = 0.06,
                       parallel: bool = True) -> np.ndarray:
        """Evaluate winding numbers via BVH (with optional exact near-threshold fallback)."""
        pts_a = np.ascontiguousarray(pts, dtype=np.float64).reshape(-1, 3)
        n_pts = pts_a.shape[0]
        if n_pts == 0 or self.n_tris == 0:
            return np.zeros(n_pts, dtype=np.float64)
        out = np.empty(n_pts, dtype=np.float64)
        if self.n_tris <= 2000:
            # Tiny mesh: exact evaluation is already sub-millisecond
            idx = np.arange(n_pts, dtype=np.int64)
            _exact_solid_angles_subset_njit(pts_a, self.bvh_soup, out, idx)
            return out
        if parallel and n_pts >= 32:
            _winding_bvh_batch_parallel_njit(
                pts_a, self.bvh_soup, self.bvh_left, self.bvh_right,
                self.bvh_start, self.bvh_count,
                self.bvh_center, self.bvh_rad, self.bvh_dipole, self.bvh_quad,
                float(beta), out,
            )
        else:
            _winding_bvh_batch_serial_njit(
                pts_a, self.bvh_soup, self.bvh_left, self.bvh_right,
                self.bvh_start, self.bvh_count,
                self.bvh_center, self.bvh_rad, self.bvh_dipole, self.bvh_quad,
                float(beta), out,
            )
        if refine_thresh is not None:
            near = np.flatnonzero(np.abs(out - float(refine_thresh)) < float(refine_band))
            if len(near):
                _exact_solid_angles_subset_njit(pts_a, self.bvh_soup, out, near)
        return out

    def winding_inside(self, pts: np.ndarray, thresh: float = 0.5, *,
                       parallel: bool = True) -> np.ndarray:
        """Inside boolean flags using BVH winding + exact refinement near ``thresh``."""
        vals = self.winding_values(
            pts, beta=2.6, refine_thresh=thresh, refine_band=0.06,
            parallel=parallel,
        )
        return vals > float(thresh)


@njit(cache=True, nogil=True)
def _has_box_aabb_overlap_njit(blo: np.ndarray, bhi: np.ndarray,
                               tmin: np.ndarray, tmax: np.ndarray,
                               glo: np.ndarray, inv_h: np.ndarray,
                               gx: int, gy: int, gz: int,
                               offsets: np.ndarray, tri_ids: np.ndarray) -> bool:
    """Short-circuit test: does any triangle AABB overlap ``[blo, bhi]``?"""
    x0, x1 = _tri_bin_bounds(blo[0], bhi[0], glo[0], inv_h[0], gx - 1)
    y0, y1 = _tri_bin_bounds(blo[1], bhi[1], glo[1], inv_h[1], gy - 1)
    z0, z1 = _tri_bin_bounds(blo[2], bhi[2], glo[2], inv_h[2], gz - 1)
    ax0, ay0, az0 = blo[0], blo[1], blo[2]
    ax1, ay1, az1 = bhi[0], bhi[1], bhi[2]
    for kz in range(z0, z1 + 1):
        base_z = kz * gy
        for jy in range(y0, y1 + 1):
            base_y = (base_z + jy) * gx
            for ix in range(x0, x1 + 1):
                v = base_y + ix
                s = offsets[v]
                e = offsets[v + 1]
                for p in range(s, e):
                    tid = tri_ids[p]
                    if (tmax[tid, 0] >= ax0 and tmin[tid, 0] <= ax1
                            and tmax[tid, 1] >= ay0 and tmin[tid, 1] <= ay1
                            and tmax[tid, 2] >= az0 and tmin[tid, 2] <= az1):
                        return True
    return False


@njit(cache=True, parallel=True)
def _has_boxes_aabb_overlap_batch_njit(blos: np.ndarray, h: np.ndarray,
                                       tmin: np.ndarray, tmax: np.ndarray,
                                       glo: np.ndarray, inv_h: np.ndarray,
                                       gx: int, gy: int, gz: int,
                                       offsets: np.ndarray, tri_ids: np.ndarray,
                                       out: np.ndarray) -> None:
    """Parallel short-circuit AABB overlap check for ``blos`` of shape (M, 3)."""
    n_boxes = blos.shape[0]
    hx, hy, hz = h[0], h[1], h[2]
    for idx in prange(n_boxes):
        ax0 = blos[idx, 0]
        ay0 = blos[idx, 1]
        az0 = blos[idx, 2]
        ax1 = ax0 + hx
        ay1 = ay0 + hy
        az1 = az0 + hz
        x0, x1 = _tri_bin_bounds(ax0, ax1, glo[0], inv_h[0], gx - 1)
        y0, y1 = _tri_bin_bounds(ay0, ay1, glo[1], inv_h[1], gy - 1)
        z0, z1 = _tri_bin_bounds(az0, az1, glo[2], inv_h[2], gz - 1)
        hit = False
        for kz in range(z0, z1 + 1):
            base_z = kz * gy
            for jy in range(y0, y1 + 1):
                base_y = (base_z + jy) * gx
                for ix in range(x0, x1 + 1):
                    v = base_y + ix
                    s = offsets[v]
                    e = offsets[v + 1]
                    for p in range(s, e):
                        tid = tri_ids[p]
                        if (tmax[tid, 0] >= ax0 and tmin[tid, 0] <= ax1
                                and tmax[tid, 1] >= ay0 and tmin[tid, 1] <= ay1
                                and tmax[tid, 2] >= az0 and tmin[tid, 2] <= az1):
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    break
            if hit:
                break
        out[idx] = hit


_ACCEL_CACHE: OrderedDict[tuple[int, int, tuple[int, ...], tuple[int, ...]], TriSoupAccel] = OrderedDict()
_ACCEL_LOCK = threading.Lock()
_ACCEL_MAX = 8


def get_accel(nodes: np.ndarray, tris: np.ndarray) -> TriSoupAccel:
    """Return a cached, mutation-checked :class:`TriSoupAccel` for ``(nodes, tris)``."""
    nodes_a = np.ascontiguousarray(nodes, dtype=np.float64)
    tris_a = np.ascontiguousarray(tris, dtype=np.int64)
    fp = int(_mesh_fingerprint_njit(nodes_a, tris_a))
    key = (fp, int(nodes_a.shape[0]), tuple(int(s) for s in nodes_a.shape),
           tuple(int(s) for s in tris_a.shape))
    with _ACCEL_LOCK:
        hit = _ACCEL_CACHE.get(key)
        if hit is not None:
            _ACCEL_CACHE.move_to_end(key)
            return hit
    acc = TriSoupAccel(nodes_a, tris_a)
    with _ACCEL_LOCK:
        _ACCEL_CACHE[key] = acc
        if len(_ACCEL_CACHE) > _ACCEL_MAX:
            _ACCEL_CACHE.popitem(last=False)
    return acc
