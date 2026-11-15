"""S3 classification: triangle→cell SAT binning + vertex ray parity (v2 §7.3).

Cell states: 0 OUT, 1 IN, 2 CUT (matches ``cell`` table, v2 §4.2).
``eta`` is left to S4 (exact volume); this module only sets states.

Numba-accelerated implementation:
- ``clip_poly_to_box`` delegates to a stack-allocated ``@njit`` kernel.
- ``find_cut_cells`` performs a race-free two-pass ``@njit(parallel=True)`` SAT
  sweep with zero Python-per-triangle loops.
- ``classify_vertices`` bins triangles onto uniform grid lines in ``O(M)`` via
  ``@njit`` CSR construction and runs all scanlines in ``prange``, resolving
  axis disagreements with :class:`wht_iga.geometry.accel.TriSoupAccel` BVH winding.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange  # type: ignore[import-untyped]

from .accel import TriSoupAccel, get_accel

OUT, IN, CUT = 0, 1, 2


@njit(cache=True, nogil=True)
def _tri_box_overlap(t0: np.ndarray, t1: np.ndarray, t2: np.ndarray,
                     blo: np.ndarray, bhi: np.ndarray) -> bool:
    """Akenine-Moller triangle/AABB SAT: 1 face + 9 edge-cross + 3 axis tests."""
    c0 = 0.5 * (blo[0] + bhi[0])
    c1 = 0.5 * (blo[1] + bhi[1])
    c2 = 0.5 * (blo[2] + bhi[2])
    e0 = 0.5 * (bhi[0] - blo[0])
    e1 = 0.5 * (bhi[1] - blo[1])
    e2 = 0.5 * (bhi[2] - blo[2])

    v0x, v0y, v0z = t0[0] - c0, t0[1] - c1, t0[2] - c2
    v1x, v1y, v1z = t1[0] - c0, t1[1] - c1, t1[2] - c2
    v2x, v2y, v2z = t2[0] - c0, t2[1] - c1, t2[2] - c2

    # 3 box axes
    if min(v0x, v1x, v2x) > e0 or max(v0x, v1x, v2x) < -e0:
        return False
    if min(v0y, v1y, v2y) > e1 or max(v0y, v1y, v2y) < -e1:
        return False
    if min(v0z, v1z, v2z) > e2 or max(v0z, v1z, v2z) < -e2:
        return False

    # Triangle edges
    f0x, f0y, f0z = v1x - v0x, v1y - v0y, v1z - v0z
    f1x, f1y, f1z = v2x - v1x, v2y - v1y, v2z - v1z
    f2x, f2y, f2z = v0x - v2x, v0y - v2y, v0z - v2z

    # Triangle face plane
    nx = f0y * (v2z - v0z) - f0z * (v2y - v0y)
    ny = f0z * (v2x - v0x) - f0x * (v2z - v0z)
    nz = f0x * (v2y - v0y) - f0y * (v2x - v0x)
    d = -(nx * v0x + ny * v0y + nz * v0z)
    r = e0 * abs(nx) + e1 * abs(ny) + e2 * abs(nz)
    if abs(d) > r + 1e-12 * (r + 1.0):
        return False

    # 9 edge-cross-axis tests (unrolled allocation-free)
    for fx, fy, fz in ((f0x, f0y, f0z), (f1x, f1y, f1z), (f2x, f2y, f2z)):
        # ax = ( -fy, fx, 0 )
        p0 = -fy * v0x + fx * v0y
        p1 = -fy * v1x + fx * v1y
        p2 = -fy * v2x + fx * v2y
        r = e0 * abs(fy) + e1 * abs(fx)
        if min(p0, p1, p2) > r or max(p0, p1, p2) < -r:
            return False
        # ax = ( 0, -fz, fy )
        p0 = -fz * v0y + fy * v0z
        p1 = -fz * v1y + fy * v1z
        p2 = -fz * v2y + fy * v2z
        r = e1 * abs(fz) + e2 * abs(fy)
        if min(p0, p1, p2) > r or max(p0, p1, p2) < -r:
            return False
        # ax = ( fz, 0, -fx )
        p0 = fz * v0x - fx * v0z
        p1 = fz * v1x - fx * v1z
        p2 = fz * v2x - fx * v2z
        r = e2 * abs(fx) + e0 * abs(fz)
        if min(p0, p1, p2) > r or max(p0, p1, p2) < -r:
            return False

    return True


@njit(cache=True, parallel=True)
def _sat_candidates(tris: np.ndarray, tri_idx: np.ndarray,
                    cell_lo: np.ndarray, cell_hi: np.ndarray,
                    out: np.ndarray) -> None:
    """Evaluate SAT for pre-binned (triangle, cell) candidate pairs."""
    for k in prange(len(tri_idx)):
        t = tris[tri_idx[k]]
        if _tri_box_overlap(t[0], t[1], t[2], cell_lo[k], cell_hi[k]):
            out[k] = True


@njit(cache=True, nogil=True)
def _clip_poly_to_box_njit(poly: np.ndarray, n_in: int,
                           blo: np.ndarray, bhi: np.ndarray,
                           out_buf: np.ndarray) -> int:
    """Sutherland–Hodgman clip of a 3D polygon vs AABB into ``out_buf`` (cap >= 32)."""
    if n_in <= 0:
        return 0
    work = np.empty((32, 3), dtype=np.float64)
    tmp = np.empty((32, 3), dtype=np.float64)
    n = n_in
    for k in range(n_in):
        work[k, 0] = poly[k, 0]
        work[k, 1] = poly[k, 1]
        work[k, 2] = poly[k, 2]

    for ax in range(3):
        for pass_id in range(2):
            if n == 0:
                return 0
            side = 1.0 if pass_id == 0 else -1.0
            bound = -blo[ax] if pass_id == 0 else bhi[ax]
            sx, sy, sz = work[n - 1, 0], work[n - 1, 1], work[n - 1, 2]
            ds = side * work[n - 1, ax] + bound
            m = 0
            for i in range(n):
                ex, ey, ez = work[i, 0], work[i, 1], work[i, 2]
                de = side * work[i, ax] + bound
                ie = de >= -1e-12
                ins = ds >= -1e-12
                if ie:
                    if not ins:
                        t = ds / (ds - de)
                        tmp[m, 0] = sx + t * (ex - sx)
                        tmp[m, 1] = sy + t * (ey - sy)
                        tmp[m, 2] = sz + t * (ez - sz)
                        m += 1
                    tmp[m, 0] = ex
                    tmp[m, 1] = ey
                    tmp[m, 2] = ez
                    m += 1
                elif ins:
                    t = ds / (ds - de)
                    tmp[m, 0] = sx + t * (ex - sx)
                    tmp[m, 1] = sy + t * (ey - sy)
                    tmp[m, 2] = sz + t * (ez - sz)
                    m += 1
                sx, sy, sz, ds = ex, ey, ez, de
            for i in range(m):
                work[i, 0] = tmp[i, 0]
                work[i, 1] = tmp[i, 1]
                work[i, 2] = tmp[i, 2]
            n = m

    for i in range(n):
        out_buf[i, 0] = work[i, 0]
        out_buf[i, 1] = work[i, 1]
        out_buf[i, 2] = work[i, 2]
    return n


def clip_poly_to_box(poly: np.ndarray, blo: np.ndarray,
                     bhi: np.ndarray) -> np.ndarray:
    """Sutherland–Hodgman clip of a 3D polygon vs AABB (convex output)."""
    p_arr = np.ascontiguousarray(poly, dtype=np.float64)
    if p_arr.shape[0] == 0:
        return np.zeros((0, 3), dtype=np.float64)
    out = np.empty((32, 3), dtype=np.float64)
    n = _clip_poly_to_box_njit(
        p_arr, p_arr.shape[0],
        np.ascontiguousarray(blo, dtype=np.float64),
        np.ascontiguousarray(bhi, dtype=np.float64),
        out,
    )
    return out[:n].copy() if n > 0 else np.zeros((0, 3), dtype=np.float64)


@njit(cache=True, parallel=True)
def _find_cut_cells_njit(soup: np.ndarray, tmin: np.ndarray, tmax: np.ndarray,
                         lo: np.ndarray, h: np.ndarray,
                         nx: int, ny: int, nz: int) -> np.ndarray:
    """Two-pass parallel SAT test mapping all overlapping triangle-cell pairs."""
    m = soup.shape[0]
    counts = np.zeros(m, dtype=np.int64)
    lo0, lo1, lo2 = lo[0], lo[1], lo[2]
    inv0, inv1, inv2 = 1.0 / h[0], 1.0 / h[1], 1.0 / h[2]
    nx1, ny1, nz1 = nx - 1, ny - 1, nz - 1

    for t in prange(m):
        i0 = int(np.floor((tmin[t, 0] - lo0) * inv0))
        i1 = int(np.floor((tmax[t, 0] - lo0) * inv0))
        j0 = int(np.floor((tmin[t, 1] - lo1) * inv1))
        j1 = int(np.floor((tmax[t, 1] - lo1) * inv1))
        k0 = int(np.floor((tmin[t, 2] - lo2) * inv2))
        k1 = int(np.floor((tmax[t, 2] - lo2) * inv2))
        if i0 < 0:
            i0 = 0
        elif i0 > nx1:
            i0 = nx1
        if i1 < 0:
            i1 = 0
        elif i1 > nx1:
            i1 = nx1
        if j0 < 0:
            j0 = 0
        elif j0 > ny1:
            j0 = ny1
        if j1 < 0:
            j1 = 0
        elif j1 > ny1:
            j1 = ny1
        if k0 < 0:
            k0 = 0
        elif k0 > nz1:
            k0 = nz1
        if k1 < 0:
            k1 = 0
        elif k1 > nz1:
            k1 = nz1
        counts[t] = (i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1)

    off = np.empty(m + 1, dtype=np.int64)
    off[0] = 0
    for t in range(m):
        off[t + 1] = off[t] + counts[t]

    hits = np.full(off[m], -1, dtype=np.int64)
    hx, hy, hz = h[0], h[1], h[2]

    for t in prange(m):
        i0 = int(np.floor((tmin[t, 0] - lo0) * inv0))
        i1 = int(np.floor((tmax[t, 0] - lo0) * inv0))
        j0 = int(np.floor((tmin[t, 1] - lo1) * inv1))
        j1 = int(np.floor((tmax[t, 1] - lo1) * inv1))
        k0 = int(np.floor((tmin[t, 2] - lo2) * inv2))
        k1 = int(np.floor((tmax[t, 2] - lo2) * inv2))
        if i0 < 0:
            i0 = 0
        elif i0 > nx1:
            i0 = nx1
        if i1 < 0:
            i1 = 0
        elif i1 > nx1:
            i1 = nx1
        if j0 < 0:
            j0 = 0
        elif j0 > ny1:
            j0 = ny1
        if j1 < 0:
            j1 = 0
        elif j1 > ny1:
            j1 = ny1
        if k0 < 0:
            k0 = 0
        elif k0 > nz1:
            k0 = nz1
        if k1 < 0:
            k1 = 0
        elif k1 > nz1:
            k1 = nz1

        t0 = soup[t, 0]
        t1 = soup[t, 1]
        t2 = soup[t, 2]
        blo_loc = np.empty(3, dtype=np.float64)
        bhi_loc = np.empty(3, dtype=np.float64)
        pos = off[t]
        for i in range(i0, i1 + 1):
            blo_loc[0] = lo0 + hx * i
            bhi_loc[0] = blo_loc[0] + hx
            for j in range(j0, j1 + 1):
                blo_loc[1] = lo1 + hy * j
                bhi_loc[1] = blo_loc[1] + hy
                for k in range(k0, k1 + 1):
                    blo_loc[2] = lo2 + hz * k
                    bhi_loc[2] = blo_loc[2] + hz
                    if _tri_box_overlap(t0, t1, t2, blo_loc, bhi_loc):
                        hits[pos] = i + nx * (j + ny * k)
                    pos += 1

    cut = np.zeros(nx * ny * nz, dtype=np.bool_)
    for q in range(hits.shape[0]):
        cidx = hits[q]
        if cidx >= 0:
            cut[cidx] = True
    return cut


def find_cut_cells(nodes: np.ndarray, tris: np.ndarray,
                   lo: np.ndarray, h: np.ndarray, n: tuple[int, int, int],
                   *, accel: TriSoupAccel | None = None) -> np.ndarray:
    """Cells overlapped by any triangle → CUT candidates (bool, C-order flat)."""
    n = tuple(int(v) for v in n)
    if len(tris) == 0:
        return np.zeros(int(np.prod(n)), dtype=bool)
    acc = get_accel(nodes, tris) if accel is None else accel
    return _find_cut_cells_njit(
        acc.soup, acc.tmin, acc.tmax,
        np.ascontiguousarray(lo, dtype=np.float64),
        np.ascontiguousarray(h, dtype=np.float64),
        int(n[0]), int(n[1]), int(n[2]),
    )


@njit(cache=True)
def _ray_hits(orig: np.ndarray, direc: np.ndarray, tris: np.ndarray,
              out: np.ndarray) -> int:
    """Moller-Trumbore hits of a ray against triangles."""
    n = 0
    eps = 1e-12
    for m in range(tris.shape[0]):
        v0 = tris[m, 0]
        e1 = tris[m, 1] - v0
        e2 = tris[m, 2] - v0
        p = np.cross(direc, e2)
        det = np.dot(e1, p)
        if -eps < det < eps:
            continue
        inv = 1.0 / det
        tv = orig - v0
        u = np.dot(tv, p) * inv
        if u < -1e-9 or u > 1.0 + 1e-9:
            continue
        q = np.cross(tv, e1)
        v = np.dot(direc, q) * inv
        if v < -1e-9 or u + v > 1.0 + 1e-9:
            continue
        t = np.dot(e2, q) * inv
        if t > 1e-9:
            out[n] = t
            n += 1
    return n


@njit(cache=True)
def _ray_hits_axis_csr(ox: float, oy: float, oz: float, axis: int,
                       soup: np.ndarray, tri_ids: np.ndarray, s: int, e: int,
                       out: np.ndarray) -> int:
    """Allocation-free Moller-Trumbore hits of an axis-aligned ray against CSR slice."""
    n = 0
    eps = 1e-12
    dx = 1.0 if axis == 0 else 0.0
    dy = 1.0 if axis == 1 else 0.0
    dz = 1.0 if axis == 2 else 0.0
    for idx in range(s, e):
        m = tri_ids[idx]
        v0x = soup[m, 0, 0]
        v0y = soup[m, 0, 1]
        v0z = soup[m, 0, 2]
        e1x = soup[m, 1, 0] - v0x
        e1y = soup[m, 1, 1] - v0y
        e1z = soup[m, 1, 2] - v0z
        e2x = soup[m, 2, 0] - v0x
        e2y = soup[m, 2, 1] - v0y
        e2z = soup[m, 2, 2] - v0z
        px = dy * e2z - dz * e2y
        py = dz * e2x - dx * e2z
        pz = dx * e2y - dy * e2x
        det = e1x * px + e1y * py + e1z * pz
        if -eps < det < eps:
            continue
        inv = 1.0 / det
        tvx = ox - v0x
        tvy = oy - v0y
        tvz = oz - v0z
        u = (tvx * px + tvy * py + tvz * pz) * inv
        if u < -1e-9 or u > 1.0 + 1e-9:
            continue
        qx = tvy * e1z - tvz * e1y
        qy = tvz * e1x - tvx * e1z
        qz = tvx * e1y - tvy * e1x
        v = (dx * qx + dy * qy + dz * qz) * inv
        if v < -1e-9 or u + v > 1.0 + 1e-9:
            continue
        t = (e2x * qx + e2y * qy + e2z * qz) * inv
        if t > 1e-9:
            out[n] = t
            n += 1
    return n


@njit(cache=True)
def _parity_scan(xs: np.ndarray, hits: np.ndarray, nh: int,
                 out: np.ndarray, merge_eps: float) -> None:
    """Inside flags: odd hit-cluster count before x → inside.

    Hits within ``merge_eps`` collapse to one: a ray through a shared
    mesh edge/vertex meets several triangles at the same t but crosses
    the surface once (topological dedup).
    """
    hs = np.sort(hits[:nh])
    m = 0
    for k in range(nh):
        if k == 0 or hs[k] - hs[k - 1] > merge_eps:
            hs[m] = hs[k]
            m += 1
    c = 0
    for i in range(xs.shape[0]):
        while c < m and hs[c] < xs[i]:
            c += 1
        out[i] = (c % 2) == 1


@njit(cache=True)
def _build_line_csr_njit(tmin: np.ndarray, tmax: np.ndarray,
                         o0: int, o1: int,
                         lo0: float, h0: float, n0: int,
                         lo1: float, h1: float, n1: int,
                         marg: float):
    """Bin triangles onto (n0+1) x (n1+1) axis scanlines in O(M + total_bins)."""
    m = tmin.shape[0]
    n_lines = (n0 + 1) * (n1 + 1)
    counts = np.zeros(n_lines, dtype=np.int64)
    inv0 = 1.0 / h0 if h0 > 0.0 else 0.0
    inv1 = 1.0 / h1 if h1 > 0.0 else 0.0

    for t in range(m):
        a0 = int(np.ceil((tmin[t, o0] - marg - lo0) * inv0))
        b0 = int(np.floor((tmax[t, o0] + marg - lo0) * inv0))
        a1 = int(np.ceil((tmin[t, o1] - marg - lo1) * inv1))
        b1 = int(np.floor((tmax[t, o1] + marg - lo1) * inv1))
        if a0 < 0:
            a0 = 0
        if b0 > n0:
            b0 = n0
        if a1 < 0:
            a1 = 0
        if b1 > n1:
            b1 = n1
        if a0 <= b0 and a1 <= b1:
            for i0 in range(a0, b0 + 1):
                base = i0 * (n1 + 1)
                for i1 in range(a1, b1 + 1):
                    counts[base + i1] += 1

    offsets = np.empty(n_lines + 1, dtype=np.int64)
    offsets[0] = 0
    for r in range(n_lines):
        offsets[r + 1] = offsets[r] + counts[r]

    tri_ids = np.empty(offsets[n_lines], dtype=np.int64)
    cursor = offsets[:-1].copy()

    for t in range(m):
        a0 = int(np.ceil((tmin[t, o0] - marg - lo0) * inv0))
        b0 = int(np.floor((tmax[t, o0] + marg - lo0) * inv0))
        a1 = int(np.ceil((tmin[t, o1] - marg - lo1) * inv1))
        b1 = int(np.floor((tmax[t, o1] + marg - lo1) * inv1))
        if a0 < 0:
            a0 = 0
        if b0 > n0:
            b0 = n0
        if a1 < 0:
            a1 = 0
        if b1 > n1:
            b1 = n1
        if a0 <= b0 and a1 <= b1:
            for i0 in range(a0, b0 + 1):
                base = i0 * (n1 + 1)
                for i1 in range(a1, b1 + 1):
                    r = base + i1
                    tri_ids[cursor[r]] = t
                    cursor[r] += 1

    return offsets, tri_ids


@njit(cache=True, parallel=True)
def _axis_votes_grid_njit(soup: np.ndarray, tmin: np.ndarray, tmax: np.ndarray,
                          lo: np.ndarray, h: np.ndarray,
                          nx: int, ny: int, nz: int,
                          axis: int, marg: float) -> np.ndarray:
    """Compute scanline parity votes on a (nx+1, ny+1, nz+1) grid along ``axis``."""
    dims = (nx, ny, nz)
    if axis == 0:
        o0, o1 = 1, 2
    elif axis == 1:
        o0, o1 = 0, 2
    else:
        o0, o1 = 0, 1

    na = dims[axis]
    n0 = dims[o0]
    n1 = dims[o1]
    offsets, tri_ids = _build_line_csr_njit(
        tmin, tmax, o0, o1, lo[o0], h[o0], n0, lo[o1], h[o1], n1, marg
    )
    xs_rel = np.empty(na + 1, dtype=np.float64)
    ha = h[axis]
    org_a = lo[axis] - 1e-6
    for q in range(na + 1):
        xs_rel[q] = (lo[axis] + q * ha) - org_a
    span = h[axis] * na + 1.0
    merge_eps = 1e-9 * span

    votes = np.zeros((nx + 1, ny + 1, nz + 1), dtype=np.bool_)
    n_lines = (n0 + 1) * (n1 + 1)

    for r in prange(n_lines):
        i0 = r // (n1 + 1)
        i1 = r % (n1 + 1)
        s = offsets[r]
        e = offsets[r + 1]
        nh_cap = e - s
        if nh_cap == 0:
            continue
        p0 = lo[o0] + i0 * h[o0]
        p1 = lo[o1] + i1 * h[o1]
        if axis == 0:
            ox, oy, oz = org_a, p0, p1
        elif axis == 1:
            ox, oy, oz = p0, org_a, p1
        else:
            ox, oy, oz = p0, p1, org_a

        hits = np.empty(nh_cap * 2 + 8, dtype=np.float64)
        nh = _ray_hits_axis_csr(ox, oy, oz, axis, soup, tri_ids, s, e, hits)
        if nh == 0:
            continue
        f = np.empty(na + 1, dtype=np.bool_)
        _parity_scan(xs_rel, hits, nh, f, merge_eps)
        for q in range(na + 1):
            if axis == 0:
                votes[q, i0, i1] = f[q]
            elif axis == 1:
                votes[i0, q, i1] = f[q]
            else:
                votes[i0, i1, q] = f[q]

    return votes.ravel()


def _axis_votes(verts: np.ndarray, soup: np.ndarray, axis: int) -> np.ndarray:
    """Scanline parity votes along one axis; returns bool per vertex."""
    nv = verts.shape[0]
    votes = np.zeros(nv, dtype=bool)
    other = [a for a in range(3) if a != axis]
    origin = verts.min(axis=0) - 1e-6
    direc = np.zeros(3)
    direc[axis] = 1.0
    ukeys, inv = np.unique(verts[:, other], axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    order_all = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order_all], np.arange(ukeys.shape[0] + 1))
    hits = np.empty(soup.shape[0] * 2 + 8)
    tmin, tmax = soup.min(axis=1), soup.max(axis=1)
    scale = float(np.max(np.abs(soup))) + 1.0 if soup.shape[0] else 1.0
    marg = 1e-7 * scale
    o0, o1 = other
    for r in range(ukeys.shape[0]):
        idx = order_all[bounds[r]:bounds[r + 1]]
        k0, k1 = ukeys[r, 0], ukeys[r, 1]
        sel = ((tmin[:, o0] <= k0 + marg) & (tmax[:, o0] >= k0 - marg)
               & (tmin[:, o1] <= k1 + marg) & (tmax[:, o1] >= k1 - marg))
        soup_r = np.ascontiguousarray(soup[sel])
        org = origin.copy()
        org[other[0]] = ukeys[r, 0]
        org[other[1]] = ukeys[r, 1]
        nh = _ray_hits(org, direc, soup_r, hits)
        xs_r = np.sort(verts[idx, axis])
        f = np.empty(idx.shape[0], dtype=np.bool_)
        span = xs_r[-1] - xs_r[0] + 1.0
        _parity_scan(xs_r - org[axis], hits, nh, f, 1e-9 * span)
        order = np.argsort(verts[idx, axis])
        for s_idx, ii in enumerate(idx[order]):
            votes[ii] = f[s_idx]
    return votes


@njit(cache=True, inline="always", nogil=True)
def _solid_angle_one(px: float, py: float, pz: float, tris: np.ndarray) -> float:
    """Winding number of one point (Van Oosterom-Strackee), allocation-free."""
    s = 0.0
    for m in range(tris.shape[0]):
        ax = tris[m, 0, 0] - px
        ay = tris[m, 0, 1] - py
        az = tris[m, 0, 2] - pz
        bx = tris[m, 1, 0] - px
        by = tris[m, 1, 1] - py
        bz = tris[m, 1, 2] - pz
        cx = tris[m, 2, 0] - px
        cy = tris[m, 2, 1] - py
        cz = tris[m, 2, 2] - pz
        la = np.sqrt(ax * ax + ay * ay + az * az)
        lb = np.sqrt(bx * bx + by * by + bz * bz)
        lc = np.sqrt(cx * cx + cy * cy + cz * cz)
        num = (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz)
               + az * (bx * cy - by * cx))
        den = (la * lb * lc + (ax * bx + ay * by + az * bz) * lc
               + (bx * cx + by * cy + bz * cz) * la
               + (cx * ax + cy * ay + cz * az) * lb)
        s += 2.0 * np.arctan2(num, den)
    return s / (4.0 * np.pi)


@njit(cache=True, parallel=True)
def _solid_angles(pts: np.ndarray, tris: np.ndarray, out: np.ndarray) -> None:
    """Generalized winding number numerator (Van Oosterom-Strackee).

    Parallel over points; use ``_solid_angles_serial`` for tiny batches.
    """
    for k in prange(pts.shape[0]):
        out[k] = _solid_angle_one(pts[k, 0], pts[k, 1], pts[k, 2], tris)


@njit(cache=True, nogil=True)
def _solid_angles_serial(pts: np.ndarray, tris: np.ndarray, out: np.ndarray) -> None:
    for k in range(pts.shape[0]):
        out[k] = _solid_angle_one(pts[k, 0], pts[k, 1], pts[k, 2], tris)


_PARALLEL_MIN_PTS = 16


def winding_inside(pts: np.ndarray, nodes: np.ndarray, tris: np.ndarray,
                   thresh: float = 0.5, *,
                   accel: TriSoupAccel | None = None,
                   parallel: bool = True) -> np.ndarray:
    """Inside test via generalized winding number (Jacobson 2013 / BVH fast winding)."""
    pts_a = np.ascontiguousarray(pts, dtype=np.float64)
    if pts_a.shape[0] == 0 or len(tris) == 0:
        return np.zeros(pts_a.shape[0], dtype=bool)
    acc = get_accel(nodes, tris) if accel is None else accel
    return acc.winding_inside(pts_a, thresh=thresh, parallel=parallel)


def classify_vertices(nodes: np.ndarray, tris: np.ndarray,
                      lo: np.ndarray, h: np.ndarray,
                      n: tuple[int, int, int], *,
                      accel: TriSoupAccel | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Grid-vertex inside flags with 3-axis majority vote.

    Returns (inside (nv,), disagree (nv,)): vertices where the three
    axis parities disagree are resolved by the winding-number fallback.
    """
    n = tuple(int(v) for v in n)
    nx, ny, nz = n
    nv = (nx + 1) * (ny + 1) * (nz + 1)
    if len(tris) == 0:
        return np.zeros(nv, dtype=bool), np.zeros(nv, dtype=bool)

    acc = get_accel(nodes, tris) if accel is None else accel
    lo_a = np.ascontiguousarray(lo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)
    scale = float(np.max(np.abs(acc.soup))) + 1.0
    marg = 1e-7 * scale

    v0 = _axis_votes_grid_njit(acc.soup, acc.tmin, acc.tmax, lo_a, h_a, nx, ny, nz, 0, marg)
    v1 = _axis_votes_grid_njit(acc.soup, acc.tmin, acc.tmax, lo_a, h_a, nx, ny, nz, 1, marg)
    v2 = _axis_votes_grid_njit(acc.soup, acc.tmin, acc.tmax, lo_a, h_a, nx, ny, nz, 2, marg)

    vsum = v0.astype(np.int8) + v1.astype(np.int8) + v2.astype(np.int8)
    inside = vsum >= 2
    disagree = (vsum == 1) | (vsum == 2)
    if disagree.any():
        idx = np.flatnonzero(disagree)
        nyz = (ny + 1) * (nz + 1)
        ii = idx // nyz
        rem = idx % nyz
        jj = rem // (nz + 1)
        kk = rem % (nz + 1)
        pts_d = np.column_stack([
            lo_a[0] + h_a[0] * ii,
            lo_a[1] + h_a[1] * jj,
            lo_a[2] + h_a[2] * kk,
        ])
        inside[disagree] = acc.winding_inside(pts_d)
    return inside, disagree


def classify_cells(nodes: np.ndarray, tris: np.ndarray,
                   lo: np.ndarray, h: np.ndarray,
                   n: tuple[int, int, int], *,
                   accel: TriSoupAccel | None = None) -> dict[str, np.ndarray]:
    """Cell states from SAT CUT set + vertex votes (eta left to S4).

    Corner test is fully vectorized (8 shifted AND/OR grids): zero
    Python-per-cell overhead on the hot path.
    """
    n = tuple(int(v) for v in n)
    nx, ny, nz = n
    acc = get_accel(nodes, tris) if accel is None else accel
    cut = find_cut_cells(nodes, tris, lo, h, n, accel=acc)
    inside, disagree = classify_vertices(nodes, tris, lo, h, n, accel=acc)
    vgrid = inside.reshape(nx + 1, ny + 1, nz + 1)
    c = [vgrid[i:i + nx, j:j + ny, k:k + nz]
         for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    allin = c[0]
    for g in c[1:]:
        allin = allin & g
    anyin = c[0]
    for g in c[1:]:
        anyin = anyin | g
    state = np.where(cut, CUT, np.where(allin.ravel(), IN,
                                        np.where(anyin.ravel(), CUT, OUT)))
    return {"state": state.astype(np.int64), "cut_sat": cut,
            "vertex_inside": inside, "vertex_disagree": disagree}
