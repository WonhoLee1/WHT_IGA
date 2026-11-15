"""S4 cut quadrature: recursive-divergence moments + NNLS fit (v2 §7.4).

m_ijk = ∫_{Q∩Ω} P_i P_j P_k dV (axis Legendre ≤ q; default q = p+1,
q = 2p for mass). F = e_x∫P dx kills x0/y/z faces; remainder = (a)
clipped-triangle integrals + (b) x = x1 face moments via 2D→1D pieces.
Weights: NNLS + backtracking prune (tables cached, candidates capped);
residual > tol → octree fallback (``QUAD_FALLBACK``).

Numba-accelerated implementation (uses :class:`wht_iga.geometry.accel.TriSoupAccel`
for O(local tris) AABB narrowing and BVH winding, plus compiled ``@njit``
kernels for moments, cross-sections, Vandermonde matrices, octree leaves, and
surface integration).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numba import njit  # type: ignore[import-untyped]
from numpy.polynomial.legendre import leggauss
from scipy.optimize import nnls

from .accel import TriSoupAccel, get_accel
from .classify import (
    _clip_poly_to_box_njit,
    _solid_angles_serial,
    _tri_box_overlap,
    clip_poly_to_box,
    winding_inside,
)
from .issues import IssueLog

_NNLS_TOL = 1e-12  # fit residual vs |cell| (G1 polynomial-accuracy gate)
_NNLS_RELAXED = 1e-6
_NNLS_MAX_CANDIDATES = 2048  # per-cell NNLS time budget (deterministic stride)
_COPLANAR_RTOL = 1e-9


@lru_cache(maxsize=None)
def legendre_monomials(q: int) -> np.ndarray:
    """L[j, a]: coeff of u^a in Legendre L_j on [-1, 1] (3-term recurrence)."""
    L = np.zeros((q + 1, q + 1))
    L[0, 0] = 1.0
    if q >= 1:
        L[1, 1] = 1.0
    for j in range(1, q):
        L[j + 1, 1:] += (2 * j + 1) / (j + 1) * L[j, :-1]
        L[j + 1, :] -= j / (j + 1) * L[j - 1, :]
    return L


def _polyval(C: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Evaluate monomial-coeff rows C[j] at x (broadcast)."""
    P = np.zeros((C.shape[0],) + x.shape)
    for a in np.nonzero((C != 0).any(axis=0))[0]:
        P += C[:, a].reshape((-1,) + (1,) * x.ndim) * x**a
    return P


@lru_cache(maxsize=None)
def _axis_pw_tuple(hx: float, q: int, antideriv: bool) -> np.ndarray:
    """Cached local-coordinate (x0=0) Legendre/antiderivative monomial table."""
    L = legendre_monomials(q)
    return _axis_pw(0.0, float(hx), L, antideriv)


def _axis_pw(x0: float, hx: float, L: np.ndarray, antideriv: bool) -> np.ndarray:
    """Mapped Legendre (or ∫_{x0}^x) as monomials in x (binomial shift)."""
    from math import comb

    q = L.shape[0] - 1
    W = np.zeros((q + 1, q + 2 if antideriv else q + 1))
    for j in range(q + 1):
        c = L[j].copy()
        if antideriv:  # F(u) = ∫L_j with F(-1) = 0
            F = np.append(0.0, c / np.arange(1, q + 2))
            F[0] = -float(np.dot(F, (-1.0) ** np.arange(q + 2)))
            c = F * hx / 2.0
        al, be = 2.0 / hx, -(2.0 * x0 / hx + 1.0)
        for m in np.nonzero(c)[0]:
            for k in range(m + 1):
                W[j, k] += c[m] * comb(m, k) * al**k * be ** (m - k)
    return W


@lru_cache(maxsize=None)
def _gauss1d(n: int) -> tuple[np.ndarray, np.ndarray]:
    z, w = leggauss(n)  # Gauss-Legendre nodes/weights on [0, 1]
    return 0.5 * (z + 1.0), 0.5 * w


@lru_cache(maxsize=None)
def _duffy(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Duffy-mapped tensor rule on ref triangle (exact to 2n-1)."""
    z, w = _gauss1d(n)
    Z1, Z2 = np.meshgrid(z, z, indexing="ij")
    W1, W2 = np.meshgrid(w, w, indexing="ij")
    pts = np.stack([Z1.ravel(), (Z2 * (1 - Z1)).ravel()], axis=1)
    return pts, (W1 * W2 * (1 - Z1)).ravel() * 2.0  # ×2: ref area 1/2


@lru_cache(maxsize=None)
def _tensor_gauss(order: int) -> tuple[np.ndarray, np.ndarray]:
    z, w = _gauss1d(order)  # tensor Gauss-Legendre pts/weights on [0, 1]^3
    g = np.meshgrid(z, z, z, indexing="ij")
    return np.stack(g, axis=-1).reshape(-1, 3), \
        np.einsum("i,j,k->ijk", w, w, w).ravel()


@dataclass
class CellQuadrature:
    """Volume rule for one cell + diagnostics."""

    x: np.ndarray  # (N, 3)
    w: np.ndarray  # (N,) nonnegative
    moments: np.ndarray  # exact (q+1)^3 reference moments
    residual: float = 0.0
    n_candidates: int = 0
    fallback: bool = False
    volume: float = 0.0
    relaxed: bool = False  # accepted at _NNLS_RELAXED, not _NNLS_TOL


@njit(cache=True, inline="always")
def _horner_vec(W: np.ndarray, x: float, out: np.ndarray) -> None:
    """Evaluate all rows W[j] at scalar x using Horner's method into ``out``."""
    r = W.shape[0]
    d = W.shape[1]
    for j in range(r):
        val = W[j, d - 1]
        for a in range(d - 2, -1, -1):
            val = val * x + W[j, a]
        out[j] = val


@njit(cache=True)
def _cell_moments_part_a_njit(tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                              Ax: np.ndarray, Py: np.ndarray, Pz: np.ndarray,
                              DP: np.ndarray, DW: np.ndarray, q: int,
                              coplanar_rtol: float) -> np.ndarray:
    """Compiled part (a) of divergence moments (zero heap allocs inside triangle loop)."""
    M = np.zeros((q + 1, q + 1, q + 1), dtype=np.float64)
    bx, by, bz = blo[0], blo[1], blo[2]
    hx, hy, hz = h[0], h[1], h[2]
    hmax = hx if hx >= hy else hy
    if hz > hmax:
        hmax = hz
    cop_tol = coplanar_rtol * hmax

    b_lo = np.zeros(3, dtype=np.float64)
    b_hi = np.array([hx, hy, hz], dtype=np.float64)
    t_loc = np.empty((3, 3), dtype=np.float64)
    poly_buf = np.empty((32, 3), dtype=np.float64)
    va = np.empty(q + 1, dtype=np.float64)
    vb = np.empty(q + 1, dtype=np.float64)
    vc = np.empty(q + 1, dtype=np.float64)
    ndp = DP.shape[0]
    qp1 = q + 1

    for m in range(tris.shape[0]):
        for k in range(3):
            t_loc[k, 0] = tris[m, k, 0] - bx
            t_loc[k, 1] = tris[m, k, 1] - by
            t_loc[k, 2] = tris[m, k, 2] - bz

        e1x = t_loc[1, 0] - t_loc[0, 0]
        e1y = t_loc[1, 1] - t_loc[0, 1]
        e1z = t_loc[1, 2] - t_loc[0, 2]
        e2x = t_loc[2, 0] - t_loc[0, 0]
        e2y = t_loc[2, 1] - t_loc[0, 1]
        e2z = t_loc[2, 2] - t_loc[0, 2]
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        ln = np.sqrt(nx * nx + ny * ny + nz * nz)
        if ln == 0.0:
            continue
        nx_unit = nx / ln
        if abs(nx_unit) < 1e-15:
            continue
        if (nx > 0.0
                and abs(t_loc[0, 0] - hx) <= cop_tol
                and abs(t_loc[1, 0] - hx) <= cop_tol
                and abs(t_loc[2, 0] - hx) <= cop_tol):
            continue

        n_poly = _clip_poly_to_box_njit(t_loc, 3, b_lo, b_hi, poly_buf)
        if n_poly < 3:
            continue

        p0x, p0y, p0z = poly_buf[0, 0], poly_buf[0, 1], poly_buf[0, 2]
        for f in range(1, n_poly - 1):
            fe1x = poly_buf[f, 0] - p0x
            fe1y = poly_buf[f, 1] - p0y
            fe1z = poly_buf[f, 2] - p0z
            fe2x = poly_buf[f + 1, 0] - p0x
            fe2y = poly_buf[f + 1, 1] - p0y
            fe2z = poly_buf[f + 1, 2] - p0z
            cx = fe1y * fe2z - fe1z * fe2y
            cy = fe1z * fe2x - fe1x * fe2z
            cz = fe1x * fe2y - fe1y * fe2x
            lf = np.sqrt(cx * cx + cy * cy + cz * cz)
            if lf == 0.0:
                continue
            w_scale = 0.5 * lf * nx_unit
            for p in range(ndp):
                u = DP[p, 0]
                v = DP[p, 1]
                px = p0x + fe1x * u + fe2x * v
                py = p0y + fe1y * u + fe2y * v
                pz = p0z + fe1z * u + fe2z * v
                wp = DW[p] * w_scale
                _horner_vec(Ax, px, va)
                _horner_vec(Py, py, vb)
                _horner_vec(Pz, pz, vc)
                for ii in range(qp1):
                    w_i = wp * va[ii]
                    for jj in range(qp1):
                        w_ij = w_i * vb[jj]
                        for kk in range(qp1):
                            M[ii, jj, kk] += w_ij * vc[kk]
    return M


@njit(cache=True, inline="always")
def _clip_seg_rect_njit(ax: float, ay: float, bx: float, by: float,
                        y0: float, y1: float, z0: float, z1: float) -> tuple[bool, float, float, float, float]:
    """Liang–Barsky clip of a 2D segment to [y0, y1] x [z0, z1]."""
    t0, t1 = 0.0, 1.0
    dx, dy = bx - ax, by - ay
    for side in range(4):
        if side == 0:
            p, qq = -dx, ax - y0
        elif side == 1:
            p, qq = dx, y1 - ax
        elif side == 2:
            p, qq = -dy, ay - z0
        else:
            p, qq = dy, z1 - ay
        if abs(p) < 1e-18:
            if qq < 0.0:
                return False, 0.0, 0.0, 0.0, 0.0
            continue
        r = qq / p
        if p < 0.0:
            if r > t0:
                t0 = r
        else:
            if r < t1:
                t1 = r
        if t0 > t1:
            return False, 0.0, 0.0, 0.0, 0.0
    return True, ax + t0 * dx, ay + t0 * dy, ax + t1 * dx, ay + t1 * dy


@njit(cache=True, inline="always")
def _emit_2d_seg_njit(a2x: float, a2y: float, b2x: float, b2y: float,
                      mx: float, my: float, from_tri: bool,
                      y0: float, y1: float, z0: float, z1: float, tol: float,
                      GZ: np.ndarray, GW: np.ndarray, Hy: np.ndarray, Pz: np.ndarray,
                      seen_keys: np.ndarray, n_seen: int,
                      Njk: np.ndarray) -> int:
    """Clip, deduplicate, and accumulate one 2D oriented segment into ``Njk``."""
    ok, p0, p1, r0, r1 = _clip_seg_rect_njit(a2x, a2y, b2x, b2y, y0, y1, z0, z1)
    if not ok:
        return n_seen
    dx, dy = r0 - p0, r1 - p1
    ls = np.sqrt(dx * dx + dy * dy)
    if ls <= tol or mx == 0.0:
        return n_seen
    if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(r0) and np.isfinite(r1)
            and np.isfinite(mx) and np.isfinite(my)):
        return n_seen
    if from_tri:
        if ((abs(p0 - y0) <= tol and abs(r0 - y0) <= tol)
                or (abs(p0 - y1) <= tol and abs(r0 - y1) <= tol)
                or (abs(p1 - z0) <= tol and abs(r1 - z0) <= tol)
                or (abs(p1 - z1) <= tol and abs(r1 - z1) <= tol)):
            return n_seen

    k0 = np.int64(np.rint((p0 - y0) / tol))
    k1 = np.int64(np.rint((p1 - z0) / tol))
    k2 = np.int64(np.rint((r0 - y0) / tol))
    k3 = np.int64(np.rint((r1 - z0) / tol))
    k4 = np.int64(np.rint(mx * 1e9))
    k5 = np.int64(np.rint(my * 1e9))
    for s in range(n_seen):
        if (seen_keys[s, 4] == k4 and seen_keys[s, 5] == k5
                and ((seen_keys[s, 0] == k0 and seen_keys[s, 1] == k1
                      and seen_keys[s, 2] == k2 and seen_keys[s, 3] == k3)
                     or (seen_keys[s, 0] == k2 and seen_keys[s, 1] == k3
                         and seen_keys[s, 2] == k0 and seen_keys[s, 3] == k1))):
            return n_seen
    seen_keys[n_seen, 0] = k0
    seen_keys[n_seen, 1] = k1
    seen_keys[n_seen, 2] = k2
    seen_keys[n_seen, 3] = k3
    seen_keys[n_seen, 4] = k4
    seen_keys[n_seen, 5] = k5
    n_seen += 1

    qp1 = Njk.shape[0]
    vh = np.empty(qp1, dtype=np.float64)
    vc = np.empty(qp1, dtype=np.float64)
    pref = mx * ls
    for p in range(GZ.shape[0]):
        t = GZ[p]
        _horner_vec(Hy, p0 + dx * t, vh)
        _horner_vec(Pz, p1 + dy * t, vc)
        wp = pref * GW[p]
        for jj in range(qp1):
            w_j = wp * vh[jj]
            for kk in range(qp1):
                Njk[jj, kk] += w_j * vc[kk]
    return n_seen


@njit(cache=True)
def _cross_section_step1_njit(tris: np.ndarray, bx: float, by: float, bz: float,
                              hx: float, hy: float, hz: float, q: int,
                              GZ: np.ndarray, GW: np.ndarray, Hy: np.ndarray, Pz: np.ndarray,
                              coplanar_rtol: float):
    """Cross-section triangle segments + rect-edge probes preparation (cell-local)."""
    x1 = hx
    y0, z0, y1, z1 = 0.0, 0.0, hy, hz
    hmax = hx if hx >= hy else hy
    if hz > hmax:
        hmax = hz
    tol = coplanar_rtol * hmax
    eps = 1e-7 * hmax

    mt = tris.shape[0]
    Njk = np.zeros((q + 1, q + 1), dtype=np.float64)
    seen_cap = max(64, 4 * mt + 16)
    seen_keys = np.empty((seen_cap, 6), dtype=np.int64)
    n_seen = 0

    param_cap = max(32, 4 * mt + 8)
    params = np.empty((4, param_cap), dtype=np.float64)
    pcounts = np.full(4, 2, dtype=np.int64)
    for e in range(4):
        params[e, 0] = 0.0
        params[e, 1] = 1.0

    # Local coords of triangle
    tx = np.empty(3, dtype=np.float64)
    ty = np.empty(3, dtype=np.float64)
    tz = np.empty(3, dtype=np.float64)
    s_arr = np.empty(3, dtype=np.float64)
    px = np.empty(3, dtype=np.float64)
    py = np.empty(3, dtype=np.float64)
    pz = np.empty(3, dtype=np.float64)

    for m in range(mt):
        for k in range(3):
            tx[k] = tris[m, k, 0] - bx
            ty[k] = tris[m, k, 1] - by
            tz[k] = tris[m, k, 2] - bz
            s_arr[k] = tx[k] - x1

        all_on = (abs(s_arr[0]) <= tol) and (abs(s_arr[1]) <= tol) and (abs(s_arr[2]) <= tol)
        if not all_on:
            e1x, e1y, e1z = tx[1] - tx[0], ty[1] - ty[0], tz[1] - tz[0]
            e2x, e2y, e2z = tx[2] - tx[0], ty[2] - ty[0], tz[2] - tz[0]
            ny = e1z * e2x - e1x * e2z
            nz = e1x * e2y - e1y * e2x
            np_pts = 0
            for edge_i in range(3):
                e0 = edge_i
                e1 = 0 if edge_i == 2 else edge_i + 1
                if (s_arr[e0] < 0.0) != (s_arr[e1] < 0.0):
                    r = s_arr[e0] / (s_arr[e0] - s_arr[e1])
                    px[np_pts] = tx[e0] + r * (tx[e1] - tx[e0])
                    py[np_pts] = ty[e0] + r * (ty[e1] - ty[e0])
                    pz[np_pts] = tz[e0] + r * (tz[e1] - tz[e0])
                    np_pts += 1
            if np_pts != 2:
                np_pts = 0
                for k in range(3):
                    if abs(s_arr[k]) <= tol and np_pts < 2:
                        px[np_pts] = tx[k]
                        py[np_pts] = ty[k]
                        pz[np_pts] = tz[k]
                        np_pts += 1
            if np_pts == 2:
                a2x, a2y = py[0], pz[0]
                b2x, b2y = py[1], pz[1]
                ld = np.sqrt((b2x - a2x) ** 2 + (b2y - a2y) ** 2)
                nn = np.sqrt(ny * ny + nz * nz)
                if ld > tol and nn > 1e-15:
                    n_seen = _emit_2d_seg_njit(
                        a2x, a2y, b2x, b2y, ny / nn, nz / nn, True,
                        y0, y1, z0, z1, tol, GZ, GW, Hy, Pz, seen_keys, n_seen, Njk,
                    )
                    # Split rect edges at a2 -> b2
                    dx, dy = b2x - a2x, b2y - a2y
                    for e in range(4):
                        d_ax = dy if (e == 0 or e == 2) else dx
                        if abs(d_ax) < 1e-15:
                            continue
                        val = z0 if e == 0 else (y1 if e == 1 else (z1 if e == 2 else y0))
                        seg_ax = a2y if (e == 0 or e == 2) else a2x
                        u = (val - seg_ax) / d_ax
                        if u < -1e-9 or u > 1.0 + 1e-9:
                            continue
                        w = (a2x + u * dx) if (e == 0 or e == 2) else (a2y + u * dy)
                        lo_o = y0 if (e == 0 or e == 2) else z0
                        hi_o = y1 if (e == 0 or e == 2) else z1
                        if w < lo_o - 1e-9 or w > hi_o + 1e-9:
                            continue
                        ln = hy if (e == 0 or e == 2) else hz
                        tt = (w - lo_o) / ln if (e == 0 or e == 1) else (hi_o - w) / ln
                        if tt < 0.0:
                            tt = 0.0
                        elif tt > 1.0:
                            tt = 1.0
                        params[e, pcounts[e]] = np.rint(tt * 1e9) * 1e-9
                        pcounts[e] += 1
        else:
            for edge_i in range(3):
                k0 = edge_i
                k1 = 0 if edge_i == 2 else edge_i + 1
                a2x, a2y = ty[k0], tz[k0]
                b2x, b2y = ty[k1], tz[k1]
                dx, dy = b2x - a2x, b2y - a2y
                for e in range(4):
                    d_ax = dy if (e == 0 or e == 2) else dx
                    if abs(d_ax) < 1e-15:
                        continue
                    val = z0 if e == 0 else (y1 if e == 1 else (z1 if e == 2 else y0))
                    seg_ax = a2y if (e == 0 or e == 2) else a2x
                    u = (val - seg_ax) / d_ax
                    if u < -1e-9 or u > 1.0 + 1e-9:
                        continue
                    w = (a2x + u * dx) if (e == 0 or e == 2) else (a2y + u * dy)
                    lo_o = y0 if (e == 0 or e == 2) else z0
                    hi_o = y1 if (e == 0 or e == 2) else z1
                    if w < lo_o - 1e-9 or w > hi_o + 1e-9:
                        continue
                    ln = hy if (e == 0 or e == 2) else hz
                    tt = (w - lo_o) / ln if (e == 0 or e == 1) else (hi_o - w) / ln
                    if tt < 0.0:
                        tt = 0.0
                    elif tt > 1.0:
                        tt = 1.0
                    params[e, pcounts[e]] = np.rint(tt * 1e9) * 1e-9
                    pcounts[e] += 1

    # Build rect pieces & probes
    rect_y = np.array([y0, y1, y1, y0], dtype=np.float64)
    rect_z = np.array([z0, z0, z1, z1], dtype=np.float64)
    en_y = np.array([0.0, 1.0, 0.0, -1.0], dtype=np.float64)
    en_z = np.array([-1.0, 0.0, 1.0, 0.0], dtype=np.float64)

    max_pieces = pcounts[0] + pcounts[1] + pcounts[2] + pcounts[3]
    piece_e = np.empty(max_pieces, dtype=np.int64)
    piece_a = np.empty((max_pieces, 2), dtype=np.float64)
    piece_b = np.empty((max_pieces, 2), dtype=np.float64)
    probes = np.empty((max_pieces * 4, 3), dtype=np.float64)
    n_pieces = 0

    for e in range(4):
        cnt = pcounts[e]
        arr = np.sort(params[e, :cnt])
        r0y, r0z = rect_y[e], rect_z[e]
        e_next = 0 if e == 3 else e + 1
        r1y, r1z = rect_y[e_next], rect_z[e_next]
        dry, drz = r1y - r0y, r1z - r0z
        eny, enz = en_y[e], en_z[e]

        # Dedup sorted parameters while skipping identical duplicates
        prev_t = arr[0]
        for idx in range(1, cnt):
            tb = arr[idx]
            if tb == prev_t:
                continue
            ta = prev_t
            prev_t = tb
            if tb - ta < 1e-12:
                continue
            tm = 0.5 * (ta + tb)
            mid_y = r0y + tm * dry
            mid_z = r0z + tm * drz
            pbase = 4 * n_pieces
            for depth_i in range(2):
                ee = eps if depth_i == 0 else eps * 1e-4
                iny = mid_y - ee * eny
                inz = mid_z - ee * enz
                probes[pbase + 2 * depth_i, 0] = bx + x1 - ee
                probes[pbase + 2 * depth_i, 1] = by + iny
                probes[pbase + 2 * depth_i, 2] = bz + inz
                probes[pbase + 2 * depth_i + 1, 0] = bx + x1 + ee
                probes[pbase + 2 * depth_i + 1, 1] = by + iny
                probes[pbase + 2 * depth_i + 1, 2] = bz + inz
            piece_e[n_pieces] = e
            piece_a[n_pieces, 0] = r0y + ta * dry
            piece_a[n_pieces, 1] = r0z + ta * drz
            piece_b[n_pieces, 0] = r0y + tb * dry
            piece_b[n_pieces, 1] = r0z + tb * drz
            n_pieces += 1

    return (Njk, seen_keys, n_seen,
            piece_e[:n_pieces], piece_a[:n_pieces], piece_b[:n_pieces],
            probes[:4 * n_pieces], tol)


@njit(cache=True)
def _cross_section_step2_njit(Njk: np.ndarray, seen_keys: np.ndarray, n_seen: int,
                              piece_e: np.ndarray, piece_a: np.ndarray, piece_b: np.ndarray,
                              flags4: np.ndarray, hy: float, hz: float, tol: float,
                              GZ: np.ndarray, GW: np.ndarray, Hy: np.ndarray, Pz: np.ndarray) -> None:
    """Emit all inside rect-edge sub-segments into ``Njk``."""
    en_y = np.array([0.0, 1.0, 0.0, -1.0], dtype=np.float64)
    en_z = np.array([-1.0, 0.0, 1.0, 0.0], dtype=np.float64)
    y0, z0, y1, z1 = 0.0, 0.0, hy, hz
    for p in range(piece_e.shape[0]):
        if flags4[p, 0] or flags4[p, 1] or flags4[p, 2] or flags4[p, 3]:
            e = piece_e[p]
            n_seen = _emit_2d_seg_njit(
                piece_a[p, 0], piece_a[p, 1],
                piece_b[p, 0], piece_b[p, 1],
                en_y[e], en_z[e], False,
                y0, y1, z0, z1, tol, GZ, GW, Hy, Pz, seen_keys, n_seen, Njk,
            )


def _cell_moments(tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                  q: int, whole: np.ndarray | None = None, *,
                  accel: TriSoupAccel | None = None) -> np.ndarray:
    """Exact moments via recursive divergence theorem (tessellation-exact)."""
    blo_a = np.ascontiguousarray(blo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)
    tris_a = np.ascontiguousarray(tris, dtype=np.float64)
    if tris_a.shape[0] == 0:
        return np.zeros((q + 1,) * 3, dtype=np.float64)

    hx, hy, hz = float(h_a[0]), float(h_a[1]), float(h_a[2])
    Ax = _axis_pw_tuple(hx, q, True)
    Hy = _axis_pw_tuple(hy, q, True)
    Py = _axis_pw_tuple(hy, q, False)
    Pz = _axis_pw_tuple(hz, q, False)
    DP, DW = _duffy((3 * q + 4) // 2)
    GZ, GW = _gauss1d(q + 1)

    M = _cell_moments_part_a_njit(tris_a, blo_a, h_a, Ax, Py, Pz, DP, DW, q, _COPLANAR_RTOL)

    Njk, seen_keys, n_seen, pe, pa, pb, probes, tol = _cross_section_step1_njit(
        tris_a, float(blo_a[0]), float(blo_a[1]), float(blo_a[2]),
        hx, hy, hz, q, GZ, GW, Hy, Pz, _COPLANAR_RTOL,
    )
    if pe.shape[0] > 0:
        if accel is not None:
            ins = accel.winding_inside(probes).reshape(-1, 4)
        else:
            whole_a = tris_a if whole is None else np.ascontiguousarray(whole, dtype=np.float64)
            if whole_a.shape[0] <= 2000:
                wvals = np.empty(probes.shape[0], dtype=np.float64)
                _solid_angles_serial(np.ascontiguousarray(probes, dtype=np.float64), whole_a, wvals)
                ins = (wvals > 0.5).reshape(-1, 4)
            else:
                SN = whole_a.reshape(-1, 3)
                ST = np.arange(len(whole_a) * 3, dtype=np.int64).reshape(-1, 3)
                ins = winding_inside(probes, SN, ST).reshape(-1, 4)
        _cross_section_step2_njit(Njk, seen_keys, n_seen, pe, pa, pb, ins,
                                  hy, hz, tol, GZ, GW, Hy, Pz)

    Axi = _polyval(Ax, np.array([hx]))[:, 0]
    return M + Axi[:, None, None] * Njk[None, :, :]


def _cross_section_moments(tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                           q: int, Pz: np.ndarray, Hy: np.ndarray,
                           whole: np.ndarray | None = None, *,
                           accel: TriSoupAccel | None = None) -> np.ndarray:
    """Njk = ∫_{Ω∩{x=x1}} P_j(y)P_k(z) dydz via 1D-exact boundary pieces."""
    blo_a = np.ascontiguousarray(blo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)
    tris_loc = np.ascontiguousarray(np.asarray(tris, dtype=np.float64) - blo_a, dtype=np.float64)
    # Note: white-box callers pass local or shifted blo/Pz/Hy directly. Use general path:
    GZ, GW = _gauss1d(q + 1)
    # Shift Pz, Hy evaluated cleanly in local coordinates by delegating to step1/step2 with blo=(0,0,0) if needed
    # Preserve exact white-box contract for global coordinates:
    Njk, seen_keys, n_seen, pe, pa, pb, probes, tol = _cross_section_step1_njit(
        np.ascontiguousarray(tris, dtype=np.float64),
        float(blo_a[0]), float(blo_a[1]), float(blo_a[2]),
        float(h_a[0]), float(h_a[1]), float(h_a[2]),
        q, GZ, GW,
        _axis_pw(0.0, float(h_a[1]), legendre_monomials(q), True),
        _axis_pw(0.0, float(h_a[2]), legendre_monomials(q), False),
        _COPLANAR_RTOL,
    )
    if pe.shape[0] > 0:
        if accel is not None:
            ins = accel.winding_inside(probes).reshape(-1, 4)
        else:
            whole_a = np.ascontiguousarray(tris if whole is None else whole, dtype=np.float64)
            SN = whole_a.reshape(-1, 3)
            ST = np.arange(len(whole_a) * 3, dtype=np.int64).reshape(-1, 3)
            ins = winding_inside(probes, SN, ST).reshape(-1, 4)
        _cross_section_step2_njit(
            Njk, seen_keys, n_seen, pe, pa, pb, ins,
            float(h_a[1]), float(h_a[2]), tol, GZ, GW,
            _axis_pw(0.0, float(h_a[1]), legendre_monomials(q), True),
            _axis_pw(0.0, float(h_a[2]), legendre_monomials(q), False),
        )
    return Njk


@njit(cache=True)
def _design_matrix_njit(pts: np.ndarray, blo: np.ndarray, h: np.ndarray, q: int) -> np.ndarray:
    """Legendre tensor-product Vandermonde matrix ((q+1)^3, N) compiled in Numba."""
    n_pts = pts.shape[0]
    qp1 = q + 1
    n_rows = qp1 * qp1 * qp1
    A = np.empty((n_rows, n_pts), dtype=np.float64)
    bx, by, bz = blo[0], blo[1], blo[2]
    invx, invy, invz = 2.0 / h[0], 2.0 / h[1], 2.0 / h[2]

    px = np.empty(qp1, dtype=np.float64)
    py = np.empty(qp1, dtype=np.float64)
    pz = np.empty(qp1, dtype=np.float64)

    for p in range(n_pts):
        ux = (pts[p, 0] - bx) * invx - 1.0
        uy = (pts[p, 1] - by) * invy - 1.0
        uz = (pts[p, 2] - bz) * invz - 1.0

        px[0] = 1.0
        py[0] = 1.0
        pz[0] = 1.0
        if q >= 1:
            px[1] = ux
            py[1] = uy
            pz[1] = uz
        for j in range(1, q):
            c1 = (2.0 * j + 1.0) / (j + 1.0)
            c2 = float(j) / (j + 1.0)
            px[j + 1] = c1 * ux * px[j] - c2 * px[j - 1]
            py[j + 1] = c1 * uy * py[j] - c2 * py[j - 1]
            pz[j + 1] = c1 * uz * pz[j] - c2 * pz[j - 1]

        r = 0
        for i in range(qp1):
            pi_v = px[i]
            for j in range(qp1):
                pij = pi_v * py[j]
                for k in range(qp1):
                    A[r, p] = pij * pz[k]
                    r += 1
    return A


def _design_matrix(pts: np.ndarray, blo: np.ndarray, h: np.ndarray,
                   q: int) -> np.ndarray:
    """Vandermonde of Legendre tensor products at points (rows = moments)."""
    pts_a = np.ascontiguousarray(pts, dtype=np.float64)
    if pts_a.shape[0] == 0:
        return np.zeros(((q + 1) ** 3, 0), dtype=np.float64)
    return _design_matrix_njit(
        pts_a,
        np.ascontiguousarray(blo, dtype=np.float64),
        np.ascontiguousarray(h, dtype=np.float64),
        int(q),
    )


def fit_moment_rule(moments: np.ndarray, cand: np.ndarray,
                    blo: np.ndarray, h: np.ndarray,
                    cell_vol: float) -> tuple[np.ndarray, float, np.ndarray]:
    """NNLS fit + backtracking prune (residual ≤ tol, target ≤ (q+1)^3 pts)."""
    if len(cand) > _NNLS_MAX_CANDIDATES:  # time budget: deterministic stride
        cand = cand[:: (len(cand) + _NNLS_MAX_CANDIDATES - 1)
                    // _NNLS_MAX_CANDIDATES]
    max_pts = moments.shape[0] ** 3
    A, b = _design_matrix(cand, blo, h, moments.shape[0] - 1), moments.ravel()
    try:
        w, rnorm = nnls(A, b, maxiter=10 * A.shape[1])
    except RuntimeError:
        return (np.zeros(len(cand)), float("inf"),
                np.ascontiguousarray(cand))
    keep = w > 0
    if rnorm > _NNLS_TOL * cell_vol or int(keep.sum()) <= max_pts:
        return w * keep, float(rnorm), np.ascontiguousarray(cand)
    idx = np.nonzero(keep)[0][np.argsort(w[keep])]
    best_w, best_r = w.copy(), float(rnorm)
    dropped = ~keep
    for ii in idx:
        dropped[ii] = True
        try:
            ww, rr = nnls(A[:, ~dropped], b)
        except RuntimeError:
            dropped[ii] = False
            break
        if rr <= _NNLS_TOL * cell_vol:
            best_w = np.zeros(len(w))
            best_w[~dropped], best_r = ww, float(rr)
            if int((ww > 0).sum()) <= max_pts:
                break
        else:
            dropped[ii] = False
            break
    return best_w * (best_w > 0), best_r, np.ascontiguousarray(cand)


@njit(cache=True)
def _octree_cut_flags_njit(tris: np.ndarray, tmin: np.ndarray, tmax: np.ndarray,
                           leaf_blos: np.ndarray, leaf_hs: np.ndarray) -> np.ndarray:
    """Compiled AABB-prefiltered SAT test across all octree leaves."""
    n_leaves = leaf_blos.shape[0]
    m_tris = tris.shape[0]
    cut = np.zeros(n_leaves, dtype=np.bool_)
    bhi = np.empty(3, dtype=np.float64)
    for l in range(n_leaves):
        bx, by, bz = leaf_blos[l, 0], leaf_blos[l, 1], leaf_blos[l, 2]
        bhi[0] = bx + leaf_hs[l, 0]
        bhi[1] = by + leaf_hs[l, 1]
        bhi[2] = bz + leaf_hs[l, 2]
        bhi0, bhi1, bhi2 = bhi[0], bhi[1], bhi[2]
        blo_l = leaf_blos[l]
        hit = False
        for m in range(m_tris):
            if (tmax[m, 0] >= bx and tmin[m, 0] <= bhi0
                    and tmax[m, 1] >= by and tmin[m, 1] <= bhi1
                    and tmax[m, 2] >= bz and tmin[m, 2] <= bhi2):
                if _tri_box_overlap(tris[m, 0], tris[m, 1], tris[m, 2], blo_l, bhi):
                    hit = True
                    break
        cut[l] = hit
    return cut


def octree_fallback(tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                    q: int, p: int, depth: int,
                    log: IssueLog, whole: np.ndarray | None = None, *,
                    accel: TriSoupAccel | None = None) -> CellQuadrature:
    """Depth-`depth` octree fallback (interior leaves get Gauss (p+1)^3, cut leaves 2^3)."""
    log.add("warning", "QUAD_FALLBACK", f"octree depth {depth} in cell {blo}",
            "quadrature")
    tris_a = np.ascontiguousarray(tris, dtype=np.float64)
    blo_a = np.ascontiguousarray(blo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)

    n_sub = 2 ** int(depth)
    hh = h_a / n_sub
    # Build leaves in exact binary-subdivision Z-order matching recursive expansion
    coords = np.zeros((1, 3), dtype=np.int64)
    for _ in range(depth):
        offs = np.array([[c & 1, (c >> 1) & 1, (c >> 2) & 1] for c in range(8)], dtype=np.int64)
        coords = (coords[:, None, :] * 2 + offs[None, :, :]).reshape(-1, 3)
    leaf_blos = blo_a + coords * hh
    leaf_hs = np.broadcast_to(hh, leaf_blos.shape).copy()

    Xo, Xc = (_tensor_gauss(o)[0] for o in (p + 1, 2))
    if tris_a.shape[0] > 0:
        tmin, tmax = tris_a.min(axis=1), tris_a.max(axis=1)
        cut = _octree_cut_flags_njit(tris_a, tmin, tmax, leaf_blos, leaf_hs)
    else:
        cut = np.zeros(leaf_blos.shape[0], dtype=np.bool_)

    queries, spans = [], []
    offset = 0
    for b, is_cut in zip(leaf_blos, cut):
        q_pts = b + Xc * hh if is_cut else (b + 0.5 * hh)[None, :]
        spans.append((offset, len(q_pts)))
        offset += len(q_pts)
        queries.append(q_pts)

    if queries:
        Q_all = np.vstack(queries)
        if accel is not None:
            flags = accel.winding_inside(Q_all)
        else:
            whole_a = tris_a if whole is None else np.ascontiguousarray(whole, dtype=np.float64)
            SN, ST = whole_a.reshape(-1, 3), np.arange(len(whole_a) * 3, dtype=np.int64).reshape(-1, 3)
            flags = winding_inside(Q_all, SN, ST)
    else:
        flags = np.zeros(0, dtype=bool)

    xs, ws = [], []
    sub_vol = float(np.prod(hh))
    for b, is_cut, (o, n_q), q_pts in zip(leaf_blos, cut, spans, queries):
        f = flags[o:o + n_q]
        if is_cut:
            if f.any():
                xs.append(q_pts[f])
                ws.append(np.full(int(f.sum()), sub_vol / 8.0))
        elif f[0]:
            xs.append(b + Xo * hh)
            ws.append(np.full(len(Xo), sub_vol / len(Xo)))
    X = np.vstack(xs) if xs else np.zeros((0, 3))
    W = np.concatenate(ws) if ws else np.zeros(0)
    M = (_cell_moments(tris_a, blo_a, h_a, q, whole=whole, accel=accel) if len(X)
         else np.zeros((q + 1,) * 3))
    return CellQuadrature(x=X, w=W, moments=M, n_candidates=len(X),
                          fallback=True, volume=float(W.sum()),
                          residual=abs(float(M.ravel()[0]) - float(W.sum())))


def _local_inside(nodes: np.ndarray, tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                  margin_cells: float = 2.0, *,
                  accel: TriSoupAccel | None = None):
    """Fast inside test for candidate filtering using ``TriSoupAccel`` BVH."""
    del blo, h, margin_cells
    acc = get_accel(nodes, tris) if accel is None else accel
    return acc.winding_inside


@njit(cache=True)
def _filter_cell_candidates_njit(soup: np.ndarray, cand_ids: np.ndarray,
                                 blo: np.ndarray, h: np.ndarray, tiny: np.ndarray):
    """Exact SAT test on spatial-grid candidates for cell and shrunken interior box."""
    n_c = cand_ids.shape[0]
    touched = np.zeros(n_c, dtype=np.bool_)
    bhi = blo + h
    alo = blo + tiny
    ahi = bhi - tiny
    short_any = False
    n_touch = 0
    for q in range(n_c):
        tid = cand_ids[q]
        t0, t1, t2 = soup[tid, 0], soup[tid, 1], soup[tid, 2]
        if _tri_box_overlap(t0, t1, t2, blo, bhi):
            touched[q] = True
            n_touch += 1
            if not short_any and _tri_box_overlap(t0, t1, t2, alo, ahi):
                short_any = True
    T = np.empty((n_touch, 3, 3), dtype=np.float64)
    idx = 0
    for q in range(n_c):
        if touched[q]:
            tid = cand_ids[q]
            T[idx] = soup[tid]
            idx += 1
    return T, short_any


@njit(cache=True)
def _subbox_cand_pts_njit(T: np.ndarray, Tmin: np.ndarray, Tmax: np.ndarray,
                          blo: np.ndarray, h: np.ndarray,
                          depth: int, Xo: np.ndarray) -> np.ndarray:
    """Compiled sub-box SAT sweep and leaf quadrature candidate generation."""
    nd = 1 << depth
    hhx, hhy, hhz = h[0] / nd, h[1] / nd, h[2] / nd
    n_sub = nd * nd * nd
    bx0, by0, bz0 = blo[0], blo[1], blo[2]
    mt = T.shape[0]
    nxo = Xo.shape[0]
    mask = (1 << depth) - 1

    # Pass 1: flag sub-boxes overlapped by T
    flags = np.zeros(n_sub, dtype=np.bool_)
    b = np.empty(3, dtype=np.float64)
    bh = np.empty(3, dtype=np.float64)
    n_hit = 0
    for c in range(n_sub):
        ox = c & mask
        oy = (c >> depth) & mask
        oz = (c >> (2 * depth)) & mask
        b[0] = bx0 + ox * hhx
        b[1] = by0 + oy * hhy
        b[2] = bz0 + oz * hhz
        bh[0] = b[0] + hhx
        bh[1] = b[1] + hhy
        bh[2] = b[2] + hhz
        hit = False
        for m in range(mt):
            if (Tmax[m, 0] >= b[0] and Tmin[m, 0] <= bh[0]
                    and Tmax[m, 1] >= b[1] and Tmin[m, 1] <= bh[1]
                    and Tmax[m, 2] >= b[2] and Tmin[m, 2] <= bh[2]):
                if _tri_box_overlap(T[m, 0], T[m, 1], T[m, 2], b, bh):
                    hit = True
                    break
        if hit:
            flags[c] = True
            n_hit += 1

    out = np.empty((n_hit * nxo, 3), dtype=np.float64)
    pos = 0
    for c in range(n_sub):
        if not flags[c]:
            continue
        ox = c & mask
        oy = (c >> depth) & mask
        oz = (c >> (2 * depth)) & mask
        px0 = bx0 + ox * hhx
        py0 = by0 + oy * hhy
        pz0 = bz0 + oz * hhz
        for k in range(nxo):
            out[pos, 0] = px0 + Xo[k, 0] * hhx
            out[pos, 1] = py0 + Xo[k, 1] * hhy
            out[pos, 2] = pz0 + Xo[k, 2] * hhz
            pos += 1
    return out


def cell_quadrature(nodes: np.ndarray, tris: np.ndarray, blo: np.ndarray,
                    h: np.ndarray, q: int = 2, p: int = 2,
                    log: IssueLog | None = None, *,
                    accel: TriSoupAccel | None = None) -> CellQuadrature:
    """Volume rule for one cell: exact moments → NNLS fit → prune/fallback."""
    log = log or IssueLog()
    blo_a = np.ascontiguousarray(blo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)
    cell_vol = float(np.prod(h_a))
    acc = get_accel(nodes, tris) if accel is None else accel
    tiny = 1e-9 * h_a

    cand_ids = acc.query_box_ids(blo_a, blo_a + h_a)
    T, short_any = _filter_cell_candidates_njit(acc.soup, cand_ids, blo_a, h_a, tiny)

    if not short_any:  # cell interior sees no surface: center decides
        inside = acc.winding_inside((blo_a + 0.5 * h_a)[None, :])[0]
        if not inside:
            return CellQuadrature(np.zeros((0, 3)), np.zeros(0),
                                  np.zeros((q + 1,) * 3))
        Xg, Wg = _tensor_gauss(q + 1)
        M = np.zeros((q + 1,) * 3)
        M[0, 0, 0] = cell_vol
        return CellQuadrature(blo_a + Xg * h_a, Wg * cell_vol, M,
                              n_candidates=len(Xg), volume=cell_vol)

    win = _local_inside(nodes, tris, blo_a, h_a, accel=acc)
    M = _cell_moments(T, blo_a, h_a, q, accel=acc)
    fb = lambda: octree_fallback(T, blo_a, h_a, q, p, 3, log, accel=acc)
    if float(M[0, 0, 0]) <= 0:
        return fb()

    Xo, _ = _tensor_gauss(2)
    Xg, _ = _tensor_gauss(q + 2)
    base = (blo_a + Xg * h_a)[win(blo_a + Xg * h_a)]
    surf: list[np.ndarray] = []
    best = {"rel": float("inf"), "w": None, "cand": None}
    Tmin = T.min(axis=1)
    Tmax = T.max(axis=1)

    def _fit(Mx, cand_x, b_x, h_x, vol_x):
        wx, rx, cx = fit_moment_rule(Mx, cand_x, b_x, h_x, vol_x)
        if rx / vol_x < best["rel"]:
            best.update(rel=rx / vol_x, w=wx.copy(), cand=np.asarray(cx).copy())
        return wx, rx, cx

    for depth in (2, 3, -1, -2):  # surface leaves, dense grid, then active sub-box
        if depth == -1:
            Xd, _ = _tensor_gauss(2 * q + 4)
            surf.append((blo_a + Xd * h_a)[win(blo_a + Xd * h_a)])
        elif depth == -2:
            pts = []
            poly_buf = np.empty((32, 3), dtype=np.float64)
            bhi_a = blo_a + h_a
            for m in range(T.shape[0]):
                n_p = _clip_poly_to_box_njit(T[m], 3, blo_a, bhi_a, poly_buf)
                if n_p:
                    pts.append(poly_buf[:n_p].copy())
            corners = blo_a + np.array([[a, b, c] for c in (0, 1) for b in (0, 1)
                                        for a in (0, 1)], float) * h_a
            cin = corners[acc.winding_inside(blo_a + h_a / 2 + (corners - blo_a - h_a / 2)
                                             * (1 - 1e-9))]
            if len(cin):
                pts.append(cin)
            if not pts:
                continue
            P = np.vstack(pts)
            alo, ahi = np.maximum(P.min(axis=0), blo_a), np.minimum(P.max(axis=0), bhi_a)
            ext = np.maximum(ahi - alo, 1e-12 * h_a)
            sub_vol = float(np.prod(ext))
            sub = []
            for o in (q + 2, 2 * q + 2):
                Xs, _ = _tensor_gauss(o)
                cs = alo + Xs * ext
                sub.append(cs[win(cs)])
            if any(len(x) for x in sub):
                cand = np.vstack([x for x in sub if len(x)])
                w, rnorm, cand_i = _fit(M, cand, blo_a, h_a, cell_vol)
                if rnorm <= _NNLS_TOL * cell_vol:
                    cand = cand_i
                    break
                Ms = _cell_moments(T, alo, ext, q, accel=acc)
                if abs(float(Ms[0, 0, 0]) - float(M[0, 0, 0])) > 1e-9 * cell_vol:
                    continue
                w, rnorm_s, cand = _fit(Ms, cand, alo, ext, sub_vol)
                if rnorm_s <= _NNLS_TOL * sub_vol:
                    rnorm = rnorm_s * cell_vol / sub_vol
                    break
            continue
        else:
            leaf_pts = _subbox_cand_pts_njit(T, Tmin, Tmax, blo_a, h_a, depth, Xo)
            if len(leaf_pts):
                surf.append(leaf_pts[win(leaf_pts)])
        parts = [x for x in [base, *surf] if len(x)]
        if parts:
            cand = np.vstack(parts)
            w, rnorm, cand = _fit(M, cand, blo_a, h_a, cell_vol)
            if rnorm <= _NNLS_TOL * cell_vol:
                break
    else:
        if best["rel"] <= _NNLS_RELAXED:
            bw, bc = best["w"], best["cand"]
            nzb = bw > 0
            return CellQuadrature(bc[nzb], bw[nzb], M, float(best["rel"] * cell_vol), len(bc),
                                  False, float(M[0, 0, 0]), relaxed=True)
        return fb()
    nz = w > 0
    return CellQuadrature(cand[nz], w[nz], M, float(rnorm), len(cand),
                          False, float(M[0, 0, 0]))


@njit(cache=True)
def _surface_quadrature_njit(soup: np.ndarray, face_ids: np.ndarray,
                             blo: np.ndarray, bhi: np.ndarray,
                             DP: np.ndarray, DW: np.ndarray):
    """Compiled surface quadrature over all clipped triangles inside [blo, bhi]."""
    mt = soup.shape[0]
    ndp = DP.shape[0]
    poly_buf = np.empty((32, 3), dtype=np.float64)

    # Upper bound on output points per tri is ~7 fans * ndp; grow safely in pass 1 count + pass 2 fill
    tri_fans = np.zeros(mt, dtype=np.int64)
    total_fans = 0
    for m in range(mt):
        t = soup[m]
        n_p = _clip_poly_to_box_njit(t, 3, blo, bhi, poly_buf)
        if n_p < 3:
            continue
        e1x, e1y, e1z = t[1, 0] - t[0, 0], t[1, 1] - t[0, 1], t[1, 2] - t[0, 2]
        e2x, e2y, e2z = t[2, 0] - t[0, 0], t[2, 1] - t[0, 1], t[2, 2] - t[0, 2]
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        if nx * nx + ny * ny + nz * nz == 0.0:
            continue
        p0x, p0y, p0z = poly_buf[0, 0], poly_buf[0, 1], poly_buf[0, 2]
        f_ok = 0
        for f in range(1, n_p - 1):
            fe1x = poly_buf[f, 0] - p0x
            fe1y = poly_buf[f, 1] - p0y
            fe1z = poly_buf[f, 2] - p0z
            fe2x = poly_buf[f + 1, 0] - p0x
            fe2y = poly_buf[f + 1, 1] - p0y
            fe2z = poly_buf[f + 1, 2] - p0z
            ax = fe1y * fe2z - fe1z * fe2y
            ay = fe1z * fe2x - fe1x * fe2z
            az = fe1x * fe2y - fe1y * fe2x
            if ax * ax + ay * ay + az * az > 0.0:
                f_ok += 1
        tri_fans[m] = f_ok
        total_fans += f_ok

    n_pts = total_fans * ndp
    xs = np.empty((n_pts, 3), dtype=np.float64)
    ws = np.empty(n_pts, dtype=np.float64)
    ns = np.empty((n_pts, 3), dtype=np.float64)
    fs = np.empty(n_pts, dtype=np.int64)
    if n_pts == 0:
        return xs, ws, ns, fs

    pos = 0
    for m in range(mt):
        if tri_fans[m] == 0:
            continue
        t = soup[m]
        n_p = _clip_poly_to_box_njit(t, 3, blo, bhi, poly_buf)
        e1x, e1y, e1z = t[1, 0] - t[0, 0], t[1, 1] - t[0, 1], t[1, 2] - t[0, 2]
        e2x, e2y, e2z = t[2, 0] - t[0, 0], t[2, 1] - t[0, 1], t[2, 2] - t[0, 2]
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        ln = np.sqrt(nx * nx + ny * ny + nz * nz)
        ux, uy, uz = nx / ln, ny / ln, nz / ln
        fid = face_ids[m]
        p0x, p0y, p0z = poly_buf[0, 0], poly_buf[0, 1], poly_buf[0, 2]
        for f in range(1, n_p - 1):
            fe1x = poly_buf[f, 0] - p0x
            fe1y = poly_buf[f, 1] - p0y
            fe1z = poly_buf[f, 2] - p0z
            fe2x = poly_buf[f + 1, 0] - p0x
            fe2y = poly_buf[f + 1, 1] - p0y
            fe2z = poly_buf[f + 1, 2] - p0z
            ax = fe1y * fe2z - fe1z * fe2y
            ay = fe1z * fe2x - fe1x * fe2z
            az = fe1x * fe2y - fe1y * fe2x
            a2 = np.sqrt(ax * ax + ay * ay + az * az)
            if a2 == 0.0:
                continue
            half_a2 = 0.5 * a2
            for k in range(ndp):
                u = DP[k, 0]
                v = DP[k, 1]
                xs[pos, 0] = p0x + fe1x * u + fe2x * v
                xs[pos, 1] = p0y + fe1y * u + fe2y * v
                xs[pos, 2] = p0z + fe1z * u + fe2z * v
                ws[pos] = DW[k] * half_a2
                ns[pos, 0] = ux
                ns[pos, 1] = uy
                ns[pos, 2] = uz
                fs[pos] = fid
                pos += 1
    return xs, ws, ns, fs


def surface_quadrature(nodes: np.ndarray, tris: np.ndarray,
                       face_ids: np.ndarray, blo: np.ndarray, h: np.ndarray,
                       p: int = 2) -> dict[str, np.ndarray]:
    """Boundary rule: clipped-triangle fan + Duffy Gauss (exact ≥ 2p+1)."""
    DP, DW = _duffy(p + 1)
    blo_a = np.ascontiguousarray(blo, dtype=np.float64)
    h_a = np.ascontiguousarray(h, dtype=np.float64)
    soup = np.ascontiguousarray(np.asarray(nodes, dtype=np.float64)[np.asarray(tris, dtype=np.int64)], dtype=np.float64)
    fids = np.ascontiguousarray(face_ids, dtype=np.int64)
    xs, ws, ns, fs = _surface_quadrature_njit(soup, fids, blo_a, blo_a + h_a, DP, DW)
    return {"x": xs, "w": ws, "n": ns, "face_id": fs}


def check_volume_rule(rule: CellQuadrature, tess_volume: float) -> float:
    """|Σw − tess volume| / tess volume (mandatory, rel 1e-12)."""
    return abs(float(rule.w.sum()) - tess_volume) / abs(tess_volume)


def check_surface_rule(surf: dict[str, np.ndarray], tess_area: float) -> float:
    """|Σw − tess area| / tess area."""
    return abs(float(surf["w"].sum()) - tess_area) / abs(tess_area)


def closed_surface_flux(surf: dict[str, np.ndarray]) -> np.ndarray:
    """Σ w·n over a closed surface (must vanish)."""
    return (surf["w"][:, None] * surf["n"]).sum(axis=0)
