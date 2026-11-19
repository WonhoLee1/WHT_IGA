# -*- coding: utf-8 -*-
"""Trimmed (immersed) mean-dilatation F-bar / P0 with cut-cell pressure groups.

Extends :mod:`wht_iga.analyses.fbar` (Simo-Taylor-Pister mean dilatation,
Total-Lagrangian Neo-Hookean, docs/spikes/S40.md) from untrimmed boxes to a
Cartesian ``BoxSpline`` background grid trimmed by
:func:`wht_iga.solver.trimmed.build_trimmed_table`:

* OUT cells (``out_cells``) are dropped. Cut cells (``qp_table`` keys) use their
  moment-fitted rules ``[x, y, z, w]`` (physical points; the box map is the
  identity, so dN/dX = dN/dxi). Every other active cell uses the tensor Gauss
  rule. Basis functions of a cell are always evaluated with the cell's own
  knot span, so rule points on a cell face stay on the cell's polynomial piece.
* Inactive control points (support only in OUT cells,
  :func:`wht_iga.solver.active.inactive_ctrl_mask`) belong to no cell, so their
  rows and columns are empty; :func:`solve_static` fixes them to 0.
* Pressure groups (macro elements). Cell-wise P0 gives a small cut cell
  (volume fraction ``eta = |cell & Omega| / |cell| < eta_thresh``) its own
  volumetric constraint on a sliver seen by few basis functions.
  :func:`pressure_groups` attaches each small active cell to a large one
  (``eta >= eta_thresh``) and the mean dilatation is taken per group::

      W   = sum_q w_q mu/2 (tr C_q - 3) + sum_g W_g U(Jbar_g)
      U(J) = lam/2 (ln J)^2 - mu ln J
      Jbar_g = sum_{q in g} w_q J_q / W_g,          W_g = sum_{q in g} w_q
      R_I = sum_q w_q mu F_q G_qI + U'(Jbar_g) V_gI,
      V_gI = sum_{q in g} w_q J_q F_q^-T G_qI
      K   = K_dev + U''(Jbar_g) / W_g * V_g (x) V_g + U'(Jbar_g) dV_g/du

  The middle term couples every pair of control points of a group (not only
  of a cell), so the CSR pattern is built from per-group node unions.
  ``eta_thresh = 0`` gives one group per cell, i.e. the cell-wise P0 of
  ``fbar`` (untrimmed box: identical residual and tangent).

  Measured 2026-11-19 on V07 (Kirsch, nu = 0.49999, p = 2, h 0.5 / 0.25,
  smallest cut fraction 0.087 / 0.042): cell-wise P0 does NOT lock there
  (e_u 1.92e-3 / 2.88e-4, e_sigma_xx 5.23e-2 / 2.50e-2 vs agglomerated
  2.78e-3 / 3.27e-4 and 5.62e-2 / 2.59e-2). The volumetric penalty of a cell
  scales with its volume, so a sliver constrains little; a radius-1.0605
  variant (fraction 2e-7, h 0.25) did not lock either (cell-wise e_u 2.90e-4,
  e_sigma_xx 2.43e-2; agglomerated 3.68e-4, 2.44e-2). What
  agglomeration does improve is the conditioning: cond(D K D) (Jacobi-scaled
  u = 0 tangent) is 0.36x / 0.63x of the cell-wise value on V07 and 0.65x at
  fraction 2e-7. It also changes the local pressure: per-cell mean pressure
  error in the small cells is lower at h 0.5 (8.8e-2 vs 2.1e-1 of the
  sigma_xx scale) and higher at h 0.25 (1.4e-1 vs 3.3e-2), because a group
  averages the pressure over a larger region where its gradient is steep.

Displacement conditioning: no basis extension and no ghost penalty; the linear
solve uses symmetric Jacobi scaling instead. Agglomeration fixes the pressure
space only; the displacement space keeps every active B-spline, and a basis
function whose support meets Omega only in slivers has a tiny diagonal.
Measured 2026-11-19 on V07 (nu = 0.49999, p = 2, u = 0 tangent on the free
DOFs): cond(K) = 1.1e11 (h 0.5) / 3.7e12 (h 0.25) against 1.8e5 / 5.9e5 for
the same box without the hole, but cond(D K D) with D = diag(K)^-1/2 is
2.3e5 / 1.0e6 (agglomerated; 6.5e5 / 1.6e6 cell-wise) against 1.7e5 / 5.2e5.
A radius-1.0605 variant (smallest basis support fraction 2e-8, unscaled K
singular in double precision) still gives cond(D K D) = 8.2e5. So the
ill-conditioning is diagonal and scaling removes it without changing the
discrete space (extension, :mod:`wht_iga.stabilize.extension`) or the weak
form (ghost penalty). What scaling cannot fix: the coefficient of a weakly
supported basis is set by the few sliver points it sees, so it can be large
(V07 h 0.25: 20 x max|u|) and depends on the pressure space, while u_h on
Omega does not (agglomerated vs cell-wise differ by < 1 % of max|u| at the
QPs). ``model["diag"]["min_basis_support_ratio"]`` reports the smallest
``|supp N_I & Omega| / |supp N_I|`` of the active basis.

Scope: ``BoxSpline`` (identity geometry, uniform C^(p-1) knots, as required by
``inactive_ctrl_mask``). No Kratos dependency.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from wht_iga.analyses.fbar import _det3, _inv3
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.solver.active import inactive_ctrl_mask

#: Default small-cell threshold (cell volume fraction) for pressure agglomeration.
DEFAULT_ETA_THRESH = 0.25

_FACE_NB = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
_FULL_NB = tuple((di, dj, dk) for dk in (-1, 0, 1) for dj in (-1, 0, 1)
                 for di in (-1, 0, 1) if (di, dj, dk) != (0, 0, 0))


# ---------------------------------------------------------------------------
# pressure groups
# ---------------------------------------------------------------------------

def pressure_groups(eta: np.ndarray, eta_thresh: float = DEFAULT_ETA_THRESH) -> np.ndarray:
    """Pressure-group id per background cell.

    ``eta``: (nx, ny, nz) inside volume fraction; ``eta <= 0`` marks an OUT
    cell (no group, id -1). Every large cell (``eta >= eta_thresh``) roots one
    group; roots are numbered in flat cell order ``i + nx*(j + ny*k)``.

    Small active cells are attached by synchronous multi-source breadth-first
    search: sweep s attaches every still unassigned small cell that has a face
    neighbour assigned in an earlier sweep to the group of the neighbour with
    the largest eta (tie: smallest flat id). So each group is one large cell
    plus small cells joined to it through faces, and its volume is at least
    ``eta_thresh`` cell volumes. If a sweep attaches nothing, the next sweep also
    accepts edge/vertex neighbours; if that fails too, the remaining small cells
    have no path to a large cell and ``ValueError`` is raised (refine the grid
    or lower ``eta_thresh``).

    ``eta_thresh = 0`` makes every active cell its own group (cell-wise P0).
    """
    eta = np.asarray(eta, dtype=np.float64)
    if eta.ndim != 3:
        raise ValueError(f"eta must be (nx, ny, nz), got shape {eta.shape}")
    thr = float(eta_thresh)
    if not 0.0 <= thr <= 1.0:
        raise ValueError(f"eta_thresh={thr} outside [0, 1]")
    nx, ny, nz = eta.shape
    active = eta > 0.0
    large = active & (eta >= thr)
    grp = np.full(eta.shape, -1, dtype=np.int64)
    # argwhere on the (k, j, i) view lists cells in flat order i + nx*(j + ny*k)
    roots = np.argwhere(large.transpose(2, 1, 0))[:, ::-1]
    if active.any() and roots.shape[0] == 0:
        raise ValueError(f"no active cell reaches eta_thresh={thr}: no pressure group root")
    for g, (i, j, k) in enumerate(roots):
        grp[i, j, k] = g
    pending = active & ~large
    use_full = False
    while pending.any():
        snap = grp.copy()
        nbs = _FULL_NB if use_full else _FACE_NB
        n_new = 0
        for i, j, k in np.argwhere(pending.transpose(2, 1, 0))[:, ::-1]:
            best_key, best_g = None, -1
            for di, dj, dk in nbs:
                a, b, c = i + di, j + dj, k + dk
                if 0 <= a < nx and 0 <= b < ny and 0 <= c < nz and snap[a, b, c] >= 0:
                    key = (-eta[a, b, c], a + nx * (b + ny * c))
                    if best_key is None or key < best_key:
                        best_key, best_g = key, int(snap[a, b, c])
            if best_g >= 0:
                grp[i, j, k] = best_g
                pending[i, j, k] = False
                n_new += 1
        if n_new == 0:
            if use_full:
                left = np.argwhere(pending)[:5].tolist()
                raise ValueError(f"small cells without a path to a large cell (eta_thresh={thr}): "
                                 f"{int(pending.sum())} cells, e.g. {left}")
            use_full = True
        else:
            use_full = False
    return grp


# ---------------------------------------------------------------------------
# per-cell basis data
# ---------------------------------------------------------------------------

@njit(cache=True)
def _cell_basis(pts: np.ndarray, sx: int, sy: int, sz: int,
                kx: np.ndarray, ky: np.ndarray, kz: np.ndarray, p: int):
    """N (nq, nen) and dN/dX (nq, nen, 3) on the cell with spans (sx, sy, sz).

    Local node ii = ai + (p+1)*(aj + (p+1)*ak) (same order as ``fbar``).
    """
    nq = pts.shape[0]
    n1 = p + 1
    nen = n1 * n1 * n1
    N = np.empty((nq, nen), dtype=np.float64)
    G = np.empty((nq, nen, 3), dtype=np.float64)
    for q in range(nq):
        dx = ders_basis_funs(sx, pts[q, 0], p, 1, kx)
        dy = ders_basis_funs(sy, pts[q, 1], p, 1, ky)
        dz = ders_basis_funs(sz, pts[q, 2], p, 1, kz)
        for ak in range(n1):
            for aj in range(n1):
                for ai in range(n1):
                    ii = ai + n1 * (aj + n1 * ak)
                    N[q, ii] = dx[0, ai] * dy[0, aj] * dz[0, ak]
                    G[q, ii, 0] = dx[1, ai] * dy[0, aj] * dz[0, ak]
                    G[q, ii, 1] = dx[0, ai] * dy[1, aj] * dz[0, ak]
                    G[q, ii, 2] = dx[0, ai] * dy[0, aj] * dz[1, ak]
    return N, G


def _axis_spans(knots: np.ndarray, p: int, n_cells: int) -> np.ndarray:
    """Knot span of every cell along one axis (``find_span`` at the cell midpoint)."""
    n_cp = len(knots) - p - 1
    br = np.unique(knots)
    if len(br) - 1 != n_cells:
        raise ValueError(f"knot vector has {len(br) - 1} cells, box says {n_cells}")
    out = np.empty(n_cells, dtype=np.int64)
    for c in range(n_cells):
        lo, hi = float(br[c]), float(br[c + 1])
        s = int(find_span(n_cp, p, 0.5 * (lo + hi), knots))
        if knots[s] != lo or knots[s + 1] != hi:
            raise ValueError(f"cell {c}: span {s} does not cover [{lo}, {hi}]")
        out[c] = s
    return out


# ---------------------------------------------------------------------------
# CSR pattern over group node unions
# ---------------------------------------------------------------------------

def _group_csr_pattern(grp_nptr: np.ndarray, grp_nodes: np.ndarray, n_ctrl: int) -> dict:
    """CSR pattern of the 3-DOF-per-node tangent; nodes are coupled iff they share a group.

    Entry (3n+a, 3m+b) sits at ``indptr[3n+a] + 3*pos + b`` where, for local
    nodes (li, lj) of group g, ``pos = grp_pos[grp_pptr[g] + li*nen_g + lj]``
    is the rank of node m among the sorted neighbours of node n.
    """
    n_groups = grp_nptr.shape[0] - 1
    keys = []
    pptr = np.zeros(n_groups + 1, dtype=np.int64)
    for g in range(n_groups):
        nd = grp_nodes[grp_nptr[g]:grp_nptr[g + 1]].astype(np.int64)
        keys.append((nd[:, None] * n_ctrl + nd[None, :]).ravel())
        pptr[g + 1] = pptr[g] + nd.size * nd.size
    key = np.concatenate(keys) if keys else np.zeros(0, dtype=np.int64)
    uniq, inv = np.unique(key, return_inverse=True)
    rows, cols = uniq // n_ctrl, uniq % n_ctrl
    deg = np.bincount(rows, minlength=n_ctrl).astype(np.int64)
    nstart = np.zeros(n_ctrl + 1, dtype=np.int64)
    np.cumsum(deg, out=nstart[1:])
    rank = np.arange(uniq.size, dtype=np.int64) - nstart[rows]
    pos = np.ascontiguousarray(rank[np.asarray(inv).ravel()])
    row_node = np.repeat(np.arange(n_ctrl, dtype=np.int64), 3)
    row_len = 3 * deg[row_node]
    indptr = np.zeros(3 * n_ctrl + 1, dtype=np.int64)
    np.cumsum(row_len, out=indptr[1:])
    nnz = int(indptr[-1])
    nb_cols = (3 * cols[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
    gather = (np.repeat(3 * nstart[row_node] - indptr[:-1], row_len)
              + np.arange(nnz, dtype=np.int64))
    return {"indptr": indptr, "indices": nb_cols[gather], "grp_pptr": pptr,
            "grp_pos": pos, "nnz": nnz}


# ---------------------------------------------------------------------------
# assembly kernel
# ---------------------------------------------------------------------------

@njit(cache=True)
def _assemble(u: np.ndarray, mu: float, lam: float,
              grp_cptr: np.ndarray, grp_cells: np.ndarray,
              grp_nptr: np.ndarray, grp_nodes: np.ndarray,
              grp_pptr: np.ndarray, grp_pos: np.ndarray,
              cell_loc: np.ndarray, cell_qptr: np.ndarray,
              qG: np.ndarray, qW: np.ndarray,
              indptr: np.ndarray, nnz: int, with_K: bool):
    """Global residual, CSR data of the tangent (if ``with_K``) and Jbar per group.

    A non-positive group Jbar gives NaN (log of a negative number) in the
    residual, which ``solve_static`` treats as an inadmissible trial state.
    """
    n_groups = grp_cptr.shape[0] - 1
    nen_c = cell_loc.shape[1]
    R = np.zeros(u.shape[0], dtype=np.float64)
    data = np.zeros(nnz if with_K else 0, dtype=np.float64)
    jbar_out = np.empty(n_groups, dtype=np.float64)
    H = np.empty((nen_c, 3), dtype=np.float64)
    for g in range(n_groups):
        n0 = grp_nptr[g]
        nen = grp_nptr[g + 1] - n0
        nq = 0
        for cc in range(grp_cptr[g], grp_cptr[g + 1]):
            c = grp_cells[cc]
            nq += cell_qptr[c + 1] - cell_qptr[c]
        F = np.zeros((nq, 3, 3), dtype=np.float64)
        dJ = np.empty(nq, dtype=np.float64)
        Wg = 0.0
        sJ = 0.0
        t = 0
        for cc in range(grp_cptr[g], grp_cptr[g + 1]):
            c = grp_cells[cc]
            for q in range(cell_qptr[c], cell_qptr[c + 1]):
                for a in range(3):
                    F[t, a, a] = 1.0
                for I in range(nen_c):
                    n = grp_nodes[n0 + cell_loc[c, I]]
                    for a in range(3):
                        ua = u[3 * n + a]
                        for m in range(3):
                            F[t, a, m] += ua * qG[q, I, m]
                dJ[t] = _det3(F[t])
                Wg += qW[q]
                sJ += qW[q] * dJ[t]
                t += 1
        Jbar = sJ / Wg
        jbar_out[g] = Jbar
        lnJ = np.log(Jbar)
        pbar = (lam * lnJ - mu) / Jbar
        ppbar = (lam - lam * lnJ + mu) / (Jbar * Jbar)
        Rg = np.zeros(3 * nen, dtype=np.float64)
        Kg = np.zeros((3 * nen, 3 * nen) if with_K else (0, 0), dtype=np.float64)
        V = np.zeros((nen, 3), dtype=np.float64)
        t = 0
        for cc in range(grp_cptr[g], grp_cptr[g + 1]):
            c = grp_cells[cc]
            for q in range(cell_qptr[c], cell_qptr[c + 1]):
                FinvT = _inv3(F[t]).T
                w = qW[q]
                wJ = w * dJ[t]
                for I in range(nen_c):
                    li = cell_loc[c, I]
                    for a in range(3):
                        hs = 0.0
                        fa = 0.0
                        for m in range(3):
                            hs += FinvT[a, m] * qG[q, I, m]
                            fa += F[t, a, m] * qG[q, I, m]
                        H[I, a] = hs
                        V[li, a] += wJ * hs
                        Rg[3 * li + a] += w * mu * fa
                if with_K:
                    cq = pbar * wJ
                    for I in range(nen_c):
                        li = cell_loc[c, I]
                        r0, r1, r2 = 3 * li, 3 * li + 1, 3 * li + 2
                        for J in range(I, nen_c):
                            lj = cell_loc[c, J]
                            s = w * mu * (qG[q, I, 0] * qG[q, J, 0] + qG[q, I, 1] * qG[q, J, 1]
                                          + qG[q, I, 2] * qG[q, J, 2])
                            for a in range(3):
                                Kg[3 * li + a, 3 * lj + a] += s
                                if J != I:
                                    Kg[3 * lj + a, 3 * li + a] += s
                        # pbar * dV/du: block (I, J) = cq*(H_I H_J^T - H_J H_I^T),
                        # antisymmetric, block (J, I) its transpose, (I, I) zero
                        for J in range(I + 1, nen_c):
                            lj = cell_loc[c, J]
                            c0, c1, c2 = 3 * lj, 3 * lj + 1, 3 * lj + 2
                            v01 = cq * (H[I, 0] * H[J, 1] - H[J, 0] * H[I, 1])
                            v02 = cq * (H[I, 0] * H[J, 2] - H[J, 0] * H[I, 2])
                            v12 = cq * (H[I, 1] * H[J, 2] - H[J, 1] * H[I, 2])
                            Kg[r0, c1] += v01
                            Kg[r1, c0] -= v01
                            Kg[r0, c2] += v02
                            Kg[r2, c0] -= v02
                            Kg[r1, c2] += v12
                            Kg[r2, c1] -= v12
                            Kg[c1, r0] += v01
                            Kg[c0, r1] -= v01
                            Kg[c2, r0] += v02
                            Kg[c0, r2] -= v02
                            Kg[c2, r1] += v12
                            Kg[c1, r2] -= v12
                t += 1
        for li in range(nen):
            for a in range(3):
                Rg[3 * li + a] += pbar * V[li, a]
        if with_K:
            cv = ppbar / Wg
            for li in range(nen):
                for a in range(3):
                    va = cv * V[li, a]
                    for lj in range(nen):
                        for b in range(3):
                            Kg[3 * li + a, 3 * lj + b] += va * V[lj, b]
        p0 = grp_pptr[g]
        for li in range(nen):
            n = grp_nodes[n0 + li]
            for a in range(3):
                R[3 * n + a] += Rg[3 * li + a]
                if with_K:
                    rs = indptr[3 * n + a]
                    for lj in range(nen):
                        pos = rs + 3 * grp_pos[p0 + li * nen + lj]
                        for b in range(3):
                            data[pos + b] += Kg[3 * li + a, 3 * lj + b]
    return R, data, jbar_out


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

def build_model(box, qp_table: dict, out_cells, E: float, nu: float, *,
                eta_thresh: float = DEFAULT_ETA_THRESH,
                n_gauss: int | None = None) -> dict:
    """Precompute QP data, pressure groups and the CSR pattern of a trimmed box.

    ``box``: ``BoxSpline`` (control points at the Greville abscissae, i.e. the
    identity geometry map; checked). ``qp_table`` / ``out_cells``: the first two
    return values of ``build_trimmed_table`` (cut-cell rules ``(nq, 4)
    [x, y, z, w]`` and cells without physical volume). ``eta_thresh``: small-cell
    threshold of :func:`pressure_groups` (0 = cell-wise P0). ``n_gauss``: points
    per axis of the uncut cells (default p+1, as ``fbar``).
    """
    p = int(box.degree)
    ng = p + 1 if n_gauss is None else int(n_gauss)
    nx, ny, nz = (int(v) for v in box.n_cells)
    ncx, ncy, ncz = (int(v) for v in box.n_cp)
    n_ctrl = ncx * ncy * ncz
    knots = [np.ascontiguousarray(k, dtype=np.float64) for k in box.knots]
    gre = [greville(k, p) for k in knots]
    W3, V3, U3 = np.meshgrid(gre[2], gre[1], gre[0], indexing="ij")
    cp_id = np.column_stack([U3.ravel(), V3.ravel(), W3.ravel()])
    cp = np.asarray(box.control_points(), dtype=np.float64)
    if cp.shape != cp_id.shape or not np.allclose(cp, cp_id, rtol=0.0,
                                                  atol=1e-12 * (1.0 + np.max(np.abs(cp_id)))):
        raise ValueError("fbar_trimmed needs the identity geometry map (BoxSpline: "
                         "control points at the Greville abscissae)")
    inactive = inactive_ctrl_mask(box, out_cells)  # raises unless C^(p-1) knots
    lo = np.asarray(box.lo, dtype=np.float64)
    hi = np.asarray(box.hi, dtype=np.float64)
    h = (hi - lo) / np.array([nx, ny, nz], dtype=np.float64)
    cell_vol = float(np.prod(h))
    out_set = {tuple(int(v) for v in c) for c in out_cells}
    cut_keys = {tuple(int(v) for v in k) for k in qp_table}
    if cut_keys & out_set:
        raise ValueError(f"cells both cut and OUT: {sorted(cut_keys & out_set)[:5]}")
    for i, j, k in cut_keys:
        if not (0 <= i < nx and 0 <= j < ny and 0 <= k < nz):
            raise ValueError(f"qp_table cell {(i, j, k)} outside the {box.n_cells} grid")
    spans = [_axis_spans(knots[d], p, n) for d, n in enumerate((nx, ny, nz))]

    xg, wg = np.polynomial.legendre.leggauss(ng)
    s01 = 0.5 * (xg + 1.0)
    ref = np.array([[s01[a], s01[b], s01[c]] for c in range(ng) for b in range(ng)
                    for a in range(ng)])  # q = a + ng*(b + ng*c), as fbar
    wref = np.array([wg[a] * wg[b] * wg[c] for c in range(ng) for b in range(ng)
                     for a in range(ng)]) * (cell_vol / 8.0)
    n1 = p + 1
    ax = np.arange(n1)
    cell_ijk, cell_nodes, cell_cut, cell_eta = [], [], [], []
    qX, qW, qN, qG, q_cell = [], [], [], [], []
    eta_grid = np.zeros((nx, ny, nz))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (i, j, k) in out_set:
                    continue
                if (i, j, k) in cut_keys:
                    rule = np.asarray(qp_table[(i, j, k)], dtype=np.float64).reshape(-1, 4)
                    X, W = np.ascontiguousarray(rule[:, :3]), np.ascontiguousarray(rule[:, 3])
                    if X.shape[0] == 0 or not np.all(W >= 0.0) or float(W.sum()) <= 0.0:
                        raise ValueError(f"cut cell {(i, j, k)}: empty or negative rule")
                    eta = float(W.sum()) / cell_vol
                    cut = True
                else:
                    X = np.ascontiguousarray(lo + h * np.array([i, j, k], dtype=np.float64) + ref * h)
                    W = wref.copy()
                    eta = 1.0
                    cut = False
                sx, sy, sz = int(spans[0][i]), int(spans[1][j]), int(spans[2][k])
                N, G = _cell_basis(X, sx, sy, sz, knots[0], knots[1], knots[2], p)
                nodes = ((sx - p + ax)[None, None, :]
                         + ncx * ((sy - p + ax)[None, :, None]
                                  + ncy * (sz - p + ax)[:, None, None])).reshape(n1 ** 3)
                a_id = len(cell_ijk)
                cell_ijk.append((i, j, k))
                cell_nodes.append(nodes.astype(np.int64))
                cell_cut.append(cut)
                cell_eta.append(eta)
                eta_grid[i, j, k] = eta
                qX.append(X)
                qW.append(W)
                qN.append(N)
                qG.append(G)
                q_cell.append(np.full(W.shape[0], a_id, dtype=np.int64))
    n_act = len(cell_ijk)
    if n_act == 0:
        raise ValueError("no active cell")
    cell_nodes_a = np.ascontiguousarray(np.stack(cell_nodes))
    used = np.zeros(n_ctrl, dtype=bool)
    used[cell_nodes_a.ravel()] = True
    if not np.array_equal(used, ~inactive):
        raise RuntimeError("inactive_ctrl_mask disagrees with the active-cell connectivity")
    counts = np.array([w.shape[0] for w in qW], dtype=np.int64)
    cell_qptr = np.zeros(n_act + 1, dtype=np.int64)
    np.cumsum(counts, out=cell_qptr[1:])

    grp_grid = pressure_groups(eta_grid, eta_thresh)
    cell_group = np.array([grp_grid[c] for c in cell_ijk], dtype=np.int64)
    n_groups = int(cell_group.max()) + 1
    grp_cells = np.argsort(cell_group, kind="stable").astype(np.int64)
    grp_cptr = np.zeros(n_groups + 1, dtype=np.int64)
    np.cumsum(np.bincount(cell_group, minlength=n_groups), out=grp_cptr[1:])
    grp_node_list, cell_loc = [], np.empty_like(cell_nodes_a)
    grp_nptr = np.zeros(n_groups + 1, dtype=np.int64)
    for g in range(n_groups):
        cs = grp_cells[grp_cptr[g]:grp_cptr[g + 1]]
        un = np.unique(cell_nodes_a[cs].ravel())
        for c in cs:
            cell_loc[c] = np.searchsorted(un, cell_nodes_a[c])
        grp_node_list.append(un)
        grp_nptr[g + 1] = grp_nptr[g] + un.size
    grp_nodes = np.ascontiguousarray(np.concatenate(grp_node_list).astype(np.int64))
    csr = _group_csr_pattern(grp_nptr, grp_nodes, n_ctrl)

    qW_a = np.ascontiguousarray(np.concatenate(qW))
    q_cell_a = np.concatenate(q_cell)
    group_volume = np.bincount(cell_group[q_cell_a], weights=qW_a, minlength=n_groups)
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    inactive_dofs = (3 * np.flatnonzero(inactive)[:, None] + np.arange(3)[None, :]).ravel()
    cell_eta_a = np.asarray(cell_eta)
    model = {
        "n_ctrl": n_ctrl, "n_dof": 3 * n_ctrl, "degree": p, "mu": float(mu), "lam": float(lam),
        "E": float(E), "nu": float(nu), "eta_thresh": float(eta_thresh), "cell_vol": cell_vol,
        "box_n_cells": (nx, ny, nz),
        "cell_ijk": np.asarray(cell_ijk, dtype=np.int64), "cell_nodes": cell_nodes_a,
        "cell_cut": np.asarray(cell_cut, dtype=bool), "cell_eta": cell_eta_a,
        "cell_group": cell_group, "cell_qptr": cell_qptr, "cell_loc": np.ascontiguousarray(cell_loc),
        "eta_grid": eta_grid, "group_grid": grp_grid,
        "n_groups": n_groups, "grp_cptr": grp_cptr, "grp_cells": grp_cells,
        "grp_nptr": grp_nptr, "grp_nodes": grp_nodes, "group_volume": group_volume,
        "qX": np.ascontiguousarray(np.concatenate(qX)), "qW": qW_a,
        "qN": np.ascontiguousarray(np.concatenate(qN)),
        "qG": np.ascontiguousarray(np.concatenate(qG)), "q_cell": q_cell_a,
        "inactive_ctrl": inactive, "inactive_dofs": inactive_dofs.astype(np.int64),
        "csr": csr,
    }
    sizes = np.diff(grp_cptr)
    model["diag"] = {
        "n_active": n_act, "n_cut": int(np.sum(cell_cut)),
        "n_small": int(np.sum(cell_eta_a < float(eta_thresh))) if eta_thresh > 0 else 0,
        "n_groups": n_groups, "max_group_cells": int(sizes.max()),
        "min_eta": float(cell_eta_a.min()),
        "min_group_fraction": float(group_volume.min() / cell_vol),
        "n_inactive_ctrl": int(inactive.sum()),
        "n_qp": int(qW_a.size),
        "min_basis_support_ratio": _min_support_ratio(cell_nodes_a, cell_eta_a, box, inactive),
    }
    return model


def _min_support_ratio(cell_nodes: np.ndarray, cell_eta: np.ndarray, box, inactive) -> float:
    """min over active basis I of |supp N_I & Omega| / |supp N_I| (cells of equal size)."""
    p = int(box.degree)
    ncx, ncy, ncz = (int(v) for v in box.n_cp)
    nx, ny, nz = (int(v) for v in box.n_cells)
    inside = np.zeros(ncx * ncy * ncz)
    np.add.at(inside, cell_nodes.ravel(), np.repeat(cell_eta, cell_nodes.shape[1]))
    I = np.flatnonzero(~inactive)
    a, b, c = I % ncx, (I // ncx) % ncy, I // (ncx * ncy)

    def n_sup(idx, n_cells):  # cells in the support of 1D basis idx (C^(p-1) clamped)
        return np.minimum(idx, n_cells - 1) - np.maximum(idx - p, 0) + 1

    tot = n_sup(a, nx) * n_sup(b, ny) * n_sup(c, nz)
    return float(np.min(inside[I] / tot)) if I.size else 0.0


# ---------------------------------------------------------------------------
# residual / tangent / post
# ---------------------------------------------------------------------------

def _call(model: dict, u: np.ndarray, with_K: bool):
    s = model["csr"]
    return _assemble(np.ascontiguousarray(u, dtype=np.float64).reshape(int(model["n_dof"])),
                     float(model["mu"]), float(model["lam"]),
                     model["grp_cptr"], model["grp_cells"], model["grp_nptr"], model["grp_nodes"],
                     s["grp_pptr"], s["grp_pos"], model["cell_loc"], model["cell_qptr"],
                     model["qG"], model["qW"], s["indptr"], int(s["nnz"]), bool(with_K))


def residual(model: dict, u: np.ndarray) -> np.ndarray:
    """Global internal-force vector R(u)."""
    return _call(model, u, False)[0]


def residual_tangent_sparse(model: dict, u: np.ndarray):
    """Global R (ndarray) and consistent tangent K (scipy CSR, sorted indices)."""
    import scipy.sparse as sp

    ndof = int(model["n_dof"])
    s = model["csr"]
    R, data, _ = _call(model, u, True)
    K = sp.csr_matrix((data, s["indices"], s["indptr"]), shape=(ndof, ndof))
    K.has_sorted_indices = True
    return R, K


def group_jbar(model: dict, u: np.ndarray) -> np.ndarray:
    """Jbar of every pressure group (from the assembly kernel)."""
    return _call(model, u, False)[2]


def _qp_F(model: dict, u: np.ndarray) -> np.ndarray:
    uc = np.asarray(u, dtype=np.float64).reshape(-1, 3)
    Uq = uc[model["cell_nodes"]][model["q_cell"]]  # (n_qp, nen, 3)
    return np.eye(3)[None] + np.einsum("qia,qim->qam", Uq, model["qG"])


def displacement_at_qps(model: dict, u: np.ndarray) -> np.ndarray:
    """u_h at every quadrature point (n_qp, 3), order of ``model["qX"]``."""
    uc = np.asarray(u, dtype=np.float64).reshape(-1, 3)
    Uq = uc[model["cell_nodes"]][model["q_cell"]]
    return np.einsum("qi,qia->qa", model["qN"], Uq)


def cauchy_stress_qps(model: dict, u: np.ndarray) -> np.ndarray:
    """Cauchy stress at every quadrature point (n_qp, 6) Voigt xx, yy, zz, xy, yz, xz.

    ``sigma = mu/J B + U'(Jbar_g) I`` with the pressure of the point's group
    (numpy path, independent of the assembly kernel).
    """
    F = _qp_F(model, u)
    J = np.linalg.det(F)
    gq = model["cell_group"][model["q_cell"]]
    W = model["qW"]
    ng = int(model["n_groups"])
    Jbar = np.bincount(gq, weights=W * J, minlength=ng) / np.bincount(gq, weights=W, minlength=ng)
    mu, lam = float(model["mu"]), float(model["lam"])
    pbar = (lam * np.log(Jbar) - mu) / Jbar
    B = np.einsum("qam,qbm->qab", F, F)
    sig = (mu / J)[:, None, None] * B + pbar[gq][:, None, None] * np.eye(3)[None]
    return np.column_stack([sig[:, 0, 0], sig[:, 1, 1], sig[:, 2, 2],
                            sig[:, 0, 1], sig[:, 1, 2], sig[:, 0, 2]])


# ---------------------------------------------------------------------------
# solver
# ---------------------------------------------------------------------------

def solve_static(model: dict, f_ext: np.ndarray, fixed, n_steps: int = 1,
                 rtol: float = 1e-6, atol: float = 1e-9, max_it: int = 25,
                 prescribed: dict | None = None, dtol: float = 1e-10,
                 min_frac: float = 1e-4, history: list | None = None,
                 linear_solver: str = "auto") -> np.ndarray:
    """Adaptive load-stepped Newton (same control as ``fbar.solve_static``). Returns full u.

    Inactive control points are fixed to 0 automatically; prescribing a value
    on one raises ``ValueError`` (its basis vanishes on Omega). ``fixed`` DOFs
    stay 0; ``prescribed`` maps DOF -> final value, ramped with the load factor.

    Per increment: converged when ``|R - f| <= atol + rtol*max(|R0 - f|, |f|)``
    (max norm, free DOFs) or, after one update, ``|du| <= dtol*|u|`` (the
    residual stalls at a round-off floor at nu ~ 0.5, S40). The line search is
    a guard only (halve while the trial state has a non-finite residual,
    i.e. a group Jbar <= 0). A failed increment is restored and halved; after
    an easy one (<= 6 iterations) it grows back toward ``1/n_steps``.
    ``history`` receives ``(load_factor, iterations)`` per converged increment.
    ``linear_solver``: ``"auto"`` (pypardiso if importable, else SuperLU),
    ``"pardiso"``, ``"superlu"`` or ``"dense"`` (numpy, small models only).
    """
    ndof = int(model["n_dof"])
    u = np.zeros(ndof)
    inact = np.zeros(ndof, dtype=bool)
    inact[model["inactive_dofs"]] = True
    prescribed = {int(k): float(v) for k, v in (prescribed or {}).items()}
    for k in prescribed:
        if k < 0 or k >= ndof:
            raise ValueError(f"prescribed dof {k} out of range")
        if inact[k]:
            raise ValueError(f"prescribed dof {k} belongs to an inactive control point")
    free = ~inact
    free[np.asarray(fixed, dtype=np.int64)] = False
    if prescribed:
        free[np.fromiter(prescribed.keys(), dtype=np.int64)] = False
    free_idx = np.flatnonzero(free)
    f = np.asarray(f_ext, dtype=float)
    kind = str(linear_solver).lower()
    pardiso = None
    if kind not in ("auto", "pardiso", "superlu", "dense"):
        raise ValueError(f"unknown linear_solver {linear_solver!r}")
    if kind in ("auto", "pardiso"):
        try:
            import pypardiso
            pardiso = pypardiso.PyPardisoSolver()
        except ImportError:
            if kind == "pardiso":
                raise

    def lin_solve(K, b):
        # symmetric Jacobi scaling D A D (D = |diag A|^-1/2): removes the
        # diagonal-type ill-conditioning of basis functions with a tiny support
        # in Omega (see module docstring); the solution itself is unchanged.
        A = K[free_idx][:, free_idx].tocsr()
        d = np.abs(A.diagonal())
        if not np.all(np.isfinite(d) & (d > 0.0)):
            raise RuntimeError("fbar_trimmed: free DOF with zero or non-finite diagonal stiffness")
        s = 1.0 / np.sqrt(d)
        rows = np.repeat(np.arange(A.shape[0]), np.diff(A.indptr))
        A.data = A.data * s[rows] * s[A.indices]
        bs = s * b
        if kind == "dense":
            return s * np.linalg.solve(A.toarray(), bs)
        if pardiso is not None:
            return s * np.asarray(pypardiso.spsolve(A, bs, solver=pardiso), dtype=np.float64)
        import scipy.sparse.linalg as spla
        return s * np.asarray(spla.spsolve(A.tocsc(), bs), dtype=np.float64)

    d_nom = 1.0 / max(1, int(n_steps))
    lam, dlam = 0.0, d_nom

    def _increment(lam_new: float) -> tuple[bool, int]:
        target = f * lam_new
        for k, v in prescribed.items():
            u[k] = v * lam_new
        R, K = residual_tangent_sparse(model, u)
        rr = R[free] - target[free]
        if not np.all(np.isfinite(rr)):
            return False, 0
        gate = atol + rtol * max(float(np.max(np.abs(rr))) if rr.size else 0.0,
                                 float(np.max(np.abs(target[free]))) if rr.size else 0.0, 1e-300)
        for it in range(1, max_it + 1):
            if rr.size == 0 or float(np.max(np.abs(rr))) < gate:
                return True, it - 1
            du = lin_solve(K, -rr)
            step = 1.0
            while True:
                u_try = u.copy()
                u_try[free] += step * du
                R, K = residual_tangent_sparse(model, u_try)
                if np.all(np.isfinite(R[free])):
                    break
                step *= 0.5
                if step < 1.0 / 1024.0:
                    return False, it
            u[:] = u_try
            rr = R[free] - target[free]
            if float(np.max(np.abs(step * du))) <= dtol * max(float(np.max(np.abs(u[free]))), 1e-300):
                return True, it
        return float(np.max(np.abs(rr))) < gate, max_it

    try:
        while lam < 1.0 - 1e-14:
            lam_new = min(1.0, lam + dlam)
            u_save = u.copy()
            ok, its = _increment(lam_new)
            if ok:
                lam = lam_new
                if history is not None:
                    history.append((lam, its))
                if its <= 6:
                    dlam = min(d_nom, 1.5 * dlam)
            else:
                u[:] = u_save
                dlam *= 0.5
                if dlam < min_frac:
                    raise RuntimeError(f"fbar_trimmed Newton failed: increment below {min_frac} "
                                       f"at load factor {lam:.6g}")
    finally:
        if pardiso is not None:
            pardiso.free_memory(everything=True)
    return u
