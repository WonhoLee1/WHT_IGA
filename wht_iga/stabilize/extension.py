"""T-P2.2: extended B-spline (AgFEM-style) extension operator.

Math source: v2.1 §7.5; Höllig, Reif & Wissmann, "Weighted extended
B-spline approximation of Dirichlet problems", SIAM J. Numer. Anal. (2001);
Höllig (2003) "Finite element methods with B-splines" (extrapolation weights
as tensor products of 1D Lagrange polynomials); Marussig, Hughes et al.,
"Improved conditioning through extended B-splines" / introduction to
immersed-support classification (2018).

For each extended (slave) basis j, the extension writes its coefficient as a
linear combination of nearby master coefficients::

    u_j = sum_i e_ij u_i,   e_ij = Π_ν w^ν_{i_ν}(j_ν)

The 1D weights w^ν reproduce the B-spline coefficients of every monomial of
degree < block size on the actual knot vector (Marsden identity: coefficient of
x^k on basis i = e_k(t_{i+1}, ..., t_{i+p}) / C(p, k), the blossom). On uniform
interior knots this equals 1D Lagrange extrapolation in index space; near a
clamped boundary (p >= 2) the Greville points are non-uniform and the index-space
formula is wrong (fixed 2026-11-14, see docs/verification/V01-ext.md).
Consequences used as acceptance tests: partition of unity (each row sums to 1)
and exact reproduction of degree-p polynomials in physical coordinates.

Fallback chain (v2.1 §7.5, risk R7): full (p+1)^d master block → degree-reduced
block → nearest-master constant (AgFEM-style) → ``THIN_FEATURE`` issue with a
refinement recommendation. Coefficient magnitudes are measured and recorded;
the acceptance cap value is TBD upstream, so this module records but does NOT
invent a cap (``cap=None`` in stats).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from numba import njit

from wht_iga.stabilize.support import EXTENDED, INACTIVE, MASTER

#: Issue code for exhausted extension stencils (v2.1 §7.5, risk R7).
THIN_FEATURE = "THIN_FEATURE"


def monomial_coefficients_1d(knots: np.ndarray, p: int, i: int, n_mono: int,
                             x0: float, scale: float) -> np.ndarray:
    """B-spline coefficients on basis ``i`` of ((x - x0)/scale)^k, k < n_mono.

    Marsden/blossom: c_k(i) = e_k(tau_1..tau_p) / C(p, k) with
    tau_r = (t_{i+r} - x0) / scale (e_k = elementary symmetric polynomial).
    """
    from math import comb

    tau = (np.asarray(knots, dtype=float)[i + 1:i + p + 1] - x0) / scale
    poly = np.poly(tau) if p > 0 else np.ones(1)  # x^p - e1 x^{p-1} + e2 ...
    out = np.empty(n_mono)
    for k in range(n_mono):
        out[k] = ((-1) ** k) * poly[k] / comb(p, k)
    return out


def extrapolation_weights_1d(knots: np.ndarray, p: int, lo: int, size: int,
                             s: int) -> np.ndarray:
    """1D weights w (len ``size``) with sum_m w_m c_k(lo+m) = c_k(s), k < size.

    Exact reproduction of all polynomials of degree < size on the true knot
    vector (clamped boundaries included). Raises if the block is singular.
    """
    t = np.asarray(knots, dtype=float)
    x0 = float(np.mean(t[s + 1:s + p + 1])) if p > 0 else float(t[s])
    brk = np.unique(t)
    scale = float((brk[-1] - brk[0]) / max(1, len(brk) - 1))
    A = np.column_stack([monomial_coefficients_1d(t, p, lo + m, size, x0, scale)
                         for m in range(size)])
    b = monomial_coefficients_1d(t, p, s, size, x0, scale)
    if np.linalg.cond(A) > 1e12:
        raise ValueError(f"singular extension block lo={lo} size={size} for slave {s}")
    return np.linalg.solve(A, b)


def block_row_weights(knots: list[np.ndarray], p: int,
                      lo0: int, lo1: int, lo2: int, size: int,
                      s0: int, s1: int, s2: int) -> np.ndarray:
    """Tensor-product extrapolation row for a cubic block.

    Block lower corner (lo0, lo1, lo2), edge length ``size``; slave index
    (s0, s1, s2). Row order is flat ``da + size*(db + size*dc)``.
    """
    w0 = extrapolation_weights_1d(knots[0], p, lo0, size, s0)
    w1 = extrapolation_weights_1d(knots[1], p, lo1, size, s1)
    w2 = extrapolation_weights_1d(knots[2], p, lo2, size, s2)
    row = np.empty(size * size * size)
    q = 0
    for dc in range(size):
        for db in range(size):
            for da in range(size):
                row[q] = w0[da] * w1[db] * w2[dc]
                q += 1
    return row


def _block_all_master(mask: np.ndarray, lo: tuple[int, int, int],
                      size: int) -> bool:
    nu, nv, nw = mask.shape
    i0, j0, k0 = lo
    if i0 < 0 or j0 < 0 or k0 < 0:
        return False
    if i0 + size > nu or j0 + size > nv or k0 + size > nw:
        return False
    return bool(np.all(mask[i0:i0 + size, j0:j0 + size, k0:k0 + size]))


@njit(cache=True)
def _find_master_block_njit(
    s0: int,
    s1: int,
    s2: int,
    mask: np.ndarray,
    size: int,
    has_inward: bool,
    inw0: float,
    inw1: float,
    inw2: float,
) -> tuple[bool, int, int, int]:
    """Compiled inner search for :func:`find_master_block`."""
    nu, nv, nw = mask.shape
    found = False
    best_dist = 0.0
    best_tie = 0.0
    bi = 0
    bj = 0
    bk = 0
    half = (size - 1) / 2.0
    for k0 in range(nw - size + 1):
        for j0 in range(nv - size + 1):
            for i0 in range(nu - size + 1):
                ok = True
                for ii in range(i0, i0 + size):
                    for jj in range(j0, j0 + size):
                        for kk in range(k0, k0 + size):
                            if not mask[ii, jj, kk]:
                                ok = False
                                break
                        if not ok:
                            break
                    if not ok:
                        break
                if not ok:
                    continue
                cx0 = i0 + half - s0
                cx1 = j0 + half - s1
                cx2 = k0 + half - s2
                dist = float(np.sqrt(cx0 ** 2 + cx1 ** 2 + cx2 ** 2))
                if has_inward:
                    tie = -float(cx0 * inw0 + cx1 * inw1 + cx2 * inw2)
                else:
                    tie = float(k0 * (nv * nu) + j0 * nu + i0)
                if not found or dist < best_dist or (dist == best_dist and tie < best_tie):
                    found = True
                    best_dist = dist
                    best_tie = tie
                    bi = i0
                    bj = j0
                    bk = k0
    return found, bi, bj, bk


def find_master_block(slave: tuple[int, int, int], master_mask: np.ndarray,
                      size: int,
                      outward: np.ndarray | None = None) -> tuple[int, int, int] | None:
    """Nearest master-only cubic block of edge ``size`` for a slave index.

    Candidates are lower corners whose ``size``-brick is all-master.
    Primary key: block-center distance to the slave. Tie-break: largest
    projection of (center - slave) onto the interior direction (``-outward``);
    with ``outward=None`` the lexicographically smallest center wins. Both
    rules are deterministic; ties are measure-zero on real cut grids.

    Note on the plan wording ("outward-normal preference", v2.1 §7.5): masters
    live strictly inside the domain, so a block can never sit in the +outward
    half-space; the implementable reading is a preference for the block most
    aligned with the interior (``-outward``) among equidistant candidates.
    """
    mask = np.ascontiguousarray(master_mask, dtype=np.bool_)
    if size <= 0:
        raise ValueError(f"size must be positive, got {size}")
    inward = None if outward is None else -np.asarray(outward, dtype=float)
    has_inward = inward is not None
    inw0 = float(inward[0]) if inward is not None else 0.0
    inw1 = float(inward[1]) if inward is not None else 0.0
    inw2 = float(inward[2]) if inward is not None else 0.0
    found, bi, bj, bk = _find_master_block_njit(
        int(slave[0]),
        int(slave[1]),
        int(slave[2]),
        mask,
        int(size),
        has_inward,
        inw0,
        inw1,
        inw2,
    )
    return (int(bi), int(bj), int(bk)) if found else None


def nearest_master(slave: tuple[int, int, int],
                   master_ids: np.ndarray, n_cp: tuple[int, int, int]) -> int:
    """Flat index (into ``master_ids`` order) of the closest master basis."""
    nu, nv, _nw = (int(n) for n in n_cp)
    ids = np.asarray(master_ids, dtype=np.int64)
    d = np.abs(ids % nu - slave[0]) + np.abs((ids // nu) % nv - slave[1]) \
        + np.abs(ids // (nu * nv) - slave[2])
    return int(np.argmin(d))


@dataclass
class ExtensionOperator:
    """Discrete extension ``u_full = E @ u_master`` (CSR, sparse-only)."""

    E: sp.csr_matrix  # (n_basis, n_master)
    master_ids: np.ndarray  # flat basis ids, reduced-column order
    slave_ids: np.ndarray  # flat basis ids with extension rows
    degree_used: np.ndarray  # per slave: block degree p..1, 0 = clean constant,
    # -1 = distant constant + THIN_FEATURE issue
    issues: list[dict] = field(default_factory=list)
    stats: dict = field(default_factory=dict)


def build_extension(n_cp: tuple[int, int, int], labels: np.ndarray, p: int,
                    outward: dict[int, np.ndarray] | None = None, *,
                    knots: list[np.ndarray]) -> ExtensionOperator:
    """Build the extension operator for one classification.

    :param n_cp: control-point counts per axis.
    :param labels: flat per-basis labels (MASTER/EXTENDED/INACTIVE).
    :param p: spline degree (full block edge = p + 1).
    :param outward: optional flat-id -> outward-normal map for tie-breaks.
    :param knots: background knot vectors per axis (required: the weights
        depend on the true knots, not on the index spacing).
    """
    nu, nv, nw = (int(n) for n in n_cp)
    n_basis = nu * nv * nw
    knots = [np.asarray(k, dtype=float) for k in knots]
    for ax, (kv, n) in enumerate(zip(knots, (nu, nv, nw))):
        if len(kv) - p - 1 != n:
            raise ValueError(f"knots[{ax}] gives {len(kv) - p - 1} basis, n_cp has {n}")
    if labels.shape != (n_basis,):
        raise ValueError(f"labels shape {labels.shape} != ({n_basis},)")
    mask = (np.asarray(labels).reshape(nw, nv, nu).transpose(2, 1, 0)
            == MASTER)
    master_ids = np.flatnonzero(np.asarray(labels) == MASTER)
    slave_ids = np.flatnonzero(np.asarray(labels) == EXTENDED)
    if len(slave_ids) and len(master_ids) == 0:
        raise ValueError("extension with slaves but no masters is invalid: "
                         "refine input classification (no master basis found)")
    n_m = len(master_ids)
    col_of = {int(g): c for c, g in enumerate(master_ids)}
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    degree_used = np.full(len(slave_ids), -2, dtype=np.int64)
    issues: list[dict] = []
    max_abs = 0.0
    for r, g in enumerate(slave_ids):
        s = (int(g % nu), int((g // nu) % nv), int(g // (nu * nv)))
        nrm = None if outward is None else outward.get(int(g))
        block = None
        deg = -2
        for size in range(p + 1, 1, -1):  # full block, then degree-reduced
            blk = find_master_block(s, mask, size, nrm)
            if blk is not None:
                block, deg = blk, size - 1
                break
        if block is not None and deg >= 1:
            w = block_row_weights(knots, p, block[0], block[1], block[2],
                                  deg + 1, s[0], s[1], s[2])
            mids = []
            for dc in range(deg + 1):
                for db in range(deg + 1):
                    for da in range(deg + 1):
                        mids.append((block[0] + da)
                                    + nu * ((block[1] + db) + nv * (block[2] + dc)))
            take = np.array([col_of[m] for m in mids], dtype=np.int64)
            degree_used[r] = deg
        else:
            # Nearest-master constant (AgFEM-style). A clean single-master
            # block (deg 0) is exact; a distant master is flagged THIN_FEATURE
            # with a refinement recommendation (v2.1 §7.5, risk R7).
            c = nearest_master(s, master_ids, (nu, nv, nw))
            mflat = int(master_ids[c])
            mijk = (mflat % nu, (mflat // nu) % nv, mflat // (nu * nv))
            cheb = max(abs(mijk[0] - s[0]), abs(mijk[1] - s[1]),
                       abs(mijk[2] - s[2]))
            w = np.ones(1)
            take = np.array([c], dtype=np.int64)
            if cheb > p + 1:
                degree_used[r] = -1
                issues.append({
                    "severity": "warning", "code": THIN_FEATURE,
                    "message": ("slave basis " + str(s) + " extrapolated by a "
                                f"constant from distance {cheb} cells (> p+1 = "
                                f"{p + 1}): refine the background grid around "
                                "this feature (R7)"),
                    "entity": f"basis:{int(g)}"})
            else:
                degree_used[r] = 0
        max_abs = max(max_abs, float(np.max(np.abs(w))))
        rows.append(np.full(len(w), r, dtype=np.int64))
        cols.append(take)
        data.append(w)
    # Full E: master rows are unit, inactive rows stay zero (no node, §5.2.3/6).
    # Assembled once via COO triplets (sparse-only, no dense temporaries).
    all_rows = [np.asarray(master_ids, dtype=np.int64)]
    all_cols = [np.arange(n_m, dtype=np.int64)]
    all_data = [np.ones(n_m)]
    if len(slave_ids):
        Es = sp.coo_matrix(
            (np.concatenate(data),
             (np.concatenate(rows), np.concatenate(cols))),
            shape=(len(slave_ids), n_m))
        all_rows.append(np.asarray(slave_ids, dtype=np.int64)[Es.row])
        all_cols.append(Es.col)
        all_data.append(Es.data)
    E = sp.coo_matrix(
        (np.concatenate(all_data),
         (np.concatenate(all_rows), np.concatenate(all_cols))),
        shape=(n_basis, n_m)).tocsr()
    stats = {"max_abs_coef": float(max_abs),
             "n_constant_fallback": int(np.sum(degree_used == -1)),
             "n_degree_reduced": int(np.sum((degree_used >= 0) & (degree_used < p))),
             "n_thin_feature": len(issues),
             "cap": None,  # TBD upstream: recorded, not enforced (no invented cap)
             "n_master": int(n_m), "n_slave": len(slave_ids)}
    return ExtensionOperator(E=E, master_ids=master_ids, slave_ids=slave_ids,
                             degree_used=degree_used, issues=issues, stats=stats)
