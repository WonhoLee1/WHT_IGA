"""T-P2.1: basis-support volume ratios + master/extended/inactive classification.

Rule source: master-plan v2.1 §7.5 (NOT element-wise eta_e — basis support,
Marussig et al. 2018: a C^{p-1} basis spans (p+1)^d cells, so element-wise
small-cell tests are mathematically wrong, pitfall K20).

Definitions (v2.1 §7.5):
* large cell: eta_e >= theta (``small_cell_volume_ratio``, default 0.2,
  range [0, 0.5)).
* master: support contains >= 1 large cell.
* extended (slave): support meets Omega but has no large cell.
* inactive: support does not meet Omega.
* multi-material: support ratios use the UNION domain (material interfaces
  are interior, never extension triggers — pitfall K23).

Inputs are S3/S4-style arrays: per-cell inside volumes (S4 exact moment-fit
volumes) in cell order ``i + nx*(j + ny*k)`` (grid convention), plus the
background knot vectors. Support ranges per axis are derived from the knot
vectors (contiguous cell bricks for any clamped knot vector, any continuity).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange

#: Classification labels (ints: njit-friendly).
INACTIVE = 0
EXTENDED = 1
MASTER = 2

#: Default small-cell threshold (v2.1 §7.5, §0.3 E6).
DEFAULT_THETA = 0.2
#: Allowed theta range (v2.1 §7.5).
THETA_MIN, THETA_MAX_EXCL = 0.0, 0.5


@njit(cache=True)
def _accumulate(n_basis: int, s0: np.ndarray, s1: np.ndarray,
                cell_inside: np.ndarray, cell_vol: np.ndarray,
                nx: int, ny: int, cell_large: np.ndarray,
                ratios: np.ndarray, has_large: np.ndarray,
                has_inside: np.ndarray) -> None:
    """Per-basis support sums over the cell brick s0[b,:]..s1[b,:] (inclusive)."""
    for b in range(n_basis):
        vin = 0.0
        vtot = 0.0
        large = False
        inside = False
        for k in range(s0[b, 2], s1[b, 2] + 1):
            for j in range(s0[b, 1], s1[b, 1] + 1):
                for i in range(s0[b, 0], s1[b, 0] + 1):
                    c = i + nx * (j + ny * k)
                    vin += cell_inside[c]
                    vtot += cell_vol[c]
                    if cell_large[c]:
                        large = True
                    if cell_inside[c] > 0.0:
                        inside = True
        ratios[b] = vin / vtot if vtot > 0.0 else 0.0
        has_large[b] = large
        has_inside[b] = inside


@njit(cache=True, parallel=True)
def _accumulate_tensor_bricks_njit(nu: int, nv: int, nw: int,
                                   e0x: np.ndarray, e1x: np.ndarray,
                                   e0y: np.ndarray, e1y: np.ndarray,
                                   e0z: np.ndarray, e1z: np.ndarray,
                                   cell_inside: np.ndarray, cell_vol: np.ndarray,
                                   nx: int, ny: int, cell_large: np.ndarray,
                                   ratios: np.ndarray, has_large: np.ndarray,
                                   has_inside: np.ndarray) -> None:
    """Parallel per-basis support sums directly from 1D axis ranges (no s0/s1 alloc)."""
    n_basis = nu * nv * nw
    for b in prange(n_basis):
        c = b // (nu * nv)
        rem = b - c * (nu * nv)
        a2 = rem // nu
        a1 = rem - a2 * nu
        i0, i1 = e0x[a1], e1x[a1]
        j0, j1 = e0y[a2], e1y[a2]
        k0, k1 = e0z[c], e1z[c]
        vin = 0.0
        vtot = 0.0
        large = False
        inside = False
        for k in range(k0, k1 + 1):
            base_k = ny * k
            for j in range(j0, j1 + 1):
                base_j = nx * (j + base_k)
                for i in range(i0, i1 + 1):
                    cell_id = i + base_j
                    vi = cell_inside[cell_id]
                    vin += vi
                    vtot += cell_vol[cell_id]
                    if cell_large[cell_id]:
                        large = True
                    if vi > 0.0:
                        inside = True
        ratios[b] = vin / vtot if vtot > 0.0 else 0.0
        has_large[b] = large
        has_inside[b] = inside


def _support_ranges_1d(knots: np.ndarray, breaks: np.ndarray, p: int,
                       n_cp: int) -> tuple[np.ndarray, np.ndarray]:
    """Contiguous support cell range [e0, e1] (inclusive) per 1D basis index.

    Cell e is in the support iff [breaks[e], breaks[e+1]] meets
    (knots[i], knots[i+p+1]) with positive length.
    """
    knots = np.asarray(knots, dtype=float)
    breaks = np.asarray(breaks, dtype=float)
    nc = len(breaks) - 1
    e0 = np.empty(n_cp, dtype=np.int64)
    e1 = np.empty(n_cp, dtype=np.int64)
    for i in range(n_cp):
        klo, khi = knots[i], knots[i + p + 1]
        a = int(np.searchsorted(breaks, klo, side="right") - 1)
        b = int(np.searchsorted(breaks, khi, side="left") - 1)
        # Clamp knot-exact endpoints onto the grid (clamped ends, roundoff).
        if a < 0:
            a = 0
        if b > nc - 1:
            b = nc - 1
        if b < a:  # degenerate (zero-measure support): pin to nearest cell
            a = min(max(a, 0), nc - 1)
            b = a
        e0[i], e1[i] = a, b
    return e0, e1


@dataclass
class SupportClassification:
    """Result of :func:`classify_basis`."""

    labels: np.ndarray  # (n_basis,) ints: INACTIVE/EXTENDED/MASTER
    ratios: np.ndarray  # (n_basis,) support inside-volume fractions
    theta: float
    multi_material: bool  # True only if per-material etas were unioned
    n_master: int = 0
    n_extended: int = 0
    n_inactive: int = 0
    notes: list[str] = field(default_factory=list)


def union_inside_volume(material_etas: list[np.ndarray],
                        cell_vol: np.ndarray) -> np.ndarray:
    """Union-domain inside volume per cell from per-material fractions.

    Materials form a disjoint partition here, so the union fraction is the
    sum clipped to [0, 1]; sums above 1 + tol are rejected (overlapping
    material input is a caller bug, not a union).
    """
    etas = [np.asarray(e, dtype=float) for e in material_etas]
    tot = np.zeros_like(etas[0])
    for e in etas:
        if e.shape != tot.shape:
            raise ValueError(f"material eta shape {e.shape} != {tot.shape}")
        if np.any(e < 0.0):
            raise ValueError("material eta fractions must be >= 0")
        tot = tot + e
    if np.any(tot > 1.0 + 1e-12):
        raise ValueError(
            f"material etas overlap: max sum = {tot.max()}, expected <= 1")
    return np.clip(tot, 0.0, 1.0) * np.asarray(cell_vol, dtype=float)


def classify_basis(n_cp: tuple[int, int, int], knots: list[np.ndarray],
                   cell_inside_vol: np.ndarray, cell_vol: np.ndarray,
                   theta: float = DEFAULT_THETA,
                   material_etas: list[np.ndarray] | None = None,
                   large_tol: float = 0.0) -> SupportClassification:
    """Classify every tensor-product basis function.

    :param n_cp: control-point counts per axis.
    :param knots: background knot vectors per axis (standard clamped).
    :param cell_inside_vol: inside volume per cell (union domain), flat
        cell order ``i + nx*(j + ny*k)``.
    :param cell_vol: geometric volume per cell (same order).
    :param theta: large-cell threshold in [0, 0.5).
    :param material_etas: optional per-material inside fractions; when given
        (length >= 2) the union branch runs, else single-material is used and
        the multi-material branch is reported as not-run (never a stub).
    """
    if not (THETA_MIN <= theta < THETA_MAX_EXCL):
        raise ValueError(
            f"theta={theta} outside [{THETA_MIN}, {THETA_MAX_EXCL})")
    n_cp_t = tuple(int(n) for n in n_cp)
    notes: list[str] = []
    inside = np.asarray(cell_inside_vol, dtype=float)
    vol = np.asarray(cell_vol, dtype=float)
    multi = material_etas is not None and len(material_etas) >= 2
    if multi:
        inside = union_inside_volume(list(material_etas), vol)
        notes.append(f"multi-material union branch: {len(material_etas)} materials")
    else:
        notes.append("multi-material branch not-run: single-material input "
                     "(no multi-material cells present)")
    nx = len(np.unique(knots[0])) - 1
    # Per-axis support ranges, then flat basis loop over the brick product.
    ranges = []
    for ax in range(3):
        br = np.unique(np.asarray(knots[ax], dtype=float))
        p_ax = len(knots[ax]) - n_cp_t[ax] - 1
        ranges.append(_support_ranges_1d(knots[ax], br, p_ax, n_cp_t[ax]))
    nu, nv, nw = n_cp_t
    n_basis = nu * nv * nw
    eta = inside / np.maximum(vol, 1e-300)
    cell_large = eta >= theta - large_tol
    ratios = np.empty(n_basis)
    has_large = np.empty(n_basis, dtype=np.bool_)
    has_inside = np.empty(n_basis, dtype=np.bool_)
    ny = len(np.unique(knots[1])) - 1
    _accumulate_tensor_bricks_njit(
        nu, nv, nw,
        ranges[0][0], ranges[0][1],
        ranges[1][0], ranges[1][1],
        ranges[2][0], ranges[2][1],
        inside, vol, nx, ny, cell_large,
        ratios, has_large, has_inside,
    )
    labels = np.where(has_large, MASTER,
                      np.where(has_inside, EXTENDED, INACTIVE)).astype(np.int64)
    return SupportClassification(
        labels=labels, ratios=ratios, theta=float(theta),
        multi_material=bool(multi),
        n_master=int(np.sum(labels == MASTER)),
        n_extended=int(np.sum(labels == EXTENDED)),
        n_inactive=int(np.sum(labels == INACTIVE)),
        notes=notes)
