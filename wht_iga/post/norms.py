"""Quadrature-based L2, H1, and energy error norm calculators.

Pure NumPy/SciPy implementation operating on quadrature point arrays or
analytic callbacks over brick (hexahedral) and tetrahedral cells.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class NormResult:
    """Absolute, reference, and relative error norms."""

    abs_error: float
    ref_norm: float
    rel_error: float


@dataclass(frozen=True)
class H1NormResult:
    """Combined L2, H1 seminorm, and full H1 error norm."""

    l2: NormResult
    h1_semi: NormResult
    h1: NormResult


def _finalize_norm(err_sq: float, ref_sq: float) -> NormResult:
    abs_err = float(np.sqrt(max(err_sq, 0.0)))
    ref_norm = float(np.sqrt(max(ref_sq, 0.0)))
    if ref_norm == 0.0:
        rel_err = 0.0 if abs_err == 0.0 else float("inf")
    else:
        rel_err = abs_err / ref_norm
    return NormResult(abs_error=abs_err, ref_norm=ref_norm, rel_error=rel_err)


def _sum_sq_components(arr: np.ndarray, n_qp: int) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    if a.shape[0] != n_qp:
        raise ValueError(f"Expected leading dimension n_qp={n_qp}, got {a.shape}")
    if a.ndim == 1:
        return a * a
    return np.sum(a * a, axis=tuple(range(1, a.ndim)))


def compute_l2_norm(
    values_num: np.ndarray,
    values_ref: np.ndarray,
    weights: np.ndarray,
    jac_det: np.ndarray | float = 1.0,
) -> NormResult:
    """Compute quadrature-based ``L2`` error norm.

    Parameters
    ----------
    values_num:
        Numerical field values at quadrature points, shape ``(n_qp, ...)``.
    values_ref:
        Reference/analytic field values at quadrature points, shape ``(n_qp, ...)``.
    weights:
        Quadrature weights of shape ``(n_qp,)``.
    jac_det:
        Jacobian determinants of shape ``(n_qp,)`` or scalar ``1.0`` if already
        folded into ``weights``.
    """
    w = np.asarray(weights, dtype=np.float64).ravel()
    n_qp = w.shape[0]
    j = np.broadcast_to(np.asarray(jac_det, dtype=np.float64), (n_qp,))
    wj = w * np.abs(j)

    v_num = np.asarray(values_num, dtype=np.float64)
    v_ref = np.asarray(values_ref, dtype=np.float64)
    diff = v_num - v_ref

    err_sq = float(np.sum(wj * _sum_sq_components(diff, n_qp)))
    ref_sq = float(np.sum(wj * _sum_sq_components(v_ref, n_qp)))
    return _finalize_norm(err_sq, ref_sq)


def compute_h1_norm(
    u_num: np.ndarray,
    grad_num: np.ndarray,
    u_ref: np.ndarray,
    grad_ref: np.ndarray,
    weights: np.ndarray,
    jac_det: np.ndarray | float = 1.0,
) -> H1NormResult:
    """Compute ``L2``, ``H1`` seminorm, and full ``H1`` error norm.

    Parameters
    ----------
    u_num, u_ref:
        Field values at quadrature points, shape ``(n_qp, ...)``.
    grad_num, grad_ref:
        Spatial gradients at quadrature points, shape ``(n_qp, ..., dim)``.
    weights:
        Quadrature weights of shape ``(n_qp,)``.
    jac_det:
        Jacobian determinants of shape ``(n_qp,)`` or scalar.
    """
    l2 = compute_l2_norm(u_num, u_ref, weights, jac_det)
    h1_semi = compute_l2_norm(grad_num, grad_ref, weights, jac_det)
    h1 = _finalize_norm(
        l2.abs_error**2 + h1_semi.abs_error**2,
        l2.ref_norm**2 + h1_semi.ref_norm**2,
    )
    return H1NormResult(l2=l2, h1_semi=h1_semi, h1=h1)


def compute_energy_norm(
    strain_num: np.ndarray,
    stress_num: np.ndarray,
    strain_ref: np.ndarray,
    stress_ref: np.ndarray,
    weights: np.ndarray,
    jac_det: np.ndarray | float = 1.0,
) -> NormResult:
    """Compute elastic energy error norm ``sqrt(int (e - e_ref):(s - s_ref) dOmega)``.

    Supports Voigt vectors ``(n_qp, 6)`` / ``(n_qp, 3)`` or full tensors
    ``(n_qp, 3, 3)`` / ``(n_qp, 2, 2)``.
    """
    w = np.asarray(weights, dtype=np.float64).ravel()
    n_qp = w.shape[0]
    j = np.broadcast_to(np.asarray(jac_det, dtype=np.float64), (n_qp,))
    wj = w * np.abs(j)

    e_num = np.asarray(strain_num, dtype=np.float64)
    s_num = np.asarray(stress_num, dtype=np.float64)
    e_ref = np.asarray(strain_ref, dtype=np.float64)
    s_ref = np.asarray(stress_ref, dtype=np.float64)

    de = e_num - e_ref
    ds = s_num - s_ref
    axes = tuple(range(1, de.ndim))
    err_density = np.sum(de * ds, axis=axes)
    ref_density = np.sum(e_ref * s_ref, axis=axes)

    err_sq = float(np.sum(wj * err_density))
    ref_sq = float(np.sum(wj * ref_density))
    return _finalize_norm(err_sq, ref_sq)


def gauss_legendre_hex_rule(
    n_gp: int | tuple[int, int, int] = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return 3D Gauss-Legendre points on ``[-1, 1]^3`` and tensor-product weights."""
    if isinstance(n_gp, int):
        counts = (n_gp, n_gp, n_gp)
    else:
        counts = (int(n_gp[0]), int(n_gp[1]), int(n_gp[2]))
    pts_1d = []
    wts_1d = []
    for n in counts:
        p, w = np.polynomial.legendre.leggauss(n)
        pts_1d.append(p)
        wts_1d.append(w)

    X, Y, Z = np.meshgrid(pts_1d[0], pts_1d[1], pts_1d[2], indexing="ij")
    WX, WY, WZ = np.meshgrid(wts_1d[0], wts_1d[1], wts_1d[2], indexing="ij")
    xi = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    weights = (WX * WY * WZ).ravel()
    return xi, weights


def integrate_error_against_analytic(
    cell_bounds: np.ndarray,
    eval_num_fn: Callable[[np.ndarray], np.ndarray],
    analytic_fn: Callable[[np.ndarray], np.ndarray],
    n_gp: int | tuple[int, int, int] = 3,
    eval_grad_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    analytic_grad_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> NormResult | H1NormResult:
    """Integrate ``L2`` (or ``H1`` when gradients are provided) error over brick cells.

    Parameters
    ----------
    cell_bounds:
        Array of shape ``(n_cells, 2, 3)`` or ``(n_cells, 6)`` giving ``[lo, hi]``.
    eval_num_fn:
        Callable mapping ``pts (n_qp, 3) -> values (n_qp, ...)``.
    analytic_fn:
        Analytic reference callable mapping ``pts (n_qp, 3) -> values (n_qp, ...)``.
    n_gp:
        Gauss-Legendre points per axis.
    eval_grad_fn, analytic_grad_fn:
        Optional gradient callables ``pts (n_qp, 3) -> grad (n_qp, ..., 3)``.
        When both are provided, returns :class:`H1NormResult`.
    """
    bounds = np.asarray(cell_bounds, dtype=np.float64)
    if bounds.ndim == 2 and bounds.shape[1] == 6:
        bounds = bounds.reshape((-1, 2, 3))
    if bounds.ndim != 3 or bounds.shape[1:] != (2, 3):
        raise ValueError(f"cell_bounds must have shape (N, 2, 3) or (N, 6), got {bounds.shape}")

    xi, w_ref = gauss_legendre_hex_rule(n_gp)
    n_cells = bounds.shape[0]
    n_local = xi.shape[0]

    lo = bounds[:, 0, :]
    hi = bounds[:, 1, :]
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    pts = mid[:, None, :] + xi[None, :, :] * half[:, None, :]
    pts_flat = pts.reshape((-1, 3))

    jac_cell = half[:, 0] * half[:, 1] * half[:, 2]
    jac_flat = np.repeat(jac_cell, n_local)
    w_flat = np.tile(w_ref, n_cells)

    u_num = eval_num_fn(pts_flat)
    u_ref = analytic_fn(pts_flat)

    if eval_grad_fn is not None and analytic_grad_fn is not None:
        g_num = eval_grad_fn(pts_flat)
        g_ref = analytic_grad_fn(pts_flat)
        return compute_h1_norm(u_num, g_num, u_ref, g_ref, w_flat, jac_flat)

    return compute_l2_norm(u_num, u_ref, w_flat, jac_flat)


def integrate_error_on_tets(
    points: np.ndarray,
    tets: np.ndarray,
    eval_num_fn: Callable[[np.ndarray], np.ndarray],
    analytic_fn: Callable[[np.ndarray], np.ndarray],
) -> NormResult:
    """Integrate ``L2`` error over a tetrahedral mesh using a degree-2 4-point rule.

    Parameters
    ----------
    points:
        Vertex coordinates of shape ``(n_pts, 3)``.
    tets:
        Tetrahedron connectivity of shape ``(n_tets, 4)``.
    eval_num_fn, analytic_fn:
        Callables mapping ``(n_qp, 3) -> (n_qp, ...)``.
    """
    pts_arr = np.asarray(points, dtype=np.float64)
    tets_arr = np.asarray(tets, dtype=np.int64)
    a = (5.0 - np.sqrt(5.0)) / 20.0
    b = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
    bary = np.array(
        [
            [b, a, a, a],
            [a, b, a, a],
            [a, a, b, a],
            [a, a, a, b],
        ],
        dtype=np.float64,
    )  # weights sum to 1 (each 1/4 of tet volume)

    verts = pts_arr[tets_arr]  # (n_tets, 4, 3)
    v0 = verts[:, 0, :]
    d1 = verts[:, 1, :] - v0
    d2 = verts[:, 2, :] - v0
    d3 = verts[:, 3, :] - v0
    vols = np.abs(np.einsum("ij,ij->i", d1, np.cross(d2, d3))) / 6.0

    qp = np.einsum("qj,tjd->tqd", bary, verts).reshape((-1, 3))
    weights = np.repeat(vols * 0.25, 4)

    return compute_l2_norm(eval_num_fn(qp), analytic_fn(qp), weights, jac_det=1.0)
