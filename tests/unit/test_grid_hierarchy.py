"""Unit tests for THB-spline hierarchical refinement (wht_iga.grid.hierarchy)."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.grid.bspline import basis_1d, ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.grid.hierarchy import (
    THBHierarchy,
    midpoint_refine_knots,
    subdivision_matrix_1d,
    subdivision_matrix_3d,
)
from wht_iga.grid.knots import n_control_points, uniform_clamped


def _eval_1d_all(knots: np.ndarray, p: int, u: float) -> np.ndarray:
    """Evaluate all 1D B-spline basis functions at ``u``."""
    n_cp = n_control_points(knots, p)
    out = np.zeros(n_cp, dtype=float)
    span, ders = basis_1d(knots, p, u, n_ders=0)
    out[span - p : span + 1] = ders[0]
    return out


def _eval_3d_all(
    ku: np.ndarray,
    kv: np.ndarray,
    kw: np.ndarray,
    pu: int,
    pv: int,
    pw: int,
    pt: tuple[float, float, float],
) -> np.ndarray:
    """Evaluate all 3D tensor-product B-spline basis functions at ``pt``."""
    bu = _eval_1d_all(ku, pu, pt[0])
    bv = _eval_1d_all(kv, pv, pt[1])
    bw = _eval_1d_all(kw, pw, pt[2])
    # Ordering: i + nu * (j + nv * k) == Fortran-order flattening.
    return (bu[:, None, None] * bv[None, :, None] * bw[None, None, :]).ravel(order="F")


def _boxes_touch(
    a_lo: tuple[float, float, float],
    a_hi: tuple[float, float, float],
    b_lo: tuple[float, float, float],
    b_hi: tuple[float, float, float],
    tol: float = 1e-12,
) -> bool:
    """Return True if two closed 3D boxes share at least a vertex, edge, or face."""
    for d in range(3):
        if a_hi[d] < b_lo[d] - tol or b_hi[d] < a_lo[d] - tol:
            return False
    return True


@pytest.mark.parametrize("degree", [1, 2, 3, 4])
def test_1d_dyadic_subdivision_accuracy(degree: int):
    """N^L(u) == S_1D^T N^{L+1}(u) to machine precision (< 1e-14)."""
    u_coarse = uniform_clamped(0.0, 2.5, n_cells=4, degree=degree)
    u_fine = midpoint_refine_knots(u_coarse)
    s_1d = subdivision_matrix_1d(u_coarse, u_fine, degree)

    pts = np.linspace(0.0, 2.5, 41)
    max_err = 0.0
    for u in pts:
        n_c = _eval_1d_all(u_coarse, degree, float(u))
        n_f = _eval_1d_all(u_fine, degree, float(u))
        err = float(np.max(np.abs(n_c - s_1d.T @ n_f)))
        max_err = max(max_err, err)

    assert max_err < 1e-14


def test_1d_subdivision_nonuniform_knots():
    """Oslo subdivision holds on non-uniform clamped knot vectors."""
    p = 3
    u_coarse = np.array([0.0, 0.0, 0.0, 0.0, 0.2, 0.75, 1.0, 1.0, 1.0, 1.0])
    u_fine = midpoint_refine_knots(u_coarse)
    s_1d = subdivision_matrix_1d(u_coarse, u_fine, p)

    for u in np.linspace(0.0, 1.0, 29):
        n_c = _eval_1d_all(u_coarse, p, float(u))
        n_f = _eval_1d_all(u_fine, p, float(u))
        assert np.max(np.abs(n_c - s_1d.T @ n_f)) < 1e-14


def test_3d_subdivision_operator_accuracy():
    """N_3D^L(u,v,w) == S_3D^T N_3D^{L+1}(u,v,w) to < 1e-14."""
    pu, pv, pw = 2, 2, 3
    ku0 = uniform_clamped(0.0, 1.0, 2, pu)
    kv0 = uniform_clamped(-1.0, 1.0, 2, pv)
    kw0 = uniform_clamped(0.0, 2.0, 2, pw)

    ku1 = midpoint_refine_knots(ku0)
    kv1 = midpoint_refine_knots(kv0)
    kw1 = midpoint_refine_knots(kw0)

    su = subdivision_matrix_1d(ku0, ku1, pu)
    sv = subdivision_matrix_1d(kv0, kv1, pv)
    sw = subdivision_matrix_1d(kw0, kw1, pw)
    s3d = subdivision_matrix_3d(su, sv, sw)

    rng = np.random.default_rng(42)
    for _ in range(15):
        u = float(rng.uniform(0.0, 1.0))
        v = float(rng.uniform(-1.0, 1.0))
        w = float(rng.uniform(0.0, 2.0))
        n_c = _eval_3d_all(ku0, kv0, kw0, pu, pv, pw, (u, v, w))
        n_f = _eval_3d_all(ku1, kv1, kw1, pu, pv, pw, (u, v, w))
        assert np.max(np.abs(n_c - s3d.T @ n_f)) < 1e-14


@pytest.mark.parametrize(
    "mode,target_level",
    [
        ("corner", 2),
        ("corner", 3),
        ("center", 3),
    ],
)
def test_2_to_1_balance_propagation(mode: str, target_level: int):
    """Every face/edge/vertex-adjacent leaf cell pair satisfies |L_a - L_b| <= 1."""
    hier = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n_cells=2, degree=2)
    if mode == "corner":
        hier.mark_box((0.0, 0.0, 0.0), (0.1, 0.1, 0.1), target_level=target_level)
    else:
        hier.mark_sphere((0.5, 0.5, 0.5), radius=0.05, target_level=target_level)

    hier.enforce_2_to_1_balance()
    leaves = hier.leaf_cells()

    # Check total volume equals 1.0.
    total_vol = sum(
        (c.bounds_u[1] - c.bounds_u[0])
        * (c.bounds_v[1] - c.bounds_v[0])
        * (c.bounds_w[1] - c.bounds_w[0])
        for c in leaves
    )
    assert abs(total_vol - 1.0) < 1e-13

    # Check all touching leaf pairs for |L_a - L_b| <= 1.
    for i, ca in enumerate(leaves):
        a_lo = (ca.bounds_u[0], ca.bounds_v[0], ca.bounds_w[0])
        a_hi = (ca.bounds_u[1], ca.bounds_v[1], ca.bounds_w[1])
        for cb in leaves[i + 1 :]:
            if abs(ca.level - cb.level) <= 1:
                continue
            b_lo = (cb.bounds_u[0], cb.bounds_v[0], cb.bounds_w[0])
            b_hi = (cb.bounds_u[1], cb.bounds_v[1], cb.bounds_w[1])
            assert not _boxes_touch(a_lo, a_hi, b_lo, b_hi), (
                f"2:1 balance violated between {ca} and {cb}"
            )


def test_mark_cells_explicit():
    """Explicit cell marking refines target cells and balances cleanly."""
    hier = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n_cells=2, degree=1)
    added = hier.mark_cells([(0, 0, 0, 0)], target_level=2)
    assert added > 0
    hier.enforce_2_to_1_balance()
    assert hier.max_level == 2
    assert len(hier.leaf_cells()) > 8


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize("levels", [2, 3])
def test_partition_of_unity(degree: int, levels: int):
    """Sum of THB basis functions == 1.0 (< 1e-12) and non-negative everywhere."""
    hier = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n_cells=2, degree=degree)
    hier.mark_box((0.0, 0.0, 0.0), (0.4, 0.4, 0.4), target_level=levels - 1)
    hier.enforce_2_to_1_balance()

    grid_1d = np.linspace(0.0, 1.0, 5)
    uu, vv, ww = np.meshgrid(grid_1d, grid_1d, grid_1d, indexing="ij")
    grid_pts = np.column_stack([uu.ravel(), vv.ravel(), ww.ravel()])

    rng = np.random.default_rng(123 + degree + levels)
    rand_pts = rng.uniform(0.0, 1.0, size=(80, 3))
    pts = np.vstack([grid_pts, rand_pts])

    vals = hier.evaluate_thb_basis(pts)
    row_sums = vals.sum(axis=1)
    assert np.max(np.abs(row_sums - 1.0)) < 1e-12
    assert np.min(vals) >= -1e-14


@pytest.mark.parametrize("degree", [2, 3])
def test_polynomial_reproduction(degree: int):
    """Reproduce f(u,v,w) = 1 + 2u - 3v + 0.5 w^p + u*v to < 1e-10 on 3-level THB."""
    hier = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n_cells=2, degree=degree)
    hier.mark_sphere((0.35, 0.35, 0.35), radius=0.25, target_level=2)
    hier.enforce_2_to_1_balance()

    def poly(u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        return 1.0 + 2.0 * u - 3.0 * v + 0.5 * (w**degree) + u * v

    def poly_grad(u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        df_du = 2.0 + v
        df_dv = -3.0 + u
        df_dw = 0.5 * degree * (w ** (degree - 1))
        return np.column_stack([df_du, df_dv, df_dw])

    # Compute exact level-0 B-spline coefficients by Greville collocation on V^0.
    ku0, kv0, kw0 = hier._knots_u[0], hier._knots_v[0], hier._knots_w[0]
    gu = greville(ku0, degree)
    gv = greville(kv0, degree)
    gw = greville(kw0, degree)
    ug, vg, wg = np.meshgrid(gu, gv, gw, indexing="ij")
    grev_pts = np.column_stack(
        [ug.ravel(order="F"), vg.ravel(order="F"), wg.ravel(order="F")]
    )

    n0 = len(grev_pts)
    colloc = np.zeros((n0, n0), dtype=float)
    for r, pt in enumerate(grev_pts):
        colloc[r, :] = _eval_3d_all(
            ku0, kv0, kw0, degree, degree, degree, (float(pt[0]), float(pt[1]), float(pt[2]))
        )
    rhs = poly(grev_pts[:, 0], grev_pts[:, 1], grev_pts[:, 2])
    c0 = np.linalg.solve(colloc, rhs)

    # Project coarse coefficients to THB space.
    c_thb = hier.project_coarse_coefficients(c0)

    rng = np.random.default_rng(77 + degree)
    eval_pts = rng.uniform(0.0, 1.0, size=(100, 3))
    vals, grads = hier.evaluate_thb_basis_and_grad(eval_pts)

    f_approx = vals @ c_thb
    f_exact = poly(eval_pts[:, 0], eval_pts[:, 1], eval_pts[:, 2])
    assert np.max(np.abs(f_approx - f_exact)) < 1e-10

    g_approx = np.einsum("pjk,j->pk", grads, c_thb)
    g_exact = poly_grad(eval_pts[:, 0], eval_pts[:, 1], eval_pts[:, 2])
    assert np.max(np.abs(g_approx - g_exact)) < 1e-10


def test_c_thb_matrix_consistency_and_independence():
    """C_THB matches direct fine-basis evaluation and has no all-zero columns."""
    hier = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n_cells=2, degree=2)
    hier.mark_box((0.5, 0.5, 0.5), (1.0, 1.0, 1.0), target_level=2)
    hier.enforce_2_to_1_balance()

    c_thb = hier.build_c_thb_matrix()
    assert sp.isspmatrix_csr(c_thb)
    assert c_thb.shape == (hier.n_fine_dofs, hier.n_thb_dofs)

    # Every active THB DOF must have a non-zero column in C_THB.
    col_norms = np.asarray(c_thb.power(2).sum(axis=0)).ravel()
    assert np.all(col_norms > 1e-14)

    # Evaluate using fine B-spline basis directly and compare against evaluate_thb_basis.
    l_max = hier.max_level
    ku, kv, kw = hier._knots_u[l_max], hier._knots_v[l_max], hier._knots_w[l_max]
    pu, pv, pw = hier.degrees

    rng = np.random.default_rng(999)
    pts = rng.uniform(0.0, 1.0, size=(30, 3))
    vals_api = hier.evaluate_thb_basis(pts)

    vals_manual = np.zeros_like(vals_api)
    for idx, pt in enumerate(pts):
        n_fine = _eval_3d_all(ku, kv, kw, pu, pv, pw, (float(pt[0]), float(pt[1]), float(pt[2])))
        vals_manual[idx, :] = c_thb.T @ n_fine

    assert np.max(np.abs(vals_api - vals_manual)) < 1e-12
