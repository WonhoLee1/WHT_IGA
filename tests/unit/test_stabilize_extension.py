"""T-P2.2 acceptance: partition of unity, p-reproduction, fallback chain (Kratos-free).

Math oracle: the extension reproduces degree-p polynomials in PHYSICAL
coordinates on the true (clamped) knot vectors and has unit row sums
(Höllig 2003; Marussig, Hughes et al. 2018 — cited in ``extension.py``).
Polynomial coefficients come from Greville collocation (``_coef1d``), not from
the extension code. Index-space polynomials are not B-spline polynomials near a
clamped boundary (V01-ext defect, 2026-11-14).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.grid.knots import uniform_clamped
from wht_iga.stabilize import extension as E
from wht_iga.stabilize import support as S


def _kn(ncp, p):
    """Uniform clamped knots with unit cells matching ``ncp`` control points."""
    return [uniform_clamped(0.0, float(n - p), n - p, p) for n in ncp]


def _coef1d(knots, p, f):
    """B-spline coefficients of f by Greville collocation (independent oracle)."""
    from wht_iga.grid.bspline import basis_1d
    from wht_iga.grid.greville import greville

    g = greville(knots, p)
    M = np.zeros((len(g), len(g)))
    for r, x in enumerate(g):
        span, N_ = basis_1d(knots, p, float(x), 0)
        M[r, span - p:span + 1] = np.atleast_2d(N_)[0]
    return np.linalg.solve(M, f(g))


def _hole_labels(p=2, theta=0.2):
    nc = (8, 8, 2)
    knots = [uniform_clamped(0.0, 8.0, n, p) for n in nc]
    ncp = tuple(len(k) - p - 1 for k in knots)
    br = [np.unique(k) for k in knots]
    vols = [np.diff(b) for b in br]
    cv = (vols[0][:, None, None] * vols[1][None, :, None]
          * vols[2][None, None, :]).ravel()
    eta = np.ones(nc[0] * nc[1] * nc[2])
    nx, ny, nz = nc
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                d = (i - 3.5) ** 2 + (j - 3.5) ** 2
                c = i + nx * (j + ny * k)
                if d < 2.5 ** 2:
                    eta[c] = 0.0
                elif d < 3.5 ** 2:
                    eta[c] = 0.1
    out = S.classify_basis(ncp, knots, eta * cv, cv, theta=theta)
    return ncp, out.labels


def test_partition_of_unity_1e12():
    ncp, labels = _hole_labels()
    op = E.build_extension(ncp, labels, p=2, knots=_kn(ncp, 2))
    assert len(op.slave_ids) > 0
    rowsums = np.asarray(op.E[op.slave_ids, :].sum(axis=1)).ravel()
    np.testing.assert_allclose(rowsums, 1.0, rtol=0, atol=1e-12)
    # master rows are unit, inactive rows are zero
    assert op.E[op.master_ids, :].nnz == len(op.master_ids)
    inactive = np.flatnonzero(labels == S.INACTIVE)
    assert op.E[inactive, :].nnz == 0


@pytest.mark.parametrize("p", [1, 2, 3])
def test_degree_reproduction_1e10(p):
    """Slaves on the a == 0 clamped face: physical degree-p polynomials exact."""
    n = 7
    ncp = (n, n, n)
    kn = _kn(ncp, p)
    labels = np.full(n ** 3, S.MASTER)
    nu = n
    for c in range(n):
        for b in range(n):
            labels[(0) + nu * (b + nu * c)] = S.EXTENDED  # a == 0 face
    op = E.build_extension(ncp, labels, p=p, knots=kn)
    rng = np.random.default_rng(100 + p)
    cx = [_coef1d(kn[0], p, lambda x, k=k: x ** k) for k in range(p + 1)]
    cy = [_coef1d(kn[1], p, lambda x, k=k: x ** k) for k in range(p + 1)]
    cz = [_coef1d(kn[2], p, lambda x, k=k: x ** k) for k in range(p + 1)]
    A, B, Cc = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    u = np.zeros((n, n, n))
    for i in range(p + 1):
        for j in range(p + 1 - i):
            for k in range(p + 1 - i - j):
                u = u + rng.uniform(-1, 1) * cx[i][A] * cy[j][B] * cz[k][Cc]
    u = u.transpose(2, 1, 0).reshape(-1)  # flat a + nu*(b + nv*c)
    pred = op.E @ u[op.master_ids]
    np.testing.assert_allclose(pred, u, rtol=0, atol=1e-10)


@pytest.mark.parametrize("p", [2, 3])
def test_clamped_boundary_affine_regression(p):
    """V01-ext regression: corner-column slaves near a clamped boundary.

    The former index-space weights gave slave errors 0.125 (p=2) / 0.083 (p=3)
    on u = y here. Physical affine coefficients = Greville points.
    """
    from wht_iga.reference.box import BoxSpline

    n = 4
    box = BoxSpline(lo=(0, 0, 0), hi=(1, 1, 1), n_cells=(n, n, n), degree=p)
    kn = [np.asarray(k, float) for k in box.knots]
    eta = np.ones(n ** 3)
    for k in range(n):
        eta[(n - 1) + n * ((n - 1) + n * k)] = 0.02
    cv = 1.0 / n ** 3
    cls = S.classify_basis(tuple(box.n_cp), kn, eta * cv, np.full(n ** 3, cv), theta=0.2)
    op = E.build_extension(tuple(box.n_cp), cls.labels, p=p, knots=kn)
    assert len(op.slave_ids) > 0
    X = box.control_points()
    for d in range(3):
        f = X[:, d]
        np.testing.assert_allclose(op.E @ f[op.master_ids], f, rtol=0, atol=1e-12)


def test_interior_matches_index_lagrange():
    """On uniform interior knots the weights equal index-space Lagrange."""
    p = 2
    kn = uniform_clamped(0.0, 10.0, 10, p)
    w = E.extrapolation_weights_1d(kn, p, 4, 3, 7)  # nodes 4,5,6 -> slave 7
    x = 3.0
    ref = np.array([(x - 1) * (x - 2) / 2, -(x - 0) * (x - 2), (x - 0) * (x - 1) / 2])
    np.testing.assert_allclose(w, ref, rtol=0, atol=1e-12)


def test_knots_mismatch_raises():
    labels = np.full(27, S.MASTER, dtype=np.int64)
    with pytest.raises(ValueError):
        E.build_extension((3, 3, 3), labels, p=1, knots=_kn((4, 4, 4), 1))


def test_degree_reduced_fallback():
    ncp = (5, 5, 5)
    labels = np.full(125, S.INACTIVE, dtype=np.int64)
    nu = 5
    for c in range(2):
        for b in range(2):
            for a in range(2):
                labels[a + nu * (b + nu * c)] = S.MASTER
    s = 2 + nu * (0 + nu * 0)
    labels[s] = S.EXTENDED
    op = E.build_extension(ncp, labels, p=2, knots=_kn(ncp, 2))
    assert op.degree_used[0] == 1, op.degree_used
    row = op.E.getrow(s)
    assert abs(row.sum() - 1.0) < 1e-12


def test_clean_constant_fallback():
    ncp = (4, 4, 4)
    labels = np.full(64, S.INACTIVE, dtype=np.int64)
    labels[0] = S.MASTER
    s = 3 + 4 * (3 + 4 * 3)
    labels[s] = S.EXTENDED
    op = E.build_extension(ncp, labels, p=2, knots=_kn(ncp, 2))
    assert op.degree_used[0] == 0
    assert op.issues == []
    row = op.E.getrow(s)
    np.testing.assert_allclose(row.toarray().ravel(), [1.0], atol=0)


def test_thin_feature_issue_and_refine_advice():
    ncp = (8, 8, 8)
    labels = np.full(8 ** 3, S.INACTIVE, dtype=np.int64)
    labels[0] = S.MASTER
    s = 7 + 8 * (7 + 8 * 7)
    labels[s] = S.EXTENDED
    op = E.build_extension(ncp, labels, p=2, knots=_kn(ncp, 2))
    assert op.degree_used[0] == -1
    assert len(op.issues) == 1
    iss = op.issues[0]
    assert iss["code"] == E.THIN_FEATURE and iss["severity"] == "warning"
    assert "refine" in iss["message"].lower() and iss["entity"] == f"basis:{s}"
    row = op.E.getrow(s)
    np.testing.assert_allclose(row.toarray().ravel(), [1.0], atol=0)


def test_no_master_raises():
    labels = np.full(27, S.EXTENDED, dtype=np.int64)
    with pytest.raises(ValueError):
        E.build_extension((3, 3, 3), labels, p=1, knots=_kn((3, 3, 3), 1))


def test_outward_tiebreak_deterministic():
    ncp = (5, 4, 4)
    labels = np.full(5 * 4 * 4, S.INACTIVE, dtype=np.int64)
    nu = 5
    for c in (2, 3):
        for b in (2, 3):
            for a in (0, 1):
                labels[a + nu * (b + 4 * c)] = S.MASTER  # block A lo=(0,2,2)
            for a in (3, 4):
                labels[a + nu * (b + 4 * c)] = S.MASTER  # block B lo=(3,2,2)
    s = 2 + nu * (2 + 4 * 2)
    labels[s] = S.EXTENDED
    op_in = E.build_extension(ncp, labels, p=1, knots=_kn(ncp, 1),
                              outward={s: np.array([-1.0, 0.0, 0.0])})
    op_out = E.build_extension(ncp, labels, p=1, knots=_kn(ncp, 1),
                               outward={s: np.array([1.0, 0.0, 0.0])})
    # outward -x (interior +x): block B; outward +x: block A. Both sum to 1.
    row_in = op_in.E.getrow(s).toarray().ravel()
    row_out = op_out.E.getrow(s).toarray().ravel()
    assert abs(row_in.sum() - 1.0) < 1e-12 and abs(row_out.sum() - 1.0) < 1e-12
    cols_in = set(op_in.E.getrow(s).indices.tolist())
    cols_out = set(op_out.E.getrow(s).indices.tolist())
    assert cols_in != cols_out  # tie actually broken by the normal


def test_theta_continuity_nesting():
    prev_slaves: set[int] | None = None
    prev_masters: set[int] | None = None
    n_inactive = None
    for th in (0.05, 0.2, 0.35):
        _ncp, labels = _hole_labels(theta=th)
        slaves = set(np.flatnonzero(labels == S.EXTENDED).tolist())
        masters = set(np.flatnonzero(labels == S.MASTER).tolist())
        nin = int(np.sum(labels == S.INACTIVE))
        if prev_slaves is not None:
            assert prev_slaves <= slaves
            assert masters <= prev_masters
            assert nin == n_inactive  # inside-set is theta-independent
        prev_slaves, prev_masters, n_inactive = slaves, masters, nin


def test_coefficient_monitoring_recorded():
    ncp, labels = _hole_labels()
    op = E.build_extension(ncp, labels, p=2, knots=_kn(ncp, 2))
    assert op.stats["cap"] is None  # TBD upstream: recorded, not enforced
    assert np.isfinite(op.stats["max_abs_coef"])
    assert op.stats["max_abs_coef"] >= 1.0
    assert op.stats["n_master"] + op.stats["n_slave"] <= len(labels)
