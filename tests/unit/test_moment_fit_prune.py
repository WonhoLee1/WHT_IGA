"""fit_moment_rule: one NNLS solve when the pool is large, rule still exact.

Regression (RESPONSE_ANALYSIS_G2_BOTTLENECKS §1, 2026-11-15): the prune loop
kept zero-weight candidates in the pool and tested the pool size instead of
the point count, so a q=6 cell with ~1760 candidates re-solved NNLS once per
active point (343 calls, 5 min). Lawson-Hanson already returns <= (q+1)^3
positive weights, so a large pool must cost exactly one solve.

Oracle: the moments come from an independent tensor Gauss rule, and the
fitted rule is checked on plain monomials against their analytic integrals.
"""

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from wht_iga.geometry import quadrature as Q


def _gauss_rule(blo, h, n):
    x, w = leggauss(n)
    pts = [blo[d] + 0.5 * h[d] * (x + 1.0) for d in range(3)]
    W = [0.5 * h[d] * w for d in range(3)]
    X, Y, Z = np.meshgrid(*pts, indexing="ij")
    WW = np.einsum("i,j,k->ijk", *W)
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]), WW.ravel()


@pytest.mark.parametrize("q", [2, 4])
def test_large_pool_single_nnls_and_exact(monkeypatch, q):
    blo = np.array([0.3, -0.2, 1.0])
    h = np.array([0.5, 0.25, 0.4])
    vol = float(np.prod(h))
    gp, gw = _gauss_rule(blo, h, q + 1)
    moments = (Q._design_matrix(gp, blo, h, q) @ gw).reshape(q + 1, q + 1, q + 1)
    max_pts = (q + 1) ** 3
    # pool = Gauss points of the 2x2x2 sub-boxes (8 * max_pts > 2 * max_pts);
    # their union is an exact positive rule, so a zero-residual fit exists
    cand = np.vstack([_gauss_rule(blo + np.array(o) * h / 2, h / 2, q + 1)[0]
                      for o in np.ndindex(2, 2, 2)])

    calls = []
    real = Q.nnls

    def counting(A, b, **kw):
        calls.append(A.shape)
        return real(A, b, **kw)

    monkeypatch.setattr(Q, "nnls", counting)
    w, res, used = Q.fit_moment_rule(moments, cand, blo, h, vol)
    assert len(calls) == 1, calls
    assert res <= Q._NNLS_TOL * vol
    assert 0 < int((w > 0).sum()) <= max_pts
    # exact for every monomial x^a y^b z^c with a, b, c <= q
    for a in range(q + 1):
        for b in range(q + 1):
            for c in range(q + 1):
                f = used[:, 0] ** a * used[:, 1] ** b * used[:, 2] ** c
                exact = np.prod([((blo[d] + h[d]) ** (e + 1) - blo[d] ** (e + 1)) / (e + 1)
                                 for d, e in enumerate((a, b, c))])
                assert abs(float(f @ w) - exact) <= 1e-9 * max(1.0, abs(exact)), (a, b, c)


def test_tolerance_miss_returns_without_refits(monkeypatch):
    """A first fit that misses the tolerance cannot improve by dropping columns."""
    blo = np.zeros(3)
    h = np.ones(3)
    q = 2
    gp, gw = _gauss_rule(blo, h, q + 1)
    moments = (Q._design_matrix(gp, blo, h, q) @ gw).reshape(q + 1, q + 1, q + 1)
    cand = np.full((5, 3), 0.5)  # 5 identical points: cannot match 27 moments
    calls = []
    real = Q.nnls
    monkeypatch.setattr(Q, "nnls", lambda A, b, **kw: (calls.append(1), real(A, b, **kw))[1])
    _w, res, _used = Q.fit_moment_rule(moments, cand, blo, h, 1.0)
    assert res > Q._NNLS_TOL
    assert len(calls) == 1
