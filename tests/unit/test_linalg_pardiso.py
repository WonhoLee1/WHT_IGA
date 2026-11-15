"""Tuned pypardiso wrapper vs scipy SuperLU (independent oracle)."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.linalg import PardisoFactor, pardiso_solve


def _laplace3d(n: int) -> sp.csr_matrix:
    e = np.ones(n)
    T = sp.diags([-e[:-1], 2 * e, -e[:-1]], [-1, 0, 1])
    I_ = sp.identity(n)
    return (sp.kron(sp.kron(T, I_), I_) + sp.kron(sp.kron(I_, T), I_)
            + sp.kron(sp.kron(I_, I_), T)).tocsr()


@pytest.mark.parametrize("kind", ["spd", "sym", "general"])
@pytest.mark.parametrize("parallel_ordering,two_level", [(True, False), (False, False), (True, True)])
def test_matches_superlu(kind, parallel_ordering, two_level):
    A = _laplace3d(12)
    b = np.random.default_rng(0).standard_normal((A.shape[0], 2))
    ref = spla.splu(A.tocsc()).solve(b)
    f = PardisoFactor(A, kind=kind, parallel_ordering=parallel_ordering, two_level=two_level)
    x = f.solve(b)
    f.free()
    np.testing.assert_allclose(x, ref, rtol=0, atol=1e-10 * np.abs(ref).max())


def test_reuse_factor_many_rhs():
    A = _laplace3d(10)
    f = PardisoFactor(A, kind="spd")
    for s in range(3):
        b = np.random.default_rng(s).standard_normal(A.shape[0])
        r = A @ f.solve(b) - b
        assert np.linalg.norm(r) <= 1e-10 * np.linalg.norm(b)
    f.free()


def test_symmetric_kinds_reject_nonsymmetric():
    A = _laplace3d(5).tolil()
    A[0, 1] = A[0, 1] + 0.5
    A = A.tocsr()
    with pytest.raises(ValueError):
        PardisoFactor(A, kind="spd")
    x = pardiso_solve(A, np.ones(A.shape[0]), kind="general")
    np.testing.assert_allclose(A @ x, np.ones(A.shape[0]), atol=1e-10)


def test_spd_on_indefinite_raises():
    A = _laplace3d(4) - 3.0 * sp.identity(64, format="csr")  # indefinite
    with pytest.raises(Exception):
        PardisoFactor(A, kind="spd")


def test_bad_kind():
    with pytest.raises(ValueError):
        PardisoFactor(_laplace3d(3), kind="lu")
