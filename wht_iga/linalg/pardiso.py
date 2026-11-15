"""Tuned MKL Pardiso factorization through pypardiso.

Why not ``pypardiso.spsolve``: it shares one module-global solver, always runs
the real-nonsymmetric path (mtype 11: scaling + weighted matching, LU of the
full matrix) with all-default iparm, and copies/compares the whole matrix on
every call. Our stiffness / reduced TᵀKT matrices are symmetric (Kratos zeroes
fixed rows AND columns; penalty and TᵀKT are symmetric), so the symmetric
paths need half the factor memory and roughly half the flops.

Settings (MKL PARDISO iparm, 1-based like the MKL docs and pypardiso):

- iparm(1)  = 1   use the values below (not the all-default table)
- iparm(2)  = 3   OpenMP-parallel nested-dissection ordering (2 = serial METIS)
- iparm(8)  = 0   automatic iterative refinement (2 steps on perturbed pivots)
- iparm(10) = 8 (symmetric) / 13 (general)  pivot perturbation 1e-8 / 1e-13
- iparm(11) = iparm(13) = 0 (spd/sym) / 1 (general)  scaling + matching
- iparm(21) = 1   Bunch-Kaufman pivoting (symmetric indefinite only)
- iparm(24) = two_level (0 classic, 1 two-level parallel factorization)
- iparm(35) = 0   one-based indices (pypardiso adds 1 to indptr/indices)

``kind``: "spd" (mtype 2, Cholesky), "sym" (mtype -2, LDLᵀ), "general"
(mtype 11, LU), "auto" (symmetric → "sym", else "general"; the choice is
stored in ``.kind``). Symmetric kinds receive the upper triangle only. A failed
Cholesky raises (PyPardisoError) instead of silently switching kind.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

_MTYPE = {"spd": 2, "sym": -2, "general": 11}


class PardisoFactor:
    """One factorization, many solves. Not thread-safe (one MKL handle)."""

    def __init__(self, A: sp.spmatrix, kind: str = "spd", parallel_ordering: bool = True,
                 two_level: bool = False, symmetry_rtol: float = 1e-12) -> None:
        import pypardiso

        if kind not in (*_MTYPE, "auto"):
            raise ValueError(f"kind must be one of {(*_MTYPE, 'auto')}, got {kind!r}")
        if not sp.issparse(A):
            raise TypeError("PardisoFactor needs a scipy sparse matrix")
        A = sp.csr_matrix(A, dtype=np.float64)
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"square matrix required, got {A.shape}")
        self.n = A.shape[0]
        asym = abs(A - A.T).max() if A.nnz else 0.0
        scale = abs(A).max() if A.nnz else 1.0
        if kind == "auto":
            kind = "sym" if asym <= symmetry_rtol * scale else "general"
        self.kind = kind
        if kind in ("spd", "sym"):
            if asym > symmetry_rtol * scale:
                raise ValueError(f"kind={kind!r} needs a symmetric matrix: "
                                 f"max|A-Aᵀ|={asym:.3e} > {symmetry_rtol:g}·max|A|")
            A = sp.triu(A, format="csr")
        A.sort_indices()
        self._A = A
        ps = pypardiso.PyPardisoSolver(mtype=_MTYPE[kind])
        sym = kind != "general"
        ps.set_iparm(1, 1)
        ps.set_iparm(2, 3 if parallel_ordering else 2)
        ps.set_iparm(8, 0)
        ps.set_iparm(10, 8 if sym else 13)
        ps.set_iparm(11, 0 if sym else 1)
        ps.set_iparm(13, 0 if sym else 1)
        if kind == "sym":
            ps.set_iparm(21, 1)
        ps.set_iparm(24, 1 if two_level else 0)
        ps.set_iparm(35, 0)
        self._ps = ps
        ps.factorize(A)
        self.n_perturbed_pivots = int(ps.get_iparm(14))
        self.factor_nnz = int(ps.get_iparm(18))

    def solve(self, b: np.ndarray) -> np.ndarray:
        b = np.asarray(b, dtype=np.float64)
        if b.shape[0] != self.n:
            raise ValueError(f"rhs has {b.shape[0]} rows, matrix has {self.n}")
        x = self._ps.solve(self._A, b)
        return x.reshape(b.shape)

    def free(self) -> None:
        self._ps.free_memory(everything=True)


def pardiso_solve(A: sp.spmatrix, b: np.ndarray, kind: str = "spd") -> np.ndarray:
    """One-shot tuned solve (factor, solve, free)."""
    f = PardisoFactor(A, kind=kind)
    try:
        return f.solve(b)
    finally:
        f.free()
