"""T-P2.5: conditioning tool + V-COND analysis helpers.

Matrix source (v2 §5.6, plan appendix A-4/A-5): Kratos ``GetSystemMatrix``
after a strategy run, converted with ``scipy_conversion_tools`` and reordered
with the reference-DOF permutation; the reduced operator is Tᵀ K T applied
here (``reduced_matrix``). This module never imports Kratos: the caller (test
or driver) passes matrices in.

V-COND (v2.1 §9 G2): hole-center sweep of the V02 Roark plate mother geometry
(documented choice — §10.4 TBD resolved as V02 plate; user confirmation
pending), eta_min ∈ [1e-8, 1] over 4 h-levels. Acceptance: extension ON →
kappa variation ≤ 10× (flat); log-log slope −2 ± 0.3 in h; unstabilized
comparison graph. The sweep driver consumes an injected matrix builder so
unit tests run Kratos-free against the reference assembler while Kratos
spikes inject Kratos-extracted matrices.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _require_sparse(A: sp.spmatrix, name: str) -> sp.csr_matrix:
    if not sp.issparse(A):
        raise TypeError(f"{name} must be scipy.sparse, got {type(A).__name__}")
    return A.tocsr()


def reduced_matrix(K: sp.spmatrix, T: sp.spmatrix | None) -> sp.csr_matrix:
    """Tᵀ K T with T = None meaning identity (unstabilized path)."""
    K = _require_sparse(K, "K")
    if T is None:
        return K
    T = _require_sparse(T, "T")
    return ((T.T @ K) @ T).tocsr()


def cond_dense(K: np.ndarray) -> float:
    """Reference condition number (2-norm, exact SVD) — test oracle."""
    s = np.linalg.svd(np.asarray(K, dtype=float), compute_uv=False)
    if s[-1] <= 0.0:
        raise ValueError("matrix is singular")
    return float(s[0] / s[-1])


def cond_eigsh(K: sp.spmatrix) -> tuple[float, float, float]:
    """kappa_2 estimate of SPD K via extreme eigenpairs (largest/smallest).

    Returns (kappa, lambda_min, lambda_max). Raises on non-positive spectra.
    """
    K = _require_sparse(K, "K")
    lmax = float(spla.eigsh(K, k=1, which="LM", tol=1e-10,
                            return_eigenvectors=False)[0])
    # Smallest via shift-invert (sigma = 0 → OPinv = K⁻¹). The inverse is one
    # tuned MKL Pardiso LDLᵀ factor (wht_iga.linalg) instead of SuperLU when
    # pypardiso is available (V-COND matrices reach ~1e5 DOF).
    opinv = None
    try:
        from wht_iga.linalg import PardisoFactor

        fac = PardisoFactor(K, kind="auto")
        opinv = spla.LinearOperator(K.shape, matvec=fac.solve, dtype=float)
    except ImportError:
        fac = None
    try:
        lmin = float(spla.eigsh(K, k=1, sigma=0.0, which="LM", tol=1e-10, OPinv=opinv,
                                return_eigenvectors=False)[0])
    finally:
        if fac is not None:
            fac.free()
    if lmin <= 0.0:
        raise ValueError(f"non-positive smallest eigenvalue {lmin:.3e}")
    return lmax / lmin, lmin, lmax


def h_slope(hs: np.ndarray, kappas: np.ndarray) -> float:
    """Log-log slope of kappa vs h (V-COND acceptance: −2 ± 0.3)."""
    hs = np.asarray(hs, dtype=float)
    ks = np.asarray(kappas, dtype=float)
    if np.any(hs <= 0.0) or np.any(ks <= 0.0):
        raise ValueError("h and kappa must be strictly positive")
    return float(np.polyfit(np.log10(hs), np.log10(ks), 1)[0])


def flatness(kappas: np.ndarray) -> float:
    """max/min ratio of a kappa series (V-COND acceptance: ≤ 10×)."""
    ks = np.asarray(kappas, dtype=float)
    return float(ks.max() / ks.min())
