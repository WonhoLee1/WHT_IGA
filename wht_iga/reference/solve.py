"""Linear solve with prescribed DOFs (strong elimination)."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def solve_with_dirichlet(K: sp.csr_matrix, f: np.ndarray, fixed: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Solve K u = f with u[fixed] = values. Uses pypardiso when available."""
    n = K.shape[0]
    fixed = np.asarray(fixed, dtype=np.int64)
    free = np.setdiff1d(np.arange(n), fixed)
    u = np.zeros(n)
    u[fixed] = values
    rhs = f[free] - K[free][:, fixed] @ u[fixed]
    Kff = K[free][:, free].tocsr()
    try:
        import pypardiso

        u[free] = pypardiso.spsolve(Kff, rhs)
    except ImportError:
        u[free] = spla.spsolve(Kff, rhs)
    return u
