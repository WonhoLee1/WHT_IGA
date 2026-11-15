"""T-P2.3: T = C_THB · E composition (v2.1 §5.2.3 rules 4-6).

Phase 2 has no THB hierarchy, so C_THB = I and T = E; the composition path
(``compose_thb``) exists so Phase 3 only plugs in a real C_THB. All operators
are ``scipy.sparse`` CSR. Dense inputs are refused (Phase-2 sparse-only rule,
phase2 §7): a dense T on production grids is a memory defect, not a fallback.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _require_sparse(A: sp.spmatrix, name: str) -> sp.csr_matrix:
    if not sp.issparse(A):
        raise TypeError(f"{name} must be scipy.sparse (sparse-only rule); "
                        f"got {type(A).__name__}")
    return A.tocsr()


def extension_matrix(n_full: int, master_ids: np.ndarray,
                     slave_ids: np.ndarray,
                     slave_rows: sp.spmatrix) -> sp.csr_matrix:
    """Assemble full E (n_full × n_master): unit rows on masters.

    :param slave_rows: CSR block (n_slave × n_master) in ``slave_ids`` order.
        Inactive rows stay zero (no node is created for them, v2 §5.2.3/6).
    """
    slave_rows = _require_sparse(slave_rows, "slave_rows").tocoo()
    master_ids = np.asarray(master_ids, dtype=np.int64)
    slave_ids = np.asarray(slave_ids, dtype=np.int64)
    if slave_rows.shape[0] != len(slave_ids):
        raise ValueError(
            f"slave_rows {slave_rows.shape} vs {len(slave_ids)} slave ids")
    n_m = slave_rows.shape[1]
    rows = np.concatenate([master_ids, slave_ids[slave_rows.row]])
    cols = np.concatenate([np.arange(len(master_ids), dtype=np.int64),
                           slave_rows.col])
    data = np.concatenate([np.ones(len(master_ids)), slave_rows.data])
    return sp.coo_matrix((data, (rows, cols)),
                         shape=(n_full, n_m)).tocsr()


def compose_thb(E: sp.spmatrix, C_thb: sp.spmatrix | None = None) -> sp.csr_matrix:
    """Return T = C_THB · E. ``C_thb=None`` selects the Phase-2 path C_THB = I.

    A non-None C_THB is validated (square, matches E rows) and multiplied in
    sparse arithmetic only.
    """
    E = _require_sparse(E, "E")
    if C_thb is None:
        return E.copy()  # Phase 2: no THB → C_THB = I (path only, v2.1 §7.5)
    C = _require_sparse(C_thb, "C_thb")
    if C.shape[0] != C.shape[1] or C.shape[1] != E.shape[0]:
        raise ValueError(f"C_thb shape {C.shape} incompatible with E {E.shape}")
    return (C @ E).tocsr()


def expand_dof_operator(T_basis: sp.spmatrix, n_comp: int = 3) -> sp.csr_matrix:
    """Block-diagonal DOF expansion: T_dof[g*n_comp+d, c*n_comp+d] = T[g, c].

    The extension acts identically on every component (v2 §5.2.3: scalar
    relations per component, shared T).
    """
    T_basis = _require_sparse(T_basis, "T_basis").tocoo()
    comp = np.arange(n_comp, dtype=np.int64)
    rows = (T_basis.row[:, None] * n_comp + comp).ravel()
    cols = (T_basis.col[:, None] * n_comp + comp).ravel()
    data = np.repeat(T_basis.data, n_comp)
    n_full, n_m = T_basis.shape
    return sp.coo_matrix((data, (rows, cols)),
                         shape=(n_full * n_comp, n_m * n_comp)).tocsr()


def apply_TKT(K: sp.spmatrix, T: sp.spmatrix) -> sp.csr_matrix:
    """Reduced stiffness Tᵀ K T (R1 route, v2 §5.2.3)."""
    K = _require_sparse(K, "K")
    T = _require_sparse(T, "T")
    if K.shape[0] != K.shape[1] or K.shape[0] != T.shape[0]:
        raise ValueError(f"shape mismatch: K {K.shape}, T {T.shape}")
    return ((T.T @ K) @ T).tocsr()


def apply_Tt_vec(f: np.ndarray, T: sp.spmatrix) -> np.ndarray:
    """Reduced load Tᵀ f."""
    T = _require_sparse(T, "T")
    f = np.asarray(f, dtype=float)
    if f.shape != (T.shape[0],):
        raise ValueError(f"f shape {f.shape} != ({T.shape[0]},)")
    return np.asarray(T.T @ f).ravel()


def solve_reduced(K_red: sp.spmatrix, f_red: np.ndarray) -> np.ndarray:
    """Direct solve of the reduced system (test/oracle path)."""
    K_red = _require_sparse(K_red, "K_red")
    try:
        from wht_iga.linalg import pardiso_solve

        return pardiso_solve(K_red, np.asarray(f_red, dtype=float), kind="auto")
    except ImportError:
        return spla.spsolve(K_red, np.asarray(f_red, dtype=float))
