"""Linear solver execution on SciPy CSR matrices and Kratos ModelPart strategies."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.linalg.pardiso import PardisoFactor
from wht_iga.solver.memory import SolverSelectionDecision, select_linear_solver


DEFAULT_AMGCL_PARAMS: dict[str, Any] = {
    "solver_type": "amgcl",
    "krylov_type": "cg",
    "smoother_type": "ilu0",
    "coarsening_type": "aggregation",
    "max_iteration": 5000,
    "tolerance": 1e-10,
    "verbosity": 0,
    "use_block_matrices_if_possible": False,
}


@dataclass(frozen=True)
class LinearSolveResult:
    """Structured outcome of a linear system solve ``K u = f``."""

    solution: np.ndarray
    solver_used: str
    iterations: int
    rel_residual: float
    wall_time_s: float
    converged: bool
    decision: SolverSelectionDecision | None = None


def compute_relative_residual(K: sp.spmatrix, u: np.ndarray, f: np.ndarray) -> float:
    """Compute exact scaled L2 residual ``||K u - f||_2 / max(||f||_2, 1e-30)``."""
    u_vec = np.asarray(u, dtype=np.float64).ravel()
    f_vec = np.asarray(f, dtype=np.float64).ravel()
    if f_vec.size == 0:
        return 0.0
    rhs_norm = float(np.linalg.norm(f_vec))
    res_vec = K @ u_vec - f_vec
    res_norm = float(np.linalg.norm(res_vec))
    return float(res_norm / max(rhs_norm, 1e-30))


def scipy_csr_to_kratos_matrix(KM: Any, K_csr: sp.csr_matrix) -> Any:
    """Convert a SciPy CSR matrix into a ``Kratos.CompressedMatrix``."""
    K_csr = sp.csr_matrix(K_csr, dtype=np.float64)
    K_csr.sort_indices()
    n_rows, n_cols = K_csr.shape
    A_km = KM.CompressedMatrix(n_rows, n_cols)
    indptr = K_csr.indptr
    indices = K_csr.indices
    data = K_csr.data
    for i in range(n_rows):
        s = int(indptr[i])
        e = int(indptr[i + 1])
        for idx in range(s, e):
            A_km[i, int(indices[idx])] = float(data[idx])
    return A_km


def build_kratos_linear_solver(
    KM: Any,
    solver_name: str,
    *,
    num_threads_mkl: int | str | None = None,
    amgcl_params: dict[str, Any] | None = None,
    LA: Any = None,
    tol: float = 1e-10,
    max_iter: int = 5000,
) -> Any:
    """Instantiate a Kratos ``LinearSolver`` object by name.

    Supports both Phase-3 canonical names (``"pardiso"``, ``"amgcl"``,
    ``"cg"``, ``"splu"``) and Phase-1/v2.1 direct names
    (``"pardiso_ldlt"``, ``"pardiso_lu"``, ``"pardiso_llt"``, ``"amgcl_cg"``,
    ``"skyline_lu"``).
    """
    if LA is None:
        from wht_iga.kratos_backend import load_kratos

        LA = load_kratos()["LinearSolversApplication"]

    pardiso_map = {
        "pardiso": "PardisoLDLTSolver",
        "pardiso_ldlt": "PardisoLDLTSolver",
        "pardiso_lu": "PardisoLUSolver",
        "pardiso_llt": "PardisoLLTSolver",
    }
    name = str(solver_name).strip()
    if name in pardiso_map:
        if not LA.HasMKL():
            raise RuntimeError(
                f"linear_solver={name!r} requested, but Kratos LinearSolversApplication.HasMKL() is False"
            )
        if num_threads_mkl == "consistent":
            os.environ["MKL_NUM_THREADS"] = "1"
            os.environ["OMP_NUM_THREADS"] = "1"
        elif isinstance(num_threads_mkl, int):
            os.environ["MKL_NUM_THREADS"] = str(num_threads_mkl)
            os.environ["OMP_NUM_THREADS"] = str(num_threads_mkl)
        cls_name = pardiso_map[name]
        if hasattr(LA, cls_name):
            return getattr(LA, cls_name)()
        return KM.LinearSolverFactory().Create(KM.Parameters(json.dumps({"solver_type": name})))

    if name in ("amgcl", "amgcl_cg"):
        merged = dict(DEFAULT_AMGCL_PARAMS)
        if amgcl_params:
            merged.update(amgcl_params)
        merged["solver_type"] = "amgcl"
        params = KM.Parameters(json.dumps(merged))
        factory = KM.LinearSolverFactory()
        if factory.Has("amgcl"):
            return factory.Create(params)
        return KM.AMGCLSolver(params)

    if name == "cg":
        precond = KM.DiagonalPreconditioner()
        return KM.CGSolver(float(tol), int(max_iter), precond)

    if name == "splu":
        if hasattr(LA, "SparseLUSolver"):
            return LA.SparseLUSolver()
        raise ValueError("Kratos LinearSolversApplication has no SparseLUSolver for strategy 'splu'")

    if name == "skyline_lu":
        return KM.SkylineLUFactorizationSolver()

    raise ValueError(
        f"Unknown linear solver {solver_name!r}; expected one of "
        "'auto', 'pardiso', 'pardiso_ldlt', 'pardiso_lu', 'pardiso_llt', "
        "'amgcl', 'amgcl_cg', 'cg', 'splu', 'skyline_lu'"
    )


def _solve_scipy_cg(
    K_csr: sp.csr_matrix,
    f_vec: np.ndarray,
    tol: float,
    max_iter: int,
) -> tuple[np.ndarray, int]:
    """Solve ``K u = f`` with preconditioned SciPy Conjugate Gradient."""
    n = K_csr.shape[0]
    precond: spla.LinearOperator | None = None
    try:
        ilu = spla.spilu(K_csr.tocsc(), drop_tol=1e-4, fill_factor=10.0)
        precond = spla.LinearOperator((n, n), matvec=ilu.solve, dtype=np.float64)
    except Exception:
        diag = K_csr.diagonal().astype(np.float64, copy=True)
        safe = np.where(np.abs(diag) > 1e-30, 1.0 / diag, 1.0)
        precond = spla.LinearOperator((n, n), matvec=lambda v: safe * v, dtype=np.float64)

    iters = 0

    def _cb(_xk: np.ndarray) -> None:
        nonlocal iters
        iters += 1

    u, _info = spla.cg(
        K_csr,
        f_vec,
        rtol=float(tol),
        atol=0.0,
        maxiter=int(max_iter),
        M=precond,
        callback=_cb,
    )
    return np.asarray(u, dtype=np.float64), int(iters)


def _solve_kratos_amgcl_csr(
    K_csr: sp.csr_matrix,
    f_vec: np.ndarray,
    amgcl_params: dict[str, Any] | None,
    tol: float,
    max_iter: int,
) -> tuple[np.ndarray, int]:
    """Solve CSR system via Kratos ``LinearSolverFactory`` / ``AMGCLSolver``."""
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    LA = mods["LinearSolversApplication"]

    merged = dict(DEFAULT_AMGCL_PARAMS)
    merged["tolerance"] = float(tol)
    merged["max_iteration"] = int(max_iter)
    if amgcl_params:
        merged.update(amgcl_params)
    merged["solver_type"] = "amgcl"

    solver = build_kratos_linear_solver(KM, "amgcl", amgcl_params=merged, LA=LA, tol=tol, max_iter=max_iter)
    A_km = scipy_csr_to_kratos_matrix(KM, K_csr)
    b_km = KM.Vector(f_vec.astype(np.float64).tolist())
    x_km = KM.Vector(K_csr.shape[0], 0.0)
    solver.Solve(A_km, x_km, b_km)
    u = np.array([x_km[i] for i in range(K_csr.shape[0])], dtype=np.float64)
    iters = int(solver.GetIterationsNumber())
    return u, iters


def solve_linear_csr(
    K: sp.spmatrix,
    f: np.ndarray,
    solver_spec: str = "auto",
    *,
    degree: int = 2,
    memory_budget_ratio: float = 0.70,
    avail_bytes: int | None = None,
    dof_threshold: int = 150_000,
    tol: float = 1e-10,
    max_iter: int = 5000,
    amgcl_params: dict[str, Any] | None = None,
    conv_tol: float | None = None,
) -> LinearSolveResult:
    """Solve sparse linear system ``K u = f`` and return verified metrics.

    Supports ``"auto"``, ``"pardiso"``, ``"pardiso_ldlt"``, ``"pardiso_lu"``,
    ``"pardiso_llt"``, ``"amgcl"``, ``"amgcl_cg"``, ``"cg"``, and ``"splu"``.
    Always evaluates the true residual ``||K u - f||_2 / max(||f||_2, 1e-30)``.
    """
    K_csr = sp.csr_matrix(K, dtype=np.float64)
    K_csr.sum_duplicates()
    K_csr.sort_indices()
    f_vec = np.asarray(f, dtype=np.float64).ravel()
    if K_csr.shape[0] != K_csr.shape[1]:
        raise ValueError(f"System matrix must be square, got shape {K_csr.shape}")
    if f_vec.shape[0] != K_csr.shape[0]:
        raise ValueError(f"RHS length {f_vec.shape[0]} != matrix dimension {K_csr.shape[0]}")

    n_dof = int(K_csr.shape[0])
    nnz = int(K_csr.nnz)
    if n_dof == 0:
        return LinearSolveResult(
            solution=np.zeros(0, dtype=np.float64),
            solver_used=str(solver_spec),
            iterations=0,
            rel_residual=0.0,
            wall_time_s=0.0,
            converged=True,
            decision=None,
        )

    decision = select_linear_solver(
        solver_spec=solver_spec,
        n_dof=n_dof,
        nnz=nnz,
        degree=degree,
        memory_budget_ratio=memory_budget_ratio,
        avail_bytes=avail_bytes,
        dof_threshold=dof_threshold,
    )
    chosen = decision.selected_solver
    eff_tol = float(conv_tol if conv_tol is not None else max(tol, 1e-8))

    t0 = time.perf_counter()
    if chosen in ("pardiso", "pardiso_ldlt", "pardiso_lu", "pardiso_llt"):
        kind = "spd" if chosen in ("pardiso_ldlt", "pardiso_llt") else "auto"
        try:
            factor = PardisoFactor(K_csr, kind=kind)
            u = np.asarray(factor.solve(f_vec), dtype=np.float64)
            iters = 1
        except Exception:
            # Explicit direct fallback to Kratos MKL Pardiso on CompressedMatrix.
            from wht_iga.kratos_backend import load_kratos

            mods = load_kratos()
            KM = mods["KratosMultiphysics"]
            LA = mods["LinearSolversApplication"]
            k_solver = build_kratos_linear_solver(KM, chosen, LA=LA)
            A_km = scipy_csr_to_kratos_matrix(KM, K_csr)
            b_km = KM.Vector(f_vec.tolist())
            x_km = KM.Vector(n_dof, 0.0)
            k_solver.Solve(A_km, x_km, b_km)
            u = np.array([x_km[i] for i in range(n_dof)], dtype=np.float64)
            iters = max(int(k_solver.GetIterationsNumber()), 1)
    elif chosen in ("amgcl", "amgcl_cg"):
        u, iters = _solve_kratos_amgcl_csr(K_csr, f_vec, amgcl_params, tol=tol, max_iter=max_iter)
    elif chosen == "cg":
        u, iters = _solve_scipy_cg(K_csr, f_vec, tol=tol, max_iter=max_iter)
    elif chosen == "splu":
        lu = spla.splu(K_csr.tocsc())
        u = np.asarray(lu.solve(f_vec), dtype=np.float64)
        iters = 1
    else:
        raise ValueError(f"Unsupported linear solver {chosen!r}")

    wall_s = time.perf_counter() - t0
    rel_res = compute_relative_residual(K_csr, u, f_vec)
    converged = bool(np.all(np.isfinite(u)) and rel_res <= eff_tol)

    return LinearSolveResult(
        solution=u,
        solver_used=chosen,
        iterations=int(iters),
        rel_residual=float(rel_res),
        wall_time_s=float(wall_s),
        converged=converged,
        decision=decision,
    )
