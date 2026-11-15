"""Unit tests for Phase-3 RAM probing, solver auto-selection, 6-stage logging, and CSR solvers."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.kratos_backend.solvers import solve_linear_csr
from wht_iga.solver.memory import (
    est_amgcl_bytes,
    est_pardiso_bytes,
    estimate_solver_memory,
    get_available_ram_bytes,
    get_total_ram_bytes,
    probe_system_memory,
    select_linear_solver,
)
from wht_iga.solver.solver_log import SolverLogger


def _spd_3d_poisson(n: int = 8) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Assemble an SPD 7-point 3D Laplacian CSR matrix with known exact solution."""
    I1 = sp.eye(n, format="csr")
    D1 = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n), format="csr")
    K = sp.kron(sp.kron(I1, I1), D1) + sp.kron(sp.kron(I1, D1), I1) + sp.kron(sp.kron(D1, I1), I1)
    K = K.tocsr()
    x = np.linspace(0.1, 1.3, K.shape[0], dtype=np.float64)
    u_true = np.sin(x) + 0.25 * x
    f = np.asarray(K @ u_true, dtype=np.float64)
    return K, f, u_true


def test_host_ram_probing_and_env_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("WHT_SOLVER_AVAIL_RAM_GB", raising=False)
    total = get_total_ram_bytes()
    avail = get_available_ram_bytes()
    info = probe_system_memory()
    assert total > 0
    assert avail > 0
    assert avail <= total
    assert info.total_bytes == total
    assert info.available_bytes == pytest.approx(avail, rel=0.05)
    assert 0.0 <= info.used_percent <= 100.0
    assert info.source

    monkeypatch.setenv("WHT_SOLVER_AVAIL_RAM_GB", "2.5")
    avail_override = get_available_ram_bytes()
    assert avail_override == pytest.approx(2.5 * (1024.0 ** 3), rel=1e-6)
    info2 = probe_system_memory()
    assert info2.available_gb == pytest.approx(2.5, rel=1e-6)
    assert "WHT_SOLVER_AVAIL_RAM_GB" in info2.source


def test_est_pardiso_bytes_scaling() -> None:
    assert est_pardiso_bytes(0, 0, degree=2) == 0
    assert est_pardiso_bytes(100, 0, degree=2) == 0
    assert est_amgcl_bytes(0, 0) == 0

    b_small = est_pardiso_bytes(1_000, 50_000, degree=2)
    b_med = est_pardiso_bytes(10_000, 500_000, degree=2)
    b_deg3 = est_pardiso_bytes(10_000, 500_000, degree=3)
    assert 0 < b_small < b_med < b_deg3

    est = estimate_solver_memory(50_000, 10_000_000, degree=2)
    assert est.pardiso_bytes > est.amgcl_bytes > 0


def test_select_linear_solver_auto_and_explicit() -> None:
    gib = 1024 ** 3
    # Fits RAM budget and below DOF threshold -> pardiso (or splu if pypardiso absent).
    dec_ok = select_linear_solver(
        "auto",
        n_dof=3_000,
        nnz=90_000,
        degree=2,
        avail_bytes=16 * gib,
        pypardiso_available=True,
    )
    assert dec_ok.selected_solver == "pardiso"
    assert dec_ok.switched is False
    assert dec_ok.switched_reason is None

    dec_no_pyp = select_linear_solver(
        "auto",
        n_dof=3_000,
        nnz=90_000,
        degree=2,
        avail_bytes=16 * gib,
        pypardiso_available=False,
    )
    assert dec_no_pyp.selected_solver == "splu"
    assert dec_no_pyp.switched is False

    # Tight RAM budget -> switches to amgcl with explicit reason.
    dec_ram = select_linear_solver(
        "auto",
        n_dof=20_000,
        nnz=4_000_000,
        degree=2,
        memory_budget_ratio=0.70,
        avail_bytes=int(0.5 * gib),
        pypardiso_available=True,
    )
    assert dec_ram.selected_solver == "amgcl"
    assert dec_ram.switched is True
    assert dec_ram.switched_reason is not None
    assert "AUTO-SOLVER" in dec_ram.switched_reason
    assert "avail" in dec_ram.switched_reason

    # DOF threshold exceeded even with massive RAM -> switches to amgcl.
    dec_dof = select_linear_solver(
        "auto",
        n_dof=200_000,
        nnz=2_000_000,
        degree=2,
        avail_bytes=512 * gib,
        dof_threshold=150_000,
        pypardiso_available=True,
    )
    assert dec_dof.selected_solver == "amgcl"
    assert dec_dof.switched is True
    assert "dof_threshold" in (dec_dof.switched_reason or "")

    # Explicit user choices are never silently changed.
    for explicit in ("pardiso", "amgcl", "cg", "splu", "pardiso_ldlt"):
        dec_exp = select_linear_solver(
            explicit,
            n_dof=200_000,
            nnz=5_000_000,
            avail_bytes=int(0.1 * gib),
        )
        assert dec_exp.selected_solver == explicit
        assert dec_exp.switched is False
        assert dec_exp.switched_reason is None


def test_solver_logger_6_stages_utf8_and_switches(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("WHT_SOLVER_LOG", raising=False)
    log_path = tmp_path / "run_solver.log"
    stream = io.StringIO()

    dec = select_linear_solver(
        "auto",
        n_dof=50_000,
        nnz=10_000_000,
        degree=2,
        avail_bytes=int(0.5 * (1024 ** 3)),
    )
    with SolverLogger(
        enable_solver_log=True,
        solver_log_file=log_path,
        solver_log_stdout=True,
        solver_log_level="detailed",
        stdout_stream=stream,
    ) as logger:
        logger.log_stage_1_environment(
            job_name="unicode_테스트_job",
            cpu_count=8,
            total_ram_gb=32.0,
            avail_ram_gb=16.0,
            mkl_threads=4,
        )
        logger.log_stage_2_model(
            n_cells=(4, 4, 4),
            degree=2,
            continuity=1,
            n_ctrl_pts=216,
            n_total_dofs=648,
            n_active_dofs=600,
            n_in_cells=50,
            n_cut_cells=10,
            n_out_cells=4,
            out_mode="drop",
            n_quad_pts=1620,
            nnz=120000,
        )
        logger.log_stage_3_constraints(
            n_dirichlet_bcs=1,
            n_fixed_faces=1,
            n_inactive_fixed=48,
            n_ties=4,
            n_slave_dofs=12,
            n_master_dofs=48,
            beta=1e6,
            c_pen=100.0,
            theta_thb=0.2,
        )
        logger.log_stage_4_solver_selection(decision=dec)
        logger.log_stage_5_header()
        logger.log_stage_5_row(
            step=1, inc=1, att=1, time=1.0, dt=1.0, lin_it=12, rel_res=1.5e-11, energy_norm=0.42, wall_s=0.125
        )
        logger.log_stage_6_summary(
            converged=True,
            max_disp=1.23e-3,
            total_force=[100.0, 0.0, 0.0],
            reaction_residual=1e-12,
            setup_s=0.1,
            solve_s=0.2,
            post_s=0.05,
        )
        text_mem = logger.get_text()

    assert log_path.exists()
    text_disk = log_path.read_text(encoding="utf-8")
    assert text_disk == text_mem
    assert stream.getvalue().strip() == text_mem.strip()
    assert "HWT-IGA SOLVER — Immersed Isogeometric Analysis Engine" in text_disk
    for stage in (
        "[1/6] JOB & SYSTEM ENVIRONMENT",
        "[2/6] MODEL & IMMERSED GRID SUMMARY",
        "[3/6] CONSTRAINT & INTERACTION SUMMARY",
        "[4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION",
        "[5/6] STEP & ITERATION PROGRESS TABLE",
        "[6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY",
    ):
        assert stage in text_disk
    assert "Step | Inc | Att |      Time |        dt | LinIt |      RelRes |  EnergyNorm | Wall(s)" in text_disk
    assert "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in text_disk
    assert "unicode_테스트_job" in text_disk

    # Disabled logger emits nothing and does not create files.
    disabled_path = tmp_path / "should_not_exist.log"
    disabled_stream = io.StringIO()
    with SolverLogger(
        enable_solver_log=False,
        solver_log_file=disabled_path,
        solver_log_stdout=True,
        stdout_stream=disabled_stream,
    ) as dlog:
        dlog.log_stage_1_environment(job_name="off")
        dlog.log_stage_6_summary(converged=False)
        assert dlog.get_text() == ""
    assert not disabled_path.exists()
    assert disabled_stream.getvalue() == ""

    # Env override WHT_SOLVER_LOG=0 also suppresses output.
    monkeypatch.setenv("WHT_SOLVER_LOG", "0")
    env_off_path = tmp_path / "env_off.log"
    with SolverLogger(enable_solver_log=True, solver_log_file=env_off_path) as elog:
        elog.log_stage_1_environment(job_name="env_off")
        assert elog.enabled is False
        assert elog.get_text() == ""
    assert not env_off_path.exists()


def test_solve_linear_csr_solvers_real_residuals() -> None:
    K, f, u_true = _spd_3d_poisson(n=7)
    for spec_name in ("pardiso", "splu", "cg"):
        res = solve_linear_csr(K, f, solver_spec=spec_name, tol=1e-11, max_iter=2000)
        assert res.converged is True
        assert res.solver_used == spec_name
        assert res.rel_residual < 1e-8
        assert res.wall_time_s >= 0.0
        rel_err = np.linalg.norm(res.solution - u_true) / np.linalg.norm(u_true)
        assert rel_err < 1e-8
        if spec_name == "cg":
            assert res.iterations >= 1
        else:
            assert res.iterations == 1

    # Auto selection on small system picks pardiso with true solution.
    res_auto = solve_linear_csr(K, f, solver_spec="auto", avail_bytes=16 * (1024 ** 3))
    assert res_auto.converged is True
    assert res_auto.solver_used == "pardiso"
    assert res_auto.decision is not None and res_auto.decision.switched is False
    assert np.linalg.norm(res_auto.solution - u_true) / np.linalg.norm(u_true) < 1e-8
