"""Kratos integration tests for AMGCL solver, RAM auto-switching, 6-stage logs, drop mode, and zero-weight probes."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.kratos_backend.solvers import solve_linear_csr
from wht_iga.post.probes import PointProbe, evaluate_zero_weight_kratos_probes
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.pipeline import run_linear
from wht_iga.solver.validate import JobSpec, TractionLoad


def _spd_3d_laplacian(n: int = 12) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Assemble an SPD 1728-DOF 3D Laplacian CSR system with exact solution."""
    I1 = sp.eye(n, format="csr")
    D1 = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n), format="csr")
    K = sp.kron(sp.kron(I1, I1), D1) + sp.kron(sp.kron(I1, D1), I1) + sp.kron(sp.kron(D1, I1), I1)
    K = K.tocsr()
    rng_grid = np.linspace(0.0, 1.0, K.shape[0], dtype=np.float64)
    u_true = np.sin(2.0 * np.pi * rng_grid) + 0.5 * rng_grid
    f = np.asarray(K @ u_true, dtype=np.float64)
    return K, f, u_true


def test_kratos_amgcl_csr_and_auto_switch() -> None:
    K, f, u_true = _spd_3d_laplacian(n=12)

    res_amgcl = solve_linear_csr(
        K,
        f,
        solver_spec="amgcl",
        tol=1e-10,
        max_iter=2000,
        amgcl_params={"krylov_type": "cg", "smoother_type": "ilu0", "tolerance": 1e-11},
    )
    assert res_amgcl.converged is True
    assert res_amgcl.solver_used == "amgcl"
    assert res_amgcl.rel_residual < 1e-8
    assert res_amgcl.wall_time_s > 0.0
    assert np.linalg.norm(res_amgcl.solution - u_true) / np.linalg.norm(u_true) < 1e-7

    # Force auto switch via tiny available RAM -> switches to Kratos amgcl.
    res_auto = solve_linear_csr(
        K,
        f,
        solver_spec="auto",
        avail_bytes=1024,
        memory_budget_ratio=0.70,
        tol=1e-10,
        max_iter=2000,
    )
    assert res_auto.decision is not None
    assert res_auto.decision.switched is True
    assert res_auto.solver_used == "amgcl"
    assert res_auto.converged is True
    assert res_auto.rel_residual < 1e-8
    assert np.linalg.norm(res_auto.solution - u_true) / np.linalg.norm(u_true) < 1e-7


def test_pipeline_amgcl_and_6stage_solver_log_file(tmp_path: Path) -> None:
    log_path = tmp_path / "amgcl_run.log"
    out_dir = tmp_path / "out_amgcl"

    base_kwargs = dict(
        name="cube_amgcl",
        lo=(0.0, 0.0, 0.0),
        hi=(2.0, 1.0, 1.0),
        n_cells=(3, 2, 2),
        degree=2,
        E=1000.0,
        nu=0.3,
        fixed_faces=[(0, 0)],
        tractions=[TractionLoad((0, 1), (10.0, 0.0, 0.0))],
        probes=[{"name": "tip", "xyz": [2.0, 0.5, 0.5]}],
    )
    res_pardiso = run_linear(JobSpec(**base_kwargs, linear_solver="pardiso_ldlt", enable_solver_log=False))
    assert res_pardiso.converged is True
    assert res_pardiso.solver_log == ""

    res_amgcl = run_linear(
        JobSpec(
            **base_kwargs,
            linear_solver="amgcl_cg",
            enable_solver_log=True,
            solver_log_file=str(log_path),
            solver_log_level="detailed",
            amgcl_params={"krylov_type": "cg", "smoother_type": "ilu0", "tolerance": 1e-11, "max_iteration": 5000},
            output_dir=str(out_dir),
        )
    )
    assert res_amgcl.converged is True
    assert res_amgcl.timings["linear_solver"] == "amgcl_cg"
    assert np.allclose(res_amgcl.U, res_pardiso.U, rtol=1e-5, atol=1e-6)
    assert np.allclose(res_amgcl.probe_values["tip"], res_pardiso.probe_values["tip"], rtol=1e-5, atol=1e-6)

    assert log_path.exists()
    disk_log = log_path.read_text(encoding="utf-8")
    assert disk_log == res_amgcl.solver_log
    assert "HWT-IGA SOLVER — Immersed Isogeometric Analysis Engine" in disk_log
    for stage in (
        "[1/6] JOB & SYSTEM ENVIRONMENT",
        "[2/6] MODEL & IMMERSED GRID SUMMARY",
        "[3/6] CONSTRAINT & INTERACTION SUMMARY",
        "[4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION",
        "[5/6] STEP & ITERATION PROGRESS TABLE",
        "[6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY",
    ):
        assert stage in disk_log
    assert "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in disk_log


def test_pipeline_auto_switch_and_log_disable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    log_auto = tmp_path / "auto_switch.log"
    monkeypatch.setenv("WHT_SOLVER_AVAIL_RAM_GB", "0.000001")
    spec_auto = JobSpec(
        name="auto_switch_job",
        lo=(0.0, 0.0, 0.0),
        hi=(1.0, 1.0, 1.0),
        n_cells=(2, 2, 2),
        degree=2,
        E=1000.0,
        nu=0.3,
        fixed_faces=[(0, 0)],
        tractions=[TractionLoad((0, 1), (1.0, 0.0, 0.0))],
        linear_solver="auto",
        enable_solver_log=True,
        solver_log_file=str(log_auto),
    )
    res_auto = run_linear(spec_auto)
    assert res_auto.converged is True
    assert res_auto.solver_decision is not None
    assert res_auto.solver_decision["switched"] is True
    assert res_auto.timings["linear_solver"] == "amgcl_cg"
    assert "[HWT-IGA SOLVER][AUTO-SOLVER]" in res_auto.solver_log
    assert log_auto.exists()

    # Env WHT_SOLVER_LOG=0 disables logging cleanly without creating file.
    monkeypatch.delenv("WHT_SOLVER_AVAIL_RAM_GB", raising=False)
    monkeypatch.setenv("WHT_SOLVER_LOG", "0")
    log_off = tmp_path / "env_disabled.log"
    res_off = run_linear(
        JobSpec(
            name="log_off_job",
            lo=(0.0, 0.0, 0.0),
            hi=(1.0, 1.0, 1.0),
            n_cells=(2, 2, 2),
            degree=2,
            E=1000.0,
            nu=0.3,
            fixed_faces=[(0, 0)],
            tractions=[TractionLoad((0, 1), (1.0, 0.0, 0.0))],
            enable_solver_log=True,
            solver_log_file=str(log_off),
        )
    )
    assert res_off.converged is True
    assert res_off.solver_log == ""
    assert res_off.solver_log_file is None
    assert not log_off.exists()


def test_pipeline_out_mode_drop_fewer_active_dofs() -> None:
    out_cells = [(i, j, 3) for i in range(3) for j in range(3)]
    common = dict(
        lo=(0.0, 0.0, 0.0),
        hi=(3.0, 1.0, 2.0),
        n_cells=(3, 3, 4),
        degree=2,
        E=1000.0,
        nu=0.25,
        density=1.0,
        gravity=(100.0, 0.0, 0.0),
        fixed_faces=[(0, 0)],
        out_cells=out_cells,
        probes=[{"name": "mid", "xyz": [2.5, 0.5, 0.5]}],
        enable_solver_log=False,
    )
    res_ghost = run_linear(JobSpec(name="m_ghost", out_mode="ghost", **common))
    res_drop = run_linear(JobSpec(name="m_drop", out_mode="drop", **common))

    assert res_ghost.converged is True
    assert res_drop.converged is True
    assert res_drop.active_dofs < res_ghost.active_dofs
    assert res_drop.active_dofs < res_drop.total_dofs
    assert res_drop.probe_values["mid"][0] > 0.0
    # OUT cells in drop mode produce empty stress arrays (0 integration points).
    nx, ny, _nz = 3, 3, 4
    top_slice_start = 3 * nx * ny
    assert all(arr.shape == (0, 6) for arr in res_drop.cauchy_stress[top_slice_start:])


def test_kratos_zero_weight_probe_invariance() -> None:
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SA = mods["StructuralMechanicsApplication"]

    E, nu = 2500.0, 0.25
    eps_x = 0.004
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)

    model = KM.Model()
    mp = model.CreateModelPart("zw_probe_check")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    volume, _keep = make_nurbs_volume(KM, mp, box)

    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, SA.LinearElastic3DLaw())

    for nd in mp.Nodes:
        nd.SetSolutionStepValue(KM.DISPLACEMENT, 0, [eps_x * nd.X, 0.0, 0.0])

    probe = PointProbe("center", (1.0, 0.5, 0.5), variables=("U", "E", "S", "MISES"))
    kratos_ctx = {
        "KM": KM,
        "model_part": mp,
        "volume": volume,
        "props": props,
        "box": box,
        "element_id": 950001,
    }
    sample = probe.evaluate(box, kratos_ctx=kratos_ctx)
    assert sample.mode == "zero_weight_kratos"

    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    expected_s = np.array([(lam + 2.0 * mu) * eps_x, lam * eps_x, lam * eps_x, 0.0, 0.0, 0.0])
    assert np.allclose(sample["U"][0], [eps_x * 1.0, 0.0, 0.0], atol=1e-12)
    assert np.allclose(sample["E"][0], [eps_x, 0.0, 0.0, 0.0, 0.0, 0.0], atol=1e-12)
    assert np.allclose(sample["S"][0], expected_s, atol=1e-9)

    # Also call explicit helper with verify_zero_contribution=True on MultiPoint positions.
    zw_vals = evaluate_zero_weight_kratos_probes(
        KM,
        mp,
        volume,
        props,
        box,
        np.array([[0.3, 0.2, 0.8], [1.7, 0.9, 0.1]]),
        element_id=950002,
        verify_zero_contribution=True,
    )
    assert zw_vals["S"].shape == (2, 6)
    assert np.allclose(zw_vals["S"][0], expected_s, atol=1e-9)
    assert np.allclose(zw_vals["S"][1], expected_s, atol=1e-9)
