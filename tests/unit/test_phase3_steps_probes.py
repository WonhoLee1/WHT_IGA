"""Unit tests for Phase-3 multi-step orchestration, amplitudes, THB retry, and probes."""

from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as la

from wht_iga.grid.hierarchy import THBHierarchy
from wht_iga.post import (
    LineProbe,
    MultiPointProbe,
    PointProbe,
    export_probes_json,
    export_probes_sqlite,
)
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.steps import (
    ABAQUS_TO_KRATOS,
    VAR_MAP,
    Amplitude,
    CancelToken,
    CancelledError,
    ExplicitDynamicStep,
    ImplicitDynamicStep,
    ModalStep,
    MultiStepRunner,
    ProgressEvent,
    StaticStep,
    ThermalStep,
    auto_refine_thin_features,
    detect_thin_features,
    resolve_var_name,
)
from wht_iga.solver.validate import JobSpec, validate


def test_amplitudes_all_four_kinds() -> None:
    amp_tab = Amplitude.tabular([(0.0, 0.0), (0.5, 2.0), (1.0, 1.0)])
    assert amp_tab.evaluate(-0.2) == pytest.approx(0.0)
    assert amp_tab.evaluate(0.25) == pytest.approx(1.0)
    assert amp_tab.evaluate(0.75) == pytest.approx(1.5)
    assert amp_tab.evaluate(1.5) == pytest.approx(1.0)

    amp_ramp = Amplitude.ramp(t0=1.0, t1=3.0, a0=2.0, a1=10.0)
    assert amp_ramp.evaluate(0.0) == pytest.approx(2.0)
    assert amp_ramp.evaluate(2.0) == pytest.approx(6.0)
    assert amp_ramp.evaluate(4.0) == pytest.approx(10.0)

    amp_sm = Amplitude.smooth_step(t0=0.0, t1=2.0, a0=0.0, a1=4.0)
    assert amp_sm.evaluate(0.0) == pytest.approx(0.0)
    assert amp_sm.evaluate(1.0) == pytest.approx(2.0)
    assert amp_sm.evaluate(2.0) == pytest.approx(4.0)
    # Zero slope at endpoints: central/one-sided finite difference near t=0 and t=2.
    eps = 1e-4
    assert (amp_sm.evaluate(eps) - amp_sm.evaluate(0.0)) / eps < 1e-6
    assert (amp_sm.evaluate(2.0) - amp_sm.evaluate(2.0 - eps)) / eps < 1e-6

    amp_per = Amplitude.periodic(omega=2.0 * math.pi, a0=1.5, terms=[(2.0, 3.0)], t0=0.0)
    assert amp_per.evaluate(0.0) == pytest.approx(3.5)
    assert amp_per.evaluate(0.25) == pytest.approx(4.5)
    assert amp_per.evaluate(0.5) == pytest.approx(-0.5)


def test_var_map_bidirectional_and_resolve() -> None:
    expected = {
        "U": "DISPLACEMENT",
        "S": "CAUCHY_STRESS_TENSOR",
        "E": "GREEN_LAGRANGE_STRAIN_TENSOR",
        "NT11": "TEMPERATURE",
        "RF": "REACTION",
        "CF": "CONTACT_FORCE",
    }
    assert ABAQUS_TO_KRATOS == expected
    for abq, kra in expected.items():
        assert VAR_MAP[abq] == kra
        assert VAR_MAP[kra] == abq
        assert VAR_MAP[VAR_MAP[abq]] == abq
        assert resolve_var_name(abq.lower(), target="kratos") == kra
        assert resolve_var_name(kra.lower(), target="abaqus") == abq
        assert resolve_var_name(abq, target="opposite") == kra
        assert resolve_var_name(kra, target="opposite") == abq

    with pytest.raises(ValueError):
        resolve_var_name("NONEXISTENT_VAR")


def test_multistep_runner_chaining_physics() -> None:
    # 2-DOF structural oscillator & 2-DOF thermal system.
    K = np.array([[200.0, -100.0], [-100.0, 100.0]], dtype=np.float64)
    M = np.array([[2.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    C = 0.02 * M + 0.0005 * K
    K_th = np.array([[4.0, -1.0], [-1.0, 2.0]], dtype=np.float64)
    C_th = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float64)

    f0 = np.array([10.0, 20.0], dtype=np.float64)
    q0 = np.array([5.0, 2.0], dtype=np.float64)
    runner = MultiStepRunner(
        n_mech_dofs=2,
        n_thermal_dofs=2,
        K=K,
        M=M,
        C=C,
        K_th=K_th,
        C_th=C_th,
        f_ext=lambda step, t_abs, amp, st: f0 * amp if isinstance(step, StaticStep) else f0,
        q_ext=lambda step, t_abs, amp, st: q0,
    )

    events: list[ProgressEvent] = []
    final = runner.run(
        [
            StaticStep("preload", duration=1.0, n_increments=2, amplitude=Amplitude.ramp(0.0, 1.0, 0.0, 1.0)),
            ModalStep("modes", n_modes=2),
            ImplicitDynamicStep("imp_dyn", duration=0.05, n_increments=5),
            ExplicitDynamicStep("exp_dyn", duration=0.01, n_increments=5),
            ThermalStep("steady_th", duration=1.0, n_increments=1, transient=False),
            ThermalStep("trans_th", duration=0.5, n_increments=5, transient=True),
        ],
        on_progress=events.append,
    )

    # Static result K * u = f0 -> exact solution was passed into dynamic steps.
    u_static_ref = np.linalg.solve(K, f0)
    assert final.history[1]["max_u"] == pytest.approx(float(np.max(np.abs(u_static_ref))), rel=1e-10)

    # Modal generalized eigenvalues against scipy.linalg.eigh reference.
    evals_ref, _ = la.eigh(K, M)
    assert final.eigenvalues is not None
    assert np.allclose(final.eigenvalues, evals_ref, rtol=1e-10)

    # Under constant load starting from static equilibrium u = K^-1 f0, v=0, a=0,
    # both Newmark and explicit central difference preserve exact equilibrium!
    assert np.allclose(final.u, u_static_ref, atol=1e-10)
    assert np.allclose(final.v, 0.0, atol=1e-10)

    # Thermal steady + transient keeps steady solution K_th * T = q0.
    T_ref = np.linalg.solve(K_th, q0)
    assert np.allclose(final.T, T_ref, atol=1e-10)
    assert len(events) > 5
    assert events[-1].progress == pytest.approx(1.0)


def test_cancel_token_stops_execution() -> None:
    token = CancelToken()
    token.cancel("user abort")
    assert token.is_cancelled is True
    with pytest.raises(CancelledError, match="user abort"):
        token.check()

    K = np.eye(2)
    runner = MultiStepRunner(n_mech_dofs=2, K=K, f_ext=np.array([1.0, 2.0]))
    mid_token = CancelToken()
    seen = []

    def _on_prog(ev: ProgressEvent) -> None:
        seen.append(ev.increment)
        if ev.increment >= 2:
            mid_token.cancel("stop at inc 2")

    with pytest.raises(CancelledError, match="stop at inc 2"):
        runner.run([StaticStep("s", duration=1.0, n_increments=5)], cancel_token=mid_token, on_progress=_on_prog)
    assert seen == [1, 2]


def test_thin_feature_detection_and_max2_refinement() -> None:
    # Thick solid block -> NO thin feature detected.
    thick_labels = np.ones((4, 4, 4), dtype=int)
    rep_thick = detect_thin_features((0.0, 0.0, 0.0), (4.0, 4.0, 4.0), thick_labels)
    assert rep_thick.detected is False
    assert rep_thick.cells == []

    # 1-cell thin wall in X (i == 1 active, rest OUT = -1) -> THIN_FEATURE triggered.
    def _eval_thin(h: THBHierarchy) -> tuple[np.ndarray, list[str]]:
        lvl = h.max_level
        nx, ny, nz = h.cell_shape(lvl)
        u0, u1 = float(h._breaks_u[lvl][0]), float(h._breaks_u[lvl][-1])
        lbl = -np.ones((nx, ny, nz), dtype=int)
        # Physical thin wall spans x in [0.25, 0.50]: only 1 coarse cell at n_cells=4,
        # becomes 2 cells at level 1 -> resolved!
        for i in range(nx):
            xc = u0 + (i + 0.5) / nx * (u1 - u0)
            if 0.25 <= xc <= 0.50:
                lbl[i, :, :] = 1
        return lbl, []

    thb = THBHierarchy.from_uniform((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (4, 4, 4), (2, 2, 2))
    thb_out, reports = auto_refine_thin_features(thb, _eval_thin, max_retries=5)
    assert reports[0].detected is True
    assert any("THIN_FEATURE" in msg for msg in reports[0].issues)
    assert len(reports) <= 3  # Initial + up to 2 refinement passes max.
    assert thb_out.max_level >= 1
    assert reports[-1].detected is False  # After 1 refinement pass, 2 cells cover thin wall!

    # Persistent thin feature (aspect ratio) respects hard cap of 2 retries.
    thb2 = THBHierarchy.from_uniform((0.0, 10.0), (0.0, 1.0), (0.0, 1.0), (2, 4, 4), (2, 2, 2))

    def _eval_always_thin(h: THBHierarchy) -> tuple[np.ndarray, list[str]]:
        lbl = np.ones(h.cell_shape(h.max_level), dtype=int)
        return lbl, []

    thb2_out, rep2 = auto_refine_thin_features(
        thb2, _eval_always_thin, max_retries=10, max_aspect_ratio=5.0
    )
    assert len(rep2) == 3  # Initial check + exact 2 retries max.
    assert thb2_out.max_level == 2


def test_probes_bspline_and_thb_analytical_and_exports(tmp_path: Path) -> None:
    E, nu = 1000.0, 0.3
    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 1.0, 1.0), (2, 2, 2), degree=2)
    cp = box.control_points()

    # Exact linear displacement + temperature field:
    # u = [0.01*x + 0.002*y, -0.003*x, 0.005*z], T = 20 + 50*z
    u_cp = np.column_stack([
        0.01 * cp[:, 0] + 0.002 * cp[:, 1],
        -0.003 * cp[:, 0],
        0.005 * cp[:, 2],
    ])
    t_cp = 20.0 + 50.0 * cp[:, 2]

    # Exact analytical strain Voigt [exx, eyy, ezz, gxy, gyz, gxz]:
    # du_x/dx=0.01, du_y/dy=0, du_z/dz=0.005, du_x/dy + du_y/dx = 0.002 - 0.003 = -0.001
    e_exact = np.array([0.01, 0.0, 0.005, -0.001, 0.0, 0.0], dtype=np.float64)
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    tr_e = e_exact[0] + e_exact[1] + e_exact[2]
    s_exact = np.array([
        lam * tr_e + 2.0 * mu * e_exact[0],
        lam * tr_e + 2.0 * mu * e_exact[1],
        lam * tr_e + 2.0 * mu * e_exact[2],
        mu * e_exact[3],
        mu * e_exact[4],
        mu * e_exact[5],
    ], dtype=np.float64)

    p_probe = PointProbe("center", (1.0, 0.5, 0.5), variables=("U", "S", "E", "MISES", "NT11"))
    l_probe = LineProbe("diag", (0.2, 0.1, 0.1), (1.8, 0.9, 0.9), n_points=5,
                        variables=("DISPLACEMENT", "CAUCHY_STRESS_TENSOR", "GREEN_LAGRANGE_STRAIN_TENSOR", "TEMPERATURE"))
    m_probe = MultiPointProbe("cloud", [[0.4, 0.3, 0.2], [1.6, 0.7, 0.8]], variables=("U", "S", "NT11"))

    s_pt = p_probe.evaluate(box, u_cp=u_cp, t_cp=t_cp, E=E, nu=nu)
    assert np.allclose(s_pt["U"][0], [0.01 * 1.0 + 0.002 * 0.5, -0.003 * 1.0, 0.005 * 0.5], atol=1e-12)
    assert np.allclose(s_pt["E"][0], e_exact, atol=1e-12)
    assert np.allclose(s_pt["S"][0], s_exact, atol=1e-9)
    assert np.allclose(s_pt["NT11"][0], 20.0 + 50.0 * 0.5, atol=1e-12)
    assert s_pt["MISES"][0] > 0.0

    s_line = l_probe.evaluate(box, u_cp=u_cp, t_cp=t_cp, E=E, nu=nu)
    assert s_line.s is not None and s_line.s.shape == (5,)
    for idx, pt in enumerate(s_line.points):
        assert np.allclose(s_line.get("U")[idx], [0.01 * pt[0] + 0.002 * pt[1], -0.003 * pt[0], 0.005 * pt[2]], atol=1e-12)
        assert np.allclose(s_line.get("S")[idx], s_exact, atol=1e-9)
        assert np.allclose(s_line.get("NT11")[idx], 20.0 + 50.0 * pt[2], atol=1e-12)

    s_mp = m_probe.evaluate(box, u_cp=u_cp, t_cp=t_cp, E=E, nu=nu)
    assert s_mp["U"].shape == (2, 3)

    # Also test THBHierarchy evaluation with active-DOF exact linear field!
    thb = THBHierarchy.from_uniform((0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (2, 2, 2), (2, 2, 2))
    thb.mark_box((0.0, 0.0, 0.0), (1.0, 0.5, 0.5), target_level=1)
    thb.enforce_2_to_1_balance()
    # Recover exact linear reproduction coefficients on active THB basis via least squares on sample pts:
    pts_fit = np.column_stack(
        [m.ravel() for m in np.meshgrid(np.linspace(0.05, 1.95, 7), np.linspace(0.05, 0.95, 6), np.linspace(0.05, 0.95, 6), indexing="ij")]
    )
    R_fit = thb.evaluate_thb_basis(pts_fit)
    u_target_fit = np.column_stack([
        0.01 * pts_fit[:, 0] + 0.002 * pts_fit[:, 1],
        -0.003 * pts_fit[:, 0],
        0.005 * pts_fit[:, 2],
    ])
    u_thb_cp, *_ = np.linalg.lstsq(R_fit, u_target_fit, rcond=None)
    t_thb_cp, *_ = np.linalg.lstsq(R_fit, 20.0 + 50.0 * pts_fit[:, 2], rcond=None)
    assert u_thb_cp.shape == (R_fit.shape[1], 3) and R_fit.shape[1] > box.n_ctrl

    s_thb = l_probe.evaluate(thb, u_cp=u_thb_cp, t_cp=t_thb_cp, E=E, nu=nu)
    assert s_thb.mode == "thb"
    assert np.allclose(s_thb["U"], s_line["U"], atol=1e-10)
    assert np.allclose(s_thb["S"], s_line["S"], atol=1e-7)
    assert np.allclose(s_thb["NT11"], s_line["NT11"], atol=1e-10)

    # Export JSON & SQLite with UTF-8 metadata.
    json_file = export_probes_json(
        tmp_path / "probes.json",
        [s_pt, s_line],
        metadata={"label": "한국어 프로브 테스트 — μm"},
    )
    db_file = export_probes_sqlite(tmp_path / "probes.sqlite", [s_pt, s_line, s_mp])

    loaded_json = json.loads(json_file.read_text(encoding="utf-8"))
    assert loaded_json["metadata"]["label"] == "한국어 프로브 테스트 — μm"
    assert len(loaded_json["probes"]) == 2
    assert loaded_json["probes"][1]["n_points"] == 5

    conn = sqlite3.connect(str(db_file))
    try:
        meta_rows = conn.execute("SELECT probe_name, n_points FROM probe_meta").fetchall()
        val_count = conn.execute("SELECT COUNT(*) FROM probe_values").fetchone()[0]
    finally:
        conn.close()
    assert len(meta_rows) == 3
    assert val_count > 0


def test_validate_phase3_solver_fields() -> None:
    clean = JobSpec(
        name="p3",
        lo=(0.0, 0.0, 0.0),
        hi=(1.0, 1.0, 1.0),
        n_cells=(2, 2, 2),
        degree=2,
        fixed_faces=[(0, 0)],
        linear_solver="amgcl",
        enable_solver_log=True,
        solver_log_level="detailed",
        amgcl_params={"krylov_type": "cg", "tolerance": 1e-10},
        memory_budget_ratio=0.70,
    )
    assert validate(clean, []) == []
    roundtrip = JobSpec.from_json(clean.to_json())
    assert roundtrip.to_json() == clean.to_json()

    bad = JobSpec(
        name="bad",
        lo=(0.0, 0.0, 0.0),
        hi=(1.0, 1.0, 1.0),
        n_cells=(2, 2, 2),
        degree=2,
        fixed_faces=[(0, 0)],
        linear_solver="nonexistent_solver",
        memory_budget_ratio=1.5,
    )
    codes = {i.code for i in validate(bad, [])}
    assert "LINEAR_SOLVER" in codes
    assert "MEMORY_BUDGET_RATIO" in codes
