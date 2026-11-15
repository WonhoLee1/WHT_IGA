"""Unit tests for the Kratos-independent ``wht_iga.verification`` package."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from wht_iga.verification import (
    BIMETALLIC_TIE_HEADERS,
    E03_CONV_HEADERS,
    MULTIPLANAR_RESPONSES_HEADERS,
    ROARK_QOI_HEADERS,
    ROARK_REFINEMENT_HEADERS,
    SPHERE_QOI_HEADERS,
    V09_EIGEN_HEADERS,
    analyze_dynamic_history,
    compare_against_tolerances,
    compare_csv_rows,
    filter_participating_modes,
    lame_sphere_exact,
    line_profile_rel_l2,
    list_available_cases,
    load_case_spec,
    load_reference_csv,
    log_log_convergence_rate,
    principal_stresses,
    render_verification_markdown,
    spherical_projections,
    time_window_rms_error,
    write_coreform_csv,
)

WHT_IGA_ROOT = Path(__file__).resolve().parents[2]
COREFORM_ROARK_CSV = (
    WHT_IGA_ROOT
    / "research-cf-IGA"
    / "user-guide"
    / "generated-downloads"
    / "verification-manual"
    / "problems"
    / "solid_mechanics"
    / "linear_elastic_stress"
    / "roark-17p1-7a"
    / "roark-17p1-7a-qoi.csv"
)


def test_list_and_load_existing_case_specs() -> None:
    cases = list_available_cases()
    expected_existing = {f"v{i:02d}" for i in range(4, 13)} | {
        f"e{i:02d}" for i in range(1, 6)
    }
    assert expected_existing.issubset(set(cases))

    for cid in sorted(cases):
        spec = load_case_spec(cid)
        assert spec.case_id.upper() == cid.upper()
        assert spec.title
        assert spec.category in ("verification", "example")
        assert isinstance(spec.geometry, dict)
        assert isinstance(spec.materials, dict)
        assert isinstance(spec.acceptance, dict)

    # Case-insensitive lookup
    v04_upper = load_case_spec("V04")
    v04_lower = load_case_spec("v04")
    assert v04_upper.case_id == v04_lower.case_id
    assert v04_upper.geometry["length_L_mm"] == 100


def test_load_standardized_schema_and_missing_case(tmp_path: Path) -> None:
    custom_json = tmp_path / "v01.json"
    custom_json.write_text(
        '{"meta": {"id": "V01", "category": "verification"}, '
        '"title": "Prescribed displacement cube", "status": "complete", '
        '"model": {"geometry": {"side": 1.0}, "materials": {"E": 1.0}}, '
        '"sweeps": {"degrees": [1, 2, 3]}, '
        '"probes": {"face_grid": [21, 21]}}',
        encoding="utf-8",
    )
    spec = load_case_spec("v01", data_dir=tmp_path)
    assert spec.case_id == "V01"
    assert spec.geometry["side"] == 1.0
    assert spec.sweeps["degrees"] == [1, 2, 3]
    assert spec.probes["face_grid"] == [21, 21]

    with pytest.raises(FileNotFoundError):
        load_case_spec("v99", data_dir=tmp_path)


def test_load_reference_csv_types_and_coreform_file(tmp_path: Path) -> None:
    csv_path = tmp_path / "sample.csv"
    csv_path.write_text(
        "step,step_name,accepted,total_time,note\n"
        "1,cool_to_room,True,0.01,ok\n"
        "2,cool_to_minus_55,False,1.5e-2,cutback\n",
        encoding="utf-8",
    )
    rows = load_reference_csv(csv_path)
    assert len(rows) == 2
    assert rows[0]["step"] == 1 and isinstance(rows[0]["step"], int)
    assert rows[0]["step_name"] == "cool_to_room"
    assert rows[0]["accepted"] is True
    assert rows[1]["accepted"] is False
    assert rows[1]["total_time"] == pytest.approx(0.015)

    if COREFORM_ROARK_CSV.is_file():
        roark_rows = load_reference_csv(COREFORM_ROARK_CSV)
        assert len(roark_rows) > 0
        assert set(ROARK_QOI_HEADERS).issubset(roark_rows[0].keys())
        assert isinstance(roark_rows[0]["h"], float)
        assert isinstance(roark_rows[0]["p"], int)


def test_principal_stresses_single_and_batch() -> None:
    # Rotate diagonal stress diag(100, 20, -40) by an orthogonal matrix
    theta = 0.43
    q = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    diag = np.diag([100.0, 20.0, -40.0])
    rot = q @ diag @ q.T
    voigt = np.array(
        [rot[0, 0], rot[1, 1], rot[2, 2], rot[0, 1], rot[1, 2], rot[0, 2]]
    )

    p_mat, s1_mat = principal_stresses(rot)
    p_vgt, s1_vgt = principal_stresses(voigt)
    np.testing.assert_allclose(p_mat, [100.0, 20.0, -40.0], atol=1e-12)
    np.testing.assert_allclose(p_vgt, [100.0, 20.0, -40.0], atol=1e-12)
    assert float(s1_mat) == pytest.approx(100.0, abs=1e-12)
    assert float(s1_vgt) == pytest.approx(100.0, abs=1e-12)

    # Purely compressive tensor -> sigma1+ == 0.0
    comp = np.array([-10.0, -25.0, -50.0, 0.0, 0.0, 0.0])
    batch = np.vstack([voigt, comp])
    p_batch, s1_batch = principal_stresses(batch)
    assert p_batch.shape == (2, 3)
    np.testing.assert_allclose(s1_batch, [100.0, 0.0], atol=1e-12)


def test_lame_sphere_exact_and_spherical_projections() -> None:
    rng = np.random.default_rng(42)
    # Points in 1/8 sphere shell a=0.9, b=1.0
    dirs = np.abs(rng.normal(size=(50, 3))) + 0.1
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    radii = np.linspace(0.9, 1.0, 50)
    pts = dirs * radii[:, np.newaxis]

    a, b, p_i, p_o, E, nu = 0.9, 1.0, 1.0e-6, 0.0, 1.0, 0.0
    exact = lame_sphere_exact(pts, a=a, b=b, p_i=p_i, p_o=p_o, E=E, nu=nu)

    u_r, s_rr, s_tt = spherical_projections(pts, exact.u_vec, exact.stress_voigt)
    np.testing.assert_allclose(u_r, exact.u_r, rtol=1e-13, atol=1e-19)
    np.testing.assert_allclose(s_rr, exact.sigma_rr, rtol=1e-13, atol=1e-19)
    np.testing.assert_allclose(s_tt, exact.sigma_tt, rtol=1e-13, atol=1e-19)

    # Exact radial stress boundary conditions at r=a and r=b
    assert s_rr[0] == pytest.approx(-p_i, rel=1e-13)
    assert abs(s_rr[-1]) < 1e-18


def test_line_profile_rel_l2_and_convergence_rate() -> None:
    s = np.linspace(0.0, 1.0, 201)
    ref_scalar = np.sin(np.pi * s) + 1.0
    assert line_profile_rel_l2(s, ref_scalar, ref_scalar) == pytest.approx(0.0)

    # O(h^2) perturbation on vector profile
    ref_vec = np.column_stack([ref_scalar, np.cos(np.pi * s)])
    h_vals = np.array([0.2, 0.1, 0.05, 0.025])
    errs = []
    for h in h_vals:
        pert = np.column_stack([(h**2) * np.cos(s), (h**2) * np.sin(s)])
        errs.append(line_profile_rel_l2(s, ref_vec + pert, ref_vec))

    rate = log_log_convergence_rate(h_vals, errs)
    assert rate == pytest.approx(2.0, abs=1e-12)


def test_time_window_rms_error_and_modal_filter() -> None:
    t = np.linspace(0.0, 1.0, 1001)
    q_ref = np.sin(2.0 * np.pi * t) + 2.0
    q_num = q_ref.copy()
    # Large discrepancy outside [0.15, 0.70] must be ignored
    q_num[t < 0.15] += 100.0
    q_num[t > 0.70] -= 50.0
    assert time_window_rms_error(t, q_num, q_ref, 0.15, 0.70) == pytest.approx(0.0)

    # Add 2% constant relative perturbation inside window
    q_num_2pct = q_ref * 1.02
    assert time_window_rms_error(t, q_num_2pct, q_ref, 0.15, 0.70) == pytest.approx(
        0.02, rel=1e-12
    )

    freqs = np.array([80.8, 80.8, 450.0, 505.0, 900.0])
    gamma = np.array([1.0, 0.0, 0.005, 0.35, -0.02])
    idx = filter_participating_modes(freqs, gamma, threshold=0.01)
    np.testing.assert_array_equal(idx, np.array([0, 3, 4]))


def test_analyze_dynamic_history_undamped_and_damped() -> None:
    f_true = 25_000.0
    dt = 1.0e-6  # 40 samples per period
    t = np.arange(0.0, 2.0e-3 + 0.5 * dt, dt)
    t_free_start = 8.0e-5

    # Undamped free vibration
    u_undamped = np.cos(2.0 * np.pi * f_true * t)
    f1, norm_rms, decay, spp = analyze_dynamic_history(
        t, u_undamped, t_free_start=t_free_start, ref_freq=f_true
    )
    assert abs(f1 - f_true) / f_true < 1e-3
    assert norm_rms == pytest.approx(1.0 / np.sqrt(2.0), rel=0.02)
    assert decay == pytest.approx(1.0, abs=0.01)
    assert spp == pytest.approx(40.0, rel=1e-6)

    # Damped free vibration
    zeta = 0.02
    omega_d = 2.0 * np.pi * f_true * np.sqrt(1.0 - zeta**2)
    u_damped = np.exp(-zeta * 2.0 * np.pi * f_true * (t - t_free_start)) * np.cos(
        omega_d * (t - t_free_start)
    )
    f1_d, _, decay_d, _ = analyze_dynamic_history(
        t, u_damped, t_free_start=t_free_start
    )
    expected_f_d = omega_d / (2.0 * np.pi)
    assert abs(f1_d - expected_f_d) / expected_f_d < 1e-3
    assert 0.0 < decay_d < 0.5


def test_comparator_and_reporter_roundtrip(tmp_path: Path) -> None:
    report = compare_against_tolerances(
        "V06",
        measured={"u_L2": 3.65e-5, "sigma_rr_L2": 2.11e-3},
        reference={"u_L2": 3.65e-5, "sigma_rr_L2": 2.10e-3},
        tolerances={"u_L2": 5.0e-5, "sigma_rr_L2": 0.01},
        units={"u_L2": "-", "sigma_rr_L2": "-"},
    )
    assert report.passed is True
    assert len(report.checks) == 2

    md = render_verification_markdown(
        report, metadata={"dofs": 57300, "phase_gate": "G2"}
    )
    assert "# Verification Report: V06" in md
    assert "**Status:** PASS" in md
    assert "| u_L2 |" in md
    assert "| dofs | 57300 |" in md

    # CSV write + read roundtrip for Coreform headers
    csv_out = tmp_path / "roark-17p1-7a-qoi.csv"
    row = {
        "test_case": "immersed_h1p2",
        "h": 3.0,
        "p": 2,
        "peak_stress": 4.32e-6,
        "peak_rel_err_vs_roark": 0.00058,
        "peak_rel_err_vs_reference": 0.00643,
        "stress_profile_rel_l2": 0.004,
        "ux_profile_rel_l2": 0.001,
        "uy_profile_rel_l2": 0.008,
    }
    write_coreform_csv(csv_out, ROARK_QOI_HEADERS, [row])
    # Ensure UTF-8 without BOM
    assert not csv_out.read_bytes().startswith(b"\xef\xbb\xbf")

    loaded = load_reference_csv(csv_out)
    csv_report = compare_csv_rows(
        loaded,
        [row],
        key_cols=["test_case", "h", "p"],
        value_cols=["peak_stress", "stress_profile_rel_l2"],
        rel_tol=1e-12,
        case_id="V02",
    )
    assert csv_report.passed is True

    # Verify all Coreform header tuples are non-empty and unique
    for headers in (
        ROARK_QOI_HEADERS,
        ROARK_REFINEMENT_HEADERS,
        BIMETALLIC_TIE_HEADERS,
        MULTIPLANAR_RESPONSES_HEADERS,
        SPHERE_QOI_HEADERS,
        V09_EIGEN_HEADERS,
        E03_CONV_HEADERS,
    ):
        assert len(headers) == len(set(headers)) and len(headers) >= 4
