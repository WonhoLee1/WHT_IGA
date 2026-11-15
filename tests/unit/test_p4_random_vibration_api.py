"""Comprehensive unit & analytical benchmark tests for ``wht_iga.dynamics.random_vibration``.

Verifies the pure-Python PSD random vibration API against closed-form analytical
solutions without mocks or stubs:
    * PSD curves (``PSDTable`` 4 interpolation scales x 3 boundary policies,
      ``WhiteNoisePSD``, ``OctaveBandPSD`` IEC 61260 1/1 & 1/3 octave dB).
    * Damping models (constant, list, Rayleigh, exact hysteretic structural).
    * Multi-input ``CrossPSDMatrix`` (uncorrelated, fully correlated, coherence +
      phase, ``C - iQ`` co/quad, magnitude/phase, Hermitian + PSD validation).
    * Frequency grids (uniform, log, Abaqus sweep, half-power tangent adaptive).
    * Analytical Miles SDOF relative & absolute acceleration equations.
    * 2-DOF damped cross-correlated analytical FRF vs ``ModalRandomVibration``
      and ``DirectFRFRandomVibration`` parity + static residual flexibility.
    * Spectral moments, Segalman-Fulcher von Mises RMS & rotation invariance,
      and Basquin S-N fatigue (``NarrowBandRayleigh``, ``Steinberg3Band``,
      ``Dirlik`` PDF/closed-form integral/narrowband limit).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.integrate as integrate
import scipy.linalg as la
from scipy.special import gamma as gamma_func

import wht_iga.dynamics as dyn
from wht_iga.dynamics.case_file import parse_case
from wht_iga.dynamics.consistent_mass import assemble_consistent_mass_from_bspline, solve_eigenmodes
from wht_iga.dynamics.random_vibration import (
    BaseExcitation,
    ConstantModalDamping,
    CrossPSDMatrix,
    DirectFRFRandomVibration,
    Dirlik,
    ModalDampingList,
    ModalRandomVibration,
    NarrowBandRayleigh,
    NodalForceExcitation,
    OctaveBandPSD,
    PSDTable,
    RayleighDamping,
    SNCurve,
    Steinberg3Band,
    StructuralDamping,
    WhiteNoisePSD,
    abaqus_frequency_grid,
    adaptive_modal_frequency_grid,
    compute_residual_flexibility,
    iec61260_band_edges,
    iec61260_center_frequencies,
    log_frequency_grid,
    random_vibration_from_case,
    segalman_von_mises_rms,
    spectral_moments,
    uniform_frequency_grid,
    von_mises_psd,
)
from wht_iga.dynamics.response_spectrum import VON_MISES_A, random_response_from_case


# ---------------------------------------------------------------------------
# 1. PSD Curves
# ---------------------------------------------------------------------------


def test_psd_table_interpolation_scales_and_extrapolation() -> None:
    """Verify all four interpolation scales and zero/clamp/extrapolate policies."""
    f_tab = np.array([10.0, 100.0])
    v_tab = np.array([1.0, 100.0])

    # log-log: v(f) = 0.01 * f^2 -> at 31.6227766 (= 10^1.5), v = 10.0
    f_mid_geom = math.sqrt(10.0 * 100.0)
    t_ll = PSDTable(f_tab, v_tab, scale="log-log", extrapolation="zero")
    assert t_ll(f_mid_geom) == pytest.approx(10.0, rel=1e-13)
    assert t_ll(5.0) == 0.0 and t_ll(200.0) == 0.0

    t_ll_clamp = PSDTable(f_tab, v_tab, scale="log-log", extrapolation="clamp")
    assert t_ll_clamp(5.0) == pytest.approx(1.0, rel=1e-14)
    assert t_ll_clamp(200.0) == pytest.approx(100.0, rel=1e-14)

    t_ll_ext = PSDTable(f_tab, v_tab, scale="log-log", extrapolation="extrapolate")
    # Slope in log-log is 2.0 -> at f=200: 100 * (2)^2 = 400; at f=5: 1 * (0.5)^2 = 0.25
    assert t_ll_ext(5.0) == pytest.approx(0.25, rel=1e-13)
    assert t_ll_ext(200.0) == pytest.approx(400.0, rel=1e-13)

    # log-lin: linear in y, log in f -> at geometric mean f_mid_geom, y = 0.5*(1+100) = 50.5
    t_loglin = PSDTable(f_tab, v_tab, scale="log-lin", extrapolation="extrapolate")
    assert t_loglin(f_mid_geom) == pytest.approx(50.5, rel=1e-13)
    assert t_loglin(1000.0) == pytest.approx(199.0, rel=1e-13)

    # lin-log: log in y, linear in f -> at arithmetic mean 55 Hz, y = sqrt(1 * 100) = 10
    t_linlog = PSDTable(f_tab, v_tab, scale="lin-log", extrapolation="extrapolate")
    assert t_linlog(55.0) == pytest.approx(10.0, rel=1e-13)
    assert t_linlog(190.0) == pytest.approx(10000.0, rel=1e-12)

    # lin-lin: linear in both
    t_linlin = PSDTable(f_tab, v_tab, scale="lin-lin", extrapolation="zero")
    assert t_linlin(55.0) == pytest.approx(50.5, rel=1e-13)
    assert t_linlin(0.0) == 0.0

    # Complex table in lin-lin
    t_c = PSDTable([10.0, 20.0], [2.0, 4.0], imag_values=[-1.0, 1.0], scale="lin-lin")
    assert t_c(15.0) == pytest.approx(3.0 + 0.0j, abs=1e-14)


def test_white_noise_psd_band_limits() -> None:
    """WhiteNoisePSD returns W0 inside [f_min, f_max] and 0 outside."""
    wn = WhiteNoisePSD(0.05, f_min=10.0, f_max=500.0)
    got = wn(np.array([5.0, 10.0, 100.0, 500.0, 501.0]))
    np.testing.assert_allclose(got, [0.0, 0.05, 0.05, 0.05, 0.0], atol=1e-15)
    with pytest.raises(ValueError, match=">="):
        WhiteNoisePSD(-1.0)
    with pytest.raises(ValueError, match="f_max > f_min"):
        WhiteNoisePSD(1.0, f_min=100.0, f_max=50.0)


def test_octave_band_psd_iec61260_and_db_integration() -> None:
    """Verify IEC 61260 1/1 and 1/3 octave center freqs, dB conversion, and exact band integration."""
    fc_1_3 = iec61260_center_frequencies(100.0, 1000.0, fraction=3)
    # From 100 Hz (band 20) to 1000 Hz (band 30) -> 11 bands inclusive
    assert fc_1_3.size == 11
    assert fc_1_3[0] == pytest.approx(100.0, rel=1e-12)
    assert fc_1_3[-1] == pytest.approx(1000.0, rel=1e-12)

    fc_1_1 = iec61260_center_frequencies(125.0, 1000.0, fraction=1)
    assert fc_1_1[-1] == pytest.approx(1000.0, rel=1e-12)
    assert fc_1_1.size == 4  # 125.9, 251.2, 501.2, 1000.0 (IEC base-10)

    # psd_db mode: 0 dB = 1 unit^2/Hz -> +20 dB = 100, +10 dB = 10; custom ref_psd=2.0 -> 200, 20
    levels = np.linspace(0.0, 20.0, fc_1_3.size)
    oct_psd = OctaveBandPSD(fc_1_3, levels, fraction=3, ref_psd=1.0, level_kind="psd_db")
    tbl = oct_psd.to_psd_table()
    assert isinstance(tbl, PSDTable)
    assert tbl(100.0) == pytest.approx(1.0, rel=1e-12)
    assert tbl(1000.0) == pytest.approx(100.0, rel=1e-12)

    oct_ref = OctaveBandPSD(fc_1_3, levels, fraction=3, ref_psd=2.5, level_kind="psd_db")
    assert oct_ref(1000.0) == pytest.approx(250.0, rel=1e-12)

    # ANSI S1.11 band numbers helper
    ansi = OctaveBandPSD.from_ansi_bands([20, 30], [10.0, 20.0], fraction=3)
    assert ansi.center_freqs_hz[0] == pytest.approx(100.0, rel=1e-12)
    assert ansi.center_freqs_hz[1] == pytest.approx(1000.0, rel=1e-12)
    assert ansi(100.0) == pytest.approx(10.0, rel=1e-12)
    assert ansi(1000.0) == pytest.approx(100.0, rel=1e-12)

    # DB octave integration test: band_db with band_constant mode integrates exactly to sum(ref * 10^(L/10))
    band_levels_db = np.array([80.0, 86.0, 90.0, 84.0])
    bands = iec61260_center_frequencies(100.0, 200.0, fraction=3)
    assert bands.size == band_levels_db.size
    ref_p2 = (20e-6) ** 2
    oct_band = OctaveBandPSD(
        bands,
        band_levels_db,
        fraction=3,
        ref_psd=ref_p2,
        level_kind="band_db",
        mode="band_constant",
    )
    f_lo, f_hi = iec61260_band_edges(bands, fraction=3)
    # Build per-band quadrature grid cleanly inside each band to avoid step-jump quadrature error
    total_ms = 0.0
    for k in range(bands.size):
        fk_grid = np.linspace(f_lo[k] * (1.0 + 1e-12), f_hi[k] * (1.0 - 1e-12), 51)
        total_ms += float(integrate.trapezoid(oct_band(fk_grid), fk_grid))
    expected_ms = float(np.sum(ref_p2 * (10.0 ** (band_levels_db / 10.0))))
    assert total_ms == pytest.approx(expected_ms, rel=1e-10)


# ---------------------------------------------------------------------------
# 2. Damping Models & Frequency Grids
# ---------------------------------------------------------------------------


def test_damping_models_and_structural_closed_form() -> None:
    """Verify Rayleigh, list, constant, and exact hysteretic structural white-noise variance."""
    w_modes = np.array([10.0, 50.0, 200.0])
    c_damp = ConstantModalDamping(0.03)
    np.testing.assert_allclose(c_damp.equivalent_viscous_ratios(w_modes), 0.03)

    l_damp = ModalDampingList([0.01, 0.02, 0.05])
    np.testing.assert_allclose(l_damp.equivalent_viscous_ratios(w_modes), [0.01, 0.02, 0.05])

    alpha, beta = 1.2, 2e-4
    r_damp = RayleighDamping(alpha=alpha, beta=beta)
    np.testing.assert_allclose(
        r_damp.equivalent_viscous_ratios(w_modes),
        0.5 * alpha / w_modes + 0.5 * beta * w_modes,
        rtol=1e-14,
    )

    # Structural damping SDOF closed-form variance:
    # int_0^inf |1 / (w0^2 - w^2 + i*g*w0^2)|^2 W0 df = W0 / (2*sqrt(2)*g*w0^3*sqrt(1+sqrt(1+g^2)))
    f0, g_loss, w0_psd = 100.0, 0.04, 2.5
    wn = 2.0 * math.pi * f0
    solver = ModalRandomVibration(
        frequencies_hz=[f0],
        phi=np.array([[1.0]]),
        damping=StructuralDamping(g_loss),
    )
    exc = NodalForceExcitation([1.0], WhiteNoisePSD(w0_psd))
    grid = adaptive_modal_frequency_grid([f0], StructuralDamping(g_loss), 1.0, 1500.0, n_base=200, half_power_multiplier=40.0, points_per_mode=201)
    res = solver.solve(exc, freqs_hz=grid)
    closed_var = w0_psd / (2.0 * math.sqrt(2.0) * g_loss * (wn**3) * math.sqrt(1.0 + math.sqrt(1.0 + g_loss**2)))
    # Correct for tiny finite upper/lower truncation tail
    assert res.ru[0] ** 2 == pytest.approx(closed_var, rel=5e-4)


def test_frequency_grid_builders_and_adaptive_efficiency() -> None:
    """Adaptive half-power clustering resolves sharp zeta=0.005 peak far better than uniform."""
    f0, zeta, w0_psd = 120.0, 0.005, 0.5
    wn = 2.0 * math.pi * f0
    exact_var = w0_psd / (8.0 * zeta * (wn**3))

    solver = ModalRandomVibration([f0], np.array([[1.0]]), damping=zeta)
    exc = NodalForceExcitation([1.0], WhiteNoisePSD(w0_psd))

    g_adapt = adaptive_modal_frequency_grid([f0], zeta, 10.0, 1000.0, n_base=81, half_power_multiplier=35.0, points_per_mode=141)
    g_unif = uniform_frequency_grid(10.0, 1000.0, g_adapt.size)
    g_log = log_frequency_grid(10.0, 1000.0, 100)
    g_abq = abaqus_frequency_grid([f0], 10.0, 1000.0, n_points=20, bias=3.0)

    assert g_log.size == 100 and np.all(np.diff(g_log) > 0.0)
    assert np.any(np.isclose(g_abq, f0, rtol=1e-12))

    err_adapt = abs(solver.solve(exc, freqs_hz=g_adapt).ru[0] ** 2 / exact_var - 1.0)
    err_unif = abs(solver.solve(exc, freqs_hz=g_unif).ru[0] ** 2 / exact_var - 1.0)
    assert err_adapt < 5e-4
    assert err_unif > 20.0 * err_adapt  # uniform misses the sharp zeta=0.005 peak badly


# ---------------------------------------------------------------------------
# 3. CrossPSDMatrix Validation & Representations
# ---------------------------------------------------------------------------


def test_cross_psd_matrix_modes_and_validation() -> None:
    """Verify uncorrelated, fully correlated, coherence, co/quad, mag/phase, and PSD checks."""
    freqs = np.array([10.0, 50.0, 100.0])
    e0 = BaseExcitation("x", WhiteNoisePSD(2.0), gravity_scale=1.0)
    e1 = BaseExcitation("y", WhiteNoisePSD(8.0), gravity_scale=1.0)

    # Uncorrelated
    s_unc = CrossPSDMatrix.uncorrelated([e0, e1]).evaluate(freqs)
    assert s_unc.shape == (2, 2, 3)
    np.testing.assert_allclose(s_unc[0, 1], 0.0)
    np.testing.assert_allclose(s_unc[0, 0], 2.0)
    np.testing.assert_allclose(s_unc[1, 1], 8.0)

    # Fully correlated with 90-deg phase shift -> S_01 = sqrt(2*8) * exp(-i*pi/2) = -4j
    s_fc = CrossPSDMatrix.fully_correlated([e0, e1], phases=[0.0, 90.0], phase_unit="deg").evaluate(freqs)
    np.testing.assert_allclose(s_fc[0, 1], -4.0j, atol=1e-13)
    np.testing.assert_allclose(s_fc[1, 0], +4.0j, atol=1e-13)

    # Coherence gamma^2 = 0.25 -> |S_01| = sqrt(0.25 * 2 * 8) = 2.0
    s_coh = CrossPSDMatrix.from_coherence([e0, e1], {(0, 1): 0.25}, {(0, 1): 60.0}, phase_unit="deg").evaluate(freqs)
    np.testing.assert_allclose(s_coh[0, 1], 2.0 * np.exp(1j * math.pi / 3.0), rtol=1e-13)

    # Co-quad C - iQ
    s_cq = CrossPSDMatrix.from_co_quad([e0, e1], {(0, 1): 1.5}, {(0, 1): 0.8}).evaluate(freqs)
    np.testing.assert_allclose(s_cq[0, 1], 1.5 - 0.8j, rtol=1e-14)
    np.testing.assert_allclose(s_cq[1, 0], 1.5 + 0.8j, rtol=1e-14)

    # Magnitude + phase
    s_mp = CrossPSDMatrix.from_magnitude_phase([e0, e1], {(0, 1): 2.0}, {(0, 1): -math.pi / 4}).evaluate(freqs)
    np.testing.assert_allclose(s_mp[0, 1], 2.0 * np.exp(-1j * math.pi / 4), rtol=1e-14)

    # Invalid coherence > 1 fails
    with pytest.raises(ValueError, match="Coherence"):
        CrossPSDMatrix.from_coherence([e0, e1], {(0, 1): 1.2}).evaluate(freqs)

    # Non-Hermitian array fails
    bad_herm = np.zeros((2, 2, 3), dtype=complex)
    bad_herm[0, 0] = 2.0
    bad_herm[1, 1] = 2.0
    bad_herm[0, 1] = 1.0 + 1.0j
    bad_herm[1, 0] = 1.0 + 1.0j  # not conjugate!
    with pytest.raises(ValueError, match="Hermitian"):
        CrossPSDMatrix.from_array(freqs, bad_herm)

    # Indefinite Hermitian matrix (negative eigenvalue) fails
    bad_psd = np.zeros((2, 2, 3), dtype=complex)
    bad_psd[0, 0] = 1.0
    bad_psd[1, 1] = 1.0
    bad_psd[0, 1] = 1.5  # determinant = 1 - 2.25 < 0 -> lambda_min = -0.5
    bad_psd[1, 0] = 1.5
    with pytest.raises(ValueError, match="positive semi-definite"):
        CrossPSDMatrix.from_array(freqs, bad_psd)


# ---------------------------------------------------------------------------
# 4. Miles SDOF & Multi-DOF Analytical Verification
# ---------------------------------------------------------------------------


def test_miles_sdof_relative_and_absolute_acceleration() -> None:
    """Verify SDOF white-noise relative disp/vel and Miles' absolute acceleration equation."""
    fn = 80.0
    zeta = 0.005  # Q = 100
    q_factor = 1.0 / (2.0 * zeta)
    w0 = 0.04  # (m/s^2)^2 / Hz
    wn = 2.0 * math.pi * fn

    solver = ModalRandomVibration(
        frequencies_hz=[fn],
        phi=np.array([[1.0]]),
        M_modal=np.array([1.0]),
        damping=zeta,
        participation_factors=np.array([[1.0, 0.0, 0.0]]),
        influence_vectors=np.array([[1.0, 0.0, 0.0]]),
    )
    exc = BaseExcitation("x", WhiteNoisePSD(w0), kind="acceleration")
    grid = adaptive_modal_frequency_grid(
        [fn], zeta, 0.5, 2500.0, n_base=400, half_power_multiplier=50.0, points_per_mode=241
    )

    # Relative response: sigma_u^2 = W0 / (8 * zeta * wn^3), sigma_v^2 = W0 / (8 * zeta * wn)
    res_rel = solver.solve(exc, freqs_hz=grid, response_frame="relative")
    assert res_rel.ru[0] ** 2 == pytest.approx(w0 / (8.0 * zeta * (wn**3)), rel=5e-4)
    assert res_rel.rv[0] ** 2 == pytest.approx(w0 / (8.0 * zeta * wn), rel=5e-4)

    # Absolute acceleration Miles' equation: sigma_a_abs = sqrt((pi/2) * fn * Q * W0 * (1 + 1/Q^2))
    res_abs = solver.solve(exc, freqs_hz=grid, response_frame="absolute")
    miles_exact = math.sqrt(0.5 * math.pi * fn * q_factor * w0 * (1.0 + 1.0 / (q_factor**2)))
    miles_classic = math.sqrt(0.5 * math.pi * fn * q_factor * w0)
    assert res_abs.ra[0] == pytest.approx(miles_exact, rel=5e-4)
    assert res_abs.ra[0] == pytest.approx(miles_classic, rel=1e-3)


def test_damped_2dof_correlated_modal_vs_direct_frf_and_analytical() -> None:
    """Damped 2-DOF under complex cross-correlated forces: Modal == DirectFRF == exact analytic."""
    m1, m2 = 2.0, 1.5
    k1, kc, k2 = 1.8e5, 4.5e4, 1.2e5
    M = np.diag([m1, m2])
    K = np.array([[k1 + kc, -kc], [-kc, kc + k2]], dtype=float)
    alpha, beta = 3.0, 1.5e-5
    C = alpha * M + beta * K

    # Generalized eigenproblem with mass-orthonormal Phi
    eigvals, eigvecs = la.eigh(K, M)
    omega_modes = np.sqrt(eigvals)
    f_modes = omega_modes / (2.0 * np.pi)
    damping = RayleighDamping(alpha=alpha, beta=beta)

    # Non-trivial frequency-dependent complex cross-correlated force excitations
    p0 = PSDTable([5.0, 200.0], [10.0, 2.0], scale="log-log", extrapolation="clamp")
    p1 = PSDTable([5.0, 200.0], [4.0, 8.0], scale="log-log", extrapolation="clamp")
    e0 = NodalForceExcitation([1.0, 0.0], p0)
    e1 = NodalForceExcitation([0.0, 1.0], p1)
    cpsd = CrossPSDMatrix.from_coherence([e0, e1], {(0, 1): 0.6}, {(0, 1): 35.0}, phase_unit="deg")

    grid = adaptive_modal_frequency_grid(f_modes, damping, 5.0, 200.0, n_base=120, points_per_mode=101)

    # 1. ModalRandomVibration
    modal_solver = ModalRandomVibration(
        frequencies_hz=f_modes,
        phi=eigvecs,
        M_full=M,
        damping=damping,
    )
    res_modal = modal_solver.solve(cpsd, freqs_hz=grid)

    # 2. DirectFRFRandomVibration using exact physical H(w) = (-w^2 M + i w C + K)^{-1}
    h_direct = np.empty((2, 2, grid.size), dtype=complex)
    for k, f in enumerate(grid):
        w = 2.0 * math.pi * f
        dyn_k = K - (w * w) * M + 1j * w * C
        h_direct[:, :, k] = np.linalg.inv(dyn_k)
    direct_solver = DirectFRFRandomVibration(grid, h_direct)
    res_direct = direct_solver.solve(cpsd)

    # 3. Independent scalar loop reference
    s_in = cpsd.evaluate(grid)
    s_uu_ref = np.empty((2, 2, grid.size), dtype=complex)
    for k in range(grid.size):
        hk = h_direct[:, :, k]
        s_uu_ref[:, :, k] = hk @ s_in[:, :, k] @ hk.conj().T

    np.testing.assert_allclose(res_modal.disp_psd, res_direct.disp_psd, rtol=1e-11, atol=1e-18)
    np.testing.assert_allclose(res_modal.vel_psd, res_direct.vel_psd, rtol=1e-11, atol=1e-16)
    np.testing.assert_allclose(res_modal.acc_psd, res_direct.acc_psd, rtol=1e-11, atol=1e-14)
    np.testing.assert_allclose(res_modal.cross_psd([0, 1]), s_uu_ref, rtol=1e-11, atol=1e-18)
    np.testing.assert_allclose(res_direct.cross_psd([0, 1]), s_uu_ref, rtol=1e-12, atol=1e-18)
    np.testing.assert_allclose(res_modal.ru, res_direct.ru, rtol=1e-11)
    assert res_modal.displacement_moments(0).rms == pytest.approx(res_direct.ru[0], rel=1e-12)


def test_static_residual_flexibility_recovers_truncated_static_limit() -> None:
    """Truncated 3-DOF modal basis with compute_residual_flexibility recovers exact quasi-static FRF."""
    M = np.diag([1.0, 1.0, 1.0])
    K = np.array([[2000.0, -1000.0, 0.0], [-1000.0, 2000.0, -1000.0], [0.0, -1000.0, 2000.0]])
    eigvals, eigvecs = la.eigh(K, M)
    omega_all = np.sqrt(eigvals)
    f_all = omega_all / (2.0 * np.pi)

    # Keep only mode 1 (truncate modes 2 and 3)
    phi_1 = eigvecs[:, :1]
    w_1 = omega_all[:1]
    g_res = compute_residual_flexibility(K, phi_1, w_1)

    solver_full = ModalRandomVibration(f_all, eigvecs, M_full=M, damping=0.02)
    solver_res = ModalRandomVibration(f_all[:1], phi_1, M_full=M, damping=0.02, residual_flexibility=g_res)
    solver_trunc = ModalRandomVibration(f_all[:1], phi_1, M_full=M, damping=0.02)

    exc = NodalForceExcitation([0.0, 0.0, 1.0], WhiteNoisePSD(1.0))
    f_low = np.array([0.01, 0.1, 0.5])  # well below f1 (~3.85 Hz) -> quasi-static
    psd_full = solver_full.solve(exc, freqs_hz=f_low).disp_psd
    psd_res = solver_res.solve(exc, freqs_hz=f_low).disp_psd
    psd_trunc = solver_trunc.solve(exc, freqs_hz=f_low).disp_psd

    # Residual flexibility matches full-basis quasi-static PSD far closer than bare truncated modal basis
    np.testing.assert_allclose(psd_res[:, 0], psd_full[:, 0], rtol=1e-4)
    assert np.max(np.abs(psd_trunc[:, 0] - psd_full[:, 0]) / psd_full[:, 0]) > 0.1


# ---------------------------------------------------------------------------
# 5. Spectral Moments, Segalman-Fulcher von Mises, & Fatigue
# ---------------------------------------------------------------------------


def test_spectral_moments_closed_form_polynomial() -> None:
    """Linear PSD S(f) = a + b*f on [f1, f2]: check exact analytical lambda_0,1,2,4 & rates."""
    f1, f2 = 10.0, 60.0
    a, b = 0.5, 0.2
    grid = np.linspace(f1, f2, 4001)
    s = a + b * grid

    mom = spectral_moments(grid, s)
    # Exact integrals of (a + b f) * f^n from f1 to f2
    exact = {}
    for n in (0, 1, 2, 4):
        exact[n] = a * (f2 ** (n + 1) - f1 ** (n + 1)) / (n + 1) + b * (f2 ** (n + 2) - f1 ** (n + 2)) / (n + 2)

    assert mom.lambda_0 == pytest.approx(exact[0], rel=1e-7)
    assert mom.lambda_1 == pytest.approx(exact[1], rel=1e-7)
    assert mom.lambda_2 == pytest.approx(exact[2], rel=1e-7)
    assert mom.lambda_4 == pytest.approx(exact[4], rel=1e-6)
    assert mom.rms == pytest.approx(math.sqrt(exact[0]), rel=1e-7)
    assert mom.peak_2sigma == pytest.approx(2.0 * mom.rms, rel=1e-14)
    assert mom.peak_3sigma == pytest.approx(3.0 * mom.rms, rel=1e-14)
    assert mom.f0 == pytest.approx(math.sqrt(exact[2] / exact[0]), rel=1e-7)
    assert mom.fp == pytest.approx(math.sqrt(exact[4] / exact[2]), rel=1e-7)
    assert mom.gamma2 == pytest.approx(exact[2] / math.sqrt(exact[0] * exact[4]), rel=1e-7)
    assert 0.0 < mom.gamma2 <= 1.0


def test_segalman_von_mises_complex_off_diagonal_and_3d_rotation_invariance() -> None:
    """Verify Segalman-Fulcher RMISES with complex cross-stress PSD and SO(3) rotation invariance."""
    freqs = np.linspace(10.0, 110.0, 201)
    # Build a valid Hermitian PSD 6x6 stress CPSD S_ss(f) = B(f) B(f)^H with complex off-diagonals
    b1 = np.array([3.0, -1.0, 0.5, 1.2 + 0.8j, -0.4j, 0.7 - 0.3j])
    b2 = np.array([0.0, 2.0, 1.0, -0.5 + 0.2j, 1.1, 0.3 + 0.9j])
    s_base = np.outer(b1, b1.conj()) + 0.5 * np.outer(b2, b2.conj())
    # Verify s_base is strictly Hermitian with non-zero imaginary skew-symmetric part
    assert np.max(np.abs(s_base.imag)) > 0.5
    s_cpsd = s_base[:, :, None] * np.ones_like(freqs)[None, None, :]

    rmises_cpsd = segalman_von_mises_rms(s_cpsd, freqs)
    cov_real = np.real(integrate.trapezoid(s_cpsd, freqs, axis=-1))
    rmises_cov = segalman_von_mises_rms(cov_real)
    # Manual trace formula: sqrt(tr(A_vM * Re(Lambda_0)))
    manual = math.sqrt(float(np.trace(VON_MISES_A @ cov_real)))
    assert rmises_cpsd == pytest.approx(manual, rel=1e-13)
    assert rmises_cov == pytest.approx(manual, rel=1e-13)

    # 3D SO(3) rotation invariance of Segalman-Fulcher von Mises RMS:
    # Rotate stress tensor sigma' = R sigma R^T via 6x6 Voigt transformation M_sigma
    ax, ay, az = 0.4, -0.7, 1.1
    rx = np.array([[1, 0, 0], [0, math.cos(ax), -math.sin(ax)], [0, math.sin(ax), math.cos(ax)]])
    ry = np.array([[math.cos(ay), 0, math.sin(ay)], [0, 1, 0], [-math.sin(ay), 0, math.cos(ay)]])
    rz = np.array([[math.cos(az), -math.sin(az), 0], [math.sin(az), math.cos(az), 0], [0, 0, 1]])
    R = rz @ ry @ rx

    # Build 6x6 Voigt rotation matrix T_voigt: [xx, yy, zz, xy, yz, xz]
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]
    T_v = np.zeros((6, 6), dtype=float)
    for row, (i, j) in enumerate(pairs):
        for col, (p, q) in enumerate(pairs):
            if p == q:
                T_v[row, col] = R[i, p] * R[j, p]
            else:
                T_v[row, col] = R[i, p] * R[j, q] + R[i, q] * R[j, p]

    s_rot = np.einsum("ai,ijf,bj->abf", T_v, s_cpsd, T_v)
    rmises_rot = segalman_von_mises_rms(s_rot, freqs)
    assert rmises_rot == pytest.approx(rmises_cpsd, rel=1e-12)


def test_fatigue_rayleigh_steinberg_and_dirlik_quadrature_and_narrowband_limit() -> None:
    """Verify Rayleigh, Steinberg 3-band, Dirlik closed-form == PDF quadrature, and narrowband convergence."""
    sn = SNCurve(C=1.0e15, b=5.0, stress_definition="amplitude")

    # 1. Bimodal wideband spectrum -> test Dirlik PDF integrates to 1 & closed-form matches quadrature
    f_grid = np.linspace(10.0, 300.0, 6001)
    psd_wide = (
        100.0 * np.exp(-0.5 * ((f_grid - 40.0) / 5.0) ** 2)
        + 40.0 * np.exp(-0.5 * ((f_grid - 180.0) / 12.0) ** 2)
    )
    mom_wide = spectral_moments(f_grid, psd_wide)
    assert 0.3 < mom_wide.gamma2 < 0.85

    s_amps = np.linspace(0.0, 12.0 * mom_wide.rms, 20001)
    pdf_dirlik = Dirlik.pdf(s_amps, mom_wide)
    assert float(integrate.trapezoid(pdf_dirlik, s_amps)) == pytest.approx(1.0, rel=1e-6)

    res_dirlik = Dirlik().evaluate(mom_wide, sn, duration_seconds=3600.0)
    quad_damage_rate = (mom_wide.fp / sn.C) * float(integrate.trapezoid((s_amps**sn.b) * pdf_dirlik, s_amps))
    assert res_dirlik.damage_rate == pytest.approx(quad_damage_rate, rel=1e-6)
    assert res_dirlik.life_seconds == pytest.approx(1.0 / res_dirlik.damage_rate, rel=1e-14)
    assert res_dirlik.damage_in_duration == pytest.approx(3600.0 * res_dirlik.damage_rate, rel=1e-14)

    # 2. NarrowBandRayleigh closed-form vs Rayleigh PDF quadrature
    res_ray = NarrowBandRayleigh().evaluate(mom_wide, sn)
    pdf_ray = NarrowBandRayleigh.pdf(s_amps, mom_wide)
    assert float(integrate.trapezoid(pdf_ray, s_amps)) == pytest.approx(1.0, rel=1e-7)
    ray_manual = (mom_wide.f0 / sn.C) * ((math.sqrt(2.0) * mom_wide.rms) ** sn.b) * float(gamma_func(1.0 + 0.5 * sn.b))
    assert res_ray.damage_rate == pytest.approx(ray_manual, rel=1e-13)

    # 3. Steinberg3Band exact hand formula (68.31% @ 1s, 27.10% @ 2s, 4.33% @ 3s)
    res_st = Steinberg3Band().evaluate(mom_wide, sn)
    st_manual = (mom_wide.f0 / sn.C) * (mom_wide.rms**sn.b) * (0.6831 + 0.2710 * (2.0**5) + 0.0433 * (3.0**5))
    assert res_st.damage_rate == pytest.approx(st_manual, rel=1e-13)

    # 4. Dirlik narrowband limit converges smoothly to NarrowBandRayleigh as bandwidth -> 0
    f0 = 100.0
    rel_diffs = []
    for bw in (5.0, 1.0, 0.2, 0.02):
        fg = np.linspace(f0 - 8.0 * bw, f0 + 8.0 * bw, 2001)
        s_nb = 50.0 * np.exp(-0.5 * ((fg - f0) / bw) ** 2)
        mom_nb = spectral_moments(fg, s_nb)
        d_dir = Dirlik().damage_rate(mom_nb, sn)
        d_ray = NarrowBandRayleigh().damage_rate(mom_nb, sn)
        rel_diffs.append(abs(d_dir / d_ray - 1.0))
    assert rel_diffs[0] > rel_diffs[1] > rel_diffs[2] > rel_diffs[3]
    assert rel_diffs[-1] < 1e-4


# ---------------------------------------------------------------------------
# 6. Case-File Bridge & Tiny Consistent-Mass Modal Bridge
# ---------------------------------------------------------------------------


def test_random_vibration_from_case_matches_existing_engine() -> None:
    """Bridge helper reproduces random_response_from_case ru/rv/ra/rmises to 1e-12."""
    case = parse_case(
        {
            "psds": {"P": {"type": "BASE", "gravity": 9.81, "rows": [[0.01, 0.0, 10.0], [0.04, 0.0, 200.0]]}},
            "base_motions": [
                {"dof": 1, "load_case": 1, "kind": "ACCELERATION"},
                {"dof": 2, "load_case": 1, "kind": "ACCELERATION"},
            ],
            "correlations": [{"psd": "P", "type": "CORRELATED", "weights": {"1": [1.0, 2.0]}}],
            "lower_frequency_hz": 10.0,
            "upper_frequency_hz": 200.0,
            "n_points": 40,
            "bias": 3.0,
            "frequency_scale": "LOGARITHMIC",
            "damping": {"mode_ranges": [[1, 1, 0.03], [2, 2, 0.05]]},
        }
    )
    w = np.array([2.0 * math.pi * 60.0, 2.0 * math.pi * 85.0])
    phi = np.array([[0.6, -0.3], [0.2, 0.8]])
    gam = np.array([[1.2, -0.5, 0.0], [0.7, 1.1, 0.0]])
    psi = np.zeros((2, 2, 6))
    psi[0, 0, 0], psi[0, 1, 1] = 2.0e6, 1.5e6
    psi[1, 0, 3], psi[1, 1, 3] = 1.0e6, -0.5e6

    ref = random_response_from_case(w, phi, gam, case, stress_modes=psi)
    new = random_vibration_from_case(w, phi, gam, case, stress_modes=psi)

    np.testing.assert_allclose(new.frequencies_hz, ref.frequencies_hz, rtol=1e-14)
    np.testing.assert_allclose(new.ru, ref.ru, rtol=1e-12)
    np.testing.assert_allclose(new.rv, ref.rv, rtol=1e-12)
    np.testing.assert_allclose(new.ra, ref.ra, rtol=1e-12)
    np.testing.assert_allclose(new.rmises, ref.rmises, rtol=1e-12)


def test_tiny_bspline_consistent_mass_modal_bridge_and_fatigue() -> None:
    """End-to-end tiny consistent-mass eigenmode -> ModalRandomVibration -> stress -> Dirlik fatigue."""
    # Tiny 3D B-spline box consistent mass + SPD spring stiffness (no Kratos required)
    ku = np.array([0.0, 0.0, 1.0, 1.0])
    kv = np.array([0.0, 0.0, 1.0, 1.0])
    kw = np.array([0.0, 0.0, 1.0, 1.0])
    M = assemble_consistent_mass_from_bspline((ku, kv, kw), degrees=1, rho=7850.0, bounds=((0, 1.0), (0, 0.1), (0, 0.1)))
    n_dofs = M.shape[0]
    free = np.arange(12, n_dofs)  # clamp x=0 face (4 control points * 3 DOFs)
    # Build a simple SPD stiffness matrix on free DOFs
    n_f = free.size
    diag = np.full(n_f, 1.0e6)
    off = np.full(n_f - 1, -2.0e5)
    K_ff = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    K = np.zeros((n_dofs, n_dofs), dtype=float)
    K[np.ix_(free, free)] = K_ff

    eig = solve_eigenmodes(K, M, n_modes=3, free_dofs=free)
    psi = np.zeros((1, 3, 6))
    psi[0, :, 0] = 5.0e7  # axial stress mode scaling

    solver = ModalRandomVibration.from_eigenmodes(eig, damping=0.02, stress_modes=psi)
    exc = BaseExcitation("y", WhiteNoisePSD(0.05), gravity_scale=9.81)
    res = solver.solve(exc, f_min=0.5 * eig.frequencies_hz[0], f_max=2.0 * eig.frequencies_hz[-1], n_points=120)

    assert res.ru.shape == (n_dofs,)
    assert np.all(res.ru[:12] == 0.0)
    assert np.max(res.ru[12:]) > 0.0
    assert res.rmises is not None and res.rmises[0] > 0.0

    fat = res.fatigue(SNCurve(C=1.0e45, b=5.0), method="dirlik", point_index=0, duration_seconds=100.0)
    assert fat.damage_rate > 0.0 and math.isfinite(fat.life_seconds)
    assert fat.damage_in_duration == pytest.approx(100.0 * fat.damage_rate, rel=1e-14)


def test_public_exports_in_dynamics_init() -> None:
    """All new random vibration symbols are cleanly exported from wht_iga.dynamics."""
    for name in (
        "PSDTable",
        "WhiteNoisePSD",
        "OctaveBandPSD",
        "ConstantModalDamping",
        "ModalDampingList",
        "RayleighDamping",
        "StructuralDamping",
        "BaseExcitation",
        "NodalForceExcitation",
        "CrossPSDMatrix",
        "ModalRandomVibration",
        "DirectFRFRandomVibration",
        "spectral_moments",
        "segalman_von_mises_rms",
        "NarrowBandRayleigh",
        "Steinberg3Band",
        "Dirlik",
        "SNCurve",
        "random_vibration_from_case",
    ):
        assert hasattr(dyn, name), f"Missing export {name} on wht_iga.dynamics"
        assert name in dyn.__all__


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
