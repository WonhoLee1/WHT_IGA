"""Unit tests for Phase 4 Kratos-independent dynamics and modal analysis engine."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as la
import scipy.sparse as sp

from wht_iga.dynamics import (
    BOSSAK_PRESETS,
    MODERATE_DISSIPATION_ALPHA_B,
    TRANSIENT_FIDELITY_ALPHA_B,
    BossakConfig,
    assemble_consistent_mass,
    assemble_consistent_mass_from_bspline,
    modal_harmonic,
    modal_random_psd,
    modal_response_spectrum,
    modal_transient,
    solve_bossak_linear,
    solve_eigenmodes,
)
from wht_iga.grid.knots import uniform_clamped
from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline


def test_consistent_mass_conservation_symmetry_and_pd() -> None:
    """Verify total mass conservation (<1e-13), symmetry, and positive definiteness."""
    lx, ly, lz = 0.2, 0.05, 0.04
    rho = 7850.0
    knots = (
        uniform_clamped(0.0, lx, 3, degree=2),
        uniform_clamped(0.0, ly, 2, degree=2),
        uniform_clamped(0.0, lz, 2, degree=2),
    )
    M = assemble_consistent_mass_from_bspline(knots, degrees=(2, 2, 2), rho=rho)

    # Symmetry: ||M - M^T||_max == 0
    sym_err = np.max(np.abs((M - M.T).toarray()))
    assert sym_err < 1e-15

    # Positive definiteness: smallest eigenvalue > 0
    eigvals = la.eigvalsh(M.toarray())
    assert np.min(eigvals) > 0.0

    # Total mass conservation in x, y, z: r_d^T M r_d == rho * V
    expected_mass = rho * lx * ly * lz
    n_dofs = M.shape[0]
    for d in range(3):
        r_d = np.zeros(n_dofs, dtype=float)
        r_d[d::3] = 1.0
        mass_d = float(r_d @ (M @ r_d))
        assert abs(mass_d - expected_mass) / expected_mass < 1e-13


def test_consistent_mass_normalized_knots_and_callable_rho() -> None:
    """Verify normalized [0,1]^3 knots with physical bounds and spatial rho(x)."""
    lx, ly, lz = 2.0, 1.0, 0.5
    knots = (
        uniform_clamped(0.0, 1.0, 2, degree=1),
        uniform_clamped(0.0, 1.0, 2, degree=1),
        uniform_clamped(0.0, 1.0, 1, degree=1),
    )
    bounds = ((0.0, lx), (0.0, ly), (0.0, lz))
    # Linear density rho(x) = 1000 * (1 + x) => integral over [0, 2] is 1000*(2 + 2) = 4000
    M = assemble_consistent_mass_from_bspline(
        knots,
        degrees=1,
        rho=lambda pts: 1000.0 * (1.0 + pts[:, 0]),
        bounds=bounds,
    )
    expected_mass = 1000.0 * (lx + 0.5 * lx * lx) * ly * lz
    r_x = np.zeros(M.shape[0], dtype=float)
    r_x[0::3] = 1.0
    assert abs(float(r_x @ (M @ r_x)) - expected_mass) / expected_mass < 1e-13


def test_assemble_consistent_mass_quadrature_formats() -> None:
    """Verify general quadrature assembly matches B-spline assembly for 1-cell linear hex."""
    # Single linear element on [0, 1]^3 with 2x2x2 Gauss-Legendre points
    xg, wg = np.polynomial.legendre.leggauss(2)
    pts_1d = 0.5 * (xg + 1.0)
    w_1d = 0.5 * wg
    gx, gy, gz = np.meshgrid(pts_1d, pts_1d, pts_1d, indexing="ij")
    pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    w_detj = (w_1d[:, None, None] * w_1d[None, :, None] * w_1d[None, None, :]).ravel()

    # Linear Lagrange/B-spline basis on [0, 1]^3
    n_vals = np.empty((8, 8), dtype=float)
    for q in range(8):
        x, y, z = pts[q]
        bx = [1.0 - x, x]
        by = [1.0 - y, y]
        bz = [1.0 - z, z]
        for c in range(2):
            for b in range(2):
                for a in range(2):
                    loc = a + 2 * (b + 2 * c)
                    n_vals[q, loc] = bx[a] * by[b] * bz[c]

    cp_dofs = np.arange(8, dtype=np.int64)[None, :]
    dofs_3d = np.arange(24, dtype=np.int64)[None, :]
    rho = 2700.0

    M_cp = assemble_consistent_mass(pts[None, ...], w_detj[None, ...], n_vals[None, ...], cp_dofs, rho=rho)
    M_3d = assemble_consistent_mass(pts[None, ...], w_detj[None, ...], n_vals[None, ...], dofs_3d, rho=rho)

    knots = (
        uniform_clamped(0.0, 1.0, 1, degree=1),
        uniform_clamped(0.0, 1.0, 1, degree=1),
        uniform_clamped(0.0, 1.0, 1, degree=1),
    )
    M_ref = assemble_consistent_mass_from_bspline(knots, degrees=1, rho=rho)

    assert np.max(np.abs((M_cp - M_ref).toarray())) < 1e-13
    assert np.max(np.abs((M_3d - M_ref).toarray())) < 1e-13


def test_eigenmodes_axial_bar_v10_and_orthonormality() -> None:
    """Verify V10 axial bar first mode (25000 Hz), M-orthonormality (<1e-11), and modal mass."""
    # V10 parameters: E = 2e11 Pa, rho = 8000 kg/m^3 => c = 5000 m/s, L = 0.05 m => f1 = 25000 Hz
    L = 0.05
    E = 2.0e11
    rho = 8000.0
    box = BoxSpline(
        lo=(0.0, 0.0, 0.0),
        hi=(L, 0.01, 0.01),
        n_cells=(8, 1, 1),
        degree=2,
    )
    K = assemble_stiffness(box, E=E, nu=0.0)
    M = assemble_consistent_mass_from_bspline(box.knots, degrees=box.degree, rho=rho)

    # Pure 1D axial mode test: fix x=0 (i_u == 0) and constrain lateral DOFs (u_y, u_z)
    free_axial: list[int] = []
    nu_cp, nv_cp, nw_cp = box.n_cp
    for iw in range(nw_cp):
        for iv in range(nv_cp):
            for iu in range(1, nu_cp):
                cp_idx = iu + nu_cp * (iv + nv_cp * iw)
                free_axial.append(3 * cp_idx)

    res_axial = solve_eigenmodes(K, M, n_modes=3, sigma=0.0, free_dofs=free_axial)
    f1_analytic = 0.25 / L * np.sqrt(E / rho)  # 25000.0 Hz
    assert abs(f1_analytic - 25000.0) < 1e-9
    rel_err = abs(res_axial.frequencies_hz[0] - f1_analytic) / f1_analytic
    assert rel_err < 1e-4

    # Full 3D cantilever eigenmodes (fix all 3 DOFs on x=0 face)
    free_3d: list[int] = []
    for iw in range(nw_cp):
        for iv in range(nv_cp):
            for iu in range(1, nu_cp):
                cp_idx = iu + nu_cp * (iv + nv_cp * iw)
                free_3d.extend([3 * cp_idx, 3 * cp_idx + 1, 3 * cp_idx + 2])

    res_3d = solve_eigenmodes(K, M, n_modes=10, sigma=0.0, free_dofs=free_3d)

    # Mass-orthonormality: Phi^T M Phi = I to < 1e-11
    gram = res_3d.modes.T @ (M @ res_3d.modes)
    ortho_err = np.max(np.abs(gram - np.eye(10)))
    assert ortho_err < 1e-11

    # Sum of effective modal masses <= total mass in each direction
    sum_eff = np.sum(res_3d.effective_modal_masses, axis=0)
    assert np.all(sum_eff <= res_3d.total_mass + 1e-12)


def test_bossak_sdof_accuracy_and_dissipation() -> None:
    """Verify Bossak-Newmark free vibration accuracy, energy conservation, and HF dissipation."""
    # Undamped SDOF: m = 1, omega_n = 2*pi*10 rad/s (10 Hz), u0 = 1.0, v0 = 0.0
    m = 1.0
    omega_n = 2.0 * np.pi * 10.0
    k = m * omega_n**2
    M = sp.csr_matrix([[m]])
    K = sp.csr_matrix([[k]])
    f0 = np.zeros(1)

    dt = 0.0005  # 200 steps per period
    n_steps = 200
    cfg_cons = BossakConfig(dt=dt, n_steps=n_steps, alpha_b=0.0)
    res_cons = solve_bossak_linear(M, None, K, f0, cfg_cons, u0=np.array([1.0]), v0=np.array([0.0]))

    u_exact = np.cos(omega_n * res_cons.t)
    assert np.max(np.abs(res_cons.u[:, 0] - u_exact)) < 5e-3
    # Trapezoidal rule (alpha_b = 0) conserves total energy
    e0 = res_cons.energy_total[0]
    assert np.max(np.abs(res_cons.energy_total - e0)) / e0 < 1e-11

    # V10 presets exist with expected alpha_b
    assert BOSSAK_PRESETS["TRANSIENT_FIDELITY"] == TRANSIENT_FIDELITY_ALPHA_B == -0.05
    assert BOSSAK_PRESETS["MODERATE_DISSIPATION"] == MODERATE_DISSIPATION_ALPHA_B == -0.30

    # 2-DOF high-frequency dissipation check: mode 2 has omega_2 * dt >> 1
    M2 = sp.diags([1.0, 1.0], format="csr")
    K2 = sp.diags([10.0**2, 10000.0**2], format="csr")
    u0_2 = np.array([0.0, 1.0])  # excite only high-frequency mode
    dt_hf = 0.01  # omega_2 * dt = 100

    res_0 = solve_bossak_linear(M2, None, K2, np.zeros(2), BossakConfig(dt=dt_hf, n_steps=20, alpha_b=0.0), u0=u0_2)
    res_fid = solve_bossak_linear(M2, None, K2, np.zeros(2), BossakConfig.transient_fidelity(dt_hf, 20), u0=u0_2)
    res_mod = solve_bossak_linear(M2, None, K2, np.zeros(2), BossakConfig.moderate_dissipation(dt_hf, 20), u0=u0_2)

    assert res_mod.energy_total[-1] < res_fid.energy_total[-1] < res_0.energy_total[-1]


def test_modal_transient_and_harmonic_sdof() -> None:
    """Verify exact Duhamel transient and SDOF harmonic resonance peak 1/(2*zeta)."""
    m = 2.0
    omega_n = 2.0 * np.pi * 25.0
    k = m * omega_n**2
    phi = np.array([[1.0 / np.sqrt(m)]])
    omega = np.array([omega_n])
    M = sp.csr_matrix([[m]])

    # 1) Undamped transient under constant step force F0: u(t) = (F0/k)*(1 - cos(omega_n t))
    f_step = np.array([10.0])
    dt = 0.001
    n_steps = 100
    res_tr = modal_transient(omega, phi, f_step, dt=dt, n_steps=n_steps, zeta=0.0, M=M)
    u_exact = (10.0 / k) * (1.0 - np.cos(omega_n * res_tr.t))
    assert np.max(np.abs(res_tr.u[:, 0] - u_exact)) < 1e-12

    # 2) Harmonic resonance peak: |u(omega_n)| / (F0 / k) == 1 / (2 * zeta)
    zeta = 0.02
    f0 = 5.0
    res_harm = modal_harmonic(omega, phi, np.array([f0]), frequencies_hz=np.array([25.0]), zeta=zeta)
    u_static = f0 / k
    amplification = float(res_harm.u_magnitude[0, 0]) / u_static
    assert abs(amplification - 1.0 / (2.0 * zeta)) < 1e-12


def test_modal_random_psd_miles_formula() -> None:
    """Verify SDOF white-noise base-acceleration PSD RMS against exact Miles' formula (<1%)."""
    m = 3.5
    fn = 100.0
    omega_n = 2.0 * np.pi * fn
    zeta = 0.02
    s0 = 0.5  # (m/s^2)^2 / Hz

    phi = np.array([[1.0 / np.sqrt(m)]])
    gamma = np.array([np.sqrt(m)])  # Phi * Gamma == 1.0
    omega = np.array([omega_n])

    # Integrate over [0, 20 * fn] to capture >99.9% of infinite-band integral
    res_psd = modal_random_psd(
        omega_modes=omega,
        modes=phi,
        participation_factors=gamma,
        psd_func=s0,
        f_min=0.0,
        f_max=20.0 * fn,
        zeta=zeta,
    )

    # Exact integral over [0, inf): int |H(f)|^2 df = 1 / (8 * zeta * omega_n^3)
    u_rms_miles = np.sqrt(s0 / (8.0 * zeta * omega_n**3))
    rel_err = abs(float(res_psd.u_rms[0]) - u_rms_miles) / u_rms_miles
    assert rel_err < 0.01


def test_modal_response_spectrum_rules() -> None:
    """Verify ABS, SRSS, and CQC ordering (SRSS <= ABS) and limiting correlation behavior."""
    # Two modes with same-sign participation at a single DOF
    phi = np.array([[0.6, 0.8]])
    gamma = np.array([1.2, 0.9])
    spec_f = np.array([1.0, 500.0])
    spec_sa = np.array([10.0, 10.0])

    # Closely spaced modes -> CQC approaches ABS; SRSS < CQC <= ABS
    omega_close = 2.0 * np.pi * np.array([100.0, 101.0])
    r_abs = modal_response_spectrum(omega_close, phi, gamma, spec_f, spec_sa, zeta=0.05, rule="ABS")
    r_srss = modal_response_spectrum(omega_close, phi, gamma, spec_f, spec_sa, zeta=0.05, rule="SRSS")
    r_cqc = modal_response_spectrum(omega_close, phi, gamma, spec_f, spec_sa, zeta=0.05, rule="CQC")

    assert r_srss.u_combined[0] <= r_cqc.u_combined[0] <= r_abs.u_combined[0] + 1e-15
    assert r_cqc.correlation_matrix[0, 1] > 0.95

    # Well-separated modes -> CQC approaches SRSS
    omega_far = 2.0 * np.pi * np.array([20.0, 200.0])
    r_srss_far = modal_response_spectrum(omega_far, phi, gamma, spec_f, spec_sa, zeta=0.05, rule="SRSS")
    r_cqc_far = modal_response_spectrum(omega_far, phi, gamma, spec_f, spec_sa, zeta=0.05, rule="CQC")
    assert r_srss_far.u_combined[0] <= r_cqc_far.u_combined[0]
    assert abs(r_cqc_far.u_combined[0] - r_srss_far.u_combined[0]) / r_srss_far.u_combined[0] < 1e-3
