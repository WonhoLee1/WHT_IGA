"""Modal superposition methods for transient, harmonic, random PSD, and shock spectrum analysis.

Implements Kratos-independent modal-coordinate solvers using mass-orthonormal
eigenmodes :math:`(\\omega_k, \\Phi_k)`:
- ``modal_transient``: Exact piecewise-linear Duhamel integration.
- ``modal_harmonic``: Steady-state complex frequency response.
- ``modal_random_psd``: Base-acceleration PSD integration with modal cross-terms (V09).
- ``modal_response_spectrum``: Shock response spectrum combination via ABS, SRSS, and CQC (V09).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
import scipy.integrate as integrate
import scipy.sparse as sp

CombinationRule = Literal["ABS", "SRSS", "CQC"]
PSDInput = float | Callable[[np.ndarray], np.ndarray | float]
ModalForceInput = np.ndarray | Callable[[float], np.ndarray]


@dataclass(frozen=True)
class ModalTransientResult:
    """Result container for modal transient Duhamel integration.

    Attributes:
        t: Time instants, shape ``(n_steps + 1,)``.
        u: Physical displacement history, shape ``(n_steps + 1, n_dofs)``.
        v: Physical velocity history, shape ``(n_steps + 1, n_dofs)``.
        q: Modal displacement coordinates, shape ``(n_steps + 1, n_modes)``.
        dq: Modal velocity coordinates, shape ``(n_steps + 1, n_modes)``.
    """

    t: np.ndarray
    u: np.ndarray
    v: np.ndarray
    q: np.ndarray
    dq: np.ndarray


@dataclass(frozen=True)
class ModalHarmonicResult:
    """Result container for steady-state modal harmonic frequency response.

    Attributes:
        frequencies_hz: Excitation frequencies in Hz, shape ``(n_freq,)``.
        omega: Excitation circular frequencies in rad/s, shape ``(n_freq,)``.
        u_complex: Complex displacement response, shape ``(n_freq, n_dofs)``.
        u_magnitude: Displacement amplitude ``|u_complex|``, shape ``(n_freq, n_dofs)``.
        u_phase_rad: Displacement phase in radians, shape ``(n_freq, n_dofs)``.
        modal_coordinates: Complex modal coordinates ``q(omega)``, shape ``(n_freq, n_modes)``.
    """

    frequencies_hz: np.ndarray
    omega: np.ndarray
    u_complex: np.ndarray
    u_magnitude: np.ndarray
    u_phase_rad: np.ndarray
    modal_coordinates: np.ndarray


@dataclass(frozen=True)
class ModalRandomPSDResult:
    """Result container for base-acceleration random vibration (PSD) analysis.

    Attributes:
        u_rms: Root-mean-square displacement per DOF, shape ``(n_dofs,)``.
        u_variance: Mean-square displacement ``u_rms**2`` per DOF, shape ``(n_dofs,)``.
        modal_cross_integral: Complex modal cross-spectral integral matrix
            ``I[j, k] = int H_j(w) H_k*(w) S_aa(f) df``, shape ``(n_modes, n_modes)``.
        frequencies_hz: Quadrature frequency grid used for integration.
    """

    u_rms: np.ndarray
    u_variance: np.ndarray
    modal_cross_integral: np.ndarray
    frequencies_hz: np.ndarray


@dataclass(frozen=True)
class ModalSpectrumResult:
    """Result container for response spectrum modal combination (ABS/SRSS/CQC).

    Attributes:
        u_combined: Combined peak displacement per DOF, shape ``(n_dofs,)``.
        modal_displacements: Signed modal displacement contributions
            ``Phi[:, k] * Gamma[k] * Sa(f_k) / omega_k**2``, shape ``(n_dofs, n_modes)``.
        sa_modes: Interpolated spectral acceleration at each mode, shape ``(n_modes,)``.
        correlation_matrix: Modal correlation matrix ``rho[j, k]``, shape ``(n_modes, n_modes)``.
        rule: Combination rule applied (``"ABS"``, ``"SRSS"``, or ``"CQC"``).
    """

    u_combined: np.ndarray
    modal_displacements: np.ndarray
    sa_modes: np.ndarray
    correlation_matrix: np.ndarray
    rule: str


def _normalize_damping(zeta: float | np.ndarray, n_modes: int) -> np.ndarray:
    """Broadcast scalar or 1D modal damping ratio to shape ``(n_modes,)``."""
    z = np.asarray(zeta, dtype=float)
    if z.ndim == 0:
        z_arr = np.full(n_modes, float(z), dtype=float)
    else:
        z_arr = z.reshape(-1)
        if z_arr.size != n_modes:
            raise ValueError(f"zeta length ({z_arr.size}) must equal n_modes ({n_modes}).")
    if np.any(z_arr < 0.0):
        raise ValueError("Modal damping ratios zeta must be non-negative.")
    return z_arr


def _select_gamma(participation_factors: np.ndarray, direction: int) -> np.ndarray:
    """Extract 1D modal participation vector for the requested excitation direction."""
    gamma = np.asarray(participation_factors, dtype=float)
    if gamma.ndim == 1:
        return gamma
    if gamma.ndim == 2:
        if not 0 <= direction < gamma.shape[1]:
            raise ValueError(f"direction must be in [0, {gamma.shape[1] - 1}], got {direction}.")
        return gamma[:, direction]
    raise ValueError("participation_factors must be 1D or 2D.")


def _duhamel_step_scalar(
    q0: float,
    dq0: float,
    p0: float,
    p1: float,
    omega: float,
    zeta: float,
    dt: float,
) -> tuple[float, float]:
    """Exact piecewise-linear Duhamel step for a single modal oscillator."""
    slope = (p1 - p0) / dt
    if omega <= 0.0:
        dq1 = dq0 + p0 * dt + 0.5 * slope * dt * dt
        q1 = q0 + dq0 * dt + 0.5 * p0 * dt * dt + (1.0 / 6.0) * slope * dt * dt * dt
        return q1, dq1

    w2 = omega * omega
    w3 = w2 * omega
    qp0 = p0 / w2 - 2.0 * zeta * slope / w3
    qp1 = p1 / w2 - 2.0 * zeta * slope / w3
    dqp = slope / w2

    y0 = q0 - qp0
    dy0 = dq0 - dqp

    if zeta < 1.0 - 1e-12:
        wd = omega * np.sqrt(1.0 - zeta * zeta)
        exp_t = np.exp(-zeta * omega * dt)
        cos_t = np.cos(wd * dt)
        sin_t = np.sin(wd * dt)
        a_coef = y0
        b_coef = (dy0 + zeta * omega * y0) / wd
        y1 = exp_t * (a_coef * cos_t + b_coef * sin_t)
        dy1 = exp_t * (
            (-zeta * omega * a_coef + wd * b_coef) * cos_t
            - (zeta * omega * b_coef + wd * a_coef) * sin_t
        )
    elif zeta <= 1.0 + 1e-12:
        exp_t = np.exp(-omega * dt)
        c1 = y0
        c2 = dy0 + omega * y0
        y1 = exp_t * (c1 + c2 * dt)
        dy1 = exp_t * (c2 - omega * (c1 + c2 * dt))
    else:
        sq = np.sqrt(zeta * zeta - 1.0)
        r1 = omega * (-zeta + sq)
        r2 = omega * (-zeta - sq)
        c1 = (dy0 - r2 * y0) / (r1 - r2)
        c2 = y0 - c1
        e1 = np.exp(r1 * dt)
        e2 = np.exp(r2 * dt)
        y1 = c1 * e1 + c2 * e2
        dy1 = c1 * r1 * e1 + c2 * r2 * e2

    return float(y1 + qp1), float(dy1 + dqp)


def modal_transient(
    omega: np.ndarray,
    modes: np.ndarray,
    f: ModalForceInput,
    dt: float,
    n_steps: int,
    zeta: float | np.ndarray = 0.0,
    u0: np.ndarray | None = None,
    v0: np.ndarray | None = None,
    M: sp.spmatrix | np.ndarray | None = None,
) -> ModalTransientResult:
    """Integrate modal equations of motion using exact piecewise-linear Duhamel.

    Args:
        omega: Modal circular frequencies in rad/s, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes ``Phi``, shape ``(n_dofs, n_modes)``.
        f: Constant force vector ``(n_dofs,)``, callable ``f(t) -> (n_dofs,)``,
            or pre-sampled array ``(n_steps + 1, n_dofs)``.
        dt: Time step size in seconds.
        n_steps: Number of time steps.
        zeta: Scalar or per-mode viscous damping ratio.
        u0: Optional initial physical displacement, shape ``(n_dofs,)``.
        v0: Optional initial physical velocity, shape ``(n_dofs,)``.
        M: Mass matrix required to project non-zero ``u0`` or ``v0`` onto
            modal coordinates via ``q0 = Phi.T @ M @ u0``.

    Returns:
        ``ModalTransientResult`` with time vector, physical displacements and
        velocities, and modal coordinates.
    """
    if dt <= 0.0:
        raise ValueError(f"dt must be positive, got {dt}")
    if n_steps < 1:
        raise ValueError(f"n_steps must be >= 1, got {n_steps}")

    omega_arr = np.asarray(omega, dtype=float).reshape(-1)
    phi = np.asarray(modes, dtype=float)
    n_dofs, n_modes = phi.shape
    if omega_arr.size != n_modes:
        raise ValueError(f"omega length ({omega_arr.size}) must match modes.shape[1] ({n_modes}).")
    zeta_arr = _normalize_damping(zeta, n_modes)

    t = np.linspace(0.0, n_steps * dt, n_steps + 1)

    # Evaluate modal force time history p[step, mode] = Phi.T @ f(t_step)
    p_hist = np.empty((n_steps + 1, n_modes), dtype=float)
    if callable(f):
        for s in range(n_steps + 1):
            f_vec = np.asarray(f(float(t[s])), dtype=float).reshape(n_dofs)
            p_hist[s, :] = phi.T @ f_vec
    else:
        f_arr = np.asarray(f, dtype=float)
        if f_arr.ndim == 1:
            p_const = phi.T @ f_arr.reshape(n_dofs)
            p_hist[:, :] = p_const[None, :]
        elif f_arr.shape == (n_steps + 1, n_dofs):
            p_hist[:, :] = f_arr @ phi
        else:
            raise ValueError(f"Force array shape {f_arr.shape} incompatible with ({n_steps + 1}, {n_dofs}).")

    q = np.zeros((n_steps + 1, n_modes), dtype=float)
    dq = np.zeros((n_steps + 1, n_modes), dtype=float)

    if u0 is not None or v0 is not None:
        if M is None:
            raise ValueError("Mass matrix M is required to project initial conditions u0/v0.")
        M_csr = sp.csr_matrix(M, dtype=float)
        if u0 is not None:
            q[0, :] = phi.T @ (M_csr @ np.asarray(u0, dtype=float).reshape(n_dofs))
        if v0 is not None:
            dq[0, :] = phi.T @ (M_csr @ np.asarray(v0, dtype=float).reshape(n_dofs))

    for k in range(n_modes):
        wk = float(omega_arr[k])
        zk = float(zeta_arr[k])
        qk = float(q[0, k])
        dqk = float(dq[0, k])
        for s in range(n_steps):
            qk, dqk = _duhamel_step_scalar(qk, dqk, float(p_hist[s, k]), float(p_hist[s + 1, k]), wk, zk, dt)
            q[s + 1, k] = qk
            dq[s + 1, k] = dqk

    u = q @ phi.T
    v = dq @ phi.T
    return ModalTransientResult(t=t, u=u, v=v, q=q, dq=dq)


def modal_harmonic(
    omega_modes: np.ndarray,
    modes: np.ndarray,
    force_amplitude: np.ndarray,
    frequencies_hz: np.ndarray,
    zeta: float | np.ndarray = 0.0,
) -> ModalHarmonicResult:
    """Compute steady-state complex frequency response via modal superposition.

    Evaluates ``H_k(w) = 1 / (w_k^2 - w^2 + 2 i zeta_k w_k w)`` for each
    excitation frequency ``w = 2 * pi * f``.

    Args:
        omega_modes: Modal circular frequencies in rad/s, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes ``Phi``, shape ``(n_dofs, n_modes)``.
        force_amplitude: Complex or real spatial force amplitude ``f0``, shape ``(n_dofs,)``.
        frequencies_hz: Excitation frequencies in Hz, shape ``(n_freq,)``.
        zeta: Scalar or per-mode viscous damping ratio.

    Returns:
        ``ModalHarmonicResult`` containing complex responses, magnitudes, and phases.
    """
    w_modes = np.asarray(omega_modes, dtype=float).reshape(-1)
    phi = np.asarray(modes, dtype=float)
    n_dofs, n_modes = phi.shape
    if w_modes.size != n_modes:
        raise ValueError(f"omega_modes length ({w_modes.size}) must match modes.shape[1] ({n_modes}).")

    zeta_arr = _normalize_damping(zeta, n_modes)
    f0 = np.asarray(force_amplitude).reshape(n_dofs)
    freqs = np.asarray(frequencies_hz, dtype=float).reshape(-1)
    w_exc = 2.0 * np.pi * freqs

    # Modal force amplitudes: (n_modes,)
    modal_f = phi.T @ f0
    # Transfer function H[freq, mode]
    denom = (
        w_modes[None, :] ** 2
        - w_exc[:, None] ** 2
        + 2j * zeta_arr[None, :] * w_modes[None, :] * w_exc[:, None]
    )
    h_mat = 1.0 / denom
    q_complex = h_mat * modal_f[None, :]
    u_complex = q_complex @ phi.T

    return ModalHarmonicResult(
        frequencies_hz=freqs,
        omega=w_exc,
        u_complex=u_complex,
        u_magnitude=np.abs(u_complex),
        u_phase_rad=np.angle(u_complex),
        modal_coordinates=q_complex,
    )


def _build_psd_frequency_grid(
    f_min: float,
    f_max: float,
    f_modes: np.ndarray,
    zeta: np.ndarray,
    n_freq: int,
) -> np.ndarray:
    """Build a frequency grid combining base sampling with resonance refinement."""
    grids = [np.linspace(f_min, f_max, max(int(n_freq), 200))]
    for fk, zk in zip(f_modes, zeta, strict=True):
        if f_min <= fk <= f_max and fk > 0.0:
            bw = max(zk, 1e-4) * fk
            lo = max(f_min, fk - 15.0 * bw)
            hi = min(f_max, fk + 15.0 * bw)
            if hi > lo:
                grids.append(np.linspace(lo, hi, 401))
            lo_tight = max(f_min, fk - 2.5 * bw)
            hi_tight = min(f_max, fk + 2.5 * bw)
            if hi_tight > lo_tight:
                grids.append(np.linspace(lo_tight, hi_tight, 401))
    return np.unique(np.concatenate(grids))


def modal_random_psd(
    omega_modes: np.ndarray,
    modes: np.ndarray,
    participation_factors: np.ndarray,
    psd_func: PSDInput,
    f_min: float,
    f_max: float,
    zeta: float | np.ndarray = 0.02,
    direction: int = 0,
    n_freq: int = 2000,
) -> ModalRandomPSDResult:
    """Compute exact RMS displacement under base acceleration PSD (V09).

    Integrates the full modal cross-spectral density over ``[f_min, f_max]``
    including all cross-modal terms ``Phi_j Phi_k Gamma_j Gamma_k H_j(w) H_k*(w)``.

    Args:
        omega_modes: Modal circular frequencies in rad/s, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes ``Phi``, shape ``(n_dofs, n_modes)``.
        participation_factors: Modal participation factors, shape ``(n_modes,)``
            or ``(n_modes, 3)``.
        psd_func: One-sided base acceleration PSD ``S_aa(f)`` (in ``(accel)^2 / Hz``),
            either a constant float or callable ``psd_func(f_array)``.
        f_min: Lower integration frequency bound in Hz (must be >= 0).
        f_max: Upper integration frequency bound in Hz (must be > f_min).
        zeta: Scalar or per-mode viscous damping ratio (must be > 0 at resonances).
        direction: Excitation axis ``0, 1, 2`` if ``participation_factors`` is 2D.
        n_freq: Base number of quadrature frequency points before resonance refinement.

    Returns:
        ``ModalRandomPSDResult`` with per-DOF RMS and variance displacements.
    """
    if f_min < 0.0 or f_max <= f_min:
        raise ValueError(f"Frequency bounds must satisfy 0 <= f_min < f_max, got [{f_min}, {f_max}].")

    w_modes = np.asarray(omega_modes, dtype=float).reshape(-1)
    phi = np.asarray(modes, dtype=float)
    n_dofs, n_modes = phi.shape
    if w_modes.size != n_modes:
        raise ValueError(f"omega_modes length ({w_modes.size}) must match modes.shape[1] ({n_modes}).")

    gamma = _select_gamma(participation_factors, direction)
    if gamma.size != n_modes:
        raise ValueError(f"participation_factors length ({gamma.size}) must match n_modes ({n_modes}).")

    zeta_arr = _normalize_damping(zeta, n_modes)
    f_modes = w_modes / (2.0 * np.pi)

    freqs = _build_psd_frequency_grid(float(f_min), float(f_max), f_modes, zeta_arr, n_freq)
    w_grid = 2.0 * np.pi * freqs

    if callable(psd_func):
        s_aa = np.asarray(psd_func(freqs), dtype=float).reshape(-1)
        if s_aa.size == 1:
            s_aa = np.full(freqs.shape, float(s_aa[0]), dtype=float)
    else:
        s_aa = np.full(freqs.shape, float(psd_func), dtype=float)

    denom = (
        w_modes[None, :] ** 2
        - w_grid[:, None] ** 2
        + 2j * zeta_arr[None, :] * w_modes[None, :] * w_grid[:, None]
    )
    h_mat = 1.0 / denom  # (n_grid, n_modes)

    # Cross-spectral integral I[j, k] = int H_j(w) * conj(H_k(w)) * S_aa(f) df
    integrand = h_mat[:, :, None] * np.conj(h_mat[:, None, :]) * s_aa[:, None, None]
    cross_integral = integrate.trapezoid(integrand, freqs, axis=0)

    weighted_phi = phi * gamma[None, :]
    u_var = np.einsum("ij,jk,ik->i", weighted_phi, np.real(cross_integral), weighted_phi)
    u_var = np.maximum(u_var, 0.0)
    u_rms = np.sqrt(u_var)

    return ModalRandomPSDResult(
        u_rms=u_rms,
        u_variance=u_var,
        modal_cross_integral=cross_integral,
        frequencies_hz=freqs,
    )


def modal_response_spectrum(
    omega_modes: np.ndarray,
    modes: np.ndarray,
    participation_factors: np.ndarray,
    spectrum_freqs_hz: np.ndarray,
    spectrum_sa: np.ndarray,
    zeta: float | np.ndarray = 0.05,
    rule: CombinationRule = "CQC",
    direction: int = 0,
) -> ModalSpectrumResult:
    """Compute peak displacement from a shock/response spectrum (V09).

    Supports ``ABS``, ``SRSS``, and Der Kiureghian ``CQC`` modal combination
    rules using linearly interpolated spectral accelerations ``S_a(f_k)``.

    Args:
        omega_modes: Modal circular frequencies in rad/s, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes ``Phi``, shape ``(n_dofs, n_modes)``.
        participation_factors: Modal participation factors, shape ``(n_modes,)``
            or ``(n_modes, 3)``.
        spectrum_freqs_hz: Sorted frequencies of the shock spectrum in Hz.
        spectrum_sa: Spectral acceleration values ``S_a(f)`` matching
            ``spectrum_freqs_hz``.
        zeta: Scalar or per-mode viscous damping ratio for CQC correlation.
        rule: Modal combination rule: ``"ABS"``, ``"SRSS"``, or ``"CQC"``.
        direction: Excitation axis ``0, 1, 2`` if ``participation_factors`` is 2D.

    Returns:
        ``ModalSpectrumResult`` with combined peak displacements and modal
        contributions.
    """
    rule_upper = rule.upper()
    if rule_upper not in ("ABS", "SRSS", "CQC"):
        raise ValueError(f"Unsupported combination rule '{rule}'. Expected 'ABS', 'SRSS', or 'CQC'.")

    w_modes = np.asarray(omega_modes, dtype=float).reshape(-1)
    phi = np.asarray(modes, dtype=float)
    n_dofs, n_modes = phi.shape
    if w_modes.size != n_modes:
        raise ValueError(f"omega_modes length ({w_modes.size}) must match modes.shape[1] ({n_modes}).")
    if np.any(w_modes <= 0.0):
        raise ValueError("All modal circular frequencies must be positive for response spectrum analysis.")

    gamma = _select_gamma(participation_factors, direction)
    if gamma.size != n_modes:
        raise ValueError(f"participation_factors length ({gamma.size}) must match n_modes ({n_modes}).")

    zeta_arr = _normalize_damping(zeta, n_modes)
    spec_f = np.asarray(spectrum_freqs_hz, dtype=float).reshape(-1)
    spec_sa = np.asarray(spectrum_sa, dtype=float).reshape(-1)
    if spec_f.size != spec_sa.size or spec_f.size < 1:
        raise ValueError("spectrum_freqs_hz and spectrum_sa must have identical non-zero length.")

    f_modes = w_modes / (2.0 * np.pi)
    sa_modes = np.interp(f_modes, spec_f, spec_sa)

    # Signed modal displacement vectors: u_mode[:, k] = Phi[:, k] * Gamma[k] * Sa(f_k) / w_k^2
    modal_scale = gamma * sa_modes / (w_modes**2)
    u_modes = phi * modal_scale[None, :]

    corr = np.eye(n_modes, dtype=float)
    if rule_upper == "ABS":
        u_combined = np.sum(np.abs(u_modes), axis=1)
    elif rule_upper == "SRSS":
        u_combined = np.sqrt(np.sum(u_modes**2, axis=1))
    else:
        # Der Kiureghian (1981) CQC correlation matrix
        for j in range(n_modes):
            for k in range(j + 1, n_modes):
                r = float(w_modes[k] / w_modes[j])
                zj = float(zeta_arr[j])
                zk = float(zeta_arr[k])
                num = 8.0 * np.sqrt(zj * zk) * (zj + r * zk) * (r**1.5)
                den = (1.0 - r * r) ** 2 + 4.0 * zj * zk * r * (1.0 + r * r) + 4.0 * (zj * zj + zk * zk) * r * r
                rho_jk = float(num / den) if den > 0.0 else 0.0
                corr[j, k] = rho_jk
                corr[k, j] = rho_jk
        u_sq = np.einsum("ij,jk,ik->i", u_modes, corr, u_modes)
        u_combined = np.sqrt(np.maximum(u_sq, 0.0))

    return ModalSpectrumResult(
        u_combined=u_combined,
        modal_displacements=u_modes,
        sa_modes=sa_modes,
        correlation_matrix=corr,
        rule=rule_upper,
    )
