"""Modal dynamics analysis facade and Bossak amplification theory (P4-B, V09/V10).

* :class:`ModalDynamicsAnalysis` bundles an eigen-solution with optional stress mode
  shapes and runs modal transient / harmonic responses, case-file response spectra
  and case-file random responses (all Kratos-free once the modes exist).
* :func:`bossak_free_vibration_theory` gives the exact frequency and decay rate a
  Bossak-Newmark step produces for an undamped single-degree-of-freedom oscillator;
  V10 uses it as an analytic reference for the numerical dissipation (the Abaqus HHT
  histories of the Coreform study are not available here).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from wht_iga.dynamics.case_file import RandomResponseCase, ResponseSpectrumCase
from wht_iga.dynamics.consistent_mass import EigenmodeResult
from wht_iga.dynamics.modal_superposition import (
    ModalHarmonicResult,
    ModalTransientResult,
    modal_harmonic,
    modal_transient,
)
from wht_iga.dynamics.response_spectrum import (
    RandomResponseResult,
    ResponseSpectrumResult,
    random_response_from_case,
    response_spectrum_from_case,
)


def bossak_amplification_matrix(alpha_b: float, omega: float, dt: float) -> np.ndarray:
    """One-step amplification matrix of the Bossak-Newmark scheme for ``a + w**2 u = 0``.

    State ``(u, v, a)``; ``gamma = 1/2 - alpha``, ``beta = (1 - alpha)**2 / 4`` and the
    equilibrium ``(1 - alpha) a_{n+1} + alpha a_n + w**2 u_{n+1} = 0``.

    Args:
        alpha_b: Bossak parameter (``<= 0``).
        omega: Oscillator circular frequency.
        dt: Time step.

    Returns:
        ``(3, 3)`` matrix ``A`` with ``state_{n+1} = A state_n``.
    """
    gamma = 0.5 - alpha_b
    beta = 0.25 * (1.0 - alpha_b) ** 2
    denom = (1.0 - alpha_b) + beta * omega**2 * dt**2
    amp = np.zeros((3, 3))
    for j in range(3):
        u0, v0, a0 = np.eye(3)[j]
        a1 = (-alpha_b * a0 - omega**2 * (u0 + dt * v0 + dt * dt * (0.5 - beta) * a0)) / denom
        amp[:, j] = [
            u0 + dt * v0 + dt * dt * ((0.5 - beta) * a0 + beta * a1),
            v0 + dt * ((1.0 - gamma) * a0 + gamma * a1),
            a1,
        ]
    return amp


def generalized_alpha_amplification_matrix(alpha_m: float, alpha_f: float, omega: float,
                                           dt: float) -> np.ndarray:
    """One-step amplification matrix of generalized-alpha for ``a + w**2 u = 0``.

    Kratos coefficient convention (``TimeDiscretization::GeneralizedAlpha``):
    ``beta = (1 - am + af)**2 / 4``, ``gamma = 1/2 - am + af``; equilibrium
    ``(1-am) a_{n+1} + am a_n + w**2 ((1-af) u_{n+1} + af u_n) = 0``. Bossak is
    ``af = 0``; HHT is ``am = 0``, ``af = -alpha_HHT``. State ``(u, v, a)``.
    """
    gamma = 0.5 - alpha_m + alpha_f
    beta = 0.25 * (1.0 - alpha_m + alpha_f) ** 2
    w2 = omega**2
    denom = (1.0 - alpha_m) + (1.0 - alpha_f) * beta * w2 * dt**2
    amp = np.zeros((3, 3))
    for j in range(3):
        u0, v0, a0 = np.eye(3)[j]
        u_pred = u0 + dt * v0 + dt * dt * (0.5 - beta) * a0
        a1 = (-alpha_m * a0 - w2 * ((1.0 - alpha_f) * u_pred + alpha_f * u0)) / denom
        amp[:, j] = [u_pred + dt * dt * beta * a1,
                     v0 + dt * ((1.0 - gamma) * a0 + gamma * a1),
                     a1]
    return amp


def generalized_alpha_free_vibration_theory(alpha_m: float, alpha_f: float, omega: float,
                                            dt: float) -> tuple[float, float]:
    """Frequency (Hz) and envelope decay rate (1/s) of free vibration under
    generalized-alpha stepping (see :func:`generalized_alpha_amplification_matrix`)."""
    lam = np.linalg.eigvals(generalized_alpha_amplification_matrix(alpha_m, alpha_f, omega, dt))
    lam = lam[np.argmax(np.abs(lam.imag))]
    return float(abs(np.angle(lam)) / (2.0 * np.pi * dt)), float(-np.log(abs(lam)) / dt)


def bossak_free_vibration_theory(alpha_b: float, omega: float, dt: float) -> tuple[float, float]:
    """Frequency (Hz) and envelope decay rate (1/s) of free vibration under Bossak stepping.

    From the complex eigenvalue pair ``rho exp(+-i theta)`` of the amplification matrix:
    ``f = theta / (2 pi dt)`` and ``s = -ln(rho) / dt``.

    Args:
        alpha_b: Bossak parameter.
        omega: Exact oscillator frequency in rad/s.
        dt: Time step.

    Returns:
        ``(frequency_hz, decay_rate)``.
    """
    lam = np.linalg.eigvals(bossak_amplification_matrix(alpha_b, omega, dt))
    lam = lam[np.argmax(np.abs(lam.imag))]
    return float(abs(np.angle(lam)) / (2.0 * np.pi * dt)), float(-np.log(abs(lam)) / dt)


def modal_coordinate(mode: np.ndarray, M, history: np.ndarray) -> np.ndarray:
    """Modal coordinate history ``q(t) = phi^T M u(t)`` of one mass-normalised mode.

    Args:
        mode: Mode shape ``(n_dofs,)``.
        M: Mass matrix ``(n_dofs, n_dofs)``.
        history: Displacement history ``(n_steps, n_dofs)``.

    Returns:
        ``(n_steps,)`` array.
    """
    return np.asarray(history, dtype=float) @ (M @ np.asarray(mode, dtype=float))


@dataclass
class ModalDynamicsAnalysis:
    """Modal-basis analyses on one eigen-solution.

    Attributes:
        eig: Mass-orthonormal eigenmodes.
        stress_modes: Optional stress mode shapes ``(n_pts, n_modes, 6)`` for RMISES / stress spectra.
        M: Optional mass matrix (needed to project initial conditions in transient runs).
    """

    eig: EigenmodeResult
    stress_modes: np.ndarray | None = None
    M: object | None = None

    def transient(self, force, dt: float, n_steps: int, zeta: float | np.ndarray = 0.0) -> ModalTransientResult:
        """Exact piecewise-linear Duhamel response to ``force`` (vector, callable or sampled)."""
        return modal_transient(self.eig.omega, self.eig.modes, force, dt, n_steps, zeta, M=self.M)

    def harmonic(self, force_amplitude: np.ndarray, frequencies_hz: np.ndarray, zeta: float | np.ndarray = 0.0) -> ModalHarmonicResult:
        """Steady-state complex frequency response."""
        return modal_harmonic(self.eig.omega, self.eig.modes, force_amplitude, frequencies_hz, zeta)

    def response_spectrum(self, case: ResponseSpectrumCase, quantity: str = "U") -> ResponseSpectrumResult:
        """Case-file response spectrum for ``quantity`` in ``U, V, A, S``."""
        if quantity.upper() == "S":
            if self.stress_modes is None:
                raise ValueError("stress spectra need stress_modes")
            shapes = self.stress_modes
        else:
            shapes = self.eig.modes
        return response_spectrum_from_case(self.eig.omega, shapes, self.eig.participation_factors, case, quantity)

    def random_response(self, case: RandomResponseCase, frequencies_hz: np.ndarray | None = None,
                        bias_model: str = "abaqus") -> RandomResponseResult:
        """Case-file random response (RMS U/V/A and, with stress modes, RMISES)."""
        return random_response_from_case(
            self.eig.omega, self.eig.modes, self.eig.participation_factors, case,
            stress_modes=self.stress_modes, frequencies_hz=frequencies_hz, bias_model=bias_model,
        )
