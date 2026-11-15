"""Linear implicit Bossak-Newmark time integration for structural dynamics.

Implements the unconditionally stable second-order Bossak-Newmark scheme for
linear structural systems :math:`M \\ddot{\\mathbf{u}} + C \\dot{\\mathbf{u}} + K \\mathbf{u} = \\mathbf{f}(t)`
with V10 preset dissipation profiles (``TRANSIENT_FIDELITY`` and
``MODERATE_DISSIPATION``).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ForceInput = np.ndarray | Callable[[float], np.ndarray]

TRANSIENT_FIDELITY_ALPHA_B: float = -0.05
MODERATE_DISSIPATION_ALPHA_B: float = -0.30


@dataclass(frozen=True)
class BossakConfig:
    """Configuration for the linear Bossak-Newmark implicit integrator.

    By default, ``gamma`` and ``beta`` are set to the unconditionally stable,
    second-order-accurate values for ``alpha_b <= 0``:
    ``gamma = 0.5 - alpha_b`` and ``beta = 0.25 * (1.0 - alpha_b) ** 2``.

    Attributes:
        dt: Time step size in seconds (must be positive).
        n_steps: Number of time steps to integrate (must be >= 1).
        alpha_b: Bossak inertial weighting parameter (must be <= 0).
        gamma: Newmark velocity integration parameter.
        beta: Newmark displacement integration parameter.
    """

    dt: float
    n_steps: int
    alpha_b: float = TRANSIENT_FIDELITY_ALPHA_B
    gamma: float | None = None
    beta: float | None = None

    def __post_init__(self) -> None:
        if self.dt <= 0.0:
            raise ValueError(f"dt must be positive, got {self.dt}")
        if self.n_steps < 1:
            raise ValueError(f"n_steps must be >= 1, got {self.n_steps}")
        if self.alpha_b > 0.0:
            raise ValueError(f"alpha_b must be <= 0 for stability, got {self.alpha_b}")

        if self.gamma is None:
            object.__setattr__(self, "gamma", 0.5 - self.alpha_b)
        if self.beta is None:
            object.__setattr__(self, "beta", 0.25 * (1.0 - self.alpha_b) ** 2)

    @classmethod
    def transient_fidelity(cls, dt: float, n_steps: int) -> BossakConfig:
        """Create the V10 ``TRANSIENT_FIDELITY`` preset (``alpha_b = -0.05``)."""
        return cls(dt=dt, n_steps=n_steps, alpha_b=TRANSIENT_FIDELITY_ALPHA_B)

    @classmethod
    def moderate_dissipation(cls, dt: float, n_steps: int) -> BossakConfig:
        """Create the V10 ``MODERATE_DISSIPATION`` preset (``alpha_b = -0.30``)."""
        return cls(dt=dt, n_steps=n_steps, alpha_b=MODERATE_DISSIPATION_ALPHA_B)


BOSSAK_PRESETS: dict[str, float] = {
    "TRANSIENT_FIDELITY": TRANSIENT_FIDELITY_ALPHA_B,
    "MODERATE_DISSIPATION": MODERATE_DISSIPATION_ALPHA_B,
}


@dataclass(frozen=True)
class DynamicHistoryResult:
    """Time-history output of a transient structural dynamics simulation.

    Attributes:
        t: Time instants, shape ``(n_steps + 1,)``.
        u: Displacement history, shape ``(n_steps + 1, n_dofs)``.
        v: Velocity history, shape ``(n_steps + 1, n_dofs)``.
        a: Acceleration history, shape ``(n_steps + 1, n_dofs)``.
        energy_kinetic: Kinetic energy ``0.5 * v.T @ M @ v``, shape ``(n_steps + 1,)``.
        energy_potential: Elastic strain energy ``0.5 * u.T @ K @ u``, shape ``(n_steps + 1,)``.
        energy_total: Total mechanical energy ``energy_kinetic + energy_potential``.
    """

    t: np.ndarray
    u: np.ndarray
    v: np.ndarray
    a: np.ndarray
    energy_kinetic: np.ndarray
    energy_potential: np.ndarray
    energy_total: np.ndarray


def _eval_force(f: ForceInput, t_val: float, n_dofs: int) -> np.ndarray:
    """Evaluate constant or time-dependent force vector at time ``t_val``."""
    vec = np.asarray(f(t_val) if callable(f) else f, dtype=float).reshape(-1)
    if vec.size != n_dofs:
        raise ValueError(f"Force vector length ({vec.size}) does not match n_dofs ({n_dofs}).")
    return vec


def solve_bossak_linear(
    M: sp.spmatrix | np.ndarray,
    C: sp.spmatrix | np.ndarray | None,
    K: sp.spmatrix | np.ndarray,
    f: ForceInput,
    config: BossakConfig,
    u0: np.ndarray | None = None,
    v0: np.ndarray | None = None,
    a0: np.ndarray | None = None,
    free_dofs: np.ndarray | Sequence[int] | None = None,
) -> DynamicHistoryResult:
    """Integrate linear structural dynamics using the Bossak-Newmark scheme.

    Solves at each step ``t_{n+1}``:
    ``(1 - alpha_b) M a_{n+1} + alpha_b M a_n + C v_{n+1} + K u_{n+1} = f_{n+1}``
    with Newmark kinematic updates and a single pre-factorized effective matrix.

    Args:
        M: Mass matrix of shape ``(n_dofs, n_dofs)``.
        C: Damping matrix of shape ``(n_dofs, n_dofs)``, or ``None`` for undamped.
        K: Stiffness matrix of shape ``(n_dofs, n_dofs)``.
        f: Constant force vector of shape ``(n_dofs,)`` or callable ``f(t)``.
        config: ``BossakConfig`` specifying ``dt``, ``n_steps``, ``alpha_b``,
            ``gamma``, and ``beta``.
        u0: Initial displacement vector, shape ``(n_dofs,)``. Defaults to zeros.
        v0: Initial velocity vector, shape ``(n_dofs,)``. Defaults to zeros.
        a0: Initial acceleration vector, shape ``(n_dofs,)``. If ``None``,
            computed from initial equilibrium ``M a0 = f(0) - C v0 - K u0``.
        free_dofs: Optional indices of unconstrained DOFs. Constrained DOFs
            remain identically zero throughout the simulation.

    Returns:
        ``DynamicHistoryResult`` with displacement, velocity, acceleration,
        and mechanical energy time histories.
    """
    M_csr = sp.csr_matrix(M, dtype=float)
    K_csr = sp.csr_matrix(K, dtype=float)
    n_dofs = M_csr.shape[0]
    if M_csr.shape != (n_dofs, n_dofs) or K_csr.shape != (n_dofs, n_dofs):
        raise ValueError("M and K must be square matrices of identical shape.")

    if C is None:
        C_csr = sp.csr_matrix((n_dofs, n_dofs), dtype=float)
    else:
        C_csr = sp.csr_matrix(C, dtype=float)
        if C_csr.shape != (n_dofs, n_dofs):
            raise ValueError("C must have the same shape as M and K.")

    if free_dofs is None:
        free = np.arange(n_dofs, dtype=np.int64)
    else:
        free = np.unique(np.asarray(free_dofs, dtype=np.int64))
    if len(free) == 0:
        raise ValueError("free_dofs must contain at least one DOF.")

    M_ff = M_csr[free, :][:, free]
    C_ff = C_csr[free, :][:, free]
    K_ff = K_csr[free, :][:, free]

    dt = float(config.dt)
    n_steps = int(config.n_steps)
    alpha_b = float(config.alpha_b)
    gamma = float(config.gamma)  # type: ignore[arg-type]
    beta = float(config.beta)  # type: ignore[arg-type]

    t = np.linspace(0.0, n_steps * dt, n_steps + 1)
    u_hist = np.zeros((n_steps + 1, n_dofs), dtype=float)
    v_hist = np.zeros((n_steps + 1, n_dofs), dtype=float)
    a_hist = np.zeros((n_steps + 1, n_dofs), dtype=float)

    if u0 is not None:
        u_hist[0, free] = np.asarray(u0, dtype=float).reshape(n_dofs)[free]
    if v0 is not None:
        v_hist[0, free] = np.asarray(v0, dtype=float).reshape(n_dofs)[free]

    if a0 is not None:
        a_hist[0, free] = np.asarray(a0, dtype=float).reshape(n_dofs)[free]
    else:
        f0_f = _eval_force(f, 0.0, n_dofs)[free]
        rhs0 = f0_f - C_ff @ v_hist[0, free] - K_ff @ u_hist[0, free]
        m_lu = spla.splu(M_ff.tocsc())
        a_hist[0, free] = m_lu.solve(rhs0)

    # Pre-factorize constant effective dynamic matrix:
    # K_eff = (1 - alpha_b) * M + gamma * dt * C + beta * dt^2 * K
    K_eff = (1.0 - alpha_b) * M_ff + (gamma * dt) * C_ff + (beta * dt * dt) * K_ff
    k_eff_lu = spla.splu(K_eff.tocsc())

    e_kin = np.zeros(n_steps + 1, dtype=float)
    e_pot = np.zeros(n_steps + 1, dtype=float)

    u_f = u_hist[0, free].copy()
    v_f = v_hist[0, free].copy()
    a_f = a_hist[0, free].copy()

    e_kin[0] = 0.5 * float(v_f @ (M_ff @ v_f))
    e_pot[0] = 0.5 * float(u_f @ (K_ff @ u_f))

    c_u_pred = 0.5 - beta
    c_v_pred = 1.0 - gamma
    c_u_corr = beta * dt * dt
    c_v_corr = gamma * dt

    for step in range(n_steps):
        t_next = float(t[step + 1])
        u_pred = u_f + dt * v_f + (dt * dt * c_u_pred) * a_f
        v_pred = v_f + (dt * c_v_pred) * a_f

        f_next_f = _eval_force(f, t_next, n_dofs)[free]
        rhs = f_next_f - alpha_b * (M_ff @ a_f) - C_ff @ v_pred - K_ff @ u_pred

        a_next = k_eff_lu.solve(rhs)
        u_next = u_pred + c_u_corr * a_next
        v_next = v_pred + c_v_corr * a_next

        u_f, v_f, a_f = u_next, v_next, a_next
        u_hist[step + 1, free] = u_f
        v_hist[step + 1, free] = v_f
        a_hist[step + 1, free] = a_f

        e_kin[step + 1] = 0.5 * float(v_f @ (M_ff @ v_f))
        e_pot[step + 1] = 0.5 * float(u_f @ (K_ff @ u_f))

    return DynamicHistoryResult(
        t=t,
        u=u_hist,
        v=v_hist,
        a=a_hist,
        energy_kinetic=e_kin,
        energy_potential=e_pot,
        energy_total=e_kin + e_pot,
    )
