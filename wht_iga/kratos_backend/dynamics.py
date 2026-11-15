"""Kratos Bossak implicit dynamics driver and free-vibration analysis (P4-A, V10).

Kratos' QP-geometry elements have no usable consistent mass (spike S06), but their
lumped mass vector is valid.  :func:`run_kratos_bossak` therefore integrates with
``COMPUTE_LUMPED_MASS_MATRIX = True`` (row-sum lumped mass) and the Kratos
``ResidualBasedBossakDisplacementScheme`` (the ``alpha_f = 0`` specialization of
Generalized-alpha, distinct from HHT-alpha which sets ``alpha_m = 0``).  The consistent-mass history comes from
``wht_iga.dynamics.solve_bossak_linear``; both use the same Newmark parameters
``gamma = 1/2 - alpha``, ``beta = (1 - alpha)**2 / 4``.

Free-vibration analysis (:func:`fit_free_vibration`) extracts the frequency and the
envelope decay by variable projection on a damped sinusoid.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from wht_iga.reference.assemble_elastic import assemble_face_traction
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.steps import Amplitude


@dataclass(frozen=True)
class KratosBossakResult:
    """Probe histories of a Kratos Bossak run.

    Attributes:
        t: Time instants ``(n_steps + 1,)`` including ``t = 0``.
        u: Displacement of the probe DOFs, ``(n_steps + 1, n_probe)``.
        v: Velocity of the probe DOFs, ``(n_steps + 1, n_probe)``.
        a: Acceleration of the probe DOFs, ``(n_steps + 1, n_probe)``.
    """

    t: np.ndarray
    u: np.ndarray
    v: np.ndarray
    a: np.ndarray


@dataclass(frozen=True)
class FreeVibrationFit:
    """Damped-sinusoid fit ``u(t) ~ C + exp(-s (t - t0)) (A cos w (t-t0) + B sin w (t-t0))``.

    Attributes:
        frequency_hz: Fitted frequency ``w / (2 pi)``.
        decay_rate: Envelope decay rate ``s`` (1/s); positive means decaying.
        amplitude: Envelope amplitude ``sqrt(A**2 + B**2)`` at the window start.
        phase_rad: Phase ``phi`` in ``amp * exp(-s tau) cos(w tau - phi)``, ``phi = atan2(B, A)``.
        offset: Fitted constant ``C``.
        residual_rms: RMS of the fit residual divided by the amplitude.
    """

    frequency_hz: float
    decay_rate: float
    amplitude: float
    phase_rad: float
    offset: float
    residual_rms: float


def smooth_pulse(period: float) -> Callable[[float], float]:
    """V10 load shape: smooth rise to 1 during ``[0, T]``, smooth fall to 0 during ``[T, 2T]``, then 0.

    Uses the Abaqus quintic smooth step (``x**3 (10 - 15 x + 6 x**2)``) on both ramps.

    Args:
        period: Analytical fundamental period ``T``.

    Returns:
        Function ``A(t)``.
    """
    rise = Amplitude.smooth_step(0.0, period, 0.0, 1.0)
    fall = Amplitude.smooth_step(period, 2.0 * period, 1.0, 0.0)

    def amp(t: float) -> float:
        return rise(t) if t < period else fall(t)

    return amp


def face_traction_force(box: BoxSpline, axis: int, side: int, traction: Sequence[float]) -> np.ndarray:
    """Consistent nodal force vector ``f = int N^T t dA`` of a uniform face traction.

    Args:
        box: Box B-spline.
        axis: Face normal axis.
        side: ``0`` low face, ``1`` high face.
        traction: Traction vector ``(tx, ty, tz)`` per unit area.

    Returns:
        Force vector ``(box.n_dofs,)`` in reference order.
    """
    return np.asarray(assemble_face_traction(box, axis, side, tuple(float(x) for x in traction)), dtype=float)


def kratos_lumped_mass(box: BoxSpline, rho: float, E: float = 1.0, nu: float = 0.0) -> np.ndarray:
    """Lumped mass vector (reference DOF order) from the Kratos elements.

    Args:
        box: Box B-spline.
        rho: Mass density.
        E: Young's modulus (irrelevant to mass).
        nu: Poisson's ratio (irrelevant to mass).

    Returns:
        Vector ``(box.n_dofs,)``; equals the row sums of the consistent mass matrix.
    """
    from wht_iga.kratos_backend.mass import _build_kratos_box

    KM, mp, _volume, elements = _build_kratos_box(box, E, nu, rho)
    lump = np.zeros(box.n_dofs)
    for el in elements:
        el.Initialize(mp.ProcessInfo)
        vec = KM.Vector(0)
        el.CalculateLumpedMassVector(vec, mp.ProcessInfo)
        v = np.asarray(vec)
        ids = [g.Id for g in el.GetGeometry()]
        for a, nid in enumerate(ids):
            lump[3 * (nid - 1) : 3 * (nid - 1) + 3] += v[3 * a : 3 * a + 3]
    return lump


def run_kratos_bossak(
    box: BoxSpline,
    E: float,
    nu: float,
    rho: float,
    clamped_dofs: Sequence[int] | np.ndarray,
    force: np.ndarray,
    amplitude: Callable[[float], float],
    alpha_b: float,
    dt: float,
    n_steps: int,
    probe_dofs: Sequence[int] | np.ndarray,
) -> KratosBossakResult:
    """Linear implicit dynamics with the Kratos Bossak scheme and lumped mass.

    The load ``amplitude(t) * force`` is applied through ``PointLoadCondition3D1N``
    conditions on the loaded control points.  Initial displacement, velocity and
    acceleration are zero.

    Args:
        box: Box B-spline.
        E: Young's modulus.
        nu: Poisson's ratio.
        rho: Mass density.
        clamped_dofs: Reference DOF indices fixed to zero.
        force: Nodal force pattern ``(box.n_dofs,)`` multiplied by ``amplitude(t)``.
        amplitude: Load factor ``A(t)``.
        alpha_b: Bossak parameter (Kratos sign convention, ``-0.3 <= alpha_b <= 0``).
        dt: Constant time step.
        n_steps: Number of steps.
        probe_dofs: Reference DOF indices whose histories are returned.

    Returns:
        :class:`KratosBossakResult`.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.mass import _build_kratos_box

    KM, mp, _volume, _elements = _build_kratos_box(box, E, nu, rho, dynamic=True)
    SM = load_kratos()["StructuralMechanicsApplication"]
    mp.ProcessInfo.SetValue(KM.COMPUTE_LUMPED_MASS_MATRIX, True)

    force = np.asarray(force, dtype=float).reshape(-1)
    loaded_nodes = sorted({int(d) // 3 for d in np.nonzero(force)[0]})
    props = mp.GetProperties()[1]
    next_cond = 1
    for I in loaded_nodes:
        mp.CreateNewCondition("PointLoadCondition3D1N", next_cond, [I + 1], props)
        next_cond += 1
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for dof in np.asarray(clamped_dofs, dtype=np.int64).reshape(-1):
        mp.GetNode(int(dof) // 3 + 1).Fix(comp_vars[int(dof) % 3])

    scheme = KM.ResidualBasedBossakDisplacementScheme(float(alpha_b))
    builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()

    probe = np.asarray(probe_dofs, dtype=np.int64).reshape(-1)
    comp_idx = [(int(d) // 3 + 1, int(d) % 3) for d in probe]

    def read(var_xyz) -> np.ndarray:
        return np.array([mp.GetNode(n).GetSolutionStepValue(var_xyz)[c] for n, c in comp_idx])

    t = np.zeros(n_steps + 1)
    u = np.zeros((n_steps + 1, probe.size))
    v = np.zeros_like(u)
    a = np.zeros_like(u)
    mp.ProcessInfo.SetValue(KM.TIME, 0.0)
    for step in range(1, n_steps + 1):
        now = step * dt
        mp.CloneTimeStep(now)
        scale = float(amplitude(now))
        for I in loaded_nodes:
            mp.GetNode(I + 1).SetSolutionStepValue(SM.POINT_LOAD, [scale * force[3 * I + c] for c in range(3)])
        strategy.Solve()
        t[step] = now
        u[step] = read(KM.DISPLACEMENT)
        v[step] = read(KM.VELOCITY)
        a[step] = read(KM.ACCELERATION)
    return KratosBossakResult(t=t, u=u, v=v, a=a)


def fit_free_vibration(t: np.ndarray, u: np.ndarray, f_guess_hz: float) -> FreeVibrationFit:
    """Fit ``C + exp(-s tau)(A cos w tau + B sin w tau)`` (``tau = t - t[0]``) by variable projection.

    ``w`` and ``s`` are the nonlinear unknowns; ``A``, ``B``, ``C`` follow from linear least squares.

    Args:
        t: Time samples of the free-vibration window.
        u: Displacement samples.
        f_guess_hz: Initial frequency guess in Hz.

    Returns:
        :class:`FreeVibrationFit`.
    """
    tau = np.asarray(t, dtype=float) - float(t[0])
    y = np.asarray(u, dtype=float)
    y_scale = float(np.max(np.abs(y - y.mean()))) or 1.0
    y = y / y_scale  # the optimiser tolerances are absolute: fit a unit-amplitude signal

    def lin(w: float, s: float):
        env = np.exp(-s * tau)
        design = np.column_stack([env * np.cos(w * tau), env * np.sin(w * tau), np.ones_like(tau)])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        return coef, design @ coef - y

    def resid(x: np.ndarray) -> np.ndarray:
        return lin(x[0] * 2.0 * np.pi, x[1])[1]

    # The residual is highly oscillatory in w: a guess that is off by more than about a quarter
    # cycle over the window lies outside the basin of the refinement.  Scan +-5 % first.
    window = float(tau[-1]) if tau.size > 1 else 1.0
    n_scan = max(401, int(0.1 * f_guess_hz * window * 40))
    grid = f_guess_hz * np.linspace(0.95, 1.05, n_scan)
    f_start = float(grid[int(np.argmin([np.sum(lin(2.0 * np.pi * f, 0.0)[1] ** 2) for f in grid]))])
    sol = least_squares(resid, x0=[f_start, 0.0], x_scale=[f_start * 1e-3, 1.0], xtol=1e-14, ftol=1e-14, gtol=1e-14)
    f_hz, s = float(sol.x[0]), float(sol.x[1])
    coef, res = lin(f_hz * 2.0 * np.pi, s)
    amp = float(np.hypot(coef[0], coef[1]))
    return FreeVibrationFit(
        frequency_hz=f_hz,
        decay_rate=s,
        amplitude=amp * y_scale,
        phase_rad=float(np.arctan2(coef[1], coef[0])),
        offset=float(coef[2]) * y_scale,
        residual_rms=float(np.sqrt(np.mean(res * res)) / amp) if amp > 0.0 else float("nan"),
    )
