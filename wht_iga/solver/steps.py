"""Multi-step analysis orchestration, amplitudes, variable mapping, THB thin-feature retry, and progress/cancellation."""

from __future__ import annotations

import math
import threading
from dataclasses import dataclass, field
from typing import Any, Callable, Literal, Sequence

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.grid.hierarchy import THBHierarchy
from wht_iga.stabilize.extension import THIN_FEATURE


# ---------------------------------------------------------------------------
# Bidirectional Abaqus <-> Kratos field variable mapping
# ---------------------------------------------------------------------------

ABAQUS_TO_KRATOS: dict[str, str] = {
    "U": "DISPLACEMENT",
    "S": "CAUCHY_STRESS_TENSOR",
    "E": "GREEN_LAGRANGE_STRAIN_TENSOR",
    "NT11": "TEMPERATURE",
    "RF": "REACTION",
    "CF": "CONTACT_FORCE",
}

KRATOS_TO_ABAQUS: dict[str, str] = {v: k for k, v in ABAQUS_TO_KRATOS.items()}

# Bidirectional lookup dictionary: VAR_MAP["U"] == "DISPLACEMENT" and
# VAR_MAP["DISPLACEMENT"] == "U".
VAR_MAP: dict[str, str] = {**ABAQUS_TO_KRATOS, **KRATOS_TO_ABAQUS}


def resolve_var_name(
    name: str,
    target: Literal["kratos", "abaqus", "opposite"] = "kratos",
) -> str:
    """Resolve an Abaqus or Kratos field variable alias.

    Parameters
    ----------
    name:
        Variable identifier such as ``"U"``, ``"DISPLACEMENT"``, ``"S"``,
        ``"CAUCHY_STRESS_TENSOR"``, ``"E"``, ``"GREEN_LAGRANGE_STRAIN_TENSOR"``,
        ``"NT11"``, ``"TEMPERATURE"``, ``"RF"``, ``"REACTION"``, ``"CF"``,
        or ``"CONTACT_FORCE"``.
    target:
        ``"kratos"`` returns the canonical Kratos name, ``"abaqus"`` returns the
        canonical Abaqus code, and ``"opposite"`` toggles via ``VAR_MAP``.
    """
    key = str(name).strip().upper()
    if key not in VAR_MAP:
        raise ValueError(f"Unknown field variable {name!r}; known aliases: {sorted(VAR_MAP)}")
    if target == "opposite":
        return VAR_MAP[key]
    if target == "kratos":
        return ABAQUS_TO_KRATOS.get(key, key)
    if target == "abaqus":
        return KRATOS_TO_ABAQUS.get(key, key)
    raise ValueError(f"target must be 'kratos', 'abaqus', or 'opposite', got {target!r}")


# ---------------------------------------------------------------------------
# Cancellation & Progress reporting
# ---------------------------------------------------------------------------

class CancelledError(RuntimeError):
    """Raised when an analysis run is aborted via :class:`CancelToken`."""


class CancelToken:
    """Thread-safe cooperative cancellation token for step execution."""

    def __init__(self) -> None:
        self._event = threading.Event()
        self._reason: str = ""

    def cancel(self, reason: str = "Cancelled by caller") -> None:
        """Signal cancellation."""
        self._reason = str(reason)
        self._event.set()

    @property
    def is_cancelled(self) -> bool:
        """Return ``True`` if cancellation has been requested."""
        return self._event.is_set()

    @property
    def reason(self) -> str:
        """Return the cancellation reason string."""
        return self._reason

    def check(self) -> None:
        """Raise :class:`CancelledError` if cancellation was triggered."""
        if self._event.is_set():
            raise CancelledError(self._reason or "Execution cancelled")


@dataclass(frozen=True)
class ProgressEvent:
    """Structured progress update emitted during step/increment loops."""

    phase: str
    step_index: int
    total_steps: int
    increment: int
    total_increments: int
    time: float
    progress: float
    message: str = ""


def emit_progress(
    on_progress: Callable[[ProgressEvent], Any] | None,
    event: ProgressEvent,
) -> None:
    """Invoke progress callback if supplied."""
    if on_progress is not None:
        on_progress(event)


# ---------------------------------------------------------------------------
# Time Amplitudes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Amplitude:
    """Time-dependent scalar multiplier ``A(t)`` supporting Abaqus-style modes.

    Supported ``kind`` values:
    - ``"tabular"``: piecewise-linear interpolation of ``(time, value)`` pairs.
    - ``"ramp"``: linear ramp from ``(t0, a0)`` to ``(t1, a1)``.
    - ``"smooth_step"``: C2 quintic Hermite transition ``x^3 (10 - 15x + 6x^2)``.
    - ``"periodic"``: Fourier series ``a0 + sum_k (A_k cos(k*w*(t-t0)) + B_k sin(k*w*(t-t0)))``.
    """

    kind: Literal["tabular", "ramp", "smooth_step", "periodic"] = "ramp"
    name: str = "default"
    data: tuple[tuple[float, float], ...] = ()
    t0: float = 0.0
    t1: float = 1.0
    a0: float = 0.0
    a1: float = 1.0
    omega: float = 2.0 * math.pi
    terms: tuple[tuple[float, float], ...] = ((0.0, 1.0),)

    def __post_init__(self) -> None:
        if self.kind not in ("tabular", "ramp", "smooth_step", "periodic"):
            raise ValueError(f"Unsupported amplitude kind {self.kind!r}")
        if self.kind == "tabular":
            pts = tuple(sorted(((float(t), float(v)) for t, v in self.data), key=lambda p: p[0]))
            if len(pts) == 0:
                raise ValueError("tabular amplitude requires at least one (time, value) pair")
            object.__setattr__(self, "data", pts)
        else:
            object.__setattr__(self, "data", tuple((float(t), float(v)) for t, v in self.data))
        object.__setattr__(self, "terms", tuple((float(a), float(b)) for a, b in self.terms))

    @classmethod
    def tabular(cls, points: Sequence[tuple[float, float]], name: str = "tabular") -> Amplitude:
        return cls(kind="tabular", name=name, data=tuple((float(t), float(v)) for t, v in points))

    @classmethod
    def ramp(
        cls,
        t0: float = 0.0,
        t1: float = 1.0,
        a0: float = 0.0,
        a1: float = 1.0,
        name: str = "ramp",
    ) -> Amplitude:
        return cls(kind="ramp", name=name, t0=float(t0), t1=float(t1), a0=float(a0), a1=float(a1))

    @classmethod
    def smooth_step(
        cls,
        t0: float = 0.0,
        t1: float = 1.0,
        a0: float = 0.0,
        a1: float = 1.0,
        name: str = "smooth_step",
    ) -> Amplitude:
        return cls(kind="smooth_step", name=name, t0=float(t0), t1=float(t1), a0=float(a0), a1=float(a1))

    @classmethod
    def periodic(
        cls,
        omega: float = 2.0 * math.pi,
        a0: float = 0.0,
        terms: Sequence[tuple[float, float]] = ((0.0, 1.0),),
        t0: float = 0.0,
        name: str = "periodic",
    ) -> Amplitude:
        return cls(
            kind="periodic",
            name=name,
            t0=float(t0),
            a0=float(a0),
            omega=float(omega),
            terms=tuple((float(c), float(s)) for c, s in terms),
        )

    def evaluate(self, t: float) -> float:
        """Evaluate scalar amplitude factor ``A(t)``."""
        tf = float(t)
        if self.kind == "tabular":
            ts = np.array([p[0] for p in self.data], dtype=np.float64)
            vs = np.array([p[1] for p in self.data], dtype=np.float64)
            return float(np.interp(tf, ts, vs))

        if self.kind == "ramp":
            if self.t1 <= self.t0:
                return float(self.a1 if tf >= self.t0 else self.a0)
            xi = max(0.0, min(1.0, (tf - self.t0) / (self.t1 - self.t0)))
            return float(self.a0 + (self.a1 - self.a0) * xi)

        if self.kind == "smooth_step":
            if self.t1 <= self.t0:
                return float(self.a1 if tf >= self.t0 else self.a0)
            xi = max(0.0, min(1.0, (tf - self.t0) / (self.t1 - self.t0)))
            s = (xi ** 3) * (10.0 - 15.0 * xi + 6.0 * (xi ** 2))
            return float(self.a0 + (self.a1 - self.a0) * s)

        # periodic
        dt = max(0.0, tf - self.t0)
        val = float(self.a0)
        for k, (ak, bk) in enumerate(self.terms, start=1):
            arg = float(k) * self.omega * dt
            val += float(ak) * math.cos(arg) + float(bk) * math.sin(arg)
        return float(val)

    def __call__(self, t: float) -> float:
        return self.evaluate(t)


# ---------------------------------------------------------------------------
# Analysis Step definitions & Multi-Step Runner
# ---------------------------------------------------------------------------

@dataclass
class StaticStep:
    """Static equilibrium step configuration."""

    name: str
    duration: float = 1.0
    n_increments: int = 1
    amplitude: Amplitude | None = None
    nlgeom: bool = False
    loads_scale: float = 1.0
    max_iter: int = 15
    tol: float = 1e-8
    step_type: str = field(default="static", init=False)


@dataclass
class ModalStep:
    """Natural frequency / eigenvalue extraction step configuration."""

    name: str
    n_modes: int = 10
    shift: float = 0.0
    duration: float = 0.0
    n_increments: int = 1
    amplitude: Amplitude | None = None
    step_type: str = field(default="modal", init=False)


@dataclass
class ImplicitDynamicStep:
    """Implicit Newmark-beta / HHT-alpha transient dynamics step.

    Note: Generalized-alpha (Chung-Hulbert 1993) encompasses both HHT-alpha
    (alpha_m=0) and Kratos' Bossak scheme (alpha_f=0) as distinct specializations.
    """

    name: str
    duration: float = 1.0
    n_increments: int = 10
    amplitude: Amplitude | None = None
    beta: float = 0.25
    gamma: float = 0.5
    alpha_hht: float = 0.0
    nlgeom: bool = False
    step_type: str = field(default="implicit_dynamic", init=False)


@dataclass
class ExplicitDynamicStep:
    """Explicit central-difference transient dynamics step."""

    name: str
    duration: float = 1.0
    n_increments: int = 20
    amplitude: Amplitude | None = None
    cfl: float = 0.9
    mass_scaling: float = 1.0
    step_type: str = field(default="explicit_dynamic", init=False)


@dataclass
class ThermalStep:
    """Steady or transient heat-transfer step."""

    name: str
    duration: float = 1.0
    n_increments: int = 1
    amplitude: Amplitude | None = None
    transient: bool = False
    theta: float = 1.0  # 1.0 = Backward Euler, 0.5 = Crank-Nicolson
    step_type: str = field(default="thermal", init=False)


AnalysisStep = StaticStep | ModalStep | ImplicitDynamicStep | ExplicitDynamicStep | ThermalStep


@dataclass
class StepState:
    """Carried physical state across consecutive steps in ``MultiStepRunner``."""

    time: float
    step_index: int
    u: np.ndarray
    v: np.ndarray
    a: np.ndarray
    T: np.ndarray
    eigenvalues: np.ndarray | None = None
    eigenvectors: np.ndarray | None = None
    history: list[dict[str, Any]] = field(default_factory=list)

    def copy(self) -> StepState:
        return StepState(
            time=float(self.time),
            step_index=int(self.step_index),
            u=self.u.copy(),
            v=self.v.copy(),
            a=self.a.copy(),
            T=self.T.copy(),
            eigenvalues=None if self.eigenvalues is None else self.eigenvalues.copy(),
            eigenvectors=None if self.eigenvectors is None else self.eigenvectors.copy(),
            history=list(self.history),
        )


def _solve_dense_or_sparse(A: np.ndarray | sp.spmatrix, b: np.ndarray) -> np.ndarray:
    """Solve a real square linear system cleanly whether sparse or dense."""
    if sp.issparse(A):
        return np.asarray(spla.spsolve(A.tocsc(), b), dtype=np.float64)
    return np.asarray(np.linalg.solve(np.asarray(A, dtype=np.float64), b), dtype=np.float64)


class MultiStepRunner:
    """Execute sequential multi-physics/multi-step analyses while chaining state ``(u, v, a, T)``.

    Supports either:
    1) A custom ``step_handler(step, step_idx, state, t_step, dt, amp, cancel_token) -> StepState``
       for full solver/Kratos orchestration, or
    2) Real assembled matrix operators ``(K, M, C, K_th, C_th)`` and load callables
       ``f_ext(step, t_abs, amp_val, state)`` / ``q_ext(step, t_abs, amp_val, state)``
       solving static, modal, implicit Newmark, explicit central difference, and
       thermal systems increment-by-increment.
    """

    def __init__(
        self,
        n_mech_dofs: int = 0,
        n_thermal_dofs: int = 0,
        *,
        K: np.ndarray | sp.spmatrix | None = None,
        M: np.ndarray | sp.spmatrix | None = None,
        C: np.ndarray | sp.spmatrix | None = None,
        K_th: np.ndarray | sp.spmatrix | None = None,
        C_th: np.ndarray | sp.spmatrix | None = None,
        f_ext: Callable[[AnalysisStep, float, float, StepState], np.ndarray] | np.ndarray | None = None,
        q_ext: Callable[[AnalysisStep, float, float, StepState], np.ndarray] | np.ndarray | None = None,
        step_handler: Callable[[AnalysisStep, int, StepState, float, float, float, CancelToken | None], StepState] | None = None,
    ) -> None:
        if K is not None and n_mech_dofs == 0:
            n_mech_dofs = int(K.shape[0])
        if K_th is not None and n_thermal_dofs == 0:
            n_thermal_dofs = int(K_th.shape[0])

        self.n_mech_dofs = int(n_mech_dofs)
        self.n_thermal_dofs = int(n_thermal_dofs)
        self.K = K
        self.M = M
        self.C = C
        self.K_th = K_th
        self.C_th = C_th
        self.f_ext = f_ext
        self.q_ext = q_ext
        self.step_handler = step_handler

    def _eval_mech_force(self, step: AnalysisStep, t_abs: float, amp_val: float, state: StepState) -> np.ndarray:
        if self.f_ext is None:
            return np.zeros(self.n_mech_dofs, dtype=np.float64)
        if callable(self.f_ext):
            return np.asarray(self.f_ext(step, t_abs, amp_val, state), dtype=np.float64).ravel()
        return np.asarray(self.f_ext, dtype=np.float64).ravel() * float(amp_val)

    def _eval_thermal_load(self, step: AnalysisStep, t_abs: float, amp_val: float, state: StepState) -> np.ndarray:
        if self.q_ext is None:
            return np.zeros(self.n_thermal_dofs, dtype=np.float64)
        if callable(self.q_ext):
            return np.asarray(self.q_ext(step, t_abs, amp_val, state), dtype=np.float64).ravel()
        return np.asarray(self.q_ext, dtype=np.float64).ravel() * float(amp_val)

    def run(
        self,
        steps: Sequence[AnalysisStep],
        initial_state: StepState | None = None,
        *,
        cancel_token: CancelToken | None = None,
        on_progress: Callable[[ProgressEvent], Any] | None = None,
    ) -> StepState:
        """Execute all steps sequentially and return the final chained state."""
        state = (
            initial_state.copy()
            if initial_state is not None
            else StepState(
                time=0.0,
                step_index=0,
                u=np.zeros(self.n_mech_dofs, dtype=np.float64),
                v=np.zeros(self.n_mech_dofs, dtype=np.float64),
                a=np.zeros(self.n_mech_dofs, dtype=np.float64),
                T=np.zeros(self.n_thermal_dofs, dtype=np.float64),
            )
        )

        total_steps = len(steps)
        for s_idx, step in enumerate(steps, start=1):
            if cancel_token is not None:
                cancel_token.check()

            state.step_index = s_idx
            n_inc = max(int(getattr(step, "n_increments", 1)), 1)
            dur = float(getattr(step, "duration", 0.0))

            if isinstance(step, ModalStep):
                amp_val = step.amplitude.evaluate(state.time) if step.amplitude is not None else 1.0
                if self.step_handler is not None:
                    state = self.step_handler(step, s_idx, state, state.time, 0.0, amp_val, cancel_token)
                else:
                    if self.K is None or self.M is None:
                        raise ValueError("ModalStep requires K and M operators or a step_handler")
                    K_dense = self.K.toarray() if sp.issparse(self.K) else np.asarray(self.K, dtype=np.float64)
                    M_dense = self.M.toarray() if sp.issparse(self.M) else np.asarray(self.M, dtype=np.float64)
                    k_req = min(int(step.n_modes), max(K_dense.shape[0], 1))
                    if K_dense.shape[0] <= 256 or k_req >= K_dense.shape[0] - 1:
                        vals, vecs = la.eigh(K_dense, M_dense)
                        idx = np.argsort(vals)
                        state.eigenvalues = vals[idx[:k_req]]
                        state.eigenvectors = vecs[:, idx[:k_req]]
                    else:
                        vals, vecs = spla.eigsh(
                            sp.csc_matrix(K_dense),
                            k=k_req,
                            M=sp.csc_matrix(M_dense),
                            sigma=float(step.shift),
                            which="LM",
                        )
                        idx = np.argsort(vals)
                        state.eigenvalues = vals[idx]
                        state.eigenvectors = vecs[:, idx]
                state.history.append(
                    {
                        "step_index": s_idx,
                        "step_name": step.name,
                        "step_type": step.step_type,
                        "increment": 1,
                        "time": float(state.time),
                        "dt": 0.0,
                        "amplitude": float(amp_val),
                        "max_u": float(np.max(np.abs(state.u))) if state.u.size else 0.0,
                        "max_T": float(np.max(np.abs(state.T))) if state.T.size else 0.0,
                        "eigenvalues": None if state.eigenvalues is None else state.eigenvalues.tolist(),
                    }
                )
                emit_progress(
                    on_progress,
                    ProgressEvent(
                        phase=step.step_type,
                        step_index=s_idx,
                        total_steps=total_steps,
                        increment=1,
                        total_increments=1,
                        time=float(state.time),
                        progress=float(s_idx) / max(total_steps, 1),
                        message=f"Completed modal step {step.name}",
                    ),
                )
                continue

            dt = dur / float(n_inc) if n_inc > 0 else dur
            step_t0 = float(state.time)

            for inc in range(1, n_inc + 1):
                if cancel_token is not None:
                    cancel_token.check()

                t_local = float(inc) * dt
                t_abs = step_t0 + t_local
                if step.amplitude is not None:
                    amp_val = step.amplitude.evaluate(t_local if step.amplitude.t0 == 0.0 and dur > 0 else t_abs)
                else:
                    amp_val = float(inc) / float(n_inc) if isinstance(step, StaticStep) else 1.0

                if self.step_handler is not None:
                    state = self.step_handler(step, s_idx, state, t_abs, dt, amp_val, cancel_token)
                else:
                    self._advance_builtin_step(step, state, t_abs, dt, amp_val)

                state.time = t_abs
                state.history.append(
                    {
                        "step_index": s_idx,
                        "step_name": step.name,
                        "step_type": step.step_type,
                        "increment": inc,
                        "time": float(t_abs),
                        "dt": float(dt),
                        "amplitude": float(amp_val),
                        "max_u": float(np.max(np.abs(state.u))) if state.u.size else 0.0,
                        "max_v": float(np.max(np.abs(state.v))) if state.v.size else 0.0,
                        "max_T": float(np.max(np.abs(state.T))) if state.T.size else 0.0,
                    }
                )

                frac = ((s_idx - 1) + float(inc) / float(n_inc)) / max(float(total_steps), 1.0)
                emit_progress(
                    on_progress,
                    ProgressEvent(
                        phase=step.step_type,
                        step_index=s_idx,
                        total_steps=total_steps,
                        increment=inc,
                        total_increments=n_inc,
                        time=float(t_abs),
                        progress=float(frac),
                        message=f"Step {step.name} inc {inc}/{n_inc}",
                    ),
                )

        return state

    def _advance_builtin_step(
        self,
        step: AnalysisStep,
        state: StepState,
        t_abs: float,
        dt: float,
        amp_val: float,
    ) -> None:
        """Execute one time/load increment using assembled operators."""
        if isinstance(step, StaticStep):
            if self.K is None:
                raise ValueError("StaticStep requires K operator or step_handler")
            f = self._eval_mech_force(step, t_abs, amp_val * step.loads_scale, state)
            state.u = _solve_dense_or_sparse(self.K, f)
            return

        if isinstance(step, ImplicitDynamicStep):
            if self.K is None or self.M is None:
                raise ValueError("ImplicitDynamicStep requires K and M operators or step_handler")
            if dt <= 0.0:
                raise ValueError(f"ImplicitDynamicStep requires dt > 0, got {dt}")
            beta = float(step.beta)
            gamma = float(step.gamma)
            C_op = self.C if self.C is not None else 0.0 * self.M
            a0 = 1.0 / (beta * dt * dt)
            a1 = gamma / (beta * dt)
            a2 = 1.0 / (beta * dt)
            a3 = 1.0 / (2.0 * beta) - 1.0
            a4 = gamma / beta - 1.0
            a5 = 0.5 * dt * (gamma / beta - 2.0)

            f_next = self._eval_mech_force(step, t_abs, amp_val, state)
            K_eff = self.K + a0 * self.M + a1 * C_op
            rhs = (
                f_next
                + self.M @ (a0 * state.u + a2 * state.v + a3 * state.a)
                + C_op @ (a1 * state.u + a4 * state.v + a5 * state.a)
            )
            u_next = _solve_dense_or_sparse(K_eff, rhs)
            a_next = a0 * (u_next - state.u) - a2 * state.v - a3 * state.a
            v_next = state.v + dt * ((1.0 - gamma) * state.a + gamma * a_next)
            state.u, state.v, state.a = u_next, v_next, a_next
            return

        if isinstance(step, ExplicitDynamicStep):
            if self.K is None or self.M is None:
                raise ValueError("ExplicitDynamicStep requires K and M operators or step_handler")
            if dt <= 0.0:
                raise ValueError(f"ExplicitDynamicStep requires dt > 0, got {dt}")
            M_eff = self.M * float(step.mass_scaling)
            # Use lumped diagonal mass vector for explicit central difference.
            if sp.issparse(M_eff):
                m_diag = np.asarray(M_eff.sum(axis=1)).ravel()
            else:
                m_diag = np.asarray(np.asarray(M_eff).sum(axis=1)).ravel()
            if np.any(m_diag <= 0.0):
                raise ValueError("ExplicitDynamicStep lumped mass diagonal has non-positive entry")

            f_cur = self._eval_mech_force(step, t_abs - dt, amp_val, state)
            damping_f = (self.C @ state.v) if self.C is not None else 0.0
            if np.allclose(state.a, 0.0) and state.u.size > 0:
                state.a = (f_cur - damping_f - self.K @ state.u) / m_diag
            v_half = state.v + 0.5 * dt * state.a
            u_next = state.u + dt * v_half
            f_next = self._eval_mech_force(step, t_abs, amp_val, state)
            damping_half = (self.C @ v_half) if self.C is not None else 0.0
            a_next = (f_next - damping_half - self.K @ u_next) / m_diag
            v_next = v_half + 0.5 * dt * a_next
            state.u, state.v, state.a = u_next, v_next, a_next
            return

        if isinstance(step, ThermalStep):
            if self.K_th is None:
                raise ValueError("ThermalStep requires K_th operator or step_handler")
            q_next = self._eval_thermal_load(step, t_abs, amp_val, state)
            if not step.transient:
                state.T = _solve_dense_or_sparse(self.K_th, q_next)
                return
            if self.C_th is None:
                raise ValueError("Transient ThermalStep requires C_th operator")
            if dt <= 0.0:
                raise ValueError(f"Transient ThermalStep requires dt > 0, got {dt}")
            theta = float(step.theta)
            A_th = (1.0 / dt) * self.C_th + theta * self.K_th
            rhs = q_next + ((1.0 / dt) * self.C_th - (1.0 - theta) * self.K_th) @ state.T
            state.T = _solve_dense_or_sparse(A_th, rhs)
            return

        raise ValueError(f"Unsupported step type: {type(step).__name__}")


# ---------------------------------------------------------------------------
# THIN_FEATURE detection & max-2 THB refinement retry loop
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ThinFeatureReport:
    """Diagnostic summary of thin geometric features or high aspect ratios."""

    detected: bool
    aspect_ratio: float
    min_feature_cells: float
    cells: list[tuple[int, int, int]]
    issues: list[str]
    recommended_level: int


def _min_contiguous_run(mask_1d: np.ndarray) -> int:
    """Shortest contiguous True run along a 1D boolean slice (inf if none)."""
    n = int(mask_1d.shape[0])
    best = n + 1
    cur = 0
    for val in mask_1d:
        if bool(val):
            cur += 1
        elif cur > 0:
            best = min(best, cur)
            cur = 0
    if cur > 0:
        best = min(best, cur)
    return best


def detect_thin_features(
    bounds_lo: tuple[float, float, float],
    bounds_hi: tuple[float, float, float],
    labels: np.ndarray,
    *,
    thin_threshold_cells: float = 1.5,
    max_aspect_ratio: float = 8.0,
    extension_issues: Sequence[Any] | None = None,
    current_max_level: int = 0,
) -> ThinFeatureReport:
    """Detect thin features from cell labels, geometry aspect ratio, or ``THIN_FEATURE`` diagnostics."""
    arr = np.asarray(labels)
    if arr.ndim != 3:
        raise ValueError(f"labels must be 3D (nx, ny, nz), got shape {arr.shape}")
    nx, ny, nz = arr.shape
    hx = (float(bounds_hi[0]) - float(bounds_lo[0])) / max(nx, 1)
    hy = (float(bounds_hi[1]) - float(bounds_lo[1])) / max(ny, 1)
    hz = (float(bounds_hi[2]) - float(bounds_lo[2])) / max(nz, 1)
    hs = [max(hx, 1e-30), max(hy, 1e-30), max(hz, 1e-30)]
    cell_ar = max(hs) / min(hs)

    solid = arr >= 0  # IN (1) or CUT (0)
    thin_cells_set: set[tuple[int, int, int]] = set()
    issues: list[str] = []
    min_run = float("inf")

    if np.any(solid):
        # Scan contiguous runs along each grid axis to find thin solid walls/bridges.
        for axis in range(3):
            for idx_other in np.ndindex(*(s for ax, s in enumerate(arr.shape) if ax != axis)):
                sl: list[Any] = list(idx_other)
                sl.insert(axis, slice(None))
                line = solid[tuple(sl)]
                if not np.any(line):
                    continue
                run = _min_contiguous_run(line)
                if run < min_run:
                    min_run = float(run)
                if float(run) <= float(thin_threshold_cells):
                    # Mark active solid cells on lines exhibiting thin thickness.
                    for pos in np.where(line)[0]:
                        full_idx = list(idx_other)
                        full_idx.insert(axis, int(pos))
                        thin_cells_set.add((int(full_idx[0]), int(full_idx[1]), int(full_idx[2])))

        if np.any(thin_cells_set):
            issues.append(
                f"{THIN_FEATURE}: minimum solid thickness {min_run:.2f} cells <= {thin_threshold_cells:.2f}"
            )
    else:
        min_run = 0.0

    has_ext_thin = False
    if extension_issues:
        for iss in extension_issues:
            code = getattr(iss, "code", str(iss))
            if THIN_FEATURE in str(code) or THIN_FEATURE in str(iss):
                has_ext_thin = True
                issues.append(str(iss))

    # High grid aspect ratio also flags cell refinement need on cut cells.
    if cell_ar >= float(max_aspect_ratio):
        issues.append(
            f"{THIN_FEATURE}: cell aspect ratio {cell_ar:.2f} >= {max_aspect_ratio:.2f}"
        )
        cut_locs = np.argwhere(arr == 0)
        for i, j, k in cut_locs:
            thin_cells_set.add((int(i), int(j), int(k)))
        if not thin_cells_set and np.any(solid):
            for i, j, k in np.argwhere(solid):
                thin_cells_set.add((int(i), int(j), int(k)))

    if has_ext_thin and not thin_cells_set:
        cut_locs = np.argwhere(arr == 0)
        target_locs = cut_locs if len(cut_locs) > 0 else np.argwhere(solid)
        for i, j, k in target_locs:
            thin_cells_set.add((int(i), int(j), int(k)))

    detected = bool(issues or thin_cells_set)
    return ThinFeatureReport(
        detected=detected,
        aspect_ratio=float(cell_ar),
        min_feature_cells=float(min_run if math.isfinite(min_run) else 0.0),
        cells=sorted(thin_cells_set),
        issues=issues,
        recommended_level=int(current_max_level + 1) if detected else int(current_max_level),
    )


def auto_refine_thin_features(
    hierarchy: THBHierarchy,
    evaluator: Callable[[THBHierarchy], tuple[np.ndarray, Sequence[Any] | None]],
    *,
    max_retries: int = 2,
    thin_threshold_cells: float = 1.5,
    max_aspect_ratio: float = 8.0,
) -> tuple[THBHierarchy, list[ThinFeatureReport]]:
    """Run up to ``max_retries`` automatic THB local refinement passes on thin features.

    Parameters
    ----------
    hierarchy:
        Initial :class:`THBHierarchy` grid.
    evaluator:
        Callable taking the current hierarchy and returning ``(labels_3d, extension_issues)``
        for the finest level (or active cell grid ``labels_3d`` matching that level's shape).
    max_retries:
        Hard cap on automatic refinement iterations (default ``2``).
    """
    cap = min(int(max_retries), 2)
    reports: list[ThinFeatureReport] = []

    for attempt in range(cap + 1):
        cur_level = int(hierarchy.max_level)
        labels, ext_issues = evaluator(hierarchy)
        bu = hierarchy._breaks_u[cur_level]
        bv = hierarchy._breaks_v[cur_level]
        bw = hierarchy._breaks_w[cur_level]
        rep = detect_thin_features(
            bounds_lo=(float(bu[0]), float(bv[0]), float(bw[0])),
            bounds_hi=(float(bu[-1]), float(bv[-1]), float(bw[-1])),
            labels=labels,
            thin_threshold_cells=thin_threshold_cells,
            max_aspect_ratio=max_aspect_ratio,
            extension_issues=ext_issues,
            current_max_level=cur_level,
        )
        reports.append(rep)
        if not rep.detected or attempt >= cap:
            break
        cells_to_refine = rep.cells
        if not cells_to_refine:
            break
        cells_4d = [(cur_level, i, j, k) for i, j, k in cells_to_refine]
        hierarchy.mark_cells(cells_4d, target_level=cur_level + 1)
        hierarchy.enforce_2_to_1_balance()

    return hierarchy, reports
