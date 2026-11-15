"""Results returned by :meth:`AnalysisSession.run` (thin layer over ``wht_iga.query.ResultView``)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from wht_iga.plan.artifacts import StepArtifacts


@dataclass
class StepResult:
    """Solution, probes and artifact paths of one step."""

    name: str
    backend: str
    U: np.ndarray
    converged: bool
    n_iter: int
    total_force: np.ndarray
    probes: dict[str, dict[str, np.ndarray]]
    artifacts: StepArtifacts
    solver_log: str = ""
    solver_log_file: str | None = None
    solver_decision: dict | None = None
    timings: dict = field(default_factory=dict)
    reaction_residual: float | None = None
    issues: list[dict] = field(default_factory=list)
    vtkhdf_written: bool = False

    @property
    def control_points_u(self) -> np.ndarray:
        """Control-point displacements as ``(n_ctrl, 3)``."""
        return np.asarray(self.U, float).reshape(-1, 3)


class StudyResults:
    """All step results of a run with Coreform-style probe/field access."""

    def __init__(self, steps: dict[str, StepResult]) -> None:
        self.steps = dict(steps)

    def step(self, name: str | None = None) -> StepResult:
        """Named step, or the only/last step when ``name`` is None."""
        if name is None:
            return list(self.steps.values())[-1]
        try:
            return self.steps[name]
        except KeyError:
            raise KeyError(f"no step result {name!r}; have {list(self.steps)}") from None

    def probe(self, name: str, step: str | None = None) -> dict[str, np.ndarray]:
        """Probe values ``{variable: array}`` (one row per evaluation point)."""
        res = self.step(step)
        if name not in res.probes:
            raise KeyError(f"no probe {name!r} in step {res.name!r}; have {list(res.probes)}")
        return res.probes[name]

    def view(self, step: str | None = None) -> Any:
        """``wht_iga.query.ResultView`` over the step's VTKHDF (+ results DB if present)."""
        from wht_iga.query.results import ResultView

        res = self.step(step)
        if not res.vtkhdf_written:
            raise FileNotFoundError(f"step {res.name!r} has no VTKHDF (write_vtkhdf=False or h5py missing)")
        db = res.artifacts.results_db
        return ResultView(res.artifacts.vtkhdf, db if db.exists() else None)

    def field(self, name: str, step: str | None = None, frame: int = 0) -> np.ndarray:
        """Field on the proxy grid points (aliases U/DISPLACEMENT, S/CAUCHY_STRESS_TENSOR)."""
        return self.view(step).field(name, frame=frame)
