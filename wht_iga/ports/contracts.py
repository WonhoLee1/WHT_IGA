"""Data contracts exchanged across the Hexagonal ports (plain dataclasses, Kratos-free).

The session (Layer 4) hands a :class:`StaticJob` to any ``SolverBackendPort`` and
receives a :class:`BackendSolution`; the mesh engine returns a :class:`MeshModel`.
Progress and cooperative cancellation use :class:`ProgressEvent` and
:class:`CancelToken`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np


class AnalysisCancelled(RuntimeError):
    """Raised at the next checkpoint after :meth:`CancelToken.cancel`."""


class CancelToken:
    """Cooperative cancellation flag checked between pipeline stages."""

    def __init__(self) -> None:
        self._cancelled = False
        self.reason = ""

    def cancel(self, reason: str = "cancelled by user") -> None:
        """Request cancellation (takes effect at the next checkpoint)."""
        self._cancelled = True
        self.reason = reason

    @property
    def cancelled(self) -> bool:
        """True once :meth:`cancel` was called."""
        return self._cancelled

    def raise_if_cancelled(self, where: str = "") -> None:
        """Raise :class:`AnalysisCancelled` if cancellation was requested."""
        if self._cancelled:
            raise AnalysisCancelled(f"{self.reason} (at {where})" if where else self.reason)


@dataclass(frozen=True)
class ProgressEvent:
    """Progress notification: ``stage`` name, overall ``fraction`` in [0, 1], message."""

    stage: str
    fraction: float
    message: str = ""


ProgressCallback = Callable[[ProgressEvent], None]


@dataclass(frozen=True)
class MeshModel:
    """Background grid + per-cell quadrature produced by ``build_mesh``."""

    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    n_cells: tuple[int, int, int]
    degree: int
    continuity: int
    trimmed: bool
    qp_table: dict | None = None  # {(i,j,k): (N,4) [x,y,z,w]} for CUT cells
    out_cells: tuple[tuple[int, int, int], ...] = ()
    cell_records: list[dict] | None = None
    geo_issues: list[dict] = field(default_factory=list)
    diag: dict = field(default_factory=dict)
    tess_nodes: np.ndarray | None = None
    tess_tris: np.ndarray | None = None

    def box(self):
        """``wht_iga.reference.box.BoxSpline`` of this grid."""
        from wht_iga.reference.box import BoxSpline

        return BoxSpline(lo=self.lo, hi=self.hi, n_cells=self.n_cells, degree=self.degree,
                         continuity=self.continuity)


@dataclass(frozen=True)
class StaticJob:
    """One static step handed to a solver backend."""

    step_name: str
    job_spec: Any  # wht_iga.solver.validate.JobSpec
    mesh: MeshModel
    extension: dict | None = None


@dataclass
class BackendSolution:
    """Backend-independent solution of one static step."""

    backend: str
    U: np.ndarray  # (n_dofs,) control-point displacements, xyz interleaved
    converged: bool
    n_iter: int
    total_force: np.ndarray  # resultant of the surface loads
    issues: list[dict] = field(default_factory=list)
    timings: dict = field(default_factory=dict)
    solver_log: str = ""
    solver_log_file: str | None = None
    artifacts: dict[str, str] = field(default_factory=dict)
    solver_decision: dict | None = None
    reaction_residual: float | None = None
    native: Any = None  # backend-specific result object (e.g. pipeline.LinearResult)
    #: ``fbar.build_model`` dict when the step was solved with F-bar mean-dilatation
    #: Neo-Hookean; post-processing then reports Cauchy stress / Green-Lagrange strain.
    fbar_model: dict | None = None


@dataclass(frozen=True)
class BackendCapabilities:
    """What a backend can run (checked by the session before solving)."""

    name: str
    trimmed: bool
    partial_component_bcs: bool
    follower_pressure: bool
    gravity: bool
    writes_solver_artifacts: bool  # mesh_db / results_db / vtu written by the backend itself
