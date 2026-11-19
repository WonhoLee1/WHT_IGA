"""Layer 4: ``AnalysisSession`` -- lifecycle of one compiled plan (W2 + W4).

::

    with AnalysisSession(plan, backend=KratosPythonBackend()) as session:
        session.build_mesh()
        results = session.run()

* Entering checks the plan (no blocking diagnostics) and the backend capabilities,
  creates the output directory and writes ``<job>.spec.json``.
* ``enable_solver_log`` switches the HWT-IGA SOLVER 6-stage log on/off for every
  step (file ``<prefix>_solver.log``; ``WHT_SOLVER_LOG=0`` also disables it).
* ``on_progress`` receives :class:`ProgressEvent` with a monotone overall fraction.
* Cancellation is cooperative: :meth:`cancel` (or the shared :class:`CancelToken`)
  takes effect at the next stage boundary (a running Kratos solve is not
  interrupted). On :class:`AnalysisCancelled` / ``KeyboardInterrupt`` the files this
  session created or modified are removed (``cleanup_on_cancel=True``); files that
  existed before and were not touched are kept.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from wht_iga.plan.compiler import Plan
from wht_iga.ports.contracts import (
    AnalysisCancelled,
    CancelToken,
    MeshModel,
    ProgressCallback,
    ProgressEvent,
    StaticJob,
)
from wht_iga.runner.results import StepResult, StudyResults
from wht_iga.spec.assembly import BoxGeometrySpec
from wht_iga.spec.physics import PressureLoadSpec, TractionLoadSpec


class SessionError(RuntimeError):
    """Session misuse or a backend/plan mismatch detected before solving."""


class LoadCoverageError(SessionError):
    """A surface load sits on a grid face that the part does not cover."""


class AnalysisSession:
    """Context manager running a compiled :class:`Plan` on a solver backend."""

    def __init__(self, plan: Plan, backend: Any = None, *, mesh_engine: Any = None,
                 enable_solver_log: bool = True, solver_log_stdout: bool = False,
                 solver_log_level: str = "detailed", on_progress: ProgressCallback | None = None,
                 cancel_token: CancelToken | None = None, cleanup_on_cancel: bool = True,
                 write_vtkhdf: bool = True) -> None:
        if backend is None:
            from wht_iga.backends.kratos_python import KratosPythonBackend

            backend = KratosPythonBackend()
        if mesh_engine is None:
            from wht_iga.backends.mesh_engine import GridMeshEngine

            mesh_engine = GridMeshEngine()
        self.plan = plan
        self.backend = backend
        self.mesh_engine = mesh_engine
        self.enable_solver_log = bool(enable_solver_log)
        self.solver_log_stdout = bool(solver_log_stdout)
        self.solver_log_level = solver_log_level
        self.on_progress = on_progress
        self.cancel_token = cancel_token or CancelToken()
        self.cleanup_on_cancel = bool(cleanup_on_cancel)
        self.write_vtkhdf = bool(write_vtkhdf)
        self.mesh: MeshModel | None = None
        self.results: StudyResults | None = None
        self.removed_artifacts: list[Path] = []
        self._active = False
        self._t_enter = 0.0
        self._preexisting: dict[Path, float] = {}
        self._last_fraction = 0.0

    # ------------------------------------------------------------------ lifecycle
    def __enter__(self) -> AnalysisSession:
        self.plan.report.raise_if_invalid()
        self._check_backend()
        out = self.plan.artifacts.output_dir
        out.mkdir(parents=True, exist_ok=True)
        self._preexisting = {p: p.stat().st_mtime for p in self.plan.artifacts.all() if p.exists()}
        self._t_enter = time.time()
        self._active = True
        self.plan.artifacts.spec_json.write_text(self.plan.spec.model_dump_json(indent=2), encoding="utf-8")
        self._emit("session", 0.0, f"session start: {self.plan.spec.name} on {self.backend.name}")
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
        self._active = False
        if exc_type is not None and issubclass(exc_type, (AnalysisCancelled, KeyboardInterrupt)) \
                and self.cleanup_on_cancel:
            self.cleanup_partial_artifacts()
        return False

    def cancel(self, reason: str = "cancelled by user") -> None:
        """Request cancellation at the next stage boundary."""
        self.cancel_token.cancel(reason)

    def cleanup_partial_artifacts(self) -> list[Path]:
        """Delete artifacts created or modified since the session started."""
        removed: list[Path] = []
        for p in self.plan.artifacts.all():
            if not p.exists():
                continue
            before = self._preexisting.get(p)
            if before is None or p.stat().st_mtime > before:
                p.unlink()
                removed.append(p)
        self.removed_artifacts.extend(removed)
        return removed

    # ------------------------------------------------------------------ stages
    def build_mesh(self) -> MeshModel:
        """W2: background grid + cut-cell quadrature (S1-S6)."""
        self._require_active()
        self.cancel_token.raise_if_cancelled("build_mesh start")
        self._emit("mesh", 0.05, "build_mesh start")

        def sub(ev: ProgressEvent) -> None:
            self._emit(ev.stage, 0.05 + 0.30 * ev.fraction, ev.message)

        mesh = self.mesh_engine.build(self.plan, on_progress=sub, cancel=self.cancel_token)
        self._check_load_coverage(mesh)
        self.mesh = mesh
        self.cancel_token.raise_if_cancelled("build_mesh end")
        self._emit("mesh", 0.35, "build_mesh done")
        return mesh

    def run(self) -> StudyResults:
        """W4: solve every step, evaluate probes, write probe JSON / VTKHDF."""
        self._require_active()
        if self.mesh is None:
            self.build_mesh()
        assert self.mesh is not None
        from wht_iga.runner.post import evaluate_point_fields, write_step_vtkhdf

        spec = self.plan.spec
        mat = spec.material_of(spec.parts[0])
        steps = spec.steps
        out: dict[str, StepResult] = {}
        for k, step in enumerate(steps):
            lo_f = 0.35 + 0.60 * k / len(steps)
            span = 0.60 / len(steps)
            self.cancel_token.raise_if_cancelled(f"step {step.name} start")
            self._emit(f"solve:{step.name}", lo_f, f"solving step {step.name}")
            js = self.plan.job_spec(step.name, out_cells=self.mesh.out_cells,
                                    enable_solver_log=self.enable_solver_log,
                                    solver_log_stdout=self.solver_log_stdout,
                                    solver_log_level=self.solver_log_level)
            ext = self.plan.extension() if self.mesh.trimmed else None
            sol = self.backend.solve_static(StaticJob(step.name, js, self.mesh, ext))
            self.cancel_token.raise_if_cancelled(f"step {step.name} solved")
            self._emit(f"post:{step.name}", lo_f + 0.8 * span, f"post-processing step {step.name}")
            box = self.mesh.box()
            fbm = getattr(sol, "fbar_model", None)
            probes = {p.name: evaluate_point_fields(box, sol.U, p.points(), mat.E, mat.nu, p.variables,
                                                    fbar_model=fbm)
                      for p in spec.probes}
            arts = self.plan.artifacts.step(step.name)
            payload = {name: {"points": [p.points().tolist() for p in spec.probes if p.name == name][0],
                              **{v: arr.tolist() for v, arr in vals.items()}}
                       for name, vals in probes.items()}
            arts.probe_json.write_text(json.dumps(payload, indent=1), encoding="utf-8")
            wrote = False
            if self.write_vtkhdf:
                write_step_vtkhdf(arts.vtkhdf, box, sol.U, self.mesh.out_cells, mat.E, mat.nu, fbar_model=fbm)
                wrote = True
            out[step.name] = StepResult(
                name=step.name, backend=sol.backend, U=sol.U, converged=sol.converged, n_iter=sol.n_iter,
                total_force=sol.total_force, probes=probes, artifacts=arts, solver_log=sol.solver_log,
                solver_log_file=sol.solver_log_file, solver_decision=sol.solver_decision,
                timings=sol.timings, reaction_residual=sol.reaction_residual, issues=sol.issues,
                vtkhdf_written=wrote,
                stress_measure="fbar_neo_hookean_cauchy" if fbm is not None else "linear_elastic")
            self._emit(f"solve:{step.name}", lo_f + span, f"step {step.name} done")
        self.results = StudyResults(out)
        self._emit("done", 1.0, "run complete")
        return self.results

    # ------------------------------------------------------------------ helpers
    def _require_active(self) -> None:
        if not self._active:
            raise SessionError("use AnalysisSession as a context manager: with AnalysisSession(...) as s:")

    def _emit(self, stage: str, fraction: float, message: str) -> None:
        self._last_fraction = max(self._last_fraction, min(1.0, float(fraction)))
        if self.on_progress is not None:
            self.on_progress(ProgressEvent(stage, self._last_fraction, message))

    def _check_backend(self) -> None:
        caps = self.backend.capabilities()
        g = self.plan.report.grid
        if g is not None and g.trimmed and not caps.trimmed:
            raise SessionError(f"backend {caps.name!r} cannot run trimmed grids (plan grid is trimmed)")
        for s in self.plan.spec.steps:
            if any(isinstance(ld, PressureLoadSpec) and ld.follower for ld in s.loads) \
                    and not caps.follower_pressure:
                raise SessionError(f"backend {caps.name!r} has no follower pressure (step {s.name!r})")
            if any(not all(b.components) for b in s.bcs) and not caps.partial_component_bcs:
                raise SessionError(f"backend {caps.name!r} has no per-component BCs (step {s.name!r})")

    def _check_load_coverage(self, mesh: MeshModel) -> None:
        """Surface loads on grid faces must not touch OUT cells (load would act on void)."""
        if not mesh.trimmed or isinstance(self.plan.spec.parts[0].geometry, BoxGeometrySpec):
            return  # untrimmed, or box part already checked exactly by the plan
        out = set(mesh.out_cells)
        n = mesh.n_cells
        for s in self.plan.spec.steps:
            for ld in s.loads:
                if not isinstance(ld, (PressureLoadSpec, TractionLoadSpec)):
                    continue
                axis, side = self.plan.spec.set_face(ld.set)
                layer = 0 if side == 0 else n[axis] - 1
                bad = [c for c in out if c[axis] == layer]
                if bad:
                    raise LoadCoverageError(
                        f"step {s.name!r}: {ld.kind} on {ld.set!r}: {len(bad)} OUT cells touch that grid face; "
                        f"run_linear would integrate the load over void")
