"""``KratosPythonBackend``: production adapter over ``wht_iga.solver.pipeline.run_linear``.

Kratos-first ergonomics (plan v2.1 §8.4):

* escape hatch -- ``kratos_model`` (``Kratos.Model``), ``root_model_part``
  (``Kratos.ModelPart``) and ``strategy`` are set as soon as the model exists
  (inside ``on_before_solve``) and stay available after the solve;
* ``custom_processes`` -- factories ``f(model, parameters) -> Kratos.Process``;
  ``parameters`` is ``Kratos.Parameters(json.dumps(kratos_parameters or {}))``. Their
  ``ExecuteInitialize / ExecuteBeforeSolutionLoop / ExecuteInitializeSolutionStep``
  run before the solve, ``ExecuteFinalizeSolutionStep / ExecuteFinalize`` after it;
* hooks -- ``on_before_solve(model_part)`` and ``on_step_end(model_part, step_index)``.

Kratos is imported only inside :meth:`solve_static` (``load_kratos`` pin).
"""

from __future__ import annotations

import dataclasses
import json
from typing import Any, Callable, Sequence

from wht_iga.ports.contracts import BackendCapabilities, BackendSolution, StaticJob


class _PipelineHooks:
    """Adapter from ``run_linear(hooks=...)`` callbacks to the backend's user hooks."""

    def __init__(self, backend: KratosPythonBackend, step_index: int) -> None:
        self.backend = backend
        self.step_index = step_index
        self.processes: list[Any] = []

    def before_solve(self, model: Any, model_part: Any, strategy: Any) -> None:
        """Expose native objects, build/initialize custom processes, call the user hook."""
        b = self.backend
        b.kratos_model, b.root_model_part, b.strategy = model, model_part, strategy
        if b.custom_processes:
            from wht_iga.kratos_backend import load_kratos

            KM = load_kratos()["KratosMultiphysics"]
            params = KM.Parameters(json.dumps(b.kratos_parameters or {}))
            self.processes = [factory(model, params) for factory in b.custom_processes]
            for p in self.processes:
                p.ExecuteInitialize()
            for p in self.processes:
                p.ExecuteBeforeSolutionLoop()
            for p in self.processes:
                p.ExecuteInitializeSolutionStep()
        if b.on_before_solve is not None:
            b.on_before_solve(model_part)

    def after_solve(self, model: Any, model_part: Any, strategy: Any) -> None:
        """Finalize custom processes and call the step-end hook."""
        for p in self.processes:
            p.ExecuteFinalizeSolutionStep()
        for p in self.processes:
            p.ExecuteFinalize()
        if self.backend.on_step_end is not None:
            self.backend.on_step_end(model_part, self.step_index)


class KratosPythonBackend:
    """Solver backend running the Kratos source build through ``run_linear``."""

    name = "kratos_python"

    def __init__(self, threads: int | str | None = None, *,
                 custom_processes: Sequence[Callable[[Any, Any], Any]] = (),
                 kratos_parameters: dict | None = None,
                 on_before_solve: Callable[[Any], None] | None = None,
                 on_step_end: Callable[[Any, int], None] | None = None) -> None:
        for f in custom_processes:
            if not callable(f):
                raise TypeError("custom_processes takes factories f(model, parameters) -> Process, "
                                f"got {type(f).__name__}")
        self.threads = threads
        self.custom_processes = tuple(custom_processes)
        self.kratos_parameters = kratos_parameters
        self.on_before_solve = on_before_solve
        self.on_step_end = on_step_end
        self.kratos_model: Any = None
        self.root_model_part: Any = None
        self.strategy: Any = None
        self.last_result: Any = None
        self._step_index = 0

    def capabilities(self) -> BackendCapabilities:
        """Everything ``run_linear`` supports."""
        return BackendCapabilities(name=self.name, trimmed=True, partial_component_bcs=True,
                                   follower_pressure=True, gravity=True, writes_solver_artifacts=True)

    def solve_static(self, job: StaticJob) -> BackendSolution:
        """Run ``run_linear`` for one step and wrap its result."""
        from wht_iga.solver.pipeline import run_linear

        js = job.job_spec
        if self.threads is not None and js.num_threads_mkl is None:
            js = dataclasses.replace(js, num_threads_mkl=self.threads)
        self._step_index += 1
        hooks = _PipelineHooks(self, self._step_index)
        m = job.mesh
        res = run_linear(js, geo_issues=m.geo_issues, qp_table=m.qp_table, cell_records=m.cell_records,
                         extension=job.extension, hooks=hooks)
        self.last_result = res
        self.kratos_model = res.keepalive[2]
        self.root_model_part = res.keepalive[3]
        self.strategy = res.strategy
        arts = {k: v for k, v in (("mesh_db", res.mesh_db), ("results_db", res.results_db),
                                  ("vtu", res.vtu_path), ("probe_json", res.probe_json),
                                  ("solver_log", res.solver_log_file)) if v}
        return BackendSolution(
            backend=self.name, U=res.U, converged=bool(res.converged), n_iter=int(res.n_iter),
            total_force=res.total_force,
            issues=[dataclasses.asdict(i) for i in res.issues], timings=dict(res.timings),
            solver_log=res.solver_log, solver_log_file=res.solver_log_file, artifacts=arts,
            solver_decision=res.solver_decision, native=res)

    def close(self) -> None:
        """Drop references to native Kratos objects."""
        self.kratos_model = self.root_model_part = self.strategy = self.last_result = None
