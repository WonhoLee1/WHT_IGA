"""Coreform IGA compatibility facade (plan v2.1 §8.3): ``import wht_iga.compat.coreform as cf``.

Names and arguments follow the Coreform *portable spec* layer
(``coreform.mesh`` / ``coreform.probe`` / ``coreform.job``: ``degree``,
``continuity < degree``, ``bc_penalty_value``, ``tessellation_size``, ``sizing``,
``small_cell_volume_ratio`` in ``[0, 0.5)``, ``padding``, ``extents``;
``LineProbe(start_location, stop_location, num_eval_points)``; ``to_dict /
to_json / from_dict / from_json / copy_with / validate``). Every object converts
to the WHT_IGA core spec with :meth:`to_core`, and :class:`Job` delegates to
``StudySpec -> compile_plan -> AnalysisSession -> StudyResults``.

Not supported (raise ``ValueError`` instead of being ignored): Abaqus objects
(``mdb``, ``customData``), non-global ``frame``, ``element_layout``, probe ``region``.
"""

from __future__ import annotations

import enum
import threading
from pathlib import Path
from typing import Annotated, Any, Literal, Self, Union

from pydantic import Field, model_validator

from wht_iga.spec import (
    BodyFitSpec,
    BoxRegionSpec,
    CadEntitiesRegionSpec,
    ImmersedGridSpec,
    LineProbeSpec,
    MultiPointProbeSpec,
    PointProbeSpec,
    RefinementSpec,
    SphereRegionSpec,
    StudySpec,
)
from wht_iga.spec.base import Name, PositiveFinite, SpecBase, Vec3


class _CoreformSpec(SpecBase):
    """Coreform method names on top of the frozen Pydantic base."""

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a versioned dictionary (JSON-safe)."""
        return self.model_dump(mode="json")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Validated instance from a mapping payload."""
        return cls.model_validate(data)

    def copy_with(self, **changes: Any) -> Self:
        """Validated copy with the supplied fields replaced."""
        return self.with_updates(**changes)

    def validate(self) -> list[str]:  # type: ignore[override]
        """Validation issues of this spec: always empty (instances are validated at construction)."""
        return []


# ----------------------------------------------------------------------------- sizing
class ElementSizeSizing(_CoreformSpec):
    """Uniform target element size along every axis."""

    schemaVersion: Literal[1] = 1
    sizing_type: Literal["element_size"] = "element_size"
    element_size: PositiveFinite


class ElementIntervalCountByAxisSizing(_CoreformSpec):
    """Element-interval count per axis."""

    schemaVersion: Literal[1] = 1
    sizing_type: Literal["element_interval_count_by_axis"] = "element_interval_count_by_axis"
    element_interval_count_by_axis: tuple[Annotated[int, Field(ge=1)], Annotated[int, Field(ge=1)],
                                          Annotated[int, Field(ge=1)]]


class ElementIntervalCountSizing(_CoreformSpec):
    """Same element-interval count on all three axes."""

    schemaVersion: Literal[1] = 1
    sizing_type: Literal["element_interval_count"] = "element_interval_count"
    element_interval_count: Annotated[int, Field(ge=1)]


RectilinearSizing = Annotated[Union[ElementSizeSizing, ElementIntervalCountByAxisSizing,
                                    ElementIntervalCountSizing], Field(discriminator="sizing_type")]


def _sizing_kwargs(sizing: Any) -> dict[str, Any]:
    """Core ``element_size`` / ``n_cells`` keyword for a Coreform sizing object."""
    if isinstance(sizing, ElementSizeSizing):
        return {"element_size": sizing.element_size}
    if isinstance(sizing, ElementIntervalCountByAxisSizing):
        return {"n_cells": sizing.element_interval_count_by_axis}
    if isinstance(sizing, ElementIntervalCountSizing):
        n = sizing.element_interval_count
        return {"n_cells": (n, n, n)}
    raise ValueError("sizing is required (ElementSizeSizing / ElementIntervalCount[ByAxis]Sizing)")


# ----------------------------------------------------------------------------- regions / refinement
class CadEntitiesRegion(_CoreformSpec):
    """Region of CAD entity pointers."""

    kind: Literal["cad_entities"] = "cad_entities"
    entities: tuple[Annotated[str, Field(min_length=1)], ...] = Field(min_length=1)

    def to_core(self) -> CadEntitiesRegionSpec:
        """Core region spec."""
        return CadEntitiesRegionSpec(entities=self.entities)


class SphereRegion(_CoreformSpec):
    """Ball region."""

    kind: Literal["sphere"] = "sphere"
    center: Vec3
    radius: PositiveFinite

    def to_core(self) -> SphereRegionSpec:
        """Core region spec."""
        return SphereRegionSpec(center=self.center, radius=self.radius)


class BoxRegion(_CoreformSpec):
    """Axis-aligned box region."""

    kind: Literal["box"] = "box"
    lo: Vec3
    hi: Vec3

    def to_core(self) -> BoxRegionSpec:
        """Core region spec (validates ``hi > lo``)."""
        return BoxRegionSpec(lo=self.lo, hi=self.hi)


class Refinement(_CoreformSpec):
    """Local refinement on a region or named sets."""

    name: Name
    level: Annotated[int, Field(ge=1, le=8)]
    region: Annotated[Union[CadEntitiesRegion, SphereRegion, BoxRegion], Field(discriminator="kind")] | None = None
    sets: tuple[Name, ...] = ()
    balance: Literal["two_to_one", "none"] = "two_to_one"

    def to_core(self) -> RefinementSpec:
        """Core refinement spec."""
        return RefinementSpec(name=self.name, level=self.level, sets=self.sets,
                              region=self.region.to_core() if self.region is not None else None,
                              balance=self.balance)


# ----------------------------------------------------------------------------- meshes
class _MeshCommon(_CoreformSpec):
    """Fields shared by the Coreform mesh specs (Coreform defaults: degree 1, continuity 0)."""

    schemaVersion: Literal[1] = 1
    degree: Annotated[int, Field(ge=1, le=3)] = 1
    continuity: Annotated[int, Field(ge=0)] = 0
    bc_penalty_value: PositiveFinite | None = None
    tessellation_size: PositiveFinite | None = None
    sizing: RectilinearSizing | None = None
    refinements: tuple[Refinement, ...] = ()

    @model_validator(mode="after")
    def _continuity(self) -> Self:
        if self.continuity >= self.degree:
            raise ValueError(f"continuity must be less than degree, got {self.continuity} >= {self.degree}")
        return self


class ImmersedRectilinearMesh(_MeshCommon):
    """Coreform ``ImmersedRectilinearMeshSpec`` -> :class:`wht_iga.spec.ImmersedGridSpec`."""

    meshType: Literal["immersed_rectilinear"] = "immersed_rectilinear"
    element_layout: tuple[str, str, str] | None = None
    small_cell_volume_ratio: Annotated[float, Field(ge=0.0, lt=0.5)] = 0.2
    padding: tuple[Annotated[int, Field(ge=0)], Annotated[int, Field(ge=0)],
                   Annotated[int, Field(ge=0)]] | None = None
    frame: Any = None
    extents: tuple[Vec3, Vec3] | None = None

    @model_validator(mode="after")
    def _unsupported(self) -> Self:
        if self.frame is not None:
            raise ValueError("frame: only the global axis-aligned frame is supported (pass frame=None)")
        if self.element_layout is not None:
            raise ValueError("element_layout is not supported by WHT_IGA (pass element_layout=None)")
        return self

    def to_core(self) -> ImmersedGridSpec:
        """Core immersed-grid spec."""
        return ImmersedGridSpec(
            degree=self.degree, continuity=self.continuity, bc_penalty=self.bc_penalty_value,
            tessellation_size=self.tessellation_size, small_cell_volume_ratio=self.small_cell_volume_ratio,
            padding=self.padding or (0, 0, 0), bounds=self.extents,
            refinements=tuple(r.to_core() for r in self.refinements), **_sizing_kwargs(self.sizing))


class BodyFitMesh(_MeshCommon):
    """Coreform ``BodyFitMeshSpec`` -> :class:`wht_iga.spec.BodyFitSpec` (box parts; ``sizing`` required)."""

    meshType: Literal["body_fit"] = "body_fit"

    def to_core(self) -> BodyFitSpec:
        """Core body-fit spec."""
        return BodyFitSpec(degree=self.degree, continuity=self.continuity, bc_penalty=self.bc_penalty_value,
                           refinements=tuple(r.to_core() for r in self.refinements),
                           **_sizing_kwargs(self.sizing))


def immersed_rectilinear(spec: ImmersedRectilinearMesh | None = None, /, **fields: Any) -> ImmersedRectilinearMesh:
    """Coreform ``coreform.mesh.immersed_rectilinear`` (spec object or spec fields)."""
    return spec if spec is not None else ImmersedRectilinearMesh(**fields)


def body_fit(spec: BodyFitMesh | None = None, /, **fields: Any) -> BodyFitMesh:
    """Coreform ``coreform.mesh.body_fit`` (spec object or spec fields)."""
    return spec if spec is not None else BodyFitMesh(**fields)


# ----------------------------------------------------------------------------- probes
class _ProbeCommon(_CoreformSpec):
    """Fields shared by Coreform probe specs (``variables`` default: empty tuple)."""

    schemaVersion: Literal[1] = 1
    name: Name
    variables: tuple[str, ...] = ()
    region: Any = None

    @model_validator(mode="after")
    def _no_region(self) -> Self:
        if self.region is not None:
            raise ValueError("probe region selection is not supported (pass region=None)")
        return self


class PointProbe(_ProbeCommon):
    """Single-point probe."""

    probeType: Literal["point"] = "point"
    location: Vec3

    def to_core(self) -> PointProbeSpec:
        """Core probe spec."""
        return PointProbeSpec(name=self.name, location=self.location, variables=self.variables)


class LineProbe(_ProbeCommon):
    """Segment probe with ``num_eval_points`` (endpoints included)."""

    probeType: Literal["line"] = "line"
    start_location: Vec3
    stop_location: Vec3
    num_eval_points: Annotated[int, Field(ge=2)] = 11

    def to_core(self) -> LineProbeSpec:
        """Core probe spec."""
        return LineProbeSpec(name=self.name, start=self.start_location, end=self.stop_location,
                             n_points=self.num_eval_points, variables=self.variables)


class MultiPointProbe(_ProbeCommon):
    """Explicit list of locations."""

    probeType: Literal["multi_point"] = "multi_point"
    locations: tuple[Vec3, ...] = Field(min_length=1)

    def to_core(self) -> MultiPointProbeSpec:
        """Core probe spec."""
        return MultiPointProbeSpec(name=self.name, locations=self.locations, variables=self.variables)


def _core(obj: Any) -> Any:
    """Coreform facade object -> core spec (core objects pass through)."""
    return obj.to_core() if hasattr(obj, "to_core") else obj


# ----------------------------------------------------------------------------- job
class JobStatus(str, enum.Enum):
    """Job lifecycle states (Abaqus-style names)."""

    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    ABORTED = "ABORTED"


class Job:
    """Coreform/Abaqus-style job delegating to the WHT_IGA 5-layer core.

    Either pass a ready ``study=StudySpec`` or the study fields (``parts``,
    ``materials``, ``mesh``, ``steps``, ``probes``, ``sets``, ``units``); facade
    objects (``ImmersedRectilinearMesh``, ``PointProbe`` ...) are converted with
    ``to_core()``. ``backend`` defaults to ``KratosPythonBackend``.
    """

    def __init__(self, name: str | None = None, *, study: StudySpec | None = None, parts=None,
                 materials=None, mesh=None, steps=None, probes=(), sets=None,
                 units: str = "mm-N-tonne-s", output_dir: str | Path | None = None, backend: Any = None,
                 mesh_engine: Any = None, enable_solver_log: bool = True, on_progress: Any = None) -> None:
        if study is None:
            study = StudySpec(name=name, units=units, parts=tuple(parts or ()), materials=dict(materials or {}),
                              mesh=_core(mesh), steps=tuple(steps or ()),
                              probes=tuple(_core(p) for p in probes), sets=dict(sets or {}),
                              output_dir=str(output_dir) if output_dir is not None else None)
        elif name is not None and name != study.name:
            study = study.with_updates(name=name)
        self._study = study
        self._output_dir = output_dir
        self.backend = backend
        self.mesh_engine = mesh_engine
        self.enable_solver_log = enable_solver_log
        self.on_progress = on_progress
        self._plan = None
        self._session = None
        self._thread: threading.Thread | None = None
        self._status = JobStatus.CREATED
        self._results = None
        self._error: BaseException | None = None
        self._lock = threading.Lock()
        from wht_iga.ports.contracts import CancelToken

        self._cancel = CancelToken()  # shared with the session so kill() works before it opens

    # -- spec access (Coreform names)
    @property
    def name(self) -> str:
        """Job (study) name."""
        return self._study.name

    @property
    def spec(self) -> StudySpec:
        """Underlying immutable :class:`StudySpec`."""
        return self._study

    def to_spec(self) -> StudySpec:
        """Underlying immutable :class:`StudySpec`."""
        return self._study

    @classmethod
    def from_spec(cls, study: StudySpec, **kwargs: Any) -> Job:
        """Job from a core study spec."""
        return cls(study=study, **kwargs)

    def copy_with(self, **changes: Any) -> Job:
        """New job whose study has ``changes`` applied (facade objects accepted)."""
        changes = {k: (tuple(_core(x) for x in v) if k == "probes" else _core(v)) for k, v in changes.items()}
        return Job(study=self._study.with_updates(**changes), output_dir=self._output_dir, backend=self.backend,
                   mesh_engine=self.mesh_engine, enable_solver_log=self.enable_solver_log,
                   on_progress=self.on_progress)

    @property
    def plan(self):
        """Compiled :class:`wht_iga.plan.Plan` (compiled on first access)."""
        if self._plan is None:
            from wht_iga.plan import compile_plan

            self._plan = compile_plan(self._study, output_dir=self._output_dir)
        return self._plan

    # -- lifecycle
    @property
    def status(self) -> JobStatus:
        """Current :class:`JobStatus`."""
        return self._status

    @property
    def error(self) -> BaseException | None:
        """Exception of a FAILED / ABORTED job."""
        return self._error

    @property
    def results(self):
        """:class:`wht_iga.runner.StudyResults` of a COMPLETED job."""
        if self._status is not JobStatus.COMPLETED:
            raise RuntimeError(f"job {self.name!r} has no results (status {self._status.value})")
        return self._results

    def _open_session(self):
        if self._session is None:
            from wht_iga.runner import AnalysisSession

            self._session = AnalysisSession(self.plan, backend=self.backend, mesh_engine=self.mesh_engine,
                                            enable_solver_log=self.enable_solver_log,
                                            on_progress=self.on_progress, cancel_token=self._cancel)
            self._session.__enter__()
        return self._session

    def _close_session(self, exc: BaseException | None) -> None:
        if self._session is not None:
            self._session.__exit__(type(exc) if exc is not None else None, exc, None)

    def write_input(self) -> Path:
        """Write ``<job>.spec.json`` (input reproduction) and return its path."""
        p = self.plan.artifacts.spec_json
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(self._study.model_dump_json(indent=2), encoding="utf-8")
        return p

    def build_mesh(self):
        """Build the mesh now (else ``submit`` builds it); returns the ``MeshModel``."""
        with self._lock:
            if self._status is not JobStatus.CREATED:
                raise RuntimeError(f"build_mesh needs status CREATED, got {self._status.value}")
            self.plan.report.raise_if_invalid()
            sess = self._open_session()
        try:
            return sess.build_mesh()
        except BaseException as exc:
            self._finish(exc)
            raise

    def _finish(self, exc: BaseException | None) -> None:
        from wht_iga.ports.contracts import AnalysisCancelled

        self._close_session(exc)
        if exc is None:
            self._status = JobStatus.COMPLETED
        elif isinstance(exc, (AnalysisCancelled, KeyboardInterrupt)):
            self._status, self._error = JobStatus.ABORTED, exc
        else:
            self._status, self._error = JobStatus.FAILED, exc

    def _execute(self) -> None:
        self._status = JobStatus.RUNNING
        try:
            sess = self._open_session()
            if sess.mesh is None:
                sess.build_mesh()
            self._results = sess.run()
        except BaseException as exc:  # recorded in status/error; wait_for_completion re-raises nothing
            self._finish(exc)
            return
        self._finish(None)

    def submit(self, *, wait: bool = False) -> None:
        """Start the analysis (background thread; ``wait=True`` runs in the calling thread)."""
        with self._lock:
            if self._status is not JobStatus.CREATED:
                raise RuntimeError(f"submit needs status CREATED, got {self._status.value}")
            self.plan.report.raise_if_invalid()
            self._status = JobStatus.SUBMITTED
        if wait:
            self._execute()
            return
        self._thread = threading.Thread(target=self._execute, name=f"wht-iga-job-{self.name}", daemon=True)
        self._thread.start()

    def wait_for_completion(self, timeout: float | None = None) -> JobStatus:
        """Block until the job finishes (or ``timeout`` seconds); returns the status."""
        if self._thread is not None:
            self._thread.join(timeout)
        return self._status

    def kill(self, reason: str = "killed by user") -> None:
        """Cooperative abort: takes effect at the next stage boundary; partial artifacts are removed."""
        if self._status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.ABORTED):
            return
        self._cancel.cancel(reason)
        if self._status is JobStatus.CREATED:
            from wht_iga.ports.contracts import AnalysisCancelled

            self._finish(AnalysisCancelled(reason))
