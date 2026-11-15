"""Layer 2: pre-flight plan (W3). ``compile_plan(spec)`` -> :class:`Plan`.

Resolves everything that can be decided before any solver runs and reports it
in a :class:`ResolveReport`:

* background grid (bounds, padding, ``n_cells``, degree/continuity, trimmed or not),
* named sets -> grid faces, and whether each BC/load face really lies on the part,
* control points / DOFs / estimated non-zeros,
* host RAM snapshot and the per-step linear-solver prediction (the same
  ``select_linear_solver`` decision and ``auto`` mapping ``run_linear`` uses),
* unsupported combinations as blocking ``error`` diagnostics (never silently ignored).

No Kratos import; STEP geometry is read only to get its bounding box.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

from wht_iga.plan.artifacts import ArtifactPaths
from wht_iga.spec.assembly import BoxGeometrySpec, StepGeometrySpec
from wht_iga.spec.mesh import BodyFitSpec, ImmersedGridSpec
from wht_iga.spec.physics import GravityLoadSpec, PressureLoadSpec, TractionLoadSpec
from wht_iga.spec.study import PROBE_VARIABLE_CANONICAL, StudySpec

Severity = Literal["error", "warning", "info"]

#: probe quantities the session can evaluate at points from the control-point solution
COMPUTABLE_PROBE_VARIABLES = frozenset({"U", "S", "E", "MISES"})
_SYMMETRIC_ONLY = ("pardiso", "pardiso_ldlt", "pardiso_llt", "amgcl", "amgcl_cg", "cg")
_ITERATIVE = ("amgcl", "amgcl_cg", "cg")


@dataclass(frozen=True)
class Diagnostic:
    """One pre-flight finding."""

    severity: Severity
    code: str
    message: str
    entity: str = ""


@dataclass(frozen=True)
class GridResolution:
    """Background grid decided by the plan."""

    kind: str  # "immersed" | "body_fit"
    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    n_cells: tuple[int, int, int]
    h: tuple[float, float, float]
    degree: int
    continuity: int
    padding: tuple[int, int, int]
    trimmed: bool
    quadrature_degree: int
    tessellation_size: float | None
    part_lo: tuple[float, float, float]
    part_hi: tuple[float, float, float]


@dataclass(frozen=True)
class StepPrediction:
    """Linear-solver prediction for one step (mirrors ``run_linear``)."""

    name: str
    requested_solver: str
    predicted_solver: str  # memory-aware decision ("pardiso"/"amgcl"/... or the explicit name)
    kratos_solver: str  # Kratos solver name run_linear would build
    switched: bool
    switched_reason: str | None
    est_pardiso_bytes: int
    est_amgcl_bytes: int
    est_peak_bytes: int
    budget_bytes: int


@dataclass(frozen=True)
class ResolveReport:
    """Everything the plan resolved plus its diagnostics."""

    diagnostics: tuple[Diagnostic, ...]
    grid: GridResolution | None = None
    n_ctrl: int = 0
    total_dofs: int = 0
    est_active_dofs: int = 0  # upper bound when trimmed (no classification before build_mesh)
    est_nnz: int = 0
    total_ram_bytes: int = 0
    avail_ram_bytes: int = 0
    steps: tuple[StepPrediction, ...] = ()
    resolved_sets: dict[str, tuple[int, int]] = field(default_factory=dict)
    assumed: tuple[str, ...] = ()
    dt_crit: float | None = None  # explicit steps only (none in the static-only spec today)

    @property
    def is_valid(self) -> bool:
        """True when no error-severity diagnostic exists."""
        return not self.errors()

    def errors(self) -> tuple[Diagnostic, ...]:
        """Blocking diagnostics."""
        return tuple(d for d in self.diagnostics if d.severity == "error")

    def warnings(self) -> tuple[Diagnostic, ...]:
        """Warning diagnostics."""
        return tuple(d for d in self.diagnostics if d.severity == "warning")

    def codes(self) -> tuple[str, ...]:
        """All diagnostic codes in order."""
        return tuple(d.code for d in self.diagnostics)

    def raise_if_invalid(self) -> None:
        """Raise :class:`PlanError` listing the errors (no-op when valid)."""
        if not self.is_valid:
            raise PlanError(self)

    def summary(self) -> str:
        """Human-readable multi-line summary."""
        g = self.grid
        lines = [f"ResolveReport: valid={self.is_valid} errors={len(self.errors())} "
                 f"warnings={len(self.warnings())}"]
        if g is not None:
            lines.append(f"  grid: {g.kind} n_cells={g.n_cells} p={g.degree} c={g.continuity} "
                         f"trimmed={g.trimmed} lo={g.lo} hi={g.hi}")
        lines.append(f"  dofs: total={self.total_dofs} est_active<={self.est_active_dofs} "
                     f"est_nnz={self.est_nnz}")
        gb = 1024.0 ** 3
        lines.append(f"  RAM: total={self.total_ram_bytes / gb:.2f} GB "
                     f"avail={self.avail_ram_bytes / gb:.2f} GB")
        for s in self.steps:
            lines.append(f"  step {s.name}: {s.requested_solver} -> {s.kratos_solver} "
                         f"(peak~{s.est_peak_bytes / gb:.3f} GB, budget {s.budget_bytes / gb:.3f} GB"
                         f"{', switched' if s.switched else ''})")
        for d in self.diagnostics:
            lines.append(f"  [{d.severity}] {d.code}: {d.message}")
        return "\n".join(lines)


class PlanError(RuntimeError):
    """A plan with blocking diagnostics was used where a valid plan is required."""

    def __init__(self, report: ResolveReport) -> None:
        self.report = report
        msgs = "; ".join(f"[{d.code}] {d.message}" for d in report.errors())
        super().__init__(f"plan has {len(report.errors())} blocking error(s): {msgs}")


@dataclass(frozen=True)
class Plan:
    """Compiled study: spec + resolve report + artifact paths (``ExecutionPlan``)."""

    spec: StudySpec
    report: ResolveReport
    artifacts: ArtifactPaths

    def job_spec(self, step_name: str, *, out_cells=(), output_dir: str | Path | None = None,
                 enable_solver_log: bool = True, solver_log_stdout: bool = False,
                 solver_log_level: str = "detailed"):
        """``wht_iga.solver.validate.JobSpec`` for one step (input of ``run_linear``)."""
        from wht_iga.solver.validate import DirichletBC, JobSpec, PressureLoad, TractionLoad

        self.report.raise_if_invalid()
        g = self.report.grid
        assert g is not None  # valid report always has a grid
        spec = self.spec
        step = spec.step(step_name)
        part = spec.parts[0]
        mat = spec.material_of(part)
        mesh = spec.mesh
        bcs = [DirichletBC(face=spec.set_face(b.set), components=b.components, value=b.value)
               for b in step.bcs]
        pressures = [PressureLoad(face=spec.set_face(ld.set), magnitude=float(ld.magnitude),
                                  follower=bool(ld.follower))
                     for ld in step.loads if isinstance(ld, PressureLoadSpec)]
        tractions = [TractionLoad(face=spec.set_face(ld.set), vector=tuple(ld.vector))
                     for ld in step.loads if isinstance(ld, TractionLoadSpec)]
        grav = next((tuple(ld.acceleration) for ld in step.loads if isinstance(ld, GravityLoadSpec)), None)
        arts = self.artifacts.step(step_name)
        out_dir = str(output_dir) if output_dir is not None else str(arts.output_dir)
        return JobSpec(
            name=arts.prefix, lo=g.lo, hi=g.hi, n_cells=g.n_cells, degree=g.degree,
            continuity=g.continuity, E=float(mat.E), nu=float(mat.nu), density=float(mat.density),
            beta=mesh.bc_penalty, bcs=bcs, pressures=pressures, tractions=tractions, gravity=grav,
            out_cells=[tuple(int(v) for v in c) for c in out_cells],
            ghost_alpha=float(mesh.ghost_alpha) if isinstance(mesh, ImmersedGridSpec) else 1e-8,
            out_mode=mesh.out_mode if isinstance(mesh, ImmersedGridSpec) else "drop",
            probes=[], assumed=list(self.report.assumed), output_dir=out_dir,
            linear_solver=step.linear_solver, num_threads_mkl=step.num_threads_mkl,
            enable_solver_log=enable_solver_log, solver_log_file=str(arts.solver_log),
            solver_log_stdout=solver_log_stdout, solver_log_level=solver_log_level,  # type: ignore[arg-type]
            memory_budget_ratio=float(step.memory_budget_ratio),
        )

    def extension(self) -> dict | None:
        """Small-cell extension settings for ``run_linear`` (None = off / untrimmed)."""
        g = self.report.grid
        mesh = self.spec.mesh
        if g is None or not g.trimmed or not isinstance(mesh, ImmersedGridSpec):
            return None
        if mesh.small_cell_volume_ratio <= 0.0:
            return None
        return {"theta": float(mesh.small_cell_volume_ratio)}


ExecutionPlan = Plan


def _read_step_bounds(geom: StepGeometrySpec, diags: list[Diagnostic]):
    """Bounding box of a STEP part via the geometry kernel (OCP), or None + error."""
    path = Path(geom.path)
    if not path.is_file():
        diags.append(Diagnostic("error", "PART_FILE_MISSING", f"STEP file not found: {path}", str(path)))
        return None
    try:
        from wht_iga.geometry import cad, occ_io

        shape = occ_io.load_cad_file(path, heal=geom.heal).shape
        lo, hi = cad.shape_bbox(shape)
    except Exception as exc:  # OCP missing or unreadable file: report, never guess bounds
        diags.append(Diagnostic("error", "GEOMETRY_UNREADABLE", f"cannot read {path}: {exc}", str(path)))
        return None
    return tuple(float(v) for v in lo), tuple(float(v) for v in hi)


def _resolve_grid(spec: StudySpec, diags: list[Diagnostic], assumed: list[str]) -> GridResolution | None:
    """Decide the background grid from the part bounds and the mesh spec."""
    part = spec.parts[0]
    mesh = spec.mesh
    geom = part.geometry
    if isinstance(mesh, BodyFitSpec) and not isinstance(geom, BoxGeometrySpec):
        diags.append(Diagnostic("error", "BODYFIT_REQUIRES_BOX",
                                "body-fit grids need an analytic box part (use ImmersedGridSpec for CAD)",
                                part.name))
        return None
    if isinstance(mesh, ImmersedGridSpec) and mesh.bounds is not None:
        plo, phi = tuple(mesh.bounds[0]), tuple(mesh.bounds[1])
        if isinstance(geom, BoxGeometrySpec):
            plo_b, phi_b = tuple(geom.lo), tuple(geom.hi)
        else:
            plo_b, phi_b = plo, phi
    elif isinstance(geom, BoxGeometrySpec):
        plo, phi = tuple(geom.lo), tuple(geom.hi)
        plo_b, phi_b = plo, phi
    else:
        bb = _read_step_bounds(geom, diags)
        if bb is None:
            return None
        plo, phi = bb
        plo_b, phi_b = plo, phi
    lo = np.asarray(plo, float)
    ext = np.asarray(phi, float) - lo
    if isinstance(mesh, BodyFitSpec):
        if mesh.element_size is not None:
            n = np.maximum(1, np.rint(ext / mesh.element_size)).astype(int)
        else:
            n = np.asarray(mesh.n_cells, int)
        h = ext / n
        hi = lo + ext
        pad = (0, 0, 0)
        trimmed = False
        q = mesh.degree + 1
        tess = None
    else:
        if mesh.element_size is not None:
            n = np.maximum(1, np.ceil(ext / mesh.element_size - 1e-9)).astype(int)
            h = np.full(3, float(mesh.element_size))
        else:
            n = np.asarray(mesh.n_cells, int)
            h = ext / n
        hi = lo + n * h
        pad = tuple(int(v) for v in mesh.padding)
        lo = lo - np.asarray(pad) * h
        hi = hi + np.asarray(pad) * h
        n = n + 2 * np.asarray(pad)
        tol = 1e-9 * float(np.max(np.abs(np.concatenate([lo, hi]))) + np.max(h))
        exact_box = (isinstance(geom, BoxGeometrySpec)
                     and np.allclose(lo, plo_b, atol=tol, rtol=0) and np.allclose(hi, phi_b, atol=tol, rtol=0))
        trimmed = not exact_box
        q = mesh.resolved_quadrature_degree
        if mesh.quadrature_degree is None and trimmed:
            assumed.append(f"quadrature_degree={q} (min(2p, 4))")
        if exact_box:
            diags.append(Diagnostic("info", "IMMERSED_EXACT_FIT",
                                    "grid coincides with the box part: no trimmed cells", part.name))
        tess = mesh.tessellation_size
        if trimmed and isinstance(geom, StepGeometrySpec) and tess is None:
            tess = float(np.min(h)) / 4.0
            assumed.append(f"tessellation_size={tess:g} (min(h)/4)")
    if mesh.continuity is None:
        assumed.append(f"continuity={mesh.resolved_continuity} (degree-1)")
    if mesh.bc_penalty is None:
        assumed.append("bc_penalty=c_pen*E*p^2/h (default_beta)")
    return GridResolution(
        kind=mesh.kind, lo=tuple(float(v) for v in lo), hi=tuple(float(v) for v in hi),
        n_cells=tuple(int(v) for v in n), h=tuple(float(v) for v in h), degree=int(mesh.degree),
        continuity=int(mesh.resolved_continuity), padding=pad, trimmed=bool(trimmed),
        quadrature_degree=int(q), tessellation_size=tess,
        part_lo=tuple(float(v) for v in plo_b), part_hi=tuple(float(v) for v in phi_b))


def _face_checks(spec: StudySpec, g: GridResolution, diags: list[Diagnostic]) -> dict[str, tuple[int, int]]:
    """Resolve referenced sets to grid faces and check they lie on (and cover) the part."""
    part = spec.parts[0]
    is_box = isinstance(part.geometry, BoxGeometrySpec)
    resolved: dict[str, tuple[int, int]] = {}
    scale = max(max(abs(v) for v in g.lo + g.hi), max(g.h))
    tol = 1e-9 * scale

    def plane_on_part(axis: int, side: int) -> bool:
        grid_x = (g.lo, g.hi)[side][axis]
        part_x = (g.part_lo, g.part_hi)[side][axis]
        return abs(grid_x - part_x) <= tol

    for step in spec.steps:
        targets = [("bc", b.set) for b in step.bcs]
        targets += [(ld.kind, ld.set) for ld in step.loads if hasattr(ld, "set")]
        for what, name in targets:
            axis, side = spec.set_face(name)
            resolved[name] = (axis, side)
            ent = f"set:{name}"
            if g.padding[axis] > 0:
                diags.append(Diagnostic("error", "SET_OFF_DOMAIN",
                                        f"step {step.name!r}: {what} on set {name!r} = grid face "
                                        f"(axis={axis}, side={side}) lies in padding (padding[{axis}]="
                                        f"{g.padding[axis]}); boundary data there would act on void", ent))
                continue
            if is_box and not plane_on_part(axis, side):
                diags.append(Diagnostic("error", "SET_OFF_DOMAIN",
                                        f"step {step.name!r}: {what} on set {name!r}: grid face "
                                        f"(axis={axis}, side={side}) does not coincide with the part face "
                                        f"(element_size does not divide the part length?)", ent))
                continue
            if what in ("pressure", "traction"):
                if is_box:
                    others = [d for d in range(3) if d != axis]
                    covered = all(plane_on_part(d, 0) and plane_on_part(d, 1) for d in others)
                    if not covered:
                        diags.append(Diagnostic("error", "LOAD_FACE_NOT_COVERED",
                                                f"step {step.name!r}: {what} on {name!r}: the part does not "
                                                f"cover the whole grid face; run_linear integrates loads "
                                                f"over the full face", ent))
                else:
                    diags.append(Diagnostic("warning", "LOAD_FACE_COVERAGE_UNVERIFIED",
                                            f"step {step.name!r}: {what} on {name!r}: face coverage of a CAD "
                                            f"part is checked after classification (build_mesh)", ent))
            if what == "bc" and g.trimmed and not is_box:
                diags.append(Diagnostic("info", "BC_WHOLE_GRID_FACE",
                                        f"BC on {name!r} acts on the whole grid-face plane "
                                        f"(repo convention, V02/V06)", ent))
    return resolved


def _overlap_cp(degree: int, n_cells: tuple[int, int, int]) -> int:
    """Control points coupled to one control point (same formula as run_linear)."""
    nx, ny, nz = n_cells
    return (min(2 * degree + 1, nx + degree) * min(2 * degree + 1, ny + degree)
            * min(2 * degree + 1, nz + degree))


def compile_plan(spec: StudySpec | dict | str, *, output_dir: str | Path | None = None,
                 avail_bytes: int | None = None) -> Plan:
    """Compile a study into a :class:`Plan` (never raises for invalid physics: see report).

    ``spec`` may be a :class:`StudySpec`, a dict or a JSON string. ``output_dir``
    overrides ``spec.output_dir`` (default ``./wiga_runs/<name>``). ``avail_bytes``
    overrides the measured available RAM (tests / what-if planning).
    """
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.memory import probe_system_memory, select_linear_solver

    if isinstance(spec, str):
        spec = StudySpec.model_validate_json(spec)
    elif isinstance(spec, dict):
        spec = StudySpec.model_validate(spec)
    diags: list[Diagnostic] = []
    assumed: list[str] = []
    out = Path(output_dir) if output_dir is not None else (
        Path(spec.output_dir) if spec.output_dir else Path.cwd() / "wiga_runs" / spec.name)
    artifacts = ArtifactPaths(out.resolve(), spec.name, tuple(s.name for s in spec.steps))

    if len(spec.parts) > 1:
        diags.append(Diagnostic("error", "UNSUPPORTED_MULTI_PART",
                                f"{len(spec.parts)} parts: the solver pipeline takes one part "
                                f"(assembly/tie support is not wired into run_linear)", "parts"))
    for r in spec.mesh.refinements:
        diags.append(Diagnostic("error", "UNSUPPORTED_REFINEMENT",
                                f"refinement {r.name!r}: THB refinement is not wired into the solver "
                                f"pipeline yet; remove it or run without refinement", f"refinement:{r.name}"))
    for p in spec.probes:
        bad = sorted({PROBE_VARIABLE_CANONICAL[v] for v in p.variables} - COMPUTABLE_PROBE_VARIABLES)
        if bad:
            diags.append(Diagnostic("error", "PROBE_VARIABLE_UNSUPPORTED",
                                    f"probe {p.name!r}: {bad} cannot be evaluated at points in a "
                                    f"static structural study (supported: U, S, E, MISES)", f"probe:{p.name}"))
    for s in spec.steps:
        follower = any(isinstance(ld, PressureLoadSpec) and ld.follower for ld in s.loads)
        if follower and s.linear_solver in _SYMMETRIC_ONLY:
            diags.append(Diagnostic("error", "FOLLOWER_NEEDS_LU",
                                    f"step {s.name!r}: follower pressure gives an unsymmetric tangent; "
                                    f"linear_solver={s.linear_solver!r} needs a symmetric matrix "
                                    f"(use 'pardiso_lu' or 'auto')", f"step:{s.name}"))
        if any(isinstance(ld, GravityLoadSpec) for ld in s.loads):
            if spec.material_of(spec.parts[0]).density <= 0.0:
                diags.append(Diagnostic("warning", "GRAVITY_NO_DENSITY",
                                        f"step {s.name!r}: gravity with density 0 gives no body force",
                                        f"step:{s.name}"))

    grid = _resolve_grid(spec, diags, assumed)
    if grid is None:
        return Plan(spec, ResolveReport(tuple(diags), assumed=tuple(assumed)), artifacts)
    if isinstance(spec.mesh, ImmersedGridSpec) and spec.mesh.tessellation_size is not None \
            and not (grid.trimmed and isinstance(spec.parts[0].geometry, StepGeometrySpec)):
        diags.append(Diagnostic("info", "TESSELLATION_SIZE_UNUSED",
                                "tessellation_size only applies to trimmed CAD parts", "mesh"))
    resolved = _face_checks(spec, grid, diags)

    box = BoxSpline(lo=grid.lo, hi=grid.hi, n_cells=grid.n_cells, degree=grid.degree,
                    continuity=grid.continuity)
    n_ctrl = int(box.n_ctrl)
    total_dofs = 3 * n_ctrl
    est_active = total_dofs  # untrimmed: exact; trimmed: upper bound until classification
    est_nnz = int(est_active * 3 * _overlap_cp(grid.degree, grid.n_cells))
    if grid.trimmed:
        diags.append(Diagnostic("info", "ACTIVE_DOFS_UPPER_BOUND",
                                f"trimmed grid: memory estimate uses all {total_dofs} DOFs (upper bound); "
                                f"inactive basis functions are known only after build_mesh", "dofs"))
    mem = probe_system_memory()
    avail = int(avail_bytes) if avail_bytes is not None else int(mem.available_bytes)
    preds: list[StepPrediction] = []
    for s in spec.steps:
        move_mesh = any(isinstance(ld, PressureLoadSpec) and ld.follower for ld in s.loads)
        dec = select_linear_solver(s.linear_solver, est_active, est_nnz, degree=grid.degree,
                                   memory_budget_ratio=float(s.memory_budget_ratio), avail_bytes=avail,
                                   pypardiso_available=True)
        if s.linear_solver == "auto":
            kratos = "amgcl_cg" if (dec.switched and not move_mesh) else (
                "pardiso_lu" if move_mesh else "pardiso_ldlt")
        else:
            kratos = s.linear_solver
        peak = dec.est_amgcl_bytes if kratos in _ITERATIVE else dec.est_pardiso_bytes
        if s.linear_solver != "auto" and kratos not in _ITERATIVE and dec.est_pardiso_bytes > dec.budget_bytes:
            diags.append(Diagnostic("warning", "MEMORY_OVER_BUDGET",
                                    f"step {s.name!r}: {kratos} needs ~{dec.est_pardiso_gb:.2f} GB > budget "
                                    f"{dec.budget_bytes / 1024.0 ** 3:.2f} GB; consider linear_solver='auto'",
                                    f"step:{s.name}"))
        if dec.switched:
            diags.append(Diagnostic("info", "AUTO_SOLVER_SWITCH", str(dec.switched_reason), f"step:{s.name}"))
        preds.append(StepPrediction(
            name=s.name, requested_solver=s.linear_solver, predicted_solver=dec.selected_solver,
            kratos_solver=kratos, switched=bool(dec.switched), switched_reason=dec.switched_reason,
            est_pardiso_bytes=int(dec.est_pardiso_bytes), est_amgcl_bytes=int(dec.est_amgcl_bytes),
            est_peak_bytes=int(peak), budget_bytes=int(dec.budget_bytes)))
    report = ResolveReport(
        diagnostics=tuple(diags), grid=grid, n_ctrl=n_ctrl, total_dofs=total_dofs,
        est_active_dofs=est_active, est_nnz=est_nnz, total_ram_bytes=int(mem.total_bytes),
        avail_ram_bytes=avail, steps=tuple(preds), resolved_sets=resolved, assumed=tuple(assumed))
    return Plan(spec, report, artifacts)
