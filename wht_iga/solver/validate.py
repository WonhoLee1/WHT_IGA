"""Job spec + validation gate for Phase-1 linear runs (Track SOLVE).

``validate()`` BLOCKS the run on error-severity issues (watertight failure,
volume-sign error, QUAD_FALLBACK on load-bearing faces, bad material, no
Dirichlet support). Track GEO issue codes are consumed as plain dicts/JSON
(``{"severity","code","message","entity"}``) — no import of geometry modules.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


LinearSolverType = Literal[
    "auto",
    "pardiso",
    "pardiso_ldlt",
    "pardiso_lu",
    "pardiso_llt",
    "amgcl",
    "amgcl_cg",
    "cg",
    "splu",
    "skyline_lu",
]

_VALID_LINEAR_SOLVERS: tuple[str, ...] = (
    "auto",
    "pardiso",
    "pardiso_ldlt",
    "pardiso_lu",
    "pardiso_llt",
    "amgcl",
    "amgcl_cg",
    "cg",
    "splu",
    "skyline_lu",
)

_VALID_LOG_LEVELS: tuple[str, ...] = (
    "summary",
    "detailed",
    "debug",
    "minimal",
    "standard",
    "verbose",
)


@dataclass
class ValidationIssue:
    severity: str  # "error" | "warning" | "info"
    code: str
    message: str
    entity: str = ""


@dataclass
class AffineField:
    """JSON-safe affine Dirichlet field g(x) = A @ x + b."""

    A: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    b: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass
class DirichletBC:
    face: tuple[int, int]  # (axis, side); side 0 = lo, 1 = hi
    components: tuple[bool, bool, bool] = (True, True, True)
    value: tuple[float, float, float] = (0.0, 0.0, 0.0)
    field: AffineField | None = None  # affine g(x); overrides value when set


@dataclass
class PressureLoad:
    face: tuple[int, int]
    magnitude: float
    follower: bool = False


@dataclass
class TractionLoad:
    face: tuple[int, int]
    vector: tuple[float, float, float]


@dataclass
class JobSpec:
    """Portable Phase-1 linear job spec (v2 §8 ``to_json/from_json/validate``)."""

    name: str
    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    n_cells: tuple[int, int, int]
    degree: int
    continuity: int | None = None
    E: float = 1.0
    nu: float = 0.0
    density: float = 0.0
    beta: float | None = None  # None -> default_beta(E, h, p)
    bcs: list[DirichletBC] = field(default_factory=list)
    fixed_faces: list[tuple[int, int]] = field(default_factory=list)  # strong Fix
    pressures: list[PressureLoad] = field(default_factory=list)
    tractions: list[TractionLoad] = field(default_factory=list)
    gravity: tuple[float, float, float] | None = None
    out_cells: list[tuple[int, int, int]] = field(default_factory=list)
    ghost_alpha: float = 1e-8  # FCM OUT-region stiffness scale (v2 §7.6)
    # OUT cells: "ghost" = FCM fill at ghost_alpha*E (default, V02 evidence);
    # "drop" = no OUT elements, inactive basis fixed (v2.1 §7.5, wht_iga.solver.active)
    out_mode: str = "ghost"
    probes: list[dict] = field(default_factory=list)
    assumed: list[str] = field(default_factory=list)  # M8: defaulted inputs
    output_dir: str | None = None
    # Kratos linear solver: "pardiso_lu" (MKL, default) | "pardiso_ldlt" | "pardiso_llt"
    # (symmetric K only; no follower loads) | "skyline_lu" (small/debug only: profile
    # storage explodes on 3D grids — V02 level 4 hit 16 GB, 2026-11-14)
    # "auto" (default): pardiso_ldlt (half the factor memory of LU — matters on a
    # 32 GB machine) unless a follower load makes the tangent unsymmetric, then pardiso_lu.
    linear_solver: LinearSolverType | str = "auto"
    num_threads_mkl: int | str | None = None  # None = MKL default; int; "consistent"
    enable_solver_log: bool = True
    solver_log_file: str | None = None
    solver_log_stdout: bool = False
    solver_log_level: Literal["summary", "detailed", "debug", "minimal", "standard", "verbose"] = "detailed"
    amgcl_params: dict[str, Any] | None = None
    memory_budget_ratio: float = 0.70

    def to_json(self) -> dict:
        data = asdict(self)
        if data.get("solver_log_file") is not None:
            data["solver_log_file"] = str(data["solver_log_file"])
        return data

    @staticmethod
    def from_json(data: dict) -> JobSpec:
        def tup(x: list, n: int = 3) -> tuple:
            return tuple(float(v) for v in x) if len(x) == n else tuple(x)

        bcs = [
            DirichletBC(
                face=tuple(b["face"]),  # type: ignore[arg-type]
                components=tuple(b.get("components", (True, True, True))),  # type: ignore[arg-type]
                value=tup(b.get("value", (0.0, 0.0, 0.0))),
                field=AffineField(A=tuple(tuple(r) for r in b["field"]["A"]), b=tup(b["field"].get("b", (0.0, 0.0, 0.0))))  # type: ignore[arg-type]
                if b.get("field") else None,
            )
            for b in data.get("bcs", [])
        ]
        return JobSpec(
            name=data["name"],
            lo=tup(data["lo"]),
            hi=tup(data["hi"]),
            n_cells=tuple(int(v) for v in data["n_cells"]),  # type: ignore[arg-type]
            degree=int(data["degree"]),
            continuity=data.get("continuity"),
            E=float(data["E"]),
            nu=float(data["nu"]),
            density=float(data.get("density", 0.0)),
            beta=data.get("beta"),
            bcs=bcs,
            fixed_faces=[tuple(f) for f in data.get("fixed_faces", [])],  # type: ignore[misc]
            pressures=[
                PressureLoad(face=tuple(p["face"]), magnitude=float(p["magnitude"]),  # type: ignore[arg-type]
                             follower=bool(p.get("follower", False)))
                for p in data.get("pressures", [])
            ],
            tractions=[
                TractionLoad(face=tuple(t["face"]), vector=tup(t["vector"]))  # type: ignore[arg-type]
                for t in data.get("tractions", [])
            ],
            gravity=tup(data["gravity"]) if data.get("gravity") else None,
            out_cells=[tuple(c) for c in data.get("out_cells", [])],  # type: ignore[misc]
            ghost_alpha=float(data.get("ghost_alpha", 1e-8)),
            out_mode=str(data.get("out_mode", "ghost")),
            probes=list(data.get("probes", [])),
            assumed=list(data.get("assumed", [])),
            linear_solver=str(data.get("linear_solver", "auto")),
            num_threads_mkl=data.get("num_threads_mkl"),
            output_dir=data.get("output_dir"),
            enable_solver_log=bool(data.get("enable_solver_log", True)),
            solver_log_file=str(data["solver_log_file"]) if data.get("solver_log_file") is not None else None,
            solver_log_stdout=bool(data.get("solver_log_stdout", False)),
            solver_log_level=str(data.get("solver_log_level", "detailed")),  # type: ignore[arg-type]
            amgcl_params=dict(data["amgcl_params"]) if data.get("amgcl_params") is not None else None,
            memory_budget_ratio=float(data.get("memory_budget_ratio", 0.70)),
        )


class ValidationError(Exception):
    """Raised by ``run_linear`` when ``validate()`` reports error issues."""

    def __init__(self, issues: list[ValidationIssue]) -> None:
        self.issues = [i for i in issues if i.severity == "error"]
        super().__init__(f"{len(self.issues)} blocking validation errors: " + "; ".join(f"[{i.code}] {i.message}" for i in self.issues))


def validate(spec: JobSpec, geo_issues: list[dict] | None = None) -> list[ValidationIssue]:
    """Check a job spec + Track GEO mesh-stage issues. Errors block the run."""
    out: list[ValidationIssue] = []
    if not spec.E > 0.0:
        out.append(ValidationIssue("error", "MAT_E", f"Young's modulus must be > 0, got {spec.E}", "material"))
    if not -1.0 < spec.nu < 0.5:
        out.append(ValidationIssue("error", "MAT_NU", f"Poisson ratio must be in (-1, 0.5), got {spec.nu}", "material"))
    if spec.degree not in (1, 2, 3):
        out.append(ValidationIssue("error", "GRID_DEGREE", f"degree must be 1..3, got {spec.degree}", "grid"))
    if any(n < 1 for n in spec.n_cells):
        out.append(ValidationIssue("error", "GRID_CELLS", f"n_cells must be >= 1, got {spec.n_cells}", "grid"))
    if any(h <= l for h, l in zip(spec.hi, spec.lo)):
        out.append(ValidationIssue("error", "BOX_SHAPE", f"box hi must exceed lo, got lo={spec.lo} hi={spec.hi}", "grid"))
    if spec.beta is not None and not spec.beta > 0.0:
        out.append(ValidationIssue("error", "BETA_VALUE", f"penalty beta must be > 0, got {spec.beta}", "bc"))
    if not spec.bcs and not spec.fixed_faces:
        out.append(ValidationIssue("error", "NO_DIRICHLET", "no Dirichlet support: rigid-body singular", "bc"))
    elif spec.bcs and not spec.fixed_faces and all(not any(b.components) for b in spec.bcs):
        out.append(ValidationIssue("error", "NO_DIRICHLET", "all Dirichlet BCs have empty components: rigid-body singular", "bc"))
    if not 0.0 < spec.ghost_alpha <= 1.0:
        out.append(ValidationIssue("error", "GHOST_ALPHA", f"ghost_alpha must be in (0, 1], got {spec.ghost_alpha}", "fcm"))
    if spec.out_mode not in ("ghost", "drop"):
        out.append(ValidationIssue("error", "OUT_MODE", f"out_mode must be 'ghost' or 'drop', got {spec.out_mode!r}", "fcm"))
    if spec.linear_solver not in _VALID_LINEAR_SOLVERS:
        out.append(
            ValidationIssue(
                "error",
                "LINEAR_SOLVER",
                f"linear_solver must be one of {_VALID_LINEAR_SOLVERS}, got {spec.linear_solver!r}",
                "solver",
            )
        )
    if str(spec.solver_log_level).lower() not in _VALID_LOG_LEVELS:
        out.append(
            ValidationIssue(
                "error",
                "SOLVER_LOG_LEVEL",
                f"solver_log_level must be one of {_VALID_LOG_LEVELS}, got {spec.solver_log_level!r}",
                "solver",
            )
        )
    if not 0.0 < float(spec.memory_budget_ratio) <= 1.0:
        out.append(
            ValidationIssue(
                "error",
                "MEMORY_BUDGET_RATIO",
                f"memory_budget_ratio must be in (0, 1], got {spec.memory_budget_ratio}",
                "solver",
            )
        )
    nx, ny, nz = spec.n_cells
    for c in spec.out_cells:
        if not (0 <= c[0] < nx and 0 <= c[1] < ny and 0 <= c[2] < nz):
            out.append(ValidationIssue("error", "OUT_CELL_RANGE", f"out cell {c} outside grid {spec.n_cells}", "fcm"))
    if spec.gravity is not None and not spec.density > 0.0:
        out.append(ValidationIssue("warning", "GRAVITY_NO_DENSITY", "gravity set but density = 0: no body force", "load"))
    loaded = {tuple(p.face) for p in spec.pressures} | {tuple(t.face) for t in spec.tractions}
    for g in geo_issues or []:
        code = str(g.get("code", ""))
        sev = str(g.get("severity", "info"))
        msg = str(g.get("message", ""))
        ent = str(g.get("entity", ""))
        if code in ("WATERTIGHT", "VOLUME_SIGN"):
            out.append(ValidationIssue("error", code, msg, ent))
        elif code == "QUAD_FALLBACK":
            ent_face = _entity_face(ent)
            ent_cell = _entity_cell(ent)
            on_loaded = ent_face is not None and ent_face in loaded
            if ent_cell is not None and not on_loaded:
                on_loaded = any(_cell_touches_face(ent_cell, spec.n_cells, f) for f in loaded)
            if on_loaded:
                where = ent_face if ent_face is not None else f"cell {ent_cell} adjacent to loaded face"
                out.append(ValidationIssue("error", "QUAD_FALLBACK_LOADED",
                                           f"fallback quadrature on load-bearing {where}: {msg}", ent))
            else:
                out.append(ValidationIssue("warning", code, msg, ent))
        elif sev == "error":
            out.append(ValidationIssue("error", code or "GEO_ERROR", msg, ent))
        else:
            out.append(ValidationIssue(sev if sev in ("warning", "info") else "info", code or "GEO_INFO", msg, ent))
    for a in spec.assumed:
        out.append(ValidationIssue("info", "ASSUMED", f"defaulted input recorded as assumed: {a}", "spec"))
    return out


def _entity_face(entity: str) -> tuple[int, int] | None:
    """Parse ``face:<axis>:<side>`` entities; None when not a face entity."""
    parts = entity.split(":")
    if len(parts) == 3 and parts[0] == "face":
        try:
            return (int(parts[1]), int(parts[2]))
        except ValueError:
            return None
    return None


def _entity_cell(entity: str) -> tuple[int, int, int] | None:
    """Parse ``cell:<i>,<j>,<k>`` entities; None otherwise."""
    if not entity.startswith("cell:"):
        return None
    try:
        i, j, k = (int(v) for v in entity[len("cell:"):].split(","))
        return (i, j, k)
    except ValueError:
        return None


def _cell_touches_face(cell: tuple[int, int, int], n_cells: tuple[int, int, int],
                       face: tuple[int, int]) -> bool:
    """True when a cell touches a box face (for fallback-on-load detection)."""
    axis, side = face
    extent = (n_cells[0], n_cells[1], n_cells[2])[axis]
    c = (cell[0], cell[1], cell[2])[axis]
    return (c == 0 and side == 0) or (c == extent - 1 and side == 1)
