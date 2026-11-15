"""Physical-space PointProbe, LineProbe, and MultiPointProbe evaluators and exporters.

Evaluates displacement (``U``), strain (``E``), Cauchy stress (``S``), von Mises
stress (``MISES``), and temperature (``NT11``) using either direct B-spline/THB
basis evaluation or Kratos zero-weight (``w = 0.0``) quadrature elements so the
global stiffness matrix and RHS vector remain strictly untouched.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from wht_iga.solver.steps import resolve_var_name
from wht_iga.solver.trimmed import _d_matrix, _strain_at_points


def _normalize_var_code(name: str) -> str:
    """Normalize variable name to canonical Abaqus short key when known."""
    upper = str(name).strip().upper()
    if upper in ("MISES", "VON_MISES", "VON_MISES_STRESS"):
        return "MISES"
    try:
        return resolve_var_name(upper, target="abaqus")
    except ValueError:
        return upper


def _von_mises_from_voigt(stress_voigt: np.ndarray) -> np.ndarray:
    """Compute scalar von Mises equivalent stress from ``(N, 6)`` Voigt stress."""
    s = np.asarray(stress_voigt, dtype=np.float64)
    sxx, syy, szz = s[:, 0], s[:, 1], s[:, 2]
    sxy, syz, sxz = s[:, 3], s[:, 4], s[:, 5]
    vm2 = 0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (
        sxy ** 2 + syz ** 2 + sxz ** 2
    )
    return np.sqrt(np.maximum(vm2, 0.0))


def _is_box_spline(obj: Any) -> bool:
    return hasattr(obj, "n_cells") and hasattr(obj, "knots") and hasattr(obj, "evaluate")


def _is_thb_hierarchy(obj: Any) -> bool:
    return hasattr(obj, "evaluate_thb_basis_and_grad") and hasattr(obj, "cell_shape")


def _physical_to_param(box: Any, points: np.ndarray, eps: float = 1e-11) -> list[list[float]]:
    """Map physical ``(x, y, z)`` points into ``box`` knot-parameter coordinates with ``w = 0.0``."""
    pts = np.asarray(points, dtype=np.float64)
    out: list[list[float]] = []
    for p in pts:
        uvw: list[float] = []
        for ax in range(3):
            x0 = float(box.lo[ax])
            x1 = float(box.hi[ax])
            k0 = float(box.knots[ax][0])
            k1 = float(box.knots[ax][-1])
            xi = (float(p[ax]) - x0) / max(x1 - x0, 1e-30)
            xi = max(eps, min(1.0 - eps, xi))
            uvw.append(k0 + xi * (k1 - k0))
        uvw.append(0.0)  # Zero quadrature weight: zero contribution to global K and f.
        out.append(uvw)
    return out


def create_zero_weight_probe_element(
    KM: Any,
    model_part: Any,
    volume: Any,
    props: Any,
    box: Any,
    points: np.ndarray,
    *,
    element_id: int = 900001,
    element_name: str = "SmallDisplacementElement3D4N",
) -> Any:
    """Attach a zero-weight quadrature probe element onto ``volume``.

    Every quadrature point receives weight ``0.0`` so ``CalculateLocalSystem``
    adds exact zero stiffness and zero load to the global assembled system.
    """
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"points must have shape (N, 3), got {pts.shape}")
    qp_list = _physical_to_param(box, pts)
    qps = KM.GeometriesVector()
    volume.CreateQuadraturePointGeometries(qps, 2, qp_list)
    if len(qps) != 1:
        raise RuntimeError(f"Expected 1 probe QuadraturePointGeometry, got {len(qps)}")
    el = model_part.CreateNewElement(element_name, int(element_id), qps[0], props)
    el.Initialize(model_part.ProcessInfo)
    return el


def evaluate_zero_weight_kratos_probes(
    KM: Any,
    model_part: Any,
    volume: Any,
    props: Any,
    box: Any,
    points: np.ndarray,
    u_cp: np.ndarray | None = None,
    *,
    element_id: int = 900001,
    verify_zero_contribution: bool = True,
) -> dict[str, np.ndarray]:
    """Evaluate displacement, strain, and Cauchy stress via a Kratos ``weight=0.0`` element."""
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    el = create_zero_weight_probe_element(
        KM, model_part, volume, props, box, pts, element_id=element_id
    )
    pi = model_part.ProcessInfo

    if verify_zero_contribution:
        lhs = KM.Matrix()
        rhs = KM.Vector()
        el.CalculateLocalSystem(lhs, rhs, pi)
        lhs_max = max((abs(lhs[i, j]) for i in range(lhs.Size1()) for j in range(lhs.Size2())), default=0.0)
        rhs_max = max((abs(rhs[i]) for i in range(len(rhs))), default=0.0)
        if lhs_max != 0.0 or rhs_max != 0.0:
            raise AssertionError(
                f"Zero-weight probe element produced non-zero LHS/RHS ({lhs_max}, {rhs_max})"
            )

    if u_cp is None:
        u_list = []
        for nd in sorted(model_part.Nodes, key=lambda n: n.Id):
            disp = nd.GetSolutionStepValue(KM.DISPLACEMENT, 0)
            u_list.append([float(disp[0]), float(disp[1]), float(disp[2])])
        u_cp_arr = np.asarray(u_list, dtype=np.float64)
    else:
        u_cp_arr = np.asarray(u_cp, dtype=np.float64).reshape(-1, 3)

    u_pts = np.asarray(box.evaluate(u_cp_arr, pts), dtype=np.float64)
    s_ip = el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, pi)
    e_ip = el.CalculateOnIntegrationPoints(KM.GREEN_LAGRANGE_STRAIN_VECTOR, pi)

    stress = np.array([[float(v) for v in row] for row in s_ip], dtype=np.float64)
    strain = np.array([[float(v) for v in row] for row in e_ip], dtype=np.float64)
    mises = _von_mises_from_voigt(stress)

    return {
        "U": u_pts,
        "E": strain,
        "S": stress,
        "MISES": mises,
    }


@dataclass
class ProbeSample:
    """Evaluated field values at one or more physical probe coordinates."""

    probe_name: str
    points: np.ndarray
    fields: dict[str, np.ndarray] = field(default_factory=dict)
    s: np.ndarray | None = None
    mode: str = "bspline"

    def get(self, var_name: str) -> np.ndarray:
        """Retrieve field array by either Abaqus short code or Kratos variable name."""
        if var_name in self.fields:
            return self.fields[var_name]
        code = _normalize_var_code(var_name)
        for k, arr in self.fields.items():
            if _normalize_var_code(k) == code:
                return arr
        raise KeyError(f"Variable {var_name!r} not present in probe {self.probe_name!r}")

    def __getitem__(self, var_name: str) -> np.ndarray:
        return self.get(var_name)

    def to_dict(self) -> dict[str, Any]:
        """Convert probe sample into a UTF-8 JSON-serializable dictionary."""
        return {
            "probe_name": self.probe_name,
            "mode": self.mode,
            "n_points": int(self.points.shape[0]),
            "points": self.points.tolist(),
            "s": None if self.s is None else self.s.tolist(),
            "fields": {k: np.asarray(v).tolist() for k, v in self.fields.items()},
        }


def _evaluate_on_source(
    field_source: Any,
    points: np.ndarray,
    variables: Sequence[str],
    *,
    u_cp: np.ndarray | None = None,
    t_cp: np.ndarray | None = None,
    E: float = 1.0,
    nu: float = 0.0,
    kratos_ctx: dict[str, Any] | None = None,
) -> tuple[dict[str, np.ndarray], str]:
    """Evaluate requested field variables at physical points."""
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    n_pts = pts.shape[0]
    norm_vars = {v: _normalize_var_code(v) for v in variables}
    need_codes = set(norm_vars.values())

    computed: dict[str, np.ndarray] = {}
    mode = "bspline"

    if kratos_ctx is not None:
        mode = "zero_weight_kratos"
        k_vals = evaluate_zero_weight_kratos_probes(
            KM=kratos_ctx["KM"],
            model_part=kratos_ctx["model_part"],
            volume=kratos_ctx["volume"],
            props=kratos_ctx["props"],
            box=kratos_ctx.get("box", field_source),
            points=pts,
            u_cp=u_cp,
            element_id=int(kratos_ctx.get("element_id", 900001)),
        )
        computed.update(k_vals)
        if "NT11" in need_codes:
            if t_cp is None:
                raise ValueError("Temperature evaluation ('NT11'/'TEMPERATURE') requires t_cp")
            src = kratos_ctx.get("box", field_source)
            computed["NT11"] = _eval_scalar_field(src, pts, t_cp)

    elif _is_thb_hierarchy(field_source):
        mode = "thb"
        R, dR = field_source.evaluate_thb_basis_and_grad(pts)
        if need_codes & {"U", "E", "S", "MISES"}:
            if u_cp is None:
                raise ValueError("Mechanical field evaluation requires u_cp")
            u_arr = np.asarray(u_cp, dtype=np.float64).reshape(-1, 3)
            if "U" in need_codes:
                computed["U"] = R @ u_arr
            if need_codes & {"E", "S", "MISES"}:
                grad_u = np.einsum("paj,ai->pij", dR, u_arr)
                strain = np.zeros((n_pts, 6), dtype=np.float64)
                strain[:, 0] = grad_u[:, 0, 0]
                strain[:, 1] = grad_u[:, 1, 1]
                strain[:, 2] = grad_u[:, 2, 2]
                strain[:, 3] = grad_u[:, 0, 1] + grad_u[:, 1, 0]
                strain[:, 4] = grad_u[:, 1, 2] + grad_u[:, 2, 1]
                strain[:, 5] = grad_u[:, 0, 2] + grad_u[:, 2, 0]
                computed["E"] = strain
                if need_codes & {"S", "MISES"}:
                    D = _d_matrix(float(E), float(nu))
                    stress = (D @ strain.T).T
                    computed["S"] = stress
                    if "MISES" in need_codes:
                        computed["MISES"] = _von_mises_from_voigt(stress)
        if "NT11" in need_codes:
            if t_cp is None:
                raise ValueError("Temperature evaluation ('NT11'/'TEMPERATURE') requires t_cp")
            t_arr = np.asarray(t_cp, dtype=np.float64).ravel()
            computed["NT11"] = R @ t_arr

    elif _is_box_spline(field_source):
        mode = "bspline"
        if need_codes & {"U", "E", "S", "MISES"}:
            if u_cp is None:
                raise ValueError("Mechanical field evaluation requires u_cp")
            u_arr = np.asarray(u_cp, dtype=np.float64).reshape(-1, 3)
            if "U" in need_codes:
                computed["U"] = np.asarray(field_source.evaluate(u_arr, pts), dtype=np.float64)
            if need_codes & {"E", "S", "MISES"}:
                strain = np.asarray(_strain_at_points(field_source, u_arr.ravel(), pts), dtype=np.float64)
                computed["E"] = strain
                if need_codes & {"S", "MISES"}:
                    D = _d_matrix(float(E), float(nu))
                    stress = (D @ strain.T).T
                    computed["S"] = stress
                    if "MISES" in need_codes:
                        computed["MISES"] = _von_mises_from_voigt(stress)
        if "NT11" in need_codes:
            if t_cp is None:
                raise ValueError("Temperature evaluation ('NT11'/'TEMPERATURE') requires t_cp")
            computed["NT11"] = _eval_scalar_field(field_source, pts, t_cp)

    else:
        raise TypeError(
            f"Unsupported field_source type {type(field_source).__name__}; "
            "expected BoxSpline, THBHierarchy, or kratos_ctx"
        )

    out_fields: dict[str, np.ndarray] = {}
    for orig_name, code in norm_vars.items():
        if code not in computed:
            raise ValueError(f"Variable {orig_name!r} (code={code!r}) could not be evaluated")
        out_fields[orig_name] = computed[code]
    return out_fields, mode


def _eval_scalar_field(field_source: Any, pts: np.ndarray, t_cp: np.ndarray) -> np.ndarray:
    """Evaluate a scalar control-point field such as temperature ``NT11``."""
    t_arr = np.asarray(t_cp, dtype=np.float64).ravel()
    if _is_thb_hierarchy(field_source):
        R = field_source.evaluate_thb_basis(pts)
        return np.asarray(R @ t_arr, dtype=np.float64)
    # BoxSpline.evaluate supports shape (n_ctrl, k)
    vals = np.asarray(field_source.evaluate(t_arr.reshape(-1, 1), pts), dtype=np.float64)
    return vals.ravel()


@dataclass
class PointProbe:
    """Single physical-coordinate field probe."""

    name: str
    point: tuple[float, float, float]
    variables: tuple[str, ...] = ("U", "S", "E")

    def __post_init__(self) -> None:
        if len(self.point) != 3:
            raise ValueError(f"PointProbe expects 3 coordinates, got {self.point!r}")
        self.point = (float(self.point[0]), float(self.point[1]), float(self.point[2]))
        self.variables = tuple(str(v) for v in self.variables)

    def points_array(self) -> np.ndarray:
        """Return ``(1, 3)`` array of physical probe coordinates."""
        return np.array([self.point], dtype=np.float64)

    def evaluate(
        self,
        field_source: Any = None,
        *,
        u_cp: np.ndarray | None = None,
        t_cp: np.ndarray | None = None,
        E: float = 1.0,
        nu: float = 0.0,
        kratos_ctx: dict[str, Any] | None = None,
    ) -> ProbeSample:
        """Evaluate requested variables at this point."""
        pts = self.points_array()
        fields, mode = _evaluate_on_source(
            field_source,
            pts,
            self.variables,
            u_cp=u_cp,
            t_cp=t_cp,
            E=E,
            nu=nu,
            kratos_ctx=kratos_ctx,
        )
        return ProbeSample(probe_name=self.name, points=pts, fields=fields, s=None, mode=mode)


@dataclass
class LineProbe:
    """Straight line segment probe sampling uniformly spaced points."""

    name: str
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    n_points: int = 21
    variables: tuple[str, ...] = ("U", "S", "E")

    def __post_init__(self) -> None:
        if len(self.start) != 3 or len(self.end) != 3:
            raise ValueError("LineProbe start and end must each have 3 coordinates")
        if int(self.n_points) < 2:
            raise ValueError(f"LineProbe n_points must be >= 2, got {self.n_points}")
        self.start = (float(self.start[0]), float(self.start[1]), float(self.start[2]))
        self.end = (float(self.end[0]), float(self.end[1]), float(self.end[2]))
        self.n_points = int(self.n_points)
        self.variables = tuple(str(v) for v in self.variables)

    def points_array(self) -> np.ndarray:
        """Return ``(n_points, 3)`` uniformly spaced coordinates from start to end."""
        t = np.linspace(0.0, 1.0, self.n_points)[:, None]
        p0 = np.asarray(self.start, dtype=np.float64)[None, :]
        p1 = np.asarray(self.end, dtype=np.float64)[None, :]
        return p0 + t * (p1 - p0)

    def arc_lengths(self) -> np.ndarray:
        """Return cumulative arc-length coordinate ``s`` along the line."""
        length = float(np.linalg.norm(np.asarray(self.end) - np.asarray(self.start)))
        return np.linspace(0.0, length, self.n_points, dtype=np.float64)

    def evaluate(
        self,
        field_source: Any = None,
        *,
        u_cp: np.ndarray | None = None,
        t_cp: np.ndarray | None = None,
        E: float = 1.0,
        nu: float = 0.0,
        kratos_ctx: dict[str, Any] | None = None,
    ) -> ProbeSample:
        """Evaluate requested variables along the line segment."""
        pts = self.points_array()
        fields, mode = _evaluate_on_source(
            field_source,
            pts,
            self.variables,
            u_cp=u_cp,
            t_cp=t_cp,
            E=E,
            nu=nu,
            kratos_ctx=kratos_ctx,
        )
        return ProbeSample(
            probe_name=self.name,
            points=pts,
            fields=fields,
            s=self.arc_lengths(),
            mode=mode,
        )


@dataclass
class MultiPointProbe:
    """Arbitrary cloud of physical-coordinate probe points."""

    name: str
    points: np.ndarray | Sequence[Sequence[float]]
    variables: tuple[str, ...] = ("U", "S", "E")

    def __post_init__(self) -> None:
        arr = np.atleast_2d(np.asarray(self.points, dtype=np.float64))
        if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] < 1:
            raise ValueError(f"MultiPointProbe points must have shape (N>=1, 3), got {arr.shape}")
        self.points = arr
        self.variables = tuple(str(v) for v in self.variables)

    def points_array(self) -> np.ndarray:
        """Return ``(N, 3)`` physical coordinate array."""
        return np.asarray(self.points, dtype=np.float64).copy()

    def evaluate(
        self,
        field_source: Any = None,
        *,
        u_cp: np.ndarray | None = None,
        t_cp: np.ndarray | None = None,
        E: float = 1.0,
        nu: float = 0.0,
        kratos_ctx: dict[str, Any] | None = None,
    ) -> ProbeSample:
        """Evaluate requested variables across all probe points."""
        pts = self.points_array()
        fields, mode = _evaluate_on_source(
            field_source,
            pts,
            self.variables,
            u_cp=u_cp,
            t_cp=t_cp,
            E=E,
            nu=nu,
            kratos_ctx=kratos_ctx,
        )
        return ProbeSample(probe_name=self.name, points=pts, fields=fields, s=None, mode=mode)


def export_probes_json(
    path: str | Path,
    samples: Sequence[ProbeSample],
    *,
    metadata: dict[str, Any] | None = None,
) -> Path:
    """Export one or more :class:`ProbeSample` objects to a UTF-8 JSON file."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "schema": "wht_iga.probes.v1",
        "metadata": metadata or {},
        "probes": [sample.to_dict() for sample in samples],
    }
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    return out_path


def export_probes_sqlite(
    db_path: str | Path,
    samples: Sequence[ProbeSample],
    *,
    step_id: int = 1,
    frame_id: int = 1,
) -> Path:
    """Persist evaluated probe samples into a SQLite database (tables ``probe_meta`` & ``probe_values``)."""
    out_path = Path(db_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(out_path))
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS probe_meta (
                probe_name TEXT PRIMARY KEY,
                mode TEXT NOT NULL,
                n_points INTEGER NOT NULL,
                variables_json TEXT NOT NULL
            )
            """
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS probe_values (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                probe_name TEXT NOT NULL,
                step_id INTEGER NOT NULL,
                frame_id INTEGER NOT NULL,
                point_index INTEGER NOT NULL,
                x REAL NOT NULL,
                y REAL NOT NULL,
                z REAL NOT NULL,
                s REAL,
                var_name TEXT NOT NULL,
                values_json TEXT NOT NULL
            )
            """
        )
        for sample in samples:
            cur.execute(
                """
                INSERT OR REPLACE INTO probe_meta (probe_name, mode, n_points, variables_json)
                VALUES (?, ?, ?, ?)
                """,
                (
                    sample.probe_name,
                    sample.mode,
                    int(sample.points.shape[0]),
                    json.dumps(list(sample.fields.keys()), ensure_ascii=False),
                ),
            )
            for idx, pt in enumerate(sample.points):
                s_val = None if sample.s is None else float(sample.s[idx])
                for var_name, arr in sample.fields.items():
                    row_val = arr[idx]
                    val_list = np.atleast_1d(np.asarray(row_val, dtype=np.float64)).tolist()
                    cur.execute(
                        """
                        INSERT INTO probe_values (
                            probe_name, step_id, frame_id, point_index,
                            x, y, z, s, var_name, values_json
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            sample.probe_name,
                            int(step_id),
                            int(frame_id),
                            int(idx),
                            float(pt[0]),
                            float(pt[1]),
                            float(pt[2]),
                            s_val,
                            str(var_name),
                            json.dumps(val_list, ensure_ascii=False),
                        ),
                    )
        conn.commit()
    finally:
        conn.close()
    return out_path
