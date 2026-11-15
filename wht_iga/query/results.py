"""Unified read-only query API over VTKHDF and SQLite results.

Implements :class:`ResultView` with Coreform <-> Kratos variable alias
mapping, KDTree/barycentric probing, hotspot extraction, static equilibrium
checking, plane section force integration, ASME Sec VIII Div 2 stress
linearization, and quadrature error norm evaluation.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.spatial import cKDTree

from wht_iga.post.norms import (
    H1NormResult,
    NormResult,
    integrate_error_against_analytic,
    integrate_error_on_tets,
)
from wht_iga.post.proxy_mesh import (
    VTK_HEXAHEDRON,
    VTK_TETRA,
    ProxyMesh,
    decompose_hex_to_6_tets,
)
from wht_iga.post.vtk_writer import read_vtkhdf

# Single-source bidirectional variable alias table (Coreform <-> Kratos)
FIELD_ALIASES: dict[str, tuple[str, ...]] = {
    "U": ("U", "DISPLACEMENT"),
    "DISPLACEMENT": ("DISPLACEMENT", "U"),
    "S": ("S", "CAUCHY_STRESS_TENSOR"),
    "CAUCHY_STRESS_TENSOR": ("CAUCHY_STRESS_TENSOR", "S"),
    "E": ("E", "STRAIN", "GREEN_LAGRANGE_STRAIN_TENSOR"),
    "STRAIN": ("STRAIN", "E", "GREEN_LAGRANGE_STRAIN_TENSOR"),
    "GREEN_LAGRANGE_STRAIN_TENSOR": ("GREEN_LAGRANGE_STRAIN_TENSOR", "E", "STRAIN"),
    "MISES": ("MISES", "VON_MISES_STRESS"),
    "VON_MISES_STRESS": ("VON_MISES_STRESS", "MISES"),
    "NT11": ("NT11", "TEMPERATURE"),
    "TEMPERATURE": ("TEMPERATURE", "NT11"),
    "RF": ("RF", "REACTION"),
    "REACTION": ("REACTION", "RF"),
}


def resolve_field_name(name: str, available: Iterable[str]) -> str:
    """Resolve ``name`` against ``available`` keys using :data:`FIELD_ALIASES`."""
    avail_list = list(available)
    if name in avail_list:
        return name

    upper = name.upper()
    candidates = FIELD_ALIASES.get(upper, (name,))
    avail_upper = {k.upper(): k for k in avail_list}
    for cand in candidates:
        if cand in avail_list:
            return cand
        if cand.upper() in avail_upper:
            return avail_upper[cand.upper()]

    raise KeyError(f"Field '{name}' not found; available fields: {avail_list}")


def stress_to_matrix(s: np.ndarray) -> np.ndarray:
    """Convert Voigt stress ``(..., 6)`` or ``(..., 3, 3)`` to ``(..., 3, 3)``."""
    arr = np.asarray(s, dtype=np.float64)
    if arr.ndim >= 2 and arr.shape[-2:] == (3, 3):
        return arr
    if arr.shape[-1] == 6:
        out = np.zeros(arr.shape[:-1] + (3, 3), dtype=np.float64)
        out[..., 0, 0] = arr[..., 0]
        out[..., 1, 1] = arr[..., 1]
        out[..., 2, 2] = arr[..., 2]
        out[..., 0, 1] = out[..., 1, 0] = arr[..., 3]
        out[..., 1, 2] = out[..., 2, 1] = arr[..., 4]
        out[..., 0, 2] = out[..., 2, 0] = arr[..., 5]
        return out
    raise ValueError(f"Unsupported stress shape {arr.shape}; expected (..., 6) or (..., 3, 3)")


def _stress_to_voigt6(s: np.ndarray) -> np.ndarray:
    arr = np.asarray(s, dtype=np.float64)
    if arr.shape[-1] == 6:
        return arr
    if arr.ndim >= 2 and arr.shape[-2:] == (3, 3):
        out = np.empty(arr.shape[:-2] + (6,), dtype=np.float64)
        out[..., 0] = arr[..., 0, 0]
        out[..., 1] = arr[..., 1, 1]
        out[..., 2] = arr[..., 2, 2]
        out[..., 3] = arr[..., 0, 1]
        out[..., 4] = arr[..., 1, 2]
        out[..., 5] = arr[..., 0, 2]
        return out
    raise ValueError(f"Unsupported stress shape {arr.shape}; expected (..., 6) or (..., 3, 3)")


def von_mises_stress(s: np.ndarray) -> np.ndarray:
    """Compute von Mises equivalent stress from Voigt 6 or 3x3 stress tensors."""
    v = _stress_to_voigt6(s)
    sxx, syy, szz = v[..., 0], v[..., 1], v[..., 2]
    sxy, syz, sxz = v[..., 3], v[..., 4], v[..., 5]
    vm_sq = 0.5 * (
        (sxx - syy) ** 2
        + (syy - szz) ** 2
        + (szz - sxx) ** 2
        + 6.0 * (sxy**2 + syz**2 + sxz**2)
    )
    return np.sqrt(np.maximum(vm_sq, 0.0))


def _field_magnitude(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim == 1:
        return np.abs(a)
    if a.ndim >= 2 and (a.shape[-1] == 6 or a.shape[-2:] == (3, 3)):
        return von_mises_stress(a)
    return np.linalg.norm(a.reshape((a.shape[0], -1)), axis=1)


def _simpson_weights(n: int, length: float) -> np.ndarray:
    """Composite Simpson weights on ``n`` uniform points over ``[0, length]``."""
    if n < 2:
        raise ValueError("n_samples must be >= 2")
    if n == 2:
        return np.array([0.5 * length, 0.5 * length], dtype=np.float64)
    h = length / float(n - 1)
    w = np.zeros(n, dtype=np.float64)
    if (n - 1) % 2 == 0:
        w[0] = 1.0
        w[-1] = 1.0
        w[1:-1:2] = 4.0
        w[2:-2:2] = 2.0
        return w * (h / 3.0)
    # Odd number of intervals: Simpson 1/3 up to n-4, Simpson 3/8 on last 3 intervals
    if n >= 5:
        m = n - 3
        w[0] = 1.0
        w[m - 1] += 1.0
        w[1 : m - 1 : 2] = 4.0
        w[2 : m - 2 : 2] = 2.0
        w[:m] *= h / 3.0
        w[m - 1 :] += np.array([1.0, 3.0, 3.0, 1.0]) * (3.0 * h / 8.0)
        return w
    # n == 4: pure Simpson 3/8
    return np.array([1.0, 3.0, 3.0, 1.0], dtype=np.float64) * (3.0 * h / 8.0)


@dataclass(frozen=True)
class EquilibriumReport:
    """Static force equilibrium check summary."""

    applied: np.ndarray
    reaction: np.ndarray
    residual: np.ndarray
    residual_norm: float
    rel_error: float
    balanced: bool


@dataclass(frozen=True)
class StressLinearizationResult:
    """ASME Sec VIII Div 2 stress linearization along a classification line."""

    p1: np.ndarray
    p2: np.ndarray
    length: float
    s_coord: np.ndarray
    t_dir: np.ndarray
    sampled_stress: np.ndarray
    membrane: np.ndarray
    bending_surface1: np.ndarray
    bending_surface2: np.ndarray
    peak_surface1: np.ndarray
    peak_surface2: np.ndarray
    linearized_stress: np.ndarray
    von_mises_membrane: float
    von_mises_membrane_plus_bending_1: float
    von_mises_membrane_plus_bending_2: float


class ResultView:
    """Unified read-only query interface over ``VTKHDF`` and ``SQLite`` outputs.

    Parameters
    ----------
    vtkhdf_path:
        Path to ``<job>.vtkhdf`` produced by :func:`wht_iga.post.write_vtkhdf`.
    sqlite_path:
        Optional path to ``<job>_iga_results.sqlite`` matching
        ``wht_iga.solver.db`` schema.
    """

    def __init__(
        self,
        vtkhdf_path: str | Path,
        sqlite_path: str | Path | None = None,
    ) -> None:
        self.vtkhdf_path = Path(vtkhdf_path)
        self.sqlite_path = Path(sqlite_path) if sqlite_path is not None else None
        self._hdf = read_vtkhdf(self.vtkhdf_path)
        self._sqlite_cache: dict | None = None

    def _load_sqlite(self) -> dict:
        if self._sqlite_cache is not None:
            return self._sqlite_cache
        if self.sqlite_path is None:
            self._sqlite_cache = {}
            return self._sqlite_cache

        con = sqlite3.connect(str(self.sqlite_path))
        try:
            con.row_factory = sqlite3.Row
            tables = {
                r[0]
                for r in con.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                ).fetchall()
            }
            meta = (
                dict(con.execute("SELECT key, value FROM meta").fetchall())
                if "meta" in tables
                else {}
            )
            steps = (
                [dict(r) for r in con.execute("SELECT * FROM step ORDER BY step_id")]
                if "step" in tables
                else []
            )
            frames = (
                [dict(r) for r in con.execute("SELECT * FROM frame ORDER BY frame_id")]
                if "frame" in tables
                else []
            )
            cp_fields: dict[tuple[int, str], np.ndarray] = {}
            if "cp_field" in tables:
                for r in con.execute("SELECT * FROM cp_field"):
                    shape = tuple(json.loads(r["shape"]))
                    cp_fields[(int(r["frame_id"]), r["var"])] = (
                        np.frombuffer(r["data"], dtype=np.float64).reshape(shape).copy()
                    )
            qp_fields: dict[tuple[int, str], np.ndarray] = {}
            if "qp_field" in tables:
                for r in con.execute("SELECT * FROM qp_field"):
                    shape = tuple(json.loads(r["shape"]))
                    qp_fields[(int(r["frame_id"]), r["var"])] = (
                        np.frombuffer(r["data"], dtype=np.float64).reshape(shape).copy()
                    )
            probes = (
                [dict(r) for r in con.execute("SELECT * FROM probe ORDER BY probe_id")]
                if "probe" in tables
                else []
            )
            self._sqlite_cache = {
                "meta": meta,
                "steps": steps,
                "frames": frames,
                "cp_fields": cp_fields,
                "qp_fields": qp_fields,
                "probes": probes,
            }
        finally:
            con.close()
        return self._sqlite_cache

    def meta(self) -> dict[str, str]:
        """Return SQLite ``meta`` dictionary (empty if no SQLite DB attached)."""
        return dict(self._load_sqlite().get("meta", {}))

    def frames(self) -> list[dict]:
        """Return per-frame metadata combining SQLite and VTKHDF time steps."""
        sq_frames = self._load_sqlite().get("frames", [])
        if sq_frames:
            return [dict(f) for f in sq_frames]
        times = self._hdf["times"]
        return [
            {"frame_id": int(i), "step_id": 0, "time": float(t)}
            for i, t in enumerate(times)
        ]

    def _get_mesh(self, frame: int) -> ProxyMesh:
        frames = self._hdf["frames"]
        if not 0 <= frame < len(frames):
            raise IndexError(f"Frame {frame} out of range [0, {len(frames)})")
        return frames[frame]

    def field(
        self,
        name: str,
        frame: int = 0,
        block: int | str | None = None,
        location: str = "auto",
    ) -> np.ndarray:
        """Retrieve a field array by name or Coreform/Kratos alias.

        Parameters
        ----------
        name:
            Field name or alias (e.g. ``"U"`` / ``"DISPLACEMENT"``).
        frame:
            Frame index (0-based).
        block:
            Optional block index or block name to filter cell/point data.
        location:
            ``"auto"``, ``"point"``, ``"cell"``, or ``"cp"``.
        """
        mesh = self._get_mesh(frame)
        if location in ("auto", "point") and mesh.point_data:
            try:
                key = resolve_field_name(name, mesh.point_data.keys())
                arr = mesh.point_data[key]
                if block is not None and "BlockId" in mesh.cell_data:
                    bid = self._resolve_block_id(block)
                    cell_mask = mesh.cell_data["BlockId"] == bid
                    pt_ids = np.unique(
                        np.concatenate(
                            [mesh.cell_nodes(i) for i in np.flatnonzero(cell_mask)]
                        )
                    )
                    return arr[pt_ids]
                return arr
            except KeyError:
                if location == "point":
                    raise

        if location in ("auto", "cell") and mesh.cell_data:
            try:
                key = resolve_field_name(name, mesh.cell_data.keys())
                arr = mesh.cell_data[key]
                if block is not None and "BlockId" in mesh.cell_data:
                    bid = self._resolve_block_id(block)
                    return arr[mesh.cell_data["BlockId"] == bid]
                return arr
            except KeyError:
                if location == "cell":
                    raise

        if location in ("auto", "cp"):
            return self.control_point_field(name, frame=frame)

        raise KeyError(f"Field '{name}' not found at frame {frame} (location={location})")

    def _resolve_block_id(self, block: int | str) -> float:
        if isinstance(block, (int, np.integer)):
            return float(block)
        names = self._hdf.get("block_names", [])
        if block in names:
            return float(names.index(block))
        raise KeyError(f"Block '{block}' not found in {names}")

    def control_point_field(self, name: str, frame: int = 0) -> np.ndarray:
        """Read a control-point array from ``/WHT_IGA/control_points`` or SQLite."""
        cp_all = self._hdf.get("control_points", {}).get(frame, {})
        if cp_all:
            try:
                key = resolve_field_name(name, cp_all.keys())
                return cp_all[key]
            except KeyError:
                pass

        sq_cp = self._load_sqlite().get("cp_fields", {})
        frame_vars = {var: arr for (fid, var), arr in sq_cp.items() if fid == frame}
        if frame_vars:
            key = resolve_field_name(name, frame_vars.keys())
            return frame_vars[key]

        raise KeyError(f"Control-point field '{name}' not found for frame {frame}")

    def _all_tets(self, mesh: ProxyMesh) -> np.ndarray:
        """Convert all hex and tet cells in ``mesh`` into a ``(N_tet, 4)`` array."""
        tets: list[np.ndarray] = []
        hex_list: list[np.ndarray] = []
        tet_list: list[np.ndarray] = []
        for idx, ctype in enumerate(mesh.cell_types):
            nodes = mesh.cell_nodes(idx)
            if int(ctype) == VTK_TETRA:
                tet_list.append(nodes)
            elif int(ctype) == VTK_HEXAHEDRON:
                hex_list.append(nodes)
        if tet_list:
            tets.append(np.asarray(tet_list, dtype=np.int64))
        if hex_list:
            tets.append(decompose_hex_to_6_tets(np.asarray(hex_list, dtype=np.int64)))
        if not tets:
            return np.zeros((0, 4), dtype=np.int64)
        return np.vstack(tets)

    def probe(
        self,
        points: np.ndarray,
        field_name: str = "U",
        frame: int = 0,
        tol: float = 1e-8,
    ) -> np.ndarray:
        """Evaluate ``field_name`` at arbitrary 3D ``points`` via barycentric/KDTree lookup."""
        pts = np.asarray(points, dtype=np.float64)
        single = pts.ndim == 1
        if single:
            pts = pts[None, :]

        mesh = self._get_mesh(frame)
        values = self.field(field_name, frame=frame, location="point")
        out_shape = (pts.shape[0],) + values.shape[1:]
        out = np.empty(out_shape, dtype=np.float64)

        tets = self._all_tets(mesh)
        pt_tree = cKDTree(mesh.points)
        tet_tree = None
        tet_centroids = None
        if len(tets) > 0:
            tet_centroids = mesh.points[tets].mean(axis=1)
            tet_tree = cKDTree(tet_centroids)

        for i, p in enumerate(pts):
            found = False
            if tet_tree is not None:
                k = min(24, len(tets))
                _, cand_ids = tet_tree.query(p, k=k)
                cand_arr = np.atleast_1d(cand_ids)
                for tid in cand_arr:
                    nodes = tets[int(tid)]
                    v = mesh.points[nodes]
                    mat = np.column_stack([v[1] - v[0], v[2] - v[0], v[3] - v[0]])
                    det = float(np.linalg.det(mat))
                    if abs(det) < 1e-20:
                        continue
                    lam123 = np.linalg.solve(mat, p - v[0])
                    lam0 = 1.0 - float(np.sum(lam123))
                    bary = np.array([lam0, lam123[0], lam123[1], lam123[2]], dtype=np.float64)
                    if np.all(bary >= -tol):
                        out[i] = np.tensordot(bary, values[nodes], axes=(0, 0))
                        found = True
                        break

            if not found:
                k_nn = min(4, mesh.n_points)
                dists, nn_ids = pt_tree.query(p, k=k_nn)
                dists_arr = np.atleast_1d(dists)
                nn_arr = np.atleast_1d(nn_ids)
                if dists_arr[0] < 1e-12:
                    out[i] = values[int(nn_arr[0])]
                else:
                    w = 1.0 / np.maximum(dists_arr, 1e-15)
                    w /= np.sum(w)
                    out[i] = np.tensordot(w, values[nn_arr], axes=(0, 0))

        return out[0] if single else out

    def hotspots(
        self,
        field_name: str = "S",
        top_n: int = 10,
        frame: int = 0,
    ) -> list[dict]:
        """Return top-N locations sorted by scalar magnitude / von Mises stress."""
        mesh = self._get_mesh(frame)
        arr = self.field(field_name, frame=frame, location="point")
        mag = _field_magnitude(arr)
        order = np.argsort(mag)[::-1][: int(top_n)]
        return [
            {
                "index": int(idx),
                "coords": mesh.points[int(idx)].copy(),
                "value": float(mag[int(idx)]),
                "raw": arr[int(idx)].copy(),
            }
            for idx in order
        ]

    def reaction_force(self, frame: int = 0, field_name: str = "RF") -> np.ndarray:
        """Compute total reaction force vector by summing nodal/CP reactions."""
        arr = self.field(field_name, frame=frame)
        return np.sum(arr, axis=0)

    def check_equilibrium(
        self,
        applied_force: np.ndarray,
        reaction_force: np.ndarray | None = None,
        frame: int = 0,
        tol: float = 1e-10,
    ) -> EquilibriumReport:
        """Check global static force balance ``applied + reaction == 0``."""
        f_app = np.asarray(applied_force, dtype=np.float64).ravel()
        f_reac = (
            self.reaction_force(frame=frame)
            if reaction_force is None
            else np.asarray(reaction_force, dtype=np.float64).ravel()
        )
        res = f_app + f_reac
        res_norm = float(np.linalg.norm(res))
        denom = max(float(np.linalg.norm(f_app)), float(np.linalg.norm(f_reac)), 1e-30)
        rel_err = res_norm / denom
        return EquilibriumReport(
            applied=f_app,
            reaction=f_reac,
            residual=res,
            residual_norm=res_norm,
            rel_error=rel_err,
            balanced=rel_err <= tol,
        )

    def section_force(
        self,
        plane_point: np.ndarray,
        plane_normal: np.ndarray,
        stress_field: str = "S",
        frame: int = 0,
    ) -> np.ndarray:
        """Integrate traction ``sigma @ n`` over the intersection of a plane with the mesh."""
        p0 = np.asarray(plane_point, dtype=np.float64).ravel()
        n = np.asarray(plane_normal, dtype=np.float64).ravel()
        n_norm = float(np.linalg.norm(n))
        if n_norm == 0.0:
            raise ValueError("plane_normal must be non-zero")
        n = n / n_norm

        mesh = self._get_mesh(frame)
        stress_pt = stress_to_matrix(self.field(stress_field, frame=frame, location="point"))
        tets = self._all_tets(mesh)
        if len(tets) == 0:
            return np.zeros(3, dtype=np.float64)

        signed_d = (mesh.points - p0) @ n
        total_force = np.zeros(3, dtype=np.float64)

        for nodes in tets:
            d = signed_d[nodes]
            # Strictly on one side -> no intersection
            if np.all(d > 1e-12) or np.all(d < -1e-12):
                continue
            # Avoid double-counting faces lying exactly on the plane unless straddling
            pos = [k for k in range(4) if d[k] >= 0.0]
            neg = [k for k in range(4) if d[k] < 0.0]
            if not pos or not neg:
                continue

            inter_pts: list[np.ndarray] = []
            inter_sig: list[np.ndarray] = []
            for ia in pos:
                for ib in neg:
                    da, db = float(d[ia]), float(d[ib])
                    t = da / (da - db)
                    pa, pb = mesh.points[nodes[ia]], mesh.points[nodes[ib]]
                    sa, sb = stress_pt[nodes[ia]], stress_pt[nodes[ib]]
                    inter_pts.append(pa + t * (pb - pa))
                    inter_sig.append(sa + t * (sb - sa))

            if len(inter_pts) < 3:
                continue

            pts_poly = np.asarray(inter_pts, dtype=np.float64)
            sig_poly = np.asarray(inter_sig, dtype=np.float64)

            if len(pts_poly) == 4:
                # Order quad vertices cyclically around centroid in plane
                c = pts_poly.mean(axis=0)
                e1 = pts_poly[0] - c
                e1_norm = np.linalg.norm(e1)
                if e1_norm < 1e-15:
                    continue
                e1 /= e1_norm
                e2 = np.cross(n, e1)
                angles = np.arctan2((pts_poly - c) @ e2, (pts_poly - c) @ e1)
                order = np.argsort(angles)
                pts_poly = pts_poly[order]
                sig_poly = sig_poly[order]
                tris = [(0, 1, 2), (0, 2, 3)]
            else:
                tris = [(0, 1, 2)]

            for i0, i1, i2 in tris:
                area = 0.5 * float(
                    np.linalg.norm(
                        np.cross(pts_poly[i1] - pts_poly[i0], pts_poly[i2] - pts_poly[i0])
                    )
                )
                sig_c = (sig_poly[i0] + sig_poly[i1] + sig_poly[i2]) / 3.0
                total_force += area * (sig_c @ n)

        return total_force

    def stress_linearize(
        self,
        p1: np.ndarray,
        p2: np.ndarray,
        n_samples: int = 41,
        field_name: str = "S",
        frame: int = 0,
    ) -> StressLinearizationResult:
        """Perform ASME Sec VIII Div 2 stress linearization along ``p1 -> p2``."""
        a = np.asarray(p1, dtype=np.float64).ravel()
        b = np.asarray(p2, dtype=np.float64).ravel()
        vec = b - a
        length = float(np.linalg.norm(vec))
        if length <= 0.0:
            raise ValueError("p1 and p2 must be distinct points")
        t_dir = vec / length

        s_coord = np.linspace(0.0, length, int(n_samples))
        line_pts = a[None, :] + (s_coord[:, None] / length) * vec[None, :]
        sampled_raw = self.probe(line_pts, field_name=field_name, frame=frame)
        sampled = _stress_to_voigt6(sampled_raw)

        w = _simpson_weights(len(s_coord), length)
        membrane = np.sum(w[:, None] * sampled, axis=0) / length
        s_rel = s_coord - 0.5 * length
        bending_2 = (6.0 / (length**2)) * np.sum((w * s_rel)[:, None] * sampled, axis=0)
        bending_1 = -bending_2

        peak_1 = sampled[0] - membrane - bending_1
        peak_2 = sampled[-1] - membrane - bending_2

        lin_factor = (2.0 * s_coord / length - 1.0)[:, None]
        linearized = membrane[None, :] + lin_factor * bending_2[None, :]

        return StressLinearizationResult(
            p1=a,
            p2=b,
            length=length,
            s_coord=s_coord,
            t_dir=t_dir,
            sampled_stress=sampled,
            membrane=membrane,
            bending_surface1=bending_1,
            bending_surface2=bending_2,
            peak_surface1=peak_1,
            peak_surface2=peak_2,
            linearized_stress=linearized,
            von_mises_membrane=float(von_mises_stress(membrane)),
            von_mises_membrane_plus_bending_1=float(von_mises_stress(membrane + bending_1)),
            von_mises_membrane_plus_bending_2=float(von_mises_stress(membrane + bending_2)),
        )

    def error_norm(
        self,
        analytic_fn: Callable[[np.ndarray], np.ndarray],
        field_name: str = "U",
        frame: int = 0,
        cell_bounds: np.ndarray | None = None,
        n_gp: int = 3,
        eval_num_fn: Callable[[np.ndarray], np.ndarray] | None = None,
        eval_grad_fn: Callable[[np.ndarray], np.ndarray] | None = None,
        analytic_grad_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> NormResult | H1NormResult:
        """Compute ``L2`` or ``H1`` error norm against ``analytic_fn``.

        Uses ``eval_num_fn`` if provided, otherwise probes the proxy mesh via
        :meth:`probe`. Integrates over ``cell_bounds`` (brick Gauss-Legendre)
        when provided, or over the tetrahedralized proxy mesh otherwise.
        """
        num_fn = eval_num_fn or (lambda pts: self.probe(pts, field_name=field_name, frame=frame))
        if cell_bounds is not None:
            return integrate_error_against_analytic(
                cell_bounds=cell_bounds,
                eval_num_fn=num_fn,
                analytic_fn=analytic_fn,
                n_gp=n_gp,
                eval_grad_fn=eval_grad_fn,
                analytic_grad_fn=analytic_grad_fn,
            )
        mesh = self._get_mesh(frame)
        tets = self._all_tets(mesh)
        return integrate_error_on_tets(mesh.points, tets, num_fn, analytic_fn)
