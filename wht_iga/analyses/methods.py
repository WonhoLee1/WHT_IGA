"""Selectable multi-method discretization options and unified solver facade.

Supports five first-class structural discretization strategies over a unified
specification and solver interface, without importing Kratos at module load time:

* ``"full_iga_body_fit"`` — Untrimmed body-fit B-spline IGA via ``wht_iga.solver.pipeline``.
* ``"full_iga_immersed"`` — Immersed/trimmed B-spline IGA (FCM / cut-cell quadrature).
* ``"thb_adaptive_iga"`` — Locally refined Truncated Hierarchical B-splines (THB).
* ``"iga_fe_hybrid"`` — Coupled Native FEM + IGA domains across a split interface.
* ``"native_fe"`` — Standard Kratos solid finite elements (hex8/hex20/tet4/tet10/wedge6).
"""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Literal

import numpy as np

from wht_iga.analyses.methods_options import (
    Degree,
    DiscretizationOption,
    FeElementType,
    FullIgaOption,
    HybridCouplingKind,
    IgaFeHybridOption,
    LockingMitigation,
    MethodDescriptor,
    MethodId,
    NativeFeOption,
    ThbAdaptiveIgaOption,
    list_discretization_methods,
    parse_discretization_option,
)
from wht_iga.analyses.methods_recommend import (
    DiscretizationRecommendation,
    PhysicsKind,
    recommend_discretization,
)
from wht_iga.analyses.methods_solve import (
    BoxMethodProblem,
    MethodSolution,
    compare_method_solutions,
    solve_fbar_box,
    solve_with_discretization_method,
)

Vec3 = tuple[float, float, float]


@dataclass(frozen=True)
class ObbFrame:
    """Right-handed orthonormal Oriented Bounding Box (OBB) frame.

    Maps between global coordinates and local OBB coordinates via
    ``x_local = (x_global - origin) @ R`` and ``x_global = x_local @ R.T + origin``,
    where ``R = np.column_stack(axes)`` has orthonormal columns (``det(R) = +1``).

    Attributes
    ----------
    origin:
        Frame origin in global coordinates ``(x0, y0, z0)``.
    axes:
        Three right-handed orthonormal unit vectors ``(e0, e1, e2)`` in global
        coordinates forming the columns of the local-to-global rotation ``R``.
    """

    origin: Vec3
    axes: tuple[Vec3, Vec3, Vec3]

    def rotation_matrix(self) -> np.ndarray:
        """Return ``3x3`` orthonormal matrix ``R`` whose columns are ``axes``."""
        return np.column_stack([np.asarray(a, dtype=np.float64) for a in self.axes])

    def to_local_points(self, points: np.ndarray) -> np.ndarray:
        """Transform global points ``(..., 3)`` into local OBB coordinates."""
        pts = np.asarray(points, dtype=np.float64)
        orig = np.asarray(self.origin, dtype=np.float64)
        return (pts - orig) @ self.rotation_matrix()

    def to_global_points(self, points: np.ndarray) -> np.ndarray:
        """Transform local OBB points ``(..., 3)`` back to global coordinates."""
        pts = np.asarray(points, dtype=np.float64)
        orig = np.asarray(self.origin, dtype=np.float64)
        return pts @ self.rotation_matrix().T + orig

    def to_local_vectors(self, vecs: np.ndarray) -> np.ndarray:
        """Rotate global vectors ``(..., 3)`` into local OBB coordinates."""
        v = np.asarray(vecs, dtype=np.float64)
        return v @ self.rotation_matrix()

    def to_global_vectors(self, vecs: np.ndarray) -> np.ndarray:
        """Rotate local OBB vectors ``(..., 3)`` back to global coordinates."""
        v = np.asarray(vecs, dtype=np.float64)
        return v @ self.rotation_matrix().T

    def to_global_tensors(self, tensors: np.ndarray) -> np.ndarray:
        """Rotate local 3x3 tensors ``(..., 3, 3)`` to global ``R @ S @ R.T``."""
        rot = self.rotation_matrix()
        s = np.asarray(tensors, dtype=np.float64)
        if s.ndim == 2:
            return rot @ s @ rot.T
        return np.einsum("ij,...jk,lk->...il", rot, s, rot)


def compute_pca_obb(
    points: np.ndarray,
    padding_frac: float = 0.02,
) -> tuple[ObbFrame, Vec3, Vec3]:
    """Compute a canonical right-handed PCA OBB frame and local bounds.

    Parameters
    ----------
    points:
        Array of shape ``(N, 3)`` with ``N >= 4`` finite 3D points.
    padding_frac:
        Non-negative fractional padding applied to each local span.

    Returns
    -------
    tuple[ObbFrame, Vec3, Vec3]
        ``(frame, lo_local, hi_local)`` with canonical eigenvector signs and
        ``det(R) = +1``.
    """
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 4:
        raise ValueError(f"points must have shape (N>=4, 3), got {pts.shape}")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must contain finite values")
    pad_frac = float(padding_frac)
    if not np.isfinite(pad_frac) or pad_frac < 0.0:
        raise ValueError(f"padding_frac must be finite and >= 0, got {padding_frac}")

    origin = pts.mean(axis=0)
    centered = pts - origin
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    if float(singular[2]) <= 1e-12 * max(float(singular[0]), 1e-30):
        raise ValueError("degenerate point cloud: rank < 3, OBB frame undefined")

    # Canonicalize eigenvector signs by largest-magnitude component, then enforce right-handedness.
    v0 = vt[0].copy()
    v1 = vt[1].copy()
    for vec in (v0, v1):
        idx = int(np.argmax(np.abs(vec)))
        if vec[idx] < 0.0:
            vec *= -1.0
    v2 = np.cross(v0, v1)
    v2 /= float(np.linalg.norm(v2))

    rot = np.column_stack([v0, v1, v2])
    local_pts = centered @ rot
    q_min = local_pts.min(axis=0)
    q_max = local_pts.max(axis=0)
    span = np.maximum(q_max - q_min, 1e-12)
    pad = pad_frac * span

    frame = ObbFrame(
        origin=(float(origin[0]), float(origin[1]), float(origin[2])),
        axes=(
            (float(v0[0]), float(v0[1]), float(v0[2])),
            (float(v1[0]), float(v1[1]), float(v1[2])),
            (float(v2[0]), float(v2[1]), float(v2[2])),
        ),
    )
    lo_local: Vec3 = (
        float(q_min[0] - pad[0]),
        float(q_min[1] - pad[1]),
        float(q_min[2] - pad[2]),
    )
    hi_local: Vec3 = (
        float(q_max[0] + pad[0]),
        float(q_max[1] + pad[1]),
        float(q_max[2] + pad[2]),
    )
    return frame, lo_local, hi_local


def solve_immersed_obb(
    points_soup: np.ndarray,
    tris_soup: np.ndarray,
    *,
    E: float,
    nu: float,
    n_cells: tuple[int, int, int] = (4, 2, 2),
    degree: int = 2,
    continuity: int | None = None,
    fixed_face_local: tuple[int, int] = (0, 0),
    load_face_local: tuple[int, int] = (0, 1),
    traction_global: Vec3 = (1.0, 0.0, 0.0),
    obb: ObbFrame | None = None,
    padding_frac: float = 0.0,
    quadrature_degree: int | None = None,
    out_mode: Literal["drop", "ghost"] = "drop",
    ghost_alpha: float = 1e-8,
    bc_penalty: float | None = None,
    use_pardiso: bool = True,
) -> MethodSolution:
    """Solve an immersed triangle-soup linear elasticity problem in its local OBB frame.

    Transforms the surface soup and global traction into the OBB local frame,
    builds a minimal local ``BoxSpline`` + trimmed quadrature table, solves via
    ``wht_iga.solver.pipeline.run_linear``, and maps evaluated displacements back
    to global coordinates (``u_global = R @ u_local``).
    """
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.trimmed import build_trimmed_table
    from wht_iga.solver.validate import DirichletBC, JobSpec, TractionLoad

    t_start = time.perf_counter()
    pts_global = np.asarray(points_soup, dtype=np.float64)
    tris = np.asarray(tris_soup, dtype=np.int64)

    if obb is None:
        frame, lo_local, hi_local = compute_pca_obb(pts_global, padding_frac=padding_frac)
    else:
        frame = obb
        local_q = frame.to_local_points(pts_global)
        q_min = local_q.min(axis=0)
        q_max = local_q.max(axis=0)
        pad = float(padding_frac) * np.maximum(q_max - q_min, 1e-12)
        lo_local = (
            float(q_min[0] - pad[0]),
            float(q_min[1] - pad[1]),
            float(q_min[2] - pad[2]),
        )
        hi_local = (
            float(q_max[0] + pad[0]),
            float(q_max[1] + pad[1]),
            float(q_max[2] + pad[2]),
        )

    pts_local = frame.to_local_points(pts_global)
    cont = (degree - 1) if continuity is None else int(continuity)
    q_deg = min(2 * int(degree), 4) if quadrature_degree is None else int(quadrature_degree)

    box = BoxSpline(
        lo=lo_local,
        hi=hi_local,
        n_cells=n_cells,
        degree=int(degree),
        continuity=cont,
    )
    qp_table, out_cells, geo_issues, cell_records, _ = build_trimmed_table(
        pts_local,
        tris,
        box,
        q=q_deg,
    )

    t_loc_arr = frame.to_local_vectors(np.asarray(traction_global, dtype=np.float64))
    traction_local: Vec3 = (float(t_loc_arr[0]), float(t_loc_arr[1]), float(t_loc_arr[2]))

    job_spec = JobSpec(
        name="full_iga_immersed_obb",
        lo=lo_local,
        hi=hi_local,
        n_cells=n_cells,
        degree=int(degree),
        continuity=cont,
        E=float(E),
        nu=float(nu),
        beta=bc_penalty,
        bcs=[DirichletBC(face=fixed_face_local, components=(True, True, True), value=(0.0, 0.0, 0.0))],
        tractions=[TractionLoad(face=load_face_local, vector=traction_local)],
        out_cells=out_cells,
        ghost_alpha=float(ghost_alpha),
        out_mode=out_mode,
        linear_solver="pardiso_ldlt" if use_pardiso else "auto",
        enable_solver_log=False,
    )
    res = run_linear(
        job_spec,
        geo_issues=geo_issues,
        qp_table=qp_table,
        cell_records=cell_records,
    )
    coeffs = res.U.reshape(-1, 3)

    def _eval_global(pts: np.ndarray) -> np.ndarray:
        pts_in = np.asarray(pts, dtype=np.float64)
        pts_loc = frame.to_local_points(pts_in)
        u_loc = box.evaluate(coeffs, pts_loc)
        return frame.to_global_vectors(u_loc)

    tip_loc = 0.5 * (np.asarray(lo_local, dtype=np.float64) + np.asarray(hi_local, dtype=np.float64))
    ax_load, side_load = load_face_local
    tip_loc[ax_load] = hi_local[ax_load] if side_load == 1 else lo_local[ax_load]
    tip_global = frame.to_global_points(tip_loc.reshape(1, 3))[0]
    tip_disp_global = _eval_global(tip_global.reshape(1, 3))[0]

    total_wall = time.perf_counter() - t_start
    timings = dict(res.timings)
    timings["total_s"] = total_wall

    return MethodSolution(
        method="full_iga_immersed",
        total_dofs=int(res.total_dofs or box.n_dofs),
        active_dofs=int(res.active_dofs or box.n_dofs),
        tip_displacement=tip_disp_global,
        max_displacement=float(np.max(np.abs(res.U))) if res.U.size else 0.0,
        converged=bool(res.converged),
        timings=timings,
        native={
            "linear_result": res,
            "box": box,
            "obb": frame,
            "lo_local": lo_local,
            "hi_local": hi_local,
            "out_cells": out_cells,
            "n_out_cells": len(out_cells),
        },
        _eval_fn=_eval_global,
    )


__all__ = [
    "BoxMethodProblem",
    "Degree",
    "DiscretizationOption",
    "DiscretizationRecommendation",
    "FeElementType",
    "FullIgaOption",
    "HybridCouplingKind",
    "IgaFeHybridOption",
    "LockingMitigation",
    "MethodDescriptor",
    "MethodId",
    "MethodSolution",
    "NativeFeOption",
    "ObbFrame",
    "PhysicsKind",
    "ThbAdaptiveIgaOption",
    "compare_method_solutions",
    "compute_pca_obb",
    "list_discretization_methods",
    "parse_discretization_option",
    "recommend_discretization",
    "solve_fbar_box",
    "solve_immersed_obb",
    "solve_with_discretization_method",
]
