# -*- coding: utf-8 -*-
"""Oriented bounding box (OBB) background-grid frame adapter (Track H).

For parts placed obliquely to the global axes, an axis-aligned background grid
wastes OUT cells and misaligns the spline. The OBB adapter works purely in the
local principal frame and reuses the existing axis-aligned pipeline unchanged:

1. ``compute_obb_frame`` finds the rigid frame ``x' = R @ (x - x0)`` (PCA/SVD).
2. Geometry points, Dirichlet/Neumann vectors are rotated into the local frame.
3. The existing ``BoxSpline`` / ``classify`` / ``quadrature`` / Kratos pipeline
   solves in the local frame (isotropic elasticity is rotation-invariant).
4. Results map back: ``u = R.T @ u'``, ``sigma = R.T @ sigma' @ R``.

Rotation invariance (isotropic linear elasticity, no body forces) guarantees the
mapped-back solution equals a direct global solve; the Kratos test pins this.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OBBFrame:
    """Rigid OBB frame: ``x_local = R @ (x_global - origin)``.

    Attributes
    ----------
    origin:
        Frame origin ``x0`` in global coordinates ``(3,)``.
    R:
        Rotation matrix ``(3, 3)`` with orthonormal rows (= principal axes);
        ``R.T`` maps local vectors back to global.
    extents:
        Full local-frame sizes ``(Lx, Ly, Lz)`` (max-min along each axis).
    """

    origin: np.ndarray
    R: np.ndarray
    extents: np.ndarray


def compute_obb_frame(points: np.ndarray, method: str = "pca") -> OBBFrame:
    """Compute the OBB rigid frame of a ``(N, 3)`` point cloud.

    ``method="pca"`` (only method): principal axes from SVD of the centered
    cloud. Degenerate (rank < 3) clouds raise ``ValueError``.
    """
    if method != "pca":
        raise ValueError(f"unknown method={method!r}; expected 'pca'")
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 4:
        raise ValueError(f"points must have shape (N>=4, 3), got {pts.shape}")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")
    origin = pts.mean(axis=0)
    _, singular, vt = np.linalg.svd(pts - origin, full_matrices=False)
    if float(singular[2]) <= 1e-12 * float(singular[0]):
        raise ValueError("degenerate point cloud: rank < 3, OBB frame undefined")
    R = np.ascontiguousarray(vt)
    if float(np.linalg.det(R)) < 0.0:
        R = R.copy()
        R[2, :] = -R[2, :]
    q = (pts - origin) @ R.T
    extents = q.max(axis=0) - q.min(axis=0)
    return OBBFrame(origin=origin, R=R, extents=extents)


def to_local(points: np.ndarray, frame: OBBFrame) -> np.ndarray:
    """Map global points ``(..., 3)`` into the OBB local frame."""
    pts = np.asarray(points, dtype=np.float64)
    return (pts - frame.origin) @ frame.R.T


def to_global_points(points: np.ndarray, frame: OBBFrame) -> np.ndarray:
    """Map local-frame points ``(..., 3)`` back to global coordinates."""
    pts = np.asarray(points, dtype=np.float64)
    return pts @ frame.R + frame.origin


def to_global_vecs(vecs: np.ndarray, frame: OBBFrame) -> np.ndarray:
    """Rotate local-frame vectors (displacement, traction, normals) to global."""
    v = np.asarray(vecs, dtype=np.float64)
    return v @ frame.R


def to_local_vecs(vecs: np.ndarray, frame: OBBFrame) -> np.ndarray:
    """Rotate global vectors into the OBB local frame."""
    v = np.asarray(vecs, dtype=np.float64)
    return v @ frame.R.T


def to_global_tensors(sigmas: np.ndarray, frame: OBBFrame) -> np.ndarray:
    """Rotate local-frame symmetric tensors ``(..., 3, 3)`` back to global."""
    s = np.asarray(sigmas, dtype=np.float64)
    return frame.R.T @ s @ frame.R if s.ndim == 2 else np.einsum("ij,...jk,lk->...il", frame.R.T, s, frame.R.T)


def local_aabb(frame: OBBFrame) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Return the local-frame ``(lo, hi)`` box corners of the OBB extents."""
    half = 0.5 * np.asarray(frame.extents, dtype=np.float64)
    return tuple((-half).tolist()), tuple(half.tolist())  # type: ignore[return-value]


def rotation_about_z(angle_rad: float) -> np.ndarray:
    """Right-handed rotation matrix for angle about +z (test helper)."""
    c, s = float(np.cos(angle_rad)), float(np.sin(angle_rad))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
