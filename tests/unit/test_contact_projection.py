"""Unit tests for DeformableBoxFaceMaster and _project_box_face_njit."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses.contact import DeformableBoxFaceMaster
from wht_iga.reference.box import BoxSpline


def test_box_face_master_all_six_orientations() -> None:
    """Outward normals and gap signs are exact on all 6 faces of a BoxSpline."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 3.0, 4.0), n_cells=(2, 2, 2), degree=2)
    U0 = np.zeros((box.n_ctrl, 3), dtype=np.float64)
    center = np.array([1.0, 1.5, 2.0])

    expected_normals = {
        (0, 0): np.array([-1.0, 0.0, 0.0]),
        (0, 1): np.array([1.0, 0.0, 0.0]),
        (1, 0): np.array([0.0, -1.0, 0.0]),
        (1, 1): np.array([0.0, 1.0, 0.0]),
        (2, 0): np.array([0.0, 0.0, -1.0]),
        (2, 1): np.array([0.0, 0.0, 1.0]),
    }

    for (axis, side), n_exp in expected_normals.items():
        face = DeformableBoxFaceMaster(box, axis=axis, side=side)
        # Point 0.05 outside the face along n_exp, and 0.05 inside (penetration)
        face_pt = center.copy()
        face_pt[axis] = box.lo[axis] if side == 0 else box.hi[axis]
        x_s = np.vstack([face_pt + 0.05 * n_exp, face_pt - 0.05 * n_exp])

        g_n, n_m, x_m, ab_m, B_m = face.evaluate(x_s, U0)
        np.testing.assert_allclose(n_m[0], n_exp, atol=1e-14)
        np.testing.assert_allclose(n_m[1], n_exp, atol=1e-14)
        assert g_n[0] == pytest.approx(0.05, abs=1e-12)
        assert g_n[1] == pytest.approx(-0.05, abs=1e-12)
        np.testing.assert_allclose(x_m, np.vstack([face_pt, face_pt]), atol=1e-12)
        np.testing.assert_allclose(B_m.sum(axis=1).A.ravel(), 1.0, atol=1e-13)
        X_def = box.control_points() + U0
        np.testing.assert_allclose(B_m @ X_def, x_m, atol=1e-13)


def test_box_face_master_sheared_interior_and_boundary_clamping() -> None:
    """Closest-point projection under shear and tangential boundary clamping."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(3, 3, 3), degree=2)
    face = DeformableBoxFaceMaster(box, axis=1, side=1)

    # Apply a smooth shear + normal tilt to the master control points:
    # u_x = 0.25 * y * (1 + 0.2 * z), u_y = 0.1 * x, u_z = 0.15 * x * y
    X0 = box.control_points()
    U_m = np.column_stack([
        0.25 * X0[:, 1] * (1.0 + 0.2 * X0[:, 2]),
        0.10 * X0[:, 0],
        0.15 * X0[:, 0] * X0[:, 1],
    ])
    X_def = X0 + U_m

    # 1. Interior point: should project strictly inside (lo, hi) with zero tangential residual
    x_int = np.array([[0.65, 1.02, 0.55]])
    g_int, n_int, xm_int, ab_int, Bm_int = face.evaluate(x_int, U_m)
    assert 0.0 < ab_int[0, 0] < 1.0
    assert 0.0 < ab_int[0, 1] < 1.0
    np.testing.assert_allclose(Bm_int @ X_def, xm_int, atol=1e-12)
    # Closest-point orthogonality: (x_s - x_m) must be parallel to n_m
    diff_int = x_int[0] - xm_int[0]
    np.testing.assert_allclose(diff_int, g_int[0] * n_int[0], atol=1e-10)

    # 2. Tangential overhang beyond hi_a=1.0 and below lo_b=0.0:
    #    ab_m must clamp to [lo, hi], B_m remains a valid partition of unity,
    #    and g_n equals the signed normal projection (x_s - x_m) . n_m at the clamped edge.
    x_over = np.array([
        [1.80, 1.05, 0.50],   # past hi_a, slightly penetrating along n_m
        [-0.40, 1.20, -0.30], # past lo_a and lo_b, separated
    ])
    g_ov, n_ov, xm_ov, ab_ov, Bm_ov = face.evaluate(x_over, U_m)

    assert ab_ov[0, 0] == pytest.approx(1.0, abs=1e-14)
    assert 0.0 <= ab_ov[0, 1] <= 1.0
    assert ab_ov[1, 0] == pytest.approx(0.0, abs=1e-14)
    assert ab_ov[1, 1] == pytest.approx(0.0, abs=1e-14)

    np.testing.assert_allclose(Bm_ov.sum(axis=1).A.ravel(), 1.0, atol=1e-13)
    np.testing.assert_allclose(Bm_ov @ X_def, xm_ov, atol=1e-12)
    for q in range(2):
        assert g_ov[q] == pytest.approx(float(np.dot(x_over[q] - xm_ov[q], n_ov[q])), abs=1e-12)
