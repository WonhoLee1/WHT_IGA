# -*- coding: utf-8 -*-
"""OBB background-grid frame adapter + rotation-invariance Kratos gate (Track H).

A bar rotated 30 deg about z is solved in its PCA local frame with the
existing axis-aligned pipeline; the mapped-back tip must equal the analytic
global tip (isotropic elasticity is rotation-invariant).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses.methods import (
    BoxMethodProblem,
    FullIgaOption,
    solve_with_discretization_method,
)
from wht_iga.reference.obb import (
    compute_obb_frame,
    local_aabb,
    rotation_about_z,
    to_global_vecs,
    to_local,
)

pytestmark = pytest.mark.kratos


def _rotated_bar_corners(L: float, h: float, R: np.ndarray) -> np.ndarray:
    base = np.array(
        [[x, y, z] for x in (-L / 2, L / 2) for y in (-h / 2, h / 2) for z in (-h / 2, h / 2)],
        dtype=np.float64,
    )
    return base @ R.T


def test_obb_frame_recovers_principal_axes():
    """PCA frame of a rotated box recovers the rotation (up to axis sign)."""
    L, h = 10.0, 1.0
    R = rotation_about_z(np.deg2rad(30.0))
    corners = _rotated_bar_corners(L, h, R)
    frame = compute_obb_frame(corners)
    assert frame.extents is not None
    assert sorted(frame.extents.tolist()) == pytest.approx(sorted([L, h, h]), rel=1e-9)
    # recovered rows must match the true body axes (columns of R) up to sign/order
    body_axes = [R[:, i] for i in range(3)]
    for a in body_axes:
        dots = [abs(float(a @ frame.R[j, :])) for j in range(3)]
        assert max(dots) == pytest.approx(1.0, abs=1e-9)


def test_obb_local_solve_maps_to_global_analytic():
    """Local-frame solve mapped back matches the analytic global tip."""
    E_mod, L, h, t = 2.0e5, 10.0, 1.0, 2.0
    R = rotation_about_z(np.deg2rad(30.0))
    corners = _rotated_bar_corners(L, h, R)
    frame = compute_obb_frame(corners)
    assert np.allclose(sorted(frame.extents.tolist()), sorted([L, h, h]), rtol=1e-9)

    lo, hi = local_aabb(frame)
    problem = BoxMethodProblem(
        lo=lo,
        hi=hi,
        E=E_mod,
        nu=0.0,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(t, 0.0, 0.0),
    )
    sol = solve_with_discretization_method(
        problem, FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(4, 2, 2))
    )
    assert sol.converged is True
    # local tip may sit on any +/- x end depending on PCA sign flips: use |x|
    tip_local = np.asarray(sol.tip_displacement, dtype=float)
    L_fit = float(frame.extents[0])
    assert abs(abs(tip_local[0]) - t * L_fit / E_mod) / (t * L_fit / E_mod) < 1e-5
    tip_global = to_global_vecs(tip_local, frame)
    # length axis in global = recovered row with the largest extent
    long_ax = np.asarray(frame.R[int(np.argmax(frame.extents)), :], dtype=float)
    assert abs(float(np.linalg.norm(tip_global)) - t * L_fit / E_mod) / (t * L_fit / E_mod) < 1e-5
    assert float(np.linalg.norm(np.cross(tip_global, long_ax))) / (t * L_fit / E_mod) < 1e-6
    # local frame really is rotated: mapped tip must differ from the naive axis answer
    assert not np.allclose(tip_global, np.array([t * L / E_mod, 0.0, 0.0]), rtol=1e-3)
