"""Reference assembler acceptance tests (implementation request T-P0.5)."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.reference.assemble_elastic import assemble_face_traction, assemble_stiffness
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.solve import solve_with_dirichlet


def boundary_dofs(box: BoxSpline) -> np.ndarray:
    cps = np.unique(np.concatenate([box.boundary_control_points(a, s) for a in range(3) for s in (0, 1)]))
    return np.sort(np.concatenate([3 * cps + d for d in range(3)]))


@pytest.mark.parametrize("p,k", [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)])
def test_patch_test_linear_field(p, k):
    """u = (x, 2y, 3z), E=1, nu=0 (Coreform prescribed-displacement cube field)."""
    box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (3, 3, 3), p, k)
    K = assemble_stiffness(box, 1.0, 0.0)
    X = box.control_points()
    exact = np.column_stack([X[:, 0], 2 * X[:, 1], 3 * X[:, 2]])  # exact: identity map + linear field
    fixed = boundary_dofs(box)
    u = solve_with_dirichlet(K, np.zeros(box.n_dofs), fixed, exact.ravel()[fixed])
    np.testing.assert_allclose(u.reshape(-1, 3), exact, atol=1e-10)

    pts = np.random.default_rng(1).uniform(0, 1, (25, 3))
    np.testing.assert_allclose(box.evaluate(u.reshape(-1, 3), pts), pts * [1, 2, 3], atol=1e-10)


def test_stiffness_symmetric_and_rigid_modes():
    box = BoxSpline((0, 0, 0), (2, 1, 1), (2, 2, 2), 2)
    K = assemble_stiffness(box, 200.0, 0.3)
    assert abs(K - K.T).max() < 1e-10 * abs(K).max()
    X = box.control_points()
    modes = [np.tile(e, (box.n_ctrl, 1)) for e in np.eye(3)]
    modes.append(np.column_stack([-X[:, 1], X[:, 0], 0 * X[:, 0]]))  # rotation about z
    for m in modes:
        assert np.abs(K @ m.ravel()).max() < 1e-9 * abs(K).max()


@pytest.mark.parametrize("p", [2, 3])
def test_cantilever_tip_deflection(p):
    """Clamped at x=0, shear traction on x=L. Timoshenko reference, nu=0."""
    L, h, E, P = 10.0, 1.0, 1000.0, 1.0
    box = BoxSpline((0, 0, 0), (L, h, h), (40, 4, 4), p)
    K = assemble_stiffness(box, E, 0.0)
    f = assemble_face_traction(box, 0, 1, (0.0, 0.0, -P / (h * h)))
    clamp = box.boundary_control_points(0, 0)
    fixed = np.sort(np.concatenate([3 * clamp + d for d in range(3)]))
    u = solve_with_dirichlet(K, f, fixed, np.zeros(len(fixed)))
    tip = box.evaluate(u.reshape(-1, 3), np.array([[L, 0.5 * h, 0.5 * h]]))[0, 2]
    I = h**4 / 12
    G, kappa = E / 2, 5 / 6
    ref = -(P * L**3 / (3 * E * I) + P * L / (kappa * G * h * h))
    assert abs(tip - ref) / abs(ref) < 0.01, (tip, ref)
