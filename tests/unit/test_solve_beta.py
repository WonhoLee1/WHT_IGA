"""Unit tests for solver.vtu + penalty helpers (Kratos-free)."""

from __future__ import annotations

import numpy as np

from wht_iga.solver.penalty import analyze_sweep, box_face_triangles, default_beta, eval_g, rule_for_degree
from wht_iga.solver.validate import DirichletBC


def test_face_triangles_outward() -> None:
    from wht_iga.kratos_backend.surface import triangle_area, triangle_normal

    lo = np.zeros(3)
    hi = np.ones(3)
    for axis in range(3):
        for side in (0, 1):
            tris = box_face_triangles(lo, hi, axis, side)
            assert len(tris) == 2
            outward = np.zeros(3)
            outward[axis] = 1.0 if side == 1 else -1.0
            for t in tris:
                n = triangle_normal(t)
                assert float(n @ outward) > 0.999, (axis, side, n)
            assert abs(sum(triangle_area(t) for t in tris) - 1.0) < 1e-15


def test_face_triangles_subdiv_exact_area() -> None:
    from wht_iga.kratos_backend.surface import triangle_area, triangle_normal

    lo, hi = np.zeros(3), np.ones(3)
    br = [np.array([0.0, 0.5, 1.0]), np.array([0.0, 0.5, 1.0])]
    tris = box_face_triangles(lo, hi, 2, 1, subdiv=br)
    assert len(tris) == 8  # 2x2 face cells x 2
    assert abs(sum(triangle_area(t) for t in tris) - 1.0) < 1e-15
    for t in tris:
        assert float(triangle_normal(t) @ np.array([0.0, 0.0, 1.0])) > 0.999


def test_default_beta_formula() -> None:
    from wht_iga.solver.penalty import DEFAULT_C_PEN

    assert DEFAULT_C_PEN == 3.0e5  # ADR-002
    assert default_beta(1.0, 0.5, 2) == 3.0e5 * 1.0 * 4 / 0.5
    assert default_beta(2.0e5, 0.1, 1, 1.0e3) == 1.0e3 * 2.0e5 / 0.1


def test_rule_for_degree() -> None:
    from wht_iga.kratos_backend.surface import triangle_quadrature

    assert rule_for_degree(1) == 2
    assert rule_for_degree(2) == 4
    # S09 correction: rank needs per-face-cell points >= (p+1)^2, so p=3
    # needs order 6 (24 pts/face-cell >= 16); order 4 floors at 1.8e-4.
    assert rule_for_degree(3) == 6
    for p, need in ((1, 4), (2, 9), (3, 16)):
        npts = 2 * len(triangle_quadrature(rule_for_degree(p))[0])
        assert npts >= need, (p, npts, need)


def test_eval_g_affine() -> None:
    from wht_iga.solver.validate import AffineField

    bc = DirichletBC(face=(0, 0), field=AffineField(A=((1, 0, 0), (0, 2, 0), (0, 0, 3)), b=(0, 0, 0)))
    np.testing.assert_allclose(eval_g(bc, np.array([0.5, 0.5, 0.5])), [0.5, 1.0, 1.5])


def test_analyze_sweep_plateau() -> None:
    betas = np.array([1e2, 1e3, 1e4, 1e5, 1e6])
    errors = 3.0 / betas + 1e-9  # O(1/beta) + discretization floor
    r = analyze_sweep(betas, errors)
    assert r["plateau_start"] == 4
    assert abs(r["pre_plateau_slope"] + 1.0) < 0.05  # slope -1 in log-log
    assert r["recommended"] == 1e7
