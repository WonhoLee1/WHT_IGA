"""S09: penalty Dirichlet BC on boundary QP geometry.

Question (v2 §5.8): does IgaApp ``SupportPenaltyCondition`` impose
``∫ β (u − g)·v dA`` on surface QP geometry?

Reference: unit cube (p=2, E=1, nu=0), bottom face fixed, top face driven to
``g = (0.1, 0, 0)`` with β=1e10. Penalty theory predicts the face-average
displacement matches g to O(E/(β·h)); asserted here to 1e-6 relative.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature, triangle_quadrature
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


@pytest.mark.parametrize("order,n_pts", [(4, 6), (6, 12)])
def test_triangle_rule_exactness(order: int, n_pts: int) -> None:
    """Self-check of the Dunavant rules: exact on x^i y^j, i+j<=order.

    Analytic: ∫∫_ref x^i y^j dx dy = i! j! / (i+j+2)! (ref area 1/2).
    """
    import math

    pts, wts = triangle_quadrature(order)
    assert len(pts) == n_pts
    assert np.all(wts > 0.0)
    # Constants are quoted to ~1e-11..1e-15; 1e-9 still discriminates a wrong rule.
    np.testing.assert_allclose(float(wts.sum()), 0.5, rtol=0.0, atol=1e-9)
    for i in range(order + 1):
        for j in range(order + 1 - i):
            approx = float(sum(w * (x**i) * (y**j) for (x, y), w in zip(pts, wts)))
            exact = math.factorial(i) * math.factorial(j) / math.factorial(i + j + 2)
            np.testing.assert_allclose(approx, exact, rtol=0.0, atol=1e-9)


def _top_face_triangles(n_split: int = 1) -> list[np.ndarray]:
    """Top face z=1 as 2 * n_split**2 triangles (n_split x n_split squares, each cut in two)."""
    tris = []
    h = 1.0 / n_split
    for i in range(n_split):
        for j in range(n_split):
            x0, y0 = i * h, j * h
            tris.append(np.array([[x0, y0, 1.0], [x0 + h, y0, 1.0], [x0 + h, y0 + h, 1.0]]))
            tris.append(np.array([[x0, y0, 1.0], [x0 + h, y0 + h, 1.0], [x0, y0 + h, 1.0]]))
    return tris


def _penalty_top_face_error(degree: int, rule_order: int, beta: float, n_split: int = 1) -> tuple[float, float]:
    """(relative error of the top-face CP average u_x against g_x = 0.1, max |u_y| on top CPs).

    Unit cube, 1 cell, E=1, nu=0, bottom fixed, top face penalty to g=(0.1,0,0).
    The exact solution (simple shear u_x = 0.1 z) is in the spline space.
    """
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S09")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=degree)
    volume, _ = make_nurbs_volume(KM, mp, box)
    vol_props = mp.CreateNewProperties(1)
    vol_props.SetValue(KM.YOUNG_MODULUS, 1.0)
    vol_props.SetValue(KM.POISSON_RATIO, 0.0)
    vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    qp_vol = KM.GeometriesVector()
    volume.CreateQuadraturePointGeometries(qp_vol, 2, cell_gauss_points(box, 0, 0, 0))
    mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qp_vol[0], vol_props)
    for node in mp.Nodes:
        if abs(node.Z0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)

    qp_geoms, wrappers, _, _ = make_surface_quadrature(KM, mp, volume, _top_face_triangles(n_split), rule_order=rule_order)
    volume_keepalive = volume
    wrappers_keepalive = wrappers
    assert volume_keepalive is not None and wrappers_keepalive is not None
    pen_props = mp.CreateNewProperties(2)
    pen_props.SetValue(IGA.PENALTY_FACTOR, beta)
    g = KM.Vector(3)
    g[0], g[1], g[2] = 0.1, 0.0, 0.0
    for i, qp in enumerate(qp_geoms):
        cond = mp.CreateNewCondition("SupportPenaltyCondition", 100 + i, qp, pen_props)
        cond.SetValue(KM.DISPLACEMENT, g)

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()

    top_ids = set(int(1 + I) for I in box.boundary_control_points(axis=2, side=1))
    ux = [node.GetSolutionStepValue(KM.DISPLACEMENT_X) for node in mp.Nodes if node.Id in top_ids]
    uy = [node.GetSolutionStepValue(KM.DISPLACEMENT_Y) for node in mp.Nodes if node.Id in top_ids]
    assert len(ux) == (degree + 1) ** 2
    return abs(float(sum(ux) / len(ux)) - 0.1) / 0.1, float(max(abs(v) for v in uy))


def test_penalty_dirichlet_top_face() -> None:
    """p=2, degree-4 rule (12 points >= 9 top-face CPs): error O(E/beta)."""
    err, max_uy = _penalty_top_face_error(2, 4, 1.0e10)
    print(f"\nS09 p=2 rule4 beta=1e10: rel err {err:.3e}, max|u_y| {max_uy:.3e}")
    assert err < 1e-6
    assert max_uy < 1e-9


# Rank rule (REVIEW_P0_P1 §1 #3). Discriminating runs on the 1-cell cube
# (2026-11-14, integration lead): the beta-independent floor tracks the
# per-face-cell QP count against the (p+1)^2 face CPs, not the rule degree.
#   p=2 order-2, 6 QPs  < 9  -> floor 4.541e-4 at beta 1e8 and 1e10
#   p=2 order-2, 24 QPs >= 9 (degree 2 only) -> 3.438e-9 / 3.438e-11 (1/beta)
#   p=3 order-4, 12 QPs < 16 -> floor 1.772e-4 at beta 1e8 and 1e10
#   p=3 order-6, 12 QPs/tri = 24 >= 16 -> 3.027e-9 / 3.027e-11
#   p=3 order-2, 24 QPs >= 16 (degree 2 only) -> 3.040e-9 / 3.040e-11
@pytest.mark.parametrize(
    "degree,rule_order,n_split",
    [(2, 2, 2), (2, 4, 1), (3, 6, 1), (3, 2, 2)],
)
def test_penalty_full_rank_scales_with_beta(degree: int, rule_order: int, n_split: int) -> None:
    e8, _ = _penalty_top_face_error(degree, rule_order, 1.0e8, n_split)
    e10, uy = _penalty_top_face_error(degree, rule_order, 1.0e10, n_split)
    print(f"\nS09 full-rank p={degree} order={rule_order} split={n_split}: {e8:.3e} -> {e10:.3e}")
    assert 90.0 < e8 / e10 < 110.0
    assert e10 < 1e-9
    assert uy < 1e-9


@pytest.mark.parametrize(
    "degree,rule_order,floor",
    [(2, 2, 4.541e-4), (3, 4, 1.772e-4)],
)
def test_penalty_rank_deficient_floor(degree: int, rule_order: int, floor: float) -> None:
    e8, _ = _penalty_top_face_error(degree, rule_order, 1.0e8)
    e10, _ = _penalty_top_face_error(degree, rule_order, 1.0e10)
    print(f"\nS09 rank-deficient p={degree} order={rule_order}: {e8:.3e} / {e10:.3e}")
    np.testing.assert_allclose(e8, floor, rtol=1e-3)
    np.testing.assert_allclose(e10, floor, rtol=1e-3)
