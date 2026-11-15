"""S07: surface quadrature-point geometry on an immersed volume boundary.

Question (v2 §5.8): does ``SurfaceInNurbsVolumeGeometry(volume, Triangle3D3)``
+ user quadrature give ``DeterminantOfJacobian`` == area element and ``Normal``
matching the triangle orientation?

Reference: analytic face area (unit cube, face z=1, area 1.0) and analytic
outward normal (0, 0, 1). Tolerance 1e-13 on areas, exact sign recorded.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _top_face_triangles() -> list[np.ndarray]:
    return [
        np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]),
        np.array([[0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]]),
    ]


@pytest.mark.parametrize("degree", [1, 2])
def test_surface_qp_area_and_normal(degree: int) -> None:
    from wht_iga.kratos_backend import load_kratos

    KM = load_kratos()["KratosMultiphysics"]
    model = KM.Model()
    mp = model.CreateModelPart("S07")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=degree)
    volume, _ = make_nurbs_volume(KM, mp, box)

    qp_geoms, _, _, _ = make_surface_quadrature(KM, mp, volume, _top_face_triangles(), rule_order=2)

    from wht_iga.kratos_backend.surface import triangle_quadrature

    _, ref_w = triangle_quadrature(2)
    assert len(qp_geoms) == 6  # one QP geometry per integration point (see surface.py)
    total_area = 0.0
    for qp, w in zip(qp_geoms, list(ref_w) * 2):
        assert qp.IntegrationPointsNumber() == 1
        detj = qp.DeterminantOfJacobian(0)
        total_area += float(w) * detj
        n = qp.Normal(0)
        assert abs(float(n[0])) < 1e-12
        assert abs(float(n[1])) < 1e-12
        assert abs(abs(float(n[2])) - 1.0) < 1e-12
    np.testing.assert_allclose(total_area, 1.0, rtol=0.0, atol=1e-13)
    print(f"\nS07 p={degree}: face area sum = {total_area:.15f}, normal_z sign recorded in S07.md")


def test_surface_qp_normal_sign_matches_winding() -> None:
    """Outward (+z) winding must give +z normals; flipped winding must flip them."""
    from wht_iga.kratos_backend import load_kratos

    KM = load_kratos()["KratosMultiphysics"]
    model = KM.Model()
    mp = model.CreateModelPart("S07sign")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=2)
    volume, _ = make_nurbs_volume(KM, mp, box)

    outward = _top_face_triangles()
    inward = [t[::-1].copy() for t in outward]
    qp_out, _, _, _ = make_surface_quadrature(KM, mp, volume, outward, rule_order=1, first_node_id=20_000)
    qp_in, _, _, _ = make_surface_quadrature(KM, mp, volume, inward, rule_order=1, first_node_id=30_000)
    assert len(qp_out) == 2 and len(qp_in) == 2  # 1-pt rule × 2 triangles

    n_out = qp_out[0].Normal(0)
    n_in = qp_in[0].Normal(0)
    print(f"\nS07 winding: outward n = {[float(v) for v in n_out]}, flipped n = {[float(v) for v in n_in]}")
    assert float(n_out[2]) > 0.0
    assert float(n_in[2]) < 0.0
