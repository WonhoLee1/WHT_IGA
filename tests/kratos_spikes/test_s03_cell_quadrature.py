"""S03: one QP geometry per cell — count, weights, partition of unity, support.

Reference: tensor-product Gauss structure (27 points, weight sum = cell
volume) and the B-spline support of cell (i, j, k) = (p+1)^3 control points.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

P = 2
BOX = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=P)


def _model(KM):
    model = KM.Model()
    mp = model.CreateModelPart("S03")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def _support_ids(box: BoxSpline, i: int, j: int, k: int) -> set[int]:
    nu, nv, _ = box.n_cp
    return {
        1 + (a + nu * (b + nv * c))
        for a in range(i, i + box.degree + 1)
        for b in range(j, j + box.degree + 1)
        for c in range(k, k + box.degree + 1)
    }


def test_gauss_cell_single_geometry_conserves_weight():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, BOX)
    pts = cell_gauss_points(BOX, 0, 0, 0)
    assert len(pts) == 27
    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, pts)
    assert len(qps) == 1, "K2: one cell -> exactly one QP geometry"
    assert qps[0].IntegrationPointsNumber() == 27
    assert qps[0].PointsNumber() == 27
    assert abs(sum(p[3] for p in pts) - 0.125) < 1e-15, "cell volume 0.5^3"
    sfv = np.array(qps[0].ShapeFunctionsValues())
    assert sfv.shape == (27, 27)
    np.testing.assert_allclose(sfv.sum(axis=1), np.ones(27), atol=1e-14)
    assert {int(nd.Id) for nd in qps[0]} == _support_ids(BOX, 0, 0, 0)
    print("S03 gauss: rowsums ok, nodes=27, wsum=0.125")


def test_arbitrary_points_single_geometry_partition_of_unity():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, BOX)
    rng = np.random.default_rng(7)
    xyz = rng.uniform(0.02, 0.48, size=(7, 3))
    w = rng.uniform(0.5, 1.5, size=7)
    pts = [[float(x), float(y), float(z), float(wi)] for (x, y, z), wi in zip(xyz, w)]
    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, pts)
    assert len(qps) == 1
    assert qps[0].IntegrationPointsNumber() == 7
    sfv = np.array(qps[0].ShapeFunctionsValues())
    np.testing.assert_allclose(sfv.sum(axis=1), np.ones(7), atol=1e-14)
    assert {int(nd.Id) for nd in qps[0]} == _support_ids(BOX, 0, 0, 0)


def test_boundary_point_support_shift_recorded():
    """A point exactly on a cell face: record which span owns it (K3 evidence)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, BOX)
    on_face = [[0.5, 0.25, 0.25, 1.0]]  # x = 0.5 is the face between cells 0 and 1
    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, on_face)
    got = {int(nd.Id) for nd in qps[0]}
    in_left = got == _support_ids(BOX, 0, 0, 0)
    in_right = got == _support_ids(BOX, 1, 0, 0)
    assert in_left or in_right, sorted(got)
    print(f"face point owned by {'left' if in_left else 'right'} cell support")
