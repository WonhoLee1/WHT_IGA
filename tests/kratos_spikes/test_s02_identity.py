"""S02: identity map — physical knots + Greville nodes give x(xi) = xi, detJ = 1.

Reference: analytic identity (partition of unity + linear precision of the
uniform B-spline under the v2 §5.2.1 identity trick).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.grid.knots import to_kratos
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _model(KM):
    model = KM.Model()
    mp = model.CreateModelPart("S02")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def test_global_coordinates_is_identity():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    rng = np.random.default_rng(0)
    pts = rng.uniform(0.05, 0.95, size=(10, 3))
    errs = []
    for p in pts:
        g = vol.GlobalCoordinates(KM.Vector([float(p[0]), float(p[1]), float(p[2])]))
        errs.append(float(np.max(np.abs(np.array([g[0], g[1], g[2]]) - p))))
    print("S02 identity max err:", max(errs))
    assert max(errs) < 1e-14, max(errs)


def test_quadrature_geometry_detj_is_one():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    worst = 0.0
    for k in range(2):
        for j in range(2):
            for i in range(2):
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, j, k))
                for q in range(qps[0].IntegrationPointsNumber()):
                    worst = max(worst, abs(float(qps[0].DeterminantOfJacobian(q)) - 1.0))
    print("S02 detJ worst dev:", worst)
    assert worst < 1e-14, worst


def test_identity_under_frame_rotation():
    """Rotate all Greville nodes by 30 deg about (1,1,1): G(xi) == R xi, detJ == 1."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    _, mp = _model(KM)
    theta = np.deg2rad(30.0)
    axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) * np.cos(theta) + np.sin(theta) * K + (1 - np.cos(theta)) * np.outer(axis, axis)
    nodes = KM.NodesVector()
    for I, gp in enumerate(box.control_points()):
        x, y, z = R @ gp
        nodes.append(mp.CreateNewNode(1 + I, float(x), float(y), float(z)))
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    p = box.degree
    vol = KM.NurbsVolumeGeometry(
        nodes, p, p, p,
        KM.Vector(to_kratos(box.knots[0], p).tolist()),
        KM.Vector(to_kratos(box.knots[1], p).tolist()),
        KM.Vector(to_kratos(box.knots[2], p).tolist()),
    )
    rng = np.random.default_rng(1)
    worst_pos, worst_det = 0.0, 0.0
    for q in rng.uniform(0.05, 0.95, size=(10, 3)):
        g = vol.GlobalCoordinates(KM.Vector([float(q[0]), float(q[1]), float(q[2])]))
        worst_pos = max(worst_pos, float(np.max(np.abs(np.array([g[0], g[1], g[2]]) - R @ q))))
    for k in range(2):
        for j in range(2):
            for i in range(2):
                pts = cell_gauss_points(box, i, j, k)
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                for n in range(qps[0].IntegrationPointsNumber()):
                    worst_det = max(worst_det, abs(float(qps[0].DeterminantOfJacobian(n)) - 1.0))
    assert worst_pos < 1e-14, worst_pos
    assert worst_det < 1e-14, worst_det
    print("S02 rotated worst_pos:", worst_pos, "worst_det:", worst_det)
