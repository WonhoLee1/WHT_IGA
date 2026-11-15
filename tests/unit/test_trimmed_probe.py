# -*- coding: utf-8 -*-
"""Kratos-free tests for trimmed.py probes: Voigt order (F2), shear (F1),
face detection + explicit face id (F6), shared V02 BC helper (T1)."""

from __future__ import annotations

import numpy as np

from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T


def _box() -> BoxSpline:
    return BoxSpline(lo=(-100.0, -20.0, -0.5), hi=(100.0, 20.0, 0.5),
                     n_cells=(20, 4, 1), degree=2)


def test_voigt_order_matches_kratos() -> None:
    # Kratos order [xx,yy,zz,xy,yz,xz]; tensor must equal analytic eigendecomp.
    v = np.array([3.0, 1.0, 0.5, 0.4, -0.2, 0.7])
    s = T._tensor_from_voigt(v)
    assert s.shape == (3, 3)
    np.testing.assert_allclose(s, s.T, rtol=0, atol=0)
    w, _ = np.linalg.eigh(s)
    assert w[-1] > v[0]  # shear raises the max principal
    # pure shear: principals are exactly +/-tau
    t = T._tensor_from_voigt(np.array([0.0, 0.0, 0.0, 2.5, 0.0, 0.0]))
    np.testing.assert_allclose(sorted(np.linalg.eigvalsh(t)), [-2.5, 0.0, 2.5],
                               rtol=1e-12)


def test_pure_shear_strain_and_stress() -> None:
    # u = (gamma*y, 0, 0): engineering shear strain gamma_xy = gamma.
    box = _box()
    X = box.control_points()
    gamma = 1e-3
    U = np.column_stack([gamma * X[:, 1], np.zeros(len(X)), np.zeros(len(X))]).ravel()
    pts = np.array([[0.0, 0.0, 0.0], [33.3, -7.7, 0.2], [-88.8, 11.1, -0.4]])
    S = T._strain_at_points(box, U, pts)
    np.testing.assert_allclose(S[:, 3], gamma, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(S[:, [0, 1, 2, 4, 5]], 0.0, atol=1e-12)
    # E=2, nu=0 -> mu = E/2 = 1.0; sigma_xy = mu*gamma = 1e-3.
    D = T._d_matrix(2.0, 0.0)
    for k in range(S.shape[0]):
        sig = D @ S[k]
        # Kratos Voigt slot check: shear lives in [xy]=3, not in yz/xz slots.
        np.testing.assert_allclose(sig, [0.0, 0.0, 0.0, 1e-3, 0.0, 0.0],
                                   rtol=1e-12, atol=1e-15)
        w = np.linalg.eigvalsh(T._tensor_from_voigt(sig))
        np.testing.assert_allclose(w, [-1e-3, 0.0, 1e-3], rtol=1e-9)


def test_uniaxial_tension_recovers_sigma() -> None:
    box = _box()
    X = box.control_points()
    s0 = 2e-6
    U = np.column_stack([s0 * X[:, 0], np.zeros(len(X)), np.zeros(len(X))]).ravel()
    pts = np.array([[10.0, 0.0, 0.0], [-40.0, 12.0, 0.4]])
    S = T._strain_at_points(box, U, pts)
    np.testing.assert_allclose(S[:, 0], s0, rtol=1e-12)
    np.testing.assert_allclose(S[:, 1:], 0.0, atol=1e-14)


def test_hole_face_detection_and_explicit_id() -> None:
    from wht_iga.geometry import cad as G

    shape = G.read_step(str(__import__("pathlib").Path(__file__).resolve().parents[2] / "data" / "roark_plate_hole.step"))
    tess = G.tessellate(shape, deflection=0.25)
    fid = T.hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids,
                                       radius=10.0)
    cents = tess.nodes[tess.tris[tess.face_ids == fid]].mean(axis=1)
    rr = np.hypot(cents[:, 0], cents[:, 1])
    assert np.all(np.abs(rr - 10.0) < 1.0), (rr.min(), rr.max())
    try:
        T.hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids,
                                     radius=1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("nonsense radius must not match")
    box = _box()
    U = np.zeros(box.n_dofs)
    a = T.boundary_peak_stress(box, U, 1.0, 0.0, tess.nodes, tess.tris,
                               tess.face_ids, face_id=fid)
    assert a["face_id"] == fid and a["peak"] == 0.0
    try:
        T.boundary_peak_stress(box, U, 1.0, 0.0, tess.nodes, tess.tris,
                               tess.face_ids)
    except ValueError:
        pass
    else:
        raise AssertionError("face_id must be required")


def test_roark_bc_helper_shape() -> None:
    from wht_iga.solver import bc_roark as B

    spec = B.roark_left_pinned()
    assert spec["fixed_faces"] == [(0, 0)]
    assert len(spec["pressures"]) == 1 and spec["pressures"][0].face == (0, 1)
    assert spec["pressures"][0].magnitude == 1e-6
