"""S3 classification tests vs analytic box/sphere references (v2 §7.3)."""

from __future__ import annotations

import numpy as np

from wht_iga.geometry import cad, classify

V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
              [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
T = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
              [0, 1, 5], [0, 5, 4], [2, 3, 7], [2, 7, 6],
              [0, 4, 7], [0, 7, 3], [1, 2, 6], [1, 6, 5]])


def test_box_states_exact() -> None:
    """Non-aligned grid: IN/OUT/CUT match analytic box membership."""
    lo = np.array([-0.3, -0.3, -0.3])
    h = np.array([0.4, 0.4, 0.4])
    n = (4, 4, 4)
    res = classify.classify_cells(V, T, lo, h, n)
    assert not res["vertex_disagree"].any()  # clean case: 3 axes agree
    for k in range(4):
        for j in range(4):
            for i in range(4):
                blo = lo + h * np.array([i, j, k])
                bhi = blo + h
                c = i + 4 * (j + 4 * k)
                if (blo > 0).all() and (bhi < 1).all():
                    assert res["state"][c] == classify.IN, (i, j, k)
                elif (bhi <= 0).any() or (blo >= 1).any():
                    assert res["state"][c] == classify.OUT, (i, j, k)
                else:
                    assert res["state"][c] == classify.CUT, (i, j, k)


def test_box_vertices_match_analytic() -> None:
    """Grid-vertex flags vs strict box membership (margin keeps it exact)."""
    lo = np.array([-0.3, -0.3, -0.3])
    h = np.array([0.4, 0.4, 0.4])
    n = (4, 4, 4)
    inside, _ = classify.classify_vertices(V, T, lo, h, n)
    axes = [np.linspace(lo[d], lo[d] + h[d] * n[d], n[d] + 1) for d in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    verts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    exact = ((verts > 0.05) & (verts < 0.95)).all(axis=1)
    far = ((verts < -0.05) | (verts > 1.05)).any(axis=1)
    np.testing.assert_array_equal(inside[exact], True)
    np.testing.assert_array_equal(inside[far], False)


def _sphere_soup(r: float, deflection: float):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
    from OCP.gp import gp_Pnt

    sph = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), r).Shape()
    tess = cad.tessellate(cad.heal(sph), deflection=deflection)
    return tess.nodes, tess.tris


def test_sphere_states_with_margin() -> None:
    """Sphere r=0.5: cells strictly inside/outside (margin 0.05) classify so."""
    nodes, tris = _sphere_soup(0.5, 0.05)
    lo = np.array([-1.0, -1.0, -1.0])
    h = np.array([0.4, 0.4, 0.4])
    n = (5, 5, 5)
    res = classify.classify_cells(nodes, tris, lo, h, n)
    checked_in, checked_out = 0, 0
    for k in range(5):
        for j in range(5):
            for i in range(5):
                blo = lo + h * np.array([i, j, k])
                corners = blo + h * np.array(
                    [[a, b, c] for a in (0, 1) for b in (0, 1) for c in (0, 1)])
                d = np.linalg.norm(corners, axis=1)
                c = i + 5 * (j + 5 * k)
                if d.max() < 0.45:
                    assert res["state"][c] == classify.IN, (i, j, k)
                    checked_in += 1
                elif d.min() > 0.55:
                    assert res["state"][c] == classify.OUT, (i, j, k)
                    checked_out += 1
    assert checked_in >= 1 and checked_out >= 1
    assert (res["state"] == classify.CUT).sum() >= 1


def test_winding_fallback_matches_analytic() -> None:
    """Generalized winding number vs sphere membership (Jacobson 2013)."""
    nodes, tris = _sphere_soup(0.5, 0.05)
    rng = np.random.default_rng(7)
    pts = rng.uniform(-0.9, 0.9, size=(40, 3))
    dist = np.linalg.norm(pts, axis=1)
    keep = np.abs(dist - 0.5) > 0.05
    got = classify.winding_inside(pts[keep], nodes, tris)
    np.testing.assert_array_equal(got, dist[keep] < 0.5)
