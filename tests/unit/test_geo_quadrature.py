"""S4 quadrature tests vs analytic refs: tensor exactness, tet cut (G1),
slab R-face path, mandatory Σw checks, fallback (v2 §7.4, §10.3)."""

from __future__ import annotations

from math import factorial

import numpy as np

from wht_iga.geometry import quadrature
from wht_iga.geometry.quadrature import (
    check_surface_rule,
    check_volume_rule,
    closed_surface_flux,
)

V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
              [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
T = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
              [0, 1, 5], [0, 5, 4], [2, 3, 7], [2, 7, 6],
              [0, 4, 7], [0, 7, 3], [1, 2, 6], [1, 6, 5]])

S_TET = 0.9
TV = np.array([[0, 0, 0], [S_TET, 0, 0], [0, S_TET, 0], [0, 0, S_TET]])
TT = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
TET_VOL = S_TET**3 / 6
TET_AREA = S_TET**2 * (3 + np.sqrt(3)) / 2


def _legendre_moment_ref(blo: np.ndarray, h: np.ndarray, q: int) -> np.ndarray:
    """Independent reference: standard-Legendre moments via analytic 1D
    monomial integrals over the box (no quadrature in the reference)."""
    L = quadrature.legendre_monomials(q)
    m1 = []
    for d in range(3):
        lo, hi = blo[d], blo[d] + h[d]
        m1.append(np.array([(hi ** (a + 1) - lo ** (a + 1)) / (a + 1)
                            for a in range(q + 1)]))
    return np.einsum("ia,jb,kc,a,b,c->ijk", L, L, L,
                     m1[0], m1[1], m1[2], optimize=True)


def test_in_cell_rule_exact() -> None:
    """Fully-inside cell: tensor rule integrates degree-≤q polys to 1e-12."""
    blo = np.array([0.2, 0.1, 0.3])
    h = np.array([0.5, 0.5, 0.5])
    for q in (1, 2):
        rule = quadrature.cell_quadrature(V, T, blo, h, q=q)
        assert not rule.fallback
        ref = _legendre_moment_ref(blo, h, q)
        got = np.einsum("p,ip,jp,kp->ijk", rule.w,
                        quadrature._polyval(quadrature.legendre_monomials(q),
                                            rule.x[:, 0]),
                        quadrature._polyval(quadrature.legendre_monomials(q),
                                            rule.x[:, 1]),
                        quadrature._polyval(quadrature.legendre_monomials(q),
                                            rule.x[:, 2]))
        scale = abs(ref).max()
        np.testing.assert_allclose(got, ref, rtol=1e-12,
                                   atol=1e-12 * scale)
        assert (rule.w >= 0).all()


def _tet_legendre_ref(q: int) -> np.ndarray:
    """Analytic tet Legendre moments: ∫x^ay^bz^c = a!b!c!s^{n+3}/(n+3)!."""
    L = quadrature.legendre_monomials(q)
    out = np.zeros((q + 1,) * 3)
    for i in range(q + 1):
        for j in range(q + 1):
            for k in range(q + 1):
                s = 0.0
                for a in range(q + 1):
                    for b in range(q + 1):
                        for c in range(q + 1):
                            n = a + b + c
                            mono = (factorial(a) * factorial(b) * factorial(c)
                                    * S_TET ** (n + 3) / factorial(n + 3))
                            s += L[i, a] * L[j, b] * L[k, c] * mono
                out[i, j, k] = s
    return out


def test_tet_cut_polynomial_accuracy() -> None:
    """G1 gate: corner-tet cut cell integrates degree-≤2 polys to 1e-12."""
    blo = np.array([-0.1, -0.1, -0.1])
    h = np.array([1.1, 1.1, 1.1])
    rule = quadrature.cell_quadrature(TV, TT, blo, h, q=2)
    assert not rule.fallback
    assert len(rule.w) <= 27 and (rule.w >= 0).all()
    ref = _tet_legendre_ref(2)
    L = quadrature.legendre_monomials(2)
    got = np.einsum("p,ip,jp,kp->ijk", rule.w,
                    quadrature._polyval(L, rule.x[:, 0]),
                    quadrature._polyval(L, rule.x[:, 1]),
                    quadrature._polyval(L, rule.x[:, 2]))
    np.testing.assert_allclose(rule.volume, TET_VOL, rtol=1e-12)
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12)
    assert check_volume_rule(rule, TET_VOL) < 1e-12


def test_slab_coplanar_r_face() -> None:
    """Slab face coplanar with x1: full-rect R path vs separable analytic.

    The cell interior sees no surface (face touches x1), so the public
    path takes the IN shortcut; the R machinery is reached here through
    ``_cell_moments`` directly (white-box).
    """
    from wht_iga.geometry import classify as _cl

    V2 = V * np.array([1.0, 2.0, 2.0])
    blo = np.array([0.5, 0.2, 0.2])
    h = np.array([0.5, 0.5, 0.5])
    soup = V2[T]
    touched = [_cl._tri_box_overlap(t[0], t[1], t[2], blo, blo + h)
               for t in soup]
    M = quadrature._cell_moments(soup[np.array(touched)], blo, h, 2,
                               whole=soup)
    # mapped Legendre is orthogonal over the cell: only M000 survives
    ref = np.zeros((3, 3, 3))
    ref[0, 0, 0] = 0.5**3
    np.testing.assert_allclose(M, ref, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(M[0, 0, 0], 0.5**3, rtol=1e-12)
    # public path: same cell resolves through the IN shortcut, same volume
    rule = quadrature.cell_quadrature(V2, T, blo, h, q=2)
    assert not rule.fallback
    np.testing.assert_allclose(rule.volume, 0.5**3, rtol=1e-12)


def test_tet_surface_checks() -> None:
    """Boundary rule on tet: area match, closed-surface flux vanishes."""
    blo = np.array([-0.1, -0.1, -0.1])
    h = np.array([1.1, 1.1, 1.1])
    surf = quadrature.surface_quadrature(TV, TT,
                                         np.arange(len(TT)), blo, h, p=2)
    assert set(surf["face_id"].tolist()) == {0, 1, 2, 3}
    assert check_surface_rule(surf, TET_AREA) < 1e-12
    np.testing.assert_allclose(closed_surface_flux(surf),
                               np.zeros(3), atol=1e-12)


def test_roark_surface_checks() -> None:
    """Real CAD: whole-Roark surface rule matches tess area, flux closes."""
    import os

    from wht_iga.geometry import cad

    step = os.path.join(os.path.dirname(__file__), "..", "..", "data",
                        "roark_plate_hole.step")
    shape = cad.read_step(step)
    tess = cad.tessellate(shape, deflection=0.25)
    lo, hi = cad.shape_bbox(shape)
    pad = np.array([1.0, 1.0, 1.0])
    surf = quadrature.surface_quadrature(tess.nodes, tess.tris,
                                         tess.face_ids, lo - pad,
                                         (hi - lo) + 2 * pad, p=2)
    assert check_surface_rule(surf, cad.tess_area(tess)) < 1e-12
    np.testing.assert_allclose(closed_surface_flux(surf), np.zeros(3),
                               atol=1e-9)


def test_octree_fallback_character() -> None:
    """Fallback: QUAD_FALLBACK logged, nonnegative, converges with depth."""
    from wht_iga.geometry.issues import IssueLog

    V2 = V * np.array([1.0, 2.0, 2.0])  # axis-aligned slab: exact at any depth
    blo = np.array([0.5, 0.2, 0.2])
    h = np.array([0.5, 0.5, 0.5])
    log = IssueLog()
    rule = quadrature.octree_fallback(V2[T], blo, h, 2, 2, 3, log)
    assert rule.fallback
    assert (rule.w >= 0).all() and len(rule.w) > 0
    np.testing.assert_allclose(rule.volume, 0.5**3, rtol=1e-9)
    assert any(i.code == "QUAD_FALLBACK" for i in log)
    # diagonal tet surface: stair-step error shrinks with depth
    tb, th = np.full(3, -0.1), np.full(3, 1.1)
    e3 = abs(quadrature.octree_fallback(
        TV[TT], tb, th, 2, 2, 3, IssueLog()).volume - TET_VOL) / TET_VOL
    e4 = abs(quadrature.octree_fallback(
        TV[TT], tb, th, 2, 2, 4, IssueLog()).volume - TET_VOL) / TET_VOL
    assert e4 < e3 and e4 < 0.10


def test_roark_tangent_ridge_x1_no_fallback() -> None:
    """D4 regression: grazing ridge on x1=10 at y=0 must not suppress the rect edge via seen."""
    import os
    from wht_iga.geometry import cad

    step = os.path.join(os.path.dirname(__file__), "..", "..", "data",
                        "roark_plate_hole.step")
    shape = cad.read_step(step)
    tess = cad.tessellate(shape, deflection=0.25)
    lo = np.array([-100.0, -20.0, -0.5])
    h = np.array([10.0, 10.0, 1.0])
    r_neg = quadrature.cell_quadrature(
        tess.nodes, tess.tris, lo + h * np.array([10, 1, 0]), h, q=4, p=2)
    r_pos = quadrature.cell_quadrature(
        tess.nodes, tess.tris, lo + h * np.array([10, 2, 0]), h, q=4, p=2)
    assert not r_neg.fallback
    assert not r_pos.fallback
    np.testing.assert_allclose(r_neg.volume, r_pos.volume, rtol=1e-12)
