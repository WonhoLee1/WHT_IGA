# -*- coding: utf-8 -*-
"""Trimmed mean-dilatation F-bar with cut-cell pressure groups (fbar_trimmed).

References (all independent of the code under test):
* ``fbar.residual_tangent``: the S40 dense cell kernel (FD-verified there),
  oracle for the untrimmed / cell-wise limit;
* central finite differences of the residual, for the tangent including the
  group coupling U''(Jbar_g)/W_g V_g (x) V_g between cells of one group;
* the closed-form trimmed volume t*(b^2 - pi*a^2/4) of the polygonal hole
  (its polygon radius is chosen so the polygon has the area of the circle);
* an affine stress patch whose hole-wall traction is integrated on the hole
  facets by a separate surface Gauss rule (not by the moment-fitted volume
  rule under test);
* the Kirsch closed form (formulas copied from
  ``tests/kratos_spikes/test_p5_v07_kirsch_precursor.py``; that module imports
  Kratos, so it is not imported here).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from wht_iga.analyses import fbar as FB
from wht_iga.analyses import fbar_trimmed as FT
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

# V07 (tests/verification/data/v07.json)
E_VAL = 100000.0
NU_VAL = 0.49999
T_VAL = 10.0
A_VAL = 1.0
B_VAL = 4.0
T_THICK = 0.25


# ---------------------------------------------------------------------------
# geometry: quarter plate with a polygonal quarter hole
# ---------------------------------------------------------------------------

def plate_hole_soup(a: float, b_ext: float, z_hi: float, n_seg: int = 64):
    """Closed outward-wound soup of {0 <= x, y <= b_ext, r >= a, 0 <= z <= z_hi}.

    Walls x = 0, y = 0, z = 0 lie on the grid's lower planes and the solid
    extends beyond the grid's upper planes (b_ext, z_hi > box hi), as
    ``thermal.quarter_disc_prism_soup`` requires for ``classify`` (ray starts
    below the grid minimum must be outside). The arc polygon radius
    ``a*sqrt(ang/sin(ang))`` gives the quarter hole the exact area pi*a^2/4.
    The caps fan from a point outside the grid (cap edges off grid lines).
    Returns (nodes, tris, arc) with ``arc`` the n_seg+1 hole vertices at
    increasing angle.
    """
    n = int(n_seg)
    ang = 0.5 * np.pi / n
    r_p = a * (ang / np.sin(ang)) ** 0.5
    th = ang * np.arange(n + 1)
    arc = np.column_stack([r_p * np.cos(th), r_p * np.sin(th)])
    poly = np.vstack([arc[:1], [[b_ext, 0.0], [b_ext, b_ext], [0.0, b_ext]], arc[:0:-1]])  # CCW
    m = poly.shape[0]
    fan = (b_ext - 0.0123, b_ext - 0.0137)
    nodes = np.vstack([np.column_stack([poly, np.zeros(m)]),
                       np.column_stack([poly, np.full(m, z_hi)]),
                       [[fan[0], fan[1], 0.0], [fan[0], fan[1], z_hi]]])
    bot_c, top_c = 2 * m, 2 * m + 1
    tris = []
    for i in range(m):
        j = (i + 1) % m
        tris.append((i, j, m + j))          # wall, outward for a CCW polygon
        tris.append((i, m + j, m + i))
        tris.append((top_c, m + i, m + j))  # top cap +z
        tris.append((bot_c, j, i))          # bottom cap -z
    return nodes, np.asarray(tris, dtype=np.int64), arc


def _trimmed_case(a: float, b: float, t: float, n_cells, b_ext: float, z_hi: float,
                  q: int = 3, exact_rules: bool = True):
    """BoxSpline p=2 on [0,b]^2 x [0,t] and its trimmed table (cache off).

    ``exact_rules``: require every cut rule to meet the 1e-12 moment-fit gate
    (no octree fallback, no relaxed fit); the volume / patch / Kirsch tests need it.
    """
    nodes, tris, arc = plate_hole_soup(a, b_ext, z_hi)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(b, b, t), n_cells=n_cells, degree=2)
    qp, out, _iss, _rec, diag = build_trimmed_table(nodes, tris, box, q=q, workers=1,
                                                    cache_dir="off")
    assert diag["n_fallback"] == 0
    if exact_rules:
        assert diag["n_relaxed"] == 0
    return box, qp, out, arc


# ---------------------------------------------------------------------------
# Kirsch closed form (copied from tests/kratos_spikes/test_p5_v07_kirsch_precursor.py)
# ---------------------------------------------------------------------------

def kirsch_displacement(x, y, a=A_VAL, E=E_VAL, nu=NU_VAL, T=T_VAL):
    """Exact plane-strain Kirsch displacement (u_x, u_y)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mu = E / (2.0 * (1.0 + nu))
    kappa = 3.0 - 4.0 * nu
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    q = a / r
    q3 = q ** 3
    fac = T * a / (8.0 * mu)
    ux = fac * ((r / a) * (kappa + 1.0) * np.cos(theta)
                + 2.0 * q * ((kappa + 1.0) * np.cos(theta) + np.cos(3.0 * theta))
                - 2.0 * q3 * np.cos(3.0 * theta))
    uy = fac * ((r / a) * (kappa - 3.0) * np.sin(theta)
                + 2.0 * q * ((1.0 - kappa) * np.sin(theta) + np.sin(3.0 * theta))
                - 2.0 * q3 * np.sin(3.0 * theta))
    return ux, uy


def kirsch_stress(x, y, a=A_VAL, T=T_VAL):
    """Exact Kirsch in-plane Cauchy stresses (sigma_xx, sigma_yy, sigma_xy)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    q2 = (a / r) ** 2
    q4 = q2 * q2
    c2, s2 = np.cos(2.0 * theta), np.sin(2.0 * theta)
    c4, s4 = np.cos(4.0 * theta), np.sin(4.0 * theta)
    sxx = T * (1.0 - q2 * (1.5 * c2 + c4) + 1.5 * q4 * c4)
    syy = T * (-q2 * (0.5 * c2 - c4) - 1.5 * q4 * c4)
    sxy = T * (-q2 * (0.5 * s2 + s4) + 1.5 * q4 * s4)
    return sxx, syy, sxy


def _collocate_1d(knots, p, params, values):
    """Spline coefficients interpolating ``values`` at ``params`` (Greville)."""
    n_cp = len(knots) - p - 1
    A = np.zeros((n_cp, n_cp))
    for m, xi in enumerate(params):
        span = int(find_span(n_cp, p, float(xi), knots))
        A[m, span - p:span + 1] = ders_basis_funs(span, float(xi), p, 0, knots)[0]
    return np.linalg.solve(A, values)


def kirsch_dirichlet(box):
    """Precursor BCs on the Cartesian box: u_z = 0, u_x = 0 on x = 0, u_y = 0 on
    y = 0, exact Kirsch u (Greville interpolation) on x = b and y = b."""
    p = box.degree
    ncx, ncy, ncz = box.n_cp
    kx, ky = box.knots[0], box.knots[1]
    gx, gy = greville(kx, p), greville(ky, p)
    b = float(box.hi[0])
    uxr, uyr = kirsch_displacement(np.full_like(gy, b), gy)
    cxr, cyr = _collocate_1d(ky, p, gy, uxr), _collocate_1d(ky, p, gy, uyr)
    uxt, uyt = kirsch_displacement(gx, np.full_like(gx, float(box.hi[1])))
    cxt, cyt = _collocate_1d(kx, p, gx, uxt), _collocate_1d(kx, p, gx, uyt)
    fixed, prescribed = set(), {}
    for k in range(ncz):
        for j in range(ncy):
            for i in range(ncx):
                I = i + ncx * (j + ncy * k)
                fixed.add(3 * I + 2)
                if i == 0:
                    fixed.add(3 * I)
                if j == 0:
                    fixed.add(3 * I + 1)
                if i == ncx - 1:
                    prescribed[3 * I], prescribed[3 * I + 1] = float(cxr[j]), float(cyr[j])
                elif j == ncy - 1:
                    prescribed[3 * I], prescribed[3 * I + 1] = float(cxt[i]), float(cyt[i])
    prescribed = {d: v for d, v in prescribed.items() if d not in fixed}
    return np.array(sorted(fixed), dtype=np.int64), prescribed


def kirsch_l2_errors(model, u):
    """(e_u, e_sigma_xx) as in the precursor's ``_l2_errors``, on the trimmed rule."""
    X, W = model["qX"], model["qW"]
    uq = FT.displacement_at_qps(model, u)
    sq = FT.cauchy_stress_qps(model, u)
    ux, uy = kirsch_displacement(X[:, 0], X[:, 1])
    sxx, _, _ = kirsch_stress(X[:, 0], X[:, 1])
    du2 = (uq[:, 0] - ux) ** 2 + (uq[:, 1] - uy) ** 2 + uq[:, 2] ** 2
    e_u = float(np.sqrt(np.sum(W * du2) / np.sum(W * (ux ** 2 + uy ** 2))))
    e_s = float(np.sqrt(np.sum(W * (sq[:, 0] - sxx) ** 2) / np.sum(W * sxx ** 2)))
    return e_u, e_s


def _smooth_state(box, scale: float) -> np.ndarray:
    """Deterministic smooth displacement at the control points (no random data)."""
    X = np.asarray(box.control_points())
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    U = np.column_stack([np.sin(1.1 * x + 0.3 * y + 0.2 * z),
                         np.cos(0.7 * x - 0.5 * y + 0.4 * z) - 1.0,
                         0.5 * np.sin(0.9 * y + 0.6 * z + 0.3 * x)])
    return scale * U.ravel()


# ---------------------------------------------------------------------------
# 1. untrimmed limit == fbar (S40 kernel)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_cells", [(3, 2, 2), (2, 3, 1)])
def test_untrimmed_box_matches_fbar_kernel(n_cells):
    """No cut and no OUT cell: every eta = 1, one group per cell -> fbar's R and K."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.5, 1.0, 0.75), n_cells=n_cells, degree=2)
    model = FT.build_model(box, {}, [], E=2.0, nu=0.45)
    ref = FB.build_model(box, 2, E=2.0, nu=0.45)
    assert model["n_groups"] == int(np.prod(n_cells))
    u = _smooth_state(box, 0.04)
    R, K = FT.residual_tangent_sparse(model, u)
    R0, K0 = FB.residual_tangent(ref, u)
    assert np.max(np.abs(R - R0)) <= 1e-12 * np.max(np.abs(R0))
    assert np.max(np.abs(K.toarray() - K0)) <= 1e-12 * np.max(np.abs(K0))


def test_rule_points_on_cell_faces_use_the_cell_span():
    """A cut rule may put points on its cell's faces (moment-fit candidates
    include clipped polygon vertices). Give one interior cell the tensor
    Gauss-Lobatto rule (4 points per axis, endpoints on all six faces, exact to
    degree 5 >= the 2p = 4 of the u = 0 tangent): K(0) and the volume must equal
    the all-Gauss model. Evaluating at the point's own span instead of the
    cell's would assign upper-face points the neighbouring cell's basis."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.5, 1.0, 0.75), n_cells=(3, 2, 2), degree=2)
    h = np.array([0.5, 0.5, 0.375])
    xl = np.array([-1.0, -1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0), 1.0])
    wl = np.array([1.0, 5.0, 5.0, 1.0]) / 6.0
    lo = np.array([0.5, 0.0, 0.0])  # cell (1, 0, 0): upper faces are interior breaks
    rule = np.array([[*(lo + 0.5 * (np.array([a, b, c]) + 1.0) * h),
                      wa * wb * wc * np.prod(h) / 8.0]
                     for c, wc in zip(xl, wl) for b, wb in zip(xl, wl) for a, wa in zip(xl, wl)])
    ref = FT.build_model(box, {}, [], E=1.0, nu=0.3)
    lob = FT.build_model(box, {(1, 0, 0): rule}, [], E=1.0, nu=0.3)
    assert lob["diag"]["n_cut"] == 1
    on_upper = np.isclose(rule[:, :3], lo + h).any(axis=1)
    assert on_upper.sum() == 64 - 27
    assert np.sum(lob["qW"]) == pytest.approx(np.sum(ref["qW"]), rel=1e-14)
    _, K0 = FT.residual_tangent_sparse(ref, np.zeros(ref["n_dof"]))
    _, K1 = FT.residual_tangent_sparse(lob, np.zeros(lob["n_dof"]))
    err = np.max(np.abs((K1 - K0).toarray())) / np.max(np.abs(K0.toarray()))
    assert err < 1e-13, err


# ---------------------------------------------------------------------------
# 2. consistent tangent vs central FD, agglomerated groups included
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("nu", [0.3, 0.49])
def test_tangent_matches_central_fd_with_groups(nu):
    box, qp, out, _ = _trimmed_case(a=0.95, b=1.5, t=0.5, n_cells=(3, 3, 2),
                                    b_ext=1.87, z_hi=0.62, exact_rules=False)
    model = FT.build_model(box, qp, out, E=1.0, nu=nu, eta_thresh=0.25)
    sizes = np.diff(model["grp_cptr"])
    assert np.any(sizes >= 2), "the FD model must contain agglomerated groups"
    u = _smooth_state(box, 0.05)
    u[model["inactive_dofs"]] = 0.0
    R, K = FT.residual_tangent_sparse(model, u)
    assert np.array_equal(R, FT.residual(model, u))
    Kd = K.toarray()
    eps = 1e-6
    Kfd = np.zeros_like(Kd)
    for j in range(model["n_dof"]):
        up, um = u.copy(), u.copy()
        up[j] += eps
        um[j] -= eps
        Kfd[:, j] = (FT.residual(model, up) - FT.residual(model, um)) / (2.0 * eps)
    err = np.max(np.abs(Kd - Kfd)) / np.max(np.abs(Kfd))
    print(f"\n[FD nu={nu}] n_dof={model['n_dof']} groups={model['n_groups']} "
          f"max_group={int(sizes.max())} rel_err={err:.3e}")
    assert err < 1e-5, err
    # entries only the group coupling creates (node pair shares a group, no cell)
    cell_pairs = set()
    for nodes in model["cell_nodes"]:
        cell_pairs.update((int(a), int(b)) for a in nodes for b in nodes)
    extra = []
    for g in range(model["n_groups"]):
        nd = model["grp_nodes"][model["grp_nptr"][g]:model["grp_nptr"][g + 1]]
        extra += [(int(a), int(b)) for a in nd for b in nd if (int(a), int(b)) not in cell_pairs]
    assert extra, "agglomeration must couple node pairs that share no cell"
    rows = np.array([3 * a + c for a, _ in extra for c in range(3) for _d in range(3)])
    cols = np.array([3 * b + d for _, b in extra for _c in range(3) for d in range(3)])
    blk = np.max(np.abs(Kfd[rows, cols]))
    assert blk > 1e-6 * np.max(np.abs(Kfd))  # far above the FD noise (~1e-10)
    err_c = np.max(np.abs(Kd[rows, cols] - Kfd[rows, cols])) / blk
    print(f"[FD nu={nu}] group-only node pairs={len(extra)} max|K|={blk:.3e} rel_err={err_c:.3e}")
    assert err_c < 1e-5


# ---------------------------------------------------------------------------
# 3. group construction invariants
# ---------------------------------------------------------------------------

def _face_connected(cells: np.ndarray) -> bool:
    todo = {tuple(c) for c in cells}
    stack = [todo.pop()]
    while stack:
        i, j, k = stack.pop()
        for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            nb = (i + d[0], j + d[1], k + d[2])
            if nb in todo:
                todo.remove(nb)
                stack.append(nb)
    return not todo


@pytest.mark.parametrize("h", [0.5, 0.25])
def test_pressure_group_invariants_v07(h):
    n = int(round(B_VAL / h))
    box, qp, out, _ = _trimmed_case(A_VAL, B_VAL, T_THICK, (n, n, 2), 4.37, 0.37)
    thr = 0.25
    model = FT.build_model(box, qp, out, E_VAL, NU_VAL, eta_thresh=thr)
    nx, ny, nz = box.n_cells
    out_set = {tuple(c) for c in out}
    active = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)
              if (i, j, k) not in out_set]
    assert sorted(map(tuple, model["cell_ijk"].tolist())) == sorted(active)
    # every active cell in exactly one group; OUT cells in none
    gg = model["group_grid"]
    assert all(gg[c] >= 0 for c in active) and all(gg[c] == -1 for c in out_set)
    members = np.concatenate([model["grp_cells"][model["grp_cptr"][g]:model["grp_cptr"][g + 1]]
                              for g in range(model["n_groups"])])
    assert sorted(members.tolist()) == list(range(len(active)))
    cv = model["cell_vol"]
    eta = model["cell_eta"]
    assert np.sum(eta < thr) >= 2, "V07 grids have small cut cells"
    for g in range(model["n_groups"]):
        cs = model["grp_cells"][model["grp_cptr"][g]:model["grp_cptr"][g + 1]]
        vol = float(np.sum(model["qW"][np.isin(model["q_cell"], cs)]))
        assert vol / cv >= thr
        assert vol == pytest.approx(model["group_volume"][g], rel=1e-14)
        assert np.sum(eta[cs] >= thr) == 1  # one large root per group
        assert _face_connected(model["cell_ijk"][cs])
    assert np.all(model["group_volume"] / cv >= thr)
    # cell-wise P0: one group per active cell
    m0 = FT.build_model(box, qp, out, E_VAL, NU_VAL, eta_thresh=0.0)
    assert m0["n_groups"] == len(active) and np.all(np.diff(m0["grp_cptr"]) == 1)
    # inactive control points: exactly the mask of solver.active; no stiffness there
    assert np.array_equal(model["inactive_ctrl"], inactive_ctrl_mask(box, out))
    assert model["inactive_ctrl"].sum() > 0
    _, K = FT.residual_tangent_sparse(model, np.zeros(model["n_dof"]))
    rows = np.diff(K.indptr)
    assert np.all(rows[model["inactive_dofs"]] == 0)
    with pytest.raises(ValueError, match="inactive"):
        FT.solve_static(model, np.zeros(model["n_dof"]), [],
                        prescribed={int(model["inactive_dofs"][0]): 1e-3})
    # kernel group Jbar == independent numpy gather by group
    u = _smooth_state(box, 0.01)
    F = np.eye(3) + np.einsum("qia,qim->qam",
                              u.reshape(-1, 3)[model["cell_nodes"]][model["q_cell"]], model["qG"])
    J = np.linalg.det(F)
    gq = model["cell_group"][model["q_cell"]]
    jb = (np.bincount(gq, weights=model["qW"] * J) / np.bincount(gq, weights=model["qW"]))
    assert np.max(np.abs(FT.group_jbar(model, u) - jb)) < 1e-13


def test_pressure_groups_bfs_chain_and_unreachable():
    """Small cells reach a root through other small cells; isolated ones raise."""
    eta = np.zeros((5, 1, 1))
    eta[:, 0, 0] = [0.05, 0.1, 0.02, 0.9, 0.0]
    g = FT.pressure_groups(eta, 0.25)
    assert g[:, 0, 0].tolist() == [0, 0, 0, 0, -1]
    eta2 = np.zeros((4, 4, 1))
    eta2[0, 0, 0] = 0.1   # small, no active neighbour at all
    eta2[3, 3, 0] = 1.0
    with pytest.raises(ValueError, match="without a path"):
        FT.pressure_groups(eta2, 0.25)
    eta3 = np.full((3, 1, 1), 0.1)
    with pytest.raises(ValueError, match="eta_thresh"):
        FT.pressure_groups(eta3, 0.25)
    # tie on eta -> smallest flat id; edge neighbour only when no face neighbour
    eta4 = np.zeros((3, 3, 1))
    eta4[1, 1, 0] = 0.1
    eta4[0, 1, 0] = eta4[2, 1, 0] = 0.8
    eta4[2, 2, 0] = 0.05      # touches the root (2,1) by a face
    eta4[0, 0, 0] = 0.6       # large, edge neighbour of (1,1) only
    g4 = FT.pressure_groups(eta4, 0.25)
    root_of = {tuple(map(int, c)): int(g4[tuple(c)]) for c in np.argwhere(eta4 >= 0.25)}
    assert g4[1, 1, 0] == root_of[(0, 1, 0)]
    assert g4[2, 2, 0] == root_of[(2, 1, 0)]
    eta6 = np.zeros((3, 1, 1))
    eta6[:, 0, 0] = [0.3, 0.1, 0.9]  # attach to the larger of two large neighbours
    assert FT.pressure_groups(eta6, 0.25)[:, 0, 0].tolist() == [0, 1, 1]
    eta5 = np.zeros((2, 2, 1))
    eta5[0, 0, 0], eta5[1, 1, 0] = 1.0, 0.1   # vertex contact only
    assert FT.pressure_groups(eta5, 0.25)[1, 1, 0] == 0


# ---------------------------------------------------------------------------
# 4. trimmed volume + affine stress patch with surface-integrated hole traction
# ---------------------------------------------------------------------------

def _hole_wall_loads(box, arc, P, z_lo, z_hi):
    """f_I = int_{hole wall} N_I (P N) dA on the planar hole facets.

    Each facet (arc chord x [z_lo, z_hi]) is split at the grid lines, where
    the spline is a polynomial of degree p per axis; Gauss rules with p+2 (in
    plane) and p+1 (in z) points are exact there. N points into the hole
    (outward for the domain).
    """
    p = int(box.degree)
    knots = [np.asarray(k, dtype=float) for k in box.knots]
    ncx, ncy, ncz = box.n_cp
    f = np.zeros(3 * ncx * ncy * ncz)
    br = [np.unique(k) for k in knots]
    xs, ws = np.polynomial.legendre.leggauss(p + 2)
    zs, wz = np.polynomial.legendre.leggauss(p + 1)
    zb = br[2][(br[2] >= z_lo) & (br[2] <= z_hi)]
    for a0, a1 in zip(arc[:-1], arc[1:]):
        d = a1 - a0
        L = float(np.hypot(*d))
        N = np.array([-d[1], d[0], 0.0]) / L
        tr = P @ N
        cuts = [0.0, 1.0]
        for ax in (0, 1):
            if d[ax] != 0.0:
                s = (br[ax] - a0[ax]) / d[ax]
                cuts += [float(v) for v in s if 0.0 < v < 1.0]
        cuts = np.unique(cuts)
        for s0, s1 in zip(cuts[:-1], cuts[1:]):
            for z0, z1 in zip(zb[:-1], zb[1:]):
                xm = a0 + 0.5 * (s0 + s1) * d
                zm = 0.5 * (z0 + z1)
                sp = [int(find_span(n, p, float(c), k)) for n, c, k in
                      zip((ncx, ncy, ncz), (xm[0], xm[1], zm), knots)]
                for si, wi in zip(xs, ws):
                    s = s0 + 0.5 * (si + 1.0) * (s1 - s0)
                    x, y = a0 + s * d
                    Nx = ders_basis_funs(sp[0], x, p, 0, knots[0])[0]
                    Ny = ders_basis_funs(sp[1], y, p, 0, knots[1])[0]
                    for zi, wzi in zip(zs, wz):
                        z = z0 + 0.5 * (zi + 1.0) * (z1 - z0)
                        Nz = ders_basis_funs(sp[2], z, p, 0, knots[2])[0]
                        w = wi * 0.5 * (s1 - s0) * L * wzi * 0.5 * (z1 - z0)
                        for c in range(p + 1):
                            for b in range(p + 1):
                                for a in range(p + 1):
                                    I = (sp[0] - p + a) + ncx * ((sp[1] - p + b)
                                                                + ncy * (sp[2] - p + c))
                                    f[3 * I:3 * I + 3] += w * Nx[a] * Ny[b] * Nz[c] * tr
    return f


@pytest.mark.parametrize("eta_thresh", [0.25, 0.0])
def test_trimmed_volume_and_affine_patch(eta_thresh):
    """Multi-cell trimmed plate: exact volume, then an affine finite-strain field
    (all box faces Dirichlet, hole loaded by P N) is reproduced to 1e-11."""
    a, b, t = 1.2, 2.0, 0.5
    box, qp, out, arc = _trimmed_case(a, b, t, (4, 4, 2), b_ext=2.41, z_hi=0.66)
    E, nu = 1.0, 0.3
    model = FT.build_model(box, qp, out, E, nu, eta_thresh=eta_thresh)
    v_exact = t * (b * b - 0.25 * np.pi * a * a)
    assert abs(float(np.sum(model["qW"])) - v_exact) <= 1e-12 * v_exact
    assert model["diag"]["n_cut"] > 0
    if eta_thresh > 0:
        assert model["diag"]["max_group_cells"] >= 2

    Agrad = np.array([[0.03, -0.012, 0.004], [0.008, -0.02, 0.006], [0.005, 0.01, 0.015]])
    c0 = np.array([0.001, -0.002, 0.0005])
    F = np.eye(3) + Agrad
    J = float(np.linalg.det(F))
    mu, lam = model["mu"], model["lam"]
    P = mu * F + (lam * np.log(J) - mu) * np.linalg.inv(F).T     # first Piola of W
    sig = mu / J * (F @ F.T - np.eye(3)) + lam * np.log(J) / J * np.eye(3)
    X = np.asarray(box.control_points())
    u_aff = (X @ Agrad.T + c0).ravel()

    ncx, ncy, ncz = box.n_cp
    act = ~model["inactive_ctrl"]
    presc = {}
    for kk in range(ncz):
        for j in range(ncy):
            for i in range(ncx):
                I = i + ncx * (j + ncy * kk)
                on_face = i in (0, ncx - 1) or j in (0, ncy - 1) or kk in (0, ncz - 1)
                if on_face and act[I]:
                    for d in range(3):
                        presc[3 * I + d] = float(u_aff[3 * I + d])
    f = _hole_wall_loads(box, arc, P, 0.0, t)
    hist = []
    u = FT.solve_static(model, f, [], n_steps=1, prescribed=presc, rtol=1e-13, atol=0.0,
                        dtol=0.0, history=hist)
    free_act = np.ones(model["n_dof"], dtype=bool)
    free_act[model["inactive_dofs"]] = False
    free_act[list(presc)] = False
    assert free_act.sum() >= 30
    err_u = np.max(np.abs(u[free_act] - u_aff[free_act])) / np.max(np.abs(u_aff))
    assert err_u < 1e-11, err_u
    S = FT.cauchy_stress_qps(model, u)
    ref = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]])
    err_s = np.max(np.abs(S - ref[None, :])) / np.max(np.abs(ref))
    print(f"\n[patch eta_thresh={eta_thresh}] vol_err={abs(np.sum(model['qW']) - v_exact) / v_exact:.1e} "
          f"groups={model['n_groups']}/{model['diag']['n_active']} free={int(free_act.sum())} "
          f"newton={hist} err_u={err_u:.2e} err_sigma={err_s:.2e}")
    assert err_s < 1e-11, err_s
    assert np.all(u[model["inactive_dofs"]] == 0.0)


# ---------------------------------------------------------------------------
# 5. V07 Kirsch on the Cartesian grid: agglomerated vs cell-wise P0
# ---------------------------------------------------------------------------

def _scaled_condition(model, fixed, prescribed) -> float:
    """cond(D K D) of the u = 0 tangent on the free DOFs, D = diag(K)^-1/2 (dense eig)."""
    free = np.ones(model["n_dof"], dtype=bool)
    free[fixed] = False
    free[list(prescribed)] = False
    free[model["inactive_dofs"]] = False
    _, K = FT.residual_tangent_sparse(model, np.zeros(model["n_dof"]))
    A = K[np.flatnonzero(free)][:, np.flatnonzero(free)].toarray()
    d = 1.0 / np.sqrt(np.diag(A))
    A = d[:, None] * A * d[None, :]
    ev = np.linalg.eigvalsh(0.5 * (A + A.T))
    assert ev[0] > 0.0
    return float(ev[-1] / ev[0])


@pytest.mark.slow
def test_v07_kirsch_cartesian_agglomeration():
    """V07 (nu = 0.49999) on BoxSpline (h, h, 0.125) grids with the hole trimmed.

    Gate (task): e_u < 0.02 and e_sigma_xx < 0.25 per level, errors decrease
    h0 -> h1 with slopes > 0.5 (v07.json source-doc criterion). Both pressure
    spaces are solved: agglomerated groups (eta_thresh 0.25) and cell-wise P0
    (eta_thresh 0).

    Effect of agglomeration (measured 2026-11-19, docs in fbar_trimmed): no
    global locking of cell-wise P0 on these grids (smallest cut fraction 0.087
    / 0.042, both variants pass the gate); the measurable difference is the
    stability of the discrete system: the Jacobi-scaled tangent condition
    number is lower with agglomeration (asserted as an ordering, no tuned value).
    """
    v_exact = T_THICK * (B_VAL ** 2 - 0.25 * np.pi * A_VAL ** 2)
    rows = []
    for h in (0.5, 0.25):
        n = int(round(B_VAL / h))
        t0 = time.perf_counter()
        box, qp, out, _ = _trimmed_case(A_VAL, B_VAL, T_THICK, (n, n, 2), 4.37, 0.37)
        t_tab = time.perf_counter() - t0
        fixed, prescribed = kirsch_dirichlet(box)
        sol = {}
        for thr in (0.25, 0.0):
            t1 = time.perf_counter()
            model = FT.build_model(box, qp, out, E_VAL, NU_VAL, eta_thresh=thr)
            vol = float(np.sum(model["qW"]))
            assert abs(vol - v_exact) <= 1e-3 * v_exact
            hist = []
            u = FT.solve_static(model, np.zeros(model["n_dof"]), fixed, n_steps=1,
                                prescribed=prescribed, history=hist)
            e_u, e_s = kirsch_l2_errors(model, u)
            dt = time.perf_counter() - t1
            cond = _scaled_condition(model, fixed, prescribed)
            sol[thr] = (u, model, cond)
            d = model["diag"]
            rows.append((h, thr, e_u, e_s))
            print(f"\n[V07 trimmed h={h} eta_thresh={thr}] n_dof={model['n_dof']} "
                  f"active={d['n_active']} cut={d['n_cut']} groups={d['n_groups']} "
                  f"max_group={d['max_group_cells']} min_eta={d['min_eta']:.4f} "
                  f"min_supp={d['min_basis_support_ratio']:.2e} "
                  f"vol_err={abs(vol - v_exact) / v_exact:.1e} newton={hist} "
                  f"e_u={e_u:.4e} e_sigma_xx={e_s:.4e} cond(DKD)={cond:.3e} "
                  f"t_table={t_tab:.2f}s t_build+solve={dt:.2f}s")
            assert e_u < 0.02, (h, thr, e_u)
            assert e_s < 0.25, (h, thr, e_s)
        ua, ma, ca = sol[0.25]
        uc, mc, cc = sol[0.0]
        assert ma["diag"]["max_group_cells"] >= 2 and ma["n_groups"] < ma["diag"]["n_active"]
        assert mc["n_groups"] == mc["diag"]["n_active"]
        X = ma["qX"]
        u_ref = np.max(np.hypot(*kirsch_displacement(X[:, 0], X[:, 1])))
        d_field = np.max(np.abs(FT.displacement_at_qps(ma, ua) - FT.displacement_at_qps(mc, uc)))
        d_coef = np.max(np.abs(ua - uc))
        print(f"[V07 trimmed h={h}] u_h(agg) - u_h(cell-wise) at QPs: {d_field / u_ref:.3e} max|u|; "
              f"coefficients: {d_coef / u_ref:.3e} max|u|; cond(DKD) agg/cell = {ca / cc:.3f}")
        assert ca < cc
    for thr in (0.25, 0.0):
        r0, r1 = [r for r in rows if r[1] == thr]
        slope_u = np.log(r0[2] / r1[2]) / np.log(2.0)
        slope_s = np.log(r0[3] / r1[3]) / np.log(2.0)
        print(f"[V07 trimmed eta_thresh={thr}] slope_u={slope_u:.3f} slope_sigma_xx={slope_s:.3f}")
        assert slope_u > 0.5 and slope_s > 0.5
