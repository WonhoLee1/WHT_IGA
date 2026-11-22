# -*- coding: utf-8 -*-
"""Generalised mean dilatation with cell-wise P1 / P2-disc pressure (mixed_p1, S49).

References (all independent of the code under test):
* ``fbar.residual_tangent`` / ``fbar.cauchy_stress`` (S40 cell kernel,
  FD-verified in ``test_fbar_mean_dilatation.py``): oracle for the P0 limit;
* central finite differences of the residual, for the consistent tangent;
* closed-form polynomials, for the cell L2 projection;
* exact solutions: an affine finite-strain field with the closed-form
  Neo-Hookean stress, and plane-strain pure bending (quadratic u, linear
  stress and pressure) on an affine mesh;
* V07 Kirsch closed form, mapped quarter-plate geometry, BCs and L2 error
  helpers of ``tests/kratos_spikes/test_p5_v07_kirsch_precursor.py``. That
  module is loaded by path: at import it only binds ``load_kratos`` (no
  Kratos load), so the precursor's own functions are the reference here;
  the gate constants are parsed from ``tests/verification/data/v07.json``;
* for the stress recovery (``cauchy_stress_linearized``, ``cauchy_stress``):
  polynomial displacement fields with analytic gradients on an affine
  (sheared) box, a tensor Gauss rule built here from ``numpy.leggauss`` and a
  cell-wise weighted least-squares L2 projection written here (``lstsq``,
  not the code's Cholesky-orthonormal basis);
* the best-approximation property of the L2 projection, as a floor for the
  P0 mean-stress error on the V07 plate.

V07 is linear elastic (v07.json: "small-strain; NOT finite strain"), so the
Kirsch runs use ``solve_linearized`` (the consistent tangent at u = 0). The
precursor's Neo-Hookean Newton differs from it by O(|grad u|) ~ 1e-4 in u,
which is above the p=3 errors measured here (S49).
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses import fbar
from wht_iga.analyses import mixed_p1 as MP
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import n_control_points, uniform_clamped
from wht_iga.reference.box import BoxSpline

_PRE_PATH = Path(__file__).resolve().parents[1] / "kratos_spikes" / "test_p5_v07_kirsch_precursor.py"
_V07_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v07.json"


def _load_precursor():
    name = "_s49_v07_precursor"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _PRE_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # @dataclass resolves its module through sys.modules
    spec.loader.exec_module(mod)
    return mod


PRE = _load_precursor()
V_EXACT = PRE.T_THICK * (PRE.B_VAL ** 2 - 0.25 * np.pi * PRE.A_VAL ** 2)
# v07.json acceptance: per-level (e_u, e_sigma_xx) bounds and the plan's
# target slopes (u 2.3, sigma_xx 2.2); domain volume within 0.1 %.
GATE = {0.5: (0.01, 0.2), 0.25: (0.005, 0.1), 0.125: (0.002, 0.05)}
VOL_TOL = 1e-3
SLOPE_U, SLOPE_S = 2.3, 2.2
LEVEL = {0.5: (6, 8, 2), 1.0 / 3.0: (9, 12, 2), 0.25: (12, 16, 2), 0.125: (24, 32, 2)}


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

class _Plate:
    """``PRE.MappedQuarterPlate`` with degree p and in-plane continuity C^k.

    Copied and adapted from ``MappedQuarterPlate`` in
    ``tests/kratos_spikes/test_p5_v07_kirsch_precursor.py`` (2026-11-20): the
    same collocated inner circle / outer square and linear radial blend, the
    v = 0.5 corner knot of multiplicity p (C^0), interior knots of
    multiplicity p - k, thickness knots C^(p-1). p=2, k=1 is the precursor's
    plate (checked bit for bit in ``test_generalised_plate_is_precursor_plate``).
    """

    def __init__(self, n_u, n_v_half, n_w, degree=2, cont=1, a=PRE.A_VAL, b=PRE.B_VAL, t=PRE.T_THICK):
        p, k = int(degree), int(cont)
        self.degree, self.a, self.b, self.t = p, a, b, t
        ku = uniform_clamped(0.0, 1.0, n_u, p, continuity=k)
        br_l = np.repeat(np.linspace(0.0, 0.5, n_v_half + 1)[1:-1], p - k)
        br_r = np.repeat(np.linspace(0.5, 1.0, n_v_half + 1)[1:-1], p - k)
        kv = np.concatenate([np.full(p + 1, 0.0), br_l, np.full(p, 0.5), br_r, np.full(p + 1, 1.0)])
        kw = uniform_clamped(0.0, t, n_w, p, continuity=p - 1)
        self.knots = [ku, kv, kw]
        self.n_cp = tuple(n_control_points(kk, p) for kk in self.knots)
        gu, gv, gw = (greville(kk, p) for kk in self.knots)
        cp_in = PRE._collocate_1d(kv, p, gv, PRE._inner_circle_xy(gv, a))
        cp_out = PRE._collocate_1d(kv, p, gv, PRE._outer_square_xy(gv, b))
        nu, nv, nw = self.n_cp
        cp = np.zeros((nu * nv * nw, 3))
        for kz, z in enumerate(gw):
            for j in range(nv):
                for i, uu in enumerate(gu):
                    xy = (1.0 - uu) * cp_in[j] + uu * cp_out[j]
                    cp[i + nu * (j + nv * kz)] = [xy[0], xy[1], z]
        self._cp = cp

    def control_points(self):
        return self._cp


S_AFFINE = np.array([[1.0, 0.3, 0.1], [0.0, 1.2, -0.2], [0.05, 0.0, 0.9]])


class _ShearedBox(BoxSpline):
    """BoxSpline under the affine map X = S_AFFINE Y: parallelepiped cells."""

    def control_points(self):
        return super().control_points() @ S_AFFINE.T


def _greville_interpolant(box, f):
    """Control values (n_ctrl, 3) interpolating f(X) at the tensor Greville
    grid; exact for every f in the spline space (unisolvent collocation)."""
    p = box.degree
    g = [greville(kk, p) for kk in box.knots]
    Y = np.stack(np.meshgrid(*g, indexing="ij"), axis=-1)
    F = f(Y @ S_AFFINE.T)
    for ax in range(3):
        n = len(g[ax])
        A = np.zeros((n, n))
        for m, xi in enumerate(g[ax]):
            span = int(find_span(n, p, float(xi), box.knots[ax]))
            A[m, span - p:span + 1] = ders_basis_funs(span, float(xi), p, 0, box.knots[ax])[0]
        Fm = np.moveaxis(F, ax, 0)
        F = np.moveaxis(np.linalg.solve(A, Fm.reshape(n, -1)).reshape(Fm.shape), 0, ax)
    return np.transpose(F, (2, 1, 0, 3)).reshape(-1, 3)  # I = i + n0 (j + n1 k)


def _surface_dofs(n_cp):
    n0, n1, n2 = n_cp
    k, j, i = np.meshgrid(np.arange(n2), np.arange(n1), np.arange(n0), indexing="ij")
    on = ((i == 0) | (i == n0 - 1) | (j == 0) | (j == n1 - 1) | (k == 0) | (k == n2 - 1)).ravel()
    return np.flatnonzero(np.repeat(on, 3))


def _smooth_state(cp, scale):
    """Deterministic smooth displacement at the control points (no random data)."""
    x, y, z = cp[:, 0], cp[:, 1], cp[:, 2]
    U = np.column_stack([np.sin(1.1 * x + 0.3 * y + 0.2 * z),
                         np.cos(0.7 * x - 0.5 * y + 0.4 * z) - 1.0,
                         0.5 * np.sin(0.9 * y + 0.6 * z + 0.3 * x)])
    return scale * U.ravel()


def _rel(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(b)))


# ---------------------------------------------------------------------------
# 1. geometry and pressure space
# ---------------------------------------------------------------------------

def test_gate_constants_are_v07_json():
    """GATE / VOL_TOL / SLOPE_* used by the V07 tests are the v07.json hard
    gates (parsed from its acceptance text), at its in-plane h sweep."""
    v07 = json.loads(_V07_JSON.read_text(encoding="utf-8"))
    acc = v07["acceptance"]
    num = r"([0-9]+(?:\.[0-9]+)?)"
    pairs = [(float(a), float(b)) for a, b in re.findall(rf"\({num}, {num}\)", acc["gate"])]
    hs = [float(h) for h in v07["mesh_sweep"]["in_plane_h"]]
    assert hs == sorted(GATE, reverse=True) and pairs == [GATE[h] for h in hs]
    slopes = re.search(rf"slopes >= {num} / {num}", acc["gate"])
    assert (float(slopes.group(1)), float(slopes.group(2))) == (SLOPE_U, SLOPE_S)
    assert float(acc["domain_volume_agreement_tol"]) == VOL_TOL


def test_generalised_plate_is_precursor_plate():
    """p=2, C^1 reproduces the precursor's MappedQuarterPlate bit for bit."""
    mine, ref = _Plate(6, 8, 2, degree=2, cont=1), PRE.MappedQuarterPlate(n_u=6, n_v_half=8, n_w=2, degree=2)
    assert mine.n_cp == ref.n_cp
    for k1, k2 in zip(mine.knots, ref.knots):
        assert np.array_equal(k1, k2)
    assert np.array_equal(mine.control_points(), ref.control_points())


def test_pressure_basis_orthonormal_and_reproduces_polynomials():
    """Cell basis is W-orthonormal; Pi_Pk reproduces degree-k polynomials in X
    on the curved plate cells (physical coords); P0 is the fbar cell mean."""
    plate = _Plate(2, 2, 1, degree=2, cont=0)
    lin = lambda X: 1.0 + 0.3 * X[..., 0] - 0.2 * X[..., 1] + 0.5 * X[..., 2]  # noqa: E731
    quad = lambda X: lin(X) + 0.07 * X[..., 0] ** 2 - 0.11 * X[..., 0] * X[..., 1] + 0.4 * X[..., 2] ** 2  # noqa: E731
    out = {}
    for space, exps_n in (("P0", 1), ("P1", 4), ("P2", 10)):
        m = MP.build_model(plate, 2, E=1.0, nu=0.3, pressure=space)
        Psi, W, X = m["Psi_all"], m["W_all"], m["X_all"]
        assert Psi.shape[2] == exps_n == m["n_pressure"]
        gram = np.einsum("cqk,cq,cql->ckl", Psi, W, Psi)
        assert np.max(np.abs(gram - np.eye(exps_n))) < 1e-12
        out[space] = (m, X, W)
    # exact reproduction (round-off level; the projection is an orthonormal sum)
    m1, X, W = out["P1"]
    assert _rel(MP.pressure_projection(m1, lin(X)), lin(X)) < 1e-12
    m2 = out["P2"][0]
    assert _rel(MP.pressure_projection(m2, quad(X)), quad(X)) < 1e-12
    # P1 cannot reproduce the quadratic: the check above has power
    assert _rel(MP.pressure_projection(m1, quad(X)), quad(X)) > 1e-3
    # P0 = cell mean with the reference weights (fbar's Jbar definition)
    m0 = out["P0"][0]
    f = quad(X)
    mean = (np.sum(W * f, axis=1) / np.sum(W, axis=1))[:, None] * np.ones_like(f)
    assert _rel(MP.pressure_projection(m0, f), mean) < 1e-14
    # parametric (mapped) P1 is a different space on these curved cells
    mpar = MP.build_model(plate, 2, E=1.0, nu=0.3, pressure="P1", coords="parametric")
    assert _rel(MP.pressure_projection(mpar, lin(X)), lin(X)) > 1e-6


# ---------------------------------------------------------------------------
# 2. P0 limit and consistent tangent
# ---------------------------------------------------------------------------

def test_p0_limit_reproduces_fbar():
    """Projection onto constants is fbar's mean dilatation: R, K, sigma agree
    with fbar.residual_tangent / cauchy_stress to 1e-12 (curved cells,
    nu = 0.49999, J varying by ~1 %)."""
    plate = PRE.MappedQuarterPlate(n_u=2, n_v_half=2, n_w=1, degree=2)
    model = MP.build_model(plate, 2, E=PRE.E_VAL, nu=PRE.NU_VAL, pressure="P0")
    u = _smooth_state(plate.control_points(), 1e-2)
    Rf, Kf = fbar.residual_tangent(model, u)
    Rm, Km = MP.residual_tangent(model, u)
    assert _rel(Rm, Rf) < 1e-12
    assert _rel(Km, Kf) < 1e-12
    Sf, Sm = np.stack(fbar.cauchy_stress(model, u)), np.stack(MP.cauchy_stress(model, u))
    assert _rel(Sm, Sf) < 1e-12
    # CSR path == dense oracle (same additions, same order)
    Rs, Ks = MP.residual_tangent_sparse(model, u)
    assert _rel(Rs, Rm) < 1e-14 and _rel(Ks.toarray(), Km) < 1e-14


@pytest.mark.parametrize("p,cont,space,coords,nu", [
    (2, 0, "P1", "physical", 0.45),
    (2, 0, "P1", "physical", 0.49999),
    (2, 0, "P1", "parametric", 0.45),
    (3, 0, "P2", "physical", 0.49999),
])
def test_consistent_tangent_central_fd(p, cont, space, coords, nu):
    """K = dR/du by central differences (rel < 1e-6), K symmetric."""
    plate = _Plate(1, 1, 1, degree=p, cont=cont)
    model = MP.build_model(plate, p, E=3.0, nu=nu, pressure=space, coords=coords)
    u = _smooth_state(plate.control_points(), 1e-2)
    R, K = MP.residual_tangent(model, u)
    Rs, Ks = MP.residual_tangent_sparse(model, u)
    assert _rel(Rs, R) < 1e-14 and _rel(Ks.toarray(), K) < 1e-14
    h = 1e-7  # |u| ~ 1e-2: truncation ~ h^2, round-off ~ eps |R| / h << 1e-6 |K|
    Kfd = np.empty_like(K)
    for j in range(u.size):
        up, um = u.copy(), u.copy()
        up[j] += h
        um[j] -= h
        Kfd[:, j] = (MP.residual_tangent_sparse(model, up)[0] - MP.residual_tangent_sparse(model, um)[0]) / (2 * h)
    err = _rel(K, Kfd)
    print(f"\n[FD p={p} C^{cont} {space}/{coords} nu={nu}] ndof={u.size} rel={err:.2e}")
    assert err < 1e-6
    assert _rel(K, K.T) < 1e-12


# ---------------------------------------------------------------------------
# 3. patch tests with exact solutions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("p,space", [(2, "P1"), (3, "P2")])
def test_affine_finite_strain_patch(p, space):
    """u = A X on every surface control point of the curved plate, nu = 0.49999,
    J = 1.02: Newton returns A X, sigma = mu/J B + U'(J) I (closed form);
    the linearised solve returns A X with sigma = 2 mu eps + lam tr(eps) I."""
    E, nu = PRE.E_VAL, PRE.NU_VAL
    mu, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    A = np.array([[0.02, 0.01, 0.0], [-0.015, 0.01, 0.005], [0.0, 0.004, -0.01]])
    plate = _Plate(2, 2, 1, degree=p, cont=0)
    model = MP.build_model(plate, p, E=E, nu=nu, pressure=space)
    U = (plate.control_points() @ A.T).ravel()
    presc = {int(d): float(U[d]) for d in _surface_dofs(plate.n_cp)}
    assert 3 * len(plate.control_points()) - len(presc) > 0  # interior DOFs are free
    hist = []
    u = MP.solve_static(model, np.zeros(model["n_dof"]), np.array([], dtype=np.int64),
                        prescribed=presc, history=hist)
    F = np.eye(3) + A
    J = np.linalg.det(F)
    sig = mu / J * F @ F.T + (lam * np.log(J) - mu) / J * np.eye(3)
    sv = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]])
    # round-off: K/mu ~ 5e4 times the mesh condition number amplifies eps
    assert _rel(u, U) < 1e-9
    assert _rel(np.stack(MP.cauchy_stress(model, u)), np.broadcast_to(sv, (model["conn"].shape[0], model["W_all"].shape[1], 6))) < 1e-10
    # one increment: the consistent predictor starts Newton in its quadratic basin
    assert len(hist) == 1
    ul = MP.solve_linearized(model, np.zeros(model["n_dof"]), np.array([], dtype=np.int64), presc)
    eps = 0.5 * (A + A.T)
    sl = 2 * mu * eps + lam * np.trace(eps) * np.eye(3)
    slv = np.array([sl[0, 0], sl[1, 1], sl[2, 2], sl[0, 1], sl[1, 2], sl[0, 2]])
    assert _rel(ul, U) < 1e-8
    assert _rel(MP.cauchy_stress_linearized(model, ul), np.broadcast_to(slv, (model["conn"].shape[0], model["W_all"].shape[1], 6))) < 1e-8


def test_linear_pressure_bending_patch_p1_exact_p0_not():
    """Plane-strain pure bending sigma_xx = c y (sigma_zz = nu c y, u quadratic,
    pressure linear), nu = 0.49999, affine (sheared) p=2 C^0 mesh, Dirichlet
    on the surface: P1-disc reproduces u and sigma (Pi div u = div u, so the
    exact field solves the discrete equations); P0 cannot represent the
    linear volumetric stress."""
    E, nu, c = PRE.E_VAL, PRE.NU_VAL, 10.0
    kk, mm = c * (1 - nu ** 2) / E, nu * (1 + nu) * c / E

    def u_ex(X):
        x, y = X[..., 0], X[..., 1]
        return np.stack([kk * x * y, -0.5 * kk * x * x - 0.5 * mm * y * y, 0.0 * x], axis=-1)

    def s_ex(X):
        S = np.zeros(X.shape[:-1] + (6,))
        S[..., 0], S[..., 2] = c * X[..., 1], nu * c * X[..., 1]
        return S

    res = {}
    for p, space in ((2, "P1"), (2, "P0"), (3, "P2")):
        box = _ShearedBox(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(3, 2, 2), degree=p, continuity=0)
        model = MP.build_model(box, p, E=E, nu=nu, pressure=space)
        U = _greville_interpolant(box, u_ex).ravel()
        presc = {int(d): float(U[d]) for d in _surface_dofs(box.n_cp)}
        u = MP.solve_linearized(model, np.zeros(model["n_dof"]), np.array([], dtype=np.int64), presc)
        e_s = float(np.max(np.abs(MP.cauchy_stress_linearized(model, u) - s_ex(model["X_all"])))) / c
        res[(p, space)] = (_rel(u, U), e_s)
        print(f"\n[bending p={p} {space}] max rel u err {res[(p, space)][0]:.2e}, max sigma err / c {e_s:.2e}")
    # exact up to round-off (lam / mu ~ 5e4 amplifies eps in the stress)
    assert res[(2, "P1")][0] < 1e-9 and res[(2, "P1")][1] < 1e-8
    assert res[(3, "P2")][0] < 1e-9 and res[(3, "P2")][1] < 1e-8
    # P0: the cell-constant volumetric stress misses the linear part by
    # O(c * cell height) (cells 0.5 high in y)
    assert res[(2, "P0")][1] > 0.05


# ---------------------------------------------------------------------------
# 3b. stress recovery against an independent projection (polynomial fields)
# ---------------------------------------------------------------------------

S_INV = np.linalg.inv(S_AFFINE)
_DEG = {"P0": 0, "P1": 1, "P2": 2}
# u(X) = g(S^-1 X), g in Q_p(Y) (degree <= p in each Y_b): it lies in every
# degree-p spline space of the box, so the Greville interpolant is exact and
# grad u is analytic. Terms (coefficient, (i, j, k)) per component; the
# mixed / high-order terms make div u leave P_k(X) for every k <= 2 used.
_FIELD = {
    2: [[(0.06, (2, 0, 0)), (-0.05, (1, 1, 0)), (0.04, (0, 0, 2)), (0.03, (2, 2, 0))],
        [(0.05, (0, 2, 0)), (0.06, (1, 0, 1)), (-0.04, (2, 0, 0)), (0.025, (1, 2, 1))],
        [(-0.06, (0, 0, 2)), (0.05, (1, 1, 0)), (0.03, (2, 0, 2))]],
    3: [[(0.03, (3, 0, 0)), (-0.05, (1, 1, 0)), (0.04, (0, 0, 2)), (0.025, (2, 1, 0)), (0.015, (3, 2, 1))],
        [(0.05, (0, 2, 0)), (0.03, (1, 0, 1)), (-0.025, (0, 3, 1)), (0.015, (2, 3, 0))],
        [(-0.04, (0, 0, 3)), (0.05, (1, 1, 0)), (0.015, (3, 1, 2))]],
}


def _poly_field(p, Y):
    """u (.., 3) and grad_X u (.., 3, 3) [a, m] = du_a/dX_m at parametric Y."""
    u = np.zeros(Y.shape[:-1] + (3,))
    gY = np.zeros(Y.shape[:-1] + (3, 3))
    for a, terms in enumerate(_FIELD[p]):
        for coef, e in terms:
            u[..., a] += coef * np.prod(Y ** np.array(e), axis=-1)
            for b in range(3):
                if e[b]:
                    e1 = np.array(e)
                    e1[b] -= 1
                    gY[..., a, b] += coef * e[b] * np.prod(Y ** e1, axis=-1)
    return u, gY @ S_INV  # dY/dX = S^-1


def _sheared_gauss(box, ng):
    """Tensor Gauss rule on the cells of a ``_ShearedBox`` from numpy's
    leggauss: parametric points Y (nc, nq, 3) and physical weights (nc, nq)
    = |det S| * (cell volume / 8) * w_x w_y w_z. Cell order x fastest, then
    y, z; QP order q = gx + ng (gy + ng gz) (the order the model documents)."""
    xg, wg = np.polynomial.legendre.leggauss(ng)
    br = [np.unique(k) for k in box.knots]
    detS = abs(float(np.linalg.det(S_AFFINE)))
    wz, wy, wx = np.meshgrid(wg, wg, wg, indexing="ij")
    w3 = (wx * wy * wz).ravel()
    Ys, Ws = [], []
    for z0, z1 in zip(br[2][:-1], br[2][1:]):
        for y0, y1 in zip(br[1][:-1], br[1][1:]):
            for x0, x1 in zip(br[0][:-1], br[0][1:]):
                Z, Yy, Xx = np.meshgrid(*(lo + 0.5 * (xg + 1.0) * (hi - lo) for lo, hi in ((z0, z1), (y0, y1), (x0, x1))),
                                        indexing="ij")
                Ys.append(np.column_stack([Xx.ravel(), Yy.ravel(), Z.ravel()]))
                Ws.append(w3 * 0.125 * (x1 - x0) * (y1 - y0) * (z1 - z0) * detS)
    return np.stack(Ys), np.stack(Ws)


def _l2_project(Xq, w, f, deg, X_eval=None):
    """Cell-wise L2 projection onto P_deg(X) in the inner product
    sum_q w_q f_q g_q, as weighted least squares (``np.linalg.lstsq`` on
    sqrt(w)-scaled monomials centred / scaled per cell), evaluated at X_eval
    (default: Xq). Written independently of ``mixed_p1.orthonormal_basis``."""
    X_eval = Xq if X_eval is None else X_eval
    exps = np.array([e for e in itertools.product(range(deg + 1), repeat=3) if sum(e) <= deg])
    out = np.empty(X_eval.shape[:2])
    for c in range(Xq.shape[0]):
        ctr, scl = Xq[c].mean(axis=0), Xq[c].std(axis=0)
        mono = lambda Z: np.prod(((Z - ctr) / scl)[:, None, :] ** exps[None], axis=2)  # noqa: E731
        sw = np.sqrt(w[c])
        coef = np.linalg.lstsq(mono(Xq[c]) * sw[:, None], f[c] * sw, rcond=None)[0]
        out[c] = mono(X_eval[c]) @ coef
    return out


def _voigt(T):
    return np.stack([T[..., 0, 0], T[..., 1, 1], T[..., 2, 2], T[..., 0, 1], T[..., 1, 2], T[..., 0, 2]], axis=-1)


def _poly_case(p, space, nu):
    """Sheared p-degree C^0 box (2 x 2 x 1 cells), model, control values of
    the polynomial field, and the independent Gauss rule (checked against
    the model's QP positions and weights)."""
    box = _ShearedBox(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(2, 2, 1), degree=p, continuity=0)
    model = MP.build_model(box, p, E=1.0, nu=nu, pressure=space)
    U = _greville_interpolant(box, lambda X: _poly_field(p, X @ S_INV.T)[0]).ravel()
    Yq, wq = _sheared_gauss(box, p + 1)
    # same points and weights as the model (round-off of the affine map)
    assert np.max(np.abs(model["X_all"] - Yq @ S_AFFINE.T)) < 1e-13
    assert _rel(model["W_all"], wq) < 1e-13
    return box, model, U, Yq, wq


@pytest.mark.parametrize("p,space", [(2, "P0"), (2, "P1"), (3, "P1"), (3, "P2")])
def test_linearized_stress_is_closed_form_with_projected_divergence(p, space):
    """cauchy_stress_linearized against its derivation, on a field with a
    non-constant div u outside the pressure space.

    Derivation (linearising sigma = mu/J B + phat I, phat = Pi U'(Pi J), at
    u = 0): J = 1 + div u, B = I + grad u + grad u^T, U'(1) = -mu,
    U''(1) = lam + mu, Pi reproduces constants, so
    phat = -mu + (lam + mu) Pi(div u) and
    sigma = mu (grad u + grad u^T) - mu div u I + (lam + mu) Pi(div u) I.
    Reference: grad u analytic, Pi(div u) by ``_l2_project`` with a
    (p+3)-point rule. Both that rule and the model's (p+1)-point rule
    integrate the projection integrands exactly (affine cells, constant
    Jacobian; div u of degree <= p and P_k of degree <= k <= p per
    parametric direction, products <= 2p + 1), so the reference is the exact
    cell L2 projection. nu = 0.3 (lam ~ mu): for a prescribed field with
    div u = O(|grad u|) the -mu (div u - Pi div u) term is O(mu / lam) of the
    stress, ~2e-5 at nu = 0.49999. (On the solved V07 plate div u_h is not
    small and 'vol = lam * Pi div u' moves e_sigma_xx by 3 - 33 %, but
    downwards, so the gate and ordering asserts there could not see it, S49
    §8.) Tolerance 1e-12 of max |sigma|: round-off of the Greville solves and
    the projection (measured <= 1e-14), while the distinguishing terms are
    checked to be >= 1e-3 below. The finite-difference consistency with
    ``cauchy_stress`` (tolerance 1e-6: truncation O(s^2) at s = 1e-4,
    measured <= 4e-9) is a cross-check between the two recoveries, not the
    reference."""
    nu = 0.3
    box, model, U, Yq, _ = _poly_case(p, space, nu)
    mu, lam = float(model["mu"]), float(model["lam"])
    Gu = _poly_field(p, Yq)[1]
    div = np.einsum("cqaa->cq", Gu)
    Yh, wh = _sheared_gauss(box, p + 3)
    div_h = np.einsum("cqaa->cq", _poly_field(p, Yh)[1])
    pdiv = _l2_project(Yh @ S_AFFINE.T, wh, div_h, _DEG[space], X_eval=Yq @ S_AFFINE.T)
    ref = _voigt(mu * (Gu + np.swapaxes(Gu, -1, -2))
                 + (-mu * div + (lam + mu) * pdiv)[..., None, None] * np.eye(3))
    got = MP.cauchy_stress_linearized(model, U)
    scale = float(np.max(np.abs(ref)))
    err = float(np.max(np.abs(got - ref))) / scale
    # power: Pi div u != div u by >= 1e-3 of the stress, so dropping the
    # projection ((lam + mu) (div - Pi div)) or the -mu div u term
    # (mu (Pi div - div)) moves sigma by >= 1e-3 * max |sigma|
    gap = mu * float(np.max(np.abs(pdiv - div))) / scale
    # consistency with the finite-strain recovery: d/ds cauchy_stress(s u) at
    # s = 0 by central differences (truncation O(s^2), s = 1e-4)
    s = 1e-4
    dS = (np.stack(MP.cauchy_stress(model, s * U)) - np.stack(MP.cauchy_stress(model, -s * U))) / (2 * s)
    err_fd = float(np.max(np.abs(dS - got))) / scale
    print(f"\n[lin stress p={p} {space}] rel err {err:.2e}, mu|Pi div - div|/|sigma| {gap:.2e}, "
          f"d/ds finite-strain {err_fd:.2e}")
    assert gap > 1e-3
    assert err < 1e-12
    assert err_fd < 1e-6


@pytest.mark.parametrize("p,space", [(2, "P1"), (3, "P1"), (3, "P2")])
def test_finite_strain_stress_projected_pressure_nonuniform_J(p, space):
    """cauchy_stress on a large, non-uniform deformation (J ~ 0.84 .. 1.87,
    varying by > 0.05 inside every cell, asserted) against
    sigma = mu/J F F^T + phat I, theta = Pi J, phat = Pi U'(theta),
    U'(t) = (lam ln t - mu) / t, computed here: F analytic, Pi by
    ``_l2_project`` with an independent Gauss rule equal to the model's
    (p+1)-point cell quadrature (Pi_Q is defined by that quadrature; the
    rule's identity with the model's X_all / W_all is asserted). nu = 0.45.
    Tolerance 1e-11 of max |sigma| (theta, phat likewise): round-off of the
    least-squares projections (measured <= 1e-13). Power: U'(theta) differs
    from Pi U'(theta) and J from theta by >= 1e-2 (asserted), so recovering
    the stress with the unprojected U'(theta) or theta = J fails."""
    nu = 0.45
    _, model, U, Yq, wq = _poly_case(p, space, nu)
    mu, lam = float(model["mu"]), float(model["lam"])
    Xq = Yq @ S_AFFINE.T
    F = np.eye(3) + _poly_field(p, Yq)[1]
    J = np.linalg.det(F)
    theta = _l2_project(Xq, wq, J, _DEG[space])
    assert np.all(theta > 0.0)
    Up = (lam * np.log(theta) - mu) / theta
    phat = _l2_project(Xq, wq, Up, _DEG[space])
    ref = _voigt(mu / J[..., None, None] * (F @ np.swapaxes(F, -1, -2)) + phat[..., None, None] * np.eye(3))
    S, th_m, ph_m = MP.cauchy_stress(model, U, with_pressure=True)
    scale = float(np.max(np.abs(ref)))
    err = float(np.max(np.abs(np.stack(S) - ref))) / scale
    jrange = float(np.min(np.ptp(J, axis=1)))
    gap_p = float(np.max(np.abs(Up - phat))) / scale
    gap_t = float(np.max(np.abs(J - theta)))
    print(f"\n[finite-strain stress p={p} {space}] J in [{J.min():.3f}, {J.max():.3f}], min cell range {jrange:.3f}; "
          f"rel err {err:.2e}; |U'(th) - Pi U'(th)|/|sigma| {gap_p:.2e}, |J - theta| {gap_t:.2e}")
    assert jrange > 0.05
    assert gap_p > 1e-2 and gap_t > 1e-2
    assert err < 1e-11
    assert float(np.max(np.abs(th_m - theta))) < 1e-11
    assert _rel(ph_m, phat) < 1e-11


# ---------------------------------------------------------------------------
# 4. V07 Kirsch (nu = 0.49999) on the precursor's body-fitted plate
# ---------------------------------------------------------------------------

def _mean_stress_floor(model, u, S, X, W, deg):
    """Mean-stress error e_p and its best-approximation floor eps - delta.

    p_ex = (1 + nu)(sigma_xx + sigma_yy)/3 (Kirsch, plane strain sigma_zz =
    nu (sigma_xx + sigma_yy)). The recovered mean stress is
    tr(sigma_h)/3 = q_h - mu/3 div u_h with q_h = (lam + mu) Pi div u_h in
    the cell-wise P_deg space, so (L2 norm with the QP weights W, in which
    Pi is the orthogonal projection, hence the best approximation)
    ||p_ex - tr(sigma_h)/3|| >= ||p_ex - Pi p_ex|| - mu/3 ||div u_h||,
    i.e. e_p >= eps - delta for every solution with pressure in that space.
    Pi p_ex by ``_l2_project``; div u_h from fbar's G_all (S40 kernel)."""
    assert _rel(model["W_all"], W) < 1e-12  # Pi is orthogonal in this inner product
    nu = PRE.NU_VAL
    sxx, syy, _ = PRE.kirsch_stress(X[..., 0], X[..., 1])
    p_ex = (1.0 + nu) * (sxx + syy) / 3.0
    nrm = float(np.sqrt(np.sum(W * p_ex ** 2)))
    e_p = float(np.sqrt(np.sum(W * (S[..., :3].sum(axis=-1) / 3.0 - p_ex) ** 2))) / nrm
    eps = float(np.sqrt(np.sum(W * (p_ex - _l2_project(X, W, p_ex, deg)) ** 2))) / nrm
    div = np.einsum("cqnm,cnm->cq", model["G_all"], np.asarray(u).reshape(-1, 3)[model["conn"]])
    delta = float(model["mu"]) / 3.0 * float(np.sqrt(np.sum(W * div ** 2))) / nrm
    return {"e_p": e_p, "eps": eps, "delta": delta}


def _kirsch(plate, p, pressure, floor=False):
    """Linearised V07 solve; errors with the precursor's helpers
    (``floor=True``: also :func:`_mean_stress_floor`)."""
    t0 = time.perf_counter()
    X, W, _ = PRE._qp_positions_and_weights(plate)
    vol_err = abs(float(W.sum()) - V_EXACT) / V_EXACT
    fixed, presc = PRE._kirsch_dirichlet_dofs(plate)
    model = MP.build_model(plate, p, E=PRE.E_VAL, nu=PRE.NU_VAL, pressure=pressure)
    assert np.max(np.abs(model["X_all"] - X)) < 1e-12  # same cells / QPs as the precursor
    u = MP.solve_linearized(model, np.zeros(model["n_dof"]), fixed, presc)
    S = MP.cauchy_stress_linearized(model, u)
    e_u, e_s = PRE._l2_errors(X, W, PRE._eval_u_at_qps(plate, u), S[..., 0])
    out = {"e_u": e_u, "e_s": e_s, "vol_err": vol_err, "ndof": int(model["n_dof"]),
           "free": int(model["n_dof"]) - len(fixed) - len(presc), "cells": int(X.shape[0])}
    if floor:
        out.update(_mean_stress_floor(model, u, S, X, W, _DEG[pressure]))
    out["t"] = time.perf_counter() - t0
    del model, u, S
    return out


def _slope(r0, r1, h0, h1, key):
    return float(np.log(r0[key] / r1[key]) / np.log(h0 / h1))


def _ols(rs, hs, key):
    return float(np.polyfit(np.log(hs), np.log([r[key] for r in rs]), 1)[0])


def test_v07_kirsch_p2c0_p1disc_gate_and_slopes_vs_p0():
    """p=2 C^0 / P1-disc (Q2/P1disc-like) at the v07.json sweep h = 0.5, 0.25,
    0.125: the v07.json per-level bounds on (e_u, e_sigma_xx) hold; P1 beats
    P0 on the same mesh and DOFs and beats the precursor's p=2 C^1 / P0
    (fbar) at about the same DOF count.

    P1 vs P0 without pinned numbers: (a) the P1 mean-stress error is below
    the floor that no solution with a cell-constant pressure can undercut on
    that mesh (best approximation of the Kirsch mean stress in P0, computed
    here; :func:`_mean_stress_floor`), and both runs respect their own floor
    (a theorem, so a check of the error evaluation); (b) the P0 / P1
    e_sigma_xx ratio grows under refinement, as the projection orders
    O(h) vs O(h^2) require. The plan's 2.2 sigma slope is reported, not
    asserted: with p=2 displacements the stress converges at best O(h^2)
    (measured 1.63 / 1.85, S49)."""
    hs = [0.5, 0.25, 0.125]
    P1, P0, C1 = [], [], []
    for h in hs:
        nu_c, nvh, nw = LEVEL[h]
        p2c0 = _Plate(nu_c, nvh, nw, degree=2, cont=0)
        P1.append(_kirsch(p2c0, 2, "P1", floor=True))
        P0.append(_kirsch(p2c0, 2, "P0", floor=True))
        C1.append(_kirsch(PRE.MappedQuarterPlate(n_u=nu_c, n_v_half=nvh, n_w=nw, degree=2), 2, "P0"))
    print("\n[V07 p=2] h | space | ndof | e_u | e_sigma_xx | vol_err | t | e_p | floor eps - delta")
    for h, a, b, c in zip(hs, P1, P0, C1):
        for name, r in (("C0/P1", a), ("C0/P0", b), ("C1/P0", c)):
            fl = f" {r['e_p']:.4e} {r['eps']:.4e} - {r['delta']:.2e}" if "e_p" in r else ""
            print(f"  {h:<6} {name} {r['ndof']:>6} {r['e_u']:.4e} {r['e_s']:.4e} {r['vol_err']:.1e} {r['t']:.1f}s{fl}")
    for name, rs in (("C0/P1", P1), ("C0/P0", P0), ("C1/P0", C1)):
        print(f"  slopes {name}: u {_slope(rs[0], rs[1], .5, .25, 'e_u'):.3f} / {_slope(rs[1], rs[2], .25, .125, 'e_u'):.3f}"
              f" (OLS {_ols(rs, hs, 'e_u'):.3f}), sigma_xx {_slope(rs[0], rs[1], .5, .25, 'e_s'):.3f} / "
              f"{_slope(rs[1], rs[2], .25, .125, 'e_s'):.3f} (OLS {_ols(rs, hs, 'e_s'):.3f}); "
              f"plan targets {SLOPE_U} / {SLOPE_S}")
    for h, r, r0 in zip(hs, P1, P0):
        tol_u, tol_s = GATE[h]
        assert r["vol_err"] <= VOL_TOL
        assert r["e_u"] <= tol_u, f"h={h}: e_u={r['e_u']:.3e}"
        assert r["e_s"] <= tol_s, f"h={h}: e_sigma_xx={r['e_s']:.3e}"
        assert r["e_s"] < r0["e_s"]  # same mesh, same displacement DOFs
        # (a) best-approximation floors (theorem) and P1 below the P0 floor
        assert r0["e_p"] >= r0["eps"] - r0["delta"]
        assert r["e_p"] >= r["eps"] - r["delta"]
        assert r["e_p"] < r0["eps"] - r0["delta"], f"h={h}: P1 e_p={r['e_p']:.3e} not below the P0 floor"
    # (b) P0 / P1 ratio ~ h^-1 (projection orders): grows at every refinement
    ratio = [r0["e_s"] / r["e_s"] for r, r0 in zip(P1, P0)]
    print("  e_sigma_xx ratio P0/P1: " + ", ".join(f"{q:.2f}" for q in ratio))
    assert ratio[0] < ratio[1] < ratio[2]
    # equal-DOF comparison: C0/P1 at h=0.25 (19,500 DOFs) vs fbar C1/P0 at h=0.125 (20,904)
    assert P1[1]["ndof"] <= C1[2]["ndof"] and P1[1]["e_s"] < C1[2]["e_s"]
    # rate: P1 projection error O(h^2) vs P0 O(h) -> at least half the gap
    assert _ols(P1, hs, "e_s") >= _ols(P0, hs, "e_s") + 0.5


def test_v07_kirsch_p3c0_p2disc_slopes_on_h_0p5_0p333_0p25_not_v07_sweep():
    """p=3 C^0 / P2-disc (Q3/P2disc-like, inf-sup stable, S49) on the levels
    h = 0.5, 1/3, 0.25 -- NOT the v07.json sweep h = 0.5, 0.25, 0.125: p=3
    at h = 0.125 needs ~4 GB (over this PC's 1.5 GB run limit), so it is not
    run and nothing is claimed for it.

    Asserted: the v07.json per-level bounds at the two v07.json levels in
    this range (0.5, 0.25); slope_u >= 2.3 and slope_sigma_xx >= 2.2 between
    the end points h = 0.5 and 0.25; the 3-point OLS sigma_xx slope over
    0.5, 1/3, 0.25 >= 2.2. Whether the slopes hold over the v07.json range
    0.5 .. 0.125 is not measured. The slope asserts rest on theory (O(h^3)
    sigma, O(h^4) u for this pair) and on the measured local sigma slopes
    rising with refinement (2.33 -> 2.53, S49)."""
    hs = [0.5, 1.0 / 3.0, 0.25]
    rs = [_kirsch(_Plate(*LEVEL[h], degree=3, cont=0), 3, "P2") for h in hs]
    for h, r in zip(hs, rs):
        print(f"\n[V07 p=3 C0/P2] h={h:.4f} ndof={r['ndof']} e_u={r['e_u']:.4e} e_sigma_xx={r['e_s']:.4e} "
              f"vol_err={r['vol_err']:.1e} t={r['t']:.1f}s")
        assert r["vol_err"] <= VOL_TOL
        if h in GATE:
            assert r["e_u"] <= GATE[h][0] and r["e_s"] <= GATE[h][1]
    s_u, s_s = _slope(rs[0], rs[2], 0.5, 0.25, "e_u"), _slope(rs[0], rs[2], 0.5, 0.25, "e_s")
    loc = [_slope(rs[i], rs[i + 1], hs[i], hs[i + 1], "e_s") for i in range(2)]
    print(f"  slopes h0->h1: u {s_u:.3f}, sigma_xx {s_s:.3f}; local sigma_xx {loc[0]:.3f}, {loc[1]:.3f}; "
          f"OLS sigma_xx {_ols(rs, hs, 'e_s'):.3f}")
    assert s_u >= SLOPE_U
    assert s_s >= SLOPE_S
    assert _ols(rs, hs, "e_s") >= SLOPE_S


def test_v07_kirsch_p2c1_p1_locks():
    """On the precursor's own space (p=2 C^1) the P1 pair is not inf-sup
    stable (about 2 in-plane DOFs per cell for 3 in-plane P1 modes; beta_h ~
    2e-5, S49): the stress error exceeds 100 % while P0 on the same mesh
    passes the gate."""
    for h in (0.5, 0.25):
        plate = PRE.MappedQuarterPlate(n_u=LEVEL[h][0], n_v_half=LEVEL[h][1], n_w=LEVEL[h][2], degree=2)
        r1, r0 = _kirsch(plate, 2, "P1"), _kirsch(plate, 2, "P0")
        print(f"\n[V07 p=2 C1] h={h}: P1 e_sigma_xx={r1['e_s']:.3e}, P0 e_sigma_xx={r0['e_s']:.3e}")
        assert r0["e_s"] <= GATE[h][1]
        assert r1["e_s"] > 1.0
