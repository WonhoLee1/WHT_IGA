# -*- coding: utf-8 -*-
"""Mean-dilatation F-bar polynomial hyperelasticity (Mooney-Rivlin, Yeoh,
full polynomial N = 3 for the energy/tangent/patch checks), S48.

Every reference below is computed in this file, independently of the kernel
in ``wht_iga.analyses.hyperelastic_fbar``:

* the strain energy W(F) is coded here directly from the Abaqus formulas
  (``_w_iso`` / ``_u_vol``); stresses come from its complex-step derivative
  P = Im W(F + i h E) / h (exact to round-off, no truncation error) and
  sigma = P F^T / J;
* the residual is checked against the complex-step gradient of the
  mean-dilatation energy Pi(u) assembled here from the model's QP data;
* the tangent is checked against central finite differences of the residual;
* homogeneous solves (uniaxial, equibiaxial, simple shear) are compared with
  the stress at the deformation solved here by bisection (traction-free
  lateral faces);
* incompressible Mooney-Rivlin uniaxial: 2 (lam - 1/lam^2)(C10 + C01/lam);
* Neo-Hookean (C01 = 0): hand-derived Abaqus closed form
  sigma = (2 C10 / J) dev(Bbar) + 2 (J - 1) / D1 I.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses import hyperelastic_fbar as HF
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.cook import CookSpline

MR = {"kind": "mr", "C10": 0.3, "C01": 0.2, "D1": 0.4}  # mu0 = 1, K0 = 5
# Yeoh with all six terms active (C20 < 0 as in typical rubber fits);
# D2, D3 small enough that (J-1)^4 and (J-1)^6 matter at |J-1| ~ 0.3.
YEOH = {"kind": "yeoh", "C10": 0.5, "C20": -0.05, "C30": 0.02, "D1": 0.4, "D2": 0.3, "D3": 0.5}
# Full Abaqus polynomial N = 3 (every C_ij, i + j <= 3): the only material
# here with W12 and W22 != 0 (mixed and (I2b-3)^j, j >= 2 terms).
POLY = {"kind": "poly", "D1": 0.4, "D2": 0.3, "D3": 0.5,
        "coeffs": {(1, 0): 0.3, (0, 1): 0.15, (2, 0): 0.04, (1, 1): -0.02, (0, 2): 0.03,
                   (3, 0): 0.01, (2, 1): 0.005, (1, 2): -0.004, (0, 3): 0.002}}
SPECS = {"mr": MR, "yeoh": YEOH, "poly": POLY}

# A general (non-symmetric) finite deformation gradient, det ~ 1.30
F_GEN = np.array([[1.25, 0.15, -0.05], [0.05, 1.10, 0.10], [0.0, -0.08, 0.95]])


# --------------------------------------------------------------------------
# independent reference energy (complex-capable, batched over leading axes)
# --------------------------------------------------------------------------

def _w_iso(F, s):
    J = np.linalg.det(F)
    C = np.swapaxes(F, -1, -2) @ F
    I1 = np.trace(C, axis1=-2, axis2=-1)
    I2 = 0.5 * (I1 * I1 - np.sum(C * C, axis=(-2, -1)))  # C symmetric: sum C_ij^2 = tr C^2
    a = J ** (-2.0 / 3.0) * I1 - 3.0
    b = J ** (-4.0 / 3.0) * I2 - 3.0
    if s["kind"] == "mr":
        return s["C10"] * a + s["C01"] * b
    if s["kind"] == "nh":
        return s["C10"] * a
    if s["kind"] == "poly":
        return sum(c * a ** i * b ** j for (i, j), c in s["coeffs"].items())
    return s["C10"] * a + s["C20"] * a ** 2 + s["C30"] * a ** 3


def _u_vol(J, s):
    U = (J - 1.0) ** 2 / s["D1"]
    if "D2" in s:
        U = U + (J - 1.0) ** 4 / s["D2"] + (J - 1.0) ** 6 / s["D3"]
    return U


def _w(F, s):
    return _w_iso(F, s) + _u_vol(np.linalg.det(F), s)


def _cauchy_ref(F, s, h=1e-30):
    """sigma = P F^T / J with P from the complex-step derivative of W."""
    F = np.asarray(F, dtype=float)
    P = np.zeros((3, 3))
    for a in range(3):
        for M in range(3):
            Fc = F.astype(complex)
            Fc[a, M] += 1j * h
            P[a, M] = float(np.imag(_w(Fc, s))) / h
    return P @ F.T / np.linalg.det(F)


def _material(s):
    if s["kind"] == "mr":
        return HF.mooney_rivlin(s["C10"], s["C01"], s["D1"])
    if s["kind"] == "nh":
        return HF.mooney_rivlin(s["C10"], 0.0, s["D1"])
    if s["kind"] == "poly":
        return HF.polynomial(s["coeffs"], (s["D1"], s["D2"], s["D3"]))
    return HF.yeoh(s["C10"], s["C20"], s["C30"], s["D1"], s["D2"], s["D3"])


def _voigt_to_mat(r):
    return np.array([[r[0], r[3], r[5]], [r[3], r[1], r[4]], [r[5], r[4], r[2]]])


def _bisect(fun, lo, hi, n=200):
    flo = fun(lo)
    assert flo * fun(hi) < 0.0, "bracket"
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        fm = fun(mid)
        if fm == 0.0 or hi - lo < 1e-16:
            return mid
        if (fm < 0.0) == (flo < 0.0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------
# models and states
# --------------------------------------------------------------------------

def _multicell_cases():
    """C^1 on the distorted Cook map (2x2x1 cells) and C^0 on a box (2x1x2)."""
    return {
        "cook_C1": CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 1), degree=2),
        "box_C0": BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 0.8, 0.6), n_cells=(2, 1, 2), degree=2,
                            continuity=0),
    }


def _nonuniform_state(box, seed=3):
    """Finite affine stretch/shear (det F ~ 1.3) plus a random field, so Jbar
    differs per cell and every term of W (incl. (I1b-3)^3, (J-1)^6) matters."""
    cp = np.asarray(box.control_points())
    ext = cp.max(axis=0) - cp.min(axis=0)
    rng = np.random.default_rng(seed)
    pert = 0.03 * rng.standard_normal(cp.shape) * ext.min()
    return (cp @ (F_GEN - np.eye(3)).T + pert).reshape(-1)


def _pi_ref(model, u, s):
    """Mean-dilatation energy sum_c [sum_q w_q W_iso(F_q) + V_c U(Jbar_c)] (complex ok)."""
    uu = np.asarray(u).reshape(-1, 3)
    total = 0.0
    for c in range(model["conn"].shape[0]):
        G, W = model["G_all"][c], model["W_all"][c]
        F = np.eye(3) + np.einsum("ia,qim->qam", uu[model["conn"][c]], G)
        Jbar = np.sum(W * np.linalg.det(F)) / np.sum(W)
        total = total + np.sum(W * _w_iso(F, s)) + np.sum(W) * _u_vol(Jbar, s)
    return total


# --------------------------------------------------------------------------
# 1. residual = gradient of the mean-dilatation energy; tangent = dR/du
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", ["cook_C1", "box_C0"])
@pytest.mark.parametrize("mat", ["mr", "yeoh", "poly"])
def test_residual_is_complex_step_gradient_of_mean_dilatation_energy(mat, case):
    """R(u) == dPi/du with Pi coded here (pointwise W_iso, cell-Jbar U).

    Complex step has no truncation error and no cancellation, so both sides
    are exact to round-off; 1e-11 of max|R| leaves ~1000x the accumulated
    double-precision error of a few hundred QP terms. A pointwise J in U
    (no mean dilatation) changes R by O(|J - Jbar| / D1) ~ 1e-2 here.
    """
    s = SPECS[mat]
    box = _multicell_cases()[case]
    model = HF.build_model(box, 2, _material(s))
    u = _nonuniform_state(box)
    Jb = HF.cell_jbar(model, u)
    assert np.ptp(Jb) > 1e-2 and np.min(Jb) > 1.2  # genuinely non-uniform, finite strain
    R, _ = HF.residual_tangent_sparse(model, u)
    h = 1e-30
    g = np.empty(u.size)
    for k in range(u.size):
        uc = u.astype(complex)
        uc[k] += 1j * h
        g[k] = float(np.imag(_pi_ref(model, uc, s))) / h
    assert float(np.max(np.abs(R - g))) <= 1e-11 * float(np.max(np.abs(g)))


@pytest.mark.parametrize("case", ["cook_C1", "box_C0"])
@pytest.mark.parametrize("mat", ["mr", "yeoh", "poly"])
def test_tangent_matches_central_fd_multicell(mat, case):
    """CSR tangent vs central FD of R on a non-uniform multi-cell state.

    h = 1e-6 * max|u|: truncation O(h^2 R''') and round-off O(eps |R| / h)
    both stay near 1e-9 of max|K| (measured 1e-10 .. 1e-9 for h = 1e-5..1e-7
    relative), so rel < 1e-6 is the required gate with margin. The tangent is
    the Hessian of Pi, hence symmetric (round-off level).
    """
    s = SPECS[mat]
    box = _multicell_cases()[case]
    model = HF.build_model(box, 2, _material(s))
    u = _nonuniform_state(box)
    _, K = HF.residual_tangent_sparse(model, u)
    K = K.toarray()
    h = 1e-6 * float(np.max(np.abs(u)))
    Kfd = np.zeros_like(K)
    for j in range(u.size):
        up, um = u.copy(), u.copy()
        up[j] += h
        um[j] -= h
        Kfd[:, j] = (HF.residual_tangent_sparse(model, up)[0]
                     - HF.residual_tangent_sparse(model, um)[0]) / (2 * h)
    scale = float(np.max(np.abs(Kfd)))
    assert float(np.max(np.abs(K - Kfd))) / scale < 1e-6
    assert float(np.max(np.abs(K - K.T))) / scale < 1e-12


# --------------------------------------------------------------------------
# 2. affine patch on a distorted map; stress output vs weak form
# --------------------------------------------------------------------------

@pytest.mark.parametrize("mat", ["mr", "yeoh", "poly"])
def test_affine_patch_cook_stress_and_interior_equilibrium(mat):
    """Affine u = (F - I) X on the tapered Cook map (C^1, 2x2x1): every QP
    stress equals the complex-step reference at F, and the interior control
    points carry zero residual (constant P is in equilibrium; p+1 Gauss
    points integrate the bilinear-map integrand exactly). Both are
    round-off-exact, hence 1e-12."""
    s = SPECS[mat]
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 1), degree=2)
    model = HF.build_model(box, 2, _material(s))
    cp = np.asarray(box.control_points())
    u = (cp @ (F_GEN - np.eye(3)).T + np.array([0.3, -0.2, 0.1])).reshape(-1)
    ref = _cauchy_ref(F_GEN, s)
    for S in HF.cauchy_stress(model, u):
        assert S.shape == (27, 6)
        for r in S:
            assert np.max(np.abs(_voigt_to_mat(r) - ref)) <= 1e-12 * np.max(np.abs(ref))
    R, _ = HF.residual_tangent_sparse(model, u)
    nu_, nv_, nw_ = box.n_cp
    interior = [a + nu_ * (b + nv_ * c) for c in range(1, nw_ - 1) for b in range(1, nv_ - 1)
                for a in range(1, nu_ - 1)]
    assert len(interior) == 2 * 2 * 1  # basis functions vanishing on the whole boundary
    Ri = R.reshape(-1, 3)[interior]
    assert float(np.max(np.abs(Ri))) <= 1e-12 * float(np.max(np.abs(R)))


@pytest.mark.parametrize("mat", ["mr", "yeoh", "poly"])
def test_qp_stress_output_is_the_stress_the_residual_equilibrates(mat):
    """Non-uniform state: R_I == sum_q w_q J_q sigma_q F_q^-T G_qI with sigma_q
    from ``cauchy_stress`` and F_q recomputed here. Catches a stress output
    that uses Jbar where the weak form uses J_q (invisible in homogeneous
    tests, where Jbar == J). Same arithmetic, so 1e-12."""
    s = SPECS[mat]
    box = _multicell_cases()["box_C0"]
    model = HF.build_model(box, 2, _material(s))
    u = _nonuniform_state(box, seed=8)
    R, _ = HF.residual_tangent_sparse(model, u)
    S = HF.cauchy_stress(model, u)
    uu = u.reshape(-1, 3)
    Rr = np.zeros_like(uu)
    for c in range(model["conn"].shape[0]):
        G, W, nodes = model["G_all"][c], model["W_all"][c], model["conn"][c]
        for q in range(W.size):
            F = np.eye(3) + uu[nodes].T @ G[q]
            P = np.linalg.det(F) * _voigt_to_mat(S[c][q]) @ np.linalg.inv(F).T
            Rr[nodes] += W[q] * G[q] @ P.T
    assert float(np.max(np.abs(Rr.reshape(-1) - R))) <= 1e-12 * float(np.max(np.abs(R)))


# --------------------------------------------------------------------------
# 3. homogeneous deformations solved by Newton vs independent references
# --------------------------------------------------------------------------

LX, LY, LZ = 1.0, 0.8, 0.6
# Gate for the Newton solves: dtol = 1e-13 bounds the strain error near 1e-13;
# times K0/mu0 <= 5 (1e4 in the incompressible test is checked separately)
# and the round-off of the assembly (~1e-15 * cond) this stays < 1e-11, and
# the references (complex step, bisection to 1e-16) are round-off-exact.
# 1e-9 of the stress scale leaves a 100x margin.
TOL_HOMOG = 1e-9


def _homog_model(s, n_cells=(2, 2, 1)):
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(LX, LY, LZ), n_cells=n_cells, degree=2)
    return HF.build_model(box, 2, _material(s)), np.asarray(box.control_points())


def _symmetry_fixed(cp):
    """u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0 (no rigid modes)."""
    fixed = []
    for I, x in enumerate(cp):
        for d in range(3):
            if abs(x[d]) < 1e-12:
                fixed.append(3 * I + d)
    return np.asarray(fixed, dtype=int)


def _solve(model, fixed, prescribed, n_steps=4):
    return HF.solve_static(model, np.zeros(model["n_dof"]), fixed, n_steps=n_steps,
                           prescribed=prescribed, rtol=1e-13, atol=1e-14, dtol=1e-13)


def _assert_all_qp(model, u, ref):
    for S in HF.cauchy_stress(model, u):
        for r in S:
            err = np.max(np.abs(_voigt_to_mat(r) - ref))
            assert err <= TOL_HOMOG * np.max(np.abs(ref)), (err, ref)


@pytest.mark.parametrize("mat", ["mr", "yeoh"])
@pytest.mark.parametrize("lam1", [1.6, 0.7])
def test_uniaxial_matches_bisection_reference(mat, lam1):
    """Uniaxial stretch lam1 (tension and compression), lateral faces traction
    free. Reference: lateral stretch from sigma_22(diag(lam1, l2, l2)) = 0 by
    bisection on the complex-step stress. Checks every QP stress tensor, the
    lateral contraction and the nominal stress from the reaction force."""
    s = SPECS[mat]
    model, cp = _homog_model(s)
    l2 = _bisect(lambda x: _cauchy_ref(np.diag([lam1, x, x]), s)[1, 1], 0.5, 1.5)
    ref = _cauchy_ref(np.diag([lam1, l2, l2]), s)
    assert abs(ref[1, 1]) <= 1e-12 * abs(ref[0, 0]) and abs(ref[0, 0]) > 0.1
    tip = [I for I, x in enumerate(cp) if abs(x[0] - LX) < 1e-12]
    u = _solve(model, _symmetry_fixed(cp), {3 * I: (lam1 - 1.0) * LX for I in tip})
    _assert_all_qp(model, u, ref)
    top = [I for I, x in enumerate(cp) if abs(x[1] - LY) < 1e-12]
    assert np.max(np.abs(u[[3 * I + 1 for I in top]] - (l2 - 1.0) * LY)) <= TOL_HOMOG * LY
    R, _ = HF.residual_tangent_sparse(model, u)
    P11 = float(np.sum(R[[3 * I for I in tip]])) / (LY * LZ)  # nominal = reaction / A0
    assert abs(P11 - ref[0, 0] * l2 * l2) <= TOL_HOMOG * abs(ref[0, 0])


def test_newton_without_dirichlet_predictor_reaches_the_same_state():
    """``predictor=False`` (the fbar.solve_static start: new boundary values on
    the old interior) must converge to the same homogeneous state as the
    default predictor path; both are compared with the bisection reference
    for the lateral stretch (same TOL_HOMOG reasoning)."""
    s, lam1 = MR, 1.6
    model, cp = _homog_model(s)
    l2 = _bisect(lambda x: _cauchy_ref(np.diag([lam1, x, x]), s)[1, 1], 0.5, 1.5)
    tip = [I for I, x in enumerate(cp) if abs(x[0] - LX) < 1e-12]
    top = [I for I, x in enumerate(cp) if abs(x[1] - LY) < 1e-12]
    us = {}
    for pred in (False, True):
        hist = []
        us[pred] = HF.solve_static(model, np.zeros(model["n_dof"]), _symmetry_fixed(cp), n_steps=4,
                                   prescribed={3 * I: (lam1 - 1.0) * LX for I in tip},
                                   rtol=1e-13, atol=1e-14, dtol=1e-13, history=hist, predictor=pred)
        assert hist[-1][0] == 1.0
        assert np.max(np.abs(us[pred][[3 * I + 1 for I in top]] - (l2 - 1.0) * LY)) <= TOL_HOMOG * LY
    assert np.max(np.abs(us[True] - us[False])) <= TOL_HOMOG * np.max(np.abs(us[False]))


@pytest.mark.parametrize("mat", ["mr", "yeoh"])
def test_equibiaxial_matches_bisection_reference(mat):
    """lam = 1.4 in x and y, z face traction free (thickness stretch by bisection)."""
    s = SPECS[mat]
    lam = 1.4
    model, cp = _homog_model(s)
    l3 = _bisect(lambda x: _cauchy_ref(np.diag([lam, lam, x]), s)[2, 2], 0.3, 1.2)
    ref = _cauchy_ref(np.diag([lam, lam, l3]), s)
    assert abs(ref[2, 2]) <= 1e-12 * abs(ref[0, 0]) and abs(ref[0, 0] - ref[1, 1]) <= 1e-12 * abs(ref[0, 0])
    pres = {}
    for I, x in enumerate(cp):
        if abs(x[0] - LX) < 1e-12:
            pres[3 * I] = (lam - 1.0) * LX
        if abs(x[1] - LY) < 1e-12:
            pres[3 * I + 1] = (lam - 1.0) * LY
    u = _solve(model, _symmetry_fixed(cp), pres)
    _assert_all_qp(model, u, ref)
    top = [I for I, x in enumerate(cp) if abs(x[2] - LZ) < 1e-12]
    assert np.max(np.abs(u[[3 * I + 2 for I in top]] - (l3 - 1.0) * LZ)) <= TOL_HOMOG * LZ


@pytest.mark.parametrize("mat", ["mr", "yeoh"])
def test_simple_shear_full_stress_tensor(mat):
    """F = I + gamma e1 (x) e2 (gamma = 0.6, J = 1): all boundary control points
    prescribed to the affine field, the 2x2x2-cell interior solved by Newton.
    Reference: the complex-step tensor, including the normal stresses
    (Poynting effect), plus the hand formulas for J = 1 simple shear with
    W1 = dW/dI1b, W2 = dW/dI2b at I1b - 3 = I2b - 3 = gamma^2:
    sigma_12 = 2 (W1 + W2) gamma, sigma_11 - sigma_22 = 2 (W1 + W2) gamma^2."""
    s = SPECS[mat]
    g = 0.6
    model, cp = _homog_model(s, n_cells=(2, 2, 2))
    F = np.eye(3)
    F[0, 1] = g
    ref = _cauchy_ref(F, s)
    assert abs(ref[0, 0]) > 0.05 and abs(ref[0, 1]) > 0.3
    hi = np.array([LX, LY, LZ])
    pres = {}
    for I, x in enumerate(cp):
        if np.any(np.abs(x) < 1e-12) or np.any(np.abs(x - hi) < 1e-12):
            pres[3 * I] = g * x[1]
            pres[3 * I + 1] = 0.0
            pres[3 * I + 2] = 0.0
    assert model["n_dof"] - len(pres) == 3 * 8  # 2x2x2 interior control points free
    u = _solve(model, np.zeros(0, dtype=int), pres)
    _assert_all_qp(model, u, ref)
    a = g * g
    if s["kind"] == "mr":
        w12 = s["C10"] + s["C01"]
    else:
        w12 = s["C10"] + 2.0 * s["C20"] * a + 3.0 * s["C30"] * a * a
    assert abs(ref[0, 1] - 2.0 * w12 * g) <= 1e-12 * abs(ref[0, 1])
    assert abs(ref[0, 0] - ref[1, 1] - 2.0 * w12 * g * g) <= 1e-12 * abs(ref[0, 1])


# --------------------------------------------------------------------------
# 4. incompressible limit and Neo-Hookean special case
# --------------------------------------------------------------------------

def test_mooney_rivlin_incompressible_limit_uniaxial():
    """Small D1: nominal uniaxial stress -> 2 (lam - 1/lam^2)(C10 + C01/lam).

    The compressible solution is a regular perturbation in D1, so the error
    must fall by 10x per decade of K0/mu0: gate ratio in [8, 12] (+-20 % for
    the O(D1^2) remainder). Size, first-order estimate: J - 1 = p/K0 with
    p = sigma_11/3 = 0.49 mu0 (lateral sigma = 0, sigma_iso deviatoric), and
    d ln P / d ln J = -1 - (1/3) d ln sigma / d ln lam + 1 = -1.04 for this MR
    at lam = 1.5, so error ~ 0.51 mu0/K0; gate 1.0 mu0/K0 (2x for the
    estimate). Nominal stress from the reaction force (weak form) and from
    the QP Cauchy stress times the lateral area ratio.
    """
    C10, C01, lam = 0.4, 0.1, 1.5
    mu0 = 2.0 * (C10 + C01)
    P_inc = 2.0 * (lam - 1.0 / lam ** 2) * (C10 + C01 / lam)
    errs = []
    for ratio in (1e2, 1e3, 1e4):
        s = {"kind": "mr", "C10": C10, "C01": C01, "D1": 2.0 / (ratio * mu0)}
        model, cp = _homog_model(s)
        tip = [I for I, x in enumerate(cp) if abs(x[0] - LX) < 1e-12]
        u = _solve(model, _symmetry_fixed(cp), {3 * I: (lam - 1.0) * LX for I in tip}, n_steps=8)
        R, _ = HF.residual_tangent_sparse(model, u)
        P11 = float(np.sum(R[[3 * I for I in tip]])) / (LY * LZ)
        top = [I for I, x in enumerate(cp) if abs(x[1] - LY) < 1e-12]
        l2 = 1.0 + float(np.mean(u[[3 * I + 1 for I in top]])) / LY
        sig11 = float(np.mean([S[:, 0].mean() for S in HF.cauchy_stress(model, u)]))
        assert abs(sig11 * l2 * l2 - P11) <= 1e-8 * P11
        e = abs(P11 - P_inc) / P_inc
        assert e < 1.0 / ratio, (ratio, e)
        errs.append(e)
    for e0, e1 in zip(errs, errs[1:]):
        assert 8.0 < e0 / e1 < 12.0, errs


def _nh_hand_cauchy(F, C10, D1):
    """Abaqus NEO HOOKE closed form (hand-derived, not via W):
    sigma = (2 C10 / J) (Bbar - tr(Bbar)/3 I) + (2/D1)(J - 1) I."""
    J = np.linalg.det(F)
    Bb = J ** (-2.0 / 3.0) * F @ F.T
    return 2.0 * C10 / J * (Bb - np.trace(Bb) / 3.0 * np.eye(3)) + 2.0 / D1 * (J - 1.0) * np.eye(3)


def test_neo_hookean_special_case_matches_hand_formula():
    """Mooney-Rivlin with C01 = 0 == ``neo_hookean`` == Abaqus NH closed form:
    (a) bit-identical R, K of the two constructors on a non-uniform state,
    (b) affine general F on the Cook map: QP stress vs the hand formula,
    (c) uniaxial Newton solve vs the hand formula with the lateral stretch by
        bisection of the hand sigma_22 (no use of W or its derivatives)."""
    C10, D1 = 0.35, 0.25
    s = {"kind": "nh", "C10": C10, "D1": D1}
    box = _multicell_cases()["box_C0"]
    m_mr = HF.build_model(box, 2, HF.mooney_rivlin(C10, 0.0, D1))
    m_nh = HF.build_model(box, 2, HF.neo_hookean(C10, D1))
    u = _nonuniform_state(box)
    (R1, K1), (R2, K2) = HF.residual_tangent_sparse(m_mr, u), HF.residual_tangent_sparse(m_nh, u)
    assert np.array_equal(R1, R2) and np.array_equal(K1.data, K2.data)

    cook = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 1, 1), degree=2)
    m = HF.build_model(cook, 2, _material(s))
    cp = np.asarray(cook.control_points())
    ref = _nh_hand_cauchy(F_GEN, C10, D1)
    for S in HF.cauchy_stress(m, (cp @ (F_GEN - np.eye(3)).T).reshape(-1)):
        for r in S:
            assert np.max(np.abs(_voigt_to_mat(r) - ref)) <= 1e-12 * np.max(np.abs(ref))

    lam1 = 1.7
    l2 = _bisect(lambda x: _nh_hand_cauchy(np.diag([lam1, x, x]), C10, D1)[1, 1], 0.5, 1.2)
    ref = _nh_hand_cauchy(np.diag([lam1, l2, l2]), C10, D1)
    model, cp = _homog_model(s)
    tip = [I for I, x in enumerate(cp) if abs(x[0] - LX) < 1e-12]
    u = _solve(model, _symmetry_fixed(cp), {3 * I: (lam1 - 1.0) * LX for I in tip})
    _assert_all_qp(model, u, ref)


def test_material_constructors_reject_invalid_input():
    """D1 is required and must be > 0 (no hybrid formulation); mu0 > 0."""
    with pytest.raises(ValueError):
        HF.mooney_rivlin(0.3, 0.2, 0.0)
    with pytest.raises(ValueError):
        HF.mooney_rivlin(-0.3, 0.2, 0.1)
    with pytest.raises(ValueError):
        HF.polynomial({(2, 2): 1.0}, (0.1,))
    m = HF.yeoh(0.5, -0.05, 0.02, 0.4)
    assert m.mu0 == 1.0 and m.K0 == 5.0 and m.inv_D[1] == 0.0 and m.inv_D[2] == 0.0


# --------------------------------------------------------------------------
# 5. non-homogeneous, near-incompressible finite deformation: Cook E04
# --------------------------------------------------------------------------

@pytest.mark.slow
def test_cook_e04_abaqus_neo_hookean_within_one_percent():
    """E04 Cook (nu = 0.49999, TL Neo-Hookean): u_y within 1 % of the Coreform
    reference 6.98946 (the G5 gate of ``fbar.py``, S40), here with the Abaqus
    NEO HOOKE law C10 = mu/2, D1 = 2/K.

    Same setup as ``test_fbar_mean_dilatation._cook_case``: p=2 C^1,
    16x16x1, P0 per cell, u_z fixed everywhere, left face clamped, consistent
    load int N ds of the 6.25 N/mm^2 shear traction on the right face. The
    NH variant differs from the reference law only in the volumetric part,
    an O(mu/K) = 3e-5 effect at nu = 0.49999, far below the 1 % gate.
    About 40 s (CSR + Pardiso).
    """
    from wht_iga.reference.cook import HR_MM, T_MM, U_Y_REF_MM

    E, nu, n, p = 240.565, 0.49999, 16, 2
    mu, K = E / (2.0 * (1.0 + nu)), E / (3.0 * (1.0 - 2.0 * nu))
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, 1), degree=p)
    model = HF.build_model(box, p, HF.neo_hookean(mu / 2.0, 2.0 / K))
    assert model["n_dof"] == 3 * 18 * 18 * 3
    cp = np.asarray(box.control_points())
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    kv, kw = (np.asarray(k, dtype=float) for k in box.knots[1:])
    a = np.array([(kv[j + p + 1] - kv[j]) / (p + 1) for j in range(nv_)])
    c = np.array([(kw[k + p + 1] - kw[k]) / (p + 1) for k in range(nw_)])
    f = np.zeros(model["n_dof"])
    fixed = []
    for I, x in enumerate(cp):
        fixed.append(3 * I + 2)
        if abs(x[0]) < 1e-9:
            fixed += [3 * I, 3 * I + 1]
    for k in range(nw_):
        for j in range(nv_):
            f[3 * ((nu_ - 1) + nu_ * (j + nv_ * k)) + 1] = 6.25 * HR_MM * T_MM * a[j] * c[k]
    assert abs(f.sum() - 6.25 * HR_MM * T_MM) < 1e-9 * 6.25 * HR_MM * T_MM
    u = HF.solve_static(model, f, np.asarray(fixed, dtype=int), n_steps=10)
    tips = [(nu_ - 1) + nu_ * ((nv_ - 1) + nv_ * k) for k in range(nw_)]
    uy = float(np.mean([u[3 * I + 1] for I in tips]))
    assert abs(uy - U_Y_REF_MM) / U_Y_REF_MM < 0.01, uy
