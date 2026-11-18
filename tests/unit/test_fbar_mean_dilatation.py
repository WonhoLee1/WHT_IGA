# -*- coding: utf-8 -*-
"""Mean-dilatation F-bar Neo-Hookean: patch, FD tangent, sparse vs dense,
reduced continuity (C^0 vs independent Q9/P0), uniaxial, Cook gate."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses import fbar as FB
from wht_iga.reference.box import BoxSpline


def _unit_model(p: int = 1):
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(1, 1, 1), degree=p)
    return FB.build_model(box, p, E=3.0, nu=0.3), box


def test_fbar_constant_strain_stress():
    """Affine displacement => constant F, Jbar == J, analytic NH stress."""
    model, box = _unit_model()
    cp = np.asarray(box.control_points())
    A = np.array([[0.02, -0.01, 0.005], [0.01, 0.03, -0.008], [0.0, 0.012, -0.02]])
    b = np.array([0.001, -0.002, 0.003])
    u = (cp @ A.T + b).reshape(-1)
    F = np.eye(3) + A
    J = float(np.linalg.det(F))
    mu, lam = model["mu"], model["lam"]
    B = F @ F.T
    # S19 closed form: (mu/J)(B - I) + (lam*lnJ/J) I  (== our W with Jbar=J)
    sig_ref = mu / J * (B - np.eye(3)) + (lam * np.log(J)) / J * np.eye(3)
    got = FB.cauchy_stress(model, u)[0]
    assert got.shape == (8, 6)
    M = np.array([[[r[0], r[3], r[5]],
                   [r[3], r[1], r[4]],
                   [r[5], r[4], r[2]]] for r in got])
    assert np.max(np.abs(M - sig_ref)) / np.max(np.abs(sig_ref)) < 1e-12


@pytest.mark.parametrize("p,n_cells,cont", [
    (1, (3, 2, 1), None), (2, (3, 2, 1), None), (2, (1, 2, 3), None),
    (2, (3, 2, 1), 0), (2, (1, 2, 3), 0), (3, (2, 1, 2), 1), (3, (2, 1, 2), 0),
])
def test_fbar_multicell_volume_and_patch(p, n_cells, cont):
    """Multi-cell, unequal cells per axis, tapered Cook map, any continuity.

    Single-cell tests cannot see a per-axis knot mix-up (all axes share one
    knot vector). Measured 2026-11-19: with the z basis built from the x knots,
    the affine stress patch still passes (isoparametric reproduction) but the
    quadrature volume is wrong by +200 % (p=1) / +30 % (p=2), and (1,2,3)
    raises ZeroDivisionError. So the exact trapezoid volume is the gate.
    ``cont`` < p-1 (repeated interior knots; Coreform E04 uses C^(p-2)) checks
    the cell -> span / first-control-point map of ``build_model``.
    """
    from wht_iga.reference.cook import CookSpline, HL_MM, HR_MM, L_MM, T_MM

    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=p,
                     continuity=cont)
    model = FB.build_model(box, p, E=3.0, nu=0.3)
    vol = sum(float(c["W"].sum()) for c in model["cells"])
    vol_ref = 0.5 * (HL_MM + HR_MM) * L_MM * T_MM
    assert abs(vol - vol_ref) / vol_ref < 1e-12, (vol, vol_ref)

    cp = np.asarray(box.control_points())
    A = np.array([[0.02, -0.01, 0.005], [0.01, 0.03, -0.008], [0.0, 0.012, -0.02]])
    u = (cp @ A.T + np.array([0.001, -0.002, 0.003])).reshape(-1)
    F = np.eye(3) + A
    J = float(np.linalg.det(F))
    mu, lam = model["mu"], model["lam"]
    ref = mu / J * (F @ F.T - np.eye(3)) + lam * np.log(J) / J * np.eye(3)
    refv = np.array([ref[0, 0], ref[1, 1], ref[2, 2], ref[0, 1], ref[1, 2], ref[0, 2]])
    for S in FB.cauchy_stress(model, u):
        assert np.max(np.abs(S - refv)) / np.max(np.abs(refv)) < 1e-12


def test_fbar_tangent_matches_finite_differences():
    """Analytic consistent tangent vs central FD (catches algebra slips)."""
    model, box = _unit_model()
    cp = np.asarray(box.control_points())
    rng = np.random.default_rng(7)
    u = 0.01 * rng.standard_normal(cp.size)
    R0, K = FB.residual_tangent(model, u)
    h = 1e-7
    Kfd = np.zeros_like(K)
    for j in range(u.size):
        up, um = u.copy(), u.copy()
        up[j] += h
        um[j] -= h
        Rp, _ = FB.residual_tangent(model, up)
        Rm, _ = FB.residual_tangent(model, um)
        Kfd[:, j] = (Rp - Rm) / (2 * h)
    denom = max(1e-30, float(np.max(np.abs(Kfd))))
    assert float(np.max(np.abs(K - Kfd))) / denom < 1e-5


def test_fbar_tangent_fd_multicell_c0_sparse():
    """Sparse CSR tangent vs central FD on a multi-cell p=2 C^0 box (nu=0.45).

    The single-cell FD test above cannot see cross-cell scatter errors or the
    shared C^0 interface control points.
    """
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 1, 2),
                    degree=2, continuity=0)
    assert tuple(box.n_cp) == (5, 3, 5)
    model = FB.build_model(box, 2, E=3.0, nu=0.45)
    rng = np.random.default_rng(11)
    u = 0.01 * rng.standard_normal(model["n_dof"])
    _, K = FB.residual_tangent_sparse(model, u)
    K = K.toarray()
    h = 1e-7
    Kfd = np.zeros_like(K)
    for j in range(u.size):
        up, um = u.copy(), u.copy()
        up[j] += h
        um[j] -= h
        Kfd[:, j] = (FB.residual_tangent_sparse(model, up)[0]
                     - FB.residual_tangent_sparse(model, um)[0]) / (2 * h)
    assert float(np.max(np.abs(K - Kfd))) / float(np.max(np.abs(Kfd))) < 1e-5


@pytest.mark.parametrize("p,n_cells,cont", [(2, (3, 2, 2), None), (2, (2, 3, 1), 0), (3, (2, 1, 2), 1)])
def test_fbar_sparse_assembly_matches_dense_oracle(p, n_cells, cont):
    """CSR scatter vs the dense ``np.ix_`` oracle at a random state.

    Same cell kernel and same cell order, so any index-map error is O(1);
    1e-13 only allows for summation round-off.
    """
    from wht_iga.reference.cook import CookSpline

    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=p,
                     continuity=cont)
    model = FB.build_model(box, p, E=240.565, nu=0.49)
    u = 0.05 * np.random.default_rng(5).standard_normal(model["n_dof"])
    Rd, Kd = FB.residual_tangent(model, u)
    Rs, Ks = FB.residual_tangent_sparse(model, u)
    assert Ks.shape == Kd.shape and Ks.has_sorted_indices
    coupled = np.zeros(Kd.shape, dtype=bool)
    for cell in model["cells"]:
        d = (3 * cell["nodes"][:, None] + np.arange(3)).ravel()
        coupled[np.ix_(d, d)] = True
    assert Ks.nnz == int(coupled.sum())  # stored = structural couplings, no padding
    assert float(np.max(np.abs(Ks.toarray() - Kd))) <= 1e-13 * float(np.max(np.abs(Kd)))
    assert float(np.max(np.abs(Rs - Rd))) <= 1e-13 * float(np.max(np.abs(Rd)))


def _cook_case(n: int, p: int, cont, E: float, nu: float, traction: float):
    """Cook model + consistent right-face load (int N ds), u_z fixed, left clamped."""
    from wht_iga.reference.cook import CookSpline, HR_MM, L_MM, T_MM

    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, 1), degree=p,
                     continuity=cont)
    model = FB.build_model(box, p, E=E, nu=nu)
    cp = np.asarray(box.control_points())
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    kv, kw = (np.asarray(k, dtype=float) for k in box.knots[1:])
    a = np.array([(kv[j + p + 1] - kv[j]) / (p + 1) for j in range(nv_)])
    c = np.array([(kw[k + p + 1] - kw[k]) / (p + 1) for k in range(nw_)])
    f = np.zeros(model["n_dof"])
    fixed = []
    for I, (x, y, z) in enumerate(cp):
        fixed.append(3 * I + 2)
        if abs(x) < 1e-9:
            fixed += [3 * I, 3 * I + 1]
    for k in range(nw_):
        for j in range(nv_):
            I = (nu_ - 1) + nu_ * (j + nv_ * k)
            assert abs(cp[I, 0] - L_MM) < 1e-9
            f[3 * I + 1] = traction * HR_MM * T_MM * a[j] * c[k]
    assert abs(f.sum() - traction * HR_MM * T_MM) < 1e-9 * traction * HR_MM * T_MM
    tips = [(nu_ - 1) + nu_ * ((nv_ - 1) + nv_ * k) for k in range(nw_)]
    return model, f, np.asarray(fixed, dtype=int), tips


def test_fbar_sparse_solve_matches_dense_solve():
    """Nonlinear Cook p=2 C^0 2x2 (nu=0.49, full 6.25 load): sparse Newton
    (auto = pardiso if installed, superlu) vs the dense oracle path.

    Same Newton equations, only the linear solver differs; dtol=1e-12 makes
    the converged states agree far below the 1e-8 bound.
    """
    model, f, fixed, tips = _cook_case(2, 2, 0, 240.565, 0.49, 6.25)
    u = {k: FB.solve_static(model, f, fixed, n_steps=4, dtol=1e-12, linear_solver=k)
         for k in ("dense", "auto", "superlu")}
    scale = float(np.max(np.abs(u["dense"])))
    assert u["dense"][3 * tips[0] + 1] > 1.0  # genuinely nonlinear (mm-scale tip)
    for k in ("auto", "superlu"):
        assert float(np.max(np.abs(u[k] - u["dense"]))) <= 1e-8 * scale, k


def _q9_mean_dilatation_cook_tip_uy(n: int, E: float, nu: float, traction: float) -> float:
    """Independent reference: 2D plane-strain Lagrange Q9 / P0 mean dilatation.

    Linear energy per element (unit thickness), the small-strain limit of the
    F-bar law in plane strain: mu|dev_2D eps|^2 + (lam+mu)/2 * A_e * mean(div u)^2.
    Q9 nodes on the Cook bilinear map; consistent edge load (1/6, 2/3, 1/6) h_e.
    """
    mu, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    m = 2 * n + 1
    s1 = np.linspace(0.0, 1.0, m)
    X = np.array([(48.0 * s1[a], 44.0 * s1[a] + 44.0 * s1[b] - 28.0 * s1[a] * s1[b])
                  for b in range(m) for a in range(m)])
    D = mu * np.array([[1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    gx, gw = np.polynomial.legendre.leggauss(3)

    def lag(s):
        return np.array([0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1)])

    def dlag(s):
        return np.array([s - 0.5, -2 * s, s + 0.5])

    K = np.zeros((2 * m * m, 2 * m * m))
    for j in range(n):
        for i in range(n):
            nodes = [(2 * i + a) + m * (2 * j + b) for b in range(3) for a in range(3)]
            dofs = np.array([[2 * I, 2 * I + 1] for I in nodes]).ravel()
            Ke, dv, area = np.zeros((18, 18)), np.zeros(18), 0.0
            for s, ws in zip(gx, gw):
                for r, wr in zip(gx, gw):
                    dN = np.array([np.outer(lag(s), dlag(r)).ravel(),
                                   np.outer(dlag(s), lag(r)).ravel()])
                    Jm = dN @ X[nodes]  # Jm[i, k] = dX_k / dxi_i
                    dNx = np.linalg.solve(Jm, dN)
                    B = np.zeros((3, 18))
                    B[0, 0::2], B[1, 1::2] = dNx[0], dNx[1]
                    B[2, 0::2], B[2, 1::2] = dNx[1], dNx[0]
                    w = ws * wr * np.linalg.det(Jm)
                    Ke += w * B.T @ D @ B
                    dv += w * (B[0] + B[1])
                    area += w
            K[np.ix_(dofs, dofs)] += Ke + (lam + mu) / area * np.outer(dv, dv)
    f = np.zeros(2 * m * m)
    for j in range(n):
        for b, wt in zip(range(3), (1 / 6, 2 / 3, 1 / 6)):
            f[2 * ((m - 1) + m * (2 * j + b)) + 1] += traction * (16.0 / n) * wt
    free = np.ones(2 * m * m, dtype=bool)
    free[[2 * m * b + d for b in range(m) for d in (0, 1)]] = False
    u = np.zeros(2 * m * m)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free])
    tip = (m - 1) + m * (m - 1)
    assert np.allclose(X[tip], (48.0, 60.0))
    return float(u[2 * tip + 1])


@pytest.mark.parametrize("n,nu", [(2, 0.3), (3, 0.49)])
def test_fbar_c0_linear_tangent_matches_q9_mean_dilatation(n, nu):
    """p=2 C^0 3D F-bar tangent at u=0, Cook (u_z fixed, one cell thick),
    vs an independent 2D Lagrange Q9/P0 mean-dilatation code.

    Bernstein C^0 p=2 and Lagrange Q9 span the same space and the solution is
    z-independent (separable cell mean, consistent z load), so the tip u_y must
    agree to solver round-off (cond ~ lam/mu * h^-2 << 1e7 -> 1e-9 bound).
    The volume/patch tests cannot see a wrong but self-consistent cell map.
    """
    import scipy.sparse.linalg as spla

    model, f, fixed, tips = _cook_case(n, 2, 0, 240.565, nu, 6.25)
    _, K = FB.residual_tangent_sparse(model, np.zeros(model["n_dof"]))
    free = np.setdiff1d(np.arange(model["n_dof"]), fixed)
    u = np.zeros(model["n_dof"])
    u[free] = spla.spsolve(K[free][:, free].tocsc(), f[free])
    uy = np.array([u[3 * I + 1] for I in tips])
    ref = _q9_mean_dilatation_cook_tip_uy(n, 240.565, nu, 6.25)
    assert np.max(np.abs(uy - ref)) <= 1e-9 * abs(ref), (uy, ref)


def test_fbar_uniaxial_matches_s19_analytic():
    """Single-cell uniaxial NH vs S19 closed form (Jbar == J here)."""
    from tests.kratos_spikes.test_s19_constitutive_models import (
        _analytic_uniaxial_neohookean,
    )

    E, nu, lam1 = 3.0, 0.3, 1.2
    _, sig_ref, _ = _analytic_uniaxial_neohookean(E, nu, lam1)
    model, box = _unit_model()
    cp = np.asarray(box.control_points())
    fixed, tip = [], []
    for I, (x, y, z) in enumerate(cp):
        if abs(x) < 1e-12:
            fixed.append(3 * I)  # clamp x only; lateral free to contract
        if abs(x - 1.0) < 1e-12:
            tip.append(I)
    assert abs(np.asarray(cp[tip])[:, 0].mean() - 1.0) < 1e-12
    prescribed = {3 * I: lam1 - 1.0 for I in tip}
    u = FB.solve_static(model, np.zeros(model["n_dof"]),
                        np.asarray(fixed, dtype=int), n_steps=4, prescribed=prescribed)
    sig = FB.cauchy_stress(model, u)[0]
    sig11 = float(sig[:, 0].mean())
    assert abs(sig11 - sig_ref) / abs(sig_ref) < 1e-9, (sig11, sig_ref)


@pytest.mark.slow
def test_fbar_cook_e04_gate():
    """E04 Cook (nu=0.49999, TL Neo-Hookean): u_y within 1 % of 6.98946 (G5 gate).

    p=2, 16x16x1 cells, P0 per cell, u_z fixed everywhere (e04.json), consistent
    load for the uniform 6.25 N/mm^2 shear traction on the right face
    (int N_i ds of each B-spline, exact). Runtime about 35 s (CSR + Pardiso;
    about 3 min with the former dense Newton, same u_y 7.03806, S40).

    Measured 2026-11-19, consistent load: p=2 4/8/12/16 = +5.55 / +2.25 /
    +1.21 / +0.70 % (converges from above, observed order about 1.1, corner
    singularity). p=1 4/8/16 (tributary load = consistent for p=1) = -11.0 /
    -3.76 / -1.29 % (from below, Richardson limit about 6.987). Stock Kratos
    TL locks at -69 %. A lumped Greville load gives p=2 8x8 +0.52 % only by
    under-loading the corner node; do not use it for p >= 2.
    """
    from wht_iga.reference.cook import CookSpline, HR_MM, L_MM, T_MM, U_Y_REF_MM

    p = 2
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                     n_cells=(16, 16, 1), degree=p)
    model = FB.build_model(box, p, E=240.565, nu=0.49999)
    cp = np.asarray(box.control_points())
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    kv, kw = (np.asarray(k, dtype=float) for k in box.knots[1:])
    a = np.array([(kv[j + p + 1] - kv[j]) / (p + 1) for j in range(nv_)])  # int N_j d(eta)
    c = np.array([(kw[k + p + 1] - kw[k]) / (p + 1) for k in range(nw_)])  # int N_k d(zeta)
    assert abs(a.sum() - 1.0) < 1e-14 and abs(c.sum() - 1.0) < 1e-14
    traction = 6.25
    f = np.zeros(model["n_dof"])
    fixed = []
    for I, (x, y, z) in enumerate(cp):
        fixed.append(3 * I + 2)  # e04.json: u_z fixed throughout the body
        if abs(x) < 1e-9:
            fixed += [3 * I, 3 * I + 1]
    for k in range(nw_):
        for j in range(nv_):
            I = (nu_ - 1) + nu_ * (j + nv_ * k)
            assert abs(cp[I, 0] - L_MM) < 1e-9
            f[3 * I + 1] = traction * HR_MM * T_MM * a[j] * c[k]
    assert abs(f.sum() - traction * HR_MM * T_MM) < 1e-9  # 191.516 N
    u = FB.solve_static(model, f, np.asarray(fixed, dtype=int), n_steps=10)
    tips = [(nu_ - 1) + nu_ * ((nv_ - 1) + nv_ * k) for k in range(nw_)]
    uy = float(np.mean([u[3 * I + 1] for I in tips]))
    assert abs(uy - U_Y_REF_MM) / U_Y_REF_MM < 0.01, uy
