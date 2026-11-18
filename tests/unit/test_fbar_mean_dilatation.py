# -*- coding: utf-8 -*-
"""Mean-dilatation F-bar Neo-Hookean: patch, FD tangent, uniaxial, Cook gate."""

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


@pytest.mark.parametrize("p,n_cells", [(1, (3, 2, 1)), (2, (3, 2, 1)), (2, (1, 2, 3))])
def test_fbar_multicell_volume_and_patch(p, n_cells):
    """Multi-cell, unequal cells per axis, tapered Cook map.

    Single-cell tests cannot see a per-axis knot mix-up (all axes share one
    knot vector). Measured 2026-11-19: with the z basis built from the x knots,
    the affine stress patch still passes (isoparametric reproduction) but the
    quadrature volume is wrong by +200 % (p=1) / +30 % (p=2), and (1,2,3)
    raises ZeroDivisionError. So the exact trapezoid volume is the gate.
    """
    from wht_iga.reference.cook import CookSpline, HL_MM, HR_MM, L_MM, T_MM

    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=p)
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
    (int N_i ds of each B-spline, exact). Runtime about 3 min (dense Newton).

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
