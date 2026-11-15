"""P4-C: thermal math without Kratos -- conduction, capacity, film/radiation faces, time integration, eigenstrain.

Oracles are closed forms or exact semi-discrete solutions:

* linear profiles and ``k V`` / ``rho c V`` / ``rho c int x^2`` energies (B-splines reproduce linears);
* 1D film ``T(x)``: ``g = -h (T_A - T_inf) / (k + h L)``; 1D radiation: scalar root of
  ``k (T_A - T_L) / L = eps sigma (T_L**4 - T_inf**4)``;
* exact semi-discrete solution ``T(t) = sum exp(-lambda t) v v^T C T0`` from the generalized
  eigenproblem ``K v = lambda C v`` (isolates the time-integration error);
* the modal series of the 1D diffusion problem with ``T_B(t) = 100 sin(pi t / 40)`` (NAFEMS T3 data);
* free thermal expansion is stress free.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.linalg as la
import scipy.sparse as sp
from scipy.optimize import brentq

from wht_iga.analyses.thermal import (
    STEFAN_BOLTZMANN,
    assemble_capacity,
    assemble_conduction,
    assemble_stiffness_tdep,
    eigenstrain_load,
    face_rule,
    plane_stress_as_plane_strain,
    recover_thermoelastic_stress_tdep,
    solve_thermoelastic_tdep,
    solve_with_dirichlet,
    steady_thermal,
    temperature_dependent_eigenstrain_load,
    transient_thermal,
)
from wht_iga.reference.assemble_elastic import assemble_stiffness, isotropic_D
from wht_iga.reference.box import BoxSpline

L_BAR, SECTION = 0.1, 0.02


def _bar(p: int = 1, n: int = 10) -> BoxSpline:
    return BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L_BAR, SECTION, SECTION), n_cells=(n, 1, 1), degree=p)


def _x_cps(box: BoxSpline) -> np.ndarray:
    return box.control_points()[:, 0]


# ---------------------------------------------------------------------------
# matrices and rules
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [1, 2, 3])
def test_conduction_matrix_energy_nullspace_symmetry(p: int) -> None:
    """K 1 = 0, symmetric, and for T = x: T^T K T = k V (Greville points reproduce linears)."""
    box = BoxSpline(lo=(0, 0, 0), hi=(0.6, 1.0, 0.05), n_cells=(3, 4, 2), degree=p)
    k = 52.0
    K = assemble_conduction(box, k)
    ones = np.ones(box.n_ctrl)
    assert np.abs(K @ ones).max() < 1e-11 * abs(K).max()
    assert abs(K - K.T).max() < 1e-13 * abs(K).max()
    x = box.control_points()[:, 0]
    assert float(x @ (K @ x)) == pytest.approx(k * 0.6 * 1.0 * 0.05, rel=1e-12)


@pytest.mark.parametrize("p", [1, 2, 3])
def test_capacity_matrix_closed_forms(p: int) -> None:
    """1^T C 1 = rho c V and x^T C x = rho c int x^2 dV."""
    box = BoxSpline(lo=(0, 0, 0), hi=(0.6, 1.0, 0.05), n_cells=(3, 4, 2), degree=p)
    rho_c = 7200.0 * 440.5
    C = assemble_capacity(box, rho_c)
    ones = np.ones(box.n_ctrl)
    x = box.control_points()[:, 0]
    vol = 0.6 * 1.0 * 0.05
    assert float(ones @ (C @ ones)) == pytest.approx(rho_c * vol, rel=1e-13)
    assert float(x @ (C @ x)) == pytest.approx(rho_c * (0.6**3 / 3.0) * 1.0 * 0.05, rel=1e-12)


@pytest.mark.parametrize("axis, side", [(0, 0), (0, 1), (1, 1), (2, 0)])
@pytest.mark.parametrize("p", [1, 2])
def test_face_rule_area_partition_of_unity_and_location(p: int, axis: int, side: int) -> None:
    """Weights sum to the face area, basis rows sum to one, points lie on the face; centroid exact."""
    box = BoxSpline(lo=(0.1, 0.2, 0.3), hi=(0.7, 1.0, 0.5), n_cells=(3, 2, 2), degree=p)
    rule = face_rule(box, axis, side)
    ext = np.asarray(box.hi) - np.asarray(box.lo)
    tang = [d for d in range(3) if d != axis]
    assert rule.weights.sum() == pytest.approx(ext[tang[0]] * ext[tang[1]], rel=1e-13)
    np.testing.assert_allclose(np.asarray(rule.basis.sum(axis=1)).ravel(), 1.0, atol=1e-13)
    assert np.allclose(rule.points[:, axis], box.hi[axis] if side else box.lo[axis])
    centroid = (rule.weights[:, None] * rule.points).sum(axis=0) / rule.weights.sum()
    for d in tang:
        assert centroid[d] == pytest.approx(0.5 * (box.lo[d] + box.hi[d]), rel=1e-13)
    # interpolating a linear field through the face basis reproduces it at the points
    t = box.control_points() @ np.array([2.0, -1.0, 0.5])
    np.testing.assert_allclose(rule.basis @ t, rule.points @ np.array([2.0, -1.0, 0.5]), atol=1e-12)


# ---------------------------------------------------------------------------
# steady problems
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [1, 2])
def test_steady_film_matches_closed_form(p: int) -> None:
    """Fixed T_A at x = 0, film (h, T_inf) at x = L: linear profile with g = -h (T_A - T_inf) / (k + h L)."""
    box = _bar(p)
    k, h, t_a, t_inf = 35.0, 750.0, 100.0, 5.0
    K = assemble_conduction(box, k)
    left = np.asarray(box.boundary_control_points(0, 0))
    t, _ = steady_thermal(K, box.n_ctrl, films=[(face_rule(box, 0, 1), h, t_inf)], fixed=left, fixed_values=np.full(left.size, t_a))
    g = -h * (t_a - t_inf) / (k + h * L_BAR)
    np.testing.assert_allclose(t, t_a + g * _x_cps(box), rtol=1e-11, atol=1e-9)


@pytest.mark.parametrize("p", [1, 2])
def test_steady_radiation_matches_scalar_newton(p: int) -> None:
    """NAFEMS T2 data: T_A = 1000 K, eps = 0.98, T_inf = 300 K, k = 55.6, L = 0.1."""
    box = _bar(p)
    k, eps, t_a, t_inf = 55.6, 0.98, 1000.0, 300.0
    K = assemble_conduction(box, k)
    left = np.asarray(box.boundary_control_points(0, 0))
    t, its = steady_thermal(K, box.n_ctrl, radiations=[(face_rule(box, 0, 1), eps, t_inf)], fixed=left,
                            fixed_values=np.full(left.size, t_a))
    root = brentq(lambda tl: k * (t_a - tl) / L_BAR - eps * STEFAN_BOLTZMANN * (tl**4 - t_inf**4), 300.0, 1000.0, xtol=1e-13)
    right = np.asarray(box.boundary_control_points(0, 1))
    assert t[right] == pytest.approx(root, abs=1e-8)
    assert its <= 12
    np.testing.assert_allclose(t, t_a + (root - t_a) * _x_cps(box) / L_BAR, atol=1e-7)  # linear conduction between


def test_steady_two_faces_and_source_balance() -> None:
    """Insulated bar with a nodal source and one film face: delivered heat equals film loss."""
    box = _bar(2, n=6)
    K = assemble_conduction(box, 20.0)
    rule = face_rule(box, 0, 1)
    rho_c = np.ones(box.n_ctrl)
    source = assemble_capacity(box, 1.0) @ rho_c * 1000.0  # uniform volumetric source 1000 W/m^3
    t, _ = steady_thermal(K, box.n_ctrl, films=[(rule, 50.0, 10.0)], source=source, fixed=np.zeros(0, dtype=int), fixed_values=np.zeros(0))
    q_in = 1000.0 * (L_BAR * SECTION * SECTION)
    q_out = float((rule.weights * 50.0 * (rule.basis @ t - 10.0)).sum())
    assert q_out == pytest.approx(q_in, rel=1e-9)


# ---------------------------------------------------------------------------
# transient
# ---------------------------------------------------------------------------
def _exact_semi_discrete(C: sp.spmatrix, K: sp.spmatrix, free: np.ndarray, t0: np.ndarray, t: float) -> np.ndarray:
    lam, v = la.eigh(K.toarray()[np.ix_(free, free)], C.toarray()[np.ix_(free, free)])
    out = np.zeros(K.shape[0])
    out[free] = v @ (np.exp(-lam * t) * (v.T @ (C.toarray()[np.ix_(free, free)] @ t0[free])))
    return out


@pytest.mark.parametrize("scheme, order", [("BE", 1.0), ("CN", 2.0), ("BDF2", 2.0)])
def test_time_integrators_converge_with_their_order(scheme: str, order: float) -> None:
    """Error against the exact semi-discrete solution halves/quarters when dt is halved."""
    box = _bar(2, n=8)
    k, rho_c = 35.0, 7200.0 * 440.5
    K, C = assemble_conduction(box, k), assemble_capacity(box, rho_c)
    ends = np.concatenate([box.boundary_control_points(0, 0), box.boundary_control_points(0, 1)])
    free = np.setdiff1d(np.arange(box.n_ctrl), ends)
    t0 = np.sin(np.pi * _x_cps(box) / L_BAR)
    t0[ends] = 0.0
    t_end = 12.0
    ref = _exact_semi_discrete(C, K, free, t0, t_end)
    errs = []
    for n in (60, 120, 240):
        _times, hist = transient_thermal(C, K, t0, t_end / n, n, ends, lambda t: np.zeros(ends.size), scheme)
        errs.append(np.abs(hist[-1] - ref).max())
    assert errs[0] > errs[1] > errs[2]
    for a, b in zip(errs[:-1], errs[1:], strict=True):
        assert math.log2(a / b) == pytest.approx(order, abs=0.25)


def _t3_series(x: float, t: float, length: float = 0.1, k: float = 35.0, rho: float = 7200.0, c: float = 440.5,
               period_half: float = 40.0, amp: float = 100.0, n_modes: int = 400) -> float:
    """Modal series of u_t = kappa u_xx, u(0) = 0, u(L, t) = amp sin(pi t / 40), u(x, 0) = 0."""
    kappa = k / (rho * c)
    om = math.pi / period_half
    total = amp * math.sin(om * t) * x / length
    for n in range(1, n_modes + 1):
        lam = kappa * (n * math.pi / length) ** 2
        cn = 2.0 / (n * math.pi) * (-1.0) ** (n + 1)
        integral = (lam * math.cos(om * t) + om * math.sin(om * t) - lam * math.exp(-lam * t)) / (lam**2 + om**2)
        total += -cn * amp * om * integral * math.sin(n * math.pi * x / length)
    return total


def test_time_dependent_dirichlet_matches_modal_series() -> None:
    """V12-T3 data on a fine p = 2 mesh: BDF2 reproduces the modal series to 1e-4 K at t = 32 s."""
    box = _bar(2, n=40)
    k, rho_c = 35.0, 7200.0 * 440.5
    K, C = assemble_conduction(box, k), assemble_capacity(box, rho_c)
    left = np.asarray(box.boundary_control_points(0, 0))
    right = np.asarray(box.boundary_control_points(0, 1))
    fixed = np.concatenate([left, right])
    fn = lambda t: np.concatenate([np.zeros(left.size), np.full(right.size, 100.0 * math.sin(math.pi * t / 40.0))])  # noqa: E731
    times, hist = transient_thermal(C, K, np.zeros(box.n_ctrl), 0.05, 640, fixed, fn, "BDF2")
    assert times[-1] == pytest.approx(32.0)
    def temperature_at(history_row: np.ndarray, x: float) -> float:
        coeffs = np.column_stack([history_row, np.zeros(box.n_ctrl), np.zeros(box.n_ctrl)])
        return float(box.evaluate(coeffs, np.array([[x, 0.5 * SECTION, 0.5 * SECTION]]))[0, 0])

    assert temperature_at(hist[-1], 0.08) == pytest.approx(_t3_series(0.08, 32.0), abs=1e-4)
    assert temperature_at(hist[0], 0.08) == 0.0


# ---------------------------------------------------------------------------
# thermal strain
# ---------------------------------------------------------------------------
def test_plane_stress_conversion_relations() -> None:
    """(E_p, nu_p, e_p) satisfy E_s = E_p/(1 - nu_p^2), nu_s = nu_p/(1 - nu_p), e_s = (1 + nu_p) e_p."""
    e_s, nu_s, eps_s = 1.0e5, 0.3, 1.0e-3
    e_p, nu_p, e_pe = plane_stress_as_plane_strain(e_s, nu_s, eps_s)
    assert e_p / (1 - nu_p**2) == pytest.approx(e_s, rel=1e-14)
    assert nu_p / (1 - nu_p) == pytest.approx(nu_s, rel=1e-14)
    assert (1 + nu_p) * e_pe == pytest.approx(eps_s, rel=1e-14)


@pytest.mark.parametrize("p", [1, 2])
def test_free_thermal_expansion_is_stress_free(p: int) -> None:
    """Uniform eigenstrain on the whole box with rigid-body constraints only: u = eps* x, sigma = 0."""
    box = BoxSpline(lo=(0, 0, 0), hi=(1.0, 0.7, 0.5), n_cells=(3, 2, 2), degree=p)
    E, nu = 2.1, 0.27
    D = isotropic_D(E, nu)
    eps_star = np.array([1.0e-3, -2.0e-3, 0.5e-3, 0.8e-3, -0.3e-3, 0.6e-3])  # includes engineering shears
    cells = [(i, j, k) for k in range(2) for j in range(2) for i in range(3)]
    f = eigenstrain_load(box, D, eps_star, {}, cells)
    K = assemble_stiffness(box, E, nu)
    # remove the six rigid-body modes by fixing three corner control points suitably
    nu_cp, nv_cp, nw_cp = box.n_cp
    cp = lambda i, j, k: i + nu_cp * (j + nv_cp * k)  # noqa: E731
    fixed = np.array([3 * cp(0, 0, 0), 3 * cp(0, 0, 0) + 1, 3 * cp(0, 0, 0) + 2,
                      3 * cp(nu_cp - 1, 0, 0) + 1, 3 * cp(nu_cp - 1, 0, 0) + 2, 3 * cp(0, nv_cp - 1, 0) + 2])
    u = solve_with_dirichlet(K, f, fixed, np.zeros(fixed.size))
    # the exact field is affine: u = G x with G = [[exx, gxy/2, gxz/2], ...] up to a rotation fixed by the constraints
    pts = box.control_points()
    grad = np.zeros((3, 3))
    design = np.column_stack([pts, np.ones(len(pts))])
    sol, *_ = np.linalg.lstsq(design, u.reshape(-1, 3), rcond=None)
    grad = sol[:3].T  # du_i/dx_j
    strain = 0.5 * (grad + grad.T)
    eps_voigt = np.array([strain[0, 0], strain[1, 1], strain[2, 2], 2 * strain[0, 1], 2 * strain[1, 2], 2 * strain[0, 2]])
    np.testing.assert_allclose(eps_voigt, eps_star, atol=1e-12)
    assert np.abs(design @ sol - u.reshape(-1, 3)).max() < 1e-12  # the displacement really is affine


def _quad_x_bspline_coeffs(box: BoxSpline, c0: float, c1: float, c2: float) -> np.ndarray:
    """Exact B-spline control-point coefficients for T(x) = c0 + c1*x + c2*x^2 (degree >= 2)."""
    p = box.degree
    assert p >= 2
    kx = box.knots[0]
    nx_cp = box.n_cp[0]
    cx = np.empty(nx_cp, dtype=float)
    for i in range(nx_cp):
        t = kx[i + 1 : i + p + 1]
        xi1 = float(t.mean())
        xi2 = (2.0 / (p * (p - 1))) * sum(float(t[r] * t[s]) for r in range(p) for s in range(r + 1, p))
        cx[i] = c0 + c1 * xi1 + c2 * xi2
    return np.tile(cx, box.n_cp[1] * box.n_cp[2])


def test_tdep_reduces_to_constant_elastic_and_eigenstrain() -> None:
    """For uniform T, E, nu, alpha, temperature-dependent builders match constant assembly."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 0.5, 0.4), n_cells=(3, 2, 2), degree=2)
    E, nu, alpha, T_val, T_ref = 210.0, 0.3, 1.2e-5, 85.0, 20.0
    D = isotropic_D(E, nu)
    eps0 = alpha * (T_val - T_ref)
    eps_star = np.array([eps0, eps0, eps0, 0.0, 0.0, 0.0])
    cells = [(i, j, k) for k in range(2) for j in range(2) for i in range(3)]

    K_const = assemble_stiffness(box, E, nu)
    K_tdep = assemble_stiffness_tdep(box, np.full(box.n_ctrl, T_val), E, nu)
    assert abs(K_tdep - K_const).max() < 1e-13 * abs(K_const).max()

    f_const = eigenstrain_load(box, D, eps_star, {}, cells)
    f_tdep = temperature_dependent_eigenstrain_load(box, np.full(box.n_ctrl, T_val), E, nu, alpha, T_ref=T_ref)
    np.testing.assert_allclose(f_tdep, f_const, rtol=1e-13, atol=1e-13)


def test_tdep_quadratic_T_affine_E_alpha_free_expansion_bar() -> None:
    """Free-end 3D bar (nu=0, u_y=u_z=0) with quadratic T(x) and affine E(T), alpha(T): sigma_xx = 0."""
    L, w = 2.0, 0.2
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L, w, w), n_cells=(8, 1, 1), degree=4)
    T_ref = 25.0
    c0, c1, c2 = 25.0, 30.0, 20.0  # T(x) = 25 + 30 x + 20 x^2
    E0, E1 = 200.0, 0.5
    alpha0, alpha1 = 1.0e-5, 2.0e-8

    T_fn = lambda x, _y, _z: c0 + c1 * x + c2 * x**2  # noqa: E731
    E_fn = lambda T: E0 + E1 * (T - T_ref)  # noqa: E731
    alpha_fn = lambda T: alpha0 + alpha1 * (T - T_ref)  # noqa: E731
    t_cp = _quad_x_bspline_coeffs(box, c0, c1, c2)

    all_cps = np.arange(box.n_ctrl, dtype=np.int64)
    left_cps = np.asarray(box.boundary_control_points(0, 0), dtype=np.int64)
    fixed = np.unique(np.concatenate([3 * left_cps + 0, 3 * all_cps + 1, 3 * all_cps + 2]))

    res_cp = solve_thermoelastic_tdep(box, t_cp, E_fn, 0.0, alpha_fn, fixed, T_ref=T_ref)
    res_fn = solve_thermoelastic_tdep(box, T_fn, E_fn, 0.0, alpha_fn, fixed, T_ref=T_ref)
    np.testing.assert_allclose(res_cp.u, res_fn.u, rtol=1e-13, atol=1e-13)

    # Exact 1D free expansion u_x(x) = int_0^x (a0 + a1*(c1*s + c2*s^2))*(c1*s + c2*s^2) ds
    # = a0*(c1/2 x^2 + c2/3 x^3) + a1*(c1^2/3 x^3 + c1*c2/2 x^4 + c2^2/5 x^5)
    xs = np.linspace(0.0, L, 11)
    pts = np.column_stack([xs, np.full_like(xs, 0.5 * w), np.full_like(xs, 0.5 * w)])
    u_eval = box.evaluate(res_cp.u.reshape(-1, 3), pts)[:, 0]
    u_exact = alpha0 * (0.5 * c1 * xs**2 + (1.0 / 3.0) * c2 * xs**3) + alpha1 * (
        (c1**2 / 3.0) * xs**3 + 0.5 * c1 * c2 * xs**4 + (c2**2 / 5.0) * xs**5
    )
    # Exact u_x is degree 5, so p=4 B-spline on 8 cells achieves O(h^5) ~ 1e-7 relative error
    assert np.max(np.abs(u_eval - u_exact)) / np.max(np.abs(u_exact)) < 1e-7

    sig = recover_thermoelastic_stress_tdep(box, res_cp.u, t_cp, E_fn, 0.0, alpha_fn, pts, T_ref=T_ref)
    dT = c1 * xs + c2 * xs**2
    sig_lat_exact = -E_fn(c0 + dT) * alpha_fn(c0 + dT) * dT
    assert np.max(np.abs(sig[:, 0])) / np.max(np.abs(sig_lat_exact)) < 1e-6
    np.testing.assert_allclose(sig[:, 1], sig_lat_exact, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(sig[:, 2], sig_lat_exact, rtol=1e-6, atol=1e-8)


def test_tdep_quadratic_T_affine_E_alpha_clamped_bar_analytic_and_high_order_gauss() -> None:
    """Clamped 3D bar (nu=0, u_y=u_z=0): quadratic T(x) + affine E(T), alpha(T) vs arctan exact & ng=8."""
    L, w = 1.5, 0.2
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L, w, w), n_cells=(16, 1, 1), degree=4)
    dT, E0, E1, alpha0, alpha1 = 100.0, 100.0, 0.5, 1.0e-5, 2.5e-7
    c2 = dT / (L * L)  # T(x) = dT * (x / L)^2, T_ref = 0

    T_fn = lambda x, _y, _z: dT * (x / L) ** 2  # noqa: E731
    E_fn = lambda T: E0 + E1 * T  # noqa: E731
    alpha_fn = lambda T: alpha0 + alpha1 * T  # noqa: E731
    t_cp = _quad_x_bspline_coeffs(box, 0.0, 0.0, c2)

    all_cps = np.arange(box.n_ctrl, dtype=np.int64)
    ends = np.concatenate([box.boundary_control_points(0, 0), box.boundary_control_points(0, 1)])
    fixed = np.unique(np.concatenate([3 * ends + 0, 3 * all_cps + 1, 3 * all_cps + 2]))

    res = solve_thermoelastic_tdep(box, t_cp, E_fn, 0.0, alpha_fn, fixed, T_ref=0.0)
    res_hi = solve_thermoelastic_tdep(box, T_fn, E_fn, 0.0, alpha_fn, fixed, T_ref=0.0, n_gauss=8)
    np.testing.assert_allclose(res.u, res_hi.u, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(res.f, res_hi.f, rtol=1e-11, atol=1e-13)

    # Closed form: E(x) = E0 * (1 + beta * (x/L)^2), beta = E1 * dT / E0
    beta = E1 * dT / E0
    sqrt_b = math.sqrt(beta)
    sig0_exact = -(E0 * sqrt_b / math.atan(sqrt_b)) * (alpha0 * dT / 3.0 + alpha1 * dT**2 / 5.0)

    xs = np.linspace(0.05 * L, 0.95 * L, 19)
    pts = np.column_stack([xs, np.full_like(xs, 0.5 * w), np.full_like(xs, 0.5 * w)])
    sig = recover_thermoelastic_stress_tdep(box, res.u, t_cp, E_fn, 0.0, alpha_fn, pts, T_ref=0.0)
    assert abs(float(np.mean(sig[:, 0])) - sig0_exact) / abs(sig0_exact) < 1e-7
    assert float(np.max(np.abs(sig[:, 0] - sig0_exact)) / abs(sig0_exact)) < 5e-6

    xi = xs / L
    u_exact = sig0_exact * (L / (E0 * sqrt_b)) * np.arctan(sqrt_b * xi) + L * (
        (alpha0 * dT / 3.0) * xi**3 + (alpha1 * dT**2 / 5.0) * xi**5
    )
    u_eval = box.evaluate(res.u.reshape(-1, 3), pts)[:, 0]
    rel_u_err = float(np.max(np.abs(u_eval - u_exact)) / np.max(np.abs(u_exact)))
    assert rel_u_err < 1e-7


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
