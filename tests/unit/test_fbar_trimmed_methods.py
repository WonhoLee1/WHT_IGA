# -*- coding: utf-8 -*-
"""Trimmed F-bar (fbar_trimmed) wired into the methods API (Kratos-free).

``kind='full_iga_immersed'`` with ``locking_mitigation='fbar'`` (or ``'auto'``
at ``nu >= 0.49``) must solve the trimmed problem with
:mod:`wht_iga.analyses.fbar_trimmed`. References, independent of the routing
code under test:

* a direct ``build_trimmed_table`` + ``fbar_trimmed.build_model`` +
  ``fbar_trimmed.solve_static`` solve, with the face load from
  ``reference.assemble_elastic.assemble_face_traction`` (not the
  ``_assemble_fbar_face_traction`` the route uses) and the clamp built from
  ``BoxSpline.boundary_control_points``;
* the Kirsch closed form for the V07-like plate with a hole (formulas copied
  from ``tests/unit/test_fbar_trimmed.py``, which copied them from
  ``tests/kratos_spikes/test_p5_v07_kirsch_precursor.py``).

Routing to the linear (Kratos) path is detected by replacing
``wht_iga.solver.pipeline.run_linear`` with a sentinel that raises; no physics
is mocked.
"""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError

from wht_iga.analyses import fbar_trimmed as FT
from wht_iga.analyses.methods_options import FullIgaOption, parse_discretization_option
from wht_iga.analyses.methods_solve import (
    BoxMethodProblem,
    solve_fbar_trimmed_box,
    solve_with_discretization_method,
)
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.reference.assemble_elastic import assemble_face_traction
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.trimmed import build_trimmed_table

# V07 (tests/verification/data/v07.json)
E_V07 = 100000.0
NU_V07 = 0.49999
T_V07 = 10.0
A_V07 = 1.0
B_V07 = 4.0
T_THICK = 0.25


@pytest.fixture(autouse=True)
def _serial_uncached_trimming(monkeypatch):
    """Serial cut-cell fits (low RAM) and no on-disk trimmed-table cache for the route."""
    monkeypatch.setenv("WHT_GEO_WORKERS", "1")
    monkeypatch.setenv("WHT_GEO_CACHE", "off")


class _LinearPathReached(Exception):
    """Raised by the run_linear sentinel; deliberately not a ValueError."""


@pytest.fixture
def linear_sentinel(monkeypatch):
    import wht_iga.solver.pipeline as pipeline

    def _sentinel(*_args, **_kwargs):
        raise _LinearPathReached

    monkeypatch.setattr(pipeline, "run_linear", _sentinel)


# ---------------------------------------------------------------------------
# geometry and Kirsch reference (copied from tests/unit/test_fbar_trimmed.py)
# ---------------------------------------------------------------------------

def plate_hole_soup(a: float, b_ext: float, z_hi: float, n_seg: int = 64):
    """Closed outward-wound soup of {0 <= x, y <= b_ext, r >= a, 0 <= z <= z_hi}.

    Copy of ``tests/unit/test_fbar_trimmed.py::plate_hole_soup``: the walls
    x = 0, y = 0, z = 0 lie on the grid's lower planes, the solid extends
    beyond the grid's upper planes, and the polygon radius a*sqrt(ang/sin(ang))
    gives the quarter hole the exact area pi*a^2/4.
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
        tris.append((i, j, m + j))
        tris.append((i, m + j, m + i))
        tris.append((top_c, m + i, m + j))
        tris.append((bot_c, j, i))
    return nodes, np.asarray(tris, dtype=np.int64)


def kirsch_displacement(x, y, a=A_V07, E=E_V07, nu=NU_V07, T=T_V07):
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


def kirsch_sigma_xx(x, y, a=A_V07, T=T_V07):
    """Exact Kirsch sigma_xx."""
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    q2 = (a / r) ** 2
    q4 = q2 * q2
    return T * (1.0 - q2 * (1.5 * np.cos(2.0 * theta) + np.cos(4.0 * theta))
                + 1.5 * q4 * np.cos(4.0 * theta))


def _collocate_1d(knots, p, params, values):
    n_cp = len(knots) - p - 1
    A = np.zeros((n_cp, n_cp))
    for m, xi in enumerate(params):
        span = int(find_span(n_cp, p, float(xi), knots))
        A[m, span - p:span + 1] = ders_basis_funs(span, float(xi), p, 0, knots)[0]
    return np.linalg.solve(A, values)


def kirsch_dirichlet(box):
    """V07 precursor BCs: u_z = 0, u_x = 0 on x = 0, u_y = 0 on y = 0, exact
    Kirsch u (Greville interpolation) on x = b and y = b."""
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


def kirsch_errors(model, u):
    """(e_u, e_sigma_xx): relative L2 errors on the trimmed rule (as test_fbar_trimmed)."""
    X, W = model["qX"], model["qW"]
    uq = FT.displacement_at_qps(model, u)
    sq = FT.cauchy_stress_qps(model, u)
    ux, uy = kirsch_displacement(X[:, 0], X[:, 1])
    sxx = kirsch_sigma_xx(X[:, 0], X[:, 1])
    du2 = (uq[:, 0] - ux) ** 2 + (uq[:, 1] - uy) ** 2 + uq[:, 2] ** 2
    e_u = float(np.sqrt(np.sum(W * du2) / np.sum(W * (ux ** 2 + uy ** 2))))
    e_s = float(np.sqrt(np.sum(W * (sq[:, 0] - sxx) ** 2) / np.sum(W * sxx ** 2)))
    return e_u, e_s


# ---------------------------------------------------------------------------
# a small trimmed plate loaded through the methods facade
# ---------------------------------------------------------------------------

_PLATE = dict(a=1.0, b=2.0, t=0.5, n_cells=(4, 4, 2), b_ext=2.41, z_hi=0.66)


def _plate_problem(nu: float, **kw) -> BoxMethodProblem:
    nodes, tris = plate_hole_soup(_PLATE["a"], _PLATE["b_ext"], _PLATE["z_hi"])
    args = dict(
        lo=(0.0, 0.0, 0.0),
        hi=(_PLATE["b"], _PLATE["b"], _PLATE["t"]),
        E=1000.0,
        nu=nu,
        fixed_face=(0, 1),          # x = b, inside the solid
        load_face=(1, 1),           # y = b, inside the solid
        traction=(0.4, 2.0, 0.0),
        tess_soup=(nodes, tris),
    )
    args.update(kw)
    return BoxMethodProblem(**args)


def _direct_fbar_trimmed_model(problem: BoxMethodProblem, n_cells, q: int):
    """Independent build of the same trimmed model, load and clamp (no methods code)."""
    box = BoxSpline(lo=problem.lo, hi=problem.hi, n_cells=n_cells, degree=2)
    qp, out, _iss, _rec, _diag = build_trimmed_table(*problem.tess_soup, box, q=q, workers=1,
                                                     cache_dir="off")
    model = FT.build_model(box, qp, out, problem.E, problem.nu)
    f = assemble_face_traction(box, *problem.load_face, np.asarray(problem.traction))
    fixed = np.array([3 * int(c) + d for c in box.boundary_control_points(*problem.fixed_face)
                      for d in range(3)], dtype=np.int64)
    return box, model, f, fixed


@pytest.mark.parametrize("mitigation", ["fbar", "auto"])
def test_methods_route_equals_direct_fbar_trimmed(mitigation, linear_sentinel):
    """nu = 0.4999 trimmed plate (max |grad u| ~ 2 %): the facade solves the
    same discrete problem as a direct fbar_trimmed solve.

    Split in two exact checks: (1) the route's load equals the independent
    ``assemble_face_traction`` load and its clamp equals the
    ``boundary_control_points`` clamp, both to round-off (1e-12); (2) the
    route's solution equals ``fbar_trimmed.solve_static`` on the independently
    built model with that load (rel 1e-12, the same arithmetic). Feeding the
    reference load into (2) instead is not a 1e-12 comparison: its last-bit
    differences pass through the Newton iterations and cond(K) of the sliver
    basis (measured 4e-12 on the coefficients, 5e-13 on u_h at the points),
    so that variant is checked at 1e-9.
    """
    problem = _plate_problem(nu=0.4999)
    opt = FullIgaOption(kind="full_iga_immersed", degree=2, n_cells=_PLATE["n_cells"],
                        locking_mitigation=mitigation)
    assert parse_discretization_option(opt.model_dump_json()) == opt
    sol = solve_with_discretization_method(problem, opt, use_pardiso=False)
    assert sol.method == "full_iga_immersed" and sol.converged
    assert sol.native["locking_mitigation"] == "fbar"
    assert sol.native["quadrature_degree"] == opt.resolved_quadrature_degree == 4
    assert sol.native["eta_thresh"] == FT.DEFAULT_ETA_THRESH

    box, model, f_ref, fixed_ref = _direct_fbar_trimmed_model(problem, _PLATE["n_cells"], q=4)
    d = model["diag"]
    assert d["n_cut"] > 0 and model["n_groups"] < d["n_active"]  # trimmed, agglomerated
    assert d["n_inactive_ctrl"] > 0
    m = sol.native["model"]
    assert m["n_groups"] == model["n_groups"] and np.array_equal(m["qW"], model["qW"])
    # (1) load and clamp
    assert np.max(np.abs(sol.native["f_ext"] - f_ref)) <= 1e-12 * np.max(np.abs(f_ref))
    assert np.array_equal(sol.native["fixed"], np.sort(fixed_ref))
    # (2) solve
    u_ref = FT.solve_static(model, sol.native["f_ext"], fixed_ref, n_steps=1,
                            linear_solver="superlu")
    scale = np.max(np.abs(u_ref))
    uq_ref = FT.displacement_at_qps(model, u_ref)
    assert np.max(np.abs(uq_ref)) > 1e-3 * _PLATE["b"]  # a real deformation, not a zero solution
    err = np.max(np.abs(sol.native["u"] - u_ref)) / scale
    tip_ref = box.evaluate(u_ref.reshape(-1, 3), problem.tip_point().reshape(1, 3))[0]
    err_tip = np.max(np.abs(sol.tip_displacement - tip_ref)) / np.max(np.abs(tip_ref))
    u_f = FT.solve_static(model, f_ref, fixed_ref, n_steps=1, linear_solver="superlu")
    err_f = np.max(np.abs(sol.native["u"] - u_f)) / scale
    print(f"\n[fbar_trimmed route {mitigation}] n_dof={model['n_dof']} groups={model['n_groups']}/"
          f"{d['n_active']} max|u_h|={np.max(np.abs(uq_ref)):.4e} rel_diff={err:.2e} "
          f"tip={err_tip:.2e} (reference load: {err_f:.2e})")
    assert err <= 1e-12, err
    assert err_tip <= 1e-12, err_tip
    assert err_f <= 1e-9, err_f
    assert np.all(sol.native["u"][model["inactive_dofs"]] == 0.0)
    not_solved = set(fixed_ref.tolist()) | set(model["inactive_dofs"].tolist())
    assert sol.total_dofs == model["n_dof"]
    assert sol.active_dofs == model["n_dof"] - len(not_solved)


def test_auto_routes_to_fbar_only_for_nu_at_least_049(linear_sentinel):
    """'auto': nu = 0.49 -> trimmed F-bar (identical to explicit 'fbar');
    nu = 0.4899 -> linear Kratos path (sentinel reached); 'none' -> linear path."""
    opt = dict(kind="full_iga_immersed", degree=2, n_cells=_PLATE["n_cells"])
    p49 = _plate_problem(nu=0.49)
    sol_auto = solve_with_discretization_method(
        p49, FullIgaOption(**opt, locking_mitigation="auto"), use_pardiso=False)
    sol_fbar = solve_with_discretization_method(
        p49, FullIgaOption(**opt, locking_mitigation="fbar"), use_pardiso=False)
    assert sol_auto.native["locking_mitigation"] == "fbar"
    assert np.array_equal(sol_auto.native["u"], sol_fbar.native["u"])
    for nu, mit in ((0.4899, "auto"), (0.49, "none"), (0.3, "auto")):
        with pytest.raises(_LinearPathReached):
            solve_with_discretization_method(
                _plate_problem(nu=nu), FullIgaOption(**opt, locking_mitigation=mit),
                use_pardiso=False)


def test_immersed_fbar_box_soup_matches_body_fit_fbar():
    """Default soup (no ``tess_soup``: the box surface). Every background cell
    touches the box faces, so all 12 are "cut" and get moment-fitted rules
    (q = 4) instead of the body-fit route's 3-point tensor Gauss rule; every
    cell has eta = 1, so both pressure spaces are cell-wise P0. Both rules
    integrate the u = 0 tangent (degree <= 4 per axis for p = 2) exactly, so
    the two discrete problems share their linearization and differ only in
    the quadrature of the O(eps^2) terms of the Neo-Hookean residual
    (eps = max |grad u|). Bound: rel difference < eps^2 (computed from the
    solution, ~3e-5 here), while a wrong load, clamp or face check would
    give O(1)."""
    problem = BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(4.0, 1.0, 1.0), E=240.565, nu=0.49999,
                               fixed_face=(0, 0), load_face=(0, 1), traction=(0.0, 0.01, 0.0))
    kw = dict(degree=2, n_cells=(3, 2, 2), locking_mitigation="fbar")
    imm = solve_with_discretization_method(problem, FullIgaOption(kind="full_iga_immersed", **kw))
    bf = solve_with_discretization_method(problem, FullIgaOption(kind="full_iga_body_fit", **kw))
    model = imm.native["model"]
    d = model["diag"]
    assert d["n_cut"] == 12 and d["n_inactive_ctrl"] == 0 and d["n_groups"] == 12
    u_i, u_b = imm.native["u"], bf.native["u"]
    Uq = u_i.reshape(-1, 3)[model["cell_nodes"]][model["q_cell"]]
    eps = float(np.max(np.abs(np.einsum("qia,qim->qam", Uq, model["qG"]))))
    rel = np.max(np.abs(u_i - u_b)) / np.max(np.abs(u_b))
    print(f"\n[box soup] max|u|={np.max(np.abs(u_b)):.4e} eps={eps:.2e} "
          f"immersed vs body-fit rel={rel:.2e}")
    assert 1e-4 < eps < 1e-2  # moderate strain: the bound is meaningful
    assert rel < eps ** 2, (rel, eps)


# ---------------------------------------------------------------------------
# unsupported configurations raise (no silent linear fallback)
# ---------------------------------------------------------------------------

def test_unsupported_immersed_fbar_configs_raise(linear_sentinel):
    imm = dict(kind="full_iga_immersed", degree=2, n_cells=_PLATE["n_cells"])
    with pytest.raises(ValidationError, match="out_mode"):
        FullIgaOption(**imm, locking_mitigation="fbar", out_mode="ghost")
    with pytest.raises(ValidationError, match="continuity"):
        FullIgaOption(**imm, locking_mitigation="fbar", continuity=0)
    with pytest.raises(ValidationError, match="bc_penalty"):
        FullIgaOption(**imm, locking_mitigation="fbar", bc_penalty=1e3)
    # 'auto' passes option validation (nu unknown) and must raise at solve time
    # for nu >= 0.49, never reach the linear path
    p = _plate_problem(nu=0.4999)
    for extra, match in (({"out_mode": "ghost"}, "out_mode"), ({"continuity": 0}, "continuity"),
                         ({"bc_penalty": 1e3}, "bc_penalty")):
        opt = FullIgaOption(**imm, locking_mitigation="auto", **extra)
        with pytest.raises(ValueError, match=match):
            solve_with_discretization_method(p, opt)
    # the same 'auto' options below nu = 0.49 keep the linear path
    with pytest.raises(_LinearPathReached):
        solve_with_discretization_method(
            _plate_problem(nu=0.3), FullIgaOption(**imm, locking_mitigation="auto", out_mode="ghost"))
    # BC faces outside the solid: x = 0 crosses the hole (clamp), z = 0.5 too (load)
    fb = FullIgaOption(**imm, locking_mitigation="fbar")
    with pytest.raises(ValueError, match="fixed_face=\\(0, 0\\).*outside the trimmed solid"):
        solve_with_discretization_method(_plate_problem(nu=0.4999, fixed_face=(0, 0)), fb)
    with pytest.raises(ValueError, match="load_face=\\(2, 1\\).*outside the trimmed solid"):
        solve_with_discretization_method(_plate_problem(nu=0.4999, load_face=(2, 1)), fb)
    # fixed and prescribed overlap; prescribing an inactive control point
    with pytest.raises(ValueError, match="both fixed and prescribed"):
        solve_fbar_trimmed_box(p, n_cells=_PLATE["n_cells"], prescribed={0: 1e-3},
                               fixed_dofs=[0], clamp_fixed_face=False, use_pardiso=False)
    with pytest.raises(ValueError, match="inactive"):
        solve_fbar_trimmed_box(p, n_cells=_PLATE["n_cells"], prescribed={0: 1e-3},
                               use_pardiso=False)  # control point 0 sits in the hole


# ---------------------------------------------------------------------------
# V07-like Kirsch case through the methods API
# ---------------------------------------------------------------------------

def test_v07_kirsch_h05_through_methods_api_reproduces_fbar_trimmed():
    """V07 (nu = 0.49999, h = 0.5, q = 3 as test_fbar_trimmed) through
    ``solve_fbar_trimmed_box`` with the precursor's Kirsch BCs.

    The facade's errors against the Kirsch closed form must equal those of a
    direct fbar_trimmed solve of the same model (rel 1e-9: same discrete
    problem, the gap is the solver round-off amplified by cond(DKD) ~ 2e5),
    and meet the v07.json source-doc gate (e_u < 0.02, e_sigma_xx < 0.25).
    Both pressure spaces (agglomerated 0.25 and cell-wise 0) are checked.
    """
    nodes, tris = plate_hole_soup(A_V07, 4.37, 0.37)
    problem = BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(B_V07, B_V07, T_THICK), E=E_V07, nu=NU_V07,
                               tess_soup=(nodes, tris))
    n_cells = (8, 8, 2)
    box = BoxSpline(lo=problem.lo, hi=problem.hi, n_cells=n_cells, degree=2)
    qp, out, _iss, _rec, diag = build_trimmed_table(nodes, tris, box, q=3, workers=1, cache_dir="off")
    assert diag["n_fallback"] == 0 and diag["n_relaxed"] == 0
    fixed, prescribed = kirsch_dirichlet(box)
    for thr in (0.25, 0.0):
        sol = solve_fbar_trimmed_box(problem, n_cells=n_cells, quadrature_degree=3, eta_thresh=thr,
                                     clamp_fixed_face=False, fixed_dofs=fixed,
                                     prescribed=prescribed, use_pardiso=False)
        model_m = sol.native["model"]
        e_u, e_s = kirsch_errors(model_m, sol.native["u"])
        model = FT.build_model(box, qp, out, E_V07, NU_V07, eta_thresh=thr)
        u_ref = FT.solve_static(model, np.zeros(model["n_dof"]), fixed, n_steps=1,
                                prescribed=prescribed, linear_solver="superlu")
        e_u0, e_s0 = kirsch_errors(model, u_ref)
        print(f"\n[V07 h=0.5 eta_thresh={thr}] groups={model['n_groups']}/{model['diag']['n_active']} "
              f"methods e_u={e_u:.6e} e_sigma_xx={e_s:.6e} | direct e_u={e_u0:.6e} "
              f"e_sigma_xx={e_s0:.6e}")
        assert abs(e_u - e_u0) <= 1e-9 * e_u0
        assert abs(e_s - e_s0) <= 1e-9 * e_s0
        assert e_u < 0.02 and e_s < 0.25
        assert np.max(np.abs(sol.native["u"] - u_ref)) <= 1e-9 * np.max(np.abs(u_ref))
