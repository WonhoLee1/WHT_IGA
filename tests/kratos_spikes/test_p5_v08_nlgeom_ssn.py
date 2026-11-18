# -*- coding: utf-8 -*-
"""V08 Cases R1 & R2 with NLGEOM (Total-Lagrangian Neo-Hookean) bulk + SSN contact.

Tests:
1. Two p=2 cubes (2x2x2) with TotalLagrangianElement3D4N + CL.HyperElastic3DLaw (nu=0):
   compression (-0.05), unilateral separation (+0.02), and Coulomb friction (mu=0.3).
2. V08 Case R1 NLGEOM at h=0.25, sym='half', n=48, eps_n=1e10, t=0.03->0.1 (8 steps).
3. V08 Case R2 NLGEOM at h=0.25, sym='half', mu=0.3, t=0->0.40, sliding window 0.15<=t<=0.40.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import DeformableBoxFaceMaster, sample_surface_triangles_gauss
from wht_iga.analyses.contact_nlgeom import KratosBodyOperator, solve_nlgeom_contact_step
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"


def _box_face_triangles(box: BoxSpline, axis: int = 1, side: int = 0):
    """Build triangulated mesh ``(nodes, tris)`` on a BoxSpline boundary face."""
    others = [d for d in range(3) if d != axis]
    b0 = box.breakpoints(others[0])
    b1 = box.breakpoints(others[1])
    x_fixed = float(box.lo[axis] if side == 0 else box.hi[axis])
    n0, n1 = len(b0), len(b1)
    nodes = np.empty((n0 * n1, 3), dtype=np.float64)
    for j, v1 in enumerate(b1):
        for i, v0 in enumerate(b0):
            idx = i + n0 * j
            nodes[idx, axis] = x_fixed
            nodes[idx, others[0]] = float(v0)
            nodes[idx, others[1]] = float(v1)
    tris = []
    for j in range(n1 - 1):
        for i in range(n0 - 1):
            a = i + n0 * j
            b = (i + 1) + n0 * j
            c = (i + 1) + n0 * (j + 1)
            d = i + n0 * (j + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    return nodes, np.asarray(tris, dtype=np.int64)


def _prescribed_vector(presc_mask: np.ndarray, ux: float, uy: float, uz: float = 0.0) -> np.ndarray:
    """Build prescribed-DOF vector matching ``presc_mask`` component indices."""
    idx = np.flatnonzero(presc_mask)
    u_p = np.zeros(idx.size, dtype=np.float64)
    comp = idx % 3
    u_p[comp == 0] = float(ux)
    u_p[comp == 1] = float(uy)
    u_p[comp == 2] = float(uz)
    return u_p


def _cut_face_reaction(body_s: KratosBodyOperator, res: dict) -> np.ndarray:
    """Recover total ``(Fx, Fy, Fz)`` reaction on the prescribed cut face."""
    pm = body_s.presc_mask
    R_p = res["f_int_s"][pm] + res["f_c_s_full"][pm]
    idx = np.flatnonzero(pm) % 3
    return np.array(
        [float(R_p[idx == 0].sum()), float(R_p[idx == 1].sum()), float(R_p[idx == 2].sum())],
        dtype=np.float64,
    )


def _setup_v08_nlgeom_half(KM, SM, CL, inp: H.V08Inputs, h: float, n_soup: int = 48):
    """Build NLGEOM trimmed deformable indenter + deformable block for V08 sym='half'."""
    lo_b, hi_b, nc_b = H.v08_block_grid(h, sym="half")
    box_blk = BoxSpline(lo_b, hi_b, nc_b, degree=2)

    def _blk_zero(x):
        if abs(x[1] - lo_b[1]) < 1e-12:
            return (0, 1, 2)
        if abs(x[2] - lo_b[2]) < 1e-12:
            return (2,)
        return ()

    body_m = KratosBodyOperator(
        KM, SM, box_blk, inp.E_block, inp.nu_block,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=_blk_zero,
    )

    nodes_ind, tris_ind, fids_ind = H.hollow_hemisphere_soup(
        inp.R, inp.Ri, inp.center, n=n_soup, sym="half"
    )
    lo_i, hi_i, nc_i = H.v08_indenter_grid(h, inp.center, inp.R, sym="half")
    box_ind = BoxSpline(lo_i, hi_i, nc_i, degree=2)

    body_s = KratosBodyOperator(
        KM, SM, box_ind, inp.E_indenter, inp.nu_indenter,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (2,) if abs(x[2] - lo_i[2]) < 1e-12 else (),
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - hi_i[1]) < 1e-12 else (),
        trim_soup=(nodes_ind, tris_ind),
        q=4,
        ghost_alpha=1e-8,
    )
    assert body_s.diag["n_fallback"] == 0

    tris_outer = tris_ind[fids_ind == H.FACE_OUTER]
    X_s, w_q, B_s = sample_surface_triangles_gauss(nodes_ind, tris_outer, box_ind, order=2)
    master_face = DeformableBoxFaceMaster(box_blk, axis=1, side=1)
    return box_blk, box_ind, body_m, body_s, X_s, w_q, B_s, master_face


def test_nlgeom_two_block_compression_separation_and_friction():
    """Test 1: two p=2 NLGEOM cubes: compression, separation, and Coulomb friction."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    E1, E2, nu = 1000.0, 2000.0, 0.0
    eps_n, eps_t = 1e6, 1e5
    delta = 0.05

    box_m = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    box_s = BoxSpline(lo=(0.0, 1.0, 0.0), hi=(1.0, 2.0, 1.0), n_cells=(2, 2, 2), degree=2)

    body_m = KratosBodyOperator(
        KM, SM, box_m, E1, nu,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
    )
    body_s = KratosBodyOperator(
        KM, SM, box_s, E2, nu,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - 2.0) < 1e-12 else (),
    )

    tri_nodes, tris = _box_face_triangles(box_s, axis=1, side=0)
    X_s, w_q, B_s = sample_surface_triangles_gauss(tri_nodes, tris, box_s, order=2)
    assert float(w_q.sum()) == pytest.approx(1.0, rel=1e-14)
    master_face = DeformableBoxFaceMaster(box_m, axis=1, side=1)

    # 1a. Compression -0.05
    um_full = np.zeros(box_m.n_dofs, dtype=np.float64)
    us_full = np.zeros(box_s.n_dofs, dtype=np.float64)
    us_full[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, 0.0, -delta, 0.0)

    res_comp = solve_nlgeom_contact_step(
        body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
        um_full, us_full, outer_tol=1e-8, use_pardiso=False,
    )
    assert res_comp["outer_rel"] <= 1e-8
    assert res_comp["n_active"] == len(w_q)
    R_comp = _cut_face_reaction(body_s, res_comp)
    F_N_comp = float(-res_comp["force_total"][1])
    assert F_N_comp > 0.0
    assert abs(R_comp[1] - float(res_comp["force_total"][1])) / F_N_comp < 1e-8
    assert abs(R_comp[0]) < 1e-8 and abs(R_comp[2]) < 1e-8

    # 1b. Lift +0.02 -> unilateral separation
    us_up = res_comp["us_full"].copy()
    us_up[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, 0.0, 0.02, 0.0)
    res_sep = solve_nlgeom_contact_step(
        body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
        res_comp["um_full"], us_up, outer_tol=1e-8, use_pardiso=False,
    )
    assert res_sep["outer_rel"] <= 1e-8
    assert res_sep["n_active"] == 0
    assert float(np.linalg.norm(res_sep["force_total"])) == 0.0

    # 1c. Press (-0.05) + tangential slide (mu = 0.3) on flat-master cubes (E1=1e9, E2=2000)
    mu = 0.3
    body_m_stiff = KratosBodyOperator(
        KM, SM, box_m, 1.0e9, nu,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
    )
    um_f = np.zeros(box_m.n_dofs, dtype=np.float64)
    us_f = np.zeros(box_s.n_dofs, dtype=np.float64)
    us_f[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, 0.0, -delta, 0.0)
    res_f = solve_nlgeom_contact_step(
        body_m_stiff, body_s, X_s, w_q, B_s, master_face, eps_n,
        um_f, us_f, mu=mu, eps_t=eps_t, outer_tol=1e-8, use_pardiso=False,
    )
    um_f, us_f, slip_state = res_f["um_full"], res_f["us_full"], res_f["slip_state"]

    for step in range(1, 26):
        ux = 0.005 * step
        us_f = us_f.copy()
        us_f[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, ux, -delta, 0.0)
        res_f = solve_nlgeom_contact_step(
            body_m_stiff, body_s, X_s, w_q, B_s, master_face, eps_n,
            um_f, us_f, mu=mu, eps_t=eps_t, slip_state=slip_state,
            outer_tol=1e-8, use_pardiso=False,
        )
        um_f, us_f, slip_state = res_f["um_full"], res_f["us_full"], res_f["slip_state"]

    fn_f = float(-res_f["force_total"][1])
    ft_f = float(np.hypot(res_f["force_total"][0], res_f["force_total"][2]))
    ratio = ft_f / fn_f
    R_fric = _cut_face_reaction(body_s, res_f)
    eq_x = abs(R_fric[0] - float(res_f["force_total"][0])) / max(ft_f, 1e-30)
    eq_y = abs(R_fric[1] - float(res_f["force_total"][1])) / max(fn_f, 1e-30)

    print(
        f"\n[Test 1 NLGEOM] comp F_N={F_N_comp:.6f} (n_outer={res_comp['n_outer']}), "
        f"sep n_act={res_sep['n_active']}, "
        f"fric |F_T|/F_N={ratio:.6f} (n_slip={res_f['n_slip']}/{res_f['n_active']}), "
        f"eq_x={eq_x:.2e}, eq_y={eq_y:.2e}"
    )
    assert res_f["n_slip"] > 0
    assert ratio == pytest.approx(mu, rel=1e-4)
    assert eq_x < 1e-8
    assert eq_y < 1e-8


def test_v08_r1_nlgeom_h1():
    """Test 2: V08 Case R1 NLGEOM at h=0.25, sym='half', n=48, eps_n=1e10, t=0.03->0.1."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    F_hertz = H.v08_reference(inp)["normal_force_F_N_lbf"]
    F_sneddon = H.v08_model_form_gaps(inp)["exact_F_N_lbf"]

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    h = 0.25
    n_soup = 48
    eps_n = 1e10
    t_steps = np.linspace(0.03, 0.1, 8)

    t0 = time.perf_counter()
    _bb, box_ind, body_m, body_s, X_s, w_q, B_s, master_face = _setup_v08_nlgeom_half(
        KM, SM, CL, inp, h, n_soup=n_soup
    )
    t_setup = time.perf_counter() - t0

    um_full = np.zeros(body_m.box.n_dofs, dtype=np.float64)
    us_full = np.zeros(box_ind.n_dofs, dtype=np.float64)

    t_s = time.perf_counter()
    rows = []
    for t in t_steps:
        ux, uy = H.v08_cut_face_displacement(float(t))
        us_full = us_full.copy()
        us_full[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, float(ux), float(uy), 0.0)
        res = solve_nlgeom_contact_step(
            body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
            um_full, us_full, outer_tol=1e-8, rtol=1e-9, use_pardiso=True,
        )
        um_full, us_full = res["um_full"], res["us_full"]
        fn_iface = 2.0 * float(-res["force_total"][1])
        R_cut = _cut_face_reaction(body_s, res)
        fn_cut = 2.0 * float(-R_cut[1])
        eq_rel = abs(fn_cut - fn_iface) / fn_iface
        rows.append((float(t), res["n_outer"], res["iterations"], res["n_active"], fn_iface, fn_cut, eq_rel))
    dt_solve = time.perf_counter() - t_s

    fn_final = rows[-1][4]
    max_eq = max(r[6] for r in rows)
    max_outer = max(r[1] for r in rows)
    err_h = (fn_final / F_hertz - 1.0) * 100.0
    err_s = (fn_final / F_sneddon - 1.0) * 100.0

    print(
        f"\n[Test 2 V08 R1 NLGEOM] h={h}, n={n_soup}, eps_n={eps_n:.0e}, "
        f"setup={t_setup:.2f}s, solve={dt_solve:.2f}s, max_outer={max_outer}: "
        f"F_N(0.1)={fn_final:,.1f} lbf, vs Hertz {err_h:+.2f}%, vs Sneddon {err_s:+.2f}%, "
        f"max_eq_rel={max_eq:.2e}"
    )
    assert body_s.diag["n_fallback"] == 0
    assert max_eq < 1e-5
    assert abs(err_h) < 4.0


def test_v08_r2_nlgeom_h1():
    """Test 3: V08 Case R2 NLGEOM at h=0.25, sym='half', mu=0.3, t=0->0.40."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    h = 0.25
    n_soup = 24
    eps_n, eps_t = 1e10, 1e9
    t_indent = np.linspace(0.03, 0.1, 6)
    t_slide = np.linspace(0.115, 0.40, 20)
    t_all = np.concatenate([t_indent, t_slide])

    t0 = time.perf_counter()
    _bb, box_ind, body_m, body_s, X_s, w_q, B_s, master_face = _setup_v08_nlgeom_half(
        KM, SM, CL, inp, h, n_soup=n_soup
    )
    t_setup = time.perf_counter() - t0

    um_full = np.zeros(body_m.box.n_dofs, dtype=np.float64)
    us_full = np.zeros(box_ind.n_dofs, dtype=np.float64)
    slip_state = None

    t_s = time.perf_counter()
    rows = []
    for t in t_all:
        ux, uy = H.v08_cut_face_displacement(float(t))
        us_full = us_full.copy()
        us_full[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, float(ux), float(uy), 0.0)
        res = solve_nlgeom_contact_step(
            body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
            um_full, us_full, mu=inp.mu, eps_t=eps_t, slip_state=slip_state,
            outer_tol=1e-8, rtol=1e-8, use_pardiso=True,
        )
        um_full, us_full, slip_state = res["um_full"], res["us_full"], res["slip_state"]
        fn_iface = 2.0 * float(-res["force_total"][1])
        ft_iface = 2.0 * float(res["force_total"][0])
        R_cut = _cut_face_reaction(body_s, res)
        fn_cut = 2.0 * float(-R_cut[1])
        ft_cut = 2.0 * float(R_cut[0])
        rows.append((
            float(t), res["n_outer"], res["iterations"], res["n_active"], res["n_slip"],
            fn_iface, ft_iface, fn_cut, ft_cut,
        ))
    dt_solve = time.perf_counter() - t_s

    t_arr = np.array([r[0] for r in rows])
    fn_arr = np.array([r[5] for r in rows])
    ft_arr = np.array([r[6] for r in rows])
    fn_cut_arr = np.array([r[7] for r in rows])
    ft_cut_arr = np.array([r[8] for r in rows])

    win = (t_arr >= 0.15) & (t_arr <= 0.40)
    ref = H.v08_reference_history(t_arr[win], inp)
    rms_fn = H.relative_rms(fn_arr[win], ref["F_N"])
    rms_ft = H.relative_rms(ft_arr[win], ref["F_T"])
    ratios_win = ft_arr[win] / fn_arr[win]
    eq_fn = float(np.max(np.abs(fn_cut_arr[win] - fn_arr[win]) / fn_arr[win]))
    eq_ft = float(np.max(np.abs(ft_cut_arr[win] - ft_arr[win]) / ft_arr[win]))
    max_outer = max(r[1] for r in rows)

    print(
        f"\n[Test 3 V08 R2 NLGEOM] h={h}, n={n_soup}, eps_n={eps_n:.0e}, eps_t={eps_t:.0e}, "
        f"setup={t_setup:.2f}s, solve={dt_solve:.2f}s, max_outer={max_outer}: "
        f"F_N_win={fn_arr[win].mean():,.1f} lbf, F_T_win={ft_arr[win].mean():,.1f} lbf, "
        f"F_T/F_N in [{ratios_win.min():.4f}, {ratios_win.max():.4f}] (mean={ratios_win.mean():.4f}), "
        f"eq_err=(N:{eq_fn:.2e}, T:{eq_ft:.2e}), "
        f"rel_RMS=(F_N:{rms_fn*100:.2f}%, F_T:{rms_ft*100:.2f}%)"
    )
    assert body_s.diag["n_fallback"] == 0
    assert eq_fn < 1e-5
    assert eq_ft < 2e-3
    assert abs(float(ratios_win.mean()) - inp.mu) < 0.015
    assert abs(float(ratios_win[-1]) - inp.mu) < 0.015
    assert rms_fn < 0.06  # h=0.25 coarse NLGEOM ~5.09%; G5 4% gate requires h=0.125
    assert rms_ft < 0.05  # h=0.25 coarse NLGEOM ~4.23% (< 4.5% G5 gate)
