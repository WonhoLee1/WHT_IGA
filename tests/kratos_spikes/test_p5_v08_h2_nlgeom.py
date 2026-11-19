# -*- coding: utf-8 -*-
"""V08 h=0.125 R1 NLGEOM (G5 gate) & h=0.25 R2 NLGEOM full sliding window (0.15<=t<=0.70).

Tests:
1. ``test_v08_r1_nlgeom_h2``: V08 Case R1 NLGEOM at ``h=0.125``, ``sym='half'``,
   ``n_soup=48``, ``eps_n=1e10``, ``t=0.03->0.10`` (8 steps), ``outer_tol=1e-7, rtol=1e-8``.
2. ``test_v08_r2_nlgeom_full_window_h1``: V08 Case R2 NLGEOM at ``h=0.25``, ``sym='half'``,
   ``n_soup=24``, ``mu=0.3``, ``eps_n=1e10``, ``eps_t=1e9``, ``t=0.03->0.70``, verifying
   both ``0.15<=t<=0.40`` and the full ``v08.json`` ``0.15<=t<=0.70`` sliding window.
"""

from __future__ import annotations

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import psutil
import pytest

from wht_iga.analyses.contact import DeformableBoxFaceMaster, sample_surface_triangles_gauss
from wht_iga.analyses.contact_nlgeom import KratosBodyOperator, solve_nlgeom_contact_step
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"
_GIB = 1024.0 ** 3


def _rss_gb() -> float:
    """Return current process RSS (working set) in GiB."""
    return float(psutil.Process(os.getpid()).memory_info().rss) / _GIB


def _peak_wset_gb() -> float:
    """Return cumulative OS process peak working set in GiB."""
    return float(psutil.Process(os.getpid()).memory_info().peak_wset) / _GIB


def _release_solver_memory() -> None:
    """Release cached MKL Pardiso factorization and collect Python/Kratos objects."""
    try:
        from pypardiso.scipy_aliases import pypardiso_solver

        pypardiso_solver.free_memory(everything=True)
    except Exception:
        pass
    gc.collect()


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


def _zero_mask_reaction(body_s: KratosBodyOperator, res: dict) -> np.ndarray:
    """Recover total ``(Fx, Fy, Fz)`` reaction on ``zero_mask`` (inactive ghost DOFs for x, y)."""
    zm = body_s.zero_mask
    R_z = res["f_int_s"][zm] + res["f_c_s_full"][zm]
    idx = np.flatnonzero(zm) % 3
    return np.array(
        [float(R_z[idx == 0].sum()), float(R_z[idx == 1].sum()), float(R_z[idx == 2].sum())],
        dtype=np.float64,
    )


def _setup_v08_nlgeom_half(KM, SM, CL, inp: H.V08Inputs, h: float, n_soup: int):
    """Build NLGEOM trimmed deformable indenter + deformable block for V08 ``sym='half'``."""
    lo_b, hi_b, nc_b = H.v08_block_grid(h, sym="half")
    box_blk = BoxSpline(lo_b, hi_b, nc_b, degree=2)

    def _blk_zero(x: np.ndarray) -> tuple[int, ...]:
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

    tris_outer = tris_ind[fids_ind == H.FACE_OUTER]
    X_s, w_q, B_s = sample_surface_triangles_gauss(nodes_ind, tris_outer, box_ind, order=2)
    master_face = DeformableBoxFaceMaster(box_blk, axis=1, side=1)
    return box_blk, box_ind, body_m, body_s, X_s, w_q, B_s, master_face


def test_v08_r1_nlgeom_h2():
    """V08 Case R1 NLGEOM at h=0.125, sym='half', n_soup=48, eps_n=1e10, t=0.03->0.10."""
    with open(_V08_JSON, encoding="utf-8") as f:
        v08_data = json.load(f)
    inp = H.V08Inputs.from_json(v08_data)

    ref_v08 = H.v08_reference(inp)
    F_hertz = float(v08_data["reference_values"]["Hertz_Coulomb_accepted_at_full_indentation"]["normal_force_F_N_lbf"])
    p0_hertz = float(ref_v08["maximum_pressure_p0_psi"])
    F_sneddon = H.v08_model_form_gaps(inp)["exact_F_N_lbf"]

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    h = 0.125
    n_soup = 48
    eps_n = 1e10
    t_steps = np.linspace(0.03, 0.10, 8)

    rss_start = _rss_gb()
    peak_step_rss = rss_start
    rows: list[tuple[float, float, int, int, int, float, float, float, float, float]] = []
    n_fallback = -1
    free_m = free_s = total_m = total_s = n_gp = 0
    t_setup = dt_solve = 0.0

    box_blk = box_ind = body_m = body_s = X_s = w_q = B_s = master_face = None
    um_full = us_full = res = None
    try:
        t0 = time.perf_counter()
        box_blk, box_ind, body_m, body_s, X_s, w_q, B_s, master_face = _setup_v08_nlgeom_half(
            KM, SM, CL, inp, h, n_soup=n_soup
        )
        t_setup = time.perf_counter() - t0
        n_fallback = int(body_s.diag["n_fallback"])
        free_m = int(body_m.free_mask.sum())
        free_s = int(body_s.free_mask.sum())
        total_m = int(box_blk.n_dofs)
        total_s = int(box_ind.n_dofs)
        n_gp = int(len(w_q))
        peak_step_rss = max(peak_step_rss, _rss_gb())

        print(
            f"\n[R1 h=0.125 setup] dofs=({free_m}+{free_s})/{total_m}+{total_s}, "
            f"n_gp={n_gp}, n_fb={n_fallback}, setup={t_setup:.2f}s, rss={_rss_gb():.3f} GiB",
            flush=True,
        )

        um_full = np.zeros(box_blk.n_dofs, dtype=np.float64)
        us_full = np.zeros(box_ind.n_dofs, dtype=np.float64)

        t_s = time.perf_counter()
        for t in t_steps:
            ux, uy = H.v08_cut_face_displacement(float(t))
            us_full = us_full.copy()
            us_full[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, float(ux), float(uy), 0.0)
            res = solve_nlgeom_contact_step(
                body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
                um_full, us_full, outer_tol=1e-7, rtol=1e-8, use_pardiso=True,
            )
            um_full, us_full = res["um_full"], res["us_full"]
            fn_iface = 2.0 * float(-res["force_total"][1])
            R_cut = _cut_face_reaction(body_s, res)
            fn_cut = 2.0 * float(-R_cut[1])
            eq_rel = abs(fn_cut - fn_iface) / fn_iface
            p_max = float(res["max_pressure"])
            rss_now = _rss_gb()
            peak_step_rss = max(peak_step_rss, rss_now)
            delta = float(H.v08_delta(float(t), inp.g0))
            rows.append((
                float(t), delta, int(res["n_outer"]), int(res["iterations"]),
                int(res["n_active"]), fn_iface, fn_cut, eq_rel, p_max, rss_now,
            ))
            print(
                f"  t={t:.4f} (delta={delta:.5f}): outer={res['n_outer']}, it={res['iterations']}, "
                f"n_act={res['n_active']}, F_N={fn_iface:,.1f} lbf, eq_rel={eq_rel:.2e}, "
                f"p_max={p_max:,.0f} psi, rss={rss_now:.3f} GiB",
                flush=True,
            )
        dt_solve = time.perf_counter() - t_s
    finally:
        del res, um_full, us_full, body_m, body_s, master_face, B_s, X_s, w_q, box_blk, box_ind
        _release_solver_memory()

    fn_final = rows[-1][5]
    fn_cut_final = rows[-1][6]
    pmax_final = rows[-1][8]
    max_eq = max(r[7] for r in rows)
    total_outer = sum(r[2] for r in rows)
    max_outer = max(r[2] for r in rows)
    err_h = (fn_final / F_hertz - 1.0) * 100.0
    err_s = (fn_final / F_sneddon - 1.0) * 100.0
    err_p0 = (pmax_final / p0_hertz - 1.0) * 100.0
    d_arr = np.array([r[1] for r in rows])
    fn_arr = np.array([r[5] for r in rows])
    slope, _ = np.polyfit(np.log(d_arr), np.log(fn_arr), 1)
    peak_wset = _peak_wset_gb()

    print(
        f"[R1 h=0.125 summary] h={h}, n_soup={n_soup}, eps_n={eps_n:.0e}, "
        f"setup={t_setup:.2f}s, solve={dt_solve:.2f}s, total_wall={t_setup + dt_solve:.2f}s, "
        f"sum_outer={total_outer} (max={max_outer}), slope={slope:.4f}, "
        f"F_N(0.1)={fn_final:,.1f} lbf (cut={fn_cut_final:,.1f} lbf), "
        f"vs Hertz {err_h:+.2f}%, vs Sneddon {err_s:+.2f}%, "
        f"p_max={pmax_final:,.1f} psi (vs Hertz p0 {err_p0:+.2f}%), "
        f"max_eq_rel={max_eq:.2e}, peak_step_rss={peak_step_rss:.3f} GiB, "
        f"peak_wset={peak_wset:.3f} GiB, post_free_rss={_rss_gb():.3f} GiB",
        flush=True,
    )

    assert n_fallback == 0
    assert max_eq < 1e-4
    assert abs(err_h) < 4.0


def test_v08_r2_nlgeom_full_window_h1():
    """V08 Case R2 NLGEOM at h=0.25, sym='half', n_soup=24, mu=0.3, t=0.03->0.70."""
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
    t_slide = np.linspace(0.115, 0.70, 40)
    t_all = np.concatenate([t_indent, t_slide])

    rss_start = _rss_gb()
    peak_step_rss = rss_start
    rows: list[tuple[float, int, int, int, int, float, float, float, float, float, float, float, float]] = []
    n_fallback = -1
    free_m = free_s = 0
    t_setup = dt_solve = 0.0

    box_blk = box_ind = body_m = body_s = X_s = w_q = B_s = master_face = None
    um_full = us_full = slip_state = res = None
    try:
        t0 = time.perf_counter()
        box_blk, box_ind, body_m, body_s, X_s, w_q, B_s, master_face = _setup_v08_nlgeom_half(
            KM, SM, CL, inp, h, n_soup=n_soup
        )
        t_setup = time.perf_counter() - t0
        n_fallback = int(body_s.diag["n_fallback"])
        free_m = int(body_m.free_mask.sum())
        free_s = int(body_s.free_mask.sum())
        peak_step_rss = max(peak_step_rss, _rss_gb())

        print(
            f"\n[R2 h=0.25 full setup] dofs=({free_m}+{free_s}), n_gp={len(w_q)}, "
            f"n_fb={n_fallback}, setup={t_setup:.2f}s, rss={_rss_gb():.3f} GiB",
            flush=True,
        )

        um_full = np.zeros(box_blk.n_dofs, dtype=np.float64)
        us_full = np.zeros(box_ind.n_dofs, dtype=np.float64)
        slip_state = None

        t_s = time.perf_counter()
        for step_idx, t in enumerate(t_all, start=1):
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
            R_zero = _zero_mask_reaction(body_s, res)
            fn_cut = 2.0 * float(-R_cut[1])
            ft_cut = 2.0 * float(R_cut[0])
            fn_tot = 2.0 * float(-(R_cut[1] + R_zero[1]))
            ft_tot = 2.0 * float(R_cut[0] + R_zero[0])
            p_max = float(res["max_pressure"])
            rss_now = _rss_gb()
            peak_step_rss = max(peak_step_rss, rss_now)
            rows.append((
                float(t), int(res["n_outer"]), int(res["iterations"]),
                int(res["n_active"]), int(res["n_slip"]),
                fn_iface, ft_iface, fn_cut, ft_cut, fn_tot, ft_tot, p_max, rss_now,
            ))
            if step_idx in (1, 6) or step_idx % 5 == 0 or abs(float(t) - 0.40) < 1e-9 or step_idx == len(t_all):
                ratio = ft_iface / fn_iface if fn_iface > 0 else 0.0
                print(
                    f"  step {step_idx:2d}/{len(t_all)} t={t:.4f}: outer={res['n_outer']}, "
                    f"it={res['iterations']}, act={res['n_active']}, slip={res['n_slip']}, "
                    f"F_N={fn_iface:,.1f}, F_T={ft_iface:,.1f}, F_T/F_N={ratio:.4f}, "
                    f"p_max={p_max:,.0f}, rss={rss_now:.3f} GiB",
                    flush=True,
                )
        dt_solve = time.perf_counter() - t_s
    finally:
        del res, um_full, us_full, slip_state, body_m, body_s, master_face, B_s, X_s, w_q, box_blk, box_ind
        _release_solver_memory()

    t_arr = np.array([r[0] for r in rows])
    fn_arr = np.array([r[5] for r in rows])
    ft_arr = np.array([r[6] for r in rows])
    fn_cut_arr = np.array([r[7] for r in rows])
    ft_cut_arr = np.array([r[8] for r in rows])
    fn_tot_arr = np.array([r[9] for r in rows])
    ft_tot_arr = np.array([r[10] for r in rows])
    pmax_arr = np.array([r[11] for r in rows])
    max_outer = max(r[1] for r in rows)
    total_outer = sum(r[1] for r in rows)
    peak_wset = _peak_wset_gb()

    win_stats = {}
    for label, lo_t, hi_t in (("0.15<=t<=0.40", 0.15, 0.40), ("0.15<=t<=0.70", 0.15, 0.70)):
        win = (t_arr >= lo_t - 1e-12) & (t_arr <= hi_t + 1e-12)
        ref = H.v08_reference_history(t_arr[win], inp)
        rms_fn = H.relative_rms(fn_arr[win], ref["F_N"])
        rms_ft = H.relative_rms(ft_arr[win], ref["F_T"])
        rms_p0 = H.relative_rms(pmax_arr[win], ref["p0"])
        ratios_win = ft_arr[win] / fn_arr[win]
        eq_fn = float(np.max(np.abs(fn_cut_arr[win] - fn_arr[win]) / fn_arr[win]))
        eq_ft = float(np.max(np.abs(ft_cut_arr[win] - ft_arr[win]) / ft_arr[win]))
        eq_fn_tot = float(np.max(np.abs(fn_tot_arr[win] - fn_arr[win]) / fn_arr[win]))
        eq_ft_tot = float(np.max(np.abs(ft_tot_arr[win] - ft_arr[win]) / ft_arr[win]))
        win_stats[label] = {
            "n_frames": int(win.sum()),
            "fn_mean": float(fn_arr[win].mean()),
            "ft_mean": float(ft_arr[win].mean()),
            "ratio_min": float(ratios_win.min()),
            "ratio_max": float(ratios_win.max()),
            "ratio_mean": float(ratios_win.mean()),
            "ratio_end": float(ratios_win[-1]),
            "eq_fn": eq_fn,
            "eq_ft": eq_ft,
            "eq_fn_tot": eq_fn_tot,
            "eq_ft_tot": eq_ft_tot,
            "rms_fn": rms_fn,
            "rms_ft": rms_ft,
            "rms_p0": rms_p0,
        }
        print(
            f"[R2 window {label}] frames={int(win.sum())}, "
            f"F_N_mean={fn_arr[win].mean():,.1f} lbf, F_T_mean={ft_arr[win].mean():,.1f} lbf, "
            f"F_T/F_N in [{ratios_win.min():.4f}, {ratios_win.max():.4f}] (mean={ratios_win.mean():.4f}), "
            f"eq_cut=(N:{eq_fn:.2e}, T:{eq_ft:.2e}), eq_tot=(N:{eq_fn_tot:.2e}, T:{eq_ft_tot:.2e}), "
            f"rel_RMS=(F_N:{rms_fn*100:.2f}%, F_T:{rms_ft*100:.2f}%, p_max:{rms_p0*100:.2f}%)",
            flush=True,
        )

    print(
        f"[R2 h=0.25 summary] setup={t_setup:.2f}s, solve={dt_solve:.2f}s, "
        f"total_wall={t_setup + dt_solve:.2f}s, sum_outer={total_outer} (max={max_outer}), "
        f"peak_step_rss={peak_step_rss:.3f} GiB, peak_wset={peak_wset:.3f} GiB, "
        f"post_free_rss={_rss_gb():.3f} GiB",
        flush=True,
    )

    assert n_fallback == 0
    for label, st in win_stats.items():
        assert st["eq_fn"] < 1e-4, f"{label} eq_fn={st['eq_fn']}"
        assert st["eq_fn_tot"] < 1e-4, f"{label} eq_fn_tot={st['eq_fn_tot']}"
        # Whole-body Newton equilibrium (cut face + pinned inactive ghost DOFs) matches interface force:
        assert st["eq_ft_tot"] < 1e-4, f"{label} eq_ft_total={st['eq_ft_tot']}"
        # Cut-face-only eq_ft includes O((u_x/h)^3 * ghost_alpha * E) reaction absorbed by
        # inactive_ctrl_mask across OUT ghost cells under TotalLagrangianElement3D4N (~1.6e-3
        # at t=0.40, growing by ~2^3 to ~1.3e-2 at t=0.70 where u_x=2.833 in).
        eq_ft_bound = 5e-3 if label == "0.15<=t<=0.40" else 2e-2
        assert st["eq_ft"] < eq_ft_bound, f"{label} eq_ft={st['eq_ft']}"
        assert abs(st["ratio_mean"] - inp.mu) < 0.015, f"{label} ratio_mean={st['ratio_mean']}"
        assert abs(st["ratio_end"] - inp.mu) < 0.015, f"{label} ratio_end={st['ratio_end']}"
        assert st["rms_fn"] < 0.045, f"{label} rms_fn={st['rms_fn']}"
        assert st["rms_ft"] < 0.045, f"{label} rms_ft={st['rms_ft']}"
