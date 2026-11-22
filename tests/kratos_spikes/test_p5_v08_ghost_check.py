# -*- coding: utf-8 -*-
"""Fast V08 R2 ghost-cell hypothesis check at h=0.25 (sym='half', n_soup=24).

Compares ghost-cell configurations on ``OUT`` cells of the trimmed indenter in
``KratosBodyOperator`` across key R2 NLGEOM steps ``t = (0.03, 0.10, 0.40, 0.70)``:

1. ``ghost_alpha=1e-8, ghost_element="TotalLagrangianElement3D4N"``
   (pre-patch St. Venant-Kirchhoff ghost-cell force leak).
2. ``ghost_alpha=1e-8, ghost_element="SmallDisplacementElement3D4N"``
   (current ``KratosBodyOperator`` default with ``LinearElastic3DLaw`` on ``OUT``).
3. ``ghost_alpha=0.0`` (checks whether ``KratosBodyOperator`` supports omitting
   ``OUT`` elements cleanly with ``inactive_ctrl_mask`` in ``zero_mask``).

Also verifies the algebraic equilibrium identity:
``R_cut_x - force_total[0] == -(f_int_s[zero_mask_x].sum())``.
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
from wht_iga.solver.active import inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

pytestmark = [pytest.mark.kratos]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"
_MIB = 1024.0 ** 2
_T_STEPS = (0.03, 0.10, 0.40, 0.70)


def _rss_mb() -> float:
    """Return current process RSS in MiB."""
    return float(psutil.Process(os.getpid()).memory_info().rss) / _MIB


def _release_solver_memory() -> None:
    """Release cached MKL Pardiso memory and collect Python/Kratos objects."""
    try:
        from pypardiso.scipy_aliases import pypardiso_solver

        pypardiso_solver.free_memory(everything=True)
    except Exception:
        pass
    gc.collect()


def _prescribed_vector(presc_mask: np.ndarray, ux: float, uy: float, uz: float = 0.0) -> np.ndarray:
    """Build prescribed-DOF displacement vector ordered by ``presc_mask`` indices."""
    idx = np.flatnonzero(presc_mask)
    u_p = np.zeros(idx.size, dtype=np.float64)
    comp = idx % 3
    u_p[comp == 0] = float(ux)
    u_p[comp == 1] = float(uy)
    u_p[comp == 2] = float(uz)
    return u_p


def _component_mask(mask: np.ndarray, comp: int) -> np.ndarray:
    """Filter a full-DOF boolean mask to a single coordinate component (0=x, 1=y, 2=z)."""
    out = np.asarray(mask, dtype=bool).copy()
    dofs = np.arange(out.size)
    out &= (dofs % 3 == int(comp))
    return out


def _print_case_summary(tag: str, out: dict) -> None:
    """Print per-step ghost-cell measurements and memory/runtime summary."""
    print(
        f"\n[{tag}] n_elem_s={out['n_elem_s']}, sys_size={out['sys_size']}, "
        f"free_s={out['free_s']}, inact_ctrl={out['n_inact_ctrl']}, "
        f"wall={out['elapsed_s']:.2f}s, peak_rss={out['peak_rss_mb']:.1f} MiB",
        flush=True,
    )
    for r in out["rows"]:
        print(
            f"  t={r['t']:.2f} ux={r['ux']:.3f} outer={r['n_outer']} it={r['it']} "
            f"F_N={r['F_N']:,.2f} F_T={r['F_T']:,.2f} "
            f"Rcut_x-ft0={r['lhs_id']:+.6e} -fint_zx={r['rhs_id']:+.6e} "
            f"id_err={r['id_err']:.2e} eq_cut_N={r['eq_cut_N']:.3e} eq_cut_T={r['eq_cut_T']:.3e}",
            flush=True,
        )


def _run_ghost_case(
    KM,
    SM,
    CL,
    inp: H.V08Inputs,
    *,
    ghost_alpha: float,
    ghost_element: str = "SmallDisplacementElement3D4N",
    h: float = 0.25,
    n_soup: int = 24,
    t_steps: tuple[float, ...] = _T_STEPS,
) -> dict:
    """Run short V08 R2 NLGEOM sequence for one ``(ghost_alpha, ghost_element)`` setting."""
    eps_n, eps_t = 1e10, 1e9
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

    _, out_cells, _, _, _ = build_trimmed_table(nodes_ind, tris_ind, box_ind, q=4, cache_dir="off")
    inact = inactive_ctrl_mask(box_ind, out_cells)

    body_s = KratosBodyOperator(
        KM, SM, box_ind, inp.E_indenter, inp.nu_indenter,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (2,) if abs(x[2] - lo_i[2]) < 1e-12 else (),
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - hi_i[1]) < 1e-12 else (),
        trim_soup=(nodes_ind, tris_ind),
        q=4,
        ghost_alpha=float(ghost_alpha),
        ghost_element=ghost_element,
    )

    tris_outer = tris_ind[fids_ind == H.FACE_OUTER]
    X_s, w_q, B_s = sample_surface_triangles_gauss(nodes_ind, tris_outer, box_ind, order=2)
    master_face = DeformableBoxFaceMaster(box_blk, axis=1, side=1)

    zero_mask_x = _component_mask(body_s.zero_mask, 0)
    inact_x = np.zeros(box_ind.n_dofs, dtype=bool)
    inact_x[3 * np.flatnonzero(inact)] = True
    assert np.array_equal(zero_mask_x, inact_x), "x-direction zero_mask must equal inactive_ctrl_mask x-DOFs"

    presc_x = _component_mask(body_s.presc_mask, 0)
    presc_y = _component_mask(body_s.presc_mask, 1)

    um_full = np.zeros(box_blk.n_dofs, dtype=np.float64)
    us_full = np.zeros(box_ind.n_dofs, dtype=np.float64)
    slip_state = None
    rows: list[dict] = []
    peak_rss_mb = _rss_mb()
    sys_size = int(body_s._A.Size1())
    n_elem_s = int(len(body_s.mp.Elements))
    free_s = int(body_s.free_mask.sum())

    t0 = time.perf_counter()
    try:
        for t in t_steps:
            ux, uy = H.v08_cut_face_displacement(float(t))
            us_full = us_full.copy()
            us_full[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, float(ux), float(uy), 0.0)
            res = solve_nlgeom_contact_step(
                body_m, body_s, X_s, w_q, B_s, master_face, eps_n,
                um_full, us_full, mu=inp.mu, eps_t=eps_t, slip_state=slip_state,
                outer_tol=1e-8, rtol=1e-8, use_pardiso=True,
            )
            um_full, us_full, slip_state = res["um_full"], res["us_full"], res["slip_state"]

            f_int_s = res["f_int_s"]
            f_c_s = res["f_c_s_full"]
            force_total = res["force_total"]

            R_cut_x = float((f_int_s[presc_x] + f_c_s[presc_x]).sum())
            R_cut_y = float((f_int_s[presc_y] + f_c_s[presc_y]).sum())
            zero_leak_x = float(f_int_s[zero_mask_x].sum())
            fc_zero_x = float(f_c_s[zero_mask_x].sum())

            lhs_id = R_cut_x - float(force_total[0])
            rhs_id = -zero_leak_x
            id_err = abs(lhs_id - rhs_id)

            fn_iface = 2.0 * float(-force_total[1])
            ft_iface = 2.0 * float(force_total[0])
            eq_cut_n = abs(-R_cut_y - float(-force_total[1])) / max(abs(float(force_total[1])), 1e-30)
            eq_cut_t = abs(R_cut_x - float(force_total[0])) / max(abs(float(force_total[0])), 1e-30)

            rss_now = _rss_mb()
            peak_rss_mb = max(peak_rss_mb, rss_now)
            rows.append({
                "t": float(t),
                "ux": float(ux),
                "n_outer": int(res["n_outer"]),
                "it": int(res["iterations"]),
                "n_active": int(res["n_active"]),
                "n_slip": int(res["n_slip"]),
                "F_N": fn_iface,
                "F_T": ft_iface,
                "R_cut_x": R_cut_x,
                "force_total_x": float(force_total[0]),
                "zero_leak_x": zero_leak_x,
                "fc_zero_x": fc_zero_x,
                "lhs_id": lhs_id,
                "rhs_id": rhs_id,
                "id_err": id_err,
                "eq_cut_N": eq_cut_n,
                "eq_cut_T": eq_cut_t,
            })
        elapsed_s = time.perf_counter() - t0
    finally:
        del body_m, body_s, master_face, B_s, X_s, w_q, box_blk, box_ind, um_full, us_full, slip_state
        _release_solver_memory()

    return {
        "ghost_alpha": float(ghost_alpha),
        "ghost_element": ghost_element,
        "n_elem_s": n_elem_s,
        "sys_size": sys_size,
        "free_s": free_s,
        "n_inact_ctrl": int(inact.sum()),
        "elapsed_s": elapsed_s,
        "peak_rss_mb": peak_rss_mb,
        "rows": rows,
    }


def test_v08_ghost_cell_hypothesis_and_zero_alpha_probe():
    """Verify ghost-cell leak identity and probe ghost_alpha=1e-8 vs ghost_alpha=0.0."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    res_tl = _run_ghost_case(
        KM, SM, CL, inp, ghost_alpha=1e-8, ghost_element="TotalLagrangianElement3D4N"
    )
    _print_case_summary("TL_ghost(1e-8)", res_tl)

    res_sd = _run_ghost_case(
        KM, SM, CL, inp, ghost_alpha=1e-8, ghost_element="SmallDisplacementElement3D4N"
    )
    _print_case_summary("SD_ghost(1e-8)", res_sd)

    # Probe ghost_alpha=0.0: verify KratosBodyOperator state when OUT elements are skipped
    nodes_ind, tris_ind, _ = H.hollow_hemisphere_soup(
        inp.R, inp.Ri, inp.center, n=24, sym="half"
    )
    lo_i, hi_i, nc_i = H.v08_indenter_grid(0.25, inp.center, inp.R, sym="half")
    box_ind = BoxSpline(lo_i, hi_i, nc_i, degree=2)
    body_zero = KratosBodyOperator(
        KM, SM, box_ind, inp.E_indenter, inp.nu_indenter,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (2,) if abs(x[2] - lo_i[2]) < 1e-12 else (),
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - hi_i[1]) < 1e-12 else (),
        trim_soup=(nodes_ind, tris_ind),
        q=4,
        ghost_alpha=0.0,
    )
    sys_size_zero = int(body_zero._A.Size1())
    n_elem_zero = int(len(body_zero.mp.Elements))
    n_dofs_ind = int(box_ind.n_dofs)
    _, out_zero, _, _, _ = build_trimmed_table(nodes_ind, tris_ind, box_ind, q=4, cache_dir="off")
    inact_dofs_zero = np.repeat(inactive_ctrl_mask(box_ind, out_zero), 3)
    outside_zero = ~body_zero.in_system_mask
    print(
        f"\n[No_ghost(0.0) probe] n_elem_s={n_elem_zero} (vs {res_sd['n_elem_s']}), "
        f"sys_size={sys_size_zero} vs box.n_dofs={n_dofs_ind}, "
        f"free_s={int(body_zero.free_mask.sum())}, outside_system={int(outside_zero.sum())}",
        flush=True,
    )
    # OUT elements skipped: Kratos drops the inactive DOFs, eval maps back to full size
    assert sys_size_zero < n_dofs_ind
    assert np.array_equal(outside_zero, inact_dofs_zero)
    assert np.all(body_zero.zero_mask[outside_zero])
    K_zero, f_zero = body_zero.eval(np.zeros(n_dofs_ind, dtype=np.float64))
    assert K_zero.shape == (n_dofs_ind, n_dofs_ind) and f_zero.shape == (n_dofs_ind,)
    out_idx = np.flatnonzero(outside_zero)
    assert K_zero[out_idx].nnz == 0 and K_zero[:, out_idx].nnz == 0
    assert np.all(f_zero[outside_zero] == 0.0)
    del body_zero, box_ind, K_zero, f_zero
    _release_solver_memory()

    # 1. Algebraic identity R_cut_x - force_total[0] == -(f_int_s[zero_mask_x].sum())
    for out in (res_tl, res_sd):
        assert out["peak_rss_mb"] < 800.0
        for r in out["rows"]:
            assert abs(r["fc_zero_x"]) < 1e-12
            rel_id = r["id_err"] / max(abs(r["force_total_x"]), 1.0)
            assert rel_id < 1e-7, f"Identity mismatch at t={r['t']}: rel_id={rel_id}"

    # 2. Verify TL ghost (u_x/h)^3 scaling (ux doubles t=0.40 -> 0.70 => ~8x leak)
    tl_040 = next(r for r in res_tl["rows"] if abs(r["t"] - 0.40) < 1e-12)
    tl_070 = next(r for r in res_tl["rows"] if abs(r["t"] - 0.70) < 1e-12)
    ratio_tl_leak = abs(tl_070["zero_leak_x"]) / abs(tl_040["zero_leak_x"])
    assert 7.5 < ratio_tl_leak < 8.5

    # 3. Verify SD ghost linear scaling (~2x leak) and negligible effect on F_N, F_T
    sd_040 = next(r for r in res_sd["rows"] if abs(r["t"] - 0.40) < 1e-12)
    sd_070 = next(r for r in res_sd["rows"] if abs(r["t"] - 0.70) < 1e-12)
    ratio_sd_leak = abs(sd_070["zero_leak_x"]) / abs(sd_040["zero_leak_x"])
    assert 1.9 < ratio_sd_leak < 2.1
    assert sd_070["eq_cut_T"] < 1e-4

    rel_fn_070 = abs(sd_070["F_N"] - tl_070["F_N"]) / sd_070["F_N"]
    rel_ft_070 = abs(sd_070["F_T"] - tl_070["F_T"]) / sd_070["F_T"]
    print(
        f"\n[Summary @ t=0.70] leak_scale=(TL:{ratio_tl_leak:.2f}x, SD:{ratio_sd_leak:.2f}x), "
        f"SD vs TL rel_diff=(F_N:{rel_fn_070:.4%}, F_T:{rel_ft_070:.4%}), "
        f"eq_cut_T=(TL:{tl_070['eq_cut_T']:.3e}, SD:{sd_070['eq_cut_T']:.3e})",
        flush=True,
    )
    assert rel_fn_070 < 5e-3
    assert rel_ft_070 < 5e-3
