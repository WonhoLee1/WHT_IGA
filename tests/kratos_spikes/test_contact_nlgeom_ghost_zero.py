# -*- coding: utf-8 -*-
"""``KratosBodyOperator(ghost_alpha=0)``: OUT cells dropped, full-size K / f_int.

With ``ghost_alpha=0`` no element is created on OUT cells, so Kratos assembles
a reduced system without the inactive control points. ``eval`` must still
return ``(n_dofs, n_dofs)`` / ``(n_dofs,)`` arrays with exact zeros on the DOFs
outside the system, all of which are in ``zero_mask``.

References (independent of the reduced-system mapping under test):

(a) The ghost cells use ``SmallDisplacementElement3D4N`` with the linear law at
    ``alpha * E``: their stiffness and force are exactly linear in ``alpha``
    (fixed state ``u``). The full-system path (``alpha > 0``, unchanged code)
    at two values gives the ``alpha = 0`` operator by linear extrapolation,
    ``K(0) = 2 K(a) - K(2a)``, ``f(0) = 2 f(a) - f(2a)``, up to round-off
    (bound 1e-13 of max|K| / max|f|; measured 4e-16). And
    ``|K(alpha) - K(0)| <= alpha * max|G|`` with ``G = (K(2a) - K(a)) / a``:
    the default ``alpha = 1e-8`` operator converges to the ``alpha = 0`` one
    linearly.
(b) A two-step V08 R2 NLGEOM contact run (h = 0.5, t = 0.03 indenting and
    t = 0.15 sliding): the ``alpha = 0`` run must satisfy the cut-face
    identity ``R_cut = F_contact`` to the solver tolerance (no ghost force
    leaks into the inactive DOFs: |R_cut - F_c| / |F_c| < 1e-8 over the x, y
    components, z also carries the symmetry-plane reactions; the inner and
    outer tolerances are 1e-8) and match the ``alpha = 1e-10`` run.
    Tolerance of that match: the ghost perturbation is first order in alpha.
    Spike 2026-11-20 (h = 0.5, same steps): |dF|/|F| = 2.0e-8 / 5.0e-7 and
    surface-displacement rel 4.9e-7 / 3.4e-6 at alpha = 1e-10, scaling by
    89-111x at alpha = 1e-8. Bounds 1e-5 (forces) and 1e-4 (surface field)
    allow 20-30x on top and still reject an O(alpha^0) error. Spline
    coefficients of weakly supported basis functions are NOT compared: they
    move by 1e-3 relative at alpha = 1e-10 while u_h on the surface does not.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import DeformableBoxFaceMaster, sample_surface_triangles_gauss
from wht_iga.analyses.contact_nlgeom import KratosBodyOperator, solve_nlgeom_contact_step
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

pytestmark = [pytest.mark.kratos]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"


def _inputs() -> H.V08Inputs:
    with open(_V08_JSON, encoding="utf-8") as f:
        return H.V08Inputs.from_json(json.load(f))


def _indenter(KM, SM, CL, inp, h: float, ghost_alpha: float):
    """Trimmed V08 half indenter (n_soup = 24) as in test_p5_v08_ghost_check."""
    nodes, tris, fids = H.hollow_hemisphere_soup(inp.R, inp.Ri, inp.center, n=24, sym="half")
    lo, hi, nc = H.v08_indenter_grid(h, inp.center, inp.R, sym="half")
    box = BoxSpline(lo, hi, nc, degree=2)
    op = KratosBodyOperator(
        KM, SM, box, inp.E_indenter, inp.nu_indenter,
        law_factory=CL.HyperElastic3DLaw,
        element="TotalLagrangianElement3D4N",
        zero_fn=lambda x: (2,) if abs(x[2] - lo[2]) < 1e-12 else (),
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - hi[1]) < 1e-12 else (),
        trim_soup=(nodes, tris),
        q=4,
        ghost_alpha=ghost_alpha,
    )
    return op, box, (nodes, tris, fids)


def _smooth_state(box: BoxSpline, scale: float) -> np.ndarray:
    """Deterministic smooth control-point displacement (no random data)."""
    X = np.asarray(box.control_points())
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    U = np.column_stack([np.sin(2.1 * x + 0.7 * y + 0.4 * z),
                         np.cos(1.3 * x - 0.9 * y + 0.5 * z) - 1.0,
                         0.5 * np.sin(1.7 * y + 0.6 * z + 0.2 * x)])
    return scale * U.ravel()


def test_ghost_alpha_zero_eval_equals_alpha_extrapolation(kratos):
    """(a) h = 0.25 indenter (175 control points, 65 inactive), TL Neo-Hookean
    at a nonzero smooth state (max |u| = 2e-3, R = 0.5)."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]
    CL = kratos["ConstitutiveLawsApplication"]
    inp = _inputs()
    a = 1e-2
    op0, box, (nodes, tris, _) = _indenter(KM, SM, CL, inp, 0.25, 0.0)
    op1, _, _ = _indenter(KM, SM, CL, inp, 0.25, a)
    op2, _, _ = _indenter(KM, SM, CL, inp, 0.25, 2.0 * a)
    op8, _, _ = _indenter(KM, SM, CL, inp, 0.25, 1e-8)

    n = box.n_dofs
    _, out_cells, _, _, _ = build_trimmed_table(nodes, tris, box, q=4, cache_dir="off")
    inactive_dofs = np.repeat(inactive_ctrl_mask(box, out_cells), 3)
    assert inactive_dofs.sum() == 3 * 65
    # reduced system: exactly the inactive control points are missing, all in zero_mask
    assert int(op0._A.Size1()) == n - int(inactive_dofs.sum())
    assert np.array_equal(~op0.in_system_mask, inactive_dofs)
    assert np.all(op0.zero_mask[~op0.in_system_mask])
    assert all(int(op._A.Size1()) == n for op in (op1, op2, op8))  # ghost path: full system

    u = _smooth_state(box, 2e-3)
    K0, f0 = op0.eval(u)
    K1, f1 = op1.eval(u)
    K2, f2 = op2.eval(u)
    K8, f8 = op8.eval(u)
    assert K0.shape == (n, n) and f0.shape == (n,)
    miss = np.flatnonzero(~op0.in_system_mask)
    assert K0[miss].nnz == 0 and K0[:, miss].nnz == 0
    assert np.all(f0[miss] == 0.0)

    kmax, fmax = abs(K0).max(), np.abs(f0).max()
    err_K = abs(K0 - (2.0 * K1 - K2)).max() / kmax
    err_f = np.abs(f0 - (2.0 * f1 - f2)).max() / fmax
    G = (K2 - K1) / a
    g = (f2 - f1) / a
    act = np.flatnonzero(op0.in_system_mask)
    g_couple = abs(G[act][:, miss]).max()
    d8_K = abs(K8 - K0).max()
    d8_f = np.abs(f8 - f0).max()
    print(f"\n[ghost_alpha=0 eval] n_dofs={n} n_sys={op0._A.Size1()} max|K0|={kmax:.4e} "
          f"max|f0|={fmax:.4e} extrapolation err K={err_K:.2e} f={err_f:.2e} "
          f"max|G| active-inactive={g_couple:.3e} |K(1e-8)-K0|={d8_K:.3e} "
          f"(bound {1e-8 * abs(G).max():.3e}) |f(1e-8)-f0|={d8_f:.3e} "
          f"(bound {1e-8 * np.abs(g).max():.3e})")
    assert err_K <= 1e-13, err_K
    assert err_f <= 1e-13, err_f
    # the removed ghost part is real: it couples active and inactive DOFs
    assert g_couple > 1e-2 * kmax
    assert np.abs(g[act]).max() > 1e-2 * fmax
    # alpha -> 0: linear convergence, |X(alpha) - X(0)| <= alpha max|G| (+ round-off)
    assert d8_K <= 1e-8 * abs(G).max() * (1.0 + 1e-6) + 1e-13 * kmax
    assert d8_f <= 1e-8 * np.abs(g).max() * (1.0 + 1e-6) + 1e-13 * fmax


def test_ghost_alpha_zero_nlgeom_contact_steps_match_tiny_ghost(kratos):
    """(b) V08 R2 NLGEOM, h = 0.5 half model, t = 0.03 (indenting) and
    t = 0.15 (sliding): ghost_alpha = 0 vs 1e-10."""
    KM = kratos["KratosMultiphysics"]
    SM = kratos["StructuralMechanicsApplication"]
    CL = kratos["ConstitutiveLawsApplication"]
    inp = _inputs()
    h = 0.5
    eps_n, eps_t = 1e10, 1e9
    lo_b, hi_b, nc_b = H.v08_block_grid(h, sym="half")
    box_blk = BoxSpline(lo_b, hi_b, nc_b, degree=2)

    def _blk_zero(x: np.ndarray) -> tuple[int, ...]:
        if abs(x[1] - lo_b[1]) < 1e-12:
            return (0, 1, 2)
        if abs(x[2] - lo_b[2]) < 1e-12:
            return (2,)
        return ()

    body_m = KratosBodyOperator(KM, SM, box_blk, inp.E_block, inp.nu_block,
                                law_factory=CL.HyperElastic3DLaw,
                                element="TotalLagrangianElement3D4N", zero_fn=_blk_zero)
    master_face = DeformableBoxFaceMaster(box_blk, axis=1, side=1)

    runs = {}
    for alpha in (0.0, 1e-10):
        body_s, box_ind, (nodes, tris, fids) = _indenter(KM, SM, CL, inp, h, alpha)
        X_s, w_q, B_s = sample_surface_triangles_gauss(nodes, tris[fids == H.FACE_OUTER], box_ind,
                                                       order=2)
        comp = np.arange(box_ind.n_dofs) % 3
        pm = body_s.presc_mask
        zero_xy = [body_s.zero_mask & (comp == c) for c in (0, 1)]  # = inactive x / y DOFs
        um = np.zeros(box_blk.n_dofs)
        us = np.zeros(box_ind.n_dofs)
        slip = None
        rows = []
        for t in (0.03, 0.15):
            ux, uy = H.v08_cut_face_displacement(t)
            us = us.copy()
            us[pm & (comp == 0)] = float(ux)
            us[pm & (comp == 1)] = float(uy)
            us[pm & (comp == 2)] = 0.0
            res = solve_nlgeom_contact_step(
                body_m, body_s, X_s, w_q, B_s, master_face, eps_n, um, us,
                mu=inp.mu, eps_t=eps_t, slip_state=slip, outer_tol=1e-8, rtol=1e-8,
                use_pardiso=True,
            )
            um, us, slip = res["um_full"], res["us_full"], res["slip_state"]
            fi, fc = res["f_int_s"], res["f_c_s_full"]
            F_c = np.asarray(res["force_total"], dtype=float)
            R_cut = np.array([float((fi + fc)[pm & (comp == c)].sum()) for c in range(3)])
            rows.append({
                "t": t, "F_c": F_c, "R_cut": R_cut, "n_active": int(res["n_active"]),
                # x, y only: z also has the reactions of the z = 0 symmetry plane
                "eq": float(np.linalg.norm(R_cut[:2] - F_c[:2]) / np.linalg.norm(F_c[:2])),
                "leak": [float(fi[m].sum()) for m in zero_xy],
                "surf": B_s @ us.reshape(-1, 3), "um": um.copy(),
            })
            print(f"\n[alpha={alpha:g} t={t}] n_outer={res['n_outer']} n_active={res['n_active']} "
                  f"F_c={F_c} cut-face eq={rows[-1]['eq']:.2e} leak_xy={rows[-1]['leak']}")
        runs[alpha] = rows
        del body_s

    r0, r1 = runs[0.0], runs[1e-10]
    for a0, a1 in zip(r0, r1):
        # alpha = 0: nothing reaches the inactive DOFs, cut face balances the contact force
        assert a0["leak"] == [0.0, 0.0]
        assert a0["eq"] < 1e-8, a0
        assert a1["leak"][0] != 0.0  # the tiny ghost does carry force (it was really dropped)
        assert a0["n_active"] == a1["n_active"] > 0
        dF = float(np.linalg.norm(a0["F_c"] - a1["F_c"]) / np.linalg.norm(a0["F_c"]))
        d_surf = float(np.abs(a0["surf"] - a1["surf"]).max() / np.abs(a0["surf"]).max())
        d_blk = float(np.abs(a0["um"] - a1["um"]).max() / np.abs(a0["um"]).max())
        print(f"[t={a0['t']}] alpha 0 vs 1e-10: |dF|/|F|={dF:.2e} surface u rel={d_surf:.2e} "
              f"block u rel={d_blk:.2e}")
        assert dF < 1e-5, dF
        assert d_surf < 1e-4, d_surf
        assert d_blk < 1e-4, d_blk
    # the sliding step has a real tangential force (eq is not trivially relative to 0)
    assert abs(r0[1]["F_c"][0]) > 0.1 * abs(r0[1]["F_c"][1])
