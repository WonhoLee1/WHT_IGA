# -*- coding: utf-8 -*-
"""V08 Case R2 NLGEOM at h=0.125 over the full v08.json sliding window (0.15 <= t <= 0.70).

Promotes the 2026-11-19 scratchpad run (``v08/nlgeom_r2.py``, TL OUT ghost
cells, 8,018 s: rel_RMS F_N 1.67 %, F_T 0.94 %) to a pytest gate. Since
2026-11-20 ``KratosBodyOperator`` defaults to small-strain ghost cells, which
removed the ghost-force leak into the zero DOFs (cut-face tangential
equilibrium 5.4e-2 -> ~1e-5 at h=0.25), so the equilibrium gate here is the same
1e-4 as for the normal force.

Setup: ``sym='half'``, soup n=48, eps_n=1e11 (penalty converged at h=0.125, S41),
eps_t=1e10, mu=0.3, 6 indentation + 40 sliding steps (same schedule as
``test_v08_r2_nlgeom_full_window_h1``). Runtime about 2-2.5 h, peak RSS < 3 GB.

Gates (v08.json / master plan G5): rel_RMS(F_T) <= 4.5 % and rel_RMS(F_N) <= 4.5 %
over both windows, mean F_T/F_N within 0.015 of mu, cut-face equilibrium
<= 1e-4, n_fallback == 0.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pytest

from tests.kratos_spikes.test_p5_v08_h2_nlgeom import (
    _V08_JSON,
    _cut_face_reaction,
    _prescribed_vector,
    _release_solver_memory,
    _setup_v08_nlgeom_half,
)
from wht_iga.analyses.contact_nlgeom import solve_nlgeom_contact_step
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference import hertz as H

pytestmark = [pytest.mark.kratos, pytest.mark.slow]


def test_v08_r2_nlgeom_full_window_h2():
    """V08 R2 NLGEOM, h=0.125, n=48, eps_n=1e11, eps_t=1e10, t = 0.03 -> 0.70."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    h, n_soup, eps_n, eps_t = 0.125, 48, 1e11, 1e10
    t_all = np.concatenate([np.linspace(0.03, 0.1, 6), np.linspace(0.115, 0.70, 40)])

    t0 = time.perf_counter()
    _bb, box_ind, body_m, body_s, X_s, w_q, B_s, master = _setup_v08_nlgeom_half(KM, SM, CL, inp, h, n_soup=n_soup)
    assert body_s.diag["n_fallback"] == 0
    print(f"\n[R2 h=0.125 setup] ghost cells, n_gp={len(w_q)}, {time.perf_counter() - t0:.1f}s", flush=True)

    um = np.zeros(body_m.box.n_dofs)
    us = np.zeros(box_ind.n_dofs)
    slip = None
    rows = []
    try:
        for step, t in enumerate(t_all, start=1):
            ux, uy = H.v08_cut_face_displacement(float(t))
            us = us.copy()
            us[body_s.presc_mask] = _prescribed_vector(body_s.presc_mask, float(ux), float(uy), 0.0)
            res = solve_nlgeom_contact_step(body_m, body_s, X_s, w_q, B_s, master, eps_n, um, us,
                                            mu=inp.mu, eps_t=eps_t, slip_state=slip,
                                            outer_tol=1e-8, rtol=1e-8, use_pardiso=True)
            um, us, slip = res["um_full"], res["us_full"], res["slip_state"]
            fn, ft = 2.0 * float(-res["force_total"][1]), 2.0 * float(res["force_total"][0])
            Rc = _cut_face_reaction(body_s, res)
            rows.append((float(t), fn, ft, 2.0 * float(-Rc[1]), 2.0 * float(Rc[0])))
            if step in (1, 6) or step % 5 == 0 or step == len(t_all):
                print(f"  step {step:2d}/{len(t_all)} t={t:.4f} outer={res['n_outer']} n_act={res['n_active']} "
                      f"F_N={fn:,.1f} F_T={ft:,.1f} ratio={ft / fn if fn else 0.0:.4f}", flush=True)
    finally:
        _release_solver_memory()

    A = np.array(rows)
    for lo, hi in ((0.15, 0.40), (0.15, 0.70)):
        w = (A[:, 0] >= lo) & (A[:, 0] <= hi + 1e-12)
        ref = H.v08_reference_history(A[w, 0], inp)
        rms_fn = H.relative_rms(A[w, 1], ref["F_N"])
        rms_ft = H.relative_rms(A[w, 2], ref["F_T"])
        eq_n = float(np.max(np.abs(A[w, 3] - A[w, 1]) / A[w, 1]))
        eq_t = float(np.max(np.abs(A[w, 4] - A[w, 2]) / A[w, 2]))
        ratio = A[w, 2] / A[w, 1]
        print(f"[R2 h=0.125 window {lo}<=t<={hi}] n={int(w.sum())} rel_RMS=(F_N {rms_fn * 100:.2f}%, "
              f"F_T {rms_ft * 100:.2f}%) ratio mean={ratio.mean():.4f} [{ratio.min():.4f}, {ratio.max():.4f}] "
              f"eq_cut=(N {eq_n:.2e}, T {eq_t:.2e})", flush=True)
        assert rms_ft <= 0.045, f"[{lo},{hi}] F_T rel_RMS {rms_ft:.4f}"
        assert rms_fn <= 0.045, f"[{lo},{hi}] F_N rel_RMS {rms_fn:.4f}"
        assert abs(float(ratio.mean()) - inp.mu) < 0.015
        assert eq_n < 1e-4 and eq_t < 1e-4, f"[{lo},{hi}] eq_cut=({eq_n:.2e}, {eq_t:.2e})"
    print(f"[R2 h=0.125 total] {time.perf_counter() - t0:.0f}s", flush=True)
