"""P2-STAB Kratos spike K-D: V02 trimmed runs, FCM-alpha baseline (T-P2.4).

h = 10 level (20x4x1, p=2) via the integration-owned trimmed bridge
(read-only reuse: geometry tessellation + ``trimmed.build_trimmed_table`` +
``pipeline.run_linear``): default beta, plus the sweep endpoints, with the
FCM ghost-alpha baseline (alpha = 1e-8, v2 §7.6).

The extension-vs-FCM comparison table (phase2 T-P2.4) is filled as far as
the current pipeline allows:
* extension-ON V02: BLOCKED — ``run_linear`` has no constraint hook, and
  ``LinearResult`` does not expose the strategy/system matrix, so neither
  the T-constrained solve nor the V02 kappa can run without modifying
  SOLVE-owned files (forbidden). Recorded as blocking issue + shared-module
  change proposal (P24 doc); the extension leg is covered Kratos-free by
  test_stabilize_extension.py and the reduced TKT agreement.
* V02 kappa (either leg): BLOCKED for the same reason.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from wht_iga.solver.validate import JobSpec

pytestmark = pytest.mark.kratos

STEP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                    "data", "roark_plate_hole.step")
LO, HI = (-100.0, -20.0, -0.5), (100.0, 20.0, 0.5)
REF_PEAK = 4.3175e-6


def _run(nc: tuple[int, int, int], beta: float | None,
         ghost_alpha: float = 1e-8) -> dict:
    from wht_iga.geometry import cad as G
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver import bc_roark as _bc
    from wht_iga.solver import trimmed as T
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.validate import JobSpec

    box = BoxSpline(lo=LO, hi=HI, n_cells=nc, degree=2)
    h = float(min((np.asarray(HI) - np.asarray(LO)) / np.asarray(nc)))
    shape = G.read_step(STEP)
    tess = G.tessellate(shape, deflection=float(h / 4))
    qp_table, out_cells, issues, cell_records, diag = T.build_trimmed_table(
        tess.nodes, tess.tris, box, q=4)
    fid = T.hole_face_by_cylinder_wall(tess.nodes, tess.tris,
                                       tess.face_ids, radius=10.0)
    spec = JobSpec(
        name=f"v02_fcm_h{h:g}", lo=LO, hi=HI, n_cells=nc, degree=2,
        E=1.0, nu=0.0, beta=beta, out_cells=out_cells,
        ghost_alpha=ghost_alpha, **_bc.roark_left_pinned())  # type: ignore[arg-type]
    res = run_linear(spec, geo_issues=issues, qp_table=qp_table,
                     cell_records=cell_records)
    vol_peak = T.peak_max_principal(res.cauchy_stress)
    bnd = T.boundary_peak_stress(box, res.U, 1.0, 0.0, tess.nodes, tess.tris,
                                 tess.face_ids, face_id=fid)
    return {"h": h, "peak": bnd["peak"], "vol_peak": vol_peak,
            "rel_err": abs(bnd["peak"] - REF_PEAK) / REF_PEAK,
            "beta": res.beta_used, "n_cut": diag["n_cut"],
            "n_fallback": diag["n_fallback"],
            "umax": float(np.max(np.abs(res.U))),
            "n_out": len(out_cells)}


def test_v02_fcm_baseline_h10():
    row = _run((20, 4, 1), None)
    print(f"K-D h=10 FCM baseline: peak={row['peak']:.6e} "
          f"err={row['rel_err']:.3%} beta={row['beta']:.3e} "
          f"cut={row['n_cut']} fb={row['n_fallback']} umax={row['umax']:.3e} "
          f"out={row['n_out']}", flush=True)
    assert np.isfinite(row["peak"]) and row["peak"] > 0.0
    assert row["n_fallback"] == 0


def test_v02_beta_endpoints_h10():
    rows = [_run((20, 4, 1), 1e2 * 1.0 / 10.0),
            _run((20, 4, 1), 1e6 * 1.0 / 10.0)]
    for row in rows:
        print(f"K-D h=10 beta={row['beta']:.3e}: peak={row['peak']:.6e} "
              f"err={row['rel_err']:.3%}", flush=True)
    # Higher beta must not worsen the boundary peak error here (penalty
    # floor, if any, sits below the cut-cell discretization error).
    assert rows[1]["rel_err"] <= rows[0]["rel_err"] * 2.0
