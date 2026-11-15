"""P2-STAB Kratos spike K-C: penalty beta sweep on the V01 cube (T-P2.4).

V01 2x2x2 p=2 affine patch (E=1, nu=0) with beta = mult·E/h for
mult in {1e2 … 1e6} (h = 0.5): verifies the O(1/β) law and the plateau
verdict on the SOLVE-owned pipeline (read-only reuse), and checks the
measured constants against the ADR-002 patch-cube values (p=2: C ≈ 3.7).
ADR-002 itself is NOT rewritten: this test extends its evidence with the
V01-case column (phase2 T-P2.4 instruction).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec
from wht_iga.stabilize import fcm as F

pytestmark = pytest.mark.kratos

AFFINE_A = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
REF_PEAK = 4.3175e-6  # V02 reference (kept for the K-D table header only)


def _v01_spec(beta: float) -> JobSpec:
    bcs = [DirichletBC(face=(a, s), field=AffineField(A=AFFINE_A, b=(0.0, 0.0, 0.0)))
           for a in range(3) for s in (0, 1)]
    return JobSpec(name="v01_beta_sweep", lo=(0.0, 0.0, 0.0),
                   hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2,
                   E=1.0, nu=0.0, beta=beta, bcs=bcs)


def test_beta_sweep_v01():
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.pipeline import run_linear

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(2, 2, 2), degree=2)
    Xc = box.control_points()
    Uex = np.column_stack([Xc[:, 0], 2 * Xc[:, 1], 3 * Xc[:, 2]]).ravel()
    betas, errs = [], []
    for mult in sorted(F.SWEEP_MULTS):
        beta = mult * 1.0 / 0.5
        res = run_linear(_v01_spec(beta))
        rel = float(np.abs(res.U - Uex).max() / np.abs(Uex).max())
        betas.append(beta)
        errs.append(rel)
        print(f"K-C: mult={mult:.0e} beta={beta:.3e} max-rel-err={rel:.3e}",
              flush=True)
    slope = F.loglog_slope(np.array(betas), np.array(errs))
    print(f"K-C: O(1/beta) slope={slope:.4f} (acceptance -1 +- 0.1)")
    assert abs(slope + 1.0) <= 0.1
    # Constants C = err_abs·beta vs ADR-002 patch-cube C ≈ 3.7 (p=2).
    # ADR-002 uses the ABSOLUTE max error; the rel-err above divides by
    # max|U| = 3, hence C_rel ≈ 1.24 = 3.72/3 — same physics, both reported.
    abs_errs = np.array(errs) * 3.0
    C = abs_errs * np.array(betas)
    print(f"K-C: C_abs = err_abs·beta = {[f'{c:.3f}' for c in C]} "
          f"(ADR-002 p=2: 3.72)")
    np.testing.assert_allclose(C, 3.72, rtol=0.2)
    # Plateau verdict: no floor in 1e2..1e6 (patch field exactly
    # representable) — same finding as ADR-002, now on the V01 case.
    assert errs[-1] < errs[0] / 100.0
