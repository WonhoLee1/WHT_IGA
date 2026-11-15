"""Beta sweep on the patch cube: O(1/beta) evidence for ADR-002 (Track SOLVE).

Setup: unit cube 2x2x2, p = 2, E = 1, nu = 0, all-face BC-A affine field
u = (x, 2y, 3z). Reference: analytic field (exact in p >= 1, so the sweep
isolates the penalty error; reference-assembler cross-check in TRACK_SOLVE).
Asserts exact 1/beta scaling (slope -1 +- 0.1); prints the table consumed by
ADR-002-penalty-beta.md.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec

pytestmark = pytest.mark.kratos


def test_beta_sweep_inverse_scaling() -> None:
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.penalty import analyze_sweep
    from wht_iga.solver.pipeline import run_linear

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    Xc = box.control_points()
    uex = np.column_stack([Xc[:, 0], 2 * Xc[:, 1], 3 * Xc[:, 2]]).ravel()
    h, E = 0.5, 1.0
    mults = [1e2, 1e3, 1e4, 1e5, 1e6]
    errs = []
    for m in mults:
        beta = m * E / h
        spec = JobSpec(name="sweep", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                       n_cells=(2, 2, 2), degree=2, E=E, nu=0.0, beta=beta,
                       bcs=[DirichletBC(face=(a, s), field=AffineField(
                           A=((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
                           b=(0.0, 0.0, 0.0))) for a in range(3) for s in (0, 1)])
        res = run_linear(spec)
        errs.append(float(np.abs(res.U - uex).max()))
    errs = np.array(errs)
    betas = np.array([m * E / h for m in mults])
    slope = float(np.polyfit(np.log(betas), np.log(errs), 1)[0])
    print("\nP1 sweep p=2 (beta, maxerr): " + ", ".join(f"({b:.1e},{e:.3e})" for b, e in zip(betas, errs)))
    print(f"slope dlog(err)/dlog(beta) = {slope:.4f} (expect -1)")
    assert abs(slope + 1.0) < 0.1
    r = analyze_sweep(betas, errs)
    assert r["min_error"] < 2.5e-6  # mult=1e6 point pins the constant C ~= 3.7
