"""V01: prescribed-displacement cube via run_linear (Track SOLVE, G1).

Reference: analytic u = (x, 2y, 3z), E = 1, nu = 0, unit cube, penalty BC-A
on all 6 faces (affine field). Criterion: max |u - u_exact| <= 1e-5 for
p = 1, 2, 3 (v2 §9 Phase 1 G1; Coreform level per §10.1 V01).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec

pytestmark = pytest.mark.kratos

AFFINE_A = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))


def _v01_spec(p: int, n: int = 2) -> JobSpec:
    bcs = [DirichletBC(face=(a, s), field=AffineField(A=AFFINE_A, b=(0.0, 0.0, 0.0)))
           for a in range(3) for s in (0, 1)]
    return JobSpec(name=f"v01_p{p}", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(n, n, n), degree=p, E=1.0, nu=0.0, bcs=bcs)


@pytest.mark.parametrize("p", [1, 2, 3])
def test_v01_prescribed_cube(p: int) -> None:
    from wht_iga.solver.pipeline import run_linear

    res = run_linear(_v01_spec(p))
    from wht_iga.reference.box import BoxSpline

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=p)
    Xc = box.control_points()
    Uex = np.column_stack([Xc[:, 0], 2 * Xc[:, 1], 3 * Xc[:, 2]]).ravel()
    rel = float(np.abs(res.U - Uex).max() / np.abs(Uex).max())
    print(f"\nV01 p={p}: beta={res.beta_used:.3e} max-rel-err={rel:.3e} (tol 1e-5)")
    assert rel <= 1e-5


def test_v01_fcm_alpha_one_reproduces_full_cube() -> None:
    """FCM option with alpha=1.0 on OUT cells must reproduce the full-cube solve."""
    from wht_iga.solver.pipeline import run_linear

    base = run_linear(_v01_spec(2))
    spec = _v01_spec(2)
    spec.name = "v01_fcm1"
    spec.out_cells = [(0, 0, 0), (1, 1, 1)]
    spec.ghost_alpha = 1.0
    ghost = run_linear(spec)
    diff = float(np.abs(base.U - ghost.U).max())
    print(f"\nFCM alpha=1: max|U_base - U_ghost| = {diff:.3e} (tol 1e-12)")
    assert diff <= 1e-12


def test_v01_fcm_alpha_small_completes() -> None:
    """FCM ghost alpha=1e-8 (v2 §7.6 baseline): solve completes, U finite."""
    from wht_iga.solver.pipeline import run_linear

    spec = _v01_spec(1)
    spec.name = "v01_fcm8"
    spec.out_cells = [(0, 0, 0)]
    spec.ghost_alpha = 1e-8
    res = run_linear(spec)
    assert np.all(np.isfinite(res.U))
    print(f"\nFCM alpha=1e-8: max|U| = {float(np.abs(res.U).max()):.3e}, finite OK")


def test_v01_bc_b_symmetry_plane() -> None:
    """BC-B per-component: ux=0 on x=0 plane, analytic field elsewhere; <= 1e-5."""
    from wht_iga.solver.pipeline import run_linear

    bcs = [DirichletBC(face=(0, 0), components=(True, False, False), value=(0.0, 0.0, 0.0))]
    bcs += [DirichletBC(face=(a, s), field=AffineField(A=AFFINE_A, b=(0.0, 0.0, 0.0)))
            for (a, s) in [(0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]]
    spec = JobSpec(name="v01_bcb", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(2, 2, 2), degree=2, E=1.0, nu=0.0, bcs=bcs)
    res = run_linear(spec)
    from wht_iga.reference.box import BoxSpline

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    Xc = box.control_points()
    Uex = np.column_stack([Xc[:, 0], 2 * Xc[:, 1], 3 * Xc[:, 2]]).ravel()
    rel = float(np.abs(res.U - Uex).max() / np.abs(Uex).max())
    print(f"\nV01 BC-B: beta={res.beta_used:.3e} max-rel-err={rel:.3e} (tol 1e-5)")
    assert rel <= 1e-5
