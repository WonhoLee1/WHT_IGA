"""S7 loads via run_linear: pressure sign/total + gravity vs analytic (Track SOLVE).

References (analytic, unit cube, E = 1000, nu = 0.3, bottom face strongly fixed):
- PRESSURE p on top (outward +z): bottom reaction sum = -p.A (S08: traction = +p.n_out).
- Gravity (0,0,-g), density rho: tip u_z = -rho.g.L^2/(2E) (uniaxial, exact for any nu).
- SURFACE_LOAD traction t: reaction = -t.A.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.solver.validate import JobSpec, PressureLoad, TractionLoad

pytestmark = pytest.mark.kratos


def _base(name: str) -> JobSpec:
    return JobSpec(name=name, lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(2, 2, 2), degree=1, E=1000.0, nu=0.3,
                   fixed_faces=[(2, 0)])


def test_pressure_sign_and_total() -> None:
    from wht_iga.solver.pipeline import run_linear

    spec = _base("p1_press")
    spec.pressures = [PressureLoad(face=(2, 1), magnitude=2.0)]
    res = run_linear(spec)
    np.testing.assert_allclose(res.total_force, [0.0, 0.0, 2.0], rtol=0.0, atol=1e-12)
    # bottom reaction balances the +z pressure load (S08 sign convention)
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mp = res.keepalive[3]
    fz = sum(nd.GetSolutionStepValue(KM.REACTION_Z) for nd in mp.Nodes if abs(nd.Z0) < 1e-12)
    print(f"\nP1 pressure: total={res.total_force}, bottom reaction z={fz:.12f} (expect -2)")
    np.testing.assert_allclose(fz, -2.0, rtol=0.0, atol=1e-9)


def test_traction_total() -> None:
    from wht_iga.solver.pipeline import run_linear

    spec = _base("p1_trac")
    spec.tractions = [TractionLoad(face=(2, 1), vector=(0.0, 0.0, -3.0))]
    res = run_linear(spec)
    np.testing.assert_allclose(res.total_force, [0.0, 0.0, -3.0], rtol=0.0, atol=1e-12)
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mp = res.keepalive[3]
    fz = sum(nd.GetSolutionStepValue(KM.REACTION_Z) for nd in mp.Nodes if abs(nd.Z0) < 1e-12)
    print(f"\nP1 traction: total={res.total_force}, bottom reaction z={fz:.12f} (expect +3)")
    np.testing.assert_allclose(fz, 3.0, rtol=0.0, atol=1e-9)


def test_gravity_vs_analytic() -> None:
    """Uniaxial gravity column: u_z(top) = -rho.g.L^2/(2E), exact for nu = 0.

    nu = 0 keeps the strong bottom fix consistent (no lateral strain); p = 2
    represents the quadratic gravity profile exactly.
    """
    from wht_iga.solver.pipeline import run_linear

    spec = _base("p1_grav")
    spec.degree = 2
    spec.nu = 0.0
    spec.density = 1.0
    spec.gravity = (0.0, 0.0, -9.81)
    res = run_linear(spec)
    from wht_iga.reference.box import BoxSpline

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    top = box.boundary_control_points(axis=2, side=1)
    uz_top = float(res.U.reshape(-1, 3)[top, 2].mean())
    exact = -1.0 * 9.81 * 1.0 ** 2 / (2 * 1000.0)
    print(f"\nP1 gravity: avg top u_z={uz_top:.12e} (analytic {exact:.12e})")
    np.testing.assert_allclose(uz_top, exact, rtol=1e-6, atol=0.0)


def test_follower_matches_dead_load() -> None:
    from wht_iga.solver.pipeline import run_linear

    spec = _base("p1_fol")
    spec.E = 1.0e6
    spec.pressures = [PressureLoad(face=(2, 1), magnitude=2.0, follower=True)]
    fol = run_linear(spec)
    spec2 = _base("p1_dead")
    spec2.E = 1.0e6
    spec2.pressures = [PressureLoad(face=(2, 1), magnitude=2.0)]
    dead = run_linear(spec2)
    rel = float(np.abs(fol.U - dead.U).max() / max(1e-30, np.abs(dead.U).max()))
    print(f"\nP1 follower vs dead: rel diff={rel:.3e} (tol 1e-5, S08 pattern)")
    assert rel < 1e-5
