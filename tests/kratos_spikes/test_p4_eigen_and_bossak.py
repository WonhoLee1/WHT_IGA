"""P4-A: Kratos lumped mass, Kratos Bossak versus Python Bossak, V10 rod versus Bossak theory.

Oracles
-------
* Kratos lumped mass == row sums of the Python consistent mass (partition of unity).
* Kratos ``ResidualBasedBossakDisplacementScheme`` (C++) == ``solve_bossak_linear``
  (scipy) for the same lumped mass: two independent implementations of one scheme.
* The modal coordinate of the axial mode of the V10 rod must decay and oscillate
  exactly as the eigenvalues of the one-step Bossak amplification matrix predict.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.analyses.modal_dynamics import bossak_free_vibration_theory, modal_coordinate
from wht_iga.dynamics import BossakConfig, solve_bossak_linear, solve_eigenmodes
from wht_iga.kratos_backend.dynamics import (
    face_traction_force,
    fit_free_vibration,
    kratos_lumped_mass,
    run_kratos_bossak,
    smooth_pulse,
)
from wht_iga.kratos_backend.mass import assemble_box_mass, face_dofs, kratos_box_stiffness
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

E, RHO = 1.0e6, 1.0e-6
PRESETS = {"TRANSIENT_FIDELITY": -0.05, "MODERATE_DISSIPATION": -0.30}


@pytest.mark.parametrize("p", [1, 2])
def test_kratos_lumped_mass_equals_consistent_row_sums(p: int) -> None:
    """Row sums of the Python consistent mass equal Kratos' lumped mass vector entry by entry."""
    box = BoxSpline(lo=(0, 0, 0), hi=(4.0, 1.0, 2.0), n_cells=(4, 1, 2), degree=p)
    lump = kratos_lumped_mass(box, RHO)
    rowsum = np.asarray(assemble_box_mass(box, RHO).sum(axis=1)).ravel()
    assert np.abs(lump - rowsum).max() < 1e-13 * np.abs(rowsum).max()


@pytest.mark.parametrize("name", list(PRESETS))
def test_kratos_bossak_matches_python_bossak_with_lumped_mass(name: str) -> None:
    """Free-end history from the Kratos scheme equals the python scheme on the same lumped M."""
    alpha = PRESETS[name]
    length = 4.0
    box = BoxSpline(lo=(0, 0, 0), hi=(length, 1.0, 1.0), n_cells=(4, 1, 1), degree=2)
    clamp = face_dofs(box, 0, 0)
    f0 = face_traction_force(box, 0, 1, (1.0, 0.0, 0.0))
    t1 = 4.0 * length / np.sqrt(E / RHO)
    amp = smooth_pulse(t1)
    dt, n = t1 / 20.0, 120
    probe = np.array([3 * int(c) for c in box.boundary_control_points(0, 1)])
    kr = run_kratos_bossak(box, E, 0.0, RHO, clamp, f0, amp, alpha, dt, n, probe)
    K = kratos_box_stiffness(box, E, 0.0, clamp)
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    py = solve_bossak_linear(
        sp.diags(kratos_lumped_mass(box, RHO)), None, K, lambda t: amp(t) * f0,
        BossakConfig(dt=dt, n_steps=n, alpha_b=alpha), free_dofs=free,
    )
    u_py = py.u[:, probe]
    scale = np.abs(u_py).max()
    assert scale > 0.0
    assert np.abs(kr.u - u_py).max() < 1e-10 * scale
    assert np.abs(kr.v - py.v[:, probe]).max() < 1e-10 * np.abs(py.v[:, probe]).max()
    assert np.abs(kr.a - py.a[:, probe]).max() < 1e-9 * np.abs(py.a[:, probe]).max()


@pytest.mark.parametrize("name", list(PRESETS))
def test_v10_rod_axial_mode_follows_bossak_amplification_theory(name: str) -> None:
    """V10 rod (10 x 2 x 2, p = 2): axial mode is 25 kHz; its free vibration obeys the Bossak eigenvalues."""
    alpha = PRESETS[name]
    length = 10.0
    box = BoxSpline(lo=(0, 0, 0), hi=(length, 2.0, 2.0), n_cells=(10, 2, 2), degree=2)
    clamp = face_dofs(box, 0, 0)
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    K = kratos_box_stiffness(box, E, 0.0, clamp)
    M = assemble_box_mass(box, RHO)
    eig = solve_eigenmodes(K, M, n_modes=12, free_dofs=free)
    ax = int(np.argmax(np.abs(eig.participation_factors[:, 0])))
    f_exact = np.sqrt(E / RHO) / (4.0 * length)
    assert eig.frequencies_hz[ax] == pytest.approx(f_exact, rel=1e-6)  # discrete axial mode, p=2, h=1

    t1 = 1.0 / f_exact
    amp = smooth_pulse(t1)
    f0 = face_traction_force(box, 0, 1, (1.0, 0.0, 0.0))
    dt, n = 2.0e-6, 1000
    res = solve_bossak_linear(M, None, K, lambda t: amp(t) * f0, BossakConfig(dt=dt, n_steps=n, alpha_b=alpha), free_dofs=free)
    q = modal_coordinate(eig.modes[:, ax], M, res.u)
    sel = res.t >= 2.0 * t1  # load fully removed
    fit = fit_free_vibration(res.t[sel], q[sel], f_exact)
    f_th, s_th = bossak_free_vibration_theory(alpha, eig.omega[ax], dt)
    assert fit.frequency_hz == pytest.approx(f_th, rel=1e-5)
    assert fit.decay_rate == pytest.approx(s_th, rel=1e-2)
    assert fit.residual_rms < 1e-3  # single mode: a pure damped sinusoid
    # the more dissipative preset must decay faster and lower the measured frequency
    other = bossak_free_vibration_theory(-0.30 if alpha > -0.1 else -0.05, eig.omega[ax], dt)
    assert (s_th < other[1]) == (alpha > -0.1)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
