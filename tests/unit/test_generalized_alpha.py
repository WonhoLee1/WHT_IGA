"""Generalized-alpha / HHT linear integrator (Kratos coefficient convention), Kratos-free."""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.analyses.modal_dynamics import (
    bossak_amplification_matrix,
    generalized_alpha_amplification_matrix,
    generalized_alpha_free_vibration_theory,
)
from wht_iga.dynamics import (
    HHT_ABAQUS_PRESETS,
    BossakConfig,
    GeneralizedAlphaConfig,
    solve_bossak_linear,
    solve_generalized_alpha_linear,
)


def test_abaqus_presets_match_documented_coefficients() -> None:
    # Abaqus *DYNAMIC application defaults (HHT-TF / HHT-MD)
    tf = GeneralizedAlphaConfig.abaqus("TRANSIENT_FIDELITY", 1e-3, 1)
    md = GeneralizedAlphaConfig.abaqus("MODERATE_DISSIPATION", 1e-3, 1)
    assert (tf.alpha_m, tf.alpha_f) == (0.0, 0.05)
    assert tf.beta == pytest.approx(0.275625, abs=1e-12) and tf.gamma == pytest.approx(0.55, abs=1e-12)
    assert md.beta == pytest.approx(0.5, abs=1e-5) and md.gamma == pytest.approx(0.91421, abs=1e-12)
    assert set(HHT_ABAQUS_PRESETS) == {"TRANSIENT_FIDELITY", "MODERATE_DISSIPATION"}


@pytest.mark.parametrize("bad", [dict(alpha_m=0.1, alpha_f=0.0), dict(alpha_m=0.0, alpha_f=0.6)])
def test_config_rejects_unstable_parameters(bad) -> None:
    with pytest.raises(ValueError):
        GeneralizedAlphaConfig(dt=1e-3, n_steps=1, **bad)
    with pytest.raises(ValueError):
        GeneralizedAlphaConfig.hht(1e-3, 1, 0.1)


def _chain(n: int = 6):
    k = sp.diags([2 * np.ones(n), -np.ones(n - 1), -np.ones(n - 1)], [0, 1, -1]).tocsr() * 1e4
    m = sp.diags(np.linspace(1.0, 2.0, n)).tocsr()
    f = np.zeros(n); f[-1] = 1.0
    return m, k, lambda t: f * math.sin(300.0 * t)


@pytest.mark.parametrize("alpha_b", [-0.05, -0.3])
def test_alpha_f_zero_reproduces_bossak(alpha_b: float) -> None:
    m, k, f = _chain()
    b = solve_bossak_linear(m, None, k, f, BossakConfig(dt=1e-3, n_steps=200, alpha_b=alpha_b))
    g = solve_generalized_alpha_linear(m, None, k, f, GeneralizedAlphaConfig(dt=1e-3, n_steps=200, alpha_m=alpha_b))
    assert np.abs(g.u - b.u).max() <= 1e-12 * np.abs(b.u).max()
    assert np.abs(bossak_amplification_matrix(alpha_b, 50.0, 1e-2)
                  - generalized_alpha_amplification_matrix(alpha_b, 0.0, 50.0, 1e-2)).max() < 1e-12


@pytest.mark.parametrize("alpha_hht", [-0.05, -0.41421])
def test_hht_sdof_free_vibration_matches_amplification_theory(alpha_hht: float) -> None:
    """Single oscillator: integrated history equals the amplification-matrix recursion,
    and its spectral radius gives the measured envelope decay."""
    w, dt, n = 2 * math.pi * 10.0, 1e-3, 400
    cfg = GeneralizedAlphaConfig.hht(dt, n, alpha_hht)
    res = solve_generalized_alpha_linear(np.array([[1.0]]), None, np.array([[w * w]]), np.zeros(1), cfg,
                                         u0=np.array([1.0]), v0=np.array([0.0]))
    amp = generalized_alpha_amplification_matrix(cfg.alpha_m, cfg.alpha_f, w, dt)
    state = np.array([1.0, 0.0, -w * w])
    for _ in range(n):
        state = amp @ state
    assert res.u[-1, 0] == pytest.approx(state[0], rel=1e-10, abs=1e-12)
    f_th, s_th = generalized_alpha_free_vibration_theory(cfg.alpha_m, cfg.alpha_f, w, dt)
    assert f_th < w / (2 * math.pi)  # period elongation
    assert s_th > 0.0  # numerical dissipation for alpha_hht < 0
    rho = np.abs(np.linalg.eigvals(amp)).max()
    assert rho <= 1.0 + 1e-12  # unconditionally stable range
