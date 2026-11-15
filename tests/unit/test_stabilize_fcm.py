"""T-P2.4 module tests: FCM ghost-alpha law + beta sweep helpers (Kratos-free).

Physics oracle: ADR-002 table (patch-cube sweep: exact 1/β, C ≈ 4–8,
slope −1.00 ± 0.01; V01 defaults p1 6.0e5 / p2 2.4e6 / p3 5.4e6 at h = 0.5).
The O(1/β) law is verified on the exact scalar spring-penalty model
u(β) = (f + βg)/(k + β), err = k|f/k − g|/(k + β) → C/β.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.stabilize import fcm as F


def test_fcm_alpha_default_and_scale():
    assert F.FCM_ALPHA == 1e-8
    assert F.fictitious_modulus(1.0) == pytest.approx(1e-8)
    assert F.fictitious_modulus(2.0e5) == pytest.approx(2.0e-3)
    with pytest.raises(ValueError):
        F.fictitious_modulus(0.0)
    with pytest.raises(ValueError):
        F.fictitious_modulus(1.0, alpha=2.0)


def test_default_beta_matches_adr002_table():
    # ADR-002 §Evidence, h = 0.5, E = 1: p1 6.0e5 / p2 2.4e6 / p3 5.4e6.
    assert F.DEFAULT_C_PEN == 3.0e5
    assert F.default_beta(1.0, 0.5, 1) == pytest.approx(6.0e5)
    assert F.default_beta(1.0, 0.5, 2) == pytest.approx(2.4e6)
    assert F.default_beta(1.0, 0.5, 3) == pytest.approx(5.4e6)
    with pytest.raises(ValueError):
        F.default_beta(1.0, 0.0, 2)


def test_beta_sweep_values():
    sw = F.beta_sweep(1.0, 0.5, 2)
    assert sorted(sw) == [1e2, 1e3, 1e4, 1e5, 1e6]
    for m, b in sw.items():
        assert b == pytest.approx(m * 1.0 / 0.5)


def test_one_over_beta_law_on_exact_model():
    # Scalar penalty model: (k + β) u = f + β g, k=2, f=3, g=0.5.
    k, f, g = 2.0, 3.0, 0.5
    betas = np.array([F.beta_sweep(1.0, 1.0, 1)[m]
                      for m in sorted(F.beta_sweep(1.0, 1.0, 1))])
    u = (f + betas * g) / (k + betas)
    err = np.abs(u - g)
    slope = F.loglog_slope(betas, err)
    # asymptotic 1/β law: β=1e2 decade carries the k=2 offset (ADR-002
    # reports slope −1.00 ± 0.01 on the patch cube for the same reason).
    assert slope == pytest.approx(-1.0, abs=0.01)
    # exact constant C = k|f/k − g| = 2.0 (asymptotic: check upper decades)
    C = err * betas
    np.testing.assert_allclose(C[2:], 2.0, rtol=1e-2)


def test_loglog_slope_rejects_nonpositive():
    with pytest.raises(ValueError):
        F.loglog_slope(np.array([1.0, 0.0]), np.array([1.0, 2.0]))
