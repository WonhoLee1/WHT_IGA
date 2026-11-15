"""Line-probe quadrature and self-convergence rate (Coreform V02 convention)."""

from __future__ import annotations

import numpy as np

from wht_iga.solver import probes as PR


def test_line_gauss_integrates_polynomial():
    P, w = PR.line_gauss_points((0.0, 10.0, 0.0), (0.0, 20.0, 0.0), n_seg=4, n_gauss=4)
    y = P[:, 1]
    assert abs(w.sum() - 10.0) < 1e-12
    assert abs(np.sum(w * y**3) - (20**4 - 10**4) / 4) < 1e-8
    assert np.all((y > 10) & (y < 20))


def test_self_convergence_recovers_known_rate():
    P, w = PR.line_gauss_points((0.0, 10.0, 0.0), (0.0, 20.0, 0.0))
    y = P[:, 1]
    exact = 1.0 + 1.0 / y**2
    hs = np.array([6.0, 3.0, 1.5, 0.75, 0.375])
    profiles = [exact + 0.3 * h**2 * np.sin(y) for h in hs]
    sc = PR.self_convergence(hs, profiles, w)
    # errors vs finest: c*(h^2 - h_f^2) -> slope -> 2 as h >> h_f
    assert len(sc["errors"]) == 4
    assert 1.9 < sc["rate_all"] < 2.3
    assert np.all(np.diff(sc["errors"]) < 0)


def test_rel_l2_zero_for_identical():
    w = np.ones(5)
    a = np.arange(1.0, 6.0)
    assert PR.rel_l2(a, a, w) == 0.0


def test_peak_max_principal_fast_matches_loop():
    from wht_iga.solver.trimmed import peak_max_principal

    rng = np.random.default_rng(3)
    data = [rng.standard_normal((n, 6)) for n in (1, 27, 8, 125)]
    assert PR.peak_max_principal_fast(data) == peak_max_principal(data)
