"""T-P2.1 acceptance: support classification vs S3/S4-style inputs (Kratos-free).

References: untrimmed grid (analytic: all-master), synthetic hole grid with
fractional CUT ring (band match), theta-sweep monotonicity, multi-material
union branch (v2.1 §7.5, pitfalls K20/K23).
"""

from __future__ import annotations

import numpy as np

from wht_iga.grid.knots import uniform_clamped
from wht_iga.stabilize import support as S


def _grid(nc, p, lo=0.0, hi=8.0):
    knots = [uniform_clamped(lo, hi, n, p) for n in nc]
    ncp = tuple(len(k) - p - 1 for k in knots)
    br = [np.unique(k) for k in knots]
    d = [np.diff(b) for b in br]
    cv = (d[0][:, None, None] * d[1][None, :, None]
          * d[2][None, None, :]).ravel()
    return knots, ncp, cv


def _ranges_1d(knots, br, p, ncp):
    """Independent reimplementation of support ranges (test oracle)."""
    e0, e1 = np.empty(ncp, dtype=int), np.empty(ncp, dtype=int)
    for i in range(ncp):
        a = int(np.searchsorted(br, knots[i], side="right") - 1)
        b = int(np.searchsorted(br, knots[i + p + 1], side="left") - 1)
        e0[i], e1[i] = max(a, 0), min(b, len(br) - 2)
    return e0, e1


def test_untrimmed_all_master():
    knots, ncp, cv = _grid((4, 3, 2), 2)
    out = S.classify_basis(ncp, knots, cv.copy(), cv, theta=0.2)
    assert out.n_master == np.prod(ncp)
    assert out.n_inactive == 0 and out.n_extended == 0
    np.testing.assert_allclose(out.ratios, 1.0, rtol=0, atol=0)
    assert out.multi_material is False
    assert any("not-run" in n for n in out.notes)


def _hole_case():
    nc, p = (8, 8, 2), 2
    knots, ncp, cv = _grid(nc, p)
    nx, ny, nz = nc
    eta = np.ones(nx * ny * nz)
    hole, ring1, ring2 = [], [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                d = (i - 3.5) ** 2 + (j - 3.5) ** 2
                c = i + nx * (j + ny * k)
                if d < 2.5 ** 2:
                    eta[c] = 0.0
                    hole.append(c)
                elif d < 3.5 ** 2:
                    eta[c] = 0.1
                    ring1.append(c)
                elif d < 4.5 ** 2:
                    eta[c] = 0.3
                    ring2.append(c)
    return knots, ncp, cv, eta, hole, ring1, ring2


def test_hole_extended_band_matches_cut_boundary():
    knots, ncp, cv, eta, hole, ring1, _ring2 = _hole_case()
    nx, ny, _nz = (8, 8, 2)
    out = S.classify_basis(ncp, knots, eta * cv, cv, theta=0.2)
    assert out.n_extended > 0 and out.n_inactive > 0 and out.n_master > 0
    r0 = [_ranges_1d(knots[ax], np.unique(knots[ax]), 2, ncp[ax])
          for ax in range(3)]
    small = set(np.flatnonzero((eta > 0.0) & (eta < 0.2)).tolist())
    assert len(small) > 0
    nu, nv, _nw = ncp
    covered = set()
    for b in np.flatnonzero(out.labels == S.EXTENDED):
        a = (int(b % nu), int((b // nu) % nv), int(b // (nu * nv)))
        cells = {(i + nx * (j + ny * k))
                 for k in range(r0[2][0][a[2]], r0[2][1][a[2]] + 1)
                 for j in range(r0[1][0][a[1]], r0[1][1][a[1]] + 1)
                 for i in range(r0[0][0][a[0]], r0[0][1][a[0]] + 1)}
        assert cells & small, f"extended basis {a} touches no small-cut cell"
        covered |= (cells & small)
    # every small-cut cell is seen by a non-inactive support
    assert small <= covered, f"uncovered small cells: {sorted(small - covered)}"
    # inactive supports sit fully inside the zero region
    zeros = set(np.flatnonzero(eta == 0.0).tolist())
    for b in np.flatnonzero(out.labels == S.INACTIVE):
        a = (int(b % nu), int((b // nu) % nv), int(b // (nu * nv)))
        cells = {(i + nx * (j + ny * k))
                 for k in range(r0[2][0][a[2]], r0[2][1][a[2]] + 1)
                 for j in range(r0[1][0][a[1]], r0[1][1][a[1]] + 1)
                 for i in range(r0[0][0][a[0]], r0[0][1][a[0]] + 1)}
        assert cells <= zeros, f"inactive basis {a} meets Omega"


def test_theta_sweep_monotone():
    knots, ncp, cv, eta, *_ = _hole_case()
    counts = []
    prev_slaves: set[int] | None = None
    prev_masters: set[int] | None = None
    for th in (0.05, 0.2, 0.35):
        out = S.classify_basis(ncp, knots, eta * cv, cv, theta=th)
        counts.append(out.n_extended)
        slaves = set(np.flatnonzero(out.labels == S.EXTENDED).tolist())
        masters = set(np.flatnonzero(out.labels == S.MASTER).tolist())
        if prev_slaves is not None:
            assert prev_slaves <= slaves, f"slave set not nested at {th}"
            assert masters <= prev_masters, f"master set not nested at {th}"
        prev_slaves, prev_masters = slaves, masters
    assert counts[0] <= counts[1] <= counts[2]
    assert counts[0] < counts[2], counts  # sweep actually exercises both rings


def test_multi_material_union_matches_single():
    knots, ncp, cv, eta, *_ = _hole_case()
    nx = 8
    idx = np.arange(len(eta))
    i = idx % nx
    eta_a = np.where(i < 4, eta, 0.0)
    eta_b = np.where(i >= 4, eta, 0.0)
    single = S.classify_basis(ncp, knots, eta * cv, cv, theta=0.2)
    union = S.classify_basis(ncp, knots, eta * cv, cv, theta=0.2,
                             material_etas=[eta_a, eta_b])
    assert union.multi_material is True
    np.testing.assert_array_equal(union.labels, single.labels)
    np.testing.assert_allclose(union.ratios, single.ratios, rtol=0, atol=1e-15)


def test_multi_material_overlap_rejected():
    knots, ncp, cv, eta, *_ = _hole_case()
    with np.testing.assert_raises(ValueError):
        S.classify_basis(ncp, knots, eta * cv, cv, theta=0.2,
                         material_etas=[eta, eta])


def test_theta_range_enforced():
    knots, ncp, cv = _grid((2, 2, 2), 1)
    for bad in (-0.1, 0.5, 0.7):
        try:
            S.classify_basis(ncp, knots, cv.copy(), cv, theta=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"theta={bad} accepted")
