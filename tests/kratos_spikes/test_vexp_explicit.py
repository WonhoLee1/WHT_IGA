# -*- coding: utf-8 -*-
"""V-EXP / V-EXP-C verification (Phase 6 gates, master-plan v2.1 section 10.1).

V-EXP: 1D elastic wave, step end load sigma0 at x = L, x = 0 fixed.
Model data come from tests/verification/data/vexp.json: p = 1..3, nx = 40,
dt = 0.005 (CFL 0.2), t_end = 2.0. The analytic (d'Alembert) solution has
front speed c = 10 and particle velocity v0 = sigma0 c / E = 0.1.
Plan gates: arrival time and peak stress within 2 % (vexp.json still says
0.03 s for arrival; the plan value is used).
- Tip history: u(L, t) = v0 t for t <= 2L/c. Gate: relative L2.
- Peak stress: 2 sigma0 at the fixed end x = 0 for L/c < t < 3L/c
  (reflection doubling, the largest stress in the bar). Field value
  E sum N_i'(0) u_i, averaged over t in [1.5, 2.0] s.
- Arrival at x = 5.0 (field value sum N_i(5) u_i), analytic 0.5 s.
  - First moment of the front: gate 2 % at nx = 40.
  - Threshold 0.5 v0 (vexp.json definition): misses 2 % at nx = 40 because of
    numerical dispersion (gap pin + Airy theory + h^(2/3) convergence below).
    The 2 % gate is met on the Airy-predicted refined mesh.
- Independent reference for the discrete behaviour: a 1D SciPy B-spline bar
  (no Kratos, no wht_iga) with the same lumped mass and leapfrog start.
V-EXP-C (free-free two-block elastic collision, direct Python penalty on
  surface layers): p = 1 and p = 2, momentum conserved (< 1e-6), rebound
  velocity exchange, total energy (KE + PE_strain + E_penalty) drift (< 2%),
  bounded penetration, dt sweep {0.0025, 0.005}.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import BSpline
from scipy.optimize import brentq
from scipy.special import itairy

from tests.kratos_spikes import test_s24b_hybrid_explicit as S24B
from tests.kratos_spikes.test_s24b_hybrid_explicit import (
    _apply_tip_traction,
    _build_wave_model,
)
from tests.kratos_spikes.test_s25_explicit_contact import (
    _assemble_block_stiffness,
    _build_two_blocks,
)
from wht_iga.dynamics import explicit as EX
from wht_iga.grid.bspline import basis_1d
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_DATA = json.loads((Path(__file__).resolve().parents[1] / "verification" / "data"
                    / "vexp.json").read_text(encoding="utf-8"))
_L = float(_DATA["model"]["L"])
_NX = int(_DATA["model"]["n_cells"][0])
_E = float(_DATA["material"]["E"])
_NU = float(_DATA["material"]["nu"])
_RHO = float(_DATA["material"]["rho"])
_SIGMA0 = float(_DATA["loads"]["traction"]["magnitude"])
_DT = float(_DATA["transient"]["dt"][0])
_T_END = float(_DATA["transient"]["t_end"])
_C = float(np.sqrt(_E / _RHO))
_V0 = _SIGMA0 * _C / _E
_CFL = _DT * _C / (_L / _NX)
_X_PROBE = 0.5 * _L
_T_ARRIVAL = float(_DATA["reference_values"]["analytic_wavefront_arrival_s"])
_T_FRONT = 0.6  # run length for arrival-only histories (0.5 v0 crossing at x = 5 by 0.53 s)
_GATE = 0.02  # master-plan v2.1 section 10.1: arrival time and peak stress within 2 %
# Airy leading-order lag: relative remainder O((h/d)^(2/3)) = 0.14 at nx = 40 with
# small constants; 10 % of the lag is the allowed asymptotic remainder.
_AIRY_TOL = 0.10
_ROUNDOFF_S = 1e-9  # Kratos vs the 1D reference: same discrete system, round-off only


def test_vexp_reference_data_matches_model_builder():
    """vexp.json values equal the constants the shared S24b builder actually uses."""
    assert (S24B._E, S24B._NU, S24B._RHO, S24B._SIGMA0) == (_E, _NU, _RHO, _SIGMA0)
    assert float(_DATA["model"]["length"]) == _L
    assert abs(_T_ARRIVAL - (_L - _X_PROBE) / _C) < 1e-12
    assert abs(float(_DATA["reference_values"]["sigma0_times_c_over_E"]) - _V0) < 1e-12
    assert abs(_T_END - 2.0 * _L / _C) < 1e-12  # tip formula valid up to t_end
    assert abs(_CFL - 0.2) < 1e-12  # vexp.json "cfl": 0.2


@functools.lru_cache(maxsize=None)
def _kratos_wave(degree: int, nx: int, dt: float, t_end: float) -> dict:
    """One Kratos run. Per-step x-field probes (rows averaged):
    tip control point (= u(L)), u(x = 5.0), sigma(x = 0) = E u'(0)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(nx, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, 1.0, False, f"VEXP_p{degree}_nx{nx}")
    mp = mh["mp"]
    _apply_tip_traction(KM, SM, mp, box, nx, 1.0, False, degree)
    nu = box.n_cp[0]
    n_rows = box.n_ctrl // nu
    first, N, dN, X, _ = box.axis_tables(0, np.array([-1.0]))  # left end of each cell
    c_probe = int(round(_X_PROBE / (_L / nx)))
    assert abs(X[c_probe, 0] - _X_PROBE) < 1e-12 and abs(X[0, 0]) < 1e-12
    w_probe = np.zeros(nu)
    w_probe[first[c_probe]:first[c_probe] + degree + 1] = N[c_probe, 0, :]
    w_dx0 = np.zeros(nu)
    w_dx0[first[0]:first[0] + degree + 1] = dN[0, 0, :]
    nodes = [mp.GetNode(I + 1) for I in range(box.n_ctrl)]
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    rec = {"tip": [], "u_probe": [], "sigma_x0": []}

    def probe(_mp):
        ux = np.array([nd.GetSolutionStepValue(KM.DISPLACEMENT_X) for nd in nodes])
        cx = ux.reshape(n_rows, nu).mean(axis=0)  # x-coefficients, averaged over (y, z) rows
        rec["tip"].append(float(cx[-1]))
        rec["u_probe"].append(float(w_probe @ cx))
        rec["sigma_x0"].append(_E * float(w_dx0 @ cx))
        return rec["tip"][-1]

    n_steps = int(round(t_end / dt))
    out = EX.run_explicit_history(KM, mp, strat, dt, n_steps, probe_fn=probe)
    assert len(out["history"]) == n_steps
    hist = {k: np.array(v) for k, v in rec.items()}
    hist["t"] = (np.arange(n_steps) + 1) * dt
    return hist


def _wave_histories(degree: int) -> dict:
    """The vexp.json run (nx = 40, dt = 0.005, t_end = 2.0)."""
    return _kratos_wave(degree, _NX, _DT, _T_END)


@functools.lru_cache(maxsize=None)
def _ref1d_operators(degree: int, nx: int):
    """Independent 1D discretisation of the bar (SciPy B-splines; no Kratos, no wht_iga).

    Open uniform knots on [0, L], unit cross-section, K by (p + 1)-point Gauss
    (exact), row-sum lumped mass m_i = rho int N_i. The 3D Kratos bar (one cell
    across, nu = 0, y/z fixed, uniform face traction) reduces to this exactly.
    """
    p, h = degree, _L / nx
    knots = np.concatenate([np.zeros(p), np.linspace(0.0, _L, nx + 1), np.full(p, _L)])
    n = len(knots) - p - 1
    basis = BSpline(knots, np.eye(n), p)
    g, w = np.polynomial.legendre.leggauss(p + 1)
    xq = ((np.arange(nx)[:, None] + 0.5 + 0.5 * g[None, :]) * h).ravel()
    wq = np.tile(0.5 * h * w, nx)
    Nq, dNq = basis(xq), basis.derivative()(xq)
    K = _E * (dNq.T * wq) @ dNq
    m = _RHO * (Nq.T @ wq)
    return basis, K, m


@functools.lru_cache(maxsize=None)
def _ref1d_history(degree: int, nx: int, dt: float, t_end: float) -> dict:
    """Leapfrog from rest, half-step start v_{1/2} = (dt/2) a_0, step end load sigma0."""
    basis, K, m = _ref1d_operators(degree, nx)
    Kf, mf = K[1:, 1:], m[1:]  # x = 0 fixed (first control point)
    f = np.zeros(len(mf))
    f[-1] = _SIGMA0  # N_n(L) = 1
    w_probe = basis(np.array([_X_PROBE]))[0, 1:]
    w_dx0 = basis.derivative()(np.array([0.0]))[0, 1:]
    u = np.zeros(len(mf))
    v = 0.5 * dt * f / mf
    n_steps = int(round(t_end / dt))
    up, s0 = np.empty(n_steps), np.empty(n_steps)
    for i in range(n_steps):
        u = u + dt * v
        up[i], s0[i] = w_probe @ u, _E * (w_dx0 @ u)
        v = v + dt * (f - Kf @ u) / mf
    return {"t": (np.arange(n_steps) + 1) * dt, "u_probe": up, "sigma_x0": s0}


def _beta(degree: int, cfl: float) -> float:
    """Leading dispersion coefficient: omega = c k (1 - beta (k h)^2) for lumped-mass leapfrog.

    Interior symbol: omega_h^2 = (c/h)^2 (2 - 2 cos th) S_{p-1}(th), S_0 = 1 and
    S_q = 1 - (q + 1) th^2 / 12 + ... (Gram symbol of degree-q B-splines). So
    omega_h = (c/h) th (1 - b_p th^2) with b_1 = 1/24 and b_p = (p + 1)/24 for p >= 2.
    Leapfrog sin(omega dt / 2) = omega_h dt / 2 subtracts CFL^2 / 24.
    """
    b_p = 1.0 / 24.0 if degree == 1 else (degree + 1.0) / 24.0
    return b_p - cfl * cfl / 24.0


def _airy_lag(degree: int, nx: int, cfl: float = _CFL) -> float:
    """Relative delay of the 0.5 v0 crossing at x = 5 predicted by Airy front theory.

    The dispersion above gives u_t + c u_x + c beta h^2 u_xxx = 0 near the front,
    so a step becomes int_{z}^{inf} Ai(s) ds with z = (s_front - c t)/(3 c beta h^2 t)^(1/3).
    The value at the nominal front (z = 0) is 1/3; the 50 % level sits at z = -a with
    int_{-a}^{0} Ai = 1/6 (a = 0.4096). Over the travel distance d = L - x:
    delay / (d / c) = a (3 beta)^(1/3) (h / d)^(2/3).
    """
    a = brentq(lambda s: itairy(s)[2] - 1.0 / 6.0, 0.0, 2.0)
    d = _L - _X_PROBE
    return a * (3.0 * _beta(degree, cfl)) ** (1.0 / 3.0) * ((_L / nx) / d) ** (2.0 / 3.0)


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wave_p_sweep(degree: int):
    """V-EXP: step-load wave tip L2 < 2% for p = 1, 2, 3 (untrimmed)."""
    h = _wave_histories(degree)
    exact = _V0 * h["t"]
    rel = float(np.linalg.norm(h["tip"] - exact) / np.linalg.norm(exact))
    assert len(h["tip"]) == int(round(_T_END / _DT))
    assert rel < _GATE, f"V-EXP p={degree} tip L2 too large: {rel}"


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_peak_stress(degree: int):
    """V-EXP peak stress 2 sigma0 = 2.0 at the fixed end x = 0 (field E u'(0)), 2 %.

    Gated value: mean of sigma_h(0, t) over t in [1.5, 2.0] s, the latter half of
    the in-run peak interval (1, 2]. Measured -0.20 / +0.27 / -0.38 % for p = 1 / 2 / 3.
    The raw maximum is not gated:
    - It is the Airy overshoot of the dispersive step front, not a solver spike.
      The limit is 2 sigma0 x max int_z^inf Ai = 2 x 1.2744 (+27.4 %).
    - Kratos at nx = 40: +25.5 / +25.4 / +24.8 % for p = 1 / 2 / 3.
    - 1D reference, p = 1: +25.5 / +26.2 / +26.7 % at nx = 40 / 80 / 160. It grows
      towards +27.4 % under refinement, so it never converges to 2 sigma0.
    """
    h = _wave_histories(degree)
    t, s = h["t"], h["sigma_x0"]
    peak_exact = 2.0 * _SIGMA0  # reflection doubling at x = 0 for L/c < t < 3L/c
    t_on = _L / _C
    window = (t >= t_on + 0.5 * (_T_END - t_on)) & (t <= _T_END + 1e-12)
    plateau = float(np.mean(s[window]))
    rel = plateau / peak_exact - 1.0
    print(f"V-EXP p={degree} peak stress x=0: plateau[{t[window][0]:.3f},{t[window][-1]:.3f}]="
          f"{plateau:.5f} ({rel:+.4f}) raw max={s.max():.4f} ({s.max() / peak_exact - 1:+.4f}, "
          f"not gated) max|sigma| before L/c={np.max(np.abs(s[t < 0.9 * t_on])):.3e}")
    assert abs(rel) <= _GATE, f"peak stress {plateau} vs {peak_exact} (rel {rel:+.4f})"


def _arrival_threshold(t: np.ndarray, u: np.ndarray, dt: float, level: float) -> float:
    """First crossing of the half-step velocity (u_n - u_{n-1})/dt at t_{n-1/2}, interpolated."""
    v = np.diff(np.concatenate([[0.0], u])) / dt
    tv = t - 0.5 * dt
    k = int(np.argmax(v > level))
    assert k >= 1 and v[k] > level, "front never crossed the threshold"
    return float(tv[k - 1] + (level - v[k - 1]) / (v[k] - v[k - 1]) * dt)


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wavefront_arrival_gap_pin(degree: int):
    """PIN, open item: the 0.5 v0 threshold arrival misses 2 % at nx = 40.

    Not a pass of the plan gate. The gap is the numerical dispersion of the
    discretisation, not a probe or solver error. Measured 2026-11-19:
    +2.85 / +3.89 / +4.25 % for p = 1 / 2 / 3.
    - Kratos equals the independent 1D reference to round-off (same discrete system).
    - The lag is the Airy prediction a (3 beta)^(1/3) (h/d)^(2/3) within 10 %
      (measured ratio 1.04 / 0.98 / 0.97).
    - It converges as h^(2/3) (next test). The first-moment arrival is within
      0.5 % at nx = 40, so the wave speed itself is right.
    - Mesh needed for 2 %: nx = 80 (p = 1) and nx = 160 (p = 2, 3), see
      test_vexp_wavefront_arrival_refined_mesh.
    If this pin fails because the lag fell below 2 %, restore the 2 % gate at nx = 40.
    User decision recorded in docs/phase6/TRACK2_COMPLETE.md (open items).
    """
    h = _wave_histories(degree)
    t_k = _arrival_threshold(h["t"], h["u_probe"], _DT, 0.5 * _V0)
    r = _ref1d_history(degree, _NX, _DT, _T_FRONT)
    t_r = _arrival_threshold(r["t"], r["u_probe"], _DT, 0.5 * _V0)
    lag, pred = t_k / _T_ARRIVAL - 1.0, _airy_lag(degree, _NX)
    print(f"V-EXP p={degree} arrival(threshold) x=5.0: Kratos {t_k:.6f} s ({lag:+.4%}) "
          f"1D ref {t_r:.6f} s |diff|={abs(t_k - t_r):.1e} Airy {pred:+.4%} ratio {lag / pred:.3f}")
    assert abs(t_k - t_r) <= _ROUNDOFF_S, f"Kratos {t_k} vs 1D reference {t_r}"
    assert abs(lag / pred - 1.0) <= _AIRY_TOL, f"lag {lag:+.4%} vs Airy {pred:+.4%}"
    assert lag > _GATE, f"gap closed ({lag:+.4%}): restore the 2 % gate at nx = {_NX}"


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wavefront_arrival_lag_converges(degree: int):
    """Threshold-arrival lag of the 1D reference at CFL 0.2: order 2/3 (+-0.05), towards Airy.

    - beta from the reference's interior stencil equals the closed form (1e-4).
    - Observed order log2(lag_h / lag_h/2) over nx = 40 -> 80 -> 160 is 2/3 +- 0.05.
    - The ratio lag / Airy prediction approaches 1 under refinement.
    Measured 2026-11-19 (lag at nx = 40 / 80 / 160):
    p = 1 +2.85 / +1.77 / +1.11 %, p = 2 +3.89 / +2.48 / +1.57 %, p = 3 +4.25 / +2.72 / +1.73 %.
    """
    _, K, m = _ref1d_operators(degree, 64)
    i0 = K.shape[0] // 2
    js = np.arange(K.shape[0]) - i0
    th = np.array([0.05, 0.1, 0.15])
    w_h = np.sqrt([np.real(np.sum(K[i0] / m[i0] * np.exp(1j * js * q))) for q in th])
    dt64 = _CFL * (_L / 64) / _C
    w = (2.0 / dt64) * np.arcsin(0.5 * w_h * dt64)
    rel = 1.0 - w / (_C * th / (_L / 64))
    beta_fit = np.linalg.lstsq(np.vstack([th ** 2, th ** 4]).T, rel, rcond=None)[0][0]
    assert abs(beta_fit / _beta(degree, _CFL) - 1.0) <= 1e-4, f"beta {beta_fit} vs {_beta(degree, _CFL)}"
    nxs = (_NX, 2 * _NX, 4 * _NX)
    lags, ratios = [], []
    for nx in nxs:
        dt = _CFL * (_L / nx) / _C
        r = _ref1d_history(degree, nx, dt, _T_FRONT)
        lags.append(_arrival_threshold(r["t"], r["u_probe"], dt, 0.5 * _V0) / _T_ARRIVAL - 1.0)
        ratios.append(lags[-1] / _airy_lag(degree, nx))
    orders = [float(np.log2(lags[i] / lags[i + 1])) for i in range(2)]
    print(f"V-EXP p={degree} threshold lag nx={nxs}: " + " ".join(f"{x:+.4%}" for x in lags)
          + f" orders={orders[0]:.3f},{orders[1]:.3f} lag/Airy=" + ",".join(f"{x:.3f}" for x in ratios))
    for q in orders:
        assert abs(q - 2.0 / 3.0) <= 0.05, f"observed order {q} vs 2/3"
    assert abs(ratios[-1] - 1.0) < abs(ratios[0] - 1.0), f"no approach to Airy: {ratios}"


def _refined_nx(degree: int) -> int:
    """Smallest nx = 40 * 2^k whose Airy lag is <= (1 - 10 %) x 2 % (a priori, CFL 0.2)."""
    nx = _NX
    while _airy_lag(degree, nx) > (1.0 - _AIRY_TOL) * _GATE:
        nx *= 2
    return nx


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wavefront_arrival_refined_mesh(degree: int):
    """Kratos, threshold 0.5 v0 at x = 5.0, 2 % gate on the Airy-predicted mesh (CFL 0.2).

    nx = 80 / 160 / 160 for p = 1 / 2 / 3 (from _refined_nx, not from a measurement).
    """
    nx = _refined_nx(degree)
    dt = _CFL * (_L / nx) / _C
    h = _kratos_wave(degree, nx, dt, _T_FRONT)
    t_k = _arrival_threshold(h["t"], h["u_probe"], dt, 0.5 * _V0)
    rel = t_k / _T_ARRIVAL - 1.0
    print(f"V-EXP p={degree} arrival(threshold) nx={nx} dt={dt:.6f}: {t_k:.5f} s ({rel:+.4%}) "
          f"Airy {_airy_lag(degree, nx):+.4%}")
    assert nx > _NX
    assert abs(rel) <= _GATE, f"refined-mesh arrival {t_k:.5f} s vs {_T_ARRIVAL} s ({rel:+.4%})"


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wavefront_arrival_first_moment(degree: int):
    """V-EXP arrival at x = 5.0 from u(5, t) = v0 (t - 0.5) (first moment of the front).

    int s Ai(s) ds = 0, so the centroid of the dispersive front moves at c.
    Fit window: middle half of the analytic linear phase. The incident front
    arrives at (L - x)/c = 0.5 s and the fixed-end reflection at (L + x)/c = 1.5 s.
    Intercept = arrival (2 %), slope = v0 (2 %).
    """
    h = _wave_histories(degree)
    t, u = h["t"], h["u_probe"]
    t_a, t_b = (_L - _X_PROBE) / _C, (_L + _X_PROBE) / _C
    window = (t >= t_a + 0.25 * (t_b - t_a)) & (t <= t_b - 0.25 * (t_b - t_a))
    slope, icpt = np.polyfit(t[window], u[window], 1)
    t_i = -icpt / slope
    print(f"V-EXP p={degree} arrival(first moment) x=5.0: {t_i:.4f} s ({t_i / _T_ARRIVAL - 1:+.4f}) "
          f"slope={slope:.5f} ({slope / _V0 - 1:+.4f})")
    assert abs(t_i / _T_ARRIVAL - 1.0) <= _GATE, f"first-moment arrival {t_i}"
    assert abs(slope / _V0 - 1.0) <= _GATE, f"particle velocity {slope} vs {_V0}"


def test_explicit_mass_floor_unit():
    """floor_lumped_mass raises only sub-floor entries (no Kratos)."""
    m = EX.floor_lumped_mass(np.array([0.0, 0.5, 2.0]), 0.1)
    assert list(m) == [0.1, 0.5, 2.0]


def _build_p2_two_blocks(KM, SM, name: str, nz: int = 4, k_total: float = 2000.0):
    """Build two free-free p=2 IGA blocks (1x1xnz) with 3x3 Gauss surface contact pairs."""
    model = KM.Model()
    mp = model.CreateModelPart(name)
    for v in (KM.DISPLACEMENT, KM.VELOCITY, KM.ACCELERATION, KM.REACTION,
              SM.MIDDLE_VELOCITY, KM.NODAL_MASS, KM.FORCE_RESIDUAL,
              KM.RESIDUAL_VECTOR, SM.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(v)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.COMPUTE_LUMPED_MASS_MATRIX, True)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_ALPHA, 0.0)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_BETA, 0.0)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, _E)
    props.SetValue(KM.POISSON_RATIO, _NU)
    props.SetValue(KM.DENSITY, _RHO)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    gap0 = 0.1
    lo_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                       n_cells=(1, 1, nz), degree=2)
    hi_box = BoxSpline(lo=(0.0, 0.0, 1.0 + gap0), hi=(1.0, 1.0, 2.0 + gap0),
                       n_cells=(1, 1, nz), degree=2)
    n_cp = lo_box.n_ctrl
    vol_lo, _ = make_nurbs_volume(KM, mp, lo_box, first_node_id=1)
    els_lo = attach_cell_elements(KM, mp, vol_lo, lo_box, props, start_id=1)
    vol_hi, _ = make_nurbs_volume(KM, mp, hi_box, first_node_id=100)
    els_hi = attach_cell_elements(KM, mp, vol_hi, hi_box, props, start_id=100)
    lo_nodes = [mp.GetNode(1 + i) for i in range(n_cp)]
    hi_nodes = [mp.GetNode(100 + i) for i in range(n_cp)]

    xg1d, wg1d = np.polynomial.legendre.leggauss(3)
    pts_x, w_x = 0.5 * (xg1d + 1.0), 0.5 * wg1d
    smp_lo = mp.CreateSubModelPart("InterfaceLo")
    smp_hi = mp.CreateSubModelPart("InterfaceHi")
    s_lo, s_hi, k_pairs, tie_specs = [], [], [], []
    nu, nv, nw = lo_box.n_cp
    sid = 0
    for qy, wy in zip(pts_x, w_x):
        _, Ny = basis_1d(lo_box.knots[1], 2, float(qy), 0)
        for qx, wx in zip(pts_x, w_x):
            _, Nx = basis_1d(lo_box.knots[0], 2, float(qx), 0)
            k_pairs.append(k_total * float(wx * wy))
            nl = smp_lo.CreateNewNode(1000 + sid, float(qx), float(qy), 1.0)
            nh = smp_hi.CreateNewNode(2000 + sid, float(qx), float(qy), 1.0 + gap0)
            s_lo.append(nl)
            s_hi.append(nh)
            sid += 1
            weights = np.array([[float(Nx[0, a] * Ny[0, b])
                                 for b in range(nv) for a in range(nu)]])
            tie_specs.append((nl, nh, weights))
            mp.CreateNewCondition("PointLoadCondition3D1N", 5000 + sid, [nl.Id], props)
            mp.CreateNewCondition("PointLoadCondition3D1N", 6000 + sid, [nh.Id], props)

    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)

    lo_face = [lo_nodes[a + nu * (b + nv * (nw - 1))]
               for b in range(nv) for a in range(nu)]
    hi_face = [hi_nodes[a + nu * (b + nv * 0)]
               for b in range(nv) for a in range(nu)]
    cid = 1
    for nl, nh, weights in tie_specs:
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            add_linear_tie(KM, mp, cid, [nl.GetDof(var)],
                           [m.GetDof(var) for m in lo_face], weights)
            cid += 1
            add_linear_tie(KM, mp, cid, [nh.GetDof(var)],
                           [m.GetDof(var) for m in hi_face], weights)
            cid += 1
    return (model, mp, s_lo, s_hi, lo_nodes, hi_nodes,
            els_lo, els_hi, np.array(k_pairs, dtype=float))


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("dt", [0.0025, 0.005])
def test_vexp_c_dt_sweep(degree: int, dt: float):
    """V-EXP-C (p=1, 2): dt sweep — momentum < 1e-6, rebound, total energy < 2%."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    n_steps = int(round(1.0 / dt))
    v0 = 1.0
    if degree == 1:
        # NOTE (2026-11-16 probe): k=2000/pair @ dt=0.005 resolves the ~0.05 s
        # contact in ~10 steps -> +8% energy growth; k=500/pair (k_total=2000)
        # resolves it in ~20 steps and keeps penetration (~0.012) < 0.02.
        model, mp, s_lo, s_hi = _build_two_blocks(KM, SM, f"VEXPC_p1")
        lo_nodes = [nd for nd in mp.Nodes if nd.Id < 100]
        hi_nodes = [nd for nd in mp.Nodes if 100 <= nd.Id < 1000]
        els_lo, els_hi = [mp.GetElement(1)], [mp.GetElement(100)]
        k_pairs = np.full(len(s_lo), 500.0, dtype=float)
    else:
        (model, mp, s_lo, s_hi, lo_nodes, hi_nodes,
         els_lo, els_hi, k_pairs) = _build_p2_two_blocks(KM, SM, f"VEXPC_p2")

    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    for nd in hi_nodes + s_hi:
        nd.SetSolutionStepValue(KM.VELOCITY_Z, -v0)
    strat.Initialize()

    K_lo = _assemble_block_stiffness(KM, mp, els_lo, lo_nodes)
    K_hi = _assemble_block_stiffness(KM, mp, els_hi, hi_nodes)
    masses = {nd.Id: float(nd.GetValue(KM.NODAL_MASS)) for nd in mp.Nodes}
    m_vec = np.array([masses[nd.Id] for nd in lo_nodes + hi_nodes for _ in range(3)])
    up_ids = [nd.Id for nd in hi_nodes]
    lo_ids = [nd.Id for nd in lo_nodes]
    m_up = sum(masses[i] for i in up_ids)
    p0 = -m_up * v0
    e0 = 0.5 * m_up * v0 * v0

    def disp_vec(nodes):
        return np.array([nd.GetSolutionStepValue(v) for nd in nodes
                         for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)])

    def vel_vec(nodes):
        return np.array([nd.GetSolutionStepValue(v) for nd in nodes
                         for v in (KM.VELOCITY_X, KM.VELOCITY_Y, KM.VELOCITY_Z)])

    def gaps_vec():
        return np.array([b.GetSolutionStepValue(KM.DISPLACEMENT_Z) + b.Z0
                         - (a.GetSolutionStepValue(KM.DISPLACEMENT_Z) + a.Z0)
                         for a, b in zip(s_lo, s_hi)], dtype=float)

    def com_vel(ids):
        m = sum(masses[i] for i in ids)
        return sum(masses[i] * float(mp.GetNode(i).GetSolutionStepValue(KM.VELOCITY_Z))
                   for i in ids) / m

    gaps, energy_devs = [], []
    state = {"u_lo": disp_vec(lo_nodes), "u_hi": disp_vec(hi_nodes), "g": gaps_vec()}

    def hook(m, _i):
        g = gaps_vec()
        gaps.append(float(np.min(g)))
        forces = EX.paired_penalty_forces(g, 1.0) * k_pairs
        for (a, b), ff in zip(zip(s_lo, s_hi), forces):
            a.SetSolutionStepValue(SM.POINT_LOAD_Z, float(-ff))
            b.SetSolutionStepValue(SM.POINT_LOAD_Z, float(ff))

    def probe(_mp):
        # Evaluate total energy (KE + PE_strain + E_penalty) at time t_n:
        # VELOCITY(0) is v_n, while state holds u_n and g_n prior to Solve().
        vn = vel_vec(lo_nodes + hi_nodes)
        ke = 0.5 * float(np.sum(m_vec * vn * vn))
        pe = (0.5 * float(state["u_lo"] @ K_lo @ state["u_lo"])
              + 0.5 * float(state["u_hi"] @ K_hi @ state["u_hi"]))
        e_pen = 0.5 * float(np.sum(k_pairs * np.maximum(0.0, -state["g"]) ** 2))
        energy_devs.append(abs((ke + pe + e_pen) - e0) / e0)
        state["u_lo"] = disp_vec(lo_nodes)
        state["u_hi"] = disp_vec(hi_nodes)
        state["g"] = gaps_vec()
        return sum(masses[nd.Id] * float(nd.GetSolutionStepValue(KM.VELOCITY_Z))
                   for nd in mp.Nodes)

    out = EX.run_explicit_history(KM, mp, strat, dt, n_steps,
                                  probe_fn=probe, step_hook=hook)
    hist = np.array(out["history"])
    assert len(hist) == n_steps and np.all(np.isfinite(hist))
    assert min(gaps) < 0.0, "contact never engaged"
    assert min(gaps) > -0.02, f"excessive penetration: {min(gaps)}"
    dev = float(np.max(np.abs(hist - p0)) / abs(p0))
    assert dev < 1e-6, f"momentum not conserved: {dev}"
    v_up_end = com_vel(up_ids)
    v_lo_end = com_vel(lo_ids)
    assert abs(v_up_end) < 0.15 * v0, f"upper block did not stop: {v_up_end}"
    assert abs(v_lo_end + v0) / v0 < 0.15, f"lower block speed wrong: {v_lo_end}"
    max_de = float(max(energy_devs))
    end_de = float(energy_devs[-1])
    assert max_de < 0.02, f"total energy drift during run: {max_de:.4%}"
    assert end_de < 0.02, f"end total energy drift: {end_de:.4%}"
    print(f"V-EXP-C p={degree} dt={dt}: mom_dev={dev:.2e} min_gap={min(gaps):.4f} "
          f"v_up={v_up_end:.4f} v_lo={v_lo_end:.4f} max_dE={max_de:.4%} end_dE={end_de:.4%}")
