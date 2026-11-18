# -*- coding: utf-8 -*-
"""V-EXP / V-EXP-C verification (Phase 6 gates, master-plan v2.1 section 10.1).

V-EXP: 1D elastic wave, step end load sigma0 at x = L, x = 0 fixed.
Model data come from tests/verification/data/vexp.json: p = 1..3, nx = 40,
dt = 0.005, t_end = 2.0. The analytic (d'Alembert) solution has front speed
c = 10 and particle velocity v0 = sigma0 c / E = 0.1. All gates are 2 %.
- Tip history: u(L, t) = v0 t for t <= 2L/c. Gate: relative L2.
- Peak stress: sigma = 2 sigma0 = 2.0 at the fixed end x = 0 for
  L/c < t < 3L/c. This is reflection doubling, and it is the largest stress
  anywhere in the bar.
  - Sampled as the field value sigma_h(0, t) = E sum N_i'(0) u_i(t), averaged over
    t in [1.5, 2.0] s. That is the latter half of the in-run peak interval (1, 2],
    which keeps the dispersive front out.
  - The raw maximum is not gated. It is about 2.51 (+25 %) and grows slightly
    under refinement (+25.5 / +26.2 / +26.7 % at nx = 40 / 80 / 160), towards
    the Airy overshoot (about 1.27x) of a dispersive step front.
- Wavefront arrival at x = 5.0 (field value sum N_i(5) u_i, not a control
  point): analytic 0.5 s, gate 0.01 s.
  - (a) Threshold: first crossing of the half-step velocity v_{n-1/2} = 0.5 v0.
    This is the vexp.json probe definition.
  - (b) First moment: intercept of a linear fit of u(5, t) over the middle half
    of the analytic linear phase (0.5, 1.5) s. The fitted slope must equal v0.
V-EXP-C (free-free two-block elastic collision, direct Python penalty on
  surface layers): momentum conserved (< 1e-6), rebound velocity exchange,
  energy drift (< 2%), bounded penetration, dt sweep {0.0025, 0.005}.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path

import numpy as np
import pytest

from tests.kratos_spikes import test_s24b_hybrid_explicit as S24B
from tests.kratos_spikes.test_s24b_hybrid_explicit import (
    _apply_tip_traction,
    _build_wave_model,
)
from tests.kratos_spikes.test_s25_explicit_contact import _build_two_blocks
from wht_iga.dynamics import explicit as EX
from wht_iga.kratos_backend import load_kratos
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
_X_PROBE = 0.5 * _L
_T_ARRIVAL = float(_DATA["reference_values"]["analytic_wavefront_arrival_s"])
_GATE = 0.02  # master-plan v2.1 section 10.1: arrival time and peak stress within 2 %


def test_vexp_reference_data_matches_model_builder():
    """vexp.json values equal the constants the shared S24b builder actually uses."""
    assert (S24B._E, S24B._NU, S24B._RHO, S24B._SIGMA0) == (_E, _NU, _RHO, _SIGMA0)
    assert float(_DATA["model"]["length"]) == _L
    assert abs(_T_ARRIVAL - (_L - _X_PROBE) / _C) < 1e-12
    assert abs(float(_DATA["reference_values"]["sigma0_times_c_over_E"]) - _V0) < 1e-12
    assert abs(_T_END - 2.0 * _L / _C) < 1e-12  # tip formula valid up to t_end


@functools.lru_cache(maxsize=None)
def _wave_histories(degree: int) -> dict:
    """One Kratos run per degree. Per-step x-field probes (rows averaged):
    tip control point (= u(L)), u(x = 5.0), sigma(x = 0) = E u'(0)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(_NX, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, 1.0, False, f"VEXP_p{degree}")
    mp = mh["mp"]
    _apply_tip_traction(KM, SM, mp, box, _NX, 1.0, False, degree)
    nu = box.n_cp[0]
    n_rows = box.n_ctrl // nu
    first, N, dN, X, _ = box.axis_tables(0, np.array([-1.0]))  # left end of each cell
    c_probe = int(round(_X_PROBE / (_L / _NX)))
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

    n_steps = int(round(_T_END / _DT))
    out = EX.run_explicit_history(KM, mp, strat, _DT, n_steps, probe_fn=probe)
    assert len(out["history"]) == n_steps
    hist = {k: np.array(v) for k, v in rec.items()}
    hist["t"] = (np.arange(n_steps) + 1) * _DT
    return hist


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
    """V-EXP peak stress 2 sigma0 = 2.0 at the fixed end x = 0 (field E u'(0)), 2 %."""
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
def test_vexp_wavefront_arrival(degree: int):
    """V-EXP arrival at x = 5.0, threshold 0.5 v0 (vexp.json definition): 0.5 s within 2 %.

    Known failure at the vexp.json mesh (nx = 40, measured 2026-11-19):
    +2.85 / +3.89 / +4.25 % for p = 1 / 2 / 3.
    - Cause: numerical dispersion of lumped-mass central differences. The
      Airy-shaped front lags at the 50 % level.
    - The lag shrinks as h^(2/3) at CFL 0.2. For p = 2 it is +3.89 / +2.48 / +1.57 %
      at nx = 40 / 80 / 160.
    - The first-moment arrival (next test) is within 0.2 %, so the wave speed is right.
    - Open item in docs/phase6/TRACK2_COMPLETE.md. The gate is not loosened.
    """
    h = _wave_histories(degree)
    t_arr = _arrival_threshold(h["t"], h["u_probe"], _DT, 0.5 * _V0)
    rel = t_arr / _T_ARRIVAL - 1.0
    print(f"V-EXP p={degree} arrival(threshold) x=5.0: {t_arr:.4f} s ({rel:+.4f})")
    assert abs(rel) <= _GATE, f"wavefront arrival {t_arr:.4f} s vs {_T_ARRIVAL} s ({rel:+.4f})"


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wavefront_arrival_first_moment(degree: int):
    """V-EXP arrival at x = 5.0 from u(5, t) = v0 (t - 0.5) (first moment of the front).

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


@pytest.mark.parametrize("dt", [0.0025, 0.005])
def test_vexp_c_dt_sweep(dt: float):
    """V-EXP-C: dt sweep — momentum < 1e-6, rebound exchange, energy < 2%."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    n_steps = int(round(1.0 / dt))
    model, mp, s_lo, s_hi = _build_two_blocks(KM, SM, "VEXPC")
    v0 = 1.0
    # NOTE (2026-11-16 probe): k=2000 @ dt=0.005 resolves the ~0.05 s contact
    # in ~10 steps -> +8% energy growth; k=500 resolves it in ~20 steps and
    # keeps penetration (~0.012) well below the 0.02 bound at both dts.
    k_pen = 500.0
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    for nd in mp.Nodes:
        if 100 <= nd.Id < 1000 or nd.Id >= 2000:
            nd.SetSolutionStepValue(KM.VELOCITY_Z, -v0)
    strat.Initialize()
    masses = {nd.Id: float(nd.GetValue(KM.NODAL_MASS)) for nd in mp.Nodes}
    up_ids = [nd.Id for nd in mp.Nodes if 100 <= nd.Id < 1000 or nd.Id >= 2000]
    lo_ids = [nd.Id for nd in mp.Nodes if nd.Id < 100 or 1000 <= nd.Id < 2000]
    m_up = sum(masses[i] for i in up_ids)
    m_lo = sum(masses[i] for i in lo_ids)
    p0 = -m_up * v0
    ke0 = 0.5 * m_up * v0 * v0

    def com_vel(ids):
        m = sum(masses[i] for i in ids)
        return sum(masses[i] * float(mp.GetNode(i).GetSolutionStepValue(KM.VELOCITY_Z))
                   for i in ids) / m

    gaps = []

    def hook(m, _i):
        g = np.array([b.GetSolutionStepValue(KM.DISPLACEMENT_Z) + b.Z0
                      - (a.GetSolutionStepValue(KM.DISPLACEMENT_Z) + a.Z0)
                      for a, b in zip(s_lo, s_hi)], dtype=float)
        gaps.extend(g.tolist())
        for (a, b), ff in zip(zip(s_lo, s_hi), EX.paired_penalty_forces(g, k_pen)):
            a.SetSolutionStepValue(SM.POINT_LOAD_Z, float(-ff))
            b.SetSolutionStepValue(SM.POINT_LOAD_Z, float(ff))

    out = EX.run_explicit_history(KM, mp, strat, dt, n_steps,
                                  probe_fn=lambda m: sum(
                                      masses[nd.Id] * float(nd.GetSolutionStepValue(KM.VELOCITY_Z))
                                      for nd in mp.Nodes),
                                  step_hook=hook)
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
    ke_end = 0.5 * sum(masses[nd.Id] * float(nd.GetSolutionStepValue(KM.VELOCITY_Z)) ** 2
                       for nd in mp.Nodes)
    assert abs(ke_end - ke0) / ke0 < 0.02, f"energy drift: {ke_end} vs {ke0}"
    print(f"V-EXP-C dt={dt}: mom_dev={dev:.2e} min_gap={min(gaps):.4f} "
          f"v_up={v_up_end:.4f} v_lo={v_lo_end:.4f} dE={abs(ke_end-ke0)/ke0:.4f}")
