# -*- coding: utf-8 -*-
"""V-EXP / V-EXP-C verification (Phase 6 gates).

V-EXP (1D elastic wave, step load): tip history vs analytic u = s*c/E*t
  (L2 < 2%) and mid-point wavefront arrival (0.5 s +/- 0.03), p = 1..3.
V-EXP-C (free-free two-block elastic collision, direct Python penalty on
  surface layers): momentum conserved (< 1e-6), rebound velocity exchange,
  energy drift (< 2%), bounded penetration, dt sweep {0.0025, 0.005}.
Reference data: tests/verification/data/vexp.json, vexp_c.json.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.kratos_spikes.test_s24b_hybrid_explicit import (
    _apply_tip_traction,
    _build_wave_model,
)
from tests.kratos_spikes.test_s25_explicit_contact import _build_two_blocks
from wht_iga.dynamics import explicit as EX
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_L, _E, _NU, _RHO = 10.0, 100.0, 0.0, 1.0
_C = float(np.sqrt(_E / _RHO))
_SIGMA0 = 1.0


def _tip_and_mid_ids(box):
    nu_cp = box.n_cp[0]
    tip = [I + 1 for I in range(box.n_ctrl) if I % nu_cp == nu_cp - 1]
    mid = [I + 1 for I in range(box.n_ctrl) if I % nu_cp == nu_cp // 2]
    return tip, mid


def _run_wave(KM, SM, mp, box, nx, dt, t_end):
    nu_cp = box.n_cp[0]
    tip_ids, mid_ids = _tip_and_mid_ids(box)
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()

    def probe(m):
        return float(np.mean([m.GetNode(i).GetSolutionStepValue(KM.DISPLACEMENT_X)
                              for i in tip_ids]))

    out = EX.run_explicit_history(KM, mp, strat, dt, int(round(t_end / dt)),
                                  probe_fn=probe)
    hist = np.array(out["history"])
    t_hist = (np.arange(len(hist)) + 1) * dt
    exact = (_SIGMA0 * _C / _E) * t_hist
    rel = float(np.linalg.norm(hist - exact) / np.linalg.norm(exact))
    return rel, len(hist)


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_vexp_wave_p_sweep(degree: int):
    """V-EXP: step-load wave tip L2 < 2% for p = 1, 2, 3 (untrimmed)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    nx, dt, t_end = 40, 0.005, 2.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(nx, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, 1.0, False, f"VEXP_p{degree}")
    mp = mh["mp"]
    _apply_tip_traction(KM, SM, mp, box, nx, 1.0, False, degree)
    rel, n = _run_wave(KM, SM, mp, box, nx, dt, t_end)
    assert n == int(round(t_end / dt))
    assert rel < 0.02, f"V-EXP p={degree} tip L2 too large: {rel}"


def test_vexp_wavefront_arrival():
    """V-EXP: mid-point wavefront arrival 0.5 s +/- 0.03 (p = 2)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    nx, degree, dt, t_end = 40, 2, 0.005, 1.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_L, 1.0, 1.0),
                    n_cells=(nx, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, 1.0, False, "VEXP_arr")
    mp = mh["mp"]
    _apply_tip_traction(KM, SM, mp, box, nx, 1.0, False, degree)
    _, mid_ids = _tip_and_mid_ids(box)
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    v0 = _SIGMA0 * _C / _E

    def mid_probe(m):
        return float(np.mean([m.GetNode(i).GetSolutionStepValue(KM.DISPLACEMENT_X)
                              for i in mid_ids]))

    out = EX.run_explicit_history(KM, mp, strat, dt, int(round(t_end / dt)),
                                  probe_fn=mid_probe)
    hist = np.array(out["history"])
    # midpoint velocity proxy via discrete differences; arrival = first passage
    # of displacement above half the eventual linear trend is ill-posed early;
    # use velocity threshold like S24 Q2 (v > 0.5 * v0).
    vel = np.diff(hist, prepend=0.0) / dt
    idx = int(np.argmax(vel > 0.5 * v0))
    t_arr = float((idx + 1) * dt)
    assert abs(t_arr - 0.5) < 0.03, f"wavefront arrival {t_arr} vs 0.5 s"


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
