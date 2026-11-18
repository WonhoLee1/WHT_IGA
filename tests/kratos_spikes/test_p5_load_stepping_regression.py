# -*- coding: utf-8 -*-
"""Regression: ``analyses.nonlinear.run_load_stepping`` on path-dependent cases.

Review 2026-11-17 (dev_log/track3_review_20261117_opus.md, H1/H2) found:
- a new strategy per load step re-initialized the elements, so J2 plastic
  history and the UL reference deformation were reset every step
  (J2 unload gave 0.5 instead of 0.1; UL 2-step NH ramp was 54 % off);
- a cutback converged at ``target * frac`` was counted as reaching ``target``.

References are analytic: uniaxial J2 (E = 100, sigma_y = 0.5, perfect
plasticity) and the S19 uniaxial Neo-Hookean solution.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.kratos_spikes.test_p5_nlgeom_contact_prep import _build_uniaxial, _cauchy_xx
from tests.kratos_spikes.test_s19_constitutive_models import _analytic_uniaxial_neohookean
from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_E, _YS = 100.0, 0.5
_E_NH, _NU_NH, _LAM = 3.0, 0.3, 1.2


class _ReportFailure:
    """Strategy proxy: runs the real Newton solve but reports the first
    ``n_fail`` attempts as non-converged (state is modified, not finalized)."""

    def __init__(self, strategy, n_fail: int):
        self._s = strategy
        self._n = int(n_fail)

    def __getattr__(self, name):
        return getattr(self._s, name)

    def SolveSolutionStep(self):  # noqa: N802 (Kratos API name)
        ok = self._s.SolveSolutionStep()
        if self._n > 0:
            self._n -= 1
            return False
        return ok


def _mods():
    mods = load_kratos()
    return (mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"],
            mods["ConstitutiveLawsApplication"])


def _unit_box():
    return BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)


def _j2_bar(KM, SM, CL):
    _, mp, els = _build_uniaxial(
        KM, SM, _unit_box(), NL.j2_law(CL), "SmallDisplacementElement3D4N", 1.0,
        lambda p: (p.SetValue(KM.YOUNG_MODULUS, _E), p.SetValue(KM.POISSON_RATIO, 0.0),
                   NL.set_j2_properties(p, KM, _YS, 0.0)))
    tip = [n for n in mp.Nodes if abs(n.X0 - 1.0) < 1e-12]

    def load(level):  # level 1.0 = 2 % strain
        for n in tip:
            n.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.02 * level)
    return mp, els, load


def _nh_bar(KM, SM, CL, element_name):
    _, mp, els = _build_uniaxial(
        KM, SM, _unit_box(), NL.hyperelastic_neohookean_law(CL), element_name, 1.0,
        lambda p: (p.SetValue(KM.YOUNG_MODULUS, _E_NH), p.SetValue(KM.POISSON_RATIO, _NU_NH)))
    tip = [n for n in mp.Nodes if abs(n.X0 - 1.0) < 1e-12]

    def load(level):
        for n in tip:
            n.SetSolutionStepValue(KM.DISPLACEMENT_X, (_LAM - 1.0) * level)
    return mp, els, tip, load


def test_j2_unload_keeps_plastic_history():
    """2 % strain (plastic) then 1.6 %: elastic unloading gives E*(0.016-0.015) = 0.1."""
    KM, SM, CL = _mods()
    mp, els, load = _j2_bar(KM, SM, CL)
    out = NL.run_load_stepping(mp, KM, lambda: NL.make_newton_strategy(KM, mp),
                               [1.0, 0.8], load)
    assert out["steps"] == 2 and out["level"] == 0.8
    sig_ref = _E * (0.016 - (0.02 - _YS / _E))
    sig = _cauchy_xx(KM, mp, els)
    assert abs(sig - sig_ref) / sig_ref < 1e-9, (sig, sig_ref)


def _lateral_error(KM, mp):
    """Max |u_y - (lam2 - 1) y0| over nodes (uniaxial NH lateral contraction)."""
    lam2, _, _ = _analytic_uniaxial_neohookean(_E_NH, _NU_NH, _LAM)
    return max(abs(n.GetSolutionStepValue(KM.DISPLACEMENT_Y) - (lam2 - 1.0) * n.Y0)
               for n in mp.Nodes)


@pytest.mark.parametrize("element_name", ["UpdatedLagrangianElement3D4N",
                                          "TotalLagrangianElement3D4N"])
def test_nh_two_step_ramp_matches_analytic(element_name: str):
    """NH uniaxial ramp [0.5, 1.0]: solution (lateral stretch) analytic for UL
    and TL (old driver: UL off by 1 %); TL Cauchy stress analytic.

    UL stress *output* is not checked: Kratos 10.4 ``UpdatedLagrangian``
    reports the last increment only after ``FinalizeSolutionStep`` when
    ``STEP <= 1`` and applies F0 twice when ``STEP > 1`` (measured 0.268 /
    0.904 vs 0.587 with the exact displacement field)."""
    KM, SM, CL = _mods()
    _, sig_ref, _ = _analytic_uniaxial_neohookean(_E_NH, _NU_NH, _LAM)
    mp, els, _, load = _nh_bar(KM, SM, CL, element_name)
    out = NL.run_load_stepping(mp, KM, lambda: NL.make_newton_strategy(KM, mp),
                               [0.5, 1.0], load)
    assert out["steps"] == 2
    assert _lateral_error(KM, mp) < 1e-9
    if element_name.startswith("Total"):
        sig = _cauchy_xx(KM, mp, els)
        assert abs(sig - sig_ref) / abs(sig_ref) < 1e-9, (sig, sig_ref)


def test_cutback_reaches_target_after_failed_attempt():
    """A failed attempt (UL, moved mesh) is rolled back; cutback levels lead
    to the full target, which is the only level counted as reached."""
    KM, SM, CL = _mods()
    mp, _, tip, load = _nh_bar(KM, SM, CL, "UpdatedLagrangianElement3D4N")
    mismatch = []

    def load_and_check(level):
        # every attempt must start from a state whose moved coordinates match
        # the displacement (a failed attempt can leave a distorted mesh)
        mismatch.append(max(abs(n.X - n.X0 - n.GetSolutionStepValue(KM.DISPLACEMENT_X))
                            for n in mp.Nodes))
        load(level)

    out = NL.run_load_stepping(
        mp, KM, lambda: _ReportFailure(NL.make_newton_strategy(KM, mp), 1), [1.0],
        load_and_check)
    levels = [(h["level"], h["converged"]) for h in out["history"]]
    assert levels == [(1.0, False), (0.5, True), (1.0, True)], levels
    assert max(mismatch) < 1e-12, mismatch
    assert out["steps"] == 1 and out["level"] == 1.0
    ux = float(np.mean([n.GetSolutionStepValue(KM.DISPLACEMENT_X) for n in tip]))
    assert abs(ux - (_LAM - 1.0)) < 1e-12
    for n in mp.Nodes:  # moved coordinates consistent with the displacement
        assert abs(n.X - (n.X0 + n.GetSolutionStepValue(KM.DISPLACEMENT_X))) < 1e-12
        assert abs(n.Y - (n.Y0 + n.GetSolutionStepValue(KM.DISPLACEMENT_Y))) < 1e-12
    assert _lateral_error(KM, mp) < 1e-9


def test_failed_attempt_does_not_commit_plastic_state():
    """J2: the rolled-back attempt at 2 % must not leave plastic strain behind.

    Path 0 -> [fail at 2 %] -> 1 % -> 2 % -> 1.6 %. The stress evaluated when
    the retry starts (displacement restored to zero) must be 0; a committed
    failed attempt (eps_p = 0.015) would give a non-zero stress there.
    The final stress is the history-correct 0.1."""
    KM, SM, CL = _mods()
    mp, els, load = _j2_bar(KM, SM, CL)
    seen = []
    calls = []

    def load_and_record(level):
        if calls:  # stress of the state left by the previous attempt
            seen.append(_cauchy_xx(KM, mp, els))
        calls.append(level)
        load(level)

    out = NL.run_load_stepping(
        mp, KM, lambda: _ReportFailure(NL.make_newton_strategy(KM, mp), 1),
        [1.0, 0.8], load_and_record)
    assert [h["converged"] for h in out["history"]] == [False, True, True, True]
    assert abs(seen[0]) < 1e-12, seen  # retry starts from the virgin state
    sig_ref = _E * (0.016 - (0.02 - _YS / _E))
    sig = _cauchy_xx(KM, mp, els)
    assert abs(sig - sig_ref) / sig_ref < 1e-9, (sig, sig_ref)


def test_persistent_failure_raises():
    """No silent substitute: if every attempt fails, RuntimeError is raised."""
    KM, SM, CL = _mods()
    mp, _, _, load = _nh_bar(KM, SM, CL, "TotalLagrangianElement3D4N")
    with pytest.raises(RuntimeError, match="load stepping failed"):
        NL.run_load_stepping(
            mp, KM, lambda: _ReportFailure(NL.make_newton_strategy(KM, mp), 10**6),
            [1.0], load, min_frac=1e-2)
