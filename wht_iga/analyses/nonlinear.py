# -*- coding: utf-8 -*-
"""Quasi-static nonlinear analysis (Track 3, Phase 5 Python precursor).

Owns the S19/S21-verified recipe so every nonlinear run gets it right:
``MoveMeshFlag=True``, ``SetBufferSize(2)``, and
``AndCriteria(ResidualCriteria, DisplacementCriteria)`` (single-criterion
runs can report false convergence — S21 negative tests).

Scope: Neo-Hookean hyperelasticity (``ConstitutiveLawsApplication``),
small-strain J2 plasticity, load stepping with cutback. Polynomial
hyperelasticity (Mooney-Rivlin/Yeoh/Ogden) and MixedUP0/VMS locking-free
elements need C++ work and are out of scope (reported only).

Known Kratos 10.4 limitation: after a multi-step run the ``UpdatedLagrangian``
element's integration-point stress *output* is wrong although the solution is
right (``STEP <= 1``: last increment only; ``STEP > 1``: F0 applied twice).
Use TotalLagrangian for stress output, or recover stress from displacements.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np


def enforce_static_settings(mp) -> None:
    """Buffer size required by TL/UL time-stepping schemes."""
    mp.SetBufferSize(2)


def make_newton_strategy(KM: Any, mp, max_iters: int = 25,
                         tol_rel: float = 1e-8, tol_abs: float = 1e-10,
                         move_mesh: bool = True):
    """ResidualBasedNewtonRaphsonStrategy with the S21-mandated settings."""
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    rc = KM.ResidualCriteria(tol_rel, tol_abs)
    rc.SetEchoLevel(0)
    dc = KM.DisplacementCriteria(tol_rel, tol_abs)
    dc.SetEchoLevel(0)
    conv = KM.AndCriteria(rc, dc)
    conv.SetEchoLevel(0)
    strat = KM.ResidualBasedNewtonRaphsonStrategy(
        mp, scheme, conv, builder, int(max_iters), True, False, bool(move_mesh))
    strat.SetEchoLevel(0)
    return strat


def hyperelastic_neohookean_law(CL: Any):
    """Compressible Neo-Hookean law (S19-verified on QP geometry)."""
    return CL.HyperElastic3DLaw()


def j2_law(CL: Any):
    """Small-strain Simo J2 plasticity law (properties set separately)."""
    return CL.SmallStrainJ2Plasticity3DLaw()


def set_j2_properties(props, KM: Any, yield_stress: float,
                      hardening_modulus: float = 0.0) -> None:
    """Perfect plasticity by default (hardening_modulus=0)."""
    props.SetValue(KM.YIELD_STRESS, float(yield_stress))
    props.SetValue(KM.ISOTROPIC_HARDENING_MODULUS, float(hardening_modulus))


def solve_static(mp, strategy) -> bool:
    """Solve one static step with an initialized strategy (single-shot).

    Protocol (S19-verified): ``CloneTimeStep`` first, then apply BCs/loads,
    then call this function with NO further cloning. Any ``CloneTimeStep``
    after BC application wipes the prescribed solution (UL returns zeros).
    Returns Kratos' converged flag (see S21: verify independently).

    Single-shot only: ``Solve()`` commits element state even when Newton
    fails, and a *new* strategy re-initializes every element on
    ``Initialize()`` (fresh constitutive-law clones: plastic history and the
    UL reference deformation are lost). Multi-step runs must use
    :func:`run_load_stepping`, which keeps one strategy (review 2026-11-17:
    J2 unload gave 0.5 instead of 0.1, UL 2-step lateral displacement 1 %
    off).
    """
    strategy.Initialize()
    strategy.Solve()
    return bool(strategy.IsConverged())


def solve_increment(strategy) -> bool:
    """One Newton increment that commits element state only on convergence.

    Runs ``InitializeSolutionStep`` -> ``Predict`` -> ``SolveSolutionStep``
    and calls ``FinalizeSolutionStep`` only when converged, so a failed
    attempt never commits constitutive-law history (plastic strain, UL
    reference deformation). After a failure the same Kratos step can be
    retried (``InitializeSolutionStep`` is a no-op until the step is
    finalized). ``strategy.Initialize()`` must have been called once.
    """
    strategy.InitializeSolutionStep()
    strategy.Predict()
    converged = bool(strategy.SolveSolutionStep())
    if converged:
        strategy.FinalizeSolutionStep()
    return converged


def snapshot_state(mp, KM) -> dict[int, tuple[list[float], list[float]]]:
    """Nodal DISPLACEMENT and current coordinates (cutback restore)."""
    return {int(nd.Id): ([float(x) for x in nd.GetSolutionStepValue(KM.DISPLACEMENT)],
                         [float(nd.X), float(nd.Y), float(nd.Z)])
            for nd in mp.Nodes}


def restore_state(mp, KM, snap: dict[int, tuple[list[float], list[float]]]) -> None:
    """Restore DISPLACEMENT (current buffer) and coordinates (MoveMesh runs)."""
    for nd in mp.Nodes:
        u, x = snap[int(nd.Id)]
        nd.SetSolutionStepValue(KM.DISPLACEMENT, np.array(u, dtype=float))
        nd.X, nd.Y, nd.Z = x


def run_load_stepping(mp, KM, make_strategy, fractions: list[float],
                      apply_load: Callable[[float], None],
                      cutback: float = 0.5,
                      min_frac: float = 1e-3,
                      start_level: float = 0.0) -> dict:
    """Load stepping with cutback over absolute load levels ``fractions``.

    ``apply_load(level)`` sets the absolute load level. ``make_strategy()`` is
    called ONCE: one strategy (one scheme) for the whole history, so element
    and constitutive-law state carries over between steps (a new strategy per
    step re-initialized the elements and silently reset plastic history and
    the UL reference deformation).

    Each converged increment is one Kratos step (``CloneTimeStep`` with
    strictly increasing TIME, continued from the caller's clock). On failure
    the state is restored (displacements and moved coordinates; uncommitted
    constitutive state is discarded because the step was not finalized) and
    the same step is retried with the increment cut to ``cutback`` times the
    remaining increment ``target - level``; after a converged sub-increment
    the cut is relaxed again until ``target`` itself has converged. A target
    is counted only when the applied level equals it. Raises RuntimeError if
    the relative cut drops below ``min_frac``.

    Returns ``{"steps": targets reached, "level": last converged level,
    "n_solves": attempts, "history": per-attempt records}``.
    """
    if not 0.0 < float(cutback) < 1.0:
        raise ValueError(f"cutback must be in (0, 1), got {cutback}")
    strat = make_strategy()
    strat.Initialize()
    history: list[dict] = []
    level = float(start_level)
    reached = 0
    # TIME must strictly increase across CloneTimeStep calls (a repeated TIME
    # wipes the prescribed solution for UL). Continue from the caller's clock.
    t_abs = float(mp.ProcessInfo[KM.TIME])
    for target in (float(f) for f in fractions):
        frac = 1.0
        new_step = True
        while True:
            trial = target if frac >= 1.0 else level + frac * (target - level)
            if new_step:
                snap = snapshot_state(mp, KM)
                t_abs += 1.0
                mp.CloneTimeStep(t_abs)
            apply_load(trial)
            ok = solve_increment(strat)
            history.append({"target": target, "frac": frac, "level": trial,
                            "converged": ok})
            if ok:
                level = trial
                new_step = True
                if frac >= 1.0:
                    break
                frac = min(1.0, frac / float(cutback))
            else:
                restore_state(mp, KM, snap)
                new_step = False
                frac *= float(cutback)
                if frac < min_frac:
                    raise RuntimeError(
                        f"load stepping failed at target {target} from level "
                        f"{level} (relative cut {frac:.3g} < {min_frac})")
        reached += 1
    return {"steps": reached, "level": level, "n_solves": len(history),
            "history": history}
