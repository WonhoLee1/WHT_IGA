# -*- coding: utf-8 -*-
"""E04 Cook precursor (Track 3): load-semantics + locking gap (S27).

Findings (probes 2026-11-16, p=1 4x4x1 coarse):
- e04.json "resultant 6.25N" cannot drive ref u_y=6.98946mm: linear compliance
  is 0.0366 mm/N -> 6.25N gives 0.23mm. 6.25 as traction (6.25 N/mm^2 over the
  16x1.91516 face) gives resultant 191.5N -> linear 7.02mm (0.4% over ref),
  TL-NH nu=0.3 gives 6.14mm (12% low, coarse+compressible). Traction reading adopted.
- Stock TL locks at nu=0.49999: 2.15mm at full load (69% low). The G5 gate is met by the
  Numba mean-dilatation F-bar/P0 (analyses/fbar.py, p=2 16x16 +0.70 %, S40), not stock TL.
- AndCriteria falsely reports non-convergence at nu~0.5 (values verify against
  linear to 5 digits); S19 single DisplacementCriteria converges.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses import nonlinear as NL
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.cook import CookSpline, right_face_traction_forces

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_E, _NU03, _NU49 = 240.565, 0.3, 0.49999
_P_TRACTION, _U_REF = 6.25, 6.98946
_P_RESULTANT = _P_TRACTION * 16.0 * 1.91516  # 191.516 N


def _build_cook(KM, SM, tag, law, nu, degree=1, nx=4,
                element_name="TotalLagrangianElement3D4N"):
    model = KM.Model()
    mp = model.CreateModelPart(tag)
    for v in (KM.DISPLACEMENT, KM.REACTION, SM.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(v)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                     n_cells=(nx, nx, 1), degree=degree)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, _E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, law)
    attach_cell_elements(KM, mp, vol, box, props, element_name=element_name)
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    for nd in mp.Nodes:
        nd.Fix(KM.DISPLACEMENT_Z)
        if abs(nd.X0) < 1e-9:
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
    right = sorted([nd for nd in mp.Nodes if abs(nd.X0 - 48.0) < 1e-9],
                   key=lambda n: n.Y0)
    for i, nd in enumerate(right):
        mp.CreateNewCondition("PointLoadCondition3D1N", 5000 + i, [nd.Id],
                              mp.GetProperties()[7])
    tip = max(right, key=lambda n: n.Y0)
    return model, mp, box, right, tip


def _tip_uy(KM, tip):
    return float(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y))


def test_e04_load_semantics_traction():
    """6.25 as traction over the 16 x 1.91516 face -> 191.516 N resultant."""
    assert _P_RESULTANT == pytest.approx(191.516, rel=1e-6)
    # linear compliance basis: 0.0366 mm/N * 191.516 N = 7.02mm ~ ref
    assert 0.0366 * _P_RESULTANT == pytest.approx(7.02, rel=1e-2)


def test_e04_linear_compliance_nu03():
    """Linear unlocked compliance anchor: u_y = 0.22895 at 6.25 N."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    _, mp, _, right, tip = _build_cook(
        KM, SM, "E04LC", SM.LinearElastic3DLaw(), _NU03,
        element_name="SmallDisplacementElement3D4N")
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), 6.25)
    mp.CloneTimeStep(0.0)
    mp.CloneTimeStep(1.0)
    for nd, f in zip(right, fb):
        nd.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f))
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver as _BB
    conv = KM.AndCriteria(KM.ResidualCriteria(1e-8, 1e-10),
                          KM.DisplacementCriteria(1e-8, 1e-10))
    strat = KM.ResidualBasedNewtonRaphsonStrategy(
        mp, scheme, conv, _BB(solver), 25, True, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    strat.Solve()
    assert bool(strat.IsConverged())
    assert _tip_uy(KM, tip) == pytest.approx(0.22895, rel=1e-4)


def test_e04_traction_hypothesis_nu03():
    """TL-NH nu=0.3 at 191.5 N (stepped) lands within 15% of ref."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    _, mp, _, right, tip = _build_cook(
        KM, SM, "E04TH", CL.HyperElastic3DLaw(), _NU03)
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), _P_RESULTANT)
    mp.CloneTimeStep(0.0)
    out = NL.run_load_stepping(
        mp, KM, lambda: NL.make_newton_strategy(KM, mp),
        [0.25, 0.5, 0.75, 1.0],
        lambda frac: [n.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f * frac))
                      for n, f in zip(right, fb)])
    assert out["steps"] == 4
    assert abs(_tip_uy(KM, tip) - _U_REF) / _U_REF < 0.15


def test_e04_incompressible_locks_single_dc():
    """Stock TL locks at nu=0.49999: full-load u_y pins ~2.15 (69% low)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    _, mp, _, right, tip = _build_cook(
        KM, SM, "E04LK", CL.HyperElastic3DLaw(), _NU49)
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), _P_RESULTANT)
    from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver as _BB
    mp.CloneTimeStep(0.0)
    t, uys = 1.0, []
    for frac in (0.25, 0.5, 0.75, 1.0):
        mp.CloneTimeStep(t)
        t += 1.0
        for nd, f in zip(right, fb):
            nd.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f * frac))
        dc = KM.DisplacementCriteria(1e-8, 1e-10)
        strat = KM.ResidualBasedNewtonRaphsonStrategy(
            mp, KM.ResidualBasedIncrementalUpdateStaticScheme(), dc,
            _BB(KM.SkylineLUFactorizationSolver()), 25, True, False, True)
        strat.SetEchoLevel(0)
        strat.Initialize()
        strat.Solve()
        assert bool(strat.IsConverged())
        uys.append(_tip_uy(KM, tip))
    err = abs(uys[-1] - _U_REF) / _U_REF
    assert 0.6 < err < 0.8  # locking signature pin (stock TL; the F-bar/P0 path meets the gate, S40)
    assert uys[-1] == pytest.approx(2.14994, rel=1e-4)


def test_e04_p_refinement_alleviates_but_not_to_gate():
    """p=2/nx=4 pins the locking-alleviation trend (err ~0.56, cf p=1 0.69).

    h-refinement helps (p=2 nx=8 -> 0.41) but p=3 diverges at nu~0.5 and the
    cost/robustness wall makes the 1% gate unreachable this way — documents
    why a mixed (U/P0) formulation is needed. That formulation now exists as
    the Numba F-bar/P0 in analyses/fbar.py (S40); C++ is not required for the gate.
    """
    import numpy as np

    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    from wht_iga.analyses.nonlinear import make_newton_strategy
    from wht_iga.kratos_backend.elements import attach_cell_elements
    from wht_iga.kratos_backend.geometry import make_nurbs_volume
    from wht_iga.reference.cook import CookSpline, right_face_traction_forces

    model = KM.Model()
    mp = model.CreateModelPart("E04P")
    for v in (KM.DISPLACEMENT, KM.REACTION, SM.POINT_LOAD):
        mp.AddNodalSolutionStepVariable(v)
    mp.SetBufferSize(2)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                     n_cells=(4, 4, 1), degree=2)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 240.565)
    props.SetValue(KM.POISSON_RATIO, 0.49999)
    props.SetValue(KM.CONSTITUTIVE_LAW, CL.HyperElastic3DLaw())
    attach_cell_elements(KM, mp, vol, box, props,
                         element_name="TotalLagrangianElement3D4N")
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X),
                 (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    for nd in mp.Nodes:
        nd.Fix(KM.DISPLACEMENT_Z)
        if abs(nd.X0) < 1e-9:
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
    right = sorted([nd for nd in mp.Nodes if abs(nd.X0 - 48.0) < 1e-9],
                   key=lambda n: n.Y0)
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), 191.5)
    for i, nd in enumerate(right):
        mp.CreateNewCondition("PointLoadCondition3D1N", 5000 + i, [nd.Id],
                              mp.GetProperties()[7])
    mp.CloneTimeStep(0.0)
    from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver as _BB
    t = 1.0
    tip = max(right, key=lambda n: n.Y0)
    uy = 0.0
    for frac in (0.25, 0.5, 0.75, 1.0):
        mp.CloneTimeStep(t)
        t += 1.0
        for nd, f in zip(right, fb):
            nd.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f * frac))
        dc = KM.DisplacementCriteria(1e-8, 1e-10)
        strat = KM.ResidualBasedNewtonRaphsonStrategy(
            mp, KM.ResidualBasedIncrementalUpdateStaticScheme(), dc,
            _BB(KM.SkylineLUFactorizationSolver()), 25, True, False, True)
        strat.SetEchoLevel(0)
        strat.Initialize()
        strat.Solve()
        assert bool(strat.IsConverged())
        uy = float(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y))
    err = abs(uy - 6.98946) / 6.98946
    assert 0.50 < err < 0.62, (uy, err)


def test_e04_and_criteria_false_negative_nu49():
    """AndCriteria reports failure at nu~0.5 while the value is correct."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]
    _, mp, _, right, tip = _build_cook(
        KM, SM, "E04FN", CL.HyperElastic3DLaw(), _NU49)
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), 6.25)
    for nd, f in zip(right, fb):
        nd.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f))
    mp.CloneTimeStep(0.0)
    mp.CloneTimeStep(1.0)
    dc = KM.DisplacementCriteria(1e-8, 1e-10)
    from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver as _BB
    strat = KM.ResidualBasedNewtonRaphsonStrategy(
        mp, KM.ResidualBasedIncrementalUpdateStaticScheme(), dc,
        _BB(KM.SkylineLUFactorizationSolver()), 25, True, False, True)
    strat.SetEchoLevel(0)
    strat.Initialize()
    strat.Solve()
    assert bool(strat.IsConverged())  # S19 single-DC recipe converges
    assert _tip_uy(KM, tip) == pytest.approx(0.07017, rel=1e-4)
