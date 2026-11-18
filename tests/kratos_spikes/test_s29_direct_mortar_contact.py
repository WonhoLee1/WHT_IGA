# -*- coding: utf-8 -*-
"""S29: Kratos penalty mortar contact on IGA blocks through surface layers.

S22/S25: the ContactStructuralMechanics processes fail on models with
QuadraturePointGeometry elements (``GenerateFaces``). Two bypasses make the
penalty process initialize:
1. ``_interface_preprocess`` overridden to a no-op (the tied Triangle3D3
   surface layers already are the contact conditions);
2. ``ProcessInfo[IS_RESTARTED] = 1`` only during ``ExecuteInitialize`` (skips
   ``NormalCheckProcess``, which needs element faces; the layer normals are
   built outward). It is reset to 0 before elements are initialized.

Compression between two stacked p=1 blocks converges to the monolithic analytic
stress E*delta/2 as O(1/penalty) (plain Newton). Separation (S29 follow-up):
the weighted gap that drives the active set is computed from the *current node
coordinates*, so it only follows the deformation with ``MoveMeshFlag=True``;
with that and Kratos' mortar active-set criterion the interface opens under
tension. The reusable recipe (quad layers, displacement criterion) lives in
``wht_iga.analyses.kratos_contact`` and is tested in ``test_s29b``.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_E = 1000.0
_DELTA = 0.01
_APPS = ("StructuralMechanicsApplication", "ContactStructuralMechanicsApplication")


def _solve_two_blocks(penalty: float, top_uz: float, mode: str = "plain") -> dict:
    """Two unit p=1 blocks (z in [0,1], [1,2]) in contact at z = 1 through tied
    surface layers; bottom fixed in z, top prescribed.

    mode: ``plain`` (S29: Newton + AndCriteria, fixed mesh), ``mortar_fixed_mesh``
    or ``mortar_moving_mesh`` (Kratos contact solver pieces: INTERACTION flag,
    MortarAndConvergenceCriteria(DisplacementContactCriteria,
    PenaltyFrictionlessMortarConvergenceCriteria), contact Newton strategy, with
    MoveMeshFlag False / True). Returns bottom reaction ``R``, layer opening
    ``gap`` (mean z of upper minus lower layer), slave ``WEIGHTED_GAP`` by node
    and the converged flag."""
    mods = load_kratos(_APPS)
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA
    import KratosMultiphysics.ContactStructuralMechanicsApplication.penalty_contact_process as pcp

    model = KM.Model()
    mp = model.CreateModelPart("Structure")
    for v in (KM.DISPLACEMENT, KM.REACTION, KM.NORMAL, KM.NODAL_H, KM.NODAL_AREA,
              CSMA.WEIGHTED_GAP):
        mp.AddNodalSolutionStepVariable(v)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.SetBufferSize(2)
    props = mp.GetProperties()[1]
    for p in (props, mp.GetProperties()[2]):
        p.SetValue(KM.YOUNG_MODULUS, _E)
        p.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    lo = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    hi = BoxSpline(lo=(0.0, 0.0, 1.0), hi=(1.0, 1.0, 2.0), n_cells=(1, 1, 1), degree=1)
    vol_lo, _ = make_nurbs_volume(KM, mp, lo, first_node_id=1)
    attach_cell_elements(KM, mp, vol_lo, lo, props, start_id=1)
    vol_hi, _ = make_nurbs_volume(KM, mp, hi, first_node_id=100)
    attach_cell_elements(KM, mp, vol_hi, hi, props, start_id=100)
    pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    layers = {}
    for name, base in (("InterfaceLo", 1000), ("InterfaceHi", 2000)):
        smp = mp.CreateSubModelPart(name)
        layers[name] = [smp.CreateNewNode(base + i, x, y, 1.0) for i, (x, y) in enumerate(pts)]
    cprop = mp.GetProperties()[2]
    a = [n.Id for n in layers["InterfaceLo"]]
    b = [n.Id for n in layers["InterfaceHi"]]
    # outward normals: lower block +z, upper block -z
    mp.GetSubModelPart("InterfaceLo").CreateNewCondition("SurfaceCondition3D3N", 1000, [a[0], a[1], a[2]], cprop)
    mp.GetSubModelPart("InterfaceLo").CreateNewCondition("SurfaceCondition3D3N", 1001, [a[0], a[2], a[3]], cprop)
    mp.GetSubModelPart("InterfaceHi").CreateNewCondition("SurfaceCondition3D3N", 2000, [b[0], b[2], b[1]], cprop)
    mp.GetSubModelPart("InterfaceHi").CreateNewCondition("SurfaceCondition3D3N", 2001, [b[0], b[3], b[2]], cprop)
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    cid = 1
    for name, ids in (("InterfaceLo", range(1, 100)), ("InterfaceHi", range(100, 1000))):
        for sn in layers[name]:
            mn = next(n for n in mp.Nodes if n.Id in ids and abs(n.X0 - sn.X0) < 1e-12
                      and abs(n.Y0 - sn.Y0) < 1e-12 and abs(n.Z0 - 1.0) < 1e-12)
            for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
                add_linear_tie(KM, mp, cid, [sn.GetDof(var)], [mn.GetDof(var)], np.array([[1.0]]))
                cid += 1
    for n in mp.Nodes:
        if n.Id >= 1000:
            continue
        if abs(n.X0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_X)
        if abs(n.Y0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Y)
        if abs(n.Z0) < 1e-12 or abs(n.Z0 - 2.0) < 1e-12:
            n.Fix(KM.DISPLACEMENT_Z)

    class _SurfaceLayerPenaltyContact(pcp.PenaltyContactProcess):
        def _interface_preprocess(self, partial_model_part, key="0"):
            pass  # layers are already Triangle3D3 conditions; GenerateFaces fails on QP geometry

    cfg = {"model_part_name": "Structure",
           "contact_model_part": {"0": ["InterfaceLo", "InterfaceHi"]},
           "assume_master_slave": {"0": ["InterfaceLo"]},
           "contact_property_ids": {"0": 3},
           "contact_type": "Frictionless",
           "advance_ALM_parameters": {"manual_ALM": True, "penalty": float(penalty)}}
    proc = _SurfaceLayerPenaltyContact(model, KM.Parameters(json.dumps(cfg)))
    mp.ProcessInfo[KM.IS_RESTARTED] = 1  # skip NormalCheckProcess (element faces)
    proc.ExecuteInitialize()
    mp.ProcessInfo[KM.IS_RESTARTED] = 0  # elements must initialize their laws
    assert mp.ProcessInfo[KM.INITIAL_PENALTY] == penalty
    proc.ExecuteBeforeSolutionLoop()

    def start_step():
        mp.CloneTimeStep(1.0)
        mp.ProcessInfo[KM.STEP] = 1
        for n in mp.Nodes:
            if n.Id < 1000 and abs(n.Z0 - 2.0) < 1e-12:
                n.SetSolutionStepValue(KM.DISPLACEMENT_Z, float(top_uz))
        proc.ExecuteInitializeSolutionStep()  # contact search -> ComputingContact conditions

    if mode == "plain":
        start_step()
        conv = KM.AndCriteria(KM.ResidualCriteria(1e-10, 1e-12), KM.DisplacementCriteria(1e-10, 1e-12))
        strat = KM.ResidualBasedNewtonRaphsonStrategy(
            mp, KM.ResidualBasedIncrementalUpdateStaticScheme(), conv,
            KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()), 30, True, False, False)
    else:
        mp.Set(KM.INTERACTION, True)
        base = CSMA.DisplacementContactCriteria(1e-8, 1e-12, 1e-8, 1e-12, False)
        mortar = CSMA.PenaltyFrictionlessMortarConvergenceCriteria(False, False, False)
        conv = CSMA.MortarAndConvergenceCriteria(base, mortar, False, None)
        conv.SetActualizeRHSFlag(True)
        linear_solver = KM.SkylineLUFactorizationSolver()
        builder = CSMA.ContactResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        strat = CSMA.ResidualBasedNewtonRaphsonContactStrategy(
            mp, scheme, conv, builder, 30, True, True,
            mode == "mortar_moving_mesh", KM.Parameters("{}"), None, None)
    strat.SetEchoLevel(0)
    strat.Initialize()
    if mode == "plain":
        strat.Solve()
        converged = bool(strat.IsConverged())
    else:  # AnalysisStage order: strategy initialized before the step and its contact search
        start_step()
        strat.InitializeSolutionStep()
        strat.Predict()
        converged = bool(strat.SolveSolutionStep())
        # not IsConverged(): with reform_dofs=True FinalizeSolutionStep clears A, Dx, b and
        # IsConverged() re-runs the criteria on them (access violation)
        strat.FinalizeSolutionStep()
    proc.ExecuteFinalizeSolutionStep()
    z = lambda name: np.mean([n.Z0 + n.GetSolutionStepValue(KM.DISPLACEMENT_Z)
                              for n in mp.GetSubModelPart(name).Nodes])
    return {"R": float(sum(n.GetSolutionStepValue(KM.REACTION_Z) for n in mp.Nodes
                           if n.Id < 100 and abs(n.Z0) < 1e-12)),
            "gap": float(z("InterfaceHi") - z("InterfaceLo")),
            "wg": {n.Id: n.GetSolutionStepValue(CSMA.WEIGHTED_GAP) for n in mp.GetSubModelPart("InterfaceLo").Nodes},
            "converged": converged}


def test_s29_compression_converges_to_monolithic_as_one_over_penalty():
    """Bottom reaction -> E*delta/2 = 5 (monolithic 1x1x2 bar, nu = 0);
    error drops 10x per 10x penalty (measured 3.4e-3, 3.4e-4, 3.4e-5)."""
    ref = _E * _DELTA / 2.0
    errs = []
    for pen in (1e6, 1e7, 1e8):
        out = _solve_two_blocks(pen, -_DELTA)
        assert out["converged"]
        errs.append(abs(out["R"] - ref) / ref)
    assert errs[-1] < 1e-4, errs
    for coarse, fine in zip(errs, errs[1:]):
        assert 8.0 < coarse / fine < 12.0, errs


def test_s29_pull_apart_separates_with_mortar_criterion_and_moving_mesh():
    """Unilateral: pulling the top up transmits no tension, the upper block lifts
    rigidly (opening = delta) and the slave weighted gap is delta * integral(N_i):
    1/3 at the two nodes shared by both triangles, 1/6 at the other two."""
    out = _solve_two_blocks(1e6, +_DELTA, mode="mortar_moving_mesh")
    assert out["converged"]
    assert abs(out["R"]) < 1e-8 * _E * _DELTA, out["R"]
    assert out["gap"] == pytest.approx(_DELTA, abs=1e-12)
    expected = {1000: 1 / 3, 1001: 1 / 6, 1002: 1 / 3, 1003: 1 / 6}
    for nid, wg in out["wg"].items():
        assert wg == pytest.approx(_DELTA * expected[nid], rel=1e-9), (nid, wg)


def test_s29_root_cause_weighted_gap_ignores_displacement_on_fixed_mesh():
    """Root cause of the S29 'bilateral' result: without mesh motion the weighted gap
    (the active-set input) stays at its initial value 0 although the layers open by
    more than half of delta; tracking the opening it would be >= 0.5*delta/6 per node."""
    out = _solve_two_blocks(1e6, +_DELTA, mode="mortar_fixed_mesh")
    assert out["gap"] > 0.5 * _DELTA, out["gap"]
    assert max(abs(w) for w in out["wg"].values()) < 1e-12, out["wg"]
