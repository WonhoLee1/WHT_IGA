# -*- coding: utf-8 -*-
"""S25: explicit penalty contact on IGA blocks via surface layers.

1. ``explicit_penalty_contact_process`` cannot initialize on a model containing
   QuadraturePointGeometry elements (same GenerateFaces root cause as S22).
2. Direct-penalty alternative: surface Triangle3D3 layers tied to the splines
   via LinearMasterSlaveConstraint + per-step Python penalty POINT_LOAD pairs
   + Kratos MechanicalExplicitStrategy stepping. Free-free impact conserves
   total momentum to roundoff and stays stable with bounded penetration.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from tests.kratos_spikes.test_s22_surface_contact import (
    _setup_two_blocks_with_surface_layers,
)
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_E, _NU, _RHO = 100.0, 0.0, 1.0
_C = float(np.sqrt(_E / _RHO))
_GAP0, _V0, _KPEN = 0.1, 1.0, 2000.0
_DT, _NSTEPS = 0.005, 200


def _build_two_blocks(KM, SM, name: str):
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
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, _E)
    props.SetValue(KM.POISSON_RATIO, _NU)
    props.SetValue(KM.DENSITY, _RHO)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    lo_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                       n_cells=(1, 1, 1), degree=1)
    hi_box = BoxSpline(lo=(0.0, 0.0, 1.0 + _GAP0), hi=(1.0, 1.0, 2.0 + _GAP0),
                       n_cells=(1, 1, 1), degree=1)
    vol_lo, _ = make_nurbs_volume(KM, mp, lo_box, first_node_id=1)
    attach_cell_elements(KM, mp, vol_lo, lo_box, props, start_id=1)
    vol_hi, _ = make_nurbs_volume(KM, mp, hi_box, first_node_id=100)
    attach_cell_elements(KM, mp, vol_hi, hi_box, props, start_id=100)
    smp_lo = mp.CreateSubModelPart("InterfaceLo")
    smp_hi = mp.CreateSubModelPart("InterfaceHi")
    z_lo, z_hi = 1.0, 1.0 + _GAP0
    pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    s_lo = [smp_lo.CreateNewNode(1000 + i, x, y, z_lo) for i, (x, y) in enumerate(pts)]
    s_hi = [smp_hi.CreateNewNode(2000 + i, x, y, z_hi) for i, (x, y) in enumerate(pts)]
    cprop = mp.GetProperties()[2]
    smp_lo.CreateNewCondition("SurfaceCondition3D3N", 1000, [1000, 1001, 1002], cprop)
    smp_lo.CreateNewCondition("SurfaceCondition3D3N", 1001, [1000, 1002, 1003], cprop)
    smp_hi.CreateNewCondition("SurfaceCondition3D3N", 2000, [2000, 2002, 2001], cprop)
    smp_hi.CreateNewCondition("SurfaceCondition3D3N", 2001, [2000, 2003, 2002], cprop)
    for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                 (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
        KM.VariableUtils().AddDof(v, r, mp)
    lo_face = [n for n in mp.Nodes if n.Id < 100 and abs(n.Z0 - 1.0) < 1e-12]
    hi_face = [n for n in mp.Nodes if 100 <= n.Id < 1000 and abs(n.Z0 - z_hi) < 1e-12]
    cid = 1
    for sn in list(s_lo) + list(s_hi):
        face = lo_face if sn.Id < 2000 else hi_face
        mn = next(n for n in face if abs(n.X0 - sn.X0) < 1e-12 and abs(n.Y0 - sn.Y0) < 1e-12)
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            add_linear_tie(KM, mp, cid, [sn.GetDof(var)], [mn.GetDof(var)], np.array([[1.0]]))
            cid += 1
    lcid = 5000
    for sn in list(s_lo) + list(s_hi):
        mp.CreateNewCondition("PointLoadCondition3D1N", lcid, [sn.Id], props)
        lcid += 1
    return model, mp, s_lo, s_hi


def test_s25_explicit_penalty_process_fails_on_qp_geometry():
    """explicit_penalty_contact_process init fails on QP models (GenerateFaces)."""
    mods, model, mp, s_lo, s_hi, props, _, _ = _setup_two_blocks_with_surface_layers()
    KM = mods["KratosMultiphysics"]
    import KratosMultiphysics.ContactStructuralMechanicsApplication.explicit_penalty_contact_process as epc

    cfg = {
        "model_part_name": "Structure",
        "interval": [0.0, "End"],
        "contact_model_part": {"0": ["InterfaceLo", "InterfaceHi"]},
        "assume_master_slave": {"0": ["InterfaceLo"]},
        "contact_property_ids": {"0": 3},
        "contact_type": "Frictionless",
    }
    proc = epc.ExplicitPenaltyContactProcess(
        model, KM.Parameters(json.dumps({"Parameters": cfg}))["Parameters"])
    with pytest.raises(RuntimeError, match="GenerateFaces"):
        proc.ExecuteInitialize()


def test_s25_direct_penalty_free_impact():
    """Free-free impact via surface layers + Python penalty: momentum + stability."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model, mp, s_lo, s_hi = _build_two_blocks(KM, SM, "S25")
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    for nd in mp.Nodes:
        if 100 <= nd.Id < 1000 or nd.Id >= 2000:
            nd.SetSolutionStepValue(KM.VELOCITY_Z, -_V0)
    strat.Initialize()
    masses = {nd.Id: float(nd.GetValue(KM.NODAL_MASS)) for nd in mp.Nodes}
    p0 = sum(m * (-_V0 if (100 <= i < 1000 or i >= 2000) else 0.0) for i, m in masses.items())
    assert p0 != 0.0
    t_abs, ok = 0.0, True
    hist_p, hist_gap = [], []
    for _s in range(_NSTEPS):
        gaps = []
        for a, b in zip(s_lo, s_hi):
            za = a.GetSolutionStepValue(KM.DISPLACEMENT_Z) + a.Z0
            zb = b.GetSolutionStepValue(KM.DISPLACEMENT_Z) + b.Z0
            gaps.append(zb - za)
        for (a, b), g in zip(zip(s_lo, s_hi), gaps):
            f = _KPEN * max(0.0, -g)
            a.SetSolutionStepValue(SM.POINT_LOAD_Z, float(-f))
            b.SetSolutionStepValue(SM.POINT_LOAD_Z, float(f))
        t_abs += _DT
        mp.CloneTimeStep(t_abs)
        mp.ProcessInfo[KM.STEP] = _s + 1
        strat.Solve()
        pz = sum(masses[nd.Id] * float(nd.GetSolutionStepValue(KM.VELOCITY_Z))
                 for nd in mp.Nodes)
        hist_p.append(pz)
        hist_gap.append(min(gaps))
        if not np.isfinite(pz):
            ok = False
            break
    assert ok, "direct-penalty impact went unstable"
    assert min(hist_gap) < 0.0, "contact never engaged"
    assert min(hist_gap) > -0.02, f"excessive penetration: {min(hist_gap)}"
    dev = float(np.max(np.abs(np.array(hist_p) - p0)) / abs(p0))
    assert dev < 1e-6, f"momentum not conserved: {dev}"
    print(f"S25 impact: steps={len(hist_p)} mom_dev={dev:.2e} min_gap={min(hist_gap):.4f}")
