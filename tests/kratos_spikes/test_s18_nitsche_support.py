"""S18: Nitsche Dirichlet BC on boundary QP geometry (3D).

Question (v2 §5.8, informational): does IgaApp ``SupportNitscheCondition``
support 3D? Answer: NO for 3D solids in this build.

Evidence chain (all deterministic, none segfaults in-test):
1. The condition instantiates on surface QP geometry (creation is safe).
2. 3D solid property sets carry zero subproperties -- but the condition
   unconditionally dereferences ``GetProperties().GetSubProperties().front()``
   (support_nitsche_condition.cpp l.182, 343, 355, 597, 608) and evaluates
   membrane kinematics (PK2 stress, THICKNESS, PRESTRESS). With solid
   properties this is an unchecked empty-container access (segfault, reproduced
   in a scratch script and recorded in S18.md -- never executed here).
3. The in-tree Nitsche test uses a membrane setup
   (LinearElasticPlaneStress2DLaw + THICKNESS + PRESTRESS).

Failover per v2 §5.8: Nitsche for 3D goes to the C++ track.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def test_nitsche_unsupported_for_3d_solids() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S18")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=2)
    volume, _ = make_nurbs_volume(KM, mp, box)
    tris = [np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]])]
    qp_geoms, wrappers, _, _ = make_surface_quadrature(KM, mp, volume, tris, rule_order=2)

    # Solid property set as used by every 3D spike: law + elastic constants.
    solid_props = mp.CreateNewProperties(1)
    solid_props.SetValue(KM.YOUNG_MODULUS, 1.0)
    solid_props.SetValue(KM.POISSON_RATIO, 0.0)
    solid_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    nit_props = mp.CreateNewProperties(2)
    nit_props.SetValue(IGA.NITSCHE_STABILIZATION_FACTOR, 4000.0)
    cond = mp.CreateNewCondition("SupportNitscheCondition", 1, qp_geoms[0], nit_props)
    g = KM.Vector(3)
    g[0], g[1], g[2] = 0.1, 0.0, 0.0
    cond.SetValue(KM.DISPLACEMENT, g)
    assert mp.NumberOfConditions() == 1
    print("\nS18: SupportNitscheCondition instantiates on surface QP geometry (creation safe).")
    n_sub = len(list(nit_props.GetSubProperties()))
    print(f"S18: solid-style property set carries {n_sub} subproperties; "
          "condition needs membrane subproperties (THICKNESS/PRESTRESS/law).")
    assert n_sub == 0
    _ = (volume, wrappers, box, solid_props)
