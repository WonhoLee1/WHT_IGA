"""Spike: Kratos TotalLagrangianQ1P0MixedElement on our CELL QP geometry (all cell Gauss points in one geometry)."""
import sys, time, numpy as np
from tests.kratos_spikes.test_e04_cook_precursor import _build_cook, _P_RESULTANT, _U_REF, _NU49
from wht_iga.kratos_backend import load_kratos
import wht_iga.kratos_backend.elements as _EL
_EL.check_qp_element_name = lambda name: name  # spike only: bypass QP guard to re-probe Q1P0 (S39)
from wht_iga.reference.cook import right_face_traction_forces
mods = load_kratos(); KM = mods["KratosMultiphysics"]; SM = mods["StructuralMechanicsApplication"]; CL = mods["ConstitutiveLawsApplication"]
from KratosMultiphysics import ResidualBasedBlockBuilderAndSolver as BB

def run(elname, nx, p, nsteps):
    _, mp, box, right, tip = _build_cook(KM, SM, f"Q{nx}{p}{elname[:6]}", CL.HyperElastic3DLaw(), 0.495, degree=p, nx=nx, element_name=elname)
    el = next(iter(mp.Elements)); g = el.GetGeometry()
    try: nip = g.IntegrationPointsNumber()
    except Exception as e: nip = f"? ({e.__class__.__name__})"
    fb = right_face_traction_forces(np.array([n.Y0 for n in right]), _P_RESULTANT)
    dc = KM.DisplacementCriteria(1e-8, 1e-10)
    strat = KM.ResidualBasedNewtonRaphsonStrategy(mp, KM.ResidualBasedIncrementalUpdateStaticScheme(), dc,
                                                  BB(KM.SkylineLUFactorizationSolver()), 30, True, False, True)
    strat.SetEchoLevel(0); strat.Initialize()
    mp.CloneTimeStep(0.0); t = 1.0; ok = True
    for s in range(1, nsteps + 1):
        mp.CloneTimeStep(t); t += 1.0
        for nd, f in zip(right, fb):
            nd.SetSolutionStepValue(SM.POINT_LOAD_Y, float(f * s / nsteps))
        strat.Solve(); ok = ok and bool(strat.IsConverged())
    uy = float(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y))
    pr = el.GetValue(KM.PRESSURE) if el.Has(KM.PRESSURE) else None
    print(f"{elname} nx={nx} p={p} nodes/el={len(g)} ip/el={nip} conv={ok} uy={uy:.5f} p_el0={pr}", flush=True)

for elname, nx, p in [("TotalLagrangianQ1P0MixedElement3D8N", 16, 1), ("TotalLagrangianQ1P0MixedElement3D8N", 16, 2),
                      ("TotalLagrangianQ1P0MixedElement3D8N", 32, 1)]:
    t0 = time.time()
    try: run(elname, nx, p, 10)
    except Exception as e: print(elname, nx, p, "EXC", str(e)[:300], flush=True)
    print(f"   {time.time()-t0:.1f}s", flush=True)
