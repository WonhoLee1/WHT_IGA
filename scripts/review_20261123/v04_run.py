import json, sys, time
from pathlib import Path
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
from wht_iga.kratos_backend import load_kratos
from wht_iga.grid.knots import n_control_points
from wht_iga.reference.box import BoxSpline
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.kratos_backend.assembly_builder import AssemblyBuilder, AssemblySpec, InstanceSpec, PartSpec
from wht_iga.kratos_backend.elements import solve_linear
from wht_iga.kratos_backend.interactions import SurfaceMesh, TieSpec
from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe
import tests.kratos_spikes.test_p3_gate_g3_benchmarks as T
mods = load_kratos(); KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/v04")
OUT.mkdir(parents=True, exist_ok=True)

def run(p):
    t0 = time.perf_counter()
    box_m = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 4, 16), degree=p)
    ky = [-5.0] * (p + 1) + [-2.5] + [0.0] * p + [2.5] + [5.0] * (p + 1)
    box_m.knots[1] = np.array(ky)
    box_m.n_cp = tuple(n_control_points(k, p) for k in box_m.knots)
    mp = KM.Model().CreateModelPart(f"V04_merge_p{p}")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT); mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD); mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box_m)
    pr_cu = mp.CreateNewProperties(1); pr_cu.SetValue(KM.YOUNG_MODULUS, 1.0e5); pr_cu.SetValue(KM.POISSON_RATIO, 0.0); pr_cu.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    pr_w = mp.CreateNewProperties(2); pr_w.SetValue(KM.YOUNG_MODULUS, 4.0e5); pr_w.SetValue(KM.POISSON_RATIO, 0.0); pr_w.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    ng = 3 if p < 3 else 4
    gx, wg = np.polynomial.legendre.leggauss(ng)
    bu, bv, bw = [box_m.breakpoints(d) for d in range(3)]
    eid = 1
    for ck in range(len(bw) - 1):
        for cj in range(len(bv) - 1):
            yc = 0.5 * (bv[cj] + bv[cj + 1]); prop = pr_cu if yc < 0.0 else pr_w
            for ci in range(len(bu) - 1):
                pts = []
                for c, wc in zip(gx, wg):
                    z = 0.5 * (bw[ck] + bw[ck + 1]) + 0.5 * (bw[ck + 1] - bw[ck]) * c
                    for b_, wb in zip(gx, wg):
                        y = yc + 0.5 * (bv[cj + 1] - bv[cj]) * b_
                        for a, wa in zip(gx, wg):
                            x = 0.5 * (bu[ci] + bu[ci + 1]) + 0.5 * (bu[ci + 1] - bu[ci]) * a
                            w = 0.125 * (bu[ci + 1] - bu[ci]) * (bv[cj + 1] - bv[cj]) * (bw[ck + 1] - bw[ck]) * wa * wb * wc
                            pts.append([float(x), float(y), float(z), float(w)])
                qps = KM.GeometriesVector(); vol.CreateQuadraturePointGeometries(qps, 2, pts)
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], prop); eid += 1
    for I in box_m.boundary_control_points(2, 0):
        nd = mp.GetNode(1 + int(I)); nd.Fix(KM.DISPLACEMENT_X); nd.Fix(KM.DISPLACEMENT_Y); nd.Fix(KM.DISPLACEMENT_Z)
    lp = mp.CreateNewProperties(10)
    T._apply_face_traction_on_instance(KM, SM, mp, box_m, list(range(1, box_m.n_ctrl + 1)), 2, 1, (0.0, -0.01, 0.0), lp, 1)
    solve_linear(KM, mp, use_pardiso=True)
    U_m = np.array([[mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_X), mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_Y), mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_Z)] for i in range(box_m.n_ctrl)])
    uy_merge = float(box_m.evaluate(U_m, np.array([[0.0, 0.0, 100.0]]))[0, 1])
    # 2-part tie
    b_cu = BoxSpline((-5.0, -5.0, 0.0), (5.0, 0.0, 100.0), (2, 2, 16), degree=p)
    b_w = BoxSpline((-5.0, 0.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 16), degree=p)
    built_t = AssemblyBuilder(KM, SM).build(AssemblySpec("v04_tie", parts={
        "cu": PartSpec("cu", "iga", box=b_cu, young_modulus=1.0e5, poisson_ratio=0.0),
        "w": PartSpec("w", "iga", box=b_w, young_modulus=4.0e5, poisson_ratio=0.0)},
        instances=[InstanceSpec("cu", "1"), InstanceSpec("w", "1")], use_pardiso=True))
    mp2 = built_t.structure; icu, iw = built_t.instances["cu_1"], built_t.instances["w_1"]
    for b_box, inst in [(b_cu, icu), (b_w, iw)]:
        for I in b_box.boundary_control_points(2, 0):
            nd = mp2.GetNode(inst.node_ids[int(I)]); nd.Fix(KM.DISPLACEMENT_X); nd.Fix(KM.DISPLACEMENT_Y); nd.Fix(KM.DISPLACEMENT_Z)
    lp2 = mp2.CreateNewProperties(99)
    cid = T._apply_face_traction_on_instance(KM, SM, mp2, b_cu, icu.node_ids, 2, 1, (0.0, -0.01, 0.0), lp2, 100)
    T._apply_face_traction_on_instance(KM, SM, mp2, b_w, iw.node_ids, 2, 1, (0.0, -0.01, 0.0), lp2, cid)
    idx_cu = b_cu.boundary_control_points(1, 1); idx_w = b_w.boundary_control_points(1, 0)
    sm_cu = SurfaceMesh("cu_if", "iga", b_cu.control_points()[idx_cu], np.array([[0, 1, 2]], dtype=np.int64), node_ids=np.array([icu.node_ids[int(I)] for I in idx_cu], dtype=np.int64))
    sm_w = SurfaceMesh("w_if", "iga", b_w.control_points()[idx_w], np.array([[0, 1, 2]], dtype=np.int64), node_ids=np.array([iw.node_ids[int(I)] for I in idx_w], dtype=np.int64))
    solve_hybrid_iga_fe(KM, SM, built_t, ties=[TieSpec("v04_tie", sm_cu, sm_w, method="mpc_matching")], use_pardiso=True)
    U_cu = np.array([[mp2.GetNode(n).GetSolutionStepValue(KM.DISPLACEMENT_X), mp2.GetNode(n).GetSolutionStepValue(KM.DISPLACEMENT_Y), mp2.GetNode(n).GetSolutionStepValue(KM.DISPLACEMENT_Z)] for n in icu.node_ids])
    uy_tie = float(b_cu.evaluate(U_cu, np.array([[0.0, 0.0, 100.0]]))[0, 1])
    return dict(p=p, uy_merge=uy_merge, uy_tie=uy_tie, rel_merge_tie=abs(uy_merge - uy_tie) / abs(uy_tie), time_s=time.perf_counter() - t0)

res = {}
for p in (1, 2, 3):
    r = run(p); res[f"p{p}"] = r
    print(f"p={p} merge={r['uy_merge']:.10e} tie={r['uy_tie']:.10e} rel={r['rel_merge_tie']:.2e} t={r['time_s']:.1f}s", flush=True)
json.dump(res, open(OUT / "v04_sweep.json", "w"), indent=1)
print("DONE", flush=True)
