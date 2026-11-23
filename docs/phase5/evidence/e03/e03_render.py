"""E03 result renders: run the p=2 RBE3 thermal cycle (Coreform schedule, tol_rel 1e-6),
capture control-point displacements at the three step ends, map them onto the
real STEP tessellation of each solid, write VTP and render x / y / z / iso views."""
import sys, os, csv
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import pyvista as pv
import wht_iga.analyses.thermal_plasticity as TP
import tests.kratos_spikes.test_p5_e03_p2_coupling as T

OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/render/out/E03"
os.makedirs(OUT, exist_ok=True)
pv.OFF_SCREEN = True

G = r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal/glass-to-metal-seal-convergence-history.csv"
levels = [float(r["total_time"]) for r in csv.DictReader(open(G, encoding="utf-8")) if r["accepted"] == "True"]

captured = {}
_orig = TP.run_thermal_steps


def _run(mp, KM, make_strategy, groups, increments=10, on_converged=None, **kw):
    def hook(level):
        if on_converged is not None:
            on_converged(level)
        if abs(level - round(level)) < 1e-9:
            n = max(int(nd.Id) for nd in mp.Nodes)
            u = np.zeros((n, 3))
            for nd in mp.Nodes:
                u[int(nd.Id) - 1] = [nd.GetSolutionStepValue(KM.DISPLACEMENT_X), nd.GetSolutionStepValue(KM.DISPLACEMENT_Y),
                                     nd.GetSolutionStepValue(KM.DISPLACEMENT_Z)]
            captured[int(round(level))] = u
    return _orig(mp, KM, make_strategy, groups, increments=increments, on_converged=hook, **kw)


TP.run_thermal_steps = _run
bundle = T.e03_p2_bundle.__wrapped__()
mats = TP.e03_materials(T._DATA_DIR)
res = T._build_and_run_e03_p2(bundle, [mats["SS304L_header"], mats["S8061_glass"], mats["Alloy52_pins"]],
                              increments=[15, 15, 15], bc_mode="rbe3", tol_rel=1e-6, tol_abs=1e-10, levels=levels)
print("steps", res["out"]["steps"], "solves", res["out"]["n_solves"], "captured", sorted(captured), flush=True)

box = bundle["box"]
groups = bundle["groups"]
part_of = {}
for name, idxs in groups.items():
    for i in idxs:
        part_of[i] = name
part_code = {"header": 0, "glass": 1, "pins": 2}
step_name = {1: "1_cool_to_room_20C", 2: "2_cool_to_minus55C", 3: "3_heat_to_200C"}
for lvl, u in sorted(captured.items()):
    coeffs = u[:box.n_ctrl]
    blocks = []
    for info in bundle["solid_info"]:
        tess = info["tess"]
        X = np.asarray(tess.nodes, float)
        U = box.evaluate(coeffs, X)
        faces = np.hstack([np.full((len(tess.tris), 1), 3), np.asarray(tess.tris)]).ravel()
        pd = pv.PolyData(X, faces)
        pd.point_data["U"] = U
        pd.point_data["|U| [mm]"] = np.linalg.norm(U, axis=1)
        pd.point_data["part"] = np.full(len(X), part_code[part_of[info["index"]]])
        blocks.append(pd)
    m = blocks[0].merge(blocks[1:])
    m.save(os.path.join(OUT, f"E03_{step_name[lvl]}.vtp"))
    umax = float(np.max(m.point_data["|U| [mm]"]))
    b = m.bounds
    r = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
    warp = 0.05 * r / umax if umax > 0 else 0.0
    w = m.warp_by_vector("U", factor=warp)
    c = w.center
    views = {"x": ((1, 0, 0), (0, 0, 1)), "y": ((0, 1, 0), (0, 0, 1)), "z": ((0, 0, 1), (0, 1, 0)), "iso": ((1, 1, 1), (0, 0, 1))}
    for vname, (d, up) in views.items():
        p = pv.Plotter(off_screen=True, window_size=(1000, 760))
        p.set_background("#f6f6f3")
        p.add_mesh(w, scalars="|U| [mm]", cmap="turbo", clim=(0, umax),
                   scalar_bar_args={"title": "|U| [mm]", "color": "#151a22", "fmt": "%.3g"})
        p.add_text(f"E03 {step_name[lvl]}  |U| max {umax:.3g} mm  deformation x{warp:.0f}  view {vname}",
                   font_size=9, color="#151a22")
        k = 1.7 if vname == "iso" else 2.4
        p.camera_position = [(c[0] + k * r * d[0], c[1] + k * r * d[1], c[2] + k * r * d[2]), tuple(c), up]
        p.reset_camera()
        fn = os.path.join(OUT, f"E03_{step_name[lvl]}_U_{vname}.png")
        p.screenshot(fn)
        p.close()
        print(fn, flush=True)
print("DONE")
