"""E03 extended result renders (Coreform schedule, tol_rel 1e-6).

At each step end (20 C, -55 C, 200 C) capture: spline displacement, and per Kratos
integration point the position, Cauchy stress and tracked PEEQ. Surface fields are
the value of the nearest integration point of the same part (header / glass / pins)
on the real STEP tessellation. Outputs (PNG + VTP):
  parts (geometry), |U|, von Mises (all), glass max principal (tension screening),
  metal PEEQ, and a half-section (cut y = 0) with integration-point von Mises.
"""
import sys, os, csv, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
import wht_iga.analyses.thermal_plasticity as TP
import tests.kratos_spikes.test_p5_e03_p2_coupling as T

OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/render/out/E03"
os.makedirs(OUT, exist_ok=True)
pv.OFF_SCREEN = True
G = r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal/glass-to-metal-seal-convergence-history.csv"
levels = [float(r["total_time"]) for r in csv.DictReader(open(G, encoding="utf-8")) if r["accepted"] == "True"]
MAT2PART = {"SS304L_header": "header", "S8061_glass": "glass", "Alloy52_pins": "pins"}

cap = {}
_orig = TP.run_thermal_steps


def _run(mp, KM, make_strategy, groups, increments=10, on_converged=None, **kw):
    def hook(level):
        if on_converged is not None:
            on_converged(level)
        if abs(level - round(level)) > 1e-9:
            return
        n = max(int(nd.Id) for nd in mp.Nodes)
        u = np.zeros((n, 3))
        for nd in mp.Nodes:
            u[int(nd.Id) - 1] = [nd.GetSolutionStepValue(KM.DISPLACEMENT_X), nd.GetSolutionStepValue(KM.DISPLACEMENT_Y),
                                 nd.GetSolutionStepValue(KM.DISPLACEMENT_Z)]
        peeq = groups.equivalent_plastic_strain()
        X, S, P, M = [], [], [], []
        for el in mp.Elements:
            g = groups.group_of_element.get(int(el.Id))
            if g is None:
                continue
            part = MAT2PART[groups.materials[groups.group_material[g]].name]
            xs = np.array([list(v) for v in el.CalculateOnIntegrationPoints(KM.INTEGRATION_COORDINATES, mp.ProcessInfo)])
            ss = np.array([list(v) for v in el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)])
            pe = peeq.get(int(el.Id), np.zeros(len(xs)))
            X.append(xs); S.append(ss); P.append(np.asarray(pe)); M += [part] * len(xs)
        cap[int(round(level))] = {"u": u, "X": np.vstack(X), "S": np.vstack(S), "P": np.concatenate(P), "part": np.array(M)}
    return _orig(mp, KM, make_strategy, groups, increments=increments, on_converged=hook, **kw)


TP.run_thermal_steps = _run
bundle = T.e03_p2_bundle.__wrapped__()
mats = TP.e03_materials(T._DATA_DIR)
res = T._build_and_run_e03_p2(bundle, [mats["SS304L_header"], mats["S8061_glass"], mats["Alloy52_pins"]],
                              increments=[15, 15, 15], bc_mode="rbe3", tol_rel=1e-6, tol_abs=1e-10, levels=levels)
print("steps", res["out"]["steps"], "solves", res["out"]["n_solves"], "captured", sorted(cap), flush=True)

box = bundle["box"]
part_of = {i: name for name, idxs in bundle["groups"].items() for i in idxs}
STEP = {1: ("1_cool_to_room_20C", "20 °C (after 450→20 °C)"), 2: ("2_cool_to_minus55C", "−55 °C"), 3: ("3_heat_to_200C", "200 °C")}
PARTC = {"header": 0, "glass": 1, "pins": 2}


def principal(s):
    sxx, syy, szz, sxy, syz, sxz = (s[:, i] for i in range(6))
    m = np.zeros((len(s), 3, 3))
    m[:, 0, 0], m[:, 1, 1], m[:, 2, 2] = sxx, syy, szz
    m[:, 0, 1] = m[:, 1, 0] = sxy; m[:, 1, 2] = m[:, 2, 1] = syz; m[:, 0, 2] = m[:, 2, 0] = sxz
    e = np.linalg.eigvalsh(m)
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + syz ** 2 + sxz ** 2))
    return vm, e[:, 2]


def surfaces(c):
    coeffs = c["u"][:box.n_ctrl]
    # reference positions of integration points (Kratos reports current coordinates; one fixed-point step)
    Xr = c["X"] - box.evaluate(coeffs, c["X"])
    vm, s1 = principal(c["S"])
    blocks = []
    for info in bundle["solid_info"]:
        part = part_of[info["index"]]
        sel = c["part"] == part
        tree = cKDTree(Xr[sel])
        X = np.asarray(info["tess"].nodes, float)
        _, j = tree.query(X)
        faces = np.hstack([np.full((len(info["tess"].tris), 1), 3), np.asarray(info["tess"].tris)]).ravel()
        pd = pv.PolyData(X, faces)
        U = box.evaluate(coeffs, X)
        pd.point_data["U"] = U
        pd.point_data["|U| [mm]"] = np.linalg.norm(U, axis=1)
        pd.point_data["von Mises [MPa]"] = vm[sel][j]
        pd.point_data["max principal [MPa]"] = s1[sel][j]
        pd.point_data["PEEQ"] = c["P"][sel][j]
        pd.point_data["part"] = np.full(len(X), PARTC[part])
        pd.field_data["part_name"] = [part]
        blocks.append((part, pd))
    qp = pv.PolyData(Xr)
    qp.point_data["von Mises [MPa]"] = vm
    qp.point_data["max principal [MPa]"] = s1
    qp.point_data["PEEQ"] = c["P"]
    qp.point_data["part"] = np.array([PARTC[p] for p in c["part"]])
    return blocks, qp


VIEWS = {"x": ((1, 0, 0), (0, 0, 1)), "y": ((0, 1, 0), (0, 0, 1)), "z": ((0, 0, 1), (0, 1, 0)), "iso": ((1, 1, 1), (0, 0, 1))}


def shot(meshes, scalar, title, fn, view, cmap="turbo", clim=None, warp=0.0, categories=False, points=None):
    p = pv.Plotter(off_screen=True, window_size=(1000, 760))
    p.set_background("#f6f6f3")
    allm = meshes[0].merge(meshes[1:]) if len(meshes) > 1 else meshes[0]
    b = allm.bounds; r = max(b[1] - b[0], b[3] - b[2], b[5] - b[4]); c = allm.center
    sb = {"title": scalar, "color": "#151a22", "fmt": "%.3g"}
    for m in meshes:
        mm = m.warp_by_vector("U", factor=warp) if warp else m
        if categories:
            p.add_mesh(mm, scalars=scalar, cmap=["#94a3b8", "#38bdf8", "#f59e0b"], clim=(0, 2), show_scalar_bar=False)
        else:
            p.add_mesh(mm, scalars=scalar, cmap=cmap, clim=clim, scalar_bar_args=sb)
            sb = None
            p.scalar_bars  # noqa
    if points is not None:
        p.add_mesh(points, scalars=scalar, cmap=cmap, clim=clim, point_size=6, render_points_as_spheres=True,
                   show_scalar_bar=False)
    if categories:
        p.add_legend([("header 304L", "#94a3b8"), ("glass S-8061", "#38bdf8"), ("pins Alloy 52", "#f59e0b")],
                     bcolor="#f6f6f3", size=(0.22, 0.12))
    p.add_text(title, font_size=9, color="#151a22")
    d, up = VIEWS[view]
    dn = np.asarray(d, float) / np.linalg.norm(d)
    p.camera_position = [tuple(np.asarray(c) + (1.7 if view == "iso" else 2.4) * r * dn), tuple(c), up]
    p.reset_camera()
    p.screenshot(fn); p.close()
    print(os.path.basename(fn), flush=True)


summary = {}
for lvl in sorted(cap):
    tag, lab = STEP[lvl]
    blocks, qp = surfaces(cap[lvl])
    meshes = [m for _, m in blocks]
    full = meshes[0].merge(meshes[1:]); full.save(os.path.join(OUT, f"E03_{tag}_fields.vtp")); qp.save(os.path.join(OUT, f"E03_{tag}_qp.vtp"))
    vm_all = qp.point_data["von Mises [MPa]"]
    summary[tag] = {k: float(v) for k, v in {
        "U_max_mm": np.max(full.point_data["|U| [mm]"]), "vm_max_MPa": vm_all.max(),
        "glass_s1_max_MPa": qp.point_data["max principal [MPa]"][qp.point_data["part"] == 1].max(),
        "peeq_max": qp.point_data["PEEQ"].max()}.items()}
    if lvl == 1:
        for v in ("iso", "x"):
            shot(meshes, "part", "E03 parts: 304L header, S-8061 glass, 11 Alloy 52 pins (STEP)", os.path.join(OUT, f"E03_parts_{v}.png"), v, categories=True)
    # von Mises, all parts, 4 views (true geometry)
    for v in VIEWS:
        shot(meshes, "von Mises [MPa]", f"E03 {lab}: von Mises (surface = nearest integration point)", os.path.join(OUT, f"E03_{tag}_vonMises_{v}.png"), v, clim=(0, float(vm_all.max())))
    # glass only: max principal (tension screening)
    gm = [m for p_, m in blocks if p_ == "glass"]
    g1 = qp.point_data["max principal [MPa]"][qp.point_data["part"] == 1]
    lim = float(np.max(np.abs(g1)))
    for v in ("iso", "x"):
        shot(gm, "max principal [MPa]", f"E03 {lab}: glass max principal stress", os.path.join(OUT, f"E03_{tag}_glass_s1_{v}.png"), v, cmap="coolwarm", clim=(-lim, lim))
    # metals: PEEQ
    mm = [m for p_, m in blocks if p_ != "glass"]
    pmax = max(float(qp.point_data["PEEQ"].max()), 1e-6)
    for v in ("iso", "x"):
        shot(mm, "PEEQ", f"E03 {lab}: equivalent plastic strain (header, pins)", os.path.join(OUT, f"E03_{tag}_PEEQ_{v}.png"), v, cmap="magma", clim=(0, pmax))
    # half section y >= 0 with integration points in the cut plane slab
    cut = [m.clip(normal=(0, -1, 0), origin=(0, 0, 0)) for m in meshes]
    cut = [m for m in cut if m.n_points]
    slab = qp.extract_points(np.abs(qp.points[:, 1]) < 0.15, adjacent_cells=False)
    for v in ("y", "iso"):
        shot(cut, "von Mises [MPa]", f"E03 {lab}: half section (y ≥ 0), integration points in y ≈ 0 slab", os.path.join(OUT, f"E03_{tag}_section_vonMises_{v}.png"), v,
             clim=(0, float(vm_all.max())), points=slab)
json.dump(summary, open(os.path.join(OUT, "E03_summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
print("DONE")
