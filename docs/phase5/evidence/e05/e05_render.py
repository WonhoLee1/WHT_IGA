"""E05 result renders: run the real-SAT door-lock gate test unchanged, capture each part's
spline displacement after every accepted contact step (wrappers in this script only),
map the final (u_z = 5 mm) and a mid step onto the SAT tessellations, render x/y/z/iso."""
import sys, os
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import pyvista as pv
import tests.kratos_spikes.test_p5_e05_sat_door_lock as T5

OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/render/out/E05"
os.makedirs(OUT, exist_ok=True)
pv.OFF_SCREEN = True

ops_info = []      # (box, [ (nodes, tris), ... ]) in construction order: hinge, pin, link
history = []       # u_full_list after each successful contact step

_Op = T5.NonlinearBodyOperator


class _CapOp(_Op):
    def __init__(self, KM, mods, box, regions, *a, **kw):
        super().__init__(KM, mods, box, regions, *a, **kw)
        ops_info.append((box, [r["trim_soup"] for r in regions]))


_step = T5.solve_multibody_nlgeom_contact_step


def _cap_step(*a, **kw):
    res = _step(*a, **kw)
    history.append([np.array(u, copy=True) for u in res["u_full_list"]])
    return res


T5.NonlinearBodyOperator = _CapOp
T5.solve_multibody_nlgeom_contact_step = _cap_step
T5.test_e05_sat_door_lock_full_nlgeom_j2_quantitative()
print("captured steps:", len(history), "operators:", len(ops_info), flush=True)

names = ["hinge", "pin", "link"]


def build(us):
    blocks = []
    for k, ((box, soups), u) in enumerate(zip(ops_info, us)):
        coeffs = np.asarray(u, float).reshape(-1, 3)[:box.n_ctrl]
        for nodes, tris in soups:
            X = np.asarray(nodes, float)
            U = box.evaluate(coeffs, X)
            faces = np.hstack([np.full((len(tris), 1), 3), np.asarray(tris)]).ravel()
            pd = pv.PolyData(X, faces)
            pd.point_data["U"] = U
            pd.point_data["|U| [mm]"] = np.linalg.norm(U, axis=1)
            pd.point_data["part"] = np.full(len(X), k)
            blocks.append(pd)
    return blocks[0].merge(blocks[1:])


picks = {"final_uz5.0mm": history[-1], "mid_step": history[len(history) // 2]}
for tag, us in picks.items():
    m = build(us)
    m.save(os.path.join(OUT, f"E05_{tag}.vtp"))
    umax = float(np.max(m.point_data["|U| [mm]"]))
    b = m.bounds
    r = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
    c = m.center
    views = {"x": ((1, 0, 0), (0, 0, 1)), "y": ((0, 1, 0), (0, 0, 1)), "z": ((0, 0, 1), (0, 1, 0)), "iso": ((1, -1, 0.8), (0, 0, 1))}
    for vname, (d, up) in views.items():
        p = pv.Plotter(off_screen=True, window_size=(1000, 760))
        p.set_background("#f6f6f3")
        p.add_mesh(m.warp_by_vector("U", factor=1.0), scalars="|U| [mm]", cmap="turbo", clim=(0, umax),
                   scalar_bar_args={"title": "|U| [mm]", "color": "#151a22", "fmt": "%.3g"})
        p.add_text(f"E05 door lock {tag}  |U| max {umax:.3g} mm  (true scale)  view {vname}", font_size=9, color="#151a22")
        dn = np.asarray(d, float) / np.linalg.norm(d)
        k = 1.7 if vname == "iso" else 2.4
        p.camera_position = [tuple(np.asarray(c) + k * r * dn), tuple(c), up]
        p.reset_camera()
        fn = os.path.join(OUT, f"E05_{tag}_U_{vname}.png")
        p.screenshot(fn)
        p.close()
        print(fn, flush=True)
print("DONE")
