"""Re-render E03 parts / section / PEEQ views from the saved VTP fields (no re-solve)."""
import os, json
import numpy as np
import pyvista as pv
pv.OFF_SCREEN = True
D = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/render/out/E03"
STEP = {"1_cool_to_room_20C": "20 °C (after 450->20 °C)", "2_cool_to_minus55C": "-55 °C", "3_heat_to_200C": "200 °C"}
COL = ["#94a3b8", "#38bdf8", "#f59e0b"]
NAMES = ["header 304L", "glass S-8061", "pins Alloy 52"]


def parts(m):
    return [m.threshold((k - 0.5, k + 0.5), scalars="part").extract_surface() for k in range(3)]


def cam(p, m, d, up, k):
    b = m.bounds; r = max(b[1] - b[0], b[3] - b[2], b[5] - b[4]); c = np.asarray(m.center)
    dn = np.asarray(d, float) / np.linalg.norm(d)
    p.camera_position = [tuple(c + k * r * dn), tuple(c), up]
    p.reset_camera()


def legend(p):
    p.add_legend([(n, c) for n, c in zip(NAMES, COL)], bcolor="#f6f6f3", face="rectangle", size=(0.32, 0.16), loc="lower left")


m1 = pv.read(os.path.join(D, "E03_1_cool_to_room_20C_fields.vtp"))
ps = parts(m1)
# parts: full iso, and cutaway (y >= 0 removed) iso showing the glass around the pins
for tag, clip in (("iso", False), ("cutaway_iso", True)):
    p = pv.Plotter(off_screen=True, window_size=(1000, 760)); p.set_background("#f6f6f3")
    for k, pm in enumerate(ps):
        mm = pm.clip(normal=(0, 1, 0), origin=(0, 0, 0)) if clip else pm
        if mm.n_points:
            p.add_mesh(mm, color=COL[k], smooth_shading=False)
    legend(p)
    p.add_text("E03 parts (STEP): 304L header, S-8061 glass, 11 Alloy 52 pins" + ("  - half removed (y > 0)" if clip else ""), font_size=9, color="#151a22")
    cam(p, m1, (0.6, 1, 0.5) if clip else (1, 1, 1), (0, 0, 1), 1.7)
    p.screenshot(os.path.join(D, f"E03_parts_{tag}.png")); p.close()
os.remove(os.path.join(D, "E03_parts_x.png")) if os.path.exists(os.path.join(D, "E03_parts_x.png")) else None

for tag, lab in STEP.items():
    m = pv.read(os.path.join(D, f"E03_{tag}_fields.vtp"))
    qp = pv.read(os.path.join(D, f"E03_{tag}_qp.vtp"))
    vmax = float(qp.point_data["von Mises [MPa]"].max())
    # section y = 0: nearest integration point on a 0.05 mm grid, kept where the nearest point is within 0.3 mm
    from scipy.spatial import cKDTree
    b = m.bounds
    xs = np.arange(b[0], b[1], 0.08); zs = np.arange(b[4], b[5], 0.08)
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    P = np.column_stack([XX.ravel(), np.zeros(XX.size), ZZ.ravel()])
    vm_all = np.asarray(qp.point_data["von Mises [MPa]"]); qpart = np.asarray(qp.point_data["part"])
    val = np.full(len(P), np.nan); lab_ = np.full(len(P), -1)
    cloud = pv.PolyData(P)
    for k, pm in enumerate(parts(m)):
        for body in pm.connectivity().split_bodies():
            body = body.extract_surface().clean()
            if body.n_cells < 20:
                continue
            body = body.compute_normals(auto_orient_normals=True, consistent_normals=True, cell_normals=False)
            sd = cloud.compute_implicit_distance(body)["implicit_distance"]
            enc = cloud.select_enclosed_points(body, tolerance=1e-6, check_surface=False)["SelectedPoints"].astype(bool)
            ins = (np.asarray(sd) < 0.0) if k == 1 else ((np.asarray(sd) < 0.0) & enc)
            ins &= lab_ < 0
            if ins.any():
                sel = qpart == k
                dd, jj = cKDTree(qp.points[sel]).query(P[ins])
                idx = np.flatnonzero(ins)[dd < 1.5]
                val[idx] = vm_all[sel][jj[dd < 1.5]]; lab_[idx] = k
    grid = pv.StructuredGrid(XX, np.zeros_like(XX), ZZ)
    grid.point_data["von Mises [MPa]"] = val.reshape(XX.shape).ravel(order="F")
    grid.point_data["part"] = lab_.reshape(XX.shape).ravel(order="F")
    sec = grid.threshold(-0.5, scalars="part", all_scalars=True)
    for vname, d in (("y", (0, -1, 0)), ("iso", (0.6, -1, 0.5))):
        p = pv.Plotter(off_screen=True, window_size=(1000, 760)); p.set_background("#f6f6f3")
        for pm in parts(m):
            cut = pm.clip(normal=(0, -1, 0), origin=(0, 0, 0))
            if cut.n_points:
                p.add_mesh(cut, color="#d6d9de", opacity=0.25)
        p.add_mesh(sec, scalars="von Mises [MPa]", cmap="turbo", clim=(0, vmax), nan_opacity=0.0,
                   scalar_bar_args={"title": "von Mises [MPa]", "color": "#151a22", "fmt": "%.3g"})
        p.add_text(f"E03 {lab}: section y = 0 (inside each solid: nearest integration point of that part)", font_size=9, color="#151a22")
        cam(p, m, d, (0, 0, 1), 1.7 if vname == "iso" else 2.2)
        p.screenshot(os.path.join(D, f"E03_{tag}_section_vonMises_{vname}.png")); p.close()
    # PEEQ, metals, readable colormap
    pmax = max(float(qp.point_data["PEEQ"].max()), 1e-6)
    for vname, d, up in (("iso", (1, 1, 1), (0, 0, 1)), ("x", (1, 0, 0), (0, 0, 1))):
        p = pv.Plotter(off_screen=True, window_size=(1000, 760)); p.set_background("#f6f6f3")
        sb = {"title": "PEEQ", "color": "#151a22", "fmt": "%.2e"}
        for k, pm in enumerate(parts(m)):
            if k == 1:
                p.add_mesh(pm, color="#e2e8f0", opacity=0.25)
            else:
                p.add_mesh(pm, scalars="PEEQ", cmap="viridis", clim=(0, pmax), scalar_bar_args=sb); sb = None
        p.add_text(f"E03 {lab}: equivalent plastic strain (header, pins; glass translucent)", font_size=9, color="#151a22")
        cam(p, m, d, up, 1.7 if vname == "iso" else 2.4)
        p.screenshot(os.path.join(D, f"E03_{tag}_PEEQ_{vname}.png")); p.close()
print("rerendered")
