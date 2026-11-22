import csv, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
REPO = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
CF = REPO / "research-cf-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/pressurized-sphere/pressurized-sphere-qoi.csv"
EV = REPO / "docs/phase2/evidence/v06_nnlsfix"
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
# Coreform: (p, level, quantity) -> (h, err, tol)
cf = {}
with open(CF, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        cf[(int(r["p"]), int(r["level"]), r["quantity"])] = (float(r["h"]), float(r["rel_l2_error_vs_lame"]), float(r["tolerance"]))
# ours (drop mode): json key -> Coreform quantity
ours = {}
for p in (1, 2, 3):
    for lv in (0, 1, 2):
        d = json.load(open(EV / f"L{lv}_p{p}_drop.json", encoding="utf-8"))
        ours[(p, lv)] = d
qmap = {"u": "displacement", "ur": "radial_displacement", "srr": "radial_stress", "stt": "tangential_stress"}
print("p lvl qty            ours         coreform     tol    ours/cf  pass_ours")
rows = []
for (p, lv), d in sorted(ours.items()):
    for k, q in qmap.items():
        key = (p, lv, q)
        if key not in cf:
            continue
        h, e_cf, tol = cf[key]
        e_u = d[k]
        ok = e_u <= tol
        rows.append((p, lv, q, h, e_u, e_cf, tol, ok))
        print(f"{p}  {lv}  {q:20s} {e_u:.4e}  {e_cf:.4e}  {tol:.3g}  {e_u/e_cf:6.2f}  {ok}")
fig, axs = plt.subplots(1, 4, figsize=(11.5, 3.6), sharey=True)
colors = {1: "#64748b", 2: "#0f766e", 3: "#7c3aed"}
for ax, (k, q) in zip(axs, qmap.items()):
    for p in (1, 2, 3):
        pts = [(cf[(p, lv, q)][0], cf[(p, lv, q)][1]) for lv in (0, 1, 2) if (p, lv, q) in cf]
        hs, es = zip(*pts)
        ax.loglog(hs, es, "--o", color=colors[p], ms=4, lw=1, alpha=0.6, label=f"Coreform p={p}")
        ou = [(ours[(p, lv)]["h"], ours[(p, lv)][k]) for lv in (0, 1, 2)]
        ax.loglog([a for a, _ in ou], [b for _, b in ou], "-s", color=colors[p], ms=5, lw=1.6, label=f"WHT_IGA p={p}")
    ax.set_title(q.replace("_", " "), fontsize=9.5)
    ax.set_xlabel("h  [mm]", fontsize=8.5)
    ax.invert_xaxis()
    ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
axs[0].set_ylabel("relative L² error vs Lamé")
axs[0].legend(fontsize=6.5, frameon=False, loc="lower left")
fig.suptitle("Fig. 5  Pressurized sphere: error–h curves, WHT_IGA (drop mode, solid) vs Coreform (dashed)", fontsize=10, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "fig5_v06_sphere.png", dpi=170)
print("saved", OUT / "fig5_v06_sphere.png")
