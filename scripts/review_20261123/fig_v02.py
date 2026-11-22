import csv
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
G = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/roark-17p1-7a/roark-17p1-7a-refinement.csv")
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
rows = list(csv.DictReader(open(G, encoding="utf-8")))
h = [float(r["finest_h"]) for r in rows]
cf_loc = [100 * float(r["local_peak_rel_err_vs_roark"]) for r in rows]
cf_uni = [100 * float(r["uniform_peak_rel_err_vs_roark"]) for r in rows]
# WHT_IGA: uniform cut grids (level-set run, point A vs Roark) at h=1.5, 0.75; THB L3 at h=0.375 (V02-local.md)
ours_h = [1.5, 0.75]; ours_e = [100 * 0.027820341817500826, 100 * 0.00912338609258481]
thb_h, thb_e = 0.375, 100 * 0.008106
fig, ax = plt.subplots(figsize=(7.4, 4.2))
ax.loglog(h, cf_loc, "--o", color="#64748b", label="Coreform local refinement (THB)")
ax.loglog(h, cf_uni, "--^", color="#94a3b8", label="Coreform uniform (immersed, p=2)")
ax.loglog(ours_h, ours_e, "-s", color="#0f766e", lw=1.6, ms=6, label="WHT_IGA cut grid, point A (uniform)")
ax.loglog([thb_h], [thb_e], "-D", color="#7c3aed", ms=7, label="WHT_IGA THB L3 (V02-local)")
ax.axhline(1.0, color="#b42318", ls=":", lw=1, label="gate 1 % (peak vs Roark)")
ax.invert_xaxis(); ax.set_xlabel("finest element size h [mm]"); ax.set_ylabel("peak σ_max relative error vs Roark [%]")
ax.set_title("Fig. 9  Notched plate (Roark 17.1-7a): peak error vs refinement", fontsize=9.8, loc="left")
ax.legend(fontsize=7, frameon=False, loc="lower left"); ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
fig.tight_layout(); fig.savefig(OUT / "fig9_v02_roark.png", dpi=170); print("saved")
