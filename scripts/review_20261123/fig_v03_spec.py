import csv, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
G = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/beam-bending/beam-bending-traction-tip-deflection.csv")
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
ours = json.load(open(Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/v03_spec/v03_sweep_spec.json")))
cf = {}
with open(G, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["load_case"] != "traction": continue
        cf[(int(r["p"]), int(r["level"]))] = dict(h=float(r["h_z"]), uy=float(r["tip_uy"]), ab=float(r["tip_uy_abaqus_reference"]), eb=float(r["tip_uy_euler_bernoulli"]))
AB = next(iter(cf.values()))["ab"]
EB_CF = next(iter(cf.values()))["eb"]
SCALE = 1.0
print("p lvl  h_z     coreform      ours_as_stated  ours   ratio_cf/ours  err_ours*1.05_vs_cf_abaqus")
for p in (1, 2, 3):
    for k in range(5):
        c = cf[(p, k)]; o = ours[f"p{p}_k{k}"]
        ratio = c["uy"] / o["tip_uy"]
        e_s = o["tip_uy"] * SCALE / AB - 1
        print(f"{p}  {k}  {c['h']:.5f}  {c['uy']:.6e}  {o['tip_uy']:.6e}  {o['tip_uy']*SCALE:.6e}  {ratio:6.3f}  {100*e_s:+.2f} %")
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0))
col = {1: "#64748b", 2: "#0f766e", 3: "#7c3aed"}
ax = axs[0]
for p in (1, 2, 3):
    hs = [cf[(p, k)]["h"] for k in range(5)]
    ax.loglog(hs, [abs(cf[(p, k)]["uy"]) for k in range(5)], "--o", color=col[p], ms=4, alpha=0.6, label=f"Coreform p={p}")
    ax.loglog(hs, [abs(ours[f"p{p}_k{k}"]["tip_uy"]) * SCALE for k in range(5)], "-s", color=col[p], ms=5, lw=1.6, label=f"WHT_IGA p={p} (spec load)")
ax.axhline(abs(AB), color="#111827", lw=1, ls=":", label="Abaqus finest C3D20 = 1.00596e-6")
ax.invert_xaxis(); ax.set_xlabel("h_z [m]"); ax.set_ylabel("|tip u_y| [m]")
ax.set_title("tip deflection vs element size", fontsize=9.5, loc="left")
ax.legend(fontsize=6.5, frameon=False, loc="lower left"); ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
ax = axs[1]
for p in (1, 2, 3):
    hs = [cf[(p, k)]["h"] for k in range(5)]
    e_cf = [cf[(p, k)]["uy"] / AB - 1 for k in range(5)]
    e_ours = [ours[f"p{p}_k{k}"]["tip_uy"] * SCALE / AB - 1 for k in range(5)]
    ax.semilogx(hs, [100 * e for e in e_cf], "--o", color=col[p], ms=4, alpha=0.6, label=f"Coreform p={p}")
    ax.semilogx(hs, [100 * e for e in e_ours], "-s", color=col[p], ms=5, lw=1.6, label=f"WHT_IGA p={p}")
ax.axhline(0, color="#111827", lw=0.8)
ax.invert_xaxis(); ax.set_xlabel("h_z [m]"); ax.set_ylabel("deviation from Abaqus finest [%]")
ax.set_title("signed deviation (WHT_IGA with spec E=2.0e11)", fontsize=9.5, loc="left")
ax.legend(fontsize=6.5, frameon=False, loc="lower left"); ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
fig.suptitle("Fig. 6  Cantilever end shear: Coreform (dashed) vs WHT_IGA (solid), grid spec (0.1,0.1,0.25)/2^k", fontsize=9.8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "fig6_v03_beam_spec.png", dpi=170); print("saved")
