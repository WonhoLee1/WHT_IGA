import csv, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
CF = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/bimetallic-beam-bending/bimetallic-beam-tie-comparison.csv")
OURS = json.load(open(Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/v04/v04_sweep2.json")))
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
cf = {}
with open(CF, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        cf[(int(r["p"]), r["interface_workflow"])] = float(r["tip_uy"])
EB_TM = -0.002191780821917808  # transformed-section analytic (v04.json)
col = {1: "#64748b", 2: "#0f766e", 3: "#7c3aed"}
print("p  k   n            ours_merge      ours_tie        coreform_merge   merge/cf-1 (%)  tie-merge rel")
for p in (1, 2, 3):
    for k in range(4):
        key = f"p{p}_k{k}"
        if key in OURS:
            o = OURS[key]
            print(f"{p}  {k}  {str(o['n']):12s} {o['uy_merge']:.8e}  {o['uy_tie']:.8e}  {cf[(p,'boolean_merge')]:.8e}  {100*(o['uy_merge']/cf[(p,'boolean_merge')]-1):+.3f}  {o['rel_merge_tie']:.1e}")
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6))
ax = axs[0]
for p in (1, 2, 3):
    ks = [k for k in range(4) if f"p{p}_k{k}" in OURS]
    nz = [OURS[f"p{p}_k{k}"]["n"][2] for k in ks]
    ax.plot(nz, [OURS[f"p{p}_k{k}"]["uy_merge"] for k in ks], "-s", color=col[p], lw=1.6, ms=5, label=f"WHT_IGA p={p} (merge)")
    ax.axhline(cf[(p, "boolean_merge")], color=col[p], ls="--", lw=1, label=f"Coreform p={p}")
ax.axhline(EB_TM, color="#111827", ls=":", lw=1, label="transformed-section EB")
ax.set_xscale("log", base=2); ax.set_xlabel("cells along z (= 10·2^k)"); ax.set_ylabel("tip u_y [mm]")
ax.set_title("tip deflection vs refinement", fontsize=9.5, loc="left")
ax.legend(fontsize=6.5, frameon=False, loc="upper right"); ax.grid(True, which="both", color="#e2e8f0", lw=0.6)
ax = axs[1]
ps = [p for p in (1, 2, 3) if f"p{p}_k3" in OURS]
x = np.arange(len(ps))
d_m = [100 * (OURS[f"p{p}_k3"]["uy_merge"] / cf[(p, "boolean_merge")] - 1) for p in ps]
d_t = [100 * (OURS[f"p{p}_k3"]["uy_tie"] / cf[(p, "tie")] - 1) for p in ps]
ax.bar(x - 0.18, d_m, 0.34, color="#0f766e", label="merge vs Coreform merge")
ax.bar(x + 0.18, d_t, 0.34, color="#94a3b8", label="tie vs Coreform tie")
ax.axhline(0, color="#111827", lw=0.8)
for xi, v in zip(x, d_m): ax.text(xi - 0.18, v + (0.6 if v >= 0 else -1.6), f"{v:+.2f}", ha="center", fontsize=7.5)
ax.set_xticks(x, [f"p={p}" for p in ps]); ax.set_ylabel("deviation at finest level [%]")
ax.set_title("finest level (k=3) vs Coreform", fontsize=9.5, loc="left")
ax.legend(fontsize=7, frameon=False, loc="upper right")
fig.suptitle("Fig. 7  Bimetallic beam tip deflection: WHT_IGA vs Coreform tie comparison", fontsize=9.8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "fig7_v04_bimetal.png", dpi=170); print("saved")
