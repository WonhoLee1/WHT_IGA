import csv
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
G = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal")
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
rows = list(csv.DictReader(open(G / "glass-to-metal-seal-convergence-history.csv", encoding="utf-8")))
steps = {}
for r in rows:
    steps.setdefault(r["step_name"], []).append(r)
print("Coreform increments per step (accepted):")
for s, rs in steps.items():
    acc = [r for r in rs if r["accepted"] == "True"]
    print(f"  {s}: {len(acc)} accepted, attempts={len(rs)-len(acc)} rejected, total_time_end={acc[-1]['total_time']}, equil_iters_total={sum(int(r['equilibrium_iterations']) for r in acc)}")
perf = {r["phase"]: float(r["wall_s"]) for r in csv.DictReader(open(G / "glass-to-metal-seal-performance.csv", encoding="utf-8"))}
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
ax = axs[0]
cols = {"cool_to_room": "#0f766e", "cool_to_minus_55": "#7c3aed", "heat_to_200": "#b45309"}
for s, rs in steps.items():
    acc = [r for r in rs if r["accepted"] == "True"]
    t0 = float(acc[0]["total_time"]) - float(acc[0]["increment_size"])
    tt = [t0] + [float(r["total_time"]) for r in acc]
    ax.step(range(len(tt)), tt, where="post", color=cols.get(s, "k"), lw=1.8, label=f"Coreform {s} ({len(acc)} inc.)")
    ax.scatter(range(len(tt)), tt, s=12, color=cols.get(s, "k"))
ax.set_xlabel("accepted increment index (cumulative)"); ax.set_ylabel("total time (Coreform units)")
ax.set_title("Coreform increment history (46 rows)", fontsize=9.5, loc="left")
ax.legend(fontsize=7, frameon=False, loc="upper left"); ax.grid(True, color="#e2e8f0", lw=0.6)
ax = axs[1]
labels = ["Coreform\nIGA mesh\n(16 MPI)", "Coreform\nAbaqus solve\n(16 thr.)", "WHT_IGA\nprep\n(1 process)", "WHT_IGA\nsolve\n(1 process)"]
vals = [perf["IGA mesh"], perf["Abaqus solve"], 24.64, 218.56]
bars = ax.bar(labels, vals, color=["#94a3b8", "#94a3b8", "#0f766e", "#0f766e"])
for b, v in zip(bars, vals): ax.text(b.get_x() + b.get_width() / 2, v + 8, f"{v:.0f} s", ha="center", fontsize=8)
ax.set_ylabel("wall time [s]"); ax.tick_params(axis="x", labelsize=7)
ax.set_title("wall time (hardware differs)", fontsize=9.5, loc="left")
fig.suptitle("Fig. 8  Glass-to-metal seal: increment count and cost, Coreform vs WHT_IGA (E03)", fontsize=9.8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "fig8_e03_seal.png", dpi=170); print("saved")
