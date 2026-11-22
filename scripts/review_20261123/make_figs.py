import re, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/fig")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# Fig 1: reported relative deviation vs gate (values copied from verification docs)
rows = [
 ("E01 crown σxx vs Roark", 2.43, 5.0),
 ("V03 curved beam tip vs T-G", 0.214, 5.0),
 ("V04 EB section vs beam theory", 0.532, 1.0),
 ("V05 EB vs 7-box tie", 1.28, 2.0),
 ("V08 R1 h=0.125 F_N vs Hertz", 0.24, 4.0),
 ("V09 f1 vs Coreform 81.8 Hz", 0.36, 1.0),
 ("V12 T2 steady vs ref", 0.0008, 1.0),
 ("V12 T4 transient T_E p=1", 1.13, 3.0),
 ("V12 T3 transient p=1 (linear grid)", 2.32, 2.0),
 ("E03 total volume vs B-Rep", 0.145, 1.0),
 ("E04 Cook u_y p=2 16x16", 0.70, 1.0),
 ("E05 RF3 vs -46.946 kN", 0.39, 5.0),
]
fig, ax = plt.subplots(figsize=(7.2, 4.6))
y = np.arange(len(rows))[::-1]
for yi, (lab, v, g) in zip(y, rows):
    ok = v <= g
    ax.plot([0, g], [yi, yi], color="#cbd5e1", lw=6, solid_capstyle="butt", zorder=1)
    ax.scatter([v], [yi], s=46, color="#0f766e" if ok else "#b91c1c", zorder=3)
    ax.text(v + 0.12, yi, f"{v:g} %", va="center", fontsize=8, color="#334155")
ax.set_yticks(y, [r[0] for r in rows], fontsize=8.5)
ax.set_xlabel("relative deviation [%]  (grey bar = gate)")
ax.set_xlim(0, 6.2)
ax.set_title("Fig. 1  Reported deviation against each case gate", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig1_gates.png", dpi=170); plt.close(fig)

# Fig 2: V08 R2 h=0.125 force history from scratch log
t, fn, ft = [], [], []
for line in open(Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/v08/nlgeom_r2.log"), encoding="utf-8"):
    m = re.match(r"t=([\d.]+) .*?F_N=([\d,]+\.\d) F_T=([\d,]+\.\d)", line)
    if m:
        t.append(float(m.group(1))); fn.append(float(m.group(2).replace(",", ""))); ft.append(float(m.group(3).replace(",", "")))
t, fn, ft = map(np.array, (t, fn, ft))
fig, ax = plt.subplots(figsize=(7.2, 3.8))
ax.axvspan(0.15, 0.70, color="#e2e8f0", alpha=0.7, lw=0, label="sliding window 0.15–0.70")
ax.plot(t, fn / 1e3, "-o", ms=3, color="#1e3a8a", label="F_N (×10³ lbf)")
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
ax2.plot(t, ft / 1e3, "-s", ms=3, color="#b45309", label="F_T (×10³ lbf)")
ax.set_xlabel("indentation / sliding parameter t [in]"); ax.set_ylabel("F_N [10³ lbf]", color="#1e3a8a"); ax2.set_ylabel("F_T [10³ lbf]", color="#b45309")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right", frameon=False)
ax.set_title(f"Fig. 2  V08 R2 NLGEOM h=0.125 force history ({len(t)} steps)", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig2_v08_history.png", dpi=170); plt.close(fig)

# Fig 3: Cook membrane convergence (values from S40 and S48 tables)
h = np.array([1/4, 1/8, 1/12, 1/16])
c1 = np.array([7.37756, 7.14649, 7.07429, 7.03806])
c0 = np.array([9.14585, 7.96336, np.nan, 7.62193])  # 4x4, 8x8, 16x16 ; 12x12 not run
p1 = np.array([6.21865, 6.72656, np.nan, 6.89904])
ref = 6.98946
fig, ax = plt.subplots(figsize=(7.2, 3.9))
ax.axhline(ref, color="#111827", lw=1.2, ls="--", label=f"Coreform reference {ref}")
ax.axhspan(ref * 0.99, ref * 1.01, color="#dcfce7", alpha=0.7, lw=0, label="±1 % gate")
ax.plot(h, c1, "-o", color="#0f766e", label="ours p=2 C¹ (F̄/P0)")
ax.plot(h[~np.isnan(c0)], c0[~np.isnan(c0)], "-s", color="#7c3aed", label="ours p=2 C⁰ (Q9/P0-equivalent)")
ax.plot(h[~np.isnan(p1)], p1[~np.isnan(p1)], "-^", color="#64748b", label="ours p=1 C⁰ (F̄/P0)")
ax.set_xlabel("element size h (= 1/n, n cells per side)"); ax.set_ylabel("tip u_y [mm]")
ax.set_xticks(h, ["1/4", "1/8", "1/12", "1/16"])
ax.legend(fontsize=7.5, frameon=False, loc="center right")
ax.set_title("Fig. 3  Cook membrane tip displacement (ν = 0.49999)", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig3_cook.png", dpi=170); plt.close(fig)

# Fig 4: V07 Kirsch stress slope vs Coreform and plan target
labels = ["Coreform (source)", "ours P0 body-fit h=0.5..0.125", "ours P2 pressure, p=3 C⁰ (S49)", "master-plan target"]
vals = [2.38, 0.852, 2.41, 2.2]
cols = ["#64748b", "#b91c1c", "#0f766e", "#111827"]
fig, ax = plt.subplots(figsize=(7.2, 3.3))
ax.barh(labels[::-1], vals[::-1], color=cols[::-1], height=0.55)
for i, v in enumerate(vals[::-1]):
    ax.text(v + 0.04, i, f"{v:g}", va="center", fontsize=9)
ax.axvline(0.5, color="#94a3b8", ls=":", lw=1); ax.text(0.5, 3.45, "source doc > 0.5", fontsize=7.5, color="#64748b")
ax.set_xlabel("log-log slope of σxx error (higher is faster)"); ax.set_xlim(0, 3)
ax.set_title("Fig. 4  V07 stress convergence slope", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig4_v07_slope.png", dpi=170); plt.close(fig)
print("points in V08 log:", len(t), "F_N range", fn.min(), fn.max())
