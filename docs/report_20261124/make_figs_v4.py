"""Report v4 (2026-11-24) figures: 1, 3, 4, 8 updated; 10, 11, 12 new."""
import csv, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

SP = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad")
OUT = SP / "report" / "fig"
ON = SP / "overnight"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
TEAL, PURPLE, AMBER, RED, GREY, INK = "#0f766e", "#7c3aed", "#b45309", "#b91c1c", "#64748b", "#111827"

# ---------------- Fig 1: deviation vs gate
rows = [
    ("E01 crown σxx vs Roark", 2.43, 5.0),
    ("V03 tip traction vs EB (spec input)", 0.595, 5.0),
    ("V04 EB section vs beam theory", 0.532, 1.0),
    ("V05 EB vs 7-box tie", 1.28, 2.0),
    ("V08 R1 h=0.125 F_N vs Hertz", 0.24, 4.0),
    ("V09 f1 p=1 vs 80.8 Hz", 0.874, 1.0),
    ("V10 HHT MD f1 vs Coreform", 0.0069, 2.0),
    ("V12 T2 steady vs ref", 0.0008, 1.0),
    ("V12 T3 p=2 (time-converged)", 0.13, 2.0),
    ("V12 T3 p=1 10 cells, BE dt 0.5", 1.56, 2.0),
    ("V12 T4 T_E p=1", 1.13, 3.0),
    ("E03 total volume vs B-Rep", 0.145, 1.0),
    ("E04 ν=0.495 C⁰ vs C¹ limits", 0.24, 1.0),
    ("E04 ν=0.495 C⁰ 64x64 vs C¹ limit", 2.42, 1.0),
    ("E05 RF3 vs -46.946 kN", 0.39, 5.0),
]
fig, ax = plt.subplots(figsize=(7.4, 5.4))
y = np.arange(len(rows))[::-1]
for yi, (lab, v, g) in zip(y, rows):
    ax.plot([0, g], [yi, yi], color="#cbd5e1", lw=6, solid_capstyle="butt", zorder=1)
    ax.scatter([v], [yi], s=46, color=TEAL if v <= g else RED, zorder=3)
    ax.text(v + 0.12, yi, f"{v:g} %", va="center", fontsize=8, color="#334155")
ax.set_yticks(y, [r[0] for r in rows], fontsize=8.5)
ax.set_xlabel("relative deviation [%]  (grey bar = gate)")
ax.set_xlim(0, 6.2)
ax.set_title("Fig. 1  Reported deviation against each case gate", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig1_gates.png", dpi=170); plt.close(fig)

# ---------------- Fig 3: Cook membrane, nu = 0.495
d = json.load(open(ON / "e04_nu0495.json"))
ser = {"C1": {}, "C0": {}, "p1": {}}
for k, v in d.items():
    ser["C0" if k.startswith("C0") else "C1" if k.startswith("C1") else "p1"][v["n"]] = v["uy"]
lim1, lim0 = 6.964, 6.947
fig, ax = plt.subplots(figsize=(7.4, 4.1))
ax.axhspan(lim1 * 0.99, lim1 * 1.01, color="#dcfce7", alpha=0.7, lw=0, label="±1 % of C¹ limit")
ax.axhline(lim1, color=TEAL, lw=1, ls="--", label=f"C¹ extrapolated limit {lim1}")
ax.axhline(lim0, color=PURPLE, lw=1, ls=":", label=f"C⁰ extrapolated limit {lim0}")
for key, st, col, lab in (("C1", "-o", TEAL, "ours p=2 C¹ F̄/P0"), ("C0", "-s", PURPLE, "ours p=2 C⁰ F̄/P0 (Coreform U(p,C^(p−2))/P0)"),
                          ("p1", "-^", GREY, "ours p=1 F̄/P0")):
    n = sorted(ser[key]); ax.plot([1 / k for k in n], [ser[key][k] for k in n], st, color=col, label=lab, ms=5)
ax.scatter([1 / 16], [5.08215], marker="x", color=RED, s=50, zorder=4, label="Kratos TL p=1 (no F̄): 5.082")
ax.scatter([1 / 16], [6.78691], marker="+", color=RED, s=70, zorder=4, label="Kratos TL p=2 (no F̄): 6.787")
ax.set_xscale("log"); ax.set_xticks([1 / k for k in (16, 32, 64)], ["1/16", "1/32", "1/64"])
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlabel("element size h (= 1/n, n cells per side)"); ax.set_ylabel("tip u_y [mm]")
ax.set_ylim(4.9, 7.8); ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("Fig. 3  Cook membrane, ν = 0.495, Neo-Hookean NLGEOM", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig3_cook.png", dpi=170); plt.close(fig)

# ---------------- Fig 4: V07 slope
labels = ["Coreform (source)", "ours p=2 P0 F̄ (Coreform-like path)", "ours p=3 C⁰/P2, standard grid 0.5–0.125", "master-plan target"]
vals = [2.38, 0.852, 2.58, 2.2]
cols = [GREY, RED, TEAL, INK]
fig, ax = plt.subplots(figsize=(7.4, 3.3))
ax.barh(labels[::-1], vals[::-1], color=cols[::-1], height=0.55)
for i, v in enumerate(vals[::-1]):
    ax.text(v + 0.04, i, f"{v:g}", va="center", fontsize=9)
ax.axvline(0.5, color="#94a3b8", ls=":", lw=1); ax.text(0.52, 3.35, "source doc > 0.5", fontsize=7.5, color=GREY)
ax.set_xlabel("log-log slope of σxx error (higher is faster)"); ax.set_xlim(0, 3)
ax.set_title("Fig. 4  V07 stress convergence slope", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig4_v07_slope.png", dpi=170); plt.close(fig)

# ---------------- Fig 8: E03 increments and Newton iterations
G = Path(r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal")
acc = [r for r in csv.DictReader(open(G / "glass-to-metal-seal-convergence-history.csv", encoding="utf-8")) if r["accepted"] == "True"]
cf_it = [int(r["equilibrium_iterations"]) for r in acc]
tol6 = json.load(open(ON / "e03_tol1e-6.json")); tol8 = json.load(open(ON / "e03_tol1e-8.json"))
it6 = [i for i, _ in tol6["iters"]]
it8 = [i for i, _ in tol8["iters"]]; ok8 = [o for _, o in tol8["iters"]]
fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.0), gridspec_kw={"width_ratios": [1.7, 1]})
ax = axs[0]
x = np.arange(1, 46)
ax.bar(x - 0.2, cf_it, width=0.4, color="#94a3b8", label="Coreform/Abaqus (default controls, ~5e-3)")
ax.bar(x + 0.2, it6, width=0.4, color=TEAL, label="ours, Kratos tol_rel = 1e-6")
for xs in (15.5, 30.5):
    ax.axvline(xs, color="#cbd5e1", ls=":")
ax.set_xlabel("accepted increment (Coreform schedule, 15 per step)"); ax.set_ylabel("Newton / equilibrium iterations")
ax.legend(fontsize=7.5, frameon=False, loc="upper left"); ax.set_ylim(0, 9)
ax.set_title("iterations per accepted increment: 0 cutbacks in both", fontsize=9.3, loc="left")
ax = axs[1]
lab = ["Coreform", "ours 1e-6", "ours 1e-8"]
solves = [45, len(it6), len(it8)]; cut = [0, 0, sum(not o for o in ok8)]
ax.bar(lab, solves, color=["#94a3b8", TEAL, RED])
for i, (s, c) in enumerate(zip(solves, cut)):
    ax.text(i, s + 0.8, f"{s} solves\n{c} cutbacks", ha="center", fontsize=8)
ax.set_ylim(0, 60); ax.set_ylabel("Newton solves (attempts)")
ax.set_title("1e-8 stalls at the J2 residual floor", fontsize=9.3, loc="left")
fig.suptitle("Fig. 8  Glass-to-metal seal (E03): replay of Coreform's 45 accepted increments", fontsize=9.8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "fig8_e03_seal.png", dpi=170); plt.close(fig)

# ---------------- Fig 10: V06 cut-cell quadrature (fraction position per quantity)
q = {r["q"]: r for r in json.load(open(ON / "v06_quad_sweep.json"))}
cf = {"Coreform h0p1": (0.584, 0.568, 28.5, 1.89), "Coreform regproj / wt2": (0.421, 0.409, 30.4, 1.98)}
keys, names = ("u", "ur", "srr", "stt"), ["u", "u_r", "σ_rr", "σ_θθ"]
fig, axs = plt.subplots(1, 2, figsize=(10.6, 3.7), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]
x = np.arange(4); w = 0.2
vals = {"ours q=2 (exact)": [100 * q[2][k] for k in keys], "ours q=1 (under-integrated)": [100 * q[1][k] for k in keys]}
vals.update(cf)
cols = {"ours q=2 (exact)": TEAL, "ours q=1 (under-integrated)": AMBER, "Coreform h0p1": PURPLE, "Coreform regproj / wt2": "#c4b5fd"}
for j, (lab, v) in enumerate(vals.items()):
    rel = [vv / vals["ours q=2 (exact)"][i] for i, vv in enumerate(v)]
    ax.bar(x + (j - 1.5) * w, rel, width=w, color=cols[lab], label=lab)
ax.axhline(1, color=INK, lw=0.7)
ax.set_xticks(x, names); ax.set_ylabel("error / ours q=2 error")
ax.set_ylim(0, 1.75); ax.legend(fontsize=7.3, frameon=False, loc="upper left", ncol=2); ax.set_title("errors relative to our exact-quadrature result", fontsize=9.3, loc="left")
ax = axs[1]
for lab, mk in (("Coreform h0p1", "D"), ("Coreform regproj / wt2", "^")):
    t = [(cf[lab][i] - 100 * q[2][k]) / (100 * q[1][k] - 100 * q[2][k]) for i, k in enumerate(keys)]
    ax.plot(x, t, mk, color=cols[lab], ms=9, label=f"{lab}: {min(t):.2f}–{max(t):.2f}")
ax.axhline(0, color=TEAL, lw=1.2, ls="--"); ax.axhline(1, color=AMBER, lw=1.2, ls="--")
ax.text(3.35, 0.02, "q=2 exact", color=TEAL, fontsize=8, ha="right", va="bottom"); ax.text(3.35, 0.98, "q=1", color=AMBER, fontsize=8, ha="right", va="top")
ax.set_xticks(x, names); ax.set_ylim(-0.1, 1.1); ax.set_ylabel("position between q=2 (0) and q=1 (1)")
ax.legend(fontsize=7.5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.08)); ax.set_title("same position for all four quantities", fontsize=9.3, loc="left")
fig.suptitle("Fig. 10  V06 L0 p=1: Coreform behaves like a partly under-integrated cut-cell rule", fontsize=9.8, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(OUT / "fig10_v06_quadrature.png", dpi=170); plt.close(fig)

# ---------------- Fig 11: V10 HHT vs Bossak vs Coreform
hh = json.load(open(ON / "v10_hht.json"))
apps = ["TRANSIENT_FIDELITY", "MODERATE_DISSIPATION"]
cfb = {"TRANSIENT_FIDELITY": 24768.6, "MODERATE_DISSIPATION": 24703.8}
cfi = {"TRANSIENT_FIDELITY": 24766.8, "MODERATE_DISSIPATION": 24702.1}
nat = {"TRANSIENT_FIDELITY": 24745.6, "MODERATE_DISSIPATION": 24681.0}
bos = {"TRANSIENT_FIDELITY": 24767.457, "MODERATE_DISSIPATION": 24602.979}
fig, ax = plt.subplots(figsize=(7.4, 3.6))
w = 0.16
for i, a in enumerate(apps):
    ref_ = cfb[a]
    series = [("Coreform immersed", cfi[a], "#94a3b8"), ("native Abaqus", nat[a], "#cbd5e1"),
              ("ours Bossak (earlier)", bos[a], RED), ("ours HHT (Kratos gen.-α, αm=0)", hh[a]["f1_hz"], TEAL)]
    for j, (lab, v, col) in enumerate(series):
        d = 100 * (v / ref_ - 1)
        ax.bar(i + (j - 1.5) * w, d, width=w, color=col, label=lab if i == 0 else None)
        ax.text(i + (j - 1.5) * w, d + (0.008 if d >= 0 else -0.008), f"{d:+.3f}", ha="center", va="bottom" if d >= 0 else "top", fontsize=6.8)
ax.axhline(0, color=INK, lw=0.8)
ax.set_xticks([0, 1], ["TRANSIENT_FIDELITY\n(HHT α=−0.05)", "MODERATE_DISSIPATION\n(HHT α=−0.41421)"])
ax.set_ylabel("f₁ deviation from Coreform body-fitted [%]"); ax.set_ylim(-0.5, 0.06)
ax.legend(fontsize=7.3, frameon=False, loc="lower left")
ax.set_title("Fig. 11  V10 fundamental frequency by integrator", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig11_v10_hht.png", dpi=170); plt.close(fig)

# ---------------- Fig 12: V12 T3 Coreform method bracket
dts = [0.1, 0.2, 0.25, 0.32, 0.4, 0.5]
tc = [37.3922, 37.3363, 37.3084, 37.2694, 37.2249, 37.1693]
fig, ax = plt.subplots(figsize=(7.4, 3.6))
ax.axhspan(36.60 * 0.98, 36.60 * 1.02, color="#dcfce7", alpha=0.7, lw=0, label="NAFEMS 36.60 °C ± 2 %")
ax.axvspan(0.1, 0.5, color="#f1f5f9", lw=0, label="Coreform increment range 0.1–0.5 s")
ax.plot(dts, tc, "-o", color=TEAL, label="ours, 10 linear cells, consistent capacity, BE")
ax.axhline(37.25446, color=PURPLE, ls="--", label="Coreform 37.25446")
ax.axhline(37.4482, color=GREY, ls=":", label="same mesh, time-converged 37.448")
ax.axhline(36.648, color=AMBER, ls="-.", label="ours p=2, time-converged 36.648")
ax.set_xlabel("constant backward-Euler increment [s]"); ax.set_ylabel("T_C at 32 s [°C]")
ax.set_ylim(35.7, 37.6)
ax.legend(fontsize=7.2, frameon=False, loc="lower left", ncol=2)
ax.set_title("Fig. 12  NAFEMS T3: Coreform value inside our reproduction of its method", loc="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT / "fig12_v12_t3.png", dpi=170); plt.close(fig)
print("figures written")
