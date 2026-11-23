import sys, json, math
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import tests.kratos_spikes.test_p4_gate_g4_benchmarks as G
from wht_iga.dynamics import (BossakConfig, GeneralizedAlphaConfig, solve_bossak_linear,
                              solve_generalized_alpha_linear, solve_eigenmodes)
from wht_iga.reference.box import BoxSpline
box = BoxSpline(lo=(0, 0, 0), hi=(G.V10_L, 2.0, 2.0), n_cells=(10, 2, 2), degree=2)
clamp = G.face_dofs(box, 0, 0); free = np.setdiff1d(np.arange(box.n_dofs), clamp)
K = G.kratos_box_stiffness(box, G.V10_E, 0.0, clamp); M = G.assemble_box_mass(box, G.V10_RHO)
f_exact = math.sqrt(G.V10_E / G.V10_RHO) / (4 * G.V10_L); t1 = 1 / f_exact
dt, n = 2e-6, 1000; amp = G.smooth_pulse(t1); f0 = G.face_traction_force(box, 0, 1, (1.0, 0, 0))
probe = np.array([3 * int(c) for c in box.boundary_control_points(0, 1)])
force = lambda t: amp(t) * f0
# regression: alpha_f = 0 == Bossak
for ab in (-0.05, -0.3):
    b = solve_bossak_linear(M, None, K, force, BossakConfig(dt=dt, n_steps=n, alpha_b=ab), free_dofs=free)
    g = solve_generalized_alpha_linear(M, None, K, force, GeneralizedAlphaConfig(dt=dt, n_steps=n, alpha_m=ab, alpha_f=0.0), free_dofs=free)
    print("bossak parity alpha", ab, "max|du|/max|u| =", float(np.abs(b.u - g.u).max() / np.abs(b.u).max()))
out = {}
CF = {"TRANSIENT_FIDELITY": (24766.8, 24768.6, 24745.6), "MODERATE_DISSIPATION": (24702.1, 24703.8, 24681.0)}
for name in ("TRANSIENT_FIDELITY", "MODERATE_DISSIPATION"):
    cfg = GeneralizedAlphaConfig.abaqus(name, dt, n)
    r = solve_generalized_alpha_linear(M, None, K, force, cfg, free_dofs=free)
    w = r.t >= 2 * t1
    fit = G.fit_free_vibration(r.t[w], r.u[:, probe].mean(axis=1)[w], f_exact)
    cf = CF[name]
    out[name] = {"alpha_hht": G and cfg.alpha_f * -1, "beta": cfg.beta, "gamma": cfg.gamma, "f1_hz": fit.frequency_hz,
                 "rel_vs_25000": fit.frequency_hz / f_exact - 1, "rel_vs_coreform_immersed": fit.frequency_hz / cf[0] - 1,
                 "rel_vs_coreform_bodyfit": fit.frequency_hz / cf[1] - 1, "decay": fit.decay_rate}
    print(name, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in out[name].items()})
json.dump(out, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v10_hht.json", "w"), indent=1)
