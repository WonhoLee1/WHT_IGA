import sys, json, time
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
from tests.unit.test_fbar_mean_dilatation import _cook_case
import wht_iga.analyses.fbar as FB
from wht_iga.reference.cook import U_Y_REF_MM
out = {}
for name, n, p, cont in [("p3C1_8", 8, 3, 1), ("p3C1_16", 16, 3, 1), ("p3C1_32", 32, 3, 1), ("p2C1_16_ref", 16, 2, None)]:
    t0 = time.perf_counter()
    model, f, fixed, tips = _cook_case(n, p, cont, 240.565, 0.49999, 6.25)
    u = FB.solve_static(model, f, fixed, n_steps=10)
    uy = float(np.mean([u[3 * I + 1] for I in tips]))
    out[name] = {"n": n, "p": p, "cont": cont, "uy": uy, "rel_pct": 100 * (uy - U_Y_REF_MM) / U_Y_REF_MM, "ndof": int(model["n_dof"]), "t": time.perf_counter() - t0}
    print(name, f"uy={uy:.6f} rel={out[name]['rel_pct']:+.3f}% ndof={model['n_dof']} t={out[name]['t']:.1f}s", flush=True)
    json.dump(out, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/e04_p3c1.json", "w"), indent=1)
print("DONE", flush=True)
