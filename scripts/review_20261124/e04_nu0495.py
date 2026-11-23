import sys, json, time
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
from tests.unit.test_fbar_mean_dilatation import _cook_case
import wht_iga.analyses.fbar as FB
NU = 0.495
OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/e04_nu0495.json"
cases = [("C1_p2_16", 16, 2, None), ("C1_p2_32", 32, 2, None), ("C0_p2_16", 16, 2, 0), ("p1_32", 32, 1, None),
         ("C0_p2_32", 32, 2, 0), ("C1_p2_64", 64, 2, None), ("p1_64", 64, 1, None), ("C0_p2_64", 64, 2, 0)]
out = {}
for name, n, p, cont in cases:
    t0 = time.perf_counter()
    model, f, fixed, tips = _cook_case(n, p, cont, 240.565, NU, 6.25)
    u = FB.solve_static(model, f, fixed, n_steps=10)
    uy = float(np.mean([u[3 * I + 1] for I in tips]))
    out[name] = {"n": n, "p": p, "cont": cont, "nu": NU, "uy": uy, "t": time.perf_counter() - t0, "ndof": int(model["n_dof"])}
    print(name, f"uy={uy:.6f} ndof={model['n_dof']} t={out[name]['t']:.1f}s", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
print("DONE", flush=True)
