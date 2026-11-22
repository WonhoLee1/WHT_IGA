import sys, time, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit")
import test_mixed_p1 as T
out = {}
for h in (0.5, 0.25, 0.125):
    if h == 0.125:
        t0 = time.perf_counter()
    r = T._kirsch(T._Plate(*T.LEVEL[h], degree=3, cont=0), 3, "P2")
    out[str(h)] = {k: (float(v) if isinstance(v, float) else v) for k, v in r.items()}
    print(f"h={h} ndof={r['ndof']} e_u={r['e_u']:.4e} e_s={r['e_s']:.4e} t={r['t']:.1f}s", flush=True)
hs = [0.5, 0.25, 0.125]; rs = [out[str(h)] for h in hs]
su = T._slope(rs[0], rs[2], 0.5, 0.125, "e_u"); ss = T._slope(rs[0], rs[2], 0.5, 0.125, "e_s")
out["slope_end_u"] = float(su); out["slope_end_sigma"] = float(ss)
out["slope_ols_sigma"] = float(T._ols(rs, hs, "e_s"))
print("slopes end-point u %.3f sigma %.3f OLS sigma %.3f" % (su, ss, out["slope_ols_sigma"]), flush=True)
json.dump(out, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v07_p3_standard.json", "w"), indent=1)
print("DONE", flush=True)
