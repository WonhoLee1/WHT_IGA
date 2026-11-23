import sys, csv, json, time
import os, dataclasses
VAR = os.environ.get("VARIANT", "base")
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import wht_iga.analyses.thermal_plasticity as TP
import tests.kratos_spikes.test_p5_e03_p2_coupling as T
G = r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal/glass-to-metal-seal-convergence-history.csv"
rows = [r for r in csv.DictReader(open(G, encoding="utf-8")) if r["accepted"] == "True"]
levels = [float(r["total_time"]) for r in rows]
steps_count = {}
for r in rows:
    steps_count[r["step_name"]] = steps_count.get(r["step_name"], 0) + 1
print("CSV accepted increments per step:", steps_count, "levels:", len(levels), flush=True)
TP.thermal_levels = lambda n_steps, inc: list(levels)  # Coreform increment history (accepted total_time)
TOLREL = float(os.environ.get("TOLREL", "1e-8"))
TOLABS = float(os.environ.get("TOLABS", "1e-10"))
TAG = os.environ.get("TAG", "tol")
ITERS = []
class _IterProxy:
    def __init__(self, s, KM, mp):
        self._s, self._KM, self._mp = s, KM, mp
    def __getattr__(self, n):
        return getattr(self._s, n)
    def SolveSolutionStep(self):
        ok = self._s.SolveSolutionStep()
        ITERS.append((int(self._mp.ProcessInfo[self._KM.NL_ITERATION_NUMBER]), bool(ok)))
        return ok
_orig = T._make_pardiso_newton_strategy
T._make_pardiso_newton_strategy = lambda KM, LA, mp, **kw: _IterProxy(_orig(KM, LA, mp, **kw), KM, mp)
bundle_fn = T.e03_p2_bundle.__wrapped__
t0 = time.perf_counter()
bundle = bundle_fn()
print(f"bundle {time.perf_counter()-t0:.1f}s", flush=True)
mats = TP.e03_materials(T._DATA_DIR)
if VAR == "t20":
    h = mats["SS304L_header"].hardening
    st = list(h.stresses)
    st[0] = st[1]
    mats["SS304L_header"] = dataclasses.replace(mats["SS304L_header"], hardening=TP.HardeningTable(h.plastic_strains, tuple(st), h.temperatures))
t1 = time.perf_counter()
res = T._build_and_run_e03_p2(bundle, [mats["SS304L_header"], mats["S8061_glass"], mats["Alloy52_pins"]],
                              increments=[15, 15, 15], bc_mode="rbe3", tol_rel=TOLREL, tol_abs=TOLABS)
out = res["out"]
summary = {"history": out.get("history"), "steps": out.get("steps"), "n_solves": out.get("n_solves"), "level": out.get("level"),
           "n_levels": len(levels), "solve_s": time.perf_counter() - t1, "peak_rss_mb": res.get("peak_rss_mb"),
           "iters": ITERS, "tol_rel": TOLREL, "tol_abs": TOLABS, "snapshots_full": res["snapshots"],
           }
print(json.dumps(summary, indent=1, default=str), flush=True)
json.dump(summary, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/e03_" + TAG + ".json", "w"), indent=1, default=str)
print("DONE", flush=True)
