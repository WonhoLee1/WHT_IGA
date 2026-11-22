import json, sys, time
from pathlib import Path
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
from wht_iga.solver.pipeline import run_linear
from wht_iga.solver.validate import DirichletBC, JobSpec, TractionLoad
E, NU = 2.0e11, 0.0
B = H = 0.1; L = 1.0; IX = B * H ** 3 / 12.0
P_N = 500.0 * B * B  # 500 Pa on 0.1 x 0.1 end face -> 5 N
EB_TEXT = P_N * L ** 3 / (3 * E * IX)
OUT = Path(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/report/v03_spec")
OUT.mkdir(parents=True, exist_ok=True)
res_all = {}
for p in (1, 2, 3):
    for k in range(5):
        n = (2 ** k, 2 ** k, 4 * 2 ** k)
        spec = JobSpec(name=f"v03sw_p{p}_k{k}", output_dir=str(OUT),
                       lo=(-0.05, -0.05, 0.0), hi=(0.05, 0.05, 1.0), n_cells=n, degree=p,
                       E=E, nu=NU, fixed_faces=[(2, 0)],
                       tractions=[TractionLoad(face=(2, 1), vector=(0.0, -500.0, 0.0))],
                       probes=[{"name": "tip", "xyz": [0.0, 0.0, 1.0]}])
        t0 = time.perf_counter()
        r = run_linear(spec)
        uy = float(np.atleast_1d(r.probe_values["tip"])[1])
        res_all[f"p{p}_k{k}"] = {"p": p, "k": k, "n": list(n), "h_z": 1.0 / n[2], "tip_uy": uy,
                                 "eb_text": -EB_TEXT, "time_s": time.perf_counter() - t0}
        print(f"p={p} k={k} n={n} tip={uy:.6e} eb_text={-EB_TEXT:.6e} t={time.perf_counter()-t0:.1f}s", flush=True)
json.dump(res_all, open(OUT / "v03_sweep_spec.json", "w"), indent=1)
print("DONE")
