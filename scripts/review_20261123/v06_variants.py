import sys, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
from tests.kratos_spikes.test_p2mesh_v06 import run_v06
import contextlib, io
res = {}
variants = [
    ("base_penalty_q2_ng24", dict(q=None, sym_bc="penalty"), 24),
    ("sym_fix", dict(q=None, sym_bc="fix"), 24),
    ("q4", dict(q=4, sym_bc="penalty"), 24),
    ("ngeo96", dict(q=None, sym_bc="penalty"), 96),
    ("ext0_ngeo96", dict(q=None, sym_bc="penalty", ext_theta=0.0), 96),
]
for name, kw, ng in variants:
    with contextlib.redirect_stdout(io.StringIO()):
        r = run_v06((10, 10, 10), 1, ng, out_mode="drop", **kw)
    res[name] = {k: r[k] for k in ("u", "ur", "srr", "stt") if k in r}
    print(name, {k: f"{v:.4e}" for k, v in res[name].items()}, flush=True)
json.dump(res, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_variants.json", "w"), indent=1)
print("DONE", flush=True)
