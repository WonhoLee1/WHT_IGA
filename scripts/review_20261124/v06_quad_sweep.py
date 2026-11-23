import sys, json, inspect
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import numpy as np
import v06_penalty_sweep as PS
from wht_iga.solver import trimmed as T
import tests.kratos_spikes.test_p2mesh_v06geo as GEO
from tests.kratos_spikes.test_p2mesh_v06geo import sphere_octant_tris, A, B
res = []
for q in (1, 2, 3):
    orig = T.build_trimmed_table
    T.build_trimmed_table = lambda nodes, tris, box, q=None, _o=orig, _q=q, **kw: _o(nodes, tris, box, q=_q, **kw)
    S = PS.setup(10, 1)
    T.build_trimmed_table = orig
    e = PS.solve(S)
    row = {"q": q, "n_qp": int(len(S["W"])), "vol": float(S["W"].sum()), **e}
    print(row, flush=True); res.append(row)
json.dump(res, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_quad_sweep.json", "w"), indent=1)
