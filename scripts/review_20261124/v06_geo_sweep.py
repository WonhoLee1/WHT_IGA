import sys, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import v06_penalty_sweep as PS
import tests.kratos_spikes.test_p2mesh_v06geo as GEO
orig = GEO.sphere_octant_tris
res = []
for ng, comp in [(8, False), (12, False), (16, False), (24, False), (24, True), (48, False)]:
    PS.sphere_octant_tris = lambda a, b, n, _ng=ng, _c=comp: orig(a, b, _ng, compensate=_c)
    S = PS.setup(10, 1)
    e = PS.solve(S)
    row = {"n_geo": ng, "compensate": comp, "n_active_cells": None, "vol": float(S["W"].sum()), **e}
    print(row, flush=True); res.append(row)
json.dump(res, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_geo_sweep.json", "w"), indent=1)
