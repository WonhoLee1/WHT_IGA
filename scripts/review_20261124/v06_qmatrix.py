import sys, json, time
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import v06_penalty_sweep as PS
from wht_iga.solver import trimmed as T
CF = {(1, 0): (0.00584, 0.00568, 0.285, 0.0189), (1, 1): (0.00199, 0.00194, 0.155, 0.0093),
      (2, 0): (0.00094, 0.00093, 0.0185, 0.0016), (2, 1): (0.000173, 0.000166, 0.00681, 0.000445),
      (3, 0): (0.000317, 0.000287, 0.00837, 0.000859)}
OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_qmatrix.json"
res = []
orig = T.build_trimmed_table
for p, lvl, qs in [(1, 0, (1, 2)), (2, 0, (1, 2, 3, 4)), (1, 1, (1, 2)), (3, 0, (2, 3, 4)), (2, 1, (2, 3, 4))]:
    n = 10 * 2 ** lvl
    for q in qs:
        t0 = time.perf_counter()
        T.build_trimmed_table = lambda nodes, tris, box, q=None, _o=orig, _q=q, **kw: _o(nodes, tris, box, q=_q, **kw)
        try:
            S = PS.setup(n, p)
            e = PS.solve(S)
        finally:
            T.build_trimmed_table = orig
        cf = CF[(p, lvl)]
        row = {"p": p, "level": lvl, "q": q, **e, "coreform": dict(zip(("u", "ur", "srr", "stt"), cf)),
               "t": time.perf_counter() - t0}
        print(f"p={p} L{lvl} q={q}: u={e['u']:.3e} ur={e['ur']:.3e} srr={e['srr']:.3e} stt={e['stt']:.3e} | CF u={cf[0]:.3e} srr={cf[2]:.3e} stt={cf[3]:.3e}  ({row['t']:.0f}s)", flush=True)
        res.append(row)
        json.dump(res, open(OUT, "w"), indent=1)
print("DONE")
