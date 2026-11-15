# -*- coding: utf-8 -*-
"""Pardiso settings benchmark on real V02 matrices (Kratos path + pypardiso path).

Part A (Kratos): run_linear on V02 levels with linear_solver in
  pardiso_lu / pardiso_ldlt / pardiso_llt and num_threads_mkl variants;
  reports setup_s / solve_s and max|U - U_lu| / max|U_lu| (must be <= 1e-10).
Part B (pypardiso): the assembled K of the largest level (from the Kratos
  strategy) solved with pypardiso.spsolve default vs wht_iga.linalg.PardisoFactor
  kinds / ordering / two-level; reports factor + solve time and relative residual.

Env: WHT_BENCH_LEVELS (default "3,4" = V02 h 2.5, 1.5). Writes
WHT_IGA/docs/phase2/evidence/bench_pardiso/results.json.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)

EVID = os.path.join(HERE, "..", "docs", "phase2", "evidence", "bench_pardiso")


def _v02_case(level: int):
    import v02_roark as V
    from wht_iga.geometry import cad as G
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver import bc_roark
    from wht_iga.solver import trimmed as T

    nc = V.LEVELS[level]
    shape = G.read_step(V.STEP)
    box = BoxSpline(lo=V.LO, hi=V.HI, n_cells=nc, degree=2)
    h = float(min((np.asarray(V.HI) - np.asarray(V.LO)) / np.asarray(nc)))
    tess = G.tessellate(shape, deflection=float(h / 4))
    qp_table, out_cells, issues, cell_records, _diag = T.build_trimmed_table(
        tess.nodes, tess.tris, box, q=4)
    return nc, V, qp_table, out_cells, issues, cell_records, bc_roark


def part_a(levels: list[int]) -> tuple[list[dict], object]:
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.validate import JobSpec

    rows, last_strategy_res = [], None
    variants = [("pardiso_lu", None), ("pardiso_ldlt", None), ("pardiso_llt", None),
                ("pardiso_llt", 8), ("pardiso_llt", 16)]
    for lv in levels:
        nc, V, qp, out, iss, rec, bc = _v02_case(lv)
        U_ref = None
        for name, thr in variants:
            spec = JobSpec(name=f"bench_l{lv}", lo=V.LO, hi=V.HI, n_cells=nc, degree=2,
                           E=1.0, nu=0.0, out_cells=out, linear_solver=name,
                           num_threads_mkl=thr, **bc.roark_left_pinned())
            t = time.perf_counter()
            res = run_linear(spec, geo_issues=iss, qp_table=qp, cell_records=rec)
            wall = time.perf_counter() - t
            if U_ref is None:
                U_ref = res.U.copy()
            dU = float(np.abs(res.U - U_ref).max() / np.abs(U_ref).max())
            row = {"level": lv, "n_cells": list(nc), "n_dofs": int(res.U.size),
                   "solver": name, "num_threads_mkl": thr, "wall_s": wall,
                   "setup_s": res.timings["setup_s"], "solve_s": res.timings["solve_s"],
                   "n_iter": res.n_iter, "rel_dU_vs_lu": dU}
            print(f"[A l{lv} {name} thr={thr}] dofs={row['n_dofs']} solve={row['solve_s']:.2f}s "
                  f"setup={row['setup_s']:.2f}s iter={res.n_iter} dU={dU:.1e}", flush=True)
            rows.append(row)
            last_strategy_res = res
    return rows, last_strategy_res


def part_b(res) -> list[dict]:
    import pypardiso

    from wht_iga.kratos_backend.elements import system_matrix_csr
    from wht_iga.linalg import PardisoFactor

    K = system_matrix_csr(res.strategy).tocsr()
    K.sort_indices()
    b = np.random.default_rng(0).standard_normal(K.shape[0])
    rows = []

    def rec(tag, fn):
        t = time.perf_counter()
        x = fn()
        dt = time.perf_counter() - t
        r = float(np.linalg.norm(K @ x - b) / np.linalg.norm(b))
        print(f"[B {tag}] n={K.shape[0]} nnz={K.nnz} {dt:.2f}s resid={r:.1e}", flush=True)
        rows.append({"tag": tag, "n": int(K.shape[0]), "nnz": int(K.nnz), "time_s": dt, "rel_resid": r})

    rec("pypardiso.spsolve default (mtype 11)", lambda: pypardiso.spsolve(K.copy(), b))
    for kind in ("general", "sym", "spd"):
        for par, two in ((False, False), (True, False), (True, True)):
            def run(kind=kind, par=par, two=two):
                f = PardisoFactor(K, kind=kind, parallel_ordering=par, two_level=two)
                try:
                    return f.solve(b)
                finally:
                    f.free()
            rec(f"PardisoFactor kind={kind} par_order={par} two_level={two}", run)
    return rows


def main() -> None:
    os.makedirs(EVID, exist_ok=True)
    levels = [int(v) for v in os.environ.get("WHT_BENCH_LEVELS", "3,4").split(",")]
    a, last = part_a(levels)
    b = part_b(last)
    with open(os.path.join(EVID, "results.json"), "w", encoding="utf-8") as f:
        json.dump({"kratos": a, "pypardiso": b}, f, indent=1)
    bad = [r for r in a if r["rel_dU_vs_lu"] > 1e-10] + [r for r in b if r["rel_resid"] > 1e-8]
    print(f"[bench] written; accuracy violations: {len(bad)}")
    if bad:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
