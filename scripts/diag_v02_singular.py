# -*- coding: utf-8 -*-
"""Diagnose the V02 level-5 "Decomposition failed!" (Kratos Eigen Pardiso hides the code).

Assembles the exact V02 system (same spec as scripts/v02_roark.py), captures
the Kratos strategy with a 1-iteration CG solver (so Solve does not throw),
extracts K, then factorizes it with pypardiso (mtype -2, LDLᵀ) to read the
real MKL error code and reports matrix health: zero/negative/tiny diagonals,
empty rows, symmetry, and which control points own the suspicious rows.

Usage: WHT_DIAG_LEVEL=5 python WHT_IGA/scripts/diag_v02_singular.py
Writes WHT_IGA/docs/phase1/evidence/v02_diag/diag_l{level}.json.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
EVID = os.path.join(HERE, "..", "docs", "phase1", "evidence", "v02_diag")


def main() -> None:
    import scipy.sparse as sp

    import v02_roark as V
    from wht_iga.geometry import cad as G
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.elements import system_matrix_csr
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver import bc_roark
    from wht_iga.solver import pipeline as PL
    from wht_iga.solver import trimmed as T
    from wht_iga.solver.validate import JobSpec, ValidationError

    lv = int(os.environ.get("WHT_DIAG_LEVEL", "5"))
    nc = V.LEVELS[lv]
    KM = load_kratos()["KratosMultiphysics"]
    shape = G.read_step(V.STEP)
    box = BoxSpline(lo=V.LO, hi=V.HI, n_cells=nc, degree=2)
    h = float(min((np.asarray(V.HI) - np.asarray(V.LO)) / np.asarray(nc)))
    tess = G.tessellate(shape, deflection=float(h / 4))
    qp, out, iss, rec, diag = T.build_trimmed_table(tess.nodes, tess.tris, box, q=4)
    print(f"[diag l{lv}] cut={diag['n_cut']} fb={diag['n_fallback']} out={len(out)}", flush=True)

    captured = {}
    real_strategy = KM.ResidualBasedNewtonRaphsonStrategy

    def _capture(*a, **k):
        s = real_strategy(*a, **k)
        captured["s"] = s
        return s

    KM.ResidualBasedNewtonRaphsonStrategy = _capture
    PL._make_linear_solver = lambda KM_, name, thr=None: KM_.CGSolver(1e-30, 1)
    spec = JobSpec(name=f"diag_l{lv}", lo=V.LO, hi=V.HI, n_cells=nc, degree=2, E=1.0, nu=0.0,
                   out_cells=out, **bc_roark.roark_left_pinned())
    try:
        PL.run_linear(spec, geo_issues=iss, qp_table=qp, cell_records=rec)
    except ValidationError:
        pass  # expected: 1-iteration CG does not converge; K is assembled
    finally:
        KM.ResidualBasedNewtonRaphsonStrategy = real_strategy
    K = system_matrix_csr(captured["s"]).tocsr()
    captured["s"].Clear()  # free Kratos' own copy of K before the Pardiso factorizations
    captured.clear()
    K.sort_indices()
    n = K.shape[0]
    d = K.diagonal()
    nnz_row = np.diff(K.indptr)
    asym = float(abs(K - K.T).max())
    kmax = float(abs(K).max())
    unit = (nnz_row == 1) & (d == 1.0)
    free = ~unit
    dfree = d[free]
    rep = {"level": lv, "n": int(n), "nnz": int(K.nnz), "n_fixed_unit_rows": int(unit.sum()),
           "empty_rows": int((nnz_row == 0).sum()), "max_abs_K": kmax,
           "max_asym": asym, "rel_asym": asym / kmax if kmax else 0.0,
           "diag_free_min": float(dfree.min()), "diag_free_max": float(dfree.max()),
           "n_diag_le_0": int((dfree <= 0).sum()),
           "n_diag_lt_1e-12_max": int((dfree < 1e-12 * dfree.max()).sum()),
           "n_diag_lt_1e-8_max": int((dfree < 1e-8 * dfree.max()).sum())}
    # rows with the smallest free diagonals -> control point ijk
    nu, nv = box.n_cp[0], box.n_cp[1]
    idx_free = np.flatnonzero(free)
    worst = idx_free[np.argsort(d[idx_free])[:10]]
    rep["worst_rows"] = [{"dof": int(r), "cp": int(r // 3), "comp": int(r % 3),
                          "ijk": [int((r // 3) % nu), int((r // 3 // nu) % nv), int((r // 3) // (nu * nv))],
                          "diag": float(d[r])} for r in worst]
    print(json.dumps({k: v for k, v in rep.items() if k != "worst_rows"}, indent=1), flush=True)

    import pypardiso

    for mtype, tag in ((-2, "ldlt"), (11, "lu")):
        ps = pypardiso.PyPardisoSolver(mtype=mtype)
        A = sp.triu(K, format="csr") if mtype == -2 else K
        A.sort_indices()
        def _iparm() -> dict:
            return {"perturbed_pivots": int(ps.get_iparm(14)),
                    "factor_nnz": int(ps.get_iparm(18)),
                    "peak_mem_kb_analysis": int(ps.get_iparm(15)),
                    "mem_kb_factor": int(ps.get_iparm(16)) + int(ps.get_iparm(17)),
                    "n_pos_eig": int(ps.get_iparm(22)), "n_neg_eig": int(ps.get_iparm(23))}

        stage = "analysis(11)"
        try:
            ps.set_phase(11)
            ps._call_pardiso(A, np.zeros((n, 1)))
            after_analysis = _iparm()  # keep: phase-11 memory/nnz estimates survive a phase-22 failure
            stage = "factorization(22)"
            ps.set_phase(22)
            ps._call_pardiso(A, np.zeros((n, 1)))
            rep[f"pardiso_{tag}"] = {"error": 0, **_iparm()}
        except pypardiso.pardiso_wrapper.PyPardisoError as e:
            rep[f"pardiso_{tag}"] = {"error": str(e), "failed_stage": stage,
                                     "after_analysis": locals().get("after_analysis")}
        finally:
            ps.free_memory(everything=True)
        print(f"[pardiso {tag}] {rep[f'pardiso_{tag}']}", flush=True)
    os.makedirs(EVID, exist_ok=True)
    with open(os.path.join(EVID, f"diag_l{lv}.json"), "w", encoding="utf-8") as f:
        json.dump(rep, f, indent=1)


if __name__ == "__main__":
    main()
