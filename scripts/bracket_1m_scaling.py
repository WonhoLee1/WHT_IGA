"""Phase 6 1M-DOF gate measurement: Generator_bracket.step, PCA-OBB background grid (11,7,4)*k,
p=2, q=2, theta=0.2 extension, SmallDisplacement linear elastic, prescribed end displacement (S45 setup).

Usage (cwd WHT_IGA): python -u run_1m.py <k> <solver>[,<solver>...] [min_avail_gb]
  solvers: PardisoLU, AMGCL_CG_ILU0, AMGCL_GMRES_ILU0 (S45 _make_amgcl_or_pardiso, tol 1e-10)
Prints one JSON line 'RESULT {...}' per solver with stage times (s), RSS/peak working set (GB),
n_eq, nnz, reaction balance and (if several solvers) the relative difference to the first one.
"""
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import psutil

from tests.kratos_spikes.test_p6_industrial_bracket_amgcl import (
    _RESOURCES,
    _build_trimmed_and_classify,
    _make_amgcl_or_pardiso,
    _populate_kratos_model,
)
from wht_iga.analyses.methods import compute_pca_obb
from wht_iga.geometry import cad
from wht_iga.geometry.occ_io import load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference.box import BoxSpline

GB = 1024.0 ** 3
P = psutil.Process()


def avail_gb():
    return psutil.virtual_memory().available / GB


def mem():
    mi = P.memory_info()
    return mi.rss / GB, getattr(mi, "peak_wset", mi.rss) / GB


k = float(sys.argv[1])
solvers = sys.argv[2].split(",")
min_avail = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
T = {}
t0 = time.perf_counter()
geom = load_cad_file(_RESOURCES / "Generator_bracket.step", heal=False)
bvol = cad.brep_volume(geom.shape)
lo, hi = cad.shape_bbox(geom.shape)
diag = float(np.linalg.norm(hi - lo))
tess = cad.tessellate(geom.shape, deflection=0.003 * diag)
frame, lo_obb, hi_obb = compute_pca_obb(tess.nodes, padding_frac=0.05)
pts_local = frame.to_local_points(tess.nodes)
n_cells = tuple(max(1, int(round(c * k))) for c in (11, 7, 4))
box = BoxSpline(lo=lo_obb, hi=hi_obb, n_cells=n_cells, degree=2)
T["cad_tess_s"] = time.perf_counter() - t0
print(f"k={k} grid={n_cells} cells={int(np.prod(n_cells))} n_ctrl={box.n_ctrl} tris={len(tess.tris)} "
      f"avail={avail_gb():.1f}GB", flush=True)

t1 = time.perf_counter()
prep = _build_trimmed_and_classify(pts_local, tess.tris, box, bvol, workers=8)
T["trim_classify_ext_s"] = time.perf_counter() - t1
rss, peak = mem()
print(f"trim: active_cells={prep['n_active_cells']} cut={prep['n_cut']} fb={prep['n_fallback']} "
      f"vol_err={prep['vol_err'] * 100:.4f}% active_dofs={prep['active_dofs']} "
      f"(trim {prep['build_trimmed_table_s']:.1f}s cls {prep['classify_basis_s']:.1f}s ext {prep['build_extension_s']:.1f}s) "
      f"rss={rss:.2f}GB avail={avail_gb():.1f}GB", flush=True)

mods = load_kratos()
KM, SM, LA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["LinearSolversApplication"]
d_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
r_vars = (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z)
u_first = None
for skey in solvers:
    if avail_gb() < min_avail:
        print(f"ABORT before {skey}: available {avail_gb():.1f} GB < {min_avail} GB", flush=True)
        break
    gc.collect()
    tb = time.perf_counter()
    built = _populate_kratos_model(KM, SM, f"GB1M_{skey}", box, prep["qp_table"], prep["out_cells"],
                                   prep["cell_records"], axis_idx=0)
    t_build = time.perf_counter() - tb
    mp = built["mp"]
    rss_b, _ = mem()
    print(f"[{skey}] kratos model: elems={built['n_elements']} ties={built['n_ties']} "
          f"(vol {built['make_vol_s']:.1f}s elem {built['attach_elem_s']:.1f}s ties {built['apply_ties_s']:.1f}s) "
          f"rss={rss_b:.2f}GB avail={avail_gb():.1f}GB", flush=True)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = _make_amgcl_or_pardiso(KM, LA, skey)
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, True, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    ts = time.perf_counter()
    strategy.Solve()
    t_solve = time.perf_counter() - ts
    A = strategy.GetSystemMatrix()
    n_eq, nnz = int(A.Size1()), int(A.nnz()) if hasattr(A, "nnz") else -1
    U = np.array([[mp.GetNode(i + 1).GetSolutionStepValue(c) for c in d_vars] for i in range(box.n_ctrl)])
    R = np.array([[mp.GetNode(i + 1).GetSolutionStepValue(c) for c in r_vars] for i in range(box.n_ctrl)])
    U, R = frame.to_global_vectors(U), frame.to_global_vectors(R)
    finite = bool(np.all(np.isfinite(U)) and np.all(np.isfinite(R)))
    rx_bal = float(np.linalg.norm(R.sum(axis=0)) / max(float(np.abs(R).sum()), 1e-300))
    diff = None
    if u_first is None:
        u_first = U.copy()
    else:
        diff = float(np.linalg.norm(U - u_first) / max(float(np.linalg.norm(u_first)), 1e-300))
    rss_a, peak_a = mem()
    out = {"k": k, "grid": n_cells, "solver": skey, "active_dofs": prep["active_dofs"], "n_eq": n_eq, "nnz": nnz,
           "n_active_cells": prep["n_active_cells"], "n_cut": prep["n_cut"], "n_fallback": prep["n_fallback"],
           "vol_err": prep["vol_err"], "cad_tess_s": T["cad_tess_s"], "trim_classify_ext_s": T["trim_classify_ext_s"],
           "kratos_build_s": t_build, "strategy_solve_s": t_solve,
           "total_s": T["cad_tess_s"] + T["trim_classify_ext_s"] + t_build + t_solve,
           "rss_gb": rss_a, "peak_wset_gb": peak_a, "finite": finite, "reaction_balance_rel": rx_bal,
           "rel_diff_to_first": diff, "max_u": float(np.max(np.abs(U)))}
    print("RESULT " + json.dumps(out), flush=True)
    del strategy, builder, solver, scheme, built, mp, A, U, R
    gc.collect()
print(f"TOTAL {time.perf_counter() - t0:.1f}s", flush=True)
