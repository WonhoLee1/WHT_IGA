# -*- coding: utf-8 -*-
"""V01-ext: prescribed-displacement cube with extension ON (G2 gate, v2.1 §9).

Reference: analytic u = (x, 2y, 3z), E = 1, nu = 0, unit cube, penalty BC-A on
all 6 faces (same spec as tests/kratos_spikes/test_p1_pipeline_cube.py).
Gate: max relative error <= 1e-5 for p = 1, 2, 3 with extension enabled.

Two cases per p:
  A "noop":      untrimmed cube (all eta = 1). Extension must create 0 ties and
                 reproduce the extension-OFF solution (|dU| <= 1e-12).
  B "synthetic": eta of a corner cell column set small (ETA_SMALL) so the support
                 classification marks slave basis functions and real MPC ties are
                 built. Integration is still the full cube (no qp_table), so the
                 affine exact solution stays in the constrained space; error must
                 stay <= 1e-5. n_slave must be > 0 (otherwise the case is void and
                 the script fails instead of passing by construction).

Env: WHT_V01_EXT_THETA (default 0.2), WHT_V01_EXT_EVID (evidence dir).
Writes <EVID>/results.json.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

EVID = os.environ.get(
    "WHT_V01_EXT_EVID",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs",
                 "phase2", "evidence", "v01_ext"))
THETA = float(os.environ.get("WHT_V01_EXT_THETA", "0.2"))
AFFINE_A = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
TOL = 1e-5
NOOP_TOL = 1e-12
N_NOOP = 2       # same grid as Phase-1 V01
N_SYN = 4        # finer grid so a small-cell column leaves other cells large
ETA_SMALL = 0.02


def _spec(p: int, n: int, name: str):
    from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec

    bcs = [DirichletBC(face=(a, s), field=AffineField(A=AFFINE_A, b=(0.0, 0.0, 0.0)))
           for a in range(3) for s in (0, 1)]
    return JobSpec(name=name, lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(n, n, n), degree=p, E=1.0, nu=0.0, bcs=bcs,
                   output_dir=os.path.join(EVID, ""))


def _exact(p: int, n: int) -> np.ndarray:
    from wht_iga.reference.box import BoxSpline

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, n), degree=p)
    Xc = box.control_points()
    return np.column_stack([Xc[:, 0], 2 * Xc[:, 1], 3 * Xc[:, 2]]).ravel()


def _rel(U: np.ndarray, Uex: np.ndarray) -> float:
    return float(np.abs(U - Uex).max() / np.abs(Uex).max())


def _records(n: int, small: set[tuple[int, int, int]]) -> list[dict]:
    return [{"i": i, "j": j, "k": k, "eta": ETA_SMALL if (i, j, k) in small else 1.0}
            for k in range(n) for j in range(n) for i in range(n)]


def run_p(p: int) -> dict:
    from wht_iga.solver.pipeline import run_linear

    # A: untrimmed no-op
    off = run_linear(_spec(p, N_NOOP, f"v01ext_p{p}_off"))
    on = run_linear(_spec(p, N_NOOP, f"v01ext_p{p}_noop"),
                    cell_records=_records(N_NOOP, set()), extension={"theta": THETA})
    Uex = _exact(p, N_NOOP)
    a = {"n": N_NOOP, "rel_err_off": _rel(off.U, Uex), "rel_err_on": _rel(on.U, Uex),
         "max_dU_on_off": float(np.abs(on.U - off.U).max()), "tie_info": on.tie_info}

    # B: synthetic small-cell column (i = n-1, j = n-1, all k)
    small = {(N_SYN - 1, N_SYN - 1, k) for k in range(N_SYN)}
    syn = run_linear(_spec(p, N_SYN, f"v01ext_p{p}_syn"),
                     cell_records=_records(N_SYN, small), extension={"theta": THETA})
    b = {"n": N_SYN, "eta_small": ETA_SMALL, "small_cells": sorted(small),
         "rel_err": _rel(syn.U, _exact(p, N_SYN)), "beta": syn.beta_used,
         "tie_info": syn.tie_info}

    n_slave_b = int((b["tie_info"] or {}).get("n_slave", 0))
    n_ties_a = int((a["tie_info"] or {}).get("n_ties", 0))
    ok = (a["rel_err_on"] <= TOL and a["max_dU_on_off"] <= NOOP_TOL and n_ties_a == 0
          and b["rel_err"] <= TOL and n_slave_b > 0)
    print(f"[V01-ext p={p}] noop: err_off={a['rel_err_off']:.3e} err_on={a['rel_err_on']:.3e} "
          f"|dU|={a['max_dU_on_off']:.3e} ties={n_ties_a} | syn: err={b['rel_err']:.3e} "
          f"ties={b['tie_info']} -> {'PASS' if ok else 'FAIL'}", flush=True)
    return {"p": p, "noop": a, "synthetic": b, "pass": bool(ok)}


def main() -> None:
    os.makedirs(EVID, exist_ok=True)
    rows = [run_p(p) for p in (1, 2, 3)]
    out = {"theta": THETA, "rows": rows,
           "gate": {"rel_err_max": TOL, "noop_dU_max": NOOP_TOL, "synthetic_n_slave_min": 1},
           "all_pass": all(r["pass"] for r in rows)}
    with open(os.path.join(EVID, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=str)
    print(f"[V01-ext] all_pass={out['all_pass']} -> {os.path.join(EVID, 'results.json')}")
    if not out["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
