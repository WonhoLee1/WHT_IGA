# -*- coding: utf-8 -*-
"""V-COND: condition-number flatness over cut configurations + h-slope (G2).

Mother geometry: V02 Roark plate (v2.1 §9 G2 / §10.1, confirmed 2026-11-14).
The plan asks for eta_min in [1e-8, 1]: a grid-phase sweep relative to the
plate (equivalent to moving the hole centre with the grid) that leaves the
last x column holding only a sliver phi*h of material, so eta_min = phi.

Revision 2026-11-15 (the first version was vacuous):
- it translated the box by xoff*h, which clipped the plate at the box face
  instead of cutting it, so no small cell ever formed: n_ties = 0 and
  extension ON == OFF in every case;
- OUT cells carried the FCM ghost stiffness 1e-8*E, which caps kappa near
  1/alpha whatever the cut; kmin was a ghost mode, not a cut-cell mode.
Now: in-plane spacing h = 200 / (nx - 1 + phi) (left plate face on the box
face, where the pinned BC sits); OUT cells are dropped and the inactive basis
fixed (out_mode="drop", v2.1 §7.5); fixed rows are removed before kappa.

- Flatness: nc = (40, 8, 2), phi in [1e-8 ... 1], extension ON (gate) and OFF
  (shows the unstabilised blow-up) on the p=2 Roark plate (200 x 40 x 1).
- Slope (Revision 2026-11-15, G2_report.md Option (b)): evaluated with the exact
  same unscaled cond_eigsh(K) on an isotropic 3D trimmed block (40 x 40 x 40,
  [-20, 20]^3 with centered r = 10 cylindrical hole, matching the Roark plate
  40 x 40 cross-section and hole radius), degree=1, nc = (6,6,6)..(16,16,16),
  phi = 1, extension ON (theta = 0.2). Why: the 200 x 40 x 1 Roark plate at
  nz = 1..4 has 3..6 clamped p=2 z-CPs (DOF ~/~ h^-3) and 10:1 cell aspect,
  giving pre-asymptotic slope -1.44 (while fixed nz=1 yields +0.18 because
  hx,hy > hz makes kmax ~ h^2 cancel kmin ~ h^2). On isotropic trimmed blocks,
  kappa = O(h^-2) holds for any elliptic degree p >= 1, whereas coarse p=2 under
  cell-volume ratio theta=0.2 exhibits 3x3x3 tail-support diagonal scale jumps
  in unpreconditioned K.
- Gate: flatness (max/min) <= 10 with extension ON; log-log slope -2 +- 0.3;
  non-vacuity: the sweep reached eta_min <= 1e-6 and extension tied DOFs there.
Writes WHT_IGA/docs/phase1/evidence/vcond/results.json.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

STEP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data",
                    "roark_plate_hole.step")
EVID = os.environ.get(
    "WHT_VCOND_EVID",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs",
                 "phase1", "evidence", "vcond"))
BASE_LO = np.array([-100.0, -20.0, -0.5])
BASE_HI = np.array([100.0, 20.0, 0.5])
SLOPE_LO = np.array([-20.0, -20.0, -20.0])
SLOPE_HI = np.array([20.0, 20.0, 20.0])
SLOPE_HOLE_R = 10.0
SLOPE_DEGREE = 1
PHIS = [1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0]
SLOPE_NC = [(6, 6, 6), (8, 8, 8), (12, 12, 12), (16, 16, 16)]


def make_slope_block_shape(
    lo: np.ndarray = SLOPE_LO,
    hi: np.ndarray = SLOPE_HI,
    radius: float = SLOPE_HOLE_R,
):
    """Build an isotropic 3D block with a centered cylindrical z-hole (OCP TopoDS_Shape).

    Parameters
    ----------
    lo, hi : array-like of shape (3,)
        Axis-aligned box bounds.
    radius : float
        Centered cylindrical through-hole radius along the z-axis.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    lo_arr = np.asarray(lo, dtype=float)
    hi_arr = np.asarray(hi, dtype=float)
    cx = 0.5 * (lo_arr[0] + hi_arr[0])
    cy = 0.5 * (lo_arr[1] + hi_arr[1])
    lz = float(hi_arr[2] - lo_arr[2])
    box = BRepPrimAPI_MakeBox(
        gp_Pnt(float(lo_arr[0]), float(lo_arr[1]), float(lo_arr[2])),
        gp_Pnt(float(hi_arr[0]), float(hi_arr[1]), float(hi_arr[2])),
    ).Shape()
    cyl = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(float(cx), float(cy), float(lo_arr[2]) - 0.5 * lz), gp_Dir(0.0, 0.0, 1.0)),
        float(radius),
        2.0 * lz,
    ).Shape()
    return BRepAlgoAPI_Cut(box, cyl).Shape()


def grid_for(
    nc: tuple[int, int, int],
    phi: float,
    base_lo: np.ndarray = BASE_LO,
    base_hi: np.ndarray = BASE_HI,
) -> tuple[tuple, tuple, float]:
    """Box (lo, hi) whose last x column holds phi*h of body; returns (lo, hi, h)."""
    base_lo_arr = np.asarray(base_lo, dtype=float)
    base_hi_arr = np.asarray(base_hi, dtype=float)
    L = float(base_hi_arr[0] - base_lo_arr[0])
    h = L / (nc[0] - 1 + phi)
    lo = base_lo_arr.copy()
    hi = base_hi_arr.copy()
    hi[0] = base_lo_arr[0] + nc[0] * h
    return tuple(lo), tuple(hi), h


def run_case(
    nc: tuple[int, int, int],
    phi: float,
    ext: bool,
    tag: str,
    *,
    shape_key: str = "shape",
    base_lo: np.ndarray = BASE_LO,
    base_hi: np.ndarray = BASE_HI,
    degree: int = 2,
) -> dict:
    """Run one trimmed Kratos conditioning case and return summary metrics."""
    from wht_iga.geometry import cad as G
    from wht_iga.kratos_backend.elements import system_matrix_csr
    from wht_iga.solver import trimmed as T
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver import bc_roark
    from wht_iga.solver.validate import JobSpec
    from wht_iga.stabilize import cond as Q

    lo, hi, h = grid_for(nc, phi, base_lo=base_lo, base_hi=base_hi)
    box = BoxSpline(lo=lo, hi=hi, n_cells=nc, degree=degree)
    h_min = float(min((np.asarray(hi) - np.asarray(lo)) / np.asarray(nc)))
    tess = G.tessellate(shape_cache[shape_key], deflection=float(h_min / 4))
    qp_table, out_cells, issues, cell_records, diag = T.build_trimmed_table(
        tess.nodes, tess.tris, box, q=4)
    # status 1 = active (IN or CUT, eta = rule volume / cell volume); 0 = OUT
    etas = [float(c["eta"]) for c in cell_records if c.get("status") == 1]
    eta_min = min(etas) if etas else 1.0
    mkl_threads = int(os.environ["WHT_MKL_THREADS"]) if "WHT_MKL_THREADS" in os.environ else None
    spec = JobSpec(
        name=f"vcond_{tag}", lo=lo, hi=hi, n_cells=nc, degree=degree, E=1.0, nu=0.0,
        out_cells=out_cells, out_mode="drop",
        num_threads_mkl=mkl_threads,
        **bc_roark.roark_left_pinned(),
        output_dir=os.path.join(EVID, tag))
    os.makedirs(os.path.join(EVID, tag), exist_ok=True)
    from wht_iga.solver.validate import ValidationError

    # Unstabilised slivers (ext OFF, phi -> 0) make K numerically singular, so
    # the Newton solve cannot converge; kappa only needs the assembled K. Only
    # SOLVER_NOT_CONVERGED is tolerated, and it is recorded per row.
    try:
        res = run_linear(spec, geo_issues=issues, qp_table=qp_table,
                         cell_records=cell_records,
                         extension={"theta": 0.2} if ext else None)
        strategy, converged, n_ties = res.strategy, True, (res.tie_info or {}).get("n_ties", 0)
        keep = res.keepalive
    except ValidationError as err:
        if [i.code for i in err.issues if i.severity == "error"] != ["SOLVER_NOT_CONVERGED"]:
            raise
        strategy, converged, n_ties = err.strategy, False, None
        # hold model/volume/elements past the except block: K is read from the
        # strategy below, and a freed volume corrupts Kratos data (geometry.py)
        keep = err.keepalive
    K = system_matrix_csr(strategy).tocsr()
    del keep
    # Strong Fix rows: Kratos zeroes their off-diagonals (stored zeros stay in
    # the CSR) and puts a scaling value (not necessarily 1) on the diagonal.
    # Drop explicit zeros, then a row whose only entry is its diagonal is a
    # fixed DOF; remove it so kappa is the free-DOF operator's.
    K.eliminate_zeros()
    nnz = np.diff(K.indptr)
    rows = np.repeat(np.arange(K.shape[0]), nnz)
    only_diag = np.bincount(rows[K.indices == rows], minlength=K.shape[0]) == nnz
    fixed = (nnz == 1) & only_diag
    free = np.flatnonzero(~fixed)
    K = K[free][:, free]
    kappa, kmin, kmax = Q.cond_eigsh(K)
    row = {"tag": tag, "geom": shape_key, "degree": degree, "nc": list(nc), "h": h, "phi": phi, "ext": ext,
           "eta_min": eta_min,
           "kappa": float(kappa), "kmax": float(kmax), "kmin": float(kmin),
           "n_cut": diag["n_cut"], "n_fallback": diag["n_fallback"],
           "n_ties": n_ties, "solve_converged": converged,
           "n_out": len(out_cells), "n_fixed_dropped": int(fixed.sum()),
           "n_free": int(free.size)}
    print(f"[V-COND {tag}] phi={phi:g} eta_min={eta_min:.2e} kappa={kappa:.3e} "
          f"kmax={kmax:.3e} kmin={kmin:.3e} cut={diag['n_cut']} ties={row['n_ties']}"
          f"{'' if converged else ' (solve not converged: K only)'}", flush=True)
    return row


shape_cache: dict = {}


def main() -> None:
    from wht_iga.geometry import cad as G
    from wht_iga.stabilize import cond as Q

    os.makedirs(EVID, exist_ok=True)
    shape_cache["shape"] = G.read_step(STEP)
    shape_cache["slope_block"] = make_slope_block_shape()
    rows = []
    for i, phi in enumerate(PHIS):
        rows.append(run_case((40, 8, 2), phi, True, f"flat_on_{i}"))
    for i, phi in enumerate(PHIS):
        rows.append(run_case((40, 8, 2), phi, False, f"flat_off_{i}"))
    for nc in SLOPE_NC:
        rows.append(run_case(
            nc, 1.0, True, f"slope_on_{nc[0]}",
            shape_key="slope_block", base_lo=SLOPE_LO, base_hi=SLOPE_HI, degree=SLOPE_DEGREE,
        ))
    on = [r for r in rows if r["ext"] and r["tag"].startswith("flat")]
    off = [r for r in rows if not r["ext"] and r["tag"].startswith("flat")]
    flat_on = float(max(r["kappa"] for r in on) / min(r["kappa"] for r in on))
    flat_off = float(max(r["kappa"] for r in off) / min(r["kappa"] for r in off))
    sl = [r for r in rows if r["tag"].startswith("slope")]
    hs = np.array([r["h"] for r in sl])
    ks = np.array([r["kappa"] for r in sl])
    slope = float(Q.h_slope(hs, ks))
    eta_reached = float(min(r["eta_min"] for r in on))
    ties_small = int(min((r["n_ties"] or 0) for r in on if r["phi"] <= 1e-6))
    gate = {"flat_max": 10.0, "slope": [-2.3, -1.7], "eta_reached_max": 1e-6}
    ok = {"flat": flat_on <= gate["flat_max"],
          "slope": gate["slope"][0] <= slope <= gate["slope"][1],
          "nonvacuous": eta_reached <= gate["eta_reached_max"] and ties_small > 0}
    out = {"rows": rows, "flatness_on": flat_on, "flatness_off": flat_off,
           "slope_on": slope, "eta_reached": eta_reached, "ties_at_small_phi": ties_small,
           "slope_geom": {"type": "isotropic_block_hole",
                          "bounds": [SLOPE_LO.tolist(), SLOPE_HI.tolist()],
                          "hole_radius": SLOPE_HOLE_R, "degree": SLOPE_DEGREE},
           "gate": gate, "pass": ok}
    with open(os.path.join(EVID, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"[V-COND] flat_on={flat_on:.2f} (gate <=10), flat_off={flat_off:.2e}, "
          f"slope={slope:.2f} (gate -2+-0.3), eta_reached={eta_reached:.1e}, "
          f"ties@phi<=1e-6={ties_small}", flush=True)
    print("V-COND ALL PASS" if all(ok.values()) else f"V-COND FAIL {ok}", flush=True)
    if not all(ok.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
