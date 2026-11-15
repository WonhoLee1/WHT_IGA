# -*- coding: utf-8 -*-
"""V02: Roark 17p1-7a finite plate with hole (G1/G2 gate, Coreform convention).

Setup mirrors Coreform verification (200 x 40 x 1, r=10, E=1, nu=0,
P=1e-6 tension on +x face, left face pinned). Grid == plate bbox (padding 0):
outer faces conforming (box-face BC/load paths), hole interior trimmed via
GEO cut QP tables (qp_table hook).

Gate (user decision 2026-11-15, Coreform roark-17p1-7a §5/§7.5 convention):
- levels ``coreform`` (default): nested h ~ 6/2^k, k = 0..4 -> n_cells =
  (33, 7, 1) * 2^k in-plane (h = 6.06 x 5.71 / 2^k; nz = 1 since nu = 0 makes
  the solution z-independent). ``WHT_V02_SET=legacy`` = old 10/5/3/2.5/1.5/0.75.
- peak (max principal on the hole face) vs Roark 4.3175e-6 (gate <= 1% on the
  finest level) and vs the body-fitted Abaqus reference 4.34797e-6 (reported).
- self-convergence of the max-principal-stress profile along the shortest span
  hole edge -> free edge (x = 0, y = 10 -> 20, z = 0) toward the FINEST level:
  fitted log-log rate (gate >= 1.8; Coreform p=2: 2.23).
Writes <EVID>/results.json.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

STEP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data",
                    "roark_plate_hole.step")
EVID = os.environ.get(
    "WHT_V02_EVID",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs",
                 "phase1", "evidence", "v02"))
REF_PEAK = 4.3175e-6
REF_ABAQUS = 4.34796947956784e-06  # Coreform roark-17p1-7a-accepted-peak.csv (body-fitted C3D20)
PROFILE_P0, PROFILE_P1 = (0.0, 10.0, 0.0), (0.0, 20.0, 0.0)
LO, HI = (-100.0, -20.0, -0.5), (100.0, 20.0, 0.5)
# Level 5 thickness split 8 -> 2 (user decision 2026-11-15): nu = 0 and an
# in-plane load make the exact solution z-independent, so nz does not change
# the in-plane resolution (h = 0.75) or the reference; (267, 53, 8) needed
# ~11.5 GB (LDLT) / 22.8 GB (LU) of factor and failed inside Kratos.
LEGACY_LEVELS = [(20, 4, 1), (40, 8, 2), (67, 13, 2), (80, 16, 4),
                 (133, 27, 4), (267, 53, 2)]
LEGACY_INPLANE = [10.0, 5.0, 3.0, 2.5, 1.5, 0.75]
CF_LEVELS = [(33 * 2**k, 7 * 2**k, 1) for k in range(5)]
CF_INPLANE = [6.0 / 2**k for k in range(5)]  # nominal Coreform h (actual 6.06 x 5.71 / 2^k)
if os.environ.get("WHT_V02_SET", "coreform") == "legacy":
    LEVELS, INPLANE = LEGACY_LEVELS, LEGACY_INPLANE
else:
    LEVELS, INPLANE = CF_LEVELS, CF_INPLANE


def run_level(nc: tuple[int, int, int], shape, outdir: str, h_ip: float) -> dict:
    import os as _os
    from wht_iga.geometry import cad as G
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver import trimmed as T
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.validate import JobSpec

    _p = int(_os.environ.get("WHT_V02_P", "2"))
    import os as _os
    _ext = _os.environ.get("WHT_V02_EXT", "")
    _theta = float(_os.environ.get("WHT_V02_THETA", "0.2"))
    _ghost = _os.environ.get("WHT_V02_GHOST_ALPHA", "")
    _out = os.path.join(EVID + ("_ext" if _ext else ""), "")
    box = BoxSpline(lo=LO, hi=HI, n_cells=nc, degree=_p)
    from wht_iga.progress import log

    h = float(min((np.asarray(HI) - np.asarray(LO)) / np.asarray(nc)))
    # Hole tessellation sag = in-plane h / 40 at every level (2026-11-15). The
    # old min(hx, hy, hz)/4 tied it to the thickness split: ratios ranged
    # h/6..h/40, and the level-5 nz 8->2 change coarsened the hole 4x
    # (deflection 0.031 -> 0.125, ~20-gon), polluting the finest peak.
    defl = float(h_ip / 40.0)
    log(f"V02 level n_cells={nc} h={h_ip:g} p={_p} ext={'theta=' + str(_theta) if _ext else 'off'}: "
        f"tessellate (deflection {defl:g})")
    _geom = _os.environ.get("WHT_V02_GEOM", "")
    if _geom:  # G2 gate: geometry from a mesh FILE (outer faces of a volume mesh)
        from types import SimpleNamespace

        from wht_iga.meshio_ext.model import outer_faces
        from wht_iga.meshio_ext.read_msh import read_msh

        _mm, _mi = read_msh(_geom, "SI")
        _mn, _mt, _notes = outer_faces(_mm)
        tess = SimpleNamespace(nodes=_mn, tris=_mt, face_ids=None)
        log(f"V02 level h={h_ip:g}: geometry from mesh file {os.path.basename(_geom)} "
            f"({len(_mt)} outer triangles; notes={_notes})")
    else:
        tess = G.tessellate(shape, deflection=defl)
    log(f"V02 level h={h_ip:g}: {len(tess.tris)} surface triangles; trimmed table")
    qp_table, out_cells, issues, cell_records, diag = T.build_trimmed_table(
        tess.nodes, tess.tris, box, q=2 * _p)
    fid = (None if tess.face_ids is None else
           T.hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids, radius=10.0))
    from wht_iga.solver import bc_roark as _bc
    _bc_kw = _bc.roark_left_pinned()
    if _ghost:
        _bc_kw["ghost_alpha"] = float(_ghost)
    # Memory knobs for h=0.375 on a ~15 GB-free machine (RESPONSE_ANALYSIS_G2_
    # BOTTLENECKS §4): WHT_V02_OUT_MODE=drop skips the hole-interior OUT cells,
    # WHT_V02_THREADS caps MKL threads (Pardiso per-thread buffers). Defaults
    # keep the G1/G2 evidence configuration (ghost, MKL default threads).
    _thr = _os.environ.get("WHT_V02_THREADS", "")
    spec = JobSpec(
        name=f"v02_h{h_ip:g}", lo=LO, hi=HI, n_cells=nc, degree=_p, E=1.0, nu=0.0,
        out_cells=out_cells, output_dir=_out,
        linear_solver=_os.environ.get("WHT_V02_SOLVER", "auto"),
        out_mode=_os.environ.get("WHT_V02_OUT_MODE", "ghost"),
        num_threads_mkl=int(_thr) if _thr else None,
        **_bc_kw)  # type: ignore[arg-type]
    os.makedirs(_out, exist_ok=True)
    res = run_linear(spec, geo_issues=issues, qp_table=qp_table,
                     cell_records=cell_records,
                     extension={"theta": _theta} if _ext else None)
    from wht_iga.solver.probes import peak_max_principal_fast

    vol_peak = peak_max_principal_fast(res.cauchy_stress)  # batched, == T.peak_max_principal
    if fid is None:  # mesh-file geometry: no CAD face ids -> facet peak not defined
        bnd = {"peak": float("nan"), "n_points": 0}
    else:
        bnd = T.boundary_peak_stress(box, res.U, 1.0, 0.0, tess.nodes, tess.tris,
                                     tess.face_ids, face_id=fid)
    peak = bnd["peak"]
    from wht_iga.solver import probes as PR

    Pp, _wp = PR.line_gauss_points(PROFILE_P0, PROFILE_P1)
    profile = PR.max_principal_profile(box, res.U, 1.0, 0.0, Pp)
    # Coreform convention: peak evaluated AT point A (hole edge, on the exact
    # circle r = 10, x = 0, mid-thickness), not max over tessellated facets
    # (inscribed-polygon points lie inside the hole -> extrapolated stress).
    peak_A = float(PR.max_principal_profile(
        box, res.U, 1.0, 0.0, np.array([[0.0, 10.0, 0.0], [0.0, -10.0, 0.0]])).max())
    np.savez_compressed(os.path.join(_out, f"v02_h{h_ip:g}_U.npz"), U=res.U,
                        lo=np.asarray(LO), hi=np.asarray(HI), n_cells=np.asarray(nc),
                        degree=_p)
    umax = float(np.max(np.abs(res.U)))
    err = abs(peak - REF_PEAK) / REF_PEAK
    # crown-zone fallback cells: fallback CUT cells whose center is within
    # 2r of the hole axis — these directly pollute the peak probe.
    hx = (HI[0] - LO[0]) / nc[0]
    hy = (HI[1] - LO[1]) / nc[1]
    crown_fb = 0
    for g in issues:
        if g.get("code") == "QUAD_FALLBACK":
            try:
                ci, cj, _ck = (int(v) for v in str(g.get("entity", "")).split(":")[1].split(","))
                cx, cy = LO[0] + (ci + 0.5) * hx, LO[1] + (cj + 0.5) * hy
                if cx * cx + cy * cy < (2.0 * 10.0) ** 2:
                    crown_fb += 1
            except (ValueError, IndexError):
                pass
    row = {"n_cells": list(nc), "h": h_ip, "peak": peak, "vol_peak": vol_peak,
           "bnd_pts": bnd["n_points"], "ref": REF_PEAK,
           "rel_err": err, "rel_err_abaqus": abs(peak - REF_ABAQUS) / REF_ABAQUS,
           "peak_A": peak_A, "rel_err_A": abs(peak_A - REF_PEAK) / REF_PEAK,
           "rel_err_A_abaqus": abs(peak_A - REF_ABAQUS) / REF_ABAQUS,
           "profile": profile.tolist(), "n_relaxed": diag.get("n_relaxed", 0),
           "timings": res.timings, "beta": res.beta_used, "n_cut": diag["n_cut"],
           "n_fallback": diag["n_fallback"], "crown_fallback": crown_fb,
           "umax": umax,
           "tie_info": res.tie_info,
           "total_force": [float(v) for v in np.atleast_1d(res.total_force)]}
    print(f"[V02 h={h_ip:g}] peak={peak:.6e} (vol {vol_peak:.3e}) ref={REF_PEAK:.6e} "
          f"err={err:.3%} errAbq={row['rel_err_abaqus']:.3%} peakA={peak_A:.6e} "
          f"errA={row['rel_err_A']:.3%} errA_Abq={row['rel_err_A_abaqus']:.3%} cut={diag['n_cut']} fb={diag['n_fallback']} "
          f"crownfb={crown_fb} umax={umax:.3e} out={len(out_cells)} ties={res.tie_info}", flush=True)
    favg = {str(k): [float(vv) for vv in np.atleast_1d(v)] for k, v in res.face_averages.items()}
    print(f"  total_force={row['total_force']} face_avgs={favg}", flush=True)
    return row


def main() -> None:
    from wht_iga.geometry import cad as G
    from wht_iga.solver import probes as PR

    import os as _os
    evd = EVID + ("_ext" if _os.environ.get("WHT_V02_EXT", "") else "")
    os.makedirs(evd, exist_ok=True)
    shape = G.read_step(STEP)
    _lv = _os.environ.get("WHT_V02_LEVELS", "")
    idx = list(range(len(LEVELS))) if not _lv else [int(v) for v in _lv.split(",")]
    rows = [run_level(LEVELS[i], shape, evd, INPLANE[i]) for i in idx]
    hs = np.array([r["h"] for r in rows])
    es = np.array([r["rel_err"] for r in rows])  # facet-max peak (info)
    esA = np.array([r["rel_err_A"] for r in rows])  # point-A peak (gate, Coreform)
    _P, wp = PR.line_gauss_points(PROFILE_P0, PROFILE_P1)
    sc = PR.self_convergence(hs, [np.asarray(r["profile"]) for r in rows], wp)
    rate_peak_roark = (float(np.polyfit(np.log(hs), np.log(es), 1)[0])
                       if len(es) >= 2 and np.all(es > 0) else float("nan"))
    out = {"level_set": _os.environ.get("WHT_V02_SET", "coreform"),
           "rows": rows, "profile_self_convergence": sc,
           "rate_peak_vs_roark_info_only": rate_peak_roark,
           "gate": {"peak_A_vs_roark_finest_max": 0.01, "profile_rate_min": 1.8,
                    "peak_definition": "max principal at point A = (0, +-10, 0) on the exact hole circle",
                    "convention": "Coreform roark-17p1-7a §7.5 (self-convergence of the "
                                  "transverse max-principal profile toward the finest level)"}}
    with open(os.path.join(evd, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"[V02] finest point-A peak err vs Roark={esA[-1]:.3%} (gate 1%), vs Abaqus="
          f"{rows[-1]['rel_err_A_abaqus']:.3%} (facet-max {es[-1]:.3%}); profile self-conv errors="
          f"{['%.2e' % e for e in sc['errors']]} rate_all={sc['rate_all']:.2f} "
          f"rate_last3={sc['rate_last3']:.2f} (gate 1.8)", flush=True)
    ok = esA[-1] <= 0.01 and sc["rate_all"] >= 1.8
    print("V02 ALL PASS" if ok else "V02 FAIL")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
