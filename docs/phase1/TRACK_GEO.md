# TRACK_GEO — S1+S3+S4 geometry pipeline (Phase 1)

- Scope: `wht_iga/geometry/*` (S1 CAD, S3 classify, S4 cut quadrature, quad DB).
- Branch: `WHTIGA` (verified `git status --short --branch` before work).
- Env: `micromamba run -n whtiga` (Python 3.12.16, OCP/numpy/scipy/numba preinstalled;
  nothing installed). No Kratos import in any module or test. No aieng/wht_topopt.
- Ownership respected: only geometry package + `tests/unit/test_geo_*` + `data/roark_*`
  + this file touched. No commits.

## 1. Results table (measured, not claimed)

| Stage | Case | Result |
|---|---|---|
| S1 read/heal | `roark_plate_hole.step` (gen script) | bbox ±(100,20,0.5) exact; BRep vol 7685.840735 = analytic L·D·T−πr²T (rtol 1e-12); BRep area 15914.513322 = analytic (rtol 1e-9) |
| S1 tessellate | deflection 0.125, 92 nodes / 184 tris / 7 faces | watertight (0 bnd edges), signed vol +7687.011 (orientation flip fix), area err 1.4e-4, edge-midpoint max_dev 0.028 < 0.125, `validate()` green, zero issues |
| S1 STEP roundtrip | write→read | bbox + area identical (rtol 0 / 1e-12); generator rerun reproduces geometry |
| S1 fuse | two 10³ boxes (general fuse) | vol 2000 (rtol 1e-9), watertight, validate green |
| S3 box | unit box, non-aligned grid 4³ | IN/OUT/CUT exact vs analytic; 3-axis votes agree everywhere |
| S3 sphere | r=0.5 OCP sphere, 5³ grid | margin-0.05 IN/OUT exact; CUT set nonempty |
| S3 winding | 40 random pts vs analytic | exact (margin 0.05 from surface) |
| S4 IN cell | q=1,2 tensor rule | degree-≤q polys exact, rtol 1e-12 |
| S4 tet cut | corner tet s=0.9, q=2 **(G1 gate)** | 27 pts, all w≥0, no fallback, volume + all 27 Legendre moments vs analytic tet formula, max rel err < 1e-12; Σw check 1e-12 |
| S4 slab R-face | coplanar x1 path, white-box `_cell_moments` | matches analytic (mapped basis, only M000 nonzero) 1e-12 |
| S4 surface | tet: area + flux≈0 (1e-12); full Roark: Σw=tess area 1e-12, flux <1e-9 | green |
| S4 fallback | slab: exact 1e-9; tet d3→d4 error shrinks (<10%) | QUAD_FALLBACK logged, w≥0 |
| S4 Roark grid | 42×10×4, q=2: 632 rules, **0 fallback** | Σw = tess vol, rel **1.2e-16**; max residual 0.0; max npts 27 |
| DB | `quad_db` write→read roundtrip | meta/rows/ΣdV/blobs lossless; Roark DB 17064 rows, ΣdV exact |
| Tests | `pytest WHT_IGA/tests/unit/test_geo_*.py` (marker-free) | **16/16 pass** (`-m "not kratos"` clean) |

## 2. G1-moment-fitting evidence

`test_tet_cut_polynomial_accuracy`: exact tessellation (flat faces) corner-tet in
one cell; rule reproduces all 27 q=2 Legendre moments of the analytic
`a!b!c!s^{n+3}/(n+3)!` reference to < 1e-12 relative, volume to 1e-12, Σw check
1e-12, ≤27 nonnegative points, no fallback. Roark end-to-end: 616 fitted CUT/IN
cells all residuals ≤ 1e-12·|cell|.

## 3. Performance (additional requirement)

Wall times, Roark grid (1680 cells, 184 tris), 16 CPUs. Cold = first process run
(numba compile one-off ~11 s, cached across runs); warm = steady state.

| Hotspot | Cold | Warm (16t) | Warm (1t) |
|---|---|---|---|
| S1 bbox / mesh / validate | – | 3 / 20 / 133 ms | – |
| S3 broadphase (bbox binning, numpy) | – | 28 ms | 28 ms |
| S3 SAT narrow (`@njit parallel`, 8768 pairs) | 8.8 s (compile) | **5.9 ms** | 10.8 ms (1.8×) |
| S3 scanline parity axes x/y/z (Python lines + njit) | – | 10 / 34 / 76 ms | – |
| S3 winding fallback, 537 pts | – | 6 ms | – |
| S3 `classify_cells` total | 11.4 s (compile) | **0.15 s** | 0.17 s |
| S4 moments / genuine-CUT cell (numpy+Duffy) | – | mean 6.8 / max 9.1 ms | 8.8 / 14.8 ms |
| S4 cell total (candidates+fit) | – | mean 78 / max 122 ms | mean 148 / max 293 ms |
| NNLS fit N=64/250/352 | – | 0.45 / 1.3 / 28 ms | same (scipy) |
| S4 Roark total (632 rules) | – | 3.5–4.5 s | – |
| S1→S4→DB end-to-end | – | **≈5 s** | – |

- NNLS budget: hard cap 2048 candidates (deterministic stride; never hit —
  max pool ≈1200); residual gate 1e-12·|cell|; fallback triggers on Roark: **0**
  (dense-grid escalation resolved all 16 cells that fail with boundary-only
  candidates). Worst measured single fit 28 ms; worst cell total 122 ms.
- Caches: `legendre_monomials`, `_gauss1d`, `_duffy`, `_tensor_gauss` are
  `lru_cache`d (hits confirmed: 103/47/47/62); OCP shape meshed once per
  geometry; no re-mesh/re-query of any face (validate does read-only traversals).
- Cell loops: SAT, raycast, parity-scan, winding are `@njit` (`parallel=True`
  where measured beneficial: winding 4.5×, SAT 1.8× on this size); cell-state
  assignment is fully vectorized (no Python-per-cell); per-cell driver stays
  Python (orchestrates scipy/NumPy, cannot be njit-compiled).
- S16 projection (52.5 s / 1e5 elems = 0.5 ms/elem): S3 ≈ 0.09 ms/cell ✓;
  S1 ≈ 0.1 ms/cell ✓; S4 genuine-CUT ≈ 80 ms/cell applies only to surface
  cells (~4k at 1e5 scale → ~330 s single-core). **Threatens Phase-6 pace
  single-core.** Proposals (not silent): (P1) batch-prune NNLS (prune 8–16 per
  refit, ~5–10× fewer scipy solves — fit dominates cell cost, moments only
  6.8 ms); (P2) multiprocessing over cell blocks per v2 §3.3 (embarrassingly
  parallel → /16 cores ≈ 25 s); (P3) keep octree fallback as the cheap path
  for low-value cells (G-gate decision); (P4) C++ only if P1+P2 insufficient
  (no single Python hotspot >30% besides scipy/njit internals — v2 §2.2 bar
  not met, revisit Phase 6). Cold-compile 11 s is one-off (cache persists).
- Grid-alignment noise (measured): face-aligned grids flip a few touch-only
  cells between runs (IN 215 vs 224) — OCP tessellation coordinates differ in
  last-ulp across processes (nodes hash differs run to run), moving coplanar
  touch cases across the 0.5 winding threshold. Quadrature stays exact per
  run (moments follow that run's tessellation). Recommendation: production
  verification on non-aligned grids (offset or non-dividing h); DB `meta`
  already carries job/grid/q for traceability — add tessellation hash there.

## 4. Algorithm notes vs v2 (deliberate, tested deviations)

- Surface points use Duffy-mapped tensor Gauss-Legendre (exact ≥ 2p+1) instead
  of Dunavant tables — same exactness, verifiable by the Σw checks.
- `_cross_section_moments` orients segments by parent-triangle outward normal
  (exact, no on-surface winding query); rect edges use an x1±eps inside OR-test
  (tangent rulings count); shared on-plane edges deduplicated by endpoints.
- Touch-only cells (surface grazes closed box, misses shrunk box) resolve by
  center → IN/OUT, no moments (exact for degenerate touch, avoids 0-volume fits).
- Candidates escalate: (q+2)³ Gauss → depth-2/3 surface leaves → dense
  (2q+4)³ grid → fallback. The dense stage is load-bearing (sliver cells fail
  without it; verified by NNLS residual 6.5 at N=64 vs 0.0 at N=352).
- Winding queries always use the full closed soup (partial patches give
  half-space values — found via slab R-face zero-moment bug).

## 5. Bugs found and fixed (all covered by tests)

1. REVERSED OCP faces need index flip (signed volume 208 vs 7686).
2. `STEPControl_Reader.ReadStep` → `ReadFile` in cadquery-ocp 7.8.
3. `clip_poly_to_box` bound signs flipped (emptied every clip).
4. Missing `q` arg in `_design_matrix` call inside `fit_moment_rule`.
5. Duffy order (3q+3)//2 under-integrates (ξ-degree d+1) → (3q+4)//2.
6. Fan triangles used parent area instead of own (52× volume error on Roark).
7. Ray double-count through shared edges → parity dedup by hit clustering.
8. OCP API shapes: `TopoDS.Face(shape)` downcast, `OCP.collections.List_TopoDS_Shape`,
   `BOPAlgo_Builder` arguments-only general fuse (matches v2 §7.1), pole-degenerate
   triangles break watertight check (filtered).
9. Test-reference bugs of mine: standard-vs-mapped Legendre basis (2×), einsum
   label collision, symmetric-interval hand expectations.

## 6. Shared-module change proposals

None — all work localized to owned files. Notes for others: Bézier-pattern
caching belongs to `wht_iga/grid/` (SOLVE track, not touched); no `grid/`,
`reference/`, `kratos_backend/`, `solver/`, `conftest.py`, `pyproject.toml`,
`uv.lock`, `docs/adr/*` changes made.

## 7. Stop conditions / blocking issues

None triggered. OCP API differs from v2 §7.1 assumptions in names only
(§5.8 documented); no new package needed; workload fits estimate.

## 8. Final state pointer (2026-11-15 G2 COMPLETE — this doc is the Phase-1 baseline record)

- GEO pipeline (S1/S3/S4/quad DB) is the quadrature foundation that G2 final runs
  built on. Fallback reduction (16–48 → 0–1) and NNLS candidate pruning are recorded in
  `docs/phase2/G2_report.md` §2 (not re-edited here).
- Final gates: `docs/phase2/G2_COMPLETE.md` (V02 0.674%, V06 L0–L2 × p=1–3, V-COND).
