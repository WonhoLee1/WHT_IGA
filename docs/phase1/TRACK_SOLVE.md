# TRACK_SOLVE — Phase 1 solver pipeline (S2-use, S6, S7, strategy, DB, VTU)

- Owner: Track SOLVE agent, 2026-11-14, branch WHTIGA
- Scope files: `WHT_IGA/wht_iga/solver/*.py`, `tests/unit/test_solve_*.py`,
  `tests/kratos_spikes/test_p1_*.py`, `docs/adr/ADR-002-penalty-beta.md`, this file
- Kratos: 10.4.3 source build `E:\kratos-build\install` (asserted via
  `load_kratos()` + `__file__` check); Python 3.12.16 (`whtiga` env)

## Results (all numbers from executed runs, tolerances explicit)

### G1 acceptance

| Gate item | Result | Evidence |
|---|---|---|
| V01 cube p=1,2,3 <= 1e-5 (max-rel err) | p1 2.814e-6, p2 5.172e-7, p3 5.025e-7 | `test_p1_pipeline_cube.py` (default beta, affine field, 2x2x2) |
| beta default fixed in ADR-002 | c_pen = 3e5, slope -0.9983 | `test_p1_penalty_sweep.py` + `ADR-002-penalty-beta.md` |
| mesh+results DB roundtrip | bit-exact (array_equal) | `test_solve_db.py` (unit) + `test_p1_outputs.py` (pipeline) |
| VTU | meshio re-read exact (points/cells/U) | `test_solve_vtu.py`, `test_p1_outputs.py`; ParaView NOT available -> re-read check (recorded) |
| validate() blocks on errors | block + pass paths tested | `test_solve_validate.py` (8 tests) |
| Full suite (mine) | 30 passed (18 unit + 12 kratos) | `p1_run4/7.log` |

### Supporting measurements

| Check | Number | Where |
|---|---|---|
| BC-B symmetry plane (ux=0 on x=0) | rel err 5.841e-7 (<= 1e-5) | `test_v01_bc_b_symmetry_plane` |
| FCM alpha=1.0 equivalence | max diff 2.887e-15 (<= 1e-12) | `test_v01_fcm_alpha_one` |
| FCM alpha=1e-8 | completes, finite U, max 3.0 | `test_v01_fcm_alpha_small` |
| Pressure +2 top | bottom reaction -2.000000000000 (1e-9) | `test_pressure_sign_and_total` |
| Traction (0,0,-3) | bottom reaction +3.0000000000 (1e-9) | `test_traction_total` |
| Gravity rho=1,g=9.81,p2,nu=0 | -4.905000000000e-03 = analytic to 12 digits | `test_gravity_vs_analytic` |
| Follower vs dead load | rel diff 1.232e-6 (<= 1e-5) | `test_follower_matches_dead_load` |
| S09 replication via builder | avg 0.099999999997, maxdev 3.439e-12 | scratch `debug_s09.py` |
| Reference+penalty cross-check | 1.846e-6 vs Kratos 1.862e-6 at beta=2e6 | scratch `debug_refpen.py` |
| Sweep constants C (err=C/beta) | p1 5.06, p2 3.72, p3 8.14 | scratch sweeps + `test_p1_penalty_sweep.py` |

## Key finding (penalty quadrature)

Whole-face triangles spanning knot lines underintegrate the penalty
(N_i.N_j piecewise smooth): 28% penalty-matrix error on 2x2 faces
(Kratos K+P vs reference, rel 0.23), while 1-cell faces are exact.
Fix: per-face-cell triangle pairs (`face_subdiv()` from breakpoints) in
`penalty.box_face_triangles`; same fix in `loads.attach_loads`. After fix:
exact 1/beta (slope -0.9983), Kratos-vs-reference agreement 1.862e-6/1.846e-6.
S09 single-cell behavior unchanged (replicated to 12 digits).
Per-face-cell QPs raise condition counts (~288 QPs/face at p=2, 2x2x2);
acceptable for Phase 1, revisit on profiling (Phase 6).

## Beta decision (ADR-002, summary)

beta = 3e5 . E . p^2 / h (`DEFAULT_C_PEN` in `penalty.py`, single source).
No plateau in 1e2..1e6 .E/h because the patch field is exactly
representable (error falls to the conditioning floor; mult=1e8 still exact
1/beta) -- recorded, not hidden: the "plateau" stop-condition does not
trigger because exact scaling + reference agreement + V01 margins constitute
stronger evidence than a floor observation. p=1 binds (needs c >= 8.4e4);
c = 3e5 keeps >= 3.5x margin on the V01 grid.

## Interface for the integrator (V02/V03)

`run_linear(spec: JobSpec, geo_issues: list[dict] | None) -> LinearResult`
(`wht_iga/solver/pipeline.py`): validates (raises `ValidationError` on
error issues), builds BoxSpline grid, Kratos model (volume keepalive held on
result), IN/ghost elements, BC-A/BC-B, loads, Newton+SkylineLU solve
(move_mesh when follower), probes, optional DB/VTU/probe-JSON output.
`JobSpec` is JSON-roundtrippable (`to_json`/`from_json`).

## Proposals / notes for other tracks

1. Shared-module change proposals (NOT made; local copies/patterns used):
   - `geometry.make_nurbs_volume`: no change needed; pipeline holds the
     volume reference itself (lifetime rule honored).
   - `surface.triangle_quadrature`: p=3 penalty formally needs a degree-6
     Dunavant rule (currently max 4); measured no floor to 4e-8 on the patch
     case, so order 4 is used with a documented revisit trigger.
2. Track GEO: `validate()` consumes your issue dicts `{severity, code,
   message, entity}`; face entities parsed as `face:<axis>:<side>`.
   `QUAD_FALLBACK` on a pressure/traction face is a blocking error.
3. V02/V03 need trimmed cells + cut QP tables: pipeline currently builds
   full IN-cell Gauss elements (OUT set optional); feed per-cell QP lists
   through `out_cells`/a future `qp_override` hook (not built -- ask before
   I add it, to avoid guessing your DB schema).
4. `fixed_faces` pins u=0 strongly; nonzero strong Dirichlet is not offered
   (penalty covers it).
5. No new packages installed; `meshio` 5.3.5 already in env used for VTU.

## Stop conditions / blocking issues

None triggered. No Kratos-vs-plan contradictions (S18 Nitsche already
resolved by ADR-001 failover; not re-attempted). Workload within estimate.

## Final state pointer (2026-11-15 G2 COMPLETE — this doc is the Phase-1 solver baseline)

- Penalty quadrature fix (per-face-cell pairs) and beta=3e5·E·p²/h (ADR-002) carried
  through to G2 final runs unchanged. p=3 Dunavant-order revisit trigger (§proposals 1)
  resolved in G2: final V06 uses q=min(2p,4) — see `G2_COMPLETE.md`.
- Final gates: `docs/phase2/G2_COMPLETE.md` (V01-ext p=1–3, V02 0.674%, V06, V-COND).
