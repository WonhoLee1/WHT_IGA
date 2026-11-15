# ADR-002: penalty default beta (c_pen)

- Status: decided (Track SOLVE, Phase 1)
- Date: 2026-11-14 / branch: WHTIGA / Kratos 10.4.3 source build (`E:\kratos-build\install`)
- Plan ref: `WHT_IGA/research-cf-IGA/master-plan_v2.md` section 5.4, 7.6
  (beta = c_pen . E . p^2 / h, default c_pen = 1e3 placeholder, fixed by sweep)

## Decision

Default `c_pen = 3.0e5`, i.e. `beta = 3e5 . E . p^2 / h`, implemented as
`wht_iga.solver.penalty.DEFAULT_C_PEN` (single source; `default_beta()`,
`run_linear()` use it when the spec leaves `beta = None`).
An explicit user `beta` (Coreform `bc_penalty_value` equivalent) always wins.

## Evidence: beta sweep on the patch cube

Setup (all runs `wht_iga.solver.pipeline.run_linear`, 2x2x2 unit cube,
E = 1, nu = 0, all-face BC-A affine field u = (x, 2y, 3z), SkylineLU):
exact field representable in p >= 1, so the sweep isolates penalty error.
Reference: analytic field; cross-checked against the reference assembler +
consistent penalty (agreement 1.846e-6 vs 1.862e-6 at beta = 2e6, p = 2).

| p | mult (beta = mult.E/h, h = 0.5) | max |u-uex| | scaling |
|---|---|---|---|
| 1 | 1e2 / 1e3 / 1e4 / 1e5 / 1e6 | 2.467e-2 / 2.526e-3 / 2.532e-4 / 2.532e-5 / 2.532e-6 | 1/beta, C ~= 5.1 |
| 2 | 1e2 / 1e3 / 1e4 / 1e5 / 1e6 | 1.827e-2 / 1.858e-3 / 1.861e-4 / 1.862e-5 / 1.862e-6 | 1/beta, C ~= 3.7 |
| 3 | 1e2 ... 1e6 / 1e7 / 1e8 | ... / 4.070e-6 / 4.070e-7 / 4.070e-8 | 1/beta, C ~= 8.1, no floor |

(log-log slope -1.00 +- 0.01 all p; `test_p1_penalty_sweep.py` asserts
slope -1 +- 0.1 for p = 2. No plateau in 1e2..1e6: the patch field is
exactly representable, so penalty error falls to the conditioning floor;
mult = 1e8 still exact 1/beta, direct solver unaffected.)

V01 check with the default (max-relative error, max|u| = 3):

| p | default beta (c = 3e5, h = 0.5) | V01 rel err (measured) | margin to 1e-5 |
|---|---|---|---|
| 1 | 6.0e5 | 2.814e-6 | 3.5x |
| 2 | 2.4e6 | 5.172e-7 | 19x |
| 3 | 5.4e6 | 5.025e-7 | 20x |

p = 1 is the binding case and sets c_pen: worst measured constant
C ~= 8.2 (p = 3) needs c >= 1.5e4 for 1e-5; p = 1 needs c >= 8.4e4;
c = 3e5 keeps >= 3x margin everywhere on this grid.

## Findings recorded alongside

1. Per-face-cell triangulation is REQUIRED (not whole-face triangles):
   B-spline products are only piecewise smooth, so a Dunavant rule over a
   triangle spanning knot lines underintegrates (measured 28% penalty-matrix
   error on 2x2 faces; Kratos P vs reference K+P rel 0.23). Fix:
   `box_face_triangles(..., subdiv=face_subdiv(box, axis))` in `penalty.py`;
   same fix applied to `loads.py`. Single-cell faces unaffected (S09
   1-cell replication through the same builder: avg 0.099999999997,
   maxdev 3.4e-12).
2. p = 3 uses Dunavant order 6 (S09 rank rule: per-face-cell points ≥
   (p+1)²; 12 < 16 at order 4). Re-verified 2026-11-14 on the ADR patch sweep
   (2x2x2 cube, p=3, β = 1e2..1e8): order 4 AND order 6 both scale exact 1/β
   (slopes −0.997; 4.07e-8 vs 3.74e-8 at the bottom) — the 1-cell order-4
   floor (1.8e-4) does not appear on the patch because neighboring face
   cells constrain the shared control points. The order-6 default is kept
   as the guarantee (single-cell safety), not because the sweep demanded it.
3. `SupportPenaltyCondition` penalizes all 3 components (K6) -- per-component
   needs use BC-B (virtual node + LinearMasterSlaveConstraint +
   NodalConcentratedElement3D1N + PointLoadCondition3D1N), verified in
   `test_v01_bc_b_symmetry_plane` (<= 1e-5).
4. Nitsche stays out (S18: 3D unsupported, segfault without Check) --
   penalty is the Phase 1 path; Nitsche remains a C++ item.

## Consequences

- `run_linear` default beta meets G1 V01 (<= 1e-5) on the 2x2x2 grid.
- Finer grids keep the margin: C ~ u.E/h scaling cancels the 1/h in beta.
- Very stiff materials (E ~ 2e5) give beta ~ 1e11 at h ~ O(1): direct
  solver verified to beta = 2e8 on E = 1; conditioning for large E is
  scale-invariant (K and P scale together), no action.
- Re-tune if V02/V03 trimmed cases show a penalty plateau above 1e-5.

## Reproduce

`... python -m pytest WHT_IGA/tests/kratos_spikes/test_p1_penalty_sweep.py
WHT_IGA/tests/kratos_spikes/test_p1_pipeline_cube.py`
