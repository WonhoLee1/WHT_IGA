# V01-ext — prescribed-displacement cube with extension ON

- Status: **PASS (p=1,2,3) after fix** — extension operator defect found and fixed here (2026-11-14, 통합 담당)
- Gate (G2, v2.1 §9 "V01, V02를 extension 켠 상태로 재통과"): V01 p ∈ {1,2,3} with
  extension enabled, max relative error ≤ 1e-5 (reference: analytic u = (x, 2y, 3z), E = 1, nu = 0)
- Script: `scripts/v01_ext.py` → `docs/phase2/evidence/v01_ext/results.json` (θ = 0.2)

## Cases

- A "noop": untrimmed 2×2×2 cube, all η = 1 → extension must build 0 ties and equal the OFF solve.
- B "synthetic": 4×4×4 cube, corner column cells (i=3, j=3, all k) given η = 0.02 so the
  support classification creates slaves. Integration is still the full cube, so the affine
  exact solution must remain representable. n_slave > 0 is required (no by-construction pass).

## Result before fix (run 2026-11-14)

| p | A err OFF | A err ON | A max\|ΔU\| | A ties | B err | B ties / slaves | 판정 |
|---|---|---|---|---|---|---|---|
| 1 | 2.81e-6 | 2.81e-6 | 2.2e-15 | 0 | 1.12e-6 | 15 / 5 | PASS |
| 2 | 5.17e-7 | 5.17e-7 | 2.7e-15 | 0 | **1.34e-2** | 18 / 6 | FAIL |
| 3 | 4.61e-7 | 4.61e-7 | 7.6e-14 | 0 | **4.24e-3** | 21 / 7 | FAIL |

## Root cause (verified in Python, no Kratos)

`wht_iga/stabilize/extension.py` extrapolates with 1D Lagrange weights in **uniform index
space** (`lagrange_weights_1d(size, s - lo)`). On clamped (open) knot vectors with p ≥ 2 the
Greville abscissae near the grid boundary are non-uniform (p=2, n=4: 0, .125, .375, .625,
.875, 1). Affine fields have coefficients = affine function of Greville points, so index-space
extrapolation does not reproduce them near the boundary.

Direct check of `op.E` on exact coefficients (same 4×4×4 setup, scratch script):

| p | x | y | z | 1 |
|---|---|---|---|---|
| 1 | 0 | 0 | 0 | 0 |
| 2 | 0 | **0.125** | 0 | 0 |
| 3 | 0 | **0.0833** | 1e-16 | 0 |

(max slave-coefficient error; constants pass = partition of unity holds.)

The Phase-2 unit tests (`tests/unit/test_stabilize_extension.py::test_degree_reproduction_1e10`,
`test_operators_compose.py`) build the manufactured polynomial **in index space**, so they pass
by construction and cannot catch this.

## Impact

- Slaves whose master block lies within p cells of a clamped grid boundary: wrong (not
  polynomial-reproducing) extension → consistency error. Interior uniform-knot regions are exact.
- V02 ext: hole is in the grid interior, outer faces conforming → likely small effect, but must be
  re-checked after a fix. V06 / any geometry cut near the bounding box is affected.

## Fix (applied 2026-11-14, user decision: fix in integration session)

- `extension.py`: 1D weights now solve `sum_m w_m c_k(lo+m) = c_k(s)`, k < block size, with
  c_k = B-spline coefficients of monomials on the TRUE knots (Marsden blossom
  `e_k(t_{i+1..i+p}) / C(p,k)`, `extrapolation_weights_1d`). `build_extension(..., knots=)` is
  now a required keyword (no silent index-space default). Uniform interior = old Lagrange
  (`test_interior_matches_index_lagrange`).
- Tests: physical-coordinate reproduction (Greville-collocation oracle) p=1..3,
  `test_clamped_boundary_affine_regression` (the failing config above), knots mismatch guard;
  compose unit + Kratos tests switched to physical-coordinate manufactured fields.
- Scratch coefficient check after fix: x, y, z, 1 errors ≤ 2.2e-16 (p=1..3), x^p·y² ≤ 3.3e-16.

## Result after fix (run 2026-11-14)

| p | A err ON | A max\|ΔU\| | A ties | B err | B ties / slaves | 판정 |
|---|---|---|---|---|---|---|
| 1 | 2.81e-6 | 2.2e-15 | 0 | 1.12e-6 | 15 / 5 | PASS |
| 2 | 5.17e-7 | 8.9e-15 | 0 | 3.10e-7 | 18 / 6 | PASS |
| 3 | 4.61e-7 | 1.1e-13 | 0 | 1.80e-7 | 21 / 7 | PASS |

Evidence: `docs/phase2/evidence/v01_ext/results.json` (all_pass = true).
