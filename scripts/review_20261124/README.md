# Review 2026-11-24 — diagnostics

Scripts were run from the session scratchpad; output paths inside them point there. Edit the paths before re-running. Python: `D:\mamba\envs\whtiga\python.exe`, cwd `WHT_IGA`.

| File | Case | What it shows |
|---|---|---|
| `e03_diag_elastic.py` / `.json` | E03 | Metals elastic (yield and hardening removed): 45/45, no cutback |
| `e03_diag_t20.py` / `.json` | E03 | 304L −55 °C row replaced by the 20 °C row: heat-start failures gone, one failure at level 1.62 |
| `e03_echo3.py`, `e03_echo3_failed_increment.log` | E03 | Kratos Newton residual history of a failed increment: stall at 4e-7–1e-6 |
| `e03_tol.py`, `e03_tol1e-6.json`, `e03_tol1e-8.json` | E03 | Coreform history at `tol_rel=1e-6` (45/45, no cutback) and `1e-8` (51 solves, 3 cutbacks): Newton iterations per solve, full step-end snapshots; step-end values agree to ≤ 3.1e-8 relative |
| `v06_energy_proj.py` / `.json` | V06 | Kratos-free Galerkin / energy projection / L2 best approximation on our trimmed space (L0 p1, L1 p1, L0 p2) |
| `v06_facet_load.py` / `.json` | V06 | Effect of facetted inner-sphere pressure (n_geo) on the p=2 L0 error |
| `v06_p2_variants.py` / `.json` | V06 | Current-code Kratos L0 p=2: extension on/off, penalty/strong symmetry |
| `e04_p3c1.py` / `.json` | E04 | Cook membrane with U(3, C¹)/P0: n=8, 16 |
