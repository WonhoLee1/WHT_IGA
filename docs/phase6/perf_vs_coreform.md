# Phase 6 G-8 — Performance Profiling, Iterative Solver Spike & Coreform Comparison

- **Date / Environment**: 2026-11-19, `WHT_IGA` (`whtiga` Python 3.12.16, Kratos 10.4.3 at `E:\kratos-build\install`), 16 logical CPUs, 31.60 GiB total RAM (14.93 GiB available during run).
- **Scripts & Tests**:
  - Benchmark harness: [`scripts/perf_vs_coreform.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/scripts/perf_vs_coreform.py)
  - Integration test suite: [`tests/kratos_spikes/test_p6_perf_vs_coreform.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_perf_vs_coreform.py) (`5 passed in 43.37s`)
- **Peak RSS Guard**: All scaling points (`1,029`–`20,904` DOFs) ran with `gc.collect()` cleanup between stages; overall measured peak RSS was **`906.9 MB` (`0.907 GB`)**, strictly `< 1.0 GB`.

---

## 1. Coreform Reference Baseline (`glass-to-metal-seal`)

Source CSVs:
- [`glass-to-metal-seal-performance.csv`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal/glass-to-metal-seal-performance.csv)
- [`glass-to-metal-seal-problem-size.csv`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/user-guide/generated-downloads/example-problems/glass-to-metal-seal/glass-to-metal-seal-problem-size.csv)

| Metric / Phase | Resources / Unit | Measured Value |
|---|---|---:|
| CAE script | 1 process | `4.98 s` |
| IGA mesh | 16 MPI ranks | `318.81 s` (`5.31 min`) |
| Interop | 1 process | `3.60 s` |
| Abaqus solve (3 NLGEOM steps) | 1 host; 16 threads | `691.00 s` (`11.52 min`) |
| **Total Coreform wall time** | — | **`1,018.39 s` (`16.97 min`)** |
| Background / Active / Cut elements | count | `146,461` / `32,098` / `3,420` |
| Abaqus elements / nodes / variables / equations | count | `56,595` / `66,180` / `198,546` / `178,278` |
| Derived IGA mesh throughput (bg / active) | elem/s (16 ranks) | `459.40` / `100.68` |
| Derived Abaqus solve throughput | eq/s (16 threads) | `258.00` |
| Phase 0 S16 WHT_IGA reference (`p=2`, `46³`) | 1 process | `97,336` elems in `52.5 s` (`1,854.0 elem/s`), `+2.92 GB` RSS |
| **Phase 6 1M DOF Gate** | 16 cores | **`< 1,800 s` (`30 min`), `< 64 GB`** |

---

## 2. Per-Stage Wall Time & Memory Scaling (`p=1`, `p=2`, 3D B-Spline)

### 2.1 Untrimmed 3D Cube (`make_nurbs_volume`, `attach_cell_elements` S16 Rate, Kratos Assembly, `PardisoLU`)

| `p` | Grid | Elems | DOFs | `nnz` | `make_vol` (s) | `attach_elem` (s) | Rate (elem/s) | `BuildRHS` (s) | `Build` (s) | `to_csr` (s) | `Solve net` (s) | Peak RSS (MB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `6×6×6` | 216 | 1,029 | 61,731 | 0.0041 | 0.0621 | 3,480.1 | 0.0135 | 0.0293 | 0.0134 | 0.1347 | 164.4 |
| 1 | `10×10×10` | 1,000 | 3,993 | 268,119 | 0.0096 | 0.2510 | 3,984.7 | 0.0340 | 0.0717 | 0.0601 | 0.0835 | 202.2 |
| 1 | `14×14×14` | 2,744 | 10,125 | 715,563 | 0.0244 | 1.1894 | 2,307.1 | 0.0401 | 0.0830 | 0.1425 | 0.1625 | 292.9 |
| 1 | `16×16×16` | 4,096 | 14,739 | 1,058,841 | 0.0447 | 1.3743 | 2,980.4 | 0.0455 | 0.0902 | 0.2194 | 0.3078 | 378.8 |
| 2 | `6×6×6` | 216 | 1,536 | 353,736 | 0.0049 | 0.1360 | 1,588.2 | 0.0693 | 0.6892 | 0.1318 | 0.0874 | 216.6 |
| 2 | `9×9×9` | 729 | 3,993 | 1,058,841 | 0.0092 | 0.4711 | 1,547.6 | 0.1714 | 1.2740 | 0.2118 | 1.1766 | 323.0 |
| 2 | `12×12×12` | 1,728 | 8,232 | 2,359,296 | 0.0768 | 1.1050 | 1,563.9 | 0.2677 | 1.4112 | 0.3764 | 0.6987 | 533.2 |
| 2 | `13×13×13` | 2,197 | 10,125 | 2,956,581 | 0.0494 | 1.0514 | 2,089.7 | 0.2248 | 1.4832 | 0.4400 | 0.6752 | 643.8 |

Observations:
- `attach_cell_elements` element creation rate stabilizes at **`~2,980 elem/s`** (`p=1`) and **`~1,550–2,090 elem/s`** (`p=2`), matching the S16 `46³` (`97,336` elems) observation (`1,854 elem/s`).
- `p=2` at `10,125` DOFs has `2.96M` nonzeros (`292 nnz/DOF`), `4.13×` denser than `p=1` at the same `10,125` DOFs (`715k` nonzeros, `70.7 nnz/DOF`), driving `PardisoLU` RSS to `643.8 MB` vs `292.9 MB` (`p=1`).

### 2.2 Trimmed 3D Sphere Octant (`r_in=0.45, r_out=0.95` in `[0,1]³`, `out_mode="drop"`, `theta=0.2`)

Note: [`trimmed.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/trimmed.py#L21) runs serial below `_PARALLEL_MIN_CELLS = 400` active cells and switches to an 8-thread pool at `>= 400` active cells (`cache_dir="off"` throughout).

| `p` (`q`) | Grid | Active / Cut / OUT | Active DOFs (`nnz`) | Master / Ext / Ties | `trimmed_table` (s) | `extension` (s) | `attach` (s) | `ties` (s) | `Build` (s) | `Solve net` (s) | Peak RSS (MB) |
|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 (`2`) | `6×6×6` | 132 / 119 / 84 | 732 (`50,094`) | 210 / 34 / 102 | 2.367 (1T) | 0.0335 | 0.0197 | 0.0339 | 0.0474 | 0.0111 | 252.6 |
| 1 (`2`) | `10×10×10` | 535 / 368 / 465 | 2,427 (`197,829`) | 707 / 102 / 306 | 0.950 (8T) | 0.0726 | 0.1569 | 0.1454 | 0.1641 | 0.0099 | 297.2 |
| 1 (`2`) | `14×14×14` | 1,353 / 752 / 1,391 | 5,556 (`462,096`) | 1,661 / 191 / 573 | 1.689 (8T) | 0.2008 | 0.4002 | 0.1894 | 0.1693 | 0.0957 | 348.5 |
| 2 (`4`) | `6×6×6` | 132 / 119 / 84 | 1,194 (`333,630`) | 355 / 43 / 129 | 8.129 (1T) | 0.0343 | 0.1303 | 0.0418 | 1.3114 | 0.0000* | 331.8 |
| 2 (`4`) | `9×9×9` | 393 / 287 / 336 | 2,697 (`852,921`) | 821 / 78 / 234 | 17.121 (1T) | 0.1067 | 0.6025 | 0.1358 | 4.1705 | 0.0000* | 423.9 |
| 2 (`4`) | `12×12×12` | 878 / 542 / 850 | 5,124 (`1,730,448`) | 1,569 / 139 / 417 | 6.051 (8T) | 0.1161 | 1.2158 | 0.2553 | 3.8280 | 1.0853 | 607.4 |
| 2 (`4`) | `13×13×13` | 1,100 / 650 / 1,097 | 6,186 (`2,179,476`) | 1,863 / 199 / 597 | 6.133 (8T) | 0.1025 | 0.6234 | 0.2307 | 4.1655 | 0.2402 | 696.6 |

*\*At small `<= 2,697` DOFs, `Solve net = max(0, strategy_s - build_s)` floors at `0.0` when `Build` dominates (`> 95%`) within run-to-run CPU jitter; pure linear solve times are isolated in §3.*

Key findings:
- **8-thread pool speedup**: Crossing `_PARALLEL_MIN_CELLS=400` drops `p=2` `build_trimmed_table` from `17.12 s` (`287` cut cells, serial) to `6.05 s` (`542` cut cells, 8 threads) — **`~5.3×` higher cut-cell throughput** (`~9.4–11.2 ms/cut cell`) with negligible extra RSS (`+36–56 MB`).
- **MPC tie assembly overhead**: Adding `LinearMasterSlaveConstraint` (`T^T K T` inside `ApplyConstraints`) increases `kratos_build_s` roughly `2.7–2.8×` vs untrimmed (`4.17 s` for `1,100` trimmed elems vs `1.48 s` for `2,197` untrimmed elems at `p=2`, `13³`).

---

## 3. Linear Solver Comparison & Iterative Solver Spike (`AMGCL` / `AdditiveSchwarz`)

### 3.1 Kratos Strategy Solver Matrix (`p=2`)

Tested on fresh `ModelPart` instances per solver (`tol=1e-10`, `max_iter=2000`):
- **Untrimmed**: `9×9×9` (`3,993` DOFs, `nnz=1,058,841`, `assembly_ref_s=0.8691 s`)
- **Trimmed + Extension**: `8×8×8` (`2,121` active eqs, `nnz=703,035`, `assembly_ref_s=1.7460 s`)
- **Trimmed without Extension**: `8×8×8` (`2,121` active eqs, `nnz=485,631`, `assembly_ref_s=1.5312 s`)

| Solver | Untrimmed `9³` Strat (s) / Iters / Diff vs LU / RSS (MB) | Trimmed+Ext `8³` Strat (s) / Iters / Diff vs LU / RSS (MB) | Trimmed No-Ext `8³` Strat (s) / Iters / Diff vs LU |
|---|---|---|---|
| `PardisoLU` | `1.352s` / — / `0.0` / `492.4` | `1.638s` / — / `0.0` / `532.7` | `1.727s` / — / `0.0` |
| `PardisoLDLT` | `2.415s` / — / `2.11e-15` / `456.9` | `1.889s` / — / `1.77e-13` / `508.9` | `1.793s` / — / **`1.000` (FAIL)** |
| `SkylineLU` | `6.817s` / — / `8.01e-15` / `410.0` | `2.180s` / — / `2.49e-13` / `479.8` | `2.665s` / — / `4.39e-13` |
| `AMGCL_CG_ILU0` | `0.844s` / *0* / `4.14e-11` / `412.9` | `1.643s` / *0* / `5.63e-11` / `480.3` | `1.678s` / *0* / `3.57e-10` |
| `AMGCL_CG_SPAI0_Block` | `1.030s` / *0* / `4.23e-10` / `414.6` | `2.756s` / *0* / **`8.10e-03` (STALL)** / `480.6` | `2.289s` / *0* / **`1.005` (FAIL)** |
| `AMGCL_GMRES_ILU0` | `1.276s` / *0* / `1.55e-10` / `417.7` | `1.715s` / *0* / `5.09e-11` / `480.4` | `1.617s` / *0* / `3.00e-10` |
| `CG_Diag` | `0.864s` / `67` / `2.49e-10` / `418.0` | `1.715s` / `151` / `1.11e-08` / `480.4` | `1.650s` / `156` / `4.92e-08` |
| `CG_ILU0` | `1.110s` / `47` / `5.22e-10` / `430.9` | `4.080s` / `2000` / **`2.950` (DIV)** / `488.8` | `3.185s` / `2000` / **`212.6` (DIV)** |
| `CG_AdditiveSchwarz` | `0.968s` / `84` / `2.04e-10` / `437.1` | `2.719s` / `2000` / **`30.17` (DIV)** / `457.4` | `2.092s` / `2000` / **`4.00e14` (DIV)** |
| `Scaling_CG_Diag` | `0.754s` / `67` / `2.49e-10` / `389.2` | `1.813s` / `152` / `1.03e-08` / `450.3` | `1.609s` / `155` / `6.30e-08` |
| `Scaling_AMGCL_CG` | `1.018s` / *0* / `3.44e-10` / `397.8` | `1.894s` / *0* / `4.62e-11` / `453.6` | `1.535s` / *0* / `2.68e-10` |

*Note: Kratos `AMGCLSolver` does not override `LinearSolver::GetIterationsNumber()`, so Kratos reports `0` with a base-method warning.*

Also observed on `SkylineLU` at larger sizes (separate untrimmed probe):
- `p=1, 16³` (`14,739` DOFs): `SkylineLU` `77.50 s` vs `PardisoLU` `1.09 s` (**`71.1×` slower**).
- `p=2, 12³` (`8,232` DOFs): `SkylineLU` `51.24 s` vs `PardisoLU` `3.29 s` (**`15.6×` slower**).

### 3.2 Isolated CSR Linear Solve (`wht_iga.kratos_backend.solvers.solve_linear_csr`)

| Case | `pardiso_lu` (s / rel_res) | `pardiso_ldlt` (s / rel_res) | `amgcl` (s / rel_res) | `cg` SciPy ILU (s / it / rel_res) | `splu` (s / rel_res) |
|---|---:|---:|---:|---:|---:|
| Untrimmed `9³ p=2` (`3,993` eq) | `0.1113s` / `6.4e-16` | `0.0935s` / `6.2e-16` | `1.7318s` / `3.4e-11` | `0.7454s` / `6` / `7.8e-12` | `0.6502s` / `1.2e-15` |
| Trimmed+Ext `8³ p=2` (`2,121` eq) | `0.0731s` / `5.1e-15` | `0.0745s` / `4.9e-15` | `1.1583s` / `4.0e-11` | `3.9597s` / `2000` / `1.7e-05` (fail) | `0.1707s` / `1.1e-14` |
| Trimmed No-Ext `8³ p=2` (`2,121` eq) | `0.0632s` / `6.0e-15` | `0.0523s` / `4.4e-15` (diff `1.60`) | `0.7481s` / `7.5e-12` (diff `1.60`) | `3.9597s` / `2000` / `1.0` (fail) | `0.1551s` / `2.5e-14` (diff `1.60`) |

Pure CSR `AMGCL` vs `PardisoLDLT` multi-point scaling (`untrimmed p=2`):

| Grid | DOFs | `nnz` | `AMGCL` (s) | `AMGCL` rel_res | `PardisoLDLT` (s) | RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| `6×6×6` | 1,536 | 353,736 | 0.4653 | `4.07e-12` | 0.0410 | 539.8 |
| `9×9×9` | 3,993 | 1,058,841 | 1.3324 | `3.38e-11` | 0.1817 | 601.8 |
| `12×12×12` | 8,232 | 2,359,296 | 2.9770 | `3.84e-11` | 0.1852 | 760.4 |
| `13×13×13` | 10,125 | 2,956,581 | 4.3378 | `7.14e-11` | 0.4867 | 906.9 |

Empirical CSR scaling fit for `AMGCL`: `t = 9.365e-5 * N^1.157`.

### 3.3 Iterative Solver Spike Findings (`AMGCL` & `IgaApplication.AdditiveSchwarzPreconditioner`)

1. **`AMGCLSolver` Production Configuration Confirmed**:
   - `wht_iga`'s `DEFAULT_AMGCL_PARAMS` (`krylov_type="cg"`, `smoother_type="ilu0"`, `coarsening_type="aggregation"`, `use_block_matrices_if_possible=False`) and `ScalingSolver(AMGCL_CG_ILU0)` converged cleanly on every tested model (`rel_l2_diff <= 4.6e-11` with extension, `<= 3.6e-10` without extension), using **`< 50%` of `PardisoLU` memory** at `10,125` DOFs (`254.1 MB` vs `568.0 MB` in standalone Kratos probe).
   - **Pitfall discovered**: Setting `use_block_matrices_if_possible=True` with `smoother_type="spai0"` and `coarsening_type="smoothed_aggregation"` (`AMGCL_CG_SPAI0_Block`) stalls at `0.81%` error (`8.10e-3`) on trimmed+extension meshes and `100.5%` error on unextended meshes due to scalar constraint diagonals (`mScaleFactor * I_slave`) breaking 3×3 block aggregation. Keep `use_block_matrices_if_possible=False`.
2. **`IgaApplication.AdditiveSchwarzPreconditioner` — Not Ready for Production (2 C++ Bugs Identified)**:
   - On **untrimmed** `p=2`, `CG + AdditiveSchwarz` converges (`2.04e-10` at `9³`), but requires **more iterations (`84`) than plain `DiagonalPreconditioner` (`67`) and `ILU0Preconditioner` (`47`)**. Inspecting Kratos C++ [`additive_schwarz_preconditioner.h:290-294`](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_solvers/additive_schwarz_preconditioner.h#L290-L294) shows why: `M_invert = omega * M_invert;` is applied to `ZeroMatrix` *before* `MathUtils<double>::InvertMatrix(M, M_invert, M_det, -1e-6);`, so the 3D relaxation weight `omega = 1/18` is never applied to the overlapping block inverses.
   - On **trimmed** `p=2`, `CG + AdditiveSchwarz` **diverges** (`2,000` iters, `rel_l2_diff = 30.17` with extension, `4.00e14` without extension). With `ext=True`, `ApplyConstraints` produces `T^T A T + mScaleFactor * I_slave` (slave CSR rows retain only the diagonal entry), while `ProvideAdditionalData` still queries `r_geometry[j].GetDofs()` (which reference slave `EquationId`s) and uses `ForwardFind` expecting full local block entries; without extension (`ext=False`), cut-cell slivers destabilize `InvertMatrix(..., -1e-6)`.
3. **`CG_Diag` / `Scaling_CG_Diag` vs `CG_ILU0`**:
   - `CG_Diag` and `Scaling_CG_Diag` reliably converge (`151–156` iters, `1.0e-8–6.3e-8`) even on trimmed meshes, making `Scaling_CG_Diag` a solid fallback when `AMGCL` is unavailable, whereas `CG_ILU0` diverges (`rel_l2_diff = 2.95–212.6`) on trimmed systems.
4. **Extension Necessity**:
   - Without `extension` (`ext=False`), `PardisoLDLT` silently failed (`rel_l2_diff = 1.000`, `max_abs_diff = 52,645.9`), and tiny sliver modes (`~1e-12`) caused `1.60` relative solution sensitivity across `PardisoLU`, `PardisoLDLT`, `AMGCL`, and `splu` despite `1e-15` linear residuals. Enabling `theta=0.2` extension eliminates all solver instability.

---

## 4. `fbar` (`residual_tangent_sparse` & `solve_static`: `pypardiso` vs `superlu`) on E04 & V07 Sizing

| Case | Grid | DOFs (Free) | `nnz` | `build` (s) | `res_tan_sparse` (s) | 1-Solve `pardiso` / `superlu` (s) | `solve_static` `pardiso` / `superlu` (s) | Rel Diff |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| E04 Cook `4×4×1` | `(4,4,1)` | 324 (180) | 46,656 | 0.0263 | 0.0020 | 0.0044 / 0.0021 | 0.550 / 0.310 (65 it) | `8.2e-16` |
| E04 Cook `8×8×1` | `(8,8,1)` | 900 (540) | 156,816 | 0.0474 | 0.0057 | 0.0070 / 0.0076 | 2.489 / 1.515 (73 it) | `9.5e-16` |
| E04 Cook `16×16×1` | `(16,16,1)` | 2,916 (1,836) | 571,536 | 0.1639 | 0.0197 | **0.0153** / 0.0535 (`3.49×`) | **12.066** / 14.654 (86 it) | `2.1e-15` |
| V07 Kirsch `h0=0.5` | `(6,16,2)` | 1,824 (1,008) | 372,708 | 0.1396 | 0.0157 | 0.0151 / 0.0167 | 3.229 / 2.607 (8 it) | `9.5e-13` |
| V07 Kirsch `h1=0.25` | `(12,32,2)` | 5,880 (3,536) | 1,346,688 | 0.5719 | 0.0848 | **0.0945** / 0.1627 (`1.72×`) | 11.074 / 10.826 (8 it) | `4.3e-13` |
| V07 Kirsch `h2=0.125` | `(24,64,2)` | 20,904 (13,200) | 5,109,048 | 1.7952 | 0.2204 | 0.1518 / — | 2.009 (4 it) / — | — |

Observations:
- Numba `residual_tangent_sparse` is remarkably fast: `19.7 ms` for E04 `16×16×1` and `220.4 ms` for V07 `h2=0.125` (`5.11M` nonzeros).
- Single-solve `PardisoFactor` beats `SuperLU` by **`3.49×`** at E04 `16×16` (`15.3 ms` vs `53.5 ms`) and **`1.72×`** at V07 `h1=0.25` (`94.5 ms` vs `162.7 ms`). Inside `fbar.solve_static`, however, `pypardiso.spsolve(...)` uses default nonsymmetric `mtype=11` + CSR copy/comparison overhead on every Newton iteration, compressing the speedup (`12.07 s` vs `14.65 s` on E04 `16×16`). Reusing `PardisoFactor(kind="auto")` or symbolic pattern reuse in `solve_static` is an immediate `~3×` linear-solve speedup.

---

## 5. Extrapolation to Coreform Baseline & Phase 6 1M DOF Gate

Empirical power-law fits ($t = c \cdot N^\alpha$ on measured `p=2` scaling points) and conservative $O(N)$ linear-bound check for stages measured over narrow ranges:

| Target Scale | `est_nnz` | `Pardiso` Est. RAM | `AMGCL` Est. RAM | Power-Law Total (`AMGCL`) | Linear-Stage Envelope (`AMGCL`) | Gate Status |
|---|---:|---:|---:|---:|---:|---|
| **Coreform GTM-seal** (`178,278` DOFs, `32,098` active, `3,420` cut) | 52.06M | `146.1 GB` (FAIL) | **`1.91 GB`** | **`146.7 s`** | `~206 s` (vs Coreform `1,018.4 s`) | Time & RAM OK (`AMGCL`) |
| **Phase 6 1M DOF** (`1,000,000` DOFs, `310,000` cells, `15,000` cut) | 292.01M | `1,258.4 GB` (FAIL) | **`10.71 GB`** | **`1,026.2 s`** (`17.1 min`) | **`~1,374 s`** (`22.9 min`) | `< 30 min`, `< 64 GB` (w/ `AMGCL` + parallel mesh) |

*Caveat*: On `build_trimmed_table`, only two 8-thread points (`542` and `650` cut cells, `6.05 s -> 6.13 s`) were available above `_PARALLEL_MIN_CELLS=400`, producing an artificially flat power-law exponent (`alpha = 0.074`, projecting `7.74 s` for 15k cut cells). The linear-stage envelope column corrects this using the measured `9.44 ms/cut cell` (`8` threads) -> `~141.5 s` for 15k cut cells and linear cell scaling for `kratos_build` (`~209.2 s`). Even under the conservative linear-stage envelope, total time (`~1,374 s = 22.9 min`) sits under `1,800 s` (`30 min`), while **Pardiso is strictly infeasible in 64 GB RAM above ~100k–130k DOFs (`146 GB` at 178k DOFs, `1.26 TB` at 1M DOFs)**.

---

## 6. Phase 6 Hotspot Ranking & Action Plan

Ranked by projected share at `1M` DOFs (`AMGCL` path):

1. **Rank 1 — Iterative Linear Solve (`AMGCL_CG_ILU0`: `819.1 s`, `~60–80%` share)**:
   - Bottleneck: Scalar `AMGCL` (`use_block_matrices_if_possible=False`) + `ilu0` on `~292M` nonzeros.
   - Action: Benchmark block-aware / geometric-IGA smoothers (`ScalingSolver` + tuned `AMGCL` smoother `chebyshev`/`damped_jacobi` with coordinate nullspace `provide_coordinates=True`, or fix `omega=1/18` and slave-row handling in `IgaApplication.AdditiveSchwarzPreconditioner` in C++).
2. **Rank 2 — Element Creation Loop (`attach_cell_elements`: `119.0–167.2 s`, `~11–12%` share; `~1,850 elem/s`)**:
   - Bottleneck: Per-cell Python loop calling `cell_gauss_points` + `CreateQuadraturePointGeometries` + `CreateNewElement`.
   - Action: Precompute tensor Gauss tables once per uniform grid (`leggauss` hoisting), batch C++ element/QP-geometry builder in `WhtIgaApplication`, and implement `ExtendedSplineSolidElement` (eliminating MPC ties).
3. **Rank 3 — Kratos Element Assembly & MPC Ties (`kratos_build` + `apply_extension_ties`: `32–232 s`, `~3–17%` share)**:
   - Bottleneck: Python per-DOF `add_linear_tie` loop and `T^T K T` triple matrix product inside `ApplyConstraints` (~`2.8×` build cost).
   - Action: Direct `ExtendedSplineSolidElement` C++ extraction operator $\tilde{\mathbf{C}}^e$ eliminates `LinearMasterSlaveConstraint` and `T^T K T` completely; also reuse dependent-instance QP data (§5.9.3).
4. **Rank 4 — Cut-Cell Moment Fitting (`build_trimmed_table`: `~7.7–141.5 s` at 8T, up to `10%+` for dense CAD)**:
   - Bottleneck: Python `_cell_moments` / `_cross_section_moments` loop + per-cell `scipy.optimize.nnls` (`~9.4 ms/cut cell` on 8 threads at `q=4`).
   - Action: Move `_cell_moments` / `_cross_section_moments` to `@njit` (plan §2.2), lower `_PARALLEL_MIN_CELLS` from `400` to `~100` for `p >= 2`, and tune `q=4` cache reuse.
5. **Rank 5 — `make_nurbs_volume` & `build_extension` (`~4–48 s`, `< 5%` share)**:
   - Action: Batch node creation in C++ when scaling past 500k nodes; `classify_basis` (`@njit parallel`) is already sub-second (`< 2 ms` at 10k DOFs).
6. **Quick Win for Phase 5/6 `fbar`**:
   - Replace `pypardiso.spsolve(..., solver=pardiso)` inside `fbar.solve_static` and `fbar_trimmed.solve_static` with `wht_iga.linalg.PardisoFactor(kind="auto")` (symmetric `mtype=-2`), capturing the measured `1.7–3.5×` single-solve speedup over SuperLU without default `mtype=11` overhead.
