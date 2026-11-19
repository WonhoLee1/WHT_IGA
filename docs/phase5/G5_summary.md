# Phase 5 Gate (G5) & Verification Summary

> **Scope & Provenance:** Summary compiled strictly from spike reports (`S31.md`–`S33.md`, `S40.md`–`S43.md`, `S46_fbar_trimmed.md`, `S47_hardening_peeq.md`), verification reports (`V07.md`, `V08.md`, `E03.md`, `E05.md`), verification specs (`e03.json`, `e04.json`, `e05.json`, `v08.json`), pytest runs (`test_fbar_pipeline_wiring.py`, `test_fbar_trimmed.py`, `test_p5_v08_h2_nlgeom.py`, `test_p5_v08_ghost_check.py`, `test_p5_e03_p2_coupling.py`, `test_p5_e05_nlgeom_hardening.py`), and external session handoff ([g5_claude_handoff_20261117.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/dev_log/g5_claude_handoff_20261117.md)).

---

## 1. Benchmark Status Overview

| Benchmark | Description | G5 Gate Target | Current Status | Key Verified Metric | Primary Open Gap |
|---|---|---|---|---|---|
| **V07** | Kirsch plate w/ hole ($\nu=0.49999$) | `e_u`, `e_sigma_xx` & convergence slopes | **PASS** (3-Level Body-Fitted & Trimmed F-bar) | Body-fitted `h2=0.125`: `e_u=8.097e-05`, `e_sigma_xx=1.379e-02`, OLS `slope_u=2.492`; Trimmed `h1=0.25`: `e_u=3.275e-04`, `e_sigma_xx=2.587e-02`, `slope_u=3.087` | `P0` stress slope (`0.852` body-fitted / `1.120` trimmed) meets doc `> 0.5` but `<` plan `2.2`; `fbar_trimmed` pipeline wiring pending |
| **E04** | Cook's membrane ($\nu=0.49999$, NLGEOM) | `u_y` within `1%` of `6.98946 mm` | **Verified Pass** | `p=2, 16x16` Numba F-bar/P0: `u_y = 7.03806 mm` (`+0.70%`, CSR+Pardiso `21.25s`) | No cut-cell or Kratos strategy wiring |
| **V08** | Spherical indenter Hertz-Coulomb contact | `F_N` within `4%` (`481,579.3 lbf`), `F_T` RMS `<= 4.5%` | **PASS** (`h=0.125` R1 & `h=0.25` R2 pytest; `h=0.125` R2 external handoff) | NLGEOM R1 `h=0.125, n=48, eps_n=1e10`: `+0.24%` (`n=96, eps_n=1e11`: `+0.95%`); R2 `0.15<=t<=0.70`: `h=0.25` `rel_RMS=(3.25%, 2.67%)` (pytest), `h=0.125` `rel_RMS=(1.67%, 0.94%)` (Claude 10:40 handoff) | `h=0.125` R2 not yet in pytest gate suite; `ghost_alpha=0.0` `IndexError` in `KratosBodyOperator`; `p_max` scatter (record-only) |
| **E03** | Glass-to-metal seal 13-solid thermal cycle | 3-step convergence, `<=2x` Coreform inc, qualitative PEEQ/stress | **PASS** (`p=2` Extended + Annular RBE3) | `(9,6,6) p=2` NLGEOM+T-hardening+RBE3: `15/15` inc (`5/step <= 2x15`), `0` cutbacks, `||W_op u_surf|| <= 8.97e-17`, `571.0 MB` peak RSS | `0.5 mm` full-resolution mesh (`146,461` background cells) deferred to Phase 6 parallel/AMGCL |
| **E05** | Door-lock 3-part frictional contact | `link RF3 = -46.946 kN` (`±5%`) | **PARAMETRIC PASS** (3D Curved Contact + NLGEOM + J2; Quant Gate Blocked) | Curved patch `gap_max_err=1.34e-4`; `Contact.STEP` (`mu=0.3`) `eq_rel=1.33e-11`; E05 5-solid NLGEOM+J2 `uz=1.0mm` `PEEQ_max=0.006017`, `RF3_j2=-6.7490 kN` (`0.9691x` elastic), `bal_z=4.16e-8` | External Coreform E05 STEP CAD file missing -> `-46.946 kN` quantitative gate **BLOCKED** |

---

## 2. Detailed Benchmark Breakdown

### 2.1 V07 — Kirsch Infinite Plate with Circular Hole ($\nu = 0.49999$)
- **Verified Pass — 3-Level Body-Fitted & Immersed/Trimmed F-bar** ([V07.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V07.md), [S42.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S42.md), [S46_fbar_trimmed.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S46_fbar_trimmed.md)):
  - **Kratos mixed element incompatibility (Check A):** `SmallDisplacementMixedVolumetricStrainElement3D4N` crashes with Access Violation (`returncode=3221225477`, `0xC0000005`) in `Initialize()`, and `SmallDisplacementMixedStrainElement3D4N` crashes in `strategy.Solve()` on `p=2` `QuadraturePointGeometry`.
  - **Body-fitted `fbar` 3-level convergence (`MappedQuarterPlate`, `p=2`, `test_p5_v07_kirsch_precursor.py` `4 passed in 113.07s`):**
    - Analytical oracle (`test_kirsch_analytical_verification`): hole boundary traction/stress `< 1e-13`, peak `sigma_xx = 30.0` at `(0, a)` (`< 1e-13`), `div(sigma) = 0` (`< 1e-8`).
    - `h0 = 0.5` (`(6, 16, 2)`, 608 CPs / 192 cells, `6.93s`): volume rel error `3.279e-08`, `e_u = 2.562389e-03 <= 0.01`, `e_sigma_xx = 4.492931e-02 <= 0.2`.
    - `h1 = 0.25` (`(12, 32, 2)`, 1,960 CPs / 768 cells, `33.76s`): volume rel error `2.117e-09`, `e_u = 4.030790e-04 <= 0.005`, `e_sigma_xx = 2.398510e-02 <= 0.1` (vs Kratos `SmallDisplacementElement3D4N`: `e_u = 2.122799e-02`, `e_sigma_xx = 2.694238e+01`; 52.7x `e_u` and 1,123x `e_sigma_xx` improvement).
    - `h2 = 0.125` (`(24, 64, 2)`, 6,968 CPs / 20,904 DOFs / 3,072 cells, `17.30s`): volume rel error `1.344e-10`, `e_u = 8.096688e-05 <= 0.002`, `e_sigma_xx = 1.379114e-02 <= 0.05`.
    - 3-level OLS slopes: `slope_u = 2.492` (`h0->h1`: `2.668`, `h1->h2`: `2.316`; meets plan `>= 2.3` and doc `> 0.5`), `slope_sigma_xx = 0.852` (`h0->h1`: `0.906`, `h1->h2`: `0.798`; meets source doc `> 0.5`).
  - **Immersed/Trimmed F-bar on orthogonal `BoxSpline` (`wht_iga/analyses/fbar_trimmed.py`, `tests/unit/test_fbar_trimmed.py` `11 passed in 18.87s`, [S46_fbar_trimmed.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S46_fbar_trimmed.md)):**
    - Pressure-cell merging (`eta_thresh = 0.25`) + symmetric Jacobi scaling (`D K D`, no ghost/extension) + FD tangent verification (`rel_err = 1.551e-10` at `nu=0.3`, `2.457e-10` at `nu=0.49`) + affine stress patch (`err_u = 4.17e-12`, `err_sigma = 4.47e-14`).
    - `h = 0.5` (`1,200` free DOFs, `126` active / `78` cut cells, `122` groups): `e_u = 2.7822e-03 <= 0.01`, `e_sigma_xx = 5.6236e-02 <= 0.2`, `cond(D K D) = 2.290e+05` (`0.355x` vs unmerged `eta_thresh=0.0` `6.453e+05`).
    - `h = 0.25` (`3,888` free DOFs, `496` active / `279` cut cells, `490` groups): `e_u = 3.2748e-04 <= 0.005`, `e_sigma_xx = 2.5867e-02 <= 0.1`, `cond(D K D) = 1.000e+06` (`0.632x` vs unmerged `1.582e+06`).
    - Convergence slopes (`h0 -> h1`): `slope_u = 3.087`, `slope_sigma_xx = 1.120` (unmerged cell-wise P0 `eta_thresh=0.0`: `h=0.25` `e_u = 2.8779e-04`, `e_sigma_xx = 2.5037e-02`, slopes `2.737 / 1.062`).
- **Remaining Scope / Open Gaps:**
  - `P0` cell/group-constant pressure limits stress $L^2$ slope to $O(h)$ (`0.852` body-fitted, `1.120` trimmed; satisfies source doc `> 0.5` but `<` master-plan `>= 2.2` / Coreform `2.385`).
  - `FullIgaOption(kind="full_iga_immersed", locking_mitigation="fbar")` pipeline/User API wiring to `fbar_trimmed.py` is pending.

### 2.2 E04 — Cook's Membrane ($\nu = 0.49999$, TL Neo-Hookean)
- **Verified Pass** ([S40.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S40.md), [e04.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/e04.json), [test_fbar_pipeline_wiring.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_pipeline_wiring.py)):
  - **Cook's Membrane Gate (`wht_iga/analyses/fbar.py`, Simo–Taylor–Pister F-bar/P0):**
    - `p=2, 16x16` (`2,916` DOF): `u_y = 7.03806 mm`, error **`+0.70%`** vs `6.98946 mm` (passes G5 `1%` gate).
    - `p=1` series: `4x4` `6.21865 mm` (`-11.0%`), `8x8` `6.72656 mm` (`-3.76%`), `16x16` `6.89904 mm` (`-1.29%`), Richardson limit `~6.987 mm` (`-0.03%`).
    - `p=2` series: `4x4` `7.37756 mm` (`+5.55%`), `8x8` `7.14649 mm` (`+2.25%`), `12x12` `7.07429 mm` (`+1.21%`).
  - **Pipeline & User API Wiring:** `recommend_discretization` (`nu >= 0.49`), `FullIgaOption`, `solve_with_discretization_method`, `solve_fbar_box` (`C^(p-2)` `CookSpline`), and `StudySpec -> compile_plan -> ReferenceNumbaBackend` all route to `fbar.solve_static` and match direct calls within `1e-12`.
- **Remaining Scope / Open Gaps:**
  - Immersed/cut cells and Kratos strategy coupling (contact, explicit dynamics, C++ UEL) are not yet implemented. (CSR sparse assembly + `C^(p-2)` support in `fbar.py` reduced `p=2, 16x16` runtime to `21.25s`, and `methods.py` / User API wiring is verified.)

### 2.3 V08 — Elastic Spherical Indenter Sliding on Block (Hertz–Coulomb)
- **Verified Pass** ([V08.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V08.md), [S31.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S31.md), [S32.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S32.md), [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md), [v08.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v08.json)):
  - **Case R0 (Rigid Sphere SSN, `wht_iga/analyses/contact.py`):** Flat punch (`49.999500 lbf`, `1e-11` vs analytic), separation (`n_active=0`), Coulomb patch (`|F_T|/F_N = 0.30000000`), and quarter-block V08 R0 (`h2=0.125`, `23,460` free DOFs, `eps_n=3e8`: `F_N = 472,601.8 lbf`, `-1.86%` vs Hertz, `+0.14%` vs Sneddon; `h1` exponent `1.477` vs `1.500`).
  - **Cases R1 & R2 Small-Strain (`sym="half"`, `h1=0.25, n=24, eps_n=1e9`):** R1 `F_N = 486,248.6 lbf` (`+0.97%` vs Hertz, eq error `5.08e-07`); R2 (`0.15 <= t <= 0.40`, `mu=0.3`) `rel_RMS = (F_N: 2.63%, F_T: 2.60%)`. S41 sweep showed `+0.97%` relied on error cancellation (small-strain extrapolated limit `~+5.0%`).
  - **Cases R1 & R2 NLGEOM (`wht_iga/analyses/contact_nlgeom.py`, `TotalLagrangianElement3D4N` + `HyperElastic3DLaw`):**
    - Patch tests (compression `F_N=33.786323`, separation `n_active=0`, friction `|F_T|/F_N=0.299999`) passed.
    - **R1 NLGEOM `h=0.25`** (`n=48, eps_n=1e10`, `303.35s`): `F_N(0.1) = 496,764.5 lbf` (**`+3.15%`** vs Hertz, passes `4%` gate; eq error `3.16e-6`).
    - **R1 NLGEOM `h=0.125` (pytest gate `test_p5_v08_h2_nlgeom.py`):** `n_soup=48, eps_n=1e10` (`1368.98s`, `peak_step_rss=1.666 GiB`, `post_free_rss=0.259 GiB`): `F_N(0.1) = 482,739.8 lbf` (**`+0.24%`** vs Hertz, `+2.29%` vs Sneddon, slope `1.4739`, `max_eq_rel = 5.14e-6`). Prior S41 `eps_n=1e11` runs: `n=48` `484,352.0 lbf` (`+0.58%`), `n=96` `486,150.8 lbf` (`+0.95%`), $n\to\infty$ estimate `~+1.07%` (Abaqus `h2` `1.07%`, Coreform `3.33%`).
    - **R2 NLGEOM `h=0.25` full window (pytest gate `test_p5_v08_h2_nlgeom.py`):** `n_soup=24, mu=0.3, eps_n=1e10, eps_t=1e9`, `46` steps (`1073.45s`, `peak_step_rss=0.430 GiB`): `0.15 <= t <= 0.40` `rel_RMS = (F_N: 3.42%, F_T: 2.75%)`, `0.15 <= t <= 0.70` **`rel_RMS = (F_N: 3.25%, F_T: 2.67%)`** (passes `4.5%` gate), `eq_tot_T = 5.94e-8`.
    - **R2 NLGEOM `h=0.125` (external Claude session handoff):** Per [g5_claude_handoff_20261117.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/dev_log/g5_claude_handoff_20261117.md) (Claude 2026-11-19 10:40 handoff, `scratchpad/v08/nlgeom_r2.py`, `n=48, eps_n=1e11, eps_t=1e10`, `46` steps, `8018s`, `n_fallback=0`, pre-patch TL ghost cells): `0.15 <= t <= 0.40` `rel_RMS = (F_N: 1.81%, F_T: 0.91%)`; `0.15 <= t <= 0.70` **`rel_RMS = (F_N: 1.67%, F_T: 0.94%)`**, `F_T/F_N` mean `0.2961` (`[0.2923, 0.2995]`), `eq_cut = (N: 3.8e-4, T: 5.4e-2)`.
    - **Ghost-cell hypothesis & `ghost_alpha` verification (`test_p5_v08_ghost_check.py`, `1 passed in 292.57s`, `peak_rss = 395.3 MiB`):**
      - Verified `R_cut_x - force_total[0] == -(f_int_s[zero_mask_x].sum())` (`id_err <= 3.00e-3`, rel `< 1e-7`, `fc_zero_x == 0.0`) on `h=0.25, n_soup=24`, `t = (0.03, 0.10, 0.40, 0.70)`.
      - Pre-patch `TL_ghost(1e-8)` (`TotalLagrangianElement3D4N`) leaks `+928.51 lbf` (`2x = 1,857.0 lbf`, `eq_cut_T = 1.301e-02`, `7.71x` cubic scaling `0.40 -> 0.70`), whereas current default `SD_ghost(1e-8)` (`SmallDisplacementElement3D4N`) leaks `+2.53 lbf` (`2x = 5.05 lbf`, `eq_cut_T = 3.541e-05`, `2.01x` linear scaling, `367x` reduction) with `F_N/F_T` changed by only `0.0120% / 0.0261%`.
      - `ghost_alpha = 0.0` skips `OUT` cells (`18` vs `45` elements) and pins `inactive_ctrl_mask` (`65` CPs) in `zero_mask`, but `KratosBodyOperator.eval()` raises `IndexError: index (524) out of range` because Kratos sizes `_A` to `330` active DOFs while `self._inv` has length `box.n_dofs = 525`.
- **Remaining Scope / Open Gaps:**
  - `h=0.125` R2 NLGEOM has not yet been promoted from the external Claude scratchpad run into `test_p5_v08_h2_nlgeom.py`.
  - Supporting `ghost_alpha = 0.0` cleanly in `KratosBodyOperator` requires mapping `self._inv` only over DOFs attached to active elements; `p_max` discretization scatter remains record-only per `v08.json`.

### 2.4 E03 — Glass-to-Metal Seal Thermal Cycling (13-Solid Multicell, T-Dependent J2)
- **Verified Pass (`p=2` Extended B-Spline + Annular RBE3 Coupling)** ([E03.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E03.md), [S33.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S33.md), [S47_hardening_peeq.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S47_hardening_peeq.md), [e03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/e03.json)):
  - **STEP Inspection (`glass_to_metal_seal_connector.stp`):** 13 solids (`SS304L_header`, `S8061_glass`, 11 `Alloy52_pins`), 429 faces, `watertight=True` (`free_edges_after=0`), B-Rep volume `2590.4895 mm³` (solid sum `2590.4877 mm³`, rel diff `< 1e-6`).
  - **Hardening & PEEQ Sync ([S47_hardening_peeq.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S47_hardening_peeq.md)):** `HardeningTable` (`HARDENING_CURVE=6` + `point_curve_margin` check, min `0.0872` across `-55..500 °C`) and incremental `PEEQ -> PLASTIC_DISSIPATION` sync verified (`test_s44_hardening_initial_strain.py` + `test_p5_thermal_plasticity.py`: `25 passed in 27.29s`).
  - **Coarse `p=1` Baseline (`(8, 6, 6)`, `S33.md`):** `n_fallback = 0`, `9/9` converged (`0` cutbacks, `485.4 MB` RSS), free thermal stress `< 5e-12 MPa`.
  - **`p=2` Extended B-Spline + Annular RBE3 (`(9, 6, 6)`, `test_p5_e03_p2_coupling.py` `3 passed in 315.96s`, Peak RSS `571.0 MB`):**
    - Trimmed model: `704` CPs (`416` master / `124` extended / `164` inactive, `max_abs_coef = 6.0`), `276` QP elements (`n_fallback = 0`), volume `2586.7213 mm³` (`0.145%` rel error), mounting face area `12.6376 mm²` (`0.250%` rel error).
    - Annular RBE3 (`build_rbe3` + QR pivot `cond(Ws) = 7.40`, `378` flattened `LinearMasterSlaveConstraint`s): free thermal expansion stress `4.223e-09 MPa < 1e-6 MPa` (`||W_op u_surf|| < 1e-11`, radial contraction `u_yz > 1e-3 mm`).
    - Full 3-step NLGEOM thermal cycle (`450 -> 20 -> -55 -> 200 °C`, `increments=[5, 5, 5]`, `218.56s`): **`15/15` converged, `0` cutbacks** (`5 <= 2 x 15` Coreform gate), `||W_op u_surf|| <= 8.97e-17`, Header `PEEQ_max: 3.897e-3 -> 4.421e-3 -> 4.421e-3`, Pins `PEEQ_max: 6.536e-3 -> 7.830e-3 -> 7.960e-3`, Glass `sigma_3,min = -490.974 -> -564.739 -> -292.392 MPa`, `sigma_1,max = +100.519 -> +116.116 -> +117.987 MPa`.
- **Remaining Scope / Open Gaps:**
  - `0.5 mm` full-resolution mesh (`146,461` background elements) sweep is deferred to Phase 6 parallel/AMGCL scale-up.

### 2.5 E05 — Door-Lock Assembly (3-Part Frictional Contact, NLGEOM & J2 Hardening)
- **Verified Parametric Pass (Quantitative `-46.946 kN` Gate Blocked by Missing External CAD)** ([E05.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E05.md), [S43.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S43.md), [e05.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/e05.json)):
  - **Contract & Parametric B-Rep (`test_p5_e05_door_lock_precursor.py` `3 passed in 4.42s`):** Extracted MMTS contract (`weak_steel`, `strong_steel`, 3 contact pairs, `RF3 = -46.946 kN ± 5%`); 5 watertight B-Rep solids (`hinge`, `link`, `pin_cyl`, `pin_cap_l`, `pin_cap_r`, all `free_edges = 0`, `watertight = True`); `u_z = +0.2 mm` SSN (`mu = 0.0` `link RF3 = -198.868 kN`, `Z_bal = 1.32e-15`; `mu = 0.3` `link RF3 = -201.145 kN`, `Z_bal = 0.00e+00`).
  - **3D Curved Contact + NLGEOM + 2-Point J2 Hardening (`wht_iga/analyses/contact_surface.py`, `test_p5_e05_nlgeom_hardening.py` `3 passed in 171.32s`):**
    - Fixed `NonlinearBodyOperator(..., ghost_alpha=0.0)` active-DOF `_P_act` injection mapping when `A.Size1() < box.n_dofs`.
    - Curved `DeformableTriangleSurfaceMaster` patch (`M=36`): `row_sum_err = 4.44e-16`, `xm_interp_err = 2.66e-15`, `gap_max_err = 1.34e-04`, open-boundary exclusion `overhang_gn = +0.8533 > 0`.
    - `wht_iga/resources/Contact.STEP` 3D curved NLGEOM contact (`p=2, q=2, ghost_alpha=0.0, mu=0.3`): `u_x = 0.7 -> 1.0 -> 1.3 mm`, `n_act = 62 -> 101 -> 112`, `F_N = 69.21 -> 301.66 -> 481.84 N`, `eq_rel <= 3.68e-11`, `RSS = 453.9 MB`.
    - E05 5-solid parametric NLGEOM + 2-point J2 hardening (`weak_steel` `1500->150000 MPa`, `strong_steel` `2000->200000 MPa`, `u_z = 0.25 -> 1.00 mm`, `mu=0.3`, `2,388` DOFs / `1,848` free): `PEEQ_max = 0.0 -> 0.0 -> 0.002657 -> 0.006017`, `RF3_j2 = -6.7490 kN` vs `RF3_el = -6.9644 kN` (`0.9691x`), `bal_z <= 4.37e-08`, `RSS = 446.5 MB`.
- **Remaining Scope / Open Gaps:**
  - **External STEP/ACIS file missing (BLOCKER):** Original Coreform E05 CAD geometry is not provided in `research-cf-IGA/user-guide/example-problems/door-lock-assembly/`, blocking quantitative `-46.946 kN` (`±5%`) validation.

---

## 3. Verified Test Suite Inventory

| Test File | Scope / Benchmark | Reported Result | Reference |
|---|---|---|---|
| `tests/unit/test_fbar_mean_dilatation.py` | E04 Cook gate (`p=2, 16x16`) & `fbar` unit tests (CSR + `C^(p-2)`) | `18 passed` (Cook gate `21.25s`) | [S40.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S40.md) |
| `tests/unit/test_fbar_trimmed.py` | V07 immersed/trimmed orthogonal `BoxSpline` F-bar/P0 & merging | `11 passed in 18.87s` | [S46_fbar_trimmed.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S46_fbar_trimmed.md) |
| `tests/unit/test_fbar_pipeline_wiring.py` | F-bar `methods.py` & `StudySpec` User API wiring (Kratos-free) | Unit test suite | [test_fbar_pipeline_wiring.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_pipeline_wiring.py) |
| `tests/kratos_spikes/test_p5_v07_kirsch_precursor.py` | V07 Kratos mixed crash check & 3-level body-fitted `fbar` Kirsch | `4 passed in 113.07s` | [V07.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V07.md) |
| `tests/kratos_spikes/test_p5_v08_r0_ssn.py` | V08 R0 rigid sphere SSN + punch/friction/separation patches | `4 passed in 47.32s` | [S31.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S31.md) |
| `tests/kratos_spikes/test_p5_v08_r1_r2_ssn.py` | V08 R1 & R2 small-strain deformable SSN (+ R0 & `test_isolation.py`) | `12 passed in 169.39s` | [S32.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S32.md) |
| `tests/kratos_spikes/test_p5_v08_nlgeom_ssn.py` | V08 R1 & R2 NLGEOM (`h=0.25`, `TotalLagrangianElement3D4N`) + patches | Verified in S41 | [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md) |
| `tests/kratos_spikes/test_p5_v08_h2_nlgeom.py` | V08 `h=0.125` R1 NLGEOM & `h=0.25` R2 NLGEOM (`0.15<=t<=0.70`) | `2 passed in 2445.31s` | [V08.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V08.md) |
| `tests/kratos_spikes/test_p5_v08_ghost_check.py` | V08 R2 ghost-cell leak identity (`TL` vs `SD` `1e-8`) & `ghost_alpha=0.0` probe | `1 passed in 292.57s` | [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md) |
| `tests/kratos_spikes/test_s44_hardening_initial_strain.py` | E03 `HardeningTable` & PEEQ sync (+ `test_p5_thermal_plasticity.py`) | `25 passed in 27.29s` | [S47_hardening_peeq.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S47_hardening_peeq.md) |
| `tests/kratos_spikes/test_p5_e03_reduced_step.py` | E03 13-solid STEP, trimmed quadrature, `p=1` 3-step thermal + NLGEOM | `4 tests` | [S33.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S33.md) |
| `tests/kratos_spikes/test_p5_e03_p2_coupling.py` | E03 `p=2` extended B-spline + annular RBE3 3-step NLGEOM | `3 passed in 315.96s` | [E03.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E03.md) |
| `tests/kratos_spikes/test_p5_e05_door_lock_precursor.py` | E05 5-solid watertight B-Rep & 3-part elastic SSN (`mu=0, 0.3`) | `3 passed in 4.42s` | [S43.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S43.md) |
| `tests/kratos_spikes/test_p5_e05_nlgeom_hardening.py` | E05 curved surface contact + `Contact.STEP` NLGEOM + 5-solid J2 hardening | `3 passed in 171.32s` | [E05.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E05.md) |
