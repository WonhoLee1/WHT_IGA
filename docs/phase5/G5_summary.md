# Phase 5 Gate (G5) & Verification Summary

> **Scope & Provenance:** Summary compiled strictly from spike reports (`S31.md`–`S33.md`, `S40.md`–`S43.md`), verification specs (`e03.json`, `e04.json`, `v08.json`), and unit test files (`test_fbar_pipeline_wiring.py`).

---

## 1. Benchmark Status Overview

| Benchmark | Description | G5 Gate Target | Current Status | Key Verified Metric | Primary Open Gap |
|---|---|---|---|---|---|
| **V07** | Kirsch plate w/ hole ($\nu=0.49999$) | `e_u`, `e_sigma_xx` & convergence slopes | **Precursor Pass** (Body-fitted) | `h1=0.25`: `e_u=4.031e-04`, `e_sigma_xx=2.399e-02`, `slope_u=2.668` | Immersed/trimmed F-bar unsupported; `slope_sigma_xx=0.906` (< plan `2.2`, > doc `0.5`) |
| **E04** | Cook's membrane ($\nu=0.49999$, NLGEOM) | `u_y` within `1%` of `6.98946 mm` | **Verified Pass** | `p=2, 16x16` Numba F-bar/P0: `u_y = 7.03806 mm` (`+0.70%`, CSR+Pardiso `21.25s`) | No cut-cell or Kratos strategy wiring |
| **V08** | Spherical indenter Hertz-Coulomb contact | `F_N` within `4%` (`481,579.3 lbf`), `F_T` RMS `<= 4.5%` | **R1 Verified Pass (h=0.125)** / R2 Coarse-Verified | NLGEOM R1 `h=0.125, n=96, eps_n=1e11`: `+0.95%` (n->inf about `+1.07%`, Abaqus h2 `1.07%`); R2 `h=0.25` `rel_RMS=(3.42%, 2.75%)` (`n=24`), `(5.09%, 4.23%)` (`n=48`) | `h=0.125` R2 NLGEOM pending; R2 tested on `0.15<=t<=0.40` (vs `0.70`) |
| **E03** | Glass-to-metal seal 13-solid thermal cycle | 3-step convergence, `<=2x` Coreform inc, qualitative PEEQ/stress | **Precursor Pass** (Coarse) | `(8,6,6) p=1` NLGEOM+hardening: `9/9` inc, `0` cutbacks, `<1 GB` RSS | Full-resolution / `p=2` mesh run pending; 3-2-1 BC vs distributing coupling |
| **E05** | Door-lock 3-part frictional contact | `link RF3 = -46.946 kN` (`±5%`) | **Precursor Pass** (Parametric Elastic) | 5 watertight solids; `u_z=+0.2 mm` SSN (`mu=0, 0.3`) global eq `<=1.32e-15` | External STEP file missing; plasticity, `5 mm` NLGEOM, 3D trimmed projection pending |

---

## 2. Detailed Benchmark Breakdown

### 2.1 V07 — Kirsch Infinite Plate with Circular Hole ($\nu = 0.49999$)
- **Verified / Precursor Pass** ([S42.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S42.md)):
  - **Kratos mixed element incompatibility (Check A):** `SmallDisplacementMixedVolumetricStrainElement3D4N` crashes with Access Violation (`returncode=3221225477`, `0xC0000005`) in `Initialize()`, and `SmallDisplacementMixedStrainElement3D4N` crashes in `strategy.Solve()` on `p=2` `QuadraturePointGeometry`.
  - **Body-fitted `fbar` (`MappedQuarterPlate`, `p=2`, Check B):**
    - Analytical oracle (`test_kirsch_analytical_verification`): hole boundary traction/stress `< 1e-13`, peak `sigma_xx = 30.0` at `(0, a)` (`< 1e-13`), `div(sigma) = 0` (`< 1e-8`).
    - `h0 = 0.5` (`(6, 16, 2)`, 608 CPs / 192 cells): volume rel error `3.279e-08`, `e_u = 2.562e-03 < 0.01`, `e_sigma_xx = 4.493e-02 < 0.2`.
    - `h1 = 0.25` (`(12, 32, 2)`, 1,960 CPs / 768 cells): volume rel error `2.117e-09`, `e_u = 4.031e-04 < 0.005`, `e_sigma_xx = 2.399e-02 < 0.1` (vs Kratos `SmallDisplacementElement3D4N`: `e_u = 2.123e-02`, `e_sigma_xx = 2.694e+01`; 52.7x `e_u` and 1,123x `e_sigma_xx` improvement).
    - Convergence slopes (`h0 -> h1`): `slope_u = 2.668` (meets plan `>= 2.3` and doc `> 0.5`), `slope_sigma_xx = 0.906` (meets doc `> 0.5`).
- **Remaining Scope / Open Gaps:**
  - `fbar` cell-constant pressure (`Jbar`) limits stress $L^2$ slope to $O(h)$ (`0.906` vs master-plan `>= 2.2` / Coreform `2.385`). Gate selection (`> 0.5` vs `>= 2.2`) needs confirmation.
  - Immersed/trimmed F-bar is pending (`fbar.py` lacks cut-cell support; `FullIgaOption(kind="full_iga_immersed", locking_mitigation="fbar")` raises `ValidationError`).

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
- **Verified Pass / Precursor Pass** ([S31.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S31.md), [S32.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S32.md), [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md), [v08.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v08.json)):
  - **Case R0 (Rigid Sphere SSN, `wht_iga/analyses/contact.py`):** Flat punch (`49.999500 lbf`, `1e-11` vs analytic), separation (`n_active=0`), Coulomb patch (`|F_T|/F_N = 0.30000000`), and quarter-block V08 R0 (`h2=0.125`, `23,460` free DOFs, `eps_n=3e8`: `F_N = 472,601.8 lbf`, `-1.86%` vs Hertz, `+0.14%` vs Sneddon; `h1` exponent `1.477` vs `1.500`).
  - **Cases R1 & R2 Small-Strain (`sym="half"`, `h1=0.25, n=24, eps_n=1e9`):** R1 `F_N = 486,248.6 lbf` (`+0.97%` vs Hertz, eq error `5.08e-07`); R2 (`0.15 <= t <= 0.40`, `mu=0.3`) `rel_RMS = (F_N: 2.63%, F_T: 2.60%)`. S41 sweep showed `+0.97%` relied on error cancellation (small-strain extrapolated limit `~+5.0%`).
  - **Cases R1 & R2 NLGEOM (`wht_iga/analyses/contact_nlgeom.py`, `TotalLagrangianElement3D4N` + `HyperElastic3DLaw`):**
    - Patch tests (compression `F_N=33.786323`, separation `n_active=0`, friction `|F_T|/F_N=0.299999`) passed.
    - **R1 NLGEOM** (`h=0.25, n=48, eps_n=1e10`, `303.35s`): `F_N(0.1) = 496,764.5 lbf` (**`+3.15%`** vs Hertz, passes `4%` gate; eq error `3.16e-6`).
    - **R2 NLGEOM** (`h=0.25, mu=0.3, eps_n=1e10, eps_t=1e9`, `0.15 <= t <= 0.40`): `n=24` (`585.97s`) gives `rel_RMS = (F_N: 3.42%, F_T: 2.75%)` (passes `4.5%` gate); `n=48` (`901.49s`) gives `rel_RMS = (F_N: 5.09%, F_T: 4.23%)`.
  - **R1 NLGEOM h=0.125** (S41, 2026-11-19): `n=48` `484,352.0 lbf` (`+0.58%`), `n=96` `486,150.8 lbf` (`+0.95%`), `eps_n=1e11` (penalty converged); small-strain at the same setting `+4.38%`. Converged estimate about `+1.07%` (Abaqus h2 `1.07%`, Coreform `3.33%`).
- **Remaining Scope / Open Gaps:**
  - `h=0.125` R2 NLGEOM fine run is pending (R1 h=0.125 done, S41).
  - Sliding window verified in S32/S41 is `0.15 <= t <= 0.40` (vs `v08.json` full window `0.15 <= t <= 0.70`).

### 2.4 E03 — Glass-to-Metal Seal Thermal Cycling (13-Solid Multicell, T-Dependent J2)
- **Verified Precursor Pass** ([S33.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S33.md), [e03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/e03.json)):
  - **STEP Inspection (`glass_to_metal_seal_connector.stp`):** 13 solids (`SS304L_header`, `S8061_glass`, 11 `Alloy52_pins`), 429 faces, `watertight=True` (`free_edges_after=0`), B-Rep volume `2590.4895 mm³` (solid sum `2590.4877 mm³`, rel diff `< 1e-6`).
  - **Reduced Coarse Model (`(8, 6, 6) p=1`, 288 cells, 180 active, `q=2`):**
    - Trimmed quadrature: `n_fallback = 0`, `6,962` CUT QPs, integrated volume `2586.7213 mm³` (`1.454e-3` rel error).
    - Free thermal expansion (`3/3` increments): max Cauchy stress `< 5e-12 MPa`, `PLASTIC_DISSIPATION == 0.0`.
    - Initial-yield small-strain (`increments=[3, 3, 3]`): `9/9` converged, `0` cutbacks, RSS `485.4 MB` (`< 1 GB`).
    - Hardening + NLGEOM (`TotalLagrangianElement3D4N`, `finite_strain=True`, `increments=[3, 3, 3]`): `9/9` converged, `0` cutbacks; monotonic `PEEQ_max` in header (`2.750e-3 -> 3.060e-3 -> 3.060e-3`) and pins (`9.314e-3 -> 1.188e-2 -> 1.230e-2`); glass `sigma_1,max` rises `+13.22 -> +38.63 -> +97.66 MPa`.
- **Remaining Scope / Open Gaps:**
  - Full-resolution mesh & `p=2` run is pending (`e03.json` Coreform baseline: `146,461` background elements / `0.5 mm` mesh, `15` increments/step).
  - Current BC uses `fix_isostatic_321` (6-DOF minimal restraint) rather than a distributing coupling on the annular mounting surface.

### 2.5 E05 — Door-Lock Assembly (3-Part Frictional Contact)
- **Verified Precursor Pass** ([S43.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S43.md)):
  - **Contract & Parametric B-Rep:** Extracted MMTS contract (`weak_steel`, `strong_steel`, 3 contact pairs, `RF3 = -46.946 kN ± 5%`); generated 5 watertight B-Rep solids (`hinge`, `link`, `pin_cyl`, `pin_cap_l`, `pin_cap_r`, all `free_edges = 0`, `watertight = True`).
  - **Reduced 3-Part Elastic SSN Contact (`p=2` BoxSpline, `u_z = +0.2 mm`):** `mu = 0.0` (`6` iters, `rel_res = 7.57e-08`, `link RF3 = -198.868 kN`, Z eq error `1.32e-15`) and `mu = 0.3` (`4` iters, `rel_res = 6.32e-09`, `link RF3 = -201.145 kN`, Z eq error `0.00e+00`).
- **Remaining Scope / Open Gaps:**
  - **External STEP/ACIS file missing:** Original CAD geometry is not provided in `door-lock-assembly/`, blocking quantitative `-46.946 kN` (`±5%`) validation.
  - **Plastic hardening, `5 mm` NLGEOM, and 3D trimmed surface contact:** Full E05 requires coupling 2-point linear-extrapolated J2 hardening with `TotalLagrangianElement3D4N` / `contact_nlgeom.py` and 3D trimmed B-Rep triangle contact projection (`SFM3D3` / `Triangle3D3`) instead of `DeformableBoxFaceMaster`.

---

## 3. Verified Test Suite Inventory

| Test File | Scope / Benchmark | Reported Result | Reference |
|---|---|---|---|
| `tests/unit/test_fbar_mean_dilatation.py` | E04 Cook gate (`p=2, 16x16`) & `fbar` unit tests (CSR + `C^(p-2)`) | `18 passed` (Cook gate `21.25s`) | [S40.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S40.md) |
| `tests/unit/test_fbar_pipeline_wiring.py` | F-bar `methods.py` & `StudySpec` User API wiring (Kratos-free) | Unit test suite | [test_fbar_pipeline_wiring.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_pipeline_wiring.py) |
| `tests/kratos_spikes/test_p5_v07_kirsch_precursor.py` | V07 Kratos mixed crash check & body-fitted `fbar` Kirsch | `3 passed in 54.99s` | [S42.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S42.md) |
| `tests/kratos_spikes/test_p5_v08_r0_ssn.py` | V08 R0 rigid sphere SSN + punch/friction/separation patches | `4 passed in 47.32s` | [S31.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S31.md) |
| `tests/kratos_spikes/test_p5_v08_r1_r2_ssn.py` | V08 R1 & R2 small-strain deformable SSN (+ R0 & `test_isolation.py`) | `12 passed in 169.39s` | [S32.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S32.md) |
| `tests/kratos_spikes/test_p5_v08_nlgeom_ssn.py` | V08 R1 & R2 NLGEOM (`TotalLagrangianElement3D4N`) + patches | Verified in S41 | [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md) |
| `tests/kratos_spikes/test_p5_e03_reduced_step.py` | E03 13-solid STEP, trimmed quadrature, 3-step thermal + NLGEOM | `4 tests` | [S33.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S33.md) |
| `tests/kratos_spikes/test_p5_e05_door_lock_precursor.py` | E05 5-solid watertight B-Rep & 3-part elastic SSN (`mu=0, 0.3`) | `3 passed in 3.13s` | [S43.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S43.md) |
