# Phase 6 S45 — Industrial 3D CAD Trimmed IGA + Extension + AMGCL Scaling

- **Date / Environment**: 2026-11-19, `WHT_IGA` (`D:\mamba\envs\whtiga\python.exe`, Python 3.12, Kratos 10.4.3 at `E:\kratos-build\install`, Intel MKL 2026.0.2).
- **Test Suite**: [`tests/kratos_spikes/test_p6_industrial_bracket_amgcl.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_industrial_bracket_amgcl.py)
- **Verification Command**:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests/kratos_spikes/test_p6_industrial_bracket_amgcl.py
  ```
- **Outcome**: **`3 passed in 52.51s`**, overall Peak RSS **`970.4 MB` (`0.970 GB`)** (`< 1.5 GB` budget while Track 1 ran in parallel).

---

## 1. Industrial CAD Inspection (`wht_iga/resources/`)

All three STEP files were loaded with `load_cad_file(..., heal=True)` and tessellated via `cad.tessellate(...)` (`angular=0.3`):

| Model | Solids | Faces | Free / Bnd Edges | Watertight | B-Rep Span `(dx, dy, dz)` | B-Rep Volume | Tess Triangles | Tess Vol Error |
|---|---:|---:|---:|:---:|---|---:|---:|---:|
| `Generator_bracket.step` | 1 | 251 | `0 / 0` | Yes | `[64.204, 95.002, 114.499]` | `104,578.98` | `11,860` | `0.0323%` |
| `Flange.STEP` | 12 | 194 | `0 / 0` | Yes | `[152.400, 132.452, 76.200]` | `385,905.69` | `7,842` | `0.4704%` |
| `Support_plate.STEP` | 1 | 81 | `0 / 0` | Yes | `[1000.000, 105.000, 1000.000]` | `10,028,443.93` | `612` | `0.0929%` |

Inspection notes:
- `Generator_bracket.step`: Single solid cast/machined bracket oriented obliquely to global axes (volume fill ratio `~15.0%` of its global AABB), making it an ideal testbed for PCA OBB alignment.
- `Flange.STEP`: 12-solid bolted flange assembly (2 main flange halves `solids[0]` / `solids[1]`, 1 central ring `solids[2]`, 9 fasteners `solids[3..11]`). `solids[0]` (`vol = 171,489.33`, `2,254` tris at `0.003 * diag`) provides a clean single-solid curved industrial flange benchmark.
- `Support_plate.STEP`: Single solid ribbed plate (`1000 × 105 × 1000`) with `5 mm` bottom plate and `5 mm` stiffener ribs (`1:200` thickness-to-span aspect ratio). On coarse uniform grids (`ny = 3..12`), thin-rib / L-junction cut cells trigger `octree_fallback` (`32–185` fallback cells), inflating per-element quadrature point counts and Kratos element assembly time unless `ny` resolves the `5 mm` plate thickness or THB local refinement is used.

---

## 2. `Generator_bracket.step`: Axis-Aligned vs PCA OBB (`compute_pca_obb`)

Configuration: `p=2`, `q=2`, `theta=0.2` extension, `E=210,000`, `nu=0.3`, `cache_dir="off"`, `workers=8`, `0.003 * diag` tessellation (`11,860` tris, `tess = 0.20s`).

### 2.1 Grid & Trimmed Classification Metrics

| Frame | Grid | Padding | Total Cells | Active Cells | Active-Cell Ratio | `n_fallback` | Quad Vol Error | `n_master / n_extended / n_inactive` | Active DOFs | `trim` (s) |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| Axis-Aligned (AA) | `(8, 12, 14)` | `0.03` | 1,344 | 386 | **`0.2872`** | **`0`** | `0.0323%` | `852 / 189 / 1,199` | 3,123 | `2.33 s` |
| PCA OBB (`compute_pca_obb`) | `(11, 7, 4)` | `0.05` | **308** | 141 | **`0.4578`** | **`0`** | `0.0323%` | `295 / 171 / 236` | 1,398 | **`0.93 s`** |

Key OBB gains:
- **Active-cell ratio**: `0.2872 -> 0.4578` (**`+59.4%`** relative improvement; **`+17.06` percentage points**).
- **Background cell & DOF reduction**: `4.36×` fewer total background cells (`308` vs `1,344`), `5.08×` fewer inactive control points (`236` vs `1,199`), `2.23×` fewer active DOFs (`1,398` vs `3,123`), and **`2.51×` faster** `build_trimmed_table` (`0.93 s` vs `2.33 s`) with `n_fallback == 0` and `0.0323%` volume error on both frames.

### 2.2 Kratos `PardisoLU` vs `AMGCL` Comparison (`Generator_bracket.step`)

| Case (`n_eq`, `nnz`, `asm`) | Solver | Strategy (s) | Solve Net (s) | Diff vs `PardisoLU` | `\|sum(R)\|` | RSS Before → After (Delta) |
|---|---|---:|---:|---:|---:|---:|
| `GB_AA` (`3,123`, `1,313,919`, `0.819s`) | `PardisoLU` | `1.405 s` | `0.587 s` | `0.0` | `1.840e-11` | `451.6 → 611.5 MB` (`+159.9 MB`) |
| `GB_AA` (`3,123`, `1,313,919`, `0.819s`) | `AMGCL_CG_ILU0` | `1.270 s` | `0.451 s` | `2.156e-09` | `2.869e-11` | `481.2 → 525.5 MB` (**`+44.3 MB`**) |
| `GB_AA` (`3,123`, `1,313,919`, `0.819s`) | `AMGCL_GMRES_ILU0` | **`1.249 s`** | **`0.430 s`** | `2.508e-09` | `1.194e-11` | `488.6 → 529.8 MB` (**`+41.2 MB`**) |
| `GB_OBB` (`1,398`, `596,700`, `0.325s`) | `PardisoLU` | `0.507 s` | `0.181 s` | `0.0` | `4.875e-12` | `492.6 → 564.4 MB` (`+71.8 MB`) |
| `GB_OBB` (`1,398`, `596,700`, `0.325s`) | `AMGCL_CG_ILU0` | **`0.406 s`** | **`0.081 s`** | `3.067e-09` | `1.979e-12` | `497.2 → 524.2 MB` (**`+27.0 MB`**) |
| `GB_OBB` (`1,398`, `596,700`, `0.325s`) | `AMGCL_GMRES_ILU0` | `0.442 s` | `0.116 s` | `3.398e-09` | `2.565e-12` | `513.7 → 527.5 MB` (**`+13.7 MB`**) |

*(In `GB_OBB`, displacements and reactions were solved in the local PCA OBB frame and mapped back to global coordinates via `ObbFrame.to_global_vectors`.)*

---

## 3. `Flange.STEP` (Solid 0): Medium-Scale Trimmed `p=2` + Extension AMGCL Scaling

Configuration: `Flange.STEP` `solids[0]` (`2,254` tris at `0.003 * diag`, `tess = 0.03s`), `p=2`, `q=2`, `theta=0.2` extension, `E=210,000`, `nu=0.3`.

### 3.1 Stage Breakdown

| Case | Grid (`pad`) | Active / Total Cells (Ratio) | `cut` / `fb` | Vol Error | `M / E / I` | Active DOFs (`nnz`) | `trim` (s) | `cls + ext` (s) | `vol + elem` (s) | `ties` (s) | `asm` (s) |
|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| `Flange_5k` | `(22, 10, 11)` (`0.05`) | `737 / 2,420` (`0.3045`) | `641 / 0` | `0.4930%` | `1,530 / 223 / 1,991` | `5,259` (`1,998,405`) | `0.577` | `0.1249` | `0.244` | `0.296` | `1.081` |
| `Flange_8k5` | `(27, 12, 13)` (`0.04`) | `1,371 / 4,212` (`0.3255`) | `1,046 / 2` | `0.4929%` | `2,428 / 407 / 3,255` | `8,505` (`3,626,127`) | `1.314` | `0.2711` | `0.517` | `0.544` | `1.296` |

### 3.2 `PardisoLU` vs `AMGCL_CG_ILU0` vs `AMGCL_GMRES_ILU0` Runtime & RSS

| Case | Solver | Strategy (s) | Solve Net (s) | Diff vs `PardisoLU` | `\|sum(R)\|` | RSS Before → After (Delta) | Peak RSS (MB) |
|---|---|---:|---:|---:|---:|---:|---:|
| `Flange_5k` | `PardisoLU` | `1.606 s` | `0.526 s` | `0.0` | `1.690e-11` | `517.4 → 748.9 MB` (`+231.4 MB`) | `748.9 MB` |
| `Flange_5k` | `AMGCL_CG_ILU0` | `1.604 s` | `0.523 s` | `2.438e-10` | `6.808e-11` | `523.6 → 595.9 MB` (**`+72.3 MB`**) | `595.9 MB` |
| `Flange_5k` | `AMGCL_GMRES_ILU0` | **`1.576 s`** | **`0.495 s`** | `5.746e-10` | `1.683e-11` | `524.9 → 600.3 MB` (**`+75.4 MB`**) | `600.3 MB` |
| `Flange_8k5` | `PardisoLU` | `1.904 s` | `0.609 s` | `0.0` | `4.077e-11` | `534.4 → 970.4 MB` (`+436.0 MB`) | **`970.4 MB`** |
| `Flange_8k5` | `AMGCL_CG_ILU0` | `1.898 s` | `0.602 s` | `1.412e-09` | `3.923e-11` | `547.1 → 664.9 MB` (**`+117.8 MB`**) | **`664.9 MB`** |
| `Flange_8k5` | `AMGCL_GMRES_ILU0` | **`1.855 s`** | **`0.559 s`** | `1.463e-09` | `3.182e-11` | `548.9 → 674.5 MB` (**`+125.6 MB`**) | **`674.5 MB`** |

---

## 4. Phase 6 Takeaways

1. **AMGCL Memory & Speed Advantage on Real Trimmed+Extension CAD**:
   - Across all four industrial CAD systems (`1,398` to `8,505` active DOFs, `p=2`, `q=2`, `theta=0.2` extension), both `AMGCL_CG_ILU0` and `AMGCL_GMRES_ILU0` (`use_block_matrices_if_possible=False`, `tol=1e-10`) converged cleanly to **`2.4e-10 – 3.4e-9`** relative L2 agreement vs `PardisoLU` (well below the `1e-7` gate) and preserved global reaction equilibrium **`<= 6.8e-11`**.
   - At `8,505` active DOFs (`3.63M` nonzeros), `AMGCL_CG_ILU0` used **`3.70×` less solver RSS growth** (`+117.8 MB` vs `+436.0 MB` for `PardisoLU`), cutting Peak RSS by **`305.5 MB`** (`664.9 MB` vs `970.4 MB`) while matching or beating `PardisoLU` wall time (`1.855–1.898 s` vs `1.904 s`).
2. **PCA OBB Efficiency on Oblique CAD (`Generator_bracket.step`)**:
   - `compute_pca_obb` raised active-cell density from **`28.72%`** to **`45.78%`** (`+59.4%` gain) and cut total background cells `4.36×` (`308` vs `1,344`) with `n_fallback == 0` and `0.0323%` volume error.
3. **Kratos Naming & Thin-Feature Caveats Observed**:
   - **Kratos `ModelPart` naming**: Dots (`.`) in a root `ModelPart` name (e.g. `Flange_8.5k`) are parsed by Kratos as SubModelPart path separators, causing `SetBufferSize(2)` to raise `RuntimeError`. Sanitize model part names with underscores.
   - **High-aspect thin-walled CAD (`Support_plate.STEP`, `1:200` plate/rib thickness)**: Uniform coarse grids (`ny=4..12`) generate `octree_fallback` cells (`depth=3`, `>10k` QPs/cell) in `5 mm` L-junction ribs; thin-plate structures require `ny >= 20` or THB local refinement around thin ribs rather than uniform coarse background grids.
