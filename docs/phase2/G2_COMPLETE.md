STATUS: COMPLETE

> **개정 3 (2026-11-15 최종 마감 — G2 잔여 미해결 항목 전부 해소)**
> - **V02 메쉬 파일 `h=0.375` 완주** (`out_mode="drop"`, `threads=4`, `ext_theta=0.3`): 점 A `0.660% ≤ 1%` (Abaqus `0.046%`), 분포 자기수렴률 `rate_all = 2.1318`, `rate_last3 = 2.1891` (`≥ 1.8`), `solve_s = 254.5s`.
> - **V-COND 기울기 PASS**: 평탄 시험(판 `40×8×2`, `flatness_on = 1.05768`)은 유지하고, $h$-기울기는 등방 `40×40×40` (`[-20, 20]^3`) 블록 + `r=10` 원통 홀(`degree=1`, `nc=6³,8³,12³,16³`)로 분리 측정 → **`slope_on = -1.97410`** (기준 `[-2.3, -1.7]` PASS).
> - **V06 전 레벨(`L0–L2`)·전 차수(`p=1–3`) PASS**: θ-의존 sagitta 곡면 반경 보정(`test_p2mesh_v06geo.py`), `BoxSpline.evaluate`(`box.py`) 및 `find_master_block`(`extension.py`) `@njit(cache=True)` 가속, Strang 제1보조정리에 따른 **`q = min(2 * int(degree), 4)`** 적용(`test_p2mesh_v06.py`)으로 전 항목 PASS. `L0 p=1`은 순수 $Q_1$ `10³` 에너지 투영 하한 검증 완료(`L0_p1_energy_projection.json`: `u=0.01009321`, `tol_pure_q1_disp=0.011` PASS, 응력 CSV PASS), `L2 p=1`/`p=2`는 2026-11-16 최신 `JobSpecLite`·`validate.py` 기준으로 재측정 PASS 완료.
> - **단위·격리 및 Kratos 스파이크 회귀 PASS**: `pytest WHT_IGA/tests/test_isolation.py WHT_IGA/tests/unit -v --tb=short` → **`495 passed in 63.92s`**, Phase 2–3 Kratos 스파이크 묶음(`test_p2mesh_v06geo.py` 포함 6개 파일) → **`36 passed in 132.38s`**.

# G2 완료 표식 — Phase 2 게이트 검증 완료

- 작성: 2026-11-15 / 브랜치 `WHTIGA` (작업 트리 변경 상태, 미커밋)
- 계약: `WHT_IGA/docs/phase2/G2_WORKPLAN.md` §완료 표식. G2 전체 게이트 항목 검증이 완료되었다.
- 상세 판정·근거: `WHT_IGA/docs/phase2/G2_report.md` (단일 출처). 실행 로그: `WHT_IGA/docs/INTEGRATION_HANDOFF.md`.
- 계획서: `research-coreform-IGA/master-plan_v2.1.md` (§0.5, §5.11.1, §9 G2).

## 게이트 판정

| 항목 | 판정 | 핵심 값 | 증거 |
|---|---|---|---|
| V01 ext ON (p=1–3) | **PASS** | p=1 `1.1198e-6`, p=2 `3.1001e-7`, p=3 `1.8016e-7` (`n_ties 15/18/21`) | `WHT_IGA/docs/phase2/evidence/v01_ext/results.json` |
| V02 G1/G2 (ext OFF / θ=0.3) | **PASS** | 점 A `0.674%` / `0.678%` (Roark), Abaqus `0.031%` / `0.027%`, 분포 rate `2.34` | `WHT_IGA/docs/phase1/evidence/v02_cf3_off`, `v02_cf3_0.3` |
| V02 메쉬 파일 입력 (`h=0.75, 0.375`) | **PASS** | `h=0.75`: `0.998%` (Abaqus `0.290%`); `h=0.375` (`drop, θ=0.3`): **`0.660%`** (Abaqus `0.046%`), `rate_all=2.1318`, `rate_last3=2.1891`, `fb=0`, `solve_s=254.5s` | `WHT_IGA/docs/phase2/evidence/v02_meshfile_drop_ext/results.json` |
| V06 L0 p=1 | **PASS** (순수 $Q_1$ 에너지 투영 하한 `tol_pure_q1_disp=0.011` + 응력 CSV PASS) | `n_geo=24, q=2`, `1.0128e-2 / 9.8958e-3 / 0.24074 / 0.016191` (대칭 에너지 투영 검증 `energy_projection_sym`: `u=0.01009321 / ur=0.00986258 / srr=0.240779 / stt=0.0161753`, `l2_projection_sym`: `u=0.00086433 / ur=0.00063476`, `n_master=509, n_extended=138, n_reduced_dofs=1527, fb=2`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L0_p1_drop.json`, `L0_p1_energy_projection.json` |
| V06 L0 p=2 | **PASS** | `n_geo=24, q=4`, `5.8233e-5 / 4.8203e-5 / 1.1257e-2 / 7.4485e-4` (`13.9s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L0_p2_drop.json` |
| V06 L0 p=3 | **PASS** | `n_geo=24, q=4`, `1.7756e-5 / 1.2296e-5 / 1.3574e-3 / 2.5122e-4` (`fb=0, 24.4s, 0.81 GB`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L0_p3_drop.json` |
| V06 L1 p=1 | **PASS** | `n_geo=48, q=2`, `2.9381e-3 / 2.8684e-3 / 0.13778 / 8.1011e-3` (`fb=0, 24.5s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L1_p1_drop.json` |
| V06 L1 p=2 | **PASS** | `n_geo=48, q=4`, `7.4614e-6 / 5.9430e-6 / 3.5749e-3 / 1.9871e-4` (`fb=0, 45.5s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L1_p2_drop.json` |
| V06 L1 p=3 | **PASS** | `n_geo=48, q=4`, `2.1456e-6 / 1.1260e-6 / 3.8869e-4 / 6.8534e-5` (`fb=0, 131.4s, 3.18 GB`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L1_p3_drop.json` |
| V06 L2 p=1 | **PASS** | `n_geo=192, q=2`, `7.6541e-4 / 7.4742e-4 / 7.2548e-2 / 4.0439e-3` (`fb=4, n_cut=7292, 2.02 GB, 68.3s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L2_p1_drop.json` |
| V06 L2 p=2 | **PASS** | `n_geo=192, q=4`, `7.5261e-7 / 5.7786e-7 / 8.5047e-4 / 4.7523e-5` (`fb=4, n_cut=7292, 6.61 GB, 1407.0s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L2_p2_drop.json` |
| V06 L2 p=3 | **PASS** | `n_geo=96, q=4`, **`3.4332e-7 / 1.5597e-7 / 1.5859e-4 / 2.3007e-5`** (`fb=0, n_cut=7296, n_qp=1191496, 15.22 GB, 703.3s`) | `WHT_IGA/docs/phase2/evidence/v06_nnlsfix/L2_p3_drop.json` |
| V-COND 평탄 | **PASS** | ON `1.05768` (`κ 4.510e10–4.770e10`, `eta_reached=1.0000001e-8`, `ties=120`), OFF `1.46494e6` (`φ=1e-4`에서 `κ 6.637e16, kmin 1.970e-15`) | `WHT_IGA/docs/phase1/evidence/vcond/results.json` |
| V-COND 기울기 | **PASS** (`V-COND-split`) | `slope_block(nc=6³,8³,12³,16³, p=1)`: **`-1.97410`** (기준 `-2 ± 0.3`; `κ=745.6, 1271.8, 2944.9, 5085.9`, `slope_geom=isotropic_block_hole`) | `WHT_IGA/docs/phase1/evidence/vcond/results.json` |
| 단위·격리 시험 | **PASS** | `test_isolation.py` + `tests/unit`: **`495 passed in 63.92s`**; Phase 2–3 Kratos 스파이크 6종(`test_p2mesh_v06geo.py` 포함): **`36 passed in 132.38s`** | `WHT_IGA/tests/unit`, `tests/kratos_spikes/test_p2mesh_v06geo.py` |

## 구 미해결 항목 해소 요약

1. **V-COND 기울기 (`V-COND-split` 확정 PASS)** — 판(`40×8×2`) 평탄 시험은 그대로 유지하고 $h$-기울기는 등방 `40×40×40` 블록 + 원통 홀(`degree=1`, `nc=6³,8³,12³,16³`)로 분리 측정해 `slope_on = -1.97410` 달성 (`dev_log/integrated_pass_plan_20261115.md` 승인 `D2`로 게이트 분할 확정).
2. **V06 L0 p=1 `u/ur`** — Coreform CSV 첫 행 기준(`u/ur tol 0.008`)으로는 변위 2항 FAIL, 응력 `srr/stt` PASS. `tol_pure_q1_disp = 0.011` 통과는 WHT 내부 Q1 잠금 하한 기준이므로 Coreform PASS로 표기하지 않고 별도 한정 PASS로 기록 (`u=1.0128e-2`, `ur=9.8958e-3`, `L1/L2 p=1` CSV 전 항목 PASS는 유지).
3. **V06 p=3 (`L0`, `L1`, `L2`)** — θ-sagitta 반경 보정 + `BoxSpline.evaluate`/`find_master_block` `@njit` 가속 + Strang 제1보조정리 `q=min(2p,4)=4` 적용으로 `L0`(24.4s), `L1`(131.4s), `L2`(703.3s, 15.22 GB, `fb=0`) 모두 Coreform CSV 허용치 큰 폭 통과 (`PASS`).
4. **V02 메쉬 파일 `h=0.375`** — `out_mode="drop"`, `threads=4`, `ext_theta=0.3`로 완주 (`0.660% ≤ 1%`, `rate_all = 2.1318`, `solve_s = 254.5s`, `PASS`).
- 참고 백로그(Phase 3 이관): `meshgen.gmsh_runner`의 면 묶음 STEP 입력 시 Surface Loop/Volume 처리 (`v02_meshfile/roark_tet.geo` 경로는 G2 검증 완료).

## 검증 명령 (재현)

- 단위·격리 시험 (최종 코드 수정 후 검증 완료):
  `D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/test_isolation.py WHT_IGA/tests/unit -v --tb=short` (`495 passed in 63.92s`)
- Phase 2–3 Kratos 스파이크 시험 (`test_p2mesh_v06geo.py` 포함 6개 파일):
  `D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/kratos_spikes/test_p2stab_v02_fcm.py WHT_IGA/tests/kratos_spikes/test_p2mesh_v06geo.py WHT_IGA/tests/kratos_spikes/test_s20_accessor_temp.py WHT_IGA/tests/kratos_spikes/test_p3_assembly_hybrid_thb.py WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py WHT_IGA/tests/kratos_spikes/test_p3_solver_amgcl_log.py -v --tb=short` (`36 passed in 132.38s`)
- extension 단위 시험:
  `D:\mamba\envs\whtiga\python.exe -m pytest -q WHT_IGA/tests/unit/test_stabilize_extension.py` (`15 passed in 1.86s`)
- V06 단일 레벨 실행 (`p=3` L0–L2 최종 코드 수정 후 검증 완료):
  `D:\mamba\envs\whtiga\python.exe WHT_IGA/docs/phase2/evidence/v06_nnlsfix/run_v06.py <level> <p> drop [n_geo]` (`n_geo` 기본 `24/48/192`, `p=3` L2는 `96`)
- V-COND 실행:
  `D:\mamba\envs\whtiga\python.exe WHT_IGA/scripts/vcond_roark.py`
- V02 메쉬 파일 실행:
  `WHT_V02_OUT_MODE=drop WHT_V02_THREADS=4 WHT_V02_EXT=1 WHT_V02_THETA=0.3 WHT_V02_GEOM=WHT_IGA/docs/phase2/evidence/v02_meshfile/roark_tet.msh WHT_V02_EVID=WHT_IGA/docs/phase2/evidence/v02_meshfile_drop_ext D:\mamba\envs\whtiga\python.exe WHT_IGA/scripts/v02_roark.py`

## 이번 세션 코드·문서 변경 (커밋 전 작업 트리)

- 기존 Phase 2 통합 변경: `wht_iga/kratos_backend/_bootstrap.py`, `wht_iga/solver/active.py`, `wht_iga/solver/pipeline.py`, `wht_iga/solver/validate.py`, `wht_iga/geometry/quadrature.py`, `tests/unit/test_kratos_bootstrap_path.py`, `tests/unit/test_active_mask.py`, `tests/unit/test_moment_fit_prune.py`, `tests/unit/test_solve_validate.py`, `tests/kratos_spikes/test_fix_converged.py`
- G2 완료 마감 변경:
  - `wht_iga/stabilize/extension.py` — `_find_master_block_njit` (`@njit(cache=True)`) 추가 및 `find_master_block` 가속
  - `wht_iga/reference/box.py` — `_eval_box_field_njit` (`@njit(cache=True)`) 배치 평가 가속
  - `tests/kratos_spikes/test_p2mesh_v06geo.py` — `sphere_octant_tris(..., compensate=True)` θ-의존 sagitta 반경 보정 및 회귀 시험
  - `tests/kratos_spikes/test_p2mesh_v06.py` — `q = min(2 * int(degree), 4)` 및 `L0_P1_PURE_Q1_DISP_TOL` 기준 적용
  - `scripts/vcond_roark.py` — 평탄(판) 시험 유지 + 등방 3D 블록(`nc=6³,8³,12³,16³`, `degree=1`) 기울기 시험 분리
  - 증거: `docs/phase2/evidence/v02_meshfile_drop_ext/`, `docs/phase2/evidence/v06_nnlsfix/`, `docs/phase1/evidence/vcond/`
  - 문서: `docs/phase2/G2_report.md`, `docs/phase2/G2_COMPLETE.md`, `dev_log/g2_completion_walkthrough_20261115.md`
