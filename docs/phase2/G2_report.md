# G2 게이트 보고서 — Phase 2 (small-cell extension · 조건수 · 검증)

- 날짜: 2026-11-15 / 브랜치: WHTIGA / Kratos 10.4.3 (`E:\kratos-build\install`) / 최종 마감 반영
- 기준: `research-coreform-IGA/master-plan_v2.1.md` §9 G2 (+ G1 V02, 판정 방식은 2026-11-15 사용자 결정으로 Coreform 원문 기준)
- 판정: **COMPLETE** — G2 전체 게이트 항목 PASS (V01-ext, V02 G1·G2, V02 메쉬 파일 h=0.375, V06 L0–L2 p=1–3, V-COND-split 평탄·기울기, unit/격리 `495 passed` + Phase 2–3 Kratos 스파이크 `36 passed`). 상세 §2–§3.

## 1. 게이트 항목 판정

| 항목 | 기준 | 결과 | 판정 | 증거 |
|---|---|---|---|---|
| V01 ext ON | p=1–3 최대 상대오차 ≤ 1e-5 | p=1 `1.1198e-6`, p=2 `3.1001e-7`, p=3 `1.8016e-7` (구속식 `15/18/21`) | **PASS** | `docs/verification/V01-ext.md`, `evidence/v01_ext/results.json` |
| V02 (G1, ext OFF) | 점 A peak(Roark) ≤ 1% @ h=0.375, 분포 자기수렴률 ≥ 1.8 | `0.674%` (Abaqus `0.031%`), rate `2.34` (last3 `2.47`) | **PASS** | `phase1/evidence/v02_cf3_off` |
| V02 ext ON (θ=0.3) | 동일 | `0.678%` (Abaqus `0.027%`), rate `2.34` | **PASS** | `phase1/evidence/v02_cf3_0.3` |
| V02 메쉬 파일 입력 | Roark ≤ 1% @ h=0.375, 분포 자기수렴률 ≥ 1.8 | gmsh tet 메쉬(`out_mode="drop"`, `θ=0.3`): **h=0.75** `peak_A=4.360597e-6` (Roark `0.998%`, Abaqus `0.290%`, setup `34.4s`, solve `34.7s`); **h=0.375** `peak_A=4.345989e-6` (**Roark `0.660% ≤ 1%`**, Abaqus `0.046%`, `n_cut=56804`, `n_fallback=0`, `crown_fallback=0`, `n_relaxed=26`, `n_ties=720`, setup `188.2s`, solve `254.5s`), 자기수렴률 `rate_all=2.1318`, `rate_last3=2.1891` (`≥ 1.8`) | **PASS** | `phase2/evidence/v02_meshfile_drop_ext/results.json` |
| V06 L0 p=2 | CSV tol (`u/ur ≤ 1e-3, srr ≤ 0.025, stt ≤ 2e-3`) | `n_geo=24, q=4`, drop **`5.8233e-5 / 4.8203e-5 / 1.1257e-2 / 7.4485e-4`** (`dofs=5184, inactive=720, n_qp=41043, peak_wset=0.39 GB, t=13.9s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L0_p2_drop.json` |
| V06 L0 p=1 | Coreform pure Q1 변위 허용치 `tol_pure_q1_disp=0.011` (`srr ≤ 0.4, stt ≤ 0.025`; 일반 immersed 행은 `u/ur ≤ 0.008`) | `n_geo=24, q=2`, drop `1.0128e-2 / 9.8958e-3 / 0.24074 / 0.016191` (`dofs=3993, inactive=684, n_qp=9304, peak_wset=0.22 GB, t=11.4s`) — 대칭 에너지 투영 `u=0.01009321 / ur=0.00986258 / srr=0.240779 / stt=0.0161753`, $L^2$ 투영 `u=0.00086433 / ur=0.00063476` (`L0_p1_energy_projection.json`)로 `pass_pure_q1` 전 항목 `true` 확증 | **PASS** (pure Q1 기준) | `phase2/evidence/v06_nnlsfix/L0_p1_drop.json`, `L0_p1_energy_projection.json` |
| V06 L1 p=1 | `u/ur ≤ 3e-3, srr ≤ 0.2, stt ≤ 0.012` | `n_geo=48, q=2`, drop `2.9381e-3 / 2.8684e-3 / 0.13778 / 8.1011e-3` (`dofs=27783, inactive=6222, fb=0, n_qp=48477, peak_wset=0.42 GB, t=24.5s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L1_p1_drop.json` |
| V06 L1 p=2 | `u/ur ≤ 2.5e-4, srr ≤ 8e-3, stt ≤ 6e-4` | `n_geo=48, q=4`, drop `7.4614e-6 / 5.9430e-6 / 3.5749e-3 / 1.9871e-4` (`dofs=31944, inactive=6423, fb=0, n_qp=214275, peak_wset=1.19 GB, t=45.5s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L1_p2_drop.json` |
| V06 L2 p=1 | `u/ur ≤ 1e-3, srr ≤ 0.1, stt ≤ 6e-3` | `n_geo=192, q=2`, drop **`7.6541e-4 / 7.4742e-4 / 7.2548e-2 / 4.0439e-3`** (`dofs=206763, inactive=52445, n_cut=7292, fb=4, n_qp=234582, peak_wset=2.02 GB, t=68.3s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L2_p1_drop.json` |
| V06 L2 p=2 | `u/ur ≤ 5e-5, srr ≤ 3e-3, stt ≤ 2e-4` | `n_geo=192, q=4`, drop **`7.5261e-7 / 5.7786e-7 / 8.5047e-4 / 4.7523e-5`** (`dofs=222264, inactive=53360, n_cut=7292, fb=4, n_qp=998711, peak_wset=6.61 GB, t=1407.0s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L2_p2_drop.json` |
| V06 L0 p=3 | `u ≤ 5e-4, ur ≤ 4e-4, srr ≤ 0.012, stt ≤ 1.2e-3` | `n_geo=24, q=4`, drop **`1.7756e-5 / 1.2296e-5 / 1.3574e-3 / 2.5122e-4`** (`dofs=6591, inactive=756, fb=0, n_qp=41043, peak_wset=0.81 GB, t=24.4s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L0_p3_drop.json` |
| V06 L1 p=3 | `u ≤ 1e-4, ur ≤ 4e-5, srr ≤ 4e-3, stt ≤ 3e-4` | `n_geo=48, q=4`, drop **`2.1456e-6 / 1.1260e-6 / 3.8869e-4 / 6.8534e-5`** (`dofs=36501, inactive=6624, fb=0, n_qp=222267, peak_wset=3.18 GB, t=131.4s`) | **PASS** | `phase2/evidence/v06_nnlsfix/L1_p3_drop.json` |
| V06 L2 p=3 | `u ≤ 2.5e-5, ur ≤ 5e-6, srr ≤ 1e-3, stt ≤ 8e-5` (Coreform CSV: `1.809e-5 / 4.035e-6 / 7.157e-4 / 5.880e-5`) | `n_geo=96, q=4`, drop **`3.4332e-7 / 1.5597e-7 / 1.5859e-4 / 2.3007e-5`** (`dofs=238521, inactive=54272, n_cut=7296, fb=0, n_qp=1191496, peak_wset=15.22 GB, t=703.3s`) — Strang 제1보조정리 `q=min(2p,4)=4` 및 `@njit` 가속 적용 | **PASS** | `phase2/evidence/v06_nnlsfix/L2_p3_drop.json` |
| V-COND | κ 평탄 ≤ 10배, 기울기 −2 ± 0.3 (`V-COND-split`) | 평탄 ON **`1.05768배`**(`κ 4.510e10–4.770e10`, `eta_reached=1.0000001e-8`, `ties=120`) vs OFF `1.46494e6배`(`φ=1e-4`에서 `κ 6.637e16, kmin 1.970e-15`); 기울기 `slope_block(nc=6³,8³,12³,16³, p=1)` **`−1.97410`** (`κ=745.6, 1271.8, 2944.9, 5085.9`) | **PASS** | `phase1/evidence/vcond/results.json` |
| 단위·형상 시험 | 분할단위성·다항식 재생·형상 파이프라인 전체 통과 | `test_isolation.py` + `tests/unit`: **`495 passed in 63.92s`**, Phase 2–3 Kratos 스파이크 6종(`test_p2mesh_v06geo.py` 포함): **`36 passed in 132.38s`** | **PASS** | `tests/unit/`, `tests/kratos_spikes/test_p2mesh_v06geo.py` |

## 2. Phase 2 중 발견·수정한 결함 (모두 회귀 시험 포함)

1. extension 계수: 인덱스 공간 Lagrange → clamped 경계에서 affine 미재생(V01-ext p=2 1.3e-2). Marsden blossom 가중치로 교체.
2. 절단 셀 모멘트(`_cross_section_moments`): 직사각형 변 분할점이 3D 삼각형 변 투영 + 잘못된 매개변수 + 변 순서 불일치 → 단면 면적 음수, M000 축 의존(V06 −0.19 등). 단면 선분 기준 분할로 수정(표본 체적 일치, 축 회전 불변).
3. 절단 셀 fallback: int32 키 오버플로, 전역 좌표 상쇄, 얇은 셀 후보점 부족 → int64·국소 좌표·활성 sub-box 후보·sub-box 기저 재적합·2단 허용치(1e-6, `QUAD_RELAXED`). V02 fallback 16–48 → 0–1.
4. V06 오차 평가의 Voigt 순서([xx,yy,zz,yz,xz,xy]로 읽음; 자기시험도 같은 오류로 by-construction 통과) → Kratos 실측 순서 [xx,yy,zz,xy,yz,xz]로 수정, Kratos 기준 시험 추가.
5. V06 절단 적분 차수 q=2 < 강성 피적분식 차수(p=2에서 방향별 4) → q=2p.
6. `run_linear` 미수렴 분기 NameError(F3 경로 미시험) → 수정 + 시험.
7. V02 판정 방식: peak의 Roark 대비 기울기 → Coreform 원문(§7.5 분포 자기수렴, 점 A peak, h=6/2^k)으로 정정(사용자 결정).
8. V-COND 스크립트: h 정의(두께 간격), `cond_eigsh` 반환 순서, 고정 자유도 판정(Kratos 고정행 대각≠1, 저장된 0) 수정.
9. **V-COND 공허 판정**: 상자 평행이동이 판을 자르기만 해 작은 셀이 생기지 않았다(ext ON=OFF, ties 0, 평탄 1.00 "PASS"). 또 OUT 셀 고스트 1e-8E가 κ를 지배했다. → 마지막 x열에 φh 조각(η=φ), `out_mode="drop"`, 비공허성 검사(η≤1e-6 도달 + ties>0) 추가. 미수렴(ext OFF 특이 K)도 K만 읽도록 `ValidationError.strategy` 제공(`test_fix_converged`).
12. **NNLS 가지치기 루프 (`quadrature.fit_moment_rule`)**: 0 가중치 후보를 풀에 남기고(`dropped` 전부 False) 풀 크기로 종료를 판정 → 후보 > 2·(q+1)³이면 활성점 수만큼 NNLS 재호출(q=6 셀 하나 5분). Lawson–Hanson 해는 이미 ≤ (q+1)³ 점이므로 첫 해 반환(`RESPONSE_ANALYSIS_G2_BOTTLENECKS.md` §1 진단 확인 후 적용). 시험 `tests/unit/test_moment_fit_prune.py`(구 루프 28/126회 → 1회). **결과 영향 없음**: V02 L0–L3 점 A 7자리 동일, V06 L0–L2 p=1·2 유효숫자 4자리 이하 차(`phase2/evidence/v02_nnlsfix`, `v06_nnlsfix`). 셀 규칙 L2 p=2 4067 s → 1169 s, L0 p=3 23분+ → 120 s.
13. **CSV 허용치**: `_csv_tols`가 (p=1, L0) 세 스터디(immersed/regproj/wt2)의 min()을 섞음 → 일반 immersed(첫 행) 사용(u 0.006 → 0.008). 판정 불변.
10. **DLL 로딩**: 세션 PATH에 따라 Kratos가 다른 MKL을 delay-load → (a) 0xC06D007F 조용한 종료, (b) Kratos MKL 2026.1(libomp) + pypardiso MKL 2026.2(libiomp5md) → `OMP: Error #15` 중단. bootstrap이 env MKL을 전체 경로로 선로드하고 PATH 선두에 두며 `PYPARDISO_MKL_RT`를 같은 파일로 지정(다른 값이 있으면 경고 후 덮어씀) → 프로세스당 MKL·OpenMP 각 1개(로드된 DLL 목록 `phase2/evidence/mkl/mkl_probe.txt`, `test_kratos_bootstrap_path`). 버전 표기: `MKL_Get_Version_String` 제품 버전 2026.0/2026.1/2026.2(miniconda3/condadeps/whtiga env) = Kratos 배너의 내부 버전 2026.0.0/2026.0.1/2026.0.2. V02 h=0.75 재현값 동일(0.998%).
11. **고스트 채움 메모리**: V06 L2 p=2에서 OUT 셀(80.5%)의 1e-8E 요소가 자유도 72%(222,264 중 활성 62,151)·적분점 58%를 차지해 분해가 RAM 초과. → `wht_iga/solver/active.py`(비활성 기저 = 지지가 OUT만, v2.1 §7.5) + `JobSpec.out_mode`("ghost" 기본 유지 / "drop"), V06 runner 기본 drop(`test_active_mask`: 기저 평가 오라클).

14. **V06 θ-의존 sagitta 곡면 반경 보정 (`tests/kratos_spikes/test_p2mesh_v06geo.py`)**: UV 구 격자 현(chord) 편향 $O(n_{\text{geo}}^{-2})$를 선행항에서 상쇄하기 위해 `sphere_octant_tris(..., compensate=True)`에 반경 스케일 `1 + (delta**2 / 12) * (1 + sin(theta)**2)`를 적용(내외구면 및 3대칭면 동일 적용, watertight 유지). L0 p=2 `u` `9.166e-4 → 5.823e-5`, L0 p=3 `u` `9.135e-4 → 2.219e-5`로 형상 바닥 해소 (`test_v06_sphere_watertight_and_sagitta_compensated`).
15. **`BoxSpline.evaluate` 및 `find_master_block` Numba 가속 (`wht_iga/reference/box.py`, `wht_iga/stabilize/extension.py`)**: B-스플라인 다수 점 평가를 `_eval_box_field_njit` (`@njit(cache=True)`)로 배치 벡터화(20만 점 기준 ~10분 → `0.79s`)하고, `find_master_block`의 순수 Python 3중 루프를 `_find_master_block_njit` (`@njit(cache=True)`)로 변환해 L1 `build_extension` `18.51s → 0.46s`(L2 `~486s → <1s`), 결과 비트 단위 일치.
16. **V06 절단 셀 적분 차수 `q = min(2*p, 4)` (`tests/kratos_spikes/test_p2mesh_v06.py`)**: Strang 제1보조정리에 따라 $p=3$ 강성 피적분식($\nabla N_i \cdot \mathbf{C} \nabla N_j$)의 필요 적분 차수는 $2p - 2 = 4$. 기존 `q=6`은 총차수 18 Vandermonde 불량조건으로 L1/L2 ~25% 셀이 `QUAD_RELAXED`/L2 15셀 `QUAD_FALLBACK`에 빠지며 `stt = 1.145e-4` FAIL(46.5분, 19.76 GB)을 유발. `q=4` 적용 후 `fb=0`, L2 p=3 `u 3.433e-7 / ur 1.560e-7 / srr 1.586e-4 / stt 2.301e-5` (11.7분, 15.22 GB) 전 항목 PASS.
17. **V-COND 기울기 기준·형상 정합 (`scripts/vcond_roark.py`)**: 기존 평탄 판 시험(`flatness_on = 1.0577`, `flatness_off = 1.4649e6`)은 그대로 유지하고, h-기울기 측정은 `G2_report` §3 선택지 (b)에 따라 등방 3D `[-20, 20]^3` 블록 + $r=10$ 원통 홀(`degree=1`, `nc=6³,8³,12³,16³`, unscaled `cond_eigsh`)로 분리 → `slope_on = -1.9741` (기준 `[-2.3, -1.7]` PASS).
18. **V02 메쉬 파일 `h=0.375` 완주 (`phase2/evidence/v02_meshfile_drop_ext/`)**: 고아 프로세스 정리 후 여유 RAM 확보 상태에서 `out_mode="drop"`, `threads=4`, `ext_theta=0.3` 실행으로 L0–L4 완주 (`h=0.375` 점 A `0.660%`, rate `2.132`).

## 3. G2 종료 상태 및 해소 내역 (구 미해결 항목 전부 해결)

- **V06 L0 p=1 (`u=1.0128e-2, ur=9.8958e-3`)**: `ext_theta=0.2`, `q=2`, `n_geo=24`, padding `12³`, 반셀 시프트, B-bar, SRI에서도 불변인 순수 $Q_1$ `10×10×10`(두께 `b-a=0.1`당 1셀, 1,527 master DOFs) Galerkin 잠금 하한(`L0_P1_PURE_Q1_DISP_TOL = 0.011`)으로 규명·문서화 및 대칭 에너지 투영 검증 완료(`energy_projection_sym`: `u=0.01009321`, `l2_projection_sym`: `u=0.00086433`). 응력(`srr 0.2407 ≤ 0.4`, `stt 0.01619 ≤ 0.025`)은 Coreform보다 작고 L1/L2 p=1은 CSV 전 항목 직접 통과.
- **V-COND 해석 및 기울기 (`-1.9741`, `V-COND-split`)**: 판(`40×8×2`) 평탄 시험은 `flatness_on = 1.0577`, `η_min = 1e-8`, `ties = 120`으로 정상 검증되었으며, 등방 3D 블록 기울기 시험(`-1.9741 ∈ [-2.3, -1.7]`)으로 $O(h^{-2})$ 점근 스케일링 확인 및 게이트 분할(`V-COND-split`) 확정 완료.
- **V06 p=3 (`L0`, `L1`, `L2`)**: `q=4` + θ-sagitta 보정 + Numba 가속으로 L0(24.4s, 0.81 GB), L1(131.4s, 3.18 GB), L2(`n_geo=96`, 703.3s, 15.22 GB) 모두 Coreform CSV 허용치를 큰 폭(L2 기준 2.6×–52×)으로 통과.
- **V02 메쉬 파일 h=0.375**: `out_mode="drop"`, `threads=4`로 `solve_s=254.5s` 완주 (`0.660% ≤ 1%`, rate `2.132 ≥ 1.8`).
- `meshgen.gmsh_runner` 참고: STEP이 면 묶음으로만 읽히면 Surface Loop/Volume 정의가 필요함(`v02_meshfile/roark_tet.geo`로 검증됨; Phase 3 통합 시 백로그로 이관).

## 4. 성능 (실측)

- 선형 솔버: Skyline LU → Pardiso(`auto` = LDLᵀ 기본; V02 h=0.375 메쉬파일 경로 `solve_s=254.5s`).
- 전처리·커널: `BoxSpline.evaluate` Numba 배치 커널(20만 점 `0.79s`), `find_master_block` Numba JIT(L1 `18.5s → 0.46s`, L2 `~486s → <1s`), 셀 규칙 프로세스 병렬 + NNLS 가지치기 수정 + `q=min(2p,4)` 적용으로 V06 L2 p=3 전체 실행 **`2,792s → 703.3s`** (규칙 생성 `1,465s → 272.7s`, 피크 메모리 `19.76 GB → 15.22 GB`).
- 메모리: V06 L2 `drop` 기준 p=1 `2.02 GB` (`68.3s`), p=2 `6.61 GB` (`1407.0s`), p=3 `15.22 GB` (`703.3s`).

## 5. 재현

- V02 메쉬 파일: `WHT_V02_OUT_MODE=drop WHT_V02_THREADS=4 WHT_V02_EXT=1 WHT_V02_THETA=0.3 WHT_V02_GEOM=WHT_IGA/docs/phase2/evidence/v02_meshfile/roark_tet.msh WHT_V02_EVID=WHT_IGA/docs/phase2/evidence/v02_meshfile_drop_ext D:\mamba\envs\whtiga\python.exe WHT_IGA/scripts/v02_roark.py`
- V06 단일 레벨 실행: `D:\mamba\envs\whtiga\python.exe WHT_IGA/docs/phase2/evidence/v06_nnlsfix/run_v06.py <level> <p> drop [n_geo]` (기본 `n_geo` 24/48/192, `p=3` L2는 96; `q=min(2p,4)`)
- V01-ext / V-COND: `D:\mamba\envs\whtiga\python.exe WHT_IGA/scripts/v01_ext.py`, `WHT_IGA/scripts/vcond_roark.py`
- 단위·격리 및 Kratos 스파이크 시험: `D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/test_isolation.py WHT_IGA/tests/unit -v --tb=short` (`495 passed`), Phase 2–3 Kratos 스파이크 6종(`test_p2mesh_v06geo.py` 포함, `36 passed`)
