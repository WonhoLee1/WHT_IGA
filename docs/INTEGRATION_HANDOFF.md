# 통합 담당 인계 작업 현황 (2026-11-14 23:30~)

Phase 2 에이전트가 넘긴 Kratos 실행 작업. 사용자 지시: 여기(통합 세션)에서 처리.

## 인계 받은 상태 (Phase 2 보고, 검증은 아직 안 함)
- 완료: F1/F2(Voigt, /2), F4(fallback 인접 + 빈 BC 거부), F5(단조 ID), F6(face_id 필수 + `hole_face_by_cylinder_wall`), BC 헬퍼 3곳 통일, AGENTS v2.1, ADR-002 p=3 재검증(1/β), penalty 규칙 p=3 → order 6. 단위 33개, AST 통과.

## 해야 할 일 (Kratos 필요)
| # | 작업 | 명령/조건 | 상태 |
|---|---|---|---|
| 1 | V02 6레벨, extension 없음 | `scripts/v02_roark.py` 기본 → `docs/phase1/evidence/v02/` (로그 `evidence/v02_run.log`) | PID 32824는 사용자 결정으로 중단(Skyline 16 GB). Pardiso 재실행(00:00–00:40): 레벨 0–4 완료 — h=10 10.55%, 5 8.24%, 3 0.622%, 2.5 6.55%, 1.5 0.575% (peak 오차, 비단조: 셀 정렬·crown fallback(16–48셀) 영향, extension OFF). **레벨 5(h=0.75, ~44만 DOF) `pardiso_lu` "Decomposition failed!"** (메모리 vs 특이성 미판별, Eigen 래퍼가 오류 코드 미전달) → `pardiso_ldlt` + 메모리 추적으로 레벨 5만 재실행 중 (`evidence/v02_l5`) |
| 2 | V02 extension on/off + θ(0.1/0.2/0.35)·α sweep | `WHT_V02_EXT=1`, `WHT_V02_THETA`, `WHT_V02_GHOST_ALPHA`. ext 결과는 `evidence/v02_ext`. sweep마다 `WHT_V02_EVID`로 폴더 분리할 것(덮어쓰기 방지) | 대기(1 종료 후) |
| 3 | T6 p=3 재확인 | `WHT_V02_P=3` (q=6). NNLS 발산·fallback 셀 수 감시 | 대기 |
| 4 | F3 Kratos 시험, STAB fcm, mesh Kratos 시험 재실행 | `test_p2stab_*`, `test_p2mesh_*`, F3 관련 | STAB: **6 passed (83.6 s, extension 수정 후)**. mesh·F3 대기 |
| 5 | V06 풀게이트 | `test_p2mesh_v06.py` (V06.md에 레벨0 FAIL 기록, extension 후 재측정). PID 32352 `test_v06_level1` 실행 중(22:39~) | 실행 중 프로세스 결과 확인 |
| 6 | V-COND sweep | `scripts/vcond_roark.py` (모체 형상 = V02 Roark 판, v2.1 확정) | 대기 |
| 7 | V01-ext 스크립트 신규 | `scripts/v01_ext.py` | **PASS p=1,2,3 (수정 후)**. 수정 전 p=2/3 FAIL — extension 연산자가 index 공간 Lagrange 외삽이라 clamped 경계 근처 affine 미재현. 사용자 결정으로 여기서 수정 (`extension.py` Marsden blossom 가중치, `knots=` 필수). `docs/verification/V01-ext.md`, `docs/phase2/P22-extension.md` 정정절 |
| 6a | V-COND 스크립트 정정 | h=면내 크기(기존 두께 간격 사용 버그), h 4단계(20/40/60/80), BC=`roark_left_pinned`, 고정 자유도 제외 κ, `cond_eigsh` 반환 순서 버그(κ 자리에 λmax) 수정, `WHT_VCOND_EVID` | 완료, 미실행 |
| 8 | G2 보고서 수치 기입 | `docs/phase2/G2_report.md` | 마지막 |

## 고속화 (2026-11-14~15, 사용자 결정: PID 32824·32352 중단 후 고속화·재실행)

진단 (py-spy 20 s 표본):
- V02 PID 32824: 레벨 4(133×27×4, ~7만 DOF) `SkylineLUFactorizationSolver` 풀이 중 16 GB (총 31.6 GB). 레벨 5(~44만 DOF)는 메모리 부족 확실. 또 구버전 스크립트로 `v02_bcnew`에 쓰고 있었음(파일명 h1/h0.5/h0.25로 레벨 덮어씀).
- V06 PID 32352: 95% `winding_inside` — octree 대체 경로가 leaf마다(512/셀) 1~8점으로 전체 표면 winding 호출. 또 extension 없이 실행 중(레벨 0 FAIL 원인 그대로 → 판정 무의미).

수정:
- `JobSpec.linear_solver` (기본 `"pardiso_lu"`, `"skyline_lu"` 선택 가능, 미지 이름·MKL 없음은 오류). `pipeline._make_linear_solver`. V02 스크립트 `WHT_V02_SOLVER`.
- `geometry/classify.py`: winding 커널 무할당 스칼라 연산 + 16점 미만은 직렬 커널. 원 커널과 20000점 값 차이 0, 판정 불일치 0.
- `geometry/quadrature.py`: octree leaf 판정·직사각형 변 조각 판정을 셀당 1회 winding 호출로 묶음, 삼각형 AABB 선필터 후 SAT. 결과 점·가중치·순서 동일 (V06 형상 418셀, 대체 12셀 포함, 차이 0).
- 실측 (V06 형상 n_geo=24, h=0.1, 같은 셀): 대체 셀 4개 6.61 s → 0.39 s (17×), 일반 셀 20개 7.62 s → 1.28 s (6×).
- V06 시험: extension θ=0.2 기본 ON (`ext_theta=None`이면 OFF), Pardiso.
- `trimmed.build_trimmed_table(workers=)`: 셀별 적분 규칙을 프로세스 풀로 병렬 (기본 물리 코어 수, `WHT_GEO_WORKERS`, 400셀 미만은 직렬, 작업 프로세스는 numba/MKL/OMP 스레드 1). V02 레벨 2(잘린 셀 1667): 직렬 76.3 s → 7 프로세스 23.6 s, 결과 5종 완전 동일. numba 범위: classify의 SAT·축 투표·winding 커널. quadrature의 모멘트 계산(Python 루프)·NNLS(scipy)·octree 구동은 numba 아님 → 다음 고속화 후보.
- 버그 수정: `pipeline.run_linear` 미수렴 분기가 `ValidationIssue`를 import하지 않아 `ValidationError` 대신 `NameError` (F3 경로 미시험). import 추가 + `test_fix_converged.py::test_nonconverged_raises_validation_error`(1회 CG로 미수렴 강제). 9 passed.
- 절단 적분 캐시: `build_trimmed_table(cache_dir=)`, 키 = 입력 + 코드 소스 SHA-256, `WHT_IGA/.cache/trimmed`, `WHT_GEO_CACHE=off`. 시험 3건 통과.
- Pardiso: Kratos `linear_solver` = pardiso_lu/ldlt/llt/skyline_lu + `num_threads_mkl` (Eigen 래퍼라 iparm은 C++ 고정, follower+ldlt/llt는 오류). Python `wht_iga/linalg/pardiso.py` `PardisoFactor` (spd/sym/general/auto, 병렬 ordering iparm2=3, 상삼각 전달) → reference/solve, compose.solve_reduced, cond_eigsh shift-invert. 단위 13 + 전체 unit 267 passed. 최적 조합은 `scripts/bench_pardiso.py`로 측정 예정.
- `test_trimmed_probe.py` STEP 경로를 `__file__` 기준으로 (cwd 무관). unit 245 passed.

## 주의
- 같은 evidence 폴더에 두 실행이 동시에 쓰지 않게 한다. 한 번에 무거운 실행 1개.
- `v02_bcnew` evidence 폴더는 다른 주체의 것 — 건드리지 않음.
- V07 gate: **최신 계획서(v2.1) 기준 사용 (사용자 결정 2026-11-14)**. v2.1 l.119/1189/1222 — Hard Gate = `slope > 0.5` + 레벨별 `(e_u, e_σxx) ≤ (0.01,0.2),(0.005,0.1),(0.002,0.05)` + 체적 0.1%; `2.3/2.2`는 Target(목표치, 미달해도 FAIL 아님).
- V-COND 모체 형상: v2.1 l.119에서 **V02 Roark 판(구멍 중심 이동)**으로 확정됨 → 6번 작업은 이 형상으로 진행.
- 전체 시험 재실행(S16 61.8 s 재확인 포함)은 무거운 실행이 없을 때.
- 통합 담당 이번 세션 수정 내역: `docs/REVIEW_P0_P1.md` §5.

## RESPONSE_ANALYSIS_V02_CONVERGENCE 검증 (2026-11-15, 통합 담당)

| 항목 | 분석 주장 | 통합 실측 |
|---|---|---|
| 2-A int32 키 오버플로 | 확인 | numpy 1.26.4 `dtype(int)=int32` 재현. int64 + 원점 기준 키로 수정 → cast 경고 0. fallback 수는 불변(16/16) |
| 2-C 전역 좌표 상쇄 | 그럴듯 | `_cell_moments` 셀 국소 좌표로 수정. fallback 수 불변 |
| 2-B 접선 프로브 M000≤0 | 주원인 주장 | **반증**: fallback 셀 M000/|cell| = 0.0907(L1)/0.05–0.064(L3), 전부 양수. 2단 깊이 프로브 추가했으나 fallback 불변 |
| 2-D 후보점 부족 | 부원인 | **실제 주원인**: 잔차가 후보 수에 따라 4e-3→7.9e-6로 줄다 1e-12 미달. 활성 sub-box Gauss 후보(셀 기저 → sub-box 기저 재적합, 체적 일치 검사) + 2단 허용치(`_NNLS_RELAXED=1e-8`, `QUAD_RELAXED` 경고) |
| 결과 | — | fallback L1 16→0, L2 16→0, L3 48→2, L4 32→0 (relaxed 0/8/30/16) |
| 4 레벨5 메모리 | -2 추정 | 행렬 정상(대칭, 양정치 n_neg=0, 빈 행 0). pypardiso 직접: LDLᵀ 성공 factor 11.5 GB(nnz 1.10e9, 섭동 9360), LU 성공 22.8 GB(nnz int32 넘침). Kratos 실패는 프로세스 14 GB + factor로 메모리 한계 추정 |
| 처방 | nz 축소 | 사용자 결정: 레벨 5 = (267,53,2) |
| 기준 문구 | "3레벨 ≤2%" | 오기. G2 = 4레벨+, ≤1%, 수렴률 ≥1.8 |

- 진단 스크립트 `scripts/diag_v02_singular.py` 4-A(phase 11 정보 유실)·4-B(행렬 이중 보유) 수정.
- `test_geo_quad_db.py` 합 비교 rtol 0 → 1e-14 (SQL 순차합 vs numpy 쌍별합, 2.8e-17 차이). unit 322 passed.
- V-COND: Kratos 고정행은 비대각이 저장된 0으로 남고 대각은 스케일값 → `vcond_roark.py`가 `eliminate_zeros()` 후 "대각만 남은 행"을 고정으로 판정하도록 수정 (실행 전).
- **V02 재실행 (ext OFF, 적분 수정 후, 레벨5=(267,53,2), LDLᵀ, 642 s)**: h=10 10.55%, 5 7.28%, 3 0.724%, 2.5 5.52%, 1.5 0.445%, 0.75 1.292% (fb 0/0/0/2/0/4). 최세 1.29% > 1%, 수렴률 1.03 < 1.8 → **G2 FAIL (ext OFF)**. fallback 제거 효과는 L1 −1.0%p, L3 −1.0%p뿐 → 주원인은 미확장 소지지 기저(분석 가설 b). 결과 `evidence/v02/results.json`, 로그 `evidence/v02_run2.log`. 다음: ext θ 0.2/0.3/0.35 sweep (`evidence/v02_theta*`).
- 메모리(사용자 PC 31.6 GB): `JobSpec.linear_solver` 기본 `auto` = follower 없으면 `pardiso_ldlt`(LU 대비 factor 약 절반), 있으면 `pardiso_lu`. V02 스크립트도 `auto`. 시험 12 passed.
- **V02 판정 방식 변경 (사용자 결정 2026-11-15)**: Coreform 원문(`roark-17p1-7a.md` §7.5, `-qoi.csv`) 기준. 격자 n_cells = (33,7,1)·2^k (h ≈ 6/2^k, nz=1), peak는 Roark·Abaqus(4.34797e-6) 둘 다, 수렴률은 x=0, y=10→20 최대 주응력 분포의 최세 해 대비 L2 자기 수렴률 ≥ 1.8. 근거: Coreform p=2 peak의 Roark 대비 오차도 비단조(24.6/5.66/0.117/0.228/0.826%), Roark–Abaqus 0.7% 차. 실측 근거: 구멍 테셀레이션을 촘촘히 하자 h=1.5 peak가 4.3367e-6(+0.45% Roark) → 4.3642e-6(+1.08% Roark, +0.37% Abaqus). 계획서 §9 G1·§10.1·§0.5 갱신. 구 격자는 `WHT_V02_SET=legacy`. 실행 `evidence/v02_cf_{off,0.2,0.3,0.35}`.
- **V02 원문 기준 결과 (ext OFF, `evidence/v02_cf_off`, 1096 s)**: peak(Roark/Abaqus) h=6 19.75/20.32%, 3 8.63/7.87%, 1.5 2.77/2.05%, 0.75 1.06/0.35%, 0.375 **1.157/0.448%**. 분포 자기 수렴 오차 1.51e-1, 3.12e-2, 9.25e-3, 1.02e-3 → **rate_all 2.34 (≥1.8 PASS, Coreform 2.23)**, rate_last3 2.47. peak(Roark) 1.157% > 1% → FAIL(0.16%p). 가설: peak를 테셀레이션 내접 다각형 위 최대값으로 잼(구멍 쪽 외삽) vs Coreform은 점 A 평가 → 스크립트에 `peak_A`(정확한 원 위 (0,±10,0)) 추가, 레벨별 U 저장(`v02_h*_U.npz`). θ=0.2: ext 효과 ≤0.01%p. h=0.375 풀이 151.6 s(54만 DOF, LDLᵀ) 정상.
- 고속화: 응력 변환 `np.array(list(vals))`(20×, 후처리 62% 부분), 체적 peak 일괄 eigvalsh(`probes.peak_max_principal_fast`, 결과 비트 동일). 남은 병목: `classify_cells` 직렬 136 s(h=0.375), Kratos 응력 계산(38%).
- **점 A peak (Coreform 방식, 정확한 원 위 (0,±10,0))**: θ=0.3 실행에서 h=6 3.207e-6(Roark −25.7%, Coreform 3.257e-6 −24.6%와 1.5% 차) / h=3 +7.9%(Coreform −5.66%) / h=1.5 +2.78%(Coreform +0.12%) / **h=0.75 0.912%(Roark), 0.205%(Abaqus)**. h=3·1.5의 과대는 crown fallback 셀(4–5개)과 격자–crown 정렬 의심.
- 남은 fallback 진단(h=3, 4셀): 셀 기저 sub-box 후보 잔차 2.4e-7에서 정체(> 완화 1e-8) → **`_NNLS_RELAXED` 1e-8 → 1e-6** (fallback 모멘트 오차 ~1e-3 대비). 별도 결함: 대칭 셀 2개에서 sub-box 기준 모멘트 M000이 셀 값과 불일치(0.240 vs 0.077) — 체적 검사로 안전 차단 중, 원인 미조사(TODO 아님: 기록된 미해결 항목).
- θ 0.35·offA 대기 실행 취소. θ=0.3(1e-8 기준) 완료 후 `evidence/v02_cf2_{off,0.3}`(1e-6 기준, 점 A 포함)으로 재실행.
- **θ=0.3 (완화 1e-8) h=0.375 점 A peak: 4.346793e-6 → Roark 0.678%, Abaqus 0.027%** (Coreform 0.83/0.12%). 분포 수렴률 2.34. 점 A 기준이면 G2 V02(ext ON) 조건 충족. 스크립트 판정을 점 A peak로 변경, 계획서 §9 G1 반영. 확정은 `v02_cf2_{off,0.3}`(완화 1e-6) 결과로.
- **V02 PASS (ext OFF, 완화 1e-6, 원문 기준, `evidence/v02_cf2_off`)**: 점 A peak h=6/3/1.5/0.75/0.375 = 25.72/7.92/2.78/0.906/**0.674%**(Roark), Abaqus 대비 최세 **0.031%**; fallback 0/0/1/0/0; 분포 자기 수렴 1.51e-1, 3.13e-2, 9.24e-3, 1.02e-3 → **rate 2.34** (last3 2.47). G1 V02 충족. G2용 ext ON(θ=0.3)은 `v02_cf2_0.3` 실행 중.
- **G2 V02 (ext ON, θ=0.3) PASS** (`evidence/v02_cf2_0.3`, 497 s): 점 A peak 최세 0.678%(Roark)/0.027%(Abaqus), fallback 0, 분포 자기 수렴 rate 2.34 (last3 2.46). ext OFF와 ≤0.005%p 차.
- **V06 Voigt 결함 수정**: `test_p2mesh_v06geo.lame_errors`가 Kratos 응력을 [xx,yy,zz,yz,xz,xy]로 읽었고, 자기시험 `_exact_voigt`도 같은 틀린 순서로 데이터를 만들어 by-construction 통과. Kratos 실측 순서 [xx,yy,zz,xy,yz,xz] 확인(`tests/kratos_spikes/test_s04_voigt_order_kratos.py`, 3 passed). 둘 다 수정.
- **V06 레벨0 재측정 (ext θ=0.2, LDLᵀ, n_geo=24)**: p=1 u 1.05e-2 / ur 9.9e-3 / srr 0.248 / stt 0.0264 (tol 0.006/0.006/0.4/0.025), p=2 u 9.7e-3 / ur 9.4e-3 / srr 0.065 / stt 0.0228 (tol 0.001/0.001/0.025/0.002). 이전 u 0.30, srr 1.87에서 크게 개선. p=1≈p=2 → 비이산화 오차 바닥 의심(구 테셀레이션 n_geo=24 내접 다면체; 추정 0.4–0.8%). fb=50(q=2). n_geo 24/96 × ext on/off 실험 실행 중.
- **모멘트 계산 결함 수정 (`quadrature._cross_section_moments`, 2026-11-15)**: x=x1 단면의 직사각형 변 분할점을 (i) 3D 삼각형 변 투영 교점으로 구하고 (ii) 삼각형 변 매개변수를 직사각형 변 매개변수로 저장, 게다가 (iii) `elines` 선 순서가 변 순서와 불일치 → 경계 닫는 변 누락, 단면 면적 음수. V06 셀 M000 −0.193/−0.056 등(실제 ≈0), 축 회전에 따라 값 변동. 수정: 단면 선분 끝점과 직사각형 변의 교점을 변 매개변수로(`_split_at`), x1 평면 위 삼각형은 2D 변으로. 결과: M000이 표본 체적과 일치(0.8550 vs 0.8551), 축 회전 불변. 회귀 시험 `test_p2mesh_v06geo.py::test_cell_moment_volume_matches_sampling_and_rotation`(4건). 분석 문서 2-B 처방 문구("단면 교선 선분의 교점을 params에")가 이 결함을 가리켰음. V02 결과도 영향 가능 → `v02_cf3_{off,0.3}` 재실행, V06 레벨0 `v06_l0_momfix.log`.
- **V06 레벨0 (모멘트 수정 후)**: p=2 n_geo=24: q=2 → u 9.44e-4 / ur 8.93e-4 / srr 0.0325 / stt 1.74e-3 (srr 불합격), **q=4 → u 9.17e-4 / ur 8.72e-4 / srr 0.0159 / stt 1.04e-3 → 4항목 PASS** (Coreform 9.40e-4 / 9.30e-4 / 0.0185 / 1.60e-3; tol 1e-3/1e-3/0.025/2e-3). n_geo=96: u 3.3e-4, srr 불변 → srr은 형상 아님. β 0.1–1000배 무영향. ext on/off 무영향. **원인: q=2가 p=2 강성 피적분식(방향별 4차) 미달 적분 → V06 기본 q=2p로 변경.** p=1 q=2: u 1.02e-2 > tol 0.006 (Coreform 0.0042–0.0058) — 조사 중.
- V02 `v02_cf3`(모멘트 수정 후) h=6..0.75 결과가 `cf2`와 동일 → V02 형상에는 결함 미발현.
- **V02 최종 (모멘트 수정 후 재확인)**: `v02_cf3_off` 0.674%/0.031%, rate 2.34 PASS; `v02_cf3_0.3` 0.678%/0.027%, rate 2.34 PASS. **G1·G2 V02 확정.**
- **V06 p=1 레벨0**: u≈1.0e-2 고정(ext 0.2/0.35, q=4, n_geo=96 모두 동일), stress는 Coreform보다 양호(srr 0.24 vs 0.285, stt 0.016 vs 0.019). Coreform p=1 u: L0 0.0042–0.0058, L1 0.0020, L2 0.00065. 레벨1(n_geo=96, q=2p, ext 0.2) 실행 중 `evidence/v06_l1_r2.log`.
- **고속화 (결과 비트 동일)**: `classify._axis_votes` 격자선 그룹화 unique/argsort(기존 O(lines×verts) isclose) + 광선별 삼각형 AABB 선필터 → V06 20³ n_geo96 42 s→3.2 s, V02 264×56 7.6→3.7 s. (스칼라 MT 커널은 FMA/순서로 경계 꼭짓점 판정 246개가 달라져 폐기.) `trimmed._cell_rules`: 셀×삼각형 ≥ 2e6이면 셀<400이라도 병렬. 캐시 키: classify.py·quadrature.py 전체 + trimmed의 `_build_trimmed_table`·`_eval_rule_chunk` 소스만 (2026-11-18 스레드-안전 정리로 모듈 전역 `_W`/`_rule_chunk` 제거). unit 326 passed.
- 다음: V06 L1 결과 → p=1 판단, V06 L2, V-COND(고정자유도 판정 수정됨), 메쉬 파일 입력 V02, G2 보고서. 남은 미조사: V02 sub-box 모멘트 불일치(모멘트 결함 수정으로 해소됐을 가능성 — 재확인 필요).

## 2026-11-15 오후 (새 세션) — G2 마무리 실행

- **DLL/MKL 결함**: 새 세션 PATH(whtiga env 경로 없음)에서 첫 Pardiso 호출 시 조용한 종료(`0xC06D007F`, delay-load 함수 없음). condadeps 선두로 고치자 Kratos MKL 2026.1(libomp) + pypardiso MKL 2026.2(libiomp5md) 동시 로드 → `OMP: Error #15` 중단. 최종: `_bootstrap.configure_paths`가 `<sys.prefix>/Library/bin/mkl_rt.3.dll`을 전체 경로 선로드 + PATH 선두 + `PYPARDISO_MKL_RT` 지정 → 프로세스당 MKL 1개(2026.2)·OpenMP 1개(libiomp5md) 확인(psutil 맵). `WHT_IGA_MKL_DIR`로 변경 가능. 새벽 실행들은 miniconda3 MKL 2026.0을 썼음(로그 배너 2026.0.0). V02 h=0.75 재현 동일(0.998%).
- **V06 메모리 원인 실측**: L2 p=2 고스트 채움 — OUT 셀 51,529/64,000(80.5%), 자유도 222,264 중 활성 62,151, 적분점 239만 중 고스트 139만. → `out_mode="drop"`(`wht_iga/solver/active.py`, JobSpec 기본 "ghost" 유지). drop 결과 = ghost (L0 p1/p2, L1 p1/p2 유효숫자 3–4자리 일치). **L2 p=2 drop PASS: u 1.429e-5 / ur 1.358e-5 / srr 8.67e-4 / stt 4.88e-5, peak 6.6 GB, 674 s.**
- **V-COND 재설계·결과**: 첫 판은 상자 평행이동이 판을 자르기만 해 작은 셀 없음(ties 0, ON=OFF, 고스트가 κ 지배) → 공허. 개정: 마지막 x열 φh 조각(φ∈[1e-8,1], 셀 η_min 1e-8 도달), drop, 미수렴(OFF 특이) 시 `ValidationError.strategy`로 K만 사용. 결과: 평탄 ON 1.06(PASS, η 1e-8 도달·ties 120), OFF는 φ=1e-4에서 κ 6.6e16(반올림 수준 kmin, 정성적)·φ=1e-8은 ON과 동일(특이 모드가 pivot 보정에 가려짐), 기울기 −1.44 FAIL(두께 방향 nz 1–4 → 점근 전; kmin = 외팔보 굽힘 모드). 기준 개정 여부 사용자 결정 필요.
- **V02 메쉬 파일**: h=0.75 점 A 0.998%/Abaqus 0.290% (STEP 0.912%/0.205%). h=0.375는 LDLᵀ ~18 GB 초과로 2회 메모리 감시 중단(여유 13 GB). AMGCL+CG(`linear_solver="amgcl_cg"` 추가)는 h=0.75에서 10분+ 결과 없이 중단(원인 미진단) → 현 설정 기각(옵션은 남김, 대형 SPD 실험용).
- 실수 기록: V06 drop 첫 실행은 runner의 `__main__` 보호 누락(spawn 재실행)과 `JobSpecLite.out_mode` 누락, L1 n_geo 24(증거는 96) 오설정 → 모두 수정 후 재실행. n_geo 24 L1 p=2 drop(u 9.08e-4)은 형상 오차 비교용으로만 보존(`L1_p2_drop_ngeo24.json`).
- 단위시험 332 passed (+ bootstrap 2, active mask 3, out_mode 1).

## 2026-11-15 15:10 — 분석 문서 반영 후 인계 (G2_COMPLETE 개정 2)

- `RESPONSE_ANALYSIS_G2_BOTTLENECKS.md` §1 NNLS 가지치기 버그 확인·수정(결과 불변, 셀 규칙 L2 p=2 4067→1169 s, L0 p=3 23분+→120 s). §3.1 CSV 허용치 정정(판정 불변). §3.2 BC 가설 실험 기각. §4 V02 옵션 추가(미실행).
- 신규 진단: V06 p≥2 오차 = 구 면 분할 한계(u/sagitta 1.71 일정). L0 p=3 u 9.135e-4 FAIL(p=2와 동일). L2 p=3 보류(사용자 결정).
- 실행 중 작업 없음(사용자 요청 중단). 다른 에이전트 인계: `docs/phase2/G2_COMPLETE.md` 개정 2 블록과 "다음 작업" 0번부터.
