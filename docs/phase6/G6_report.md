# G6 게이트 보고서 (Phase 6 검증 — 성능, 반복 솔버, Explicit 동역학, SBM 비교)

- **날짜 / 브랜치 / 환경**: 2026-11-19 / `WHTIGA` / Python 3.12 (`D:\mamba\envs\whtiga\python.exe`), Kratos 10.4.3 (`E:\kratos-build\install`), 16 논리 코어, 31.60 GiB RAM
- **기준 문서**: [master-plan_v2.1](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/master-plan_v2.1.md) §9–§10.5
- **종합 판정**: **Python·Kratos 선행 트랙 PASS (`AMGCL` 스케일링, 산업용 STEP OBB+AMGCL, SBM 비교, V-EXP / V-EXP-C); Full 1M DOF 실측 및 C++ `ExtendedSplineSolidElement` 승격은 C++ 트랙 이관**

---

## 1. Phase 6 워크스트림별 판정 요약

| 워크스트림 / 게이트 | 범위 및 검증 대상 | 판정 | 핵심 실측 근거 | 문서 / 테스트 |
|---|---|---|---|---|
| **G-8 성능·Coreform 비교** | `make_nurbs_volume`, `attach_cell_elements`, `trimmed_table`, `AMGCL`/`Pardiso` 스케일링 및 Coreform E03(`1,018.4 s`) 대비 외삽 | **PASS (프로파일 & 모델 기반 1M 게이트 확인; 1M 실측 미실행)** | Peak RSS `906.9 MB`; `8T` `build_trimmed_table` `~5.3x` 가속; `AMGCL_CG_ILU0`가 `PardisoLU` 대비 `< 50%` 메모리(`254 vs 568 MB` @ `10,125` DOF); `1M DOF` 외삽 `~1,026 s` (`17.1 min`, 선형 상한 `~1,374 s = 22.9 min < 30 min`), **`10.71 GB < 64 GB`** (Pardiso `1,258 GB` 불가) | [perf_vs_coreform.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/perf_vs_coreform.md), [test_p6_perf_vs_coreform.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_perf_vs_coreform.py) (`5 passed in 39.51s`) |
| **S45 산업용 STEP OBB + AMGCL** | `Generator_bracket.step` OBB 정렬 비교 + `Flange.STEP` (`8,505` DOFs) Trimmed+Ext `AMGCL_CG_ILU0` vs `PardisoLU` | **PASS** | `Generator_bracket.step` OBB `BoxSpline`: active ratio `+59.4%`, 배경 셀 `4.36x` 절감(`p=2` 확장 OBB=AABB 동일 안정성 `5.83e-4`); `Flange.STEP` `8,505` DOFs: `AMGCL_CG_ILU0` vs `PardisoLU` 상대차 **`1.41e-9`**, 솔버 RSS 증가 **`3.70x` 절감 (`105 MB vs 389 MB`)** | [S45_bracket_amgcl.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/S45_bracket_amgcl.md), [test_p6_industrial_bracket_amgcl.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_industrial_bracket_amgcl.py) (`3 passed in 55.72s`) |
| **§10.5 SBM 비교 트랙** | Kratos `NurbsGeometryModelerSbm` vs 동일 DOF `cut+ext`·`FCM alpha` 비교 (V02/V06) | **PASS (주 경로 `cut+ext` 우위 확정)** | Kratos SBM은 **2D 전용**(3D V06 N/A)이며 비-dyadic knot `1e-15` 교차 결함 보유; 2D Lamé 동일 DOF 비교에서 `cut+ext`가 SBM 대비 **6.5배~58배 높은 정밀도** 달성 | [sbm_comparison.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/sbm_comparison.md), [test_p6_sbm_comparison.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_sbm_comparison.py) (`4 passed in 57.88s`) |
| **Explicit (V-EXP, V-EXP-C)** | `MechanicalExplicitStrategy` + HRZ lumped `M` + 확장 `dt_crit` + penalty 충돌 (`p=1..3`) | **COMPLETE_WITH_OPEN_ITEMS** | **V-EXP**: 평탄부 피크 응력 `-0.20%..+0.27% <= 2%`, 1차 모멘트 도달 시각 `+0.04%..+0.43%`, Airy 세분 격자 `0.5 v0` `<= 1.77%`; **V-EXP-C**: `p=1` (`0.31%, 1.23%`) & `p=2` (`0.35%, 1.50%`) 전체 에너지 `< 2%`, `mom_dev <= 2.33e-15`; **CFL**: 확장 `dt_ext/dt_U = 1.015..1.110` | [TRACK2_COMPLETE.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/TRACK2_COMPLETE.md), [test_vexp_explicit.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_vexp_explicit.py) |
| **Track H 방법론 & 번들 STEP** | `recommend_discretization`, OBB, 프리셋, `wht_iga/resources/` PrePoMax STEP 14종 벤치마크 | **PASS** | `L_bracket.step`, `Wrench.stp`, `Hertz.STEP`, `Flange.STEP`, `Contact.STEP` 등 14종 `wht_iga/resources/` 번들 완료; 4종 STEP 벤치마크 통과 | [TRACKH_COMPLETE.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/TRACKH_COMPLETE.md), [S44.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S44.md) (`4 passed in 29.77s`) |

---

## 2. 반복 솔버(`AMGCL` vs `AdditiveSchwarz`) 및 1M DOF 병목 분석

1. **선정된 반복 솔버 기본값 (`DEFAULT_AMGCL_PARAMS`)**:
   - `krylov_type="cg"`, `smoother_type="ilu0"`, `coarsening_type="aggregation"`, **`use_block_matrices_if_possible=False`**.
   - `use_block_matrices_if_possible=True` + `spai0` 조합은 Trimmed+Extension 모델의 스칼라 slave 구속식 대각(`mScaleFactor * I_slave`)으로 3x3 블록 집계가 깨져 `8.10e-3` 정체 / `1.005` 발산을 유발하므로 사용 금지([perf_vs_coreform.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/perf_vs_coreform.md) §3.3).
   - `IgaApplication.AdditiveSchwarzPreconditioner`는 2건의 Kratos C++ 결함(`M_invert = omega * M_invert` 순서 오류 및 `ProvideAdditionalData`의 slave row 부재)으로 Trimmed+Ext에서 발산(`30.17` / `4.00e14`)해 프로덕션에서 제외했다.
2. **1M DOF 스케일링 병목 (Coreform E03 178k DOFs = `1,018.39 s` 대비)**:
   - 직접 솔버(`PardisoLU`/`PardisoLDLT`)는 `178k` DOFs에서 `~146 GB`, `1M` DOFs에서 `~1,258 GB`가 필요해 64 GB 워크스테이션에서 **실행 불가**하며, **`AMGCL_CG_ILU0`** 경로만 `~10.71 GB`로 64 GB 게이트를 통과한다.
   - `build_trimmed_table`은 `>= 400` 활성 셀에서 8스레드 풀 가동 시 `~5.3x` 빠르다(`9.44 ms/cut cell`).
   - 남은 3대 핫스팟(1M DOF `AMGCL` 전체 시간 기준):
     1. `AMGCL_CG_ILU0` 선형 해법 (`~60–80%`, `819 s`)
     2. `attach_cell_elements` Python 루프 (`~11–12%`, `~1,850 elem/s`)
     3. MPC `LinearMasterSlaveConstraint` `T^T K T` 조립 오버헤드 (`~2.8x` build 비용) → C++ `ExtendedSplineSolidElement` ($\tilde{\mathbf{C}}^e$) 승격 대상.

---

## 3. 정직한 미결 항목 (Open Items Requiring User Decision or C++ Work)

1. **1M DOF 3D 브래킷 실제 완주 & C++ 최적화 미실행**:
   - 본 Phase 6에서는 `1,029`–`20,904` DOFs (`Flange.STEP` `8,505` DOFs 포함) 실측 + 법칙 외삽으로 1M DOF 예산(`< 30 min`, `< 64 GB`)과 C++ 우선순위(`ExtendedSplineSolidElement`, Numba `@njit` moment fitting, `PardisoFactor` 재사용)를 도출했다. 실제 `1M` DOF end-to-end 실행과 Kratos C++ (`WhtIgaApplication` / `AdditiveSchwarzPreconditioner`) 수정은 후속 C++ 통합 단계 과제로 남긴다.
2. **Explicit 동역학 정책 결정 대기 ([TRACK2_COMPLETE.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase6/TRACK2_COMPLETE.md))**:
   - `V-EXP` 파면 도달 시각 정의(`1차 모멘트` 채택 vs `vexp.json` 격자 `nx=80/160` 상향 vs 허용오차 조정) 및 `V-EXP` 피크 응력의 평탄부(`t ∈ [1.5, 2.0] s`) 평균 정의 승인.
   - `vexp.json` / `vexp_c.json` 정합화, 계획서 대비 HRZ 질량 + 직접 penalty 접촉 경로 승인, `dynamics/explicit.py`의 `build_x_extension` 범용화(`nx=1` / `-x` 확장).
3. **`aieng-ui` F-bar 브리지 / `WHT_IGA` 경계**:
   - F-bar 프로브 Cauchy 응력([test_fbar_probe_stress.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_probe_stress.py))과 `aieng-ui/backend/app/iga_runner.py` 연계는 완료되었으나, `fbar_trimmed.py`의 `StudySpec`/`methods.py` (`full_iga_immersed`) 배선은 후속 통합 대상이다.

---

## 4. 직접 검증된 대표 스위트 목록

- `tests/kratos_spikes/test_p6_perf_vs_coreform.py`: `5 passed in 39.51s`
- `tests/kratos_spikes/test_p6_industrial_bracket_amgcl.py`: `3 passed in 55.72s`
- `tests/kratos_spikes/test_p6_sbm_comparison.py`: `4 passed in 57.88s`
- `tests/kratos_spikes/test_p5_prepomax_step_benchmarks.py`: `4 passed in 29.77s`
- `tests/kratos_spikes/test_vexp_explicit.py` (`test_s25_explicit_contact.py`, `test_vexp_c_dt_sweep` 등): 주요 부분집합 `12 passed in 46.33s`, 서브에이전트 `a130df46` 6개 파일 전체 회귀 `64 passed in 570.56s`
