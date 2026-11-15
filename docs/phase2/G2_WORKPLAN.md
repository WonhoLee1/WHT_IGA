# G2 마무리 작업 계획 (통합 담당, 2026-11-15)

사용자 지시: 아래 1–6을 끝내고, 테스트 과정의 미적용 고속화 점검을 한 뒤 **완료 표식**을 남긴다.
다른 에이전트는 표식 파일이 생기면 후속 작업(Phase 3 등)을 시작한다.

## 완료 표식 (다른 에이전트용 계약)

- 경로: `WHT_IGA/docs/phase2/G2_COMPLETE.md`
- 생성 조건: 아래 체크리스트 1–7이 모두 끝났을 때만 생성 (도중 생성 금지).
- 내용: 첫 줄 `STATUS: COMPLETE` (게이트 미통과 항목이 있으면 `STATUS: COMPLETE_WITH_OPEN_ITEMS`),
  게이트 항목별 판정표, 증거 파일 경로, 미해결 항목, 다음 에이전트 작업 목록.
- 표식이 없으면 = 통합 담당 작업 진행 중. 다른 에이전트는 시작하지 않는다.

## 체크리스트

| # | 작업 | 상태 |
|---|---|---|
| 1 | V06 레벨0 p=1 (u 0.0102 > tol 0.006) | 원인 조사: locking 기각(반경 편향 −0.15%, 접선 0.2%) → 잔여는 평균 0 진동(절단 셀 제어점). ext/q4/n_geo96 무영향. L1부터 PASS. **판정 보류 항목으로 보고서 기재** |
| 2 | V06 레벨1 | **PASS** p=1 (2.94e-3/2.87e-3/0.138/8.1e-3), p=2 (5.76e-5/5.47e-5/3.37e-3/1.94e-4) (`evidence/v06_l1_r2.log`) |
| 3 | V06 레벨2 | p=1 PASS(ghost). **p=2 PASS (drop: 1.429e-5/1.358e-5/8.67e-4/4.88e-5, 6.6 GB)**. p=3 미실행(q=6 NNLS 비용) → 미해결 |
| 4 | V-COND | 재설계 후 실행: 평탄 ON 1.06 **PASS**(η 1e-8 도달, ties 120), OFF 폭증은 φ=1e-4에서만 뚜렷(정성적); 기울기 −1.44 **FAIL** → 기준 결정 대기 |
| 5 | V02 메쉬 파일 입력 | h=0.75 0.998%/Abaqus 0.290% **조건부 PASS**; h=0.375 LDLᵀ ≥18 GB로 미실행 |
| 6 | `docs/phase2/G2_report.md` | **확정** (COMPLETE_WITH_OPEN_ITEMS, PENDING 0) |
| 7 | 미적용 고속화 점검 (완료) | 발견: 셀 적분 규칙의 winding이 후보점마다 전체 표면(V06 n_geo96 4만 삼각형) → L1 p=2 1172 s. 처방: 근거리 삼각형 정확 + 원거리 기여 셀 중심 1회(결과 동일성 검증 후 적용). 이미 적용: classify 그룹화/선필터, 병렬 기준, 응력 변환, peak 일괄 |

## 확보된 결과 (상세 `docs/INTEGRATION_HANDOFF.md`)

- V01-ext PASS (p=1–3 ≤ 3.1e-7). V02 PASS (원문 기준, ext OFF 0.674%/0.031%, ext ON θ=0.3 0.678%/0.027%, 분포 rate 2.34; `evidence/v02_cf3_*`).
- V06 L0 p=2 q=4 PASS (u 9.17e-4, ur 8.72e-4, srr 0.0159, stt 1.04e-3). V06 L1 p=1 PASS (u 2.94e-3, ur 2.87e-3, srr 0.138, stt 8.1e-3).
- 결함 수정: extension Marsden 계수, 모멘트 단면 분할, V06 Voigt 순서, cut-cell fallback(sub-box+relaxed 1e-6), int64 키, 국소 좌표, pipeline NameError.

## 진행 메모 (2026-11-15 12:00, 새 세션)

- **환경 결함 수정**: 새 세션 PATH에 Kratos 의존 DLL 폴더가 없어 첫 Pardiso 호출에서 조용히 종료(0xC06D007F, delay-load 함수 없음). `kratos_backend/_bootstrap.configure_paths`가 DLL 폴더를 PATH 앞에 넣음. 시험 `tests/unit/test_kratos_bootstrap_path.py`.
- **V06 L2 메모리 원인(실측)**: 고스트 채움. 64,000셀 중 OUT 51,529(80.5%)에 1e-8E 요소. p=2 L2: 자유도 222,264 중 활성 62,151(28%), 적분점 239만 중 고스트 139만(58%). → `out_mode="drop"` 구현(`wht_iga/solver/active.py`, pipeline/JobSpec `out_mode`, 기본 ghost 유지, V06 runner 기본 drop). 시험 `tests/unit/test_active_mask.py`, `test_solve_validate.py::test_out_mode_*`.
- **V-COND 설계 결함**: 첫 판은 상자 평행이동으로 판을 잘라내기만 해 작은 셀이 생기지 않음(ties 0, ON=OFF), 고스트 1e-8E가 κ를 지배. 개정: 마지막 x열에 φh 조각(φ∈[1e-8,1], 셀 η_min 1e-8 도달), drop, 비공허성 검사(η≤1e-6 도달 + ties>0). 첫 판 수동 결과(참고): flat 1.00, slope −2.68.
- 실행 체인: V02 메쉬 파일(L3 errA 0.998%/Abaqus 0.290%) → 개정 V-COND → V06 drop(L0/L1 p=1,2 재검증 → L2 p=2, `evidence/v06_drop/`).

## 진행 메모 (compact 대비, 2026-11-15 04:50)

- 실행 중: V06 L2 p=1→2→3 (`scripts` 아님: scratch `v06_l2.py` = `run_v06((40,40,40), p, 192, ext_theta=0.2)`, 로그 `phase1/evidence/v06_l2.log`, 메모리 `v06_l2_mem.log`, 2 GB 미만이면 자동 종료).
- 대기 체인: V06 L2 종료 후 → V02 메쉬 파일(`WHT_V02_GEOM=phase2/evidence/v02_meshfile/roark_tet.msh`, levels 3,4, θ=0.3, 로그 `phase2/evidence/v02_meshfile/run.log`; rate는 2레벨이라 nan → peak만 판정) → V-COND(`phase1/evidence/vcond_run.log`).
- 보고서 초안 `docs/phase2/G2_report.md` 작성 완료 — PENDING 3행(메쉬 파일, V06 L2, V-COND)만 채우면 됨. 그 뒤 `G2_COMPLETE.md` 생성.
- 고속화 점검(7) 완료 기록: 보고서 §4.

## 종료 (2026-11-15 13:30)

- 체크리스트 1–7 종료. 독립 검증 1차 REJECT(증거 불일치·재현 명령·V-COND 과장·수명 결함) → 수정 → 2차 APPROVE.
- **완료 표식 생성: `WHT_IGA/docs/phase2/G2_COMPLETE.md` (`STATUS: COMPLETE_WITH_OPEN_ITEMS`)**. 다음 에이전트는 그 파일의 "다음 작업"부터.
