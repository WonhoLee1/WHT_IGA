# Track H 실행 체크 (계획 대비 실행, 상용급 목표)

> 작업 원칙: 모킹 금지 (전부 Kratos 실측), TDD Red→Green, 증거 기반 허용오차.

## 계획 항목별 상태 (2026-11-19 현재)

| # | 계획 | 상태 | 증거 |
|---|---|---|---|
| 1 | MPC 정석화 (auto→mortar_l2) | 완료 | Red 3.3e-5 → Green 15/15, `test_discretization_methods_kratos.py` |
| 2 | THB 본선 연결 | 완료 (Go) | `test_thb_pipeline.py` 4/4 (단축 1e-9, 굽힘 3%) |
| 3 | 비교 매트릭스 | 완료 | `test_method_comparison.py` 2/2, S35 수치표 |
| 4 | numba 가속 | 완료 (변경 없음) | 웜업 Python 0.104 s / wall 0.7 s, 5% 컷 미달 → leggauss 호이스팅 integrator 이관 |
| 5 | 문서 (S34-S37, TRACKH) | 완료 | 비교표·THB·리포트·실증 기록 |
| 6 | 상용화 (CI·견고성·G-8·사용자 가이드) | 부분 진행 | verify_report + evidence 세트 + presets + 기법 선택 가이드 + API 분할 완료; CI·G-8 잔여 |
| - | recommend_discretization + OBB | 완료 (전기간) | 단위 16 + OBB 2/2 |

## 테스트 결과 원장 (요약; 전체: `.opencode/plan/test_results_ledger_20261118.md`)
- 교차 검증 15/15 (7.65 s), 비교 매트릭스 2/2, THB 파이프라인 5/5 (3.0 s),
  단위+API+격리 47/47, 분할 후 전체 91/91, 광역 unit 604 passed.
- 광역 unit 4 failed = 외부 gmsh 바이너 wrapper의 `gmsh` 파이썬 모듈 부재
  (miniconda3 환경 문제, 당사 코드 무관·미수정 파일, 선재함).
- THB: 단축 exact (1e-15), 굽힘 p=2 3.78e-4 / p=3 5.98e-3, p=1 굽힘 제외 (47.6%).
- barycentric: O(h) 정합 + 게이트 2e-3 (근본 원인 미특정, mortar 정석).
- penalty 1/beta 법칙 (3 decade 실측).

## 검증 커버리지 (전체: `.opencode/plan/verification_coverage_20261119.md`)
- PASS: V01, V02, V03, V06, V-COND-split, V-EXP/C, V08 R1/R2(h=0.25),
  E03 축소, V04, V05, E01, E02 (G3 9/9), V09, V10, V11 (G4 10/10, 378 s),
  V12 (인용, G4_report 수치 일치).
- UNKNOWN: 없음 (G3/G4로 해소).
- NOT STARTED: V07, E04/E05 풀, V08 풀 3레벨, E03 풀 13셀.
- 실증 자료: `docs/verification/report/` (HTML+PNG 4종+data.json),
  `evidence_bulk_bar/` 6종 (vtkhdf/vtu/msh/fem/stl/png).

## 남은 결정 사항
- V-COND gate closure, e04.json 표기, push 승인, C++ 마일스톤, V08 스윕/E03 풀 착수 판단.
