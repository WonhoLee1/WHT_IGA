# 검증 리뷰 스크립트와 결과 (2026-11-23)

이 폴더는 V03, V04, V06, V07, E03, E04, V02 검증 판정 개선 작업의 실행 스크립트, 결과, 그림을 묶은 것이다. 스크립트는 작업 당시 scratchpad 경로를 그대로 쓴다. 저장소 기준으로 다시 돌리려면 경로를 바꿔야 한다.

## 스크립트
- `v03_sweep_spec.py`: V03 끝단 전단, 스펙 E=2.0e11, 격자 (0.1,0.1,0.25)/2^k, p=1–3 k=0–4.
- `v07_p3_standard.py`: V07 p=3 C⁰/P2, 표준 격자 0.5/0.25/0.125 (`mixed_p1` 헬퍼 사용).
- `v04_run.py`, `v04_sweep2.py`: V04 이종 보, 병합 vs tie, 명세 격자 (5,5,10)/2^k, p=1–3, k=0–3.
- `e04_convergence.py`, `e04_c0_64.py`: E04 Cook 막 C¹/C⁰, p=1,2, 격자 16–64.
- `e03_csv_levels_hist.py`: E03, Coreform 증분 이력(15/단계)을 하중 레벨로 사용, 시도 기록 저장.
- `e03_csv_levels_2x.py`: E03, 증분을 2배로 나눈 실행. 레벨 1.0002에서 수렴 실패로 중단(JSON 없음).
- `v06_variants.py`, `run_v06_theta.py`: V06 p=1 L0 진단(확장 θ, 경계 방식, 적분 차수, 기하 분할).
- `t3_1d.py`: V12 T3 1차원 시간 적분 비교(lumped, consistent, 혼합).
- `fig_*.py`, `make_figs.py`: 보고서 그림 생성.

## 결과 (`results/`)
- `v03_sweep_spec.json`, `v04_sweep2.json`, `v07_p3_standard.json`, `e04_convergence.json`, `e04_c0_64.json`, `v06_variants.json`.
- `e03_csv_levels_hist.json`: max_iters=60으로 돌린 실행(기본 30에서 같은 컷백 3회, 같은 위치). 기본 30 실행의 기록은 E03.md에 요약.

## 주의
- E03 컷백 3회는 `cool_to_minus_55`가 아니라 `heat_to_200` 첫 증분 근처(레벨 2.00–2.01)에서 난다.
- 원인은 아직 확정되지 않았다. Kratos Newton 잔차는 Python 표준 출력으로 나오지 않아 추적하지 못했다.
