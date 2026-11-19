# S47: HardeningTable 및 PEEQ 동기화 검증 (`InitialState` 경로)

**Date:** 2026-11-19  
**Module:** `wht_iga/analyses/thermal_plasticity.py` (`HardeningTable`, `point_curve_margin`, `CellThermalGroups`)  
**Verification Suites:**
- `tests/kratos_spikes/test_s44_hardening_initial_strain.py` (`HardeningTable` 검증 스위트; 파일명의 `S44`는 기존 명명을 유지한 것이며 `docs/spikes/S44.md`는 PrePoMax 벤치마크 문서)
- `tests/kratos_spikes/test_p5_thermal_plasticity.py` (완전소성 및 열탄성 기준 스위트)

---

## 핵심 요약

1. **구현 확인 (`wht_iga/analyses/thermal_plasticity.py`)**
   - `HARDENING_CURVE = 6` (`CurveDefinedByPoints`) + 끝단 평탄 연장(`FLAT_EXTENSION`)으로 Abaqus `*PLASTIC` 다중선형·온도의존 등방경화 `sigma_y(PEEQ, T)`를 구현.
   - `point_curve_margin`: Kratos 점 곡선의 세그먼트별 표현 가능성 조건 $a_j^2 = S_j^2 - 2 h_j W_j > 0$을 `HardeningTable` 생성 시 및 그룹 온도 보간 곡선 적용 시 검증.
   - PEEQ 동기화(`CellThermalGroups.commit` / `apply`): Kratos의 반복말 응력 기반 $W_p$ 적분 과잉 누적 및 온도의존 곡선에서의 일(work) 기반 조회 결함을 제거하고, 매 증분 시작 시 `SetValuesOnIntegrationPoints(PLASTIC_DISSIPATION)`으로 `kappa = l_c W_T(PEEQ) / G_f`를 기록.

2. **테스트 기준·증분 결함 수정 (`tests/kratos_spikes/test_s44_hardening_initial_strain.py`)**
   - `test_c_multilinear_curve_uniaxial_stress`: 70증분(`deps = 0.001`)에서는 첫 구간(`H = 80`)에서 증분당 해석 상계 `_bound` 자체가 `~0.0045` 및 꺾임점 부근 `0.010478`로 `max(tols) < 2e-3` 조건을 초과함. 증분 수를 `200`으로 세분하여 `max(tols) = 0.001078 < 2e-3`과 Kratos 해석 결과 검증을 모두 통과.
   - `test_f_three_cells_two_hardening_groups_residual`: 가열 종료 시점(`i_heat = 15`) 기준응력 검증에서 9자리 반올림 리터럴(`-0.182411067 / 0.065573123`)로 인한 `2.29e-8` 절단 오차를 정확한 닫힌형 유리식 `-(0.07 + 3.6 / 55.0 + 2.7 / 57.5) / (0.03 + 1.0 / 55.0 + 1.0 / 57.5)` (`-4615 / 1659`, `abs=1e-12`)로 교체. Kratos 잔차·응력·PEEQ 허용오차는 무수정 유지.

---

## 재현 및 실행 결과

```powershell
D:\mamba\envs\whtiga\python.exe -m pytest -q -p no:cacheprovider tests/kratos_spikes/test_s44_hardening_initial_strain.py tests/kratos_spikes/test_p5_thermal_plasticity.py
```

- **결과:** `25 passed in 27.29s` (15 passed in `test_s44_hardening_initial_strain.py` + 10 passed in `test_p5_thermal_plasticity.py`).
