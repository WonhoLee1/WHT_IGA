# P21 support — basis-support classification

- 작업: T-P2.1
- 결과: 완료
- 근거:
  - `wht_iga/stabilize/support.py` (신규): S3/S4 입력(셀 내부체적, 셀 체적, 백그라운드 knot) → 기저별 지지 체적비(@njit 누적) → master/extended/inactive 분류. θ 기본 0.2, 범위 [0, 0.5) 강제.
  - 수락 기준:
    - untrimmed 균일 격자: 전부 master, inactive 0 (ratios 정확히 1.0) — `test_untrimmed_all_master` PASS.
    - 구멍 격자: extended 집합이 절단 경계 밴드와 일치 — 각 extended 지지가 소절단 셀(0<η<θ)과 만나고, 모든 소절단 셀이 non-inactive 지지에 덮이며, inactive 지지는 전부 영영역 내부 — `test_hole_extended_band_matches_cut_boundary` PASS.
    - θ sweep [0.05, 0.2, 0.35]: slave 집합 단조 증가, master 집합 단조 감소, inactive 불변 — `test_theta_sweep_monotone` PASS.
    - 다중재료 = 합집합 영역(과lapped 입력 거부): 분할 재료 eta 합집합이 단일재료 분류와 정확히 일치 — `test_multi_material_union_matches_single` PASS. 본 Phase에 다중재료 셀 없음 → 단일재료 경로 + "multi-material branch not-run" 명시 보고 (TODO 스텁 아님).
- 환경: Python 3.12.16 (whtiga), Kratos 불필요, 2026-11-14
- 결론과 영향: T-P2.2 입력 준비 완료. G2 extension 검증 경로의 분류 Leg.
- 재현: `python -m pytest WHT_IGA/tests/unit/test_stabilize_support.py`
