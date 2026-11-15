# P22 extension — small-cell extension 연산자

- 작업: T-P2.2
- 결과: 완료
- 근거:
  - `wht_iga/stabilize/extension.py` (신규): slave별 최근접 master-only (p+1)^d 블록(블록 중심 거리 최소, 동점 시 interior 방향 우선 — §7.5 "바깥 법선 우선"의 구현 가능 해석을 주석에 명시), e_ij = Π Lagrange 외삽(@njit), 폴백 체인(차수 저하 → 최근접 master 상수 → THIN_FEATURE + 세분 권고), 계수 감시 기록(cap TBD — 기록만, 강제 없음).
  - 수락 기준:
    - 분할단위성: slave 행합 = 1 (1e-12) — `test_partition_of_unity_1e12` PASS.
    - p차 재생 (p=1,2,3, 1e-10, Höllig 2003 / Marussig 2018) — `test_degree_reproduction_1e10` PASS.
    - θ 연속성: slave/master 중첩 단조 + inactive 불변 — `test_theta_continuity_nesting` PASS.
    - 폴백: 차수저하 deg=1 확인, clean 상수 확인, 원거리 상수 THIN_FEATURE issue(코드·심각도·개체·세분 권고 포함) 확인, master 0개 입력 거부 — PASS.
    - 계수 상한: `stats["cap"] is None` (TBD 기록), `max_abs_coef` 유한 — PASS.
- 환경: Python 3.12.16 (whtiga), Kratos 불필요, 2026-11-14
- 결론과 영향: T-P2.3 합성 입력 준비 완료. 계수 상한값 확정은 후속 ADR 과제 (측정치: 구멍 격자 p=2 max|e|=3.0, n_master=352, n_slave=32, 상수 폴백 0, THIN_FEATURE 0).
- 재현: `python -m pytest WHT_IGA/tests/unit/test_stabilize_extension.py`

## 정정 (2026-11-14, 통합 담당 — V01-ext)

- e_ij를 index 공간 Lagrange로 구하던 방식은 clamped 매듭(p≥2) 경계 근처에서 Greville 점이 비균일해 affine도 재현하지 못했다 (u=y slave 계수 오차 p=2 0.125, p=3 0.083; V01-ext B케이스 1.3e-2/4.2e-3). 기존 `test_degree_reproduction_1e10`은 다항식을 index 공간에서 만들어 by-construction으로 통과했다.
- 수정: 1D 가중치를 실제 매듭의 단항식 B-spline 계수(Marsden blossom, `e_k(t_{i+1..i+p})/C(p,k)`)를 재현하는 연립식으로 구함 (`extrapolation_weights_1d`). `build_extension(..., knots=)` 필수 인자. 균일 내부에서는 기존 결과와 같음 (`test_interior_matches_index_lagrange`).
- 시험: 재생 시험을 물리 좌표 다항식(Greville collocation 오라클)으로 교체, `test_clamped_boundary_affine_regression`, `test_knots_mismatch_raises` 추가. compose 시험 2개(단위·Kratos)도 물리 좌표 필드로 교체. 상세 `docs/verification/V01-ext.md`.
