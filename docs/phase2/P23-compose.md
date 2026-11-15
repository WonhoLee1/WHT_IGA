# P23 compose — 연산자 합성 + 구속식 배치 생성

- 작업: T-P2.3
- 결과: 완료 (단위) / Kratos 확인 진행 중
- 근거:
  - `wht_iga/operators/compose.py` (신규): T = C_THB·E 합성 경로 (Phase 2 C_THB=I), `apply_TKT`/`apply_Tt_vec`/`solve_reduced`, `expand_dof_operator`(성분 블록대각 확장) — scipy.sparse 전용 (dense 입력 거부).
  - `wht_iga/operators/numbering.py` (신규): `ReducedNumbering`, 스칼라 스펙 배치 생성, c≠0 거부 (`require_zero_constant` — kratos_backend/constraints.py 미수정, STAB 자체 검사), 연쇄구속 거부 (`validate_no_chains`).
  - 수락 기준:
    - TᵀKT 일치 ≤ 1e-10 (V01 축소판, 참조 조립기 vs KKT 독립 경로) — `test_tkt_agreement_vs_kkt_1e10` PASS.
    - 1e5 스칼라 스펙 생성 2.85 s (≤ 30 s; S05 Kratos 생성 2.7 s와 구분 — 본 측정은 스펙 생성, Kratos 생성은 K-B에서 측정) — PASS.
    - 1e5 Kratos `LinearMasterSlaveConstraint` 생성 4.42 s (S05 2.7 s 대비, ≤ 30 s) — K-B PASS.
    - c≠0 거부, 연쇄 거부 단위시험 — PASS.
    - Kratos 경로 TᵀKT (K-A): twin 무구속 모델 추출 K == 참조 K (1e-10, 자유 DOF), tied 모델 Kratos 해 == T-축소 scipy 해 err=7.461e-14 (≤ 1e-10), slave recovery=0.0 (≤ 1e-12) — PASS.
- 환경: Python 3.12.16 (whtiga), 2026-11-14
- 결론과 영향: R1 구속식 경로의 STAB측 합성/배치 준비 완료. K-A/K-B 결과로 G2 증거 보강.
- 재현: `python -m pytest WHT_IGA/tests/unit/test_operators_compose.py`
  `python -m pytest WHT_IGA/tests/kratos_spikes/test_p2stab_compose_kratos.py -m kratos`
