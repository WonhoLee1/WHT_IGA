# P25 cond — 조건수 도구 + V-COND

- 작업: T-P2.5
- 결과: 부분 완료 (도구 + Kratos 앵커 + analog study 완료 / V-COND sweep은 차단됨)
- 근거:
  - `wht_iga/stabilize/cond.py` (신규): `reduced_matrix`(TᵀKT, T=None 항등), `cond_eigsh`(shift-invert 극소 + LM 극대), `cond_dense`(SVD oracle), `h_slope`, `flatness` — sparse 전용.
  - 도구 검증:
    - 1D Laplacian 해석 스펙트럼과 1% 일치, dense 기준과 5% 일치 — PASS.
    - Kratos GetSystemMatrix 앵커: 무구속 추출 K == 참조 조립기 (1e-10), κ 동일 (K-A2: kratos 9.3472e+01 = reference, 오차 0%) — PASS.
    - h-slope: full-Dirichlet p=1에서 −2.02 ± 0.3 달성 (교과서 −2 재현, 도구+게이트 동시 검증) — PASS.
  - η sweep analog study (4³ p=2, 내부 1기저 η∈[1,1e-8]):
    - 미안정화 κ: 2.8e2 → 5.8e8 (×2e6 폭증) — `test_eta_sweep…` PASS.
    - extension κ: 281 → 269 (flatness 1.04 ≤ 10) — PASS (V-COND 평탄 패턴의 analog 증거).
    - FCM-α κ: 2.03e8 (미안정화 5.83e8 대비 낮음) — PASS.
  - 중요 관측 (V-COND 설계용, 계획서 §9와의 관계): 단일면 고정 p=2 박스의 κ(h) 기울기는 nc≤12에서 −1.35 (raw) / −1.67 (Jacobi) — 전형적인 coarse-grid pre-asymptotic이며, full 구속 시 −2.02로 회복. V-COND의 −2±0.3 G2 바는 V-COND 자체 레벨에서 판정해야 하며, 본 analog는 그 바를 미리 검증하지 않음 (§10.7 해당 시 보고 — 현阶段은 관측 기록이며 모순 아님).
  - V-COND sweep (구멍중심 이동 × 4 h-levels, extension ON κ ≤ 10×): 차단됨 — extension-ON Kratos 해와 행렬 추출 모두 파이프라인 훅 필요 (P24와 동일 사유). 모체 형상: V02 Roark 판 (phase2 §10.4 TBD 해소안 — 사용자 확인 요청).
- 환경: Python 3.12.16 (whtiga), Kratos 10.4.3 (E:\kratos-build\install), 2026-11-14
- 결론과 영향: κ 도구는 Kratos 앵커 완료. V-COND 본 sweep은 파이프라인 구속식 훅 + strategy 노출 후 실행 (G2 V-COND 항목).
- 재현: `python -m pytest WHT_IGA/tests/unit/test_stabilize_cond.py`
  `python -m pytest WHT_IGA/tests/kratos_spikes/test_p2stab_kappa_kratos.py -m kratos`
