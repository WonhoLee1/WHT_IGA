# P24 fcm-beta — FCM α 기준선 + β sweep

- 작업: T-P2.4
- 결과: 부분 완료 (단위 + V01 sweep 완료 / V02 FCM 완료 / extension 비교표는 차단됨)
- 근거:
  - `wht_iga/stabilize/fcm.py` (신규): FCM_ALPHA=1e-8, `fictitious_modulus`, `default_beta`(solver.penalty 단일원천의 STAB측 사본 — 값 일치 단위시험), `beta_sweep` {1e2…1e6}, `loglog_slope`.
  - O(1/β) 법칙: 정확 스프링-페널티 모델에서 기울기 −1 ± 0.01, C=2.0 확인 (상위 3 decade) — PASS.
  - ADR-002 표 대조: default_beta(1, 0.5, p) = 6.0e5/2.4e6/5.4e6 정확 일치 — PASS (ADR-002 미수정, 증거만 확장).
  - K-C (V01 β sweep, Kratos): slope −0.9983 (±0.1 PASS); C_abs = 3.654/3.716/3.723/3.723/3.723 = ADR-002 p=2 값 3.72와 정확 일치 (rel-err 정규화로 3 division — 양쪽 metric 보고). plateau 없음 (1e2…1e6 단조 1/β, 100배 이상 감소) — PASS, ADR-002 "no plateau" 관찰과 일치.
  - K-D (V02 FCM-α, h=10, 20x4x1, p=2): peak=3.299184e-6, err=23.586%, beta=1.2e6(기본값), cut=80, fb=0, umax=2.238e-4 — G1 V02 h=10행(3.30e-6/23.6%)과 일치 (FCM-α 경로의 G1 재현). β 끝점: β=10 → err 23.543%, β=1e5 → err 23.586% — 4 decade에 걸쳐 penalty 기여 무시 수준 (cut-이산화 오차 지배, G1 원인분석과 일치: extension 부재가 주범, penalty 바닥 아님).
  - V02 extension-vs-FCMα 비교표: extension-ON 열 + 조건수 열은 차단됨 — `run_linear`에 구속식 훅이 없고 `LinearResult`가 strategy/행렬을 노출하지 않아 SOLVE 파일 수정 없이는 불가 (금지됨). FCM-α 열만 측정, extension 열은 "미수행(차단)"으로 기록. 공유모듈 변경 제안: `run_linear`에 선택적 `constraints` 훅 + `LinearResult.strategy` 노출.
- 환경: Python 3.12.16 (whtiga), Kratos 10.4.3 (E:\kratos-build\install), 2026-11-14
- 결론과 영향: β 기본값(ADR-002 c_pen=3e5) 유지 — V01 sweep이 독립 증거 추가. V02 extension 재측정은 파이프라인 훅 후 가능 (G2 V02 항목과 연결).
- 재현: `python -m pytest WHT_IGA/tests/unit/test_stabilize_fcm.py`
  `python -m pytest WHT_IGA/tests/kratos_spikes/test_p2stab_beta_sweep.py WHT_IGA/tests/kratos_spikes/test_p2stab_v02_fcm.py -m kratos`
