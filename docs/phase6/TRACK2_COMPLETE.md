STATUS: COMPLETE_WITH_OPEN_ITEMS

# Track 2 Phase 6 Explicit 동역학 (2026-11-16, 리뷰 반영 2026-11-17)

> 2026-11-17 리뷰(`dev_log/track2_review_20261116_opus.md`)에 따라 상태를 `COMPLETE`에서 낮췄습니다.
> - 이유: 계획서 §10.1의 V-EXP·V-EXP-C 기준 일부를 시험하지 않았습니다.
> - 수정한 항목과 미결 항목은 문서 끝 "리뷰 반영" 절에 있습니다.
> - 첫 줄 형식은 저장소 규약(`STATUS: ...`)으로 맞췄습니다.

- 최종 회귀: `test_s24_explicit_dynamics.py` + `test_s24b_hybrid_explicit.py` +
  `test_s25_explicit_contact.py` + `test_vexp_explicit.py` +
  `tests/unit/test_explicit_alpha.py` + `tests/test_isolation.py` →
  **27 passed in 17.91s** (0 failed). 개별: S24 6/6, S24b 4/4, S25 2/2, V-EXP 7/7.
- 환경: Python 3.12.16 (`whtiga`), Kratos 10.4.3 (`E:\kratos-build\install`).

## S24 혼합안 (정정 §33-44 검증 완료)
- Kratos `MechanicalExplicitStrategy` + `Initialize()` 후 `NODAL_MASS`를
  Python HRZ(`T^T M_c T`) 대각 값으로 덮어쓰면 그대로 Kratos 경로 사용 가능.
- 근거 (`test_s24b`): untrimmed tip L2 0.75%, trimmed+T (η=1e-4/1e-6) **1.2%**
  (2026-11-17 정정: 이전 4.36%는 slave 제어점 프로브 + 반사 무시 기준식의
  인공물. 물리 끝단 장 값 + 반사 포함 기준식으로 바꾸고 게이트 5% → 2%),
  overwrite persist 확인 (Kratos는 매 스텝 질량 미갱신), 전 구간 안정.
- HRZ 전부 양수 (p=1/2 × η=1/1e-2/1e-4/1e-6) — 음수 가중치 소멸 우려 해소.

## S25 (직접 penalty 접촉)
- `explicit_penalty_contact_process` 초기화 실패 (GenerateFaces, S22와 동일).
- 대안 확정: 표면층 + MPC tie + Python penalty 내력쌍 + Kratos explicit 적분.
  자유-자유 충돌 200스텝 안정, 운동량 편차 4.4e-16, 침투 0.006.

## V-EXP / V-EXP-C (게이트 통과)
- V-EXP: p=1/2/3 파면 tip L2 < 2%, 파면 도달 0.5s ± 0.03.
- V-EXP-C (dt 0.0025/0.005): 운동량 < 1e-6 (실측 1e-15 수준),
  반발 속도교환 (상부 ~0, 하부 ~v0), 에너지 드리프트 0.18%/0.26% (< 2%),
  침투 ≈ 0.011 (< 0.02). k_pen=500 (k=2000은 dt=0.005에서 접촉 스프링
  분해능 부족으로 +8% 성장 — CFL 이유로 문서화).

## 산출물
- 신규: `wht_iga/dynamics/explicit.py` (HRZ/extension T/overwrite/dt_crit/
  mass floor/paired penalty/`run_explicit_history`),
  `docs/spikes/S25.md`, `tests/kratos_spikes/test_s24b_hybrid_explicit.py`,
  `test_s25_explicit_contact.py`, `test_vexp_explicit.py`,
  `tests/unit/test_explicit_alpha.py`, `tests/verification/data/vexp{,_c}.json`,
  본 파일.
- 기존 파일 수정: 없음 (프레임워크 변경 0건 — 별도 리뷰 불필요).
- 금지 준수: commit/push 없음, L2/L3 벤치 없음.

## 리뷰 반영 (2026-11-17, Claude Opus)

### 수정
- `tests/kratos_spikes/test_s24b_hybrid_explicit.py` 절단 막대 시험 (H2):
  - 프로브: slave 제어점 → 물리 끝단 장 값 `Σ N_i(x_tip) u_i`.
  - 기준식: 고정-자유 막대 반사 포함(2L_trim/c = 1.95 s 이후 기울기 −v0).
  - 게이트: 5 % → 2 %. 실측 1.2 %.
  - untrimmed 시험 docstring 정정(HRZ 덮어쓰기 없음).
- 회귀(2026-11-17): S24 + S24b + S25 + V-EXP + `test_explicit_alpha` → 24 passed.

### 미결 (계획서 §10.1과 차이, 리뷰 H1·M1–M5)
1. V-EXP 피크 응력 2 % 시험 없음.
2. V-EXP 파면 도달 허용치 ±0.03 s(6 %)를 2 %로 강화해야 합니다. 프로브도 x = 5.125가 아니라 x = 5.0 장 값이어야 합니다.
3. CFL·mass scaling 실행 시험 없음. `estimate_dt_crit` 호출처가 없고, S24의 "소형 셀 Δt_crit ≥ 0.5 × 비절단"을 측정하지 않았습니다.
4. V-EXP-C의 차이:
   - 블록이 p=2가 아니라 p=1입니다.
   - 에너지를 운동에너지만으로 판정합니다(변형·penalty 에너지 빠짐).
   - S25에 에너지 검사가 없습니다.
5. 기준 데이터와 시험이 연결되지 않았습니다.
   - `vexp.json`·`vexp_c.json`을 시험이 읽지 않습니다.
   - `vexp_c.json`은 p=2, k=2000이고 시험은 p=1, k=500입니다.
   - `title` 누락으로 `test_verification_core`가 실패하던 문제는 2026-11-17에 해소했습니다(`id`·`title`·`status` 추가, 값 무수정). `vexp_c.json`은 `status: draft`이고 불일치는 `review_note`에 적었습니다.
6. 계획서와 다른 경로 2건은 사용자 승인과 계획서 반영이 필요합니다.
   - 접촉: ContactSM 프로세스 대신 직접 penalty.
   - 질량: Kratos lumped 대신 HRZ 덮어쓰기.
7. `dynamics/explicit.py`의 문제:
   - `build_x_extension`이 dense 행렬이고 +x 끝 전용입니다.
   - x 방향 셀이 1개면 오류 없이 잘못된 T를 만듭니다.
   - `run_explicit_history`는 시간을 0에서 다시 시작합니다.
