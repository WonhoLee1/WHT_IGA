STATUS: COMPLETE_WITH_OPEN_ITEMS

# Track 3 Phase 5 순수 Python 선행 (2026-11-16, 리뷰 반영 2026-11-17)

> 2026-11-17 리뷰(`dev_log/track3_review_20261117_opus.md`)에 따라 상태를 `COMPLETE`에서 낮췄습니다.
> - 이유: `run_load_stepping` 결함 2건(이력 초기화, 컷백 오보고)을 확인했습니다. 두 건은 수정했습니다.
> - G5 기준(V08 4 %, E03, E05)에는 아직 도달하지 않았습니다.
> - 상세는 문서 끝 "리뷰 반영" 절에 있습니다.

- 최종 회귀: `tests/kratos_spikes/test_p5_nlgeom_contact_prep.py` →
  **9 passed in 3.27s** (0 failed).
  내역: NH UL/TL 해석해 2건, J2 완전소성 1건, J2 선형경화 1건,
  load stepping ramp 1건, tied surface layer 1건,
  바이메탈 열응력 1건, Hertz SSN 1건, 마찰 plateau 1건.
- 환경: Python 3.12 (`whtiga`), Kratos 10.4.3 (`E:\kratos-build\install`).

## 비선형 (S19/S21 레시피 고정)
- `wht_iga/analyses/nonlinear.py`: MoveMeshFlag + SetBufferSize(2) +
  AndCriteria(Residual + Displacement). S21 단독 criterion 오수렴 회피.
- S19 보완 발견 2건: (1) BC 적용 후 `CloneTimeStep` 금지 (UL 0 붕괴),
  (2) `Solve()` 반환값 대신 `IsConverged()`로 판정.
- Neo-Hookean UL/TL 단축 해석해 대비 1e-12, J2 완전소성 plateau 1e-9,
  J2 선형경화 접선기울기 E·H/(E+H) 2%, 2-스텝 ramp 최종값 해석해 1e-9.

## 접촉 (S22/S25 패턴, Python precursor)
- `wht_iga/analyses/contact.py`: tied Triangle3D3 표면층 builder
  (일치/외향/tie 검증) + `paired_gap_z` helper.
- Staggered penalty 한계 확정: k·C >> 1이면 ω와 무관하게 발산
  (ω_crit ≈ 0.05). 정답 = 반평활 Newton
  (Kratos 추출 K + nodal penalty + 일치 tangent): 3-5회 반복, 잔차 1e-12 수준.
- Hertz 무마찰 (1/4 대칭, k=1e4): NX=12 err 26% → NX=24 err 16.8% →
  NX=48 err 7.1% (48/2401 active, 5회, 855s). 잔차 = 유한 깊이 + penalty + 성김.
- Coulomb 마찰: stick predictor + radial return mapping.
  슬라이딩 probe 40스텝에서 step 30부터 F_T/F_N = 0.3000 고정 (μ=0.3).
  ux_slip = μ·k_n·depth/k_t = 0.06 일치.
- E03 축소판 (2-재료 고정단, 균일 dT): 해석해와 5.7e-16 일치 (S20 경로).

## 범위 외 (C++ 필요, 보고만)
- 다항 초탄성 (Mooney-Rivlin/Yeoh/Ogden), MixedUP0/VMS locking-free 요소,
  thermal J2 (Kratos 10.4에 법칙 없음 — J2+열연성 결합 불가로 문서화).

## 산출물
- 신규: `wht_iga/analyses/nonlinear.py`, `wht_iga/analyses/contact.py`,
  `tests/kratos_spikes/test_p5_nlgeom_contact_prep.py` (9 tests),
  `docs/spikes/S26.md`, 본 파일.
- 프레임워크 변경 1건 (가산만, `dev_log/track3_framework_changes_20261116.md` 기록):
  `wht_iga/dynamics/explicit.py`에 `coulomb_friction_update` 순수 함수 추가.
  기존 함수 무수정, Track 2 회귀 5 passed 유지.
- 금지 준수: commit/push 없음, L2/L3 벤치 없음.

## 리뷰 반영 (2026-11-17, Claude Opus)

### 수정: `wht_iga/analyses/nonlinear.py`
- **H1 이력 초기화**
  - 원인: 하중 단계마다 새 전략을 만들고 `Initialize()`를 불러 요소를 재초기화했습니다. 그래서 구성법칙 이력과 UL 기준 변형이 매번 초기화됐습니다.
  - 수정 전 실측: J2를 2 %까지 당긴 뒤 1.6 %로 풀면 σ = 0.5였습니다(정답 0.1). UL 2단계 측면 변위는 1 % 틀렸습니다.
  - 수정: `run_load_stepping`이 전략을 한 번만 만들고 `Initialize()`도 한 번만 부릅니다.
  - 새 `solve_increment()`는 수렴했을 때만 `FinalizeSolutionStep`을 호출합니다. 실패한 시도는 구성법칙 상태를 확정하지 않습니다.
- **H2 컷백 오보고**
  - 원인: 실패하면 `target × frac`, 즉 절대 목표값의 절반을 풀고 그것을 목표 도달로 셌습니다.
  - 수정: 마지막 수렴 수준 λ_c에서 λ_c + frac·(target − λ_c)로 진행합니다. 목표값 자체가 수렴할 때까지 반복하고, 그때만 도달로 셉니다.
  - 실패 후에는 변위와 이동 좌표를 함께 복원합니다(`snapshot_state` / `restore_state`).
  - 반환값에 `level`과 `n_solves`를 추가했습니다.
- **Kratos 10.4 한계 (기록만)**
  - UL 요소는 다단계 실행 후 적분점 응력 **출력**이 틀립니다. 해(변위)는 정확합니다.
  - `STEP ≤ 1`이면 마지막 증분만 반영되고, `STEP > 1`이면 F0이 두 번 곱해집니다.
  - 실측: 0.268 / 0.904 (해석해 0.587).
  - 응력 출력에는 TL을 쓰거나 변위에서 응력을 다시 계산해야 합니다.
- 새 시험: `tests/kratos_spikes/test_p5_load_stepping_regression.py`(6건, 해석해 대조).
  - 대상: J2 제하, UL·TL 2단계 ramp, 실패 주입 컷백(UL), 실패 시도 미확정(J2), 지속 실패 시 예외.
  - 변이 검출:
    - 옛 드라이버: 이력 초기화·컷백 오보고·예외 미발생 검출.
    - 확정 강제, 좌표 미복원, 복원 생략 변이 3종 모두 검출.
- 회귀(2026-11-17): `test_p5_*` + `test_e04_cook_precursor` + S19 + S21 → 55 passed.

### 미결 (리뷰 H3·M1–M4)
1. G5 기준에 아직 도달하지 않았습니다.

   | 항목 | 기준 | 현재 |
   |---|---|---|
   | V08 Hertz | 4 % | 시험 게이트 35 %, 최선 7.1 %는 probe에만 있음 |
   | V08 마찰 F_T | 경계값 검증 | 검증 없음 |
   | E03 | T 의존 J2 + NLGEOM | 탄성 축소판 |
   | E05 | 5 % | 없음 |

2. "Kratos 10.4에 thermal J2 없음"과 Python SSN 접촉 경로 채택은 사용자 결정이 필요합니다(AGENTS.md Must 6).
3. 접촉 해법 위치와 모델:
   - SSN 해법이 시험 코드에만 있습니다.
   - penalty가 절점당 상수(면적 비가중)입니다.
   - p≥2에서는 실제 표면이 아니라 control net으로 접촉을 판정합니다.
4. 마찰 함수가 `dynamics/explicit.py`(Track 2 모듈)에 있고, 시험은 구성상 성립합니다(by-construction).
