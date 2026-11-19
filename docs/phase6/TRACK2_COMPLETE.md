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

## V-EXP / V-EXP-C
- V-EXP: p=1/2/3 tip L2 < 2 %. 피크 응력과 파면 도달은 아래 미결 1·2를 보세요
  (2026-11-19 갱신; 이전 "도달 0.5 s ± 0.03"은 6 % 허용치에 Greville 프로브였습니다).
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

### 미결 (계획서 §10.1과 차이, 리뷰 H1·M1–M5; 1–3은 2026-11-19 갱신)
1. V-EXP 피크 응력: **시험 추가, "피크" 해석은 사용자 승인 필요.**
   - `test_vexp_peak_stress`: x = 0 장 응력 E u'(0)을 t ∈ [1.5, 2.0] s에서 평균해 해석해 2σ0와 비교합니다(2 %).
   - 실측 −0.20 / +0.27 / −0.38 % (p = 1 / 2 / 3, nx = 40).
   - 원시 최댓값은 게이트하지 않습니다. 실측 +25.5 / +25.4 / +24.8 %입니다.
     - 원인: 분산 파면의 Airy 과도(극한 2σ0 × 1.2744, +27.4 %)이고, 솔버 결함이 아닙니다.
     - 세분해도 줄지 않고 극한 쪽으로 커집니다(p = 1: nx = 40 / 80 / 160에서 +25.5 / +26.2 / +26.7 %).
   - 결정 필요: 계획서의 "피크 응력 2 %"를 평탄부 값으로 읽어도 되는지.
2. V-EXP 파면 도달 2 %: **원인 규명 완료(수치 분산), 기준 정의는 사용자 결정 필요.**
   - 프로브는 x = 5.0 장 값 Σ N_i(5) u_i로 고쳤습니다.
   - 도달 시각 실측(nx = 40):

     | 정의 | p = 1 | p = 2 | p = 3 | 2 % 게이트 |
     |---|---:|---:|---:|---|
     | 1차 모멘트(u(5, t) 직선 절편) | +0.04 % | +0.22 % | +0.43 % | 통과 |
     | 임계값 0.5 v0 (`vexp.json` 정의) | +2.85 % | +3.89 % | +4.25 % | **미달** |

   - 원인: lumped mass + 중앙차분의 수치 분산입니다. 계단 파면이 Airy 적분 모양으로 퍼져 50 % 지점이 늦게 옵니다. 시험 결함이 아닌 근거:
     - Kratos와 독립 1D 기준(SciPy B-spline, 같은 이산계)이 반올림 수준으로 같습니다(|Δt| ≤ 6e-16 s).
     - Airy 예측 지연 a (3β)^(1/3) (h/d)^(2/3)와 맞습니다(실측/예측 1.04 / 0.98 / 0.97).
       a = 0.4096, β = (p + 1 − CFL²)/24 (p = 1은 (1 − CFL²)/24), d = 5.
     - h^(2/3)로 수렴합니다. nx = 40 → 80 → 160의 관측 차수는 0.64–0.69입니다.
     - 1차 모멘트는 분산과 무관하게 c로 움직입니다(∫ s Ai(s) ds = 0).
   - 세분 격자에서는 임계값 정의도 2 %를 만족합니다(Kratos, CFL 0.2): p = 1 nx = 80 +1.77 %, p = 2 nx = 160 +1.57 %, p = 3 nx = 160 +1.73 %.
     격자는 Airy 예측으로 미리 정했습니다.
   - 시험(skip/xfail 없음):
     - `test_vexp_wavefront_arrival_gap_pin`: 게이트가 아닙니다. 간극을 고정해 두고, 간극이 2 % 밑으로 줄면 실패해서 게이트 복원을 알립니다.
     - `..._lag_converges`, `..._refined_mesh`, `..._first_moment`
   - 선택지:
     - (a) `vexp.json` 격자를 nx = 160으로 올립니다(p = 1은 80).
     - (b) nx = 40을 유지하고 1차 모멘트 정의를 게이트로 씁니다.
     - (c) 허용치를 완화합니다. 계획서 2 %와 어긋나므로 권하지 않습니다.
   - `vexp.json`의 acceptance는 아직 0.03 s(6 %)이고 계획서는 2 %입니다. 시험은 계획서 값을 씁니다. 데이터 파일은 이 트랙에서 고치지 않았습니다.
3. CFL·mass scaling 실행 시험: **해소.** `test_explicit_cfl_mass_scaling.py`, 23개.
   - `estimate_dt_crit`을 Kratos가 실제로 적분하는 계에 적용합니다.
     - K: 요소 LHS(정적 builder는 0에 가까운 대각을 바꿉니다. p = 3, η = 1e-6에서 3.7e-30 → 19.0).
     - M: Kratos `NODAL_MASS` 또는 덮어쓴 HRZ. 확장 시 TᵀKT를 씁니다.
   - p = 1 비절단 dt_crit은 1D 고정-자유 lumped 막대 해석해 (h/c)/cos(π/16) = 0.1019591158과 상대 1e-9 이내로 같습니다.
   - 실행 bracket 12건(p 1–3 × Kratos 질량·HRZ·절단+T η = 1e-4·1e-6):
     - 0.95 dt_crit에서 max E/E0 = 1.27–1.65입니다(leapfrog 상한 2).
     - 1.05·1.5 dt_crit에서는 100스텝 안에 ≥ 1e49로 발산합니다.
   - S24 주장 "절단 셀 dt_crit ≥ 0.5 × 비절단"을 측정했습니다.
     - 절단+T / 비절단 = 1.015 / 1.110 / 1.032 (p = 1 / 2 / 3, η = 1e-4·1e-6 같음).
     - 0.525 dt_U에서 안정합니다.
     - 확장이 없으면 (1.2–1.4) η 배로 붕괴합니다.
   - mass floor(η = 1e-2, 하한 = 비절단 최소 HRZ 질량의 1 %):
     - dt_F / dt_U = 1.018 / 1.118 / 1.037, 추가 질량 ≤ 0.17 %.
     - floor를 빼면 0.95 dt_F에서 발산합니다.
   - Kratos는 `NODAL_MASS` ≤ eps인 자유도의 가속도를 0으로 둡니다.
     - 확장 없는 η = 1e-4, p = 2(최소 질량 2.3e-22)는 dt_crit의 약 4000배에서도 유한합니다.
     - 따라서 유계 실행만으로는 dt_crit 복원의 증거가 되지 않습니다. 이 동작은 시험으로 고정했습니다.
4. V-EXP-C의 차이: **해소 (2026-11-19).**
   - `test_vexp_explicit.py::test_vexp_c_dt_sweep`를 `p = 1`과 `p = 2` 모두로 확장했습니다(`p = 2`는 `1x1x4` 블록, `3x3` Gauss 접촉점 `add_linear_tie` 결합, `k_total = 2000`).
   - 전 스텝에서 시간 정렬(`t_n`) 전체 에너지(`KE + PE_strain + E_penalty`)를 판정해 `≤ 2%` 게이트를 통과했습니다(`p = 1`: `0.3088%` / `1.2345%`, `p = 2`: `0.3549%` / `1.5000%` at `dt = 0.0025` / `0.005`).
   - `test_s25_explicit_contact.py::test_s25_direct_penalty_free_impact`에도 전 스텝 전체 에너지(`KE + PE_strain + E_penalty`) 검사를 추가했습니다(실측 `max_dE = 7.9645% < 0.10`, `_KPEN = 2000` / `dt = 0.005`에서의 기존 문서화된 +8% 접촉 스프링 분해능 성장과 일치).
5. 기준 데이터와 시험 연결이 일부만 되어 있습니다.
   - `test_vexp_explicit.py`는 이제 `vexp.json`을 읽습니다(모델·재료·하중·dt·t_end·도달 기준값). 빌더 상수와 같은지도 검사합니다.
   - `vexp.json`의 acceptance(0.03 s)는 갱신해야 합니다(미결 2).
   - `vexp_c.json`은 아직 어떤 시험도 읽지 않습니다.
   - `vexp_c.json`은 `degree = 2, penalty_stiffness_per_pair = 2000.0`이고, V-EXP-C 시험(`p = 1, 2`)은 에너지 2% 게이트를 만족하기 위해 `k_total = 2000.0`을 사용합니다(`p = 1` 4쌍 `k_pen = 500`, `p = 2` 3x3 Gauss 가중합 `2000`).
   - `title` 누락으로 `test_verification_core`가 실패하던 문제는 2026-11-17에 해소했습니다(`id`·`title`·`status` 추가, 값 무수정). `vexp_c.json`은 `status: draft`이고 불일치는 `review_note`에 적었습니다.
6. 계획서와 다른 경로 2건은 사용자 승인과 계획서 반영이 필요합니다.
   - 접촉: ContactSM 프로세스 대신 직접 penalty.
   - 질량: Kratos lumped 대신 HRZ 덮어쓰기.
7. `dynamics/explicit.py`의 문제:
   - `build_x_extension`이 dense 행렬이고 +x 끝 전용입니다.
   - x 방향 셀이 1개면 오류 없이 잘못된 T를 만듭니다.
   - `run_explicit_history`는 시간을 0에서 다시 시작합니다.

## 2026-11-19 검증 보강 (미결 1–4 해소 요약)
- 전체 회귀 명령:
  `D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests/kratos_spikes/test_vexp_explicit.py tests/kratos_spikes/test_explicit_cfl_mass_scaling.py tests/kratos_spikes/test_s24_explicit_dynamics.py tests/kratos_spikes/test_s24b_hybrid_explicit.py tests/kratos_spikes/test_s25_explicit_contact.py tests/unit/test_explicit_alpha.py`
  → **64 passed in 570.56s (0:09:30)** (0 failed).
  - 개별 시험 수: `test_vexp_explicit.py` **24**, `test_explicit_cfl_mass_scaling.py` **23**, `test_s24_explicit_dynamics.py` **6**, `test_s24b_hybrid_explicit.py` **4**, `test_s25_explicit_contact.py` **2**, `test_explicit_alpha.py` **5**.
- **항목 1 (V-EXP 피크 응력 ≤ 2%)**: `x = 0` 장 응력 `E sum N_i'(0) u_i` (`t ∈ [1.5, 2.0]` s 평균) vs `2σ0 = 2.0`:
  - `p = 1`: `1.99601` (`-0.20%`, raw max `2.5110` `+25.55%`)
  - `p = 2`: `2.00550` (`+0.27%`, raw max `2.5073` `+25.36%`)
  - `p = 3`: `1.99238` (`-0.38%`, raw max `2.4964` `+24.82%`)
- **항목 2 (V-EXP 파면 도달 `x = 5.0` 장 프로브 `Σ N_i(5) u_i`, 해석해 `0.5 s`, 게이트 ≤ 2%)**:
  - 1차 모멘트 (`nx = 40`): `p = 1` `0.5002 s` (`+0.04%`, slope `+0.02%`), `p = 2` `0.5011 s` (`+0.22%`, slope `+0.17%`), `p = 3` `0.5021 s` (`+0.43%`, slope `+0.37%`).
  - `0.5 v0` 임계값 (`nx = 40`): `p = 1` `0.514237 s` (`+2.85%`), `p = 2` `0.519456 s` (`+3.89%`), `p = 3` `0.521270 s` (`+4.25%`) — Kratos vs 1D SciPy 기준 차이 `≤ 4.4e-16 s`, Airy 예측 대비 비율 `1.038 / 0.975 / 0.967`, `h^(2/3)` 수렴 차수 `0.644–0.685`.
  - `0.5 v0` 임계값 (Airy 예측 세분 격자, CFL 0.2): `p = 1 (nx = 80)` `0.50886 s` (`+1.77%`), `p = 2 (nx = 160)` `0.50786 s` (`+1.57%`), `p = 3 (nx = 160)` `0.50864 s` (`+1.73%`) — 모두 `≤ 2%` 통과.
- **항목 3 (CFL·mass scaling 실행 검증)**: `test_explicit_cfl_mass_scaling.py` 23/23 통과.
  - `p = 1` 비절단 `dt_crit = 0.1019591158` (해석해와 일치), 실행 bracket 12건 (`0.95 dt_crit`에서 `max E/E0 = 1.266–1.645 ≤ 2`, `1.05 dt_crit`에서 `≥ 2.47e+49` 발산).
  - 절단+T (`η = 1e-4, 1e-6`) `dt_ext / dt_U = 1.0154 / 1.1105 / 1.0322` (`p = 1 / 2 / 3`), 미확장 시 `~1.21e-4 / 1.21e-6` 배로 붕괴; mass floor (`η = 1e-2`) `dt_F / dt_U = 1.0176 / 1.1182 / 1.0373`, 추가 질량 `5.74e-4 – 1.66e-3`.
- **항목 4 (V-EXP-C `p = 1, 2` 충돌 및 전체 에너지 `KE + PE_strain + E_penalty`)**:
  - `p = 1` (`1x1x1`, 4쌍 `k_pen = 500`):
    - `dt = 0.0025`: `mom_dev = 2.33e-15`, `min_gap = -0.0113`, `v_up = 0.0005`, `v_lo = -1.0005`, `max_dE = 0.3088%`, `end_dE = 0.1922%` (`< 2%`).
    - `dt = 0.005`: `mom_dev = 4.44e-16`, `min_gap = -0.0112`, `v_up = 0.0009`, `v_lo = -1.0009`, `max_dE = 1.2345%`, `end_dE = 0.2757%` (`< 2%`).
  - `p = 2` (`1x1x4`, 3x3 Gauss 접촉점 `k_total = 2000`):
    - `dt = 0.0025`: `mom_dev = 2.22e-16`, `min_gap = -0.0070`, `v_up = -0.0743`, `v_lo = -0.9257`, `max_dE = 0.3549%`, `end_dE = 0.0528%` (`< 2%`).
    - `dt = 0.005`: `mom_dev = 4.44e-16`, `min_gap = -0.0077`, `v_up = -0.0781`, `v_lo = -0.9219`, `max_dE = 1.5000%`, `end_dE = 0.2090%` (`< 2%`).
  - `S25` (`test_s25_direct_penalty_free_impact`, `k_pen = 2000`, `dt = 0.005`): `mom_dev = 4.44e-16`, `min_gap = -0.0060`, `max_dE = 7.9645% < 10%` (전체 에너지 검사 추가).
- **미결 유지**: 항목 5(`vexp_c.json` 시험 미연결 및 `k_pen` 차이), 항목 6(계획서 대비 경로 승인), 항목 7(`dynamics/explicit.py` 제한사항) 및 항목 1–2의 기준 정의에 대한 사용자 최종 결정은 미결로 유지하므로 `STATUS: COMPLETE_WITH_OPEN_ITEMS`를 유지합니다.
