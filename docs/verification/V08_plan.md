# V08 계획 — 탄성 구형 압자 압입·미끄럼 (Hertz–Coulomb)

- 상태 (2026-11-17): **기준값 재현 완료, 해석 미실행.** 접촉 경로(Kratos mortar / Python SSN) 결정 전 단계입니다.
- 기준: `tests/verification/data/v08.json`, 원문 `research-cf-IGA/user-guide/verification-manual/problems/solid_mechanics/contact/contact-indenter/contact-indenter.md` (§2–§7, §8.3 Abaqus 입력)
- 코드: `wht_iga/reference/hertz.py` (Hertz·정확한 구 해·V08 이력·격자·압자 삼각형 soup)
- 시험: `tests/unit/test_hertz_reference.py` → **29 passed (5.6 s)**, Kratos 무의존. 변이 11종(p0 계수, E* 식, 정확해 F, χ₂ 부호, 외피 법선, 링 정점 용접, 간극 차감, RMS 정의, 대칭면 법선, Dundurs, nnz 띠폭)을 모두 잡습니다.
- 게이트 G5: F_N(t = 0.1) 4 % 이내, 미끄럼 F_T RMS ≤ 4.5 %. p_max는 기록·비교만 합니다.

## 1. 기준값 재현

입력(IPS: in, lbf, psi)은 json에서 읽습니다. 블록 E = 2.9e7, ν = 0.3. 압자 E = 3.0e7, ν = 0.29. R = 0.5, δ = 0.1, μ = 0.3(`boundary_conditions.contact` 문자열에서 추출).

가정(원문 §6):
- 두 물체를 모두 선형 탄성 반무한체로 보고 `1/E* = (1−ν_b²)/E_b + (1−ν_i²)/E_i`로 합칩니다.
- 압자 곡률 반경은 외경 R이고, 블록 윗면은 평면(R₂ = ∞)입니다.
- δ는 절단면 강제 변위 0.125에서 초기 간극 g₀ = 0.025를 뺀 값입니다.
- 마찰은 수직 해에 영향을 주지 않습니다. Dundurs β = −1.04e-3로 거의 0이기 때문입니다.

| 양 | 계산값 | v08.json | 상대차 |
|---|---:|---:|---:|
| E* (psi) | 16,152,659 | 16,152,660 | 원문 Table 3 1.615266e7 대비 −4.1e-8 |
| a (in) | 0.2236068 | 0.223607 | −9.0e-7 |
| F_N (lbf) | 481,579.26 | 481,579.3 | −8.9e-8 |
| F_T = μF_N (lbf) | 144,473.78 | 144,473.8 | −1.6e-7 |
| p₀ (psi) | 4,598,743.2 | 4,598,743.0 | +3.4e-8 |
| p_mean = ⅔p₀ (psi) | 3,065,828.8 | — | — |

- 시험 허용차는 5e-6입니다. json이 6–7자리로 반올림돼 있어 a 하나만 해도 최대 2.2e-6 차이가 납니다.
- **json 결함 2건 수정 완료**:
  1. `reduced_modulus_E_star_psi`를 두 곳 모두 ksi 절삭값 `16152.0`에서 원문 Table 3의 `16152660.0` psi로 정정했습니다.
  2. `block_Lx_Ly_Lz_in`을 `[5, 2, 1]`에서 원문 Table 1 및 Abaqus `*Node`와 일치하는 `[5, 2, 3]`으로 정정했습니다(`V08_BLOCK_LO/HI`와 일치).
  - 시험 `test_v08_json_E_star_and_block_dimensions`가 두 항목(`1.615266e7` psi, `[5, 2, 3]`)을 직접 검증합니다.

## 2. Hertz가 빼놓는 것 (모형 오차)

json 기준값은 Hertz 식 그대로라 정확히 재현됩니다. 문제는 V08 문제 자체가 Hertz 가정을 벗어난다는 점입니다.

| 항목 | 크기 | 근거 |
|---|---|---|
| a/R | 0.447 (Hertz 가정 a ≪ R 위반) | `v08_model_form_gaps` |
| 정확한 구 간극(포물면 아님), 같은 반무한체·같은 δ = 0.1 | a −3.34 %, **F_N −2.00 %**, p₀ +1.14 % | Segedin/Sneddon 닫힌식. 시험에서 Sneddon 적분을 수치로 다시 계산해 1e-9로 일치 |
| 블록 바닥 고정(깊이 2 in) | Boussinesq 축 변위가 δ의 4.1 % | 크기 추정만. 측면이 자유이고 블록이 유한해 부호와 크기 모두 불확실 |
| 압자 하중면(접촉점에서 0.5 in = 2.2a) | Boussinesq 값으로 δ의 16 % | 압자는 반무한체가 아니라 중공 반구(벽 두께 0.25 = 1.1a)라 이 추정은 쓸 수 없음. 닫힌 해 없음 |
| NLGEOM (δ/R = 0.2) | 미정 | 원문은 기하 비선형 + 소변형 탄성 |

- 결론: 물체가 반무한체였다고 해도, 이산화 오차가 0인 해는 Hertz보다 약 2 % 낮습니다. 4 % 게이트의 절반입니다.
- 유한 블록·중공 압자 효과는 닫힌 해로 계산할 수 없습니다. 원문 Coreform IGA(1.84–3.74 %)와 Abaqus h₂(1.07 %)를 보면 모형 오차와 이산화 오차를 합친 차이가 수 % 수준입니다.
- 권고: F_N은 Hertz와 정확한 구 해 **둘 다**와 비교해 기록합니다. 게이트 판정은 원문대로 Hertz 기준입니다.

## 3. 케이스 정의

공통 조건:
- 블록 바닥 y = −1은 u = 0(Abaqus pinned)입니다.
- 압자 절단면 y = y_c에 강제 변위를 줍니다. `u_y = −0.125·min(t/0.1, 1)`, `u_x = 4.25·max(0, (t−0.1)/0.9)`, `u_z = 0` (`v08_cut_face_displacement`).
- 접촉: 압자 외구(soup face 1)가 slave, 블록 윗면 y = 1이 master입니다. Abaqus `*CONTACT PAIR` 순서와 같습니다.
- 출력:
  - F_N = −ΣRF_y, F_T = |ΣRF_x|. 절단면 제어점 반력의 합입니다.
  - p_max
  - 상대 RMS(원문 §7.1, `relative_rms`)
- 대칭 모델은 힘을 2배(half) 또는 4배(quarter)로 환산합니다.

| 케이스 | 모델 | 물성 | 이력 | 판정 |
|---|---|---|---|---|
| **R0** 축소·압입 | 1/4 블록(x ≥ −1.5, z ≥ 0), **강체 해석 구** | 블록 E′ = (1−ν²)E* = 1.46989e7, ν = 0.3 → E* 동일 | δ 0 → 0.1, 8–10 증분, μ = 0 | F_N(δ) 곡선 vs Hertz·정확한 구. 4 % 이내 |
| **R1** 축소·압입 | 전체 모델(또는 §4 버그 수정 후 반모델), 두 물체 모두 탄성 | 원문 그대로 | t 0 → 0.1, Δt = 0.01, μ = 0 (β ≈ 0이므로 F_N은 μ와 무관) | F_N(0.1) 4 %, p₀ 기록 |
| **R2** 축소·미끄럼 | R1과 같음 | 원문 그대로, μ = 0.3 | t 0 → 0.40 (x_c = −0.083), 압입 10 + 미끄럼 26 증분 | F_T/F_N → μ, RMS 창 0.15 ≤ t ≤ 0.40 |
| **F** 전체 V08 | 전체(또는 버그 수정 후 z = 0 반모델, 전 이력에서 정확) | 원문 그대로 | t 0 → 1, 미끄럼 Δx ≤ a/4 = 0.056 → 약 76 증분 | G5: F_N(0.1) 4 %, F_T RMS 4.5 % (창 0.15–0.70) |

- R0는 블록이 body-fitted라 trimming이 없습니다. 따라서 §4 버그와 무관합니다. 강체 압자로 바꾸면 Dundurs β가 0.286이 되므로 μ = 0으로만 씁니다.
- 압자 격자 규칙(`v08_indenter_grid`):
  - 셀 크기는 정확히 h입니다. 윗평면은 절단면 y = y_c입니다.
  - 그 평면의 제어점 전부에 강제 변위를 주면, 경계 trace가 경계 제어점만으로 정해지므로 환형 하중면 BC가 정확히 들어갑니다(구멍 위쪽은 가상 영역).
  - 나머지 면은 구와 0.1h 이상 떨어뜨려, 접촉점이 격자면 위에 오지 않게 합니다.
  - 대칭면은 격자의 아래쪽 면에 둡니다(V06 octant 관례).
- 압자 soup(`hollow_hemisphere_soup`):
  - 극을 절단면 테두리(z = ±r)에 두었습니다. 그래서 접촉점은 규칙 격자 정점이고, 2a 이내 삼각형의 최소각은 n = 12, 24, 48 모두 30° 이상입니다. 접촉점에 극을 두면 7.5°(n = 24)가 됩니다.
  - 정점은 정확한 곡면 위에 있고, 면 사이 정점은 비트 단위로 같아 위치로 용접됩니다. 닫혀 있고 바깥 방향입니다.
  - 체적 결손은 n = 8, 16, 32에서 6.3 / 1.6 / 0.40 %로, O(n⁻²)입니다.

## 4. 격자와 메모리 (p = 2, C¹, Coreform과 같음)

| h | 블록 셀 (full / half / quarter) | 압자 셀 (full / half / quarter) | 접촉 지름당 셀 2a/h |
|---|---|---|---:|
| 0.5 | 10×4×6 / 10×4×3 / 8×4×3 | 3×2×3 / 3×2×2 / 2×2×2 | 0.89 |
| 0.25 | 20×8×12 / 20×8×6 / 16×8×6 | 5×3×5 / 5×3×3 / 3×3×3 | 1.79 |
| 0.125 | 40×16×24 / 40×16×12 / 32×16×12 | 9×5×9 / 9×5×5 / 5×5×5 | 3.58 |

전체 격자 DOF 수입니다(비활성 제어점 포함, 상한). 계수 메모리는 Pardiso 기호 분석(phase 11)의 예측값입니다. 같은 희소 패턴의 합성 SPD 행렬로 실측했고 수치 분해는 하지 않았습니다(스크래치 `T3/probe_fill.py`).

| h | 모델 | DOF (블록 + 압자) | nnz(K) | LU 계수 nnz | LU 메모리 | LDLᵀ 메모리 |
|---|---|---:|---:|---:|---:|---:|
| 0.5 | half | 1,080 + 240 | 0.26 M | 0.56 M | 25 MB | 22 MB |
| 0.5 | full | 1,728 + 300 | 0.44 M | 1.18 M | 32 MB | 26 MB |
| 0.25 | half | 5,280 + 525 | 1.49 M | 6.2 M | 84 MB | 55 MB |
| 0.25 | full | 9,240 + 735 | 2.78 M | 15.5 M | 168 MB | 99 MB |
| 0.125 | half | 31,752 + 1,617 | 10.2 M | 94 M | 846 MB | 450 MB |
| 0.125 | full | 58,968 + 2,541 | 19.8 M | 245 M | 2,084 MB | 1,080 MB |

- R0 quarter 블록 DOF: 900 / 4,320 / 25,704.
- 이 PC 기준(31.6 GB, 측정 시 여유 16 GB, 에이전트 3개 병행, 실행당 2 GB 미만 목표):
  - h₀·h₁: 전부 가능합니다.
  - h₂ half: 버그 수정 후 가능합니다. LU 0.85 GB + K·Kratos 모델.
  - h₂ full: 계수만 2.1 GB입니다. PC가 비어 있을 때만 돌립니다.
  - 마찰이 있으면 접선 행렬이 비대칭이라 LU입니다.
- 저장소 추정식 `est_pardiso_bytes`는 0.41 / 3.8 / 42.5 GB(full)를 냅니다. Pardiso 예측의 약 20배입니다. `auto`가 h₂를 불필요하게 AMGCL로 넘길 수 있습니다.

### 해결: `classify` 셀 순서 버그 (`wht_iga/geometry/classify.py`)

- 원인(수정 완료):
  - 기존 `_find_cut_cells_njit`가 SAT 절단 플래그를 `i + nx·(j + ny·k)` 순서(i가 가장 빠름)로 기록한 반면, `classify_cells`의 꼭짓점 판정(`allin.ravel()`)과 `trimmed._build_trimmed_table`(`state.reshape(n)[i, j, k]`)은 C 순서 `(i·ny + j)·nz + k`로 읽어 절단 플래그가 전치된 셀에 매핑되는 문제가 있었습니다.
- 수정 및 검증:
  - `_find_cut_cells_njit`의 기록식을 `hits[pos] = (i * ny + j) * nz + k`로 C 순서에 맞춰 수정했습니다.
  - `T3/repro_classify_order.py`(상자 `[0.01, 0.24] × [0.01, 0.49] × [0.01, 0.74]`, 4×2×3 및 3×3×3 격자)에서 `V_quad == V_exact = 0.080592` 일치를 확인했습니다.
  - 회귀 시험 `tests/unit/test_geo_classify_order.py`에서 비정방 격자(`(4, 2, 3)`, `(2, 3, 5)`) SAT-only 상자 인덱스·체적과 비대칭 중공 반구(`sym="half"`) 및 편심 상자 체적을 검증합니다.
  - 이제 V08 half 압자 격자(h = 0.5, 0.25, 0.125) 등 비대칭·비정방 trimmed 모델을 정상 사용할 수 있습니다.

## 5. 접촉 경로 두 가지와 남은 증명

### A. Kratos penalty mortar + S29 표면층
- 증명된 것 (S29):
  - 우회 2건으로 IGA 블록과 tied `Triangle3D3` 표면층에서 초기화됩니다.
  - p = 1 블록 2개 압축이 일체 해석해로 O(1/penalty) 수렴합니다. 실효 강성은 약 0.148·`INITIAL_PENALTY`입니다.
- 남은 것:
  1. **분리(단방향)**: 미해결입니다. 당기면 인장을 전달하고, contact strategy를 쓰면 부호 대칭인 오답이 나옵니다. V08 미끄럼 후단(앞쪽 분리)과 R0 하중 감소에 필수입니다.
  2. Coulomb 마찰: frictional mortar 조건을 같은 우회로 쓸 수 있는지 시험한 적이 없습니다.
  3. 곡면 slave:
     - S29 표면층은 경계 제어점에 1:1로 묶였습니다(p = 1 평면).
     - trimmed 압자 외구에는 제어점이 없습니다. 표면 정점을 `u = Σ N_i(ξ) u_i` 다점 구속(MPC)으로 묶는 연결이 필요합니다.
     - p = 2 블록 윗면도 control net 보간과 B-spline trace가 다릅니다.
  4. 유한 미끄럼 4.25 in: 매 증분 mortar 짝 재탐색, NLGEOM(TL) 요소와 결합.
  5. penalty 오차 예산: O(1/k) 오차를 F_N의 0.5 % 미만으로 맞출 k와 조건수의 균형.

### B. Python 반평활 Newton (Track 3, S26)
- 증명된 것:
  - 강체 해석 구 + p = 1 body-fitted 블록 무마찰 접촉이 3–5회 반복으로 수렴합니다.
  - F_N 오차는 NX = 12 / 24 / 48에서 26 / 16.8 / 7.1 %입니다(1/4, k = 1e4, 깊이 2). penalty·깊이·이산화 기여는 분리되지 않았습니다.
  - 마찰 plateau F_T/F_N = 0.3000(장난감 모델, return map).
  - `run_load_stepping`: 전략 1개, 수렴했을 때만 확정.
- 남은 것:
  1. SSN 해법이 아직 시험 코드 안에만 있습니다. 모듈로 옮겨야 합니다(Track 3 미결 3).
  2. penalty를 절점 상수에서 표면 적분점 면적 가중으로 바꿔야 합니다. p ≥ 2에서는 control net이 아니라 실제 B-spline 표면으로 간극을 계산해야 합니다.
  3. 변형체–변형체:
     - slave 적분점은 trimmed 압자 외구에 둡니다(soup face 1 → 표면 적분).
     - master는 블록 윗면(untrimmed B-spline 면)이고, 변형 후 최근접점 투영이 필요합니다.
  4. Coulomb return map의 일치 접선(비대칭)과 LU 경로.
  5. NLGEOM: 반복마다 Kratos TL에서 K를 꺼내고 접촉 기여를 더해야 합니다(UL 응력 출력 오류 회피).
  6. R0로 7.1 % 잔차를 분해해야 합니다: penalty 수렴 → 정확한 구 해(−2.0 %) 대비 → 격자 수렴.
- 장점: 활성 집합으로 분리가 기본 지원됩니다(경로 A의 미해결 1번이 없음).
- 사용자 결정(AGENTS.md Must 6)이 필요합니다. TRACK3_COMPLETE 미결 2와 같습니다.

## 6. 다음 단계
1. 완료 / 결정:
   - `classify` 셀 순서 수정(geometry) — **완료** (`tests/unit/test_geo_classify_order.py`)
   - v08.json의 E* 단위와 Lz 정정 — **완료** (`tests/unit/test_hertz_reference.py`)
   - 접촉 경로: **B. Python SSN 주력 진행** (`wht_iga/analyses/contact.py`, `docs/spikes/S31.md`); Kratos C++ 마찰 포크는 심층 리팩토링 규모로 보류
2. R0(경로 B, 강체 구, quarter body-fitted) — **h₀–h₂ 검증 완료** (`test_p5_v08_r0_ssn.py`).
3. R1 → R2 → F 순서로 진행합니다(§4 `classify` 수정 완료로 half 모델 사용 가능).

## 재현
- 시험: `D:\mamba\envs\whtiga\python.exe -m pytest -q tests/unit/test_hertz_reference.py` (`WHT_IGA`에서 실행)
- 수치 조사 스크립트(저장소 밖, 스크래치 `T3/`):
  - `probe1.py`: 기준값·모형 오차
  - `probe_fill.py`: Pardiso 기호 분석
  - `probe_soup.py`: soup·trimmed 체적
  - `repro_classify_order.py`: 버그 재현
- 환경: Python 3.12 (`whtiga`), Kratos 미사용, 2026-11-17.
