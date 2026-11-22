# S48 다항식 초탄성(Mooney-Rivlin·Yeoh)을 Numba 평균 팽창 F̄/P0로 구현 (2026-11-20)

## 질문
master plan Phase 5 재료 사다리의 다항식 초탄성(Mooney-Rivlin, Yeoh, Ogden)은 C++를 가정했다.
Kratos `ConstitutiveModelsApplication`이 Python에 노출되지 않고 `QuadraturePointGeometry::Volume()`이 없기 때문이다(S19).
사용자 승인에 따라 S40 F̄ Neo-Hookean(`fbar.py`)과 같은 Numba 경로로 만들 수 있는지 확인했다.

## 결과: 가능. `wht_iga/analyses/hyperelastic_fbar.py`
- 재료(Abaqus 규약, 압축성, 등체적/체적 분리):
  - W = W_iso(Ī1, Ī2) + U(J), Ī1 = J^(−2/3) I1, Ī2 = J^(−4/3) I2
  - W_iso = Σ_{1≤i+j≤3} C_ij (Ī1−3)^i (Ī2−3)^j, U = Σ_{k=1..3} (1/D_k)(J−1)^(2k)
  - `mooney_rivlin(C10, C01, D1)`, `yeoh(C10, C20, C30, D1, D2, D3)`, `neo_hookean(C10, D1)`(= MR, C01 = 0),
    `polynomial({(i, j): C_ij}, (D1, D2, D3))`(Abaqus POLYNOMIAL N ≤ 3). μ0 = 2(C10 + C01), K0 = 2/D1.
  - D1은 유한 양수만 허용한다(hybrid·완전 비압축 미지원, 작은 D1로 근사). μ0 ≤ 0이면 거부한다.
- 정식화(Simo–Taylor–Pister 평균 팽창, S40과 같은 구조):
  - Π_cell = Σ_q w_q W_iso(F_q) + V_c U(J̄), J̄ = (1/V_c) Σ_q w_q J_q
  - 등체적 부분은 점별 적분한다. Ī1, Ī2는 F̄ 치환 F → (J̄/J)^(1/3) F에 불변이다.
  - 체적 부분은 셀 J̄만 본다. 압력 P0 = U′(J̄)는 닫힌 식으로 응축한다.
  - R_I = Σ_q w_q P_iso(F_q) G_qI + U′(J̄) V_I, V_I = Σ_q w_q J_q F_q^−T G_qI
  - K = Σ_q w_q G_qI·A_iso·G_qJ + U″(J̄) V_I⊗V_J / V_c + U′(J̄) ∂V_I/∂u_J. Π의 정확한 헤시안이라 대칭이다.
  - A_iso = ∂²W_iso/∂F∂F는 Ī1·Ī2의 1·2차 도함수 닫힌 식으로 만든다(모듈 docstring `_iso_PA`).
- 응력 출력(QP): σ_q = P_iso(F_q) F_q^T / J_q + U′(J̄) I. 약형식이 평형시키는 응력이다(`fbar.cauchy_stress`와 같은 규약).
  균일 변형에서는 J̄ = J라 Abaqus 초탄성 Cauchy 응력과 같다.
- 기하·적분점·CSR 패턴은 `fbar.build_model`을 import해서 쓴다(`fbar.py` 무수정). CSR 산포는 `fbar._assemble_csr`와 같고 셀 커널만 다르다.
- `solve_static`은 `fbar.solve_static`의 적응 Newton을 복사·수정했다. 바뀐 점은 처방 변위의 선형화 Dirichlet 예측자 하나다:
  K_ff Δu_f = −(R_f − f_f) − K_fp Δu_p. 잔차 게이트는 이 예측자 우변 기준이다. `predictor=False`이면 fbar 방식으로 시작한다.
- 범위: 비절단 단일 패치 박스, 매듭 연속성 C^(p−1)..C^0. Kratos 의존 없음.

## 검증 (`tests/unit/test_hyperelastic_fbar.py`, 31 passed)
모든 기준값은 시험 파일 안에서 독립적으로 계산한다. W(F)는 Abaqus 식을 그대로 다시 코딩했다.
응력은 그 W의 complex-step 미분(P = Im W(F + ihE)/h, 절단오차 없음)으로 구한 뒤 σ = P F^T / J로 바꾼다.

| 시험 | 기준 | 측정 최대 상대오차 | 게이트(근거) |
|---|---|---:|---:|
| 잔차 = 평균 팽창 에너지 Π의 기울기 (MR·Yeoh·N=3 다항식 × Cook C^1 2×2×1·박스 C^0 2×1×2, 비균일 상태 J̄ 1.28–1.36) | 시험에서 조립한 Π의 complex-step 기울기 | 1.6e-15 – 2.8e-15 | 1e-11 (양쪽 모두 반올림 수준이며, QP 수백 항 누적의 1000배 여유) |
| 접선 vs 중앙 FD (같은 6경우) | FD(R), h = 1e-6·max\|u\| | 8.2e-11 – 6.0e-10 | 1e-6 (요구 기준. h = 1e-5..1e-7에서 1e-10–6e-9) |
| 접선 대칭 | K^T | ≤ 2.4e-16 | 1e-12 |
| 아핀 패치 (MR·Yeoh·N=3, Cook 사다리꼴 사상, 일반 비대칭 F, det F ≈ 1.30) | complex-step σ(F) · 내부 제어점 잔차 0 | σ 2.8e-15 – 3.1e-15, 내부 잔차/max\|R\| ≤ 2.0e-16 | 1e-12 (둘 다 반올림 수준이어야 함) |
| QP 응력 출력 = 잔차가 평형시키는 응력 (MR·Yeoh·N=3, 비균일) | Σ w J σ F^−T G를 시험에서 재조립 | ≤ 5.0e-16 | 1e-12 (같은 산술) |
| 단축 λ = 1.6 / 0.7 (MR, Yeoh), 측면 무견인 | 이분법 측면 늘림 + complex-step σ, 반력/A0 공칭응력 | QP 3.0e-15 – 6.5e-15, 반력 ≤ 2.2e-15 | 1e-9 (dtol 1e-13 × K0/μ0 ≤ 5, 100배 여유) |
| 등이축 λ = 1.4 (MR, Yeoh) | 이분법 두께 늘림 + complex-step σ | 4.0e-15 / 5.1e-15 | 1e-9 |
| 단순전단 γ = 0.6 (MR, Yeoh), 2×2×2 내부 Newton | complex-step σ 전 성분(Poynting 수직응력 포함) + 손 공식 σ12 = 2(W1+W2)γ, σ11−σ22 = 2(W1+W2)γ² | 4.2e-15 / 5.1e-15 | 1e-9 |
| 비압축 극한 MR (C10 = 0.4, C01 = 0.1, λ = 1.5) | 2(λ − 1/λ²)(C10 + C01/λ) = 0.985185 | K0/μ0 = 1e2/1e3/1e4: 5.05e-3 / 5.12e-4 / 5.12e-5 (오차×K0/μ0 = 0.506, 0.512, 0.512) | 오차 < μ0/K0 (1차 추정 0.51 μ0/K0의 2배), 10년당 감소비 8–12 |
| Neo-Hookean 특수형 (C01 = 0) | 손으로 유도한 Abaqus NH 닫힌 식 σ = (2C10/J) dev(B̄) + 2(J−1)/D1 I, 측면 늘림도 이 식으로 이분법 | 단축 λ = 1.7: 3.4e-15. `neo_hookean` = `mooney_rivlin(C01=0)` R, K 비트 동일 | 1e-9 / 1e-12 |
| Newton 예측자 끔 = 켬 | 같은 이분법 기준 | 상태 차 ≤ 1.7e-13 | 1e-9 |
| Cook E04 (ν = 0.49999, p=2 C^1 16×16, Abaqus NH C10 = μ/2, D1 = 2/K) | Coreform 6.98946 | −0.066 % | 1 % (G5 게이트, S40과 동일) |

- 비압축 극한 1차 추정: 측면 σ = 0이므로 p = σ11/3 ≈ 0.49 μ0이고 J − 1 = p/K0다.
  d ln P / d ln J = −1 − (1/3) d ln σ/d ln λ + 1 ≈ −1.04(MR, λ = 1.5). 따라서 오차는 약 0.51 μ0/K0이고, 측정 0.506–0.512와 맞는다.
- 비압축 극한 반력 공칭응력과 QP Cauchy × 측면 면적비는 1e-8 이내에서 같다.

### 변이 시험 (스크래치 사본을 pytest 플러그인으로 주입, 저장소 무변경)
| 변이 | 실패한 시험 | 통과 |
|---|---|---|
| M1 ∂²Ī2 교차항 4/3 → 2/3 (접선만, C01) | FD MR 2건, 단축 λ=0.7 MR·예측자 끔(Newton 실패) 2건 | 나머지 20 |
| M2 Yeoh U″의 D2·D3 항 계수 (접선만) | FD Yeoh 2건 | 나머지 22 (응력 시험은 못 봄 → FD 시험이 필요한 이유) |
| M3 응력 출력의 등체적 1/J_q → 1/J̄ | "QP 응력 = 약형식 응력" 2건만 | 균일 변형 시험 전부(J̄ = J라 못 봄) |
| M4 dĪ1의 −2/3 I1 H → −1/3 | 21건 | 생성자·QP 출력 일관성 3건 |
| M5 평균 팽창 제거(U′, U″를 첫 QP의 J로) | 에너지 기울기 4, FD 4, QP 출력 2, 예측자 끔 1 | 나머지 13 (균일 변형·아핀 패치: J 균일이라 못 봄) |
| M6 W22 계수 j(j−1) → j² | FD 다항식 2건 | 나머지 28 (MR·Yeoh는 W22 = 0) |
| M7 W12 ×0.5 | FD 다항식 2건 | 나머지 28 |

M1–M5는 24건 시점(N=3 다항식·Cook 시험 추가 전), M6·M7은 30건 시점(Cook 시험 추가 전)에 실행했다.

## 측정: Cook E04, 두 Neo-Hookean 이산화 (p=2 C^1, `n_steps=10`, Pardiso, 다른 에이전트와 CPU 공유)
| 격자 | `fbar.py` NH (ln J 법칙) | 오차 | Abaqus NH (본 모듈) | 오차 | MR 0.8/0.2 (μ0 동일) | 오차 |
|---|---:|---:|---:|---:|---:|---:|
| 4×4 | 7.377563 | +5.553 % | 7.168385 | +2.560 % | 7.169395 | +2.574 % |
| 8×8 | 7.146492 | +2.247 % | 7.047734 | +0.834 % | 7.048273 | +0.841 % |
| 16×16 | 7.038057 | +0.695 % | 6.984879 | −0.066 % | 6.984804 | −0.067 % |

- fbar 16×16 값은 S40의 7.0380574와 같다.
- 시간(16×16): fbar 33.3 s, Abaqus NH 41.9 s, MR 37.1 s. 11–12 증분, 64–86 반복.
- ν = 0.49999에서 두 NH 법칙의 연속체 차이는 O(μ/K) = 3e-5뿐이다. 따라서 위 차이는 이산화에서 온다.
  fbar는 μ/2(I1 − 3)를 점별로 두고 −μ ln J̄ + λ/2 ln² J̄만 셀 평균한다. 그래서 점별 부분이 등체적이 아니다.
  본 모듈은 체적 응답 전체를 P0로 본다. 8×8 이하에서 본 모듈이 기준에 더 가깝다.
  16×16에서는 기준 아래로 내려갔다(−0.066 %). 어느 극한으로 수렴하는지는 32×32를 측정하지 않아 모른다.
- Newton 예측자 효과(2×2×1 p=2, 균일 단축 MR):
  - λ = 0.7: 15 → 4 증분(88 → 12 반복)
  - K0/μ0 = 1e4: 80 → 8 증분(623 → 24 반복). 1e5: 215 → 8 증분
  - 수렴 상태 차는 ≤ 1.7e-13이다.

## 남은 것
1. **Ogden 미구현.** 주늘림 고유분해와 중근(λ_a = λ_b) 근방 접선의 극한 공식이 필요하다.
   같은 기준(FD 1e-6, 에너지 complex-step 기울기)으로 검증하려면 중근을 포함한 상태를 따로 시험해야 해서 이번 범위에서 뺐다.
2. **StudySpec / `methods` / User API 연결 안 됨.** `methods_solve.py` 등은 다른 에이전트가 수정 중이라 건드리지 않았다.
   재료 스펙(`hyperelastic: {model, C_ij, D_k}`)과 `ReferenceNumbaBackend` 경로 추가가 다음 단계다.
3. **`runner.post` 탐침 응력 미연결.** 현재는 QP 응력(`cauchy_stress`)과 셀 J̄(`cell_jbar`)만 있다. `fbar_point_stress_strain`에 대응하는 함수가 필요하다.
4. **Kratos 결합(접촉·동역학·열) 없음.** S40과 같다. 필요해지면 C++ UEL 승격을 결정한다.
5. **절단(immersed) 셀 미지원.** `fbar_trimmed`의 압력 그룹 방식을 이 커널에 붙여야 한다.
6. **완전 비압축(D1 = 0, hybrid) 미지원.** Drucker 안정성 검사도 없고 μ0 > 0만 확인한다.
7. **성능.** QP마다 4계 A(81성분)를 전부 만들고 (I, J) 블록 대칭을 쓰지 않는다. Cook 16×16에서 fbar보다 1.1–1.3배 느리다(측정 잡음 범위 큼).
8. **fbar 반영 결정 필요.** Dirichlet 예측자는 `fbar.solve_static`에 없다. `fbar.py`는 이번 작업에서 수정 금지라 반영할지 사용자가 결정해야 한다.

## 재현
- `D:\mamba\envs\whtiga\python.exe -m pytest -q -p no:cacheprovider tests/unit/test_hyperelastic_fbar.py tests/test_isolation.py` (cwd `WHT_IGA`)
  → 35 passed, 80.2 s(Cook 시험 36 s, Numba 컴파일 포함).
- 스크래치(저장소 밖) `C-hyperelastic/`:
  - 측정: `probe1.py`, `measure.py`, `measure2.py`, `measure_poly.py`
  - Newton 예측자: `predictor.py`
  - Cook: `cook_nh.py`
  - 변이: `make_mutants.py`, `make_mutants2.py`, `mutplug.py`, `mutants/results.txt`
