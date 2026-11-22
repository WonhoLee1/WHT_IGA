# §10.4 Kratos 교차검증 — Kratos 경로 대 참조 조립기 (V 케이스 축소판)

- 계획 기준 (master-plan v2.1 §10.4, §5.5): "모든 V 케이스의 축소판에서 Kratos 경로와 참조 조립기 결과 차이 ≤ 1e-10".
- 날짜 / 환경: 2026-11-20 (같은 날 검증자 지적을 반영해 개정), WHTIGA 브랜치, `D:\mamba\envs\whtiga\python.exe` (Python 3.12), Kratos 10.4.3 (`E:\kratos-build\install`).
- 시험: `tests/kratos_spikes/test_p6_cross_validation_10_4.py` (marker `kratos`, 38개 시험)
- 재현 (cwd `WHT_IGA`):
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests/kratos_spikes/test_p6_cross_validation_10_4.py tests/test_isolation.py
  ```
  결과: **`42 passed in 31.99s`** (이 파일 38개 + `test_isolation.py` 4개). `-s`로 실행하면 사례마다 `XVAL ...` 한 줄에 측정값이 찍힌다. 아래 수치는 모두 이 실행의 출력이다. 실행마다 바뀌는 값은 Kratos 쪽 Pardiso 스레딩과 V09 ARPACK 시작 벡터에서 오는 반올림 수준 변동뿐이다(예: V09 7.98e-13–8.98e-13).

## 1. 결론

| 구분 | 내용 |
|---|---|
| 1e-10 충족: 행렬·하중 벡터·해(+응력·반력)를 모두 비교 | V01(BC-B), V03(강체 고정, 3하중, p 1–3), V02 절단판(drop p 1–3, **FCM ghost α = 1e-8 p 1–2**), V07(ν = 0.3). 해 최대 1.3e-12, 행렬 ≤ 1.3e-15, 하중 벡터 ≤ 7.6e-15 |
| 1e-10 충족: 해만 비교 (하중 벡터는 참조에서 넣음) | V04(2부품 tie) 7.7e-12, V05(7 슬랩) 7.9e-12. 두 시험은 참조 면 견인력 벡터를 Kratos에 `PointLoadCondition3D1N`으로 넣는다. 그래서 강성·tie·셀 영역 경로만 비교하고 Kratos 표면 하중 조립은 비교하지 않는다 |
| 1e-10 충족: 기타 | V09(모드), V10(Bossak), V11(영역 K), V12(전도), V-EXP(explicit) |
| 1e-10 충족, 단 **같은 표면 적분 규칙 기준** | V01(BC-A), V03(BC-A 고정), 대칭면 BC-B. 정확 적분 기준과는 행렬 3.5e-3–6.9e-3, 해 8.8e-10–2.3e-8 차이 → **미해결 이슈 1** |
| 전방 1e-10 미충족 (기록, 허용치는 낮추지 않음) | V07 ν = 0.49999 (κ = 2.39e9): 행렬 5.3e-16, 하중 벡터 1.1e-15, 역방향 오차 6.0e-18로 조립된 이산계는 같다. 그러나 Kratos U는 정밀 참조해와 **3.5e-9** 다르다. Kratos 자신의 배정밀도 선형 해법 오차다(이슈 3) |
| 비교 불가 (참조 경로 없음) | V08, V-EXP-C (접촉), V03 follower·모멘트(RBE3), V06 내부 구면 압력, V07 F-bar, V11 Nitsche·고유변형 하중, V12 필름·복사·용량 (모두 Python 단일 구현이거나 참조 조립기 미지원) |

§10.4 항목은 **조건부로 닫는다**. 참조 경로가 있는 모든 V 케이스 축소판에서 Kratos가 조립한 행렬과 하중 벡터가 참조와 반올림 수준(행렬 ≤ 1.3e-15, 하중 ≤ 7.6e-15)으로 같고, V07 ν = 0.49999 하나를 빼면 해도 1e-10 안에서 같다. 남는 것은 두 가지다. (a) 벌칙(penalty) 경계면 적분 규칙이 정확 적분과 다르다(이슈 1). (b) V07 ν = 0.49999에서 Kratos 선형 해법의 전방 오차가 3.5e-9다(이슈 3).

**이전 판 정정.** 이전 판은 V02 ghost α = 1e-8 (p = 1, 2)과 V02 drop p = 3을 "배정밀도로 1e-10 도달 불가(조건수)"로 분류하고 전방 단언을 껐다. 이 판정은 틀렸다. 차이(6.4e-8, 7.8e-8, 1.09e-10)는 **참조 쪽** 해법(스케일링 없는 scipy `spsolve`/SuperLU)의 오차였다. 정밀 참조해(§2)와 비교하면 Kratos U는 6.5e-15, 2.8e-14, 6.3e-13만큼 다르다. 이 세 사례는 이제 1e-10 전방 단언을 한다. p = 3의 1.09e-10은 허용치 바로 위의 측정값이었고, 단언을 끈 것은 측정값에 맞춘 게이트였다.

**"같은 이산계"의 근거.** 이 판은 그 근거를 행렬 비교(≤ 1e-12)와 하중 벡터 비교(≤ 1e-10)에 둔다. 역방향 오차만으로는 같은 이산계라고 주장하지 않는다. 역방향 오차는 하중 오차를 거의 보지 못한다. Kratos 하중을 1 + 1e-9배 하는 변이(§6 N2a)에서 역방향 오차는 1.2e-16에 그쳤다. 하중 벡터 비교를 끄면 이 변이는 살아남는다(N2b).

## 2. 참조 측 구성

Kratos를 import하지 않는다.

- `wht_iga.reference.assemble_elastic` (`assemble_stiffness`, `assemble_face_traction`, `isotropic_D`), `backends.reference_numba.body_force_vector`, `analyses.thermal.assemble_conduction`, `kratos_backend.mass.assemble_box_mass`(Python 일관 질량, Kratos 미사용), scipy `eigsh`/Bossak (`wht_iga.dynamics`).
- 참조 조립기(`wht_iga.reference`)에는 **절단 셀(QP 표)·벌칙 경계·셀별 재료 경로가 없다**. 그래서 시험 파일 안에 numpy oracle 3개를 두었다.
  - `_qp_stiffness`: 임의의 `[x, y, z, w]` 점에서 Σ w·s·BᵀDB를 계산한다. 전체 Gauss 규칙에서 `assemble_stiffness`와 **8.8e-16 / 1.0e-15 / 1.1e-15** (p = 1/2/3)로 일치한다(`test_oracle_matches_reference_assembler`).
  - `_penalty_exact`: β∫N Nᵀ dA를 면 셀마다 텐서 Gauss p+1점으로 적분한다(정확 적분).
  - `_penalty_same_rule`: 같은 피적분 함수를 파이프라인과 같은 면 삼각형(`box_face_triangles`)과 Dunavant 규칙(`rule_for_degree`, `triangle_quadrature`, 순수 Python 표)으로 적분한다. 벌칙 합(분할 단위성) 검사는 이제 1e-13으로 단언한다. 측정값은 exact 3.3e-16–4.4e-16, same 3.3e-16–1.8e-15다. 1e-13은 15자리 표와 10³개 이하 항의 반올림에서 정한 상한이다.
- **정밀 참조 해법 `_ref_solve`** (이번 개정의 핵심):
  1. 대칭 Jacobi 스케일링 D K D (D = diag(K)^-1/2) 후 dense Cholesky로 푼다.
  2. 반복 개선(iterative refinement)을 한다. 잔차 f − K u는 `_residual`로 **항마다 올바르게 반올림해** 계산한다. 곱은 Dekker TwoProduct로 오차 없이 두 double로 나누고, 행 합은 `math.fsum`(정확 반올림)으로 한다.
  3. 보정량이 더 이상 절반으로 줄지 않거나 eps 이하가 되면 멈춘다. 마지막 보정량의 상대 크기 `ref_dx`가 REF_ACC = 1e-13(기준의 1/1000)을 넘으면 시험이 실패한다. 그래서 참조가 수렴하지 않은 채 통과할 수 없다. 모든 정적 사례에서 개선 단계 `ref_ir` = 2, `ref_dx` = 5.5e-17–1.1e-16이었다.
  - 정지 조건 메모: 작업 지시는 "잔차가 더 이상 줄지 않을 때까지"였다. 처음 구현은 그대로 했지만, 조건수가 큰 계에서는 첫 Cholesky 해의 잔차가 이미 반올림 바닥(~eps·|K||u|)에 있다. 그래서 해가 크게 좋아져도 잔차는 줄지 않는다. pascal(12)에서 첫 보정량이 4.7e-7인데 잔차는 그대로였고, 잔차 기준 정지는 그 자리에서 멈췄다(REF_ACC 단언이 잡았다). 그래서 보정량 기준(Higham, *Accuracy and Stability of Numerical Algorithms* §12.1의 정지 판정)으로 바꾸고 잔차 이력은 함께 기록한다.
  - 독립 검증 (`test_exact_residual_and_refined_solve`): (a) 상쇄가 심한 40×40 계에서 `_residual`이 `fractions.Fraction` 정확 유리수 잔차를 한 번 반올림한 값과 **비트 단위로 같다**. (b) pascal(12)(SPD, 정수 성분, κ = 8.8e11, f = P u_true가 binary64로 정확)에서 `_ref_solve`는 u_true를 **오차 0**으로 복원한다. 스케일링 없는 `spsolve` 오차는 1.9e-5다.
  - 추가 독립 검증 (저장소 밖, 1회): V07 ν = 0.49999 축약계(n = 270)를 50자리 `decimal` 가우스 소거(부분 피벗)로 풀었다. IR 참조와 **6.2e-25**, Kratos U와 **3.51e-9** 차이다(`C:\Users\GOODMAN\AppData\Local\Temp\claude\694f846d5550b343\57a213e2-2165-46f7-af77-ce5cc827bb7c\scratchpad\B-xval-fix\mp_check.py`, 2.1 s).
- Gauss 점은 시험 안에서 직접 만든다(`kratos_backend.geometry.cell_gauss_points`는 쓰지 않는다. 변이 M2 참조).
- 확장(extension) 연산자 T는 파이프라인과 같은 `stabilize.support/extension`으로 만든다. T 자체는 비교 대상이 아니라 같은 입력이다. Kratos는 T를 MPC로 처리하고(block builder 소거) 참조는 PᵀKP 축약계를 `_ref_solve`로 푼다.

## 3. Kratos 측 추출

- 행렬 `_kratos_matrix`: Kratos가 실제로 인수분해한 strategy 시스템 행렬을 참조 DOF 순서로 바꾼 것이다(MPC 적용 후 TᵀKT, Dirichlet 행·열 0). 비교는 풀린 계의 DOF에서만 한다. 강체 고정 사례는 자유 DOF, V02는 slave도 고정도 아닌 DOF(`keep`), BC-B는 spline DOF 전체(가상 노드는 slave라 소거됨)다.
- 하중 벡터 `_kratos_rhs`: 풀이가 끝난 model part를 새 `ResidualBasedBlockBuilderAndSolver`로 u = 0에서 다시 조립한다(`Build`). 선형 요소의 RHS = 외력이 된다(−K u 항이 0). 이어서 model part의 MPC를 적용한다(`ApplyConstraints`, Tᵀb). 그러면 BC-B 가상 노드 하중 k·g와 확장 slave 하중이 master로 옮겨진다. Dirichlet 행은 0으로 만들지 않으므로 고정면 하중도 비교된다. 변위는 끝난 뒤 원래 값으로 되돌린다. 방정식 번호가 원래 풀이와 같은지 단언한다. 고정 DOF의 처방값이 0인지도 단언한다(아니면 b(0)이 외력이 아니다). 대상은 `SURFACE_LOAD`·압력 조건, 체적력, `SupportPenaltyCondition`의 g 항, BC-B 가상 노드 점하중이다. 하중 벡터는 Kratos 자신의 조립 결과이고 참조 벡터가 아니다.
- V04/V05는 `run_linear`을 거치지 않는다. 참조 면 견인력 벡터를 그대로 점하중으로 넣으므로 하중 비교는 의미가 없다(같은 벡터). 그래서 해만 비교한다.

판정량: `max|a − ref| / max|ref|` (참조 장의 최댓값 기준 상대). 행렬은 `max|ΔK| / max|K|`다. 역방향 오차는 `‖K u_k − f‖∞ / (‖K‖∞‖u_k‖∞ + ‖f‖∞)` (u_k = Kratos 해, K·f = 참조계, 잔차는 올바르게 반올림)다.

## 4. V 케이스별 표

t_K = 해당 사례의 Kratos 경로 실행 시간(초, 위 실행). 전체 파일 실행 시간은 약 32 초다.

| V | 물리 | 참조 경로 | 기존 비교 시험 (파일::시험, 허용치) | 이 파일의 시험 | 측정 최대 차이 | 판정 | t_K [s] |
|---|---|---|---|---|---|---|---|
| V01 | 선형 정적, 강제 변위 정육면체 (E = 1, ν = 0) | 있음 (box + 벌칙 oracle) | `test_s04_sm_elements.py::test_global_stiffness_matches_reference` (K 1e-12, 파이프라인 아님), `test_p2stab_compose_kratos.py::test_kratos_tkt_matches_reference_1e10` (TᵀKT 해 1e-10, 파이프라인 아님), `test_p1_pipeline_cube.py` (해석해 1e-5) | `test_v01_bca_affine_cube[p1-3]` (BC-A 6면 affine), `test_v01_bcb_normal_components[p1-3]` (면마다 법선 성분 BC-B) | BC-A 같은 규칙: U 1.0e-13 / K 8.5e-16 / f 6.9e-16; **정확 적분: U 2.3e-8 / K 6.9e-3** (f 3.8e-15). BC-B: U 8.0e-14 (정확 규칙 3.5e-13) / K 8.8e-16 / f 0 (비트 단위 일치, 변이 N3이 1e-9 검출) | BC-A: 같은 규칙 기준 통과, 정확 적분 기준 **미해결 이슈 1**; BC-B 통과 | 0.07–0.83 |
| V02 | 선형 정적, 원공 판 (절단 IGA) | QP oracle (참조 조립기에는 절단 경로 없음) | 없음 | `test_v02_trimmed_notched_plate[p, mode]`: 모서리 1/4 원 노치 판 (r = 1.09, 4×4×1, η = 0.004 sliver, OUT 1셀, 절단 15셀, 확장 θ = 0.2, slave 4/6/8 CP), drop p 1–3 + ghost α = 1e-8 p 1–2 | U: drop 5.8e-15 / 3.0e-14 / 6.3e-13, ghost 6.5e-15 / 2.8e-14. σ ≤ 4.7e-14 (ghost 셀은 α·E 응력으로 비교). K(keep) ≤ 1.3e-15, f(keep) ≤ 3.6e-15, 역방향 ≤ 8.1e-17, slave 복원 0 | 5개 모두 1e-10 통과 | 0.04–0.69 |
| V03 | 선형 정적 외팔보 (E = 200 GPa, ν = 0) | 있음 | `test_user_api_kratos.py::test_kratos_bar_matches_closed_form_reference_and_exposes_native` (벌칙 대 강체 고정, 다른 BC라 atol 1e-6·scale) | `test_v03_cantilever_strong_clamp[p1-3 × traction/pressure/gravity]` (하중 벡터, U, 모든 QP의 Cauchy 응력, 고정면 반력), `test_v03_cantilever_penalty_clamp` | 강체 고정: U 1.3e-12, σ 1.3e-12, 반력 1.4e-12, f 7.6e-15. BC-A 고정 p = 2: U 1.6e-12 (같은 규칙·정확 모두), K 5.2e-16, f 1.4e-15 | 통과 | 0.01–0.16 |
| V03 follower / 모멘트 | follower 압력, RBE3 모멘트 | **없음** (`ReferenceNumbaBackend`는 follower에서 `BackendUnsupported`, RBE3 참조 없음) | `test_p3_gate_g3_benchmarks.py::test_g3_01_...` (해석해) | — | — | 비교 불가 | — |
| V04 | 이종 재료 2부품 tie | 있음 (참조 box 2개를 계면 CP 공유로 접합 = C⁰ merge와 수학적으로 같음) | `test_p3_gate_g3_benchmarks.py::test_g3_02_...` (merge 대 tie 1e-6, 해석해 1%) | `test_v04_two_part_tie` (`AssemblyBuilder` + `mpc_matching`, 2×2×8, p = 2) | U 7.7e-12 | 통과. **참조 면 견인력 벡터를 Kratos에 점하중으로 넣는다 → 강성·tie 경로만 비교**. merge(이중 knot) Kratos 모델은 기존 시험 안에서 손으로 만든 것이라 라이브러리 진입점이 없어 따로 비교하지 않음 | 0.21 |
| V05 | 7 단면 재료 경계, 공유 spline | 셀별 E oracle | `test_g3_03_...` (모든 영역 같은 E로 monolithic 대비 3.3e-6) | `test_v05_cell_domains_shared_spline` (`CellDomainSpec` 7개, E_k = 1e11·(1 + 0.5k)) | U 7.9e-12 | 통과. **참조 면 견인력 벡터를 점하중으로 넣는다 → 강성·셀 영역 경로만 비교**. 기존 시험은 영역 E가 모두 같아 셀별 재료 경로를 검증하지 못했다 | 0.10 |
| V06 | 내압 구 octant (절단, BC-B 대칭) | 체적·대칭면: 있음 (V02 절단판, 롤러 시험이 같은 코드 경로) / 내부 구면 압력: **없음** | `test_p2mesh_v06.py` (Lamé·CSV) | 대리: `test_v02_trimmed_notched_plate`, `test_symmetry_rollers_bcb_bending` | 롤러 BC-B p = 2: U 4.4e-15, K 5.2e-16, f 1.7e-15 (같은 규칙), **U 8.8e-10 (정확)** | 체적·BC-B는 위 판정과 같음. 구면 압력은 `test_p2mesh_v06.py` 안의 시험 국소 workaround(`make_surface_quadrature` + `LoadCondition`)라 진입점·참조 없음 | 0.20 |
| V07 | ν = 0.49999 Kirsch (주 경로 F-bar Python) | 소변형: 있음 / F-bar: Python만 (Kratos F-bar QP 요소 없음) | `test_p5_v07_kirsch_precursor.py` (해석해) | `test_v07_near_incompressible_block[0.3, 0.49999]` (4×4×0.25 블록, p = 2, 비절단) | ν = 0.3: U 3.2e-13, K 8.0e-16, f 1.1e-15. ν = 0.49999: **U 3.5e-9**, K 5.3e-16, f 1.1e-15, 역방향 6.0e-18, ‖Δu‖₂/‖u‖₂ 2.5e-9 ≤ 잔차 상한 1.6e-8 | ν = 0.3 통과. 0.49999는 행렬·하중·역방향·잔차 상한 통과, 전방은 기록(이슈 3) | 0.09 |
| V08 | Hertz 접촉 + 마찰 | **없음** (접촉 SSN/mortar의 참조 조립기 없음) | `test_p5_v08_*` (Hertz 해석해) | — | — | 비교 불가 | — |
| V09 | 모드 (K Kratos + M Python) | 있음 | `test_p4_gate_g4_benchmarks.py::test_v09_...` (해석해), `test_p4_eigen_and_bossak.py` (양쪽 모두 Kratos K) | `test_v09_cantilever_modes` (0.1×0.1×1 보, 2×2×10, p = 2, 10모드) | f 8.2e-13 (f1 = 80.27 Hz) | 통과 | 0.18 |
| V10 | Bossak 과도 | 있음 | `test_p4_eigen_and_bossak.py::test_kratos_bossak_matches_python_bossak_with_lumped_mass` (u, v 1e-10, a 1e-9 — **양쪽 모두 Kratos K**, 스킴 비교만) | `test_v10_bossak[-0.05, -0.30]` (Kratos K·질량 대 참조 K + Python 행합 질량) | u 5.0e-15, v 1.5e-14, a 5.7e-14 | 통과 | 1.1 |
| V11 | 열변형 응력 (이중 재료 Nitsche) | 영역 K: QP oracle / Nitsche·고유변형: Python만 | `test_p4_bimaterial.py` (Kratos 없음), `test_g4::test_v11_...` (해석해) | `test_v11_region_stiffness` (`kratos_region_stiffness`, 11×11×1, p = 1, 실제 1/4 원 분할) | K 7.5e-15 | Kratos 부분(영역 K) 통과. 나머지는 Python 단일 구현 | 0.24 |
| V12 | 열전도 (K_th Kratos, 나머지 Python) | 있음 (`assemble_conduction`) | `test_p4_thermal_kratos.py::test_kratos_conductivity_equals_python_assembly` (K 1e-12) | `test_v12_conduction[p1-2]` (행렬 + 2D 정상 온도장) | K 5.2e-16, T 1.4e-15 | 통과. 필름·복사·용량은 Python만 | 0.02 |
| V-COND | 조건수 | 있음 | `test_p2stab_kappa_kratos.py::test_kratos_kappa_matches_reference` (K atol 1e-10, κ 5%) | — (기존으로 충분) | — | 기존 시험으로 충족 | — |
| V-EXP | explicit 1D 파동 | 있음 | `test_vexp_explicit.py::test_vexp_wavefront_arrival_gap_pin` (1D SciPy 참조, 도달 시각 ≤ 1e-9 s; 3D 참조 조립기 아님) | `test_vexp_explicit_bar[p1-2]` (S24b 파동 모델, 10셀, CFL 0.2, 100 스텝, 모든 CP 이력) | 이력 1.4e-14 | 통과 | 0.04–0.10 |
| V-EXP-C | explicit 벌칙 접촉 | **없음** | `test_vexp_explicit.py::test_vexp_c_dt_sweep` (운동량·에너지) | — | — | 비교 불가 | — |

### 4.1 수정 전후 최대 차이 (해 U, 사례별)

- (A) 초판: Dunavant 오타가 있었고, 참조는 `spsolve`였다(초판 문서 값).
- (B) `surface.py` 표 수정 뒤, 수정 전 시험 파일(`spsolve` 참조)로 이번에 다시 잰 값(`37 passed in 29.06s`).
- (C) 이번 개정(정밀 참조 + 행렬·하중 비교).

| 사례 | (A) | (B) | (C) | (C) 단언 |
|---|---|---|---|---|
| V01 BC-A (같은 규칙, p 1–3 최대) | 1.0e-13 | 1.1e-13 | 1.0e-13 | U·K·f |
| V01 BC-B (p 1–3 최대) | 9.3e-14 | 7.5e-14 | 8.0e-14 | U·K·f |
| V03 강체 고정 (9 사례 최대 U / σ / 반력) | 3.0e-11 / 3.1e-11 / 3.1e-11 | 2.3e-12 / 2.8e-12 / 2.8e-12 | 1.3e-12 / 1.3e-12 / 1.4e-12 | U·σ·반력·f |
| V03 BC-A 고정 p = 2 | 3.2e-11 | 1.0e-12 | 1.6e-12 | U·K·f |
| 롤러 BC-B p = 2 | 3.2e-11 | 4.9e-15 | 4.4e-15 | U·K·f |
| V07 ν = 0.3 | 3.5e-11 | 4.6e-13 | 3.2e-13 | U·K·f·역방향 |
| V07 ν = 0.49999 | 6.2e-9 | 6.1e-9 | **3.5e-9** | K·f·역방향·잔차 상한 (U 기록) |
| V02 drop p = 1 | 6.9e-15 | 6.9e-15 | 5.8e-15 | U·σ·K·f·역방향 |
| V02 drop p = 2 | 3.2e-11 | 8.4e-13 | 3.0e-14 | U·σ·K·f·역방향 |
| V02 drop p = 3 | 1.09e-10 (미단언) | 1.09e-10 (미단언) | **6.3e-13** | U·σ·K·f·역방향 |
| V02 ghost p = 1 | 6.4e-8 (미단언) | 6.4e-8 (미단언) | **6.5e-15** | U·σ·K·f·역방향 |
| V02 ghost p = 2 | 7.8e-8 (미단언) | 7.8e-8 (미단언) | **2.8e-14** | U·σ·K·f·역방향 |
| V04 tie | 1.5e-12 | 1.5e-12 | 7.7e-12 | U |
| V05 7 슬랩 | 5.3e-12 | 5.3e-12 | 7.9e-12 | U |

(A)→(B)의 감소는 Dunavant 표 수정 때문이다(이슈 2). (B)→(C)의 V02·V07 감소는 참조 해법 교체 때문이다. 시험 출력의 `u_plain_spsolve_ref`가 같은 Kratos U를 `spsolve` 참조와 비교한 값이며 (B)와 같다. V04/V05는 (C)에서 조금 커졌다. `spsolve`가 우연히 Kratos와 같은 쪽으로 반올림했던 것이고, 정밀 참조 기준 7.7e-12 / 7.9e-12가 Kratos의 실제 차이다.

## 5. 미해결·해결 이슈 (근거 포함)

### 이슈 1 — 벌칙 경계면 행렬이 정확 적분이 아니다 (BC-A, BC-B) — 미해결

- 위치: `wht_iga/solver/penalty.py::rule_for_degree` (p = 1/2/3 → Dunavant 차수 2/4/6), 면 셀마다 삼각형 2개.
- 피적분 함수 β N_i N_j는 삼각형 위에서 총차수 4p다. 이 규칙은 총차수 2p(= ∫N_i dA)까지만 정확하다. Kratos `SupportPenaltyCondition`의 식은 β·w·detJ·HᵀH(`support_penalty_condition.cpp` l.53–90)로 맞고, 차이는 적분 규칙에서만 생긴다.
- 측정 (`test_v01_bca_affine_cube`): Kratos 조립 행렬은 같은 규칙 참조와 **4.4e-16 – 8.5e-16**으로 같다. 정확 적분 참조와는 **6.9e-3 / 3.5e-3 / 4.8e-3** (p = 1/2/3, 최대 항 기준) 다르다. 하중 벡터(β∫g N dA)는 정확 적분과도 ≤ 3.8e-15로 맞는다. 변위 차이는 O(1/β)로 줄어든 **2.1e-8 / 4.9e-9 / 2.3e-8**이다(V01 BC-A, 기본 β = 3e5·E·p²/h). BC-B 롤러 굽힘에서는 **8.8e-10**이다. 경계 편차가 균일하면 ∫N_i dA만 들어가므로 차이가 사라진다(V01 BC-B ≤ 3.5e-13, V03 BC-A g = 0: 1.6e-12).
- 영향: V01 1e-5 게이트 등 정확도 게이트에는 영향이 없다(O(1/β) 항). 그러나 정확 적분 참조 조립기 기준 §10.4 1e-10은 벌칙 경계가 있는 케이스에서 성립하지 않는다.
- 제안: 면 셀을 사각형 텐서 Gauss (p+1)²점으로 적분하거나, 총차수 4p 이상의 삼각형 규칙을 쓴다(`triangle_quadrature`에는 차수 6까지만 있다).

### 이슈 2 — Dunavant 차수 4 표 가중치 전사 오타 — **해결 (2026-11-20)**

- 원인: `kratos_backend/surface.py::triangle_quadrature(4)`의 b 가중치가 `0.0549758718227661`로 적혀 있었다. 바른 값은 `0.054975871827661`(Dunavant 0.109951743655322의 절반)이다. 숫자 중간에 `8`이 하나 더 끼어든 전사 오타다. 초판은 이를 "14자리 표"라고 썼는데 틀린 설명이었다. 자릿수가 아니라 오타가 원인이다.
- 수정 전 가중치 합 3(wa + wb) = 0.49999999998544836 (0.5 대비 −1.46e-11, 상대 −2.9e-11). 수정 후 (`wa = 0.1116907948390055`, `wb = 0.054975871827661`) 합 = 0.4999999999999995 (상대 −1.0e-15).
- 재측정 (§4.1 (A)→(B)): V03 p = 2 표면 하중 3.0e-11 → 4.6e-13, V07 ν = 0.3 3.5e-11 → 4.6e-13, V02 drop p = 2 3.2e-11 → 8.4e-13, BC-A 고정·롤러 3.2e-11 → 1.0e-12 / 4.9e-15. 벌칙 합 oracle 검사는 1e-10에서 1e-13으로 좁혔다(측정 ≤ 1.8e-15).
- `test_user_api_kratos.py`의 "~3e-11 relative (measured)" 주석도 같은 원인이었다. 이 파일은 이 작업의 수정 범위 밖이라 손대지 않았다.

### 이슈 3 — V07 ν = 0.49999: Kratos 선형 해법의 전방 오차 3.5e-9 — 기록

| 양 | 값 |
|---|---|
| κ(축약 K), eps·κ | 2.39e9, 5.3e-7 |
| κ(D K D), Jacobi 스케일링 후 | 2.51e9 (스케일링으로 줄지 않음) |
| 행렬 차이 / 하중 벡터 차이 | 5.3e-16 / 1.1e-15 |
| 역방향 오차 (올바르게 반올림한 잔차) | 6.0e-18 |
| Kratos U − 정밀 참조 (max, 상대) | **3.51e-9** |
| 정밀 참조 − 50자리 decimal 해 | 6.2e-25 |
| `spsolve` 해 − 정밀 참조 | 2.7e-9 (`spsolve` 대비 Kratos 6.1e-9) |
| ‖u_k − u‖₂ / ‖u_k‖₂ 대 상한 ‖K u_k − f‖₂ / (λ_min ‖u_k‖₂) | 2.5e-9 ≤ 1.6e-8 (단언) |

- 해석: Kratos가 조립한 행렬과 하중 벡터는 참조와 반올림 수준으로 같다. 참조 해는 50자리 해와 6.2e-25까지 맞는다. 따라서 3.5e-9는 Kratos 자신의 배정밀도 선형 해법(`pardiso_ldlt`) 오차다. 이 오차는 Kratos 해의 잔차가 허용하는 범위(엄밀한 상한 1.6e-8) 안에 있다. ν = 0.49999의 나쁜 조건은 체적/전단 강성 비(≈ 1/(1 − 2ν) = 5e4)에서 오는 물리적 조건수라서 대각 스케일링으로 없어지지 않는다(κ(DKD) = 2.51e9). 같은 배정밀도 직접 해법인 `spsolve`도 2.7e-9 틀린다. 이 사례에서 1e-10 전방 일치를 얻으려면 Kratos 쪽에 확장 정밀도 반복 개선 같은 해법 변경이 필요하다. 라이브러리 해법 설정은 이 작업 범위가 아니어서 바꾸지 않았다.
- 시험: 행렬(1e-12)·하중(1e-10)·역방향(1e-10)과 잔차 상한(허용치 없는 부등식)을 단언하고, 전방 3.5e-9는 출력·기록만 한다. ν = 0.3 쌍은 전방 1e-10을 단언한다(3.2e-13). 허용치는 낮추지 않았다.
- 정정: 초판은 V02 ghost α = 1e-8 (κ 3.2e9, 5.1e9)과 V02 drop p = 3 (κ 5.5e8)도 이 이슈로 분류했다. 실제로는 참조 측 `spsolve` 오차였고, Kratos 해는 정밀 참조와 ≤ 6.3e-13으로 맞는다(§1). 그 사례들의 나쁜 조건은 축약계 대각의 크기 차이(ghost α·E, sliver)에서 왔다. Jacobi 스케일링 뒤 κ(DKD)는 ghost p = 1 8.8e2, ghost p = 2 2.0e3, drop p = 3 3.8e3이다(스케일링 전 3.2e9, 5.1e9, 5.5e8). 그래서 스케일링 + 반복 개선으로 참조 오차가 사라졌다. 측정 스크립트: `C:\Users\GOODMAN\AppData\Local\Temp\claude\694f846d5550b343\57a213e2-2165-46f7-af77-ce5cc827bb7c\scratchpad\B-xval-fix\kappa_scaled.py`.

### 참고 — ghost + 확장에서 master도 slave도 아닌 CP

OUT 셀에만 지지되는 CP(예: 노치 모서리 CP 0, 25)는 확장 연산자에서 master도 slave도 아니다. Kratos 경로에서는 MPC 없이 α·E ghost 강성만 가진 **자유 DOF**로 남는다. 참조 측이 확장 행렬 T(이 CP의 행이 0)를 그대로 쓰면 그 CP에서 상대 차이가 0.8까지 난다(초기 실험에서 확인). 이 시험의 참조는 Kratos MPC 소거와 같게 "slave가 아닌 모든 DOF = 항등"으로 만든다. 이번 개정에서 이 DOF들을 포함한 `keep` 블록의 행렬이 Kratos와 1.3e-15 이하로 같음을 단언한다. 결함은 아니지만 다른 참조 구현이 같은 함정에 빠지지 않도록 기록한다.

## 6. 시험 검출력 (변이, 저장소 파일은 수정하지 않음)

`monkeypatch`만 쓰는 임시 스크립트로 실행했다(`C:\Users\GOODMAN\AppData\Local\Temp\claude\694f846d5550b343\57a213e2-2165-46f7-af77-ce5cc827bb7c\scratchpad\B-xval-fix\mutations.py`, 저장소 밖, 로그 `C:\Users\GOODMAN\AppData\Local\Temp\claude\694f846d5550b343\57a213e2-2165-46f7-af77-ce5cc827bb7c\scratchpad\B-xval-fix\mutations.log`). N은 이번 개정용 변이다. M은 초판 변이를 새 파일에 다시 돌린 것이다.

| 변이 | 대상 시험 | 결과 |
|---|---|---|
| N1a 참조 해법을 예전 `spsolve`로 되돌림 | `test_v02_...[1-ghost]` | 검출 (U 6.4e-8) |
| N1b 같음 | `test_v02_...[3-drop]` | 검출 (U 1.09e-10) |
| N2a Kratos 견인력 ×(1 + 1e-9) (전방 미단언 사례) | `test_v07_...[0.49999]` | 검출, **하중 벡터 비교만** (f 1.0e-9). 역방향 오차는 1.2e-16으로 통과 수준 |
| N2b N2a + 하중 벡터 비교 끔 | 같음 | **생존** → 이 사례에서 하중 오류를 잡는 검사는 하중 벡터 비교뿐이다 |
| N3 Kratos BC-B 처방값 ×(1 + 1e-9) | `test_v01_bcb_normal_components[2]` | 검출 (f 1.0e-9; f 비교 0.0이 공허한 일치가 아님을 확인) |
| N4 ghost α 1e-8 → 1.1e-8 (Kratos 측) | `test_v02_...[2-ghost]` | 검출 (K 2.6e-9, U 1.1e-5) |
| N5 절단 셀 QP 하나 제거 (Kratos 측) | `test_v02_...[3-drop]` (초판 미단언 사례) | 검출 (QP 수 1794 ≠ 1795) |
| N6 `_residual`을 일반 double f − K u로 | `test_exact_residual_and_refined_solve` | 검출 (비트 일치 실패, 유리수 정확값과 최대 절대 차이 9.2e-10) |
| M1 벌칙 규칙 차수를 모든 p에서 6으로 (Kratos 측) | `test_v01_bca_affine_cube[2]` | 검출 (K 3.4e-3) |
| M2 Kratos 체적 Gauss 가중치 ×(1 + 1e-9) | `test_v03_cantilever_strong_clamp[2-gravity]` | 검출 (f 1.0e-9, 반력 1.0e-9; U는 1.0e-12라 U만으로는 못 잡는다) |
| M3 표면 하중 규칙 차수 2 (p = 2) | `test_v03_...[2-traction]` | 검출 (f 2.8e-2, U 2.0e-5) |
| M4 절단 셀 QP 하나 제거 | `test_v02_...[2-drop]` | 검출 (QP 수 1794 ≠ 1795) |
| M5 확장 θ 0.2 → 0.48 (Kratos 측) | `test_v02_...[1-drop]` | 검출 (slave 수 불일치, U 9.9e-2) |
| M6 ghost α 1e-8 → 1e-7 (Kratos 측) | `test_v02_...[1-ghost]` | 검출 (K 1.3e-7, 역방향 3.9e-9, U 6.9e-7) |
| M7 ν 0.49999 → 0.4999 (Kratos 측) | `test_v07_...[0.49999]` | 검출 (K 0.90) |
| M8 참조 측 W 재료 E 4e5 → 4.4e5 | `test_v04_two_part_tie` | 검출 (U 4.0e-2) |

## 7. 범위와 한계

- V02 축소판은 Roark 판의 실제 형상(STEP/OCP 테셀레이션)이 아니다. 모서리 1/4 원 노치 판의 해석적 삼각형 soup이다. 비교 대상인 코드 경로(`build_trimmed_table`의 절단 셀 QP 표, OUT 셀 ghost/drop, 확장 MPC, 박스 면 하중·고정)는 V02 파이프라인과 같다. soup 부피는 해석값과 1e-3 안에서 맞는지 단언한다.
- 축소판 크기: V01 2×2×2, V03 1×1×4 (계획의 가장 성긴 h), V07 4×4×1, V02 4×4×1, V04 2×(2×2×8), V05 2×2×7, V09 2×2×10, V10 4×1×1, V11 11×11×1, V12 10×2×1, V-EXP 10×1×1.
- 하중 벡터 비교는 `run_linear` 경로(V01, V02, V03, V07, 롤러)에만 있다. V04/V05는 참조 하중을 점하중으로 넣고, V09–V12·V-EXP는 하중 벡터를 같은 입력으로 쓰거나(V10, V-EXP) 행렬만 비교한다.
- `_kratos_rhs`는 풀이 뒤 새 builder로 다시 조립한 값이다. 풀이 중 strategy가 쓴 RHS 벡터 자체는 아니다. strategy의 `GetSystemVector()`는 풀이 끝에 반력 계산(`CalculateReactions`)이 f − K u로 덮어쓰고, Newton strategy는 builder를 Python에 내보내지 않는다. 같은 요소·조건·MPC를 같은 방정식 번호로 조립한 값이며, 번호 일치는 단언한다.
- 열·접촉·F-bar·follower·RBE3·Nitsche는 비교하지 않았다. Python 단일 구현이거나 참조 조립기가 지원하지 않기 때문이다(§4 표).
