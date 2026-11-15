# P301 — THB-Splines 계층 구조 및 Multi-Cell 공유 B-Spline 필드

## 1. 개요 및 구현 범위

Phase 3의 국부 세분화(T-P3.1a)와 다중 재료/다중 영역 공유 스플라인(T-P3.1b)을 구현한 모듈 기술 문서다.  
두 경로 모두 단일 `NurbsVolumeGeometry` 또는 계층적 세분화 트리를 기반으로 Kratos `Structure` 모델에 정밀하게 통합되며, 기존 안정화/확장 연산자(`E`, `T = C_THB @ E`)와 직접 결합된다.

### 관련 소스 파일
- `WHT_IGA/wht_iga/grid/hierarchy.py` — 1D/3D Oslo 세분화 행렬, 2:1 균형, 활성 DOF 판정, `THBHierarchy`, `C_THB` 조립
- `WHT_IGA/wht_iga/operators/compose.py` — 희소 행렬 기반 연산자 합성(`compose_thb`, `expand_dof_operator`, `apply_TKT`, `apply_Tt_vec`, `solve_reduced`)
- `WHT_IGA/wht_iga/kratos_backend/assembly_builder.py` — `CellDomainSpec`, `_subdivide_box_by_domains`, `_build_iga_instance`, `_build_thb_instance`, `attach_thb_master_constraints`, `assemble_and_solve_thb_reduced`, `compose_thb_operator`

---

## 2. Path A: THB-Splines (`THBHierarchy`)

### 2.1 공개 클래스 및 함수 (`wht_iga/grid/hierarchy.py`)
| 심볼 | 역할 |
|---|---|
| `midpoint_refine_knots(knots, degree)` | 비어 있지 않은 각 매듭 구간(knot span) 중점을 삽입해 이분(dyadic) 세분화 매듭 벡터 생성 |
| `subdivision_matrix_1d(coarse_knots, degree)` | Oslo Algorithm 1 기반 1D 세분화 행렬 $S \in \mathbb{R}^{n_f \times n_c}$ (또는 전치 규약 $N^L = S^T N^{L+1}$) 생성 |
| `subdivision_matrix_3d(su, sv, sw)` | 텐서곱 세분화 행렬 `kron(sw, kron(sv, su))` (`I = i + nu * (j + nv * k)` DOF 순서 일치) |
| `LeafCell` | 계층 트리의 리프 셀(`level`, `ijk`, `lo`, `hi`, `bounds`) 표현 데이터클래스 |
| `THBHierarchy` | 다단계 이분 B-spline 계층 관리자 |

### 2.2 수학적 정식화
1. **Oslo Algorithm 1 (1D/3D 세분화)**
   - 거친 매듭 벡터 $\Xi^l$ 를 중점 세분화한 미세 매듭 벡터 $\Xi^{l+1}$ 에 대해 $N^l_j(x) = \sum_i S_{i,j} N^{l+1}_i(x)$ 가 기계 정밀도로 성립하도록 $S$ 를 구성한다.
   - 3D 텐서곱 제어점 평탄화 인덱스 $I = i + n_u (j + n_v k)$ 에 맞추어 $S_{3\text{D}} = S_w \otimes S_v \otimes S_u$ 를 구성한다.
2. **2:1 균형 강제 (`enforce_2_to_1_balance`)**
   - 26-이웃(면·모서리·꼭짓점 공유) 거친 셀을 기준으로 BFS 전파를 수행하여 인접 리프 셀 간 레벨 차가 최대 1($|l_a - l_b| \le 1$)이 되도록 보증한다.
3. **활성 DOF 판정 (`active_dofs_at_level`)**
   - 각 레벨 $l$ 의 영역 마스크 $\Omega^l$ 에 대해 3D 누적합 테이블(`_prefix_sum_3d`, `_box_all_true`)을 사용해 기저 함수 지지 영역(support)을 $O(1)$ 로 검사한다.
   - Giannelli et al. (2012) 규정에 따라 $\operatorname{supp}(N^l_j) \subseteq \Omega^l$ 이고 $\operatorname{supp}(N^l_j) \not\subseteq \Omega^{l+1}$ 인 기저만 레벨 $l$ 의 활성 THB 기저로 선택한다.
4. **절단 행렬 조립 (`build_c_thb_matrix`)**
   - 최상위(가장 거친) 레벨에서 항등/활성 선택으로 시작해 미세 레벨로 내려가며 다음 점화식을 적용한다:
     $$C^{(l+1)} = M_{l+1}\left(S_{3\text{D}}^{(l)} C^{(l)}\right) \;\Big|\; P_{l+1}$$
   - 여기서 $M_{l+1}$ 은 $\Omega^{l+1}$ 내부 미세 기저에 해당하는 행을 0으로 절단(truncation)하여 자식 기저의 기여를 제거하고, $P_{l+1}$ 은 레벨 $l+1$ 에서 새로 활성화된 THB DOF 열을 병합한다.
   - 최종 출력은 형상 `(n_fine_dofs, n_thb_dofs)` 인 SciPy CSR 행렬이며, 임의의 2:1 균형 메쉬에서 **단위 분할($\sum_i N^{\text{THB}}_i(\mathbf{x}) = 1$)** 과 선형 다항식 재현을 만족한다.
5. **기저/구배 평가 및 투영**
   - `evaluate_thb_basis(points)` / `evaluate_thb_basis_and_grad(points)`: 최미세 B-spline 격자에서 기저 및 공간 구배를 계산한 뒤 `C_THB` 로 축소한다.
   - `project_coarse_coefficients(u_coarse)`: 거친 격자 계수를 THB 공간으로 투영한다.

---

## 3. 연산자 합성 (`wht_iga/operators/compose.py`)

THB 절단 행렬 $C_{\text{THB}}$ 와 안정화/확장 연산자 $E$ 를 희소 행렬 연산만으로 합성한다(밀집 행렬 변환 금지: `_require_sparse`).

| 심볼 | 설명 |
|---|---|
| `extension_matrix(n_total, active_indices, ties)` | 비활성 제어점을 활성 제어점의 선형 결합으로 확장하는 스칼라 행렬 $E$ 조립 |
| `compose_thb(E, C_thb=None)` | $C_{\text{THB}}$ 가 주어지면 `C_thb @ E` (직사각 `(n_fine_dofs, n_thb_dofs)` 허용) 반환, 없으면 `E` 반환 |
| `expand_dof_operator(T_basis, n_comp=3)` | 스칼라 기저 연산자 $T_{\text{basis}}$ 를 `sp.kron(T_basis, sp.eye(n_comp))` 로 벡터 DOF 연산자 $T_{\text{dof}}$ 로 확장 |
| `apply_TKT(K, T)` | 축소 강성 행렬 $K_{\text{red}} = T^T K T$ 계산 |
| `apply_Tt_vec(T, v)` | 축소 하중/벡터 $f_{\text{red}} = T^T f$ 계산 |
| `solve_reduced(K, f, T, solver="pardiso")` | $K_{\text{red}} u_{\text{red}} = f_{\text{red}}$ 해석 후 $u_{\text{full}} = T u_{\text{red}}$ 복원 (`pardiso_solve` 또는 `spla.spsolve`) |

`assembly_builder.py` 의 `compose_thb_operator(E, C_thb)` 역시 동일한 규약으로 직사각 `C_thb` 에 대해 `C_thb @ E` 를 계산한다.

---

## 4. Kratos 통합 (`assembly_builder.py`)

### 4.1 THB 인스턴스 조립 (`_build_thb_instance`) 및 두 가지 해석 경로
- `_build_thb_instance` 는 등방성 차수(`pu == pv == pw`)를 요구하며, 최미세 레벨 $L_{\max}$ (분할수 $2^{L_{\max}}$) 의 `BoxSpline` 및 `NurbsVolumeGeometry` 를 Kratos `ModelPart` 에 생성한다.
- 각 `LeafCell` 은 레벨 차 $L_{\max} - \text{leaf.level}$ 에 따라 최미세 격자의 부구간(sub-span) 집합으로 매핑되며, 각 미세 셀 구간마다 가우스 구적점(`SmallDisplacementElement3D4N`)을 생성해 적분 정밀도를 보장한다.
- **경로 1 — Kratos MPC 구속 (`attach_thb_master_constraints`)**:
  - `thb_masters` 서브모델파트에 가상 마스터 절점을 생성하고, 미세 제어점(slave)에 대해 `T = C_thb @ E` 행별 가중치로 `LinearMasterSlaveConstraint` 를 부착한다.
  - **S04 제약 준수**: Kratos `LinearMasterSlaveConstraint` 의 `constant` 는 반드시 `0.0` 으로 강제(`constant=None` 전달)한다.
- **경로 2 — 축소계 직접 조립/해석 (`assemble_and_solve_thb_reduced`)**:
  - 최미세 계층 모델에서 $K_{\text{hier}}, f_{\text{hier}}$ 를 추출하고 $T_{\text{dof}}^T K_{\text{hier}} T_{\text{dof}} u_{\text{red}} = T_{\text{dof}}^T f_{\text{hier}}$ 를 풀어 디리클레(Dirichlet) 구속과 함께 `u_red`, `u_fine` 을 복원한다.

### 4.2 Path B: Multi-Cell 공유 B-Spline 필드 (`CellDomainSpec`)
- V04 이종 금속 보(bimetallic beam), V05 분할 큐브, E03 다중 영역 봉지 구조처럼 **물성이 다른 인접 셀이 단일 연속 스플라인 제어점 격자를 공유**하는 경우에 사용한다.
- `CellDomainSpec(domain_id, E, nu, density, constitutive_law, bounds, cells, predicate)` 를 통해 영역별 `Properties` ID를 분리 지정한다.
- `_build_iga_instance` 와 `_subdivide_box_by_domains` 는 단일 `NurbsVolumeGeometry` 를 유지하므로 영역 경계에서 추가 인터페이스 구속 없이도 $C^{p-1}$ 연속성을 자연 보존하며, 각 셀이 속한 도메인에 따라 서로 다른 재료 `Properties` 를 `SmallDisplacementElement3D4N` 에 바인딩한다.

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_grid_hierarchy.py` | Oslo 1D/3D, 2:1 균형, 활성 DOF, `C_THB` 단위분할·선형재현·구배 | **19 passed** (1.16s) | Kratos 비의존 단위 테스트 |
| `WHT_IGA/tests/unit/test_operators_compose.py` | `compose_thb`, `expand_dof_operator`, `apply_TKT`, `solve_reduced` | **6 passed** (3.39s) | Kratos 비의존 단위 테스트 |
| `WHT_IGA/tests/kratos_spikes/test_p3_assembly_hybrid_thb.py` | `test_p3_1_thb_hierarchical_solve_and_cantilever`, `test_p3_1_multicell_v04_bimetallic_and_v05_partitioned` | **2 passed** (해당 2건 포함 전체 5건 모두 통과) | Kratos 10.4.3 스파이크 검증 |
| `WHT_IGA/tests/kratos_spikes/test_p2stab_compose_kratos.py` | `test_kratos_tkt_matches_reference_1e10`, `test_kratos_1e5_creation_timing` | **2 passed** | $T^T K T$ 및 $10^5$ 생성 시간 검증 |
