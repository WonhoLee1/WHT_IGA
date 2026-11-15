# P302 — Kratos 어셈블리 빌더, ID 오프셋 테이블 및 상호작용 구속

## 1. 개요 및 구현 범위

Phase 3 다중 파트 어셈블리(T-P3.2)는 파트를 `Structure` 루트 `ModelPart` 하나로 통합 컴파일하는 `AssemblyBuilder`, 컨테이너별 비충돌 1-based ID 할당기 `IdOffsetTable`, 그리고 비정합/정합 인터페이스 구속(`Tie`, `RBE2`, `RBE3`)을 제공하는 `interactions.py` 로 구성된다.

### 관련 소스 파일
- `WHT_IGA/wht_iga/kratos_backend/assembly_builder.py` — `PartSpec`, `InstanceSpec`, `AssemblySpec`, `IdRangeEntry`, `IdOffsetTable`, `AssemblyBuilder`, `validate_assembly_spec`, `validate_rotation_matrix`, `transform_stress_voigt`
- `WHT_IGA/wht_iga/kratos_backend/interactions.py` — 구속 그래프 검증/평탄화, `SurfaceMesh`, `TieSpec`/`build_tie`/`apply_tie`, `RBE2Spec`/`build_rbe2`/`apply_rbe2`, `RBE3Spec`/`build_rbe3`/`distribute_rbe3_loads`/`apply_rbe3`

---

## 2. 단일 루트 `Structure` 어셈블리 및 ID 테이블 (`assembly_builder.py`)

### 2.1 설계 제약 (S15)
- **단일 루트 `ModelPart("Structure")`**: Kratos `Structure` 단일 루트 모델파트 하위에 `Parts.<part>_<inst>`, `Surfaces`, `Couplings` 서브모델파트 계층을 구성하여 다중 파트 간 구속 및 일괄 해석이 가능하게 한다.
- **인스턴스 독립성/종속성 규칙 (`validate_assembly_spec`)**: 동일 `PartSpec` 을 참조하는 모든 `InstanceSpec` 은 `dependent` 플래그가 일치해야 한다. 종속(`dependent=True`) 및 독립(`dependent=False`) 인스턴스 모두 아핀 강체 변환 $\mathbf{x}' = R\mathbf{x} + \mathbf{t}$ 와 직교 회전 행렬 검증(`validate_rotation_matrix`: $R^T R \approx I, \det R \approx +1$)을 거친다.
- **응력 텐서 회전 (`transform_stress_voigt`)**: Voigt 6성분 응력 `[sxx, syy, szz, sxy, syz, sxz]` 을 대칭 행렬로 변환해 $\boldsymbol{\sigma}' = R \boldsymbol{\sigma} R^T$ 회전 후 다시 Voigt 6성분으로 복원한다.

### 2.2 `IdOffsetTable` 과 감사 추적
- **모노토닉 1-based 연속 범위 할당**: 5개 컨테이너(`nodes`, `elements`, `conditions`, `properties`, `constraints`) 각각에 대해 인스턴스별·상호작용 블록별 비충돌(disjoint) `[start_id, end_id]` 범위를 기록한다(`IdRangeEntry`).
- **직렬화 스키마**: `IdOffsetTable.to_dict()` / `write_resolved_json()` 는 스키마 버전 `"wht_iga.resolved_assembly.v1"` 으로 `<job>.resolved.json` 을 생성해 모든 엔티티 ID 출처를 재현·감사할 수 있다.

---

## 3. 구속 그래프 검증 및 컴파일 (`interactions.py`)

### 3.1 `LinearRelation` 및 S04 제약 검증
- `DofKey(node_id, var_name)` 과 `LinearRelation(slave, masters, weights, constant=0.0)` 로 다중점 구속을 추상화한다.
- `validate_constraint_graph(relations)` 는 다음 4가지 결함을 검사해 `ConstraintGraphReport` 를 생성한다:
  1. **S04 비영 상수항 (`nonzero_constants`)**: `constant != 0.0` 인 경우 차단(Kratos `LinearMasterSlaveConstraint` 와 `ResidualBasedBlockBuilderAndSolver` 의 호환성 규칙).
  2. **이중 종속절점 (`double_slaves`)**: 동일 `DofKey` 가 둘 이상의 관계에서 slave 로 지정된 경우 차단.
  3. **연쇄 구속 (`chains`)**: 한 관계의 slave 가 다른 관계의 master 로 등장하는 경우 탐지.
  4. **순환 구속 (`cycles`)**: 3색(white/gray/black) DFS 로 순환 종속성 탐지.
- `compile_constraint_graph(relations, chain_policy="reject"|"flatten")`:
  - `chain_policy="flatten"` 선택 시 순환이 없는 DAG 연쇄 구속을 치환 평탄화(flattening)하여 모든 master 를 루트 독립 DOF 로 환원한다.
  - `apply_relations_to_model_part` 는 Kratos `MasterSlaveConstraints` 컨테이너에 `LinearMasterSlaveConstraint` 를 단조 ID 로 생성한다.

---

## 4. 표면 매칭 Tie, RBE2 및 RBE3 (`interactions.py`)

### 4.1 표면 메쉬 추상화 (`SurfaceMesh`)
- `SurfaceMesh.from_fe(...)`: `hex8`, `hex20`, `tet4`, `tet10`, `wedge6` 외면 삼각화 지원. 2차 요소(`hex20`, `tet10`)는 중간절점을 포함한 4분할/8분할 서브-삼각형 분할과 형상함수 기반 `_build_fe_surface_evaluator` 를 구축한다.
- `SurfaceMesh.from_iga_face(...)`: IGA 경계면 샘플링 및 B-spline 면 평가기 생성.

### 4.2 Tie 구속 (`TieSpec`, `build_tie`, `apply_tie`)
| 기법 (`method`) | 수학적 정식화 및 동작 |
|---|---|
| `auto` | 절점-대-절점 정합 여부(허용오차 내 1:1 대응)에 따라 `mpc_matching` 또는 `mortar_l2` 자동 선택 |
| `mpc_matching` | 1:1 최근접 매칭 또는 삼각형 바리센트릭(barycentric) 보간 가중치 기반 `LinearRelation` 생성 |
| `mortar_l2` | 슬레이브·마스터 투영 삼각형을 2D 접평면에서 Sutherland-Hodgman 클리핑(`_clip_polygon_2d`) 후 양의 가중치 가우스 구적점으로 적분하여 $M_{ss} P = C_{sm}$ 을 풀어 L2 투영 연산자 $P$ (행 합 = 1) 및 MPC 관계 생성 |
| `penalty_spring` | 강성 스케일 $\alpha = \beta_{\text{pen}} E_{\max} / h$ 로 차분 행렬 $D$ 를 구성하고, 대칭 반양정치 행렬 $K_{\text{scalar}} = D^T W_\alpha D$ 를 `kron(K_scalar, I_3)` 으로 3D 벡터 DOF 로 확장 |

### 4.3 강체 커플링 RBE2 (`RBE2Spec`, `build_rbe2`, `apply_rbe2`)
- 기준절점 $\mathbf{x}_{\text{ref}}$ 의 병진 $\mathbf{u}_{\text{ref}}$ 과 회전 $\boldsymbol{\theta}_{\text{ref}}$ 로부터 각 종속절점 $\mathbf{x}_i$ 의 변위를 선형화 강체 운동학으로 구속한다:
  $$\mathbf{u}_i = \mathbf{u}_{\text{ref}} + \boldsymbol{\theta}_{\text{ref}} \times (\mathbf{x}_i - \mathbf{x}_{\text{ref}})$$
- 각 슬레이브 성분 $u_{i,k}$ 를 `u_ref[k]` 및 `ROTATION` 성분들의 `LinearRelation` (`constant = 0.0`)으로 변환한다.

### 4.4 분포 커플링 RBE3 (`RBE3Spec`, `build_rbe3`, `distribute_rbe3_loads`, `apply_rbe3`)
- 상대 위치 벡터 $\mathbf{d}_i = \mathbf{x}_i - \mathbf{x}_{\text{ref}}$ 와 가중치 $w_i$ (`weight_mode`: `"uniform"`, `"area"`, `"inverse_distance"`)를 사용해 최소제곱 관성 행렬을 구성한다:
  $$J = \sum_i w_i \left((\mathbf{d}_i \cdot \mathbf{d}_i) I_3 - \mathbf{d}_i \mathbf{d}_i^T\right)$$
- **두 가지 모드**:
  - `"mpc"` 모드: $\mathbf{u}_{\text{ref}}, \boldsymbol{\theta}_{\text{ref}}$ 를 slave, 커플링 면 절점을 master 로 두는 가중 최소제곱 `LinearRelation` 생성(선택적 `slave_basis_mat` 통해 B-spline 제어점 공간으로 역변환 가능).
  - `"load_distribution"` 모드: 기준절점에 작용하는 등가 힘/모멘트 $(\mathbf{F}, \mathbf{M})$ 를 정역학적으로 등가인 절점 하중 $\mathbf{f}_i = w_i \left(\frac{\mathbf{F}}{\sum w_k} + \mathbf{d}_i \times (J^{-1}\mathbf{M}_{\text{eff}})\right)$ 로 분배한다(IGA 제어점 하중 변환 지원).

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_phase3_assembly_interactions.py` | ID 오프셋 비중복, 구속 그래프(순환/이중슬레이브/평탄화), `mortar_l2`, `penalty_spring`, `RBE2`, `RBE3` | **7 passed** (1.50s) | `native_fe` 단위 테스트 포함 |
| `WHT_IGA/tests/unit/test_constraints_guard.py` | S04 구속 조건 위반 가드 검증 | 실행 예정 (unit 스위트 내 포함) | 기존 구속 가드 회귀 검증 |
| `WHT_IGA/tests/kratos_spikes/test_p3_assembly_hybrid_thb.py` | `test_p3_2_dependent_instances_rotation_translation_invariance`, `test_p3_2_nonmatching_mortar_rbe2_and_rbe3_v03_moment` | **2 passed** (스파이크 전체 5건 통과) | V03 굽힘 모멘트 RBE3 포함 |
