# P303 — Native FE 지원 및 IGA–FE 하이브리드 해석

## 1. 개요 및 구현 범위

Phase 3의 Native 유한요소(FE) 지원 및 IGA–FE 하이브리드 결합 해석(T-P3.3)을 구현한 `wht_iga/kratos_backend/native_fe.py` 기술 문서다.  
1차·2차 연속체 솔리드 요소 5종을 Kratos StructuralMechanics 요소로 매핑하고, 음의 야코비안(inverted element) 자동 플립 교정 및 순수 Kratos/페널티 결합 하이브리드 솔버를 제공한다.

### 관련 소스 파일
- `WHT_IGA/wht_iga/kratos_backend/native_fe.py`

---

## 2. 요소 매핑 및 절점 순서 규약

### 2.1 `KRATOS_FE_ELEMENT_MAP`
| 정규화 요소 토큰 | 절점 수 | Kratos 요소 클래스 |
|---|---|---|
| `hex8` | 8 | `SmallDisplacementElement3D8N` |
| `tet4` | 4 | `SmallDisplacementElement3D4N` |
| `tet10` | 10 | `SmallDisplacementElement3D10N` |
| `hex20` | 20 | `SmallDisplacementElement3D20N` |
| `wedge6` | 6 | `SmallDisplacementElement3D6N` |

### 2.2 Hex20 절점 재정렬 (`to_kratos_connectivity`)
- 표준/Abaqus 규약의 20절점 육면체(bottom 0–3, top 4–7, bottom midside 8–11, top midside 12–15, vertical midside 16–19)를 Kratos 3D20N 절점 순서로 재배열한다:
  - `_HEX20_TO_KRATOS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15]`
- 나머지 요소(`hex8`, `tet4`, `tet10`, `wedge6`)는 항등 순서를 유지한다.

---

## 3. 야코비안 부호 검증 및 자동 교정

| 함수 / 심볼 | 동작 |
|---|---|
| `check_element_jacobian(points, connectivity, elem_type)` | 요소별 대표 3-모서리 행렬식 $\det([\mathbf{x}_a - \mathbf{x}_0, \mathbf{x}_b - \mathbf{x}_0, \mathbf{x}_c - \mathbf{x}_0])$ 을 계산하여 양수 여부(`bool`)와 최소 행렬식 반환 |
| `ensure_positive_jacobian(points, connectivity, elem_type, auto_flip=True)` | 음의 야코비안 요소 탐지 시 `auto_flip=True` 이면 요소별 뒤집기 순열(`_FLIP_HEX8`, `_FLIP_TET4`, `_FLIP_TET10`, `_FLIP_HEX20`, `_FLIP_WEDGE6`)을 적용해 양의 체적 방향으로 교정하고 교정된 요소 인덱스 목록 반환; `auto_flip=False` 이면 예외 발생 |

---

## 4. Kratos 파트 생성 및 하이브리드 솔버 (`solve_hybrid_iga_fe`)

### 4.1 `populate_native_fe_part` 및 `NativeFEPartHandle`
- `AssemblyBuilder` 의 `IdOffsetTable` 또는 명시적 시작 ID를 기반으로 Kratos `ModelPart` 에 절점, 재료 `Properties`(`YOUNG_MODULUS`, `POISSON_RATIO`, `DENSITY`, `LinearElastic3DLaw`), 및 `KRATOS_FE_ELEMENT_MAP` 기반 솔리드 요소를 등록한다.
- 반환 객체 `NativeFEPartHandle` 은 `part_name`, `model_part`, `node_ids`, `element_ids`, `property_id`, `element_types`, `flipped_elements` 를 유지한다.

### 4.2 하이브리드 IGA–FE 해석 경로 (`solve_hybrid_iga_fe` -> `HybridSolveResult`)
- **경로 1 — 순수 Kratos MPC / 무페널티 경로 (`penalty_ties=None`)**:
  - `LinearMasterSlaveConstraint` (`mortar_l2`, `mpc_matching`, `RBE2`, `RBE3` 등)만 존재하는 경우 Kratos `ResidualBasedBlockBuilderAndSolver` + `ResidualBasedLinearStrategy` (기본 `PardisoLDLTSolver` / `PardisoLUSolver` / `SkylineLUFactorizationSolver`)를 직접 호출해 해석한다.
- **경로 2 — 페널티 결합 하이브리드 경로 (`penalty_ties` 존재 시)**:
  - Kratos 시스템 행렬/우변 벡터 `(K_base, f_base)` 와 구속 변환 행렬 `T_mpc` 를 추출하고, `penalty_ties` 로부터 생성된 대칭 반양정치 페널티 행렬 `K_pen` 을 합산해 $K_{\text{full}} = K_{\text{base}} + \sum K_{\text{pen}}$ 을 구성한다.
  - MPC 구속 축소 $K_{\text{red}} = T_{\text{mpc}}^T K_{\text{full}} T_{\text{mpc}}$, $f_{\text{red}} = T_{\text{mpc}}^T f_{\text{base}}$ 수행 후 Kratos `IsFixed()` 디리클레 경계조건 파티셔닝(`free_dofs` vs `fixed_dofs`)을 적용해 `scipy.sparse.linalg.spsolve` 로 해를 구하고, $u_{\text{full}} = T_{\text{mpc}} u_{\text{red}}$ 로 전체 자유도를 복원해 Kratos 절점 `DISPLACEMENT` 에 기록한다.
- `HybridSolveResult` 는 전체 변위 벡터 `u_full`, 절점별 `displacements`, `mpc_count`, `penalty_tie_count`, `linear_solver`, `timings` 를 반환한다.

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_phase3_assembly_interactions.py` | `hex20` Kratos 순열 매핑 및 5종 요소 음의 야코비안 자동 플립 검증 포함 | **7 passed** (1.50s) | 어셈블리·상호작용 단위 테스트와 함께 실행 |
| `WHT_IGA/tests/kratos_spikes/test_p3_assembly_hybrid_thb.py` | `test_p3_3_native_fe_all_types_and_hybrid_iga_fe_v04` (5종 FE 요소 단일 패치 + V04 IGA/FE 하이브리드 보 해석) | **1 passed** (전체 5건 통과) | Kratos 10.4.3 스파이크 검증 |
