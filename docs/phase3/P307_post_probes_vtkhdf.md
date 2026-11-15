# P307 — 후처리 오차 노름, Proxy 메쉬, VTKHDF 2.1, 무가중치 프로브 및 결과 쿼리

## 1. 개요 및 구현 범위

Phase 3의 후처리(`wht_iga/post/`)와 통합 읽기 전용 결과 조회 API(`wht_iga/query/results.py`) 기술 문서다(T-P3.7).  
두 패키지 모두 임포트 시 Kratos 의존성 없이 동작하며(Kratos 프로브 평가 시에만 컨텍스트 주입), VTKHDF 2.1 표준과 B-spline 제어점 보존 확장을 결합해 ParaView 시각화 및 정량 검증을 동시에 지원한다.

### 관련 소스 파일
- `WHT_IGA/wht_iga/post/__init__.py` — 후처리 공개 API 재수출
- `WHT_IGA/wht_iga/post/norms.py` — $L_2$, $H^1$, 에너지 노름 및 해석해 적분 오차
- `WHT_IGA/wht_iga/post/proxy_mesh.py` — `IN` 육면체 세분화 + `CUT` 사면체 클리핑 Proxy 메쉬 생성기
- `WHT_IGA/wht_iga/post/vtk_writer.py` — VTKHDF 2.1 (`UnstructuredGrid` + `Steps` + `/WHT_IGA` 그룹) 라이터/리더
- `WHT_IGA/wht_iga/post/probes.py` — `PointProbe`, `LineProbe`, `MultiPointProbe`, `w=0.0` Kratos 프로브 요소, JSON/SQLite 내보내기
- `WHT_IGA/wht_iga/query/results.py` — `ResultView`, 변수 별칭 매핑, 평형/단면력/ASME 응력 선형화

---

## 2. 오차 노름 (`wht_iga/post/norms.py`)

| 심볼 | 설명 |
|---|---|
| `NormResult`, `H1NormResult` | 절대 오차, 기준 노름, 상대 오차, 적분 체적 저장 데이터클래스 |
| `gauss_legendre_hex_rule(n_gp)` | $[-1, 1]^3$ 상의 3D 텐서곱 Gauss-Legendre 구적점 및 가중치 반환 |
| `compute_l2_norm`, `compute_h1_norm`, `compute_energy_norm` | 수치장/해석장 차이의 $L_2$, $H^1$ (값+구배), 탄성 에너지 노름 계산 |
| `integrate_error_against_analytic` | 셀 경계 상자(`cell_bounds`) 기반 3D Gauss-Legendre 오차 적분 |
| `integrate_error_on_tets` | 사면체(`VTK_TETRA`) 메쉬 위에서 수치장 vs 해석해 오차 적분 |

---

## 3. Proxy 메쉬 및 VTKHDF 2.1 (`proxy_mesh.py`, `vtk_writer.py`)

### 3.1 Proxy 메쉬 생성 (`proxy_mesh.py`)
- `ProxyMesh` (`VTK_HEXAHEDRON`, `VTK_TETRA` 지원), `BSplineFieldEvaluator`.
- `generate_hex_proxy_for_cells(cell_bounds, n_sub=2, deduplicate_tolerance=1e-12)`: `IN` 셀을 `n_sub` 세분화한 육면체 블록 생성(인접 셀 공유 절점 중복 제거 지원).
- `generate_cut_tet_proxy(...)`: `CUT` 셀을 6-사면체 Kuhn 분할(`decompose_hex_to_6_tets`) 후 `inside_fn` 필터링 또는 `phi_fn <= 0` 음반공간 클리핑(`_clip_tet_negative_halfspace`)으로 경계 적합 사면체 블록 생성(또는 사전 계산된 `cut_points`/`cut_tets` 직접 수용).
- `merge_proxy_meshes`, `build_proxy_mesh`: `IN`/`CUT` 블록 병합(`BlockId` 셀 데이터 부여) 및 B-spline 장 평가.

### 3.2 VTKHDF 2.1 저장 및 제어점 보존 (`vtk_writer.py`)
- `write_vtkhdf`, `read_vtkhdf`:
  - 표준 `/VTKHDF` 그룹(`Type="UnstructuredGrid"`, `Version=[2, 1]`, `Steps` 과도/멀티프레임 그룹)을 생성하여 ParaView 에서 즉시 로드 가능.
  - 동시에 확장 그룹 `/WHT_IGA/block_names` 및 `/WHT_IGA/control_points/<frame>/<var>` 에 원시 B-spline 제어점 장을 보존하여 사후 재평가·변환 손실을 방지.

---

## 4. 물리 좌표 프로브 및 무가중치(`w=0`) Kratos 평가 (`probes.py`)

### 4.1 프로브 클래스 및 내보내기
- `PointProbe(name, point, variables)`: 단일 물리 좌표 평가.
- `LineProbe(name, start, end, n_points=21, variables)`: 균등 선분 좌표 및 누적 호장(arc-length) `s` 제공.
- `MultiPointProbe(name, points, variables)`: 임의 3D 점군 일괄 평가.
- `ProbeSample`: `get(var_name)` / `__getitem__` 으로 Coreform/Abaqus 코드와 Kratos 변수명 모두 조회 가능.
- `export_probes_json` (스키마 `"wht_iga.probes.v1"`), `export_probes_sqlite` (`probe_meta`, `probe_values` 테이블).

### 4.2 평가 백엔드 3종
1. **`BoxSpline` 직접 평가** (`mode="bspline"`)
2. **`THBHierarchy` 계층 기저 평가** (`mode="thb"`): `evaluate_thb_basis_and_grad` 기반 `U, E, S, MISES, NT11` 계산
3. **Kratos 무가중치(`w=0.0`) 프로브 요소 (`mode="zero_weight_kratos"`)**:
   - `create_zero_weight_probe_element` / `evaluate_zero_weight_kratos_probes`:
   - 물리 좌표를 모수 좌표로 역매핑하고 구적 가중치를 **`w = 0.0`** 으로 설정해 `SmallDisplacementElement3D4N` 을 생성한다.
   - `verify_zero_contribution=True` 일 때 `CalculateLocalSystem` 의 강성 `lhs` 와 하중 `rhs` 최댓값이 정확히 `0.0` 임을 검증하여 전역 계 행렬에 영향을 주지 않고 Kratos 구성방정식 기반 응력/변형률을 직접 추출한다.

---

## 5. 통합 결과 쿼리 (`wht_iga/query/results.py`)

### 5.1 변수 별칭 테이블 (`FIELD_ALIASES`, `resolve_field_name`)
- `U` ↔ `DISPLACEMENT`
- `S` ↔ `CAUCHY_STRESS_TENSOR`
- `E` / `STRAIN` ↔ `GREEN_LAGRANGE_STRAIN_TENSOR`
- `MISES` ↔ `VON_MISES_STRESS`
- `NT11` ↔ `TEMPERATURE`
- `RF` ↔ `REACTION`

### 5.2 `ResultView(vtkhdf_path, sqlite_path=None)` 주요 기능
| 메서드 | 동작 및 수학적 기법 |
|---|---|
| `field()`, `control_point_field()` | 점/셀/제어점 데이터(`VTKHDF` + `SQLite`) 통합 조회 및 `block` 필터링 |
| `probe(points, field_name)` | 사면체 분할(`_all_tets`) + 24-최근접 무게중심 `cKDTree` 탐색 후 바리센트릭(barycentric) 보간; 외부점은 4-NN 역거리 가중 보간 |
| `hotspots(field_name, top_n)` | 스칼라 크기 또는 `von_mises_stress` 기준 상위 `top_n` 위치 추출 |
| `reaction_force()`, `check_equilibrium()` | 전역 정역학 힘 평형 $\mathbf{F}_{\text{app}} + \mathbf{R}_{\text{sum}} = \mathbf{0}$ 검사(`EquilibriumReport`) |
| `section_force(plane_point, plane_normal)` | 절단 평면과 사면체 교차 다각형(삼각/사각) 추출 후 $\int_{\Gamma} \boldsymbol{\sigma} \cdot \mathbf{n}\,dA$ 단면력 적분 |
| `stress_linearize(p1, p2, n_samples=41)` | 복합 Simpson 1/3·3/8 구적 가중치(`_simpson_weights`)로 **ASME Sec VIII Div 2** 막(membrane)·굽힘(bending)·피크(peak) 응력 선형화 계산(`StressLinearizationResult`) |
| `error_norm(...)` | 해석해 대비 $L_2$/$H^1$ 오차 노름 계산 |

---

## 6. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_post_and_query.py` | Proxy 메쉬, VTKHDF 2.1 라운드트립, `ResultView` 프로브·평형·단면력·ASME 선형화·오차노름 | **5 passed** (1.26s) | Kratos 무의존 단위 테스트 |
| `WHT_IGA/tests/unit/test_phase3_steps_probes.py` | `PointProbe`/`LineProbe`/`MultiPointProbe`, THB/B-spline 평가, JSON/SQLite 내보내기 포함 | **7 passed** (1.10s) | 스텝 단위 테스트와 함께 실행 |
| `WHT_IGA/tests/unit/test_probes_selfconv.py` | 자기수렴(self-convergence) 프로브 검증 | **4 passed** (0.35s) | — |
| `WHT_IGA/tests/kratos_spikes/test_p3_solver_amgcl_log.py` | `test_kratos_zero_weight_probe_invariance` (`w=0` 요소 전역 행렬 불변 및 응력 일치) | **1 passed** (전체 5건 통과) | Kratos 10.4.3 스파이크 검증 |
