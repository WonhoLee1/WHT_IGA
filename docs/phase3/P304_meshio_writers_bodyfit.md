# P304 — 다중 포맷 Mesh I/O (`meshio_ext`) 및 Body-Fitted B-Spline 변환

## 1. 개요 및 구현 범위

Phase 3의 범용 메쉬 입출력 계층(`wht_iga/meshio_ext/`)과 구조화된 육면체 메쉬를 Body-Fitted B-Spline/NURBS 패치로 변환하는 `bodyfit.py` 기술 문서다(T-P3.4).  
패키지 임포트 시 **Kratos 무의존(Kratos-free)** 설계를 유지하여 전처리·변환 단계에서 Kratos 런타임 없이 동작하며, K26 규약에 따라 모든 리더는 `unit_system` 인자를 필수로 요구한다.

### 관련 소스 파일
- `WHT_IGA/wht_iga/meshio_ext/__init__.py` — 통합 디스패처 `read_mesh`, `write_mesh`
- `WHT_IGA/wht_iga/meshio_ext/model.py` — `MeshModel`, `ElemBlock`, 정규화/절점순서 테이블, 외면 추출(`outer_faces`)
- `WHT_IGA/wht_iga/meshio_ext/read_inp.py`, `write_inp.py` — Abaqus `.inp`
- `WHT_IGA/wht_iga/meshio_ext/read_nastran.py`, `write_nastran.py` — Nastran `.bdf/.nas` 및 OptiStruct `.fem`
- `WHT_IGA/wht_iga/meshio_ext/read_k.py`, `write_k.py` — LS-DYNA `.k/.key/.dyn`
- `WHT_IGA/wht_iga/meshio_ext/read_rad.py`, `write_rad.py` — Altair RADIOSS `.rad`
- `WHT_IGA/wht_iga/meshio_ext/read_msh.py`, `write_msh.py` — Gmsh MSH 2.2 `.msh`
- `WHT_IGA/wht_iga/meshio_ext/bodyfit.py` — Body-Fitted 멀티패치 NURBS 변환 및 검증

---

## 2. 정규화 모델 및 디스패처 (`model.py`, `__init__.py`)

### 2.1 `MeshModel` 및 정규화 규약
- **지원 3D 솔리드 정규 토큰 (`CANONICAL_COUNTS`)**: `hex8` (8), `tet4` (4), `tet10` (10), `hex20` (20), `wedge6` (6).
- **K26 단위계 강제 (`require_unit_system`)**: 단위계 미지정 시 예외 발생, `MeshModel` 메타데이터에 단위계를 보존.
- **절점 순서 변환 (`NODE_ORDER`, `apply_order`)**: 포맷별 요소 절점 인덱싱을 표준 규약으로 통일하여 라운드트립 중 요소 뒤집힘을 방지.
- **저차원/비솔리드 필터링 (`NON_SOLID_TOKENS`)**: 쉘/빔/트러스 등 비솔리드 요소는 제외하고 `NON_SOLID_EXCLUDED` 경고 이슈를 기록.
- **외면 추출 (`outer_faces`)**: `HEX8_FACES`, `TET4_FACES`, `WEDGE6_FACES` 정의를 이용해 체적 메쉬의 경계 외면 삼각/사각형을 추출.
- **직렬화**: `MeshModel.to_dict()`, `from_dict()`, `roundtrip()`.

### 2.2 통합 읽기/쓰기 디스패치 (`__init__.py`)
| 확장자 | 리더 | 라이터 | 주요 특징 |
|---|---|---|---|
| `.inp` | `read_inp` | `write_inp` | `*NODE`, `*ELEMENT`, `*NSET`/`*ELSET` (`GENERATE` 포함), `*SURFACE`, `*PART`/`*END PART`, `*SOLID SECTION` |
| `.bdf`, `.fem`, `.nas` | `read_nastran` | `write_nastran` | `pyNastran` 기반 파싱 + OptiStruct `SET`/`SURF` + `$$HMNAME COMP/SETS`; `format_dialect="auto"\|"nastran"\|"optistruct"`, `GRID*`, `PSOLID`/`MAT1`, `THRU` 압축 |
| `.k`, `.key`, `.dyn` | `read_k` | `write_k` | `*INCLUDE` 순환 참조 가드(`INCLUDE_CYCLE`), `*PART`, `*SET_NODE_LIST`, `*SET_SOLID`, `*SET_SEGMENT` |
| `.rad` | `read_rad` | `write_rad` | `/NODE`, `/BRICK`, `/TETRA4`, `/TETRA10`, `/PENTA6`, `/GRNOD`, `/GRBRIC`, `/SURF` |
| `.msh` | `read_msh` | `write_msh` | `meshio` 기반 MSH 2.2 ASCII + `$Comments` JSON(`wht_iga_v=1`)에 파트/셋/메타데이터 무손실 저장 |

---

## 3. Body-Fitted B-Spline/NURBS 변환 (`bodyfit.py`)

### 3.1 공개 심볼
- `BodyFitValidationError`
- `PatchInterface`
- `BodyFitSplinePatch` (`control_points_flat`, `weights_flat`, `kratos_knots`, `evaluate`)
- `BodyFitConversionResult`
- `mesh_to_bodyfit_patches(model, continuity="c0"|"smooth")`
- `build_kratos_nurbs_volumes(result, KM, model_part, ...)` — `KM` 을 인자로 주입받아 임포트 시 Kratos 무의존 유지

### 3.2 위상 검증 및 구조화 패치 추출 알고리즘
`mesh_to_bodyfit_patches` 는 Hex 메쉬가 텐서곱 B-Spline 패치 집합으로 표현 가능한지 엄격히 검증한다:
1. **요소/절점 무결성**: 혼합 요소 타입 거부, 중복 절점 거부.
2. **다양체 및 내점/내변 밸런스(valence) 검사**:
   - 공유면 3개 이상인 비다양체 면(non-manifold face) 거부.
   - 내부 절점 밸런스 $\ne 8$, 내부 모서리 밸런스 $\ne 4$ 거부.
   - **T-junction / hanging node 탐지 (`_check_t_junctions`)**: 비정합 면 분할 불일치 차단.
3. **연결 성분별 구조화 격자 전파 (`_propagate_component_grid`)**:
   - 각 연결 성분에 대해 $(n_x, n_y, n_z)$ 직교 인덱스 좌표계를 BFS 로 전파한다.
   - 비구조적 육면체 메쉬인 경우 `ValidationIssue(code="UNSTRUCTURED_HEX_UNSUPPORTED")` 를 발행하고 `BodyFitValidationError` 로 안전하게 차단한다.
4. **연속성 및 멀티패치 인터페이스**:
   - `continuity="c0"`: 요소 경계마다 $C^0$ 보간 매듭 구성.
   - `continuity="smooth"`: 내부 매듭 중복도를 줄여 $C^{p-1}$ 매끄러운 B-Spline 패치 구성.
   - 인접 패치 간 경계 제어점 일치 쌍을 `PatchInterface` 로 추출하며, `build_kratos_nurbs_volumes` 를 통해 Kratos `NurbsVolumeGeometry` 로 조립할 수 있다.

---

## 4. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_meshio_model.py` | `MeshModel` 정규화, 외면, 직렬화 | **6 passed** (1.33s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_meshio_inp.py` | Abaqus `.inp` 읽기/라운드트립 | **3 passed** (1.38s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_meshio_nastran.py` | Nastran/OptiStruct `.bdf/.fem` | **3 passed** (1.43s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_meshio_k.py` | LS-DYNA `.k` 및 `*INCLUDE` 가드 | **3 passed** (1.16s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_meshio_rad.py` | RADIOSS `.rad` 읽기/셋 | **3 passed** (0.81s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_meshio_msh.py` | Gmsh MSH 2.2 및 `$Comments` 메타 | **3 passed** (0.92s) | Kratos 무의존 |
| `WHT_IGA/tests/unit/test_phase3_meshio_roundtrip_bodyfit_iges.py` | 전 포맷 라운드트립, `bodyfit` 유효/비구조 거부, IGES/Gmsh | **22 passed, 1 failed** (3.03s) | 실패 1건(`test_gmsh_surface_loop_safeguard_meshes_iges_mode0_3d`)은 호스트 `PATH` 의 `miniconda3\Scripts\gmsh.bat` 내 `No module named 'gmsh'` 외부 CLI 환경 오류에 기인 |
