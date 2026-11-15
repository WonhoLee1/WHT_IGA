# P305 — CAD I/O (IGES 치유), Gmsh CLI 래퍼, 엔티티 선택 및 GPL 격리

## 1. 개요 및 구현 범위

Phase 3의 CAD 파일 로드·치유(`occ_io.py`), GPL 라이선스 격리형 외부 Gmsh 메쉬 생성기(`gmsh_runner.py`), CAD/FEM 공통 기하·위상 선택 엔진(`selection/`), 그리고 패키지 격리 검증(`test_isolation.py`)을 다루는 문서다(T-P3.5).

### 관련 소스 파일
- `WHT_IGA/wht_iga/geometry/occ_io.py` — OpenCASCADE(OCP) STEP/IGES/BREP 로드 및 IGES Mode 0 치유 파이프라인
- `WHT_IGA/wht_iga/meshgen/gmsh_runner.py` — GPL 격리 외부 `gmsh` CLI 실행기 및 IGES Surface Loop 가드
- `WHT_IGA/wht_iga/selection/{sets.py, selector.py, measure.py, contact_check.py}` — 집합 연산, 기하/위상 쿼리, CAD↔FEM 매핑, 계측, 접촉 간극 검사
- `WHT_IGA/tests/test_isolation.py` — 패키지·GPL·지연 임포트 격리 테스트

---

## 2. CAD I/O 및 IGES Mode 0 치유 (`occ_io.py`)

### 2.1 지원 포맷 및 공개 API
- 확장자 그룹: `STEP_EXTENSIONS`, `IGES_EXTENSIONS`, `BREP_EXTENSIONS`, `SUPPORTED_CAD_EXTENSIONS`
- `HealingReport`: 수정 전/후 쉘·솔리드·면 개수, 수행 단계, 미봉합/개방 쉘 이슈 기록
- `CADGeometry`: OCP `TopoDS_Shape`, 치유 리포트, 체적/경계상자 등 메타데이터 래퍼
- `load_cad_file(path, heal=True, sewing_tol=1e-6, unit_system=None) -> CADGeometry`

### 2.2 IGES Mode 0 치유 파이프라인
1. **수밀 솔리드 고속 경로(Fast Path)**: 입력 형상이 이미 유효한 닫힌 솔리드이면 불필요한 재구성을 생략한다.
2. **열린 트림면(IGES Mode 0) 치유 경로**:
   - `ShapeFix_Shape` (기하/토폴로지 보정) → `BRepBuilderAPI_Sewing` (`sewing_tol` 기반 면 봉합) → `ShapeFix_Shell` → `BRepBuilderAPI_MakeSolid` → `ShapeFix_Solid`.
   - 생성된 솔리드의 부피(`brep_volume < 0`)가 음수인 경우 `.Reverse()` 로 외향 법선 방향을 교정한다.
   - 치유 후에도 열린 쉘이나 미봉합 면이 남으면 `OPEN_SHELL_REMAINS` 또는 `UNSEWN_FACE_REMAINS` 경고/오류 이슈를 기록한다.

---

## 3. GPL 격리 Gmsh CLI 래퍼 (`gmsh_runner.py`)

### 3.1 라이선스 격리 원칙
- **Python `gmsh` 바인딩 임포트 전면 금지**: `wht_iga` 코드베이스 어디에서도 `import gmsh` 를 사용하지 않으며, `find_gmsh_executable()` 로 외부 바이너리를 탐색해 `subprocess.run(["gmsh", script_geo, "-"], ...)` 로만 호출한다.

### 3.2 주요 API 및 IGES Surface Loop 세이프가드
- `GmshOptions`: 요소 타입(`tet4`, `tet10`, `hex8`), 출력 포맷(`msh2`, `msh4`), 메쉬 크기(`lc_min`, `lc_max`), `hex8` 선택 시 `Mesh.SubdivisionAlgorithm = 2` 적용.
- `GmshResult`: 생성된 `.msh` 경로, `.geo` 스크립트, 표준출력/오류, 구조화 이슈 목록(`GMSH_NOT_FOUND`, `CAD_NOT_FOUND`, `GMSH_TIMEOUT`, `GMSH_EXEC_ERROR`, `GMSH_FAILED`, `MSH_NOT_CREATED`).
- **IGES Mode 0 Surface Loop 세이프가드 (`build_geo_script`)**:
  - 트림면으로만 구성된 IGES 파일은 OpenCASCADE 병합 직후 볼륨(`#_wht_vols()`)이 `0` 일 수 있다.
  - 스크립트 내에 `If (#_wht_vols() == 0 && #_wht_surfs() > 0)` 분기를 자동 삽입하여 `Coherence; Surface Loop(...) = {...}; Volume(1) = {...};` 를 실행함으로써 3D 체적 메쉬 생성이 실패하지 않도록 보호한다.

---

## 4. 엔티티 선택, 계측 및 접촉 검사 (`wht_iga/selection/`)

### 4.1 `SelectionSet` 및 `SetRegistry` (`sets.py`)
- 7가지 `SetKind` 를 지원하며 Python 집합 연산자 `|` (합집합), `&` (교집합), `-` (차집합), `^` (대칭차), `~` (여집합)를 오버로드한다.
- `SetRegistry` 는 이름 기반 셋 등록 및 `to_dict()` / `from_dict()` JSON 직렬화를 제공한다.

### 4.2 `EntityQuery` 및 어댑터 (`selector.py`)
- `CadAdapter` (OCP 형상) 및 `MeshAdapter` (`MeshModel`) 에 대해 단일 체이닝 쿼리 인터페이스를 제공한다:
  - **필터**: `by_geom_type`, `within_box`, `within_sphere`, `within_cylinder`, `on_plane`, `between_planes`, `by_normal`, `by_area`, `by_length`, `nearest_to`, `flood_fill`.
  - **CAD ↔ FEM 변환**: `to_node_set()`, `to_element_set()`, `to_surface_set()`, `to_surface_tris()`, `measure()`.

### 4.3 기하 계측 및 접촉 간극 검사 (`measure.py`, `contact_check.py`)
- `measure_cad`, `measure_mesh` → `BBoxResult`, `MeasureResult` (체적, 표면적, 모서리 길이, 경계상자 치수 계산).
- `inspect_contact` → `ContactGapReport` (두 표면/바디 쌍 간 최소·최대 간극, 관통(penetration) 및 법선 정렬 상태 진단).

---

## 5. 격리 규칙 및 테스트 커버리지 (`test_isolation.py` 외)

`tests/test_isolation.py` 는 정규식, Python AST 정적 검사 및 서브프로세스 런타임 검사를 통해 `wht_iga` 패키지가 다음을 위반하지 않음을 검증한다:
1. `wht_topopt`, `aieng` 워크스페이스 패키지 임포트 금지.
2. GPL 라이브러리 `gmsh` Python 모듈 임포트 금지.
3. 패키지 최상위 임포트 시 `KratosMultiphysics` 즉시(eager) 임포트 금지(지연 `load_kratos()` 만 허용).

### 직접 실행 관찰 결과
| 테스트 파일 | 실행 결과 | 비고 |
|---|---|---|
| `WHT_IGA/tests/test_isolation.py` | **4 passed** (1.68s) | AST + 서브프로세스 격리 전항목 통과 |
| `WHT_IGA/tests/unit/test_geo_cad.py` | **5 passed** (1.69s) | STEP/IGES/BREP 로드 및 IGES 치유 통과 |
| `WHT_IGA/tests/unit/test_selection.py` | **5 passed** (2.77s) | OCP CAD↔FEM 선택·계측·접촉검사 통과 |
| `WHT_IGA/tests/unit/test_gmsh_runner.py` | **5 passed, 3 failed** (2.30s) | `.geo` 스크립트/옵션/에러처리 단위 테스트 5건 통과; 실제 외부 `gmsh` 실행 3건(`tet4`, `tet10`, `hex8`)은 현 호스트 `PATH` 의 `C:\Users\GOODMAN\miniconda3\Scripts\gmsh.bat` 이 `ModuleNotFoundError: No module named 'gmsh'` 로 종료되어 실패 |
| `WHT_IGA/tests/unit/test_phase3_meshio_roundtrip_bodyfit_iges.py` | **22 passed, 1 failed** (3.03s) | `test_gmsh_surface_loop_safeguard_meshes_iges_mode0_3d` 역시 위와 동일한 호스트 `gmsh.bat` 환경 미설치로 실패 |
