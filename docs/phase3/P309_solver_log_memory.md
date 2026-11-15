# P309 — 메모리 기반 솔버 자동 선택, 6단계 솔버 로그 및 검증 게이트

## 1. 개요 및 구현 범위

Phase 3의 호스트 RAM 탐색·솔버 메모리 추정·자동 전환(`memory.py`), 상용 CAE 수준의 6단계 ASCII/UTF-8 솔버 진단 로그(`solver_log.py`), 실행 전 차단 검증 게이트(`validate.py`), CSR/Kratos 선형 솔버 백엔드(`kratos_backend/solvers.py`), 그리고 Phase-1/통합 파이프라인 연동(`pipeline.py`)을 다루는 문서다(T-P3.9).

### 관련 소스 파일
- `WHT_IGA/wht_iga/solver/memory.py` — RAM 프로빙, Pardiso/AMGCL 메모리 추정식, `select_linear_solver`
- `WHT_IGA/wht_iga/solver/solver_log.py` — 6단계 `SolverLogger` 및 `WHT_SOLVER_LOG` 환경 제어
- `WHT_IGA/wht_iga/solver/validate.py` — `JobSpec`, `ValidationIssue`, `ValidationError`, `validate`
- `WHT_IGA/wht_iga/kratos_backend/solvers.py` — `build_kratos_linear_solver`, `solve_linear_csr`, `compute_relative_residual`
- `WHT_IGA/wht_iga/solver/pipeline.py` — `run_linear` 6단계 로깅 및 자동 솔버 선택 통합

---

## 2. 호스트 RAM 탐색 및 메모리 기반 솔버 선택 (`memory.py`)

### 2.1 RAM 탐색 계층 (`probe_system_memory`, `get_available_ram_bytes`)
1. `psutil.virtual_memory()` 우선 사용
2. 미설치/실패 시 Windows `GlobalMemoryStatusEx` (`_WindowsMemoryStatusEx` ctypes) 또는 Linux `/proc/meminfo` 폴백
3. 환경 변수 `WHT_SOLVER_AVAIL_RAM_GB` 설정 시 가용 메모리를 명시적으로 오버라이드(테스트 및 컨테이너 메모리 제한 대응)

### 2.2 피크 메모리 추정 공식
- **Pardiso 직접 해법 (`est_pardiso_bytes(n_dof, nnz, degree=2)`)**:
  $$\text{fill\_factor} = \max\left(8.0,\; (p + 1)^2 \cdot \max(n_{\text{dof}}, 1)^{0.25}\right)$$
  $$\text{bytes} = \left\lceil 1.35 \cdot \left(12\,\text{nnz}\cdot\text{fill\_factor} + 12\,\text{nnz} + 64\,n_{\text{dof}}\right) \right\rceil$$
- **AMGCL 반복 해법 (`est_amgcl_bytes(n_dof, nnz)`)**:
  $$\text{csr\_bytes} = 12\,\text{nnz} + 8\max(n_{\text{dof}}, 1), \qquad \text{bytes} = \left\lceil 3.2\,\text{csr\_bytes} + 256\,n_{\text{dof}} \right\rceil$$

### 2.3 `select_linear_solver` 자동 전환 정책
- 기본 예산 비율 `memory_budget_ratio = 0.70`, 기본 자유도 임계치 `dof_threshold = 150_000`.
- `solver_spec == "auto"` 일 때 `est_pardiso_bytes > memory_budget_ratio * avail_bytes` 또는 `n_dof > dof_threshold` 이면 직접 해법 대신 **`"amgcl"`** 로 자동 전환하고 `[HWT-IGA SOLVER][AUTO-SOLVER]` 사유 메시지를 `SolverSelectionDecision.switched_reason` 에 기록한다.
- 예산 이내인 경우 `pypardiso` 가용 여부에 따라 `"pardiso"` 또는 `"splu"` 를 선택한다(`pipeline.py` 에서는 비대칭 팔로워 압력 하중 `move_mesh=True` 유무에 따라 `"pardiso_lu"` 또는 `"pardiso_ldlt"` 로 분기).

---

## 3. 6단계 상용 스타일 솔버 로그 (`solver_log.py`)

### 3.1 로그 제어 및 등급
- `resolve_solver_log_enabled`: `enable_solver_log=False` 또는 환경변수 `WHT_SOLVER_LOG=0|false|no|off` 시 파일 생성 및 표준출력을 비활성화.
- 로그 등급(`SolverLogLevel`): `"summary"`/`"minimal"` (1), `"detailed"`/`"standard"` (2), `"debug"`/`"verbose"` (3).
- UTF-8 (`encoding="utf-8"`) 파일 기록 및 선택적 `stdout` 미러링 지원.

### 3.2 6단계 구조
1. **`[1/6] JOB & SYSTEM ENVIRONMENT`**: 작업명, 호스트/OS, Python/Kratos 버전, 논리 코어·MKL 스레드, 호스트 총/가용 RAM.
2. **`[2/6] MODEL & IMMERSED GRID SUMMARY`**: 격자 분할수, 차수 `p`, 연속성, `IN`/`CUT`/`OUT` 셀 개수, `out_mode`, 제어점·총/활성 DOF, 구적점 및 추정 `nnz`.
3. **`[3/6] CONSTRAINT & INTERACTION SUMMARY`**: 디리클레/고정면, 비활성 고정 DOF 수, 페널티 $\beta$/$C_{\text{pen}}$, 확장 tie 통계.
4. **`[4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION`**: 요청→선택 솔버, Pardiso vs AMGCL 추정 GB vs 가용 예산, 자동 전환 사유, 솔버 파라미터.
5. **`[5/6] STEP & ITERATION PROGRESS TABLE`**: `Step | Inc | Att | Time | dt | LinIt | RelRes | EnergyNorm | Wall(s)` 표 형식 진행 이력.
6. **`[6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY`**: 최대 변위 `max|U|`, 전체 하중, 반력 평형 잔차 $\|\mathbf{R}_{\text{sum}} + \mathbf{F}_{\text{ext}}\|$, `setup/solve/post/total` 소요 시간, 경고 및 최종 수렴 판정 배너.

---

## 4. 입력 검증 게이트 및 Kratos/CSR 솔버 (`validate.py`, `solvers.py`)

### 4.1 `validate.py` 차단 규칙
- `JobSpec` 의 `linear_solver`, `solver_log_level`, `memory_budget_ratio` 를 JSON 직렬화/역직렬화로 보존한다.
- `validate(spec, geo_issues=None)` 는 오류 등급 이슈 발생 시 `ValidationError` 로 실행을 차단한다:
  - 물성/격자/BC 오류: `MAT_E`, `MAT_NU`, `GRID_DEGREE`, `GRID_CELLS`, `BOX_SHAPE`, `BETA_VALUE`, `NO_DIRICHLET`, `GHOST_ALPHA`, `OUT_MODE`, `LINEAR_SOLVER`, `SOLVER_LOG_LEVEL`, `MEMORY_BUDGET_RATIO`, `OUT_CELL_RANGE`.
  - 기하 이슈 연동: `WATERTIGHT`, `VOLUME_SIGN`, 하중면 위 `QUAD_FALLBACK` → `QUAD_FALLBACK_LOADED` 승격 차단.

### 4.2 `kratos_backend/solvers.py`
- `DEFAULT_AMGCL_PARAMS`: `solver_type="amgcl"`, `krylov_type="cg"`, `smoother_type="ilu0"`, `coarsening_type="aggregation"`, `max_iteration=5000`, `tolerance=1e-10`.
- `build_kratos_linear_solver`: `pardiso`, `pardiso_ldlt`, `pardiso_lu`, `pardiso_llt`, `amgcl`, `amgcl_cg`, `cg`, `splu`, `skyline_lu` 생성 및 MKL/스레드 수 설정.
- `solve_linear_csr(K, f, solver_spec="auto", ...)`: SciPy CSR 행렬에 대해 메모리 기반 선택 후 `PardisoFactor` (실패 시 Kratos MKL Pardiso `CompressedMatrix` 폴백), Kratos AMGCL CSR(`_solve_kratos_amgcl_csr`), SciPy preconditioned CG(`_solve_scipy_cg`), `splu` 로 풀고 실제 상대 잔차 `||K u - f||_2 / max(||f||_2, 1e-30)` 를 검증해 `LinearSolveResult` 반환.
- **S04/대칭 제약 주의(`pipeline.py`)**: 팔로워 압력(`move_mesh=True`)으로 접선 강성이 비대칭이 되면 대칭 전용 솔버(`pardiso`, `pardiso_ldlt`, `pardiso_llt`, `amgcl`, `amgcl_cg`, `cg`) 사용을 명시적으로 거부한다.

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_phase3_solver_memory_log.py` | 메모리 추정식, 예산 초과 시 `auto`→`amgcl` 전환, `WHT_SOLVER_AVAIL_RAM_GB`/`WHT_SOLVER_LOG`, 6단계 로그 형식 | **5 passed** (1.56s) | Kratos 무의존 단위 테스트 |
| `WHT_IGA/tests/unit/test_solve_validate.py` | `JobSpec` 직렬화 및 `validate()` 에러/경고 차단 게이트 | **11 passed** (0.55s) | Kratos 무의존 단위 테스트 |
| `WHT_IGA/tests/kratos_spikes/test_p3_solver_amgcl_log.py` | Kratos `AMGCL` CSR/파이프라인 해법, `out_mode="drop"`, 6단계 로그 파일 생성, `WHT_SOLVER_LOG=0` 비활성화 (무가중치 프로브 포함) | **5 passed** (전체 12건 스파이크 13.41s 내 포함) | Kratos 10.4.3 스파이크 검증 |
