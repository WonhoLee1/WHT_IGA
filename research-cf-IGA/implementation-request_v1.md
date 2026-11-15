# WHT_IGA 구현 요청서 v1 (구현 담당 에이전트용)

- 작성일: 2026-11-14
- 기준 계획서: `research-coreform-IGA/master-plan_v2.md` (이하 "v2"). 이 요청서와 v2가 충돌하면 **작업을 멈추고 사용자에게 묻는다.**
- 브랜치: `WHTIGA`
- 대상: 이 저장소에서 WHT_IGA를 구현할 모든 AI 에이전트 (Claude Code, Codex, OpenCode, Cursor 등)
- 요청 범위: **Step E0–E1(환경) + Phase 0(spike와 골격)**. Phase 1 이후는 G0 게이트를 사용자가 승인한 다음에 별도 요청으로 연다 (§9에 미리보기만 있음).

---

## 진행 현황 (2026-11-14 갱신) — 인수인계 받는 에이전트는 여기부터

| 작업 | 상태 | 비고 |
|---|---|---|
| T-E0.1 백업 | 완료 | `C:\Users\GOODMAN\env-backup-2026-11-14\` (해시 일치) |
| T-E0.2 micromamba | 완료 | 2.10.0, SHA256 검증, `MAMBA_ROOT_PREFIX=D:\mamba` |
| T-E0.3 `.mambarc` | 완료 | 기존 env도 목록에 보임 |
| T-E0.4 셸 프로파일 | 완료 | PS7/PS5/cmd/Git Bash(`~/.bash_profile`) 모두 `micromamba activate vdmc` 확인. OMC 줄 보존. `where conda` 유지 |
| T-E0.5 `whtiga` env | 완료 | Python 3.12.16 |
| T-E0.6 패키지(uv) | 완료 | `WHT_IGA/uv.lock` 생성, Kratos wheel 없음 확인 |
| T-E0.7 Kratos 연결 | 완료 (**변경**) | `.pth`는 쓰지 않음. `wht_iga/kratos_backend/_bootstrap.py`가 경로·DLL 설정 (Kratos가 필요 없는 코드를 깨끗하게 유지하려고). 문서: `WHT_IGA/docs/ENVIRONMENT.md` §4 |
| T-E0.8 전역 지침·메모리 | 완료 | |
| T-E0b 빌드 재현 문서 | 완료 | `WHT_IGA/docs/KRATOS_BUILD.md`, `docs/kratos_build/` (CMakeCache, local_changes, local_status, `build_kratos_full.bat.draft`) |
| T-E1 공유 문서 | 완료 | `WHT_IGA/AGENTS.md`, `CLAUDE.md`, `docs/ENVIRONMENT.md`, 루트 `AGENTS.md` 한 줄 |
| T-P0.1 골격 | 완료 | pyproject, README, conftest, test_isolation, `docs/{adr,spikes,verification}/`, `tests/kratos_spikes/` |
| T-P0.2 bootstrap | 완료 | `load_kratos()` |
| T-P0.3 check_env | 완료 | **13/13 통과** |
| T-P0.4 grid 기초 | 완료 | `wht_iga/grid/{knots,bspline,greville,bezier_extraction}.py`, 단위시험 94개 통과 |
| T-P0.5 참조 조립기 | 완료 | `wht_iga/reference/{box,assemble_elastic,solve}.py`. 패치시험 p1–3 1e-10, 대칭·강체모드, 외팔보 Timoshenko 1% (9개 통과) |
| **T-P0.6 spike S01–S18** | **미착수** | 다음 작업. S01은 `check_env.py`로 사실상 확인됨 — pytest로 옮기기만 하면 됨 |
| **T-P0.7 ADR-001** | **미착수** | |
| **T-P0.8 G0 보고서** | **미착수** | |

현재 시험 상태: `micromamba run -n whtiga python -m pytest WHT_IGA/tests -m "not kratos"` → 103개 통과.

**기준 문서 변경 (2026-11-14)**: 계획서 기준은 이제 `research-coreform-IGA/master-plan_v2.1.md`(v2 + add_1 + add_2 + 트랙 A 결과 통합)이다. 이 요청서의 "v2 §…" 참조는 v2.1의 같은 절을 뜻한다(절 번호 유지, 신규 §3.4 워크플로, §5.9 멀티파트·어셈블리, §5.10 메쉬 입출력, spike S19–S23 추가). 트랙 B 범위는 바뀌지 않았다.

재사용할 것 (spike에서 그대로 쓰기):
- `wht_iga.reference.box.BoxSpline`: 직육면체 B-spline (knot, Greville 제어점, `axis_tables`, `evaluate`, `boundary_control_points`). 제어점 번호 `a + n_u*(b + n_v*c)`, 자유도 `3*I + d`.
- `wht_iga.reference.assemble_elastic.assemble_stiffness / assemble_face_traction`: S04의 K 비교 기준.
- `wht_iga.grid.knots.to_kratos`: Kratos knot 규약 변환 (K4).

### 병렬 진행 구성 (에이전트 2개)

T-P0.6 spike를 두 트랙으로 나눈다. **파일 소유권이 겹치지 않게** 나눴으므로 두 에이전트는 서로 기다리지 않고 동시에 진행한다.

| | 트랙 A — 체적·요소·구속 | 트랙 B — 경계·하중·물리 |
|---|---|---|
| spike | S02, S03, S04, S05, S06, S13, S15, S16 | S01, S07, S08, S09, S10, S11, S12, S14, S17, S18 |
| 핵심 질문 | 셀 단위 적분점 기하 + SM 요소 + 구속식(MPC) 경로가 성립하나 | 경계 적분점 기하에 penalty·하중·스프링·열·접촉 조건이 붙나 |
| 쓰기 허용 파일 | `tests/kratos_spikes/test_s02_*`…`test_s16_*` (자기 spike만), `docs/spikes/S02.md`… (자기 spike만), 새 모듈 `wht_iga/kratos_backend/elements.py`, `wht_iga/kratos_backend/constraints.py` | `tests/kratos_spikes/test_s01_*`, `test_s07_*`…`test_s18_*` (자기 spike만), `docs/spikes/S01.md`… (자기 spike만), 새 모듈 `wht_iga/kratos_backend/surface.py`, `wht_iga/kratos_backend/thermal.py` |
| 읽기만 | 그 외 전부 | 그 외 전부 |

공통 규칙:
- **공유 모듈은 둘 다 수정 금지**: `wht_iga/kratos_backend/geometry.py`, `_bootstrap.py`, `wht_iga/grid/*`, `wht_iga/reference/*`, `tests/conftest.py`, `pyproject.toml`, `uv.lock`. 변경이 필요하면 작업을 멈추지 말고 자기 트랙 파일 안에 임시 구현을 두고, 보고서에 "공유 모듈 변경 제안"으로 적는다.
- 공유 헬퍼 (이미 Kratos에서 동작 확인됨, 2026-11-14):
  - `make_nurbs_volume(KM, model_part, box)`: Greville 제어점 노드 + `NurbsVolumeGeometry`. 노드 id = 1 + 제어점 번호.
  - `cell_gauss_points(box, i, j, k)`: 셀 Gauss 점 `[[x,y,z,w],…]`.
  - 확인 결과: 셀 1개 → 기하 1개, 적분점 27, 노드 27, 분할단위성 오차 2e-16, detJ = 1, x(ξ) = ξ 오차 2e-16 (p=2). → S02·S03 일부의 선행 증거. 트랙 A가 정식 spike로 완성한다.
- 시험 실행은 **자기 spike 파일만** 지정해서 돌린다 (`pytest WHT_IGA/tests/kratos_spikes/test_s04_*.py`). 전체 시험은 통합 단계에서 한다.
- 새 패키지 설치가 필요하면 설치하지 말고 보고한다 (`uv.lock` 충돌 방지).
- 트랙 B의 S07–S09는 트랙 A의 결과를 기다리지 않는다. 볼륨은 공유 헬퍼로 만든다.

통합 (T-P0.7 ADR-001, T-P0.8 G0 보고서): 두 트랙 보고가 모두 끝난 뒤 **통합 담당**이 한다. 통합 담당은 이 문서를 만든 메인 세션(Claude Code) 또는 사용자가 지정한 에이전트. 트랙 에이전트는 ADR과 G0 보고서를 쓰지 않는다.

트랙 완료 조건:
1. 담당 spike 각각에 시험 파일 + `docs/spikes/SNN.md` (결과 예/아니오/부분 + 근거 수치)
2. 담당 spike 시험 실행 결과 첨부 (실패도 그대로)
3. 트랙 요약 보고 `docs/spikes/TRACK_A.md` 또는 `TRACK_B.md`: 표 형식 결과, 경로 결정(R1/R2/R3)에 미치는 영향, 공유 모듈 변경 제안

---

## 0. 먼저 읽을 것 (순서대로)

1. 이 문서 전체
2. v2 §0.4 (사용자 결정), §2 (원칙), §5 (Kratos 연동), §9 Phase 0, §11 (위험), §12 (환경)
3. 저장소 루트 `AGENTS.md` (aieng 작업공간 규칙. WHT_IGA는 aieng MCP 도구를 쓰지 않는다)
4. 작업 시작 전 확인: `git status --short --branch`가 `WHTIGA` 브랜치인지 본다. 아니면 멈추고 묻는다.

---

## 1. 목표 한 줄

Coreform IGA의 기능을 Kratos(소스빌드 10.4.3) + Python + Numba로 재현하는 `WHT_IGA/` 패키지를 만든다. 이번 요청에서는 **환경을 세우고, Kratos Python 바인딩으로 핵심 경로가 성립하는지 spike로 판정**한다.

---

## 2. 절대 규칙 (위반 시 작업 무효)

### 2.1 해야 하는 것

| # | 규칙 |
|---|---|
| M1 | 모든 Python 실행은 `micromamba run -n whtiga python ...` 또는 `D:\mamba\envs\whtiga\python.exe`로 한다 |
| M2 | Python 패키지 설치는 **uv pip만** 쓴다: `uv pip install --python D:\mamba\envs\whtiga\python.exe ...` |
| M3 | Kratos는 `E:\kratos-build\install`만 쓴다. import 후 `KratosMultiphysics.__file__`이 이 경로 아래인지 반드시 확인한다 |
| M4 | 수치 결과를 주장할 때는 실제 실행 출력(숫자)을 보고서에 붙인다 |
| M5 | 시험은 **기준값과 비교**한다(해석해, 참조 조립기, Coreform CSV). 실행 성공 여부만 보는 시험은 금지 |
| M6 | 각 spike 결과를 `WHT_IGA/docs/spikes/SNN.md`에 기록한다 (형식 §8) |
| M7 | 파일을 덮어쓰거나 지우기 전에 내용을 확인하고, 셸 프로파일·설정 파일은 먼저 백업한다 |
| M8 | 기본값으로 대체한 입력은 결과에 `assumed`로 표시한다 |
| M9 | 막히면 추측으로 우회하지 말고 §10 "멈추고 물을 조건"에 따라 사용자에게 묻는다 |

### 2.2 하면 안 되는 것

| # | 금지 |
|---|---|
| X1 | `conda install`, `conda create`, `pip install`(uv 없이) 사용 |
| X2 | `whtiga` env에 pip Kratos wheel(`KratosMultiphysics*`) 설치 |
| X3 | `wht_topopt` 또는 `aieng` 패키지 import (복사·적응은 허용, 출처 주석 필수) |
| X4 | `D:\PythonCodeStudy\kratos-src` 수정, Kratos 재빌드, `E:\kratos-build` 안의 파일 수정 |
| X5 | 기존 conda env(`vdmc`, `aieng311`, `calculix-env` 등) 수정·삭제, miniconda 제거, PATH에서 miniconda 항목 제거 |
| X6 | 사용자 요청 없이 git commit / push |
| X7 | `test.skip`, `xfail`로 실패 숨기기, TODO 스텁을 완료로 보고 |
| X8 | `research-coreform-IGA/master-plan_v1.md`, `master-plan_v2.md` 수정 (변경 제안은 보고서로) |
| X9 | 루트 `AGENTS.md`의 §4 표 한 줄 추가 외 수정 (§9의 `conda activate` 줄은 사용자 승인 후에만) |
| X10 | Rust 도구체인 설치, C++ 커스텀 앱 착수 (Phase 0 범위 밖) |

---

## 3. 확정된 사실 (다시 조사하지 말 것)

| 항목 | 값 |
|---|---|
| Kratos 설치 | `E:\kratos-build\install` (10.4.3, 커밋 96bad120, Release, MSVC 1941) |
| Kratos 빌드에 쓴 Python | `C:\Users\GOODMAN\miniconda3\envs\vdmc\python.exe` (3.12) |
| Kratos 런타임 DLL | `E:\kratos-build\install\libs`, `D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin` |
| Kratos 소스 | `D:\PythonCodeStudy\kratos-src` (읽기 전용 참고) |
| Kratos 빌드 트리 | `E:\kratos-build\Release\CMakeCache.txt` |
| 빌드 스크립트 | `D:\PythonCodeStudy\kratos-deps\build_kratos_full.bat` — **1바이트, 내용 소실** |
| 이미 import 성공한 조합 | vdmc python + `PYTHONPATH=E:\kratos-build\install;E:\kratos-build\install\libs` + PATH 앞에 condadeps `Library\bin` |
| 셸 프로파일 | `C:\Users\GOODMAN\OneDrive\문서\PowerShell\profile.ps1`, `C:\Users\GOODMAN\OneDrive\문서\WindowsPowerShell\profile.ps1` (각 1–6행 conda init 블록, 이후 OMC 환경변수 줄은 **보존**) |
| cmd AutoRun | `HKCU\Software\Microsoft\Command Processor` → `conda_hook.bat` |
| `~/.condarc` | conda-forge + defaults, flexible |
| Git Bash | `~/.bashrc` 없음 |
| aieng CalculiX | `conda run -n calculix-env ccx`로 자동 탐색 → miniconda와 PATH 유지 필요 |
| Rust | 미설치 |
| MSVC / CMake | 14.41 / 4.4 (`D:\PythonCodeStudy\kratos-deps\buildenv\Scripts`) — Phase 0에서는 쓰지 않음 |

---

## 4. Kratos 사용 시 알려진 함정 (반드시 숙지)

| # | 함정 | 올바른 처리 |
|---|---|---|
| K1 | Python에서 `QuadraturePointGeometry`를 직접 만들 수 없다 | `NurbsVolumeGeometry.CreateQuadraturePointGeometries(결과, 2, [[x,y,z,w],…])`를 쓴다 |
| K2 | 위 호출은 기본 GAUSS 정보를 써서 **Option 1**로 간다: 주어진 점 **전부를 담은 기하 1개**가 나온다. `nderiv`는 **반드시 2** | 셀 하나(또는 셀×재료 하나)당 1회 호출. 다른 셀의 점을 섞지 않는다 |
| K3 | 제어점 집합은 점들의 **중심**이 속한 knot span으로 정해진다 | 셀 면 위의 점은 ε(예: 1e-12·h)만큼 셀 안쪽으로 이동 |
| K4 | Kratos knot 규약은 끝 knot을 **p번** 반복한다 (길이 = 제어점 수 + p − 1). 단일 span p=2면 `[0,0,1,1]` | `grid/knots.py`에서 표준 clamped(p+1번 반복) ↔ Kratos 규약 변환. 왕복 단위시험 필수 |
| K5 | 요소 이름의 "3D4N"은 prototype일 뿐, 실제 기하는 생성 시 교체된다 | `model_part.CreateNewElement("SmallDisplacementElement3D4N", id, qp_geom, props)` 형태. `test_nurbs_volume_element.py` 참고 |
| K6 | IgaApp `SupportPenaltyCondition`은 **3성분을 동시에** penalty한다 | 일부 성분 구속(대칭면)은 v2 §5.4 BC-B 방식 |
| K7 | IgaApp `SolidElement`, `SupportSolidCondition`은 **2D 전용** | 3D 솔리드는 StructuralMechanics 요소 사용 |
| K8 | Kratos core에 HHT scheme이 없다 | **HHT를 구현하지 않는다** (사용자 결정 D5). 동역학은 Kratos Bossak(기본)/Newmark/BDF scheme 사용 |
| K9 | `AssignIntegrationPointsToBackgroundElementsProcess`는 결과 매핑용(output 단계)이다 | 해석 경로에 쓰지 않는다 |
| K10 | pardiso 솔버는 `LinearSolversApplication.HasMKL()`이 True일 때만 등록된다 | S01에서 확인하고 분기 |
| K11 | master 자유도는 다른 구속의 slave가 될 수 없다 (연쇄 구속 금지) | 모든 합성은 Python에서 미리 계산 (v2 §5.2.3) |
| K12 | 요소가 하나도 붙지 않은 노드의 자유도는 DofSet에 안 들어갈 수 있다 | S05에서 확인. 필요하면 강성 0 `NodalConcentratedElement3D1N` 부착 |
| K13 | 노드 변수는 노드 생성 **전에** `AddNodalSolutionStepVariable`로 등록해야 한다 | 빌더에서 순서 고정 |
| K14 | 자유도는 `VariableUtils().AddDof(...)` 또는 `node.AddDof(...)`로 추가해야 `GetDof`가 동작한다 | 빌더에서 노드 생성 직후 |

---

## 5. 작업 목록 — Step E0–E1 (환경)

각 작업은 **수락 기준을 전부 만족**해야 완료다. 완료 보고에는 실행한 명령과 출력 요약을 붙인다.

### T-E0.1 백업

- 할 일: 아래를 `C:\Users\GOODMAN\env-backup-2026-11-14\`에 복사
  - 두 `profile.ps1`
  - cmd AutoRun 값: `reg export "HKCU\Software\Microsoft\Command Processor" cmd_autorun.reg`
  - `~/.condarc`
- 수락 기준: 백업 폴더에 4개 파일 존재, 원본과 해시 일치

### T-E0.2 micromamba 설치

- 할 일: micromamba 단일 실행 파일을 `%LOCALAPPDATA%\micromamba\micromamba.exe`에 설치. 사용자 환경변수 `MAMBA_ROOT_PREFIX=D:\mamba`.
- 출처: 공식 배포(`https://micro.mamba.pm/api/micromamba/win-64/latest` 또는 GitHub `mamba-org/micromamba-releases`). 다운로드 파일은 빈 새 폴더에 받는다.
- 수락 기준: `micromamba --version` 출력, `D:\mamba` 생성

### T-E0.3 `~/.mambarc`

```yaml
channels:
  - conda-forge
channel_priority: strict
envs_dirs:
  - D:\mamba\envs
  - C:\Users\GOODMAN\miniconda3\envs
```

- 수락 기준: `micromamba env list`에 기존 env(vdmc 등)도 보인다

### T-E0.4 셸 프로파일 교체

- PowerShell 7/5 프로파일: 1–6행 conda init 블록만 `micromamba shell init -s powershell -r D:\mamba`가 생성하는 블록으로 교체. **OMC 환경변수 줄은 그대로 둔다.**
- cmd AutoRun: `micromamba shell init -s cmd.exe -r D:\mamba`
- Git Bash: `micromamba shell init -s bash -r D:\mamba` (없던 `~/.bashrc` 생성)
- base 자동 활성화 없음
- 수락 기준:
  - 새 PowerShell, cmd, Git Bash 각각에서 `micromamba activate vdmc` 성공
  - 프로파일 diff에서 OMC 줄 보존 확인
  - `where conda`가 여전히 miniconda를 가리킴 (aieng CalculiX용)

### T-E0.5 `whtiga` env 생성

```
micromamba create -n whtiga -c conda-forge python=3.12 mkl hdf5 -y
```

- Python은 3.12.x. 버전을 기록한다.
- uv가 없으면 `micromamba install -n whtiga -c conda-forge uv`로 넣는다 (이것만 예외적으로 micromamba 설치 허용).
- 수락 기준: `micromamba run -n whtiga python -c "import sys; print(sys.version)"` → 3.12.x

### T-E0.6 패키지 설치 (uv)

- `WHT_IGA/pyproject.toml`을 먼저 만든다 (T-P0.1 참고). 그다음:
  ```
  uv pip install --python D:\mamba\envs\whtiga\python.exe -e "WHT_IGA[cad,viz,solve,test]"
  ```
- extras 정의:
  - 필수: `numpy`, `scipy`, `numba`
  - `cad`: `cadquery-ocp`
  - `viz`: `pyvista`, `meshio`, `h5py`
  - `solve`: `pypardiso`
  - `test`: `pytest`, `hypothesis`
- `uv.lock` 생성 (`uv lock` 또는 `uv pip compile`로 고정 파일). 방식은 보고서에 적는다.
- 수락 기준: `import numpy, scipy, numba, OCP, pyvista, meshio, h5py, pypardiso` 성공. **Kratos wheel이 설치되지 않았음** 확인(`uv pip list`에 `KratosMultiphysics` 없음)

### T-E0.7 Kratos 연결

- `D:\mamba\envs\whtiga\Lib\site-packages\kratos_ebuild.pth`:
  ```
  E:\kratos-build\install
  E:\kratos-build\install\libs
  ```
- DLL 경로는 `wht_iga/kratos_backend/_bootstrap.py`에서 `os.add_dll_directory`로 등록 (T-P0.2).
- 수락 기준: §6 S01 통과

### T-E0.8 메모리·전역 지침 갱신 (Claude Code 에이전트인 경우만)

- `~/.claude/CLAUDE.md` 끝에 한 줄: `- Python 환경: conda 사용 금지 — micromamba(root D:\mamba) + uv pip.`
- 다른 에이전트는 이 작업을 건너뛰고 보고서에 "미수행"으로 적는다.

### T-E0b Kratos 빌드 재현 문서

- 할 일:
  1. `E:\kratos-build\Release\CMakeCache.txt`를 `WHT_IGA/docs/kratos_build/CMakeCache.txt`로 복사
  2. 캐시에서 `KRATOS_APPLICATIONS`(또는 앱 목록 관련 키), `CMAKE_BUILD_TYPE`, `PYTHON_EXECUTABLE`, MKL·HDF5·SuiteSparse 관련 키를 추출
  3. `git -C D:\PythonCodeStudy\kratos-src diff --stat > WHT_IGA/docs/kratos_build/local_changes.txt`
  4. `WHT_IGA/docs/KRATOS_BUILD.md` 작성: 재현 절차(vcvars64 → cmake 구성 → ninja `-j 6` → install), 추출한 옵션, 로컬 수정 요약, 주의(`-j 6` 초과 시 cl.exe 충돌 기록)
  5. 복원한 빌드 스크립트 초안은 `WHT_IGA/docs/kratos_build/build_kratos_full.bat.draft`로 둔다. 원래 위치의 1바이트 파일은 **건드리지 않는다**
- 금지: 실제 빌드 실행
- 수락 기준: 문서와 3개 파일 존재, 옵션 목록이 CMakeCache와 일치

### T-E1 다른 에이전트와 공유

| 파일 | 내용 |
|---|---|
| `WHT_IGA/docs/ENVIRONMENT.md` | micromamba 설치·활성화, uv 명령, Kratos 연결(.pth, DLL, wheel 금지), `check_env.py` 사용법, 자주 나는 오류와 해결 |
| `WHT_IGA/AGENTS.md` | §2 절대 규칙 요약(10–20줄) + ENVIRONMENT.md와 이 요청서 링크 |
| `WHT_IGA/CLAUDE.md` | 한 줄: `@AGENTS.md` |
| 루트 `AGENTS.md` | §4 Workspace Layout 표에 한 줄만 추가: `` | `WHT_IGA/` | **Active** | Kratos 기반 Immersed IGA 솔버 (자체 env: micromamba `whtiga`, 규칙은 `WHT_IGA/AGENTS.md`) | `` |

- 수락 기준: 루트 `AGENTS.md`의 diff가 정확히 한 줄 추가

---

## 6. 작업 목록 — Phase 0

### T-P0.1 패키지 골격

```
WHT_IGA/
├── pyproject.toml        # name="wht_iga", version="0.0.1", requires-python=">=3.12,<3.13"
├── env/whtiga.yml        # T-E0.5와 같은 내용
├── README.md
├── docs/{adr,spikes,verification}/
├── scripts/check_env.py
├── wht_iga/__init__.py
├── wht_iga/kratos_backend/__init__.py, _bootstrap.py
├── wht_iga/grid/__init__.py, knots.py, bspline.py, greville.py, bezier_extraction.py
├── wht_iga/reference/__init__.py, assemble_elastic.py, solve.py
└── tests/{unit,kratos_spikes}/ , tests/test_isolation.py, tests/conftest.py
```

- `conftest.py`: marker `kratos`(Kratos 필요), `slow` 등록. Kratos import 실패 시 `kratos` 시험은 **실패**로 처리한다(skip 금지). 단, `-m "not kratos"`로 실행하면 제외.
- `test_isolation.py`: (1) `wht_iga` 소스에 `import wht_topopt`, `from wht_topopt`, `import aieng`, `from aieng` 문자열 없음 (2) `import wht_iga` 후 `sys.modules`에 해당 패키지 없음.
- 수락 기준: `micromamba run -n whtiga python -m pytest WHT_IGA/tests -m "not kratos"` 통과

### T-P0.2 Kratos bootstrap

- `_bootstrap.py`:
  - 환경변수 `WHT_IGA_KRATOS_PATH`(기본 `E:\kratos-build\install`), `WHT_IGA_KRATOS_DLL_PATH`(기본 `E:\kratos-build\install\libs;D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin`)
  - `sys.path` 앞에 삽입(이미 .pth로 있으면 중복 삽입 안 함), 각 DLL 경로 `os.add_dll_directory`
  - `import KratosMultiphysics` 후 `__file__` 경로 검사. 다르면 `RuntimeError("wrong Kratos: ...")`
  - 함수 `load_kratos(apps=("IgaApplication","StructuralMechanicsApplication", ...)) -> dict[name, module]`
- Kratos import는 이 모듈 안에서만. `wht_iga/__init__.py`는 Kratos를 import하지 않는다.
- 수락 기준: S01

### T-P0.3 `check_env.py`

- 출력 항목: Python 버전, numpy/scipy/numba/OCP 버전, Kratos 경로·버전, 앱별 import 결과(IgaApplication, StructuralMechanics, LinearSolvers, ConstitutiveLaws, ConvectionDiffusion, ContactStructuralMechanics), `HasMKL()`, numba `@njit` 함수 1개 컴파일·실행 결과.
- 하나라도 실패하면 종료 코드 1.

### T-P0.4 grid 기초 (Kratos 불필요)

- `knots.py`: `uniform_clamped(lo, hi, n_cells, p, continuity) -> ndarray`(표준 규약), `to_kratos(knots, p)`, `from_kratos(knots, p)`
- `bspline.py`: `@njit` Cox–de Boor 값·1차·2차 미분 (`basis_funs(span, u, p, knots)`, `find_span`)
- `greville.py`: Greville 점
- `bezier_extraction.py`: Borden 2011 1D 추출 + 텐서곱
- 단위시험(수락 기준):
  - 분할단위성 1e-14, 미분의 유한차분 일치 1e-7
  - Kratos 규약 왕복 변환 일치
  - 항등 사상: Greville 제어점으로 x(ξ) = ξ (p=1,2,3, 1e-14)
  - Cᵉ·B = N (1e-14)

### T-P0.5 참조 조립기 최소판

- 범위: untrimmed 직육면체, 균일 텐서 B-spline, 선형 탄성(등방), 변위 Dirichlet(강한 부과: 해당 면 제어점 고정, 또는 penalty), 체적력.
- `assemble_elastic.py`: 셀 루프 @njit, (p+1)³ Gauss, scipy.sparse COO → CSR
- `solve.py`: pypardiso(있으면) / scipy spsolve
- 수락 기준:
  - 지정 변위 패치시험 u = (x, 2y, 3z) (E=1, ν=0)를 p=1,2,3에서 기계 정밀도(1e-10)로 재현
  - 외팔보 끝단 처짐이 Euler–Bernoulli 대비 격자 수렴 경향을 보임 (p=3 세밀 격자 1% 이내)

### T-P0.6 Spike S01–S18

- 파일: `WHT_IGA/tests/kratos_spikes/test_sNN_<이름>.py` (marker `kratos`)
- 각 spike는 v2 §5.8 표의 질문·통과 기준·실패 시 대응을 그대로 따른다. 아래는 핵심 spike의 구체 절차다.

**S01 환경**
- `check_env.py`와 같은 검사를 pytest로. 통과 기준: 전부 import, `HasMKL()` 값 기록.

**S02 항등 사상**
- 절차: 2×2×2 셀, p=2. 표준 knot(물리 좌표) → Kratos 규약 변환 → Greville 점에 노드 생성 → `KM.NurbsVolumeGeometry(nodes, 2, 2, 2, kU, kV, kW)`.
- 셀 내부 임의 점 10개에서 `geom.GlobalCoordinates(...)` 결과가 입력 좌표와 1e-14 이내.
- `CreateQuadraturePointGeometries(qp, 2, [[x,y,z,w]])` 후 `qp[0].DeterminantOfJacobian(0)` ≈ 1 (1e-14).
- Frame 회전(30°, 임의 축) 변형도 같은 검사.

**S03 셀 단위 적분점 기하**
- 절차: 한 셀 안의 점 n=27(Gauss)과 n=7(임의, 가중치 임의)을 넣는다.
- 확인: 결과 기하 개수 1, `IntegrationPointsNumber()==n`, 가중치 합 보존, 각 점에서 `ShapeFunctionsValues()` 행 합 = 1, 기하의 노드 집합 = 해당 셀 지지 제어점 (p+1)³개.
- 점 하나를 셀 면 위(경계)에 두었을 때 노드 집합이 바뀌는지 확인하고, K3 처리의 필요성을 기록.

**S04 SM 요소 부착**
- 절차: T-P0.5와 같은 큐브. 셀마다 S03 방식 기하 + `SmallDisplacementElement3D4N`, `LinearElastic3DLaw`.
- 확인 1: Kratos 전역 K(`GetSystemMatrix` 또는 B&S Build 후 scipy 변환)와 참조 조립기 K가 같은 자유도 번호 대응에서 1e-12 상대오차.
- 확인 2: 패치시험 변위 1e-10.
- 확인 3: `CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, …)`, von Mises가 값 반환.
- `TotalLagrangianElement3D4N`, `UpdatedLagrangianElement3D4N`도 초기화 가능 여부만 확인.

**S05 구속식 성능·정확도**
- 절차 A(정확도): 작은 문제에서 임의 T(slave 20개, 각 master 27개)를 `LinearMasterSlaveConstraint`로 만들고 풀어서, scipy로 직접 푼 Tᵀ K T 결과와 1e-10 비교. slave 값 복원 확인.
- 절차 B(성능): 3 slave × 81 master 구속 1e5개 생성 시간 측정(목표 < 30 s), 조립·풀이 시간 기록.
- 절차 C: 요소가 없는 master 전용 노드의 자유도가 DofSet에 포함되는지.

**S06–S18**
- v2 §5.8 표대로. 각 결과를 "예/아니오 + 근거 수치"로 기록.

### T-P0.7 ADR-001 경로 결정

- `WHT_IGA/docs/adr/ADR-001-route.md`: S03/S04/S05 결과로 v2 §5.3 결정 트리 적용 → R1/R2/R3 중 하나 **제안**. 확정은 사용자가 한다.

### T-P0.8 G0 게이트 보고서

- `WHT_IGA/docs/G0_report.md`:
  - S01–S18 결과 요약 표
  - untrimmed 큐브 p = 1, 2, 3에서 Kratos 경로와 참조 조립기 변위 최대 차이 (목표 ≤ 1e-10)
  - 실패·보류 항목과 제안
  - v2 수정 제안 (있으면)

---

## 7. 코딩 규칙

- Python 3.12, 타입 힌트, `from __future__ import annotations`.
- Numba 커널: 순수 함수, NumPy 배열 입출력, Python 객체·딕셔너리 금지, `@njit(cache=True)` (병렬은 측정 후 `parallel=True`).
- 단위: 내부는 사용자 입력 단위 그대로(무차원 일관 단위계). 단위 변환 코드 금지.
- 좌표·인덱스 규약: 셀 인덱스 (i, j, k) → 선형 `i + nx*(j + ny*k)`. 이 규약을 `grid/` 모듈 docstring에 명시.
- 부동소수 비교: 시험에서 허용오차를 반드시 명시(`np.testing.assert_allclose(rtol=..., atol=...)`).
- 출처 표기: 알고리즘 구현에 논문 번호를 주석으로(예: `# Borden et al. 2011, Alg. 1`).
- 파일 하나가 400줄을 넘으면 책임별로 나눈다.

---

## 8. 보고 형식

### 8.1 spike 기록 (`docs/spikes/SNN.md`)

```markdown
# SNN <제목>
- 질문: (v2 §5.8 그대로)
- 결과: 예 / 아니오 / 부분
- 근거: (실행 명령, 핵심 출력 숫자, 허용오차 대비)
- 환경: Python x.y.z, Kratos 10.4.3 (경로), 날짜
- 결론과 영향: (경로 선택에 미치는 영향, 실패 시 대응 적용 여부)
- 재현: pytest WHT_IGA/tests/kratos_spikes/test_sNN_*.py
```

### 8.2 작업 완료 보고 (작업마다)

- 작업 ID, 변경·생성 파일 목록, 실행 명령, 수락 기준별 통과/실패와 근거 숫자, 남은 문제.
- 실패를 성공처럼 쓰지 않는다. 미수행 항목은 "미수행"으로 적는다.

---

## 9. 정의: 이번 요청의 완료 (Definition of Done)

1. T-E0.1 – T-E1 모든 수락 기준 통과
2. T-P0.1 – T-P0.8 모든 수락 기준 통과, 또는 실패 항목이 근거와 함께 보고됨
3. `micromamba run -n whtiga python WHT_IGA/scripts/check_env.py` 종료 코드 0
4. `micromamba run -n whtiga python -m pytest WHT_IGA/tests` 결과 첨부 (실패가 있으면 그대로 보고)
5. 변경 파일에 TODO 스텁, skip, xfail 없음 (있으면 사유와 함께 보고)
6. G0 보고서와 ADR-001 제안 제출 → **사용자 승인 대기**. 승인 전 Phase 1 착수 금지.

### 9.1 Phase 1 이후 미리보기 (이번 요청 범위 아님)

| Phase | 핵심 산출물 | 게이트 (v2 §9) |
|---|---|---|
| 1 | S1 형상, S2 격자, S3 분류, S4 절단 적분, penalty BC, 선형 정적, DB | V01 ≤ 1e-5, V02 peak ≤ 1% / 수렴률 ≥ 1.8, V03 5% |
| 2 | small-cell extension, 구속식 배치, 조건수 도구 | V-COND, V06 CSV 허용오차 |
| 3 | 다중재료, tie, coupling, probe, THB, VTKHDF | V03 moment, V04, V05, V02 국소, E01, E02 |
| 4 | 고유치, 암시적 동역학(Kratos Bossak), 스펙트럼, 열, 열변형 | V09, V10, V11, V12 |
| 5 | 비선형, 잠김 방지(C++), 접촉 | V07, V08, E03, E04, E05 |
| 6 | 성능 (C++ 이전, 반복 솔버, 병렬) | 1M DOF < 30분 |

---

## 10. 멈추고 사용자에게 물을 조건

다음 중 하나라도 생기면 **추측으로 진행하지 말고** 상황, 시도한 것, 선택지를 정리해 묻는다.

1. `whtiga` env(uv Python 3.12)에서 Kratos import가 실패하고, vdmc Python으로 바꾸는 것 외에 방법이 없을 때 (런타임 변경은 사용자 결정)
2. 셸 프로파일 교체 후 기존 기능(aieng CalculiX 탐색, `micromamba activate vdmc`)이 깨질 때 → 즉시 백업으로 복원 후 보고
3. S03, S04, S05 중 하나라도 실패해 경로(R1/R2/R3)가 바뀔 때
4. v2 내용이 Kratos 실제 동작과 다를 때 (v2는 수정하지 말고 보고)
5. 라이선스가 불분명한 코드·데이터를 가져와야 할 때
6. 사용자 승인이 필요한 파일(루트 `AGENTS.md` §9, 프로파일 외 시스템 설정)을 바꿔야 할 때
7. 작업량이 예상(Phase 0 약 2주)의 2배를 넘을 것으로 보일 때

---

## 11. 참고 위치

| 내용 | 위치 |
|---|---|
| 계획서 | `research-coreform-IGA/master-plan_v2.md` |
| Coreform 문서 (목표 기능·검증값) | `research-coreform-IGA/user-guide/` |
| 검증 CSV | `research-coreform-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/` |
| 문헌 요약 | `research-coreform-IGA/papers/`, `문헌종합.md`, `papers/implementation_synthesis.md` |
| Kratos 3D IGA 예시 | `D:\PythonCodeStudy\kratos-src\applications\IgaApplication\tests\test_nurbs_volume_element.py` |
| 적분점 기하 Python 바인딩 | `D:\PythonCodeStudy\kratos-src\kratos\python\add_geometries_to_python.cpp` l.166–181, l.313–330, l.504–505 |
| Option 1 동작 | `D:\PythonCodeStudy\kratos-src\kratos\geometries\nurbs_volume_geometry.h` l.820–840 |
| 구속식 Python 사용 예 | `D:\PythonCodeStudy\kratos-src\kratos\tests\test_constraint_restart.py` |
| penalty 조건 | `D:\PythonCodeStudy\kratos-src\applications\IgaApplication\custom_conditions\support_penalty_condition.cpp` |
| 기존 Kratos 실행 패턴(참고, import 금지) | `aieng/src/aieng/topopt/adapters/kratos.py`, `aieng/docs/topopt_solver_adapters.md` |
