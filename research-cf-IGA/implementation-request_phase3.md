# WHT_IGA 구현 요청서 Phase 3 (구현 담당 에이전트용)

- 작성일: 2026-11-14
- 기준 계획서: `research-coreform-IGA/master-plan_v2.1.md` (이하 "v2.1"). 이 요청서와 v2.1이 충돌하면 **작업을 멈추고 사용자에게 묻는다.**
- 브랜치: `WHTIGA`
- 대상: 이 저장소에서 WHT_IGA Phase 3을 구현할 모든 AI 에이전트
- 요청 범위: **Phase 3 (다중재료, 멀티파트·어셈블리, IGA–FE 혼합 Beta, 메쉬 입출력/Gmsh/Body-fit, `wht_iga.selection`, coupling, 다단계/amplitude, THB 세분화, probe, VTKHDF, L2 노름, `HWT-IGA SOLVER` 상용 로깅 및 실시간 RAM 기반 지능형 솔버 전환/`out_mode`/`AMGCL`)**. Phase 1–2 산출물 및 G2 게이트 사용자 승인이 진입 전제이며, Phase 4 이후는 G3 게이트를 사용자가 승인한 다음에 별도 요청으로 연다 (§9에 미리보기만 있음).

---

## 진입 상태 (Phase 3 착수 전 확인) — 미확인 항목은 TBD로 둔다

| # | 전제 | 상태 |
|---|---|---|
| P3-0 | Phase 1–2 산출물 존재 | 전제: `wht_iga/geometry/`, `grid/`, `classify/`, `quadrature/`, `stabilize/` (`support.py`, `extension.py`, `fcm.py`, 조건수), `operators/` (`compose.py`, `numbering.py`), `meshio_ext/` (리더 6종 + `MeshModel`), `kratos_backend/`, `solver/`, `docs/phase1/`, `docs/phase2/` |
| P3-1 | G2 게이트 | **사용자 승인 대기 전제** (V-COND 평탄성, V06 CSV 전량 통과, V01/V02 extension 재통과, V02 메쉬 파일 입력 ≤ 1%). **G2 승인 전 착수 금지** |
| P3-2 | 트랙 C 실측 반영 | 확인 (`WHT_IGA/docs/spikes/TRACK_C.md`, 40/40 통과): S21 `MoveMeshFlag=True`+`AndCriteria`, S22 `ContactSM` `GenerateFaces()` 실패 및 `Triangle3D3`+MPC+`SpringDamperElement3D` tie(`5.220e-11`) 통과, S23 `nearest_element` 운동학적 패치(`4.44e-16`) 통과 vs `AssignNeighboursUtility`(`7.94e-02`) 미달 |
| P3-3 | `gmsh` / `h5py` / `pyNastran` / `meshio` | `whtiga` env 내 설치 여부 착수 전 `check_env.py`로 확인 (Gmsh는 별도 실행 파일/프로세스로만 호출, `import gmsh` 금지) |

---

## 0. 먼저 읽을 것 (순서대로)

1. 이 문서 전체
2. v2.1 §5.9 (멀티파트·어셈블리), §5.10 (메쉬 입출력·Body-fit·Gmsh·Set 키워드 출력)
3. v2.1 §7.2 THB, §7.7 (다중재료와 결합), §7.12 (후처리와 시각화)
4. v2.1 §8.5 (`wht_iga.selection` 통합 선택·계측·접촉 사전 점검·`GeomType` 가이드)
5. v2.1 §9 Phase 3 (산출물·게이트 G3), §10.1–10.2 (V02-local, V03-moment, V04, V05, E01, E02), §4.3/§4.5 (`VTKHDF` + `SQLite` 하이브리드, `probe_data.json`)
6. `WHT_IGA/docs/spikes/TRACK_C.md` (S22, S23 실측)
7. 작업 시작 전 확인: `git status --short --branch`가 `WHTIGA` 브랜치인지 본다. 아니면 멈추고 묻는다.

---

## 1. 목표 한 줄

멀티셀 공유 장과 THB 국소 계층 세분화(2:1 균형·truncation)를 완성하고, Kratos `ModelPart` 계층 위에서 멀티파트·어셈블리(dependent/independent instance, tie 자동 선택, RBE2/RBE3, IGA–FE 혼합 Beta), 5종 메쉬+Set 키워드 왕복 출력, `wht_iga.selection`, 다단계·amplitude·가중치 0 probe·`VTKHDF`·$L^2$ 노름, **`HWT-IGA SOLVER` 상용 6단 로깅(On/Off)·실시간 메모리 기반 `auto` 솔버 전환(`pardiso_ldlt`↔`AMGCLSolver`)·`out_mode`**를 통합하여 게이트 G3(V02-L3, V03-moment, V04, V05, E01, E02)를 전량 통과한다.

---

## 2. 절대 규칙 (위반 시 작업 무효)

### 2.1 해야 하는 것

| # | 규칙 |
|---|---|
| M1 | 모든 Python 실행은 `micromamba run -n whtiga python ...` 또는 `D:\mamba\envs\whtiga\python.exe`로 한다 |
| M2 | Python 패키지 설치는 **uv pip만** 쓴다: `uv pip install --python D:\mamba\envs\whtiga\python.exe ...` |
| M3 | Kratos는 `E:\kratos-build\install`만 쓴다. import 후 `KratosMultiphysics.__file__`이 이 경로 아래인지 반드시 확인한다 |
| M4 | 수치 결과를 주장할 때는 실제 실행 출력(숫자)을 보고서에 붙인다 |
| M5 | 시험은 **기준값과 비교**한다(V02-refinement/V03/V04/V05 CSV, E01 Roark Kt=2.23568, E02 정성 추세, dependent instance 1e-10 일치, 메쉬/Set 왕복 무손실). 실행 성공 여부만 보는 시험은 금지 |
| M6 | 각 작업 결과를 `WHT_IGA/docs/phase3/` 아래에 기록한다 (형식 §8) |
| M7 | 파일을 덮어쓰거나 지우기 전에 내용을 확인하고, 설정 파일은 먼저 백업한다 |
| M8 | 기본값으로 대체한 입력(재료, penalty, search radius 등)은 결과에 `assumed`로 표시한다 |
| M9 | 막히면 추측으로 우회하지 말고 §10 "멈추고 물을 조건"에 따라 사용자에게 묻는다 |

### 2.2 하면 안 되는 것

| # | 금지 |
|---|---|
| X1 | `conda install`, `conda create`, `pip install`(uv 없이) 사용 |
| X2 | `whtiga` env에 pip Kratos wheel(`KratosMultiphysics*`) 설치 |
| X3 | `wht_topopt` 또는 `aieng` 패키지 import (복사·적응은 허용, 출처 주석 필수) |
| X4 | `D:\PythonCodeStudy\kratos-src` 수정, Kratos 재빌드, `E:\kratos-build` 안의 파일 수정 |
| X5 | 기존 conda env 수정·삭제, miniconda 제거, PATH에서 miniconda 항목 제거 |
| X6 | 사용자 요청 없이 git commit / push |
| X7 | `test.skip`, `xfail`로 실패 숨기기, TODO 스텁을 완료로 보고 |
| X8 | `research-coreform-IGA/master-plan_v1.md`, `master-plan_v2.md`, `master-plan_v2.1.md` 수정 (변경 제안은 보고서로) |
| X9 | Phase 1–2 확정 소유 파일(`WHT_IGA/docs/phase1/`, `WHT_IGA/docs/phase2/`, `WHT_IGA/tests/verification/data/v01*`–`v03*`, `v06*`) 임의 수정 — 공유 모듈(`grid/`, `operators/`, `kratos_backend/`, `geometry/`) 확장이 필요하면 기존 API 하위호환을 유지하고, 불가피한 변경은 보고서에 명시한다 |
| X10 | Rust 도구체인 설치, C++ 커스텀 앱 착수 (Phase 5 범위) |

---

## 3. 확정된 사실 (다시 조사하지 말 것)

| 항목 | 값 (출처) |
|---|---|
| 멀티파트 개념 모델 | `Part` 1개 = 형상 + 이산화 1개(`Immersed`/`BodyFit`/`Native`/`Rigid`) + 재료 셀. 한 패밀리 내 dependent/independent instance 혼용 금지 (`validate()` 오류) (v2.1 §5.9.1) |
| Kratos `ModelPart` 계층 | 루트 `Structure` 1개 아래 `Parts/<part>.<instance>`(`Cells.<mat>`, `Layers.<surf>`, `Virtual.<surf>`, `Sets.<name>`), `Interactions/{Tie_n, Contact_n, Coupling_n}`, `Resolved/{BC_n, Load_n}` (v2.1 §5.9.2) |
| 구속 상수·연쇄 금지 | **c = 0 고정**, 동일 DOF 이중 slave 및 slave→master 연쇄 구속 금지 (`compile_plan()` / `job.resolve()`에서 거부) (v2.1 §5.2.3, §5.9.4) |
| S22 실측 (tie/접촉) | `ContactSM` 고수준 프로세스(`mesh_tying_process`, `penalty_contact_process`)는 `InterfacePreprocessCondition`이 `QuadraturePointGeometry::GenerateFaces()` 호출로 `RuntimeError` 실패. 대안 **`Triangle3D3` 표면층 + MPC + `SpringDamperElement3D` tie**는 단일 블록 대비 **`5.220e-11`** 일치 (`TRACK_C.md` S22) |
| S23 실측 (비정합 계면) | `MappingApplication` `nearest_element` mapper: 단위분할 `0.0`, 좌표/일정/운동학적 선형 패치 `5.55e-17`~`4.44e-16` 통과. 단 자유 계면 응력 패치는 비변분 오차 `0.086`. `AssignMasterSlaveConstraintsToNeighboursUtility`는 `7.94e-02`로 미달(제외) (`TRACK_C.md` S23) |
| S21 실측 (Newton 규칙) | 비선형 Newton 시 `MoveMeshFlag=True` + `AndCriteria(ResidualCriteria, DisplacementCriteria)` 필수 (`TRACK_C.md` S21) |
| THB 규약 | 레벨 $L$에서 $h_L = h/2^L$, 2:1 균형(이웃 레벨 차 ≤ 1), 활성 함수 $\text{supp} \subseteq \Omega^L \land \text{supp} \not\subseteq \Omega^{L+1}$, truncation(Giannelli 2012), $T = C_{\text{THB}} \cdot E$ Python 사전 합성 (v2.1 §7.2, §5.2.3) |
| Gmsh 격리 (D13) | WHT_IGA 프로세스 내 `import gmsh` **금지**(GPL 격리). `gmsh_runner.py`가 외부 프로세스 실행 → `.msh` → meshio 리더 → `MeshModel` (v2.1 §5.10.5) |
| 메쉬 쓰기 (D12) | `inp`, `fem`, `bdf`, `rad`, `k` 5종으로 메쉬 + **모든 Set 키워드** + 노드 결과(`U`, `NT11`) + 요소 평균 응력 출력. IGA 파트는 proxy hex-tet 격자 사용 (v2.1 §5.10.6) |
| Body-fit 범위 | all-hex 메쉬에서 정렬(mapped) 블록만 감지해 `NurbsVolumeGeometry` 변환 (sweep/bottom-up 불가) (v2.1 §5.10.4) |
| IGA–FE 혼합 (D14) | Phase 3 **Beta** 편입. `NativeMesh` 파트는 Kratos SM 일반 요소(`SmallDisplacementElement3D8N`, `3D10N` 등)로 동일 `Structure` 루트에 조립 (v2.1 §0.4 D14, §5.10.3) |
| 저장소·Probe (D17) | `SQLite` (`_iga_results.sqlite`, 메타·요약·Probe·인덱스) + `VTKHDF` (`<job>.vtkhdf`, 대용량 필드·시계열·proxy 격자·`/WHT_IGA/control_points`), `<job>_probe_data.json`, 가중치 0 적분점 probe (v2.1 §4.3, §4.5, §7.12) |

---

## 4. Phase 3 고유 함정 (반드시 숙지)

| # | 함정 | 올바른 처리 (근거) |
|---|---|---|
| K30 | 멀티파트 조립 시 노드·요소·Properties ID 오프셋 충돌 | `assembly_builder.py`에서 instance마다 연속 노드/요소 ID 구간을 배정하고 `<job>.resolved.json`에 `(instance, 구간, 출처: 제어점/FE노드/표면층/가상/RP)` 오프셋 표를 기록·추적한다 (v2.1 §5.9.2) |
| K31 | 비정합 계면에 `AssignNeighboursUtility`를 쓰거나 자유 계면 응력 패치에 collocation을 씀 | S23 실측대로 `AssignNeighboursUtility`(`7.94e-02` 미달)는 영구 제외한다. 운동학적 변위 패치는 `nearest_element`(`4.44e-16`), 자유 응력 계면은 Python Numba **mortar $L^2$ projection**으로 변분 일관 가중치를 계산한다 (`TRACK_C.md` S23, v2.1 §5.9.4) |
| K32 | IGA 조립체에서 `ContactSM` `mesh_tying_process`를 그대로 호출함 | `InterfacePreprocessCondition`이 `QuadraturePointGeometry::GenerateFaces()`를 호출해 `RuntimeError`로 실패한다. IGA가 낀 tie는 `Triangle3D3` 표면층 + 합성 MPC 또는 `SpringDamperElement3D` tie / Python mortar $L^2$를 사용한다 (`TRACK_C.md` S22, v2.1 §5.9.4) |
| K33 | IGA–IGA 또는 IGA–FE tie에서 가상 노드(스플라인 slave)를 다른 구속의 master로 씀 | 가상 노드는 제어점 자유도의 slave이므로 master가 될 수 없다(K11). IGA 면 가중치는 **제어점 자유도까지 Python에서 합성**해 한 줄의 `LinearMasterSlaveConstraint`로 만들거나, IGA–IGA는 penalty 스프링(BC-C) / 표면층 MPC를 쓴다 (v2.1 §5.9.4) |
| K34 | THB 세분화 시 2:1 균형 미준수 또는 truncation 누락으로 분할단위성 붕괴 | 세분 셀 마킹 후 이웃 셀 레벨 차 ≤ 1이 되도록 2:1 균형 전파를 강제하고, Giannelli 2012 truncation을 적용해 `Σ N_THB = 1` (1e-12)을 단위시험으로 검증한다. 레벨 $L$ `NurbsVolume`은 세분 영역 bbox로 잘라 노드 폭증을 막는다 (v2.1 §7.2) |
| K35 | `wht_iga` 코드에서 `import gmsh`를 직접 수행함 | GPL 라이선스 오염 방지(D13)를 위해 `gmsh` Python 바인딩 직접 import를 금지하고, `subprocess`로 외부 `gmsh` 프로세스를 실행해 `.msh`를 생성한 뒤 `meshio` 리더로 읽는다 (v2.1 §5.10.5) |
| K36 | 메쉬 쓰기 시 Set 키워드 누락 또는 1-based/로컬 면 번호 불일치로 왕복 손실 | `.inp`(`S1`–`S6`/`S1`–`S4` 변환), `.bdf`/`.fem`(`SET`/`THRU` 압축 + `SURF`/`BSURF` + `$$HMNAME SETS`), `.k`(`*SET_*_TITLE`), `.rad`(`GRNOD`/`GRBRIC`/`SURF`) 1:1 매핑 규약을 지키고 읽기→쓰기→읽기 왕복 시험으로 검증한다 (v2.1 §5.10.6) |
| K37 | `THIN_FEATURE` 자동 세분(G-5) 무한 루프 | 세분 재시도 상한(기본 최대 2회)과 최소 셀 크기 가드를 두고, 상한 초과 시 실패 사유와 함께 `issue`로 정직하게 보고한다 (v2.1 §9, add_1 G-5) |

---

## 5. 작업 목록 — Phase 3

각 작업은 **수락 기준을 전부 만족**해야 완료다. 완료 보고에는 실행한 명령과 출력 요약을 붙인다.

### T-P3.1 멀티셀 공유 장 + THB 계층 세분화 + 2:1 균형 (`wht_iga/grid/hierarchy.py`, `wht_iga/operators/compose.py`)

- 할 일:
  1. 멀티셀 공유 장: 한 파트 내 여러 재료 셀이 동일 스플라인 장(제어점)을 공유하도록 셀별 적분점·재료 마스크·`Properties` ID 매핑을 구현한다 (v2.1 §5.9.1, §7.7).
  2. `grid/hierarchy.py` (신규): 레벨 $L$ ($h_L = h/2^L$) 세분 영역 마킹(bbox / `CadEntitiesRegion` / offset 기본 $p \cdot h / 2^L$), 2:1 균형 전파(이웃 레벨 차 ≤ 1), 활성 함수 판정($\text{supp} \subseteq \Omega^L \land \text{supp} \not\subseteq \Omega^{L+1}$), Giannelli 2012 truncation, 1D 세분 행렬 텐서곱 기반 $C_{\text{THB}}$ 생성, 세분 영역 bbox 절단 레벨 $L$ `NurbsVolume` 생성 (v2.1 §7.2).
  3. `operators/compose.py` 확장: $T = C_{\text{THB}} \cdot E$ 합성(`scipy.sparse` CSR)으로 THB + small-cell extension을 단일 구속 행렬로 결합한다 (Phase 2의 `C_THB = I` 경로 확장).
- 수락 기준:
  - THB 기저 분할단위성 $\sum N_{\text{THB}} = 1$ (`1e-12`), $p$차 다항식 재생 (`1e-10`), 모든 인접 셀 레벨 차 ≤ 1.
  - $T = C_{\text{THB}} \cdot E$ 연쇄 구속 0건 단위시험 통과.
  - 멀티셀 공유 장 큐브/보 패치 테스트 해석해 대비 `1e-10` 이내.

### T-P3.2 멀티파트·어셈블리 + Coupling + 구속 그래프 검사 (`wht_iga/kratos_backend/assembly_builder.py`, `interactions.py`)

- 할 일:
  1. `assembly_builder.py` (신규): `Structure` 루트 아래 `Parts/<part>.<instance>`(`Cells`, `Layers`, `Virtual`, `Sets`), `Interactions/`, `Resolved/` SubModelPart 계층을 구축하고, instance별 연속 ID 구간 오프셋 표(`<job>.resolved.json`)를 생성한다 (v2.1 §5.9.2).
  2. Instance 구현(v2.1 §5.9.3): dependent instance는 파트 메쉬 1회 생성 후 W4에서 노드 좌표 $(R, t)$ 변환(IGA 적분점은 파라미터 좌표 유지), independent instance는 자기 위치에서 메쉬 생성. 한 패밀리 내 혼용 시 `validate()` 오류.
  3. `interactions.py` (신규): Tie 자동 선택(노드 일치 → 항등 MPC, 비정합 FE–FE → `nearest_element` 가중 MPC 또는 Python Numba mortar $L^2$ projection, IGA 포함 → 제어점 합성 MPC / `Triangle3D3` 표면층 + `SpringDamperElement3D` tie / BC-C penalty 스프링), RBE2(`ImposeRigidMovementProcess` 또는 MPC $\mathbf{u} = \mathbf{u}_{\text{RP}} + \boldsymbol{\theta} \times \mathbf{r}$), RBE3(가중 최소자승 Numba 계산 + MPC) 구현 (v2.1 §5.9.4, `TRACK_C.md` S22/S23).
  4. `job.resolve()` / `compile_plan()` 구속 그래프 검사: $c \neq 0$, 동일 DOF 이중 slave, slave→master 연쇄, 미지정 맞닿음 계면, 체적 겹침, 메쉬 소유 규칙 위반 검출 (v2.1 §5.9.7).
- 수락 기준:
  - dependent instance 3개 배치 모델: instance별 변위·응력 결과가 단일 instance 해석과 **`1e-10`** 이내 일치.
  - 비정합 계면 운동학적 패치 테스트 ≤ `1e-12`, 자유 응력 계면 mortar $L^2$ 패치 통과.
  - `job.resolve()` 연쇄 구속·이중 slave·$c \neq 0$ 거부 단위시험 통과.

### T-P3.3 일반 FE 파트 + IGA–FE 혼합 Beta(D14) + Gmsh 외부 생성 (`wht_iga/kratos_backend/native_fe.py`, `wht_iga/meshgen/gmsh_runner.py`)

- 할 일:
  1. `native_fe.py` (신규): `MeshModel`(`NativeMesh` 파트)의 `hex8`/`tet4`/`tet10`/`hex20`/`wedge6` 블록을 Kratos SM 일반 요소(`SmallDisplacementElement3D8N`, `3D10N` 등)로 루트 `Structure` ModelPart에 조립한다. IGA 파트와 같은 `Structure`에서 tie 결합되어 단일 시스템으로 풀리도록 IGA–FE 혼합 Beta 경로를 완성한다 (v2.1 §5.10.3, D14).
  2. `meshgen/gmsh_runner.py` (신규): STEP/BREP + 크기·차수(`tet4`/`tet10`/transfinite `hex`)·크기 필드 옵션을 받아 **외부 프로세스**(`subprocess`로 `gmsh` 실행 파일 또는 격리된 별도 Python 프로세스)를 호출해 `.msh`를 생성하고, `read_msh.py`를 통해 `MeshModel`로 적재한다. Physical group으로 CAD 면·솔리드 ID를 보존한다 (v2.1 §5.10.5, D13).
- 수락 기준:
  - `test_isolation.py`에 `import gmsh` 금지 검사 통과(WHT_IGA 주 프로세스 내 미사용).
  - Gmsh로 생성한 `tet4`/`tet10` 메쉬를 `NativeMesh`로 해석해 외팔보/패치 기준해 1% 이내.
  - IGA–FE 혼합 패치 및 V04 FE(Cu) + IGA(W) tie 사전 시험 통과.

### T-P3.4 Body-fit 변환 + 메쉬 쓰기 5종 + Set 키워드 왕복 무손실 (`wht_iga/meshio_ext/bodyfit.py`, `write_*.py`)

- 할 일:
  1. `meshio_ext/bodyfit.py` (신규, v2.1 §5.10.4): all-hex `MeshModel`에서 노드 차수·면 인접 그래프로 정렬(mapped) $(i, j, k)$ 블록을 감지하고 `NurbsVolumeGeometry`($p, C^{p-1}$)로 변환한다(비정렬/sweep 입력은 명확한 `issue`로 거부). 블록 간 경계 제어점 일치 시 공유, 불일치 시 tie.
  2. `meshio_ext/write_*.py` (`write_inp.py`, `write_nastran.py` 확장, `write_rad.py`, `write_k.py`, v2.1 §5.10.6 D12): `MeshModel` 및 `SelectionSet` 레지스트리로부터 5개 포맷(`.inp`, `.fem`, `.bdf`, `.rad`, `.k`)으로 메쉬 + **모든 Set 키워드** + 노드 결과(`U`, `NT11`) + 요소 평균 응력을 쓴다. IGA 파트는 proxy hex-tet 격자로 출력한다.
     - `.inp`: `*NSET`, `*ELSET`, `*SURFACE, TYPE=ELEMENT` (`S1`–`S6`/`S1`–`S4` 자동 변환)
     - `.fem`/`.bdf`: `SET <id>` (`THRU` 압축, 8/16컬럼/free-field), `SURF`/`BSURF`, `$$HMNAME SETS <id> "<name>"`
     - `.k`: `*SET_NODE_LIST_TITLE`, `*SET_SEGMENT_TITLE`, `*SET_SOLID_TITLE`
     - `.rad`: `/GRNOD/NODE`, `/GRBRIC`, `/SURF/SEG`
- 수락 기준:
  - 5개 포맷 모두 **읽기→쓰기→읽기 왕복 시험**에서 노드 좌표·요소 연결·모든 Set(이름, 종류, ID 목록) 무손실 일치.
  - Body-fit 변환 큐브/보 블록 해석해 대비 `1e-8` 이내.

### T-P3.5 IGES 입력 + `wht_iga.selection` 통합 API (`wht_iga/geometry/occ_io.py`, `wht_iga/selection/`)

- 할 일:
  1. `geometry/occ_io.py` 확장: IGES(`.igs`/`.iges`) 입력 + `BRepBuilderAPI_Sewing` / `MakeSolid` 치유 지원 (v2.1 §9).
  2. `wht_iga/selection/` 패키지 (신규, v2.1 §8.5):
     - `selector.py`: CAD(`faces/edges/vertices/solids`)와 Mesh(`nodes/elements/surface_tris/cells`) 이중 도메인 통합 `EntityQuery`, `free_surface_faces/elements/nodes()`, 공간/조건 필터(`within_box`, `within_sphere`, `within_cylinder`, `on_plane`, `between_planes`, `by_normal`, `by_area`, `by_length`, `nearest_to`, `flood_fill(angle_deg=30)`), OCP `GeomType` 필터(`by_geom_type("cylinder"|"plane"|"sphere"|"cone"|"torus"|"bspline", radius=...)`).
     - `sets.py`: `SelectionSet` (`.to_node_set()`, `.to_element_set()`, `.to_surface_set()`), 불리안 대수(`|`, `&`, `-`, `^`, `~`) 및 이종 변환 후 `materialize()`, 내보내기 레지스트리 연동.
     - `measure.py`: `bbox()`, `length()`, `area()`, `volume()`, `centroid()`, `inertia_tensor()`, `distance_to()`.
     - `contact_check.py`: `inspect_contact(master, slave)` (`KDTree`/`BVH` + 삼각형 투영 → `ContactGapReport`: `max_penetration`, `min_clearance`, `max_gap`, `normal_opposition_deg`, 접촉 면적 비율, `suggested_search_radius`) 및 `assy.detect_interfaces()`.
     - `GeomType` 자동 가이드 연결: 소형 `Cylinder`/`Torus`·고곡률 `BSpline` 감지 → `tessellation_size ≤ r/4` 및 THB `CadEntitiesRegion` 추천, `Cylinder` → `Frame.cylindrical` 추천, 6면 `Plane` → `BodyFitMesh` 추천.
- 수락 기준:
  - STEP/IGES 모델 + `MeshModel`에서 `EntityQuery`, `free_surface_*`, `SelectionSet` 불리안 차집합(`free_surf - tied_surf`) 후 `.inp`/`.fem` export 시 네이티브 Set 키워드로 정확히 출력됨.
  - `inspect_contact()`가 알려진 gap/침투 테스트 쌍에서 기하 해석값과 `1e-6` 이내 일치.

### T-P3.6 다단계 Step + Amplitude + 진행 콜백(G-2) + `THIN_FEATURE` 자동 세분(G-5) + 출력 이름 매핑(G-7)

- 할 일:
  1. 다단계 `Step` 및 `Amplitude`: 이전 step 변위·내부 상태 승계와 Python 시간함수 `Amplitude` 기반 하중·BC 갱신 경로를 구현한다 (v2.1 §9).
  2. Job 진행 콜백·취소(add_1 G-2): `AnalysisSession.run(on_progress=..., cancel_token=...)`에서 단계/스텝/수렴 진행률 콜백 호출 및 취소 토큰 감지 시 안전 중단 + 부분 결과 정리.
  3. `THIN_FEATURE` 자동 세분(add_1 G-5): Phase 2 `extension.py`가 `THIN_FEATURE` issue를 반환하면 해당 셀/영역을 THB 세분 영역으로 등록해 자동 재격자화한다(최대 재시도 2회 가드).
  4. 출력 이름 매핑표(add_1 G-7): Coreform 별칭(`U`, `S`, `E`, `NT11`, `RF`) ↔ Kratos 네이티브 변수명(`DISPLACEMENT`, `CAUCHY_STRESS_TENSOR`, `VON_MISES_STRESS`, `TEMPERATURE`, `REACTION`) 양방향 매핑표를 단일 소스로 통합해 `ResultView`·`Probe`·내보내기에서 공통 사용한다 (v2.1 §8.4).
- 수락 기준:
  - 2단계 하중-제하/다중 하중 step 승계 단위/통합 시험 통과.
  - 취소 토큰 발동 시 콜백 이후 즉시 중단되고 좀비 프로세스/손상 DB 없음.
  - 얇은 판 형상에서 `THIN_FEATURE` 감지 → 자동 THB 세분 1회 후 extension 성공 확인.

### T-P3.7 가중치 0 Probe + `VTKHDF` Proxy 격자 + $L^2$ 노름 (`wht_iga/post/`, `wht_iga/query/results.py`)

- 할 일:
  1. `post/probes.py`: `PointProbe`, `LineProbe`, `MultiPointProbe`를 Kratos 가중치 0 적분점으로 평가하여 `<job>_probe_data.json`과 `_iga_results.sqlite`에 기록한다 (v2.1 §4.3, §4.5, §7.12).
  2. `post/proxy_mesh.py` + `post/vtk_writer.py`: IN 셀 hex / CUT 셀 tet proxy 격자 생성 및 `h5py` 기반 `<job>.vtkhdf` 작성(instance별 블록 + IGA proxy + FE 요소 + 표면층 통합 파일 + `/WHT_IGA/control_points/<frame>/{U, NT11, state}` 저장).
  3. `post/norms.py`: $L^2$, $H^1$, 에너지 오차 노름 계산기(참조해/해석해 콜백 대비).
  4. `query/results.py`: `ResultView` (`probe`, `field`, `reaction_force`, `error_norm`, `check_equilibrium`, `section_force`, `hotspots`, `stress_linearize`)로 `SQLite + VTKHDF` 통합 조회 구현 (v2.1 §8.1, §8.5.6).
- 수락 기준:
  - `VTKHDF` 파일 `h5py` 스키마 및 ParaView 호환 구조 검증, proxy 격자 노드/셀 필드와 제어점 재평가 일치 (`1e-12`).
  - `norms.py` $L^2$/$H^1$ 노름이 해석해 다항식 테스트에서 `1e-12` 정확도.
  - `check_equilibrium()` 외력 합력 vs 반력 합력 평형 오차 `1e-10` 이내.

### T-P3.8 게이트 G3 통합 검증 (검증, 코드 변경 최소화)

- 할 일: v2.1 §9 Phase 3 및 §10.1–10.2 기준 전량 실행 및 리포트 자동 생성:
  1. **V03 moment**: 보 휨 moment 케이스 끝단 처짐 Euler–Bernoulli 대비 **5% 이내**.
  2. **V04 bimetallic Cu/W**: merge vs tie 끝단 처짐 상대차 **≤ 1e-6**, 변환단면 해 대비 오차 `bimetallic-beam-tie-comparison.csv` 수준(p = 3 ≈ `0.54%`) ± `0.2%p`. 추가로 **FE 파트(Cu) + IGA 파트(W) tie** 해석을 수행해 merge 대비 **0.5% 이내** 확인.
  3. **V05 multiplanar 계면 7단면 (84 해석)**: monolithic 대비 차이가 Coreform 보고 수준(p = 2 ≤ **`3.3e-6`**) 이하.
  4. **V02 국소 세분 (THB)**: L3에서 peak 오차 **≤ 0.2%**, DOF **≤ 80k** (`roark-17p1-7a-refinement.csv`, Coreform `0.105%` @ `57,300` DOF 참고).
  5. **E01 박판 구멍 (IPS, 10×5×0.1, r=1)**: Roark $K_t = 2.23568$ 기준 $\sigma_{\max} = 3.7261 \times 10^{-6}\text{ psi}$ 대비 probe 오차 ≤ Coreform(`−8.7%`, `3.4031e-6`) 수준, 세분 시 개선 추세.
  6. **E02 plate-with-hole-adaptivity**: 문서 정성 결과와 세분/오차 감소 추세 일치.
  7. **추가 게이트 2종**: dependent instance 3개 배치 모델 `1e-10` 일치, 메쉬 리더/라이터(5종 Set 키워드 포함) 왕복 무손실.
- 수락 기준: 위 7개 항목 전량 통과. 결과는 `WHT_IGA/docs/verification/` (`V03-moment.md`, `V04.md`, `V05.md`, `V02-local.md`, `E01.md`, `E02.md`)로 자동 생성 (형식 §8).

### T-P3.9 `HWT-IGA SOLVER` 상용 로깅 · 실시간 메모리 기반 지능형 솔버 전환 · `out_mode`/`AMGCL` 통합 (`wht_iga/solver/memory.py`, `solver_log.py`, `pipeline.py`, `validate.py`, `kratos_backend/solvers.py`)

- 할 일:
  1. **실시간 시스템 메모리 계측 & 지능형 솔버 전환 (`wht_iga/solver/memory.py`, 신규, v2.1 §5.11.3)**:
     - 추가 외부 패키지 없이 Windows(`ctypes.windll.kernel32.GlobalMemoryStatusEx`) / POSIX(`/proc/meminfo`)에서 `total_ram_gb`, `avail_ram_gb`를 실시간 조회한다(`WHT_SOLVER_AVAIL_RAM_GB` 테스트 오버라이드 지원).
     - 활성 DOF $N_{\text{act}}$, 추정 $\text{nnz}(K)$, 단면 프론트 폭 $F_{\text{eff}}$ 기반으로 직접 솔버(`pardiso_ldlt`/`pardiso_lu`) vs 반복 솔버(`amgcl`) 예상 피크 RAM을 산출한다.
     - `linear_solver="auto"`일 때 `safe_budget_gb = 0.65 * avail_ram_gb`를 기준으로 MKL 스레드 축소(`num_threads_mkl=4`) → 안전 예산 또는 `MKL_INT` 한계 초과 시 **`KratosMultiphysics.AMGCLSolver`로 자동 전환**한다.
     - **전환 발생 시 반드시 정량 수치 근거(총/가용 RAM, 안전 예산, 활성 DOF, 직접 vs 반복 추정 RAM, 선택된 솔버/전처리기)를 포함한 `[HWT-IGA SOLVER][AUTO-SOLVER]` 전환 메시지를 콘솔·로그·`ValidationIssue`/`ResolveReport`에 출력한다.** (`extension=None` + 절단 셀 상태에서 전환 시 조건수 경고 병기).
  2. **`out_mode="drop"|"ghost"` 및 `AMGCLSolver`/`CGSolver` 파이프라인 통합 (`validate.py`, `pipeline.py`, `kratos_backend/solvers.py`, v2.1 §5.11.1)**:
     - `JobSpec`에 `out_mode: Literal["drop", "ghost"] = "drop"`, `linear_solver` 선택지(`"auto"`, `"pardiso_lu"`, `"pardiso_ldlt"`, `"pardiso_llt"`, `"amgcl"`, `"cg"`, `"skyline_lu"`), `amgcl_params` 옵션을 추가한다.
     - `out_mode="drop"` 시 `OUT` 셀 요소를 생성하지 않고 `INACTIVE` 제어점 DOF를 고정(`Fix`)해 불필요 DOF·적분점을 제거한다(`"ghost"`는 기존 FCM 기준선 유지).
  3. **`HWT-IGA SOLVER` 상용 수준 6단 로깅 & On/Off 제어 (`wht_iga/solver/solver_log.py`, 신규, v2.1 §5.12)**:
     - `JobSpec` / `StudySpec` / `AnalysisSession`에 `enable_solver_log: bool = True`, `solver_log_file`, `solver_log_stdout: bool = True`, `solver_log_level: Literal["minimal", "standard", "verbose"] = "standard"` (및 환경변수 `WHT_SOLVER_LOG=0|1`)를 탑재한다.
     - `enable_solver_log=False` 시 배너/로그 파일 출력을 완전히 생략하고, `True` 시 공식 솔버 명칭 **`HWT-IGA SOLVER`**로 다음 6단 정보를 타임스탬프·정렬 표 형식으로 출력한다:
       1. Banner & Host/RAM/Kratos/MKL Environment
       2. Problem & Discretization Summary (파트/차수/격자/THB/`IN·CUT·OUT`/`out_mode`/활성·고정·MPC DOF/물성/하중·모멘트 사전 검산)
       3. Immersed IGA & Assembly Diagnostics (테셀레이션·`NNLS`/`QUAD_RELAXED`/`QUAD_FALLBACK`·`extension` `MASTER/EXTENDED/INACTIVE`·`THIN_FEATURE`·구속 그래프)
       4. Memory Audit & Auto-Solver Decision (가용 RAM 대비 피크 추정치 표 + `[AUTO-SOLVER]` 전환 메시지 + 최종 솔버/스레드 설정)
       5. `[STAGE 1/5]`–`[STAGE 5/5]` 진행률·RSS 메모리·Newton/Step/Krylov 잔차 수렴 이력 모니터
       6. Equilibrium Check($\sum F^{\text{ext}}$ vs $\sum R$) + Max/Peak 응답 + Wall-Time Breakdown + Peak RSS Footer (`COMPLETED` / `COMPLETED_WITH_WARNINGS` / `FAILED`)
- 수락 기준:
  - 단위시험(`test_phase3_solver_memory_log.py`): `probe_system_memory()` / 피크 추정 모형 / `WHT_SOLVER_AVAIL_RAM_GB` 압박 시 `auto` → `amgcl` 전환 판단 및 전환 메시지 검증 / `enable_solver_log=True|False` 파일·stdout On/Off 동작 검증.
  - Kratos 통합 시험: `out_mode="drop"` vs `"ghost"` 일치 확인 및 `amgcl` 솔버로 V02/V06 해석해 대비 `pardiso_ldlt` 대비 `max|ΔU|/max|U| ≤ 1e-7`·힘 평형 ≤ `1e-8` 통과.

---

## 6. 산출 모듈 목록 (신규/확장 파일)

| 모듈 | 작업 | 비고 |
|---|---|---|
| `wht_iga/grid/hierarchy.py`, `operators/compose.py` | T-P3.1 | THB, 2:1 균형, truncation, $T = C_{\text{THB}} \cdot E$ |
| `wht_iga/kratos_backend/assembly_builder.py`, `interactions.py` | T-P3.2 | `ModelPart` 계층, instance, tie/RBE2/RBE3, `job.resolve()` |
| `wht_iga/kratos_backend/native_fe.py`, `wht_iga/meshgen/gmsh_runner.py` | T-P3.3 | IGA–FE 혼합 Beta(D14), Gmsh 외부 프로세스(D13) |
| `wht_iga/meshio_ext/bodyfit.py`, `write_*.py` | T-P3.4 | Body-fit 변환, 5종 메쉬+Set 키워드 라이터(D12) |
| `wht_iga/geometry/occ_io.py`, `wht_iga/selection/{selector,sets,measure,contact_check}.py` | T-P3.5 | IGES 입력, `EntityQuery`, `SelectionSet` 불리안, `inspect_contact`, `GeomType` |
| `wht_iga/solver/` (step/amplitude/callback/auto-refine/var-map) | T-P3.6 | 다단계, G-2 콜백·취소, G-5 `THIN_FEATURE` 세분, G-7 매핑표 |
| `wht_iga/solver/{memory,solver_log,validate,pipeline}.py`, `kratos_backend/solvers.py` | T-P3.9 | **`HWT-IGA SOLVER` 6단 로거(On/Off)**, **실시간 RAM 기반 `auto`→`amgcl` 전환·메시지**, `out_mode` |
| `wht_iga/post/{probes,proxy_mesh,vtk_writer,norms}.py`, `wht_iga/query/results.py` | T-P3.7 | 가중치 0 probe, `VTKHDF`, $L^2$/$H^1$ 노름, `ResultView` |
| `WHT_IGA/tests/unit/test_phase3_*.py` | T-P3.1–P3.7 | Kratos 불필요 단위시험 (THB, 불리안 Set, 왕복 등) |
| `WHT_IGA/docs/verification/{V03-moment,V04,V05,V02-local,E01,E02}.md` | T-P3.8 | 형식 §8 |

---

## 7. 코딩 규칙 (Phase 1–2와 동일 — Phase 3 추가만 아래에)

- Python 3.12, 타입 힌트, `from __future__ import annotations`.
- Numba `@njit(cache=True)` 순수 함수, 일관 무차원 단위계, 셀 선형 인덱스 규약 유지, 파일 400줄 초과 시 분할.
- **추가 (Phase 3)**:
  - Kratos `import`는 `wht_iga/kratos_backend/` 내부로만 제한한다.
  - `gmsh`는 절대 `import`하지 않고 `subprocess`로만 실행한다 (`test_isolation.py`에 검증 추가).
  - 대용량 필드 시계열은 `VTKHDF`(`h5py`)에 연속 배열로 쓰고, `SQLite`에는 메타/Probe/인덱스만 기록한다.
  - 모든 어셈블리 구속은 `compile_plan()` / `job.resolve()`에서 그래프 검사(연쇄·이중 slave·$c \neq 0$)를 통과한 뒤에만 Kratos에 전달한다.

---

## 8. 보고 형식

### 8.1 작업 기록 (`docs/phase3/P3NN.md`)

```markdown
# P3NN <제목>
- 작업: (T-P3.x)
- 결과: 완료 / 부분 / 실패
- 근거: (실행 명령, 핵심 출력 숫자, 수락 기준별 통과/실패)
- 환경: Python x.y.z, Kratos 10.4.3 (경로), 날짜
- 결론과 영향: (G3 게이트에 미치는 영향, Phase 4로 넘기는 항목)
- 재현: pytest 경로 또는 실행 명령
```

### 8.2 검증 리포트 (`docs/verification/<ID>.md`, v2.1 §10.6)

- 케이스별 자동 생성: 입력 요약, 기하 오차, 수렴 표·그래프, 기준값 대비 오차, 통과 여부, 실행 시간, DOF 수.
- `V04`는 merge vs tie 및 FE(Cu)+IGA(W) 비교표를, `V05`는 84해석 요약표를, `V02-local`은 L1–L3 DOF·peak 오차표를 포함한다.

### 8.3 작업 완료 보고 (작업마다)

- 작업 ID, 변경·생성 파일 목록, 실행 명령, 수락 기준별 통과/실패와 근거 숫자, 남은 문제.
- 실패를 성공처럼 쓰지 않는다. 미수행 항목은 "미수행"으로 적는다.

---

## 9. 정의: 이번 요청의 완료 (Definition of Done)

1. `T-P3.1` – `T-P3.9` 모든 수락 기준 통과, 또는 실패 항목이 근거와 함께 보고됨
2. `micromamba run -n whtiga python -m pytest WHT_IGA/tests -m "not kratos"` 통과 (신규 단위시험 및 `test_isolation.py` 포함)
3. Kratos 필요 시험 (`-m kratos`) 결과 첨부 (실패가 있으면 그대로 보고)
4. 변경 파일에 TODO 스텁, `skip`, `xfail` 없음 (있으면 사유와 함께 보고)
5. G3 게이트 증거 제출 → **사용자 승인 대기**. 승인 전 Phase 4 착수 금지.

### 9.1 Phase 4 이후 미리보기 (이번 요청 범위 아님)

| Phase | 핵심 산출물 | 게이트 (v2.1 §9) |
|---|---|---|
| 4 | Python consistent M, scipy 고유치, Python Bossak 선형 과도, 모드 중첩 과도/조화, 스펙트럼 case-file 리더, 열·열변형 | V09, V10, V11, V12 |
| 5 | NLGEOM, 초탄성, J2, C++ `WhtIgaApplication` MixedUP0, 표면층 접촉 | V07, V08, E03, E04, E05 |
| 6 | 성능 (C++ 이전, `AdditiveSchwarzPreconditioner` 고도화, 병렬), Explicit 트랙 | 1M DOF < 30분, V-EXP, V-EXP-C |

---

## 10. 멈추고 사용자에게 물을 조건

다음 중 하나라도 생기면 **추측으로 진행하지 말고** 상황, 시도한 것, 선택지를 정리해 묻는다.

1. THB 세분화 후 2:1 균형 전파로 인해 DOF가 G3 한계(V02 L3 ≤ 80k)를 초과할 때
2. 비정합 계면 mortar $L^2$ projection과 `Triangle3D3`+`SpringDamperElement3D` tie 간 수렴·정확도 차이가 V04(`1e-6`) 또는 V05(`3.3e-6`) 기준을 위협할 때
3. `whtiga` env에 `gmsh` 실행 파일 또는 `h5py`가 없어 설치/경로 지정이 필요할 때
4. 메쉬 쓰기 5종(`.inp`, `.fem`, `.bdf`, `.rad`, `.k`)에서 특정 솔버 방언의 Set 키워드 왕복 무손실이 깨질 때
5. `THIN_FEATURE` 자동 세분 재시도 상한(2회) 내에서도 유효한 extension 블록을 찾지 못할 때
6. v2.1 내용이 Kratos 실제 동작과 다를 때 (v2.1은 수정하지 말고 보고)
7. 작업량이 예상(8주)의 2배를 넘을 것으로 보일 때

---

## 11. 참고 위치

| 내용 | 위치 |
|---|---|
| 계획서 | `research-coreform-IGA/master-plan_v2.1.md` (§5.9, §5.10, §7.2, §8.5, §9 Phase 3, §10.1–10.2) |
| 트랙 C 실측 (S22, S23 등) | `WHT_IGA/docs/spikes/TRACK_C.md` |
| Phase 2 구현 요청서 (선행 양식) | `research-coreform-IGA/implementation-request_phase2.md` |
| Phase 1–2 산출물 (하위호환 유지) | `WHT_IGA/wht_iga/`, `WHT_IGA/docs/phase1/`, `WHT_IGA/docs/phase2/` |
| 검증 기준값 JSON | `WHT_IGA/tests/verification/data/` (`v02*`–`v05*`, `e01*`, `e02*`) |
