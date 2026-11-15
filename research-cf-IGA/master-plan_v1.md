Kratos 기반 고차 잠김 등기하해석(Immersed IGA) 통합 솔버 종합 R&D 마스터플랜 및 아키텍처 설계서

1. 프로젝트 비전 및 시스템 엔지니어링 철학

전산역학(Computational Mechanics) 및 컴퓨터 지원 공학(CAE) 분야에서 컴퓨터 가공 설계(CAD)의 경계 표현형(B-Rep) 모델을 유한요소해석(FEA)용 3차원 체적 격자(Conformal Hex/Tet Mesh)로 변환하는 전처리 작업은 전체 해석 공정 시간의 80% 이상을 소모하는 주된 병목 구간으로 지속적인 한계가 지적되어 왔습니다. 이러한 기하-해석 간의 단절을 극복하기 위해 제안된 등기하해석(Isogeometric Analysis, IGA)은 CAD의 B-스플라인 및 NURBS 기저 함수를 해석의 시가지 함수로 직접 활용하여 기하학적 정밀도를 이산화 단계까지 완벽히 보존하는 패러다임을 확립하였습니다. 그러나 3차원 복잡 형상에 대해 위상적으로 정합하는 체적 NURBS 패치를 자동으로 구성하는 과정은 여전히 수학적·토폴로지적 난제로 남아 있습니다.

본 프로젝트의 목표는 3차원 CAD B-Rep 형상을 정합 체적 격자 생성 없이 직교 B-스플라인 및 계층적 THB-스플라인 배경 격자에 직접 침강(Immerse)시켜 해석하는 Kratos Multiphysics 기반 고차 잠김 등기하해석(Immersed IGA) 완전 통합 솔버를 개발하는 것입니다. 상용 Coreform IGA 수준의 해석 정밀도, 대수적 수치 안정성, 그리고 엔터프라이즈 CAE 연동성을 달성하기 위해 다음과 같은 3대 핵심 엔지니어링 원칙을 견지합니다.

1.1 Zero C++ Build Overhead (솔버 커널 무재컴파일 원칙)

수치 솔버 개발 시 C++ 코어의 직접 수정 및 재컴파일은 개발 파이프라인의 복잡성을 가중시키고 플랫폼 종속성을 초래합니다. 본 시스템은 글로벌 오픈소스 전산역학 프레임워크인 KratosMultiphysics의 공식 배포 Wheel 바이너리(KratosMultiphysics, IgaApplication, StructuralMechanicsApplication)를 단 한 줄의 C++ 소스 재컴파일 없이 그대로 탑재합니다. 비적합 격자 해석에서 필연적으로 발생하는 소형 절단 요소(Small Cut-Cell)의 수치 불안정성과 조건수 폭발 문제는 솔버 내부 커널을 수정하는 대신, **전처리 단계에서 베지에 추출 연산자를 대수적으로 직접 변환(\mathbf{C}^e \to \mathbf{\tilde{C}}^e)**하여 완전히 제어합니다.

1.2 Polyglot Performance (이기종 고속 컴퓨팅 아키텍처)

각 연산 단계의 특성에 최적화된 프로그래밍 언어와 가속 라이브러리를 유기적으로 결합합니다.

* Python + Numba + splinepy: 수학적 알고리즘의 신속한 프로토타이핑, 고차 스플라인 그리드 생성, 워크플로우 오케스트레이션 및 데이터 버스 제어.
* Rust (PyO3 + Rayon + sprs + rusqlite + ryu): 3차원 공간 트리 탐색, 대규모 배경 격자의 내외부 고속 판정, 소형 셀 그래프 병합 및 \mathbf{T} 행렬 조립, SQLite 바이너리 BLOB I/O 스트리밍, 대용량 CAE 결과 파일 고속 직렬화.
* Intel MKL PARDISO (PyPardiso): 고차 스플라인 요소로 인해 대역폭이 넓어진 대규모 희소 행렬(Sparse CSR Matrix)을 멀티스레드 직접 분해(Direct Factorization)하여 솔버 풀이 성능을 극대화.

1.3 Dual Enterprise Postprocessing (개방형·상용 듀얼 후처리 파이프라인)

* 연구 및 오픈소스 표준 (Track A): HDF5 기반 ParaView 고속 시계열 시각화 포맷(cad_evaluated.vtkhdf)을 생성하여 대용량 물리장의 분산 렌더링 지원.
* 상용 엔터프라이즈 표준 (Track B): OptiStruct/Nastran 벌크 데이터 포맷(.fem)과 Altair 결과 포맷(.hvascii)을 작성한 후, Altair hvtrans CLI 유틸리티를 백그라운드 호출하여 고압축 바이너리인 HyperView 포맷(cad_evaluated.h3d)으로 자동 변환.
* 정량적 설계 엔지니어링 (Track C): CAD 표면 및 내부 관심 지점에 대한 고정밀 B-스플라인 역매핑(Point Inversion)을 수행하여 응력 집중 계수 및 변위 데이터를 담은 probe_results.json 출력.

2. 엔드투엔드 시스템 4단계 계층 구조 및 중간 데이터 버스

전체 시스템은 모듈 간의 결합도를 최소화하고 데이터 교환을 표준화하기 위해 명확한 4단계 계층 구조와 바이너리 캐시 버스로 구성됩니다.

                     [ 원본 3D CAD (.step / .iges / STL) ]
                                       │
                                       ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Stage 1: Model Setup (Python, PythonOCC, Gmsh API)                      │
 │  • Dirty CAD 치유(Shape Healing, Sewing) & Watertight 검증                │
 │  • 경계면(Face ID) 태깅: Dirichlet(변위), Neumann(하중), Weak BC         │
 └─────────────────────────────────────┬────────────────────────────────────┘
                                       │ iga_config.json & ProjectParameters.json
                                       ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Stage 2-A: IGA Preprocessor Engine (Rust PyO3 + Algoim + splinepy)      │
 │  • [Knot / THB Grid] 직교 배경 격자 구축 및 베지에 추출 (C^e)            │
 │  • [Algoim / QuESo] 레벨셋 절단 고차 수치 적분점/가중치/법선 산출         │
 │  • [AgFEM / Ext-B-spline] 소형 셀 그래프 병합: C^e ──▶ C̃^e 변환         │
 └─────────────────────────────────────┬────────────────────────────────────┘
                                       │
         ┌─────────────────────────────┴─────────────────────────────┐
         ▼                                                           ▼
   [mesh_data.sqlite]                                      [bezier_grid.mdpa]
   • 요소별 가변 C̃^e (BLOB)                                • Kratos IgaApp 표준 포맷
   • 절단 적분점/가중치/법선벡터 (BLOB)                    • Nodal Coordinates & BC Elements
         │                                                           │
         └─────────────────────────────┬─────────────────────────────┘
                                       │
                                       ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Stage 2-B: Simulation Execution (Kratos Wheel + PyPardiso)              │
 │  • Kratos IgaApplication (C̃^e 주입) + StructuralMechanicsApplication     │
 │  • Kratos Sparse CSR 행렬 추출 ──▶ PyPardiso (Intel MKL 직접 솔버)       │
 └─────────────────────────────────────┬────────────────────────────────────┘
                                       │ control_points.hdf5 (해석 변위 벡터)
                                       ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Stage 3: Dual Postprocessor Engine (Rust PyO3 + ryu + h5py)             │
 │  • 평가 격자(Sub-mesh) 생성 & Bernstein 다항식 기반 물리장 내삽/투영     │
 ├──────────────────────────┬───────────────────────────┬───────────────────┤
 │ Track A: Open-Source     │ Track B: Enterprise CAE   │ Track C: Probing  │
 │ • h5py / VTKHDF 생성     │ • .fem + .hvascii 생성    │ • 점 역매핑 계산  │
 │ • ParaView 시계열 시각화 │ • Altair hvtrans ──▶ .h3d │ • 응력/변형 보고  │
 └──────────────────────────┴───────────────────────────┴───────────────────┘


2.1 중간 데이터 버스 명세

mesh_data.sqlite 데이터베이스 스키마

전처리 엔진과 Kratos 솔버, 그리고 후처리 엔진 간의 데이터 결합을 분리하기 위해 SQLite 데이터 버스를 채택합니다. 가변 크기의 적분점 세트와 대수적 추출 행렬은 원시 부동소수점 바이트 스트림(BLOB)으로 패킹되어 무복사(Zero-Copy) 직렬화를 실현합니다.

-- 요소 메타데이터 및 수정된 베지에 추출 행렬 캐시
CREATE TABLE elements (
    elem_id INTEGER PRIMARY KEY,
    status INTEGER NOT NULL,     -- 0: Exterior (Void), 1: Interior, 2: Cut
    num_qp INTEGER NOT NULL,      -- 할당된 총 적분점 개수
    c_tilde BLOB NOT NULL         -- Float64 2D Array [n_bernstein, n_bspline] BLOB
);

-- 절단 요소 내부 체적 및 경계 적분점 좌표/가중치
CREATE TABLE quadrature_points (
    elem_id INTEGER NOT NULL,
    qp_id INTEGER NOT NULL,
    x REAL NOT NULL, y REAL NOT NULL, z REAL NOT NULL, -- 물리 좌표
    weight REAL NOT NULL,                              -- 적분 가중치
    nx REAL NOT NULL, ny REAL NOT NULL, nz REAL NOT NULL, -- 법선 벡터 (경계점: 단위법선, 체적점: 0.0)
    PRIMARY KEY (elem_id, qp_id)
);
CREATE INDEX idx_qp_elem ON quadrature_points(elem_id);

-- 모델 전역 파라미터 및 스플라인 메타정보
CREATE TABLE model_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);


bezier_grid.mdpa (Kratos IgaApplication 인터페이스)

Kratos의 IgaApplication이 인식하는 모델 부분(ModelPart) 포맷으로, 베지에 제어점과 연결성, 대수적으로 안정화된 \mathbf{\tilde{C}}^e 연산자를 포함합니다.

Begin ModelPartData
    DENSITY 7850.0
    YOUNG_MODULUS 2.1e11
    POISSON_RATIO 0.3
End ModelPartData

Begin Nodes
    1  0.0  0.0  0.0
    2  1.0  0.0  0.0
    // 베지에 추출 노드 좌표 목록
End Nodes

Begin Elements Element3D8N // 3차원 솔리드 스플라인 베지에 요소
    1  1  1 2 3 4 5 6 7 8
End Elements

Begin NodalData EXTRACTION_OPERATOR_MCSR
    // C̃^e 희소 행렬 인덱스 및 가중치 데이터
End NodalData


3. 모듈식 디렉토리 및 파일 아키텍처 명세

본 솔버의 완전한 Zero C++ 빌드 파이프라인, Rust 가속 코어, Kratos 어댑터, 듀얼 후처리 엔진, 그리고 이원화된 테스트 스위트를 구현하는 전체 소스 코드 디렉토리 및 모듈 배치 구조는 다음과 같습니다:

WHT_IGA/
├── Cargo.toml                       # Rust 확장 모듈(rust_iga_core) 빌드 설정 (maturin / PyO3)
├── pyproject.toml                   # Python 패키지 배포 및 의존성 정의 (kratos, pypardiso, splinepy, gmsh, pythonocc)
├── setup.py / Makefile              # 로컬 개발 환경 원클릭 설치 스크립트
├── README.md                        # 설치 가이드 및 아키텍처 개요
│
├── rust_iga_core/                   # ─────────────────────────────────────────────────────────────
│   │                                # [Rust 가속 코어: PyO3 + Rayon + sprs + rusqlite + ryu]
│   ├── Cargo.toml                   # Rust 크레이트 의존성 설정
│   └── src/                         #
│       ├── lib.rs                   # Python 바인딩 모듈 엔트리포인트 (PyO3)
│       ├── spatial_culling.rs       # Rayon 병렬 AABB/BVH 배경 격자 내외부 판정 (In/Out/Cut)
│       ├── algoim_bridge.rs         # Algoim C-API 바인딩 및 가우스 적분 규칙 생성기
│       ├── aggregation_graph.rs     # 소형 절단 셀 그래프 병합 & T 행렬 조립 (C^e -> C̃^e 변환)
│       ├── sqlite_storage.rs        # rusqlite 기반 mesh_data.sqlite 고속 BLOB 스트리밍 입출력
│       └── hvascii_writer.rs        # ryu 고속 부동소수점 포맷팅 기반 Altair .hvascii 생성기
│
├── wht_iga/             # ─────────────────────────────────────────────────────────────
│   │                                # [Coreform & Kratos 스타일 고수준 통합 Python 패키지]
│   ├── __init__.py                  # ImmersedModel, Analysis, Material, WeakBC 등 최상위 API 노출
│   ├── config.py                    # 경로, 임계값(eta_crit=1e-2), 전역 파라미터 정의
│   │
│   ├── cad/                         # Stage 1: CAD 처리 및 기하 질의 엔진
│   │   ├── __init__.py
│   │   ├── occ_importer.py          # pythonOCC 기반 STEP/IGES 로드, ShapeFix_Shape(치유) 및 Watertight 검증
│   │   ├── boundary_tagger.py       # Face ID별 물리적 경계(Dirichlet, Neumann, Interface) 태깅
│   │   └── surface_mesher.py        # Gmsh Python API 연동: CAD B-Rep 곡률 준수 2D 표면 삼각망 생성
│   │
│   ├── preprocessor/                # Stage 2-A: IGA 배경 격자 및 절단 수치적분 전처리
│   │   ├── __init__.py
│   │   ├── background_grid.py       # splinepy 기반 텐서곱 B-스플라인 / THB 배경 격자 구축 및 기본 C^e 추출
│   │   ├── cut_quadrature.py        # Algoim 및 QuESo 연동 절단 체적/경계 고차 적분점(QP) 산출
│   │   ├── stabilizer.py            # Rust aggregation_graph 호출: C^e -> C̃^e 대수적 변환 오케스트레이션
│   │   ├── db_manager.py            # mesh_data.sqlite BLOB 패킹, 압축, 조회 관리자
│   │   └── mdpa_writer.py           # Kratos IgaApplication 표준 bezier_grid.mdpa 직렬화 모듈
│   │
│   ├── solver/                      # Stage 2-B: Kratos 솔버 및 직접 해법 연동
│   │   ├── __init__.py
│   │   ├── kratos_adapter.py        # Kratos ModelPart 생성, C̃^e 주입 및 StructuralMechanicsApp 연동
│   │   ├── weak_bcs.py              # Nitsche 기법 및 페널티 기반 약형식 경계조건(Weak Dirichlet BC) 부과
│   │   ├── pardiso_solver.py        # Kratos Sparse CSR 행렬 추출 <-> PyPardiso (Intel MKL) 고속 직접 해법
│   │   └── analysis.py              # User-facing 해석 실행기 (Run, Step loop, HDF5 결과 제어)
│   │
│   └── postprocessor/               # Stage 3: 이중(Dual) 포스트프로세싱 파이프라인
│       ├── __init__.py
│       ├── submesh_evaluator.py     # Rust 가속: 서브메시 평가점에 번스타인 다항식 물리장(변위/응력) 투영
│       ├── vtkhdf_writer.py         # Track A: h5py 기반 ParaView 시계열 .vtkhdf 생성 모듈
│       ├── altair_bridge.py         # Track B: .fem + .hvascii 작성 및 Altair hvtrans CLI 호출 (.h3d 변환)
│       └── probe.py                 # Track C: 특정 공간 좌표 B-스플라인 역매핑 및 probe_results.json 출력
│
└── tests/                           # ─────────────────────────────────────────────────────────────
    │                                # [Coreform IGA 벤치마크 기반 이원화 테스트 스위트]
    ├── 1_verification/              # 엄밀 검증군 (Coreform Verification Manual 기반, 해석해 비교)
    │   ├── test_v01_patch_test.py       # Prescribed displacement of unit cube (기초 패치 테스트)
    │   ├── test_v02_kirsch_plate.py     # Infinite plate with a hole (Kirsch 해석해, Kt=3.0)
    │   ├── test_v03_roark_hole.py       # Roark 17.1-7a: Finite thin plate with hole (응력집중계수)
    │   ├── test_v04_timoshenko_beam.py  # Cantilever beam bending (Timoshenko 휨 해석해)
    │   ├── test_v05_thick_cylinder.py   # Thick-walled cylinder under pressure (Lamé 3D 해석해)
    │   ├── test_v06_sliver_cut_cond.py  # Extreme small cut-cell (η -> 10^-8 조건수 불변성 검증)
    │   ├── test_v07_beam_dynamics.py    # Cantilever shock & vibration (모드 고유치 및 Newmark 동역학)
    │   └── test_v08_nafems_thermal.py   # NAFEMS P18 T1: Circular hot spot membrane (열탄성 연성)
    │
    └── 2_examples/                  # 실무 응용군 (Coreform Example Problems 기반, 복잡 기하/비선형)
        ├── test_e01_scordelis_lo_roof.py # Scordelis-Lo 셸 지붕 휨 해석 (곡면 셸 거동)
        ├── test_e02_free_vibe_cube.py   # 30° 기울어진 3D 단위 큐브 자유 진동 (Mindlin-Lamé 모드)
        ├── test_e03_hyperelastic_block.py# 초탄성(Mooney-Rivlin/Neo-Hookean) 대변형 인장/압축
        ├── test_e04_bimetallic_strip.py # 다중 재질 계면(Bimetallic Interface) 열-구조 해석
        └── test_e05_cad_bracket_e2e.py  # 3D CAD B-Rep 브래킷/너클 엔드투엔드 파이프라인 (VTKHDF + H3D)


4. 오픈소스 CAD 및 메쉬 처리 생태계 통합 계획

잠김 등기하해석(Immersed IGA)의 성패는 CAD 형상의 해석 적합성 확보와 고속의 공간 이산화 엔진에 달려 있습니다. 본 솔버 파이프라인은 Zero C++ 빌드 원칙 하에 파이썬 바인딩 및 고성능 C++/Rust 백엔드를 제공하는 최상위 오픈소스 라이브러리를 유기적으로 통합합니다.

3.1 3D CAD B-Rep 처리 및 기하 커널

CAD 엔진은 볼륨 메쉬를 생성하지 않고, 표면 형상의 위상 무결성을 보장하며 경계면 태깅 및 기하 질의(Geometry Query)를 수행합니다.

* OpenCASCADE (pythonOCC-core):
  * 산업 표준 3차원 형상 포맷(STEP, IGES, BREP)의 로드 및 지오메트리 파싱.
  * ShapeAnalysis 및 ShapeFix_Shape 모듈을 통한 결함 CAD 치유: 미소 틈새(Micro-gap) 결합(Sewing), 자기 교차(Self-intersection) 제거, 면 방향성 정렬(Orientation Healing)을 수행하여 완전한 닫힌 다양체(Watertight 2-manifold) 확보.
  * B-Rep Face ID별 디리클레(변위), 노이만(트랙션/압력), 약형식 계면 조건 속성 태깅.
* trimesh + embree (Intel 고속 레이 트레이싱 엔진):
  * B-Rep 표면으로부터 고밀도 삼각망을 임시 구성하여 부호 거리 함수(Signed Distance Function, SDF) 및 광선 추적(Ray Casting) 가속.
  * 수백만 개 배경 격자 노드의 내부(Interior) 및 외부(Void) 포함 판정(In/Out Test)을 멀티스레드로 고속 처리.

3.2 메쉬 생성 및 스플라인 그리드 생태계

* Gmsh (Python API):
  * OpenCASCADE 커널과 네이티브로 연동하여 CAD Face 경계 곡률을 정밀하게 보존하는 2차원 표면 비정형 삼각망 메쉬 생성.
  * 생성된 표면 메쉬는 Algoim 및 QuESo의 경계 다항식 적분 및 Nitsche 약형식 적분 패치로 직접 투입.
* splinepy (TU Wien):
  * 텐서곱(Tensor-product) 2D/3D B-스플라인 및 NURBS 볼륨 배경 격자 생성.
  * Knot Insertion 기반의 고속 베지에 추출 연산자(\mathbf{C}^e) 도출.
  * 파이썬 환경에서 C++ 코어를 직접 호출하는 고속 메모리 뷰 제공.
* thbsplines:
  * CAD 경계 곡면 근방의 국소 기하 피처 해상도를 선택적으로 높이기 위한 절단 계층형 B-스플라인(Truncated Hierarchical B-splines, THB-splines) 옥트리 관리.
  * 단위의 분할(Partition of Unity) 및 선형 독립성을 만족하는 활성 기저 함수 식별.

3.3 절단 다면체 기하 알고리즘 엔진

* Algoim (Robert Saye, LBNL):
  * 다변수 다항식 암묵적 레벨셋 곡면에 대해 다가지 높이 함수(Height Function)와 재귀적 차원 축소(Dimension Reduction)를 적용.
  * 절단 요소에 대해 기계 정밀도(Machine Precision) 수준의 1차원 가우스-르장드르 적분 가중치 및 위치 산출.
* QuESo (Quadrature for Embedded Solids, TUM / Kratos 에코시스템):
  * 복합 STL 삼각망 및 B-Rep 곡면과 배경 육면체 셀 간의 고속 계층 적분(Octree/Spacetree Sub-cell Integration) 및 모멘트 피팅(Moment Fitting) 적분점 산출.

5. 핵심 수치역학 메커니즘 및 대수적 수식 체계

4.1 베지에 추출 (Bézier Extraction)

임의의 차수 p, q, r과 연속성을 갖는 전역 스플라인 기저 함수 R_i(\boldsymbol{\xi})는 개별 요소 단위의 표준 번스타인 다항식 기저 B_j(\boldsymbol{\xi})와 국소 베지에 추출 행렬 \mathbf{C}^e의 곱으로 완전하게 분해됩니다 (Borden et al., 2011): R_i(\boldsymbol{\xi}) = \sum_{j=1}^{(p+1)(q+1)(r+1)} C_{ij}^e B_j(\boldsymbol{\xi}) 이를 통해 솔버는 대역적 스플라인 위상을 직접 처리할 필요 없이, 표준 유한요소의 요소별 강성 행렬 조립 루틴을 그대로 활용하여 전역 강성 행렬을 조립합니다: \mathbf{K}^e = (\mathbf{C}^e)^T \left( \int_{Q_e \cap \Omega} \mathbf{B}_B^T \mathbf{D} \mathbf{B}_B \, d\Omega \right) \mathbf{C}^e 여기서 \mathbf{B}_B는 번스타인 다항식의 미분 연산자 행렬이며, \mathbf{D}는 탄성 구성방정식 재료 행렬입니다.

4.2 Algoim 고차 절단 셀 적분 메커니즘

배경 요소 Q_e가 물리 영역 경계 \Gamma에 의해 절단될 때, 직선이나 평면 근사로 인한 형상 포화 오차(Geometric Saturation Error)를 방지하기 위해 Saye (2015, 2022)의 알고리즘을 적용합니다. 레벨셋 함수 \phi(\mathbf{x}) = 0에 의해 정의된 활성 체적 K = \{ \mathbf{x} \in Q_e : \phi(\mathbf{x}) \le 0 \}에 대해, 높이 함수 x_n = h(x_1, \dots, x_{n-1})의 다가지 그래프로 국소 변환한 후 1차원 가우스 적분을 재귀적으로 적용합니다: \int_{K} f(\mathbf{x}) \, d\Omega = \int_{H} \left( \int_{a(x')}^{b(x')} f(x', x_n) \, dx_n \right) dx' \approx \sum_{k=1}^{n_{QP}} w_k f(\mathbf{x}_k) 이 적분 체계는 1차원 적분점 수 q에 대해 2q차의 최적 대수적 수렴률을 보장하며, 복잡한 곡면 절단 셀에서도 최적 수렴 속도 O(h^{p+1})를 회복합니다.

4.3 대수적 자유도 집합 안정화 (\mathbf{C}^e \to \mathbf{\tilde{C}}^e 변환)

물리 영역 내부 체적비 \eta_e = \frac{|Q_e \cap \Omega|}{|Q_e|}가 극단적으로 작은 미세 절단 셀(Small Cut-Cell, \eta_e < 10^{-2})은 대역 강성 행렬의 고유값을 0으로 붕괴시켜 조건수를 무제한 증가시킵니다 (\kappa \sim h^{-2p} \cdot \eta^{-1}). 본 아키텍처는 Badia et al. (2018)의 AgFEM 및 Marussig et al. (2018)의 스플라인 확장 이론을 융합하여, 불량 절단 기저 함수(Ill-posed DOFs)를 기하학적 체적이 충분한 인접 내부 요소의 건전한 기저 함수(Well-posed Root DOFs)로 대수적 구속 결합하는 전역 변환 행렬 \mathbf{T}를 구성합니다: \mathbf{R}_{\text{ext}} = \mathbf{T} \, \mathbf{R}_{\text{act}} 이 변환 관계를 요소 단위의 베지에 추출 연산자에 전치 곱하여 대수적으로 수정된 추출 행렬 \mathbf{\tilde{C}}^e를 도출합니다: \mathbf{\tilde{C}}^e = \mathbf{C}^e \mathbf{T}_e 수정된 \mathbf{\tilde{C}}^e가 주입된 요소 강성 행렬은 다음과 같이 평가됩니다: \mathbf{\tilde{K}}^e = (\mathbf{\tilde{C}}^e)^T \mathbf{K}_B^e \mathbf{\tilde{C}}^e 이로써 솔버 내부의 C++ 유한요소 변분식이나 조립 루틴을 전혀 수정하지 않고도, 인위적 페널티 튜닝 없이 대역 행렬의 조건수 상한(\kappa \le C h^{-2p})을 완벽하게 보장합니다.

4.4 약형식 경계 조건 부과 (Nitsche 기법)

절단 경계면 \Gamma_D 상에서 디리클레 경계 조건(\mathbf{u} = \mathbf{g})을 부과하기 위해 변분 일관성을 만족하는 대칭 Nitsche 약형식을 적용합니다: a_{\text{Nitsche}}(\mathbf{u}, \mathbf{v}) = a(\mathbf{u}, \mathbf{v}) - \int_{\Gamma_D} (\boldsymbol{\sigma}(\mathbf{u}) \mathbf{n}) \cdot \mathbf{v} \, d\Gamma - \int_{\Gamma_D} \mathbf{u} \cdot (\boldsymbol{\sigma}(\mathbf{v}) \mathbf{n}) \, d\Gamma + \frac{\gamma}{h_F} \int_{\Gamma_D} \mathbf{u} \cdot \mathbf{v} \, d\Gamma 여기서 \gamma는 안정화 파라미터이며, 고차 미분 연속성(C^{p-1}) 덕분에 표준 유한요소 대비 경계면 근방의 트랙션 응력 진동이 현저히 억제됩니다.

6. 현 아키텍처의 한계점 및 심층 비판적 기술 리뷰 (Critical Review)

성공적인 솔버 완성을 위해, 현재 마스터플랜에서 간과되었거나 수치역학적·소프트웨어 공학적으로 병목이 될 수 있는 6대 핵심 기술적 결함과 위험 요소를 선제적으로 분석하고 공학적 대응책을 수립합니다.

5.1 3D CAD B-Rep의 Trimming Curve 위상 결함 및 수치 공차 불일치

* 문제점 분석: 상용 CAD(B-Rep) 모델은 곡면 간 교차 연산 시 매개변수 곡선(PCurve)과 3차원 공간 곡선 간에 수치 공차(Tolerance, 통상 10^{-3} \sim 10^{-6}) 불일치로 인해 미세 간극(Micro-gap)이나 비다양체 에지(Non-manifold edge)를 내포합니다. 형상 치유(Sewing)가 불완전한 상태에서 절단 요소 판정을 수행하면, 특정 요소에서 체적 누수(Volume Leak)가 발생하여 적분점이 물리 영역 외부로 소실되거나 내부 요소가 외각 보이드로 오판정될 위험이 존재합니다.
* 공학적 대응책: PythonOCC의 ShapeFix_Shape 파이프라인에서 공차 강제 보정(Tolerance Tightening) 루틴을 의무화하고, 광선 추적(Ray casting) 기반의 기수(Parity) 판정 엔진에 3축 다중 광선 투표제(Majority Voting) 알고리즘을 도입하여 위상 결함에 대한 기하학적 강건성을 확보합니다.

5.2 유리(Rational) NURBS와 다항식 적분 엔진(Algoim) 간의 수학적 괴리

* 문제점 분석: CAD의 원통, 구, 원뿔 등 원추 곡면은 분모에 가중치 다항식이 존재하는 유리 다항식(Rational NURBS)으로 표현됩니다. 반면 Algoim의 높이 함수 및 차원 축소 적분 엔진은 순수 다변수 다항식(Polynomial Level-Set)을 전제로 작동합니다. 또한 Antolin et al.의 무적분(Quadrature-Free) 기법 역시 비유리 베지에 곡면/곡선에 한정됩니다. 유리 NURBS를 무리하게 다항식 공간으로 변환할 경우, 테일러 급수 근사 오차가 발생하거나 텐서곱 다항식 차수가 기하급수적으로 팽창(p=4일 때 유효 피적분 다항식 차수 w > 800)하여 연산 과부하를 유발할 수 있습니다.
* 공학적 대응책: 원추 곡면 B-Rep에 대해서는 요소 단위로 적응형 베지에 다항식 투영(Bézier Re-parameterization of degree p)을 수행하여 기하학적 오차를 O(h^{p+1}) 수준으로 제한하거나, QuESo의 다면체 모멘트 피팅(Moment Fitting) 방식을 하이브리드로 채택하여 유리 표현형의 적분 안정성을 확보합니다.

5.3 3D 복잡 형상에서의 AgFEM 셀 집합(Aggregation) 그래프 위상 분기

* 문제점 분석: 2D에서는 절단 셀의 인접 내부 셀을 찾는 그래프 탐색이 단순하지만, 3D 박판 구조(Thin-shell), 날카로운 모서리(Cusp), 복합 리브(Rib) 형상에서는 단일 소형 절단 셀이 다수의 내부 셀과 면 접촉하지 못하고 능선(Edge) 또는 꼭짓점(Vertex)으로만 접촉하는 고립 절단 셀(Isolated Sliver Cell)이 발생합니다. 이 경우 불량 자유도가 엉뚱한 내부 요소로 병합되어 비국소적 제약(Non-local constraint)이 전파되거나, 연립방정식 그래프의 대역폭(Bandwidth)이 과도하게 팽창하여 스텐실 구조가 왜곡될 수 있습니다.
* 공학적 대응책: Badia et al. (2022)의 계면-면 연결성 기반 다단계 집합 알고리즘(Facet-connected Aggregation)을 도입하고, 면 접촉이 불가능한 경우 2차 근접 내부 요소를 체계적으로 탐색하는 위상 안정화 그래프 탐색기를 Rust 코어(aggregation_graph.rs)에 구현합니다.

5.4 Kratos Python 바인딩 수준의 조립(Assembly) 오버헤드 및 CSR 변환 병목

* 문제점 분석: Kratos의 IgaApplication은 대역 스플라인 모델파트를 C++ 코어로 구축하도록 설계되어 있으나, 전처리 단계에서 동적으로 생성된 수십만 개의 가변 \mathbf{\tilde{C}}^e 행렬과 비정형 적분점 세트를 Python-C++ 바인딩을 통해 주입할 경우 심각한 인터프리터 마샬링 오버헤드가 발생합니다.
* 공학적 대응책: 가변 추출 행렬과 적분점 데이터를 중간 텍스트 변환 없이 Rust에서 직접 bezier_grid.mdpa 블록으로 일괄 고속 스트리밍 직렬화하고, Kratos 솔버 구동 시 내장된 고속 바이너리 파서(MPI-ready Custom Importer)를 경유하도록 데이터 파이프라인을 최적화합니다.

5.5 직접 솔버(PyPardiso)의 메모리 한계와 대규모 3D 모델 확장성 부재

* 문제점 분석: p=3, 4의 3차원 고차 스플라인 요소는 기저 함수의 중첩 지지 영역(Support)이 넓어, 1차 유한요소 대비 전역 강성 행렬의 0이 아닌 성분(Non-zero entries)이 요소당 수십 배로 급증합니다. 이로 인해 수백만 자유도 이상의 실무 3D CAD 해석 시 Intel MKL PARDISO의 직접 분해(LU/Cholesky) 과정에서 메모리(RAM) 요구량이 테라바이트 수준으로 치솟아 메모리 부족 충돌이 불가피합니다.
* 공학적 대응책: 1단계(중소규모 검증 및 수치 엄밀성 확인)에서는 PyPardiso의 무조건적 수렴성을 활용하되, 2단계 확장 모듈로서 AgFEM 스케일링 특성을 고려한 Schwarz 영역 분할 전처리자(Subspace Splitting / Additive Schwarz Preconditioner) 기반의 Krylov 부공간 반복 솔버(Kratos C++ Conjugate Gradient / GMRES) 인터페이스를 마스터플랜에 병행 설계합니다.

5.6 Nitsche 경계 조건의 안정화 파라미터 민감성

* 문제점 분석: Nitsche 기법의 페널티 파라미터 \gamma가 임계값보다 작으면 행렬의 양의 정정성(Coercivity)이 상실되어 해가 발산하고, 너무 크면 시스템 행렬의 조건수가 악화되어 해의 품질이 저하됩니다. 절단 비율 \eta_e \to 0인 셀에서 국소 역 부등식(Inverse Inequality) 상수가 요소 기하에 따라 변동하므로 전역 고정 \gamma를 사용하는 것은 위험합니다.
* 공학적 대응책: 절단 셀의 국소 경계 적분점 및 체적 비율에 기반하여 국소 고유치 문제를 풀거나 메쉬 특성 길이를 자동 산정하는 요소별 적응형 Nitsche 파라미터(\gamma_e = C \cdot \text{trace}(\mathbf{K}^e)) 자동 계산 루틴을 구축합니다.

7. 사용자 중심 API 아키텍처 설계: Coreform IGA와 Kratos Multiphysics 패러다임 융합

실제 엔지니어링 실무자가 복잡한 수치 알고리즘의 세부 사항을 몰라도 직관적이고 생산적으로 해석을 수행할 수 있도록, **Coreform IGA의 직관적인 사용자 스크립팅 사용성(Ergonomics)**과 Kratos Multiphysics의 모듈식 모델-프로세스(Model-Process) 프레임워크를 융합한 차세대 Pythonic API를 설계합니다.

6.1 벤치마크 대상 API 아키텍처 비교 분석

* Coreform IGA API 패러다임:
  * coreform.mesh: CAD 파일 임포트, B-스플라인 배경 바운딩 박스 메쉬 생성, 적응형 세분화 제어.
  * coreform.job: 재료 물성 정의, 경계조건 부과, 시뮬레이션 실행 및 모니터링.
  * coreform.probe: 특정 관심 지점 또는 곡면/라인 상의 고정밀 물리량 샘플링 및 XY 데이터 플로팅.
  * 특징: 철저히 사용자 친화적이며 단 몇 줄의 직관적인 체이닝 스크립트로 3D Immersed 해석 완료 가능.
* Kratos Multiphysics API 패러다임:
  * Model & ModelPart: 계층적 노드, 요소, 조건(Condition), 프로퍼티 컨테이너.
  * ProjectParameters.json: 솔버 유형, 시간 적분 스킴, 수렴 기준 등을 JSON으로 엄격히 관리.
  * Process: 전처리, 경계조건 적용, 솔버 루프, 후처리가 독립된 프로세스 파이프라인으로 체계화.
  * 특징: 다중물리(Multiphysics) 확장이 용이하고 객체지향 아키텍처가 견고하나 진입 장벽이 높음.

6.2 융합 하이브리드 Python API 계층 구조 (pyimmersed_iga)

우리 솔버는 최상위 레벨에서 Coreform IGA 스타일의 유려한 파이썬 스크립팅 인터페이스를 제공하고, 백엔드에서는 Kratos의 표준 데이터 구조와 프로세스를 자동으로 제어하는 3계층 아키텍처를 취합니다.

┌────────────────────────────────────────────────────────────────────────┐
│  High-Level User Scripting API (pyimmersed_iga)                        │
│  • Model(name)                                                         │
│  • model.import_cad(file_path, heal=True)                              │
│  • model.create_background_mesh(degree, cell_sizes, bbox_margin)       │
│  • model.add_material(ElasticIsotropic / Hyperelastic)                 │
│  • model.add_boundary_condition(Dirichlet / Traction / WeakNitsche)    │
│  • model.add_probe(PointProbe / LineProbe / SurfaceProbe)              │
│  • model.solve(solver_type="pypardiso")                                │
│  • model.export(vtkhdf=True, h3d=True)                                 │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│  Mid-Level Pipeline Orchestrator & Kratos Bridge                       │
│  • Preprocessing Pipeline: splinepy + Rust Core (Algoim + AgFEM)       │
│  • Kratos ProjectParameters Auto-Synthesizer                           │
│  • Kratos ModelPart & Bézier Element Manager                           │
│  • PyPardiso Direct Sparse Solver Adapter                              │
│  • Postprocessing Dispatcher (VTKHDF / Altair Bridge)                  │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│  Low-Level Execution Kernels                                           │
│  • Rust Engine (rust_iga_core) ── sqlite_storage, aggregation_graph    │
│  • Kratos Native Wheel ── StructuralMechanicsApplication               │
└────────────────────────────────────────────────────────────────────────┘


6.3 사용자 API 클래스 인터페이스 상세 명세

# pyimmersed_iga/api.py
from typing import List, Union, Tuple, Optional
import numpy as np

class ImmersedModel:
    """사용자가 직접 대면하는 최상위 Immersed IGA 시뮬레이션 모델 클래스"""
    def __init__(self, name: str = "Immersed_Analysis"):
        self.name = name
        self.cad_geometry = None
        self.mesh_config = {}
        self.materials = {}
        self.boundary_conditions = []
        self.probes = []
        self.solver_settings = {}
        
    def import_cad(self, cad_file: str, heal: bool = True, tolerance: float = 1e-4) -> "CADGeometry":
        """STEP/IGES 형상을 로드하고 Watertight B-Rep 치유를 수행"""
        self.cad_geometry = CADGeometry(cad_file, heal=heal, tolerance=tolerance)
        return self.cad_geometry

    def create_background_mesh(
        self, 
        degree: Union[int, Tuple[int, int, int]] = 2,
        cell_size: Union[float, Tuple[float, float, float]] = 0.05,
        margin: float = 0.1
    ) -> "BackgroundMesh":
        """CAD 형상을 감싸는 최적 직교 B-스플라인 배경 격자 생성"""
        self.mesh_config = {
            "degree": degree,
            "cell_size": cell_size,
            "margin": margin,
            "stabilization": "agfem_discrete_extension",
            "eta_crit": 1e-2
        }
        return BackgroundMesh(self.cad_geometry, self.mesh_config)

    def add_material(self, name: str, material_type: str = "LinearElastic", **properties) -> None:
        """Kratos ConstitutiveLaw와 연동되는 재료 물성 정의"""
        # properties: young_modulus, poisson_ratio, density 등
        self.materials[name] = {"type": material_type, "properties": properties}

    def assign_material(self, material_name: str, region_id: Optional[int] = None) -> None:
        """형상 전체 또는 특정 영역에 재료 할당"""
        pass

    def add_displacement_bc(self, face_id: int, u: Tuple[Optional[float], Optional[float], Optional[float]]) -> None:
        """CAD Face ID에 디리클레 구속 조건(변위) 부여 (약형식 Nitsche 적용)"""
        self.boundary_conditions.append({
            "type": "Displacement", "face_id": face_id, "value": u, "method": "Nitsche"
        })

    def add_pressure_bc(self, face_id: int, pressure: float) -> None:
        """CAD Face ID에 면압 노이만 하중 조건 부여"""
        self.boundary_conditions.append({
            "type": "Pressure", "face_id": face_id, "value": pressure
        })

    def add_point_probe(self, name: str, coordinates: Tuple[float, float, float], fields: List[str] = ["STRESS", "DISPLACEMENT"]) -> None:
        """B-스플라인 고정밀 역매핑 기반 관심 지점 프로브 생성"""
        self.probes.append({"name": name, "type": "point", "coords": coordinates, "fields": fields})

    def solve(self, solver_type: str = "pypardiso", threads: int = 8) -> "AnalysisResults":
        """Stage 2-A 전처리(Rust) -> Stage 2-B Kratos 솔버 파이프라인 자동 실행"""
        # 1. iga_config.json 및 ProjectParameters.json 자동 합성
        # 2. Rust 모듈 구동: 절단 적분점 생성 및 C^e -> C̃^e 변환 후 mesh_data.sqlite 저장
        # 3. Kratos IgaApplication에 C̃^e 및 하중 주입 후 PyPardiso 해석 수행
        # 4. control_points.hdf5 출력 및 프로브 평가
        return AnalysisResults(self)

    def export_results(self, export_vtkhdf: bool = True, export_h3d: bool = True, submesh_density: int = 4) -> None:
        """Dual Postprocessor(ParaView VTKHDF 및 Altair H3D) 일괄 출력"""
        pass


6.4 엔지니어링 실무 스크립트 작성 시나리오 예시

다음은 엔지니어가 원형 구멍 평판 인장 해석(Kirsch 문제)을 수행할 때 작성하게 되는 최종 사용자 스크립트의 완전한 형태입니다. Coreform IGA의 간결함과 Kratos의 고성능 엔진이 완벽히 융합된 형태를 보여줍니다:

# run_plate_with_hole.py
from pyimmersed_iga import ImmersedModel

# 1. 모델 인스턴스 초기화 및 CAD 형상 로드
model = ImmersedModel("Plate_With_Hole_Kirsch")
cad = model.import_cad("plate_with_hole.step", heal=True)

# 2. 직교 B-스플라인 배경 격자 정의 (2차 스플라인, 격자 크기 0.05)
model.create_background_mesh(degree=(2, 2, 2), cell_size=(0.05, 0.05, 0.05), margin=0.05)

# 3. 등방성 탄성 재료 물성 정의 (강철)
model.add_material("Structural_Steel", material_type="LinearElastic", 
                   young_modulus=2.1e11, poisson_ratio=0.3, density=7850.0)
model.assign_material("Structural_Steel")

# 4. 경계 조건 부여 (대칭면 구속 및 원거리 인장 하중)
# Face 1: X-대칭면 (ux = 0), Face 2: Y-대칭면 (uy = 0)
model.add_displacement_bc(face_id=1, u=(0.0, None, None))
model.add_displacement_bc(face_id=2, u=(None, 0.0, None))
# Face 3: 우측 외곽 인장 하중 10 MPa (Traction X = 1.0e7)
model.add_pressure_bc(face_id=3, pressure=-1.0e7)

# 5. 응력 집중 관심 지점 프로브 설정 (구멍 정점: Theoretical Stress Concentration Kt = 3.0)
model.add_point_probe("Hole_Crown_Point", coordinates=(0.0, 1.0, 0.0), fields=["STRESS", "DISPLACEMENT"])

# 6. 해석 실행 (Intel MKL PARDISO 16스레드 직접 솔버)
results = model.solve(solver_type="pypardiso", threads=16)

# 7. 프로브 결과 확인 및 엔터프라이즈 듀얼 후처리 파일 저장
crown_stress = results.get_probe_data("Hole_Crown_Point")
print(f"Computed Peak Sigma_XX: {crown_stress['VON_MISES'] / 1.0e6:.2f} MPa (Target: 30.00 MPa)")

# ParaView용 .vtkhdf 및 HyperView용 .h3d 자동 생성
model.export_results(export_vtkhdf=True, export_h3d=True, submesh_density=3)


8. Coreform IGA 기반 이원화 테스트 스위트: Verification vs Example

테스트 스위트는 Coreform IGA의 검증 체계를 철저히 벤치마킹하여, 수치적 엄밀성과 수렴률을 정량 검증하는 **1_verification (검증군)**과 실제 엔터프라이즈 CAD 및 다중물리 해석 역량을 실증하는 **2_examples (실무 예제군)**으로 엄격히 이원화하여 구축합니다.

tests/
├── 1_verification/                  # [수치적 정밀도 및 해석해 일치성 엄밀 검증]
│   ├── test_v01_patch_test.py       # 단위 큐브 고차 변위 패치 테스트 (기계 정밀도 검증)
│   ├── test_v02_kirsch_plate.py     # Infinite plate with a hole (Kirsch 해석해, Kt=3.0)
│   ├── test_v03_roark_hole.py       # Roark 17.1-7a: 유한 평판 원형 홀 응력 집중
│   ├── test_v04_timoshenko_beam.py  # 경사 외팔보 휨 (Timoshenko 보 해석해)
│   ├── test_v05_thick_cylinder.py   # 내압을 받는 3D 후육 원통/구형 압력용기 (Lamé 해석해)
│   ├── test_v06_sliver_cut_cond.py  # 극단 미세 절단 셀(η -> 10^-8) 조건수 바운딩 불변성
│   ├── test_v07_beam_dynamics.py    # 외팔보 충격/진동 및 모드 고유치 해석
│   └── test_v08_nafems_thermal.py   # NAFEMS P18 T1: 원형 고온부 막 열탄성 연성 해석
│
└── 2_examples/                      # [실무 형상, 복합 기하 및 비선형/다중물리 응용]
    ├── test_e01_scordelis_lo_roof.py # Scordelis-Lo 셸 지붕 휨 (곡면 셸 거동 검증)
    ├── test_e02_free_vibe_cube.py   # 30° 기울어진 3D 단위 큐브 자유 진동 (Mindlin-Lamé 모드)
    ├── test_e03_hyperelastic_block.py# 초탄성(Mooney-Rivlin/Neo-Hookean) 대변형 인장/압축
    ├── test_e04_bimetallic_strip.py # 다중 재질 계면(Bimetallic Interface) 열-구조 변형
    └── test_e05_cad_bracket_e2e.py  # 산업용 3D CAD B-Rep 브래킷 엔드투엔드 파이프라인


7.1 Tier 1: 1_verification 상세 사양 및 정량 통과 기준

테스트 ID 및 명칭	대상 기하 및 조건	주요 수치 검증 항목	합격 기준 (Acceptance Criteria)
V-01: Prescribed Displacement Cube	임의로 회전된 단위 큐브에 고차 변위장 인가	기저 함수 일관성(Consistency) 및 완전성	|oldsymbol{\sigma}_h - oldsymbol{\sigma}_{\text{exact}}|_{L^\infty} < 10^{-12}
V-02: Infinite Plate with Hole	원형 구멍 1/4 대칭 평판 (Kirsch 문제)	Algoim 곡면 적분 수렴성, K_t=3.0 복원	피크 응력 오차 0.5% 미만, L^2 수렴률 O(h^{p+1})
V-03: Roark 17.1-7a Finite Plate	유한 폭 평판 내 원형 홀 인장 하중	Roark 공식집 응력 집중 계수 정밀도 비교	이론값 대비 오차 1.0% 이내
V-04: Timoshenko Beam Bending	배경 격자와 비정합하는 경사 외팔보 휨	전단 휨 거동 및 Nitsche 약형식 구속 오차	끝단 처짐량 이론해 대비 1.0% 이내
V-05: Thick-Walled Vessel	내압 p를 받는 3D 후육 구형/원통 압력용기	완전한 3D 곡면 절단 셀 적분 및 Lamé 해석해	반경/후프 응력 분포 오차 0.8% 미만
V-06: Extreme Sliver Cut Test	구멍 중심을 미세 이동시켜 \eta_e \to 10^{-8} 유발	\mathbf{C}^e \to \mathbf{\tilde{C}}^e 대수적 안정화 조건수 제어	조건수 \kappa(\mathbf{K}) \le C h^{-2p} (상한 불변)
V-07: Beam Shock & Dynamics	외팔보 모드 해석 및 암시적 동역학 시간 적분	일관 질량 행렬(\mathbf{M}) 고유치 및 과도 응답	1~5차 고유진동수 해석해 대비 0.5% 이내
V-08: NAFEMS P18 T1 Thermal	원형 핫스팟이 존재하는 박막 열탄성 해석	정상 열전도 및 열변형률 연성 해석 정밀도	NAFEMS 표준 벤치마크 기준치 대비 1% 이내

7.2 Tier 2: 2_examples 상세 사양 및 실무 검증 기준

테스트 ID 및 명칭	대상 기하 및 엔지니어링 문제	핵심 기술 실증 항목	통과 기준 (Evaluation Criteria)
E-01: Scordelis-Lo Roof	원통형 셸 지붕의 자중 휨 해석	곡면 막-휨(Membrane-Bending) 잠김 방지	중앙부 수직 처짐량 표준 해(0.3024) 대비 1% 이내
E-02: Free Vibrating Tilted Cube	배경 격자를 30° 대각 절단하는 3차원 큐브	복합 3차원 절단 셀 상의 질량 행렬 정밀도	Mindlin-Lamé 고유진동수 스펙트럼 완벽 일치
E-03: Hyperelastic Block	원통형 보이드를 포함한 고무 블록 대변형	Neo-Hookean/Mooney-Rivlin 비선형 수렴성	Newton-Raphson 2차 수렴성 유지, 절단 셀 안정
E-04: Bimetallic Strip	열팽창 계수가 상이한 이종 금속 접합 보	다중 재질 계면(Multi-material Interface) 결합	열하중 하 곡률 및 계면 전단 응력 이론해 일치
E-05: 3D CAD Bracket E2E	필렛, 볼트 홀이 포함된 복합 B-Rep 브래킷	PythonOCC 치유 \to PyPardiso \to Dual Post	상용 고밀도 사면체 메쉬 대비 에너지 오차 1% 이내

9. 단계별 연구개발 마일스톤 및 세부 실행 계획

전체 R&D 공정은 수치 수학적 최소 기능 제품(MVP) 관통부터 시작하여 산업용 3D CAD 대응 단계까지 4단계로 점진적으로 전개됩니다.

Phase 1: 2D 프로토타입 완전 관통 및 수치 엄밀성 입증 (주초 착수)

* 목표: 3D 기하의 복잡성을 배제하고, 2D 텐서곱 B-스플라인 및 원형 구멍 평판(Kirsch) 문제에서 파이프라인 전체를 끝단까지(End-to-End) 연결.
* 핵심 마일스톤:
  1. python_pipeline/s2_preprocessor/db_manager.py 완성: mesh_data.sqlite 테이블 생성 및 바이너리 패킹/언패킹 I/O 구현.
  2. python_pipeline/s2_preprocessor/bezier_stabilizer.py 프로토타이핑: 2D B-스플라인에서 \eta_e < 10^{-2}인 소형 절단 기저 함수를 건전한 인접 기저로 병합하는 변환 행렬 \mathbf{T} 및 \mathbf{\tilde{C}}^e = \mathbf{C}^e \mathbf{T}_e 알고리즘 작성.
  3. tests/1_verification/test_v02_kirsch_plate.py 구축: 미안정화 시 조건수 폭발과 안정화 적용 시 조건수 바운딩(\kappa \approx 10^4) 및 K_t=3.0 응력 복원 정밀도 검증.
  4. python_pipeline/s3_kratos_solver/pardiso_solver.py 작성: Kratos IgaApplication CSR 행렬을 가로채 PyPardiso로 초고속 직접 해석 수행.

Phase 2: Rust 고속 전처리 코어 및 3D 기하 절단 엔진 통합

* 목표: 대규모 3D 배경 격자 처리 속도 병목을 해소하기 위한 Rust 확장 모듈 완성 및 Algoim 고차 적분 연동.
* 핵심 마일스톤:
  1. rust_iga_core 패키지 빌드: Rayon 기반 병렬 AABB/BVH 트리 공간 탐색기(spatial_culling.rs) 및 rusqlite 고속 BLOB 스트리밍 모듈(sqlite_storage.rs) 구축.
  2. 3D 절단 셀에 대한 Algoim C-API 바인딩 연동(algoim_bridge.rs) 및 기계 정밀도 적분점 산출.
  3. 3D 복합 위상 그래프 병합 알고리즘(aggregation_graph.rs) 구현 및 V-05(구형 압력용기), V-06(슬라이버 컷 조건수 불변성) 통과.

Phase 3: 사용자 친화적 통합 API (pyimmersed_iga) 및 듀얼 후처리 파이프라인

* 목표: Coreform IGA 스타일의 고수준 사용자 API 완성 및 ParaView/HyperView 듀얼 익스포터 가동.
* 핵심 마일스톤:
  1. pyimmersed_iga 최상위 클래스(ImmersedModel, CADGeometry, BackgroundMesh) 구현 및 스크립팅 인터페이스 완성.
  2. python_pipeline/s4_dual_postprocessor/vtkhdf_writer.py 작성: h5py 기반 시간-공간 청크 데이터셋 직렬화로 ParaView 시계열 시각화 지원.
  3. python_pipeline/s4_dual_postprocessor/altair_bridge.py 완성: Rust ryu 기반 고속 .fem + .hvascii 작성 및 Altair hvtrans CLI 백그라운드 호출을 통한 .h3d 파일 자동 생성.
  4. B-스플라인 점 역매핑 프로브 엔진(probe_results.json) 구축.

Phase 4: Coreform급 실무 3D CAD B-Rep 엔드투엔드 벤치마크 및 다중물리 확장

* 목표: 복잡한 산업용 3D CAD 부품에 대한 완전 자동화 해석 입증 및 비선형/진동 벤치마크 스위트 완성.
* 핵심 마일스톤:
  1. PythonOCC 형상 치유(Sewing/Healing) 파이프라인과 Gmsh 표면 이산화 모듈의 완전 자동화 연동.
  2. Coreform Verification Manual (V-01 ~ V-08) 및 Example Problems (E-01 ~ E-05) 전체 자동화 테스트 리그 실행 및 100% 회귀 통과.
  3. Kratos StructuralMechanicsApplication 기반 대변형 초탄성(E-03) 및 구조 동역학 모드 해석(V-07) 기능 배포.

10. 종합 참고 문헌 및 오픈소스 소프트웨어 레지스트리 (References & Registry)

9.1 핵심 학술 문헌 (Academic References)

1. Borden, M. J., Scott, M. A., Evans, J. A., & Hughes, T. J. R. (2011).
Isogeometric finite element data structures based on Bézier extraction of NURBS.
International Journal for Numerical Methods in Engineering, 87(1–5), 15–47.
DOI: 10.1002/nme.2968
2. Badia, S., Verdugo, F., & Martín, A. F. (2018).
The aggregated unfitted finite element method for elliptic problems.
Computer Methods in Applied Mechanics and Engineering, 336, 533–553.
DOI: 10.1016/j.cma.2018.03.022
3. Badia, S., Neiva, E., & Verdugo, F. (2022).
Linking ghost penalty and aggregated unfitted methods.
Computer Methods in Applied Mechanics and Engineering, 388, 114232.
DOI: 10.1016/j.cma.2021.114232
4. Badia, S., Neiva, E., & Verdugo, F. (2022).
Robust high-order unfitted finite elements by interpolation-based discrete extension.
Computers & Mathematics with Applications, 127, 105–126.
DOI: 10.1016/j.camwa.2022.09.027
5. Saye, R. I. (2015).
High-order quadrature methods for implicitly defined surfaces and volumes in hyperrectangles.
SIAM Journal on Scientific Computing, 37(2), A993–A1019.
DOI: 10.1137/140966290
6. Saye, R. I. (2022).
High-order quadrature on multi-component domains implicitly defined by multivariate polynomials.
Journal of Computational Physics, 448, 110720.
DOI: 10.1016/j.jcp.2021.110720
7. Antolin, P., Buffa, A., & Hirschler, M. (2022).
Quadrature-free immersed isogeometric analysis.
Engineering with Computers, 38, 4475–4499.
DOI: 10.1007/s00366-021-01439-0 (arXiv: 2107.09024)
8. Burman, E. (2010).
Ghost penalty.
Comptes Rendus Mathematique, 348(21–22), 1217–1220.
DOI: 10.1016/j.crma.2010.10.027
9. Burman, E., Claus, S., Hansbo, P., Larson, M. G., & Massing, A. (2015).
CutFEM: Discretizing geometry and partial differential equations.
International Journal for Numerical Methods in Engineering, 104(7), 472–501.
DOI: 10.1002/nme.4823
10. Burman, E., Hansbo, P., & Larson, M. G. (2023).
On the design of locking free ghost penalty stabilization and the relation to CutFEM with discrete extension.
Computational Mechanics, 71, 1007–1019.
DOI: 10.1007/s00466-022-02263-4
11. Schillinger, D., Dedè, L., Scott, M. A., Evans, J. A., Borden, M. J., Rank, E., & Hughes, T. J. R. (2012).
An isogeometric design-through-analysis methodology based on adaptive hierarchical refinement of NURBS, immersed boundary methods, and T-spline CAD surfaces.
Computer Methods in Applied Mechanics and Engineering, 249–252, 116–150.
DOI: 10.1016/j.cma.2012.03.017
12. Schillinger, D., & Ruess, M. (2015).
The Finite Cell Method: A review in the context of higher-order structural analysis of CAD and image-based geometric models.
Archives of Computational Methods in Engineering, 22(3), 391–455.
DOI: 10.1007/s11831-014-9115-y
13. Marussig, B., & Hughes, T. J. R. (2018).
A review of trimming in isogeometric analysis: Challenges, data exchange and simulation aspects.
Archives of Computational Methods in Engineering, 25(4), 1059–1127.
DOI: 10.1007/s11831-017-9220-9
14. Meßmer, M., Teschemacher, T., Leidinger, L. F., Wüchner, R., & Bletzinger, K.-U. (2022).
Efficient CAD-integrated isogeometric analysis of trimmed solids.
Computer Methods in Applied Mechanics and Engineering, 400, 115584.
DOI: 10.1016/j.cma.2022.115584
15. Divi, S. C. (2021).
Scan-based immersed isogeometric analysis.
Ph.D. Thesis, Eindhoven University of Technology & University of Pavia.
ISBN: 978-90-386-5469-0.
16. Thomas, D. C., Scott, M. A., Evans, J. A., Tew, K., & Benzley, S. E. (2022).
U-splines: splines over unstructured meshes.
Coreform Technical Monograph / arXiv preprint.
17. de Prenter, F., Verhoosel, C. V., van Brummelen, H., & Badia, S. (2023).
Stability and conditioning of immersed finite element methods: Analysis and remedies.
Archives of Computational Methods in Engineering, 30(6), 3617–3656.
DOI: 10.1007/s11831-023-09913-0
18. Badia, S., Dilip, H., & Verdugo, F. (2023).
Space-time AgFEM.
arXiv preprint arXiv:2301.09355.
19. Ruiz-Baier, R., et al. (2026).
Divergence-free unfitted finite element discretisations for the Darcy problem.
Preprint / Zenodo DOI: 10.5281/zenodo.19199177

9.2 오픈소스 소프트웨어 리포지토리 및 기술 레지스트리

소프트웨어 명칭	개발 기관 / 저자	라이선스	기술 개요 및 모듈 위치
KratosMultiphysics	CIMNE / 글로벌 컨소시엄	BSD-4-Clause	IgaApplication, StructuralMechanicsApplication, Core Kratos Framework
Algoim	Robert Saye (LBNL)	BSD-3-Clause	고차 다항식 레벨셋 및 암묵적 기하 1D 차원축소 적분 라이브러리 (C++17)
QuESo	Manuel Meßmer (TUM)	BSD-4-Clause	Kratos 연동 임베디드 솔리드 계층 적분 및 모멘트 피팅 엔진 (C++/Python)
splinepy	TU Wien	MIT	텐서곱 B-스플라인, NURBS 및 고속 베지에 추출 연산자 추출 라이브러리
PythonOCC (OpenCASCADE)	Thomas Paviot	LGPL-3.0	OpenCASCADE 3D CAD B-Rep 지오메트리 파싱, Shape Healing 및 토폴로지 검증
Gmsh	Christophe Geuzaine, J.-F. Remacle	GPL-2.0+	CAD B-Rep 연동 고품질 2D/3D 비정형 표면 및 체적 메쉬 생성 엔진 (Python API)
PyPardiso	Adrian Haas / Intel MKL	BSD-3-Clause	Intel oneMKL PARDISO 멀티스레드 대규모 희소 직접 연립방정식 솔버 파이썬 래퍼
Gridap.jl / GridapEmbedded	Santiago Badia, Francesc Verdugo	MIT	Julia 기반 비적합 유한요소 및 집합 유한요소(AgFEM) 오픈소스 프레임워크
ngsxfem	Christoph Lehrenfeld et al.	LGPL-3.0	NGSolve 기반 고차 CutFEM, FCM 및 Ghost Penalty 연동 오픈소스 라이브러리
Nutils	Evalf group	MIT	파이썬 기반 B-스플라인/NURBS 및 파라미터 트리밍 잠김 해석 프레임워크
QUGaR	Pablo Antolin	BSD-3-Clause	비적합 절단 고차 요소를 위한 고차 수치 적분 C++ 라이브러리
QuaHOG	David Gunderman	BSD-3-Clause	절단 2D/3D 기하의 스토크스 정리 기반 무격자 고차 적분 알고리즘
Coreform IGA Documentation	Coreform LLC	상용/문서공개	Coreform IGA 공식 Example Problems 및 Verification Manual 표준 문서군
