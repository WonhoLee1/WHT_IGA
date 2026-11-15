## 핵심 오픈소스 생태계 및 기술 스택 매핑표

| 파이프라인 계층 | 오픈소스 프로젝트 | 주력 언어 | 라이선스 | GitHub / 공식 URL | 프로젝트 내 핵심 역할 및 활용 포인트 |
| :--- | :--- | :---: | :---: | :--- | :--- |
| **솔버 코어** | **Kratos Multiphysics** (`IgaApplication`, `StructuralMechanics`) | C++ / Python | BSD-3 | https://github.com | • Python 배포판(Wheel) 기반 C++ 빌드 없는 솔버 코어<br>• 베지에 추출 연산자(`C̃^e`) 및 커스텀 적분점 수용<br>• `ConstitutiveLaw`(소성·초탄성·열응력) 및 비선형·접촉 엔진 재활용 |
| **선형 솔버** | **PyPardiso** | Python / C (MKL) | BSD-3 | https://github.com | • Intel MKL PARDISO 기반 초고속 멀티스레드 희소 직접 해법<br>• 고차 스플라인 및 미세 절단 요소로 인한 행렬 조건수 악화 돌파 |
| **스플라인 연산** | **splinepy** | C++ / Python | MIT | https://github.com | • B-스플라인, NURBS, 베지에 멀티패치 고속 연산 및 역매핑(Point Inversion)<br>• 초기 베지에 추출 연산자(`C^e`) 산출 및 후처리 물리장 평가 |
| **스플라인 참조** | **G+Smo** (Geometry + Simulation Modules) | C++ | MPL-2.0 | https://github.com | • 계층적 THB-스플라인(`gsHSplines`) 및 비정형 스플라인 C++ 알고리즘 레퍼런스 |
| **절단 셀 적분** | **Algoim** | C++17 (Header) | BSD / MIT | https://github.io | • 암시적 레벨셋 및 다중 다항식 영역의 고차 가우스 적분점 기계 정밀도 생성<br>• 절단 체적 적분점 및 Nitsche 경계면 표면 적분점 산출 |
| **비적합 적분 래퍼** | **QUGaR** | C++ / Python | MIT | https://github.com | • Algoim을 내장한 비적합 격자(Unfitted/Immersed) 전용 고차 적분 및 재매개변수화 파이썬 인터페이스 |
| **CutFEM 안정화** | **ngsxfem** | C++ / Python | LGPL-2.1 | https://github.com | • 고차 레벨셋 매핑(`Θ_h`-mapping), 니체(Nitsche) 약형식 및 유령 페널티(Ghost Penalty) 구현 레퍼런스 |
| **요소 병합 안정화** | **GridapEmbedded.jl** | Julia | MIT | https://github.com | • 미세 절단 셀 자유도를 인접 안정 셀로 병합하는 AgFEM / Extended B-spline(`C^e → C̃^e`) 알고리즘 레퍼런스 |
| **CAD B-Rep 커널** | **PythonOCC** (OpenCASCADE) | C++ / Python | LGPL-2.1 | https://github.com | • 3D STEP/IGES 파일 로드, Dirty CAD 형상 치유(Shape Healing) 및 면(Face)별 경계조건 태깅 |
| **네이티브 가속** | **PyO3 + maturin** (Rust) | Rust / Python | Apache-2.0 / MIT | https://github.com | • 전·후처리 병목 구간(옥트리 분할, `C̃^e` 병합 그래프, 대용량 `.hvascii` 출력)의 Rust 네이티브 병렬(`Rayon`) 바인딩 |
| **JIT 가속기** | **Numba** | Python | BSD-2 | https://github.com | • C++ 컴파일 없이 NumPy 배열 기반 수치 적분·투영 루프를 `@njit`로 C++ 수준 고속화 |
| **중간 데이터 버스** | **SQLite** (`rusqlite` / `sqlite3`) | C / Rust / Python | Public Domain | https://sqlite.org | • 가변 크기의 `C̃^e` 행렬 및 요소별 적분점을 `BLOB` 바이너리로 저장해 전·후처리 및 솔버 간 제로 카피 수준 데이터 교환 |
| **후처리 (오픈소스)** | **h5py / VTK** (`VTKHDF`) | Python / C++ | BSD-3 | https://h5py.org | • 시계열 비선형·동역학 평가 격자 결과를 단일 HDF5 바이너리(`cad_evaluated.vtkhdf`)로 압축 저장 (ParaView 연동) |
| **후처리 (상용 CAE)** | **Altair `hvtrans` CLI** | 상용 유틸리티 | Proprietary | Altair HyperWorks 내장 | • 평가 격자(`.fem`/`.bdf`)와 절점 응력(`.hvascii`)을 압축 바이너리(`.h3d`)로 백그라운드 변환 (HyperView 연동) |
