# WHT_IGA — Kratos 기반 잠김 등기하해석(Immersed IGA) 솔버 마스터플랜 v2

- 문서 버전: v2.0 (2026-11-14)
- 작업 스트림 / 브랜치: `WHTIGA`
- 코드 위치: 저장소 루트 `WHT_IGA/`
- 선행 문서: `master-plan_v1.md` (보존, 수정하지 않음)
- 근거 자료: `나의조사.md`, `문헌종합.md`, `관련opensourceprj.md`, `papers/P01–P30`, `papers/implementation_synthesis.md`, `user-guide/` (Coreform IGA 공식 문서 크롤 111쪽 + 검증 CSV), Kratos 소스 `D:\PythonCodeStudy\kratos-src` (v10.4.3)
- 상태: **계획 단계. 이 문서가 확정될 때까지 솔버 구현은 시작하지 않는다.**

## 목차

0. 개정 이력 (v1 → v2)
1. 비전과 범위
2. 설계 원칙
3. 시스템 아키텍처
4. 데이터 버스
5. Kratos 연동 설계
6. 디렉토리·모듈 구조
7. 알고리즘과 파라미터
8. 사용자 Python API
9. 단계별 로드맵 (게이트·종료 기준)
10. 검증·시험 계획
11. 위험 등록부
12. 환경과 재현성
13. 참고문헌·라이선스 레지스트리
14. 열린 결정 사항
- 부록 A. Kratos API 근거 (파일:라인)
- 부록 B. 용어

---

## 0. 개정 이력 (v1 → v2)

### 0.1 v2를 쓴 이유

v1의 방향(직교 B-스플라인 배경 격자 + 절단 적분 + 소형 셀 안정화 + Kratos)은 맞았다. 그러나 세 가지를 확인하지 않은 채 쓰였다.

1. **목표 제품이 실제로 무엇을 하는가.** v1은 일반 문헌(THB + Algoim + AgFEM + Nitsche)에서 Coreform 방식을 추정했다.
2. **Kratos가 Python에서 실제로 무엇을 노출하는가.** v1은 존재하지 않는 데이터 경로(`.mdpa`의 `EXTRACTION_OPERATOR` NodalData)를 설계했다.
3. **개발 환경.** Kratos 설치 형태, Python 버전, Rust 설치 여부.

v2는 세 가지를 모두 조사한 결과로 다시 썼다.

### 0.2 확인된 사실

| 항목 | v1 가정 | 확인된 사실 | 출처 |
|---|---|---|---|
| Coreform IGA의 정체 | 독립 솔버 | **Coreform IGA for Abaqus (CIFA)**. Abaqus/Standard의 UEL/UELMAT 플러그인이다. 행렬 풀이, 시간적분, 접촉은 Abaqus가 맡는다. 비공개 `NativeTarget` 솔버는 이름만 보인다 | `user-guide/user-guide/index.md`, `capabilities-and-limitations.md`, `python-api/coreform/job` |
| 기저 | THB-spline | 직교 배경 격자 위의 **trimmed U-spline**. p=1–3, C⁰…C^{p−1}. 국소 세분화는 "hierarchical spline basis", 2:1 균형 | `verification-manual/.../pressurized-sphere.md`, `example-problems/plate-with-hole-adaptivity` |
| 절단 적분 | Algoim | CAD를 **삼각형으로 분할**(`tessellation_size`)해서 절단한다. 경계 적분도 그 삼각형 위에서 한다. 체적 적분 알고리즘은 공개되지 않았다 | `python-api/coreform/mesh/index.md` |
| 소형 셀 처리 | AgFEM, `eta_crit=1e-2` | **small-cell extension**(extension operator + function splitting). `small_cell_volume_ratio` 기본값 **0.2**, 허용 범위 [0, 0.5) | `python-api/coreform/mesh/index.md` l.207, 252–258 |
| 필수 경계조건 | Nitsche | **penalty**(`bc_penalty_value`, 검증에서는 1e10). CAD 면에 삼각형 표면층(SFM3D3)을 만들고 그 절점을 스플라인 장에 penalty로 결합한다 | `release-notes-2026-3.md` l.42, `infinite-plate-with-hole.md` l.75 |
| 데이터 버스 | SQLite (제안) | Coreform도 **SQLite**를 쓴다: `<job>.sql`, `_iga_results.sql`, `_iga_quadrature.sql` | `workflow-overview.md` |
| Kratos | 공식 wheel, C++ 재컴파일 없음 | 공식 pip wheel에는 **IgaApplication이 없다**. 이 PC에는 **Kratos 10.4.3 소스빌드**가 있다(`E:\kratos-build\install`, IgaApplication 포함, Python 3.12, MSVC 1941) | 환경 조사 |
| Python→Kratos 기하 주입 | `.mdpa`로 Cᵉ 전달 | Python에서 `QuadraturePointGeometry`를 **직접 만들 수 없다**. 대신 `NurbsVolumeGeometry.CreateQuadraturePointGeometries(결과, 2, [[u,v,w,가중치],…])`로 **셀 하나의 임의 적분점 묶음**을 기하 1개로 만들 수 있다 | 부록 A-1, A-2 |
| 다점 구속 | 고려 안 함 | `LinearMasterSlaveConstraint`를 Python에서 만들 수 있다. 고유치 전략도 구속을 처리한다 | 부록 A-4, A-5 |
| 가속 | Rust PyO3 | **Rust 미설치**. MSVC 14.41 + CMake/Ninja는 있다 | 환경 조사 |

### 0.3 v1 오류 정정 목록

| # | v1 내용 | 문제 | v2 처리 | 위치 |
|---|---|---|---|---|
| E1 | "Zero C++ Build Overhead" | `나의조사.md` 자체가 C++ BezierGeometry/NitscheCondition을 가정했고, 공식 wheel 전제도 사실과 다르다 | "Python-first, C++는 승격 기준을 넘을 때만"으로 재정의. Phase 0 spike로 Python 경로부터 검증 | §2, §5 |
| E2 | `.mdpa`에 `Element3D8N` + `EXTRACTION_OPERATOR_MCSR` NodalData | `Element3D8N`은 물리가 없는 8절점 placeholder다. 추출 연산자는 요소가 소유하며, 등록되지 않은 변수이고, 어떤 기본 요소도 읽지 않는다 | `.mdpa` 경로 폐기. Python에서 `NurbsVolumeGeometry` → `CreateQuadraturePointGeometries` → StructuralMechanics 요소로 구축 | §5 |
| E3 | 조건수 κ ≤ C·h^{−2p} | 2계 타원형 문제에서 안정화 후 κ = O(h⁻²)이고 p는 상수에만 영향을 준다. 미안정화일 때는 η에 따라 폭발한다 | V-COND 합격 기준: η ∈ [1e-8, 1]에서 κ가 평탄, h에 대한 기울기 −2 | §10 |
| E4 | Algoim은 유리 NURBS를 다룰 수 없다 | Algoim general 엔진은 임의의 매끄러운 level set을 처리한다. 진짜 문제는 parametric B-Rep을 implicit로 바꾸는 단계다 | 주 경로는 tessellation 기반(Coreform과 동등). Algoim은 Phase 6 선택 과제 | §7.4 |
| E5 | Algoim이 2q차 수렴을 "보장" | q점 Gauss는 2q−1차 다항식까지 정확할 뿐이다. h-수렴 2q차는 조건부 관측값이다. O(h^{p+1})는 적분 특성이 아니라 이산화 L² 수렴률이다 | 표현을 고치고, 적분 정확도 단위시험으로 대체 | §7.4, §10.3 |
| E6 | `eta_crit = 1e-2` | Coreform 기본값 0.2보다 20배 느슨하다 | 기본 0.2, 범위 [0, 0.5), sweep 가능 | §7.5 |
| E7 | 소형 셀을 요소 단위 η_e로 판정 | C^{p−1} 스플라인 기저는 (p+1)^d 셀에 걸친다. 문제가 되는 것은 **기저함수의 물리 지지영역**이 작은 경우다 | 기저 지지영역 기준 extended B-spline(Höllig WEB, Marussig 2018) | §7.5 |
| E8 | Nitsche가 기본 | Coreform 기본은 penalty다 | 기본 = penalty, Nitsche = 선택 | §7.6 |
| E9 | Nitsche γ_e = C·trace(Kᵉ) | 임시방편이다 | β = C·E·p²/h 또는 국소 고유치(Kratos `NitscheStabilizationModelPartProcess`). h는 extension root 셀 크기 | §7.6 |
| E10 | QuESo, Algoim, 자체 Rust octree를 동시에 계획 | 중복이고, QuESo는 이 PC에 없다 | 자체 Numba octree + 모멘트 피팅 하나로 시작 | §7.4 |
| E11 | 라이선스 표 | pythonocc는 LGPL-3.0이다. QuESo는 BSD-4(광고 조항), Gmsh는 GPL이다 | OCP(cadquery-ocp) 사용, Gmsh는 외부 프로세스로만 | §13 |
| E12 | 누락 | QP 수 가중 분할, 비음수 모멘트 피팅 가중치, 혼합 U/P 정적 응축, FCM α 기준선, probe = 가중치 0 적분점, 기하 오차와 근사 오차 분리 | 모두 추가 | §7, §10 |
| E13 | Altair H3D가 주 출력 트랙 | 상용 도구에 의존한다 | 선택 사항으로 격하. 주 출력은 VTKHDF/VTU + results SQLite + probe JSON | §3, §7.12 |
| E14 | Rust 코어가 Phase 2 핵심 | Rust 미설치, 사용자 결정은 Numba 우선 | Numba 우선. 전환 대상은 C++ 기본(§2.2) | §2, §9 |
| E15 | splinepy, thbsplines 의존 | 설치되어 있지 않다. 균일 직교 격자라 직접 구현이 단순하다 | Cox–de Boor, Bézier 추출, THB를 자체 Numba로 구현 | §6, §7.2 |
| E16 | Newton 점 역매핑 | 불필요하게 복잡하다 | 항등 사상 트릭: knot을 물리 좌표로 두고 제어점을 Greville 점에 놓으면 x(ξ)=ξ, detJ=1. 역매핑이 자명해진다 | §5.2 |
| E17 | 독립 검증 수단 없음 | Kratos 경로를 검증할 기준이 없다 | Numba 참조 조립기 `wht_iga.reference`를 oracle 겸 fallback으로 둔다 | §5.5 |
| E18 | 자체 구성 V01–V08, E01–E05(Scordelis-Lo, 기울어진 큐브 등) | Coreform과 정량 비교가 불가능하다 | Coreform 검증 매뉴얼 12건 + 예제 5건에 1:1 매핑, CSV 기준값 사용 | §10 |

### 0.4 사용자 결정 사항 (2026-11-14, 구속력 있음)

1. Coreform IGA와 Abaqus는 **없다고 가정**한다. 제품은 **Kratos + Python** 조합만으로 구성한다. 자체 Python API를 제공하고, Abaqus `.inp` 호환은 요구하지 않는다.
2. **Kratos Python 패키지와 바인딩**을 기준으로 개발한다. C++ 직접 구현은 이후에 검토한다.
3. 속도는 **Numba**로 확보한다(필수, 최우선). 이후 일부 코드를 C++ 또는 Rust로 바꿀 수 있다.
4. 코드 위치는 저장소 루트 **`WHT_IGA/`**. `wht_topopt`(다른 에이전트 소유)와는 서로 import하지 않는다.
5. Kratos는 **소스빌드본** `E:\kratos-build\install`(10.4.3 + IgaApplication)을 쓴다.
6. 당분간 **Zero-C++를 유지**한다. MPC(MasterSlaveConstraint), Python 조립, SBM 우회 경로가 먼저다.
7. **SBM**(Shifted Boundary Method)은 **비교 트랙**으로만 운영한다.
8. 환경은 **micromamba + uv**. 새 env `whtiga`를 만들고, 셸 프로파일의 conda 훅을 micromamba로 교체한다.
9. 동적 시간적분은 **Kratos scheme**(Bossak 기본, Newmark/BDF 선택)을 쓴다. Kratos 적분법을 HHT의 상위 호환으로 보고, HHT는 구현하지 않는다 (D5 확정).

---

## 1. 비전과 범위

### 1.1 목표

3차원 CAD(STEP/BREP/STL)를 체적 격자 생성 없이 직교 B-스플라인 배경 격자에 잠기게 해서 해석하는 솔버를 Kratos + Python만으로 만든다. "완벽한 구현"은 다음 두 조건으로 정의한다.

1. **기능 동등성**: §1.3 매트릭스의 각 기능이 목표 성숙도에 도달한다.
2. **수치 동등성**: Coreform 검증 매뉴얼 12건과 예제 5건(§10)의 기준값을 지정 허용오차 안에서 재현한다.

### 1.2 성숙도 등급

| 등급 | 의미 | 진입 조건 |
|---|---|---|
| Alpha | 동작은 하지만 검증이 불완전 | 단위시험 통과 + 예제 1개 이상 실행 |
| Beta | 검증 진행 중, 알려진 제약을 문서화 | 매핑된 검증 케이스 1개 이상 허용오차 통과 + 제약 목록 |
| Production | 검증 완료 | 매핑된 검증·예제 **전부** 통과 + 수렴률·성능 리포트 + CI 회귀시험 포함 |

### 1.3 기능 동등성 매트릭스

| 분야 | Coreform 기능 | WHT_IGA 구현 | Phase | 목표 등급 | 검증 |
|---|---|---|---|---|---|
| 격자 | 직교 배경 격자, degree 1–3, continuity C0..C^{p−1} | 균일 텐서 B-spline. 낮은 연속성은 knot 중복 | 1 | Production | V01, V02 |
| 격자 | element size / interval count / 축별 지정 | `RectilinearSizing` | 1 | Production | V01 |
| 격자 | element_layout edge/node_centered, padding | 격자 원점 정렬 규칙 | 1 | Beta | 단위시험 |
| 격자 | Frame(원점, 기저, 범위) | 회전된 Greville 제어점 (detJ=1 유지) | 1 | Beta | V01 회전 변형 |
| 격자 | tessellation_size | OCP BRepMesh 선형 편차 | 1 | Production | V02, V06 (기하 오차 분리) |
| 격자 | Body-fit mesh (정렬 hex) | 정렬 블록 형상은 trim 없이 NurbsVolume 직접 사용. 일반 body-fit은 범위 밖 | 3 | Alpha | V03 body-fit 비교 |
| 격자 | U-spline | **범위 밖**(특허) | – | – | – |
| 세분화 | 계층 세분화, 2:1 균형, CadEntities/Sphere/Box region, offset p·h/2^L | THB + 구속식 | 3 | Beta | V02 local, E02 |
| 소형 셀 | small_cell_volume_ratio=0.2 extension | 기저 지지 기준 extended B-spline (MPC) | 2 | Production | V-COND, V06 |
| 소형 셀 | 압력 셀 aggregation (혼합) | P0 셀 aggregation (C++ 요소) | 5 | Beta | V07, E04 |
| 경계 | penalty 필수 BC (bc_penalty_value) | QP penalty / 표면층 penalty | 1 | Production | V01, V03 |
| 경계 | 변위 지정, 핀, 대칭, 사용자 함수 변위 | BC-A/BC-B (§5.4) | 1 | Production | V01, V07 |
| 경계 | (Nitsche 없음) | Nitsche 선택 | 6+ | Alpha | 단위시험 |
| 하중 | 압력, 총하중 압력, 트랙션, 중력, 체적력 | IgaApp `LoadCondition`, `VOLUME_ACCELERATION` | 1 | Production | V03 |
| 하중 | follower / resultant 트랙션 | `PRESSURE_FOLLOWER_LOAD` | 5 | Beta | V03 변형 |
| 하중 | RP 모멘트 (coupling) | RP 노드 + 분배 coupling 구속 | 3 | Beta | V03 moment |
| 하중 | amplitude | Python 시간함수로 조건값 갱신 | 3 | Beta | V10 |
| 결합 | multi-cell 불리언 병합 (공유 스플라인 장) | 같은 NurbsVolume, 재료별 적분점 요소 | 3 | Production | V04, V05, E03 |
| 결합 | surface tie | 가상노드 쌍 + 스프링 penalty | 3 | Beta | V04 tie |
| 결합 | kinematic / distributing coupling | 강체 구속 + 스프링 / RBE3형 구속 | 3 | Beta | V03 moment |
| 결합 | IGA + 일반 FE 혼합 | 같은 ModelPart에 일반 FE 요소 + tie | 5 | Alpha | 자체 |
| 접촉 | surface-to-surface, penalty Coulomb μ (자기접촉 제외) | 표면층 삼각망 + ContactStructuralMechanics, 실패 시 C++ | 5 | Beta | V08, E05 |
| 해석 | 선형 정적, 다단계 | ResidualBasedLinearStrategy / step 루프 | 1/3 | Production | V01–V06 |
| 해석 | NLGEOM | TotalLagrangian / UpdatedLagrangian + Newton | 5 | Beta | E03, E04 |
| 해석 | 고유진동수 (Lanczos) | SM `EigensolverStrategy` + spectra/FEAST | 4 | Production | V09, V10 |
| 해석 | 암시적 동역학 (Coreform: HHT) | **Kratos 시간적분 scheme** (기본 Bossak, 선택 Newmark/BDF). HHT는 구현하지 않음 (사용자 결정 D5) | 4 | Beta | V10 |
| 해석 | 응답 스펙트럼 (ABS/SRSS/NRL/TENP/CQC/GRP/DSC, ALGEBRAIC/SRSS/R40/R30) | 모드 기반 numpy 후처리 | 4 | Beta | V09 |
| 해석 | 랜덤 응답 (PSD BASE/FORCE/DB, Segalman RMS von Mises) | 모드 FRF 후처리 | 4 | Beta | V09 |
| 해석 | 정상/과도 열전달 (전도, film, radiation, 열용량) | IgaApp `LaplacianElement` + ConvDiff `ThermalFace3D3N` (spike) | 4 | Beta | V12 |
| 해석 | 열변형률 (지정 온도장) | ConstitutiveLaws `ThermalElasticIsotropic3D` + 제어점 투영 T | 4 | Beta | V11 |
| 해석 | explicit 동역학 | **범위 밖**(Coreform도 미지원) | – | – | – |
| 재료 | 선형 탄성 (T 의존 표) | `LinearElastic3DLaw` + TableAccessor | 1/4 | Production | V01, E03 |
| 재료 | 초탄성, Neo-Hookean | `HyperElastic…NeoHookean3DLaw` 등 | 5 | Beta | E04 |
| 재료 | J2 소성 (다선형, T 의존) | `SmallStrainJ2Plasticity3DLaw` 등 + 테이블 (spike) | 5 | Beta | E03 |
| 재료 | 열팽창 (T 의존), 전도도, 비열, 밀도 | Properties + TableAccessor | 4 | Beta | V11, V12 |
| 잠김 방지 | U(p,C^{p−2})/P0 정적 응축, VMS, macro, U2C1/P2C1, B-bar/F-bar | **C++ `WhtIgaApplication`** MixedUP0 요소 (C++ 승격 1순위) | 5 | Beta | V07, E04 |
| 출력 | U, S, PEEQ, E, NE, LE, NT11 | `CalculateOnIntegrationPoints` → results DB | 1–5 | Production | 전체 |
| 출력 | point/line/multi-point probe (가중치 0 적분점) | 가중치 0 적분점 요소 | 3 | Production | V02, V03 |
| 출력 | 시각화 proxy 격자 (untrimmed hex / trimmed tet) | VTKHDF/VTU | 3 | Production | 단위시험 |
| 출력 | results/quadrature SQLite, probe JSON | 같은 구조 | 1/3 | Production | L2 norm 시험 |
| API | portable spec (to_json/from_json/validate → ValidationError) | dataclass spec + validator | 1 | Production | 단위시험 |
| API | Job (build_mesh, write_input, submit, wait_for_completion) | `Job` | 1 | Production | – |
| 병렬 | MPI 전처리, 스레드/MPI 풀이 | 셀 블록 multiprocessing (QP 수 가중), Kratos OpenMP | 6 | Beta | 성능 리포트 |

### 1.4 범위 밖

- **U-spline**: Coreform 특허 핵심이다. 균일 텐서곱 B-스플라인(C^{p−1})과 THB 세분화로 대체한다. Coreform의 Production 경로도 rectilinear이므로 직교 격자에서 동등한 정확도를 기대할 수 있다.
- **Abaqus 호환**: `.inp` 번역기, UEL/UELMAT, `.odb` 출력.
- **Explicit 동역학, 자기접촉**: Coreform도 지원하지 않는다.
- **셸·보 IGA**: Coreform은 3D 솔리드만 제공한다. IgaApplication Shell3p 등은 향후 확장 후보로만 둔다.
- **Altair H3D/hvtrans**: 선택 플러그인(상용 도구 보유 시).
- **Algoim, 정확 형상 적분, quadrature-free**: Phase 6 이후 선택 과제.
- **일반 body-fit 다중패치 체적 매개변수화**: 정렬 블록 외에는 하지 않는다.
- **CFD, 전자기, GPU 실행**.

---

## 2. 설계 원칙

### 2.1 Python-first, C++는 승격 기준을 넘을 때만

1. 1차 구현은 **Kratos Python 바인딩 + 기본 제공 요소·조건·전략**만 쓴다.
2. Kratos가 직접 지원하지 않는 수학은 Kratos가 이미 아는 형태로 **번역해서** 넣는다.
   - 절단 셀 적분점: 셀마다 `NurbsVolumeGeometry.CreateQuadraturePointGeometries(결과, 2, [[x,y,z,w],…])` 1회 호출 → 적분점 기하 1개 → 요소 1개 (부록 A-1, A-2).
   - small-cell extension, THB truncation, tie, coupling: 모두 `LinearMasterSlaveConstraint`로 `u_geom = T·u_master`를 부과한다 (부록 A-4).
   - 경계 적분: `SurfaceInNurbsVolumeGeometry(볼륨, Triangle3D3)` 위의 적분점 기하에 IgaApplication `SupportPenaltyCondition`과 `LoadCondition`을 붙인다 (부록 A-3).
3. **C++ 승격 기준** (§5.7): 다음 중 하나에 해당하면 그 부분만 커스텀 C++ 앱 `WhtIgaApplication`으로 옮긴다.
   - Phase 0 spike에서 Python 경로가 수학적으로 성립하지 않음이 확인됨
   - 필요한 요소가 기본 요소로 표현되지 않음 (예: 혼합 U/P0 잠김 방지 요소)
   - 측정 결과 구속식 처리나 바인딩 오버헤드가 조립 시간의 50%를 넘음

### 2.2 가속 정책: Numba 우선, 전환 대상은 C++ 기본

1. **Numba 우선**: 전처리 핫 루프(분류, 클리핑, 모멘트, 기저 평가, extension, 후처리 샘플링)는 NumPy 벡터화 → `@njit(cache=True, parallel=True)` 순으로 최적화한다.
2. **이식 가능한 형태로 작성**: 커널은 순수 함수(배열 입력 → 배열 출력)로 쓰고, Python 객체를 커널 안에 넣지 않는다. 이후 1:1 이식이 쉽다.
3. **Kratos 안에서 도는 것은 C++만 가능**: 요소, 조건, 시간적분 scheme은 Kratos C++ 클래스이므로 Rust로 대체할 수 없다.
4. **독립 커널의 전환 대상도 C++(pybind11)를 기본으로 한다.**
   - MSVC 14.41과 CMake가 이미 설치되어 있고, Kratos 빌드와 같은 컴파일러다.
   - 언어 하나로 유지보수하고, 커스텀 Kratos 앱과 같은 빌드에 넣을 수 있다.
   - Rust는 미설치라 도구체인이 하나 더 늘어난다. 별도 이유가 생길 때만 선택한다.
5. **전환 시점**: 프로파일링으로 측정한다. Numba로 이미 최적화한 함수가 전체 실행 시간의 30%를 넘을 때만 그 함수를 옮긴다. Numba 버전은 결과 비교용 기준 구현으로 남긴다.

### 2.3 고정 런타임

- Kratos: `E:\kratos-build\install` (10.4.3 소스빌드). pip Kratos wheel과 **섞지 않는다**.
- Python: 3.12.x (Kratos 빌드 ABI와 일치). env는 micromamba `whtiga` (§12).
- 단일 프로세스: 전처리, Kratos 풀이, 후처리를 같은 Python 프로세스에서 실행한다. 메모리 안의 배열을 그대로 Kratos에 넘긴다.
- Kratos 없이도 동작하는 부분(S1–S6 전처리, 후처리, 참조 조립기)은 Kratos를 import하지 않는다. Kratos는 `wht_iga.kratos_backend`에서만 지연 import한다.

### 2.4 정직성 원칙

- 기본값으로 대체된 입력(재료, penalty 등)은 결과 메타데이터에 `assumed`로 기록한다.
- 검증은 **실제 기준값과 비교**한다. `status == completed`만 보는 시험은 금지한다.
- 기하 오차(tessellation)와 근사 오차(이산화)를 따로 측정하고 따로 보고한다.
- Production 등급은 검증 증거 없이 부여하지 않는다.

### 2.5 독립성

- `wht_topopt`와 `aieng`을 import하지 않는다. 필요한 패턴은 출처를 주석에 남기고 복사·적응한다.
- `tests/test_isolation.py`가 `sys.modules`와 소스 grep으로 이를 강제한다.

---

## 3. 시스템 아키텍처

### 3.1 단계 파이프라인

```
[CAD 입력: STEP / BREP / STL]
        │
S1  형상 처리 (OCP)       읽기, 치유(ShapeFix, Sewing), watertight 검사, 셀·면 집합 태깅,
        │                  BRepMesh(tessellation_size) → 셀별 닫힌 삼각망 + face_id
S2  배경 격자 (Numba)     Frame/extents/padding/sizing → 레벨별 knot(물리 좌표), Greville 제어점,
        │                  Bézier 추출 Cᵉ (후처리·참조 조립기용), THB 계층 (Phase 3)
S3  분류 (Numba)          scanline ray parity(3축 투표) + 일반화 winding number fallback
        │                  → 셀 상태 {OUT, IN, CUT} × 재료 셀, 체적비 η_e
S4  절단 적분 (Numba)     삼각형–박스 클리핑 → 재귀 발산정리 모멘트 → 비음수 모멘트 피팅
        │                  (fallback: octree) → 체적 QP(x, w), 경계 QP(x, w, n, face_id)
S5  extension            기저 지지 체적비 → extension 계수 e_ij (Höllig) → T = C_THB · E
        │
S6  경계·표면층           BC 집합별 penalty 적분점, 표면층 삼각망, RP·coupling·tie
        │
S7  Kratos 모델 구축      노드, NurbsVolume, 적분점 기하, 요소, 조건, 구속식, Properties, 구성법칙
S8  해석 (Kratos)         선형/Newton/고유치/동적/열 전략, PARDISO 또는 대체 솔버
S9  후처리                results DB, quadrature DB(X, dV, U, S), probe JSON, VTKHDF proxy 격자,
                          응답 스펙트럼/랜덤, L2 노름
```

S1–S6은 Kratos 없이 실행된다(Coreform의 "mesh stage"에 해당). 결과는 `<job>.wiga.sqlite`에 저장된다.

### 3.2 Coreform 단계와의 대응

| Coreform | WHT_IGA |
|---|---|
| Abaqus/CAE 모델 설정 | Python API `wht_iga.api` |
| `<job>.cf` export | `ModelSpec` JSON (`Job.to_json`) |
| `coreform_iga_mesh` → `<job>.sql` | S1–S6 → `<job>.wiga.sqlite` |
| interop (SFM3D3 표면층, `*USER ELEMENT`) | S6–S7 (Kratos 조건 + MPC + 표면층 노드) |
| Abaqus/Standard 풀이 | S8 Kratos |
| `_iga.odb`, `_iga_results.sql`, `probe_data.json`, `_iga_quadrature.sql` | `<job>.vtkhdf`, `<job>_iga_results.sqlite`, `<job>_probe_data.json`, `<job>_iga_quadrature.sqlite` |

### 3.3 실행 형태

- 기본은 **`whtiga` env 단일 프로세스**다.
- Kratos 로딩: `wht_iga/kratos_backend/_bootstrap.py`가 `WHT_IGA_KRATOS_PATH`(기본 `E:\kratos-build\install`)와 `WHT_IGA_KRATOS_DLL_PATH`(기본 `E:\kratos-build\install\libs;D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin`)를 읽어 `sys.path`와 `os.add_dll_directory`를 설정한다. import 후 `KratosMultiphysics.__file__`이 지정 경로 아래인지 검증하고, 아니면 오류로 중단한다.
- 캐시: `Job.build_mesh()`가 S1–S6 결과를 저장한다. 캐시 키는 형상 해시 + 격자 파라미터 해시이고, 같으면 재사용한다.
- 외부(예: Python 3.11 aieng) 연동용 `wht_iga/runner.py`는 `subprocess`로 `micromamba run -n whtiga python -m wht_iga run job.json`을 실행하고 산출물 경로를 돌려준다. 도입 시점은 §14 결정 사항이다.
- 전처리 병렬화: `multiprocessing`(spawn)으로 셀 블록 단위 작업을 나눈다. 분할 가중치는 예상 적분점 수(CUT 셀은 IN 셀의 10–30배)로 정한다.

---

## 4. 데이터 버스

v1의 `.mdpa` 경로는 폐기한다. Kratos 모델은 §5 방식으로 Python에서 직접 구축하고, SQLite는 **캐시, 재현성, 후처리** 용도로 쓴다.

### 4.1 공통 규칙

- BLOB: little-endian float64/int64, `ndarray.tobytes()`. 배열 형상은 별도 열에 둔다. 큰 배열은 셀 블록 단위로 나눈다.
- 모든 DB에 `meta` 표: `schema_version`, `wht_iga_version`, `kratos_version`, `git_hash`, `units`, `spec_json`, `created_at`.

### 4.2 `<job>.wiga.sqlite` (mesh stage)

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

-- 형상과 기하 오차 기록
CREATE TABLE cad_cell (cell_id INTEGER PRIMARY KEY, name TEXT, material TEXT,
                       volume_tess REAL, volume_brep REAL);
CREATE TABLE cad_face (face_id INTEGER PRIMARY KEY, cell_id INTEGER, name TEXT,
                       area_tess REAL, area_brep REAL, max_dev REAL);
CREATE TABLE cad_set  (set_name TEXT, kind TEXT, entity_id INTEGER);   -- face|cell|edge|vertex
CREATE TABLE tess_tri (face_id INTEGER, block INTEGER, n INTEGER, xyz BLOB, tri BLOB);

-- 격자
CREATE TABLE grid_level (level INTEGER PRIMARY KEY, frame BLOB,
                         knots_u BLOB, knots_v BLOB, knots_w BLOB,
                         degree INTEGER, continuity INTEGER, box_lo BLOB, box_hi BLOB);
CREATE TABLE cell (cell_id INTEGER PRIMARY KEY, level INTEGER, i INTEGER, j INTEGER, k INTEGER,
                   status INTEGER,          -- 0 OUT, 1 IN, 2 CUT
                   eta REAL,                -- |Q∩Ω|/|Q|
                   mat_mask INTEGER, n_qp INTEGER, qp_offset INTEGER);
CREATE TABLE basis (basis_id INTEGER PRIMARY KEY, level INTEGER, i INTEGER, j INTEGER, k INTEGER,
                    greville BLOB, support_ratio REAL,
                    role INTEGER);           -- 0 inactive, 1 master, 2 extended, 3 thb_passive

-- 연산자
CREATE TABLE operator_T (row_geom_dof INTEGER, col_master INTEGER, coeff REAL);   -- 희소 COO, THB∘extension 합성
CREATE TABLE extension  (slave_id INTEGER, master_id INTEGER, coeff REAL);         -- 진단용 원본
CREATE TABLE extraction (cell_id INTEGER PRIMARY KEY, C BLOB);                     -- Bézier Cᵉ (참조 조립기·후처리)

-- 적분점 (모두 물리 좌표; 항등 사상이므로 파라미터 좌표와 같음)
CREATE TABLE qp_volume  (cell_id INTEGER, material INTEGER, block INTEGER, xw BLOB);   -- [x,y,z,w]×n
CREATE TABLE qp_surface (face_id INTEGER, cell_id INTEGER, block INTEGER, xwn BLOB);  -- [x,y,z,w,nx,ny,nz]×n

-- 표면층 (Coreform SFM3D3 대응)
CREATE TABLE layer_node (layer_id INTEGER, node_id INTEGER, x REAL, y REAL, z REAL, area_w REAL, face_id INTEGER);
CREATE TABLE layer_tri  (layer_id INTEGER, tri_id INTEGER, n1 INTEGER, n2 INTEGER, n3 INTEGER);

-- 프로브, 검사 결과
CREATE TABLE probe (probe_id INTEGER PRIMARY KEY, name TEXT, kind TEXT, xyz BLOB, variables TEXT);
CREATE TABLE issue (severity TEXT, code TEXT, message TEXT, entity TEXT);   -- validate() 결과
```

### 4.3 `<job>_iga_results.sqlite`

```sql
CREATE TABLE step  (step_id INTEGER PRIMARY KEY, name TEXT, kind TEXT, t_end REAL);
CREATE TABLE frame (frame_id INTEGER PRIMARY KEY, step_id INTEGER, time REAL, increment INTEGER,
                    converged INTEGER, n_iter INTEGER);
CREATE TABLE cp_field (frame_id INTEGER, var TEXT, data BLOB);      -- 제어점 계수 (U, NT11, 압력)
CREATE TABLE qp_field (frame_id INTEGER, var TEXT, data BLOB);      -- 적분점 상태 (S, PEEQ, E, LE …)
CREATE TABLE mode  (mode_id INTEGER PRIMARY KEY, freq_hz REAL, gen_mass REAL,
                    part_x REAL, part_y REAL, part_z REAL, eff_mass_x REAL, eff_mass_y REAL, eff_mass_z REAL);
CREATE TABLE energy (frame_id INTEGER, kinetic REAL, strain REAL, plastic REAL);
```

### 4.4 `<job>_iga_quadrature.sqlite` (선택, Coreform `--quadrature-output` 대응)

```sql
CREATE TABLE qp     (qp_id INTEGER PRIMARY KEY, cell_id INTEGER, material INTEGER,
                     x REAL, y REAL, z REAL, dV REAL);
CREATE TABLE qp_val (frame_id INTEGER, qp_id INTEGER, U BLOB, S BLOB);   -- L2 오차 노름 계산용
```

### 4.5 기타 출력

- `<job>_probe_data.json`: Coreform `probe_data.json` 구조를 따른다.
- `<job>.vtkhdf` / `.vtu`: 시각화 proxy 격자 (§7.12).
- `<job>.spec.json`: 입력 재현용 `ModelSpec`.

---

## 5. Kratos 연동 설계

### 5.1 확인된 Python 노출 API (요약, 상세는 부록 A)

| 기능 | Python 노출 | 비고 |
|---|---|---|
| `NurbsVolumeGeometry(points, pU, pV, pW, kU, kV, kW)` (비유리) | 예 | |
| `geom.CreateQuadraturePointGeometries(res, nderiv)` | 예 | 기본 Gauss 전체 |
| `geom.CreateQuadraturePointGeometries(res, nderiv, [[u,v,w,wt],…])` | 예 | **Option 1**: 주어진 점 전부를 담은 적분점 기하 **1개**. 제어점은 점들의 중심이 속한 knot span 기준. `nderiv == 2` 필수 |
| `IntegrationInfo` 선택(CUSTOM), `QuadraturePointGeometry` 생성자, 형상함수 컨테이너 | **아니오** | Python에서 Cᵉ·B 값 직접 주입 불가 |
| `SurfaceInNurbsVolumeGeometry(volume, surface)` | 예 | 표면 기하의 전역 좌표를 볼륨 국소 좌표로 해석 |
| `BrepVolume`, `CouplingGeometry` | 아니오 | |
| `ShapeFunctionsValues`, `ShapeFunctionDerivatives`, `DeterminantOfJacobian(i)`, `Normal(i)` | 예 | 검증용 |
| `LinearMasterSlaveConstraint`, `ModelPart.CreateNewMasterSlaveConstraint`, `node.GetDof(var)` | 예 | |
| Builder-and-solver `Build`, `ApplyConstraints`, `GetDofSet`, strategy `GetSystemMatrix`, `scipy_conversion_tools` | 예 | 행렬 추출 경로 |
| `LinearSolversApplication.HasMKL()`, `pardiso_lu/ldlt/llt`, `sparse_lu`, `spectra_sym_g_eigs_shift`, FEAST | 예 | pardiso는 MKL 빌드일 때만 등록 |
| 시간적분 scheme: Bossak(`ResidualBasedBossakDisplacementScheme`, SM `StructuralMechanicsBossakScheme`), Newmark, BDF(displacement), Static | 예 | HHT 없음 → **Kratos scheme을 그대로 사용** (D5) |
| IgaApp `SupportPenaltyCondition` | 예 | 3성분 동시 penalty (성분별 마스크 없음) |
| IgaApp `LoadCondition` | 예 | `SURFACE_LOAD`, `PRESSURE`(기하 법선), `PRESSURE_FOLLOWER_LOAD`, `DEAD_LOAD` |
| IgaApp `SolidElement`, `SupportSolidCondition` | 2D 전용 | 3D 솔리드는 StructuralMechanics 요소 사용 |
| IgaApp `AssignIntegrationPointsToBackgroundElementsProcess` | 예 | **결과 매핑용**(output 단계). 해석 경로가 아니다 |
| SM `NodalConcentratedElement3D1N`, `SpringDamperElement3D` | 예 (등록 확인) | 표면층·부분 BC용 |
| ConstitutiveLaws `ThermalElasticIsotropic3D` | 예 (등록 확인) | 열변형 |

### 5.2 핵심 기법

#### 5.2.1 항등 사상 트릭

- knot 벡터를 물리 좌표로 정의하고 제어점을 Greville 점에 놓는다. 균일 B-spline은 선형 정밀도를 가지므로 x(ξ) = ξ, detJ = 1이다. Frame 회전이 있어도 회전 행렬이므로 detJ = 1이 유지된다.
- 결과:
  - 적분점의 국소 좌표 = 물리 좌표, 가중치 = 물리 dV. 점 역매핑이 필요 없다.
  - CAD 삼각형의 물리 좌표를 그대로 `SurfaceInNurbsVolumeGeometry`에 넣을 수 있다.
  - 경계 적분점에서 `weight × detJ`가 정확히 면적 요소가 된다.
- 확인은 spike S02.

#### 5.2.2 셀 단위 적분점 기하 (Option 1)

- 각 (셀, 재료) 쌍마다 `volume.CreateQuadraturePointGeometries(g, 2, pts)`를 1회 호출한다. 결과는 n개 적분점을 담은 기하 1개다.
  - IN 셀: (p+1)³ Gauss 점
  - CUT 셀: 모멘트 피팅 점
- 제어점 집합은 점들의 중심이 속한 knot span에서 정해지므로, 셀 경계에 정확히 걸친 점은 ε만큼 셀 안쪽으로 옮긴다 (spike S03).
- 이 기하에 StructuralMechanics 요소(선형: `SmallDisplacementElement3D4N`, NLGEOM: `TotalLagrangianElement3D4N` / `UpdatedLagrangianElement3D4N`)를 붙인다. 이름의 "4N"은 prototype이고, 실제 기하는 생성 시 교체된다(`test_nurbs_volume_element.py`가 실증).

#### 5.2.3 구속식으로 기저 변환 (`u_geom = T · u_master`)

- 원래 B-spline 기저로 요소를 만들고, 기저 변환은 자유도 사이의 선형 구속으로 부과한다. 조립 결과는 K_m = Tᵀ K T이고, 요소에 C̃ᵉ = Cᵉ Tₑ를 직접 넣는 것과 수학적으로 같다.
- 하나의 틀로 다음을 모두 처리한다.
  - small-cell extension (Höllig extended B-spline: `u_slave = Σ e_ij u_master`)
  - THB truncation
  - tie, coupling, 가상노드 보간
- 규칙:
  1. **노드**: 레벨 L의 제어점마다 Kratos 노드 1개("geometry node"). CUT/IN 셀이 참조하는 노드만 생성한다.
  2. **knot 규약**: Kratos는 끝 knot을 p번 반복하는 축약 규약을 쓴다(길이 = 제어점 수 + p − 1). `grid/knots.py`에서 변환한다.
  3. **균일 격자 + extension만**: master = 지지가 충분한 geometry node 자신, slave = extension 대상 노드뿐.
  4. **THB 포함**: 모든 geometry node가 slave이고, master는 별도의 THB 자유도 노드다. 연쇄 구속(slave가 다른 구속의 master가 되는 경우)을 만들지 않기 위해서다.
  5. master는 절대 slave가 될 수 없으므로, 모든 합성(THB ∘ extension ∘ 가상노드 보간)은 Python에서 미리 계산한다.
  6. 지지가 OUT 셀에만 있는 기저는 노드를 만들지 않는다.
- 구속식 1개 = slave 노드 1개 × 성분 3개(열해석은 1개). 관계 행렬은 3×3m 블록 대각(m = 해당 T 행의 비영 개수). 열과 구조가 같은 T를 공유한다.
- 이후 C++ `ExtendedSplineSolidElement`(Phase 6)는 같은 T를 C̃ᵉ로 직접 쓴다. 데이터 경로는 바뀌지 않는다.

### 5.3 경로 선택

| 경로 | 내용 | 판정 |
|---|---|---|
| **R1** = (a) + (c) | 셀 단위 적분점 기하 + SM 요소 + 구속식 T | **주 경로** |
| (b) | Python에서 적분점 기하에 Cᵉ·B 직접 주입 | **불가**(바인딩 없음) → (e)로 이관 |
| **R2** | Kratos가 geometry node 기준으로 조립 → `GetSystemMatrix`/`scipy_conversion_tools`로 추출 → scipy에서 Tᵀ K T, Tᵀ f → pypardiso. 비선형은 Python Newton 루프에서 `Build` 반복 호출 | S05(구속식) 실패 시 |
| **R3** | Numba 참조 조립기가 주 경로. Kratos는 구성법칙 평가에만 쓰거나 쓰지 않음. C++ (e) 즉시 착수 | S04(요소 부착) 실패 시 |
| (e) | C++ `WhtIgaApplication`: `ExtendedSplineSolidElement`, `MixedUP0SplineElement`, Nitsche/ghost | 승격 기준 충족 시 (§5.7) |
| SBM | IgaApp `SbmSolidCondition` + `NurbsGeometryModelerSbm` | 비교 트랙 (§10.5) |

**ADR-001 결정 트리**: S03·S04·S05 통과 → R1. S05만 실패 → R2. S04 실패 → R3 + C++ 조기 착수.

### 5.4 경계조건, 하중, 표면층

경계 적분점 기하: 셀 안에서 클리핑된 CAD 삼각형 조각을 독립 노드의 `Triangle3D3`으로 만들고, `SurfaceInNurbsVolumeGeometry(volume_L, tri)` → `CreateQuadraturePointGeometries(g, 2, [[ξ,η,0,w],…])`로 표면 적분점 기하를 만든다. 이 기하는 표면 Jacobian과 법선을 가진다. (S07 실패 시 대안: 체적 적분점 기하에 가중치 = dA를 넣고, 법선은 조건 변수로 전달.)

| BC/하중 | 방식 | Kratos 객체 | 비고 |
|---|---|---|---|
| 3성분 지정 변위 (고정, V01) | **BC-A**: 적분점 penalty 조건 | IgaApp `SupportPenaltyCondition`. Properties `PENALTY_FACTOR`=β, 조건 `DISPLACEMENT`=g | S09 확인 |
| 일부 성분, 임의 방향 (대칭면 등) | **BC-B**: 가상노드 스프링 | 적분점마다 가상노드 V를 만들고 `LinearMasterSlaveConstraint`로 V의 자유도를 Σ N_i T u_m(또는 n·u)에 연결. `NodalConcentratedElement3D1N`(성분별 강성 β·w_q)과 `PointLoadCondition3D1N`(β·w_q·g) 부착 | 회전 좌표계 성분은 관계 행에 n을 곱해 처리 |
| Coreform식 표면층 (접촉, tie, coupling, film) | **BC-C**: 표면층 노드 + 스프링 penalty | 표면층 삼각망 정점 L(독립 자유도)과 가상노드 V(스플라인 slave)를 `SpringDamperElement3D`로 연결, 강성 β·A_vertex. L은 일반 노드이므로 `Fix`, 접촉, film을 그대로 쓸 수 있다 | Coreform SFM3D3 penalty 층과 수학적으로 같다 |
| 압력 / follower 압력 | `LoadCondition` `PRESSURE` / `PRESSURE_FOLLOWER_LOAD` | 표면 적분점 기하 | 법선 부호는 S08로 확인 |
| 트랙션 (resultant) | `SURFACE_LOAD` | 적분점별 또는 면 단위 상수 | |
| 총하중 압력 | p = F / Σ w_q (tessellation 면적) | 위와 같음 | |
| 중력·체적력 | Properties `VOLUME_ACCELERATION` | | |
| RP 모멘트, distributing coupling | RP 노드(변위+회전)를 slave, master = 해당 면 가상점들의 스플라인 자유도(합성 T). 관계식은 RBE3 가중 최소자승(선형). `PointMomentCondition3D1N` | | |
| kinematic coupling | 면 위 점 r마다 강체 노드 R(u = u_RP + θ×r)과 가상노드 V를 스프링으로 연결 | penalty kinematic | |
| tie (별개 부품) | 종속면 적분점마다 V_A(부품 A slave), V_B(부품 B slave), 스프링 β·w_q | | |

**penalty 크기**
- 사용자가 `bc_penalty_value`를 주면 그 절대값을 쓴다(Coreform 동등, 검증에서는 1e10).
- 주지 않으면 β = c_pen · E_max · p² / h_min, c_pen = 1e3. Phase 1의 β sweep으로 확정한다(ADR-002).
- 오차는 O(1/β)이고 조건수는 β·h에 비례해 커진다. 이 관계를 검증 리포트에 기록한다.

**Nitsche (선택)**: IgaApp `SupportNitscheCondition`이 3D를 지원하는지 S18로 확인한다. 지원하지 않으면 C++로 넘긴다. γ는 β = C·E·p²/h 또는 `NitscheStabilizationModelPartProcess`(국소 고유치). h는 extension root 셀 크기를 쓴다.

### 5.5 Numba 참조 조립기 (`wht_iga.reference`)

- 선형 탄성·열 문제의 요소 행렬을 Bézier 추출(Cᵉ)로 직접 조립하고, Tᵀ K T를 scipy에서 만들어 pypardiso로 푼다.
- 용도:
  1. **oracle**: 모든 검증 케이스 축소판에서 Kratos 경로와 1e-10까지 비교한다.
  2. **fallback**: R3 경로.
  3. **조건수 연구**: κ 계산용 행렬을 바로 얻는다.

### 5.6 해석 종류별 매핑

| 해석 | Kratos | 대안 |
|---|---|---|
| 선형 정적 | `ResidualBasedLinearStrategy` + `ResidualBasedIncrementalUpdateStaticScheme` + `ResidualBasedBlockBuilderAndSolver` + `pardiso_ldlt` | `sparse_lu`, 또는 추출 후 pypardiso |
| 비선형 | `ResidualBasedNewtonRaphsonStrategy` + `DisplacementCriteria`/`ResidualCriteria`, step 루프 | line search |
| 고유치 | SM `EigensolverStrategy` + `spectra_sym_g_eigs_shift` 또는 FEAST | 추출 K, M + scipy `eigsh(sigma)` |
| 선형·비선형 동적 | `ResidualBasedNewtonRaphsonStrategy`(비선형) / `ResidualBasedLinearStrategy`(선형) + **Bossak scheme** (`damp_factor_m` α_B, 기본 −0.3). Rayleigh 감쇠는 Properties | Newmark(α=0), BDF2 |
| 열 정상/과도 | IgaApp `LaplacianElement` + `SupportLaplacianCondition`(penalty T) + BDF. film·radiation은 표면층 노드(TEMPERATURE slave) 위 ConvDiff `ThermalFace3D3N` | 참조 조립기 열해석 |
| 열변형 | `ThermalElasticIsotropic3D` (노드 TEMPERATURE 보간, `REFERENCE_TEMPERATURE`). 온도장은 제어점에 L2 투영 | SM 초기 변형률 주입 |
| 응답 스펙트럼 / 랜덤 | 모드 + 참여계수 Γ = φᵀ f_r. f_r은 `VOLUME_ACCELERATION = e_k`로 `BuildRHS` → numpy 후처리 | – |

### 5.7 C++ 승격 기준 (`WHT_IGA/cpp/WhtIgaApplication`)

다음 중 하나라도 해당하면 승격한다.

1. **잠김 방지 혼합 요소** (V07, E04 필수). 1순위.
2. S05에서 구속식 오버헤드가 조립 시간의 50%를 넘음.
3. Nitsche 또는 ghost penalty가 필요함.
4. 접촉 S14 실패.

빌드: MSVC 14.41 + CMake/Ninja, Kratos `application_generator` 템플릿, 외부 앱으로 `E:\kratos-build\install`과 `KratosIgaCore.lib`에 링크한다. 요소는 C̃ᵉ와 적분점을 SQLite 또는 Python 컨테이너로 받는다. Kratos 재빌드 절차는 §12.4.

### 5.8 Phase 0 spike 목록

각 spike는 `WHT_IGA/tests/kratos_spikes/test_sNN_*.py`로 만들고 결과를 `docs/spikes/SNN.md`와 ADR로 남긴다.

| ID | 질문 (예/아니오) | 통과 기준 | 실패 시 |
|---|---|---|---|
| S01 | `whtiga` env에서 Kratos, Iga, SM, LinearSolvers, ConstitutiveLaws, ConvDiff, ContactSM import 성공? `HasMKL()` True? | import 성공, pardiso 등록 | vdmc Python을 런타임으로. MKL 없으면 `sparse_lu`, 대형은 추출 후 pypardiso |
| S02 | 물리 knot + Greville 제어점 NurbsVolume에서 `GlobalCoordinates(ξ) == ξ`, 적분점 기하 `DeterminantOfJacobian == 1`? (Frame 회전 포함) | 오차 1e-14 | 가중치를 detJ로 나눠 보정 |
| S03 | Python `CreateQuadraturePointGeometries(g, 2, pts)`가 셀당 기하 1개(n 점), 가중치 보존, 분할단위성, 올바른 제어점 집합을 만드나? 점이 셀 면 위에 있을 때도? | 모두 예 | 경계 점은 ε 이동. 안 되면 점마다 1회 호출 |
| S04 | SM `SmallDisplacement`/`TotalLagrangian`/`UpdatedLagrangian` 3D 요소가 다점 적분점 기하로 초기화되고, K가 참조 조립기와 1e-12까지 일치하며, `CalculateOnIntegrationPoints`(응력, 변형률, von Mises)가 동작하나? | 예 | 다점 실패 → 점당 요소. 전면 실패 → R3 + C++ 조기 착수 |
| S05 | `LinearMasterSlaveConstraint`(3 slave × 81 master) 1e5개를 30 s 안에 생성, Block B&S로 풀이, slave 값 복원? 요소가 없는 master 전용 노드의 자유도가 DofSet에 들어가나? | 모두 예 | 느리면 배치 생성. DofSet 누락 → master 노드에 강성 0 `NodalConcentratedElement`. 실패 → R2 |
| S06 | 구속이 있을 때 `EigensolverStrategy`가 가짜 모드 없이 참조 스펙트럼과 일치? | 앞 10모드 1e-8 | 추출 행렬 + scipy `eigsh` |
| S07 | `SurfaceInNurbsVolumeGeometry(volume, Triangle3D3)` + 사용자 적분점의 `DeterminantOfJacobian`이 면적 요소이고 `Normal`이 삼각형 방향과 일치? | 면적 합 1e-13 | 체적 적분점 기하 + 가중치 dA, 법선은 조건 변수 |
| S08 | `LoadCondition` PRESSURE/SURFACE_LOAD/FOLLOWER가 위 기하에서 총하중을 정확히 적분? 부호 규약은? | 1e-12 | SURFACE_LOAD = −p·n 사전 계산 |
| S09 | `SupportPenaltyCondition`이 위 기하에서 ∫ β (u − g)·v dA를 주나? | 예 | BC-B로 통일 |
| S10 | slave 노드 위 `NodalConcentratedElement3D1N`과 `SpringDamperElement3D`가 정상 동작? | 예 | 표면층을 slave로 두는 BC-B 변형, 접촉은 C++ |
| S11 | IgaApp `LaplacianElement` + `SupportLaplacianCondition` 3D, BDF 과도, slave 표면층 위 `ThermalFace3D3N` film/radiation 동작? | 1D 해석해 1e-6 | 참조 조립기 열해석이 주 경로 |
| S12 | `ThermalElasticIsotropic3D`가 적분점 기하 노드 TEMPERATURE로 열변형률 계산? | 균일 ΔT 자유팽창 1e-12 | 초기 변형률 주입 |
| S13 | `SmallDisplacementMixedVolumetricStrainElement*`가 적분점 기하를 받나? (실패 예상) | – | C++ MixedUP0 (Phase 5) |
| S14 | ContactSM penalty/ALM mortar가 독립 자유도 표면층(스프링 결합)에서 동작? | Hertz 축소판 수렴 | C++ 또는 Python 활성집합 penalty |
| S15 | `TableAccessor`로 E(T), α(T), σ_y(T)를 Python에서 설정 가능? | 예 | 온도 구간마다 요소 Properties 갱신 |
| S16 | 적분점 기하 + 요소 1e5개 생성 시간·메모리 (p=2, 3)? | 1e5개 < 60 s, < 8 GB | 셀 병합, C++ 생성기 앞당김 |
| S17 | `OutputQuadratureDomainProcess`, `iga_vtk_output_process` 재사용 가능? | 정보성 | 자체 writer |
| S18 | `SupportNitscheCondition` 3D 지원? | 정보성 | Nitsche는 C++로 |

---

## 6. 디렉토리·모듈 구조 (`WHT_IGA/`)

```
WHT_IGA/
├── pyproject.toml             # name=wht_iga, requires-python ">=3.12,<3.13"
│                              # deps: numpy, scipy, numba
│                              # extras: cad=[cadquery-ocp], viz=[pyvista, meshio, h5py], solve=[pypardiso],
│                              #         test=[pytest, hypothesis]
├── uv.lock
├── env/whtiga.yml             # micromamba spec (python=3.12, mkl, hdf5 …)
├── AGENTS.md                  # 다른 에이전트용 강제 규칙
├── CLAUDE.md                  # "@AGENTS.md"
├── README.md
├── docs/
│   ├── ENVIRONMENT.md         # 환경 단일 원본 (§12)
│   ├── KRATOS_BUILD.md        # Kratos 소스빌드 재현 절차 (§12.4)
│   ├── adr/                   # ADR-001 경로, ADR-002 penalty 기본값 …
│   ├── spikes/                # S01–S18 결과
│   └── verification/          # 케이스별 리포트 (자동 생성)
├── scripts/check_env.py       # Python 버전, Kratos 경로, import, numba 동작 검사
├── wht_iga/
│   ├── __init__.py            # 공개 API
│   ├── api/                   # Coreform형 사용자 API
│   │   ├── spec.py            # SpecBase(to_json/from_json/validate → ValidationError(issues))
│   │   ├── mesh.py            # ImmersedRectilinearMesh, BodyFitMesh, RectilinearSizing, Frame
│   │   ├── refinement.py      # Refinement, CadEntitiesRegion, SphereRegion, BoxRegion
│   │   ├── model.py           # Model, Part, Cell, Material, Assembly
│   │   ├── steps.py           # StaticStep, NonlinearStaticStep, FrequencyStep, DynamicStep, HeatTransferStep
│   │   ├── loads.py, bcs.py, interactions.py, amplitudes.py
│   │   ├── probe.py           # PointProbe, LineProbe, MultiPointProbe
│   │   └── job.py             # Job: build_mesh, write_input, submit, wait_for_completion, results
│   ├── geometry/              # S1
│   │   ├── occ_io.py          # OCP STEP/BREP 읽기, ShapeFix, Sewing, 셀/면 이름
│   │   ├── tessellate.py      # BRepMesh_IncrementalMesh, face_id 보존, 정점 병합
│   │   ├── mesh_io.py         # STL/OBJ 입력
│   │   ├── watertight.py      # 매니폴드, 오일러수, 체적 부호, 틈 검사
│   │   ├── selectors.py       # 이름/박스/평면/법선으로 면 선택
│   │   └── primitives.py      # 검증 형상 생성기 (Roark 판, 구 octant, 보, Cook, Hertz …)
│   ├── grid/                  # S2
│   │   ├── frame.py, sizing.py, knots.py (Kratos 축약 규약 변환), greville.py
│   │   ├── bspline.py         # @njit Cox–de Boor 값·미분
│   │   ├── bernstein.py, bezier_extraction.py   # Borden 2011
│   │   └── hierarchy.py       # THB 레벨, 2:1 균형, 세분 행렬, truncation
│   ├── classify/              # S3
│   │   ├── binning.py         # 삼각형→셀 SAT (@njit)
│   │   ├── raycast.py         # 격자선 scanline parity, 3축 투표
│   │   ├── winding.py         # 일반화 winding number (fallback)
│   │   └── cells.py           # 상태, η, 재료 마스크
│   ├── quadrature/            # S4
│   │   ├── gauss.py, triangle_rules.py (Dunavant)
│   │   ├── clip.py            # 삼각형∩AABB, Sutherland–Hodgman (@njit)
│   │   ├── moments.py         # 재귀 발산정리 다항식 모멘트
│   │   ├── momentfit.py       # 비음수 모멘트 피팅 + 점 제거
│   │   ├── octree.py          # fallback
│   │   ├── surface.py         # 경계 적분점 (x, w, n, face_id)
│   │   └── checks.py          # 체적/면적/모멘트 잔차
│   ├── stabilize/             # S5
│   │   ├── support.py         # 기저 지지 체적비
│   │   ├── extension.py       # Höllig extension 계수
│   │   └── fcm.py             # FCM α 기준선
│   ├── operators/
│   │   ├── compose.py         # T = C_THB · E, 가상노드 보간 합성 (scipy.sparse)
│   │   └── numbering.py
│   ├── bus/                   # SQLite 스키마와 입출력
│   │   ├── schema_mesh.sql, schema_results.sql, schema_quad.sql
│   │   └── blob.py, mesh_db.py, results_db.py, quad_db.py
│   ├── kratos_backend/        # S7–S8 (Kratos는 여기서만 import)
│   │   ├── _bootstrap.py      # 경로·DLL 설정, 빌드 위치 검증
│   │   ├── builder.py         # 노드, NurbsVolume, 적분점 기하, 요소
│   │   ├── constraints.py     # LinearMasterSlaveConstraint 배치 생성
│   │   ├── bcs.py             # BC-A/B/C, RP, tie, coupling
│   │   ├── loads.py, materials.py (구성법칙 매핑, TableAccessor)
│   │   ├── solvers.py         # HasMKL 분기
│   │   ├── analyses/          # static.py, nonlinear.py, modal.py, hht_linear.py, thermal.py
│   │   └── extract.py         # 적분점 결과, 제어점 계수 → bus
│   ├── reference/             # Numba oracle / fallback
│   │   └── assemble_elastic.py, assemble_thermal.py, penalty.py, solve.py
│   ├── post/                  # S9
│   │   ├── sampling.py        # 임의 점 장 평가 (@njit)
│   │   ├── proxy_mesh.py      # IN 셀 hex / CUT 셀 tet
│   │   ├── vtk_writer.py      # VTKHDF(h5py) / VTU(meshio)
│   │   ├── probes.py, norms.py (L2, H1, 에너지)
│   │   ├── spectrum.py        # 응답 스펙트럼 조합
│   │   └── random.py          # PSD, Segalman RMS von Mises
│   ├── cli.py                 # python -m wht_iga {mesh, run, post, verify}
│   └── runner.py              # 외부 프로세스 실행기
├── cpp/                       # (Phase 5+) WhtIgaApplication
├── tests/
│   ├── unit/                  # Kratos 불필요
│   ├── kratos_spikes/         # S01–S18 (marker: kratos)
│   ├── verification/          # V01–V12, V-COND (marker: kratos, slow)
│   │   └── data/*.json        # 문서에서 옮긴 치수·하중 (출처 기록)
│   ├── examples/              # E01–E05
│   └── test_isolation.py      # wht_topopt/aieng import 금지
├── examples/                  # 사용자 스크립트
└── benchmarks/
```

---

## 7. 알고리즘과 파라미터

### 7.1 S1 형상 처리

- OCP `STEPControl_Reader` → `ShapeFix_Shape` → `BRepBuilderAPI_Sewing`(허용오차 = 1e-6 × bbox 대각).
- 셀(solid)마다 `BRepMesh_IncrementalMesh(shape, 선형편차 = tessellation_size, 절대값, 각편차 = 0.3)`.
- 면별 `BRep_Tool.Triangulation_s`로 정점·삼각형을 얻고, 공유 모서리 정점은 용접 허용오차로 병합한다. `face_id`를 보존한다.
- 기본 `tessellation_size` = min(h)/4. 기하 오차(면적 오차, BRep 대비 최대 거리를 `BRepExtrema`로 샘플링)는 검증 리포트에 기록한다.
- **멀티셀**(불리언 병합): `BOPAlgo_Builder`(general fuse)로 공유 면을 가진 셀들을 만든다. 내부 계면은 양쪽 셀 삼각망이 **같은 삼각형**을 공유하게 해서, 재료별 적분점이 틈이나 겹침 없이 셀을 나누게 한다.
- watertight 검사 실패, 체적 부호 오류 등은 `issue` 표에 기록하고, 심각도 `error`면 `validate()`가 해석을 막는다.

### 7.2 S2 배경 격자

- extents = bbox + padding·h (기본 padding = 1, Coreform 검증은 padding = p).
- `element_layout`: `edge_centered`(셀 경계가 원점) 또는 `node_centered`.
- knot: 균일 clamped, 물리 좌표. continuity k < p−1이면 내부 knot을 p−k번 반복한다.
- 제어점: Greville 점. Frame이 있으면 회전을 적용한다.
- Bézier 추출: Borden 2011의 knot 삽입 알고리즘으로 축별 1D Cᵉ를 구하고, 텐서곱 Cᵉ = C_w ⊗ C_v ⊗ C_u. 균일 내부 셀은 Cᵉ가 같으므로 고유 패턴만 캐시한다.
- **THB** (Phase 3):
  - 레벨 L에서 h_L = h/2^L. 세분 영역 셀 표시(offset 기본 p·h/2^L), 2:1 균형(이웃 셀 레벨 차 ≤ 1).
  - 활성 함수: supp ⊆ Ω^L 이고 supp ⊄ Ω^{L+1}.
  - truncation: 상위 레벨 함수에서 하위 레벨 활성/통과 함수 성분을 제거한다(Giannelli 2012).
  - 1D 세분 행렬의 텐서곱으로 C_THB 행을 만든다.
  - 레벨 L NurbsVolume은 세분 영역 bbox로 잘라 생성해 노드 수 폭증을 막는다.

### 7.3 S3 분류

- 삼각형 bbox → 셀 범위 → SAT(Akenine-Möller) 판정으로 CUT 셀 집합 (@njit parallel).
- 격자 정점 분류: x/y/z 축 격자선마다 광선을 쏴서 교차점을 정렬하는 scanline parity. 3축 다수결. 비용 O(격자선 수 × 교차 수).
- 3축 결과가 갈리는 정점은 일반화 winding number(Jacobson 2013, BVH 근사)로 판정한다(임계 0.5).
- 셀 상태: CUT이 아니면 정점 상태로 IN/OUT. η_e는 S4의 정확한 체적으로 확정한다.

### 7.4 S4 절단 적분

- 체적 모멘트 m_α = ∫_{cell∩Ω} P_α dV. P_α는 축별 Legendre 텐서곱, 차수 ≤ q. 기본 q = p+1, 질량 행렬 정확도가 필요하면 q = 2p.
- **재귀 발산정리**: F = e_x ∫_{x0}^{x} P_α ds로 두면 x = x0 면과 y/z 면의 기여가 0이다. 남는 항은 (a) 클리핑된 삼각형 위 적분(Dunavant, 정확)과 (b) x = x1 면 위 Ω 단면 다각형의 2D 모멘트다. (b)는 같은 방식으로 2D → 1D로 줄인다. tessellation 기하 기준으로 정확하다.
- **비음수 모멘트 피팅**:
  - 후보점: 셀의 (q+2)³ Gauss 점 중 Ω 안의 점 + 경계 근처 octree 잎 Gauss 점.
  - `scipy.optimize.nnls`로 비음수 가중치를 구한다. 비음수 가중치는 소성·대변형 Newton 수렴에 필요하다.
  - 가중치가 작은 점부터 제거하고 다시 피팅하는 과정을 점 수 ≤ (q+1)³이 될 때까지 반복한다.
  - 잔차가 1e-10·|cell|을 넘으면 fallback.
- **fallback** (비수밀, dirty 형상): octree 깊이 3–4, 잎 중심을 winding number로 분류. 내부 잎은 Gauss (p+1)³, 절단 잎은 2×2×2 Gauss + 점별 내부 판정. `issue` 코드 `QUAD_FALLBACK`.
- 경계 적분점: 셀∩삼각형 다각형 → fan 분할 → Dunavant 차수 2p → (x, w = dA, n, face_id).
- 검사: Σw = tessellation 체적(상대 1e-12), 경계 Σw = tessellation 면적, 닫힌 면의 Σ w·n = 0.
- Algoim(정확 형상)은 Phase 6 선택 과제. 쓰게 되면 parametric B-Rep → SDF 변환이 선결 과제다.

### 7.5 S5 small-cell extension

- **large 셀**: η_e ≥ θ (`small_cell_volume_ratio`, 기본 0.2, 범위 [0, 0.5)).
- 기저 분류:
  - **master**: 지지 안에 large 셀이 1개 이상
  - **extended (slave)**: 지지가 Ω와 겹치지만 large 셀이 없음
  - **inactive**: 지지가 Ω와 겹치지 않음
- slave j마다, 인덱스 공간에서 master 함수로만 구성된 (p+1)^d 블록 I(j) 중 가장 가까운 것을 고른다(블록 중심과 j의 거리 최소, 바깥 법선 방향 우선).
- 계수 e_ij = Π_ν ℓ_{i_ν}(j_ν) (균일 knot 인덱스에서 Lagrange 외삽, Höllig 2003). p차 다항식 재생과 분할단위성이 보존된다.
- 얇은 형상으로 (p+1) 블록을 못 찾으면: 차수를 낮춘 블록 → 최근접 master 상수 외삽(AgFEM형) 순서로 처리한다. `issue` 코드 `THIN_FEATURE`, 세분 권고.
- 다중재료: 지지 체적비는 **합집합 영역**으로 계산한다(재료 계면은 내부이므로 extension 대상이 아니다).
- ghost penalty(P12/P13 locking-free, C^{p−1}에서는 k = p 항만 필요)는 C++ 단계의 대안으로만 둔다.
- 검증: 분할단위성(Σ extended 기저 = 1), 선형·이차 재생, κ sweep(V-COND).

### 7.6 S6 penalty / Nitsche

- 기본 penalty (§5.4).
- β sweep {1e2 … 1e6}·E/h로 수렴 평탄 구간을 확인하고 ADR-002로 기본값을 확정한다.
- Nitsche는 선택 사항이다(§5.4).
- FCM α 기준선: OUT 영역에 가상 강성 α = 1e-8·E 요소를 둔다. Phase 1의 조건수 완화와 extension 비교용이다.

### 7.7 다중재료와 결합

- 멀티셀: 같은 NurbsVolume, 같은 노드. 재료 셀마다 적분점을 따로 만들고 Properties만 다르게 한다. 공유 스플라인 장이므로 계면은 완전 접합이다(Coreform과 동일).
- tie, coupling, RP: §5.4 BC-B/BC-C.
- 일반 FE 혼합(Phase 5): 같은 ModelPart에 일반 FE 요소를 두고 tie로 연결한다.

### 7.8 잠김 방지 (Phase 5, C++)

- U(p, C^{p−2}) / P(0, C^{−1}): 셀별 상수 압력, 요소 단위 정적 응축 K_uu − K_up K_pp⁻¹ K_pu.
- 소형 셀 압력은 이웃 large 셀과 aggregation(같은 압력 자유도 공유).
- 선택지: B-bar(투영 체적 변형률), F-bar, equal-order U2C1/P2C1 VMS(τ = c·h²/(2μ)).
- 출력 `PRESS_IGA`.

### 7.9 고유치와 동역학

- 고유치: `spectra_sym_g_eigs_shift`(shift σ, nev) 또는 FEAST. 질량 정규화를 확인한다.
- 동적 시간적분은 **Kratos scheme을 사용**한다 (사용자 결정 D5: Kratos 적분법을 HHT의 상위 호환으로 채택).
  - 기본: Bossak(WBZ-α), `damp_factor_m` α_B ∈ [−1/3, 0], 기본 −0.3. β = (1−α_B)²/4, γ = 1/2 − α_B로 2차 정확도와 무조건 안정을 유지하면서 고주파 감쇠를 조절한다.
  - 감쇠 없는 비교가 필요하면 Newmark(α_B = 0), 강한 고주파 감쇠가 필요하면 BDF2.
  - HHT는 구현하지 않는다. Coreform 검증(V10)은 HHT 결과 자체가 아니라 물리량 기준(f1, 위상, 응답 이력 RMS)으로 비교한다.
- 초기 가속도 계산, amplitude, restart(상태 저장)는 Coreform 동작을 따른다.

### 7.10 응답 스펙트럼 / 랜덤 응답

- 참여계수 Γ_ik = φ_iᵀ f_k (f_k = ∫ ρ N e_k), 유효질량 Γ².
- 모드 응답 R_i = Γ_i S_a(ω_i, ζ)/ω_i² · φ_i. 모드 응력은 φ_i를 노드에 넣고 `CalculateOnIntegrationPoints`로 얻는다.
- 모드 조합: ABS, SRSS, NRL, TENP, CQC(Der Kiureghian ρ_ij), GRP(10% 근접 그룹), DSC(`duration` 필요).
- 방향 조합: ALGEBRAIC(모드 조합 전), SRSS, R40, R30(ASCE 4-98).
- 스펙트럼 보간: 주파수 log-log, 감쇠 곡선 사이는 감쇠 선형·크기 로그, 표 밖은 고정. S_D = S_V/ω = S_A/ω².
- 랜덤: H_i(ω) = 1/(ω_i² − ω² + 2iζω_iω), 응답 PSD = Σ_ij Γ_i Γ_j H_i H_j* S_in. RMS는 사다리꼴 적분(구간당 20점, bias 3, LINEAR/LOG).
- Segalman RMS von Mises = sqrt(Σ V_αβ T_αβ). 입력 단위 BASE(g²/Hz), FORCE, DB(대역 1–15).
- 범위 밖(Coreform도 없음): missing mass, Gupta, Lindley-Yow.

### 7.11 열

- 전도: `LaplacianElement`(CONDUCTIVITY), 열용량은 ρc 질량 행렬, 시간적분 BDF2.
- film q = h(T − T∞), radiation εσ(T⁴ − T∞⁴): 표면층 위 `ThermalFace3D3N`, Newton.
- 지정 온도: `SupportLaplacianCondition` penalty.
- NAFEMS T1용 온도장: 제어점 L2 투영(질량 행렬, 참조 조립기). 투영 오차를 따로 보고한다.

### 7.12 후처리와 시각화

- 장 평가: 임의 점 x에서 셀 찾기(항등 사상이라 인덱스 계산만) → 기저 평가 → 제어점 계수와 곱 (@njit).
- 응력: 적분점 값 직접 사용(probe = 가중치 0 적분점, Coreform 방식). 비선형 상태 변수도 정확히 얻는다.
- proxy 격자: IN 셀은 hex(서브디비전 선택), CUT 셀은 셀 내 클리핑 다면체를 tet으로 분할. 값은 장 평가로 채운다.
- 형식: VTKHDF(h5py, 시계열) 또는 VTU(meshio). ParaView에서 바로 연다.
- L2/H1/에너지 노름: quadrature DB의 X, dV, U, S로 계산한다.

---

## 8. 사용자 Python API

Coreform API의 용어와 구조를 따르되, 백엔드는 Kratos다.

```python
import wht_iga as wi

model = wi.Model("roark")
part = model.part_from_step("plate_hole.step", heal=True)
model.material("steel", E=2.0e5, nu=0.3, density=7.85e-9).assign(part)

mesh = wi.ImmersedRectilinearMesh(
    degree=2, continuity=1,
    sizing=wi.RectilinearSizing(element_size=2.5),
    tessellation_size=0.25,
    bc_penalty_value=1e10,
    small_cell_volume_ratio=0.2,
    padding=(2, 2, 2),
)
mesh.refine(wi.Refinement("hole_edge", level=2,
                          region=wi.CadEntitiesRegion(sets=["hole"]), balance="two_to_one"))

step = model.static_step("load")
step.bc(wi.Displacement(set="sym_x", ux=0.0))
step.bc(wi.Displacement(set="sym_y", uy=0.0))
step.load(wi.Pressure(set="right", magnitude=-1.0))
model.probe(wi.PointProbe("peak", location=(10, 0, 0), variables=("S", "U")))

job = wi.Job("roark", model=model, mesh=mesh, num_threads=16)
job.build_mesh()            # S1–S6 → roark.wiga.sqlite
job.submit()                # S7–S9
job.wait_for_completion()
peak = job.results.probe("peak")["S"]

spec = job.to_json()
wi.Job.from_json(spec).validate()      # 실패 시 ValidationError(issues=[ValidationIssue(path, code, message)])
```

| Coreform API | WHT_IGA |
|---|---|
| `ImmersedRectilinearMesh(degree, continuity, bcPenaltyValue, tessellationSize, elementSize/IntervalCount/ByAxis, elementLayout, smallCellVolumeRatio, padding, frame…)` | `wi.ImmersedRectilinearMesh(...)` (snake_case) |
| `BodyFitMesh` | `wi.BodyFitMesh` (정렬 블록 한정) |
| `mesh.Refinement`, `CadEntitiesRegion`, `SphereRegion`, `BoxRegion`, `refinementBalance` | `wi.Refinement`, 같은 region 3종, `balance` |
| `PointProbe`, `LineProbe(numEvalPoints)`, `MultiPointProbe` | 같은 이름 |
| `Job(buildMesh, writeInput, submit, waitForCompletion, kill, getRuntimeStatus)` | `wi.Job(build_mesh, write_input, submit, wait_for_completion, kill, status)` |
| portable spec `to_dict/to_json/from_json/copy_with/validate` | 같은 메서드 |

---

## 9. 단계별 로드맵

기간은 1인 기준 추정이다. 각 Phase는 게이트를 통과해야 다음으로 넘어간다.

### Phase 0 — 환경과 Kratos API spike (2주)

- **산출물**
  - micromamba `whtiga` env, `pyproject.toml` + `uv.lock`, `check_env.py`, `docs/ENVIRONMENT.md`, `docs/KRATOS_BUILD.md`
  - `WHT_IGA` 골격 (bootstrap, isolation test)
  - S01–S18 spike 시험과 리포트, ADR-001 (R1/R2/R3)
  - 참조 조립기 최소판 (정육면체 untrimmed 선형 탄성)
- **게이트 G0**
  - S01–S05 결론 확정
  - untrimmed 큐브 p = 1, 2, 3에서 Kratos 경로와 참조 조립기 변위 차이 ≤ 1e-10
- **위험**: S05 성능 → R2

### Phase 1 — 선형 정적, trimmed, penalty, extension 없음 (5주)

- **산출물**: S1–S4, S6, S7(BC-A/B, LoadCondition, 중력), 선형 전략, mesh/results/quadrature DB, VTU 최소 출력, FCM α 옵션, spec/validate
- **게이트 G1**
  - **V01** 지정 변위 큐브 u = (x, 2y, 3z): 모든 p ∈ {1,2,3}, C에서 최대 상대오차 ≤ 1e-5 (Coreform 수준)
  - **V02** Roark: p = 2 균일 최세 격자 peak 상대오차 ≤ 1% (Coreform 0.54%), 관측 수렴률 ≥ 1.8 (Coreform 1.86/2.23/2.38). 기준 4.3175e-6 (`roark-17p1-7a-accepted-peak.csv`)
  - **V03** 보 휨(moment 제외 5종): 끝단 처짐 Euler–Bernoulli 대비 5% 이내
  - 모멘트 피팅 다항식 정확도 1e-12
- **위험**: 소형 셀 조건수(Phase 2 전까지 FCM α로 완화), tessellation 기하 포화

### Phase 2 — small-cell extension과 조건수 (3주)

- **산출물**: S5 extension, `operators/compose`, 구속식 배치 생성, 조건수 도구, THIN_FEATURE 처리
- **게이트 G2**
  - **V-COND**(자체): 구멍 중심을 이동시켜 η_min ∈ [1e-8, 1]을 만든다. κ 변동 ≤ 10배(평탄), h에 대한 log-log 기울기 −2 ± 0.3. 미안정화 대비 개선을 그래프로 제시
  - **V06** 압력 구 octant (a = 0.9, b = 1, Lamé): `pressurized-sphere-qoi.csv`의 p별·level별 허용오차 열 전부 통과. 목표 p = 2 최세 u L2 ≤ 1e-4, σrr ≤ 0.5% (Coreform 0.00365% / 0.211%)
  - V01, V02를 extension 켠 상태로 재통과
  - extension 분할단위성·재생 단위시험
- **위험**: 얇은 형상에서 블록 부재, 구속식 성능

### Phase 3 — 다중재료, tie, coupling, probe, 세분화, 출력 (6주)

- **산출물**: 멀티셀 공유 장, tie(BC-C), RP/분배/kinematic coupling, amplitude, 다단계, 가중치 0 probe, THB 세분화, VTKHDF proxy 격자, probe JSON, L2 노름 도구
- **게이트 G3**
  - **V03** moment 케이스 5% 이내
  - **V04** bimetallic Cu/W: merge와 tie의 끝단 처짐 상대차 ≤ 1e-6. 변환단면 해 대비 오차는 `bimetallic-beam-tie-comparison.csv` 수준(p = 3 ≈ 0.54%) ± 0.2%p
  - **V05** multiplanar 계면 7단면(84 해석): monolithic 대비 차이가 Coreform 보고 수준(p = 2 ≤ 3.3e-6) 이하
  - **V02 국소 세분**: L3에서 peak 오차 ≤ 0.2%, DOF ≤ 80k (Coreform 0.105% @ 57,300 DOF, `roark-17p1-7a-refinement.csv`)
  - **E01** 박판 구멍(IPS): Roark Kt = 2.23568 기준 σ_max = 3.7261e-6 psi 대비 probe 오차 ≤ Coreform(−8.7%) 수준, 세분 시 개선 추세
  - **E02** 적응 세분 예제 재현

### Phase 4 — 고유치, 암시적 동역학(Kratos Bossak), 스펙트럼/랜덤, 열, 열변형 (6주)

- **게이트 G4**
  - **V09** 외팔보 충격·진동: f1 = 80.8 Hz 대비 1% 이내(Coreform 81.8), 랜덤 RMS 오차 ≤ 1.79% + 0.5%p, 스펙트럼 조합 규칙별 결과를 수기 계산과 1e-10까지 일치
  - **V10** 진동 보 (Kratos Bossak): f1 = 25,000 기준. 문서의 판정 기준 26개 중 물리량 기준(f1 오차, 위상, 응답 이력 RMS)을 전부 통과. HHT 고유 파라미터에 묶인 기준은 Bossak 대응값으로 바꾸고 차이를 리포트에 기록
  - **V11** NAFEMS T1: σyy(D) = 50 대비 4% 이내 (Coreform 2.99%)
  - **V12** NAFEMS T2 927 K, T3 36.60 °C @ 32 s, T4 18.3 °C: 문서 한계(1% / 2% / 3%) 이내. 2D 벤치마크는 얇은 3D 슬래브 + 대칭 penalty로 구성

### Phase 5 — 비선형, 잠김 방지, 접촉 (8주, C++ 승격 포함)

- **산출물**: NLGEOM(TL/UL), 초탄성, J2(T 의존 표), `WhtIgaApplication` MixedUP0/VMS 요소, 접촉(표면층 + ContactSM 또는 C++)
- **게이트 G5**
  - **V07** Kirsch ν = 0.49999 혼합: 문서 한계 (e_u, e_σxx) ≤ (0.01, 0.2), (0.005, 0.1), (0.002, 0.05), 수렴률 ≥ 2.3 / 2.2 (Coreform 2.79 / 2.38)
  - **E04** Cook 막: u_y = 6.98946 mm 대비 1% 이내
  - **V08** Hertz: F_N = 481,579 lbf 대비 4% 이내(Coreform 1.84–3.74%), 미끄럼 F_T RMS ≤ 4.5%
  - **E05** door-lock: RF3 = −46.946 kN 대비 5% 이내
  - **E03** glass-to-metal seal: 13셀, 3 NLGEOM 열 step 수렴. 증분 수가 Coreform(각 step 15/15, cutback 없음)의 2배 이내

### Phase 6 — 성능과 확장 (지속)

- 프로파일 기반으로 Numba 핫스팟을 C++(pybind11)로 이전 (§2.2 기준).
- C++ `ExtendedSplineSolidElement`로 구속식 제거.
- 반복 솔버: AMGCL, IgaApp `AdditiveSchwarzPreconditioner` (spike 필요).
- 전처리 병렬화(적분점 수 가중 분할). MPI는 Kratos MPI 빌드 여부를 확인한 뒤 판단.
- 선택: Nitsche/ghost penalty, Algoim 정확 형상, H3D.
- **게이트**: 1M DOF 3D 브래킷 end-to-end < 30분(16코어), 메모리 < 64 GB. 참고로 Coreform glass-to-metal seal은 16 CPU에서 mesh 318.8 s + solve 691.0 s.

---

## 10. 검증·시험 계획

### 10.1 Coreform 검증 매뉴얼 매핑

CSV 경로 접두어: `user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/`

| ID | Coreform 케이스 | 기준값 / CSV | sweep | 허용오차 | Phase |
|---|---|---|---|---|---|
| V01 | prescribed-displacement-cube (E=1, ν=0, 격자 간격 0.6) | 해석해 u = (x,2y,3z), σ = diag(1,2,3) | p 1–3 × C 0..p−1, 회전 Frame | 상대오차 ≤ 1e-5 (Coreform 최대 변위 1.5e-5, 응력 5.5e-5) | 1 |
| V02 | roark-17p1-7a (L=200, D=40, r=10, t=1, P=1e-6) | `roark-17p1-7a/` `-accepted-peak.csv`(4.3175e-6), `-qoi.csv`, `-refinement.csv` | p 1–3, h = 6/2^k (k=0–4), 국소 L1–L3 | peak ≤ 1%, 수렴률 ≥ 1.8, 국소 ≤ 0.2% | 1/3 |
| V03 | beam-bending (b=h=0.1, L=1, E=200 GPa) | `beam-bending/` `-accepted.csv`, `*-tip-deflection.csv`, `-traction-variants.csv` | 6개 하중, p 1–3, h = (0.1,0.1,0.25)/2^k | 5% (느린 케이스 10%) | 1/3 |
| V04 | bimetallic-beam-bending (Cu/W, b=10, L=100) | `bimetallic-beam-bending/` `*-accepted.csv`, `-tie-comparison.csv` | merge/tie × p × 하중 4종 | 5%, merge–tie ≤ 1e-6 | 3 |
| V05 | multiplanar-material-interfaces (7단면, 84 해석) | `multiplanar-material-interfaces/*.csv` | p 1–3 × h 4단계 | monolithic 대비 Coreform 수준 | 3 |
| V06 | pressurized-sphere (a=0.9, b=1, p_i=1e-6) | `pressurized-sphere/` `-qoi.csv`(허용오차 열) | p 1–3 × h = 0.1/2^k | CSV 허용오차 | 2 |
| V07 | infinite-plate-with-hole (ν=0.49999, T=10) | 문서 오차·한계·수렴률 | p=2 h = 1/2, 1/4, 1/8 | 문서 한계 | 5 |
| V08 | contact-indenter (Hertz, μ=0.3) | F_N 481,579 lbf, F_T 144,474 lbf, p0 4.599e6 psi | h = 0.5, 0.25, 0.125 | F_N ≤ 4% | 5 |
| V09 | cantilever-shock-vibration | f1 80.8 Hz, 랜덤 RMS u 340 µm, shock peak 965 µm | 모드 20개, ζ=0.02 | 1% / RMS +0.5%p | 4 |
| V10 | vibrating-beam (Coreform은 HHT, WHT_IGA는 Kratos Bossak; L=10, 2×2, E=1e6, ρ=1e-6) | f1 = 25,000, 26개 기준(물리량 기준) | Δt 1e-6–2e-6, α_B sweep | 문서 기준 | 4 |
| V11 | nafems-p18-t1 (열변형) | σyy(D) = 50.0 | p 1–3 | 4% | 4 |
| V12 | nafems-p18-t2/t3/t4 (열전달) | 927 K / 36.60 °C @ 32 s / 18.3 °C | h, Δt | 1% / 2% / 3% | 4 |
| V-COND | (자체) 조건수 | η sweep | η ∈ [1e-8, 1], h 4단계 | κ 평탄(≤ 10배), 기울기 −2 ± 0.3 | 2 |

### 10.2 Coreform 예제 매핑

| ID | 예제 | 기준 | 허용오차 | Phase |
|---|---|---|---|---|
| E01 | thin-plate-with-hole (IPS, 10×5×0.1, r=1) | σ_max = 3.7261e-6 psi (Kt 2.23568), Coreform probe 3.4031e-6 | Coreform 수준 이상, 세분 시 개선 | 3 |
| E02 | plate-with-hole-adaptivity | 문서 정성 결과 | 추세 일치 | 3 |
| E03 | glass-to-metal-seal (13셀, T 의존 소성, 3 열 step) | `generated-downloads/example-problems/glass-to-metal-seal/*.csv`, 재료 TSV/CSV(`example-problems/glass-to-metal-seal/data/`) | 수렴 + 정성 | 5 |
| E04 | Cook's membrane (ν=0.49999) | u_y = 6.98946 mm | 1% | 5 |
| E05 | door-lock-assembly (3 부품, μ=0.3) | RF3 = −46.946 kN | 5% | 5 |

주의: `user-guide/downloads/`의 CAD·입력 파일은 비어 있다. 모든 형상은 `geometry/primitives.py`(OCP)로 문서 치수에서 다시 만들고, 치수와 하중은 `tests/verification/data/*.json`에 출처(문서 파일·줄)와 함께 기록한다. Abaqus 참조값보다 해석해를 우선한다.

### 10.3 자체 단위시험 (`tests/unit`, Kratos 불필요)

- B-spline 분할단위성, 미분의 유한차분 일치, Kratos knot 규약 왕복 변환
- Bézier 추출: Cᵉ·B = N (1e-14). THB 분할단위성·선형독립·truncation 성질
- 분류: scanline과 winding 결과 일치(구, 토러스). 비수밀 메쉬에서 fallback 동작
- 절단 셀 적분: 평면·구 절단에서 차수 ≤ q 단항식 적분 정확도 1e-12, 가중치 비음수
- 경계 적분점: 면적 일치, 닫힌 면 Σ w·n = 0
- extension: 분할단위성, p차 재생, 계수 크기 상한, θ 변화에 대한 연속성
- penalty: β sweep에서 O(1/β) 거동과 평탄 구간
- SQLite 왕복 무손실, spec JSON 왕복, validate issue 코드
- 격리 시험(`wht_topopt`/`aieng` 미사용)

### 10.4 Kratos 교차검증 (`marker: kratos`)

- 모든 V 케이스의 축소판에서 Kratos 경로와 참조 조립기 결과 차이 ≤ 1e-10.

### 10.5 SBM 비교 트랙

- IgaApp `SbmSolidCondition`/`SbmLoadSolidCondition` + `NurbsGeometryModelerSbm`으로 V02, V06을 푼다.
- 같은 DOF 수에서 cut + extension 경로, FCM α 경로와 정확도·조건수를 비교한다.
- 주 경로 의사결정에는 쓰지 않는다. 결과는 `docs/verification/sbm_comparison.md`.

### 10.6 보고 형식

- 케이스별 `docs/verification/<ID>.md` 자동 생성: 입력 요약, 기하 오차, 수렴 표·그래프, 기준값 대비 오차, 통과 여부, 실행 시간, DOF 수.

---

## 11. 위험 등록부

| # | 위험 | 영향 | 대응 | 게이트 |
|---|---|---|---|---|
| R1 | Python으로 적분점 기하를 직접 만들 수 없음 (확인됨) | 높음 | 셀 단위 Option 1 + 구속식 경로(R1). C++는 Phase 5 이후 | G0 |
| R2 | `LinearMasterSlaveConstraint` 대량 생성·조립 성능, DofSet 누락 | 높음 | 배치 생성, R2 추출 경로, C++ C̃ 요소 | S05 |
| R3 | 다점 적분점 기하와 SM 요소 비호환 | 높음 | 점당 요소, R3 참조 조립기 | S04 |
| R4 | Python 루프 모델 구성 시간 (1e5–1e6 요소) | 중 | Numba 사전 계산 + 배치 호출, 셀 병합, Phase 6 C++ 생성기 | S16 |
| R5 | dirty CAD에서 적분 실패 (비수밀, 자기교차) | 높음 | ShapeFix + Sewing, 3축 투표, winding fallback, octree fallback, `validate` 차단 | G1 |
| R6 | 멀티셀 공유 장 구성 (계면 삼각형 불일치) | 중 | general fuse 공유 면, 셀 체적 합 = 합집합 체적 검사(1e-10) | G3 |
| R7 | 얇은 형상에서 extension 블록 부재, 계수 폭증 | 중 | 차수 저하 → 상수 외삽 → 세분 권고, 계수 상한 감시 | G2 |
| R8 | penalty 값 민감도 (정확도 대 조건수) | 중 | β sweep ADR, 자동 스케일, Nitsche 선택(C++) | G1 |
| R9 | 잠김 방지 혼합 요소는 C++ 필수 | 높음 | Phase 5 초반 C++ 착수, B-bar 대안 | S13 |
| R10 | 접촉: 표면층과 ContactSM 호환성 | 높음 | 독립 자유도 스프링 표면층, 실패 시 C++ penalty NTS | S14 |
| R11 | Coreform 동역학 검증값이 HHT 기준 | 낮음 | Kratos Bossak 사용(D5 확정). 물리량 기준으로 비교, HHT 파라미터 의존 기준은 대응값으로 변환해 기록 | G4 |
| R12 | 라이선스·IP (U-spline 특허, Gmsh GPL, QuESo BSD-4) | 중 | U-spline 제외, Gmsh 외부 프로세스, QuESo 미사용 | 상시 |
| R13 | 단일 소스빌드 Kratos 의존, **빌드 스크립트 소실** (`kratos-deps\build_kratos_full.bat` 1바이트) | 중 | `CMakeCache.txt`에서 스크립트 복원, `docs/KRATOS_BUILD.md` 작성(§12.4). 핵심 경로는 core + SM 위주로 설계 | Phase 0 |
| R14 | Python ABI 불일치 (Kratos는 vdmc Python 3.12로 빌드) | 중 | `whtiga` env Python 3.12.x 고정, S01에서 확인, 실패 시 vdmc 런타임 | S01 |
| R15 | 검증 입력 파일 부재로 형상·하중 재현 오차 | 중 | 문서 수치로 재작성, 불확실 항목 기록, 해석해 우선 | 각 V |

---

## 12. 환경과 재현성

### 12.1 env 구성

- **micromamba** 단일 실행 파일(`%LOCALAPPDATA%\micromamba\micromamba.exe`), root prefix `D:\mamba`(`MAMBA_ROOT_PREFIX`).
- `~/.mambarc`: channels [conda-forge], `channel_priority: strict`, `envs_dirs: [D:\mamba\envs, C:\Users\GOODMAN\miniconda3\envs]` (기존 env도 `micromamba activate vdmc`로 접근).
- `whtiga` env: `micromamba create -n whtiga -c conda-forge python=3.12 mkl hdf5` + 필요한 네이티브 라이브러리.
- Python 패키지: **uv pip**로 설치하고 `pyproject.toml` + `uv.lock`으로 고정.
  ```
  uv pip install --python D:\mamba\envs\whtiga\python.exe -e WHT_IGA[cad,viz,solve,test]
  ```
- 셸 프로파일: PowerShell 7/5, cmd AutoRun, Git Bash의 conda 훅을 micromamba 훅으로 교체(백업 후). miniconda 설치본과 PATH는 남긴다(aieng CalculiX 탐색 유지).

### 12.2 Kratos 연결

- venv에 pip Kratos wheel을 **설치하지 않는다**.
- `site-packages/kratos_ebuild.pth`: `E:\kratos-build\install`, `E:\kratos-build\install\libs`.
- DLL 디렉토리(`E:\kratos-build\install\libs`, `D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin`)는 `_bootstrap.py`가 `os.add_dll_directory`로 등록한다.
- import 후 `KratosMultiphysics.__file__`이 `E:\kratos-build\install` 아래인지 검증한다.

### 12.3 검증 스크립트

- `micromamba run -n whtiga python WHT_IGA/scripts/check_env.py`
  - Python 3.12.x 확인
  - Kratos 경로와 버전(10.4.3) 확인
  - IgaApplication, StructuralMechanics, LinearSolvers, ConstitutiveLaws, ConvectionDiffusion import
  - `HasMKL()` 보고
  - numba `@njit` 함수 컴파일·실행

### 12.4 Kratos 빌드 재현

- 사실: 로컬 소스빌드(Release, MSVC `cl`, `PYTHON_EXECUTABLE = C:\Users\GOODMAN\miniconda3\envs\vdmc\python.exe`, prefix `E:/kratos-build/install`, `KratosCore.dll` 2026-09-28). 소스 `D:\PythonCodeStudy\kratos-src` v10.4.3(96bad120) + 로컬 수정 약 294개 파일(IgaApplication은 수정 없음).
- 문제: 빌드 스크립트 `D:\PythonCodeStudy\kratos-deps\build_kratos_full.bat`가 1바이트(내용 소실).
- 조치(Phase 0):
  1. `E:\kratos-build\Release\CMakeCache.txt` 백업
  2. 애플리케이션 목록과 옵션을 추출해 빌드 스크립트 복원
  3. `docs/KRATOS_BUILD.md`에 재현 절차와 `git -C kratos-src diff --stat` 요약 기록
  4. 재빌드는 하지 않는다(필요해지는 시점은 C++ 승격 때)
- 주의: 병렬 빌드 `-j 6` 상한(32 GB RAM에서 그 이상은 `cl.exe` 충돌 기록).

### 12.5 다른 에이전트와 공유

- `WHT_IGA/docs/ENVIRONMENT.md`: 환경 단일 원본.
- `WHT_IGA/AGENTS.md`: 강제 규칙 (`micromamba run -n whtiga`, conda install 금지, uv pip만, Python 3.12, `wht_topopt` import 금지, pip Kratos wheel 금지) + ENVIRONMENT.md 링크.
- `WHT_IGA/CLAUDE.md`: `@AGENTS.md`.
- 루트 `AGENTS.md` §4 Workspace Layout 표에 `WHT_IGA/` 행 추가.
- `~/.claude/CLAUDE.md`: "conda 사용 금지 — micromamba(root D:\mamba) + uv pip" 한 줄.

---

## 13. 참고문헌·라이선스 레지스트리

### 13.1 핵심 문헌

1. Borden, Scott, Evans, Hughes (2011). Isogeometric finite element data structures based on Bézier extraction of NURBS. IJNME 87:15–47. DOI 10.1002/nme.2968 — Bézier 추출 (§7.2)
2. Höllig, Reif, Wipper (2001); Höllig (2003). Weighted extended B-spline approximation; Finite Element Methods with B-Splines. SIAM — extension 계수 (§7.5)
3. Marussig, Hughes (2018). A review of trimming in isogeometric analysis. ACME 25:1059–1127. DOI 10.1007/s11831-017-9220-9 — 지지 기준 안정화
4. Badia, Verdugo, Martín (2018). The aggregated unfitted finite element method for elliptic problems. CMAME 336:533–553. DOI 10.1016/j.cma.2018.03.022
5. Badia, Neiva, Verdugo (2022). Linking ghost penalty and aggregated unfitted methods. CMAME 388:114232. DOI 10.1016/j.cma.2021.114232
6. Burman, Hansbo, Larson (2023). On the design of locking free ghost penalty stabilization. Comput Mech 71:1007–1019. DOI 10.1007/s00466-022-02263-4
7. de Prenter, Verhoosel, van Brummelen, Badia (2023). Stability and conditioning of immersed finite element methods. ACME 30:3617–3656. DOI 10.1007/s11831-023-09913-0 — κ 거동 (V-COND)
8. Giannelli, Jüttler, Speleers (2012). THB-splines: the truncated basis for hierarchical splines. CAGD 29:485–498 — THB (§7.2)
9. Schillinger, Ruess (2015). The Finite Cell Method: a review. ACME 22:391–455. DOI 10.1007/s11831-014-9115-y — FCM α
10. Schillinger et al. (2012). An isogeometric design-through-analysis methodology… CMAME 249–252:116–150. DOI 10.1016/j.cma.2012.03.017
11. Meßmer, Teschemacher, Leidinger, Wüchner, Bletzinger (2022). Efficient CAD-integrated isogeometric analysis of trimmed solids. CMAME 400:115584 — QuESo + Kratos 선례
12. Loibl et al. (2026). 절단 적분 기법 비교 (`papers/P05`) — 다면체 입력에는 octree/모멘트 피팅
13. Saye (2015, 2022). High-order quadrature on implicitly defined domains. SISC 37:A993; JCP 448:110720 — Algoim (선택)
14. Antolin, Buffa, Hirschler (2022). Quadrature-free immersed isogeometric analysis. Eng Comput 38:4475–4499 — 선택
15. Jacobson, Kavan, Sorkine-Hornung (2013). Robust inside-outside segmentation using generalized winding numbers. ACM TOG 32(4) — §7.3
16. Hilber, Hughes, Taylor (1977). Improved numerical dissipation for time integration algorithms in structural dynamics. EESD 5:283–292 — HHT
17. Der Kiureghian (1981). A response spectrum method for random vibration analysis of MDF systems. EESD 9:419–435 — CQC
18. Segalman, Fulcher, Reese, Field (2000). An efficient method for calculating RMS von Mises stress in a random vibration environment. JSV 230:393–410
19. Thomas et al. (2022). U-splines: splines over unstructured meshes — **참고만, 범위 밖(특허)**
20. Coreform IGA for Abaqus 문서 (2026.11): Verification Manual, Example Problems, Python API — `user-guide/`

### 13.2 소프트웨어·라이선스

| 구성요소 | 라이선스 | 사용 방식 |
|---|---|---|
| Kratos Multiphysics | BSD-4-Clause | 소스빌드 바인딩, (선택) C++ 외부 앱 링크 |
| OCP (cadquery-ocp) / OpenCASCADE | LGPL-2.1 + OCCT 예외 | 동적 import |
| numpy, scipy, numba, llvmlite | BSD | 필수 의존 |
| pypardiso / Intel MKL | BSD-3 / Intel Simplified Software License | 선택 |
| pyvista, meshio, h5py | MIT / BSD | 선택 |
| trimesh | MIT | 선택 (STL 입력) |
| Gmsh | GPL-2.0+ | **외부 프로세스로만** (선택) |
| QuESo | BSD-4-Clause | **미사용**, 비교 참고 |
| Algoim | BSD-3 | Phase 6 선택 |
| pythonOCC | LGPL-3.0 | **미사용** (OCP로 대체) |
| U-splines (Coreform) | 특허 | **범위 밖** |
| Altair hvtrans | 상용 | 선택 플러그인 |

---

## 14. 열린 결정 사항

| # | 결정 | 선택지 | 기본안(결정 전까지 적용) |
|---|---|---|---|
| D1 | R2/R3 전환 방식 | spike 결과에 따라 자동 / 사용자 확인 후 | 사용자 확인 후 (ADR 제출) |
| D2 | C++ 승격 시점 | Phase 5 초반 / Phase 2 직후 조기 착수 | Phase 5 초반 |
| D3 | penalty 기본값 | 사용자 필수 입력(Coreform식) / 자동 스케일 | 자동 스케일 + 사용자 지정 시 우선 |
| D4 | IgaApplication 의존 범위 | IgaApp 조건 사용 / SM 객체만 쓰는 "wheel-safe 모드" 기본 | IgaApp 사용 (소스빌드 고정이므로) |
| D5 | 동적 시간적분 | **확정 (2026-11-14)**: Kratos scheme(Bossak 기본, Newmark/BDF 선택). HHT 구현 안 함 | – |
| D6 | 2D 모드 | 3D 전용(얇은 슬래브) / 별도 2D 모드 | 3D 전용 |
| D7 | Phase 6 성능 목표 | 1M DOF / 30분 / 64 GB 유지 / 조정 | 유지 |
| D8 | aieng 연동(runner) 시점 | Phase 3 / Phase 6 | Phase 6 |
| D9 | 루트 `AGENTS.md` §9 `conda activate aieng311` 수정 | `micromamba activate aieng311`로 교체 / 유지 | 사용자 확인 후 교체 |
| D10 | Gmsh 선택 기능 포함 | 포함(외부 프로세스) / 제외 | 제외 (OCP BRepMesh로 충분) |

---

## 부록 A. Kratos API 근거 (`D:\PythonCodeStudy\kratos-src`, v10.4.3)

| # | 내용 | 위치 |
|---|---|---|
| A-1 | `NurbsVolumeGeometry` Python 클래스, 생성자, 차수/knot 조회 | `kratos/python/add_geometries_to_python.cpp` l.313–330 |
| A-2 | `CreateQuadraturePointGeometries(res, nderiv)` / `(res, nderiv, [[u,v,w,weight]])` | 같은 파일 l.166–181 |
| A-2b | Option 1(GAUSS 기본): 점 전부를 담은 기하 1개, `nderiv == 2` 필수, 같은 요소 안의 점 전제 | `kratos/geometries/nurbs_volume_geometry.h` l.820–840 |
| A-3 | `SurfaceInNurbsVolumeGeometry(volume, surface)` | `add_geometries_to_python.cpp` l.504–505, `kratos/geometries/surface_in_nurbs_volume_geometry.h` |
| A-4 | `ModelPart.CreateNewMasterSlaveConstraint` 오버로드 | `kratos/python/add_model_part_to_python.cpp` l.540–584 |
| A-5 | `LinearMasterSlaveConstraint` Python 사용 예 | `kratos/tests/test_constraint_restart.py`, `kratos/python/add_constraint_to_python.cpp` |
| A-6 | 3D 솔리드 IGA 패턴 (NurbsVolume + SM `UpdatedLagrangianElement3D4N`) | `applications/IgaApplication/tests/test_nurbs_volume_element.py` |
| A-7 | `SupportPenaltyCondition`: 적분 = `PENALTY_FACTOR × weight × detJ`, 3성분 동시 | `applications/IgaApplication/custom_conditions/support_penalty_condition.cpp` l.32–67 |
| A-8 | IgaApp 조건 목록 (Support/Coupling Penalty·Lagrange·Nitsche, Load, LoadSolid, Sbm*) | `applications/IgaApplication/custom_conditions/` |
| A-9 | `IntegrationInfo::CUSTOM` ("integration points coming from QuESo") | `kratos/integration/integration_info.h` l.65 |
| A-10 | IgaApp 프로세스 (AssignIntegrationPointsToBackgroundElements, MapNurbsVolumeResults…, NitscheStabilization…, OutputQuadratureDomain, OutputEigenValues, Snake SBM) | `applications/IgaApplication/custom_processes/` |
| A-11 | SM `SpringDamperElement3D`, `NodalConcentratedElement3D1N` 등록 | `applications/StructuralMechanicsApplication/structural_mechanics_application.cpp` l.757, l.630 |
| A-12 | ConstitutiveLaws `ThermalElasticIsotropic3D` 등록 | `applications/ConstitutiveLawsApplication/constitutive_laws_application.cpp` l.404 |
| A-13 | 설치본 | `E:\kratos-build\install` (`libs\KratosIgaApplication.pyd`, `KratosIgaCore.dll/.lib`) |

## 부록 B. 용어

| 용어 | 뜻 |
|---|---|
| 배경 격자 | CAD를 감싸는 직교 B-스플라인 격자 |
| 셀 상태 OUT/IN/CUT | 셀이 형상 밖 / 완전히 안 / 경계에 걸림 |
| η (체적비) | 셀 안의 형상 체적 ÷ 셀 체적 |
| extension | 지지가 너무 작은 기저함수를 이웃 기저들의 선형 결합으로 대체하는 것 |
| master / slave | 구속식에서 독립 자유도 / 종속 자유도 |
| 적분점 기하 | Kratos `QuadraturePointGeometry`. 적분점 묶음과 그 점에서의 형상함수를 가진 기하 |
| 표면층 | CAD 면을 덮는 삼각망. 그 절점이 스플라인 장에 penalty로 결합됨 (Coreform SFM3D3 대응) |
| FCM α | 형상 밖 영역에 아주 작은 강성을 주는 Finite Cell Method 안정화 |
| tessellation | CAD 곡면을 삼각형으로 근사하는 것 |
| 항등 사상 트릭 | knot을 물리 좌표로, 제어점을 Greville 점에 두어 매개변수 좌표 = 물리 좌표(detJ = 1)로 만드는 기법 |
