# WHT_IGA — Kratos 기반 잠김 등기하해석(Immersed IGA) 솔버 마스터플랜 v2.1

- 문서 버전: **v2.1** (2026-11-14). v2.0을 대체하는 **단일 기준 문서**. 이 문서 하나로 작업한다.
- 작업 스트림 / 브랜치: `WHTIGA`
- 코드 위치: 저장소 루트 `WHT_IGA/`
- 통합한 문서 (원본은 보존, 수정하지 않음): `master-plan_v1.md`, `master-plan_v2.md`, `master-plan_v2-add_1.md`(갭 분석), `master-plan_v2-add_2.md`(add_1 검토·워크플로·메쉬 설계), `WHT_IGA/docs/spikes/TRACK_A.md` + `TRACK_A_REVIEW.md`(Phase 0 트랙 A 실측)
- 근거 자료: `나의조사.md`, `문헌종합.md`, `관련opensourceprj.md`, `papers/P01–P30`, `papers/implementation_synthesis.md`, `poc_codes/`, `user-guide/` (Coreform IGA 공식 문서 크롤 111쪽 + 검증 CSV), Kratos 소스 `D:\PythonCodeStudy\kratos-src` (v10.4.3)
- 상태: **계획 단계. Phase 0(환경·spike) 진행 중. 이 문서가 확정될 때까지 Phase 1 솔버 구현은 시작하지 않는다.**

## 목차

0. 개정 이력 (v1 → v2, v2 → v2.1 은 §0.5)
1. 비전과 범위
2. 설계 원칙
3. 시스템 아키텍처
4. 데이터 버스
5. Kratos 연동 설계 (5.9 멀티파트·어셈블리, 5.11 선형 솔버·메모리 지능화, 5.12 HWT-IGA SOLVER 로깅 신규)
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
10. **구조: Kratos Python 연결 모듈 + Python 측 대폭 조정.** 잠김 IGA 특화 작업(격자, 분류, 절단 적분, extension T, 경계 자료, consistent 질량, 후처리)은 Python(Numba)에서, 범용 FE 작업(요소 적분, 재료, 조립, 구속, 솔버, 시간적분)은 Kratos에서 한다. 필요하면 Kratos에서 행렬을 추출해 Python에서 직접 조정·결합한다(R2 경로). C++는 기본 구성이 아니라 대안이다.
11. 메쉬 리더: bdf/fem은 pyNastran, inp·rad·k는 자체 경량 파서, msh는 meshio (D11).
12. 메쉬 내보내기: 5개 포맷(inp, fem, bdf, rad, k)으로 메쉬 + 노드 결과 + 요소 평균 응력 (D12).
13. 일반(native) FE 파트 메쉬는 **Gmsh로 생성 지원**. GPL 격리를 위해 외부 프로세스로만 호출하고 WHT_IGA는 gmsh를 import하지 않는다. 외부 메쉬 파일 입력도 지원 (D13, D10 변경).
14. IGA–일반 FE 혼합 모델은 **Phase 3 Beta**로 앞당긴다 (D14).

### 0.5 v2 → v2.1 개정 이력

| 출처 | 내용 | 반영 위치 |
|---|---|---|
| Phase 0 트랙 A 실측 | S02–S05 통과 → **경로 R1(셀 단위 적분점 기하 + SM 요소 + 구속식) 확정 제안** | §5.3, §5.8 |
| 트랙 A F1 | 적분점 기하에서 SM 요소 consistent 질량 = 0 (원인: `BaseSolidElement::CalculateMassMatrix`가 한 단계 높은 Gauss 방법을 사용, Python 우회 불가). **C++ 없이 해결**: Kratos 추출 K + Python consistent M + scipy. 1차 모드 lumped 대 consistent 차이 0.07%(p=2, 40×4×4), 상위 20모드 최대 7.9%→2.3%(격자 세밀화) | §5.6, §5.8, §7.9 |
| 트랙 A F2 | Block builder가 구속 상수 c를 매 Build마다 다시 조립 → 증분 scheme에서 누적. 규칙: **c는 항상 0, 오프셋은 고정 master 자유도로 표현** + 코드 안전장치 + Newton 검증(S21) | §5.2.3, §5.8 |
| 트랙 A F3 | `SetTable` 온도 테이블은 강성에 반영 안 됨 → Accessor 확인 spike(S20), 실패 시 요소별 Properties | §5.8, §7.11 |
| 트랙 A S13 | 혼합 체적 요소는 적분점 기하에서 access violation → C++ MixedUP0 유지 (Phase 5) | §5.7 |
| 트랙 A S06 | FEAST 미빌드, spectra는 shift 0에서 실패 → 고유치는 scipy `eigsh` 기본 | §5.6 |
| add_1 (G-1–G-8) | 갭 목록. G-1 `.inp` 읽기는 메쉬 리더(§5.10)로 흡수, G-2 진행 콜백·G-5 자동 세분·G-7 출력 이름 매핑(Phase 3), G-3 스펙트럼 case-file 리더(Phase 4), G-4 재료 래더(Phase 5), G-6 단위계 가드(Phase 1), G-8 성능 비교(Phase 6) | §9 |
| add_2 A1–A9 | add_1 정정: SRI 미기재→Phase 5 옵션 추가, 온도 테이블·고유치 판정 정정, body-fit·IGA–FE 혼합·instance는 갭, **ConstitutiveModelsApplication에 Mooney-Rivlin/Yeoh/Ogden 존재**(C++ 승격 사유 아님, spike S19) | §1.3, §5.8 |
| add_2 W0–W5 | Coreform 워크플로 대응 단계, `compile_plan()`(대상 분해·진단) 신설 | §3.4 |
| add_2 M1–M10, P-1–P-10 | 파트 단위 이산화(Immersed/BodyFit/Native), 파트 분리·병합 전략, 메쉬 리더/라이터 6종, instance, body-fit 변환, 혼합, 미리보기 | §5.9, §5.10, §9 |
| 사용자 요청 | **멀티파트·어셈블리를 Kratos 기능(SubModelPart, MPC, mesh tying, 접촉 프로세스, mapper, rigid movement)과 엮어 설계** | §5.9 (신규) |
| 사용자 결정 | D11–D14, 결정 10 | §0.4, §14 |
| 사용자 추가 요구 (2026-11-14) | **모드 중첩 과도/조화 동역학** Phase 4 편입 및 **Explicit 동역학** Phase 6+ 정식 편입 (`MechanicalExplicitStrategy`, `ExplicitCentralDifferencesScheme`, `ExplicitMultiStageKimScheme` 바인딩 실측 확인; 적분점 기하/CFL 검증 S24 추가) | §1.3, §1.4, §5.6, §5.8, §7.9, §7.10, §9, §10.1 |
| 아키텍처·저장소 결정 (2026-11-14) | **D17** `SQLite + VTKHDF` 하이브리드 결과 저장소·통합 쿼리 확정, **D18** Pydantic v2 불변 Spec + Hexagonal(Ports & Adapters) + CQRS 독자 API 설계 및 `wht_iga.compat.coreform` 호환 계층 분리 | §4, §6, §8, §14 |
| 사용자 확장 요구 (2026-11-14) | **IGES(`.igs/.iges`)** 입력 통합, **CAD+Mesh 통합 선택기·기하 계측·접촉 사전 점검 API(`wht_iga.selection`)** 설계, **솔버별 Set 키워드 입출력**(OptiStruct `SET`/`BSURF`, Abaqus `*NSET`/`*ELSET`/`*SURFACE` 등) 계획 추가 | §5.10, §6, §7.1, §8.5, §9 |
| Phase 0 트랙 C 실측 (2026-11-14) | S19–S24 **40/40 통과** (`TRACK_C.md`): S19(`HyperElastic3DLaw`+`UpdatedLagrangian` `6.66e-16`~`3.33e-15`, 다항 초탄성은 C++ 보완 필요), S20($E(T)$ `TableAccessor` `6.013e-09`, $\alpha(T)$·`LinearElastic3DLaw`는 셀별 `Properties`), S21(`c=0`+`MoveMeshFlag=True`+`AndCriteria` `1e-15`), S22(`ContactSM` 고수준 실패→표면층+MPC+`SpringDamper` `5.220e-11`, Phase 5는 직접 `PairedCondition`/활성집합), S23(`nearest_element` `4.44e-16` 통과, 자유 응력 계면은 mortar $L^2$), S24(`MechanicalExplicitStrategy` 질량 소멸 `2.78e-17`→Python Explicit+HRZ lumped $M$+$T^T f^{\text{int}}$로 $\Delta t_{\text{crit}}$ **714,131배** 복원) | §1.3, §5.7, §5.8, §5.9.4, §7.11, §9 |
| 통합 실측: 성능·솔버·병렬화 (2026-11-14/15) | **선형 솔버**: Skyline LU는 3D 격자에서 profile 저장량 폭증(V02 7만 DOF에 16 GB, 44만 DOF 불가) → Kratos 기본 `pardiso_lu`, `pardiso_ldlt/llt`·`num_threads_mkl` 선택(§5.11). Python 쪽 pypardiso는 튜닝 래퍼 `wht_iga.linalg.PardisoFactor`(§5.11.2). **전처리 병렬화**: 셀별 적분 규칙을 프로세스 풀로 분할(결과 동일, V02 레벨 2 76.3 s→23.6 s, §3.3). **Numba 커널**: winding number 무할당 재작성 + 셀 단위 일괄 판정(V06 형상 대체 셀 17×, 일반 셀 6×, 결과 동일, §2.2) | §2.2, §3.3, §5.6, §5.11 |
| extension 연산자 정정 (2026-11-14) | 인덱스 공간 Lagrange 외삽은 clamped 경계 근처(p≥2, 비균일 Greville)에서 1차 함수도 재생 못함(V01-ext p=2 1.3e-2, p=3 4.2e-3). **실제 knot의 단항식 B-spline 계수(Marsden blossom) 재생 연립식**으로 교체 → V01-ext p=1–3 PASS(≤3.1e-7) | §7.5 |
| V02 판정 방식 정정 (2026-11-15, 사용자 결정) | 수렴률은 peak의 Roark 대비 오차가 아니라 최대 주응력 분포의 자기 수렴률(원문 §7.5)로, peak는 Roark·Abaqus 두 기준 보고. 구현 격자를 계획 표기 h = 6/2^k로 맞춤(이전 구현은 10/5/3/2.5/1.5/0.75로 계획과 불일치) | §9 G1, §10.1 |
| V07 Gate·V-COND 모체 확정 (2026-11-14) | **V07**(`v07.json`): Hard Gate = 원문 §6 `slope > 0.5` + 레벨별 오차 한계 `(0.01, 0.2), (0.005, 0.1), (0.002, 0.05)`, Target(목표치) = `2.3 / 2.2` (Coreform 실측 `2.79 / 2.38` 기준). **V-COND** 모체 형상 = **V02 Roark 판(구멍 중심 이동)** 확정 | §9 G2/G5, §10.1 |
| **HWT-IGA SOLVER 로깅 · 메모리 지능형 솔버 선택 편입 (2026-11-15, 사용자 요청)** | V06 L2 `p=2` 20 GB 메모리 부족 실측 반영: ① 공식 솔버 명칭 **`HWT-IGA SOLVER`** 확정 및 On/Off 가능한 상용 수준 **문제·솔버·진행/수렴·평형·리소스 로그** 규격 신설, ② **실시간 시스템 RAM 실측 + 추정 피크 메모리** 기반 `auto` 솔버 지능형 전환(`pardiso_ldlt` ↔ Kratos `AMGCLSolver`) + 명시적 전환 메시지 출력, ③ `out_mode="drop"\|"ghost"` 정식 옵션화를 **Phase 3으로 앞당겨 편입** | §5.11, §5.12, §6, §8, §9 Phase 3 |
| G2 마무리 실측 (2026-11-15 오후, 통합 담당) | **V-COND 설계 정정**: 상자 평행이동은 판을 자르기만 해 작은 셀이 생기지 않았다(ties 0, ON=OFF → 공허 PASS). 격자 위상 sweep을 "마지막 x열에 φh 조각(φ∈[1e-8,1], 셀 η_min 1e-8 도달)"으로 구현하고 `out_mode="drop"`·비공허성 검사(η≤1e-6 도달, ties>0)를 추가. 결과: 평탄 ON 1.06(PASS; κ 바닥은 외팔보 굽힘 모드), OFF 폭증은 φ=1e-4에서만 뚜렷(정성적), **기울기 −1.44(기준 −2±0.3 FAIL)** — 판 두께 방향 nz=1–4(p=2 제어점 3–6층), 셀 10:1로 점근 전. 기준(형상·세분 방식) 결정 대기. V06 L2 p=2 drop PASS(6.6 GB), L2 p=3은 q=6 셀 규칙 비용(꼬리 분포, L0도 20분+)으로 미실행 | §9 G2, §10.1, `WHT_IGA/docs/phase2/G2_report.md` |

> **Phase 1 누적 실무 규약 메모 (반영 대기 목록 — Phase 2 G2 보고 후 일괄 확정)**
> 1. **압력 부호 규약**: 외향 pull(+) 기준 통일.
> 2. **ID 공간 규칙**: 노드·요소·조건·SubModelPart ID 구간 분리 및 충돌 방지 규약.
> 3. **경계피크 probe 방식**: 경계 피크 응력 평가 시 가중치 0 적분점 probe 추출 절차.
> 4. **`qp_table` hook 제안**: 적분점 메타데이터(`qp_table`) 조회·연결 hook 인터페이스.

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
| 격자 | Body-fit mesh: 기존 정렬(mapped) all-hex 메쉬 → 스플라인 변환 (Coreform 2026.3과 같은 범위, sweep·bottom-up 불가) | 정렬 블록 감지 → 블록별 NurbsVolume, 블록 사이 C⁰ 공유 또는 tie (§5.10.4) | 3 | Beta | E04, V10 body-fit 비교 |
| 격자 | 메쉬 미리보기 (해석 전 이산화 확인) | `AnalysisSession.build_mesh(preview=True)` → 격자·셀 상태·소형 셀·표면층 VTKHDF | 1 | Production | 단위시험 |
| 메쉬 | 일반(native) FE 파트 | 메쉬 파일 6종 읽기 또는 Gmsh 외부 프로세스 생성 → Kratos SM 일반 요소 (§5.10) | 3 | Beta | 자체 + X05 |
| 메쉬 | 메쉬 파일 읽기: inp, fem, bdf, rad, k, msh | pyNastran(bdf/fem), 자체 파서(inp/rad/k), meshio(msh) → 내부 메쉬 모델 | 2–3 | Beta | 왕복 시험 |
| 메쉬 | 메쉬 파일 쓰기: inp, fem, bdf, rad, k | 메쉬 + 노드 결과 + 요소 평균 응력, IGA 파트는 proxy hex-tet | 3 | Beta | 왕복 시험 |
| 메쉬 | 메쉬 외곽 면 → immersed 형상 입력 (CAD 없는 모델) | 솔리드 메쉬 경계 면을 tessellation으로 사용 | 2 | Beta | V02 메쉬 입력판 |
| 어셈블리 | Part / Instance(dependent·independent) / Assembly | §5.9. dependent instance는 파트 메쉬 재사용 + 강체 변환 | 1(모델)/3(instance) | Production | V04, X05 |
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
| 결합 | IGA + 일반 FE 혼합 (matching / nonmatching 계면) | 같은 ModelPart에 일반 FE 요소 + MPC tie 또는 mortar $L^2$ projection (`nearest_element` mapper는 운동학적 선형 패치 `4.44e-16` 통과, 자유 응력 계면은 Python mortar $L^2$, S23), `compile_plan()` 대상 분해 (§5.9) | **3** | **Beta** | 자체 + V04 FE판 |
| 접촉 | surface-to-surface, penalty Coulomb μ (자기접촉 제외) | `Triangle3D3` 표면층 + 직접 `PairedCondition` 생성 또는 C++/Python 활성집합 경로 (`ContactSM` 고수준 프로세스는 `GenerateFaces()` 미구현으로 실패, S22) | 5 | Beta | V08, E05 |
| 해석 | 선형 정적, 다단계 | ResidualBasedLinearStrategy / step 루프 | 1/3 | Production | V01–V06 |
| 해석 | NLGEOM | TotalLagrangian / UpdatedLagrangian + Newton | 5 | Beta | E03, E04 |
| 해석 | 고유진동수 (Lanczos) | Kratos 추출 K + Python consistent M + scipy `eigsh`(shift-invert). Kratos `EigensolverStrategy`는 적분점 기하에서 질량 0이라 사용 안 함 (트랙 A F1) | 4 | Production | V09, V10 |
| 해석 | 암시적 동역학 (Coreform: HHT) | **Kratos 시간적분 scheme** (기본 Bossak, 선택 Newmark/BDF). HHT는 구현하지 않음 (사용자 결정 D5) | 4 | Beta | V10 |
| 해석 | 응답 스펙트럼 (ABS/SRSS/NRL/TENP/CQC/GRP/DSC, ALGEBRAIC/SRSS/R40/R30) | 모드 기반 numpy 후처리 | 4 | Beta | V09 |
| 해석 | 랜덤 응답 (PSD BASE/FORCE/DB, Segalman RMS von Mises) | 모드 FRF 후처리 | 4 | Beta | V09 |
| 해석 | 모드 중첩 과도 / 조화 응답 (Modal Transient / Harmonic, WHT_IGA 확장) | 추출 K + Python consistent M 고유모드 → SDOF 시간적분 / FRF 후처리 | 4 | Beta | V10 modal / 해석해 |
| 해석 | 정상/과도 열전달 (전도, film, radiation, 열용량) | IgaApp `LaplacianElement` + ConvDiff `ThermalFace3D3N` (spike) | 4 | Beta | V12 |
| 해석 | 열변형률 (지정 온도장) | ConstitutiveLaws `ThermalElasticIsotropic3D` + 제어점 투영 T | 4 | Beta | V11 |
| 해석 | explicit 동역학 (WHT_IGA 추가 목표, Coreform 미지원) | Python Explicit 루프 + HRZ lumped $M$ + $T^T f^{\text{int}}$ 축약 (또는 C++ $\tilde{\mathbf{C}}^e$ 직접 반영 요소). Kratos `MechanicalExplicitStrategy`는 비변분 $T^T m_{\text{slave}}$로 질량 소멸(`2.78e-17`) 확인 (S24) | 6 | Alpha | V-EXP |
| 재료 | 선형 탄성 (T 의존 표) | `LinearElastic3DLaw` / `CL.ThermalElasticIsotropic3D`. $E(T)$는 `TableAccessor`(`6.013e-09`, S20), $\alpha(T)$·`LinearElastic3DLaw`는 셀별 `Properties` 복제 (`SetTable`은 강성 미반영, 트랙 A S15) | 1/4 | Production | V01, E03 |
| 재료 | 초탄성, Neo-Hookean | `CL.HyperElastic3DLaw` + `UpdatedLagrangianElement3D4N` (`6.66e-16`~`3.33e-15`, S19; `TotalLagrangian`은 `DETJ0` 결함) | 5 | Beta | E04 |
| 재료 | J2 소성 (다선형, T 의존) | `SmallStrainJ2Plasticity3DLaw` 등 + 온도 처리(S20) | 5 | Beta | E03 |
| 재료 | 추가 초탄성 (Mooney-Rivlin, Yeoh, Ogden, 비압축 Neo-Hookean) | `ConstitutiveModelsApplication`은 `MooneyRivlin`/`Yeoh` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 실패(S19) → C++ 보완(`QuadraturePointGeometry::Volume()` 또는 `WhtIgaApplication`) 필요 | 5 | Alpha | patch 시험 + Cook 막 |
| 재료 | 열팽창 (T 의존), 전도도, 비열, 밀도 | $E(T)$는 `TableAccessor`, $\alpha(T)$·`LinearElastic3DLaw`는 셀별 `Properties` 복제(S20) | 4 | Beta | V11, V12 |
| 잠김 방지 | U(p,C^{p−2})/P0 정적 응축, VMS, macro, U2C1/P2C1, B-bar/F-bar, selective integration(SRI) | **C++ `WhtIgaApplication`** MixedUP0 요소 (C++ 승격 1순위, S13에서 Kratos 혼합 요소 크래시 확인) | 5 | Beta | V07, E04 |
| 출력 | U, S, PEEQ, E, NE, LE, NT11 | `CalculateOnIntegrationPoints` → results DB | 1–5 | Production | 전체 |
| 출력 | point/line/multi-point probe (가중치 0 적분점) | 가중치 0 적분점 요소 | 3 | Production | V02, V03 |
| 출력 | 시각화 proxy 격자 (untrimmed hex / trimmed tet) | VTKHDF/VTU | 3 | Production | 단위시험 |
| 출력 | results/quadrature SQLite, probe JSON | 같은 구조 | 1/3 | Production | L2 norm 시험 |
| API | portable spec (to_json/from_json/validate → ValidationError) | dataclass spec + validator | 1 | Production | 단위시험 |
| API | Job (build_mesh, write_input, submit, wait_for_completion) | `Job` | 1 | Production | – |
| 병렬 | MPI 전처리, 스레드/MPI 풀이 | 셀 블록 multiprocessing (QP 수 가중), Kratos OpenMP | 6 | Beta | 성능 리포트 |

### 1.4 범위 밖

- **U-spline**: Coreform 특허 핵심이다. 균일 텐서곱 B-스플라인(C^{p−1})과 THB 세분화로 대체한다. Coreform의 Production 경로도 rectilinear이므로 직교 격자에서 동등한 정확도를 기대할 수 있다.
- **Abaqus 해석 덱 호환**: `.inp` 키워드 번역(하중·step·재료), UEL/UELMAT, `.odb` 출력. 단 `.inp`의 **메쉬·집합 읽기와 메쉬 쓰기**는 범위 안(§5.10).
- **솔버 결과 포맷 쓰기**: `.h3d`, `.op2`, `d3plot`, Radioss animation. (메쉬 + 노드 결과 텍스트 포맷 쓰기는 범위 안)
- **자기접촉**: Coreform도 지원하지 않는다. (Explicit 동역학은 Phase 6 확장 트랙으로 편입, S24)
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
5. **현재 Numba 적용 범위 (2026-11-15 실측)**: `geometry/classify.py`의 삼각형–AABB SAT, 3축 parity 투표(`parallel=True`), winding number(무할당 스칼라 커널, 16점 이상 `prange`, 미만 직렬). 호출 측은 셀 하나의 조회점을 모아 **한 번에** 넘긴다(점 1–8개씩 반복 호출 금지 — prange 디스패치 비용만 남음). 아직 Python인 핫스팟: `quadrature._cell_moments`/`_cross_section_moments`(삼각형 루프), octree 구동부, NNLS(scipy 컴파일 코드, 셀마다 호출). 다음 Numba 대상은 모멘트 계산이며, 결과 동일성(점·가중치 bit 비교)을 수락 조건으로 한다.
6. **전환 시점**: 프로파일링으로 측정한다. Numba로 이미 최적화한 함수가 전체 실행 시간의 30%를 넘을 때만 그 함수를 옮긴다. Numba 버전은 결과 비교용 기준 구현으로 남긴다.

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
- 외부(예: Python 3.11 aieng) 연동용 `wht_iga/runner/external.py`는 `subprocess`로 `micromamba run -n whtiga python -m wht_iga run study.json`을 실행하고 산출물 경로를 돌려준다. 도입 시점은 §14 결정 사항이다.
- 전처리 병렬화: `multiprocessing`(spawn)으로 셀 블록 단위 작업을 나눈다. 분할 가중치는 예상 적분점 수(CUT 셀은 IN 셀의 10–30배)로 정한다.
  - **구현 (2026-11-15)**: `solver/trimmed.build_trimmed_table(workers=)` — OUT이 아닌 셀의 `cell_quadrature`를 `ProcessPoolExecutor`로 분할(작업 수 × 8 청크, 셀 순서 유지 → 결과·순서가 직렬과 bit 단위 동일). 기본 작업 수 = 물리 코어(`os.cpu_count()//2`), `WHT_GEO_WORKERS`로 지정, 400셀 미만은 직렬(spawn + numba 캐시 로드 비용). 작업 프로세스는 `NUMBA/OMP/MKL/OPENBLAS_NUM_THREADS=1`(중첩 스레드 과다 방지).
  - 실측: V02 레벨 2(잘린 셀 1,667) 직렬 76.3 s → 7 프로세스 23.6 s(다른 해석 1코어 점유 중 측정). QP 가중 분할은 미구현(균등 청크) — 불균형이 측정되면 도입.
  - Windows spawn이므로 스크립트는 `if __name__ == "__main__":` 보호가 필수다.
- **절단 적분 캐시 (구현 2026-11-15)**: `build_trimmed_table(cache_dir=)`는 (테셀레이션 노드·삼각형, 격자 lo/hi/n_cells/degree, q) **+ `classify.py`·`quadrature.py`·`trimmed.py` 소스**의 SHA-256을 키로 결과를 `WHT_IGA/.cache/trimmed/<key>.pkl`에 저장한다(원자적 교체, git 제외). 코드가 바뀌면 키가 바뀌어 오래된 결과를 읽을 수 없다. BC·extension θ·ghost α·솔버만 바꾸는 sweep은 두 번째 실행부터 절단 적분 0 s. `WHT_GEO_CACHE=off`로 끈다. 시험 `tests/unit/test_trimmed_cache_parallel.py`(병렬=직렬, 캐시 적중=무캐시, 키 분리).

---

### 3.4 사용자 워크플로 W0–W5 (Coreform 워크플로 대응, add_2 §3)

Coreform 단계(`user-guide/user-guide/workflow-overview.md`): ① Abaqus/CAE 모델 설정 → Job이 `<job>.inp` + `<job>.cf` 생성 → ② `coreform_iga_mesh`가 절단 메쉬·적분·추출/확장 연산자·function splitting을 `<job>.sql`로 → ③ interop가 키워드 번역, `*USER ELEMENT` 덱, 사용자 요소 자료 → ④ Abaqus/Standard가 매 반복 Coreform 사용자 요소에 K·M·f 요청, 해석 중 `_iga_results.sql`·`_iga.odb`·`probe_data.json` 생성 → ⑤ Abaqus/CAE·ParaView 후처리.

| # | Coreform | WHT_IGA | 산출물 | 모듈 |
|---|---|---|---|---|
| W0 | Abaqus/CAE 모델 설정 | Python API로 Part(형상 출처: CAD/STL/메쉬 파일)·파트별 이산화·Instance·재료·Step·하중·BC·상호작용·Probe 정의 | 메모리 `Model` | `spec/`, `plan/` |
| W1 | Job → `.inp` + `.cf` | `StudySpec.model_dump_json()`: `validate()` + 직렬화, 참조 파일 해시 기록 | `<job>.spec.json` | `spec/study.py`, `spec/base.py` |
| W2 | `coreform_iga_mesh` → `.sql` | `AnalysisSession.build_mesh()`: **파트별 병렬**. IGA 파트 S1–S6, 일반 FE 파트는 메쉬 읽기·Gmsh 생성·검사. `preview=True`면 미리보기만 | `<job>.wiga.sqlite`, (미리보기) `<job>_preview.vtkhdf` | `geometry/`, `grid/`, `classify/`, `quadrature/`, `stabilize/`, `meshio_ext/`, `meshgen/` |
| W3 | interop (키워드 번역) | `compile_plan()`: 하중·BC·출력 대상 세트를 IGA/일반 FE 파트로 분해, 상호작용(tie/contact/coupling)을 Kratos 객체 명세로 변환(§5.9), 표면층 생성, 해석 불가 조합은 거부 + 진단 | `<job>.resolved.json` | `plan/compiler.py` |
| W4 | Abaqus/Standard + 사용자 요소 | `AnalysisSession.run()`: Kratos ModelPart 구축(§5.9.2 계층) → 전략 실행 → 프레임 단위 결과 기록, 진행 콜백 | `<job>_iga_results.sqlite`, `<job>_probe_data.json`, `<job>.log` | `kratos_backend/`, `runner/session.py` |
| W5 | Abaqus/CAE, ParaView alpha reader | `job.results` API, **IGA proxy 격자 + 일반 FE 요소 + 표면층을 한 파일에**(unified output), probe JSON, 메쉬 포맷 5종 내보내기 | `<job>.vtkhdf`, 내보내기 파일 | `post/`, `meshio_ext/` |

Coreform과의 차이: 라이선스 대기 없음, 키워드 번역 대신 spec 수준 대상 분해(W3), Kratos가 직접 조립(사용자 요소 호출 없음), 캐시가 있으면 W2 생략.

**API 이름 (D18 반영)**: 이 문서의 단계 이름은 다음 API에 대응한다 — W1 직렬화 = `StudySpec.model_dump_json()`, W2 = `AnalysisSession.build_mesh()`, W3 = `compile_plan()` → `ResolveReport`(대상 분해·구속 그래프 검사·진단), W4 = `AnalysisSession.run()`, W5 = `ResultView`. 모듈은 §6(`spec/`, `plan/`, `runner/`, `query/`)을 따른다.

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

### 4.3 `<job>_iga_results.sqlite` (메타·요약·Probe·인덱스 카탈로그, D17)

대용량 다중 타임스텝 3D 필드(Implicit/Explicit 동역학)는 SQLite 행 삽입 병목을 피하기 위해 **VTKHDF(`§4.5`)에 연속 배열로 저장**하고, SQLite는 메타·요약·수렴·모드·Probe·HDF5 슬라이스 인덱스를 담당한다(소규모 정적 해석만 `cp_field`/`qp_field` 선택 병행). `ResultView`(§8)가 두 저장소를 투명하게 통합 쿼리한다.

```sql
CREATE TABLE step  (step_id INTEGER PRIMARY KEY, name TEXT, kind TEXT, t_end REAL);
CREATE TABLE frame (frame_id INTEGER PRIMARY KEY, step_id INTEGER, time REAL, increment INTEGER,
                    converged INTEGER, n_iter INTEGER, hdf_group TEXT);
CREATE TABLE cp_field (frame_id INTEGER, var TEXT, data BLOB);      -- 소규모 정적용 (대용량은 VTKHDF)
CREATE TABLE qp_field (frame_id INTEGER, var TEXT, data BLOB);      -- 소규모 정적용 (대용량은 VTKHDF)
CREATE TABLE mode  (mode_id INTEGER PRIMARY KEY, freq_hz REAL, gen_mass REAL,
                    part_x REAL, part_y REAL, part_z REAL, eff_mass_x REAL, eff_mass_y REAL, eff_mass_z REAL);
CREATE TABLE energy (frame_id INTEGER, kinetic REAL, strain REAL, plastic REAL);
CREATE TABLE reaction (frame_id INTEGER, set_name TEXT, fx REAL, fy REAL, fz REAL, mx REAL, my REAL, mz REAL);
```

### 4.4 `<job>_iga_quadrature.sqlite` (선택, Coreform `--quadrature-output` 대응)

```sql
CREATE TABLE qp     (qp_id INTEGER PRIMARY KEY, cell_id INTEGER, material INTEGER,
                     x REAL, y REAL, z REAL, dV REAL);
CREATE TABLE qp_val (frame_id INTEGER, qp_id INTEGER, U BLOB, S BLOB);   -- L2 오차 노름 계산용
```

### 4.5 기타 출력

- `<job>.vtkhdf`: **대용량 필드·시계열 주 저장소** (ParaView 네이티브 열기 + `h5py` 타임스텝/영역별 $O(1)$ 슬라이싱, §7.12). `.vtu`는 단발 내보내기용.
- **제어점 계수 저장 (재시작·재평가용)**: VTKHDF 파일 안 별도 그룹 `/WHT_IGA/control_points/<frame>/{U, NT11, state}`에 제어점 자유도와 적분점 상태 변수를 저장한다. proxy 격자 필드는 이 계수로 언제든 다시 만들 수 있다. SQLite `cp_field`/`qp_field`는 소규모 정적 해석에서만 선택 사용.
- `<job>_probe_data.json`: Coreform `probe_data.json` 구조를 따른다.
- `<job>.spec.json`: 입력 재현용 `StudySpec` (Pydantic JSON Schema 호환).

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
- **구속 상수 규칙 (트랙 A F2)**: Kratos Block builder는 상수 c를 매 `Build`마다 다시 조립하므로 증분 scheme에서 누적된다(선형 전략에서 slave = Ru + 2c 관측). 따라서 **c는 항상 0**이다. 고정값이 있는 master는 관계식에 그대로 두면 Kratos가 처리한다. 오프셋이 필요하면 고정된 보조 master 자유도로 표현한다. `constraints.add_linear_tie`는 0이 아닌 c를 오류로 거부한다(Phase 1 보강). Newton 전략 검증은 S21.
- 구속식 생성 바인딩 순서: `CreateNewMasterSlaveConstraint(name, id, master_dofs, slave_dofs, R, c)` — master 목록이 먼저. 크기 0 관계 행렬은 크래시.
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

**트랙 A 결과 (2026-11-14)**: S03·S04·S05 모두 통과 → **R1 확정 제안** (ADR-001은 트랙 B 수령 후 통합 담당이 작성). R2(행렬 추출·Python 조정)는 대체 경로가 아니라 **질량 행렬·고유치·선형 과도에 상시 병용**한다(§5.6).

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
| 선형 정적 | `ResidualBasedLinearStrategy` + `ResidualBasedIncrementalUpdateStaticScheme` + `ResidualBasedBlockBuilderAndSolver` + `pardiso_ldlt` (현재 구현 `pipeline.run_linear`: 수렴 기록(F3)을 위해 `ResidualBasedNewtonRaphsonStrategy` + `ResidualCriteria`, 솔버 기본 `pardiso_lu` — §5.11) | `sparse_lu`, 또는 추출 후 pypardiso. **`SkylineLUFactorizationSolver`는 3D 격자에 사용 금지**(소형 시험·디버그만) |
| 비선형 | `ResidualBasedNewtonRaphsonStrategy` + `DisplacementCriteria`/`ResidualCriteria`, step 루프 | line search |
| 질량 행렬 | **Python(Numba) consistent M** (적분점 기하에서 SM 요소 consistent 질량 = 0, 트랙 A F1). Kratos lumped 질량(`COMPUTE_LUMPED_MASS_MATRIX`)은 정상 | 소형 C++ 요소(정확도·성능 기준 미달 시에만) |
| 고유치·응답 스펙트럼·랜덤 | Kratos 추출 K(`GetSystemMatrix` + `scipy_conversion_tools`, 무구속 상태에서 추출 후 Python에서 Tᵀ K T) + Python consistent M + scipy `eigsh`(shift-invert) | – |
| 선형 과도 | 같은 행렬 + **Python Bossak 시간 루프**(D5의 Bossak 알고리즘을 Python에서 실행) | Kratos Bossak + lumped M |
| 비선형 동적 | `ResidualBasedNewtonRaphsonStrategy` + **Kratos Bossak scheme** (`damp_factor_m` α_B, 기본 −0.3) + lumped M (1차 모드 오차 ~0.1%) | 상위 모드가 중요하면 Python Newton 루프(매 반복 Kratos `Build`) + Python consistent M |
| Explicit 동적 (Phase 6) | SM `MechanicalExplicitStrategy` + `ExplicitCentralDifferencesScheme` (`ExplicitMultiStageKimScheme` 선택) + **Kratos lumped M**. 구속식(extension T)은 Kratos가 처리한다: 초기화 시 `ConstraintUtilities::PreComputeExplicitConstraintMassAndInertia`로 질량을 master에 모으고, 매 단계 `PreComputeExplicitConstraintConstribution`으로 slave 잔차를 `master += Tᵀ·slave`로 더한 뒤 slave 값을 복원한다(`mechanical_explicit_strategy.hpp` l.260–366, `constraint_utilities.cpp`) | S24 실패 시 Python explicit 루프(Kratos 내력 + Python Tᵀ 축약 + lumped M) |
| 열 정상/과도 | IgaApp `LaplacianElement` + `SupportLaplacianCondition`(penalty T) + BDF. film·radiation은 표면층 노드(TEMPERATURE slave) 위 ConvDiff `ThermalFace3D3N` | 참조 조립기 열해석 |
| 열변형 | `ThermalElasticIsotropic3D` (노드 TEMPERATURE 보간, `REFERENCE_TEMPERATURE`). 온도장은 제어점에 L2 투영 | SM 초기 변형률 주입 |
| 응답 스펙트럼 / 랜덤 / 모드 중첩 | 모드 + 참여계수 Γ = φᵀ f_r (`VOLUME_ACCELERATION` 또는 일반 외력 `BuildRHS`) → numpy SDOF 시간적분(과도) / FRF(조화·랜덤) / 스펙트럼 조합 | – |

### 5.7 C++ 승격 기준 (`WHT_IGA/cpp/WhtIgaApplication`)

다음 중 하나라도 해당하면 승격한다.

1. **잠김 방지 혼합 요소** (V07, E04 필수). 1순위.
2. S05에서 구속식 오버헤드가 조립 시간의 50%를 넘음.
3. Nitsche 또는 ghost penalty가 필요함.
4. 접촉 S14/S22 경로: `ContactSM` 고수준 프로세스(`mesh_tying_process`, `penalty_contact_process`)는 IGA 요소 `GenerateFaces()` 미구현으로 실패(S22) → Phase 5 접촉은 직접 `PairedCondition` 생성 또는 C++/Python 활성집합 경로 채택.
5. 다항 초탄성(Mooney-Rivlin/Yeoh/Ogden, S19): `ConstitutiveModelsApplication` Python 미노출 및 `QuadraturePointGeometry::Volume()` 미구현으로 실패 → C++ 보완(`QuadraturePointGeometry::Volume()` 또는 `WhtIgaApplication`) 필요.
6. Explicit 동역학 고성능화(S24): `MechanicalExplicitStrategy`가 비변분 $T^T m_{\text{slave}}$로 질량 소멸(`2.78e-17`) → Python Explicit + HRZ lumped $M$ + $T^T f^{\text{int}}$가 주 경로이며, 대규모 모델에서는 C++ `ExtendedSplineSolidElement`($\tilde{\mathbf{C}}^e$ 직접 반영) 승격 검토.

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
| S19 | ConstitutiveModelsApplication 재료(Mooney-Rivlin, Yeoh, Ogden)가 SM 요소 + 적분점 기하에서 동작? (로컬 복원·수정분 앱) | 단축 인장 해석해 1e-8 | ConstitutiveLaws 보유 재료만 사용, 필요 시 C++ |
| S20 | Accessor(위치·온도 기반 물성 접근자)로 적분점별 E(T), α(T), σ_y(T) 반영? | 온도 구배 막대 해석해 1e-8 | 요소(셀)별 Properties 복제 |
| S21 | 구속식(c = 0 규칙) + Newton 전략(비선형)에서 slave 값·수렴 정상? | 대변형 블록, 참조 Tᵀ K T Newton과 1e-8 | Python Newton 루프(R2) |
| S22 | ContactSM `mesh_tying_process`(mortar)와 `penalty_contact_process`가 **스프링/MPC로 스플라인에 결합된 표면층 노드**에서 동작? | 두 블록 tie: 단일 블록 대비 1e-6, Hertz 축소판 수렴 | MPC tie(§5.9.4) / C++ |
| S23 | `AssignMasterSlaveConstraintsToNeighboursUtility`와 MappingApplication `nearest_element` mapper로 비정합 계면 MPC 가중치를 만들 수 있는가? | 패치 시험(일정 변위 전달) 1e-12 | Python 투영(Numba)으로 가중치 직접 계산 |
| S24 | Kratos `MechanicalExplicitStrategy` + 적분점 기하 SM 요소 + extension MPC(T) + Kratos lumped M으로 1D 탄성파를 재현하고, 소형 절단 셀이 있어도 Δt_crit이 배경 격자 수준(~h/c)으로 유지되는가? HRZ lumping과 비교 | 1D 파동 해석해 1% 이내, 소형 셀 Δt_crit ≥ 0.5 × 비절단 격자 값 | Python explicit 루프, selective mass scaling, 필요 시 C++ `ExtendedSplineSolidElement` |
| S25 | ContactSM `explicit_penalty_contact_process`가 스프링으로 스플라인에 결합된 표면층에서 동작하는가? | 두 블록 충돌: 운동량 보존 1e-6, 에너지 오차 < 2% | Python 활성집합 penalty (explicit) |

**트랙 A 실측 (2026-11-14)**: S02·S03·S04·S05·S16 **예**, S06 **부분**(Kratos 고유치 경로 불가, Python consistent M 경로 검증), S13 **아니오**(C++ 필요, 예상대로), S15 **부분**(테이블 부착 OK, 강성 미반영). 상세: `WHT_IGA/docs/spikes/`. 트랙 B(S01, S07–S12, S14, S17, S18) 진행 중.

**트랙 C 실측 (2026-11-14, 40/40 통과 — `WHT_IGA/docs/spikes/TRACK_C.md`)**:
- **S19 (초탄성, 16/16)**: `ConstitutiveModelsApplication`은 `MooneyRivlin`/`Yeoh` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 `RuntimeError` 실패. 대안 `CL.HyperElastic3DLaw` + `UpdatedLagrangianElement3D4N`은 `6.66e-16`~`3.33e-15` 일치(`TotalLagrangian`은 `DETJ0` 결함). 다항 초탄성(Mooney-Rivlin/Ogden 등)은 C++ 보완(`QuadraturePointGeometry::Volume()` 또는 `WhtIgaApplication`) 필요.
- **S20 (온도 의존 물성, 6/6)**: `TableAccessor` + `CL.ThermalElasticIsotropic3D`로 $E(T)$ 반영 성공(상대오차 `6.013e-09` < `1e-8`). 단 $\alpha(T)$ 및 `LinearElastic3DLaw`는 `TableAccessor` 무시 → 셀별 `Properties` 복제로 $O(h^2)$ 수렴(`1.049e-05` @ `n_x=32`) 채택.
- **S21 (구속식 + 비선형 Newton, 3/3)**: `LinearMasterSlaveConstraint`($c=0$) + `TotalLagrangian` Newton에서 참조 $T^T K_t T$ 대비 `1.44e-15`~`2.12e-16` (`1e-15` 수준) 일치. **`MoveMeshFlag=True`** + **`AndCriteria(ResidualCriteria, DisplacementCriteria)`** 필수.
- **S22 (ContactSM + 표면층, 4/4)**: `ContactSM` 고수준 프로세스(`mesh_tying_process`, `penalty_contact_process`)는 IGA 요소 `GenerateFaces()` 미구현으로 실패. 대안 `Triangle3D3` 표면층 + MPC + `SpringDamperElement3D` tie는 단일 블록 대비 `5.220e-11` 통과. Phase 5 접촉은 직접 `PairedCondition` 생성 또는 C++/Python 활성집합 경로 채택.
- **S23 (비정합 계면 mapper, 5/5)**: `nearest_element` mapper로 단위분할 `0.0`, 운동학적 선형 패치 `5.55e-17`~`4.44e-16` 통과(`AssignMasterSlaveConstraintsToNeighboursUtility`는 `7.94e-02` 미달). 자유 응력 계면은 비변분 오차(`0.086`)로 Python mortar $L^2$ projection 필요.
- **S24 (Explicit 동역학, 6/6)**: `MechanicalExplicitStrategy`는 비변분 $T^T m_{\text{slave}}$를 더해 음수 확장 가중치 시 질량 소멸(`2.78e-17`) 확인. Python Explicit + HRZ lumped $M$ + $T^T f^{\text{int}}$ 축약으로 1D 파동 전파($L^2$ 오차 `1.16e-02`) 및 $\eta=10^{-6}$ 절단 셀 $\Delta t_{\text{crit}}$ **714,131배** 복원 입증.

### 5.9 멀티파트·어셈블리 설계 (Kratos 기능 연계)

Coreform은 Abaqus의 part/instance/assembly, tie, coupling, 접촉 기능 위에 IGA 사용자 요소를 얹는다. WHT_IGA는 같은 개념을 **Kratos의 기존 기능**에 대응시킨다. 원칙: 새로 만드는 것은 잠김 IGA에 특화된 부분(스플라인 장, 표면층, 계면 가중치)뿐이고, 결합·접촉·풀이는 Kratos 기능을 쓴다.

#### 5.9.1 개념 모델

| 개념 | 정의 | Coreform 대응 |
|---|---|---|
| **Part** | 형상(CAD/STL/메쉬 파일) + 이산화 1개(Immersed / BodyFit / Native / 강체) + 재료 셀 구성 | Abaqus part + Coreform Mesh |
| **Material cell** | 파트 안의 재료 영역. 같은 파트의 셀들은 **스플라인 장(또는 FE 노드)을 공유** → 계면이 완전 접합 | Boolean merge multi-cell |
| **Instance** | Part 참조 + 강체 변환(회전 R, 이동 t). dependent(파트 메쉬 재사용) / independent(자기 메쉬) | Abaqus instance |
| **Assembly** | Instance 목록 + Interaction 목록 | Abaqus assembly |
| **Interaction** | Tie, Contact, Coupling(kinematic/distributing), (이후) Connector | Abaqus interaction/constraint |
| **Surface** | Instance의 면 집합(CAD 면 ID 또는 메쉬 면). IGA 파트는 표면층으로 실체화 | Abaqus surface / SFM3D3 층 |

메쉬 소유 규칙(Coreform `scripting-guide/meshes.md`와 동일): 파트의 instance가 모두 dependent면 파트가 메쉬 1개 소유, 모두 independent면 instance마다 메쉬 1개, 한 패밀리에서 혼용 금지(`validate()` 오류).

#### 5.9.2 Kratos ModelPart 계층

전체 어셈블리는 **하나의 루트 ModelPart = 하나의 연립방정식**이다. 파트와 상호작용은 SubModelPart로 구분한다.

```
Structure                              (루트: 모든 노드·요소·조건·구속, 단일 Builder-and-Solver)
├── Parts
│   └── <part>.<instance>             (instance별 SubModelPart: 노드 id 구간 = 오프셋 표로 추적)
│       ├── Cells.<material>          (재료 셀별 요소. Properties id = 재료×파트)
│       ├── Layers.<surface>          (IGA 표면층 노드 + Triangle3D3 조건 / FE 면 조건)
│       ├── Virtual.<surface>         (표면 적분점 가상 노드: 스플라인 slave, §5.4 BC-B/C)
│       └── Sets.<name>               (사용자 집합: 노드·면)
├── Interactions
│   ├── Tie_<n>                       (MPC 또는 mortar tying 대상 노드·조건)
│   ├── Contact_<n>                   (master/slave 표면층 조건 → ContactSM 프로세스 입력)
│   └── Coupling_<n>                  (RP 노드 + 결합 구속)
└── Resolved                          (compile_plan()가 만든 하중·BC 대상 SubModelPart)
    ├── BC_<n>, Load_<n>
```

- 노드 id: instance마다 연속 구간을 배정하고 `<job>.resolved.json`에 (instance, 구간, 출처: 제어점/FE 노드/표면층/가상/RP) 표를 남긴다. 결과·진단을 instance 단위로 되찾는 데 쓴다.
- Properties: (파트, 재료 셀)마다 1개. 온도 의존이 위치별로 다르면 셀(요소)별 Properties 복제(S20 실패 시).
- DOF: 모든 instance가 같은 변수(DISPLACEMENT, 필요 시 TEMPERATURE, RP는 ROTATION 추가)를 쓴다.

#### 5.9.3 Instance 구현

- **dependent instance**: W2 메쉬 단계는 파트당 1회(파트 좌표계). W4에서 instance마다 노드를 새로 만들고 좌표에 (R, t)를 적용한다.
  - IGA 파트: 제어점 좌표에 변환을 적용하면 기하 사상은 x = R·ξ + t가 된다. 회전은 detJ = 1을 유지한다(S02에서 Frame 회전으로 확인). 적분점 입력은 계속 **파라미터(파트) 좌표**로 넣는다(S02 확인 사항).
  - 일반 FE 파트: 노드 좌표만 변환.
- **independent instance**: instance 위치에서 메쉬를 따로 만든다(전역 축 정렬 격자가 필요할 때 등).
- 같은 파트의 dependent instance 여러 개는 노드·자유도는 별개지만 요소 적분점 자료(적분점, 형상함수)는 같다 → Phase 6에서 재사용 최적화 후보.

#### 5.9.4 상호작용 → Kratos 기능 대응

| 상호작용 | Kratos 기능 | WHT_IGA가 준비하는 것 | Phase | spike |
|---|---|---|---|---|
| 파트 내부 재료 계면 (multi-cell) | 없음(같은 노드 공유) | 재료 셀별 적분점·Properties | 3 | – |
| Tie, 노드 일치 (FE–FE, FE–표면층) | `LinearMasterSlaveConstraint` 항등 | 일치 노드 쌍 탐색 | 3 | S05 |
| Tie, 비정합 (FE–FE, FE–IGA) | ① MPC: slave 노드 = 상대 면 형상함수 가중합. 가중치는 MappingApplication `nearest_element` mapper(운동학적 선형 패치 `4.44e-16` 통과; `AssignNeighboursUtility`는 `7.94e-02` 미달로 제외, S23), 자유 응력 계면은 Python mortar $L^2$ projection ② `Triangle3D3` 표면층 + MPC + `SpringDamperElement3D` tie (`5.220e-11` 통과, S22; `ContactSM` `mesh_tying_process`는 `GenerateFaces()` 미구현 실패) | 투영·가중치, IGA 쪽은 가상 노드(스플라인 slave)를 master로 쓸 수 없으므로(K11) **IGA 면의 가중치는 제어점 자유도까지 합성해서** 한 줄의 관계식으로 만든다 | 3 | S22, S23 |
| Tie, IGA–IGA | 양쪽 모두 스플라인 → 한쪽을 slave로 둘 수 없음(제어점은 이미 다른 구속의 master일 수 있음) → **penalty 스프링**(§5.4 BC-C: V_A–V_B, 강성 β·w_q) 또는 `Triangle3D3` 표면층+MPC+`SpringDamperElement3D` tie / Python mortar $L^2$ (S22/S23) | 가상 노드 쌍, β | 3 | S10, S22 |
| 접촉 (surface-to-surface, penalty Coulomb μ, 자기접촉 제외) | `Triangle3D3` 표면층 + 직접 `PairedCondition` 생성 또는 C++/Python 활성집합 경로 + 접촉 전용 수렴 기준 (`ContactSM` 고수준 프로세스는 `GenerateFaces()` 미구현으로 실패, S22) | IGA 파트 표면층(독립 자유도 노드 + Triangle3D3 조건)을 스프링으로 스플라인에 결합(Coreform SFM3D3 방식과 동일) | 5 | S14, S22 |
| Kinematic coupling (RBE2형) | SM `ImposeRigidMovementProcess` 또는 MPC u = u_RP + θ × r | RP 노드(DISPLACEMENT + ROTATION), 면 점 → 가상 노드 | 3 | – |
| Distributing coupling (RBE3형) | MPC (RP 변위·회전을 면 점들의 가중 최소자승으로) | 가중치 계산(Numba) | 3 | – |
| 거리 유지 링크 | SM `LinkConstraint` (`link_constraint_process`) | 노드 쌍 | 5+ | – |
| 볼트 프리로드, 커넥터 | – | 범위 밖 (Coreform 문서에도 IGA 검증 없음) | – | – |

규칙:
- **구속 상수 c = 0** (§5.2.3)과 **연쇄 구속 금지**(K11)를 상호작용 전체에 적용한다. `compile_plan()`가 같은 자유도가 두 구속의 slave가 되거나, slave가 다른 구속의 master가 되는 경우를 찾아 거부한다. 비선형 Newton 해석 시 **`MoveMeshFlag=True`** 및 **`AndCriteria(ResidualCriteria, DisplacementCriteria)`**를 필수 적용한다(S21 `1e-15` 일치).
- tie 방식 자동 선택: 노드 일치 → 항등 MPC, 비정합 FE 면 → 가중 MPC, IGA가 낀 경우 → 합성 MPC 또는 penalty, 사용자가 mortar 지정 시 → mesh tying.

#### 5.9.5 풀이

- 모든 instance와 상호작용을 한 루트 ModelPart에서 단일 시스템으로 푼다(Kratos Block builder가 MPC 처리).
- 접촉이 있으면 ContactSM의 정적·동적 solver 래퍼 대신, 같은 프로세스·수렴 기준을 우리 전략 구성(`kratos_backend/analyses/`)에 직접 넣는다(질량·고유치 경로가 Python 쪽이므로).
- 파트별 병렬은 W2 메쉬 단계에서만. 풀이는 Kratos OpenMP. MPI는 Phase 6에서 Kratos MPI 빌드 여부 확인 후.

#### 5.9.6 출력

- 결과 DB(`<job>_iga_results.sqlite`)에 instance 열을 둔다. VTKHDF는 instance별 블록(IGA proxy hex-tet, FE 요소, 표면층)을 한 파일에 기록한다.
- 반력·접촉력·계면력은 상호작용 SubModelPart 단위로 합산해 보고(예: door-lock RF3).

#### 5.9.7 검사 규칙 (`validate()` / `compile_plan()`)

- instance 메쉬 소유 규칙 위반, 서로 다른 파트의 면이 맞닿는데 상호작용 미지정, 파트 간 체적 겹침, 하중·BC 대상이 IGA와 FE에 걸쳐 분해 불가, 연쇄 구속, 0이 아닌 구속 상수 → 모두 issue(오류).

### 5.10 메쉬 입출력과 일반 FE 파트

#### 5.10.1 파트별 이산화 선택과 분리·병합 전략

| 상황 | 권장 | 이유 |
|---|---|---|
| 다른 재료가 완전 접합된 한 물체 (glass-to-metal seal, bimetallic) | 한 파트 + multi-cell | 계면이 장 안에 있어 정확, tie 오차 없음 (Coreform V04: merge–tie 차 1e-7) |
| 접촉·미끄럼·분리되는 부품 (door-lock) | 파트 분리 + 접촉 | 접촉은 서로 다른 장 사이에서만 의미 |
| 크기·스케일이 크게 다른 부품 | 파트 분리, 파트별 격자 크기·Frame | 한 격자면 작은 부품 해상도에 맞춰 DOF 폭증 |
| 기울어진 길쭉한 부품 | 파트별 object-aligned Frame | 절단·소형 셀 감소 |
| 체결된 별개 부품 | 파트 분리 + tie | 조립 구조 유지 |
| 같은 부품 여러 개 | 파트 1개 + dependent instance | 메쉬 단계 1회 |
| 이미 FE 메쉬가 있는 부품, 볼트·패드 등 | `NativeMesh` 파트 + tie | 기존 자산 재사용 |
| 경계 정확도가 중요한 단순 블록 형상 | `BodyFitMesh` | 절단 없음 |

자동 제안: 입력의 물체(solid / 메쉬 파트 ID)마다 파트 1개 → 재료가 달라도 면을 공유하면 multi-cell 병합을 **제안**(확정은 사용자).

#### 5.10.2 내부 메쉬 모델과 리더

- `MeshModel{nodes, blocks[elem_type, conn, part_id], parts{id: name}, node_sets, elem_sets, face_sets(elem, local_face), frames}`. 요소 유형 정규화(C3D8/HEXA8/CHEXA/BRICK/ELEMENT_SOLID → `hex8`, `tet4`, `tet10`, `hex20`, `wedge6`), 포맷별 노드 순서 변환표.
- 범위: 메쉬와 집합만. 해석 덱의 하중·재료·step 키워드는 읽지 않는다(재료·속성 ID는 파트 구분용).

| 포맷 | 파트 대응 | Set 키워드 (읽기 · 쓰기 공통 매핑) | 구현 (D11) |
|---|---|---|---|
| `.inp` (Abaqus/CalculiX) | `*PART`/`*INSTANCE`(변환 포함) 또는 `*SOLID SECTION` ELSET | `*NSET`, `*ELSET`, `*SURFACE, TYPE=ELEMENT`/`NODE` (`S1`–`S6` 매핑) | 자체 경량 파서·라이터 |
| `.bdf`/`.nas` (Nastran) | `PSOLID` PID | `SET <id>` (GRID/ELEM, `THRU` 지원), `BSURF` | pyNastran (BSD-3, `extras: mesh`) |
| `.fem` (OptiStruct) | PID + `$$HMNAME COMP` 주석 | `SET <id>` (GRID/ELEM, `THRU` 압축), `SURF`/`BSURF`, `$$HMNAME SETS` | pyNastran + 주석 파싱·라이터 |
| `.rad` (Radioss starter) | `/PART` | `/GRNOD/NODE`, `/GRBRIC`, `/SURF/SEG` | 자체 파서·라이터 |
| `.k`/`.key` (LS-DYNA) | `*PART` PID, `*INCLUDE` 추종 | `*SET_NODE_LIST`, `*SET_SEGMENT`, `*SET_SOLID` | 자체 파서·라이터 |
| `.msh` (Gmsh) | physical group | `Physical Group` (점/선/면/체적 태그 ↔ 이름) | meshio |

- 셸·빔·스프링 요소는 읽되 issue로 보고하고 해석에서 제외(3D 솔리드 전용).
- 단위는 파일에 없으므로 `unit_system` 필수(add_1 G-6). Nastran `CORD2R` 등 좌표계는 전역으로 변환.

#### 5.10.3 메쉬의 세 가지 쓰임

1. `NativeMesh` 파트 → Kratos SM 일반 요소(`SmallDisplacementElement3D8N`, `3D10N` 등)
2. `BodyFitMesh` 입력 → 정렬 블록을 스플라인으로 변환(§5.10.4)
3. `ImmersedRectilinearMesh` 형상 입력 → 솔리드 메쉬 외곽 면을 tessellation으로 사용 (CAD 없이 기존 FE 모델을 IGA로 재해석)

#### 5.10.4 Body-fit 변환

- all-hex 메쉬에서 정렬 블록 감지(노드 차수, 면 인접 그래프로 (i, j, k) 위상 확인). mapped 블록만, sweep·bottom-up 불가(Coreform 2026.3과 동일 범위).
- 블록마다 노드 격자를 제어점 망으로 보간·근사 → `NurbsVolumeGeometry`(p, C^{p−1}). 블록 사이 경계 제어점이 일치하면 같은 노드(C⁰), 아니면 tie.
- Kratos 쪽 경로는 immersed와 같다(적분점 기하 + SM 요소, 절단·extension 없음).

#### 5.10.5 Gmsh 메쉬 생성 (D13)

- `wht_iga/meshgen/gmsh_runner.py`: STEP/BREP + 크기·차수 옵션 → **Gmsh 외부 프로세스**(`gmsh` 실행 파일, 또는 별도 Python 프로세스의 gmsh API) → `.msh` → meshio 리더 → `MeshModel`.
- WHT_IGA 프로세스는 gmsh를 import하지 않는다(GPL 격리). gmsh는 `whtiga` env의 별도 실행 파일로 설치하고, 없으면 issue로 안내.
- 옵션: tet4 / tet10 / hex(transfinite 가능한 블록), 크기 필드, physical group = CAD 면·솔리드 ID 보존.

#### 5.10.6 메쉬 쓰기 및 Set 키워드 출력 (D12)

- 같은 `MeshModel` 및 `SelectionSet` 레지스트리(§8.5)에서 inp, fem, bdf, rad, k로 **메쉬 + 모든 Set 키워드(Node/Element/Surface Set) + 노드 결과(U, NT11) + 요소 평균 응력**을 쓴다. IGA 파트는 proxy hex-tet 격자로 쓴다.
- **Set 키워드 출력 규약 (§5.10.2 표와 1:1 대응):**
  - OptiStruct/Nastran(`.fem`/`.bdf`): `SET <id>`(연속 ID는 `THRU` 압축, 8/16컬럼 및 free-field 지원), 표면 세트는 `SURF`/`BSURF` + `$$HMNAME SETS <id> "<name>"` 주석 보존.
  - Abaqus/CalculiX(`.inp`): `*NSET, NSET=<name>`, `*ELSET, ELSET=<name>`, `*SURFACE, TYPE=ELEMENT, NAME=<name>` (로컬 면 번호 `S1`–`S6`/`S1`–`S4` 자동 변환).
  - LS-DYNA(`.k`): `*SET_NODE_LIST_TITLE`, `*SET_SEGMENT_TITLE`, `*SET_SOLID_TITLE`.
  - Radioss(`.rad`): `/GRNOD/NODE`, `/GRBRIC`, `/SURF/SEG`.
- 용도: HyperMesh, ANSA, LS-PrePost, OptiStruct, Abaqus 등에서 **경계조건·접촉 면 재정의 없이** 즉시 재사용. 솔버 결과 바이너리 포맷은 범위 밖.
- 시험: 읽기→쓰기→읽기 왕복에서 노드·요소·모든 Set(이름, 종류, ID 목록) 무손실 일치.

---

### 5.11 선형 솔버 및 메모리 지능화 설정 (`HWT-IGA SOLVER`, 2026-11-15 실측·확장 반영)

#### 5.11.1 Kratos 경로 (`JobSpec.linear_solver`, `JobSpec.num_threads_mkl`, `JobSpec.out_mode`)

| 값 | Kratos 클래스 | 행렬 조건 | 용도 |
|---|---|---|---|
| `auto` (기본) | **실시간 시스템 가용 RAM + 피크 메모리 추정** 기반 자동 선택 (`pardiso_ldlt` / `pardiso_lu` / `amgcl`) | — | §5.11.3 지능형 전환 규약: 여유 메모리 이내이면 `pardiso_ldlt`(follower 시 `pardiso_lu`), 안전 예산 초과 시 **명시적 전환 메시지 출력 후 `amgcl` 반복 솔버로 자동 전환** |
| `pardiso_lu` | `LinearSolversApplication.PardisoLUSolver` | 일반 | follower 하중(비대칭 접선) 포함 모든 경우 |
| `pardiso_ldlt` | `LinearSolversApplication.PardisoLDLTSolver` | 대칭 | 선형 탄성 + penalty + MPC(`TᵀKT`); LU 대비 factor 메모리 약 절반 |
| `pardiso_llt` | `LinearSolversApplication.PardisoLLTSolver` | 대칭 양정치 | 위와 같음, 가장 빠름(실패 시 오류, 묵시적 폴백 금지) |
| `amgcl` | `KratosMultiphysics.AMGCLSolver` | 일반 / 대칭 | **Phase 3 편입**: 3D $p \ge 2$ 대규모 격자·어셈블리용 저메모리 반복 솔버 (`gmres`/`cg` + `amg`/`ilu0` precond, fill-in 메모리 0) |
| `cg` | `KratosMultiphysics.CGSolver` | 대칭 양정치 | 경량/디버그·사전조건 비교용 (`tolerance`, `max_iterations` 설정 가능) |
| `skyline_lu` | `KratosMultiphysics.SkylineLUFactorizationSolver` | 일반 | 소형 시험·디버그만 |

- **`JobSpec.out_mode` (`"drop"` | `"ghost"`, 기본 `"drop"` 추천)**:
  - `"drop"`: `OUT` 셀에 대해 Kratos 적분점/요소를 생성하지 않고 `INACTIVE` 제어점 자유도만 0 고정한다 (V06 L2 `p=2` 실측: `"ghost"` 대비 전체 DOF `222,264 → 62,151`, 72% 감소).
  - `"ghost"`: FCM 가상 강성(`ghost_alpha = 1e-8`) 요소 채움. Phase 1/2 기준선 비교 및 FCM 조건수 검증용으로 유지.
- **대칭성 근거**: `ResidualBasedBlockBuilderAndSolver::ApplyDirichletConditions`가 고정 행과 **고정 열을 함께** 0으로 만든다(`residualbased_block_builder_and_solver.h` l.876–905). 요소·penalty·`TᵀKT`는 대칭. follower 압력(`MoveMeshFlag`)은 비대칭이므로 `ldlt`/`llt`/`cg`를 명시적으로 고르면 `run_linear`가 오류를 낸다(`auto`는 자동으로 `pardiso_lu` 또는 비대칭 `amgcl` BiCGStab/GMRES 선택).
- **바꿀 수 없는 것**: Kratos Pardiso는 Eigen `PardisoSupport` 래퍼라 iparm이 C++에 고정(iparm2=2 METIS, iparm8=2 반복 보정, iparm10=13 pivot 섭동, iparm11/13은 LU에서만 1). Kratos 재빌드 금지이므로 조정 불가.
- **바꿀 수 있는 것**: 분해 종류(`pardiso_*` ↔ `amgcl`), `num_threads_mkl`(정수, `"consistent"`=Kratos 스레드 수, `None`=MKL 기본 = 물리 코어), `out_mode`, 반복 솔버 파라미터(`solver_type`, `preconditioner_type`, `tolerance`, `max_iterations`, `coarsening_type`, `krylov_type`).
- **선택 기준**: `scripts/bench_pardiso.py`(V02 레벨 3·4 실행렬, 분해 종류 × 스레드) 및 `AMGCL` 정확도 검증(`max|ΔU|/max|U| ≤ 1e-7`, 평형 오차 ≤ `1e-8`)으로 확인한다.
- **실측**: V02 레벨 1에서 Skyline과 Pardiso LU 결과가 출력 자릿수까지 동일(peak 4.673141e-06). V02 레벨 0–4가 Pardiso로 19분(이전 Skyline: 45분 이상, 레벨 4에서 16 GB). 반면 V06 L2 `p=2`(`40³`, `out_mode="ghost"`)에서 PARDISO 분해 시 프로세스 20 GB로 OOM 발생(`2026-11-15`).
- **실측 (2026-11-15 오후, G2 마무리 — §5.11.3 설계의 입력값)**:
  - `out_mode="drop"`은 Phase 2에서 **구현 완료**(`wht_iga/solver/active.py`, `JobSpec.out_mode`; 코드 기본값은 V02 증거와의 일치를 위해 `"ghost"`, V06 runner·V-COND는 `"drop"`). drop = ghost 결과(V06 L0 p1/p2·L1 p1/p2 유효숫자 3–4자리 일치), **V06 L2 p=2 peak 6.6 GB**(ghost는 20 GB 초과로 중단) → PASS.
  - `linear_solver="amgcl_cg"`(Kratos `AMGCLSolver`, CG + ILU0 + aggregation, tol 1e-10) 구현했으나 **V02 h=0.75(14만 DOF, ghost 1e-8E, penalty, MPC)에서 10분 넘게 결과 없이 중단**(LDLᵀ 33 s, 원인 미진단) → 현 상태로는 `auto`의 2단계 자동 전환 대상이 될 수 없다. 전환 규칙 도입 전에 drop + 대각 scaling + block(3×3) AMG 조합의 수렴 spike가 선결이다.
  - V02 h=0.375(54.6만 DOF, nz=1, p=2) LDLᵀ 프로세스 피크 ≥ 18 GB(여유 13 GB에서 감시 종료 2회) — 추정 모델 보정용 실측점. 같은 격자는 여유 메모리가 더 컸던 새벽에 173 s로 풀렸다.
  - **MKL 단일화**: Kratos(delay-load, PATH 순서 의존)와 pypardiso(`sys.prefix` glob)가 서로 다른 MKL·OpenMP를 올리면 `0xC06D007F` 조용한 종료 또는 `OMP: Error #15` 중단. `kratos_backend._bootstrap`이 env `mkl_rt.3.dll`을 전체 경로로 선로드하고 `PYPARDISO_MKL_RT`를 같은 파일로 지정한다(`WHT_IGA_MKL_DIR`로 변경).
- `run_linear` 결과 `timings = {setup_s, solve_s, linear_solver, linear_solver_requested, auto_switched, peak_rss_gb}`로 준비·풀이 시간과 실제 선택된 솔버를 분리 기록한다.

#### 5.11.2 Python 경로 (`wht_iga.linalg.PardisoFactor`, pypardiso)

`pypardiso.spsolve`는 모듈 전역 solver 하나를 공유하고, 항상 비대칭 경로(mtype 11, scaling + matching, 전체 LU)를 기본 iparm으로 쓰며, 호출마다 행렬 전체를 복사·비교한다. 대칭 행렬에는 손해이므로 래퍼를 쓴다.

| 설정 | 값 | 이유 |
|---|---|---|
| `kind` | `spd`(mtype 2, Cholesky) / `sym`(−2, LDLᵀ) / `general`(11, LU) / `auto`(대칭이면 sym) | 대칭이면 분해 메모리·연산 약 절반 |
| 입력 | 대칭 종류는 **상삼각만** 전달, 정렬된 CSR | MKL 규약 |
| iparm(1) | 1 | 아래 값 사용 |
| iparm(2) | 3 (병렬 nested dissection), 옵션 2 | 대형 행렬 ordering 병렬화 |
| iparm(8) | 0 | 섭동 pivot 시 자동 반복 보정 |
| iparm(10) | 8 (대칭) / 13 (일반) | pivot 섭동 크기 |
| iparm(11), iparm(13) | 0 (대칭) / 1 (일반) | scaling·matching은 비대칭에만 |
| iparm(21) | 1 (sym만) | Bunch–Kaufman pivoting |
| iparm(24) | 0, 옵션 1 | 2단계 병렬 분해(코어 많을 때) |
| iparm(35) | 0 | pypardiso가 1-기반 인덱스 전달 |

- 대칭 종류에 비대칭 행렬(상대 1e-12 초과)을 주면 오류, SPD 분해 실패도 오류 — 다른 종류로 몰래 바꾸지 않는다.
- 적용처: `reference/solve.solve_with_dirichlet`, `operators/compose.solve_reduced`(`kind="auto"`), `stabilize/cond.cond_eigsh`의 shift-invert(`OPinv` = 한 번 분해한 K⁻¹). 시험: `tests/unit/test_linalg_pardiso.py`(SuperLU 대비 1e-10, 13건).
- 한 MKL 핸들은 스레드 안전하지 않다 — 인스턴스를 스레드 간 공유하지 않는다. 분해 재사용은 `PardisoFactor.solve`를 여러 번 호출한다.

#### 5.11.3 실시간 시스템 메모리 기반 지능형 솔버 전환 규칙 (`wht_iga/solver/memory.py`)

고정 DOF 상수(예: 10만 DOF)로 "대규모"를 판정하지 않고, **실행 시점의 현재 PC 물리 메모리 여유량과 조립 행렬의 추정 피크 메모리를 비교해 지능적으로 결정**한다.

1. **추가 외부 의존성 없는 시스템 메모리 계측 (`probe_system_memory()`)**:
   - Windows: `ctypes.windll.kernel32.GlobalMemoryStatusEx`로 `total_ram_gb`, `avail_ram_gb`, `memory_load_pct` 조회.
   - Linux/POSIX: `/proc/meminfo` (`MemTotal`, `MemAvailable`) 조회 (`psutil`이 설치돼 있으면 프로세스 RSS 계측에 병행).
   - 환경변수 `WHT_SOLVER_AVAIL_RAM_GB`(선택)로 테스트 시 가용 메모리 한도를 결정론적으로 시뮬레이션 가능.
2. **솔버별 피크 메모리 추정 모델 (`estimate_solver_memory()`)**:
   - 활성 자유도 $N_{\text{act}}$, 추정 비영원소 수 $\text{nnz}(K) \approx N_{\text{act}} \cdot 3(2p+1)^3$ (또는 실제 Kratos 행렬 추출/조립 통계), 단면 프론트 폭 $F_{\text{eff}}$로부터:
     - **공통 오버헤드 $M_{\text{base}}$**: Kratos `ublas` $K$ + Eigen CSR 복사본 + Python 기하/테셀레이션 상주량 $\approx 2.2 \times 16\,\text{B} \cdot \text{nnz}(K) + M_{\text{rss,current}}$.
     - **직접 솔버(`pardiso_ldlt` / `pardiso_lu`) 피크 추정 $M_{\text{direct}}$**: fill-in 계수 $\alpha_{\text{fill}} \cdot N_{\text{act}} \cdot F_{\text{eff}} \cdot 12\,\text{B}$ (`pardiso_lu`는 `pardiso_ldlt`의 약 2배) + MKL 스레드별 프론트 작업 버퍼.
     - **반복 솔버(`amgcl`) 피크 추정 $M_{\text{iter}}$**: $K$ + AMG 계층 희소 행렬/전처리기 + Krylov 벡터 $\approx 4.5 \times 16\,\text{B} \cdot \text{nnz}(K) + M_{\text{rss,current}}$ (fill-in 없음).
3. **`linear_solver="auto"` 지능형 전환 및 묵시적 폴백 금지 원칙**:
   - **안전 메모리 예산**: `safe_budget_gb = 0.65 * avail_ram_gb` (OS·타 프로세스 보호용 35% 여유 확보; `memory_budget_ratio`로 조정 가능).
   - **결정 흐름**:
     1. **1단계 (MKL 스레드 조절)**: $M_{\text{direct}}(N_{\text{threads}}) > \text{safe\_budget\_gb}$ 이고 스레드 축소(`num_threads_mkl = min(4, n_threads)`)만으로 예산 내에 들어오면 스레드 수를 우선 제한.
     2. **2단계 (`AMGCLSolver` 자동 전환)**: 스레드 제한 후에도 $M_{\text{direct}} > \text{safe\_budget\_gb}$ 이거나, 직접 분해 예상 $\text{nnz}(L+U) > 1.8 \times 10^9$ (32비트 `MKL_INT` 상한 근접)이면 **대규모 모델로 판정하고 `AMGCLSolver`로 자동 전환**한다. (follower 하중 등 비대칭이면 `amgcl` `gmres`/`bicgstab`, 대칭이면 `cg` + `amg` / `spai0` 또는 `ilu0` 설정, 필요 시 `tolerance=1e-9, max_iterations=2000`).
     3. **명시적 전환 메시지 출력 (필수)**: 사용자가 명시적으로 `pardiso_lu` 등을 강제한 경우에도 메모리 위험 경고를 출력하고, `auto` 모드에서 전환이 발생하면 콘솔·로그·`ResolveReport`에 정량 근거를 포함한 전환 메시지를 반드시 출력한다:
        ```text
        [HWT-IGA SOLVER][AUTO-SOLVER] Large-scale problem detected relative to current system RAM:
          System RAM: total=31.6 GB, available=18.7 GB, safe_budget=12.16 GB (65%)
          Model size: active_DOFs=222,264 (total_DOFs=222,264), est_nnz(K)=6.67e+07
          Estimated peak: pardiso_ldlt=18.40 GB (> budget) vs amgcl=2.85 GB (<= budget)
          Decision: Switching linear solver 'pardiso_ldlt' -> 'amgcl' (CG + AMG preconditioner, tol=1e-9).
        ```
   - `extension` 상태 연계: `amgcl`로 전환되었는데 `extension=None`이고 절단 셀($\eta < 0.2$)이 존재해 조건수 악화로 반복법 수렴이 위협받는 경우, 전환 메시지에 `WARNING: extension=None with cut cells may slow AMGCL convergence; recommend extension theta>=0.2` 경고를 함께 남긴다.

---

### 5.12 `HWT-IGA SOLVER` 로깅 및 상용 수준 문제·진행 정보 출력 (Phase 3 편입)

솔버의 공식 명칭을 **`HWT-IGA SOLVER`**로 지정하고, 상용 CAE 솔버(Abaqus `.dat`/`.sta`, ANSYS MAPDL, OptiStruct `.out`) 수준의 문제 진단·솔버 정보·단계별 진행 상황·평형 검산·리소스 요약을 제공하는 전용 로깅 모듈(`wht_iga/solver/solver_log.py`)을 구현한다.

#### 5.12.1 On/Off 및 출력 대상 제어 (`JobSpec` / `StudySpec` / `AnalysisSession`)

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `enable_solver_log: bool` | `True` | 솔빙 호출 시 `HWT-IGA SOLVER` 상세 로그 생성 **On/Off 마스터 스위치** (`False` 또는 환경변수 `WHT_SOLVER_LOG=0` 시 로그 파일 생성·배너 출력을 완전히 생략하고 최소 오류만 전달) |
| `solver_log_file: str \| Path \| None` | `None` | 로그 파일 경로 (`None`이면 `output_dir` 지정 시 `<output_dir>/<job_name>.hwt.log` 자동 생성, `output_dir=None`이면 메모리 버퍼 또는 지정 경로만 기록) |
| `solver_log_stdout: bool` | `True` | `enable_solver_log=True`일 때 표준 출력(터미널/노트북) 실시간 스트리밍 여부 |
| `solver_log_level: Literal["minimal", "standard", "verbose"]` | `"standard"` | 출력 상세도 (`minimal`: 배너+요약+오류만, `standard`: 상용 표준 6단 정보, `verbose`: Kratos 내부 echo·반복 이력 매 스텝 포함) |

#### 5.12.2 상용 소프트웨어 수준 6단 출력 구조 (`solver_log.py`)

`enable_solver_log=True`로 솔빙을 호출하면 아래 6개 블록이 일관된 타임스탬프와 정렬 표 형식으로 기록된다:

1. **Header & Environment Banner**:
   - 공식 이름 **`HWT-IGA SOLVER`**, 패키지 버전, Git/브랜치 식별자, 시작 시각(ISO 8601), OS, CPU 물리/논리 코어 수, Kratos·MKL 버전, **호스트 총 RAM / 시작 시점 가용 RAM**.
2. **Problem & Discretization Summary (문제 정의 요약)**:
   - Job 이름, 단위계, 해석 유형(`Static` / `Modal` / `Transient` / `Thermal` / `NLGEOM`), 파트 및 인스턴스 수, 이산화 유형(`Immersed IGA` / `BodyFit` / `Native FE` / `Hybrid`).
   - 스플라인 차수 $p$ / 연속성 $C^{p-1}$, 배경 격자(`lo`, `hi`, `n_cells`, $h$), THB 레벨 분포(레벨별 셀 수), 셀 분류 통계(`IN` / `CUT` / `OUT` 비율, `out_mode`).
   - DOF 통계: 전체 제어점/노드 수, **`OUT` 제거(`drop`) 후 활성 DOF**, Dirichlet 고정 DOF, MPC(`LinearMasterSlaveConstraint`) Slave/Master DOF, 최종 자유방정식 수 $N_{\text{eq}}$.
   - 재료 물성 요약(`E`, `nu`, `rho`, 구성법칙), 하중·BC 요약 및 **사전 외력 합력 검산 $\sum F^{\text{ext}}, \sum M^{\text{ext}}$**.
3. **Immersed IGA & Assembly Diagnostics (침지 IGA 특화 진단)**:
   - CAD/테셀레이션 상태: 삼각형 수, `watertight` 여부, 체적 일치 오차.
   - 절단 셀 적분 요약(`S4`): `NNLS` 정밀 수렴 셀 수, `QUAD_RELAXED` 셀 수, `QUAD_FALLBACK` 셀 수(전체 및 임계 구역 비율), 적분 차수 $q$, 총 적분점 수.
   - 안정화·계층 요약(`S5`): `extension` 활성화 여부 및 `theta`, `MASTER`/`EXTENDED`/`INACTIVE` 기저 수, 최대 외삽 계수 `max|e_ij|`, `THIN_FEATURE` 경고 건수, 구속 그래프 검사 결과(연쇄 0 / 이중 slave 0 / $c=0$).
4. **Memory Audit & Linear Solver Decision (메모리 분석 및 솔버 선택)**:
   - §5.11.3의 `probe_system_memory()` + `estimate_solver_memory()` 결과 표(안전 예산 대비 `pardiso_*` vs `amgcl` 추정 피크 RAM).
   - `auto` 전환이 일어난 경우 **`[AUTO-SOLVER]` 전환 사유 및 메시지**를 강조 블록으로 출력.
   - 최종 결정된 선형 솔버, 전처리기, 수렴 허용오차, OpenMP/MKL 스레드 수.
5. **Stage Progress & Convergence Monitor (단계별 진행률 및 수렴 모니터)**:
   - `[STAGE 1/5]` ~ `[STAGE 5/5]` 경과 시간, 진행률(`10%..100%`), 현재 프로세스 RSS 메모리.
   - 비선형 Newton / 다단계 Step / 반복 솔버(`AMGCL`) 실행 시: Iteration/Increment 별 Residual Norm, Displacement Norm, Krylov 반복 횟수, 수렴 여부 및 Step 상태 요약 테이블(Abaqus `.sta` 유사 형식).
   - `on_progress` 콜백(`T-P3.6`)과 동일 이벤트를 공유해 GUI/CLI/로그가 일관되게 동기화됨.
6. **Post-Processing Verification & Resource Summary Footer (검증·종료 요약)**:
   - 전역 힘 평형 검증: 외력 합력 $\sum F^{\text{ext}}$ vs 반력 합력 $\sum R$, 평형 상대오차 $\|\sum F^{\text{ext}} + \sum R\| / \|\sum F^{\text{ext}}\|$.
   - 주요 응답 요약: 최대 변위 $\|\mathbf{u}\|_{\max}$, 최대 주응력 / von Mises 응력 피크 및 위치 좌표, 등록된 Probe 값 요약.
   - **Wall-Clock Time Breakdown** (`Geo/Classify`, `Quadrature`, `Extension/THB`, `Assembly`, `Solve`, `Post`, `Total`), **Peak Process RSS 메모리(GB)**, Warning/Error 총 건수, 최종 종료 상태(`COMPLETED` / `COMPLETED_WITH_WARNINGS` / `FAILED`).

---

## 6. 디렉토리·모듈 구조 (`WHT_IGA/`)

```
WHT_IGA/
├── pyproject.toml             # name=wht_iga, requires-python ">=3.12,<3.13"
│                              # deps: numpy, scipy, numba, pydantic>=2.7
│                              # extras: cad=[cadquery-ocp], viz=[pyvista, meshio, h5py], solve=[pypardiso],
│                              #         mesh=[pyNastran, meshio], test=[pytest, hypothesis]
│                              # gmsh: whtiga env의 별도 실행 파일 (import 금지, GPL 격리)
├── uv.lock
├── env/whtiga.yml             # micromamba spec (python=3.12, mkl, hdf5 …)
├── AGENTS.md                  # 다른 에이전트용 강제 규칙
├── CLAUDE.md                  # "@AGENTS.md"
├── README.md
├── docs/
│   ├── ENVIRONMENT.md         # 환경 단일 원본 (§12)
│   ├── KRATOS_BUILD.md        # Kratos 소스빌드 재현 절차 (§12.4)
│   ├── adr/                   # ADR-001 경로, ADR-002 penalty 기본값 …
│   ├── spikes/                # S01–S23 결과, TRACK_A/B.md, TRACK_A_REVIEW.md
│   └── verification/          # 케이스별 리포트 (자동 생성)
├── scripts/check_env.py       # Python 버전, Kratos 경로, import, numba 동작 검사
├── wht_iga/
│   ├── __init__.py            # 공개 독자 API (StudySpec, AnalysisSession, ResultView …)
│   ├── spec/                  # Layer 1: 불변 도메인 스펙 (Pydantic v2, frozen=True, JSON Schema)
│   │   ├── base.py            # SpecBase(model_dump_json/validate_json/with_updates/schema)
│   │   ├── mesh.py            # ImmersedGridSpec, BodyFitSpec, NativeMeshSpec, RefinementSpec
│   │   ├── assembly.py        # PartSpec, MaterialCellSpec, InstanceSpec, AssemblySpec (§5.9)
│   │   ├── physics.py         # MaterialSpec, StepSpec(Static/Modal/Implicit/Explicit/Thermal/Spectrum), BCSpec, LoadSpec, InteractionSpec
│   │   └── study.py           # StudySpec, ProbeSpec
│   ├── plan/                  # Layer 2: 실행 계획 & 사전 진단 (Command/Plan)
│   │   ├── artifacts.py       # ArtifactPaths (.wiga.sqlite, _results.sqlite, .vtkhdf, .spec.json)
│   │   └── compiler.py        # compile_plan(StudySpec) → ExecutionPlan + ResolveReport (W3)
│   ├── ports/                 # Layer 3: 추상 인터페이스 (Hexagonal Ports, typing.Protocol)
│   │   ├── mesh_engine.py     # MeshEnginePort
│   │   ├── solver_backend.py  # SolverBackendPort (Kratos / Reference / C++ 교체 가능)
│   │   └── result_store.py    # ResultStorePort (SQLite + VTKHDF 하이브리드)
│   ├── runner/                # Layer 4: 실행 오케스트레이터 & 생애주기
│   │   ├── session.py         # AnalysisSession (with 문, on_progress, cancel, cross_verify, HWT-IGA SOLVER 로그 연동)
│   │   ├── sweep.py           # ParametricSweep (with_updates 기반 배치 실행)
│   │   └── external.py        # 외부 프로세스 실행기 (Python 3.11 aieng 등에서 호출, §3.3)
│   ├── solver/                # 실행 파이프라인 · 메모리 지능화 · HWT-IGA SOLVER 로거 (§5.11, §5.12)
│   │   ├── memory.py          # 시스템 RAM 실측(GlobalMemoryStatusEx//proc/meminfo), 피크 메모리 추정, auto 솔버 전환 결정
│   │   ├── solver_log.py      # HWT-IGA SOLVER 상용 수준 6단 로거 (On/Off, 파일+stdout, 진행/수렴/평형/리소스)
│   │   ├── validate.py        # JobSpec (linear_solver, out_mode, enable_solver_log, solver_log_level 등)
│   │   └── pipeline.py        # run_linear (auto 지능 전환, out_mode="drop"|"ghost", 단계별 계측)
│   ├── query/                 # Layer 5: 통합 결과 쿼리 (CQRS, SQLite+VTKHDF 은닉)
│   │   └── results.py         # ResultView(summary, probe, field, reaction_force, modal_table, error_norm)
│   ├── compat/                # 외부 API 호환 퍼사드 (코어와 격리)
│   │   └── coreform.py        # Coreform IGA 호환 어댑터 (ImmersedRectilinearMesh, Model, Job, PointProbe → StudySpec/AnalysisSession 위임)
│   ├── geometry/              # S1
│   │   ├── occ_io.py          # OCP STEP/IGES/BREP 읽기, ShapeFix, Sewing/MakeSolid, 셀/면 이름
│   │   ├── tessellate.py      # BRepMesh_IncrementalMesh, face_id 보존, 정점 병합
│   │   ├── mesh_io.py         # STL/OBJ 입력
│   │   ├── watertight.py      # 매니폴드, 오일러수, 체적 부호, 틈 검사
│   │   ├── selectors.py       # (하위 유틸) 이름/박스/평면/법선 기본 필터
│   │   └── primitives.py      # 검증 형상 생성기 (Roark 판, 구 octant, 보, Cook, Hertz …)
│   ├── selection/             # CAD+Mesh 통합 선택기·계측·접촉 점검 (§8.5)
│   │   ├── selector.py        # EntityQuery(within_box/sphere/cylinder, on_plane, by_normal, flood_fill, | & -)
│   │   ├── sets.py            # SelectionSet (CAD @face/@edge/@vertex ↔ FEM Node/Elem/Surface Set 브릿지)
│   │   ├── measure.py         # bbox, area, volume, length, centroid, inertia_tensor, distance_to
│   │   └── contact_check.py   # inspect_contact(): 침투 깊이, 이격(gap), 법선 정렬, search_radius 추천
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
│   ├── meshio_ext/            # 메쉬 입출력 (§5.10)
│   │   ├── model.py           # MeshModel, 요소 유형 정규화, 노드 순서 변환표
│   │   ├── read_inp.py, read_rad.py, read_k.py      # 자체 경량 파서
│   │   ├── read_nastran.py    # bdf/fem (pyNastran + $$HMNAME 주석)
│   │   ├── read_msh.py        # meshio
│   │   ├── write_*.py         # inp, fem, bdf, rad, k (메쉬 + 노드 결과 + 요소 평균 응력)
│   │   ├── skin.py            # 솔리드 메쉬 외곽 면 → tessellation (immersed 입력)
│   │   └── bodyfit.py         # 정렬 hex 블록 감지 → NurbsVolume 제어점 망
│   ├── meshgen/
│   │   └── gmsh_runner.py     # Gmsh 외부 프로세스 호출 → .msh
│   ├── bus/                   # SQLite 스키마와 입출력
│   │   ├── schema_mesh.sql, schema_results.sql, schema_quad.sql
│   │   └── blob.py, mesh_db.py, results_db.py, quad_db.py
│   ├── kratos_backend/        # S7–S8 (Kratos는 여기서만 import)
│   │   ├── _bootstrap.py      # 경로·DLL 설정, 빌드 위치 검증
│   │   ├── geometry.py        # make_nurbs_volume, cell_gauss_points (공유, Phase 0 작성)
│   │   ├── elements.py        # 셀 적분점 기하 + SM 요소 부착, K 추출 (트랙 A)
│   │   ├── constraints.py     # LinearMasterSlaveConstraint 배치 생성, c=0 안전장치 (트랙 A)
│   │   ├── surface.py, thermal.py   # 경계 적분점 기하·조건, 열 (트랙 B)
│   │   ├── assembly_builder.py # ModelPart 계층(§5.9.2), instance 노드 생성, 노드 id 오프셋 표
│   │   ├── interactions.py    # tie(MPC/mortar), contact 프로세스, RBE2/RBE3 결합 (§5.9.4)
│   │   ├── native_fe.py       # 일반 FE 파트 → SM 일반 요소
│   │   ├── loads.py, materials.py (구성법칙 매핑, Accessor/요소별 Properties)
│   │   ├── mass.py            # Python consistent 질량 행렬 (트랙 A F1)
│   │   ├── solvers.py         # HasMKL 분기 + AMGCLSolver/CGSolver 반복 솔버 팩토리 + auto 전환 실행 (§5.11.1, §5.11.3)
│   │   ├── analyses/          # static.py, nonlinear.py, modal.py(scipy eigsh), bossak_linear.py, thermal.py
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
│   └── cli.py                 # python -m wht_iga {mesh, run, post, verify}
├── cpp/                       # (Phase 5+) WhtIgaApplication
├── tests/
│   ├── unit/                  # Kratos 불필요
│   ├── kratos_spikes/         # S01–S25 (marker: kratos)
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

- OCP `STEPControl_Reader` / `IGESControl_Reader` / `BRepTools` → `ShapeUpgrade_UnifySameDomain` → `BRepBuilderAPI_Sewing`(허용오차 = 1e-6 × bbox 대각) → (IGES 등 열린 셸인 경우) `BRepBuilderAPI_MakeSolid` + `ShapeFix_Solid`. IGES는 단위 메타데이터를 검증해 `unit_system`과 불일치 시 경고/변환하며, 봉합 후에도 비수밀인 경우 S3 Winding Fallback으로 연계한다.
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
- 계수 e_ij = Π_ν w^ν_{i_ν}(j_ν). 1D 가중치 w는 **실제 knot 벡터에서** 블록 크기 미만 차수의 모든 단항식 B-spline 계수를 재생하는 연립식 Σ_m w_m c_k(lo+m) = c_k(j)로 구한다. c_k(i) = e_k(t_{i+1}, …, t_{i+p}) / C(p, k) (Marsden 항등식의 blossom, e_k = 기본 대칭 다항식). 균일 내부 knot에서는 인덱스 공간 Lagrange 외삽(Höllig 2003)과 같다. p차 다항식 재생과 분할단위성이 보존된다(`stabilize/extension.extrapolation_weights_1d`, `build_extension(..., knots=)` 필수).
- **정정 (2026-11-14)**: 인덱스 공간 Lagrange 외삽은 clamped 경계 근처(p≥2)에서 Greville 점이 비균일해 1차 함수도 재생하지 못한다(slave 계수 오차 p=2 0.125, p=3 0.083). 재생 시험은 다항식을 **물리 좌표**에서 만들고 계수는 Greville collocation 등 독립 오라클로 구해야 한다(인덱스 공간 시험은 by-construction).
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

- 질량 행렬 (v2.1, 트랙 A F1): 적분점 기하에서 SM 요소의 consistent 질량은 0이 된다. **consistent M은 Python(Numba)에서 조립**한다(셀 적분점 × Nᵢ Nⱼ ρ, 절단 셀은 같은 모멘트 피팅 점 사용, 구속식 T는 Tᵀ M T로 적용). Kratos lumped 질량은 정상이며, V09 보에서 1차 모드 차이 0.07%(p=2, 40×4×4)·0.02%(80×8×8), 상위 20모드 최대 7.9%→2.3%로 측정됐다.
- 고유치: Kratos 추출 K + Python consistent M → scipy `eigsh`(shift-invert, σ 지정). 질량 정규화를 확인한다. (FEAST 미빌드, spectra는 shift 0에서 실패 — 트랙 A S06)
- 동적 시간적분은 **Bossak**을 쓴다 (사용자 결정 D5: Kratos 적분법을 HHT의 상위 호환으로 채택). 선형 과도는 같은 행렬로 Python Bossak 루프, 비선형은 Kratos Bossak scheme + lumped M.
  - 기본: Bossak(WBZ-α), `damp_factor_m` α_B ∈ [−1/3, 0], 기본 −0.3. β = (1−α_B)²/4, γ = 1/2 − α_B로 2차 정확도와 무조건 안정을 유지하면서 고주파 감쇠를 조절한다.
  - 감쇠 없는 비교가 필요하면 Newmark(α_B = 0), 강한 고주파 감쇠가 필요하면 BDF2.
  - HHT는 구현하지 않는다. Coreform 검증(V10)은 HHT 결과 자체가 아니라 물리량 기준(f1, 위상, 응답 이력 RMS)으로 비교한다.
- 초기 가속도 계산, amplitude, restart(상태 저장)는 Coreform 동작을 따른다.
- **Explicit 동역학 (Phase 6)**: Kratos `StructuralMechanicsApplication`의 `ExplicitCentralDifferencesScheme` / `MechanicalExplicitStrategy` 또는 Python Explicit 중앙차분 루프를 사용한다. IGA Immersed 환경의 2대 핵심 난제를 다음 전략으로 해결한다:
  1. **질량 행렬 집중화 (lumping)**: B-스플라인 기저는 음수가 없으므로 row-sum lumping이 **항상 양의 대각 질량**을 주고, IGA explicit에서 표준으로 쓰인다(정확도는 2차로 낮아짐; Cottrell et al. 2006은 row-sum lumping의 정확도 분석). 기본은 **Kratos lumped 질량**(적분점 기하에서 정상 동작, 트랙 A 실측: V09 보 1차 모드 consistent 대비 0.07%, 상위 20모드 최대 7.9%→2.3%). HRZ는 S24에서 비교 옵션으로 측정한다.
  2. **소형 절단 셀의 임계 시간 간격**: 작은 절단 셀은 질량이 거의 0이라 최대 고유진동수가 치솟는다. extension(지지 작은 기저 제거, §7.5)이 이를 막아 Δt_crit ~ h/c를 유지하는지 S24로 확인하고, 부족하면 selective mass scaling을 쓴다(참고: Stoter et al. 2023, Leidinger 2020). extension 구속식은 Kratos `MechanicalExplicitStrategy`가 질량·잔차를 master로 모아 처리한다(§5.6). Python explicit 루프는 S24 실패 시 대안.
  3. **Explicit 접촉**: ContactSM `explicit_penalty_contact_process`를 표면층(§5.9.4와 같은 스프링 결합)에 적용한다(S25).

### 7.10 응답 스펙트럼 / 랜덤 응답 / 모드 중첩 동역학

- 참여계수 $\Gamma_{ik} = \phi_i^T f_k$ ($f_k = \int \rho N e_k$), 유효질량 $\Gamma^2$.
- **모드 중첩 과도 / 조화 응답 (Phase 4)**: 추출된 고유모드 $\phi_i$, 고유진동수 $\omega_i$, 모드 감쇠비 $\zeta_i$에 대해 일반 하중 $f(t)$ 또는 $F(\omega)$를 모드 투영하여 SDOF 방정식($\ddot{q}_i + 2\zeta_i\omega_i\dot{q}_i + \omega_i^2 q_i = \phi_i^T f(t)$)을 적분·합성($u \approx \sum \phi_i q_i$)하고, 응력은 복원된 변위 $u$를 노드에 주입해 `CalculateOnIntegrationPoints`로 산출한다.
  - 모드 절단 오차 보정: 하중 공간의 정적 보정(static correction / residual mode)을 선택 옵션으로 둔다. V10 직접법(Bossak)과 비교할 때 차이 원인을 분리하기 위함.
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
- 온도 의존 물성(S20 실측): $E(T)$ 단독 의존은 `TableAccessor` + `CL.ThermalElasticIsotropic3D`(`6.013e-09` 통과) 사용, $\alpha(T)$·`LinearElastic3DLaw`는 `TableAccessor`를 무시하므로 셀별 `Properties` 복제($O(h^2)$ 수렴, `1.049e-05` @ `n_x=32`)를 채택한다.

### 7.12 후처리와 시각화

- 장 평가: 임의 점 x에서 셀 찾기(항등 사상이라 인덱스 계산만) → 기저 평가 → 제어점 계수와 곱 (@njit).
- 응력: 적분점 값 직접 사용(probe = 가중치 0 적분점, Coreform 방식). 비선형 상태 변수도 정확히 얻는다.
- proxy 격자: IN 셀은 hex(서브디비전 선택), CUT 셀은 셀 내 클리핑 다면체를 tet으로 분할. 값은 장 평가로 채운다.
- 형식: VTKHDF(h5py, 시계열) 또는 VTU(meshio). ParaView에서 바로 연다.
- L2/H1/에너지 노름: quadrature DB의 X, dV, U, S로 계산한다.

---

## 8. 사용자 Python API (독자적 5계층 Hexagonal + Pydantic v2, D18)

Coreform API를 단순히 복제하지 않고, **불변 도메인 스펙(Pydantic v2) · 사전 계획(Plan) · 포트/어댑터(Hexagonal) · 생애주기 세션 · 통합 쿼리(CQRS)**로 독자 설계한다. Coreform 스크립트 호환은 코어와 격리된 `wht_iga.compat.coreform` 퍼사드로 제공한다.

### 8.1 계층 원칙
1. **Immutable Spec (`wht_iga.spec`):** Pydantic v2 `frozen=True` 모델. 생성 즉시 타입·물리 범위·단위계를 검증하고, `model_json_schema()`로 AI/MCP/UI 자동 연동, `with_updates(...)`로 원본 훼손 없는 파라메트릭 복제를 지원한다.
2. **Pre-flight Plan (`wht_iga.plan`):** `compile_plan(spec)`이 솔버 실행 전 `ResolveReport`(DOF·실시간 시스템 RAM 대비 예상 피크 메모리·`auto` 솔버 전환 예측·$\Delta t_{\text{crit}}$·미지원 조합 진단)와 `ArtifactPaths`를 산출한다.
3. **Ports & Adapters (`wht_iga.ports` / `backends`):** `SolverBackendPort` 프로토콜을 통해 동일 `StudySpec`을 `KratosPythonBackend`, `ReferenceNumbaBackend`(Oracle), `KratosCxxBackend`로 코드 수정 없이 교체·교차검증한다.
4. **Lifecycle Session (`wht_iga.runner`):** `with AnalysisSession(..., enable_solver_log=True)` 컨텍스트 매니저, **`HWT-IGA SOLVER` 상용 수준 6단 로그 출력(On/Off 가능)**, 진행 콜백, 취소 시 부분 산출물 정리를 보장한다.
5. **Hybrid Query (`wht_iga.query.ResultView`):** `SQLite`(메타·요약·Probe)와 `VTKHDF`(대용량 필드·시계열)를 은닉한 단일 조회 API.

### 8.2 기본 사용 예 (독자 API)

```python
import wht_iga as wi

study = wi.StudySpec(
    name="roark",
    units="mm-N-tonne-s",
    parts=[wi.PartSpec.from_step("plate_hole.step", material="steel", heal=True)],
    materials={"steel": wi.LinearElasticSpec(E=2.0e5, nu=0.3, density=7.85e-9)},
    mesh=wi.ImmersedGridSpec(
        degree=2, continuity=1, element_size=2.5, tessellation_size=0.25,
        bc_penalty=1e10, small_cell_volume_ratio=0.2, padding=(2, 2, 2),
        refinements=[wi.RefinementSpec(name="hole_edge", level=2, sets=["hole"], balance="two_to_one")],
    ),
    steps=[
        wi.StaticStepSpec(
            name="load",
            bcs=[wi.DisplacementBCSpec(set="sym_x", ux=0.0), wi.DisplacementBCSpec(set="sym_y", uy=0.0)],
            loads=[wi.PressureLoadSpec(set="right", magnitude=-1.0)],
        )
    ],
    probes=[wi.PointProbeSpec(name="peak", location=(10.0, 0.0, 0.0), variables=("S", "U"))],
)

plan = wi.compile_plan(study)              # 사전 진단(ResolveReport) 및 아티팩트 경로 확정
assert plan.report.is_valid

with wi.AnalysisSession(plan, backend=wi.KratosPythonBackend(threads=16)) as session:
    session.build_mesh()                   # S1–S6 → roark.wiga.sqlite
    results = session.run()                # S7–S9 → SQLite(메타/Probe) + VTKHDF(필드)

peak_s = results.probe("peak")["S"]
s_field = results.field("S", step="load")  # VTKHDF에서 고속 슬라이싱

# 파라메트릭 스윕 및 Oracle 교차검증 (불변 스펙 복제)
fine_study = study.with_updates(mesh=study.mesh.with_updates(element_size=1.25))
json_str = study.model_dump_json()         # Pydantic v2 직렬화 / model_json_schema() 지원
```

### 8.3 Coreform 호환 퍼사드 매핑 (`wht_iga.compat.coreform`)

외부 Coreform 스크립트 중 **포터블 spec 계층**(`coreform.mesh` / `coreform.probe` / `coreform.job`의 spec과 `to_json/from_json/validate`)과 같은 이름·인자를 `import wht_iga.compat.coreform as cf`로 받는다. Abaqus 의존 부분(`mdb`, `part.customData.CoreformIGA`, `coreform.abaqus.init_iga_mdb()`, CAE 기록)은 지원하지 않는다. 내부적으로는 5계층(`StudySpec` → `compile_plan` → `AnalysisSession` → `ResultView`)으로 변환·위임한다.

| Coreform API (`compat.coreform`) | WHT_IGA 코어 대응 | 비고 |
|---|---|---|
| `ImmersedRectilinearMesh(...)`, `BodyFitMesh(...)` | `wi.ImmersedGridSpec`, `wi.BodyFitSpec` | Pydantic v2 검증 + `with_updates` |
| `Refinement`, `CadEntitiesRegion`, `SphereRegion`, `BoxRegion` | `wi.RefinementSpec` + region union | JSON Schema 자동 생성 |
| `PointProbe`, `LineProbe`, `MultiPointProbe` | `wi.PointProbeSpec`, `LineProbeSpec`, `MultiPointProbeSpec` | `ResultView.probe()`로 통합 조회 |
| `Job(build_mesh, write_input, submit, wait_for_completion, kill, status)` | `wi.compile_plan` + `wi.AnalysisSession` | Context Manager + `SolverBackendPort` 분리 |
| `to_dict / to_json / from_json / copy_with / validate` | `model_dump / model_dump_json / model_validate_json / with_updates` | Pydantic v2 네이티브 |
| `job.results.probe(...)` | `wi.ResultView` (`probe`, `field`, `modal_table`, `reaction_force`, `error_norm`) | `SQLite + VTKHDF` 하이브리드(D17) |

### 8.4 Kratos-Native 친화 설계 (과도한 추상화 방지)

추상 계층(`SolverBackendPort`)이 Kratos의 고유 기능을 가리지 않도록, `KratosPythonBackend`는 **Kratos-First Native Ergonomics**를 기본 제공한다:

1. **네이티브 객체 직접 노출 (Escape Hatch):** `session.backend.kratos_model` (`Kratos.Model`), `session.backend.root_model_part` (`Kratos.ModelPart`), `session.backend.strategy`에 직접 접근 가능.
2. **Kratos `Process` / `Parameters` 직결 훅:** `KratosPythonBackend(..., custom_processes=[...], kratos_parameters=...)` 또는 `on_before_solve(model_part)` / `on_step_end(model_part, step)` 훅으로 Kratos 기존 Process·Utility·Mapper를 그대로 주입.
3. **1:1 개념·변수명 공용:**
   - 어셈블리 파트·인스턴스·경계 세트는 `ModelPart`/`SubModelPart` 계층 이름(§5.9.2)과 동일하게 매핑.
   - `ResultView`와 `Probe`는 Coreform 별칭(`U`, `S`, `E`, `NT11`)과 Kratos 네이티브 변수명(`DISPLACEMENT`, `CAUCHY_STRESS_TENSOR`, `VON_MISES_STRESS`, `TEMPERATURE`)을 모두 지원.
4. **Kratos 솔버·구성법칙 Passthrough:** `LinearSolver`, `Scheme`(Bossak/Newmark/BDF2/ExplicitCentralDifferences), `ConstitutiveLaw`/`ConstitutiveModel` 이름을 스펙 내 `kratos_overrides`로 직접 지정 가능.

### 8.5 통합 기하·메쉬 선택, 계측 및 접촉 사전 점검 API (`wht_iga.selection`)

상용/오픈소스 CAE·CAD 프레임워크(Abaqus `getByBounding*`/`findAt`, HyperMesh `by_angle`/인접 전파, CadQuery/build123d Fluent 필터, PyVista/SciPy `KDTree`/`BVH`)의 강점을 통합하여, **CAD(STEP/IGES/BREP)와 Mesh(`MeshModel`/표면층)를 단일 Fluent API로 다룬다.**

1. **이중 도메인 통합 선택 (`EntityQuery`):**
   - 대상: `faces()`, `edges()`, `vertices()`, `solids()` (CAD 위상) 및 `nodes()`, `elements()`, `surface_tris()`, `cells()` (FEM/IGA).
   - **자유 표면(Free Surface) 자동 추출:** 솔리드 메쉬에서 내부 공유 면을 소거한 `free_surface_faces()`, 외곽 요소 `free_surface_elements()`, 외곽 노드 `free_surface_nodes()` 지원 (IGA에서는 경계 표면층·미구속 자유면 필터링, `ALL EXTERIOR` 접촉 후보 추출에 활용).
   - 공간/조건 필터: `within_box(lo, hi)`, `within_sphere(center, r)`, `within_cylinder(p0, p1, r)`, `on_plane(origin, normal, tol)`, `between_planes(...)`, `by_normal(dir, tol_deg)`, `by_area(min, max)`, `by_length(...)`, `nearest_to(pt)`, `flood_fill(angle_deg=30)` (특징각 인접 전파).
2. **CAD ↔ FEM Set 양방향 브릿지 & Set 불리안 대수 (`SelectionSet`):**
   - CAD/Query 결과를 `.to_node_set("bolt_hole")`, `.to_element_set(...)`, `.to_surface_set(...)`로 즉시 등록.
   - **불리안 대수 (Set Algebra):** 이름/객체 기반 Set 간 합집합(`|`, `union`), 교집합(`&`, `intersection`), 차집합(`-`, `difference`), 대칭차집합(`^`, `symmetric_difference`), 여집합(`~`, `complement`) 지원. 이종 변환 후 연산(예: `free_surf - tied_surf`, `surf_set.to_nodes() - fixed_nodes`) 및 파생 Set 실체화(`materialize()`) 지원.
   - 등록/파생된 Set은 **OptiStruct(`SET`/`BSURF`), Abaqus/CalculiX(`*NSET`/`*ELSET`/`*SURFACE`), LS-DYNA(`*SET_*`), Radioss(`/GRNOD`/`/SURF`)** Export 시 네이티브 키워드로 자동 출력(§5.10.6).
3. **기하 계측 (`measure`):**
   - `bbox()`(`min`, `max`, `center`, `extents`, `diagonal`), `length()`, `area()`, `volume()`, `centroid()`, `inertia_tensor()`, `distance_to(other)`(최단거리 및 최근접 쌍).
4. **접촉 상태 사전 점검 (`inspect_contact`) & 자동 계면 탐지 (`detect_interfaces`):**
   - 해석 전 두 표면(`master_surf`, `slave_surf`) 간 `KDTree`/`BVH` + 삼각형 투영으로 **최대 침투 깊이(`max_penetration`)**, **최소/최대 이격(`min_clearance`, `max_gap`)**, **법선 대향각(`normal_opposition_deg`)**, **예상 접촉 면적 비율**, **권장 검색 반경(`suggested_search_radius`)**을 `ContactGapReport`로 반환한다.
   - `assy.detect_interfaces(max_gap, max_angle_deg)`: 멀티파트 어셈블리에서 인접·접촉 후보 면 쌍을 자동 탐지해 `Tie` / `Contact` 후보를 제안한다.
5. **국부 좌표계(`Frame`) & 사전 물성·하중 검산:**
   - `Frame.from_3points(origin, x_pt, xy_pt)`, `Frame.cylindrical(origin, axis)`: 원주/반경 방향 BC(`ur`, `utheta`), 국부 응력/변위 변환 조회(`results.field("S", frame=csys)`), `CORD2R`/`*TRANSFORM` Export와 1:1 대응.
   - `model.integrate(load, over=surf_set)`(총 하중 검산) 및 `assy.mass_properties()`(총 질량·무게중심·관성모멘트 사전 확인).
6. **공학 검증용 후처리 (`ResultView` 고수준 명령):**
   - `check_equilibrium(step)`: 외력 합력 $\sum F^{\text{ext}}$ vs 반력 합력 $\sum R$ 자동 평형 점검.
   - `section_force(origin, normal)`: 임의 절단면(Section Cut)을 통과하는 합력·모멘트 적분.
   - `hotspots(var, n=5, min_distance=...)`: 특이점 편향을 피해 일정 반경 이상 떨어진 국소 피크 Top-N 좌표·값 추출.
   - `stress_linearize(p0, p1)`: 두 점 경로를 따라 Membrane + Bending + Peak 응력 자동 분해.
7. **OCP 기하 타입(`GeomType`) 자동 판별 & 모델링 가이드 연결고리 (`BRepAdaptor` 연동):**
   - OpenCASCADE(`OCP.BRepAdaptor_Surface` / `BRepAdaptor_Curve`, `GeomAbs_SurfaceType` / `GeomAbs_CurveType`)를 직접 활용해 기하 타입별 선택 및 파라미터 추출을 지원한다.
   - **선택 필터:** `faces().by_geom_type("cylinder", radius=(4.8, 5.2))`, `"plane"`, `"sphere"`, `"cone"`, `"torus"`, `"bspline"` 및 `edges().by_geom_type("circle" | "line" | "bspline")` (예: 반지름 조건만으로 모든 볼트 홀·필렛 자동 선택).
   - **자동 가이드 연결 (`occ_io` → `plan`/`selection`):**
     - `Cylinder` 감지 → 원통 좌표계(`Frame.cylindrical`) 자동 생성 및 핀/볼트 구속(`ur=0`) 가이드.
     - 소형 `Cylinder`/`Torus`(필렛·홀 반경 $r$) 및 `BRepLProp_SLProps` 고곡률 `BSpline` 감지 → `tessellation_size ≤ r/4` 및 국소 THB 세분 영역(`CadEntitiesRegion`) 자동 추천.
     - 6면 모두 `Plane`인 정렬 블록 감지 → `BodyFitMesh` 또는 객체 정렬 `Frame` 추천.

```python
part = wi.load_cad("bracket.igs", heal=True, unit_system="mm-N-tonne-s")   # STEP/IGES 통합
hole_faces = part.faces().within_cylinder((0, 0, 0), (0, 0, 20), r=5.2).by_normal(perp=(0, 0, 1))
load_nodes = part.nodes().within_box((-10, 49.9, 0), (10, 50.1, 20))

part.define_set("bolt_hole", hole_faces.to_surface_set())   # Export 시 *SURFACE / BSURF 출력
free_outer = part.free_surface_faces() - part.sets["bolt_hole"]  # 자유 외곽면 불리안 차집합
part.define_set("contact_ext", free_outer)                  # 파생 Set 등록 및 Export 반영
box = hole_faces.measure().bbox()                           # extents, diagonal, area 등
gap = wi.inspect_contact(master=assy["pin"].faces(), slave=hole_faces)
print(gap.max_penetration, gap.min_clearance, gap.suggested_search_radius)
```

---

## 9. 단계별 로드맵

기간은 1인 기준 추정이다. 각 Phase는 게이트를 통과해야 다음으로 넘어간다.

### Phase 0 — 환경과 Kratos API spike (2주)

- **산출물**
  - micromamba `whtiga` env, `pyproject.toml` + `uv.lock`, `check_env.py`, `docs/ENVIRONMENT.md`, `docs/KRATOS_BUILD.md`
  - `WHT_IGA` 골격 (bootstrap, isolation test)
  - S01–S18 spike 시험과 리포트, ADR-001 (R1/R2/R3)
  - 참조 조립기 최소판 (정육면체 untrimmed 선형 탄성)
  - (v2.1 추가) S19–S24: ConstitutiveModels 재료, Accessor, 구속 상수 Newton, 표면층 mesh tying/접촉, 비정합 계면 가중치, Explicit 질량/CFL
- **게이트 G0**
  - S01–S05 결론 확정
  - untrimmed 큐브 p = 1, 2, 3에서 Kratos 경로와 참조 조립기 변위 차이 ≤ 1e-10
- **위험**: S05 성능 → R2
- **진행 현황 (2026-11-14)**: 환경·골격·참조 조립기 완료, 트랙 A(S02–S06, S13, S15, S16) 완료·리뷰 수락, **트랙 C(S19–S24) 40/40 통과**, 트랙 B 진행 중.

### Phase 1 — 선형 정적, trimmed, penalty, extension 없음 (5주)

- **산출물**: S1–S4, S6, S7(BC-A/B, LoadCondition, 중력), 선형 전략, mesh/results/quadrature DB, VTU 최소 출력, FCM α 옵션, spec/validate
  - (v2.1) Part/MaterialCell/Instance/Assembly **데이터 모델**(단일 파트·단일 instance로 사용, P-1), `compile_plan()` 골격과 진단 리포트(P-8), 메쉬 미리보기(P-4), 단위계 가드(add_1 G-6), `constraints.add_linear_tie` c≠0 거부(F2)
- **게이트 G1**
  - **V01** 지정 변위 큐브 u = (x, 2y, 3z): 모든 p ∈ {1,2,3}, C에서 최대 상대오차 ≤ 1e-5 (Coreform 수준)
  - **V02** Roark (Coreform 원문 판정 방식, 사용자 결정 2026-11-15): 격자 h = 6/2^k (k=0–4, 2배 중첩 세분). ① p = 2 최세 격자 peak(**점 A = 정확한 원 위 (0, ±10, 0)의 최대 주응력**, 원문 Figure 1의 평가점. 테셀레이션 면 위 최대값은 내접 다각형 외삽으로 과대 → 참고값)의 Roark(4.3175e-6) 대비 상대오차 ≤ 1%, body-fitted Abaqus 기준(4.34797e-6) 대비 오차는 함께 보고 (Coreform p=2 h=0.375: 0.83% / 0.12%). ② **수렴률 ≥ 1.8은 peak가 아니라 hole edge→free edge 최단 구간(x = 0, y = 10→20)의 최대 주응력 분포를 최세 IGA 해에 대해 잰 상대 L2 자기 수렴률**(Coreform 원문 §7.5: 1.86/2.23/2.38). 근거: Roark와 수렴 해가 약 0.7% 달라 peak의 Roark 대비 오차는 세분할수록 0.7% 근처로 수렴하므로 기울기 판정이 성립하지 않음(Coreform p=2 peak 대비 Roark 오차도 24.6/5.66/0.12/0.23/0.83%로 비단조). 구현 `scripts/v02_roark.py`, `wht_iga/solver/probes.py`
  - **V03** 보 휨(moment 제외 5종): 끝단 처짐 Euler–Bernoulli 대비 5% 이내
  - 모멘트 피팅 다항식 정확도 1e-12
- **위험**: 소형 셀 조건수(Phase 2 전까지 FCM α로 완화), tessellation 기하 포화

### Phase 2 — small-cell extension과 조건수 (3주)

- **산출물**: S5 extension, `operators/compose`, 구속식 배치 생성, 조건수 도구, THIN_FEATURE 처리
  - (v2.1) 메쉬 리더 6종 + `MeshModel`(P-2), 솔리드 메쉬 외곽 면 → immersed 형상 입력(P-3). 게이트 추가: V02를 메쉬 파일 입력(외곽 면)으로도 Roark 1% 이내
- **게이트 G2**
  - **V-COND**(자체): **V02 Roark 판(구멍 중심 이동)**을 모체 형상으로 하여 η_min ∈ [1e-8, 1]을 만든다. κ 변동 ≤ 10배(평탄), h에 대한 log-log 기울기 −2 ± 0.3. 미안정화 대비 개선을 그래프로 제시
  - **V06** 압력 구 octant (a = 0.9, b = 1, Lamé): `pressurized-sphere-qoi.csv`의 p별·level별 허용오차 열 전부 통과. 목표 p = 2 최세 u L2 ≤ 1e-4, σrr ≤ 0.5% (Coreform 0.00365% / 0.211%)
  - V01, V02를 extension 켠 상태로 재통과
  - extension 분할단위성·재생 단위시험
- **위험**: 얇은 형상에서 블록 부재, 구속식 성능

### Phase 3 — 다중재료, 멀티파트·어셈블리, IGA–FE 혼합, 메쉬 입출력, coupling, probe, 세분화, 출력 (8주)

- **산출물**: 멀티셀 공유 장, tie(BC-C), RP/분배/kinematic coupling, amplitude, 다단계, 가중치 0 probe, THB 세분화, VTKHDF proxy 격자, probe JSON, L2 노름 도구
  - (v2.1) **멀티파트·어셈블리(§5.9)**: ModelPart 계층, dependent/independent instance, tie 자동 선택(항등 MPC/가중 MPC/합성 MPC/penalty/mortar), RBE2·RBE3 결합, 구속 그래프 검사
  - (v2.1) **일반 FE 파트 + IGA–FE 혼합 Beta**(P-6, D14), Gmsh 외부 프로세스 메쉬 생성(D13), body-fit 변환(P-5), 메쉬 쓰기 5종(P-7, **솔버별 Set 키워드 `SET`/`BSURF`, `*NSET`/`*ELSET`/`*SURFACE` 등 포함**)
  - (v2.1 확장) **IGES(`.igs/.iges`) 입력** 통합(`occ_io.py` + Sewing/MakeSolid) 및 **`wht_iga.selection`** (CAD+Mesh 통합 선택기, 기하 계측 `bbox`/`area`/`volume`, 접촉 사전 점검 `inspect_contact`)
  - (v2.1) add_1 G-2 진행 콜백·취소, G-5 THIN_FEATURE 자동 세분, G-7 출력 이름 매핑표
  - (v2.1 확장, 2026-11-15) **`HWT-IGA SOLVER` 상용 수준 6단 로깅**(On/Off, 문제·기하·적분·확장·솔버·진행·평형·리소스 출력, §5.12), **실시간 PC 가용 RAM 기반 지능형 `auto` 솔버 전환(`pardiso_ldlt` ↔ Kratos `AMGCLSolver`) 및 전환 메시지 출력**(§5.11.3), **`out_mode="drop"|"ghost"`** 정식 파이프라인 편입(§5.11.1)
  - (v2.1) 게이트 추가: V04를 FE 파트(Cu) + IGA 파트(W) tie로도 해석해 merge 대비 0.5% 이내, dependent instance 3개 배치 모델에서 instance별 결과가 단일 instance 해석과 1e-10 일치, 메쉬 리더/라이터(Set 키워드 포함) 왕복 무손실, 메모리 압박 시뮬레이션(`WHT_SOLVER_AVAIL_RAM_GB`)에서 `auto` 전환 메시지 및 `AMGCLSolver` 수렴 검증 통과, `HWT-IGA SOLVER` 로그 On/Off 동작 확인
  - 기간: 6주 → **8주**(어셈블리·메쉬 입출력 추가)
- **게이트 G3**
  - **V03** moment 케이스 5% 이내
  - **V04** bimetallic Cu/W: merge와 tie의 끝단 처짐 상대차 ≤ 1e-6. 변환단면 해 대비 오차는 `bimetallic-beam-tie-comparison.csv` 수준(p = 3 ≈ 0.54%) ± 0.2%p
  - **V05** multiplanar 계면 7단면(84 해석): monolithic 대비 차이가 Coreform 보고 수준(p = 2 ≤ 3.3e-6) 이하
  - **V02 국소 세분**: L3에서 peak 오차 ≤ 0.2%, DOF ≤ 80k (Coreform 0.105% @ 57,300 DOF, `roark-17p1-7a-refinement.csv`)
  - **E01** 박판 구멍(IPS): Roark Kt = 2.23568 기준 σ_max = 3.7261e-6 psi 대비 probe 오차 ≤ Coreform(−8.7%) 수준, 세분 시 개선 추세
  - **E02** 적응 세분 예제 재현

### Phase 4 — 고유치, 암시적·모드기반 동역학, 스펙트럼/랜덤, 열, 열변형 (6주)

- **산출물 (v2.1)**: Python consistent M(`kratos_backend/mass.py`), scipy 고유치, Python Bossak 선형 과도, Kratos Bossak + lumped 비선형 동역학, **모드 중첩 과도/조화 응답(`analyses/modal_dynamics.py`)**, add_1 G-3 스펙트럼 case-file 리더, 온도 의존 물성(S20 결과 경로)
- **게이트 G4**
  - **V09** 외팔보 충격·진동: f1 = 80.8 Hz 대비 1% 이내(Coreform 81.8), 랜덤 RMS 오차 ≤ 1.79% + 0.5%p, 스펙트럼 조합 규칙별 결과를 수기 계산과 1e-10까지 일치
  - **V10** 진동 보 (Kratos Bossak): f1 = 25,000 기준. 문서의 판정 기준 26개 중 물리량 기준(f1 오차, 위상, 응답 이력 RMS)을 전부 통과. HHT 고유 파라미터에 묶인 기준은 Bossak 대응값으로 바꾸고 차이를 리포트에 기록
  - **V11** NAFEMS T1: σyy(D) = 50 대비 4% 이내 (Coreform 2.99%)
  - **V12** NAFEMS T2 927 K, T3 36.60 °C @ 32 s, T4 18.3 °C: 문서 한계(1% / 2% / 3%) 이내. 2D 벤치마크는 얇은 3D 슬래브 + 대칭 penalty로 구성

### Phase 5 — 비선형, 잠김 방지, 접촉 (8주, C++ 승격 포함)

- **산출물**: NLGEOM(TL/UL), 초탄성, J2(T 의존 표), `WhtIgaApplication` MixedUP0/VMS 요소, 접촉(표면층 + 직접 `PairedCondition` 또는 C++/Python 활성집합 경로, S22)
  - (v2.1) 재료 래더(add_1 G-4 정정): ConstitutiveLaws(`CL.HyperElastic3DLaw`+`UpdatedLagrangianElement3D4N`, J2) → 다항 초탄성(Mooney-Rivlin, Yeoh, Ogden)은 `ConstitutiveModelsApplication` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 C++ 보완 필요(S19), SRI 옵션
- **게이트 G5**
  - **V07** Kirsch ν = 0.49999 혼합 (`v07.json`): **Hard Gate** = 원문 §6 `slope > 0.5` + 레벨별 오차 한계 `(e_u, e_σxx) ≤ (0.01, 0.2), (0.005, 0.1), (0.002, 0.05)` (+ 도메인 체적 일치 0.1%), **Target(목표치)** = 수렴률 `2.3 / 2.2` (Coreform 실측 `2.79 / 2.38` 기준)
  - **E04** Cook 막: u_y = 6.98946 mm 대비 1% 이내
  - **V08** Hertz: F_N = 481,579 lbf 대비 4% 이내(Coreform 1.84–3.74%), 미끄럼 F_T RMS ≤ 4.5%
  - **E05** door-lock: RF3 = −46.946 kN 대비 5% 이내
  - **E03** glass-to-metal seal: 13셀, 3 NLGEOM 열 step 수렴. 증분 수가 Coreform(각 step 15/15, cutback 없음)의 2배 이내

### Phase 6 — 성능과 확장 (지속)

- 프로파일 기반으로 Numba 핫스팟을 C++(pybind11)로 이전 (§2.2 기준).
- (v2.1) add_1 G-8 성능 비교 하네스(`perf_vs_coreform` 리포트 절), dependent instance 적분점 자료 재사용 최적화(§5.9.3), 1e5 요소 생성 52.5초(S16) 개선.
- C++ `ExtendedSplineSolidElement`로 구속식 제거.
- 반복 솔버 고도화: `IgaApp` `AdditiveSchwarzPreconditioner` (spike 필요; 기본 Kratos `AMGCLSolver` 및 메모리 기반 `auto` 전환은 Phase 3 조기 편입 완료).
- 전처리 병렬화(적분점 수 가중 분할). MPI는 Kratos MPI 빌드 여부를 확인한 뒤 판단.
- **Explicit 동역학 트랙 (Alpha, S24·S25 통과 기반)**: Kratos `MechanicalExplicitStrategy` + `ExplicitCentralDifferencesScheme` + Kratos lumped M (+ selective mass scaling), extension 구속식은 Kratos explicit 구속 처리 사용, **explicit penalty 접촉**(표면층). 검증: V-EXP, V-EXP-C.
- 선택: Nitsche/ghost penalty, Algoim 정확 형상, H3D.
- **게이트**: 1M DOF 3D 브래킷 end-to-end < 30분(16코어), 메모리 < 64 GB. 참고로 Coreform glass-to-metal seal은 16 CPU에서 mesh 318.8 s + solve 691.0 s.

---

## 10. 검증·시험 계획

### 10.1 Coreform 검증 매뉴얼 매핑

CSV 경로 접두어: `user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/`

| ID | Coreform 케이스 | 기준값 / CSV | sweep | 허용오차 | Phase |
|---|---|---|---|---|---|
| V01 | prescribed-displacement-cube (E=1, ν=0, 격자 간격 0.6) | 해석해 u = (x,2y,3z), σ = diag(1,2,3) | p 1–3 × C 0..p−1, 회전 Frame | 상대오차 ≤ 1e-5 (Coreform 최대 변위 1.5e-5, 응력 5.5e-5) | 1 |
| V02 | roark-17p1-7a (L=200, D=40, r=10, t=1, P=1e-6) | `roark-17p1-7a/` `-accepted-peak.csv`(4.3175e-6), `-qoi.csv`, `-refinement.csv` | p 1–3, h = 6/2^k (k=0–4), 국소 L1–L3 | peak(Roark) ≤ 1% (+ Abaqus 대비 보고), **응력 분포 L2 자기 수렴률 ≥ 1.8**(원문 §7.5), 국소 ≤ 0.2% | 1/3 |
| V03 | beam-bending (b=h=0.1, L=1, E=200 GPa) | `beam-bending/` `-accepted.csv`, `*-tip-deflection.csv`, `-traction-variants.csv` | 6개 하중, p 1–3, h = (0.1,0.1,0.25)/2^k | 5% (느린 케이스 10%) | 1/3 |
| V04 | bimetallic-beam-bending (Cu/W, b=10, L=100) | `bimetallic-beam-bending/` `*-accepted.csv`, `-tie-comparison.csv` | merge/tie × p × 하중 4종 | 5%, merge–tie ≤ 1e-6 | 3 |
| V05 | multiplanar-material-interfaces (7단면, 84 해석) | `multiplanar-material-interfaces/*.csv` | p 1–3 × h 4단계 | monolithic 대비 Coreform 수준 | 3 |
| V06 | pressurized-sphere (a=0.9, b=1, p_i=1e-6) | `pressurized-sphere/` `-qoi.csv`(허용오차 열) | p 1–3 × h = 0.1/2^k | CSV 허용오차 | 2 |
| V07 | infinite-plate-with-hole (ν=0.49999, T=10) | `v07.json` (원문 §6/§7) | p=2 h = 1/2, 1/4, 1/8 | Hard Gate: `slope > 0.5` + `(e_u, e_σxx) ≤ (0.01, 0.2), (0.005, 0.1), (0.002, 0.05)` (체적 0.1%); Target: `2.3 / 2.2` (Coreform `2.79 / 2.38`) | 5 |
| V08 | contact-indenter (Hertz, μ=0.3) | F_N 481,579 lbf, F_T 144,474 lbf, p0 4.599e6 psi | h = 0.5, 0.25, 0.125 | F_N ≤ 4% | 5 |
| V09 | cantilever-shock-vibration | f1 80.8 Hz, 랜덤 RMS u 340 µm, shock peak 965 µm | 모드 20개, ζ=0.02 | 1% / RMS +0.5%p | 4 |
| V10 | vibrating-beam (Coreform은 HHT, WHT_IGA는 Kratos Bossak; L=10, 2×2, E=1e6, ρ=1e-6) | f1 = 25,000, 26개 기준(물리량 기준) | Δt 1e-6–2e-6, α_B sweep | 문서 기준 | 4 |
| V11 | nafems-p18-t1 (열변형) | σyy(D) = 50.0 | p 1–3 | 4% | 4 |
| V12 | nafems-p18-t2/t3/t4 (열전달) | 927 K / 36.60 °C @ 32 s / 18.3 °C | h, Δt | 1% / 2% / 3% | 4 |
| V-COND | (자체) 조건수 (모체: **V02 Roark 판, 구멍 중심 이동**) | η sweep | η ∈ [1e-8, 1], h 4단계 | κ 평탄(≤ 10배), 기울기 −2 ± 0.3 | 2 |
| V-EXP | (자체) 1D 탄성파 전파 / 스텝 하중 (Explicit) | 1D 파동 전파 해석해 | p 1–3, lumped M, CFL/mass-scaling | 파면 도달 시각·피크 응력 2% 이내 | 6 |
| V-EXP-C | (자체) 두 블록 explicit 충돌 (penalty 접촉) | 운동량 보존, 반발 후 속도 해석해(탄성 충돌) | p 2, Δt sweep | 운동량 1e-6, 에너지 2% | 6 |

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
| R16 | 적분점 기하에서 SM consistent 질량 = 0 (확인됨, 트랙 A F1) | 중 | Python consistent M + Kratos 추출 K (R2 상시 병용), 비선형 동역학은 lumped | G4 |
| R17 | 구속 상수 누적 (확인됨, F2), 연쇄 구속(K11) | 중 | c = 0 규칙 + 코드 거부 + `compile_plan()` 구속 그래프 검사 + S21 | G1, G3 |
| R18 | 어셈블리 상호작용이 표면층·가상 노드에서 Kratos 접촉/mesh tying과 호환되지 않을 가능성 | 높음 | S22로 조기 확인, 실패 시 MPC tie·Python 활성집합 penalty·C++ 순 | S22 |
| R19 | 메쉬 포맷 방언(Nastran 필드 형식, LS-DYNA `*INCLUDE`, OptiStruct 확장, Radioss 버전차) | 중 | 왕복 시험 + 실제 파일 샘플 수집, 지원 외 카드는 issue로 보고 | G3 |
| R20 | ConstitutiveModelsApplication이 로컬 복원·수정분(166파일)이라 안정성 불확실 | 중 | S19 통과 재료만 공개, 실패 시 ConstitutiveLaws 재료만 | S19 |
| R21 | Gmsh GPL 의존 | 낮음 | 외부 프로세스 호출만, import 금지, 미설치 시 issue | 상시 |

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
| Gmsh | GPL-2.0+ | **외부 프로세스로만** — 일반 FE 파트 메쉬 생성 (D13). import 금지 |
| pyNastran | BSD-3 | bdf/fem 리더 (선택 의존 `mesh`) |
| meshio | MIT | msh 리더, VTU 쓰기 |
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
| D10 | Gmsh 선택 기능 포함 | **확정 (2026-11-14, D13으로 대체)**: 일반 FE 파트 메쉬 생성에 Gmsh 외부 프로세스 사용. immersed tessellation은 계속 OCP BRepMesh | – |
| D11 | 메쉬 리더 구현 | **확정**: bdf/fem = pyNastran, inp/rad/k = 자체 파서, msh = meshio | – |
| D12 | 메쉬 내보내기 범위 | **확정**: inp/fem/bdf/rad/k로 메쉬 + 노드 결과 + 요소 평균 응력 | – |
| D13 | 일반 FE 메쉬 생성 | **확정**: Gmsh 외부 프로세스 생성 지원 + 외부 파일 입력 | – |
| D14 | IGA–FE 혼합 시기 | **확정**: Phase 3 Beta | – |
| D15 | 파트 자동 병합 | 제안만 / 자동 병합 | 제안만 (확정은 사용자) |
| D16 | IGA–IGA tie 기본 방식 | penalty 스프링 / mortar tying(S22 통과 시) | penalty 스프링 (Coreform과 같은 penalty 계열) |
| D17 | 결과 저장·쿼리 아키텍처 | **확정 (2026-11-14)**: `SQLite`(메타·요약·Probe·인덱스) + `VTKHDF`(대용량 필드·시계열) 하이브리드 + `ResultView` 통합 쿼리 | – |
| D18 | Python API 아키텍처 | **확정 (2026-11-14)**: Pydantic v2 불변 Spec + Hexagonal(Ports & Adapters) + CQRS(`ResultView`) 독자 설계, Coreform 호환은 `wht_iga.compat.coreform` 어댑터로 분리 | – |

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

## 부록 C. Track H 후속 통합 기록 (2026-11-19, 상용급 목표)

| 항목 | 상태 | 근거 |
|---|---|---|
| 5종 기법 교차 검증 | 완료 (Kratos 15/15 + 비교 매트릭스 2/2) | `test_discretization_methods_kratos.py`, `test_method_comparison.py`, S34/S35 |
| `auto` → `mortar_l2` (Kratos 네이티브 MPC) | 완료 (TDD Red→Green) | methods.py 레벨, 기본 beta 그대로 exact |
| `recommend_discretization` + OBB 프레임 | 완료 (단위·Kratos 게이트) | S34, `test_obb_frame.py` 2/2 |
| THB 본선 연결 | 완료 (`test_thb_pipeline.py` 5/5: 단축·굽힘·구체·거부) | `JobSpec.thb_refinement` + pipeline 분기; 컴파일러 가드·User API 무수정 |
| API 분할 (options/solve/recommend + shim) | 완료 (91/91 green) | S38; backward-compat 유지 |
| numba 가속 | 결론: 변경 없음 (5% 컷 미달) | 웜업 Python 0.1 s / wall 0.7 s; leggauss 호이스팅은 integrator 요청으로 이관 |
| 검증 리포트 생성 (HTML+PNG) | 완료 (게이트표+비교+듀얼+증거셋) | `verify_report.py`, 산출물 `docs/verification/report/` |
| 실증 자료 세트 (vtkhdf/vtu/msh/fem/stl/png) | 완료 (roundtrip 검증) | `post/evidence.py`, bulk bar 6종; hvascii 미지원 명시 |
| 프리셋 (contact_iga/bending_thb/bulk_linear) | 완료 (단위 4/4) | `analyses/presets.py` |
| Locking C++ (MixedUP0/VMS) | 미결 (차주 결정) | 미결 시 E04/V07 게이트 불가 |
| 상용화 (CI·견고성·G-8·사용자 가이드) | 계획됨 (Phase 6) | 본 문서 §9 게이트와 별도 트래킹 |
