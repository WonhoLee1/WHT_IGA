# Kratos 연결성 사전감사 (연결 기술이슈 검증)
작성: 2026-11-14 / 근거: Kratos master 소스 직접 확인 (URL 포함)
결론 먼저: **C̃^e/Algoim 주입-원안(Zero-C++ 유지)은 불가**, 우회로 3종이 있음.

## 1. 확인된 팩트 (소스 근거)

- IgaApplication 존재 + 요소 등록 (`iga_application.cpp`):
  `SolidElement`(2D 한정 — `KRATOS_ERROR_IF(dim != 2)`), `LaplacianElement`,
  `IgaMembraneElement`, Shell 3p/5p/6p, Stokes/NavierStokes, Truss/Beam.
- 조건: `SupportNitscheCondition`, `CouplingNitscheCondition`,
  `SbmLaplacian/SbmFluid/SbmSolid` Dirichlet/Neumann, Penalty/Lagrange 계열.
- 변수: `NITSCHE_STABILIZATION_FACTOR`, `EIGENVALUE_NITSCHE_STABILIZATION_*`
  (고유값 기반 자동 γ — master-plan §5.6 대응책이 이미 구현돼 있음),
  `INTEGRATION_POINTS` / `INTEGRATION_WEIGHTS` (SBM 조건용 커스텀 적분점),
  `KNOT_VECTOR_U/V/W`, `BREP_ID`, `PARAMETER_SPACE_CORNERS`.
- 모델러: `NurbsGeometryModeler`(+Sbm/GapSbm), `IgaModeler(+Sbm)`,
  `RefinementModeler`, `ImportNurbsSbmModeler`.
- **치명적 1**: `SolidElement::GetIntegrationMethod()` = 하드코딩 `GI_GAUSS_1`
  (단일 적분점). `CalculateLeftHandSide`는 `r_DN_De[0]`·단일 `INTEGRATION_WEIGHT`
  만 사용 (`custom_elements/solid_element.cpp`).
- **치명적 2**: Kratos IGA의 잠김(immersed) 경로는 FCM/CutFEM이 아니라
  **SBM(Shifted Boundary Method)** — 배경셀을 자르지 않고 대리경계에서 보정
  (`iga_modeler_sbm.cpp`, `nurbs_geometry_modeler_sbm.h`의 RegularGrid 생성).
- 긍정: `SolidElement`은 표준 `CONSTITUTIVE_LAW` 사용 → 재료 라이브러리 재활용 가능.
- 포장: `applications/IgaApplication/pyproject.toml` 존재 → pip wheel 포함 추정
  (설치 후 `import KratosMultiphysics.IgaApplication`으로 확정 필요).

## 2. 연결점별 판정

| # | 연결 | 판정 | 근거/대안 |
|---|---|---|---|
| L1 | splinepy C^e → Kratos 주입 | 불가(주입 API 없음) → **우회: 격자 복제** | Kratos가 knot/order/span JSON으로 자체 생성. 동일 격자를 `NurbsGeometryModeler` 파라미터로 재현하면 C^e 일치. 모델러 JSON 스키마 확인이 다음 단계 |
| L2 | Algoim 절단QP → 요소 적분 | **불가** (Zero-C++ 하) | 요소가 단일 Gauss점 하드코딩. 대안: (a) Python+numba 자체조립(PoC 01–04가 수학 검증済), (b) SBM 경로(절단적분 불필요), (c) C++ 커스텀 요소(Zero-C++ 파기) |
| L3 | AgFEM T행렬(C̃^e) → Kratos | **가능 예상 (MPC 경로)** | Kratos 범용 `MasterSlaveConstraint`로 이산확장 구속을 Python에서 부과. 실설치 후 PoC 필요 |
| L4 | bezier_grid.mdpa → Kratos | **부적합** | IGA는 mdpa Bézier 블록이 아니라 모델러 JSON으로 NURBS 지오메트리 생성. mdpa는 节点·하중용으로만 쓰고 격자 정의는 JSON으로 전환 |
| L5 | CSR 추출 → PyPardiso | 가능 예상 | `BuilderAndSolver.GetSystemMatrix()` + `ScipyConversionTools.ToCsr` (설치 후 API명 확정) |
| L6 | Nitsche 약형식 BC | **가능 (API 확인)** | `SupportNitscheCondition` + 고유값 자동 γ. PoC-04의 γ≳1/η 문제를 Kratos가 자체 해결 |
| L7 | THB 세분화 | 가능 예상 | `RefinementModeler` (스키마 확인 필요) |
| L8 | SBM 네이티브 경로 | 가능 (별도 트랙 권장) | 배경격자+SBM조건만으로 Zero-C++ 완결. FCM 트랙과 정확도·구현비 비교 대상 |

## 3. 권장 아키텍처 수정 (master-plan 반영용)

1. Stage 2-A 출력물을 `C̃^e 주입`에서 **`MasterSlaveConstraint 목록 + Nitsche 조건 + JSON 격자정의`** 로 변경 (L1/L3/L4).
2. 조립은 **Python+numba 자체조립을 1순위**로 격상, Kratos 조립은 SBM 트랙 전용으로 분리 (L2).
3. SBM 트랙을 병행 PoC로 추가 — Kratos 네이티브 테스트(`IgaApplication/tests`)에 SBM 케이스가 있어 회귀 기준으로 활용 가능.
4. Kratos 설치 후 즉시 확인할 3종: `IgaApplication` import, `GetSystemMatrix`/`ToCsr` 명칭, `NurbsGeometryModeler` JSON 스키마.

## 4. 출처
- 요소 목록/등록명: `applications/IgaApplication/iga_application.cpp` (master)
- 단일적분점: `.../custom_elements/solid_element.cpp` (`GetIntegrationMethod`, `r_DN_De[0]`)
- SBM 모델러: `.../custom_modelers/nurbs_geometry_modeler_sbm.h`, `iga_modeler_sbm.cpp`
- 변수: `.../iga_application.cpp` (`KRATOS_REGISTER_VARIABLE` 블록)
- 본 보고서의 PoC 수학 근거: `poc_codes/01–08` (전수 PASS)
