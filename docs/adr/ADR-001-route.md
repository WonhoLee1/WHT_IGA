# ADR-001: 해석 경로 선택 (R1/R2/R3)

- 상태: **확정** (사용자 승인 2026-11-14 — v2 D1)
- 날짜: 2026-11-14 / 브랜치: WHTIGA / Kratos 10.4.3 소스빌드 (`E:\kratos-build\install`, 96bad120)
- 결정 트리 근거: `research-coreform-IGA/master-plan_v2.md` §5.3
  ("S03·S04·S05 통과 → R1. S05만 실패 → R2. S04 실패 → R3 + C++ 조기 착수")

## 결정

**R1 = (a) + (c)** — 셀 단위 적분점 기하 + StructuralMechanics 요소 + 구속식 T.
R2/R3 전환 사유 없음. (b)는 바인딩 부재로 불가 유지 → (e)로 이관 (v2 원안대로).
(e) C++ `WhtIgaApplication`은 Phase 5 승격 기준 충족 시에만 (조기 착수 안 함).
SBM은 비교 트랙 유지 (§10.5, Roark·구면 베이스라인은 Phase 1).

## 근거 (실행 수치, 허용 기준 대비)

| 조건 | 결과 | 핵심 수치 |
|---|---|---|
| S03 셀 QP 기하 | 예 | 셀당 기하 1개·27점, 가중치합=셀체적, SFV행합 1, 지지 27노드 |
| S04 SM 요소 | 예 | K 상대오차 p=1,2 8.7e-17/1.1e-16/2.3e-16 + p=3 3.1e-16/4.2e-16 (≤1e-12, S04 파일 parametrize로 확인), 패치 3.1e-15 (≤1e-10), TL/UL 초기화 OK |
| S05 구속식 | 예 | TᵀKT 대비 4.2e-17, 1e5 생성 2.7s (<30s), master DofSet 포함 |
| S06–S18 | R1 변경 없음 | 부분/아니오는 아래 failover로 흡수 (상세 `docs/spikes/TRACK_*.md`) |

검증: `WHT_IGA/tests/kratos_spikes` 전수 61/61 통과 (통합 담당 재실행 확인,
A 27 + B 34). isolation 2건은 unit 런(`-m "not kratos"`)에 포함 — 결과는 G0 보고서 참조.

## 채택하는 failover와 v2 정정 제안 (v2 본문은 수정하지 않음, X8)

1. 표면 QP 기하는 **적분점당 1개씩** 생성 (S07) — §5.4 "표면 적분점 기하" 가정 수정.
2. penalty 표면 규칙은 **면-셀당 적분점 수 ≥ (p+1)²** (S09 정정: p=3 규칙 4는 12점 < 16점으로 1.8e-4 고착; 규칙 6 = 면-셀당 24점으로 1/β 회복. 차수 ≥2p는 정확도 권장일 뿐 고착 원인이 아님). 규칙: `rule_for_degree` p≤1→2, p=2→4, p=3→6. ADR-002의 p=3 sweep 무고착 결과와 모순 없어 보임: sweep 격자(2×2×2)는 면-셀 겹침으로 rank가 보충되나, 보장 조건은 면-셀당 개수이므로 기본값을 6으로 올림.
3. `SupportLaplacianCondition` 삭제 → 열 Dirichlet는 강한 부과; 과도 열은 참조 조립기 (S11).
4. E(T)는 `ThermalElasticIsotropic3D` + TableAccessor(`ReadMaterialsUtility` JSON "accessors")로
   동작 확인 (S20: clamped bar 해석해 대비 6.0e-9, p=4 n_x=24, 본 작업 실행). 온도 구간마다
   Properties 직접 갱신 failover는 폐기. alpha(T) accessor는 해당 법칙이 우회(const-α와 1e-15 일치) — Phase 4 과제.
5. Dirichlet 상수 처리: 고정 master는 관계식에 유지 + c=0 (Kratos가 자체 접음). 일반 규칙은
   **"c는 항상 0, 오프셋은 고정된 master 자유도로 표현"** (`add_linear_tie`가 0이 아닌 constant를
   받으면 오류). Newton 전략에서는 구속식 결과 검증 시험 추가 (REVIEW F2).
6. K 비교는 무구속 상태에서 (BlockBuilder가 고정행을 1로 치환, S04).
7. Volume 객체 수명: 파이썬 참조를 해석 끝까지 보관 (RHS 0-수렴 방지). `geometry.make_nurbs_volume`에 보관 로직 추가 제안 (공유 모듈 변경 제안 #1).
8. Thermal settings 기본값 (TEMPERATURE_GRADIENT, FACE_HEAT_FLUX) — `thermal.py` 포함됨.
9. 스프링 사용 노드의 미사용 ROTATION은 Fix (S10/S14).
10. 열 Dirichlet penalty·3D Nitsche는 C++ 트랙으로 (S11/S18). Phase 1은 penalty 기본 (S09).
11. QP consistent 질량행렬 = 0 → 원인 확정: `UseGeometryIntegrationMethod()=true` 고정 때문에
    높은 차수 Gauss 규칙을 타는데 QP 기하가 기본 규칙에만 점을 가져서 점 목록이 빔
    (TRACK_A_REVIEW F1, `base_solid_element.cpp` l.629). C++ 없이 해결:
    Kratos 추출 K + Python(Numba) consistent M + scipy eigsh.
    lumped 대 consistent 실측: 1차 모드 0.64%→0.02% 수렴, 상위 모드 최대 7.9%→2.3%.
    Phase 4 운용안 (V09 eigsh / V10 Python Bossak 루프 / 비선형 Bossak+lumped).
12. MixedVolumetric QP 크래시 → V07/E04 잠김 방지는 C++ MixedUP0 Phase 5 1순위 (S13).
13. 트랙 C(S19–S24)는 아래 "트랙 C" 절 참조 — R1 변경 없음 (S22 접촉 프로세스·S19 일부 법칙은 C++ 트랙으로).

## 트랙 C (S19–S24)

시험 파일(`tests/kratos_spikes/test_s*.py`)과 본 작업 실행 수치 기준. 상세는 `docs/spikes/SNN.md`.

| spike | 내용 | 상태 |
|---|---|---|
| S19 구성법칙 | `HyperElastic3DLaw`는 TL/UL 요소에 부착·풀이 OK. `ConstitutiveModelsApplication.LargeStrain3DLaw`은 QP 기하에서 실패 (`QuadraturePointGeometry`가 `Volume()` 미구현) → C++ 트랙 | 시험 파일 기록 기준 (본 작업 미실행) |
| S20 TableAccessor | `TableAccessor`는 Python 직접 노출 없음 → `ReadMaterialsUtility` JSON "accessors"로 부착. `ThermalElasticIsotropic3D`가 E(T)를 따름 (clamped bar 해석해 대비 rel_err 6.0e-9 <1e-8, p=4 n_x=24, 본 작업 실행). alpha(T) accessor는 우회됨 | E(T) 확인 (본 작업 실행) |
| S21 Newton+구속식 | `LinearMasterSlaveConstraint`(c=0) + `ResidualBasedNewtonRaphsonStrategy` + TL, Python TL 대조 (목표 1e-8). S05 c=0 규칙과 일관 | 본 작업 미실행, 시험·문서 참조 |
| S22 표면 접촉 | mortar/`penalty_contact` 프로세스는 IGA QP 요소 포함 시 `GenerateFaces()` 미구현으로 초기화 실패 → Phase 5 | 시험 파일 기록 기준 (본 작업 미실행) |
| S23 비적합 매퍼 | `nearest_element` 매퍼 가중치는 분할단위성 0.0·선형 재현 1.1e-16, constant/kinematic patch 1e-12. 결합 patch는 비변분 간극 0.085987 (대역 단언). `AssignMasterSlaveConstraintsToNeighboursUtility`는 48 구속·선형 오차 7.942e-02 (대역 단언, 본 작업 실행) → 매퍼 경로 채택 | 예 (경로 결정) |
| S24 explicit | Kratos explicit + HRZ lumped mass 조사. REVIEW #7: `Initialize()` 후 NODAL_MASS를 HRZ로 덮어쓰는 혼합안은 미시험으로 명시 | 시험 파일 기록 기준 (본 작업 미실행) |

트랙 C는 R1을 바꾸지 않는다 (S19 일부·S22는 C++ 트랙으로 이관, S20 E(T)는 Phase 1 penalty 경로에서 사용 가능).

## 결과

- Phase 1은 R1로 간다 (선형 정적, trimmed, penalty, extension 없음).
- SBM 비교 트랙 범위 변경 없음.
- G0 게이트의 p=3 교차비교 잔여 항목은 `docs/G0_report.md` 참조 (조건부 항목).
