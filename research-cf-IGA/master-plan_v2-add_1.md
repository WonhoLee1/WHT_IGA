# master-plan v2 Addendum 1 — Coreform IGA 기능 ↔ WHT_IGA 대응 추적표 (v1.1)

- 문서 버전: **add_1 v1.1** (2026-11-14). v1.0 원본은 `master-plan_v2-add_1_v1.0.md`에 보존.
- 역할: **추적표(감사용).** Coreform 기능마다 WHT_IGA의 담당 절, 판정, 검증 케이스를 기록한다. 계획 내용(설계·일정·결정)은 **`master-plan_v2.1.md`에만** 둔다. 이 문서와 v2.1이 다르면 v2.1이 맞다.
- 상태: v1.0의 갭·제안은 모두 v2.1에 편입되었다(§4). 판정 오류 9건을 정정했다(§6).
- 근거:
  1. Coreform 문서: `user-guide_combined.md` / `.html` (md 82개 전부 결합, `generated-downloads/` CSV는 별도), 원본 `user-guide/`
  2. 계획: `master-plan_v2.1.md` (v2 + add_1 + add_2 + 트랙 A 통합), `master-plan_v2-add_2.md` (add_1 검토·워크플로·메쉬 설계)
  3. 실측: Phase 0 트랙 A (`WHT_IGA/docs/spikes/TRACK_A.md`, `TRACK_A_REVIEW.md`, S02–S06, S13, S15, S16)
  4. Kratos 소스: `D:\PythonCodeStudy\kratos-src` (ConstitutiveLaws, **ConstitutiveModels**, ContactStructuralMechanics, Mapping, StructuralMechanics)

---

## 0. 대조 방법

1. 결합 문서 82개 섹션을 기능군으로 분류: 설치/구성(6) · 이론(2) · 사용자 가이드 핵심(7) · Python API(11) · 릴리스 노트(5) · 응답 스펙트럼(4) · 스크립팅 가이드(5) · 검증 매뉴얼(12 문제 + 인덱스) · 예제(5) · 키워드 지원 가이드(49 케이스) · CIFA 역량(JS 렌더, 내용 없음).
2. 각 Coreform 기능을 `master-plan_v2.1.md`와 1:1 대조한다. 판정은 다음 5단계:
   - `동등`: 같은 기능을 같은 성숙도 목표로 계획
   - `동등(방법 대체)`: 같은 결과를 다른 방법으로 (예: HHT → Bossak)
   - `부분`: 일부만 계획
   - `미포함(의도적)`: 범위 밖으로 결정
   - `갭`: 계획 누락 (v1.1 시점에는 모두 v2.1에 편입되어 남은 갭 없음)
3. 각 행에 **v2.1 반영 위치**와 **검증 케이스**를 적는다. 실측으로 판정이 바뀐 행은 근거(spike 번호)를 적는다.

---

## 1. 커버리지 추적표

### 1.1 격자·세분화·소형 셀

| Coreform 기능 (출처) | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| 직교 배경 격자, p=1–3, C0–C^{p−1}, element size/interval count/축별, edge/node_centered, padding, Frame(회전 포함) (`python-api/coreform/mesh`) | §1.3, §7.2, §5.2.1(항등 사상, S02 실측 detJ=1) | 동등 | V01–V03 |
| tessellation_size, 삼각형 분할 절단 | §7.1, §7.4 | 동등 | V02, V06 |
| 계층 세분화 + 2:1 균형 + CadEntities/Sphere/Box region + offset p·h/2^L | §7.2 THB + §5.2.3 구속식 | 동등 | V02-local, E02 |
| small_cell_volume_ratio 기본 0.2, [0, 0.5) | §7.5 (기저 지지 기준 extension, MPC) | 동등 | V-COND, V06 |
| U-spline | §1.4 | 미포함(의도적, 특허) | – |
| **Body-fit: 기존 정렬(mapped) all-hex 메쉬 → 스플라인 변환** (2026.3, sweep·bottom-up 불가) | §1.3, §5.10.4 | 동등 (v1.0 판정 정정: A4) | E04, V10 |
| Selective integration (nearly-incompressible, 2026.3 §3.1) | §1.3 잠김 방지 행, Phase 5 | 동등 (v1.0 판정 정정: A1) | V07, E04 |
| 메쉬 미리보기 (해석 전 이산화 확인, 2026.3 §2) | §1.3, Phase 1 | 동등 (신규) | 단위시험 |

### 1.2 메쉬 종류·파트·어셈블리 (신규 절)

| Coreform 기능 (출처) | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| IGA 메쉬와 일반 Abaqus 메쉬 공존 (`workflow-overview` Mesh settings) | §5.10 (`NativeMesh` 파트, 메쉬 파일 입력 또는 Gmsh 외부 생성) | 동등 (신규, 원래 갭: G-11) | 자체 |
| IGA 파트 + 일반 FE 파트 혼합, matching/nonmatching 계면 tie (2026.3 §3.3, 2026.9 경계 세트 분리) | §5.9.4, Phase 3 Beta (D14) | 동등 (v1.0 판정 정정: A5) | V04 FE판 |
| Instance-aware 메쉬: part 또는 independent instance 소유, dependent instance는 part 메쉬 재사용 (`scripting-guide/meshes.md`, 2026.7 §2.1, dependent는 Production) | §5.9.1, §5.9.3 | 동등 (신규, 원래 갭: G-10) | dependent 3개 배치 시험 |
| 다중 IGA 파트 어셈블리 (2026.3 §3.3) | §5.9 (루트 ModelPart + SubModelPart 계층) | 동등 (신규, G-9) | V04, E05 |
| Multi-cell(불리언 병합) 단일 스플라인 장, 셀별 재료 (2026.7 Production) | §7.7, §5.9.4 | 동등 | V04, V05, E03 |
| STL·격자 구조 입력 (Beta) | §5.10.3 (메쉬 외곽 면 → immersed 형상), STL 입력 | 동등 | V02 메쉬 입력판 |
| (WHT_IGA 추가 요구) 메쉬 파일 inp/fem/bdf/rad/k/msh 읽기, inp/fem/bdf/rad/k 쓰기 | §5.10.2, §5.10.6 (D11, D12) | Coreform에 없음 (우위 S-11) | 왕복 시험 |

### 1.3 경계·하중·결합

| Coreform 기능 (출처) | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| penalty 필수 BC (`bc_penalty_value`, 검증 1e10) + 변위/핀/대칭/사용자 함수 변위 | §5.4 BC-A/B/C, §7.6 | 동등 | V01, V07 |
| 압력·총하중 압력·트랙션·중력/체적력·RP 모멘트 (키워드 가이드 *DLOAD 10건·*DSLOAD 7건) | §5.4 | 동등 | V03 |
| amplitude, 다단계 step, 하중 변경/리셋 | §5.6, Phase 3 | 동등 | V10 |
| surface tie | §5.9.4 (항등 MPC / 가중 MPC / 합성 MPC / penalty / mortar `mesh_tying_process`) | 동등 | V04 tie |
| kinematic / distributing coupling | §5.9.4 (SM `ImposeRigidMovementProcess` 또는 MPC / RBE3형 가중 MPC) | 동등 | V03 moment |
| follower / resultant 트랙션 | §5.4 `PRESSURE_FOLLOWER_LOAD` | 동등 | V03 변형 |
| film, radiation (경계 삼각망 위 적분) | §7.11 (`ThermalFace3D3N`) | 동등 | V12 |
| surface-to-surface 접촉, penalty Coulomb μ, ALL EXTERIOR (자기접촉 제외) | §5.9.4 (ContactSM `penalty_contact_process`/`alm_contact_process` + 표면층), Phase 5 | 동등 | V08, E05 (S14, S22) |

### 1.4 해석 종류

| Coreform 기능 (출처) | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| 선형 정적·다단계, NLGEOM | §5.6, Phase 1/5 | 동등 | V01–V06, E03, E04 |
| 고유진동수 Lanczos·AMS, 질량 정규화 | §5.6, §7.9: Kratos 추출 K + **Python consistent M** + scipy `eigsh` | 동등(방법 대체) (v1.0 판정 정정: A3; 근거 S06, F1) | V09 |
| 암시적 동역학 HHT | §7.9: 선형 = Python Bossak, 비선형 = Kratos Bossak + lumped M | 동등(방법 대체, D5) | V10 |
| 응답 스펙트럼 ABS/SRSS/NRL/TENP/CQC/GRP/DSC + ALGEBRAIC/SRSS/R40/R30 | §7.10 | 동등 | V09 |
| 랜덤 응답 PSD BASE/FORCE/DB + Segalman RMS von Mises | §7.10 | 동등 | V09 |
| (WHT_IGA 확장) 모드 중첩 과도 / 조화 응답 | §7.10, Phase 4 | Coreform에 없음 (우위 S-13) | V10 modal / 해석해 |
| 정상/과도 열전달, 열변형 | §7.11, Phase 4 | 동등 | V11, V12 |
| Explicit 동역학, Abaqus/Aqua | §1.3, §5.6, §7.9, Phase 6 | Explicit은 **Coreform 미지원(WHT_IGA 추가 목표, 우위 S-13)** / Aqua는 미포함(의도적) | V-EXP (S24) |
| 좌굴 고유치 | – | 해당 없음 (양쪽 모두 없음) | – |
| 에너지 이력 (2026.11 Alpha) | §4.3 `energy` 표 | 동등 | – |

### 1.5 재료

| Coreform 기능 | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| 선형 탄성 | §1.3 | 동등 | V01 |
| 온도 의존 선형 탄성·열팽창 | §1.3, §7.11: $E(T)$는 `TableAccessor`+`CL.ThermalElasticIsotropic3D`(`6.013e-09`, S20), $\alpha(T)$·`LinearElastic3DLaw`는 셀별 `Properties` 복제 채택 | 부분 → 대안 동등 (S20 실측; `SetTable` 및 `α(T)` `TableAccessor` 미반영) | V11, E03 |
| Neo-Hookean, J2 다선형 소성(온도 의존) | §1.3, Phase 5 (`CL.HyperElastic3DLaw`+`UpdatedLagrangianElement3D4N`, S19) | 동등 | E03, E04, E05 |
| 그 밖의 초탄성 (Mooney-Rivlin, Yeoh, Ogden) | §1.3: `ConstitutiveModelsApplication`은 `MooneyRivlin`/`Yeoh` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 실패(S19) → C++ 보완(`QuadraturePointGeometry::Volume()` 또는 `WhtIgaApplication`) 필요 | 부분 (S19 실측: C++ 보완 후 동등) | patch + Cook 막 |
| 잠김 방지 혼합 요소 (U/P0, VMS, B-bar/F-bar) | §1.3, §5.7: C++ `WhtIgaApplication` (S13: Kratos 혼합 요소 크래시 확인) | 동등 (C++ 경로) | V07, E04 |

### 1.6 워크플로·출력·API·검증 방식

| Coreform 기능 | v2.1 위치 | 판정 | 검증 |
|---|---|---|---|
| 5단계 워크플로 (모델 설정 → trimmed mesh → interop → Abaqus 해석 → 후처리) | §3.4 W0–W5 (`job.resolve()`가 interop 역할) | 동등(방법 대체) (신규, G-14) | E01 end-to-end |
| 키워드 번역 시 mixed native+IGA 대상 분해, 불가 시 거부+진단 | §3.4 W3, §5.9.7 | 동등 (신규) | 단위시험 |
| U, S, PEEQ, E, NE, LE, NT11, PRESS_IGA | §1.3 출력 행 | 동등 | 전체 |
| point/line/multi-point probe + probe JSON | §5.4, §7.12 | 동등 | V02, E01 |
| conforming 시각화 proxy 격자 + ParaView 읽기 | §7.12 | 동등 | – |
| unified ODB (IGA + native 결과 한 파일, 2026.11) | §3.4 W5, §5.9.6 (VTKHDF instance 블록) | 동등 (신규) | – |
| results/quadrature SQLite | §4.2–4.4 | 동등 | V06 |
| portable spec + validate + journaling + Job Monitor·백그라운드 메쉬 | `Job.to_json/validate`, G-2 진행 콜백 | 동등 | – |
| 검증 매뉴얼 12 + 예제 5 + CSV 기준값 + 문서 자동 생성 | §10, §10.6 | 동등 이상 (S-7) | – |
| 일관 단위계 | G-6 단위계 가드 (Phase 1), 메쉬 파일 `unit_system` 필수(§5.10.2) | 동등 | 단위시험 |

---

## 2. 갭 처리 결과 (G-1–G-14)

| ID | 내용 | 처리 | v2.1 위치 | 수용 기준 |
|---|---|---|---|---|
| G-1 | `.inp` 읽기 서브셋 | **흡수·축소**: 메쉬·집합 읽기는 §5.10 메쉬 리더로. 해석 덱 키워드 번역은 범위 밖(Coreform 검증 입력 파일이 크롤에 없어 재사용 대상도 없음) | §5.10.2, §1.4 | 리더 왕복 시험 |
| G-2 | Job 진행 콜백·취소 | 채택 | Phase 3 | E03 규모 진행률 단조 증가, cancel 시 부분 산출물 정리 |
| G-3 | 스펙트럼 case-file 리더 | 채택 | Phase 4 | V09 값을 case 파일 입력으로 재현 |
| G-4 | 재료 래더 | **S19 실측 재정**: Neo-Hookean은 `CL.HyperElastic3DLaw`+`UpdatedLagrangianElement3D4N` 사용. 다항 초탄성(Mooney-Rivlin/Yeoh/Ogden)은 `ConstitutiveModels` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 C++ 보완(`QuadraturePointGeometry::Volume()` 또는 `WhtIgaApplication`) 필요 | Phase 5, S19 | 재료별 patch 1건 + Cook 막 1% (Neo-Hookean 외 1종 이상) |
| G-5 | THIN_FEATURE 자동 세분 루프 | 채택 | Phase 3 | 수동 대비 동일 정확도, 개입 0회 |
| G-6 | 일관 단위계 가드 | 채택 | Phase 1 | 단위 불일치 모델 3종 차단 |
| G-7 | 출력 이름 매핑표 | 채택 | Phase 3 | V02/V06 CSV 직접 비교 |
| G-8 | 성능 비교 하네스 | 채택 | Phase 6 | E03 1회 측정 + Phase 6 게이트 |
| G-9 | 멀티파트·어셈블리 (Kratos SubModelPart, MPC, mesh tying, 접촉, mapper, rigid movement 연계) | 채택 (add_2 + 사용자 요청) | §5.9, Phase 3 | V04 FE판, dependent instance 시험, S22/S23 |
| G-10 | Instance(dependent/independent) | 채택 | §5.9.1, §5.9.3, Phase 1(모델)/3 | instance별 결과 = 단일 해석 1e-10 |
| G-11 | 메쉬 입출력 6종 + 일반 FE 파트 + Gmsh 생성 | 채택 (D11–D13) | §5.10, Phase 2–3 | 왕복 무손실 |
| G-12 | Body-fit 변환 | 채택 | §5.10.4, Phase 3 | E04, V10 |
| G-13 | IGA–FE 혼합 | 채택, Phase 3로 앞당김 (D14) | §5.9.4 | V04 FE판 merge 대비 0.5% |
| G-14 | 워크플로 W0–W5, `job.resolve()`, 메쉬 미리보기 | 채택 | §3.4, Phase 1 | E01 end-to-end, 진단 리포트 |

---

## 3. 설계상 우위 항목 (S-1–S-12)

비공개 정보에 대한 비교는 "공개 문서에 없음"으로 표현한다(v1.0 정정: A9).

- **S-1. 솔버 라이선스 독립**: Coreform은 Abaqus/Standard + HF5 + CAE 플러그인이 필요하다. WHT_IGA는 Kratos 소스빌드 + Python만으로 같은 물리를 재현한다.
- **S-2. Oracle 교차검증**: Numba 참조 조립기와 Kratos 경로 비교(트랙 A S04 실측: K 상대오차 1e-16 수준). 이런 독립 교차검증 경로는 Coreform **공개 문서에 없다**.
- **S-3. 조건수 게이트**: V-COND (η ∈ [1e-8, 1]에서 κ 평탄, h 기울기 −2 ± 0.3). Coreform 공개 기준에 없는 정량 안정성 게이트.
- **S-4. Dirty-CAD 내성**: 3축 투표 + winding fallback + octree fallback + `QUAD_FALLBACK` issue.
- **S-5. 개방형 산출물**: VTKHDF/VTU + SQLite + probe JSON. Coreform은 ODB + 전용 ParaView 플러그인.
- **S-6. FCM α 기준선 + SBM 비교 트랙**.
- **S-7. 문서-증거 동기**: 케이스별 자동 리포트 + CSV 기준값 + oracle 차이.
- **S-8. 시간적분 선택권**: Bossak 기본 + Newmark/BDF2. (Coreform 검증은 Abaqus HHT 계열 옵션 — TRANSIENT_FIDELITY, MODERATE_DISSIPATION — 을 사용.) V10에서 물리량 기준 동등성으로 입증.
- **S-9. Nitsche 선택권 + β sweep 확정 절차**.
- **S-10. C++ 승격의 측정 기반 규칙**: 50% 조립 오버헤드·30% 핫스팟 같은 수치 기준. 트랙 A에서 질량 행렬 문제를 C++ 없이 해결한 것이 첫 적용 사례.
- **S-11. 메쉬 포맷 범위 (신규)**: inp/fem/bdf/rad/k/msh 읽기 + 5종 쓰기. Coreform은 Abaqus 메쉬만 다룬다.
- **S-12. Kratos 기능 재사용 (신규)**: 접촉(penalty/ALM/MPC), mortar mesh tying, mapper, rigid movement, link constraint를 그대로 쓴다. Coreform이 Abaqus 기능을 쓰는 것과 같은 구조를 오픈소스로 구성.
- **S-13. 모드 중첩 동역학 및 Explicit 동역학 확장 (신규)**: Phase 4에서 모드 중첩 과도/조화 응답을 지원하고, Phase 6에서 Coreform이 미지원하는 Explicit 동역학(`MechanicalExplicitStrategy` + `ExplicitCentralDifferencesScheme`, S24/V-EXP)까지 확장한다.

---

## 4. v1.0 수정 제안의 v2.1 편입 결과

| v1.0 제안 | 결과 | v2.1 위치 |
|---|---|---|
| §5.7 C++ 승격 기준에 "G-4 재료" 추가 | **채택 안 함**: ConstitutiveModels에 해당 재료가 있음(A7). 대신 S19 spike 추가 | §5.8 S19 |
| §10.6 리포트에 `perf_vs_coreform`·`unit_check`·`mapping_alias` 절 | 채택 (G-8, G-6, G-7) | §9 Phase 1/3/6 |
| Phase 3에 G-2·G-5·G-7, Phase 1에 G-6, Phase 4에 G-3, Phase 5에 G-4, Phase 6+에 G-1 | 채택. 단 G-1은 메쉬 리더로 흡수(Phase 2–3), G-4는 정정 후 | §9 |
| §14에 D11(`.inp` 읽기 채택 여부) | 번호 재배정: v2.1의 D11은 "메쉬 리더 구현"(확정). `.inp` 메쉬 읽기는 D11에 포함, 덱 번역은 범위 밖 | §14 |

---

## 5. 의도적 제외 재확인

- U-spline (특허)
- Abaqus **해석 덱** 키워드 번역(하중·step·재료), UEL/UELMAT, `.odb` 출력 — 단 `.inp` **메쉬·집합 읽기와 메쉬 쓰기는 범위 안**(§5.10)
- 솔버 결과 바이너리 포맷 쓰기(.h3d, .op2, d3plot)
- 자기접촉 (Coreform도 제외; Explicit 동역학은 Phase 6 확장 트랙으로 편입)
- 셸/보 IGA (Coreform 3D 솔리드 전용)
- missing mass, Gupta, Lindley-Yow (양쪽 범위 밖)
- Algoim 정확 형상 적분 (Phase 6 선택)
- HHT 직접 구현 (D5, Bossak)
- 볼트 프리로드·커넥터 (Coreform IGA 검증 없음)
- 키워드 가이드의 "Not supported for IGA" 16행 중 IGA 관련 핵심(U994–U999 세트, A/RF 출력, RE/PE 변수)은 우리 출력에도 두지 않는다.
- (v1.0에서 변경) **Gmsh는 사용한다**: 일반 FE 파트 메쉬 생성용, 외부 프로세스로만 (D13).

---

## 6. 정정 이력 (v1.0 → v1.1)

| ID | v1.0 내용 | 정정 | 근거 |
|---|---|---|---|
| A1 | Selective integration → "v2 §7.8 reduced/SRI" | v2 §7.8에는 SRI가 없었다. v2.1 Phase 5 잠김 방지 옵션에 추가 | `master-plan_v2.md` §7.8 원문 |
| A2 | 온도 의존 선형탄성 "동등, TableAccessor(S15)" | `부분`. 테이블은 부착·보간되지만 `LinearElastic3DLaw`·`ThermalElasticIsotropic3D`가 강성에 반영하지 않음 → S20 | 트랙 A S15 (T=100 테이블 K == 기본값 K, rtol 1e-12) |
| A3 | 고유진동수 "동등, spectra/FEAST" | `동등(방법 대체)`. FEAST 미빌드(`HasFEAST=False`), spectra는 shift 0에서 실패, 적분점 기하에서 SM consistent 질량 = 0 → Python consistent M + scipy | 트랙 A S06, 리뷰 F1 (`base_solid_element.cpp` l.629) |
| A4 | Body-fit "동등(범위 일치)" | Coreform은 다블록 정렬 hex 메쉬를 스플라인으로 변환. v2는 단일 정렬 블록뿐이었음 → v2.1 §5.10.4 | `release-notes-2026-3.md` §2.2 |
| A5 | IGA + 일반 FE 혼합 "동등" | v2는 Phase 5 Alpha 한 줄, FE 메쉬 생성·읽기 수단 없음 → `갭` → v2.1 §5.9, Phase 3 Beta | `release-notes-2026-3.md` §3.3, v2 §1.3 |
| A6 | (누락) instance-aware 메쉬 | 추가 → v2.1 §5.9.1, §5.9.3 | `scripting-guide/meshes.md`, `release-notes-2026-7.md` §2.1 |
| A7 | "Mooney-Rivlin/Yeoh/Ogden-3D는 Kratos에 없음 → C++ 승격" | ConstitutiveLawsApplication만 확인한 결과. 이 PC 빌드에 포함된 ConstitutiveModelsApplication에 Mooney-Rivlin, Yeoh, Ogden, 비압축 Neo-Hookean 모델 존재 → C++ 사유 아님, 호환성은 S19 | `kratos-src/applications/ConstitutiveModelsApplication/custom_models/elasticity_models/` (`mooney_rivlin_model.hpp`, `yeoh_model.hpp`, `ogden_model.hpp`) |
| A8 | G-1 `.inp` 읽기 서브셋 (해석 덱) | Coreform 검증 입력 파일이 크롤에 없음(`downloads/` 비어 있음, `.jnl` 404) → 재사용 원본 없음. 메쉬 읽기로 흡수, 덱 번역은 범위 밖 | `README.md` 메모, `user-guide/downloads/` |
| A9 | "Coreform에는 독립 교차검증 경로가 없다" 등 | 비공개 정보 추정 → "공개 문서에 없음"으로 | – |
| C1 (트랙 C S19) | A7 "ConstitutiveModelsApplication 존재 → C++ 불필요" | `MooneyRivlin`/`Yeoh` Python 미노출 + `QuadraturePointGeometry::Volume()` 미구현으로 실패. Neo-Hookean은 `CL.HyperElastic3DLaw`+`UpdatedLagrangian` (`6.66e-16`~`3.33e-15`), 다항 초탄성은 C++ 보완 필요 | `WHT_IGA/docs/spikes/TRACK_C.md` (16/16) |
| C2 (트랙 C S20) | A2 "TableAccessor 확인 필요" | `TableAccessor`+`CL.ThermalElasticIsotropic3D`로 $E(T)$ `6.013e-09` 성공, $\alpha(T)$·`LinearElastic3DLaw`는 셀별 `Properties` 복제 채택(`1.049e-05` @ `n_x=32`) | `WHT_IGA/docs/spikes/TRACK_C.md` (6/6) |
| C3 (트랙 C S21–S24) | 구속 Newton, ContactSM, mapper, Explicit 가정 | S21(`c=0`+`MoveMeshFlag=True`+`AndCriteria` `1e-15`), S22(`ContactSM` `GenerateFaces()` 실패→`Triangle3D3`+MPC+`SpringDamper` `5.220e-11`, Phase 5는 직접 `PairedCondition`/활성집합), S23(`nearest_element` `4.44e-16` 통과, 자유 응력은 mortar $L^2$), S24(`MechanicalExplicitStrategy` 질량 소멸 `2.78e-17`→Python Explicit+HRZ lumped $M$+$T^T f^{\text{int}}$로 $\Delta t_{\text{crit}}$ **714,131배** 복원) | `WHT_IGA/docs/spikes/TRACK_C.md` (40/40) |
