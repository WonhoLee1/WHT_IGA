# master-plan v2 Addendum 2 — 결합 문서·add_1 검토, Coreform 워크플로 대응 설계, 메쉬 전략

- 문서 버전: v2-add_2 (2026-11-14)
- 모문서: `master-plan_v2.md`, `master-plan_v2-add_1.md` (둘 다 수정하지 않음)
- 근거: `user-guide/user-guide/workflow-overview.md`, `scripting-guide/meshes.md`, `release-notes-2026-3/7/9/11.md`, `immersed-iga-overview.md`, 트랙 A spike 결과(`WHT_IGA/docs/spikes/`), Kratos 소스(`D:\PythonCodeStudy\kratos-src`)
- 상태: **계획 보완안.** §5의 편입 제안은 사용자 승인 후 v2.1에 반영한다.

---

## 1. 결합 문서 검토 (`user-guide_combined.md/.html`)

| 항목 | 결과 |
|---|---|
| 원천 범위 | `user-guide/`의 md 82개 **전부 포함** (디스크 목록과 1:1 대조, 누락 0). 제외된 29개는 `generated-downloads/`의 CSV 변환 md(의도적 제외) |
| 목차·링크 | md: 앵커 83, 목차 링크 82, 깨진 링크 0. html: `<section>` 82, 내부 링크 164, 깨진 링크 0 |
| **그림** | **0개.** 원본 html에는 `<img>` 928개가 있었으나 결합본에는 그림이 없다. 예: 워크플로 개요의 "Figure 1"(전체 흐름도), 검증 매뉴얼 수렴 그래프 |
| **수식** | LaTeX가 원문 그대로(`\(\sigma_{rr}\)`) 들어 있고 html에 MathJax/KaTeX가 없다 → 브라우저에서 수식이 렌더링되지 않음 |
| 검증 CSV 기준값 | 결합본에 없음(제외 범위). 수치 검증은 계속 `generated-downloads/`의 CSV를 직접 쓴다 |

권고: 텍스트 검색·에이전트 입력용으로는 충분하다. 사람이 읽는 용도라면 (1) 그림을 `user-guide/assets`에서 상대경로로 연결하고, (2) html에 KaTeX(cdnjs) 한 줄을 추가하면 된다.

---

## 2. add_1 검토

전반적으로 대조 방법과 갭 목록(G-1–G-8)은 쓸모 있다. 다만 다음은 사실과 다르거나 트랙 A 결과 이전 정보다.

| # | add_1 내용 | 문제 | 수정 |
|---|---|---|---|
| A1 | §1.1 Selective integration → "v2 §7.8 reduced/SRI" | v2 §7.8에는 SRI가 없다(B-bar/F-bar/VMS만) | 판정 `부분`. Phase 5 잠김 방지 옵션에 SRI 명시 추가 |
| A2 | §1.4 온도 의존 선형탄성 "동등, TableAccessor(S15)" | 트랙 A S15: 테이블은 붙지만 **강성에 반영되지 않음** | 판정 `부분`. 리뷰 F3(Accessor spike 또는 요소별 Properties) |
| A3 | §1.3 고유진동수 "동등, spectra/FEAST" | 트랙 A S06: FEAST 미빌드, spectra는 shift 0에서 실패, **consistent 질량이 0** | Python consistent M + Kratos 추출 K + scipy `eigsh` (리뷰 F1, C++ 불필요) |
| A4 | §1.1 Body-fit "동등(범위 일치)" | Coreform은 **기존 정렬 all-hex 메쉬(다블록 mapped)를 스플라인으로 변환**한다. v2는 단일 정렬 블록뿐 | 판정 `부분` → §4 M2 |
| A5 | §1.2 IGA + 일반 FE 혼합 "동등" | v2는 Phase 5 Alpha 한 줄뿐이고, 일반 FE 메쉬를 **만들거나 읽는 수단이 없다**(D10은 Gmsh 제외가 기본) | 판정 `갭` → §4 M3, M4 |
| A6 | (누락) instance-aware 메쉬 | Coreform 2026.7: dependent/independent instance별 메쉬 소유 규칙, dependent instance는 Production | `갭` → §4 M5 |
| A7 | G-4 "Mooney-Rivlin/Yeoh/Ogden-3D는 Kratos에 없음 → C++ 승격" | ConstitutiveLawsApplication만 확인했다. 이 PC의 빌드에 포함된 **ConstitutiveModelsApplication에 Mooney-Rivlin, Yeoh, Ogden, 비압축 Neo-Hookean 모델이 있다** (`custom_models/elasticity_models/`) | C++ 승격 사유에서 제외. 단 이 앱은 로컬 복원·수정분(166파일)이므로 SM 요소 + 적분점 기하와의 호환성을 spike로 확인 |
| A8 | G-1 `.inp` 읽기 서브셋 | 사용자 결정 1과 직접 충돌하지는 않지만, Coreform 검증 입력 파일은 크롤에 없다(`downloads/` 비어 있음, `.jnl` 404) → 재사용할 원본이 없다 | §4 M3의 **메쉬 파일 읽기**(inp 포함 5종)로 흡수. 해석 덱 키워드 번역은 계속 범위 밖 |
| A9 | S-항목 일부 | "Coreform에는 독립 교차검증 경로가 없다" 등은 비공개 정보에 대한 추정 | 표현을 "공개 문서에 없음"으로 |

---

## 3. Coreform 워크플로 → WHT_IGA 워크플로 대응

### 3.1 Coreform 단계 (`workflow-overview.md`)

1. **Model setup** (Abaqus/CAE): 형상, 재료, step, 하중. 메쉬는 Abaqus 일반 메쉬 또는 Coreform IGA 메쉬. IGA가 있으면 Coreform Job 도구로 `<job>.inp` + `<job>.cf` 생성.
2. **Compute trimmed mesh** (`coreform_iga_mesh`, 라이선스): 절단 스플라인 메쉬, 적분, 추출·확장 연산자, function splitting → `<job>.sql`.
3. **Initialize user element data & translate keywords** (interop): `<job>.inp` → `<job>_source.inp` 보존, IGA 대상 키워드 번역, `*USER ELEMENT` 덱 `<job>_iga_geom.inp`, 사용자 요소 자료를 `<job>.sql`에 추가.
4. **Abaqus/Pre → Abaqus/Standard** (Coreform 사용자 요소, `uniquelibs`): 매 반복 Abaqus가 Coreform에 K, M, f 요청. 해석 중 `<job>_iga_results.sql`, `<job>_iga.odb`(선형 hex-tet), `probe_data.json` 생성.
5. **Postprocessing**: Abaqus/CAE(`.odb` = 일반 요소 + IGA 노드 + 경계 표면 메쉬, `_iga.odb`), ParaView alpha reader(`_iga_results.sql`), JSON probe.

### 3.2 WHT_IGA 대응

| # | Coreform 단계 | WHT_IGA 단계 | 산출물 | 담당 모듈 |
|---|---|---|---|---|
| W0 | Abaqus/CAE 모델 설정 | **Python API로 모델 정의**: Part(형상 출처: CAD/STL/메쉬 파일), 파트별 Discretization(§4.2), 재료, Step, 하중·BC, 상호작용(tie/contact/coupling), Probe | 메모리 내 `Model` | `wht_iga.api` |
| W1 | Job 생성 → `.inp` + `.cf` | `job.write_input()`: spec 검증(`validate()`) + 직렬화. 참조한 CAD/메쉬 파일의 해시와 경로 기록 | `<job>.spec.json` | `api/job.py`, `api/spec.py` |
| W2 | `coreform_iga_mesh` → `.sql` | `job.build_mesh()`: **파트별 병렬** S1–S6. IGA 파트는 격자·분류·절단 적분·extension·표면층, 일반 FE 파트는 메쉬 읽기·검사. 캐시 키 = 형상 해시 + 파트 이산화 spec | `<job>.wiga.sqlite` (파트별 블록) | `geometry/`, `grid/`, `classify/`, `quadrature/`, `stabilize/`, `meshio_ext/`(§4.4) |
| W3 | interop (키워드 번역, 사용자 요소 덱) | `job.resolve()`: **대상 해석(resolution)**. 하중·BC·출력 요청의 대상 세트를 IGA 파트/일반 FE 파트로 분해(Coreform의 "mixed native+IGA → IGA-only/native-only 세트 분리" 규칙과 동일), tie/contact 표면층 생성, 해석 불가 조합은 거부+진단 | `<job>.resolved.json`(진단 리포트 포함) | `api/resolve.py`(신규) |
| W4 | Abaqus/Pre + Standard + 사용자 요소 | `job.submit()`: Kratos ModelPart 구축(IGA 파트 = 적분점 기하 + SM 요소 + 구속식, 일반 FE 파트 = Kratos 일반 요소), Kratos 전략 실행. 해석 중 프레임 단위로 결과 기록, 진행 콜백(add_1 G-2) | `<job>_iga_results.sqlite`, `<job>_probe_data.json`, `<job>.log` | `kratos_backend/` |
| W5 | Abaqus/CAE, ParaView alpha reader | `job.results` API + 내보내기: **IGA proxy 격자 + 일반 FE 요소 + 경계 표면층을 한 파일에**(Coreform 2026.11 unified ODB에 대응), probe JSON, 선택 시 외부 포맷 내보내기(§4.4) | `<job>.vtkhdf` / `.vtu`, 내보내기 파일 | `post/` |

차이점 (우리 쪽이 단순한 부분):
- 라이선스 대기 단계가 없다. W2–W4가 한 프로세스에서 이어질 수 있다(캐시가 있으면 W2 생략).
- 키워드 번역이 없다. 대신 W3에서 같은 의미 규칙(대상 분해, 거부+진단)을 spec 수준에서 적용한다.
- Coreform은 매 반복 Abaqus가 사용자 요소를 호출하지만, 우리는 Kratos가 직접 조립한다.

---

## 4. 메쉬 기능 설계

### 4.1 Coreform이 제공하는 메쉬 기능과 v2 반영 상태

| # | Coreform 기능 | 근거 | v2 반영 | 판정 |
|---|---|---|---|---|
| M1 | Immersed rectilinear IGA 메쉬 | mesh API, 검증 전체 | Phase 1, §7.2 | 반영 |
| M2 | Body-fit IGA: **기존 정렬 all-hex 메쉬를 스플라인으로 변환** (2026.3: fully structured, mapped all-hex만, bottom-up·sweep 불가) | `release-notes-2026-3` §2.2, Cook's membrane(E04)·vibrating beam(V10) 비교 | 단일 정렬 블록만 (Phase 3 Alpha) | **부분** |
| M3 | 일반(native) FE 메쉬 파트 | `workflow-overview` "Mesh settings" | 없음 | **갭** |
| M4 | IGA 파트 + 일반 FE 파트 혼합, matching/nonmatching 계면 tie | `release-notes-2026-3` §3.3, 2026.9 "IGA 경계 세트와 일반 FE 세트 분리" | Phase 5 Alpha 한 줄 | **갭** |
| M5 | Instance-aware 메쉬 소유 (part 또는 independent instance, dependent instance는 part 메쉬 재사용, 위치·방향 추종) | `scripting-guide/meshes.md`, `release-notes-2026-7` §2.1 | 없음 | **갭** |
| M6 | Multi-cell(불리언 병합) 단일 스플라인 장 | 2026.7 Production | Phase 3 | 반영 |
| M7 | 국소 세분화 | 2026.11 Beta | Phase 3 | 반영 |
| M8 | 메쉬 미리보기 (해석 전 이산화 확인) | 2026.3 §2 "preview the generated discretization" | 결과 proxy만, 사전 미리보기 없음 | **부분** |
| M9 | 병렬 메쉬 단계 | `num_iga_cpus` | Phase 6 | 반영 |
| M10 | STL·격자 구조 입력 (Beta) | 나의조사·릴리스 | STL 입력만 언급 | 부분 |

### 4.2 이산화 단위: 파트(Part)

- **원칙: 이산화는 파트 단위로 정한다.** Coreform의 "Mesh owner = part 또는 independent instance, 소유자당 메쉬 1개" 규칙을 따른다.
- 파트마다 하나를 고른다:

| 이산화 종류 | 입력 | 용도 |
|---|---|---|
| `ImmersedRectilinearMesh` | CAD(STEP/BREP), STL, **또는 메쉬 파일의 외곽 표면** | 기본. 복잡 형상, 메쉬 생성이 어려운 부품 |
| `BodyFitMesh` | 정렬 all-hex 메쉬 파일, 또는 블록 형상 | 경계 정확도가 중요한 단순 형상, Coreform 비교 케이스(E04, V10) |
| `NativeMesh` | 메쉬 파일(inp/fem/bdf/rad/k/msh) | 이미 메쉬가 있는 부품, 볼트·패드 등 표준 FE가 편한 부품 |
| 없음(강체/참조) | – | RP 결합 대상 등 |

### 4.3 파트를 나눌지 합칠지 — 전략

| 상황 | 권장 | 이유 |
|---|---|---|
| 서로 다른 재료가 **완전 접합**된 한 물체 (예: glass-to-metal seal의 유리·핀·헤더, bimetallic) | **한 파트 + multi-cell** (공유 배경 격자, 재료 셀별 적분) | 계면이 스플라인 장 안에 있어 연속성이 정확하고 tie 오차가 없다(Coreform V04: merge와 tie 차 1e-7). 파트 하나라 DOF가 적다 |
| 접촉·미끄럼·분리되는 부품 (door-lock의 hinge/link/pin) | **파트 분리** | 접촉은 서로 다른 장 사이에서만 의미가 있다 |
| 크기·형상 스케일이 크게 다른 부품 (큰 하우징 + 작은 핀) | **파트 분리**, 파트별 격자 크기·Frame | 한 격자로 묶으면 작은 부품 해상도에 맞춰 큰 부품 DOF가 폭증 |
| 축이 기울어진 길쭉한 부품 | 파트별 **object-aligned Frame** | 셀이 형상 방향과 맞아 절단 셀·소형 셀이 줄어든다 |
| 볼트 체결 등 강하게 묶였지만 별개 부품 | 파트 분리 + tie (또는 Phase 5 이후 볼트 프록시) | 실제 조립 구조 유지, 결과 해석이 쉬움 |
| 같은 부품이 여러 개 (패턴 배치) | **파트 1개 + dependent instance 여러 개** | 메쉬 단계는 한 번만 하고 강체 변환만 적용(§4.5) |
| 이미 FE 메쉬가 있는 부품 | `NativeMesh` 파트, IGA 파트와 tie | 기존 자산 재사용 |

기본 규칙(자동 제안): 입력 파일의 물체(solid / 메쉬 파트 ID)마다 파트 1개 → 같은 재료 그룹이 아닌데 면을 공유하면 "multi-cell로 합칠지" 제안 → 사용자가 확정. `validate()`가 다음을 issue로 낸다: 파트 간 공유 면에 tie/contact 미지정, 파트 내부 재료 셀 체적 합 ≠ 파트 체적.

### 4.4 메쉬 파일 입출력 (inp, fem, bdf, rad, k + msh)

**범위: 메쉬와 집합만 읽는다.** 노드, 요소(유형별), 파트/속성 ID, 노드·요소·면 집합, 좌표계. 해석 덱의 하중·재료·step 키워드는 읽지 않는다(사용자 결정 1, 단 재료·속성 ID는 파트 구분용으로만 사용).

| 포맷 | 구조 | 파트 대응 | 구현 방안 | 라이선스 |
|---|---|---|---|---|
| `.inp` (Abaqus) | `*NODE`, `*ELEMENT, TYPE=`, `*NSET/*ELSET`, `*PART/*INSTANCE/*ASSEMBLY` | `*PART`(있으면) 또는 `*SOLID SECTION`의 ELSET | 자체 경량 파서 (meshio는 `*PART/*INSTANCE` 변환·instance 변환을 충분히 다루지 않음) | 자체 |
| `.bdf` / `.nas` (Nastran) | `GRID`, `CHEXA/CTETRA/CPENTA`, `PSOLID`, `SET`, 고정·자유·대형 필드 | `PID`(`PSOLID`) | **pyNastran** BDF 리더(성숙, OptiStruct 일부 카드 포함) — 선택 의존 | BSD-3 |
| `.fem` (OptiStruct) | Nastran 벌크 카드 + OptiStruct 확장(`$$HM` 주석, `PSOLID` 등) | `PID`, HyperMesh 컴포넌트 주석 | pyNastran(`punch`/관대 모드) + `$$HMNAME COMP` 주석 파싱으로 컴포넌트 이름 복원 | BSD-3 |
| `.rad` (Radioss starter) | `/NODE`, `/BRICK/<part_id>`, `/TETRA4`, `/TETRA10`, `/PART`, `/GRNOD`, `/SURF` | `/PART` | 자체 파서 (블록 키워드 + 고정 폭 필드) | 자체 |
| `.k` / `.key` (LS-DYNA) | `*NODE`, `*ELEMENT_SOLID(_TET4TOTET10)`, `*PART`, `*SET_NODE/SEGMENT`, 고정 폭 10/16·`,` 구분 | `*PART` PID | 자체 파서 (`*INCLUDE` 추종) | 자체 |
| `.msh` (Gmsh) | `$Nodes`, `$Elements`, `$PhysicalNames` | physical group | meshio | MIT |

공통 설계:
- 모든 리더는 하나의 **내부 메쉬 모델**로 변환한다: `MeshModel{nodes, blocks[elem_type, conn, part_id], parts{id: name}, node_sets, elem_sets, face_sets(elem, local_face), frames}`. 요소 유형은 정규화한다(C3D8/HEXA8/CHEXA/BRICK/ELEMENT_SOLID 8절점 → `hex8`, 2차 tet → `tet10` 등, 노드 순서 변환표 포함).
- 단위·좌표계: 파일에는 단위가 없다 → `unit_system`을 사용자가 지정(add_1 G-6). Nastran `CORD2R` 등 좌표계는 전역 좌표로 변환.
- 지원 외 요소(셸, 빔, 스프링)는 읽되 `issue`로 보고하고 해석에서는 제외(Coreform도 3D 솔리드 전용).
- 쓰기(내보내기): 같은 내부 모델에서 5개 포맷으로 **메쉬 + 결과 노드값**을 내보낸다. IGA 파트는 proxy hex-tet 격자로 쓴다. 용도는 HyperMesh, ANSA, LS-PrePost 등 외부 도구에서의 확인·재사용. 결과 포맷(.h3d, .op2, d3plot)은 범위 밖(H3D는 v2 선택 항목 유지).

메쉬 파일의 세 가지 쓰임:
1. `NativeMesh` 파트 (M3)
2. `BodyFitMesh` 입력: 정렬 all-hex 블록을 감지해 스플라인으로 변환 (M2)
3. **`ImmersedRectilinearMesh`의 형상 입력**: 솔리드 메쉬의 외곽 면(경계 삼각형/사각형 분할)을 tessellation으로 바로 사용 → CAD 없이 기존 FE 모델을 IGA로 재해석. 품질 검사(watertight, 자기교차)는 S1과 같은 규칙.

### 4.5 Instance 모델 (M5)

- `Part`(형상 + 이산화) / `Instance`(Part 참조 + 강체 변환) / `Assembly`(Instance 목록).
- 소유 규칙은 Coreform과 같다: 파트의 instance가 모두 dependent면 **파트가 메쉬 1개 소유**, 모두 independent면 **instance마다 메쉬**. 한 패밀리에서 두 방식 혼용은 `validate()` 오류.
- dependent instance: W2 메쉬 단계는 파트당 1회. W4에서 노드 좌표에 변환을 적용해 instance마다 별도 자유도로 생성. Immersed 격자의 Frame은 파트 좌표계 기준이므로 강체 변환과 함께 이동(재계산 없음).
- independent instance: instance 위치에서 메쉬를 따로 만든다(전역 축 정렬 격자가 필요한 경우 등).

### 4.6 Body-fit 변환 (M2)

- 입력: all-hex 메쉬. 정렬(structured) 블록 감지 = 각 블록이 (i, j, k) 격자 위상을 가지는지 확인(노드 차수, 면 인접 그래프). 2026.3의 Coreform과 같이 **mapped 블록만** 지원, sweep·bottom-up 불가.
- 블록마다 노드 격자를 제어점 망으로 보간·근사해 `NurbsVolumeGeometry` 생성(차수 p, 연속성 C^{p−1}). 블록 사이는 공유 경계 제어점이 일치하면 C⁰ 결합(같은 노드), 아니면 tie.
- Kratos 쪽은 immersed와 같은 경로(적분점 기하 + SM 요소)라 절단 적분·extension이 필요 없다.

### 4.7 혼합 모델 (M4)

- 일반 FE 파트: Kratos SM 일반 요소(`SmallDisplacementElement3D8N`, `3D10N` 등)를 같은 ModelPart에 생성. 재료·step·전략 공유.
- IGA–FE 계면: v2 §5.4 BC-B/BC-C 방식. FE 계면 노드를 slave로, IGA 가상노드(스플라인 slave)와 MPC 또는 penalty 스프링으로 결합. matching(노드 위치 일치)과 nonmatching(투영 후 보간) 모두 같은 틀.
- 대상 분해 규칙(W3): 하중·BC·출력 대상이 IGA와 FE에 걸치면 둘로 나눠 각각 적용, 나눌 수 없으면 거부 + 진단(Coreform 키워드 가이드의 정책과 동일).
- 결과: 일반 FE 요소 결과와 IGA proxy 결과를 같은 VTKHDF 파일에 기록(unified output).

### 4.8 메쉬 미리보기 (M8)

- `job.build_mesh(preview=True)` 또는 `mesh.preview(part)`: 해석 없이 배경 격자, 셀 상태(OUT/IN/CUT), 소형 셀(extension 대상), 표면층, 국소 세분 영역을 VTKHDF로 기록. 사용자가 ParaView에서 확인 후 해석 진행.

---

## 5. master-plan v2 편입 제안

| 항목 | 내용 | Phase | 비고 |
|---|---|---|---|
| P-1 | 파트 단위 이산화 + Part/Instance/Assembly 모델(§4.2, §4.5) | 1 (데이터 모델), 3 (instance) | API 기초라 일찍 |
| P-2 | 메쉬 파일 리더 6종(inp, bdf, fem, rad, k, msh) + 내부 메쉬 모델 + 요소 정규화 (§4.4) | 2–3 | pyNastran은 선택 의존(`extras: mesh`) |
| P-3 | 메쉬 외곽 면 → immersed 형상 입력 (§4.4 쓰임 3) | 2 | 리더 완료 직후. CAD 없는 모델 지원 |
| P-4 | 메쉬 미리보기 (§4.8) | 1 | 디버깅에도 필요 |
| P-5 | Body-fit 변환 (정렬 hex 블록 → 스플라인) (§4.6) | 3 | E04, V10 비교 케이스 대상 |
| P-6 | 일반 FE 파트 + IGA–FE tie + 대상 분해 규칙 (§4.7, W3) | 3 (Phase 5에서 앞당김) | v2 Phase 5 Alpha를 Phase 3 Beta로 |
| P-7 | 메쉬 + 노드 결과 내보내기 5종 (§4.4) | 3 | |
| P-8 | `job.resolve()` 단계와 진단 리포트 (W3) | 1 | |
| P-9 | add_1 수정 A1–A9 반영 (특히 A7: ConstitutiveModels 재료 spike) | 0–1 | |
| P-10 | 결합 문서 보완(그림 연결, KaTeX) | – | 선택 |

### 5.1 열린 결정 (추가)

| # | 결정 | 기본안 |
|---|---|---|
| D11 | 메쉬 리더 의존: pyNastran(BSD-3) 사용 vs 전부 자체 파서 | bdf/fem은 pyNastran, 나머지 자체 |
| D12 | 내보내기 범위: 메쉬 + 노드 결과만 vs 요소 결과·응력까지 | 메쉬 + 노드 결과(U, NT11) + 요소 평균 응력 |
| D13 | Gmsh 사용(D10 기본 제외) — 일반 FE 파트를 우리가 생성할지 | 제외 유지. 일반 FE 메쉬는 외부 도구에서 만들어 파일로 받는다 |
| D14 | 파트 자동 병합 제안을 기본으로 켤지 | 제안만 하고 확정은 사용자 |
