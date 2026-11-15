# WHT_IGA 구현 요청서 Phase 2 (구현 담당 에이전트용)

- 작성일: 2026-11-14
- 기준 계획서: `research-coreform-IGA/master-plan_v2.1.md` (이하 "v2.1"). 이 요청서와 v2.1이 충돌하면 **작업을 멈추고 사용자에게 묻는다.**
- 브랜치: `WHTIGA`
- 대상: 이 저장소에서 WHT_IGA Phase 2를 구현할 모든 AI 에이전트
- 요청 범위: **Phase 2 (small-cell extension + 조건수 + 메쉬 리더)**. Phase 1 산출물(S1–S4, S6, S7 선형 정적, DB)이 진입 전제이며, G1 게이트 통과는 사용자가 승인한다. Phase 3 이후는 G2 게이트를 사용자가 승인한 다음에 별도 요청으로 연다 (§9에 미리보기만 있음).

---

## 진입 상태 (Phase 2 착수 전 확인) — 미확인 항목은 TBD로 둔다

| # | 전제 | 상태 |
|---|---|---|
| P2-0 | Phase 1 산출물 존재 | 확인 (통합 담당): `wht_iga/geometry/`, `grid/`, `kratos_backend/` (+`surface/thermal.py`), `solver/` (pipeline/db/vtu/validate/trimmed), mesh/results/quad DB writer, FCM α, penalty BC-A/B, `docs/phase1/`, v01–v03 JSON |
| P2-1 | G1 게이트 | 조건부 통과 (사용자 승인 2026-11-14): V01✓ V03✓ moment-fitting✓, V02 부분 23% (extension 없음이 원인, Phase 2에서 재측정) |
| P2-2 | c≠0 거부 규칙 | `constraints.add_linear_tie` 구현 확인 필요 — T-P2.3 착수 전 해당 파일에서 확인할 것 (없으면 T-P2.3에 포함) |
| P2-3 | β 기본값 (ADR-002) | 확정: c_pen=3e5 (`docs/adr/ADR-002-penalty-beta.md`, SOLVE). T-P2.4는 V02 비교표에 집중 |

---

## 0. 먼저 읽을 것 (순서대로)

1. 이 문서 전체
2. v2.1 §7.5 (S5 small-cell extension), §7.6 (S6 penalty/Nitsche), §9 Phase 2 (산출물·게이트 G2), §10.1 V-COND/V06 행
3. v2.1 §5.2.3 (구속식 기저 변환, c = 0 규칙, 연쇄 구속 금지), §5.4 (penalty BC-A/B/C), §5.10 (메쉬 입출력, 내부 메쉬 모델)
4. `WHT_IGA/tests/verification/data/v06.json` (V06 기준값·허용오차) — V-COND는 자체 정의이며 CSV가 없다
5. 작업 시작 전 확인: `git status --short --branch`가 `WHTIGA` 브랜치인지 본다. 아니면 멈추고 묻는다.

---

## 1. 목표 한 줄

소형 셀 불안정을 기저-지지 기준 extended B-spline extension(MPC)으로 제거하고, 조건수 평탄화(V-COND)와 압력 구 octant 정확도(V06 CSV 허용오차 전량 통과)로 증명한다. 동시에 외부 메쉬 6종 리더와 내부 메쉬 모델을 갖춰 V02의 메쉬-파일 입력판을 연다.

---

## 2. 절대 규칙 (위반 시 작업 무효)

### 2.1 해야 하는 것

| # | 규칙 |
|---|---|
| M1 | 모든 Python 실행은 `micromamba run -n whtiga python ...` 또는 `D:\mamba\envs\whtiga\python.exe`로 한다 |
| M2 | Python 패키지 설치는 **uv pip만** 쓴다: `uv pip install --python D:\mamba\envs\whtiga\python.exe ...` |
| M3 | Kratos는 `E:\kratos-build\install`만 쓴다. import 후 `KratosMultiphysics.__file__`이 이 경로 아래인지 반드시 확인한다 |
| M4 | 수치 결과를 주장할 때는 실제 실행 출력(숫자)을 보고서에 붙인다 |
| M5 | 시험은 **기준값과 비교**한다(V06 CSV 허용오차, V01/V02 재통과, 참조 조립기 1e-10). 실행 성공 여부만 보는 시험은 금지 |
| M6 | 각 작업 결과를 `WHT_IGA/docs/phase2/` 아래에 기록한다 (형식 §8) |
| M7 | 파일을 덮어쓰거나 지우기 전에 내용을 확인하고, 설정 파일은 먼저 백업한다 |
| M8 | 기본값으로 대체한 입력(재료, penalty 등)은 결과에 `assumed`로 표시한다 |
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
| X9 | Phase 1 소유 파일(`WHT_IGA/wht_iga/geometry/`, `grid/`, `classify/`, `quadrature/`, `kratos_backend/` 기존 모듈, `WHT_IGA/docs/phase1/`, `WHT_IGA/tests/verification/data/v01*`–`v03*`) 수정 — 변경이 필요하면 자기 작업 파일 안에 임시 구현을 두고 보고서에 "Phase 1 공유 모듈 변경 제안"으로 적는다 |
| X10 | Rust 도구체인 설치, C++ 커스텀 앱 착수 (Phase 5 범위) |

---

## 3. 확정된 사실 (다시 조사하지 말 것)

| 항목 | 값 (v2.1 출처) |
|---|---|
| extension 판정 기준 | **기저 지지 체적비** (요소 단위 η_e 아님). large 셀: η_e ≥ θ (§7.5) |
| `small_cell_volume_ratio` | 기본 **0.2**, 허용 범위 [0, 0.5). sweep 가능 (§0.3 E6, §7.5) |
| master / extended / inactive | 지지 안에 large 셀 1개 이상 → master; 지지 ∩ Ω ≠ ∅ 이나 large 셀 없음 → extended (slave); 지지 ∩ Ω = ∅ → inactive (§7.5) |
| extension 계수 | e_ij = Π_ν ℓ_{i_ν}(j_ν), 균일 knot 인덱스에서 Lagrange 외삽 (Höllig 2003). p차 재생 + 분할단위성 보존 (§7.5) |
| THB 합성 | T = C_THB · E — 구속 그래프에 올리기 전 Python에서 미리 합성, 연쇄 구속 금지 (§5.2.3) |
| 구속 상수 | **c는 항상 0**. 0이 아닌 c는 `add_linear_tie`가 오류로 거부 (§5.2.3) |
| penalty 기본 | 사용자 `bc_penalty_value` 지정 시 그 절대값 (검증에서는 1e10). 미지정 시 β = c_pen · E_max · p² / h_min, c_pen = 1e3 (T-P2.4 sweep으로 확정, ADR-002) (§5.4) |
| FCM α 기준선 | OUT 영역 가상 강성 α = 1e-8 · E. 조건수 완화 vs extension 비교용 (§7.6) |
| 질량 행렬 | 적분점 기하에서 SM 요소 consistent 질량 = 0 (트랙 A F1). consistent M은 Python(Numba) 조립 + Tᵀ M T (§5.6, §7.9) |
| V06 기준 | `pressurized-sphere-qoi.csv`의 p별·level별 허용오차 열 전부 통과 (§9 G2) |
| 메쉬 리더 매핑 (D11) | bdf/fem = pyNastran, inp/rad/k = 자체 경량 파서, msh = meshio (§5.10.2) |
| Gmsh | **Phase 2 범위 밖** (Phase 3, D13). 외부 프로세스 호출만, import 금지 |

---

## 4. Phase 2 고유 함정 (반드시 숙지)

| # | 함정 | 올바른 처리 (v2.1 출처) |
|---|---|---|
| K20 | 요소 단위 η_e로 small-cell을 판정함 | 기저함수의 물리 지지영역 기준 (Marussig 2018). C^{p−1} 기저는 (p+1)^d 셀에 걸치므로 요소 판정은 수학적으로 틀리다 (§7.5, §0.3 E7) |
| K21 | extension 블록을 찾지 못해 계수가 폭증함 (얇은 형상) | 차수 저하 블록 → 최근접 master 상수 외삽(AgFEM형) → `issue` 코드 `THIN_FEATURE` + 세분 권고. 계수 상한 감시를 둔다 (§7.5, §11 R7) |
| K22 | slave 노드를 다른 구속의 master로 씀 (연쇄 구속) | 모든 합성(THB ∘ extension ∘ 가상노드 보간)은 Python에서 미리 계산해 한 줄의 관계식으로 만든다. `job.resolve()` 수준의 구속 그래프 검사는 Phase 3이지만, Phase 2 배치 생성기가 연쇄를 만들지 않는지 단위시험으로 강제한다 (§5.2.3, §5.9.7) |
| K23 | 재료 계면을 가로질러 extension을 적용함 | 다중재료 지지 체적비는 **합집합 영역**으로 계산한다 (재료 계면은 내부이므로 extension 대상이 아니다) (§7.5) |
| K24 | β를 고정값으로 박아두고 정확도·조건수를 동시에 잃음 | β sweep {1e2 … 1e6}·E/h로 수렴 평탄 구간을 확인하고 ADR-002로 기본값을 확정한다. 오차 O(1/β), 조건수 ∝ β·h 관계를 리포트에 기록한다 (§5.4, §7.6) |
| K25 | 메쉬 리더가 해석 덱 키워드까지 해석함 | 리더 범위는 **메쉬와 집합만** (하중·재료·step 키워드는 읽지 않음). 재료·속성 ID는 파트 구분용으로만 쓴다 (§5.10.2) |
| K26 | 단위 없는 메쉬 파일의 단위를 추측함 | `unit_system` 필수 (add_1 G-6). Nastran `CORD2R` 등 좌표계는 전역으로 변환한다 (§5.10.2) |
| K27 | 셸·빔·스프링 요소를 해석에 몰래 포함함 | 읽되 `issue`로 보고하고 해석에서 제외한다 (3D 솔리드 전용) (§5.10.2) |

---

## 5. 작업 목록 — Phase 2

각 작업은 **수락 기준을 전부 만족**해야 완료다. 완료 보고에는 실행한 명령과 출력 요약을 붙인다.

### T-P2.1 기저 지지 분석 (`wht_iga/stabilize/support.py`, 신규)

- 할 일:
  1. S3/S4 결과(셀 상태, η_e, 재료 마스크)를 입력으로, 레벨별 각 기저함수의 물리 지지영역 ∩ Ω 체적비를 계산한다 (@njit).
  2. θ (`small_cell_volume_ratio`, 기본 0.2, 범위 [0, 0.5)) 입력으로 master / extended / inactive 분류를 낸다.
  3. 다중재료: 합집합 영역 기준 (TBD — Phase 2에 다중재료 셀이 없으면 단일재료로 구현하고 다중재료 분기를 TODO가 아닌 "Phase 3" 명시와 함께 미수행으로 보고).
- 수락 기준:
  - untrimmed 균일 격자에서 전부 master, inactive 0.
  - 구멍 포함 격자에서 extended 집합이 절단 경계 밴드와 일치 (V02 형상으로 육안+VTI 확인).
  - θ sweep [0.05, 0.2, 0.35]에서 분류 변화가 단조 (θ 증가 → extended 비감소).

### T-P2.2 extension 연산자 (`wht_iga/stabilize/extension.py`, 신규)

- 할 일:
  1. 각 extended slave j마다 인덱스 공간에서 master 함수만으로 구성된 (p+1)^d 블록 I(j) 중 가장 가까운 것을 선택 (블록 중심–j 거리 최소, 바깥 법선 방향 우선).
  2. e_ij = Π_ν ℓ_{i_ν}(j_ν) Lagrange 외삽 계수 계산 (@njit).
  3. 폴백 체인: (p+1) 블록 실패 → 차수 저하 블록 → 최근접 master 상수 외삽 → 그래도 실패 시 `issue THIN_FEATURE` + 세분 권고.
  4. 계수 상한 감시 (상한값 TBD — 측정 후 ADR-002 후속으로 확정하고, 그전까지는 기록만).
- 수락 기준 (v2.1 §7.5 검증):
  - 분할단위성: Σ extended 기저 = 1 (1e-12).
  - p차 다항식 재생 (1e-10).
  - θ 변화에 대한 연속성 단위시험.

### T-P2.3 연산자 합성 + 구속식 배치 생성 (`wht_iga/operators/compose.py`, `wht_iga/operators/numbering.py`, 신규)

- 할 일:
  1. T = C_THB · E 합성 (Phase 2는 THB 없음 → C_THB = I; 합성 경로만 마련). scipy.sparse.
  2. slave 노드 × 성분 3의 `LinearMasterSlaveConstraint` 배치 생성 (c = 0 강제, 0이 아니면 오류).
  3. 연쇄 구속 검사: master가 다른 구속의 slave가 되는 경우를 생성 전에 거부하는 단위시험.
- 수락 기준:
  - T-P0.5 참조 조립기와 Kratos 경로의 Tᵀ K T 일치 ≤ 1e-10 (V01 축소판).
  - S05 성능 기준 유지 (1e5개 구속 30 s 이내 목표, Phase 0 측정값과 비교 보고).

### T-P2.4 penalty β sweep + FCM α 기준선 (`wht_iga/stabilize/fcm.py` 신규, sweep 스크립트)

- 할 일:
  1. FCM α = 1e-8 · E 기준선 구현 (OUT 영역 가상 강성).
  2. V01+V02에서 β sweep {1e2 … 1e6}·E/h: 오차 O(1/β) 거동과 평탄 구간 확인.
  3. ADR-002 초안: 기본값 확정 (자동 스케일 vs 사용자 지정 우선 — v2.1 D3 기본안은 "자동 스케일 + 사용자 지정 시 우선").
- 수락 기준:
  - sweep 그래프 + 평탄 구간 식별 + ADR-002 제안값.
  - extension 경로와 FCM α 경로의 V02 peak 오차·조건수 비교표.

### T-P2.5 조건수 도구 + V-COND (신규, 위치 TBD)

- 할 일:
  1. 행렬 추출 (Kratos `GetSystemMatrix` + `scipy_conversion_tools`, Tᵀ K T 적용) → 조건수 κ 측정 도구.
  2. V-COND 자체 시험: 구멍 중심을 이동시켜 η_min ∈ [1e-8, 1] sweep (h 4단계).
- 수락 기준 (v2.1 §9 G2):
  - extension 켠 상태에서 κ 변동 ≤ 10배 (평탄), h에 대한 log-log 기울기 −2 ± 0.3.
  - 미안정화 대비 개선 그래프 제시.
- TBD: V-COND의 모체 형상 (후보: V02 Roark 판의 구멍 중심 이동 — v2.1 §9는 "구멍 중심을 이동"이라고만 하며 구체 형상을 지정하지 않음). 형상 선택을 보고서에 명시하고 사용자에게 확인받는다. 도구 모듈 위치도 TBD (후보: `wht_iga/stabilize/cond.py`).

### T-P2.6 메쉬 리더 6종 + 내부 메쉬 모델 (`wht_iga/meshio_ext/`, 신규)

- 할 일 (v2.1 §5.10.2):
  1. `model.py`: `MeshModel{nodes, blocks[elem_type, conn, part_id], parts, node_sets, elem_sets, face_sets, frames}`. 요소 유형 정규화 (C3D8/HEXA8/CHEXA/BRICK/ELEMENT_SOLID → `hex8`, `tet4`, `tet10`, `hex20`, `wedge6`), 포맷별 노드 순서 변환표.
  2. 리더: `read_inp.py`, `read_rad.py`, `read_k.py` (자체 경량 파서), `read_nastran.py` (bdf/fem, pyNastran + `$$HMNAME` 주석), `read_msh.py` (meshio).
  3. 범위 준수: 메쉬와 집합만. 셸·빔·스프링은 읽고 `issue` 보고 후 제외. `unit_system` 필수. 파트 대응: inp `*PART`/`*INSTANCE`(변환 포함) 또는 `*SOLID SECTION` ELSET; bdf `PSOLID` PID; fem PID + 주석; rad `/PART`; k `*PART` PID + `*INCLUDE` 추종; msh physical group.
- 수락 기준:
  - 읽기→쓰기→읽기 왕복 시험 (읽기만 있어도 왕복은 `MeshModel` 직렬화로 수행): 노드·요소·집합 무손실.
  - V02를 메쉬 파일 입력(외곽 면 → immersed 형상 입력, §5.10.3 용도 3)으로 Roark 1% 이내 (G2 게이트 추가 항목).
  - 지원 외 카드·방언은 `issue`로 보고 (실패 숨기기 금지).

### T-P2.7 V06 압력 구 + V01/V02 extension 재통과 (검증, 코드 변경 없음)

- 할 일:
  1. extension 켠 상태로 V01 (p ∈ {1,2,3}, C) 재통과: 상대오차 ≤ 1e-5.
  2. V02 extension 재통과: peak ≤ 1%.
  3. V06 압력 구 octant (a = 0.9, b = 1, p_i = 1e-6): `WHT_IGA/tests/verification/data/v06.json`의 CSV 허용오차 전 셀 통과. 목표 p = 2 최세 u L2 ≤ 1e-4, σrr ≤ 0.5%.
- 수락 기준: G2 게이트 3개 항목 전부. 결과는 `WHT_IGA/docs/verification/V06.md` + `V-COND.md` + `V01-ext.md`/`V02-ext.md`로 자동 생성 (형식 §8).

---

## 6. 산출 모듈 목록 (신규 파일만)

| 모듈 | 작업 | 비고 |
|---|---|---|
| `wht_iga/stabilize/__init__.py`, `support.py` | T-P2.1 | @njit 커널, 순수 함수 |
| `wht_iga/stabilize/extension.py` | T-P2.2 | @njit 커널 + THIN_FEATURE issue |
| `wht_iga/stabilize/fcm.py` | T-P2.4 | α = 1e-8·E 기준선 |
| `wht_iga/operators/__init__.py`, `compose.py`, `numbering.py` | T-P2.3 | scipy.sparse, c=0 강제 |
| 조건수 도구 (위치 TBD) | T-P2.5 | 후보 `wht_iga/stabilize/cond.py` |
| `wht_iga/meshio_ext/__init__.py`, `model.py`, `read_inp.py`, `read_rad.py`, `read_k.py`, `read_nastran.py`, `read_msh.py` | T-P2.6 | pyNastran/meshio는 `mesh` extra |
| `WHT_IGA/tests/unit/test_stabilize_*.py` | T-P2.1–P2.4 | Kratos 불필요 |
| `WHT_IGA/tests/verification/data/` 기존 JSON 사용 | T-P2.7 | v06.json이 기준값 원천; 새 JSON 작성 금지 |
| `WHT_IGA/docs/verification/V06.md`, `V-COND.md`, `V01-ext.md`, `V02-ext.md` | T-P2.7 | 형식 §8 |

---

## 7. 코딩 규칙 (v1 §7과 동일 — Phase 2 추가만 아래에)

- Python 3.12, 타입 힌트, `from __future__ import annotations`.
- Numba 커널: 순수 함수, NumPy 배열 입출력, Python 객체·딕셔너리 금지, `@njit(cache=True)` (병렬은 측정 후 `parallel=True`).
- 단위는 입력 단위 그대로 (무차원 일관 단위계). 단위 변환 코드 금지.
- 좌표·인덱스 규약: 셀 인덱스 (i, j, k) → 선형 `i + nx*(j + ny*k)` (`grid/` 모듈 docstring 규약 준수).
- 부동소수 비교는 허용오차 명시 (`np.testing.assert_allclose(rtol=..., atol=...)`).
- 출처 표기: 알고리즘에 논문 번호 주석 (Höllig 2003, MarussigHughes 2018 등).
- 파일 400줄 초과 시 책임별로 분할.
- **추가 (Phase 2)**: extension/합성 연산자는 희소 COO → CSR로만 다루고密행렬을 만들지 않는다. 1e5 구속 배치 생성 30 s 목표를 측정·보고한다.

---

## 8. 보고 형식

### 8.1 작업 기록 (`docs/phase2/P2NN.md`)

```markdown
# P2NN <제목>
- 작업: (T-P2.x)
- 결과: 완료 / 부분 / 실패
- 근거: (실행 명령, 핵심 출력 숫자, 수락 기준별 통과/실패)
- 환경: Python x.y.z, Kratos 10.4.3 (경로), 날짜
- 결론과 영향: (G2 게이트에 미치는 영향, Phase 3으로 넘기는 항목)
- 재현: pytest 경로 또는 실행 명령
```

### 8.2 검증 리포트 (`docs/verification/<ID>.md`, v2.1 §10.6)

- 케이스별 자동 생성: 입력 요약, 기하 오차, 수렴 표·그래프, 기준값 대비 오차, 통과 여부, 실행 시간, DOF 수.
- V06은 CSV 허용오차 셀별 통과표를 포함한다. V-COND는 κ 평탄 그래프 + 미안정화 대비 그래프를 포함한다.

### 8.3 작업 완료 보고 (작업마다)

- 작업 ID, 변경·생성 파일 목록, 실행 명령, 수락 기준별 통과/실패와 근거 숫자, 남은 문제.
- 실패를 성공처럼 쓰지 않는다. 미수행 항목은 "미수행"으로 적는다.

---

## 9. 정의: 이번 요청의 완료 (Definition of Done)

1. T-P2.1 – T-P2.7 모든 수락 기준 통과, 또는 실패 항목이 근거와 함께 보고됨
2. `micromamba run -n whtiga python -m pytest WHT_IGA/tests -m "not kratos"` 통과 (신규 단위시험 포함)
3. Kratos 필요 시험 (`-m kratos`) 결과 첨부 (실패가 있으면 그대로 보고)
4. 변경 파일에 TODO 스텁, skip, xfail 없음 (있으면 사유와 함께 보고)
5. G2 게이트 증거 제출 → **사용자 승인 대기**. 승인 전 Phase 3 착수 금지.

### 9.1 Phase 3 이후 미리보기 (이번 요청 범위 아님)

| Phase | 핵심 산출물 | 게이트 (v2.1 §9) |
|---|---|---|
| 3 | 멀티셀 공유 장, tie (BC-C), RP/분배/kinematic coupling, amplitude, 다단계, 가중치 0 probe, THB 세분화, VTKHDF proxy 격자, probe JSON, L2 노름 도구 + 어셈블리/메쉬 입출력/Gmsh/body-fit | V03 moment 5%, V04 (merge–tie ≤ 1e-6), V05, V02 국소 ≤ 0.2%, E01, E02 |
| 4 | Python consistent M, scipy 고유치, Python Bossak 선형 과도, 모드 중첩 과동/조화, 스펙트럼 case-file 리더 | V09, V10, V11, V12 |
| 5 | NLGEOM, 초탄성, J2, C++ MixedUP0, 접촉 | V07, V08, E03, E04, E05 |
| 6 | 성능 (C++ 이전, 반복 솔버, 병렬), Explicit 트랙 | 1M DOF < 30분 |

---

## 10. 멈추고 사용자에게 물을 조건

다음 중 하나라도 생기면 **추측으로 진행하지 말고** 상황, 시도한 것, 선택지를 정리해 묻는다.

1. extension 블록 탐색이 얇은 형상에서 대량으로 실패해 `THIN_FEATURE`가 지배적이 될 때 (세분 권고 임계는 사용자 결정)
2. T-P2.3 합성 경로에서 THB가 없다는 가정이 깨질 때 (Phase 3 THB 조기 착수 여부)
3. β sweep에 평탄 구간이 없을 때 (penalty 경로 자체의 재검토 — Nitsche 조기 착수는 범위 밖이므로 반드시 묻는다)
4. V-COND 모체 형상 선택이 필요할 때 (TBD — 후보는 V02 판, 결정은 사용자)
5. V06 CSV 허용오차와 v2.1 §9 G2 목표치가 다를 때 (CSV를 우선하고 보고 — 계획서 수정은 사용자 승인 후)
6. 메쉬 파일 방언(Nastran 필드 형식, LS-DYNA `*INCLUDE`, OptiStruct 확장, Radioss 버전차)에서 왕복 무손실이 깨질 때 (지원 외 카드 보고 vs 파서 확장 — 범위 판단 필요)
7. v2.1 내용이 Kratos 실제 동작과 다를 때 (v2.1은 수정하지 말고 보고)
8. 라이선스가 불분명한 코드·데이터를 가져와야 할 때 (pyNastran BSD-3, meshio MIT 외)
9. 작업량이 예상(3주)의 2배를 넘을 것으로 보일 때

---

## 11. 참고 위치

| 내용 | 위치 |
|---|---|
| 계획서 | `research-coreform-IGA/master-plan_v2.1.md` (§7.5, §7.6, §9 Phase 2, §10.1, §11 R7, 부록 A) |
| V06 기준값 JSON (직접 읽을 것) | `WHT_IGA/tests/verification/data/v06.json` |
| V-COND 정의 (CSV 없음, 계획서 §9·§10.1이 원천) | v2.1 §9 G2 + §10.1 V-COND 행 |
| Phase 1 산출물 (읽기 전용) | `WHT_IGA/wht_iga/geometry/`, `grid/`, `classify/`, `quadrature/`, `kratos_backend/`, `WHT_IGA/docs/phase1/` |
| Phase 1 기준값 JSON (읽기 전용) | `WHT_IGA/tests/verification/data/v01*`–`v03*` |
| 메쉬 리더 선행 패턴 (참고만, import 금지) | `aieng/src/aieng/topopt/adapters/kratos.py` |
| Kratos 행렬 추출 바인딩 | `kratos/python/add_model_part_to_python.cpp` (부록 A-4/A-5), `scipy_conversion_tools` |
