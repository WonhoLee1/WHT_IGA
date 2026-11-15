# 트랙 A 리뷰 (통합 담당, 2026-11-14)

## 판정

**수락.** 정적 경로 R1(셀 단위 적분점 기하 + SM 요소 + 구속식)은 성립한다. G0 전에 질량 행렬 운용안 확정(C++ 불필요)과 코드 보강 2건이 필요하다.

## 재현 확인

- 트랙 A spike 시험 8개 파일 독립 재실행: **27개 통과, 실패 0, skip 0** (99.3초, JUnit 요약).
- 공유 모듈(`geometry.py`, `_bootstrap.py`, `grid/`, `reference/`, `conftest.py`, `pyproject.toml`, `uv.lock`) 수정 없음 확인. 소유권 규칙 준수.
- 수치는 기준값(참조 조립기, 해석해)과 비교하는 형태. 실패를 고정하는 특성화 시험(S06, S13, S15)도 동작이 바뀌면 깨지도록 작성됨.

## 발견 사항 (심각도순)

### F1 [높음, 사용자 결정 필요] S06 consistent 질량 행렬 = 0 — 근본 원인 확인

- 원인: `BaseSolidElement::CalculateMassMatrix`(SM `base_solid_element.cpp` l.629)가 `UseGeometryIntegrationMethod()`가 true이면 `IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation`으로 **한 단계 높은 Gauss 방법**을 고른다. `UseGeometryIntegrationMethod()`는 기본 클래스에서 `true` 고정(가상 함수, `base_solid_element.h` l.282). 적분점 기하는 기본 방법에만 점을 갖고 있어 높은 방법의 점 목록이 비어 → 합 0.
- Python에서 우회 불가 (가상 함수 오버라이드 필요). lumped 질량은 기본 방법을 쓰므로 정상.
- 영향: **D5(동역학 = Kratos Bossak)와 Kratos 고유치 전략이 consistent 질량을 쓸 수 없다.**
- **C++ 없이 해결 가능** (수정, 2026-11-14): 행렬 직접 경로(v2 R2)를 쓴다. K는 Kratos에서 `GetSystemMatrix()`로 추출(S04에서 참조 조립기와 1e-16 일치 확인), consistent M은 Python(Numba)으로 조립(S06에서 검증), 풀이는 scipy.
- lumped 대 consistent 측정 (V09 보, b=h=0.1, L=1, E=200 GPa, ρ=8000, ν=0; 참조 조립기 + scipy `eigsh`, 스크립트 scratchpad `lumped_vs_consistent.py`):

  | 격자 | consistent f1 | lumped f1 | 차이 |
  |---|---|---|---|
  | p=2, 10×2×2 | 80.266 Hz | 79.755 Hz | −0.64% |
  | p=2, 40×4×4 | 80.246 Hz | 80.194 Hz | −0.07% |
  | p=2, 80×8×8 | – | – | −0.02% |
  | p=3, 20×2×2 | 80.246 Hz | 80.036 Hz | −0.26% |

  - consistent 수렴값 80.246 Hz = Euler–Bernoulli 80.770 Hz 대비 −0.65% (3D 효과). Coreform 보고 81.8 Hz(+1.2%).
  - 상위 모드(20개, 최고 약 6.5 kHz): lumped 오차 최대 7.93%(40×4×4) → 2.33%(80×8×8).
- 운용 결정안 (C++ 불필요):

  | 해석 | 방법 |
  |---|---|
  | 고유치, 응답 스펙트럼, 랜덤 응답 (V09) | Kratos 추출 K + Python consistent M + scipy `eigsh` |
  | 선형 과도 (V10) | 같은 행렬 + Python Bossak 시간 루프 (D5의 Bossak 알고리즘을 Python에서 실행). 대안: Kratos Bossak + lumped |
  | 비선형 동역학 | Kratos Bossak + lumped (1차 모드 오차 0.1% 수준). 상위 모드가 중요하면 Python Newton 루프(매 반복 Kratos `Build`) + Python M |
  | C++ 요소 | 위 방법이 정확도·성능 기준을 못 맞출 때만 |

### F2 [중간] S05 구속식 상수 규칙의 일반화와 안전장치

- 원인: Block builder가 매 `Build`마다 구속 상수 c를 상수 벡터에 다시 조립(`residualbased_block_builder_and_solver.h` `BuildMasterSlaveConstraints`). 증분형 scheme에서는 해가 증분이므로 c가 반복마다 더해진다. 선형 전략에서 관측된 "slave = Ru + 2c"가 그 결과. Newton 반복에서는 c가 반복 횟수만큼 누적될 수 있다.
- 트랙 A 규칙("고정 master는 관계식에 두고 c = 0")은 맞지만, 일반 규칙은 **"c는 항상 0, 오프셋은 고정된 master 자유도로 표현"**이어야 한다.
- 조치: `constraints.add_linear_tie`가 0이 아닌 `constant`를 받으면 오류를 내도록 보강. Phase 1 또는 Phase 5 시작 시 Newton 전략에서 구속식 결과 검증 시험 추가.

### F3 [중간] S15 대안의 한계 — 공간적으로 변하는 온도

- `Properties.SetValue(YOUNG_MODULUS)` 갱신은 Properties 단위(요소 묶음 전체)로만 바뀐다. glass-to-metal seal처럼 온도가 위치마다 다르면 적분점마다 다른 E가 필요하다.
- 트랙 A는 `SetTable`만 시험했고 Kratos Accessor(`Properties.SetAccessor` + 테이블/위치 기반 접근자) 경로는 시험하지 않았다.
- 조치: Phase 4 전에 Accessor 지원 여부를 확인하는 짧은 spike 추가. 안 되면 셀(요소)별 Properties 복제 또는 C++.

### F4 [낮음] 코드 품질

- `elements.attach_cell_elements`: 런타임 검사에 `assert` 사용 (`python -O`에서 사라짐) → `RuntimeError`로.
- `constraints._as_dof_lists`: 쓰지 않는 `KM` 인자, 의미 없는 간접 호출 → 제거.

### F5 [낮음] 시험 신뢰도

- S06 빔 이론 비교 허용 3배는 느슨함. 실질 근거는 "hex-Kratos 대비 1e-9 일치"이므로 spike로는 수용. 검증 단계(V09)에서는 수렴 연구로 대체.
- S16의 `dt < 60` 시간 단언은 기계 상태에 따라 흔들릴 수 있음(slow 마커로 분리되어 있음, 수용).

### F6 [정보]

- S03: 셀 면 위 점은 왼쪽 셀로 귀속 → Phase 1에서 K3(ε 이동) 정책 적용 필요 (이미 알려진 사항).
- S06: spectra는 shift = 0에서 실패, FEAST 미빌드(`HasFEAST = False`). 고유치는 shift 지정 또는 scipy.
- S13: 혼합 체적 요소는 적분점 기하에서 access violation → C++ MixedUP0 필요 (v2 예상대로).
- S16: 1e5 요소 생성 52.5초 — 60초 기준에 근접.

## G0/ADR-001 반영 사항

- ADR-001: R1 채택 제안 (S03·S04·S05 통과).
- G0 미결: F1(질량 — Python consistent M 경로로 C++ 없이 해결, 운용안 확정만 필요), F2(상수 안전장치 + Newton 검증), F3(Accessor spike).
- 트랙 B 결과 수령 후 통합.
