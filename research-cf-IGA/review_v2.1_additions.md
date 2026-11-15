# v2.1 추가분 검토 (2026-11-14)

대상: 다른 에이전트가 `master-plan_v2.1.md`에 추가한 내용 — 모드 중첩 동역학, Explicit 동역학(Phase 6, S24, V-EXP), D17(SQLite + VTKHDF 결과 저장소), D18(Pydantic v2 + Hexagonal + CQRS API, Coreform 호환 퍼사드).

## 판정

| 항목 | 판정 | 요약 |
|---|---|---|
| 모드 중첩 과도·조화 응답 (§7.10, Phase 4) | **수락** | 고유모드 경로(Kratos 추출 K + Python consistent M)를 그대로 쓰므로 무리 없음 |
| Explicit 동역학 (§5.6, §7.9, S24, V-EXP) | **수정 필요** | Kratos 구속식 처리에 대한 핵심 전제가 사실과 다름(R1), lumping 설명 과장(R2), explicit 접촉 누락(R3) |
| D17 결과 저장소 | **수락** (보완 1건) | 대용량 시계열은 VTKHDF, 메타·Probe는 SQLite. 재시작용 제어점 계수 저장 위치 명시 필요(R7) |
| D18 API 아키텍처 | **조건부 수락** | 방향은 좋으나 문서 내 용어·디렉토리가 옛 설계와 섞여 충돌(R4, R5), 호환 퍼사드 "100%"는 과장(R6) |

## 수정 요청

### R1 [높음] Explicit + 구속식: "Kratos가 slave 내력을 master로 축약하지 않는다"는 틀림

- 문서(§5.6 Explicit 행, §7.9 항목 2)는 `explicit_builder.h`/`mechanical_explicit_strategy.hpp`가 `ApplyConstraints`로 slave 변위만 복원하고 Tᵀ f 축약을 하지 않는다고 적었다.
- 실제 소스(`StructuralMechanicsApplication/custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp` l.260, 320, 348–366):
  - 초기화 시 `ConstraintUtilities::PreComputeExplicitConstraintMassAndInertia(r_model_part)` — 구속식에 따라 질량·관성을 master로 모음
  - 매 단계 `ConstraintUtilities::PreComputeExplicitConstraintConstribution(r_model_part, dof_variable_names, residual_variable_names)` — slave의 잔차(FORCE_RESIDUAL 등)를 **`master += Tᵀ · slave`**로 master에 더함 (`kratos/utilities/constraint_utilities.cpp`의 `noalias(master_solution_vector) = prod(trans(transformation_matrix), slave_solution_vector)` 후 `AtomicAdd`)
  - 그 뒤 `ComputeExplicitConstraintConstribution` → `ResetSlaveDofs` → slave 값 복원
- 결론: **Kratos `MechanicalExplicitStrategy`가 extension 구속식(T)을 explicit에서 이미 처리한다.** Python explicit 루프는 기본 경로가 아니라 대안이다.
- 수정: §5.6 Explicit 행 대안 열과 §7.9 항목 2를 위 내용으로 고치고, S24 질문을 "Kratos `MechanicalExplicitStrategy` + MPC(T) + lumped M으로 1D 파동 해석해를 재현하는가"로 바꾼다. Python 루프는 실패 시 대안.

### R2 [중간] Lumping 설명 정정

- 문서는 row-sum lumping이 "고유진동수 왜곡·고차 모드 에너지 비보존"을 일으키므로 HRZ 또는 "Optimal Scaled Mass Lumping (Cottrell 2006)"을 쓴다고 적었다.
- B-스플라인은 음수가 없어 row-sum lumping이 **항상 양의 대각 질량**을 주며, IGA explicit에서 표준적으로 쓰인다. 정확도는 2차로 떨어지지만 "에너지 비보존"은 아니다(중앙차분은 lumped M에 대해 에너지 일관). Cottrell 2006은 row-sum lumping의 정확도를 분석한 논문이지 "optimal scaled" 기법의 출처가 아니다.
- 실측도 이미 있다: Kratos lumped(적분점 기하에서 정상 동작) vs consistent — V09 보 1차 모드 0.07%(p=2, 40×4×4), 상위 20모드 최대 7.9%→2.3%.
- 수정: 기본 = Kratos lumped(row-sum 계열, 실측 근거), HRZ는 비교 옵션으로 S24에서 측정. 출처 표기 정정.

### R3 [중간] Explicit 접촉 누락

- 사용자가 기대한 범위(implicit dynamic, 비선형 접촉, 비선형 정적, explicit)에 explicit 접촉이 자연스럽게 포함된다.
- Kratos에 `ContactStructuralMechanicsApplication/python_scripts/explicit_penalty_contact_process.py`와 `contact_structural_mechanics_explicit_dynamic_solver.py`가 있다.
- 수정: Phase 6 Explicit 트랙에 "explicit penalty 접촉(표면층, §5.9.4와 같은 결합)" 추가, spike S25(explicit 접촉 + 스프링 결합 표면층) 추가, 검증은 두 블록 충돌(운동량·에너지 보존) 1건.

### R4 [중간] 용어·이름 충돌 (D18 반영 누락)

- §3.4, §1.3, §5.9.2, §5.9.4, §5.9.7, R17, Phase 1이 아직 옛 API(`job.write_input/build_mesh/resolve/submit`, `api/job.py`, `api/spec.py`, `api/resolve.py`, `api/assembly.py`)를 쓴다. §6·§8은 새 구조(`spec/`, `plan/compiler.py`의 `compile_plan`, `runner/session.py`의 `AnalysisSession`)다.
- 수정: 대응표 한 줄로 통일 — `job.write_input` → `StudySpec.model_dump_json`, `job.build_mesh` → `AnalysisSession.build_mesh`, `job.resolve` → `compile_plan`(→ `ResolveReport`), `job.submit` → `AnalysisSession.run`, `api/*` → `spec/*`, `plan/*`. W0–W5 표의 모듈 열 갱신.

### R5 [중간] 실제 충돌: `wht_iga/runner.py`와 `wht_iga/runner/` 패키지

- §6 트리에 `runner/`(Layer 4 패키지)와 `runner.py`(외부 프로세스 실행기, §3.3)가 같은 수준에 함께 있다. Python에서 같은 이름의 모듈과 패키지는 공존할 수 없다.
- 수정: 외부 실행기를 `wht_iga/runner/external.py`로 옮기고 §3.3 참조를 고친다.
- 같은 맥락: §8.1이 `backends`를 언급하지만 트리에 없다. `KratosPythonBackend`는 `kratos_backend/`, `ReferenceNumbaBackend`는 `reference/`에 둔다고 명시.

### R6 [낮음] 호환 퍼사드 "100% 수용"은 과장

- Coreform 스크립트는 Abaqus/CAE 안에서 `mdb`, `part.customData.CoreformIGA`, `coreform.abaqus.init_iga_mdb()` 등을 쓴다. Abaqus 객체 모델 없이 100% 수용은 불가능하다.
- 수정: "Coreform **포터블 spec 계층**(`coreform.mesh` / `coreform.probe` / `coreform.job`의 spec, `to_json/from_json/validate`)과 같은 이름·인자를 받는다. Abaqus 의존 부분(`mdb`, `customData`, CAE 기록)은 지원하지 않는다."

### R7 [낮음] D17 보완

- VTKHDF는 proxy 격자 위 필드를 담는다. 재시작·후처리 재평가에 필요한 **제어점 계수(U, 상태 변수)**를 어디 둘지 명시가 없다.
- 수정: VTKHDF 안에 별도 그룹(`/WHT_IGA/control_points/<frame>`)으로 제어점 계수를 저장하거나, 소규모는 SQLite `cp_field` 유지 — 둘 중 하나로 확정.

### R8 [낮음] 남은 옛 문장

- §5.6 l.515 "고유치 | SM `EigensolverStrategy` + spectra/FEAST" 행이 바로 아래 행(Kratos 추출 K + Python consistent M + scipy)과 모순된다 → 삭제. (v2 시절 행이 남은 것)
- §6 `tests/kratos_spikes/` 주석 "S01–S18" → "S01–S25".
- 모드 중첩 과도: 모드 절단 오차 보정(static correction / residual mode)을 선택 옵션으로 명시 권장(V10 직접법 비교 시 차이 원인).

## 우선순위

1. R1 (경로 판단에 영향) → 2. R5 (구현 시 바로 오류) → 3. R4 → 4. R2, R3 → 5. R6–R8

## 적용 결과 (2026-11-14)

R1–R8 모두 `master-plan_v2.1.md`에 반영함 (메인 세션). R1: §5.6 Explicit 행·§7.9·S24를 Kratos `MechanicalExplicitStrategy` 구속 처리 기준으로 수정. R2: lumping 설명·출처 정정. R3: explicit 접촉(§7.9 항목 3, S25, Phase 6, V-EXP-C) 추가. R4: `job.*`/`api/*` → D18 이름으로 통일 + §3.4 대응 문단. R5: `runner/external.py`로 이동. R6: 호환 퍼사드 범위 정정. R7: VTKHDF 제어점 계수 그룹 명시. R8: 옛 고유치 행 삭제, spike 주석, 모드 절단 보정 옵션 추가. 편집 전 사본: scratchpad `v2.1_before_review_fix.md`.
