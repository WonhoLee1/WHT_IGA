# WHT_IGA 사용자 매뉴얼

이 문서는 WHT_IGA를 쓰는 사용자를 위한 모듈별 사용법이다. 예제 코드는 저장소의 시험(`tests/`)이 실제로 실행하는 경로를 따랐다.
내부 구현은 `docs/verification/`와 `docs/spikes/`의 기록을 본다.

## 0. 환경

- Python: `D:\mamba\envs\whtiga\python.exe` (Python 3.12). 모든 명령은 저장소 루트 `WHT_IGA`에서 실행한다.
- Kratos: `E:\kratos-build\install` 빌드만 쓴다. `from wht_iga.kratos_backend import load_kratos` 로만 불러온다.
- 패키지 설치: `uv pip install --python D:\mamba\envs\whtiga\python.exe <패키지>`. conda 명령은 쓰지 않는다.
- 메모리: 대형 해석은 한 번에 하나씩 돌린다. 이 PC는 32 GB, 8코어다.

## 1. 입력 스펙 (`wht_iga.spec`)

해석 입력은 pydantic 스펙 객체로 만든다. 주요 객체:

| 객체 | 용도 |
|---|---|
| `StudySpec(name, units, parts, materials, mesh, steps, probes)` | 전체 해석 정의 |
| `PartSpec.from_box(lo, hi, material)` | 상자형 파트 |
| `LinearElasticSpec(E, nu)` | 선형 탄성 재료 |
| `BodyFitSpec(degree, continuity, n_cells)` | 몸체 적합 격자 |
| `ImmersedGridSpec(...)` | 몰입형 격자(배경 격자와 절단) |
| `StaticStepSpec(name, bcs, loads, locking_mitigation)` | 정적 단계 |
| `DisplacementBCSpec(set, ux, uy, uz)`, `TractionLoadSpec(set, vector)`, `PressureLoadSpec`, `GravityLoadSpec` | 경계와 하중 |
| `PointProbeSpec`, `LineProbeSpec`, `MultiPointProbeSpec` | 결과 탐침 |

`locking_mitigation`은 `"none"`(기본, 선형 탄성), `"fbar"`(F̄/P0 강제), `"auto"`(ν ≥ 0.49면 F̄)다.

## 2. 실행 (`wht_iga.plan`, `wht_iga.runner`, `wht_iga.backends`)

```python
import wht_iga.spec as s
from wht_iga.backends import ReferenceNumbaBackend   # Kratos 없이 Numba로 실행
from wht_iga.plan import compile_plan
from wht_iga.runner import AnalysisSession

plan = compile_plan(study, output_dir="out/case")    # 스펙 검증과 산출물 경로 결정
with AnalysisSession(plan, backend=ReferenceNumbaBackend()) as ses:
    res = ses.run()
step = res.step("bend")
step.U                # 제어점 변위
step.probes           # 탐침 결과 (U, S, E, MISES)
step.converged, step.n_iter, step.stress_measure
```

- `backend`는 두 가지다. `ReferenceNumbaBackend()`는 Kratos가 필요 없다. `KratosPythonBackend()`는 Kratos 요소와 풀이기를 쓴다.
- `res.step(name)`에서 단계 결과를 꺼낸다. `stress_measure`는 탐침 응력의 종류다(`"linear_elastic"` 또는 `"fbar_neo_hookean_cauchy"`).
- 진행 상황은 `AnalysisSession(..., on_progress=callback)`으로 받는다.

## 3. 방법 선택과 사용자 API (`wht_iga.analyses.methods`)

격자 방법(몸체 적합, 몰입형, F̄ 등)을 고르고 싶을 때 쓴다. 검증된 호출은 다음과 같다.

```python
from wht_iga.analyses.methods import (
    BoxMethodProblem, parse_discretization_option, solve_with_discretization_method,
    recommend_discretization,
)
problem = BoxMethodProblem(lo=(0, 0, 0), hi=(4, 1, 1), E=240.565, nu=0.49999,
                           fixed_face=(0, 0), load_face=(0, 1), traction=(0.0, 0.5, 0.0))
option = parse_discretization_option({"kind": "full_iga_body_fit", "degree": 2,
                                      "n_cells": [3, 2, 2], "locking_mitigation": "fbar"})
sol = solve_with_discretization_method(problem, option, use_pardiso=False)
sol.tip_displacement
rec = recommend_discretization("bending", nu=0.49999, n_cells=[3, 2, 2])  # 추천 옵션 패치
```

- `kind`: `"full_iga_body_fit"`(몸체 적합), `"full_iga_immersed"`(몰입형). 몰입형 F̄는 `out_mode="drop"`, 균일 C^(p-1), `bc_penalty` 없음일 때만 된다. 그 밖의 조합은 `ValueError`로 거부된다.
- `recommend_discretization`은 ν와 하중 종류로 차수와 잠김 완화를 제안한다. 규칙 기반이므로 신뢰도는 `heuristic`이다.

## 4. 고급 해석 모듈

| 모듈 | 용도 | 주요 함수/클래스 |
|---|---|---|
| `wht_iga.analyses.fbar` | 비압축 탄성의 F̄/P0 (Neo-Hookean, 몸체 적합) | `build_model`, `solve_static`, `residual_tangent_sparse`, `cauchy_stress` |
| `wht_iga.analyses.fbar_trimmed` | 절단 격자 F̄/P0, 작은 절단 셀 압력 병합 | `build_model(box, qp_table, out_cells, E, nu, ...)`, `pressure_groups`, `residual_tangent_sparse`, `group_jbar` |
| `wht_iga.analyses.hyperelastic_fbar` | 다항식 초탄성(Mooney-Rivlin, Yeoh, N≤3)의 F̄/P0, Numba | `mooney_rivlin`, `neo_hookean`, `yeoh`, `polynomial`, `build_model(box, degree, material, ...)` |
| `wht_iga.analyses.mixed_p1` | 셀별 P1/P2 불연속 압력 실험(V07 기울기 연구) | `build_model(box, degree, E, nu, pressure="P1")`, `residual_tangent`, `pressure_projection` |
| `wht_iga.analyses.contact_nlgeom` | NLGEOM 몸체 + SSN 접촉(Kratos TL 강성 재조립) | `KratosBodyOperator`, `solve_nlgeom_contact_step` |
| `wht_iga.analyses.contact` | 변형체 SSN 접촉(표면 가우스 점 penalty, 마찰) | `solve_deformable_ssn_contact` |
| `wht_iga.analyses.nonlinear` | Kratos 비선형 하중 단계와 Newton | `run_load_stepping`, `make_newton_strategy`, `j2_law`, `hyperelastic_neohookean_law` |
| `wht_iga.analyses.thermal` | 열 전도, 열변형률 | (모듈 내 함수 참고) |
| `wht_iga.analyses.explicit` | 명시적 동역학(중앙차분, 집중 질량) | (모듈 내 함수 참고) |

F̄ 해의 응력 탐침(`S`, `MISES`)은 Neo-Hookean F̄ Cauchy 응력으로 계산된다. 선형 해석의 탐침은 `D:ε`를 쓴다. 두 경우 모두 `StepResult.stress_measure`에 표시된다.

## 5. 기하와 격자 (`wht_iga.geometry`, `wht_iga.grid`, `wht_iga.solver`)

- CAD 읽기: `wht_iga.geometry.occ_io.load_cad_file(path, heal=False)` (STEP, SAT 지원). 테셀레이션: `wht_iga.geometry.cad.tessellate(shape, deflection=...)`.
- B-spline 격자: `wht_iga.reference.box.BoxSpline(lo, hi, n_cells, degree=p)`. 매듭 편집은 `knots[d]`와 `n_cp`를 맞춰서 한다.
- 절단 적분 표: `wht_iga.solver.trimmed.build_trimmed_table(nodes, tris, box, q=..., cache_dir="off")`. 반환은 `(qp_table, out_cells, issues, cell_records, diag)`다. 절단 셀 적분 규칙은 `wht_iga.geometry.quadrature`에 있다.
- 작은 절단 셀 안정화: `wht_iga.stabilize.support.classify_basis`, `wht_iga.stabilize.extension.build_extension`.
- 선형 풀이 파이프라인: `wht_iga.solver.pipeline.run_linear(JobSpec(...))`.

## 6. 결과 후처리 (`wht_iga.post`, `wht_iga.runner.post`)

- `evaluate_point_fields(box, U, points, E, nu, variables, fbar_model=None)`: 점에서 U, E, S, MISES를 계산한다. `fbar_model`을 주면 F̄ 응력을 쓴다.
- `write_step_vtkhdf(path, box, U, out_cells, E, nu, fbar_model=None)`: 결과를 VTK HDF로 쓴다(셀 대리 메쉬, 절단 셀은 전체 육면체로 근사).
- 필드 그림: 저장된 `.vtu`는 `pyvista`로 읽어 x, y, z, iso 뷰로 렌더링할 수 있다(오프스크린 가능).

## 7. 동역학 (`wht_iga.dynamics`)

- 진동 응답 스펙트럼: `response_spectrum.py`(Abaqus 주파수 스윕, 정의는 모듈 문서 참고).
- 랜덤 진동: `random_vibration.py`(PSD 테이블, 옥타브 대역, 모드 중첩 TENP/GRP/DSC). 정의는 Coreform 원문과 일치하지 않는 가정이 있으므로 V09 문서를 먼저 읽는다.
- 시간 적분: Kratos 경로는 Bossak(`ResidualBasedBossakDisplacementScheme`, 기본 `damp_factor_m = -0.3`, `newmark_beta = 0.25`)을 쓴다. Python 경로(`solve_bossak_linear`)는 집중·일관 질량 모두 지원한다. HHT는 구현되어 있지 않다(V10 문서).

## 8. 검증 실행

```powershell
D:\mamba\envs\whtiga\python.exe -m pytest -q -p no:cacheprovider tests\test_isolation.py tests\unit   # 단위(빠름)
D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests\kratos_spikes\test_p5_v08_ghost_check.py  # V08 ghost (약 5분)
D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests\kratos_spikes\test_p5_v08_r2_h2_nlgeom.py  # V08 R2 h=0.125 (약 1.9 h)
```

- `tests/test_isolation.py`는 `wht_topopt`, `aieng`, `gmsh`를 `wht_iga`에서 import하지 않는지 검사한다.
- 검증 명세는 `tests/verification/data/*.json`에 있다. 게이트는 명세의 `acceptance`에 적혀 있다.
- Coreform 기준 CSV는 `research-cf-IGA/user-guide/generated-downloads/`에 있다.

## 9. aieng-ui 연계 (MCP 도구)

WHT_IGA의 해석을 aieng 작업 공간에서 호출하는 도구다. aieng 환경은 `wht_iga`를 직접 import하지 않고, 프로세스 경계를 거쳐 JSON으로 주고받는다.

| 도구 | 성격 | 설명 |
|---|---|---|
| `iga.runner_status` | 읽기 | whtiga 실행 파일과 상태 |
| `iga.recommend_discretization` | 읽기 | 방법 추천 |
| `iga.run_study` | 승인 필요 | 스펙 해석 실행, 결과는 `<data_root>/iga_runs/<run_id>` |
| `iga.run_method` | 승인 필요 | 방법 선택 해석 실행 |

환경 변수: `AIENG_WHTIGA_PYTHON`(whtiga 실행 파일), `AIENG_WHT_IGA_ROOT`(WHT_IGA 경로), `AIENG_WHTIGA_TIMEOUT_S`(시간 제한).
MCP 서버를 재시작해야 `iga.*` 도구가 보인다.

## 10. 알려진 제약과 주의

- 응력 결과의 종류를 항상 확인한다(`stress_measure`). 선형 해의 응력을 F̄ 해에 쓰면 틀린다.
- 몰입형 F̄는 지원 범위가 좁다(4절 참고). 범위 밖이면 조용히 다른 경로로 바꾸지 않고 오류를 낸다.
- 메쉬 격자 스펙(요소 크기, 패딩)이 결과를 크게 좌우한다. 검증 비교 시 명세의 격자 스펙을 그대로 따른다(V04 사례).
- 대형 해석(수십만 자유도 이상)은 메모리를 먼저 확인한다. 측정 방법은 `docs/phase6/S50_bracket_1m_scaling.md`.
- 커밋과 결과 파일: 100 MB를 넘는 결과 파일은 저장소에 넣지 않는다.
