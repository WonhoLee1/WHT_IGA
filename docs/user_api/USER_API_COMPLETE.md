STATUS: COMPLETE_WITH_OPEN_ITEMS

# 사용자 Python API (D18) 완료 표식

- 작성: 2026-11-16T03:03+09:00 / 브랜치 `WHTIGA` (작업 트리, 미커밋)
- 계약: `WHT_IGA/research-cf-IGA/master-plan_v2.1.md` §8 (8.1–8.4), D18
- 사용 설명: `WHT_IGA/docs/user_api/USER_API.md`
- 기계 판독 상태: `WHT_IGA/docs/user_api/user_api_status.json`

## 범위 판정

| # | 요청 항목 | 판정 | 구현 |
|---|---|---|---|
| 1 | `wht_iga/spec/` Pydantic v2 frozen 모델, `with_updates()`, 단위계·물리 범위 검증, `model_json_schema()` | **완료** | `spec/base.py`, `units.py`, `mesh.py`, `assembly.py`, `physics.py`, `study.py` |
| 2 | `wht_iga/plan/` `compile_plan(spec) -> Plan(ResolveReport + ArtifactPaths)`, DOF·메모리·`auto` 솔버 전환 예측 | **완료** | `plan/compiler.py`, `plan/artifacts.py` |
| 3 | `wht_iga/ports/` `SolverBackendPort` + `backends/` `KratosPythonBackend`(escape hatch, `custom_processes`, 훅) + `ReferenceNumbaBackend` | **완료** | `ports/*.py`, `backends/*.py` (+ `MeshEnginePort` / `GridMeshEngine`) |
| 4 | `wht_iga/runner/` `AnalysisSession` (build_mesh, run, 6단 로그 On/Off, 진행 콜백, 취소 시 부분 산출물 정리) | **완료** | `runner/session.py`, `results.py`, `post.py` |
| 5 | `wht_iga/compat/coreform.py` (`ImmersedRectilinearMesh`, `BodyFitMesh`, `Refinement`, 영역 3종, 프로브 3종, `Job`) | **완료** | `compat/coreform.py` |
| 6 | `tests/unit/test_user_api.py`, `tests/test_isolation.py` 통과 | **완료** | 아래 실측 |

## 실측 (D:\mamba\envs\whtiga\python.exe, 2026-11-16)

| 명령 | 결과 |
|---|---|
| `python -m pytest -q -p no:cacheprovider tests/unit/test_user_api.py` | **28 passed** (12.1 s) |
| `python -m pytest -q -p no:cacheprovider tests/test_isolation.py tests/unit` | **542 passed** (36.4 s), 실패 0 |
| `python -m pytest -q -p no:cacheprovider tests/kratos_spikes/test_user_api_kratos.py` (소형 Kratos smoke) | **2 passed** (5.5 s) |
| 기존 `run_linear` Kratos 시험: `test_p1_loads.py`, `test_p1_outputs.py`, `test_p1_pipeline_cube.py`, `test_fix_converged.py` | 13 passed, **1 failed** (아래 미결 2, 기존 결함) |
| 변이 시험 5건 (압력 부호, 기울기 커널, 취소 정리, 체적력 계수, padding 검사) | **5/5 검출** (padding 건은 시험 강화 후 검출) |

해석해 대비 오차 (프로브, 최대 상대오차):

| 경우 | 백엔드 | `ux` | `σxx` |
|---|---|---:|---:|
| 균일 인장 막대, p=2, 8×2×2 (body-fit) | `ReferenceNumbaBackend` | 7.9e-15 | 1.1e-14 |
| 같은 경우 | `KratosPythonBackend` (penalty BC) | 4.2e-8 | 2.9e-11 |
| 절단 격자 자중 기둥, p=2, 10×3×2 (OUT 20셀, `out_mode="drop"`) | `KratosPythonBackend` | 8.8e-8 | 8.3e-8 |

## 수정한 기존 파일 (최소 변경)
- `wht_iga/__init__.py`: PEP 562 지연 export.
  - `import wht_iga`는 pydantic과 Kratos를 로드하지 않습니다(시험으로 확인).
- `wht_iga/solver/pipeline.py`: `run_linear(..., hooks=None)` 추가(+13줄).
  - `strategy.Solve()` 직전에 `hooks.before_solve`, 직후에 `hooks.after_solve`를 호출합니다.
  - 기본값 `None`이면 기존 동작과 같습니다.
- `pyproject.toml`: 의존성에 `pydantic>=2.7` 추가(plan §6).
- 환경: `uv pip install --python D:\mamba\envs\whtiga\python.exe "pydantic>=2.7"`.
  - 설치된 패키지: `pydantic 2.14.3`, `pydantic-core 2.50.4`, `annotated-types 0.8.0`, `typing-inspection 0.4.5`.

## 손대지 않은 것
- `wht_iga/geometry/*`, `wht_iga/solver/trimmed.py`: 공개 함수 호출만 했습니다.
  - 호출한 함수: `occ_io.load_cad_file`, `cad.shape_bbox`, `cad.tessellate`, `build_trimmed_table`.
- `E:\kratos-build`, `D:\PythonCodeStudy\kratos-src`.
- L2/L3 대형 벤치마크는 실행하지 않았습니다.

## 미결 항목
1. **`uv.lock` 미갱신.**
   - `pyproject.toml`에 pydantic을 추가했지만 `uv lock`은 실행하지 않았습니다. 잠금 파일 전체가 다시 풀릴 수 있기 때문입니다.
   - 갱신 여부는 사용자 결정 사항입니다.
2. **(해소됨) 기존 결함 `tests/kratos_spikes/test_fix_converged.py::test_nonconverged_raises_validation_error`.**
   - 원인: HEAD의 monkeypatch 람다가 `amgcl_params` 인자를 받지 않았습니다.
   - 다른 에이전트가 람다에 `amgcl_params=None`을 추가했습니다.
   - 2026-11-16T08:35 재실행: `test_user_api.py` + `test_isolation.py` + `test_fix_converged.py` → **35 passed**.
3. **설계상 제한.** plan이 error 진단을 내고, 침묵 대체는 없습니다.
   - THB 세분(`UNSUPPORTED_REFINEMENT`): spec은 정의되지만 `run_linear`에 연결되지 않았습니다.
   - 다중 파트.
   - CAD 면 세트: 세트는 격자 면만 지원합니다.
   - 점 프로브 `NT11`/`RF`.
   - Coreform `frame`/`element_layout`/프로브 `region`.
4. **plan §6에 있으나 이번 요청 범위 밖이라 만들지 않은 것.**
   - `ports/result_store.py`(`ResultStorePort`)
   - `runner/sweep.py`
   - `runner/external.py`
   - `KratosCxxBackend`
   - `kratos_overrides`(scheme·구성법칙 passthrough)
