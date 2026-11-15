# WHT_IGA 사용자 Python API (D18, plan v2.1 §8)

5계층 Hexagonal 구조입니다. 각 계층은 아래 계층만 사용합니다. `import wht_iga`는 가볍게 유지되며(PEP 562 지연 export), Kratos는 `KratosPythonBackend.solve_static` 안에서만 로드됩니다.

| 계층 | 패키지 | 핵심 | 비고 |
|---|---|---|---|
| 1 Spec | `wht_iga.spec` | `StudySpec`, `PartSpec`, `LinearElasticSpec`, `ImmersedGridSpec`, `BodyFitSpec`, `RefinementSpec`, `*RegionSpec`, `StaticStepSpec`, `DisplacementBCSpec`, `Pressure/Traction/GravityLoadSpec`, `Point/Line/MultiPointProbeSpec` | Pydantic v2 `frozen=True`, `extra="forbid"`, `with_updates()`, `model_json_schema()`, 단위계 타당성 검사 |
| 2 Plan | `wht_iga.plan` | `compile_plan(spec) -> Plan(spec, report: ResolveReport, artifacts: ArtifactPaths)` | 격자 결정, 세트→격자 면, DOF/nnz, RAM 스냅샷, `auto` 솔버 전환 예측, 미지원 조합 = error 진단 |
| 3 Ports/Adapters | `wht_iga.ports`, `wht_iga.backends` | `SolverBackendPort`, `MeshEnginePort` / `KratosPythonBackend`, `ReferenceNumbaBackend`, `GridMeshEngine` | 같은 `Plan`을 코드 수정 없이 백엔드 교체·교차검증 |
| 4 Runner | `wht_iga.runner` | `AnalysisSession` (with 문), `StudyResults` | 6단 솔버 로그 On/Off, 진행 콜백, 협조적 취소 + 부분 산출물 정리 |
| 5 Query | `wht_iga.query` | `ResultView` (기존) | `StudyResults.field()` / `.view()`가 VTKHDF(+SQLite)로 위임 |
| compat | `wht_iga.compat.coreform` | `ImmersedRectilinearMesh`, `BodyFitMesh`, `Element*Sizing`, `Refinement`, `CadEntitiesRegion`, `SphereRegion`, `BoxRegion`, `PointProbe`, `LineProbe`, `MultiPointProbe`, `Job` | Coreform 포터블 spec 이름·인자, `to_dict/to_json/from_dict/from_json/copy_with/validate`, `to_core()` |

## 사용 예

```python
import wht_iga as wi

study = wi.StudySpec(
    name="bar", units="mm-N-tonne-s",
    parts=[wi.PartSpec.from_box((0, 0, 0), (10, 1, 1), "steel")],   # 또는 PartSpec.from_step("plate.step", "steel")
    materials={"steel": wi.LinearElasticSpec(E=2.0e5, nu=0.3, density=7.85e-9)},
    mesh=wi.BodyFitSpec(degree=2, n_cells=(8, 2, 2)),                # 또는 ImmersedGridSpec(degree=2, element_size=0.5, padding=...)
    steps=[wi.StaticStepSpec(name="pull",
        bcs=[wi.DisplacementBCSpec(set="xmin", ux=0.0), wi.DisplacementBCSpec(set="ymin", uy=0.0),
             wi.DisplacementBCSpec(set="zmin", uz=0.0)],
        loads=[wi.TractionLoadSpec(set="xmax", vector=(2.0, 0, 0))])],
    probes=[wi.LineProbeSpec(name="axis", start=(0, .5, .5), end=(10, .5, .5), n_points=5, variables=("U", "S"))],
)
plan = wi.compile_plan(study, output_dir="runs/bar")
print(plan.report.summary()); plan.report.raise_if_invalid()

with wi.AnalysisSession(plan, backend=wi.KratosPythonBackend(threads=8), enable_solver_log=True) as session:
    session.build_mesh()
    results = session.run()
    model_part = session.backend.root_model_part      # escape hatch: Kratos.ModelPart

results.probe("axis")["S"]            # (5, 6) Voigt [xx, yy, zz, xy, yz, xz]
results.field("U")                    # VTKHDF 프록시 점 값 (ResultView)
fine = study.with_updates(mesh=study.mesh.with_updates(n_cells=(16, 4, 4)))
```

Coreform 스타일:

```python
import wht_iga.compat.coreform as cf
job = cf.Job("bar", parts=[...], materials={...},
             mesh=cf.ImmersedRectilinearMesh(degree=2, continuity=1, sizing=cf.ElementSizeSizing(element_size=0.5)),
             steps=[...], probes=[cf.PointProbe(name="tip", location=(10, .5, .5), variables=("U", "S"))])
job.write_input(); job.submit(); job.wait_for_completion(); job.status; job.results.probe("tip")
```

## 규약과 제약 (침묵 대체 없음: 모두 `compile_plan`의 error 진단이나 예외로 보고)

**세트와 경계·하중**
- 세트는 배경 격자 면입니다.
  - 기본 이름: `xmin xmax ymin ymax zmin zmax`.
  - 사용자 별칭: `StudySpec.sets = {"right": GridFaceSetSpec(axis=0, side=1)}`.
- BC·하중이 padding 안의 면에 있으면 `SET_OFF_DOMAIN` 오류입니다.
- 상자 파트의 면과 일치하지 않는 면도 `SET_OFF_DOMAIN` 오류입니다.
- 하중 면을 파트가 다 덮지 않으면 다음과 같이 처리합니다.
  - 상자 파트: `LOAD_FACE_NOT_COVERED` 오류.
  - CAD 파트: `build_mesh` 후 OUT 셀이 그 면에 닿으면 `LoadCoverageError`.

**하중 부호**
- 압력은 traction `+magnitude·n_out`입니다(저장소 S08 규약).

**단위계와 물성 범위**
- 단위계: `mm-N-tonne-s`, `m-N-kg-s`, `mm-kN-kg-ms`.
- E와 밀도를 SI로 환산해 물리 범위 밖이면 생성 시 오류입니다. SI 값을 mm 계에 넣는 실수 등을 잡습니다.

**미지원 조합 (error 진단)**
- 다중 파트(`UNSUPPORTED_MULTI_PART`).
- THB 세분(`UNSUPPORTED_REFINEMENT`): spec은 정의되지만 `run_linear`에 아직 연결되지 않았습니다.
- 점 값이 아닌 프로브 변수 `NT11`, `RF`(`PROBE_VARIABLE_UNSUPPORTED`).
- 비대칭 솔버가 필요한 follower 압력에 대칭 솔버 지정(`FOLLOWER_NEEDS_LU`).
- CAD 파트의 body-fit(`BODYFIT_REQUIRES_BOX`).

**메모리 예측**
- 절단 격자는 분류 전이라 전체 DOF를 상한으로 써서 보수적으로 예측합니다(`ACTIVE_DOFS_UPPER_BOUND`).

**백엔드**
- `ReferenceNumbaBackend`는 비절단 격자만 다룹니다. strong Dirichlet을 쓰며 follower 압력은 없습니다.
- 지원 범위를 벗어나면 세션 진입 시 `SessionError`입니다.

**프로브**
- U, E, S, MISES를 제어점 해에서 B-스플라인 값·기울기로 계산합니다. 두 백엔드가 같은 경로를 씁니다.
- S는 등방 선형탄성 `D·ε`입니다.

**취소**
- 단계 경계에서만 적용되는 협조적 취소입니다. 실행 중인 Kratos 해는 중단하지 않습니다.
- 취소 시 세션이 만들었거나 고친 산출물만 지우고, 기존 파일은 그대로 둡니다.

**Kratos 훅**
- `KratosPythonBackend(custom_processes=[f(model, parameters)], kratos_parameters={...}, on_before_solve=..., on_step_end=...)`.
- 내부적으로 `pipeline.run_linear(..., hooks=...)`의 `before_solve` / `after_solve`를 씁니다. 기본값 `None`이면 기존 동작과 같습니다.

**Coreform 퍼사드**
- `frame`(전역 축 외), `element_layout`, 프로브 `region`, Abaqus 의존부(`mdb`, `customData`)는 미지원이며 `ValidationError`가 납니다.
