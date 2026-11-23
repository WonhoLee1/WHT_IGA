# WHT_IGA

Kratos Multiphysics 기반 **몰입형(immersed) 등기하 해석(IGA)** 프레임워크.
목표는 Coreform IGA 검증 매뉴얼과 예제의 결과를 같은 조건에서 재현하고, 그 일치 정도를 게이트로 판정하는 것이다.
Kratos는 요소·풀이기 등 일부 경로에만 쓰고, 공식(B-spline 기저, 절단 적분, F̄/P0, SSN 접촉 등)은 Python(NumPy, Numba)으로 구현했다.

- 작업 규칙(에이전트용): [`AGENTS.md`](AGENTS.md)
- 환경과 Kratos 연결: [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md), [`docs/KRATOS_BUILD.md`](docs/KRATOS_BUILD.md)
- 사용자 매뉴얼(모듈별 사용법): [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md)
- 마스터 계획: [`research-cf-IGA/master-plan_v2.1.md`](research-cf-IGA/master-plan_v2.1.md)
- 검증 기록: [`docs/verification/`](docs/verification/), 스파이크 기록: [`docs/spikes/`](docs/spikes/)

## 1. 현재 상태 (2026-11-23 기준)

Coreform 검증 매뉴얼 14건과 예제 5건을 19개 케이스(V01–V12, E01–E05)로 대응시켜 판정했다.
판정 기준은 각 케이스의 기록된 게이트이며, 게이트를 조정하지 않았다.

| 판정 | 건수 | 케이스 |
|---|---:|---|
| 충족 | 15 | E01, E02, E03, E05, V01, V02, V03, V04, V05, V08, V10(HHT), V11, V12 T2, V12 T3(방법 재현), V12 T4 |
| 충족(변형) | 1 | V07 (p=3 C⁰/P2 변형, σ 기울기 2.58 ≥ 2.2; Coreform과 같은 p=2 P0 경로는 0.85) |
| 충족(단서) | 1 | V09 (p=1 게이트 +0.33 %; RMS 정의는 가정) |
| 조건부 | 2 | E04(ν = 0.495로 전환: C⁰·C¹ 극한 0.24 % 일치, C⁰ 64×64 +2.4 %), V06(p=1 L0 셀: Coreform 절단 셀 적분이 p=1 강성에 정확하지 않음) |
| 불일치 | 0 | — |

정정·해소 이력은 `docs/HANDOFF_20261124.md`. 막힌 항목(정보 필요): V06 Coreform 절단 셀 적분 규칙, E04 Coreform 격자, V09 랜덤 응답 정의, V10 Abaqus 이력(변위·위상 직접 비교).
상세 내용과 그림은 보고서 "WHT_IGA 검증 보고서"를 본다.

## 2. 구성 (모듈 맵)

| 패키지 | 역할 |
|---|---|
| `wht_iga/spec` | 사용자 입력 스펙(`StudySpec`, `PartSpec`, `StaticStepSpec` 등, pydantic) |
| `wht_iga/plan` | 스펙을 실행 계획으로 컴파일(`compile_plan`), 산출물 경로 결정 |
| `wht_iga/runner` | 실행 세션(`AnalysisSession`), 결과 객체(`StudyResults`, `StepResult`), 프로브 후처리 |
| `wht_iga/backends` | 실행 백엔드: `ReferenceNumbaBackend`(Kratos 없이 Numba), `KratosPythonBackend`(Kratos 경로) |
| `wht_iga/analyses` | 해석 모듈: F̄/P0(`fbar`, `fbar_trimmed`), 초탄성(`hyperelastic_fbar`), P1/P2 압력(`mixed_p1`), 접촉(`contact`, `contact_surface`), NLGEOM 접촉(`contact_nlgeom`), 방법 선택·사용자 API(`methods*`), 비선형 하중 단계(`nonlinear`), 열(`thermal`, `thermal_plasticity`), 명시적 동역학(`explicit`) |
| `wht_iga/solver` | 선형 풀이 파이프라인(`pipeline.run_linear`), 절단 셀 표(`trimmed`), 벌칙 경계(`penalty`), 하중(`loads`) |
| `wht_iga/grid` | B-spline 기저, 매듭(`knots`), 계층 격자(`hierarchy`) |
| `wht_iga/geometry` | CAD 입출력(STEP, SAT), 테셀레이션, 절단 분류(`classify`), 적분 규칙(`quadrature`) |
| `wht_iga/kratos_backend` | Kratos 모델 구성(요소, 경계, 표면 하중, tie, 블록 조립)과 `load_kratos()` |
| `wht_iga/stabilize` | 작은 절단 셀 안정화(지지 분류, extension 연산자, 조건수) |
| `wht_iga/reference` | 해석해와 참조 조립기(Hertz, Cook, Roark, 상자 B-spline 기준 해) |
| `wht_iga/dynamics` | 응답 스펙트럼, 랜덤 진동, 모달 동역학 |
| `wht_iga/post` | VTK/VTKHDF 출력, 프로브 필드(U, S, E, MISES) |
| `wht_iga/verification` | 검증 명세 로딩, 보고서 생성(V/E 케이스 명세는 `tests/verification/data/*.json`) |
| `wht_iga/meshio_ext` | 메쉬 파일 읽기(INP 등) |
| `wht_iga/resources` | 예제 CAD 파일(STEP/SAT: Door_Lock_Assembly.sat 등) |

aieng-ui 연계(외부 저장소 `aieng-ui/backend/app`): `iga_bridge.py`(프로세스 경계 클라이언트), `iga_runner.py`(whtiga 파이썬에서 실행), `runtime_registry/iga.py`(MCP 도구 `iga.*`).
aieng311 환경은 `wht_iga`를 직접 import하지 않는다.

## 3. 빠른 시작

```powershell
# 환경: micromamba 환경 whtiga (Python 3.12). conda 명령 사용 금지.
D:\mamba\envs\whtiga\python.exe -m pytest -q -p no:cacheprovider tests\test_isolation.py tests\unit   # 단위 시험
D:\mamba\envs\whtiga\python.exe -m pytest -q -p no:cacheprovider tests\kratos_spikes\test_p6_cross_validation_10_4.py  # Kratos 교차검증
```

Kratos는 `wht_iga.kratos_backend.load_kratos()`로만 불러오며, `E:\kratos-build\install`의 빌드만 허용한다.
Python 패키지 설치는 `uv pip install --python D:\mamba\envs\whtiga\python.exe ...`만 쓴다.

## 4. 사용 예 (정적 해석, 검증된 코드 경로)

```python
import wht_iga.spec as s
from wht_iga.backends import ReferenceNumbaBackend
from wht_iga.plan import compile_plan
from wht_iga.runner import AnalysisSession

study = s.StudySpec(
    name="bend", units="mm-N-tonne-s",
    parts=[s.PartSpec.from_box((0, 0, 0), (4, 1, 1), "rubber")],
    materials={"rubber": s.LinearElasticSpec(E=240.565, nu=0.49999)},
    mesh=s.BodyFitSpec(degree=2, continuity=1, n_cells=(3, 2, 2)),
    steps=[s.StaticStepSpec(
        name="bend",
        bcs=[s.DisplacementBCSpec(set="xmin", ux=0.0, uy=0.0, uz=0.0)],
        loads=[s.TractionLoadSpec(set="xmax", vector=(0.0, 0.5, 0.0))],
        locking_mitigation="fbar",   # "none" | "fbar" | "auto"(nu >= 0.49)
    )],
)
plan = compile_plan(study, output_dir="out/bend")
with AnalysisSession(plan, backend=ReferenceNumbaBackend()) as ses:
    res = ses.run()
U = res.step("bend").U     # 제어점 변위, (n_dof,) xyz 인터리브
```

방법 선택(직접 호출): `recommend_discretization`, `solve_with_discretization_method`(`wht_iga.analyses.methods`).
자세한 사용법과 주의점은 [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md).

## 5. 알려진 한계

- 절단(immersed) F̄는 `out_mode="drop"`, 균일 C^(p-1), `bc_penalty` 없음 조건만 지원한다.
- Ogden 초탄성은 미구현. 다항식 초탄성은 Numba 평균 팽창 경로만 있다(Kratos C++ 아님).
- 응력 탐침: F̄ 해의 S, MISES는 Neo-Hookean F̄ Cauchy 응력이다. 선형 해에는 D:ε를 쓴다.
- 1M DOF 실행은 64 GB 장비 필요(이 PC 32 GB에서는 27.8만 DOF까지 실측, `docs/phase6/S50_bracket_1m_scaling.md`).
- 명시적 동역학 정책과 HHT 시간 적분은 미결(V10, `docs/verification/V10.md`). Kratos 구조 동역학의 기본 적분기는 Bossak이다.

## 6. 기여 규칙

- 게이트와 허용치는 측정값에 맞춰 바꾸지 않는다. 실패는 기록하고 원인을 분해한다.
- 저장소에 결과 파일(evidence, 100 MB 이상)을 커밋하지 않는다. `.omc/`, `.cache/`는 제외한다.
- 커밋과 푸시는 사용자 요청이 있을 때만 한다.
