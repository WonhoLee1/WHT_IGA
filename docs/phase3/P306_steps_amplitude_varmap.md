# P306 — 다중 스텝 실행기, 진폭 곡선, 변수 매핑 및 박판 자동 세분화

## 1. 개요 및 구현 범위

Phase 3의 시간/하중 이력 해석 오케스트레이션, Coreform/Abaqus ↔ Kratos 변수명 변환, 취소/진행률 훅, 그리고 얇은 기하(thin feature) 자동 세분화 루프를 구현한 `wht_iga/solver/steps.py` 기술 문서다(T-P3.6).

### 관련 소스 파일
- `WHT_IGA/wht_iga/solver/steps.py`
- 연동 안정화 모듈: `WHT_IGA/wht_iga/stabilize/extension.py` (`THIN_FEATURE`)

---

## 2. 변수 매핑 (`VAR_MAP`) 및 제어·진행률 훅

### 2.1 `ABAQUS_TO_KRATOS` / `KRATOS_TO_ABAQUS` / `VAR_MAP`
`resolve_var_name(name, target="kratos"|"abaqus")` 를 통해 입력·출력 계층 전역에서 단일 소스 양방향 매핑을 제공한다:

| Abaqus / Coreform 코드 | Kratos 변수명 | 물리량 |
|---|---|---|
| `U` | `DISPLACEMENT` | 변위 벡터 |
| `S` | `CAUCHY_STRESS_TENSOR` | 코시 응력 텐서 |
| `E` | `GREEN_LAGRANGE_STRAIN_TENSOR` | 변형률 텐서 |
| `NT11` | `TEMPERATURE` | 절점/제어점 온도 |
| `RF` | `REACTION` | 반력 벡터 |
| `CF` | `CONTACT_FORCE` | 접촉력 벡터 |

### 2.2 취소 토큰 및 진행 이벤트
- `CancelToken` / `CancelledError`: 장시간 다중 스텝 해석 루프 중 스레드 안전한 중단 요청(`cancel()`, `is_cancelled`, `raise_if_cancelled()`)을 처리한다.
- `ProgressEvent` / `emit_progress`: 스텝 인덱스, 스텝명, 현재 시간, 진행률, 메시지 페이로드를 등록된 콜백으로 전달한다.

---

## 3. 진폭 곡선 (`Amplitude`) 및 5대 해석 스텝 (`MultiStepRunner`)

### 3.1 `Amplitude` 시간 함수
정규화된 시간 $t \in [t_0, t_1]$ 또는 절대 시간 $t$ 에 대해 4종의 하중/경계조건 스케일링 함수를 평가한다:
- `tabular`: 구간별 선형 보간 테이블 `(t_k, a_k)`
- `ramp`: 선형 램프 보간
- `smooth_step`: 5차 $C^2$-연속 에르미트(Hermite) 다항식 $\xi^3(10 - 15\xi + 6\xi^2)$ ($\xi = \operatorname{clamp}((t - t_0)/(t_1 - t_0), 0, 1)$)
- `periodic`: 푸리에 급수(cosine/sine 계수 및 원형 주파수 $\omega$) 주기 진폭

### 3.2 스텝 명세 클래스 및 `MultiStepRunner`
| 스텝 클래스 | 대상 물리 및 수학적 적분기 |
|---|---|
| `StaticStep` | 준정적 선형/증분 평형 $K u = \alpha(t) f$ |
| `ModalStep` | 일반화 고유치 문제 $K \phi = \omega^2 M \phi$ (`n_modes`, 고유진동수 Hz 및 모드 형상) |
| `ImplicitDynamicStep` | Newmark-$\beta$ ($\beta, \gamma$) 또는 HHT-$\alpha$ 음해법 적분 $(K_{\text{eff}} \Delta u = r_{\text{eff}})$ |
| `ExplicitDynamicStep` | 중앙차분(central-difference) 양해법 적분 + 임계 시간 간격 $\Delta t_{\text{crit}}$ 가드 |
| `ThermalStep` | 열전도·열용량 방정식 $C_{\text{th}} \dot{T} + K_{\text{th}} T = q(t)$ (또는 정상상태 $K_{\text{th}} T = q$) |

- `MultiStepRunner`: 각 스텝 완료 후 `StepState` (시간, 누적 변위/속도/가속도/온도, 수렴 이력)를 다음 스텝으로 전달한다.
- 사용자 정의 `step_handler` 콜백을 연결하거나, 내장 연산자(`K`, `M`, `C`, `K_th`, `C_th`) 기반 솔버 경로를 직접 구동할 수 있다.

---

## 4. 박판/박육 형상 자동 세분화 (`THIN_FEATURE`)

### 4.1 `detect_thin_features` 및 `auto_refine_thin_features`
- `detect_thin_features` 는 `wht_iga.stabilize.extension` 의 연결성분 기반 스팬(span) 검사와 연계해 셀 두께 방향으로 최소 스팬이 부족한 영역을 탐지하고 `ThinFeatureReport` 및 `THIN_FEATURE` 이슈를 생성한다.
- `auto_refine_thin_features` 는 `THIN_FEATURE` 이슈가 감지되면 셀 분할수(`n_cells`)를 세분화하여 재검증하며, 무한 루프를 방지하기 위해 **최대 2회 재시도(max 2 retries)** 제한을 강제한다.

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_phase3_steps_probes.py` | `VAR_MAP`, `Amplitude` 4종, `MultiStepRunner` 5종 스텝 연쇄, `CancelToken`, `auto_refine_thin_features` (및 프로브/VTKHDF 포함) | **7 passed** (1.10s) | Kratos 비의존 단위 테스트 |
| `WHT_IGA/tests/unit/test_stabilize_extension.py` | `THIN_FEATURE` 탐지 및 비자명/연결성분 확장 연산자 | **15 passed** (1.07s) | 안정화·확장 회귀 검증 |
