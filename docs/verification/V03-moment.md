# V03-moment — 외팔보 굽힘 모멘트 검증 (곡선 보 및 직선 보 RBE3 커플링)

> **정정 (2026-11-23):** `scripts/v03_beam.py`가 E = 2.1e11과 끝단 하중 5e4 Pa(P = 500 N)를 썼다. 스펙(`v03.json`)은 E = 2.0e11, 끝단 하중 500 Pa(P = 5 N)이다. 앞서 기록한 "Coreform 하중 1.05배 불일치"는 이 입력 오류 때문이었다. 스펙 값으로 고쳐 다시 돌렸다(`docs/phase1/evidence/v03/results.json`).
> 결과(스펙 입력): 5개 하중 케이스 모두 EB 대비 0.80 %(분포 하중) 또는 0.60 %(끝단 전단)로 게이트 5 % 안. 끝단 전단의 격자 스윕(`(0.1,0.1,0.25)/2^k`, k=4 = 16×16×64, p=2) 값 `1.005958e-6`은 Coreform Abaqus 최세 `1.00595798e-6`과 일치한다. 스크립트 격자(4×4×40)의 값은 `1.005952e-6`이다.

- 상태: **PASS** (Gate G3 기준 충족)
- 게이트 (G3):
  - 곡선 외팔보(NURBS body-fit + `RBE3` 분포 모멘트 커플링): Timoshenko–Goodier 해석해 대비 끝단 반경방향 변위 상대오차 **≤ 5%**
  - 직선 외팔보([v03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v03.json) 모멘트 하중): Euler–Bernoulli 해석해 대비 상대오차 **≤ 5%** (게이트 G3 강화 기준 ≤ 0.5%)
- 날짜 / 환경: 2026-11-15, WHTIGA 브랜치, Python 3.12.16, Kratos 10.4.3 (`PARDISO`)
- 참조 명세: [v03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v03.json)
- 실측 데이터: [g3_measurements.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/g3_measurements.json)
- 재현 명령:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py -k test_g3_01_v03_moment_curved_and_straight_cantilever -v
  ```

## 1. 입력 요약

### 1.1 곡선 외팔보 (Curved Cantilever, Body-Fit NURBS + RBE3)
- 기하 형상: 1/4 원호 평면 곡선 보, 내경 $R_i = 4.12$, 외경 $R_o = 4.32$ (두께 $h = 0.2$), 폭 $b = 0.1$, 각도 $\theta \in [0, \pi/2]$
- 재료 물성: 등방 선형탄성 $E = 1.0 \times 10^7$, $\nu = 0.25$
- 이산화: 구조화 `hex8` 격자 `[4, 32, 2]` 를 [`mesh_to_bodyfit_patches`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/meshio_ext/bodyfit.py) (`degree=2`, `continuity="smooth"`)로 2차 NURBS 볼륨 변환 후 Kratos `NurbsVolumeGeometry` + `SmallDisplacementElement3D4N` Gauss 구적점 조립
- 경계조건 및 하중:
  - 고정단($\theta = 0$): 전 변위 성분 구속 ($u_x = u_y = u_z = 0$)
  - 자유단($\theta = \pi/2$): 중심점 $(0,\, (R_i+R_o)/2,\, 0)$ 기준 절점에 [`RBE3Spec`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/interactions.py) (`mode="mpc"`) 분포 모멘트 커플링 적용, 순수 굽힘 모멘트 $M_z = 100.0$ 부가

### 1.2 직선 외팔보 (`v03.json` 모멘트 케이스)
- 기하 형상: `[-0.05, -0.05, 0.0]` ~ `[0.05, 0.05, 1.0]` ($b = 0.1\text{ m}$, $h = 0.1\text{ m}$, $L = 1.0\text{ m}$, $I_x = 8.333333333333334 \times 10^{-6}\text{ m}^4$)
- 재료 물성: $E = 2.0 \times 10^{11}\text{ Pa}$, $\nu = 0$ (Euler–Bernoulli 비교 정합용)
- 이산화 및 해석: `BoxSpline` `(2, 2, 16)`, `degree=2`, [`AssemblyBuilder`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/assembly_builder.py) + [`solve_hybrid_iga_fe`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/native_fe.py) (`use_pardiso=True`)
- 경계조건 및 하중: $z = 0$ 전면 고정, $z = 1.0\text{ m}$ 자유단 중심점 `(0.0, 0.0, 1.0)`에 `RBE3Spec`으로 $x$축 굽힘 모멘트 $M_0 = 5.0\text{ N}\cdot\text{m}$ 부가

## 2. 실측 결과 및 기준 비교

| 케이스 | Kratos 실측값 | 해석해 / 기준값 | 상대오차 | 게이트 허용치 | 판정 |
|---|---:|---:|---:|---:|---|
| 곡선 보 끝단 $u_r$ (`[4, 32, 2]`) | `-2.6650569048944446` | `-2.670772649284193` (Timoshenko–Goodier) | `0.002140` (**`0.2140%`**) | `≤ 5.0%` | **PASS** |
| 곡선 보 끝단 $u_x$ (참고 기록) | `-1.5198803768216522` | — | — | 참고치 | — |
| 직선 보 끝단 $u_y$ ([v03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v03.json)) | `-1.500374204689463e-06` m | `-1.4999999999999998e-06` m (EB) | `0.000249` (**`0.0249%`**) | `≤ 0.5%` (문서 `≤ 5.0%`) | **PASS** |

### 기준값 정합성 주석 (`prompt_ref_note`)
- 프롬프트에서 인용된 `u_r_ref = -0.0211622` 수치는 본 문제의 입력 파라미터($M = 100$, $E = 10^7$, $R_i = 4.12$, $R_o = 4.32$, $b = 0.1$)와 물리적·수치적으로 정합하지 않습니다.
- 동일 파라미터에 대한 엄밀한 Timoshenko–Goodier Airy 응력함수 해석해([`test_p3_gate_g3_benchmarks.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L120-L195) 내 `_curved_beam_timoshenko_ur`)는 **`-2.670772649284193`**(주석 기록상 `-2.67146011`)이며, Kratos body-fit NURBS + `RBE3` 실측값 `-2.6650569048944446`은 이 해석해와 **`0.2140%`** 오차로 일치합니다.

## 3. 판정 및 구현 근거

- **PASS**: 곡선 외팔보(`0.2140% ≤ 5%`)와 `v03.json` 직선 외팔보(`0.0249% ≤ 0.5%`) 모두 목표 허용오차를 큰 폭으로 통과했습니다.
- 관련 모듈: [bodyfit.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/meshio_ext/bodyfit.py), [interactions.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/interactions.py), [assembly_builder.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/assembly_builder.py), [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L553-L725).
