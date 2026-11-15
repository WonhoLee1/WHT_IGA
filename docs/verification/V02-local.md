# V02-local — Roark 원공 평판 THB-Splines 3단계 국부 세분화 (1/4 대칭)

- 상태: **PASS** (Gate G3 기준 전 항목 충족)
- 게이트 (G3):
  - `roark_plate_hole.step` 대상 THB Level 3 국부 세분화 + 2:1 균형 + AgFEM 확장 적용
  - 활성 자유도(active DOFs) **≤ 80,000** (균일 미세 격자 대비 대폭 DOF 절감)
  - C3D20 기준해(`4.34796947956784e-06`) 대비 최대응력 상대오차 **≤ 0.2%**
  - Roark 공식해(`4.3175e-06`) 대비 상대오차 **≤ 1.0%**
- 날짜 / 환경: 2026-11-15, WHTIGA 브랜치, Python 3.12.16, Kratos 10.4.3 (`PARDISO`)
- 참조 명세: [v02.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v02.json)
- 실측 데이터: [g3_measurements.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/g3_measurements.json)
- 재현 명령:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py -k test_g3_04_v02_local_thb_l3_roark_plate_hole -v
  ```

## 1. 입력 요약 및 대칭성 명시

- CAD 입력: [roark_plate_hole.step](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/data/roark_plate_hole.step) (전체 판 `[-100, 100] x [-20, 20] x [-0.5, 0.5]`, 원형 구멍 반경 $r = 10$), `tessellate(..., deflection=0.002, angular=0.02)`
- **1/4 대칭(Quarter Symmetry) 모델링 명시**:
  - 해석 영역: `[0, 100] x [0, 20] x [-0.5, 0.5]`
  - 대칭 경계조건: $x = 0$ 면 $u_x = 0$, $y = 0$ 면 $u_y = 0$, $x = 100$ 면 $u_z = 0$ (강체 Z 억제)
  - 하중: $x = 100$ 우측면에 단축 인장 트랙션 $t_x = 1.0 \times 10^{-6}$
- 재료 물성: $E = 1.0$, $\nu = 0.0$
- THB 계층 및 안정화:
  - 기본 격자: `base_cells = [33, 7, 1]`, `degree = 2`
  - 원공($r = 10$) 주변 반경 대역 `[(6.0, 22.0), (7.5, 16.0), (8.8, 13.0)]` 3단계 세분화(`thb_levels = 3`) + `enforce_2_to_1_balance()` (`balance_added == 0`)
  - 절단 셀 구적: `quadrature.cell_quadrature(..., q=3, p=2)`
  - AgFEM 확장: 최미세 활성 레벨에서 [`classify_basis`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/stabilize/support.py) + [`build_extension`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/stabilize/extension.py) (`theta = 0.2`) 적용

## 2. 메쉬 · 자유도 · 응력 실측 결과

### 2.1 계층 격자 및 DOF 절감 통계
| 항목 | 실측값 | 비고 |
|---|---:|---|
| 최미세 셀 크기 `(hx, hy)` | `(0.3787878787878788, 0.35714285714285715)` | `[100/(33*8), 20/(7*8)]` (`v02.json` L3 `0.375` 급) |
| 총 리프 셀 수 (`n_leaves`) | `6,391` | L0: `183`, L1: `208`, L2: `752`, L3: `5,248` |
| 절단 셀 수 (`n_cut_cells`) | `440` | — |
| AgFEM 확장 슬레이브 수 (`n_slave_extended`) | `190` | $\theta = 0.2$ |
| **THB L3 활성 자유도 (`active_dofs`)** | **`20,566`** | 게이트 `≤ 80,000` **PASS** |
| 동일 1/4 영역 최미세 균일 격자 DOF | `462,840` | THB 대비 **약 22.51× DOF 절감** (4.44%만 사용) |
| 전체 판 균일 `h=0.375` DOF (`v02.json`) | `516,858` | — |
| 전체 해석 소요 시간 (`elapsed_s`) | `25.78 s` (`25.7799954 s`) | 조립 + AgFEM + PARDISO 해석 포함 |

### 2.2 크라운 점(Point A) 및 원공 경계 피크 응력 비교
| 평가 위치 | Kratos THB L3 실측값 | 기준값 (`C3D20` / `Roark`) | 상대오차 vs `C3D20` | 상대오차 vs `Roark` | 판정 |
|---|---:|---:|---:|---:|---|
| 크라운 점 `(0.0, 10.0, 0.5)` $\sigma_{xx}$ | `4.352498956558145e-06` | `4.34796947956784e-06` (`C3D20`) / `4.3175e-06` (`Roark`) | `0.001042` (**`0.1042%`**) | `0.008106` (**`0.8106%`**) | **PASS** |
| 원공 경계 피크 응력 (`boundary_peak_stress`) | `4.350650572761627e-06` | `4.34796947956784e-06` (`C3D20`) / `4.3175e-06` (`Roark`) | `0.000617` (**`0.0617%`**) | `0.007678` (**`0.7678%`**) | **PASS** |

- 허용오차 대조: `C3D20` 대비 `0.1042% / 0.0617% ≤ 0.2%` (**PASS**), `Roark` 대비 `0.8106% / 0.7678% ≤ 1.0%` (**PASS**).

## 3. 판정 및 구현 근거

- **PASS**: 정확한 Giannelli 절단(truncation) 행렬 $C_{\text{lvl}}$, 2:1 균형, 그리고 AgFEM 지지 영역 확장을 결합한 3단계 THB 계층 해석으로 **20,566 활성 자유도(`25.78 s`)** 만에 균일 미세 격자(462,840 DOFs) 수준인 **C3D20 대비 `0.0617%` 오차**를 달성했습니다.
- 관련 모듈: [hierarchy.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/grid/hierarchy.py), [compose.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/operators/compose.py), [extension.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/stabilize/extension.py), [trimmed.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/trimmed.py), [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L1062-L1126).
