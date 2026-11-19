# S46: V07 Kirsch 침지(Trimmed) F-bar / P0 압력 셀 병합 검증 (2026-11-19)

- **모듈**: [`wht_iga/analyses/fbar_trimmed.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/fbar_trimmed.py)
- **시험**: [`tests/unit/test_fbar_trimmed.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_trimmed.py) (Kratos 미사용, 순수 NumPy/Numba/SciPy + `pypardiso`)
- **재현 명령** (`cwd`: `WHT_IGA`):
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest -q -s -p no:cacheprovider tests/unit/test_fbar_trimmed.py
  ```
- **결과**: **`11 passed in 18.87s`**

---

## 1. 목적 및 정식화

사상 격자(`MappedQuarterPlate`, [S42.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S42.md))가 아닌 **직각 배경 `BoxSpline` 격자 + 원공 트림(`build_trimmed_table`)**에서 V07 Kirsch ($\nu = 0.49999$) 문제를 풀기 위한 Simo–Taylor–Pister 평균 팽창 `F̄/P0` 구현:
- `OUT` 셀(`out_cells`) 제거, 절단 셀(`qp_table`)은 moment-fit 적분점, 나머지 활성 셀은 텐서 Gauss 적분점 사용 (`inactive_ctrl_mask` 고정).
- `pressure_groups(eta, eta_thresh=0.25)`: 체적분율 $\eta < 0.25$인 소형 절단 셀을 이웃 대형 활성 셀과 압력 그룹 $g$로 병합해 그룹 평균 팽창 $\bar{J}_g = \frac{\sum_{q \in g} w_q J_q}{\sum_{q \in g} w_q}$ 계산.
- 기저 확장(`extension`) 및 유령 셀(`ghost`) 없이 대칭 Jacobi 스케일링($D = \operatorname{diag}(K)^{-1/2}$) 적용.

---

## 2. 실측 결과 (`test_fbar_trimmed.py`)

### 2.1 접선행렬 FD 정합성 및 아핀 응력 패치
| 시험 항목 | 조건 | 실측 오차 / 수렴 |
|---|---|---|
| 중심 유한차분 접선 (`nu=0.3`) | `n_dof=300`, `12`그룹, `max_group=2`, 그룹결합 `648`쌍 | 전체 `rel_err = 1.551e-10`, 그룹쌍 `1.794e-08` |
| 중심 유한차분 접선 (`nu=0.49`) | `n_dof=300`, `12`그룹, `max_group=2`, 그룹결합 `648`쌍 | 전체 `rel_err = 2.457e-10`, 그룹쌍 `2.002e-08` |
| 아핀 응력 패치 (`eta_thresh=0.25`) | `vol_err = 1.5e-16`, `24/26`그룹, Newton `4`회 | `err_u = 4.17e-12`, `err_sigma = 4.47e-14` |
| 아핀 응력 패치 (`eta_thresh=0.0`) | `vol_err = 1.5e-16`, `26/26`그룹, Newton `4`회 | `err_u = 4.15e-12`, `err_sigma = 5.50e-14` |

### 2.2 V07 1/4 원공 직각 격자 Kirsch (`nu = 0.49999`, `p = 2`)

| 격자 | `eta_thresh` | 자유 DOF | 활성 / 절단 셀 | 압력 그룹 (`min_eta`) | Newton | `e_u` (Hard Gate) | `e_sigma_xx` (Hard Gate) | `cond(D K D)` | `t_table + t_solve` |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| `h=0.5` | `0.25` (병합) | 1,200 | `126 / 78` | `122` (`0.0868`) | `3` | `2.7822e-03` (`<= 0.01`) | `5.6236e-02` (`<= 0.2`) | **`2.290e+05`** | `0.14s + 0.78s` |
| `h=0.5` | `0.0` (셀별 P0) | 1,200 | `126 / 78` | `126` (`0.0868`) | `3` | `1.9181e-03` (`<= 0.01`) | `5.2273e-02` (`<= 0.2`) | `6.453e+05` | `0.14s + 0.91s` |
| `h=0.25` | `0.25` (병합) | 3,888 | `496 / 279` | `490` (`0.0421`) | `5` | `3.2748e-04` (`<= 0.005`) | `2.5867e-02` (`<= 0.1`) | **`1.000e+06`** | `0.24s + 2.95s` |
| `h=0.25` | `0.0` (셀별 P0) | 3,888 | `496 / 279` | `496` (`0.0421`) | `5` | `2.8779e-04` (`<= 0.005`) | `2.5037e-02` (`<= 0.1`) | `1.582e+06` | `0.24s + 1.50s` |

- **수렴 기울기 (`h0=0.5 -> h1=0.25`)**:
  - 병합 (`eta_thresh=0.25`): `slope_u = 3.087`, `slope_sigma_xx = 1.120`
  - 셀별 (`eta_thresh=0.0`): `slope_u = 2.737`, `slope_sigma_xx = 1.062`
- **핵심 관찰**:
  1. 셀별 체적 페널티가 셀 체적에 비례해 스케일되므로, 최소 분율 `0.0421` (`h=0.25`)에서도 셀별 P0(`eta_thresh=0.0`)가 잠김 없이 통과한다(`u_h` 차이 `<= 0.97% max|u|`).
  2. 압력 셀 병합(`eta_thresh=0.25`)의 주된 실측 효과는 대칭 Jacobi 조건수 `cond(D K D)` 개선이다 (`h=0.5` **`0.355×`**, `h=0.25` **`0.632×`**).
  3. 비스케일 `K`는 약지지 B-스플라인(`min_supp = 6.68e-3`) 때문에 `cond ~ 1e11–1e12`이지만 Jacobi 스케일링으로 `1e6` 이하로 해소되며, 약지지 제어점 계수 자체는 `12.6–16.1 × max|u|`까지 커질 수 있으나 물리 영역 $\Omega$ 위 해 `u_h`는 정상 수렴한다.
  4. `P0` 압력 정식화로 응력 수렴 차수는 `~1.1` (`> 0.5` source doc 기준 충족, `>= 2.2` Target 미만)이다.
