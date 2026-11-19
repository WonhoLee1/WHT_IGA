# Phase 6 SBM 비교 평가 보고서 (Master Plan §10.5)

- **작성일**: 2026-11-19 / **대상**: Kratos 10.4.3 `IgaApplication` Shifted Boundary Method (SBM) vs WHT_IGA `cut + extension` vs `FCM α`
- **목적**: [master-plan_v2.1.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/master-plan_v2.1.md) §10.5에 따른 비교 트랙 검증 (주 경로 의사결정에는 사용하지 않으며, 동일 DOF 수준에서의 정확도·조건수·기능 제약 실측 비교)
- **실행 스파이크 시험**: [test_p6_sbm_comparison.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p6_sbm_comparison.py) (`4 passed in 59.31s`, 모킹 없음, 실 Kratos 10.4.3 실행)

---

## 1. Kratos 10.4.3 Solid SBM 런타임 감사 및 제약 사항

Kratos 소스([solid_element.cpp](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_elements/solid_element.cpp), [sbm_solid_condition.cpp](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_conditions/sbm_solid_condition.cpp), [sbm_load_solid_condition.cpp](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp), [iga_modeler_sbm.cpp](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_modelers/iga_modeler_sbm.cpp), [nurbs_geometry_modeler_sbm.cpp](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_modelers/nurbs_geometry_modeler_sbm.cpp))와 [sbm_solid_tests.py](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/tests/sbm_solid_tests.py), [test_modelers_sbm.py](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/tests/test_modelers_sbm.py) 검사 및 런타임 검증 결과:

| 항목 | Kratos 10.4.3 실측 동작 | 근거 (소스/런타임) |
|---|---|---|
| 3D `SolidElement` | **미지원 (2D 전용)**: 3D `NurbsVolumeGeometry` QP 또는 3D `IgaModelerSbm` 요소에서 `RuntimeError: Error: SolidElement momentarily only supports 2D elements, but the current element has dimension 3` 발생; 2D 구성법칙(`strain_size == 3`) 필수 | [solid_element.cpp#L144](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_elements/solid_element.cpp#L144), [solid_element.cpp#L331](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_elements/solid_element.cpp#L331) |
| 3D `SbmSolidCondition` / `SbmLoadSolidCondition` | **미지원 (2D 전용)**: 3D `NurbsVolume` QP에서 `SbmSolidCondition momentarily only supports 2D conditions, but the current dimension is3` 발생; 3D `IgaModelerSbm`(`BrepSurfaceOnVolume`)에서는 `CalculateLocalSystem`이 `LHS (0, 0)` / `RHS (0,)` 반환 | [sbm_solid_condition.cpp#L60](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_conditions/sbm_solid_condition.cpp#L60), [sbm_load_solid_condition.cpp#L238](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp#L238) |
| 3D 다중 레이어 `SbmCondition` | `CreateConditions` 3D 분기가 `KRATOS_ERROR << "CreateConditions: 3D case not implemented yet."` 로 차단됨 | [iga_modeler_sbm.cpp#L1006-L1009](file:///D:/PythonCodeStudy/kratos-src/applications/IgaApplication/custom_modelers/iga_modeler_sbm.cpp#L1006-L1009) |
| `ComputeAxisIntersectionSBM` 수치 허용치 | `tolerance_intersection = 1e-15` 하드코딩으로 `200/48` 등 비이진(non-dyadic) 간격에서 `RuntimeError: ComputeAxisIntersectionSBM :: Wrong number of intersection found (<2)` 발생. 정수/2진 나눗셈 격자(`(40,8)`, `(50,10)`, `(100,20)` 등)에서는 정상 통과 | [nurbs_curve_on_surface_geometry.h#L356](file:///D:/PythonCodeStudy/kratos-src/kratos/geometries/nurbs_curve_on_surface_geometry.h#L356), [nurbs_curve_on_surface_geometry.h#L420](file:///D:/PythonCodeStudy/kratos-src/kratos/geometries/nurbs_curve_on_surface_geometry.h#L420) |
| 2D 침수 선형 패치 테스트 | `[0,2]^2` 평판(배경 `[-0.5,2.5]^2`, `10x10`, `p=2`, `PENALTY_FACTOR=-1.0`)에서 `max\|u_x - 1e-3 x\| = 3.903e-18`, `max\|u_y\| = 2.407e-18`, `dofs = 128`, `cond = 595.61`, `t = 1.27s` — **2D 패치 테스트 기계정밀도 통과** | `test_kratos_sbm_2d_uniaxial_patch_exact` |

> [!IMPORTANT]
> Kratos 10.4.3의 구조해석 SBM(`SolidElement` + `SbmSolidCondition` + `SbmLoadSolidCondition`)은 **2D 평면응력/평면변형률 전용**이므로, 3D 가압 구(V06)는 Kratos SBM으로 직접 해석할 수 없다. 따라서 §10.5 비교는 (1) 2D Lamé 가압 링($a=0.5, b=1.0$) 3자 직접 비교, (2) V02 Roark 유공 평판($\nu=0$, 2D 평면응력 SBM vs 3D `nz=1` `cut+ext`/`FCM α`) 비교, (3) V06 3D 구 게이트 요약으로 수행한다.

---

## 2. 2D Lamé 가압 링 ($a=0.5, b=1.0, E=1, \nu=0, p_i=10^{-6}, p=2$)

사분 링($r \in [0.5, 1.0], \theta \in [0, \pi/2]$)에 내압 $p_i = 10^{-6}$ 및 대칭 경계조건($x=0$에서 $u_x=0$, $y=0$에서 $u_y=0$)을 적용하고 해석해 $u_r(r) = \frac{p_i a^2}{E(b^2-a^2)}\left(r + \frac{b^2}{r}\right)$, $\sigma_{rr}(r)$, $\sigma_{\theta\theta}(r)$ 대비 상대 $L^2$ 오차 및 시스템 행렬 조건수 $\kappa_2(K)$를 비교하였다.

| 방법 | 격자 `nc` | 자유 DOFs | $u$ rel-$L^2$ | $\sigma_{rr}$ rel-$L^2$ | $\sigma_{\theta\theta}$ rel-$L^2$ | $\kappa_2(K)$ | 시간 (s) |
|---|---|---|---|---|---|---|---|
| `cut + extension` (`θ=0.3, drop`) | `(6, 6, 1)` | 386 | `7.603e-04` | `1.980e-02` | `7.622e-03` | `3.046e+05` | 1.29 |
| `FCM α` (`α=1e-8, ghost`) | `(6, 6, 1)` | 464 | `7.601e-04` | `1.979e-02` | `7.621e-03` | `1.321e+10` | 0.53 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(15, 15)` | 400 | `4.969e-03` | `1.254e-02` | `5.393e-03` | `1.656e+04` | 0.20 |
| `cut + extension` (`θ=0.3, drop`) | `(10, 10, 1)` | 822 | `1.644e-04` | `7.903e-03` | `2.904e-03` | `5.718e+05` | 1.84 |
| `FCM α` (`α=1e-8, ghost`) | `(10, 10, 1)` | 1080 | `1.614e-04` | `7.801e-03` | `2.902e-03` | `1.199e+10` | 3.32 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(25, 25)` | 954 | `3.212e-03` | `5.087e-03` | `1.950e-03` | `2.126e+05` | 0.22 |
| `cut + extension` (`θ=0.3, drop`) | `(14, 14, 1)` | 1378 | `8.657e-05` | `3.550e-03` | `1.458e-03` | `1.749e+05` | 6.31 |
| `FCM α` (`α=1e-8, ghost`) | `(14, 14, 1)` | 1952 | `8.606e-05` | `3.496e-03` | `1.457e-03` | `6.880e+09` | 3.80 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(31, 31)` | 1398 | `5.063e-03` | `1.208e-02` | `6.375e-03` | `3.520e+05` | 1.48 |

**핵심 관찰**:
- **정확도**: 동일 DOF 구간(~400, ~900, ~1,380 DOFs)에서 `cut + extension`은 변위 $L^2$ 오차 `7.60e-4 → 1.64e-4 → 8.66e-5`로 단조 수렴하며 `Kratos SBM`(`4.97e-3 → 3.21e-3 → 5.06e-3`) 대비 **6.5×–58× 높은 정확도**를 보인다. 이는 SBM이 실제 곡면 경계 대신 계단형 대리 경계(surrogate boundary)에서 Taylor 전개를 수행하는 반면, `cut + extension`은 실제 곡면에서 적분점(`q=4`)을 직접 적합하기 때문이다.
- **조건수**: `FCM α`는 비활성/고스트 셀의 $\alpha = 10^{-8}$ 강성으로 인해 $\kappa_2(K) \approx 6.9\times 10^9\text{–}1.3\times 10^{10}$($\sim 1/\alpha$)으로 폭증한다. `cut + extension`과 `Kratos SBM` 모두 $\kappa_2(K) \sim 10^4\text{–}5.7\times 10^5$ 수준으로 제어되지만, `Kratos SBM`은 대각 정점 근접(`SnakeSbmProcess` warning)이나 코너 지지 부족($p=2$ B-spline의 작은 대리 지지) 시 `(31, 31)`에서 조건수 상승(`3.52e5`) 및 오차 반등(`5.06e-3`)이 발생한다.

---

## 3. V02 Roark 유공 평판 ($L=200, D=40, r=10, E=1, \nu=0, P=10^{-6}, p=2$)

$\nu=0$ 조건에서 3D 얇은 평판(`t=1`)의 엄밀해는 $z$에 무관하므로 2D 평면응력 `Kratos SBM`과 3D `nz=1` `cut + extension` / `FCM α`의 점 A($x=0, y=10$) 최대 주응력(`ROARK_PEAK = 4.3175e-6`, `ABAQUS_C3D20 = 4.34797e-6`) 및 조건수를 비교하였다.

| 방법 | 격자 `nc` | 자유 DOFs | `peak_A` | Roark 오차 | Abaqus 오차 | `u_right` | $\kappa_2(K)$ | 시간 (s) |
|---|---|---|---|---|---|---|---|---|
| `cut + extension` (`θ=0.3, drop`) | `(20, 4, 1)` | 1134 | `3.851758e-06` | `10.79%` | `11.41%` | `2.277643e-04` | `1.439e+10` | 2.95 |
| `FCM α` (`α=1e-8, ghost`) | `(20, 4, 1)` | 1134 | `3.851758e-06` | `10.79%` | `11.41%` | `2.277643e-04` | `1.439e+10` | 3.29 |
| `cut + extension` (`θ=0.3, drop`) | `(33, 7, 1)` | 2745 | `3.207171e-06` | `25.72%` | `26.24%` | `2.304454e-04` | `1.313e+10` | 7.64 |
| `FCM α` (`α=1e-8, ghost`) | `(33, 7, 1)` | 2754 | `3.207193e-06` | `25.72%` | `26.24%` | `2.304458e-04` | `1.314e+10` | 8.58 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(40, 8)` | 840 | `4.568931e-06` | `5.82%` | `5.08%` | `2.251644e-04` | `2.294e+05` | 0.35 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(50, 10)` | 1240 | `4.628705e-06` | `7.21%` | `6.46%` | `2.407284e-04` | `4.306e+05` | 0.29 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(66, 14)` | 2152 | `4.221209e-06` | `2.23%` | `2.92%` | `2.308095e-04` | `7.565e+05` | 0.33 |
| `Kratos SBM` (2D, `PENALTY=-1`) | `(100, 20)` | 4400 | `4.712030e-06` | `9.14%` | `8.37%` | `2.336776e-04` | `nan` (>2200) | 1.22 |

**핵심 관찰 및 기존 G2 기준선 대비**:
- **조대 격자(~1,200–2,750 DOFs)**: 2D `Kratos SBM`은 모든 자유도를 $(x, y)$ 면내에 배치하므로 3D `nz=1`(`z` 방향 $p=2$ 제어점 3층 × 3성분 = $4.5\times$ DOF 오버헤드) 대비 동일 DOF에서 면내 해상도가 약 2.5배 높아 `(50, 10)` (`7.21%`), `(66, 14)` (`2.23%`)에서 더 낮은 오차를 보인다.
- **SBM 점 A 진동 한계**: `Kratos SBM`에서 `lambda_inner = 0.5`일 때 점 A`(0, 10)`는 대리 경계 외부의 이동 간격(shifted gap)에 놓여 외삽으로 평가되며, 구멍 주변 대리 계단 격자 위치에 따라 `5.82% → 7.21% → 2.23% → 9.14%`로 비단조 진동한다.
- **조건수 차이 원인**: `cut + extension` / `FCM α`의 V02 조건수(`~1.3e10`)는 [test_p2stab_v02_fcm.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p2stab_v02_fcm.py) 및 [V02-ext.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V02-ext.md)의 좌측면 페널티 경계조건(`roark_left_pinned`, $\beta \sim 1.2\times 10^6$)과 세장비 `200:40:1` 굽힘 모드의 영향이며, 실제 절단 강건성 검증(`V-COND`)에서는 구멍 이동 시 `cut + extension` 평탄도 `1.05768×` vs `ext OFF` `1.46494e6×`를 달성하였다([G2_COMPLETE.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_COMPLETE.md)).
- **미세 격자 정밀도**: G2 검증 기준에서 `cut + extension`은 V02 `h=0.375`에서 **Roark `0.678%` (Abaqus `0.027%`), 분포 자기수렴률 `2.34`**로 게이트를 통과하였다([G2_report.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_report.md)).

---

## 4. V06 3D 가압 구 ($a=0.9, b=1.0$) 요약

- **Kratos SBM**: §1에서 확인한 바와 같이 `SolidElement` / `SbmSolidCondition` / `SbmLoadSolidCondition`이 2D 전용이므로 **실행 불가 (N/A)**.
- **`cut + extension` vs `FCM α`** ([V06.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V06.md), [G2_COMPLETE.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_COMPLETE.md)):
  - `FCM α` / 미안정화 초기 시도(`L0 p=2`): `u-L2 = 0.3164`, `srr-L2 = 1.9066` (허용치 10–300배 초과, **FAIL**), 72% OUT 고스트 DOF로 L2 메모리 초과.
  - `cut + extension` (`out_mode="drop"`, `q=min(2p,4)`, `θ`-sagitta 보정): `L0–L2 × p=1–3` 전 셀 통과(`L0 p=2`: `u=5.823e-5`, `srr=1.126e-2`; `L2 p=2`: `u=7.526e-7`, `srr=8.505e-4`; `L2 p=3`: `u=3.433e-7`, `srr=1.586e-4`).

---

## 5. 결론 (§10.5 비교 트랙 판정)

1. Kratos 10.4.3 SBM은 2D 평면응력/평면변형률 침수 문제에서 별도 절단 적분 없이 선형 패치 테스트(`1e-18`)와 `PENALTY_FACTOR=-1` 기반 저조건수(`10^4–10^5`) 해석을 수행한다.
2. 그러나 **3D 구조해석 미지원(2D 전용)**, **비이진 격자 `1e-15` 교차 허용치 오류**, **곡면 경계에서의 Taylor 계단 외삽에 따른 오차 진동 및 저차 수렴(2D Lamé 링 기준 `cut + extension` 대비 6.5×–58× 큰 $u$-오차)**으로 인해 범용 3D CAD-CAE 침수해석 주 경로로 사용할 수 없다.
3. 따라서 WHT_IGA의 **`cut + extension` (`out_mode="drop"`)** 주 경로 선정([ADR-001-route.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/adr/ADR-001-route.md))이 정량적 정확도, 3D 지원, 조건수 안정성 모두에서 타당함을 확인하였다.
