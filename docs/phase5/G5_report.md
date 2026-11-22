# G5 게이트 보고서 (Phase 5 검증 — 비선형, 잠김 방지, 접촉)

- **날짜 / 브랜치 / 환경**: 2026-11-20 / `WHTIGA` / Python 3.12 (`D:\mamba\envs\whtiga\python.exe`), Kratos 10.4.3 (`E:\kratos-build\install`), `pypardiso` + Intel MKL
- **기준 문서**: [master-plan_v2.1](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/research-cf-IGA/master-plan_v2.1.md) §9–§10, [G5_summary.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase5/G5_summary.md)
- **종합 판정**: **G5 PASS (전체 5개 정량 벤치마크 V07·E04·V08·E03·E05 통과; Kratos C++ 승격 및 Phase 6 고해상도/배선 항목만 이관)**

---

## 1. Gate G5 벤치마크 판정 요약

| ID | 대상 / 게이트 기준 | 판정 | 핵심 실측 수치 | 상세 보고서 / 테스트 |
|---|---|---|---|---|
| **V07** | Kirsch 원공 ($\nu=0.49999$), 3레벨 오차 한계 + `slope > 0.5` (`plan Target: 2.3/2.2`) | **PASS (Hard Gate & `slope_u` 충족; `slope_sigma` P0 한계)** | 사상 `h2=0.125`: `e_u=8.097e-5 <= 0.002`, `e_sigma_xx=1.379e-2 <= 0.05`, OLS `slope_u=2.492`, `slope_sigma_xx=0.852 > 0.5`; 트림 `h1=0.25`: `e_u=3.275e-4`, `e_sigma_xx=2.587e-2`, `slope=(3.087, 1.120)` | [V07.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V07.md), [S46_fbar_trimmed.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S46_fbar_trimmed.md) |
| **E04** | Cook's membrane ($\nu=0.49999$, TL NH), `u_y` vs `6.98946 mm` `<= 1%` | **PASS** | `wht_iga.analyses.fbar` (`p=2, 16x16`): `u_y = 7.03806 mm` (**`+0.70%`**); User API·`StudySpec`·F-bar 프로브 Cauchy 응력 배선 완료 | [E04.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E04.md), [S40.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S40.md) |
| **V08** | 구형 압자 Hertz–Coulomb 접촉 (`F_N <= 4%`, `F_T RMS <= 4.5%`) | **PASS** | R1 NLGEOM `h=0.125`: `F_N(0.1)=482,739.8 lbf` (**`+0.24%`** vs `481,579.3 lbf`); R2 NLGEOM `0.15<=t<=0.70`: `h=0.25` `rel_RMS=(3.25%, 2.67%)` (pytest), `h=0.125` `rel_RMS=(1.67%, 0.94%)` (외부 Claude 인계); `SD_ghost(1e-8)` 검증 (`eq_cut_T: 1.301e-2 -> 3.541e-5`) | [V08.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V08.md), [S41.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S41.md) |
| **E03** | 유리-금속 밀봉 13셀 3-step 열사이클 (증분 `<= 2x` Coreform `15`, 무컷백) | **PASS** | `(9,6,6) p=2` 확장 B-스플라인 + 환형 RBE3(`cond=7.40`): **`15/15` 수렴(`5/step`, 컷백 `0`)**, 자유열팽창 `4.223e-9 MPa`, `||W_op u_surf|| <= 8.97e-17`, Header `PEEQ_max=4.421e-3`, Pins `7.960e-3`, `571.0 MB` RSS | [E03.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E03.md), [S47_hardening_peeq.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/spikes/S47_hardening_peeq.md) |
| **E05** | 도어 락 3-part 마찰 접촉 (`link RF3 = -46.946 kN ± 5%`) | **PASS** | 실제 `Door_Lock_Assembly.sat` (`h=3.0mm, p=2, mu=0.3`, NLGEOM+J2, [sat_io.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/sat_io.py)): `u_z=5.0mm` `RF3_link = -46.763 kN` vs `-46.946 kN` (**`0.39%`**), `bal_z=2.89e-5`, `PEEQ=0.1310`, `act=(624,535,5)`, `RSS=937 MB`; SVG 전 체크포인트 `<= 3.71%` | [E05.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E05.md), [test_p5_e05_sat_door_lock.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p5_e05_sat_door_lock.py) |

---

## 2. 핵심 구현 및 기술적 발견

1. **잠김 방지 (`F̄/P0`) — 사상 및 침지 트림 격자 모두 확보**:
   - Kratos 내장 `SmallDisplacementMixedVolumetricStrainElement3D4N` / `SmallDisplacementMixedStrainElement3D4N`은 `p>=2` `QuadraturePointGeometry`에서 Access Violation(`returncode=3221225477`)으로 사용 불가를 확정([test_p5_v07_kirsch_precursor.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p5_v07_kirsch_precursor.py)).
   - Numba + CSR 희소 조립 기반 Simo–Taylor–Pister [fbar.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/fbar.py)와 직각 트림 압력 병합([fbar_trimmed.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/fbar_trimmed.py), `eta_thresh=0.25`)으로 **V07** 3레벨·**E04** 1% 게이트를 모두 통과했다.
   - `ReferenceNumbaBackend` + `StudySpec` + `runner/post.py` + `aieng-ui` 브리지(`iga_runner.py`)까지 연결하여 F-bar 해석 시 프로브 `S`/`MISES`/`E`가 선형 `D:eps`가 아닌 Neo-Hookean F-bar Cauchy 응력/Green–Lagrange 변형률로 계산되도록 마감했다([test_fbar_probe_stress.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_fbar_probe_stress.py), [test_aieng_iga_bridge.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_aieng_iga_bridge.py)).

2. **대변형(`NLGEOM`) 표면 접촉 & 유령 셀 안정화**:
   - `TotalLagrangianElement3D4N + CL.HyperElastic3DLaw` 결합 외부 Newton + 내부 Semi-Smooth Newton(SSN) 접촉([contact_nlgeom.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/contact_nlgeom.py))으로 소변형 St. Venant–Kirchhoff 압축 연화(`-8.2%`)를 제거하고 V08 R1 `h=0.125` **`+0.24%`**, R2 `h=0.25` **`(3.25%, 2.67%)`**를 달성했다.
   - `OUT` 유령 셀을 `SmallDisplacementElement3D4N` 기본값으로 전환해 `t=0.70` 기준 `eq_cut_T`를 `1.301e-2 -> 3.541e-5`(**367배 개선**)로 줄이면서도 접촉력 변화를 `< 0.03%`로 억제함을 [test_p5_v08_ghost_check.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p5_v08_ghost_check.py)로 검증했다.
   - 일반 3D 곡면 삼각형 접촉(`DeformableTriangleSurfaceMaster`) + Conforming Red-Green `subdivide_long_edges`(내부 T-junction 균열 `0`개) + CAD 표면 페널티 BC + `ghost_alpha=0.0` 활성 DOF 주입(`_P_act`) + Anderson(`m=4`)·이중 잔차 수렴 + `abort_step()` 적응 이분 서브스텝을 [contact_surface.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/contact_surface.py)에 구현하고, 자체 ACIS `.sat` 임포터([sat_io.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/sat_io.py))로 실제 [Door_Lock_Assembly.sat](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/resources/Door_Lock_Assembly.sat)을 직접 해석해 E05 `-46.946 kN ± 5%` 게이트를 **`0.39%` 오차**(`-46.763 kN`)로 통과했다.

3. **온도 의존 다선형 J2 경화 + `p=2` 확장 + 환형 RBE3 결합 (E03)**:
   - `HARDENING_CURVE = 6` + `PEEQ -> PLASTIC_DISSIPATION` 증분 동기화([thermal_plasticity.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/analyses/thermal_plasticity.py))와 `p=2` 확장 B-스플라인 + QR 피벗 환형 RBE3 구속([test_p5_e03_p2_coupling.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p5_e03_p2_coupling.py))으로 E03 13셀 전체 열사이클(`15/15`, 무컷백)을 단일 프로세스 `571 MB` RSS에서 통과했다.

---

## 3. 정직한 미결/이관 항목 (Open Deferrals)

1. **Kratos C++ 승격/부재 항목 (`WhtIgaApplication` / `ConstitutiveModelsApplication`)**:
   - 다항 초탄성(`Mooney-Rivlin`/`Yeoh`/`Ogden`) 및 `MixedUP0`/`VMS` 요소는 Kratos C++ `QuadraturePointGeometry::Volume()` 미구현·Python 미노출 상태이며, Kratos 10.4.3 내장 열-J2 법칙 부재로 `CellThermalGroups` + `SetInitialStateProcess3D` 우회 경로를 유지한다.
   - `V07`의 `slope_sigma_xx` (`0.852` / `1.120`)는 셀 상수 압력 `P0` 특성에 의한 $O(h)$ 한계로 source doc(`> 0.5`)은 통과하나 master-plan 목표(`2.2`)에는 C++ 고차 혼합(Q2Q1 등) 요소가 필요하다.
2. **시험·배선·고해상도 이관 항목**:
   - `V08` `h=0.125` R2 NLGEOM: pytest 승격 완료 (2026-11-20, `1 passed in 6703 s`, `rel_RMS = 1.66% / 0.93%`), `KratosBodyOperator(..., ghost_alpha=0.0)` 지원 (해소 2026-11-20), `fbar_trimmed.py` → `methods.py` (`full_iga_immersed`) 연결 (해소 2026-11-20), `E03` `0.5 mm` 조밀 메쉬(`146,461` 배경 셀) 및 `E05` `h=2.0 mm` 추가 세분화는 Phase 6 병렬/AMGCL 트랙으로 이관한다.

---

## 4. 직접 검증된 대표 스위트 목록

- `tests/unit/test_fbar_trimmed.py`: `11 passed in 18.87s`
- `tests/unit/test_fbar_probe_stress.py tests/unit/test_aieng_iga_bridge.py tests/unit/test_fbar_pipeline_wiring.py`: `31 passed in 87.43s` (+ `test_user_api.py`, `test_discretization_methods.py` 포함 `89 passed in 108.18s`)
- `tests/unit/test_fbar_mean_dilatation.py`: `20 passed in 25.77s`
- `tests/kratos_spikes/test_p5_v07_kirsch_precursor.py`: `4 passed in 98.57s`
- `tests/kratos_spikes/test_p5_v08_ghost_check.py`: `1 passed in 406.88s`
- `tests/kratos_spikes/test_s44_hardening_initial_strain.py tests/kratos_spikes/test_p5_thermal_plasticity.py`: `25 passed in 30.83s`
- `tests/kratos_spikes/test_p5_e03_p2_coupling.py`: `3 passed in 181.63s`
- `tests/kratos_spikes/test_p5_e05_door_lock_precursor.py`: `3 passed in 4.42s`
- `tests/kratos_spikes/test_p5_e05_nlgeom_hardening.py`: `3 passed in 171.32s`
- `tests/unit/test_sat_io.py`: `2 passed` (+ `test_phase3_meshio_roundtrip_bodyfit_iges.py` `25 passed in 21.30s`)
- `tests/kratos_spikes/test_p5_e05_sat_door_lock.py`: `1 passed in 1504.38s` (`RF3_link = -46.763 kN` vs `-46.946 kN`, **`0.39%`**, `RSS = 937 MB`)
- `tests/kratos_spikes/test_p5_prepomax_step_benchmarks.py`: `4 passed in 29.77s`
