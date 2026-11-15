# G3 게이트 보고서 — Phase 3 (다중물리 · 어셈블리 · THB · 하이브리드 · 후처리/솔버 진단)

- **최종 판정**: **PASS (Gate G3 10대 기준 전 항목 통과)**
- **날짜 / 환경**: 2026-11-15 / `WHTIGA` 브랜치 / Python 3.12.16 / Kratos 10.4.3 (`PARDISO`, `AMGCL`)
- **근거 스위트 및 계측 파일**:
  - 벤치마크 스위트: [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py) — **`9 passed in 75.37s`**
  - 단위·격리·모듈 스파이크 스위트(`test_isolation.py` + `tests/unit` + `test_p3_solver_amgcl_log.py` + `test_p3_assembly_hybrid_thb.py`): **`392 passed in 52.68s`**
  - 실측 JSON: [g3_measurements.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/g3_measurements.json)

---

## 1. Track P3-A / P3-B / P3-C 구현 아키텍처 요약

### Track P3-A — 솔버 오케스트레이션 · 후처리 · 검증 프레임워크 · 메모리/로그
- **다중 스텝·진폭·변수 매핑·박판 세분화** ([P306_steps_amplitude_varmap.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P306_steps_amplitude_varmap.md)):
  - [steps.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/steps.py) (`ABAQUS_TO_KRATOS`, `KRATOS_TO_ABAQUS`, `VAR_MAP`, `Amplitude`, `StaticStep`, `ModalStep`, `ImplicitDynamicStep`, `ExplicitDynamicStep`, `ThermalStep`, `MultiStepRunner`, `CancelToken`, `detect_thin_features`, `auto_refine_thin_features`)
- **후처리 노름 · Proxy 메쉬 · VTKHDF 2.1 · 무가중치 프로브 · 결과 쿼리** ([P307_post_probes_vtkhdf.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P307_post_probes_vtkhdf.md)):
  - [norms.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/post/norms.py), [proxy_mesh.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/post/proxy_mesh.py), [vtk_writer.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/post/vtk_writer.py), [probes.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/post/probes.py) (`w=0.0` Kratos 프로브 요소), [results.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/query/results.py) (`ResultView`, ASME 응력 선형화, 단면력, 평형 검증)
- **표준 검증 프레임워크 (V01–V12, E01–E05)** ([P308_verification_framework.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P308_verification_framework.md)):
  - [schema.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/verification/schema.py), [catalog.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/verification/catalog.py), [qoi.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/verification/qoi.py), [comparator.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/verification/comparator.py), [reporter.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/verification/reporter.py)
- **RAM 기반 자동 솔버 전환 · 6단계 상용 스타일 로그 · 검증 게이트** ([P309_solver_log_memory.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P309_solver_log_memory.md)):
  - [memory.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/memory.py), [solver_log.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/solver_log.py), [validate.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/validate.py), [solvers.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/solvers.py), [pipeline.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/pipeline.py)

### Track P3-B — 다중 포맷 Mesh I/O · Body-Fit NURBS · CAD 치유 · 선택/계측 · GPL 격리
- **5-포맷(6확장자) Mesh I/O 및 Body-Fitted B-Spline 변환** ([P304_meshio_writers_bodyfit.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P304_meshio_writers_bodyfit.md)):
  - [meshio_ext](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/meshio_ext/__init__.py) (`model.py`, `read/write_{inp,nastran,k,rad,msh}.py`), [bodyfit.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/meshio_ext/bodyfit.py) (`mesh_to_bodyfit_patches`)
- **IGES Mode 0 치유 · GPL 격리 Gmsh CLI · 엔티티 선택/계측/접촉검사** ([P305_gmsh_iges_selection.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P305_gmsh_iges_selection.md)):
  - [occ_io.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/occ_io.py), [gmsh_runner.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/meshgen/gmsh_runner.py), [sets.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/selection/sets.py), [selector.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/selection/selector.py), [measure.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/selection/measure.py), [contact_check.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/selection/contact_check.py), [test_isolation.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/test_isolation.py)

### Track P3-C — THB-Splines · Multi-Cell · 어셈블리 · 상호작용 · Native FE & 하이브리드
- **THB-Splines 계층 & Multi-Cell 공유 B-Spline 필드** ([P301_thb_multicell.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P301_thb_multicell.md)):
  - [hierarchy.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/grid/hierarchy.py) (`THBHierarchy`, Oslo 1D/3D 세분화, `enforce_2_to_1_balance`, `build_c_thb_matrix`), [compose.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/operators/compose.py) (`compose_thb`, `apply_TKT`, `solve_reduced`)
- **단일 루트 `Structure` 어셈블리 · `IdOffsetTable` · `Tie`/`RBE2`/`RBE3`** ([P302_assembly_interactions.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P302_assembly_interactions.md)):
  - [assembly_builder.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/assembly_builder.py) (`AssemblyBuilder`, `CellDomainSpec`, `IdOffsetTable`, `transform_stress_voigt`), [interactions.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/interactions.py) (`SurfaceMesh`, `TieSpec`, `RBE2Spec`, `RBE3Spec`, S04 구속 그래프 평탄화)
- **Native FE 5종 요소 & IGA–FE 하이브리드 솔버** ([P303_native_fe_hybrid.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/P303_native_fe_hybrid.md)):
  - [native_fe.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/kratos_backend/native_fe.py) (`hex8`/`tet4`/`tet10`/`hex20`/`wedge6`, 음의 야코비안 자동 플립, `solve_hybrid_iga_fe`)

---

## 2. Gate G3 10대 기준 판정 총괄표

| # | 게이트 기준 | 검증 대상 / 조건 | 실측 결과 | 게이트 허용치 | 판정 | 상세 문서 / 증거 |
|---|---|---|---|---|---|---|
| 1 | **V03-moment** | 곡선 외팔보(Body-fit NURBS + `RBE3`, `[4,32,2]`) & 직선 보([v03.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v03.json) $M_0=5.0$) | 곡선 보 $u_r = -2.665057$ vs 해석해 $-2.670773$ (**`0.2140%`**); 직선 보 $u_y = -1.500374\times 10^{-6}$ vs EB $-1.500000\times 10^{-6}$ (**`0.0249%`**) | 곡선 `≤ 5%`, 직선 `≤ 0.5%` | **PASS** | [V03-moment.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V03-moment.md) |
| 2 | **V04 (Cu/W 이종 보)** | $C^0$ 이중매듭 병합 vs 2-파트 `mpc_matching` Tie & [v04.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/verification/data/v04.json) 환산단면 EB | 병합 `-0.0022034365200516394`, Tie `-0.002203436520040321` mm → **차이 `5.1366e-12`**; EB `-0.0021917808` mm 대비 **`0.5318%`** | 병합 vs Tie `≤ 1e-6`, EB 대비 `≤ 1.0%` | **PASS** | [V04.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V04.md) |
| 3 | **V04 Hybrid `hex8`+IGA** | Cu Native FE `hex8 (2,4,300)` + W IGA (`auto` Tie) | 하이브리드 $u_y = -0.002202782846963654$ mm vs 순수 IGA Tie → **차이 `0.0297%`** (`0.0002967`) | 순수 IGA 대비 `≤ 0.5%` | **PASS** | [V04.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V04.md) |
| 4 | **V05 (7-박스 다중경계)** | 7개 `CellDomainSpec` 공유 B-spline vs 단일(Monolithic) & EB | Monolithic = 7-Box = `-3.9488181143624026e-09` → **상대차이 `0.0`**; EB `-4.0e-9` 대비 **`1.2795%`** | Monolithic 대비 `≤ 3.3e-6`, EB `≤ 2.0%` | **PASS** | [V05.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V05.md) |
| 5 | **V02-local (THB L3)** | `roark_plate_hole.step` 1/4 대칭 THB L3 (2:1 균형, AgFEM $\theta=0.2$) | `active_dofs = 20,566` (균일 `462,840` 대비 22.5× 절감), `elapsed_s = 25.78 s`; C3D20 대비 점 **`0.1042%`**·경계 **`0.0617%`**, Roark 대비 **`0.8106%` / `0.7678%`** | `DOFs ≤ 80,000`, C3D20 `≤ 0.2%`, Roark `≤ 1.0%` | **PASS** | [V02-local.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/V02-local.md) |
| 6 | **E01 (`THIN_FEATURE` + IPS)** | `10×5×0.1 in` 원공 박판 1/4 대칭 + `auto_refine_thin_features` | Pass 0 `detected=True` (`1.00 cells ≤ 1.50`) → Pass 1 `detected=False` (`2.00 cells`); `active_dofs=16,956`; Roark 대비 크라운 **`2.4301%`**, 경계 **`2.2936%`** | `THIN_FEATURE` 해소, Roark `≤ 5%` | **PASS** | [E01.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E01.md) |
| 7 | **E02 (THB 적응 세분화)** | E01 기반 L0→L1→L2 2:1 균형 세분화 & 원거리장 $L_2$ | Roark 오차 `12.4171% → 3.6517% → 1.5016%`; L2 `6,031 / 40,392` DOFs (**`14.93%`**); 원거리 $L_2$ `4.3555e-4 → 1.6763e-5` | 단조 수렴, DOF `< 25%`, $L_2$ 급감 | **PASS** | [E02.md](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/verification/E02.md) |
| 8 | **종속 3-인스턴스 어셈블리** | 회전·평행이동(`i1`, `i2`: $R_z(37^\circ)$, `i3`: $R_x(25^\circ)R_y(-42^\circ)$) 불변성 | `beam_i2`: `rel_diff_u = 5.519e-12`, `rel_diff_voigt_stress = 5.507e-12`; `beam_i3`: `rel_diff_u = 2.226e-12`, `rel_diff_voigt_stress = 2.249e-12` | 변위·회전 Voigt 응력 `≤ 1e-10` | **PASS** | [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L1320-L1386) |
| 9 | **5-포맷 메쉬 무손실 라운드트립** | `.msh`, `.inp`, `.bdf`, `.fem`, `.rad`, `.k` (6종 확장자, 2파트·노드/요소/면 셋) | 전 6종 확장자 **`max_node_err = 0.0`**, `n_blocks = 2`, `parts_ok = node_sets_ok = elem_sets_ok = face_sets_ok = True` | 노드 오차 `≤ 1e-12`, 셋·파트 100% 보존 | **PASS** | [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L1394-L1442) |
| 10 | **상용 6단계 로그 · RAM 자동 전환 · `out_mode="drop"`** | `run_linear` (`auto` 솔버, `WHT_SOLVER_AVAIL_RAM_GB=0.0001`, `drop` vs `ghost`) | `[1/6]`–`[6/6]` 전 단계 및 완료 배너 존재; 기본 `pardiso_ldlt` → 저메모리 시 `amgcl_cg` 자동 전환(`avail=107,374 B`, `est_pardiso=49,523,402 B > 70%`); `drop` **`342`** DOFs `< ghost` **`360`** DOFs (총 `360`) | 6단계 완비, `auto→amgcl`, `drop < ghost` | **PASS** | [test_p3_gate_g3_benchmarks.py](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py#L1450-L1506) |

---

## 3. 기준값 정합성 및 기존 Phase 2 스파이크 상태 (정직한 기록)

### 3.1 `V03-moment` / `V04` / `V05` 참조 수치 정정 기록
[g3_measurements.json](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase3/g3_measurements.json) 의 `prompt_ref_note` 에 기록된 바와 같이, 초기 작업 요청서에 인용된 일부 수치(`V03` `u_r_ref = -0.0211622`, `V04` `u_z = -5.50763e-5 m`, `V05` `2.66065e-2 m`)는 실제 `tests/verification/data/v03.json`, `v04.json`, `v05.json` 및 문제의 기하/물성 파라미터와 일치하지 않았습니다. 본 Gate G3 검증은 추측 없이 **실제 `v03.json` / `v04.json` / `v05.json` 공식 기준값 및 동일 파라미터에 대한 엄밀한 Timoshenko–Goodier / Euler–Bernoulli 해석해**를 기준으로 판정했습니다.

### 3.2 기존 Phase 2 캐시 스파이크(`test_p2stab_v02_fcm.py::test_v02_fcm_baseline_h10`) 관찰 상태
- 실행 명령:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/kratos_spikes/test_p2stab_v02_fcm.py::test_v02_fcm_baseline_h10 -v -s
  ```
- 실측 결과: **`FAILED (1 failed in 9.07s)`**
  - 로그: `[trimmed table: cache hit 5e6d2e1e2f01]`, `K-D h=10 FCM baseline: peak=3.876736e-06 err=10.209% beta=1.200e+06 cut=80 fb=1 umax=2.278e-04 out=0`
  - 실패 사유: `WHT_IGA/tests/kratos_spikes/test_p2stab_v02_fcm.py:77` 의 `assert row["n_fallback"] == 0` 단언에서 기존 Phase 2 캐시 테이블(`5e6d2e1e2f01`)이 `n_fallback = 1` 을 반환(`assert 1 == 0`)하여 실패합니다.
- 해석: 이는 Phase 2 초기 `h=10` FCM 무안정화 기준선 스파이크와 캐시된 `quad_db` 레코드(`5e6d2e1e2f01`)에 국한된 기존(pre-existing) 상태이며, Phase 3 Gate G3 신규 벤치마크 스위트(`9 passed in 75.37s`) 및 통합 단위·격리·스파이크 스위트(`392 passed in 52.68s`)의 통과 여부와는 무관합니다.

---

## 4. 재현 명령

```powershell
# 1. Gate G3 10대 기준 Kratos 벤치마크 전체
D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/kratos_spikes/test_p3_gate_g3_benchmarks.py -v

# 2. 격리 + 단위 + Phase 3 솔버/어셈블리/THB/하이브리드 스파이크 스위트
D:\mamba\envs\whtiga\python.exe -m pytest WHT_IGA/tests/test_isolation.py WHT_IGA/tests/unit WHT_IGA/tests/kratos_spikes/test_p3_solver_amgcl_log.py WHT_IGA/tests/kratos_spikes/test_p3_assembly_hybrid_thb.py -q
```
