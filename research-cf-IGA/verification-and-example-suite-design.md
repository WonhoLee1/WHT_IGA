# WHT_IGA Verification Manual (V01–V12) & Example Problems (E01–E05) 1:1 재현·비교 통합 설계서

- 작성일: 2026-11-14
- 기준 문서:
  - `research-coreform-IGA/master-plan_v2.1.md` (§9, §10.1–10.2)
  - `research-coreform-IGA/user-guide/verification-manual/` (V01–V12 원문 및 `generated-downloads/` CSV 60+종)
  - `research-coreform-IGA/user-guide/example-problems/` (E01–E05 원문 및 `glass-to-metal-seal/data/` CSV)
  - `WHT_IGA/tests/verification/data/` (`v04.json`–`v12.json`, `e01.json`–`e05.json`)

---

## 1. 설계 목표 및 현재 갭(Gap) 분석

### 1.1 목표
Coreform IGA for Abaqus의 **Verification Manual 12건(V01–V12)** 과 **Example Problems 5건(E01–E05)** 에 정의된:
1. **모델 구성·경계·하중·파라미터 스윕 절차**
2. **QoI(관심 물리량) 추출 위치·적분·프로브 경로**
3. **비교 수식·CSV 컬럼 스키마·수락 기준**

을 WHT_IGA에서 **동일한 절차와 동일한 출력 포맷으로 자동 재현·비교**할 수 있는 통합 검증 아키텍처를 설계한다.

### 1.2 현재 갭(Gap) 진단
| 항목 | 현재 상태 | 설계 보완 사항 |
|---|---|---|
| 기준 JSON (`tests/verification/data/`) | `v04`–`v12`, `e01`–`e05`만 존재 (**`v01`–`v03` 누락**) | `v01.json`, `v02.json`, `v03.json` 추가 + 전 17건 JSON에 `procedure`, `probes`, `sweeps`, `reference_csvs`, `criteria` 표준 섹션 통합 |
| 실행 스크립트 | `scripts/v02_roark.py`, `v03_beam.py`, `vcond_roark.py` 개별 스크립트 위주 | 공통 `wht_iga/verification/` 패키지로 스윕·추출·비교·리포팅 표준화 |
| CSV 비교 포맷 | Coreform `generated-downloads/` CSV와 1:1 컬럼 호환 출력기 부재 | Coreform과 **동일한 CSV 파일명·헤더** 출력 및 자동 diff/오버레이 플롯 생성 |
| Example E01–E05 절차 | `e01`–`e05` JSON에 일부 미확인(`unknown`) 스윕·압력 정보 잔존 | 원문 `.md` 및 `glass-to-metal-seal/data/` 전수 반영하여 실행 절차 확정 (`E01` $P = -10^{-6}\text{ psi}$, $h=(0.3,0.3,0.1)$, `Probe-Stress-Profile` 100점 등) |

---

## 2. 통합 검증 프레임워크 아키텍처 (`wht_iga/verification/`)

```text
tests/verification/data/{v01..v12,e01..e05}.json   (Case Spec + Reference Values)
research-coreform-IGA/user-guide/generated-downloads/.../*.csv   (Coreform Reference CSV)
                    │
                    ▼
        ┌───────────────────────┐
        │  runner.py            │  1. 형상/재료/Step/BC/Load 조립 (CaseBuilder)
        │  (Sweep Orchestrator) │  2. (p, C, h, dt, workflow) 매트릭스 실행
        └───────────┬───────────┘  3. trim/solve 성능(wall_s, DOF, active/trimmed) 수집
                    │
                    ▼
        ┌───────────────────────┐
        │  qoi.py               │  7종 표준 QoI 추출기:
        │  (QoI Extractors)     │   - Point / Face-Grid / Interior-Grid Probe
        └───────────┬───────────┘   - Line Probe + L2 Profile / Principal Stress
                    │               - Volume Quadrature L2 + Log-Log Slope
                    │               - Contact Force/Pressure & Sliding RMS
                    │               - Modal Participation + PSD RMS + SRSS Peak
                    │               - Transient History (f1, Phase, Decay, Samples)
                    │               - Multi-Step Nonlinear (PEEQ, Sigma1+, Sensitivity)
                    ▼
        ┌───────────────────────┐
        │  comparator.py        │  1. Coreform CSV / Analytical / Native FE 비교
        │  & reporter.py        │  2. Coreform 동일 헤더 CSV 출력 + PASS/FAIL 판정
        └───────────────────────┘  3. docs/verification/<ID>.md + PNG 그래프 생성
```

### 2.1 표준 JSON 스키마 확장 (`tests/verification/data/<id>.json`)
기존 `v04`–`v12`, `e01`–`e05` JSON의 하위 호환을 유지하면서 다음 5개 블록을 표준화한다:
1. `meta`: `id`, `category` (`verification` | `example`), `phase_gate` (`G1`–`G6`), `unit_system`, `reference_doc`
2. `model`: `geometry`, `materials`, `steps`, `boundary_conditions`, `loads`, `interactions`
3. `sweeps`: `degrees` (`p`), `continuities` (`C`), `mesh_sizes` (`h`), `padding`, `workflows` (`immersed`, `body_fit`, `native_fe`, `boolean_merge`, `tie`, `local_refinement`), `time_steps`
4. `probes`: `point_probes`, `face_grids`, `line_probes` (`start`, `end`, `n_points`, `fields`), `quadrature_qois`, `history_outputs`
5. `acceptance`: `criteria` (지표별 수식/허용치/레벨별 상한/제외 규칙), `reference_csvs` (Coreform 기준 CSV 상대 경로 매핑)

---

## 3. 7대 표준 QoI 추출기 설계 (`wht_iga/verification/qoi.py`)

기존 `wht_iga/post/` (`norms.py`, `vtk_writer.py`) 및 `wht_iga/query/results.py` (`ResultView`, `probe`, `check_equilibrium`) 위에 구축한다.

| 추출기 클래스 / 함수 | 적용 스터디 | 추출 절차 및 수식 |
|---|---|---|
| `FaceAndGridProbeExtractor` | **V01**, V11, V12, E04 | 지정 면 위 $N\times N$ 균일 격자(`21×21`) + 내부 3D 격자(`{0.1,0.5,0.9}^3`)를 가중치 0 적분점 `Probe`로 평가. `probe_coord_err`, `max_abs_err` 계산 |
| `LineProfileExtractor` | **V02**, **V03**, **V04**, **V05**, **E01**, **E02** | `start -> end` 100점 라인 프로브에서 $U, S, \sigma_{\max}=\max\operatorname{eig}(S), \sigma_{\text{vM}}$ 추출. 선분 사다리꼴/Gauss 적분으로 프로파일 상대 $L^2$ 오차 계산 |
| `VolumeQuadratureL2Extractor` | **V06**, **V07** | 적분점 $(\mathbf{x}_q, w_q \det J_q, \mathbf{u}_h, \boldsymbol{\sigma}_h)$에서 해석해 대비 상대 $L^2$ 적분($\sqrt{\sum wJ(q_h-q_*)^2 / \sum wJ q_*^2}$), 체적 오차, $\log e$ vs $\log h$ 최소제곱 기울기 산출 |
| `ContactHistoryExtractor` | **V08**, **E05** | `RF1`, `RF2`, `RF3`, 최대 접촉압력 $p_0$ 시계열 추출. 지정 시간창($0.15\le t\le 0.70$) 상대 RMS 오차 $\sqrt{\sum(q-q_*)^2}/\sqrt{\sum q_*^2}$ 및 합산 반력 계산 |
| `ModalAndSpectrumExtractor` | **V09** | 20개 고유진동수·질량정규화 모드·모달참여율 $\Gamma_y$ 추출. $|\Gamma_y| \ge 0.01 \max|\Gamma_y|$ 참여 모드만 필터링, 랜덤 PSD RMS 및 SRSS 충격 피크 변위 산출 |
| `DynamicHistoryExtractor` | **V10** | 자유단 평균 축변위/속도/가속도 시계열에서 자유진동 구간 $f_1$ 추출, 정규화 변위 RMS, 위상 RMS, 진폭·에너지 감쇠비 차이, 허용 증분 수 비율, 주기당 샘플 수 계산 |
| `NonlinearThermoExtractor` | **E03** | 다단계 Step별 `PEEQ` 최대값, 유리 양의 최대 주응력 $\sigma_1^+ = \max(\sigma_1, 0)$, 단면 뷰컷($y=0$), 계면 경로 프로브, 수렴 증분 이력(`convergence-history`) 수집 |

---

## 4. Verification Manual (V01–V12) 1:1 재현·비교 명세

### 4.1 V01: `prescribed-displacement-cube` (Phase 1 / `v01.json` 신규)
- **모델·스윕**: 단위 큐브 $[0,1]^3$, $E=1, \nu=0$, 배경 격자 중심 $(0.5,0.5,0.5)$, $h=0.6$ (격자면 $-0.4, 0.2, 0.8, 1.4$ → 내부 $[0.2,0.8]^3$ 1개만 untrimmed). 6개 기저 조합 $(p,C) \in \{(1,0),(2,0),(2,1),(3,0),(3,1),(3,2)\}$.
- **BC**: 6개 면 직교 변위($x=0\to u_x=0$, $x=1\to u_x=1$, $y=0\to u_y=0$, $y=1\to u_y=2$, $z=0\to u_z=0$, $z=1\to u_z=3$; 나머지 성분 자유).
- **비교 항목·절차**:
  1. 6개 면 각각 $21\times 21 = 441$점 Face-probe → $u$ 잔차 및 $\sigma_{nn}$ 잔차 (`<= 1.51e-5`)
  2. 표면층(`Triangle3D3`, Coreform `SFM3D3` 대응) 노드 변위 추출 → `<= 1.42e-5`
  3. 내부 27개 점 $\{0.1,0.5,0.9\}^3$ 프로브 → 6개 응력 성분 전체 최대 절대 오차(`<= 5.46e-5`, 전단응력 $0$ 확인)
  4. 프로브 좌표 정합 오차 `0` 검증
- **출력 테이블 헤더**: `Basis, Probe-location max. abs. error, Face-probe U max. abs. error, Surface U max. abs. error, Interior S max. abs. error`

### 4.2 V02: `roark-17p1-7a` (Phase 1 전역 + Phase 3 국소 THB / `v02.json` 신규)
- **모델·스윕**: $L=200, D=40, t=1, r=10$ ($2r/D=0.5$), $E=1, \nu=0$, 좌측 `PINNED`, 우측 $P=10^{-6}$.
  - 전역 스윕: $p \in \{1,2,3\} \times h = 6/2^k$ ($k=0..4$: $6.0, 3.0, 1.5, 0.75, 0.375$).
  - 국소 스윕 (Phase 3): $h=3.0, p=2$ 베이스에서 `hole` 면 주변 THB $L=1,2,3$ (2:1 균형, offset $p\cdot h/2^L$).
- **비교 항목·절차**:
  1. Point A $(L/2, r, t/2)$에서 최대 주응력 $\sigma_{\max}$ → Roark($4.3175\times 10^{-6}$) 및 Reference FE($4.34797\times 10^{-6}$) 대비 상대 오차
  2. 횡단 라인 프로브($x=L/2, y\in[10,20], z=t/2$) → $\sigma_{\max}(y), u_x(y), u_y(y)$의 최고 해상도($h=0.375$) 대비 자기 수렴 상대 $L^2$ 오차 및 수렴 기울기($\ge 1.8$)
  3. 국소 L3에서 peak 오차 $\le 0.2\%$, DOF $\le 80\text{k}$ (Coreform 참고: $0.105\%$ @ $57,300$ DOF)
- **생성 CSV (Coreform 동일 헤더)**:
  - `roark-17p1-7a-qoi.csv`: `test_case,h,p,peak_stress,peak_rel_err_vs_roark,peak_rel_err_vs_reference,stress_profile_rel_l2,ux_profile_rel_l2,uy_profile_rel_l2`
  - `roark-17p1-7a-refinement.csv`: `level,finest_h,local_test_case,uniform_test_case,local_dofs,uniform_dofs,local_peak_rel_err_vs_reference,uniform_peak_rel_err_vs_reference,local_peak_rel_err_vs_roark,uniform_peak_rel_err_vs_roark,local_stress_profile_rel_l2,uniform_stress_profile_rel_l2,num_cpus,local_trim_wall_s,uniform_trim_wall_s,local_solve_wall_s,uniform_solve_wall_s`
  - `roark-17p1-7a-trim-performance.csv`, `roark-17p1-7a-solve-performance.csv`

### 4.3 V03: `beam-bending` (Phase 1 기본 + Phase 3 moment / `v03.json` 신규)
- **모델·스윕**: $b=0.1, h=0.1, L=1.0\text{ m}$, $E=200\text{ GPa}, \nu=0, \rho=8000\text{ kg/m}^3$. 격자 $(0.1,0.1,0.25)/2^k$ ($k=0..4$), $p \in \{1,2,3\}$, $C^{p-1}$, 패딩 $p$.
- **6종 하중 + 2종 트랙션 변형**: `gravity`, `body_force`, `pressure_uniform`, `pressure_total_force`, `traction` (+ `follower`, `resultant` — Phase 5), `moment` ($M_0=5\text{ N}\cdot\text{m}$ RP coupling — Phase 3).
- **비교 항목·절차**:
  1. 끝단 처짐 $u_y(0,0,L)$ vs Euler-Bernoulli($-wL^4/8EI$, $-PL^3/3EI$, $-M_0L^2/2EI$) & Native FE 참조해 (`beam-bending-accepted.csv`)
  2. 3개 라인 프로브: 중심축 $u_y(z)$, 상단 섬유 $\sigma_{zz}(z)$, 중앙 단면 $\sigma_{zz}(y), \tau_{zy}(y)$ → $h_4$ 대비 자기 수렴 $L^2$
  3. 점근 레벨($h_3, h_4$) 오차 $\le 5\%$ ($p=1$ 조대 격자 shear locking 구간은 기록 전용)
- **생성 CSV**: `beam-bending-<load>-tip-deflection.csv`, `beam-bending-<load>-trim-performance.csv`, `beam-bending-<load>-solve-performance.csv`, `beam-bending-traction-variants.csv`

### 4.4 V04: `bimetallic-beam-bending` (Phase 3 G3)
- **모델·스윕**: Cu($E_f=100\text{ GPa}, h_f=5$) + W($E_s=400\text{ GPa}, h_s=5$), $\nu=0, b=10, L=100\text{ mm}$, $I_t = 1520\text{ mm}^4$. 격자 $(5,5,10)/2^k$ ($k=0..3$), $p\in\{1,2,3\}$. 5종 하중(`gravity`, `body_force`, `pressure`, `traction`, `moment`).
- **비교 워크플로우 3종**:
  1. `boolean_merge` (공유 스플라인 장 / 멀티셀)
  2. `tie` (독립 파트 2개 + 계면 tie)
  3. **IGA–FE 혼합(D14)**: FE 파트(Cu) + IGA 파트(W) tie
- **비교 항목·절차**:
  - 끝단 $u_y(0,0,L)$ vs 변환단면 해 (`bimetallic-beam-accepted.csv`, $\le 5\%$), `|merge - tie|` 상대차 $\le 10^{-6}$
  - 3개 라인 프로브: 계면 $u_y(z)$, 중앙 단면($z=50$) $\sigma_{zz}(y)$($n=4$ 불연속 도약) 및 $\tau_{zy}(y)$(연속), 고정단($z=0$) 응력
- **생성 CSV**: `bimetallic-beam-tie-comparison.csv`, `bimetallic-beam-<load>-tip-deflection.csv`, trim/solve performance CSV

### 4.5 V05: `multiplanar-material-interfaces` (Phase 3 G3)
- **모델·스윕**: 7종 단면(`monolithic`, `bottom-center-third`, `bottom-left-half`, `centered-third`, `zig-zag`, `nested-l`, `tetris`) $\times 4$ 격자($k=0..3$) $\times 3$ 차수($p=1..3$) = **84 해석**. 동일 물성($E=100\text{ GPa}, \nu=0.3$), 끝단 $P=1\text{ N}$.
- **비교 항목·절차**: 동일 $(h,p)$의 `monolithic` 대비 6개 상대 오차 QoI (`tip_displacement_rel`, `axial_displacement_l2_rel`, `vertical/horizontal_bending_stress_l2_rel`, `vertical/horizontal_shear_stress_l2_rel`) 추출 → $k=3$에서 Coreform 수준(p=2 $\le 3.3\times 10^{-6}$) 확인.
- **생성 CSV**: `multiplanar-material-interface-<case>-responses.csv`, trim/solve performance CSV

### 4.6 V06: `pressurized-sphere` (Phase 2 G2)
- **모델·스윕**: 1/8 구 ($a=0.9, b=1.0, E=1, \nu=0, p_i=10^{-6}$), 대칭면 3개, $h = 0.1/2^k$ ($k=0,1,2$), $p\in\{1,2,3\}$.
- **비교 항목·절차**: 적분점에서 구좌표 투영($u_r, \sigma_{rr}, \sigma_{\theta\theta}$) 후 Lamé 엄밀해 대비 전역 상대 $L^2$ 오차 4종(`displacement`, `radial_displacement`, `radial_stress`, `tangential_stress`)을 `pressurized-sphere-qoi.csv`의 `(p, level, h, quantity)`별 `tolerance` 행렬과 1:1 비교.

### 4.7 V07: `infinite-plate-with-hole` (Phase 5 G5)
- **모델·스윕**: 1/4 대칭 $[0,4]^2 \times 0.25$, $a=1, E=10^5, \nu=0.49999, T=10$. 외부 경계에 Kirsch 해석해 변위 직접 부여(penalty $10^{10}$). Mixed $U(2,C^1)/P(2,C^1)$ + VMS (Phase 5 C++ `WhtIgaApplication`). $h \in \{0.5, 0.25, 0.125\}$.
- **비교 항목·절차**: 체적 적분 상대 $L^2$ 오차 $(e_u, e_{\sigma_{xx}})$, 적분 체적 오차($< 0.1\%$), log-log 수렴 기울기(Hard Gate $> 0.5$, Target $2.3/2.2$), 레벨별 상한 `(0.01,0.2), (0.005,0.1), (0.002,0.05)` 검증 후 `global_l2_metrics.json` 출력.

### 4.8 V08: `contact-indenter` (Phase 5 G5)
- **모델·스윕**: $5\times 2\times 3\text{ in}$ 블록 + 중공 반구($R=0.5, R_i=0.25$, 초기 간극 $0.025\text{ in}$), `NLGEOM=YES`, $\mu=0.3$. Step 1 ($0\le t\le 0.1$: $-0.125\text{ in}$ 압입), Step 2 ($0.1\le t\le 1.0$: $+4.25\text{ in}$ 미끄럼). $h \in \{0.5, 0.25, 0.125\}$, $p=2, C^1$ + 표면층 접촉.
- **비교 항목·절차**: Hertz 이론해($F_N=481{,}579.3\text{ lbf}, p_0=4{,}598{,}743\text{ psi}, F_T=144{,}473.8\text{ lbf}$) 대비 $t=0.1$ 시점 $F_N, p_0$ 오차 및 미끄럼 구간($0.15\le t\le 0.70$) $F_N, F_T, p_0$ 상대 RMS 오차 5개 지표($F_N \le 4\%$, RMS $\le 4.5\%$) 검증.

### 4.9 V09: `cantilever-shock-vibration` (Phase 4 G4)
- **모델·스윕**: $0.1\times 0.1\times 1.0\text{ m}$, $E=200\text{ GPa}, \nu=0, \rho=8000\text{ kg/m}^3$, `ENCASTRE` at $z=0$. Body-fit Native FE($50\text{ mm}$) vs Immersed IGA($25\text{ mm}, p=1$).
- **비교 항목·절차**:
  1. 고유치 20모드 추출 + $|\Gamma_y| \ge 0.01\max|\Gamma_y|$ 참여 모드 상대차 $< 1\%$ (`eigenfrequencies.csv`)
  2. 랜덤 진동($10\sim 2000\text{ Hz}$, PSD $0.01\text{ g}^2/\text{Hz}$, $\xi=0.02$) $u_y$ RMS & 응답 스펙트럼(SRSS, $\xi=0.02$) 피크 변위: 구현 오차 $< 5\times 10^{-6}$ (`implementation.csv`), 이산화 차이 $< 5\%$ (`response.csv`)

### 4.10 V10: `vibrating-beam` (Phase 4 G4)
- **모델·스윕**: $10\times 2\times 2$ 봉, $E=10^6, \rho=10^{-6}, \nu=0$. 3종 이산화(Native FE vs Body-fit IGA vs Immersed IGA, $p=2, C^1, h=1$) $\times$ 2종 시간적분 프로파일(`TRANSIENT_FIDELITY`, `MODERATE_DISSIPATION` — Kratos Bossak $\alpha_B$ 스윕, D5). $0\le t\le 2\times 10^{-3}, dt_{\text{init}}=10^{-6}, dt_{\max}=2\times 10^{-6}$, 1주기($4\times 10^{-5}\text{ s}$) 가압 후 자유진동.
- **비교 항목·절차**: 6개 케이스 $\times$ 핵심 지표(총 26개 기준): $f_1$ 상대오차($<2\%$ vs $25{,}000\text{ Hz}$), 정규화 변위 RMS($<0.01$), 위상 RMS($<0.2\text{ rad}$), 진폭 감쇠비 차($<0.03$), 에너지 감쇠비 차($<0.05$), 증분 수 비($\le 1.10$), 주기당 샘플($\ge 19.9$).

### 4.11 V11: `nafems-p18-t1` (Phase 4 G4)
- **모델·스윕**: 1/4 대칭 $10\times 10\times 0.25\text{ mm}$, 반경 $r=1\text{ mm}$ 원형 핫스팟 멀티셀(2셀 공유 장). $44\times 44\times 1$ ($h=0.227\text{ mm}$), $p=1, C^0$, $E=10^5\text{ MPa}, \nu=0.3$, 내측 $\epsilon_{\text{th}}=10^{-3}$, 외측 $0$, 평면응력 등가 구속.
- **비교 항목·절차**: 계면 외측 점 $D(1,0)\text{ mm}$에서 $\sigma_{yy}(D)$ 프로브 → NAFEMS 기준 $50.0\text{ MPa}$ 대비 상대 오차 $< 4\%$.

### 4.12 V12: `nafems-p18-t2/t3/t4` (Phase 4 G4)
- **T2 (정상 복사)**: $L=0.1\text{ m}, k=55.6, T_A=1000\text{ K}$, B단 복사($T_\infty=300\text{ K}, \epsilon=0.98$), $10\times 1\times 1$ $p=1$. $T_B$ vs $927\text{ K}$ ($< 1\%$).
- **T3 (과도 전도)**: $L=0.1\text{ m}, k=35, \rho=7200, c=440.5, T_0=0, T_A=0, T_B(t)=100\sin(\pi t/40)$, $dt=0.1\sim 0.5\text{ s}$. $T_C(x=0.08, t=32\text{ s})$ vs $36.60^\circ\text{C}$ ($< 2\%$).
- **T4 (정상 대류)**: $0.6\times 1.0\times 0.05\text{ m}, k=52$, 하단 $100^\circ\text{C}$, 상/우측 대류($T_\infty=0, h=750$), $12\times 20\times 1$ $p=1$. $T_E(0.6,0.2)$ vs $18.3^\circ\text{C}$ ($< 3\%$).

---

## 5. Example Problems (E01–E05) 1:1 재현·비교 명세

| ID | Phase | 모델·재료·물리 절차 (Coreform 동일 설정) | 프로브 / QoI 추출 및 비교 기준 |
|---|---|---|---|
| **E01** | **Phase 3 (G3)** | $10\times 5\times 0.1\text{ in}$ 강판 + 중심 홀 $r=1\text{ in}$ (IPS, $E=30\times 10^6\text{ psi}, \nu=0.29$). `X_min` PINNED, `X_max` $P=-10^{-6}\text{ psi}$ ($F=5\times 10^{-7}\text{ lbf}$). 직교 격자 $(0.3,0.3,0.1)\text{ in}$, $p=2$, Smooth($C^1$). | `Probe-Stress-Profile`: $(0, 1.0, 0.05) \to (0, 2.5, 0.05)$ **100점 라인 프로브**. $\sigma_1 = \max\operatorname{eig}(S)$ 산출. 구멍 에지($y=1.0$)에서 Roark($K_t=2.23568 \Rightarrow 3.7261\times 10^{-6}\text{ psi}$) 대비 가중치 0 Probe 및 Proxy 격자 노드 값 병행 비교 (Coreform: $-8.7\%$ / $-10.7\%$). |
| **E02** | **Phase 3 (G3)** | E01 모델에 `CadEntitiesRegion(sets=["hole"], offset=0.1)` (기본 $p\cdot h/2^L = 0.15\text{ in}$), `level=2` ($0.3 \to 0.15 \to 0.075\text{ in}$), `two_to_one` 2:1 균형 THB 세분화 적용. (`wht_iga.grid.hierarchy` + `wht_iga.selection` 활용) | 동일 100점 `Probe-Stress-Profile` 상에서 $\sigma_1, \sigma_{\text{vM}}$ 추출. 균일 격자(E01) 대비 최소 DOF로 Roark $3.7261\times 10^{-6}\text{ psi}$ 수렴 개선 정량 확인. |
| **E03** | **Phase 5 (G5)** | 13개 솔리드 불리언 머지(S-8061 유리 + Alloy 52 핀 11개 + 304L VAR 헤더, MMTS), `example-problems/glass-to-metal-seal/data/` 온도 의존 $E(T), \alpha(T)$, J2 다선형 소성 CSV 그대로 로딩. 환형 장착면 `DISTRIBUTED_COUPLING` 고정. $450^\circ\text{C} \to 20^\circ\text{C} \to -55^\circ\text{C} \to 200^\circ\text{C}$ 3단계 `NLGEOM=ON` 열응력. 격자 $(0.5,0.5,0.5)\text{ mm}, p=2$. | 유리 $\sigma_1^+ = \max(\sigma_1, 0)$ 최대값/위치, 헤더·핀 `PEEQ` 최대값, $y=0$ 종단면 뷰컷, 유리-핀/유리-헤더 계면 경로, 최대 변위, + P0~P4 민감도 케이스 비교. Coreform 동일 CSV 3종(`convergence-history`, `performance`, `problem-size`) 출력 및 증분 수 비교(Coreform 15/15/15 대비 2배 이내). |
| **E04** | **Phase 5 (G5)** | Cook's membrane ($L=48, H_L=44, H_R=16, \Delta y=16, t=1.91516\text{ mm}$), 압축성 Neo-Hookean ($E=240.565\text{ MPa}, \nu=0.49999$), 좌측 `ENCASTRE`, 평면변형률 $u_z=0$, 우측 $P_y=6.25\text{ N}$, `NLGEOM=ON`. $U(p=2, C^0)/P_0(C^{-1})$ 정적 응축 C++ 요소. | 상단 우측 모서리 점 프로브의 비선형 하중 단계별 $u_y$ 이력 추출 → 최종 $u_y = 6.98946\text{ mm}$ 대비 상대 오차 $\le 1\%$. |
| **E05** | **Phase 5 (G5)** | Door-lock 3부품(Hinge, Link, Pin) 침체 $p=2, C^1, h=2\text{ mm}$, 패딩 $(2,2,2)$, 2종 탄소성 강재, `NLGEOM=ON`. Link 고정면 `PINNED`, Hinge $(0,0,5)\text{ mm}$ 강제 변위. 3쌍 표면층(`Triangle3D3`) 마찰 접촉($\mu=0.3$). | 구속면 `hold_link_surfaces` 전체 반력 합산 `RF3` 이력 추출 → 최종 $RF3 = -46.946\text{ kN}$ 대비 상대 오차 $\le 5\%$. |

---

## 6. 단계별 구현 로드맵 (기존 로드맵과의 정렬)

1. **즉시/Phase 2 마무리 직후**:
   - 누락된 `v01.json`, `v02.json`, `v03.json` 생성 및 `v04`–`v12`, `e01`–`e05` JSON에 `probes`/`sweeps`/`reference_csvs` 필드 보강
   - `wht_iga/verification/{__init__.py, schema.py, qoi.py, comparator.py, reporter.py}` 골격 구축 (기존 `post/`·`query/`·`selection/` 재사용)
2. **Phase 3 (G3)**:
   - `V02-local`, `V03-moment`, `V04`, `V05`, **`E01`**, **`E02`** 자동 실행기 및 Coreform 호환 CSV/리포트 출력 완성
3. **Phase 4 (G4)**:
   - `V09` (3종 CSV 포함), `V10` (26개 시계열 지표), `V11`, `V12` (`T2/T3/T4`) 자동화
4. **Phase 5 (G5)**:
   - `V07` (`global_l2_metrics.json`), `V08`, **`E03`** (`data/` CSV 로딩 + 3종 이력/성능 CSV + P0–P4), **`E04`**, **`E05`** 자동화
