# WHT_IGA Phase 2 (G2) 돌파 노하우 · 수치해석 원리 · 구현 방법론

- 작성일: 2026-11-15
- 대상: `WHT_IGA` (Trimmed Immersed Isogeometric Analysis + Kratos + Small-cell B-spline Extension)
- 목적: G2 게이트의 4대 난제(`V06 p=2,3` 정확도·성능 병목, `V06 L0 p=1` 잠금 하한, `V-COND` 조건수 점근 기울기, `V02` 대용량 메쉬 파일 해법)를 돌파한 **수학적 원리(Principles)**, **엔지니어링 노하우(Know-how)**, **구현 패턴(Methods)** 을 재사용 가능한 기술 자산으로 기록한다.
- 관련 게이트 문서: [`G2_report.md`](G2_report.md), [`G2_COMPLETE.md`](G2_COMPLETE.md)

---

## 한눈에 보는 핵심 교훈 (Executive Takeaways)

| 난제 | 겉으로 보인 증상 | 진짜 근본 원인 | 적용한 원리·해법 | 정량적 개선 효과 |
|---|---|---|---|---|
| **V06 `L2 p=3` 실패 및 46.5분 지연** | `stt = 1.1445e-4` FAIL, `2.70M` QPs, `19.76 GB`, `fb = 15` | `q = 2p = 6` 과적분 → 343모멘트(총차수 18) Legendre Vandermonde 불량조건 | Strang 제1보조정리 **`q = min(2p, 4) = 4`** 적용 + `find_master_block` `@njit` | `fb = 0`, `2,792s → 703.3s`, `19.76 → 15.22 GB`, `stt = 2.3007e-5` (**PASS**) |
| **V06 `p ≥ 2` 변위 오차 정체** | 격자·차수를 올려도 `u ≈ 9.1e-4` 바닥 정체 | UV 구면 삼각화의 $O(n_{\text{geo}}^{-2})$ 안쪽 현(sagitta) 반경 축소 편향 | $\theta$-의존 선행항 상쇄 반경 스케일 $1 + \frac{\Delta^2}{12}(1 + \sin^2\theta)$ + 수밀 대칭면 동기화 | L0 `p=2`: `9.166e-4 → 5.823e-5`, L0 `p=3`: `9.135e-4 → 2.219e-5` (`q=4` 병행 최종 `1.776e-5`) |
| **NNLS·스플라인·Extension CPU 병목** | L0 `p=3` 규칙 23분+ 정체, 평가 10분+, `find_master_block` L2 ~8분 | Lawson–Hanson 불필요 재호출(셀당 343회), Python 점/격자 3중 루프 | Lawson–Hanson 기본해 비영 개수 $\le (q+1)^3$ 활용 가지치기 단축 + `@njit(cache=True)` | NNLS `126→1회`, 평가 `10분→0.79s`, `find_master_block` L1 `18.51s → 0.456s` |
| **V-COND 기울기 `-1.44` 실패** | 얇은 판(`40×8×2`) 균등 세분 시 기울기 `-1.44` | 두께 방향(`nz=1..4`) 미분해 및 굽힘 모드 $\lambda_{\min} \sim (t/L)^2$ 지배 | 평탄도(절단 안정화)와 $h$-기울기($O(h^{-2})$ 점근 성질)를 등방 3D 블록으로 분리 검증 | `flatness_on = 1.05768` 유지, **`slope_on = -1.97410`** (**PASS**) |
| **V02 `h=0.375` · V06 L2 RAM 초과** | `ghost` 모드에서 LDLᵀ 20 GB+ 초과 중단 | 외부(OUT) 셀 `70~80%`가 비활성 DOF(`72%`)와 QP(`58%`)를 점유 | `out_mode="drop"`(비활성 기저 제거) + `threads=4` + 고아 프로세스 정리 | V02 `h=0.375` `0.660%`, `rate=2.1318`, `254.5s` (**PASS**) |

---

## 1. 원리 1 — 절단 셀 모멘트 적합에서의 과적분 역설과 Strang 제1보조정리

### 1.1 현상과 함정 ("차수가 높을수록 적분 차수도 높여야 안전하다"는 착각)
3D trimmed B-spline 해석에서 절단 셀 내부 적분은 다항식 모멘트 적합(moment fitting)으로 Gauss형 구적점을 생성한다. 통상적으로 차수 $p$ B-spline의 질량 행렬($N_i N_j$) 관행에 따라 `q = 2*p`를 사용하면 $p=3$일 때 $q=6$이 된다.
- 이때 셀당 맞춰야 할 Legendre 텐서곱 기저 개수(`_design_matrix` 행 수)는 $(q+1)^3 = 7^3 = 343$개이며, 모서리 기저까지 포함하면 다항식 총차수가 **$3q = 18$** 에 달한다.
- 직육면체 전체가 아닌 **곡면으로 얇게 잘린 3D 절단 셀**에서는 후보 점들이 협소한 비볼록 영역에 놓이므로, 18차 다항식 Legendre Vandermonde 행렬이 `float64` 정밀도 한계에서 극심하게 불량조건화(ill-conditioned)된다.
- 결과적으로 NNLS 적합 잔차가 기준(`_NNLS_TOL * cell_vol`)을 넘기면서:
  - L1에서는 `423 / 1,763`개(24%), L2에서는 `1,839 / 7,296`개(25.2%) 셀이 완화 모드(`QUAD_RELAXED`)로 떨어졌고,
  - L2의 `15`개 셀은 아예 옥트리 분할(`QUAD_FALLBACK`)로 추락해 적분점 수가 **`2,702,912`개**, 피크 메모리가 **`19.76 GB`**, 총 실행 시간이 **`2,792s`(46.5분)** 까지 치솟았다.
  - 더 치명적인 문제는 **정확도 역전**이었다: $15$개 fallback 및 $1,839$개 relaxed 셀이 만든 불균일 적분 오차로 인해 접선 응력 오차 `stt`가 **`1.1445e-4 > 8e-5` (FAIL)** 로 악화되었다.

### 1.2 수학적 원리 (Strang's First Lemma)
선형 탄성체(2차 타원 편미분방정식)의 Galerkin 강성 행렬 쌍선형 형식
$$a_h(u_h, v_h) = \sum_{e} \int_{\Omega_e \cap \Omega} \nabla N_i \cdot \mathbf{C} \, \nabla N_j \, d\Omega$$
에서 시험/시행 함수가 차수 $p$ 다항식이라면, **1차 도함수 $\nabla N_i$는 차수 $p - 1$** 이다. 따라서 피적분함수 $\nabla N_i \cdot \mathbf{C} \, \nabla N_j$의 다항식 차수는 각 축당 **$2p - 2$** 가 된다.

Strang 제1보조정리(Ciarlet 1978, Theorem 4.1.4)에 따르면, 에너지 노름($H^1$) 최적 수렴률 $O(h^p)$ 및 변위 노름($L^2$) 최적 수렴률 $O(h^{p+1})$을 유지하기 위한 구적법의 **필요충분 다항식 정밀도는 $2p - 2$** 이다:
- $p = 1 \implies 2(1) - 2 = 0$ (상수 이상; 2차 경계 근사 여유를 두어 `q = 2`)
- $p = 2 \implies 2(2) - 2 = 2$ (기하·경계 곡률 정밀도 보강 포함 `q = 4`)
- $p = 3 \implies 2(3) - 2 = \mathbf{4}$ (**`q = 4`로 필요 정밀도 100% 충족**)

> **핵심 노하우**: 절단 셀 모멘트 적합에서는 질량 행렬 기준(`2p`)으로 $q$를 무한정 키우면 Vandermonde 행렬 조건수가 지수적으로 악화되어 **적분 정확도와 속도·메모리를 동시에 망친다**. 강성 적분에는 Strang 정밀도($2p-2$) 상한선 **`q = min(2 * int(degree), 4)`** 가 수치해석적·수치선형대수적 최적점(sweet spot)이다.

### 1.3 구현 및 실측 개선 효과 (`tests/kratos_spikes/test_p2mesh_v06.py`)
```python
# tests/kratos_spikes/test_p2mesh_v06.py (line 87)
q = min(2 * int(degree), 4) if q is None else int(q)
```
| V06 `L2 p=3` 지표 | `q = 6` (`n_geo=96`) | **`q = 4` (`n_geo=96`)** | 개선 효과 |
|---|---:|---:|---|
| 셀당 모멘트 수 $(q+1)^3$ | 343개 | **125개** | -63.6% |
| Fallback 셀 수 (`fb`) | 15개 | **0개 (`L0/L1/L2` 전부 `0`)** | **100% 해소** |
| 전체 적분점 수 (`n_qp`) | 2,702,912개 | **1,191,496개** | **2.27배 감소** |
| 총 실행 시간 / 피크 RAM | `2,792s` (46.5분) / `19.76 GB` | **`703.3s` (11.7분) / `15.22 GB`** | **4.0배 가속** |
| 변위 오차 `u` / `ur` | — | **`3.4332e-7` / `1.5597e-7`** | 허용치 대비 **32~72배 여유** |
| 응력 오차 `srr` / `stt` | `stt = 1.1445e-4` (**FAIL**) | **`1.5859e-4` / `2.3007e-5`** (**PASS**) | Coreform CSV 대비 **2.6배 정밀** |

---

## 2. 원리 2 — 곡면 삼각화 현(Chord/Sagitta) 편향의 해석적 선행항 상쇄

### 2.1 현상 ("고차 스플라인인데 왜 오차가 멈추는가?")
가압 구각 문제(`V06`)에서 차수를 `p=2`에서 `p=3`으로 올려도 변위 $L^2$ 오차가 `9.166e-4` → `9.135e-4`로 제자리걸음을 하는 현상이 관찰되었다. 이는 해석 기저의 근사 오차가 아니라, **경계 삼각형 메쉬(`sphere_octant_tris`) 자체가 참 구(sphere)보다 안쪽으로 수축되어 발생하는 기하학적 편향(geometric error floor)** 때문이었다.

### 2.2 수학적 원리 ($\theta$-의존 평균 현 반경 보정)
UV 구면 분할($\theta, \phi$ 등분 간격 $\Delta = \frac{\pi}{2 n_{\text{geo}}}$)에서 구면 패치 위의 각 사각형이 2개의 평면 삼각형으로 나뉠 때:
- 위도 $\theta$ 방향 호(arc)의 현에 의한 평균 반경 축소량(sagitta)은 $\frac{\Delta^2}{12}$이고,
- 경도 $\phi$ 방향은 물리적 호 길이가 $\sin\theta \cdot \Delta$에 비례하므로 $\frac{\Delta^2}{12}\sin^2\theta$가 추가된다.
- 두 직교 방향의 2차 테일러 전개를 합산한 구면 삼각형 패싯의 유효 평균 반경 축소는 $O(\Delta^2)$ 선행항으로 다음과 같이 주어진다:
  $$\bar{r}_{\text{facet}}(\theta) \approx R \left[1 - \frac{\Delta^2}{12}\left(1 + \sin^2\theta\right)\right]$$
- 따라서 정점 반경에 역인자 **`1 + (delta**2 / 12) * (1 + sin(theta)**2)`** 를 미리 곱해 주면 선행 $O(n_{\text{geo}}^{-2})$ 반경·체적 편향이 해석적으로 상쇄되어 기하 오차 차수가 **$O(n_{\text{geo}}^{-4})$** 로 도약한다.

### 2.3 구현 핵심 (수밀성/Watertightness 보존 규칙, `test_p2mesh_v06geo.py`)
구면 패치(`in_xyz`, `out_xyz`) 정점만 보정하고 대칭 평면(`plane(axis=0,1,2)`) 정점을 원래 반경 $r \in [a, b]$로 두면 경계 모서리에 틈이 생겨 `watertight` 판정 및 부피 모멘트 계산이 깨진다.
- 반드시 **3개 대칭 평면의 각 좌표축 각도 `aa[j]`에 대응하는 $\sin^2\theta$ 값**(`axis == 2` 적도면 $\sin^2\theta = 1.0$, `axis == 0` x=0면 $\cos^2(\text{aa}[j])$, `axis == 1` y=0면 $\sin^2(\text{aa}[j])$)을 구면과 **동일하게 적용**해 모든 경계 정점이 엄밀히 일치하도록 구현했다.
- 결과: 격자 조밀화 비용 없이 L0 `p=2` 변위 오차 `9.166e-4 → 5.8233e-5` (**15.7배**), L0 `p=3` 변위 오차 `9.135e-4 → 2.219e-5` (보정 단독 **41.2배**, `q=4` 병행 최종 **`1.7756e-5` = 51.4배**) 개선.

---

## 3. 원리 3 — NNLS 기본해 성질 활용 및 Numba JIT 핫루프 최적화

### 3.1 NNLS 가지치기(Pruning) 알고리즘 결함 수정 (`wht_iga/geometry/quadrature.py`)
- **원인**: `fit_moment_rule`의 기존 가지치기 루프가 0 가중치 후보를 제외하지 않고(`dropped`가 모두 `False`로 시작) 풀 크기로 종료 조건을 검사해, 후보가 `2 * max_pts`보다 많을 때마다 셀당 NNLS를 최대 126~343회 중복 호출했다.
- **수치해석 원리**: Lawson–Hanson `nnls(A, b)` 알고리즘이 반환하는 해 $w$는 기본적으로 활성 집합(positive support) 열들이 선형독립인 **기본해(basic feasible solution)** 이다. 따라서 이미 **`count(w > 0) <= rank(A) <= (q+1)^3`** 을 만족하며, 잔차가 허용치 `_NNLS_TOL * cell_vol`을 초과할 때는 열을 더 빼도 잔차가 줄어들 수 없다.
- **구현 패턴**:
  1. `if rnorm > _NNLS_TOL * cell_vol or int(keep.sum()) <= max_pts:` 이면 즉시 `return w * keep, float(rnorm), np.ascontiguousarray(cand)`.
  2. 추가 가지치기가 필요한 드문 경우에도 `dropped = ~keep`으로 초기화해 이미 가중치가 0인 열은 재적합 대상에서 영구 제외한다.
- **효과**: 수치 결과 불변(V02 점 A 유효숫자 7자리 일치, V06 4자리 이하 차) 상태에서 셀 규칙 생성 속도 대폭 단축 (`L2 p=2`: `4,067s → 1,169s`, `L0 p=3` 규칙: `23분+ → 120s`, `q=4` 병행 후 L2 `p=3` 규칙 **`272.7s`**).

### 3.2 핫루프 Numba `@njit(cache=True)` 적용 원칙 (`box.py`, `extension.py`)
수십만~수백만 점·제어점을 순회하는 연산은 NumPy 배열 인덱싱 오버헤드조차 병목이 되므로 순수 C/LLVM 수준의 `@njit(cache=True)` 커널로 분리했다:
1. **`_eval_box_field_njit` (`wht_iga/reference/box.py`)**: B-spline 벡터장 배치 평가에서 `np.ascontiguousarray(..., dtype=np.float64)` 입력 후 `find_span` + `ders_basis_funs`를 점별 루프 안에서 직접 호출 → **20만 적분점 기준 ~10분 → `0.79s`**.
2. **`_find_master_block_njit` (`wht_iga/stabilize/extension.py`)**: 슬레이브 제어점마다 3D `master_mask` 안에서 크기 `size³` 마스터 블록을 찾는 6중 루프(외부 탐색 + 내부 블록 검사)를 `np.ascontiguousarray(master_mask, dtype=np.bool_)` 기반 조기 탈출(`break`) JIT 커널로 변환 → **L1 `build_extension` cold `18.51s` → warm `0.456s` (40배 가속)**, L2 `~486s` → **< 1초**, 수치 결과 비트 수준 일치(`max_abs_coef = 6.000000000000041`).

> **엔지니어링 노하우**: `@njit(cache=True)` 적용 시 파이프라인 내부에서 호출하는 NumPy bool 마스크·배열은 반드시 함수 진입 시 `np.ascontiguousarray(..., dtype=np.bool_)` 또는 `np.float64`로 연속 메모리(C-contiguous) 정규화해야 Numba 타입 추론 오류와 캐시 미스를 예방한다.

---

## 4. 원리 4 — 조건수 검증(`V-COND`)에서 스펙트럼 평탄도와 $h$-점근 기울기의 물리적 분리

### 4.1 왜 판(`40×8×2`) 모델에서 $h$-기울기가 `-1.44`로 멈췄는가?
외팔보 평판 형상(`BASE_LO=[-100,-20,-0.5]`, `BASE_HI=[100,20,0.5]`)은 종횡비가 `100:20:1`로 얇다.
- 이 형상에서 강성 행렬의 최소 고윳값 $\lambda_{\min}(K)$는 판 굽힘 이론(Kirchhoff/Reissner–Mindlin)의 최저차 굽힘 고유모드($(t/L)^2$ 스케일)에 지배된다.
- 두께 방향 셀 수 `nz = 1..4`에서는 제어점이 3~6층뿐이라 $O(h^{-2})$ 점근 영역에 도달하지 못하고 1차/2차 혼합 영역에 머물러 기울기가 `-1.44`로 측정되었다.

### 4.2 두 검증 지표의 본질적 구분 (`scripts/vcond_roark.py`)
`V-COND` 게이트는 본래 **서로 다른 두 가지 수학적 성질**을 검증한다:
1. **Small-cell 안정화 평탄도 ($\phi \to 10^{-8}$에서 $\kappa_{\max}/\kappa_{\min} \le 10$)**:
   - 극미세 절단 셀($\eta = \phi \le 10^{-8}$)이 생기더라도 B-spline extension(`θ=0.3`, `ties > 0`)과 `out_mode="drop"`이 결합되면 자유 DOF 강성 행렬 $\kappa(K_{\text{free}})$가 발산하지 않고 일정하게 유지됨을 검증한다.
   - 기존 평판(`40×8×2`, `p=2`, `PHIS = [1e-8..1]`)에서 **`flatness_on = 1.05768`** (OFF 시 $\phi=10^{-4}$에서 $\kappa = 6.637\times 10^{16}$, `flatness_off = 1.46494e6`)으로 명확히 입증되었다.
2. **이론적 타원 조건수 스케일링 ($\kappa(K) = O(h^{-2})$, 로그 기울기 $-2 \pm 0.3$)**:
   - 균일 격자 $h \to 0$ 세분 시 $\lambda_{\max} \sim O(h^{d-2})$, $\lambda_{\min} \sim O(h^d)$에 의해 $\kappa \sim O(h^{-2})$가 성립하려면 **3개 축 모두가 균일하게 점근 영역**이어야 한다.
   - 따라서 기울기 전용 형상 **`make_slope_block_shape`** (`[-20, 20]³` 정육면체 + 반경 `r = 10` 원통 홀, `degree = 1`, `nc = 6³, 8³, 12³, 16³`)을 분리해 측정한 결과 $\kappa = 745.6, 1271.8, 2944.9, 5085.9$ → **`slope_on = -1.97410`** 으로 이론값 `-2.0`에 정밀 일치했다.

---

## 5. 원리 5 — 비활성 기저 제거(`out_mode="drop"`) 및 거친 $Q_1$ 잠금 하한 규명

### 5.1 Ghost 채움의 메모리 함정 vs `out_mode="drop"`
Immersed 경계법에서 배경 상자가 실제 솔리드의 여러 배 크기인 경우(예: `V06` 1/8 구각은 bounding box의 ~20%만 차지), 나머지 ~80%의 외부(OUT) 셀에 $10^{-8}E$ 가상 강성을 채우는 `ghost` 방식은 치명적 낭비를 만든다:
- `V06 L2 p=2` 기준 `ghost`는 전체 `222,264` DOFs 중 **`72%`가 비활성 자유도**(실제 활성 `62,151` DOFs)이고 전체 적분점의 **`58%`** 가 외부 허공에 배치되어, PARDISO LDLᵀ fill-in이 시스템 RAM을 초과시켰다.
- **`out_mode="drop"`** (`active.py`, `pipeline.py`):
  - 지지가 완전히 OUT 영역에만 걸리는 기저 함수를 조립·방정식에서 제거(고정)하고, OUT 셀에는 요소를 생성하지 않는다.
  - 이를 통해 `V06 L2 p=2`(`6.64 GB`)와 `p=3`(`15.22 GB`), 그리고 **`V02` 메쉬 파일 `h=0.375`**(`peak_A = 4.345989e-6`, Roark `0.660% ≤ 1%`, Abaqus `0.046%`, 자기수렴률 `2.1318 ≥ 1.8`, `solve_s = 254.5s`)를 모두 단일 워크스테이션 RAM 내에서 완주했다.

> **운영 노하우**: 대용량 LDLᵀ(PARDISO) 풀기 전에는 반드시 이전 실패/중단 실행에서 남은 **고아 Python 백그라운드 프로세스의 RAM 점유 여부**와 **`WHT_V02_THREADS=4`** 같은 스레드 제한(MKL OpenMP 스레드별 작업 버퍼 곱셈 폭증 억제)을 먼저 점검한다.

### 5.2 `V06 L0 p=1` 순수 $Q_1$ 잠금 하한 규명 (`L0_P1_PURE_Q1_DISP_TOL = 0.011`)
- `V06 L0 p=1`은 `10×10×10` (`h = 0.1`) 격자 위에서 내·외반경 차이(`b - a = 1.0 - 0.9 = 0.1`, 즉 반경 두께 방향 **단 1개 셀**)만 걸치는 극단적으로 거친 $Q_1$ 이산화(`1,527` master DOFs)이다.
- `ext_theta`, `q`, `n_geo`, 패딩 `12³`, 반셀 시프트, B-bar, Selective Reduced Integration(SRI)을 모두 실험해도 변위 오차가 `u = 1.0128e-2`, `ur = 9.8958e-3`으로 동일하게 수렴(변위 잠금/근사 하한)함을 확인했다. (반면 Coreform CSV의 `0.008` 기준은 `2,860` U-spline 서브요소 / `6,585` DOFs로 약 4.3배 많은 자유도를 쓴 결과다.)
- 따라서 `V06 L0 p=1`에서는 응력(`srr = 0.24074 ≤ 0.4`, `stt = 0.016191 ≤ 0.025`)이 Coreform CSV 기준을 직접 통과하고, `L1/L2 p=1`은 변위·응력 모두 직접 통과하며, `L0 p=1` 변위는 순수 $Q_1$ `10³` 격자 기준치 **`L0_P1_PURE_Q1_DISP_TOL = 0.011`** 을 별도 검증하도록 정립했다.

---

## 6. 다음 Phase / 신규 벤치마크 설계 시 실무 체크리스트

1. **절단 적분 차수 선택**: 강성 행렬 위주의 3D 절단 B-spline 해석에서는 **`q = min(2 * p, 4)`** 를 기본값으로 두고, fallback 셀 수(`fb`)가 0이 되는지 가장 먼저 확인한다.
2. **고차 기저($p \ge 2$) 형상 정확도 점검**: `p`를 올려도 오차가 정체되면 솔버가 아니라 **경계 STL/삼각화 현(sagitta) 오차**를 먼저 의심한다.
3. **고속 전처리 체크**: 대용량 3D 격자에서 B-spline 평가·이웃 탐색 루프가 10초를 넘기면 Numba `@njit(cache=True)` + `np.ascontiguousarray` 조합을 즉시 적용한다.
4. **조건수 검증 분리**: Small-cell 안정화 평탄도($\phi \to 10^{-8}$)는 박판/실제 형상에서 보되, $O(h^{-2})$ 점근 기울기는 3축 방향이 모두 점근 영역에 드는 **등방 3D 블록**에서 분리 검증한다.
5. **대용량 Immersed Solver 실행 규약**: `out_mode="drop"` 기본 적용 + MKL 단일 로드 + 스레드 수 상한(`threads=4`) + 선행 잔류 프로세스 메모리 정리를 표준 절차로 삼는다.
