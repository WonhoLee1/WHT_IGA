# V02 (Roark Plate-with-Hole) 비단조 수렴 · Crown NNLS Fallback · Level 5 Pardiso 분해 실패 원인 분석 및 해결 처방

- **작성일**: 2026-11-15
- **대상 요청서**: `WHT_IGA/docs/REQUEST_ANALYSIS_V02_CONVERGENCE.md`
- **분석 원칙**: 기존 소스 코드 무수정(Read-Only) 상태에서 코드 라인(`파일:줄`), 수식 전개, 로그(`docs/phase1/evidence/v02_run.log`) 및 재현 가능한 수치만을 근거로 원인 규명 및 최소 처방 제시.

---

## 요약 (Executive Summary)

1. **질문 1 (비단조 수렴의 1차 원인)**: `h=5`(Level 1, 오차 `8.237%`)와 `h=2.5`(Level 3, 오차 `6.549%`)에서만 오차가 급증하고 `h=3`(Level 2, `0.622%`)과 `h=1.5`(Level 4, `0.575%`)에서는 `< 0.65%`로 수렴하는 현상은 **(d) 배경 격자 노트 평면(`x = 0`, `y = ±10`)과 원공 정점(Crown, `(0, ±10)`)의 접선 정렬(Grid-Crown Tangency)** 이 촉발한 **(b) 미확장(`extension=None`) 소지지(small-support) 기저함수 불안정성 + (a) Crown 쐐기 셀의 `octree_fallback` + (c) `boundary_peak_stress`의 경계 외삽 미분 증폭**의 결합 현상이다. 실제로 체적 적분점 최대응력(`vol_peak`)은 Level 1에서 `4.302e-6`(오차 **`-0.36%`**), Level 3에서 `4.375e-6`(오차 **`+1.33%`**)로 이미 정해에 근접해 있으며, 오차의 대부분(`+8.6%p`, `+5.2%p`)은 `boundary_peak_stress`가 정점 쐐기(cusp) 끝단에서 불안정한 소지지 자유도의 미분을 직접 평가할 때 발생한다.
2. **질문 2 (`q=4` NNLS 모멘트 피팅 실패 원인)**: `quadrature.py`에서 `crown_fb > 0`이 발생하는 원인은 단일 파라미터 부족이 아니라 코드 내 **4가지 수치/기하 결함**의 중첩이다:
   - **(원인 2-A, 치명적)** `_cross_section_moments`(`quadrature.py:178-180`)에서 `np.round(p2 / tol).astype(int)`가 Windows 32비트 `int`(`max = 2,147,483,647`)를 오버플로하여 $|y| > 2.147\,\max(h)$인 모든 Crown 셀에서 `-2147483648`(`INT_MIN`)로 붕괴됨 (`v02_run.log:5-8`의 `RuntimeWarning`으로 100% 실증됨).
   - **(원인 2-B, 치명적)** `_cross_section_moments`(`quadrature.py:227-235`)에서 직사각형 경계 안쪽 판정 시 `eps = 1e-7 * max(h)`만큼 안쪽으로 이동해 판정하므로, 원공 정점 접선 셀($\delta y < \epsilon$)에서 상단 경계 변(`y = 10`, $m_y = +1$)이 누락되어 $N_{00} \approx -h_y h_z$의 부호 반전이 발생하고 `M[0,0,0] <= 0`(`quadrature.py:385`)으로 즉시 `octree_fallback`으로 빠짐.
   - **(원인 2-C)** `_axis_pw`(`quadrature.py:52-68`)가 셀 국소 좌표 $u \in [-1, 1]$ 대신 **전역 좌표 $(x, y, z)$의 단항식 거듭제곱**으로 $(q+1)=5$차 원시함수를 이항 전개하여 $|\beta_y|^5 = |-(2y_0/h_y + 1)|^5 \sim 10^4 \sim 1.4 \times 10^7$ 배의 부동소수점 상쇄 오차(Catastrophic Cancellation)를 유발함.
   - **(원인 2-D)** 얇은 쐐기(sliver) 절단 셀에서 전체 셀 $[blo, blo+h]$ 기준 텐서 가우스 후보점 중 $\Omega$ 내부에 남는 층(layer) 수가 $q+1=5$개 미만이어서 설계 행렬 $A \in \mathbb{R}^{125 \times N}$가 구조적으로 계수 부족($\text{rank}(A) < 125$)이 되며, `fit_moment_rule`(`quadrature.py:264-266`)의 균등 스트라이딩(`cand[::stride]`)이 세밀한 표면 후보점들을 추가로 탈락시킴.
3. **질문 3 (Extension 재활성화 조건)**: 현재의 올바른 경계조건(`bc_roark.roark_left_pinned()`) 하에서 `extension={"theta": 0.2}`를 재활성화하면 Level 1(`h=5`)과 Level 3(`h=2.5`)의 정점 쐐기 셀($\eta = 0.0868, 0.0419 < 0.2$)에만 걸치는 소지지 기저함수들이 내부 마스터 자유도에 구속되어 경계 응력 오버슈트가 해소된다. 단, `support.py:58, 182`의 분류 규칙이 "지지 영역 내 **단 하나의 셀이라도** $\eta_c \ge \theta$이면 `MASTER`"로 판정하므로, 인접 셀이 $\eta \approx 0.298$인 Crown 주변 기저까지 안정화하려면 `theta = 0.3` (또는 지지 전체 체적비 `ratios[b] >= theta_supp` 조건 병행) 및 절단 셀 외부 영역에 대한 FCM 유령 강성(`ghost_alpha`) 보강이 가장 안전하다.
4. **질문 4 (Level 5 Pardiso 분해 실패 원인)**: Level 5(`267x53x8`, $N_{\text{dof}} = 443,850$, $\text{nnz}(K) \approx 1.33 \times 10^8$)는 $p=2$ 3D 스플라인의 단면 프론트 폭 $F = 2 \times (55 \times 10 \times 3) = 3,300$으로 인해 **3D 직접법 분해 시 인수 비영 원소 수 $\text{nnz}(L+U) \sim 1.5\times 10^9 \sim 2.2\times 10^9$가 32비트 `MKL_INT` 한계($2^{31}-1 = 2.147\times 10^9$) 및 31.6 GB 물리 RAM(Kratos `ublas` $K$ + Eigen 복사본 + 16스레드 MKL 작업 버퍼)을 초과하여 발생하는 `MKL error = -2`(메모리/32비트 정수 부족)** 가 주원인이며, 동시에 `diag_v02_singular.py:103-113`에는 `phase 22` 예외 발생 시 `iparm(14..23)` 진단 정보를 덮어써 버리는 버그가 존재한다. 또한 물리적으로 V02는 $\nu = 0.0$인 평면응력/무포아송 문제이므로 해가 $z$에 완벽히 무관($\partial/\partial z \equiv 0$)하여 Level 5를 `(267, 53, 2)` 또는 `(267, 53, 1)`로 두면 **동일한 면내 해상도($h_x = h_y = 0.75$)와 동일한 이론해를 유지하면서 메모리와 계산 시간을 $1/10$ 이하로 줄여 즉시 해결**된다.

---

## 질문 1. 비단조 수렴의 1차 원인

### 1.1 결론
`h=5`(Level 1)와 `h=2.5`(Level 3)에서만 `boundary_peak_stress` 오차가 각각 `8.237%`, `6.549%`로 치솟는 1차 원인은 **`(d) 격자 노트 평면(x=0, y=±10)과 원공 정점(0, ±10)의 접선 정렬`이 만든 퇴화 쐐기(cusp) 셀에서 `(b) 미확장(extension=None) 소지지 기저함수`가 `(a) QUAD_FALLBACK`으로 부정확하게 적분된 채 `(c) boundary_peak_stress`에 의해 경계 끝단에서 미분 평가되기 때문**이다. 영향도 순위는 **`(d) 격자 접선 정렬 [촉발 조건] ➔ (b) 미확장 소지지 기저함수 [증폭 메커니즘] ➔ (c) 경계 미분 외삽 [관측 증폭] ➔ (a) Crown 적분 Fallback [국소 강성 오차]`** 순이다.

### 1.2 근거 (`파일:줄`, 수식, 재현 가능한 수치)

1. **짝수/4배수 격자(Level 0, 1, 3) vs 홀수 격자(Level 2, 4)의 기하학적 차이 (`scripts/v02_roark.py:29-32`)**:
   - 배경 영역은 `LO = (-100, -20, -0.5)`, `HI = (100, 20, 0.5)`이고 원공 반지름은 $R = 10$이므로 최대 응력이 발생하는 정점(Crown)은 $(x, y) = (0, \pm 10)$이다.
   - `LEVELS`(`v02_roark.py:32`)의 $(n_x, n_y)$를 대조하면 노트 평면과 정점 $(0, \pm 10)$의 위치 관계가 정확히 두 그룹으로 갈린다:

   | Level | `(nx, ny, nz)` | $h_x, h_y$ | $x=0$ 노트 평면 여부 | $y=\pm 10$ 노트 평면 여부 | Crown $(0, 10)$ 인접 최소 양수 $\eta$ | `vol_peak` 오차 | `boundary_peak` 오차 |
   | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
   | **0** | `(20, 4, 1)` | `10.0, 10.0` | **Yes** ($i=10$) | **Yes** ($j=1,3$) | `0.2146` | `-19.63%` | `10.552%` |
   | **1** | `(40, 8, 2)` | `5.0, 5.0` | **Yes** ($i=20$) | **Yes** ($j=2,6$) | **`0.0868`** (쐐기) | **`-0.36%`** | **`8.237%`** |
   | **2** | `(67, 13, 2)` | `2.985, 3.077` | **No** (홀수, 셀 중앙) | **No** ($y \in [7.69, 10.77]$) | **`0.2500`** (정칙) | `-1.96%` | **`0.622%`** |
   | **3** | `(80, 16, 4)` | `2.5, 2.5` | **Yes** ($i=40$) | **Yes** ($j=4,12$) | **`0.0419`** (쐐기) | **`+1.33%`** | **`6.549%`** |
   | **4** | `(133, 27, 4)` | `1.504, 1.481` | **No** (홀수, 셀 중앙) | **No** ($y \in [9.63, 11.11]$) | **`0.2500`** (정칙) | **`-0.71%`** | **`0.575%`** |

2. **정점 접선 정렬 시 쐐기(Cusp) 셀의 체적비 수식**:
   - $x = 0$과 $y = +10$이 동시에 격자 평면일 때(Level 1, Level 3), 정점 바로 아래 셀 $[0, h_x] \times [10 - h_y, 10]$에서 물리 영역 $\Omega$($x^2 + y^2 \ge 10^2$)의 $y$-두께는 다음과 같이 $x \to 0$에서 2차 함수로 소멸한다:
     $$\Delta y(x) = 10 - \sqrt{100 - x^2} = \frac{x^2}{20} + O(x^4)$$
   - 이 셀의 활성 체적비 $\eta$는 적분으로 정확히 계산된다:
     $$\eta(h_x, h_y) = \frac{1}{h_x h_y} \int_0^{h_x} \left(10 - \sqrt{100 - x^2}\right) dx$$
     - **Level 1 ($h_x = h_y = 5$)**: $\eta = \frac{1}{25}\left[50 - 50\left(\frac{\sqrt{3}}{4} + \frac{\pi}{6}\right)\right] = \mathbf{0.08682}$ ($8.68\% < \theta = 0.2$).
     - **Level 3 ($h_x = h_y = 2.5$)**: $\eta \approx \frac{2.5^3}{60 \times 6.25} = \mathbf{0.04193}$ ($4.19\% \ll \theta = 0.2$).
   - 반면 **Level 2 (`67x13x2`)** 와 **Level 4 (`133x27x4`)** 는 $n_x, n_y$가 모두 홀수여서 정점 $(0, 10)$이 단일 셀의 내부($y_9 = 7.692 < 10 < y_{10} = 10.769$ 및 $y_{20} = 9.630 < 10 < y_{21} = 11.111$)를 가로지르며, 최소 두께가 $\Delta y \ge 0.25 h_y$ 이상인 정칙 절단 셀이 된다.

3. **`vol_peak`와 `boundary_peak`의 괴리가 증명하는 국소 경계 오버슈트 (`docs/phase1/evidence/v02_run.log:1-14`, `solver/trimmed.py:346-420`)**:
   - `v02_run.log`에서 `vol_peak`(Kratos 적분점에서의 최대 주응력)와 `boundary_peak`(`trimmed.py:383`의 CAD 원공 표면 삼각형 Duffy 점에서의 최대 주응력)의 차이를 비교하면:
     - **Level 1 (`h=5`)**: `vol_peak = 4.302e-6` (**정해 `4.3175e-6` 대비 오차 `-0.36%`**) $\to$ `boundary_peak = 4.673141e-6` (**`+8.237%`**, 경계 평가 시 **`+8.60%p` 급증**).
     - **Level 2 (`h=3`)**: `vol_peak = 4.233e-6` (`-1.96%`) $\to$ `boundary_peak = 4.290654e-6` (**`-0.622%`**, 경계 평가 시 `+1.34%p`의 정상 물리적 증가).
     - **Level 3 (`h=2.5`)**: `vol_peak = 4.375e-6` (**`+1.33%`**) $\to$ `boundary_peak = 4.600234e-6` (**`+6.549%`**, 경계 평가 시 **`+5.22%p` 급증**).
     - **Level 4 (`h=1.5`)**: `vol_peak = 4.287e-6` (`-0.71%`) $\to$ `boundary_peak = 4.342333e-6` (**`+0.575%`**, 경계 평가 시 `+1.28%p`의 정상 물리적 증가).
   - 또한 `pipeline.py:173` (`props = ghost_props if (i, j, k) in out_set else vol_props`)에서 FCM 유령 강성(`ghost_alpha = 1e-8`)은 **오직 `OUT` 셀에만 부여되고 `CUT` 셀의 빈 영역($C \setminus \Omega$)에는 전혀 부여되지 않는다**. 따라서 `extension=None`일 때 정점 쐐기 셀($\eta = 0.0419$)에만 걸치는 스플라인 기저함수는 유령 강성조차 없이 미세한 `octree_fallback` 적분점들에만 의존하여 $(0, \pm 10)$ 경계에서 도함수($\nabla \mathbf{u}_h$) 오버슈트를 일으킨다.

### 1.3 검증 실험
- **명령**:
  ```powershell
  $env:WHT_V02_LEVELS="1,2,3"; D:\mamba\envs\whtiga\python.exe -c "
  import v02_roark as V
  from wht_iga.geometry import cad as G
  shape = G.read_step(V.STEP)
  for nc in [(39, 7, 2), (40, 8, 2), (41, 9, 2), (79, 15, 4), (80, 16, 4), (81, 17, 4)]:
      r = V.run_level(nc, shape, None, 200.0/nc[0])
  "
  ```
  *(실행 시 `PYTHONPATH="WHT_IGA;WHT_IGA/scripts"` 설정)*
- **기대 결과**: 동일한 $h \approx 5$ 및 $h \approx 2.5$ 해상도에서도 $(n_x, n_y)$를 홀수인 `(39, 7, 2)` / `(41, 9, 2)` 및 `(79, 15, 4)` / `(81, 17, 4)`로 바꾸어 $x=0, y=\pm 10$ 접선 정렬만 피해도 `err_pct`가 `8.24% / 6.55%`에서 즉시 `~1%`대로 하락한다.

### 1.4 불확실성
- **확인된 사실**: Level 1, 3에서만 $x=0$ 및 $y=\pm 10$이 격자 평면과 정확히 일치하여 $\eta = 0.0868, 0.0419$의 정점 쐐기 셀이 생성된다는 점과, `vol_peak`는 이미 오차 `0.36%`, `1.33%` 이내라는 점은 수식과 `v02_run.log` 수치로 100% 확정되었다.
- **추정(실험으로 분리 확인 필요)**: `8.237%`(Level 1) 및 `6.549%`(Level 3) 중 순수하게 `QUAD_FALLBACK` 적분 오차가 기여한 비중(약 1~2%p 추정)과 `extension=None` 소지지 자유도 진동이 기여한 비중(약 5~6%p 추정)의 정확한 분리 수치는 질문 2의 적분 수정과 질문 3의 Extension 활성화를 각각 단독 실행하여 계측해야 한다.

---

## 질문 2. NNLS 4차 모멘트(`q=4`) 피팅 실패 원인

### 2.1 결론
`cell_quadrature(..., q=4)`가 Crown 절단 셀에서 `octree_fallback`으로 빠지는 원인은 **`quadrature.py` 내부의 4가지 구체적 수치/알고리즘 버그**(① `_cross_section_moments`의 Windows 32비트 `int` 오버플로, ② 접선 경계에서 `eps` 안쪽 프로브로 인한 상단 직사각형 변 누락 및 `M[0,0,0] <= 0` 발생, ③ `_axis_pw`의 전역 좌표 5차 단항식 이항 전개에 따른 상쇄 오차, ④ 얇은 쐐기 영역의 후보점 층 부족 및 `fit_moment_rule`의 균등 스트라이딩)이다.

### 2.2 근거 (`파일:줄`, 수식, 재현 가능한 수치)

1. **[버그 2-A] `_cross_section_moments` 중복 제거 키의 32비트 정수 오버플로 (`wht_iga/geometry/quadrature.py:163, 178-180`, `docs/phase1/evidence/v02_run.log:5-8`)**:
   - `quadrature.py:163`에서 `tol = _COPLANAR_RTOL * max(h) = 1e-9 * max(h)`로 정의하고, `178-180`줄에서 다음과 같이 정수 튜플 키를 만든다:
     ```python
     key = (tuple(np.round(p2 / tol).astype(int)),
            tuple(np.round(r2 / tol).astype(int)),
            tuple(np.round(m / tol).astype(int)))
     ```
   - Windows 64비트(LLP64) 환경의 NumPy에서 `.astype(int)`의 기본 정수형이 32비트 부호 있는 정수(`int32`, 최댓값 $2^{31}-1 = 2,147,483,647$)로 동작할 때, 전역 좌표 $|y|$가 다음 임계값을 넘으면 오버플로가 발생한다:
     $$\frac{|y|}{10^{-9}\max(h)} > 2,147,483,647 \iff |y| > 2.147483647 \times \max(h)$$
   - 각 Level에서 원공 정점 $|y| \le 10$에 대한 오버플로 여부를 계산하면 `v02_run.log:5-8`의 경고 발생 시점과 소수점까지 일치한다:
     - **Level 0 ($\max(h)=10$)**: 임계 $|y| > 21.47$ $\implies$ 오버플로 없음.
     - **Level 1 ($\max(h)=5$)**: 임계 $|y| > 10.74$ $\implies$ $|y| \le 10$이므로 오버플로 없음.
     - **Level 2 ($\max(h)=3.077$)**: 임계 $|y| > 6.608$ $\implies$ **Crown 부근 $|y| \in [6.61, 10]$에서 최초로 오버플로 발생!** 실제로 `v02_run.log:5-8`에서 `[V02 h=3]` 직전에 정확히 `quadrature.py:178`과 `179`줄에서만 `RuntimeWarning: invalid value encountered in cast`가 출력되었다(`180`줄의 `m`은 $|m|\le 1 < 6.61$이므로 경고가 뜨지 않음).
     - x86_64 `cvttsd2si` 명령은 오버플로된 모든 부동소수점을 동일한 상수 `-2147483648`(`0x80000000`)로 변환하므로, $|y| > 2.147\max(h)$인 셀에서 서로 다른 $y$ 위치를 갖는 선분들이 동일한 `key`로 충돌하여 `seen` 집합(`quadrature.py:181`)에서 누락될 수 있다.

2. **[버그 2-B] 접선 Crown 셀($x_1 = 0, y_1 = 10$)에서 `eps` 프로브가 원공 내부로 넘어가며 $M_{000} \le 0$ 유발 (`wht_iga/geometry/quadrature.py:189-235, 385`)**:
   - 정점 바로 왼쪽 아래 셀 $[-h_x, 0] \times [10 - h_y, 10] \times [z_0, z_1]$에서 오른쪽 단면은 $x_1 = 0$, 상단 변은 $y_1 = 10$이다.
   - $x_1 = 0$ 단면에서 원공 테셀레이션의 수직 모선(또는 현)이 $y = 10 - \delta y$ ($\delta y < \epsilon = 10^{-7}\max(h)$)에 위치하면:
     1. 삼각형 루프(`quadrature.py:189-204`)는 $y \approx 10$에서 원공 안쪽을 향하는 법선 $m_y = -1$을 가진 수직 선분을 방출하여 $N_{00}^{\text{tri}} = -H_0(y_1)(z_1 - z_0) = -h_y h_z$를 더한다.
     2. 2차원 발산 정리가 성립하려면 직사각형 상단 변 $e=1$ ($y = y_1 = 10$, 외향 법선 `en[1] = [+1, 0]`)이 방출되어 $+H_0(y_1)(z_1 - z_0) = +h_y h_z$로 이를 상쇄해야 한다.
     3. 그러나 `quadrature.py:227`의 안쪽 프로브 점 `inward = ... - eps * en[1]`은 $y = 10 - \epsilon$에서 `winding_inside`를 평가한다. $\delta y < \epsilon$이면 이 프로브 점은 얇은 고체 층을 건너뛰어 **원공 내부(고체 외부, `winding_inside == False`)** 에 찍히므로 `if f.any():`(`234`줄)가 거짓이 되어 상단 변 $e=1$이 방출되지 않는다!
     4. 그 결과 $N_{00} \approx -h_y h_z$가 그대로 남아 `quadrature.py:127`(`M + Axi * Njk`)에서 $M_{000} \approx \text{vol}(\Omega_c) - h_x h_y h_z < 0$이 되며, `quadrature.py:385`(`if float(M[0, 0, 0]) <= 0: return fb()`)에 의해 **NNLS를 호출조차 하지 못하고 즉시 `octree_fallback`으로 직행**한다.

3. **[버그 2-C] `_axis_pw`의 전역 좌표 단항식 전개로 인한 부동소수점 상쇄 오차 (`wht_iga/geometry/quadrature.py:52-68, 103-105, 184-186`)**:
   - `_axis_pw(x0, hx, L, integrate)`는 $L_j\left(\frac{2(x-x_0)}{h_x}-1\right)$ 및 그 원시함수를 **전역 좌표 $x$의 거듭제곱** $\sum_{k=0}^{q+1} W_{j,k} x^k$로 전개한다(`al = 2.0 / hx`, `be = -(2.0 * x0 / hx + 1.0)`).
   - $q = 4$일 때 원시함수 $A_4(x), H_4(y)$는 5차 다항식이며, 이항 계수 전개 항의 크기는 $\frac{35}{8}\binom{4}{k} |\beta|^{5-k} \frac{h}{2}$에 비례한다.
   - Level 4($h_y = 1.481, y_0 = 9.63 \implies \beta_y = -14$)에서 $|\beta_y|^5 = 537,824$로 인해 Horner 평가(`_polyval`) 시 개별 항이 $\sim 6 \times 10^6 h_y$에 달해 약 $7$자리의 유효숫자가 상쇄되며, $x_0 \approx \pm 100$ 부근에서는 $|\beta_x| \approx 133 \implies |\beta_x|^5 \approx 4.2 \times 10^{10}$으로 **$10 \sim 11$자리의 정밀도가 손실**된다.
   - `_NNLS_TOL = 1e-12`(`quadrature.py:22`)는 모멘트 벡터 `M`이 최소 $10^{-13} \cdot \text{cell\_vol}$ 이하의 정확도를 가질 때만 달성 가능하므로, 전역 좌표 전개로 오염된 모멘트 벡터 $b$는 양수 가중치 원뿔 $\text{cone}(A)$ 밖으로 벗어나 NNLS가 수렴할 수 없다.

4. **[버그 2-D] 얇은 쐐기 셀의 후보점 층 부족(Rank Deficiency) 및 `fit_moment_rule`의 균등 스트라이딩 (`wht_iga/geometry/quadrature.py:264-266, 388-408`)**:
   - Level 3 Crown 셀 $[0, 2.5] \times [7.5, 10.0]$에서 고체 영역은 $y \in [\sqrt{100-x^2}, 10.0] \subset [9.6825, 10.0]$ (셀 높이의 상위 $12.7\%$ 이하)에만 존재한다.
   - 전체 셀 $[7.5, 10.0]$에 배치된 `base = _tensor_gauss(6)` 중 $y > 9.6825$인 점은 최상단 **단 1개 층**($y = 9.9156$)뿐이고, `depth = -1`의 `_tensor_gauss(12)`도 상위 **단 2개 층**($y = 9.7603, 9.9539$)뿐이다.
   - 게다가 `depth = 2, 3, -1`을 누적할 때 후보점 수가 `_NNLS_MAX_CANDIDATES = 2048`을 넘으면 `fit_moment_rule`(`quadrature.py:265-266`)이 `cand = cand[::stride]`로 **단순 균등 추출**하여 후반부에 추가된 고밀도 후보점의 절반 이상을 버린다.

### 2.3 최소 수정 코드 요지 (`wht_iga/geometry/quadrature.py`)
1. **셀 국소 좌표 평행이동 (`_cell_moments`, `100-127`줄)**:
   - `_cell_moments` 진입 시 `blo`를 원점으로 평행이동(`tris_loc = tris - blo`, `blo_loc = np.zeros(3)`)하여 `_axis_pw(0.0, h[d], L, ...)`를 호출하면 항상 $\beta_d = -1.0$이 되어 전역 좌표 거듭제곱 상쇄 오차(버그 2-C)가 **완전히 제거**된다.
2. **`_cross_section_moments` 정수 오버플로 및 접선 프로브 수정 (`163-235`줄)**:
   - `178-180`줄의 `key` 생성 시 `(p2 - [y0, z0]) / tol`을 사용하고 `.astype(np.int64)`로 캐스팅하여 32비트 오버플로(버그 2-A)를 원천 차단한다.
   - `205-217`줄에서 3D 삼각형 변 전체가 아니라 **$x = x_1$ 단면 교선 선분 `(a2, b2)`** 가 직사각형 변과 만나는 교점을 `params[e]`에 추가하도록 고치고, `227`줄의 `inward` 오프셋을 셀 내 절단 조각 크기에 맞춘 적응형 `eps_in = min(eps, 0.1 * piece_len)`으로 줄이거나 `M[0,0,0]`을 클리핑된 3D 다면체 발산 적분(6면 중 $x=x_1$ 면도 `clip_poly_to_box`로 닫힌 단면 다각형 triangulation 사용)으로 계산하여 접선 누락(버그 2-B)을 제거한다.
3. **절단 활성 영역 AABB 내부에 텐서 가우스 후보점 배치 (`cell_quadrature`, `388-408`줄)**:
   - 내부로 판정된 후보점들(및 클리핑된 표면 다각형 정점들)의 실제 바운딩 박스 `[alo, ahi]`를 구해 그 sub-box `[alo, ahi]` 내부에 `_tensor_gauss(q + 2)` 후보점을 생성하면, 아무리 얇은 Crown 쐐기 셀($\eta \sim 0.04$)이라도 두께 방향으로 최소 $q+2=6$개의 독립적인 가우스 층을 확보하여 $\text{rank}(A) = 125$가 보장된다.
   - 또한 `fit_moment_rule`(`265`줄)에서 `len(cand) > _NNLS_MAX_CANDIDATES`일 때 앞쪽(`depth=2`) 대신 **뒤쪽 고해상도 점들을 우선 보존**(`cand = cand[-_NNLS_MAX_CANDIDATES:]` 또는 QR 열 피벗 선택)하도록 수정한다.

### 2.4 검증 실험
- **명령**:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -c "
  import numpy as np, sys
  sys.path.insert(0, 'WHT_IGA'); sys.path.insert(0, 'WHT_IGA/scripts')
  import v02_roark as V
  from wht_iga.geometry import cad as G
  from wht_iga.reference.box import BoxSpline
  from wht_iga.solver import trimmed as T
  shape = G.read_step(V.STEP)
  for lv in (1, 2, 3, 4):
      nc = V.LEVELS[lv]
      box = BoxSpline(lo=V.LO, hi=V.HI, n_cells=nc, degree=2)
      h = float(min((np.asarray(V.HI) - np.asarray(V.LO)) / np.asarray(nc)))
      tess = G.tessellate(shape, deflection=float(h / 4))
      _, _, _, rec, diag = T.build_trimmed_table(tess.nodes, tess.tris, box, q=4, cache_dir='0')
      cfb = sum(1 for r in rec if r['n_qp'] > 125)
      print(f'Level {lv} {nc}: cut={diag[\"n_cut\"]} fb={diag[\"n_fallback\"]}')
  "
  ```
- **기대 결과**: 위 3가지 수정을 적용하고 캐시를 무효화(`cache_dir='0'`)한 뒤 실행하면, Level 1~4에서 `RuntimeWarning`이 사라지고 `fb`(`n_fallback`) 및 `crown_fb`가 `0`으로 감소한다.

### 2.5 불확실성
- **확인된 사실**: `quadrature.py:178-179`의 `int32` 오버플로(`v02_run.log:5-8`), `_axis_pw`의 전역 좌표 이항 전개 구조, 그리고 Level 3 Crown 셀 $[0, 2.5]\times[7.5, 10]$에서 `_tensor_gauss(6)` 및 `_tensor_gauss(12)`가 $y \ge 9.6825$ 영역에 각각 1개, 2개 층만 배치한다는 기하학적 사실은 100% 확정되었다.
- **추정**: 테셀레이션 현(chord)의 미세한 비평면성 때문에 `_NNLS_TOL = 1e-12`가 일부 극단적 마이크로 슬리버($\eta < 10^{-4}$)에서 여전히 엄격할 수 있으며, 이 경우 허용오차를 `max(_NNLS_TOL * cell_vol, 1e-10 * float(M[0, 0, 0]))` (활성 체적 기준 상대오차 $10^{-10}$ 병행)으로 두면 불필요한 fallback이 완전히 사라진다.

---

## 질문 3. 확장(Extension) 재활성화 조건

### 3.1 결론
현재의 `roark_left_pinned()` 경계조건에서 **`extension={"theta": 0.2}` (권장: `theta = 0.25 ~ 0.30`)를 재활성화하면 Level 1(`h=5`)과 Level 3(`h=2.5`)의 정점 쐐기 셀($\eta = 0.0868, 0.0419$)에 걸친 소지지 자유도가 내부 마스터 자유도로 외삽 구속되어 비단조 오버슈트가 억제된다.** 과거 G1 단계에서 나타났던 `-23%` 과강성(under-prediction)은 Extension 때문이 아니라 자유변/대칭면의 잘못된 성분별 페널티 구속(`bc_roark.py:5-8` 주석에 명시됨) 때문이었으므로, 현재 경계조건에서는 Extension을 켜는 것이 수학적·물리적으로 타당하다.

### 3.2 근거 (`파일:줄`, 수식, 재현 가능한 수치)

1. **`support.py`와 `extension.py`의 동작 구조 검증 (`wht_iga/stabilize/support.py:47-64, 182-190`, `wht_iga/stabilize/extension.py:61-98`)**:
   - `support.py:182`는 각 셀 $c$에 대해 `cell_large = eta >= theta - large_tol`을 판정하고, 기저함수 $b$의 지지 영역 $(p+1)^3 = 27$개 셀 중 **하나라도** `cell_large[c] == True`이면 `MASTER`(`labels == 2`), 하나도 없으면서 `cell_inside > 0`이면 `EXTENDED`(`labels == 1`, 슬레이브), 전부 `cell_inside == 0`이면 `INACTIVE`(`labels == 0`)로 분류한다.
   - `trimmed.py:254-264`(`apply_extension_ties`)는 `op.slave_ids`(`EXTENDED` 기저들)에 대해서만 Kratos `LinearMasterSlaveConstraint` ($u_{\text{slave}} = \sum_{m \in \text{master}} e_{s,m} u_m + 0$)를 생성한다.
   - `extension.py:61-77`(`extrapolation_weights_1d`)은 Marsden 항등식을 이용해 실제 노트 벡터 위에서 $p=2$ 다항식을 정확히 재현($\sum_m w_m = 1$)하므로, 균일 인장 및 선형/2차 변형률장을 인위적 잠김(locking) 없이 보존한다.

2. **Crown 정점 $(0, \pm 10)$에서 $\theta$ 임계값과 기저 분류의 정량적 관계**:
   - **Level 1 (`40x8x2`, $h=5$)**: 정점 아래 셀 $(19, 5), (20, 5)$의 체적비는 $\eta = 0.0868$이고, 그 옆 셀 $(18, 5), (21, 5)$의 체적비는 $\eta \approx 0.614$이다.
     - $x \in [-5, 5]$ 및 $y \in [-5, 10]$에 지지되는 정점 기저함수는 지지 내 최대 셀 체적비가 $\max \eta_c = 0.0868$이므로 **$\theta \in [0.10, 0.49]$ 전 구간에서 `EXTENDED`(슬레이브)로 분류되어 안정화**된다.
   - **Level 3 (`80x16x4`, $h=2.5$)**: 정점 주변 $y \in [7.5, 10.0]$ ($j=11$) 열의 셀 체적비는 다음과 같다:
     - $i = 39, 40$ ($x \in [-2.5, 2.5]$): $\eta = \mathbf{0.0419}$
     - $i = 38, 41$ ($x \in [-5.0, -2.5], [2.5, 5.0]$): $\eta = \mathbf{0.2980}$
     - $i = 37, 42$ ($x \in [-7.5, -5.0], [5.0, 7.5]$): $\eta = \mathbf{0.8190}$
     - 따라서 $y \in [2.5, 10.0]$ ($b=11$, 지지 셀 $j \in \{9, 10, 11\}$)이면서 $x \in [-5.0, 2.5]$ ($a=40$, 지지 셀 $i \in \{38, 39, 40\}$)인 기저함수 $B_{40, 11, c}$는 **전체 9개 $(x,y)$ 셀 중 6개가 $\eta=0$, 2개가 $\eta=0.0419$, 단 1개(`i=38, j=11`)가 $\eta=0.2980$** 이어서 전체 지지 체적비(`ratios[b]`)는 불과 **$4.24\%$** 에 불과하지만:
       - $\theta = 0.20$이면 단 하나의 코너 셀 $\eta = 0.2980 \ge 0.20$ 때문에 `MASTER`로 남아 구속되지 않는다!
       - **$\theta = 0.30$ (또는 $0.35$)으로 설정하면** $\max \eta_c = 0.2980 < 0.30$이 되어 이 기저함수까지 정확히 `EXTENDED`로 포섭된다!

3. **주의할 결합 지점 (`trimmed.py:254`의 `INACTIVE` 기저 처리)**:
   - `op.slave_ids`는 `EXTENDED` 기저만 포함하고 `INACTIVE`($\eta=0$, `OUT` 셀에만 걸치는 기저)는 포함하지 않는다(`extension.py:198`). `INACTIVE` 기저들은 `pipeline.py:173`의 `ghost_props`(`ghost_alpha = 1e-8`)에 의해서만 지지된다.

### 3.3 검증 실험
- **명령**:
  ```powershell
  D:\mamba\envs\whtiga\python.exe -c "
  import sys, numpy as np
  sys.path.insert(0, 'WHT_IGA'); sys.path.insert(0, 'WHT_IGA/scripts')
  import v02_roark as V
  from wht_iga.geometry import cad as G
  from wht_iga.reference.box import BoxSpline
  from wht_iga.solver import bc_roark, trimmed as T
  from wht_iga.solver.pipeline import run_linear
  from wht_iga.solver.validate import JobSpec
  shape = G.read_step(V.STEP)
  for theta in (0.1, 0.2, 0.3, 0.4):
      print(f'=== theta={theta} ===')
      for lv in (1, 2, 3, 4):
          nc = V.LEVELS[lv]
          box = BoxSpline(lo=V.LO, hi=V.HI, n_cells=nc, degree=2)
          h = float(min((np.asarray(V.HI) - np.asarray(V.LO)) / np.asarray(nc)))
          tess = G.tessellate(shape, deflection=float(h / 4))
          fid = T.hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids, radius=V.HOLE_R)
          qp, out, iss, rec, diag = T.build_trimmed_table(tess.nodes, tess.tris, box, q=4)
          spec = JobSpec(name=f'v02_t{theta}_l{lv}', lo=V.LO, hi=V.HI, n_cells=nc, degree=2,
                         E=1.0, nu=0.0, out_cells=out, **bc_roark.roark_left_pinned())
          res = run_linear(spec, geo_issues=iss, qp_table=qp, cell_records=rec, extension={'theta': theta})
          bp = T.boundary_peak_stress(box, res.U, 1.0, 0.0, tess.nodes, tess.tris, tess.face_ids, face_id=fid)['peak']
          err = abs(bp - V.REF_PEAK) / V.REF_PEAK * 100
          print(f'  L{lv} {nc}: peak={bp:.6e} err={err:.3f}% ties={res.tie_info}')
  "
  ```
- **기대 결과**:
  - 기존 캐시된 `qp_table`을 그대로 재사용하므로 4개 $\theta \times 4$개 Level 전체 스윕이 **약 60~90초 내**에 완료된다.
  - $\theta = 0.30$에서 Level 3(`80x16x4`)의 $B_{40,11,c}, B_{41,11,c}$ 기저가 `EXTENDED`로 편입되면서 Level 1, 2, 3, 4 전 구간에서 오차가 단조 감소한다.

### 3.4 불확실성
- **확인된 사실**: Level 3에서 정점 기저 $B_{40,11,c}$의 인접 셀 체적비가 $\eta = 0.2980$이어서 $\theta = 0.20$일 때는 `MASTER`로 남고 $\theta \ge 0.30$일 때 `EXTENDED`로 전환된다는 점은 기하 적분과 `support.py:58, 182` 코드로부터 100% 수학적으로 확정되었다.
- **추정**: `Kratos`의 `ResidualBasedBlockBuilderAndSolver`가 `LinearMasterSlaveConstraint`($T^T K T$)를 조립할 때 슬레이브 수가 많아지면 희소 행렬 대역폭이 약간 증가할 수 있으나, Level 1~4 규모($N_{\text{dof}} \le 70,470$)에서는 풀이 시간 차이가 1초 미만이다.

---

## 질문 4. Level 5 (`267x53x8`, $h=0.75$) Pardiso 분해 실패 원인

### 4.1 결론
Level 5(`267x53x8`, $N_{\text{dof}} = 443,850$)에서 `pardiso_lu`와 `pardiso_ldlt`가 모두 `"Decomposition failed!"`로 실패한 원인은 **① 3차원 $p=2$ B-스플라인 격자의 큰 단면 프론트 폭($F = 3,300$)으로 인해 직접법 분해 시 인수 행렬 비영 원소 수 $\text{nnz}(L+U) \approx 1.5\times 10^9 \sim 2.2\times 10^9$가 32비트 `MKL_INT` 상한($2^{31}-1$) 및 31.6 GB 시스템 메모리(Kratos `ublas` $K$ + Eigen 내부 복사본 + 16스레드 MKL 작업 버퍼)를 초과하는 `MKL error = -2`**, 또는 **② 미세 슬리버/유령 자유도 스케일링($10^{-8}$)에 의한 수치적 피벗 붕괴(`MKL error = -4`)** 이다. 특히 **현재 작성된 `scripts/diag_v02_singular.py:103-113`에는 `phase 22` 예외 발생 시 `iparm` 정보를 덮어써 버리고 Kratos C++ 전략 객체를 메모리에 중복 보유하는 결함**이 있어, 이를 먼저 수정해야 정확한 MKL 에러 코드를 읽을 수 있다.

### 4.2 근거 (`파일:줄`, 수식, 재현 가능한 수치)

1. **Level 4 vs Level 5의 3D $p=2$ 스플라인 직접법 메모리/Fill-in 스케일링**:
   - $p=2$ $C^1$ 3D B-스플라인 내부 제어점은 $(2p+1)^3 = 125$개 제어점($375$ DOFs/row)과 결합되어 1차 육면체 유한요소(27개 노드, $81$ DOFs/row)보다 행당 비영 원소가 **$4.63\times$ 조밀**하고, 단면 분리자(separator) 두께도 $p=2$개 층이다.
   - Level 4와 Level 5의 행렬 규모를 정량 비교하면:

   | 항목 | Level 4 (`133x27x4`, $p=2$) | Level 5 (`267x53x8`, $p=2$) | 증가 배율 |
   | :--- | :---: | :---: | :---: |
   | 제어점 수 $N_{cp}$ | $135 \times 29 \times 6 = 23,490$ | $269 \times 55 \times 10 = 147,950$ | **$6.30\times$** |
   | 총 자유도 $N_{\text{dof}}$ | $70,470$ | **$443,850$** | **$6.30\times$** |
   | 시스템 행렬 $\text{nnz}(K)$ | $\approx 1.8 \times 10^7$ (~0.22 GB) | $\approx \mathbf{1.33 \times 10^8}$ (~1.6 GB) | **$7.4\times$** |
   | $(y,z)$ 단면 프론트 폭 $F = 2 \cdot (N_y N_z \cdot 3)$ | $2 \times (29 \times 6 \times 3) = 1,044$ | $2 \times (55 \times 10 \times 3) = \mathbf{3,300}$ | **$3.16\times$** |
   | LU 인수 비영 원소 $\text{nnz}(L+U) \sim O(N_{\text{dof}} F)$ | $\approx 7.3 \times 10^7$ (~0.9 GB) | $\approx \mathbf{1.4 \times 10^9 \sim 2.2 \times 10^9}$ (**17~26 GB**) | **$\sim 20\times$** |

   - Kratos의 `EigenDirectSolver`(`v02_run.log:64`, `eigen_direct_solver.h:101`)가 실행될 때 프로세스 내에는 ① Kratos C++ `ublas::compressed_matrix` $K$(~1.8 GB), ② `Eigen::PardisoImpl` 내부의 `m_matrix` 복사본(~1.8 GB), ③ Python/테셀레이션/QP 기하 객체(~3 GB), ④ **16개 OpenMP 스레드(`v02_run.log:22`)에 대한 MKL PARDISO 프론트 작업 버퍼**가 동시에 상주한다.
   - 따라서 `pardiso_lu`에서는 $\text{nnz}(L+U)$가 32비트 부호 있는 정수(`MKL_INT` = `int32`, $2,147,483,647$)에 육박하거나 초과하여 `iparm[17]` 오버플로(`error = -2`)가 발생하고, `pardiso_ldlt`에서도 약 $10 \sim 14$ GB의 추가 연속 메모리 할당이 31.6 GB 시스템의 가용 커밋 한계(이미 RSS 14.4 GB 사용 중)에 걸려 `error = -2`를 반환한다.

2. **`scripts/diag_v02_singular.py`의 2가지 버그 (`scripts/diag_v02_singular.py:69, 97-115`)**:
   - **버그 4-A (`except` 블록의 `iparm` 유실)**: `diag_v02_singular.py:103-112`에서 `ps.set_phase(11)`이 성공한 뒤 `ps.set_phase(22)`에서 `PyPardisoError`가 발생하면 `except` 블록(`111-112`줄)이 `rep[f"pardiso_{tag}"] = {"error": str(e)}`만 기록하고 **`phase 11`이 이미 계산해 놓은 `iparm(15)`(분석 메모리 KB), `iparm(16)+iparm(17)`(인수 메모리 KB), `iparm(18)`(인수 `nnz`), `iparm(14)`(교란 피벗 수)를 모두 버린다!**
   - **버그 4-B (중복 메모리 점유)**: `69`줄에서 `K = system_matrix_csr(captured["s"]).tocsr()`로 추출한 뒤 `captured.clear()` / `del res`를 호출하지 않아 C++ `ublas` 행렬(~1.8 GB)과 SciPy CSR 행렬(~1.8 GB)과 `sp.triu(K)`(~0.9 GB)가 동시에 메모리를 점유한다.

3. **물리적 차원 축소를 통한 근본적 해결 ($\nu = 0.0$ 활용)**:
   - `v02_roark.py:64`에서 포아송 비는 **`nu = 0.0`** 이고 하중과 원공 형상은 모두 $z \in [-0.5, 0.5]$에 대해 균일하다.
   - $\nu = 0.0$인 3차원 선형 탄성체에서 평면 내 하중만 작용하면 두께 방향 변형률 $\varepsilon_{zz} = \frac{1}{E}[\sigma_{zz} - \nu(\sigma_{xx} + \sigma_{yy})] \equiv 0$이므로 **3차원 엄밀해 자체가 $z$에 완전히 무관($u_z \equiv 0$, $\partial \mathbf{u}/\partial z \equiv 0$)** 하다.
   - 따라서 Level 5의 두께 분할을 `nz = 8` 대신 `nz = 2`(또는 `nz = 1`), 즉 **`(267, 53, 2)`** 로 설정하면:
     - 면내 $(x, y)$ 격자 크기 $h_x = 0.749, h_y = 0.755$와 해의 정확도는 `(267, 53, 8)`과 **수학적으로 100% 동일**하면서,
     - 제어점 수는 $269 \times 55 \times 4 = 59,180$ ($N_{\text{dof}} = 177,540$, $2.5\times$ 감소), 단면 프론트 폭은 $F = 2 \times (55 \times 4 \times 3) = 1,320$ ($2.5\times$ 감소)이 되어 **인수 메모리(`nnz(L)`)가 $6.25\times$ 감소하고 연산량이 $15.6\times$ 감소하며, 적분해야 할 절단 셀 수도 $1/4$로 줄어들어 31.6 GB PC에서 약 60초 내에 여유 있게 풀린다!**
   - 만약 `(267, 53, 8)` 원본 격자를 그대로 풀어야 한다면 `JobSpec`에서 `linear_solver="pardiso_ldlt"`, `num_threads_mkl=4` (`pipeline.py:78-79`에서 이미 지원됨)로 스레드별 작업 버퍼를 줄여 실행하면 된다.

### 4.3 검증 실험
- **명령 1 (`diag_v02_singular.py` 보완 후 Level 5 진단 실행)**:
  `diag_v02_singular.py`의 `103-115`줄에서 `phase 11` 직후와 `except` 블록 내부에서도 `iparm(14..23)`을 `rep`에 저장하고, `69`줄 직후 `captured["s"].Clear(); captured.clear()`로 C++ 행렬 메모리를 해제한 뒤 실행:
  ```powershell
  $env:WHT_DIAG_LEVEL="5"; D:\mamba\envs\whtiga\python.exe WHT_IGA\scripts\diag_v02_singular.py
  ```
- **명령 2 (Level 5에서 `num_threads_mkl=4`, `pardiso_ldlt` 및 `(267, 53, 2)` 비교)**:
  - `diag_l5.json`에서 `empty_rows == 0`, `n_diag_le_0 == 0`인지 확인하고, `pardiso_lu` / `pardiso_ldlt`의 `iparm(16)+iparm(17)`(요구 메모리 KB)과 `error` 코드(`-2` vs `-4`)를 직접 확인한다.

### 4.4 불확실성
- **확인된 사실**: Level 5(`267x53x8`)에서 $N_{\text{dof}} = 443,850$, 단면 프론트 폭 $F = 3,300$이라는 점, `v02_run.log:22`에서 OpenMP 최대 스레드가 16개라는 점, `diag_v02_singular.py:111-112`가 예외 시 `iparm`을 기록하지 않는다는 점, 그리고 $\nu = 0$에서 해가 $z$에 무관하다는 점은 100% 확정된 사실이다.
- **추정(`diag_l5.json`으로 즉시 판별 가능)**: Level 5에서 Kratos `EigenPardisoLDLTSolver`가 던진 `"Decomposition failed!"`가 순수하게 `error = -2`(메모리 부족)인지, 아니면 버그 2-A/2-B로 깨진 슬리버 셀의 `error = -4`(수치적 영 피벗)가 동반되었는지는 위 명령 1의 `diag_l5.json` 출력으로 즉시 확정된다.

---

## 질문 5. V02 수용 기준 통과를 위한 최소 수정 처방 및 실행 순서

### 5.1 결론
V02 수용 기준(**연속 3개 이상 레벨 단조 수렴** 및 **최종 상대 오차 $\le 2\%$**)을 가장 적은 코드 변경으로 확실하게 달성하는 최소 처방은 다음 **3단계(우선순위 순)** 이다:
1. **[P1 — 필수, 솔버/모델 안정화]** `scripts/v02_roark.py`에서 `extension={"theta": 0.3}`(또는 `0.2`) 활성화 + `linear_solver="pardiso_ldlt"`, `num_threads_mkl=4` 지정 (필요 시 Level 5를 `(267, 53, 2)`로 설정하거나 Level 0~4의 5개 레벨에서 단조 수렴 달성).
2. **[P2 — 필수, 적분 무결성 수정]** `wht_iga/geometry/quadrature.py`에서 ① `_cell_moments`의 셀 국소 좌표 평행이동(`tris - blo`), ② `_cross_section_moments:178-180`의 `(p2 - [y0, z0]) / tol` 및 `np.int64` 캐스팅(`RuntimeWarning` 제거), ③ `cell_quadrature:388-408`의 절단 활성 바운딩 박스 `[alo, ahi]` 내 텐서 가우스 후보점 생성 및 `fit_moment_rule:265`의 후반부 고밀도 후보 보존.
3. **[P3 — 권장, FCM 유령 강성 일관성]** `wht_iga/solver/pipeline.py:166-175` 또는 `trimmed.py:186`에서 극소 체적비 절단 셀($\eta < 10^{-6}$)에도 최소 유령 강성(`ghost_alpha`)이 보장되도록 처리.

### 5.2 파일별 최소 수정 명세 (Surgical Diff Blueprint)

#### (1) `WHT_IGA/wht_iga/geometry/quadrature.py`
- **`_cell_moments` (`100-127`줄)**: 전역 좌표 거듭제곱 상쇄 오차 제거를 위해 셀 원점 `blo`로 평행이동:
  ```python
  def _cell_moments(tris: np.ndarray, blo: np.ndarray, h: np.ndarray,
                    q: int, whole: np.ndarray | None = None) -> np.ndarray:
      """Exact moments via recursive divergence theorem (tessellation-exact)."""
      blo = np.asarray(blo, float)
      tris_loc = np.asarray(tris, float) - blo
      whole_loc = tris_loc if whole is None else (np.asarray(whole, float) - blo)
      b0 = np.zeros(3)
      L = legendre_monomials(q)
      Px, Py, Pz = (_axis_pw(0.0, h[d], L, False) for d in range(3))
      Ax, Hy = _axis_pw(0.0, h[0], L, True), _axis_pw(0.0, h[1], L, True)
      x1, bhi = float(h[0]), np.asarray(h, float)
      # 이하 tris -> tris_loc, blo -> b0, whole -> whole_loc 동일
  ```
- **`_cross_section_moments` (`178-180`줄)**: Windows `int32` 오버플로 제거:
  ```python
      origin2 = np.array([y0, z0])
      key = (tuple(np.round((p2 - origin2) / tol).astype(np.int64)),
             tuple(np.round((r2 - origin2) / tol).astype(np.int64)),
             tuple(np.round(m * 1e6).astype(np.int64)))
  ```
- **`fit_moment_rule` (`265-266`줄)** 및 **`cell_quadrature` (`388-408`줄)**:
  - `fit_moment_rule`에서 후보점 수 초과 시 세밀한 표면/부영역 후보점들이 앞쪽 거친 점들보다 우선 보존되도록 `cand = cand[-_NNLS_MAX_CANDIDATES:]` (또는 고밀도 우선 추출) 적용.
  - `cell_quadrature`에서 `depth < 0` 단계 진입 시 이미 찾은 내부 점들(`base` + `surf`)의 바운딩 박스 `[alo, ahi]`에 맞춰 `(alo + Xg * (ahi - alo))` 후보점을 추가하여 얇은 Crown 쐐기 셀에서도 두께 방향 가우스 층 $q+2=6$개를 보장.

#### (2) `WHT_IGA/scripts/v02_roark.py`
- **`run_level` (`63-67`줄)**:
  - `extension` 파라미터를 환경변수 `WHT_V02_THETA` (기본값 `0.3` 또는 `0.2`, `"none"`이면 `None`)로 제어 가능하게 연결하고, `JobSpec`에 `linear_solver=os.environ.get("WHT_V02_SOLVER", "pardiso_ldlt")`, `num_threads_mkl=int(os.environ.get("WHT_MKL_THREADS", "4"))` 전달.

---

## 6. 통합자(Integrator) 실행 순서 제안

통합자는 아래 4단계를 순서대로 실행하여 **기존 캐시를 활용한 즉각적인 가설 검증 ➔ 적분 수정 검증 ➔ 전 레벨 수렴표 생성**을 수행할 수 있다.

| 단계 | 목적 | 실행 명령 (PowerShell, `D:\PythonCodeStudy\CAD-CAE-COPILOT` 기준) | 예상 소요 시간 | 성공 판정 기준 |
| :---: | :--- | :--- | :---: | :--- |
| **Step 1** | **Extension 단독 효과 즉시 검증** (`quadrature.py` 수정 전, 기존 `.cache/trimmed` 재사용) | `D:\mamba\envs\whtiga\python.exe`로 §3.3의 `theta in (0.2, 0.3)` Level 1~4 스윕 스크립트 실행 | **~60초** (캐시 히트) | Level 1(`h=5`) 및 Level 3(`h=2.5`)의 `err`가 `8.24% / 6.55%`에서 대폭 감소하고 Level 1→2→3→4 단조 수렴 경향 확인 |
| **Step 2** | **Level 5 Pardiso 실패 코드 확정** (`diag_v02_singular.py` 보완 후 실행) | `diag_v02_singular.py`에서 `captured["s"].Clear()` 및 `except` 내 `iparm(14..23)` 저장 추가 후:<br>`$env:WHT_DIAG_LEVEL="5"; D:\mamba\envs\whtiga\python.exe WHT_IGA\scripts\diag_v02_singular.py` | **~180초** (캐시 히트 시 조립+분석만 수행) | `docs/phase1/evidence/v02_diag/diag_l5.json` 생성, `empty_rows`, `diag_free_min`, `pardiso_ldlt.error` (`-2` vs `-4`) 및 `mem_kb_factor` 수치 확정 |
| **Step 3** | **`quadrature.py` 최소 수정 적용 및 Crown Fallback 해소 검증** | §5.2 (1)의 `quadrature.py` 국소 좌표화 + `int64` 키 + 슬리버 sub-box 후보점 적용 후:<br>`$env:WHT_V02_LEVELS="0,1,2,3,4"; D:\mamba\envs\whtiga\python.exe WHT_IGA\scripts\v02_roark.py` | **~150초** (새 해시로 절단 셀 적분 재계산) | `RuntimeWarning` 0건, Level 1~4에서 `crown_fb`가 0(또는 대폭 감소)이 되며 4개 이상 연속 레벨 단조 수렴 및 Level 4 `err < 1.0%` 달성 |
| **Step 4** | **Level 5 포함 전체 수렴표 완성** | `v02_roark.py`에서 `linear_solver="pardiso_ldlt"`, `num_threads_mkl=4` (및 필요 시 Level 5 `(267, 53, 2)`) 설정 후:<br>`$env:WHT_V02_LEVELS="0,1,2,3,4,5"; D:\mamba\envs\whtiga\python.exe WHT_IGA\scripts\v02_roark.py` | **~180초** (`267x53x2` 기준) / **~450초** (`267x53x8` 기준) | `docs/phase1/evidence/v02_roark.json`에 Level 0~5 전 구간 수렴 기록 완료, 최종 오차 `< 1.0%`, `PASS` 출력 |
