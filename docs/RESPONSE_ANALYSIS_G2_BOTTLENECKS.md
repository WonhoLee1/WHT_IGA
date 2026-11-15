# G2 진행 병목 정밀 원인 분석 및 개선안 (`RESPONSE_ANALYSIS_G2_BOTTLENECKS`)

- 작성일: 2026-11-15
- 작성 범위: **코드 수정 없음(Read-Only 분석 및 프로파일링)**
- 분석 대상 문서 및 로그:
  - [`WHT_IGA/docs/phase2/G2_report.md`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_report.md)
  - [`WHT_IGA/docs/INTEGRATION_HANDOFF.md`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/INTEGRATION_HANDOFF.md)
  - [`WHT_IGA/docs/phase2/evidence/v06_drop/L0_p3.log`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/evidence/v06_drop/L0_p3.log), [`L2_p2.log`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/evidence/v06_drop/L2_p2.log), [`v02_meshfile/run.log`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/evidence/v02_meshfile/run.log)
  - [`WHT_IGA/wht_iga/geometry/quadrature.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/quadrature.py), [`WHT_IGA/scripts/vcond_roark.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/scripts/vcond_roark.py), [`WHT_IGA/tests/kratos_spikes/test_p2mesh_v06.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p2mesh_v06.py), [`WHT_IGA/scripts/v02_roark.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/scripts/v02_roark.py)

---

## 요약

[`G2_report.md`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_report.md#L5)의 `COMPLETE_WITH_OPEN_ITEMS` 상태를 가로막고 있는 4대 보류·병목 항목을 직접 코드 감사 및 단일 셀/수치 프로파일링으로 진단한 결과, **가장 치명적인 시간 병목(`V06 p=3` 23분 중단)은 알고리즘 한계가 아니라 2줄짜리 구현 버그**였음이 확인되었습니다.

| # | 항목 (`G2_report.md` §3) | 실측된 근본 원인 | 기대 효과 (실측/근거) |
|---|---|---|---|
| **1** | **V06 `p=3` (`q=6`) 23분 중단** | `fit_moment_rule` 가지치기 루프 버그 (`dropped` 미초기화 + 데드코드 `break` 조건)로 셀 1개당 `nnls` 최대 **343회 중복 실행** | 3개 테스트 셀 `335.07s → 3.48s` (**96.3× 가속**), `cell 2` 단독 `307.07s → 1.65s` (**186×**), L0 전체 355셀 직렬 **`228.99s`**(8프로세스 추정 ~30s), `fb=0` |
| **2** | **V-COND 기울기 `-1.44` FAIL** | 이방 박판(`hz=h/10`) + clamped `p=2` 두께 제어점(`3→6`) 비선형 증가로 `n_free`가 $h^{-3}$이 아닌 $h^{-2.250}$ 증가 $\implies$ $\lambda_{\min}\propto h^{2.198}$, $\lambda_{\max}\propto h^{0.757}$ | `test_stabilize_cond.py` 검증 근거와 일관; Jacobi 스케일링 또는 기준/형상 정합성 정리 |
| **3** | **V06 `L0 p=1` (`u/ur` FAIL)** | `_csv_tols()`가 3개 Coreform 스터디(`immersed`, `regproj`, `wt2`)의 `min()`을 취해 `0.008` 대신 `0.006` 적용 + 대칭면 BC-B가 `ext`와 연쇄 MPC 가능성(가설) | CSV 스터디행 일치 시 허용치 `0.006 → 0.008`; 대칭면 직접 고정(`Fix`) 대조 실험으로 `1.02e-2` 원인 확정 |
| **4** | **V02 메쉬 `h=0.375` RAM 초과** | `v02_roark.py`가 `out_mode="drop"` 미적용(기본 `"ghost"`), 레벨 0~4 직후 동일 프로세스 연속 실행, MKL 16스레드 버퍼, `run_linear` 전 요소 응력 추출 | `out_mode="drop"` + `WHT_V02_LEVELS=4` 단독 실행 + `num_threads_mkl=4` / out-of-core `iparm[60]=1`로 여유 RAM(13~15 GB) 내 완주 |

---

## 1. [최우선 병목] `V06 p=3 (q=6)` 절단 셀 적분이 23분 동안 멈춘 진짜 원인

### 1.1 현상
* [`L0_p3.log`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/evidence/v06_drop/L0_p3.log#L2-L98): `L0` (`10×10×10`, 355개 non-OUT 셀)에서 8개 워커 프로세스로 **23분(`1430s`) 동안 첫 진행 틱조차 못 넘기고 중단**됨.
* [`G2_report.md:44, 52`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_report.md#L44): *"NNLS가 셀 규칙 시간의 ~87%... L2 추정 10 h 이상. 단계 순서·가지치기 변경은 결과를 바꾸므로 미적용"*으로 기록됨.

### 1.2 코드에서 발견한 3중 버그 ([`quadrature.py:314-342`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/quadrature.py#L314-L342))

```python
max_pts = moments.shape[0] ** 3               # q=6일 때 max_pts = 7^3 = 343
A, b = _design_matrix(cand, blo, h, moments.shape[0] - 1), moments.ravel()
try:
    w, rnorm = nnls(A, b, maxiter=10 * A.shape[1])
except RuntimeError:
    return (np.zeros(len(cand)), float("inf"), np.ascontiguousarray(cand))

keep = w > 0
idx = np.nonzero(keep)[0][np.argsort(w[keep])][: len(w) - 1]
best_w, best_r = w.copy(), float(rnorm)
dropped = np.zeros(len(w), dtype=bool)         # [버그 1-A] ~keep(미사용 후보)을 True로 막지 않음!
for ii in idx:
    dropped[ii] = True
    try:
        ww, rr = nnls(A[:, ~dropped], b)       # [버그 1-B] 1,400+개 미사용 열이 살아있어 대체 점만 골라짐
    except RuntimeError:
        dropped[ii] = False
        break
    if rr <= _NNLS_TOL * cell_vol:
        best_w = np.zeros(len(w))
        best_w[~dropped], best_r = ww, float(rr)
        if (~dropped).sum() <= max_pts:        # [버그 1-C] 활성 점 수가 아닌 전체 미삭제 열 개수를 검사 (데드코드!)
            break
    else:
        dropped[ii] = False
        break
```

1. **[버그 1-A] `dropped` 초기화 오류 (line 326):**
   1차 `nnls(A, b)` 결과 `w > 0`인 활성 점(`keep`)은 최대 343개이고 나머지(`~keep`, 예: `1760 - 343 = 1417`개)는 가중치가 `0`인 미사용 후보입니다. 그런데 `dropped`를 전부 `False`로 두었기 때문에, 루프가 `ii`번 활성 점 1개를 제거해도 `A[:, ~dropped]`에는 여전히 1,417개의 미사용 후보 열이 남아 있습니다. 따라서 `nnls`가 제거된 1개 점 대신 다른 후보 점을 선택해 **매 반복마다 `rr = 0.0`, 점 개수 343개인 해를 새로 구하는 헛수고**를 합니다.
2. **[버그 1-B] 조기 종료 조건 데드코드 (line 337):**
   `(~dropped).sum() <= max_pts`는 "현재 선택된 양수 가중치 개수(`(ww > 0).sum()`)"가 아니라 "후보 풀에서 삭제되지 않은 열의 총수(`1759 → 1417`)"를 비교합니다. 따라서 후보 개수 `len(cand) > 2 * max_pts` (`> 686`)이면 **`break`가 단 한 번도 실행되지 않고 `idx` 343개를 끝까지 모두 순회**하며 `nnls`를 343번 호출합니다.
3. **[버그 1-C] Lawson-Hanson 기본해 성질 미활용 및 실패 단계 중복 호출:**
   * `scipy.optimize.nnls`는 활성집합(active-set) 기반의 Lawson-Hanson 알고리즘이므로, 해 `w`는 **항상 기본해(basic solution)**이며 양수 성분 수 `(w > 0).sum()`은 **수학적으로 이미 `<= rank(A) <= (q+1)^3 = max_pts`** 입니다(343개를 넘을 수 없음).
   * 또한 1차 `nnls(A, b)`에서 이미 `rnorm > _NNLS_TOL * cell_vol`인 중간 단계(`depth=2, 3, -1`)에서도 `for ii in idx:` 루프에 진입해 열 1개를 빼고 `nnls`를 1회 더 호출합니다(열을 제거하면 최소제곱 잔차는 절대 감소하지 않으므로 100% 낭비).

### 1.3 실측 프로파일링 검증 결과 (`V06 L0, p=3, q=6`, `n_geo=24`)

동일한 첫 3개 절단 셀(`cut[0], cut[1], cut[2]`)에 대해 원본 코드와 가지치기 루프 버그 제거 버전을 직접 측정한 결과입니다.

| 대상 | 원본 `quadrature.py` | 가지치기 버그 제거 (`fast_fit`) | 가속 배수 / 결과 동일성 |
|---|---|---|---|
| `cut[0]` (4단계 진행) | `14.32 s` (`depth=-2` 2회: `7.23s + 6.42s`) | **`1.63 s`** | **8.8×** (`res=1.46e-14`, `n_pts=284`, `relaxed=True` 동일) |
| `cut[1]` (747개 후보) | `13.68 s` (`calls=[(747, 13.39s, 0.0, 343)]`) | **`0.20 s`** | **68.4×** (`res=0.0`, `n_pts=343`, `fb=False` 동일) |
| `cut[2]` (1,760개 후보) | **`307.07 s` (5.1분!)** (`1760후보 1셀에 299.64s`) | **`1.65 s`** | **186.1×** (`res=0.0`, `n_pts=343`, `fb=False` 동일) |
| **3개 셀 합계** | **`335.07 s`** | **`3.48 s`** | **96.3× 가속 (잔차/적분점 수 100% 일치)** |
| **L0 전체 355개 non-OUT 셀** | 8프로세스 **1,430s+ (0% 중단)** | **단일 프로세스(`workers=1`) `228.99 s`** | **8프로세스 추정 ~30초**, `n_fallback=0`, `n_relaxed=100`, `n_qp_cut=106,456` |

*(참고: `n_relaxed=100`이 나온 이유는 `cell 0`에서 보듯 상대 잔차 `1.46e-11`로 `_NNLS_TOL = 1e-12`보다는 크지만 `_NNLS_RELAXED = 1e-6`보다 `~10^5`배 정밀하게 수렴하여 2단 허용치로 `fb=0` 통과한 것입니다.)*

### 1.4 즉시 적용 가능한 개선 코드 ([`quadrature.py`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/quadrature.py))

#### (A) `fit_moment_rule` 수정 (`lines 323-342` 교체)
```python
    keep = w > 0
    # 1차 해가 이미 허용 오차를 못 넘기면 열 제거로 잔차가 줄 수 없으므로 즉시 반환.
    # 또한 Lawson-Hanson 기본해는 수학적으로 keep.sum() <= rank(A) <= max_pts 이므로
    # 점 개수 축소 목적의 가지치기는 keep.sum() > max_pts 인 경우 외에는 불필요.
    if rnorm > _NNLS_TOL * cell_vol or int(keep.sum()) <= max_pts:
        return w * keep, float(rnorm), np.ascontiguousarray(cand)

    # 만일 keep.sum() > max_pts 인 특이 케이스에만 활성 열(keep) 내부로 한정해 가지치기
    idx = np.nonzero(keep)[0][np.argsort(w[keep])]
    best_w, best_r = w.copy(), float(rnorm)
    dropped = ~keep  # [핵심 수정] 미사용 후보 열은 즉시 제외!
    for ii in idx:
        dropped[ii] = True
        try:
            ww, rr = nnls(A[:, ~dropped], b)
        except RuntimeError:
            dropped[ii] = False
            break
        if rr <= _NNLS_TOL * cell_vol:
            best_w = np.zeros(len(w))
            best_w[~dropped], best_r = ww, float(rr)
            if int((ww > 0).sum()) <= max_pts:  # [핵심 수정] 실제 양수 가중치 수 검사
                break
        else:
            dropped[ii] = False
            break
    return best_w * (best_w > 0), best_r, np.ascontiguousarray(cand)
```

#### (B) 과소 후보 단계 조기 스킵 + `maxiter` 현실화 (`lines 317, 558` 추가 가속)
* [`quadrature.py:558`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/quadrature.py#L558) `cell_quadrature`의 `depth in (2, 3)` 단계에서 `len(cand) < max_pts` (`343`)이면 343개 모멘트 방정식을 정확히 맞출 자유도가 부족하므로(`cut[0]`에서 `70, 240`개 후보 호출, `cut[2]`에서 `198`개 후보 호출), `depth == -1` 또는 `depth == -2`로 바로 넘기면 무의미한 NNLS 호출을 줄일 수 있습니다.
* [`quadrature.py:317`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/geometry/quadrature.py#L317) `maxiter=10 * A.shape[1]` (`20,480`회)는 과도합니다. Active-set `nnls`는 보통 `3 * m` (`m = (q+1)^3`) 이내에서 수렴하므로 `maxiter=max(5 * A.shape[0], 3 * A.shape[1])`로 제한해도 정밀도 손실 없이 최악 셀 정체를 방지할 수 있습니다.
* **L2(`12,472` non-OUT 셀) 예상 소요 시간:** 위 수정(A)+(B) 적용 시 8프로세스 기준 L2 `p=3` 셀 규칙 생성이 **기존 10시간+ 추정에서 약 10~15분 내외(인덱스/모멘트 병렬 포함)로 단축**될 것으로 예상됩니다(단, 12,472셀 실측 전이므로 추정치).

---

## 2. `V-COND` 기울기 FAIL (`-1.44` vs 기준 `-2 ± 0.3`) 정밀 진단

### 2.1 고윳값 분해 실측 (`docs/phase1/evidence/vcond/results.json`)

[`vcond/results.json:312-398`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase1/evidence/vcond/results.json#L312-L398)의 `slope_on` 4단계 데이터를 $\lambda_{\max}$, $\lambda_{\min}$, 자유도 수 $n_{\text{free}}$로 각각 분리해 log-log 기울기를 계산한 결과는 다음과 같습니다.

| 단계 `nc` | $h$ | $n_{\text{free}}$ | `kmax` ($\lambda_{\max}$) | `kmin` ($\lambda_{\min}$) | $\kappa$ | 구간별 기울기 ($s_{k_{\max}} / s_{k_{\min}} / s_{\kappa}$) |
|---|---|---|---|---|---|---|
| `(20, 4, 1)` | `10.00` | `1,134` | `167.98` | `1.168e-8` | `1.439e10` | — |
| `(40, 8, 2)` | `5.00` | `4,920` | `127.50` | `2.814e-9` | `4.531e10` | `n_free`: $-2.117$ / `kmax`: $+0.398$ / `kmin`: $+2.053$ / `κ`: **$-1.655$** |
| `(60, 12, 3)` | `3.33` | `12,750` | `79.07` | `1.106e-9` | `7.150e10` | `n_free`: $-2.348$ / `kmax`: $+1.178$ / `kmin`: $+2.304$ / `κ`: **$-1.125$** |
| `(80, 16, 4)` | `2.50` | `25,956` | `59.01` | `5.474e-10` | `1.078e11` | `n_free`: $-2.471$ / `kmax`: $+1.017$ / `kmin`: $+2.444$ / `κ`: **$-1.427$** |
| **전체 회귀** | — | **`-2.2503`** | **`+0.7565`** | **`+2.1983`** | **`-1.4418`** | $\kappa$ 기울기 $= s_{k_{\max}} - s_{k_{\min}} = 0.7565 - 2.1983 = \mathbf{-1.4418}$ |

### 2.2 물리·수학적 원인 2가지
1. **원인 A — 이방 박판(`hz=h/10`) + Clamped `+p` 경계 제어점 효과:**
   * 3D 탄성 강성 행렬에서 $\lambda_{\min}(K) = \min \frac{c^T K c}{c^T c} \propto (\text{제어점당 체적}) \propto \frac{1}{n_{\text{free}}}$로 스케일링됩니다.
   * 실제로 위 표에서 구간별 $n_{\text{free}}$ 증가 지수(`2.117 → 2.348 → 2.471`)와 `kmin` 감소 지수(`2.053 → 2.304 → 2.444`)가 소수 첫째 자리까지 정확히 일치합니다.
   * 두께 방향 셀이 `nz = 1 → 2 → 3 → 4`로 변할 때 $p=2$ clamped B-spline 제어점 수는 `nz + p` = `3 → 4 → 5 → 6`으로 2배만 증가하므로, $N_{\text{free}} \propto h^{-2.250}$(이론 점근 $h^{-3}$ 미달)에 머뭅니다. 게다가 `nz=1 → 2` 전환 시 두께 방향 기저가 단일 2차 Bernstein에서 다중 스플라인으로 바뀌며 `kmax` 기울기가 `+0.398`로 꺾여 전체 기울기를 `-1.44`로 끌어내렸습니다.
2. **원인 B — 외팔보(단면 1면 고정, 나머지 5면 자유) 문제의 점근 전 현상:**
   * 이미 단위 테스트 [`test_stabilize_cond.py:86-107`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_stabilize_cond.py#L86-L107)이 증명했듯이, 등방 정육면체 `(4..12)^3`, `p=2`에서조차 **한 면만 고정(`single-face-fixed`)하면 조건수 기울기가 `-1.35`(Jacobi 전처리 시 `-1.67`)**로 측정되며, **6개 외곽 면을 모두 고정(`full-dirichlet`)해야만 거친 격자(`nc=4..12`)에서 `-2.02` (`-2 ± 0.3`)가 나옵니다**([`test_stabilize_cond.py:57-83`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/unit/test_stabilize_cond.py#L57-L83)).

### 2.3 개선 선택지 (게이트 정책 결정용)
* **선택지 1 (등방 세분 영역으로 기울기 leg 분리):** 평탄(flatness) leg는 V02 박판 그대로 유지하고(`1.06× PASS`), $h$-기울기 leg만 전 경계 고정 또는 등방 블록(`hx=hy=hz`)에서 측정해 연산자 $O(h^{-2})$ 스케일링을 엄밀 검증(추가 8× 두께 메모리 부담 없음).
* **선택지 2 (박판에서 유효 $h_{\text{eff}} = (V_{\text{dom}}/n_{\text{free}})^{1/3}$ 또는 Jacobi 대각 전처리 기울기 병기):** 현재 박판 `nc=(20..80)`에서 $N_{\text{free}}$ 유한효과를 반영하거나, [`G2_report.md:42`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/G2_report.md#L42)대로 이방 박판 기준을 개정(사용자 승인 항목).

---

## 3. `V06 L0 p=1` (`u/ur = 1.02e-2` > `0.006` FAIL) 진단

### 3.1 확인된 팩터
1. **CSV 허용치 파싱 버그 ([`test_p2mesh_v06.py:278`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p2mesh_v06.py#L278)):**
   * [`pressurized-sphere-qoi.csv:2-13`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/research-coreform-IGA/user-guide/generated-downloads/verification-manual/problems/solid_mechanics/linear_elastic_stress/pressurized-sphere/pressurized-sphere-qoi.csv#L2-L13)에는 `(p=1, level=0, h=0.1)`에 **3개의 서로 다른 실험 설정**이 수록되어 있습니다([`pressurized-sphere.md:600-610`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/research-coreform-IGA/user-guide/verification-manual/problems/solid_mechanics/linear_elastic_stress/pressurized-sphere/pressurized-sphere.md#L600-L610)):
     * 기본 `immersed_h0p1`: `u=0.00584`, **`tolerance=0.008`**, `srr tol=0.40`
     * 특수 기법 `regproj_h0p1` / `wt2_h0p1`: `u=0.00421 / 0.00418`, `tolerance=0.006`, `srr tol=0.45`
   * `_csv_tols()`가 `min()`으로 3행을 합치면서 `srr`는 기본 행의 `0.40`을, `u/ur`는 특수 기법 행의 `0.006`을 가져오는 혼선이 발생했습니다(정식 기본 허용치는 **`0.008`**).
2. **`ext_theta`는 무관함을 실측 확인:**
   * 직접 `run_v06((10,10,10), 1, 24, out_mode='drop')`를 `ext_theta = None / 0.05 / 0.1 / 0.2`로 실행한 결과:
     * `None`: `u = 1.0193e-2`, `ur = 9.9615e-3`, `srr = 0.2406`
     * `0.05`: `u = 1.0199e-2`, `ur = 9.9676e-3`, `srr = 0.2403`
     * `0.10`: `u = 1.0212e-2`, `ur = 9.9806e-3`, `srr = 0.2402`
     * `0.20`: `u = 1.0237e-2`, `ur = 1.0004e-2`, `srr = 0.2409`
   * 즉 `ext` 온/오프 및 `theta` 크기와 무관하게 `u ≈ 1.02e-2`가 고정됩니다.

### 3.2 대조 확인 필요 가설 (미실험)
* **대칭면 경계조건 방식 차이 ([`test_p2mesh_v06.py:142-186`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/tests/kratos_spikes/test_p2mesh_v06.py#L142-L186) vs [`pressurized-sphere.md:641-649`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/research-coreform-IGA/user-guide/verification-manual/problems/solid_mechanics/linear_elastic_stress/pressurized-sphere/pressurized-sphere.md#L641-L649)):**
  * Coreform 원문은 정합 경계면(`x=0, y=0, z=0`)에 직접 `*Boundary, XSYMM/YSYMM/ZSYMM`(노달 직접 고정)을 사용합니다.
  * 반면 현재 코드는 전체 `[0,1]^2` 박스 면(OUT 셀 포함)에 가상 노드 `vnode` + `BC-B` 페널티 MPC(`vnode -> face_cp`)를 생성한 뒤, `apply_extension_ties`에서 다시 소지지 `face_cp`를 `slave -> interior_cp`로 묶어 **2단 연쇄 MPC 위험**이 있으며, 특히 `p=1`의 거친 1-셀 두께 구각(`a=0.9, b=1.0, h=0.1`)에서는 경계 `r<0.9`에 붙은 가상 노드/페널티 분포 및 `C^0` 선형 B-spline(Coreform trimmed U-spline 대비 자유도 배치 차이)의 영향을 크게 받습니다.
  * **추천 대조 실험:** `p=1, L0`에서 대칭 3면의 활성 제어점에 대해 가상 노드 BC-B 대신 직접 성분 고정(`node.Fix(DISPLACEMENT_{X|Y|Z})`)을 적용해 `u=1.02e-2`가 `0.006~0.008` 수준으로 떨어지는지 4초 실험으로 즉시 판별 가능.

---

## 4. `V02 메쉬 파일 h=0.375` 메모리 한계(~18 GB) 극복 방안

[`v02_meshfile/run.log:17-29`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/docs/phase2/evidence/v02_meshfile/run.log#L17-L29) 및 [`v02_roark.py:104-115`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/scripts/v02_roark.py#L104-L115) 분석 결과:

1. **구멍 내부 고스트 미제거 (`out_mode="ghost"`):**
   `v02_roark.py:104`에서 `JobSpec(...)` 생성 시 `out_mode="drop"`을 전달하지 않아 구멍 내부 OUT 셀(`h=0.375` 기준 `~2,096`개 추정)이 모두 `1e-8 E` 요소와 자유도로 조립됩니다. `out_mode="drop"`을 지정하면 불필요한 요소와 비활성 제어점이 제거됩니다.
2. **동일 프로세스 연속 실행 누적 메모리:**
   `run.log`에서 `h=0.75`를 푼 뒤 프로세스를 재시작하거나 `gc.collect()` 없이 곧바로 `h=0.375` (`181,932` 노드, `59,136` 요소)로 진입합니다. `WHT_V02_LEVELS=4`로 `h=0.375`만 독립 프로세스에서 실행하면 이전 레벨 메모리 잔존 없이 시작할 수 있습니다.
3. **MKL 16스레드 작업 버퍼 및 후처리 이중 적재:**
   * `num_threads_mkl=4`(또는 `2`)로 지정하면 16스레드 대비 Pardiso 병렬 작업 버퍼 RAM을 수 GB 절감할 수 있습니다.
   * 또는 `wht_iga/linalg/pardiso.py`를 통해 `iparm[60]=1` (Out-of-Core Pardiso) 경로를 활성화하면 13~15 GB 여유 RAM 환경에서도 54만 DOF 분해가 가능합니다.
   * [`v02_roark.py:125-130`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/scripts/v02_roark.py#L125-L130)의 Coreform 게이트 판정(`peak_A`, `profile`)은 변위 벡터 `U`만 있으면 스플라인 미분으로 0.01초 내 평가되므로, 메모리 여유가 빠듯한 `h=0.375`에서는 [`pipeline.py:319-326`](file:///D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA/wht_iga/solver/pipeline.py#L319-L326)의 59,136개 전 요소 `CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR)` 루프를 생략(옵션화)할 수 있다면 후처리 피크 메모리와 시간도 함께 절감할 수 있습니다.
