# P308 — 검증 프레임워크 (`wht_iga/verification/`)

## 1. 개요 및 구현 범위

Coreform IGA for Abaqus Verification Manual(`V01`–`V12`) 및 Example Problems(`E01`–`E05`) 총 17개 벤치마크를 표준화·자동화하는 검증 엔진 기술 문서다(T-P3.8).  
전체 검증 계층은 Kratos 의존성이 없는 순수 NumPy/SciPy 구조로 설계되어 솔버 백엔드와 독립적으로 QoI 추출, CSV 비교, 마크다운 리포트 생성을 수행한다.

### 관련 소스 파일 및 데이터
- `WHT_IGA/wht_iga/verification/__init__.py` — 공개 API 재수출
- `WHT_IGA/wht_iga/verification/schema.py` — 케이스 JSON 로더(`VerificationCaseSpec`) 및 기준 CSV 파서(`load_reference_csv`)
- `WHT_IGA/wht_iga/verification/catalog.py` — 17개 케이스(`V01`–`V12`, `E01`–`E05`) 카탈로그(`CASE_CATALOG`)
- `WHT_IGA/wht_iga/verification/qoi.py` — Kratos 비의존 물리 QoI 및 해석해 계산기
- `WHT_IGA/wht_iga/verification/comparator.py` — 스칼라 허용오차 비교 및 Coreform CSV 1:1 대조 엔진
- `WHT_IGA/wht_iga/verification/reporter.py` — Coreform 호환 CSV 라이터 및 마크다운 검증 리포터
- `WHT_IGA/tests/verification/data/` — `v01.json`–`v12.json` (12건), `e01.json`–`e05.json` (5건), 총 17개 JSON 명세

---

## 2. 스키마 및 17개 케이스 카탈로그 (`schema.py`, `catalog.py`)

### 2.1 `schema.py`
- `VerificationCaseSpec`: `case_id`, `title`, `category`, `status`, `geometry`, `materials`, `boundary_conditions`, `loads`, `mesh_sweep`, `acceptance`, `reference_values`, `provenance`, `meta`, `probes`, `sweeps`, `raw` 보관.
- `load_case_spec(case_id, data_dir=None)`: 레거시 플랫 JSON(`v04`–`v12`, `e01`–`e05`) 과 표준화된 `meta`/`model`/`sweeps`/`probes` 구조 JSON(`v01`–`v03`) 모두 투명하게 파싱.
- `list_available_cases(data_dir=None)`: 디렉터리 내 `*.json` 목록 반환.
- `load_reference_csv(path)`: `utf-8-sig` 로 CSV 를 읽어 셀 값을 `bool`/`int`/`float`/`str` 로 자동 타입 변환.

### 2.2 `catalog.py`
- 데이터클래스: `ProbeSpec` (종류: `"point"`, `"face_grid"`, `"line"`, `"volume_quadrature"`, `"history"`), `SweepSpec` (`degrees`, `continuities`, `mesh_sizes`, `workflows`, `load_cases`, `time_steps`, `padding`), `CaseCatalogEntry` (Phase Gate `"G1"`–`"G5"` 매핑 포함).
- `CASE_CATALOG`: 17개 전체 케이스(`V01`–`V12`, `E01`–`E05`) 등록.
- 헬퍼 함수: `get_catalog_entry(case_id)`, `list_cases_by_phase(phase_gate)`, `resolve_reference_csvs(case_id, user_guide_root=None)` (`WHT_IGA/research-cf-IGA/user-guide/` 하위 CSV 경로 해석).

---

## 3. QoI 추출기 및 해석해 (`qoi.py`)

| 심볼 | 대상 케이스 및 수학적 정의 |
|---|---|
| `principal_stresses(stress)` | Voigt `(6,)`/`(N,6)` 또는 `(3,3)` 대칭 응력 텐서 고유치 분해 → 내림차순 주응력 $(\sigma_1 \ge \sigma_2 \ge \sigma_3)$ 및 인장 최대주응력 $\sigma_1^+ = \max(\sigma_1, 0)$ (E03 유리 응력 등) |
| `line_profile_rel_l2(s, v_num, v_ref)` | 사다리꼴 적분 기반 1D 프로파일 상대 $L_2$ 오차 $\sqrt{\int \|v_{\text{num}} - v_{\text{ref}}\|^2 ds / \int \|v_{\text{ref}}\|^2 ds}$ |
| `spherical_projections(points, u, stress)` | 구좌표 반경 투영 $(u_r, \sigma_{rr}, \sigma_{tt})$, $\sigma_{tt} = \frac{1}{2}(\operatorname{tr}\boldsymbol{\sigma} - \sigma_{rr})$ (V07) |
| `lame_sphere_exact(points, a, b, p_i, p_o, E, nu)` | 내외압을 받는 두꺼운 구각(V07)의 정확한 Lamé 해석해(`LameSphereExact`) |
| `log_log_convergence_rate(h_values, errors)` | $\log(e) = \alpha \log(h) + c$ 최소제곱 수렴 차수 $\alpha$ 추정 |
| `time_window_rms_error(t, q_num, q_ref, t_min, t_max)` | 닫힌 시간 구간 $[t_{\min}, t_{\max}]$ 내 상대 RMS 오차 |
| `filter_participating_modes(freqs, participation_factors, threshold=0.01)` | $|\Gamma_i| \ge \text{threshold} \cdot \max|\Gamma|$ 참여 모드 필터링 (V09) |
| `analyze_dynamic_history(t, u, t_free_start, ref_freq=None)` | 자유진동 구간 상향 영교차/FFT 기본진동수 $f_1$, 정규화 RMS, 진폭 감쇠비, 주기당 샘플 수 계산 (V10) |

---

## 4. 비교 엔진 및 리포터 (`comparator.py`, `reporter.py`)

### 4.1 `comparator.py`
- `CheckItem`, `ComparisonReport`
- `compare_against_tolerances(case_id, measured, reference, tolerances, mode="rel"|"abs", units=None)`: 상대/절대 오차 판정 및 `CheckItem` 목록 생성.
- `compare_csv_rows(measured_rows, reference_rows, key_cols, value_cols, rel_tol, case_id=...)`: 복합 키(`key_cols`) 기준 Coreform CSV 와 1:1 행·열 비교.

### 4.2 `reporter.py`
- 표준 Coreform CSV 헤더 상수: `ROARK_QOI_HEADERS`, `ROARK_REFINEMENT_HEADERS`, `BIMETALLIC_TIE_HEADERS`, `MULTIPLANAR_RESPONSES_HEADERS`, `SPHERE_QOI_HEADERS`, `V09_EIGEN_HEADERS`, `E03_CONV_HEADERS`.
- `write_coreform_csv(path, headers, rows)`: UTF-8(BOM 없음) CSV 출력.
- `render_verification_markdown(report, metadata=None)`: 실행 메타데이터 및 비교 체크 표를 포함한 GitHub-flavored Markdown 리포트 생성.

---

## 5. 검증 및 테스트 커버리지 (직접 실행 관찰 결과)

| 테스트 파일 | 범위 | 실행 결과 | 비고 |
|---|---|---|---|
| `WHT_IGA/tests/unit/test_verification_catalog.py` | 17개 케이스 카탈로그, JSON 파싱(`V01`–`V12`, `E01`–`E05`), Gate/CSV 경로 해석 | **44 passed** (0.78s) | Kratos 무의존 단위 테스트 |
| `WHT_IGA/tests/unit/test_verification_core.py` | QoI 계산기(Lamé, 주응력, 수렴차수, 모드/동적 이력), 비교기, CSV/Markdown 리포터 | **9 passed** (0.42s) | Kratos 무의존 단위 테스트 |
