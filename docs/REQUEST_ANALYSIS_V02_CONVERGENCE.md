# 분석 요청: V02 Roark 수렴 비단조 + 절단 셀 적분 규칙 대체(fallback) — 근본 원인과 처방

작성: 통합 담당, 2026-11-15. 이 문서 하나로 작업할 수 있게 썼다. 코드는 읽기만 하고, 결론은 근거(파일:줄, 수식, 재현 가능한 수치)로 제시해 달라.

## 0. 프로젝트 맥락 (최소)

- 저장소: `D:\PythonCodeStudy\CAD-CAE-COPILOT\WHT_IGA` (Kratos 10.4.3 + Python immersed IGA). 계획서 `WHT_IGA/research-cf-IGA/master-plan_v2.1.md` (§7.4 절단 적분, §7.5 extension, §10.1 V02).
- 경로 R1: 배경 직교 B-spline 격자(p=2), 셀마다 Kratos 적분점 기하 + `SmallDisplacementElement3D4N`, 절단 셀은 Python에서 만든 적분 규칙(`qp_table`)을 넘긴다. OUT 셀은 ghost 강성 α=1e-8.
- 실행 환경: `D:\mamba\envs\whtiga\python.exe`. Kratos 소스(`D:\PythonCodeStudy\kratos-src`)·빌드(`E:\kratos-build`) 수정 금지, 재빌드 금지. skip/xfail로 실패 숨기기 금지.

## 1. 문제 (V02)

Roark 17.1-7a: 200×40×1 판, 중심 구멍 r=10, E=1, ν=0, 왼쪽 면 강고정, +x 면 인장 압력 1e-6. 기준 peak(구멍 가장자리 최대 주응력) 4.3175e-6. G2 합격선: 최세 격자 peak 오차 ≤1%, 관측 수렴률 ≥1.8, 4레벨 이상.

스크립트 `scripts/v02_roark.py` (레벨 `LEVELS`, 면내 h `INPLANE`), BC는 `wht_iga/solver/bc_roark.py::roark_left_pinned`. 절단 적분 `q = 2p = 4`. extension OFF.

### 1.1 실측 (2026-11-15 00:00–00:40, `docs/phase1/evidence/v02_run.log`)

| 레벨 | n_cells | 면내 h | peak (경계 probe) | vol peak (체적 QP) | 오차 | 잘린 셀 | fallback 셀 | 구멍 주변(2r 이내) fallback | OUT 셀 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 20×4×1 | 10 | 3.861924e-6 | 3.470e-6 | 10.552% | 80 | 0 | 0 | 0 |
| 1 | 40×8×2 | 5 | 4.673141e-6 | 4.302e-6 | 8.237% | 626 | 16 | 16 | 8 |
| 2 | 67×13×2 | ≈3 | 4.290654e-6 | 4.233e-6 | 0.622% | 1667 | 16 | 16 | 42 |
| 3 | 80×16×4 | 2.5 | 4.600234e-6 | 4.375e-6 | 6.549% | 3298 | 48 | 48 | 128 |
| 4 | 133×27×4 | 1.5 | 4.342333e-6 | 4.287e-6 | 0.575% | 8521 | 32 | 32 | 468 |
| 5 | 267×53×8 | 0.75 | — | — | — | — | — | — | — |

- 레벨 5는 `pardiso_lu`에서 "Decomposition failed!". **`pardiso_ldlt`로 재실행해도 같은 실패**(전체 719 s, 프로세스 메모리 최대 14.4 GB / 시스템 31.6 GB). 메모리 한계 전이므로 특이성(또는 Pardiso 내부 오류) 쪽이 유력하나 오류 코드는 미확인 — Kratos Eigen 래퍼가 코드를 버린다. 행렬을 꺼내 pypardiso로 오류 코드·대각 상태를 보는 진단 `scripts/diag_v02_singular.py`를 실행 중(결과 `docs/phase1/evidence/v02_diag/diag_l5.json`, 나오면 덧붙임).
- 오차가 비단조다. 레벨 2·4(격자가 구멍과 비정렬, 셀이 정사각형 아님: 200/67=2.985, 40/13=3.077)는 작고, 레벨 1·3(정렬, h가 r의 약수)은 크다. 레벨 1·3의 경계 peak가 기준보다 **크다**(과대).
- 모든 fallback 셀이 구멍 주변(crown)에 있다. 반력 합은 정확(3.99999999988e-5 = P·A).

## 2. 절단 셀 적분 규칙 (의심 지점)

`wht_iga/geometry/quadrature.py`:
- `cell_quadrature` (l.353–): 정확한 모멘트 M(Legendre 텐서, 차수 ≤ q) → 후보점(셀 Gauss (q+2)^3 내부점 + 표면 leaf 깊이 2·3 + 조밀 격자) → `fit_moment_rule`(NNLS + 역추적 가지치기) → 잔차 `rnorm ≤ _NNLS_TOL·|cell|`(`_NNLS_TOL=1e-12`, 절대·체적 비례)이면 채택, 아니면 `octree_fallback`(깊이 3, 내부 leaf Gauss (p+1)^3, 절단 leaf 2^3점 + 안팎 판정 — 모멘트를 맞추지 않음, 정확도 낮음).
- `fit_moment_rule` (l.256–): 후보 2048개 초과 시 등간격 추출(`_NNLS_MAX_CANDIDATES`), scipy `nnls(maxiter=10·n)`; RuntimeError면 잔차 +inf → fallback.
- 모멘트: `_cell_moments` / `_cross_section_moments` (발산정리 + 단면 1D 정확 적분).
- 관찰: q=4(V02)에서 fallback이 생기고, V06 시험은 "q=4 NNLS가 얇은 sliver 셀에서 발산"이라며 q=2를 쓴다(`tests/kratos_spikes/test_p2mesh_v06.py` run_v06 docstring).

경계 peak probe: `wht_iga/solver/trimmed.py::boundary_peak_stress` (구멍 원통면 삼각형 위 Duffy 점에서 이산해 변형률 → 응력 → 최대 주응력). 체적 peak: `peak_max_principal` (적분점 응력).

## 3. 질문

1. **비단조 원인 분리**: 레벨 1·3의 과대 오차가 (a) crown fallback 셀의 부정확한 적분, (b) extension OFF 상태의 작은 절단 셀(조건수·국부 진동), (c) 경계 probe 방식(이산해를 경계에서 외삽), (d) 격자–구멍 정렬에 따른 특정 cut 형상 중 무엇이 지배적인가? 각 가설을 가르는 **최소 실험**(무엇을 바꾸고 무엇을 보면 되는지, 비용)을 설계해 달라. 가능하면 코드를 읽고 사전 판정까지.
2. **fallback의 근본 원인**: q=4에서 NNLS가 실패/잔차 미달하는 기하 조건은? 허용치 `1e-12·|cell|`(절대·체적 비례)가 작은 체적비 셀에서 과하게 엄격한가(모멘트 자체의 반올림 오차 수준 대비)? Legendre 기저 조건수, 후보점 배치(구멍 근처 얇은 sliver), 등간격 추출이 원인일 수 있는가?
3. **처방 비교**: (i) 허용치를 상대값(‖M‖ 또는 M₀₀ 기준)으로, (ii) 후보점 생성 개선(표면 근처 Gauss 점, 셀 내부 부분영역 Gauss), (iii) fallback을 모멘트 맞춤 octree(leaf마다 저차 모멘트 맞춤)로 교체, (iv) Saye(2015) 고차 implicit quadrature 또는 다른 공개 알고리즘의 적합성. 각각의 정확도 영향, 구현량, 위험을 비교하고 권고안 1개를 골라 달라.
4. **G2 판정 전망**: extension ON(θ=0.2, `stabilize/extension.py`, 2026-11-14 Marsden 계수로 정정됨)과 fallback 처방 후, 4레벨 이상에서 ≤1%·수렴률 ≥1.8이 현실적인가? 수렴률 계산에서 비단조 레벨을 어떻게 다뤄야 정직한가(계획서 G2 문구: "단조 감소 구간 log-log 최소자승" — `docs/phase2/G2_report.md` §4)?
5. **레벨 5 분해 실패**(결과가 덧붙여져 있으면): 메모리인지 특이성인지 판정과, extension OFF에서 지지가 아주 작은 기저가 Pardiso LU pivot에 미치는 영향.

## 4. 산출물 형식

- 질문별: 결론 1–2문장 → 근거(파일:줄, 수식) → 확인 실험(명령·기대값) → 불확실성.
- 마지막에 **실행 순서 제안**(통합 담당이 바로 돌릴 수 있게: 스크립트, 환경변수, 예상 소요).
- 추정과 확인을 구분해 표시. 숫자를 만들어내지 말 것.

## 5. 참고 자료

- `docs/INTEGRATION_HANDOFF.md` (현재 작업 상태, 고속화 기록)
- `docs/verification/V01-ext.md` (extension 결함과 정정)
- `docs/phase1/TRACK_GEO.md` (절단 적분 설계 원래 근거)
- `docs/phase2/G2_report.md`, `docs/REQUEST_P2_V02_FIX.md`
- 계획서 §7.4, §7.5, §10.1 V02, §0.5 개정 이력
