# Phase 0 나머지 · Phase 1 리뷰 종합 (2026-11-14, 통합 담당)

범위: Phase 0 트랙 B·C(S01, S07–S24), ADR-001/002, G0, Phase 1 코드·시험. Phase 2 진행분(stabilize, operators, meshio_ext, test_p2*)은 제외. 트랙 A는 `spikes/TRACK_A_REVIEW.md`.
V02 23% 원인 검증(3번 리뷰)은 진행 중 — 결과는 아래 §3에 추가 예정.

## 1. Phase 0 트랙 B·C, ADR, G0 — 수정 후 수락

재실행: 16개 파일, **74개 통과**, 실패·skip 0 (14.0 s). 공유 모듈 수정 없음.

| # | 심각도 | 위치 | 내용 | 수정 |
|---|---|---|---|---|
| 1 | 높음 | `kratos_backend/constraints.py:47-55`, ADR-001:33-35, G0 체크리스트 3 | ADR/G0는 "c≠0 거부" 완료라 하지만 코드는 여전히 상수를 받아 Kratos에 넘김 | `if np.any(constant != 0): raise ValueError` |
| 2 | 중상 | G0_report.md:12 | p=3 Kratos 대 참조 ≤1e-10 근거가 임시 스크립트뿐. S04는 p=1,2만 | S04 parametrize에 p=3 추가 |
| 3 | 중 | S09.md, ADR-001:30, ADR-002:56 | penalty 오차 4.5e-5(β 무관)의 원인은 적분 차수가 아니라 **rank 부족**(1셀 p=2 윗면 제어점 9개 대 적분점 6개)일 가능성. p=3 order 4는 12점 대 16 자유도 | 규칙: 면-셀당 적분점 ≥ (p+1)² 그리고 총차수 4p 정확. p=3 1셀 penalty 시험 + 4.5e-5 회귀 시험 |
| 4 | 중 | ADR-002:41-47,56-58,70 | β sweep이 정확히 표현되는 장(patch)·격자 1개뿐 → 적분/rank 오차 검출 불가 | 비표현 장(굽힘), 격자 2개로 재측정 전까지 c_pen=3e5는 가설 |
| 5 | 중 | test_s08:90,124 | 압력 부호 규약을 문서는 주장하지만 시험은 `abs`만 | `fz == -2.0` 단언 |
| 6 | 중 | test_s23:266-271 | 결정 근거 수치(lin_err 7.94e-2)에 단언 없음 | 대역 단언 추가 |
| 7 | 중 | S24.md, TRACK_C | Kratos explicit 질량 이전은 맞음. "비변분적"은 과장 — 실제 문제는 음수 extension 가중치에서 row-sum 질량이 0 되는 것. `Initialize()` 후 NODAL_MASS를 HRZ로 덮어쓰는 혼합안 미시험 | 표현 수정 + 혼합안 시험 또는 미시험 명시 |
| 8 | 중 | S19.md:13, S21 | TotalLagrangian + 적분점 기하에서 간헐적 DETJ0=-0.5/NaN, 그런데 S21이 TL을 표준으로 채택. MoveMeshFlag·ResidualCriteria 주장에 음성 시험 없음 | 음성 시험 추가, Phase 5 전 DETJ0 원인 조사 |
| 9 | 중 | ADR-001:24,32, G0 | 트랙 C(S19–S24) 미반영. ADR 항목 4(S15 Properties 갱신)는 S20(TableAccessor + ThermalElasticIsotropic3D로 E(T) 동작)으로 대체됨 | 트랙 C 절 추가 |
| 10–18 | 낮음 | S18 시험이 Kratos를 검사 안 함, S14/S17 존재 확인만, `surface.py` assert·Nitsche 문구·volume 수명, `thermal.py` 하드코딩 CONDUCTIVITY/DENSITY=1·degree 6 규칙 없음, ADR 항목 12 중복, ADR "승인" vs G0 "승인 대기" 불일치 | 각각 정리 |

판정: 트랙 B·C, ADR-001(R1 선택 자체는 근거 충분), ADR-002, G0 모두 **수정 후 수락**.

## 2. Phase 1 코드 — 수정 후 수락

재실행: **46개 통과**, 실패·skip 0. 모듈 import 시 Kratos 미로드, 금지 import 없음.

| 심각도 | 위치 | 내용 | 수정 |
|---|---|---|---|
| **높음** | `solver/trimmed.py:228-237` `boundary_peak_stress` | D·ε(공학 전단변형률)로 이미 실제 전단응력인 값을 다시 /2 → 전단 절반 | `/=2` 6줄 삭제 + 순수 전단 시험 |
| **높음** | `trimmed.py:141-150` `peak_max_principal` | Voigt 순서를 [xx,yy,zz,yz,xz,xy]로 가정, Kratos는 [xx,yy,zz,xy,yz,xz] | 순서 매핑 수정 + 시험 |
| 중 | `solver/validate.py:166-170` | `QUAD_FALLBACK_LOADED`가 `face:` 엔티티만 보는데 브리지는 `cell:` 보고 → 실제로 절대 발동 안 함 | 셀-면 인접 검사 |
| 중 | `pipeline.py:171-180`, `db.py:165` | 수렴 확인 없음, 결과 DB `converged=1,n_iter=1` 고정(꾸며낸 값) | 실제 값 기록·미수렴 경고 |
| 중 | `pipeline.py:134-157`, `trimmed.py:81` | 고정 ID 오프셋(+10000, +40000, 50000) → extension 구속 약 4만 개 넘으면 충돌 | 단조 ID 할당기 |
| 중 | `pipeline.py:123-157`, `penalty.py:146` | BC·하중이 배경 박스 면에만. CAD 면 적분점 미사용, BC-C 없음 (V02는 격자=판 bbox라 우연히 동작) | Phase 1 한계로 명시 또는 연결 |
| 중 | `pipeline.py:218-246`, `geometry/quad_db.py` | 데이터 버스 qp 표 비어 있음, quad_db 미사용·0 값 기록·재실행 시 오류, 스키마 버전 2종 | 정리 |
| 중 | `trimmed.py:179-184` | 구멍 면을 V02 전용 경험식으로 찾음 | 면 ID 명시 전달 |
| 낮음 | 여러 곳 | 공차 불일치(1e-12·h vs 1e-9·max h), fallback Gauss 차수 p 미전달, assert, p>3 하중 규칙 묵시 대체, Phase 2 로직이 Phase 1 파일에, 빈 BC 통과, `+PRESSURE` 바깥 방향 규약 | 각각 정리 |

시험 약점: VTU 값 비교 없음, V02/V03 회귀 시험 없음(스크립트뿐), **응력 후처리 시험 없음**(높음 2건이 숨은 이유).

V02 관련: 높음 2건은 구멍 정점에 전단이 없어 V02 오차(20–23%)를 설명하지 못함. "소형 셀 탓"은 미증명. 재측정 조건: 응력 버그 수정, octree fallback 셀 제거(h=2.5에서 48개), `boundary_peak_stress` 해석해 보간 검사.

## 3. G1 조건부 종료 · V02 23% 원인 — **조건부 종료 근거 불인정**

- 기준값(4.3175e-6, 순단면 Kt 2.15875 × σ_nom 2e-6), 측정 위치(구멍 CAD 면 1656점, 스플라인 미분), tessellation, 적분(모멘트 1e-12, 합력 정확)은 문제 아님.
- **주원인(신뢰도 ~85%): 경계조건 오류.** V02 스크립트(`scripts/v02_roark.py`)와 `test_p2stab_v02_fcm.py`, `test_p2mesh_v02.py`가 **자유 모서리 y=±20 전체에 uy=0**과 양 z면 uz=0을 걸고 있음. Coreform 정의는 "왼쪽 면 pinned"뿐. 자유 모서리 구속이 횡방향 강성을 더해 정점 응력을 낮춤.
- 진단(같은 코드, extension 없음, 왼쪽 면만 고정):

  | 면내 h | 기존 BC | 오차 | 왼쪽 면 고정 | 오차 |
  |---|---|---|---|---|
  | 10 | 3.299e-6 | −23.6% | 3.814e-6 | −11.7% |
  | 5 | 3.482e-6 | −19.4% | 4.586e-6 | +6.2% |
  | 2.5 | 3.330e-6 | −22.9% | 4.502e-6 | +4.3% |

  Coreform 자체 immersed p=2: h=6 −24.6%, h=3 −5.7%, h=1.5 +0.1%, h=0.375 +0.8% (`roark-17p1-7a-qoi.csv`, 소형 셀 처리 언급 없음).
- 남은 오차(+4.3%, 마지막 단계 수렴률 ~0.5): 약지지 절단 셀 기저에서의 probe 평가, 최세 격자 fallback 셀 48개 등이 후보이나 대조 실험 없음 → 소형 셀 기여분은 아직 분리 불가.
- G1 보고서의 "수렴 정체 = 소형 셀 징후" 주장은 대조 실험 없이 나온 것. `test_p2mesh_v02.py`의 `0.5×G1 < peak < 1.5×G1` 단언과 `V02-ext.md`의 "~23% 재현 기대"는 **틀린 값을 고정**함.

**G2 전 필수**: ① V02 BC 수정(왼쪽 면 pinned 또는 sym_x/sym_y 대칭 모델) 후 extension 없이 Phase 1 재기준선 ② 같은 격자·BC에서 extension 유무 대조 + θ/α sweep — G2 보고에 BC 수정 효과와 extension 효과를 분리 기재 ③ 면내 h 3, 1.5, 0.75 포함 4단계, 최세 ≤1%, 단조 수열에서 수렴률 ≥1.8 ④ 체적/경계 probe 별도 보고, 정점 주변 fallback 셀 보고 ⑤ G1 고정 단언을 기준값 기반 허용오차로 교체 ⑥ BC 수정 후 β·p=3 재확인.

## 4. 우선 조치 (Phase 2와 병행 가능)

0. **[최우선, Phase 2 진행 측에 즉시 전달] V02 경계조건 수정 + 재기준선** (§3) — 지금 STAB/MESH가 V02를 틀린 BC로 검증 중
1. `constraints.add_linear_tie` c≠0 거부 (트랙 B·C #1)
2. 응력 후처리 버그 2건 + 시험 (Phase 1 높음)
3. 결과 DB 수렴 값 실제 기록
4. penalty 적분 rank 규칙 + p=3 1셀 시험 (S09/ADR-002)
5. ADR-001·G0에 트랙 C 반영, 승인 상태 일치

## 5. 조치 현황 (2026-11-14, 통합 담당 처리분)

| 항목 | 조치 | 시험 |
|---|---|---|
| §1 #1 c≠0 | `add_linear_tie`가 0이 아니거나 길이가 틀린 상수를 거부 | `tests/unit/test_constraints_guard.py` (Phase 2 쪽 `test_p0_constraints_guard.py`도 있음) |
| §1 #2 S04 p=3 | parametrize에 (3,1셀), (3,2×2×2) 추가 | `test_s04_sm_elements.py` |
| §1 #3 penalty rank | 판별 실험: 원인은 rank 부족(면-셀 적분점 < (p+1)²), 적분 차수 아님. p=3 + 규칙 4는 1.772e-4 고착. 6차 12점 규칙 추가 | `test_s09_penalty_support.py` 13개, `S09.md` 정정. `rule_for_degree(3)` 변경은 Phase 2 쪽에 전달 |
| §1 #5 S08 부호 | `fz == -2.0` 단언, follower 부호 | `test_s08_pressure_load.py` |
| §1 #6 S23 | 0.086, 7.94e-2 대역 단언 | `test_s23_nonmatching_mapper.py` |
| §1 #7 S24 | "비변분적" 정정: 행합 1인 T에서 Kratos 질량 이전 = TᵀMT row-sum lumping. 실제 문제는 음수 가중치. HRZ 덮어쓰기 혼합안은 미시험으로 명시 | `S24.md` 정정 절 |
| §1 #8 DETJ0 | **근본 원인: NurbsVolume 수명.** TL 기준 야코비안이 해제된 parent를 읽음. `make_nurbs_volume`에서 `model_part.AddGeometry(volume)`로 수정(TL NaN 6/12 → 0/12). S21 음성 시험 2개: 두 설정 모두 Kratos가 `converged=True`로 틀린 답을 냄 | `test_s19` TL 추가(반복 4회), `test_s21` 음성 2개, `S19.md`·`S21.md` |
| 낮음 S14 | 두 mortar 조건이 모두 생성되는지 단언 | `test_s14_mortar_contact.py` |
| 낮음 S17 | "Brep 필요" 주장은 검증 불가(최소 설정이 무관한 경로 오류로 먼저 실패) → 주석으로 명시, 단언 추가 안 함 | `test_s17_output_processes.py` |
| 낮음 surface.py | `assert` → `RuntimeError`, Nitsche 문구 정정(3D 불가, S18), 표면 래퍼도 `model_part.AddGeometry`로 수명 고정 | S07·S08·S09·S11·P1 loads·P1 cube 30개 통과 |
| 낮음 thermal | 면 CONDUCTIVITY, DENSITY 하드코딩을 인자로(기본값 유지) | S11·S12 통과 |

미처리: "degree 6 규칙 없음"(thermal 쪽 위치 미확인), S18(Phase 2 쪽 처리 중). ADR·G0 갱신은 Phase 2 쪽이 진행함.
