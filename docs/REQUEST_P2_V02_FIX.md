# Phase 2 수정 지시 — V02 경계조건 오류 + Phase 1 후처리 결함 (2026-11-14)

근거: `WHT_IGA/docs/REVIEW_P0_P1.md` §2, §3. 기준 문서: `research-coreform-IGA/master-plan_v2.1.md`, 규칙: `research-coreform-IGA/implementation-request_v1.md` §2, §4(K1–K14), `WHT_IGA/AGENTS.md`.

## 0. 요약

V02(Roark 17.1-7a, 구멍 판 인장) 오차 −23%의 주원인은 소형 셀이 아니라 **경계조건 오류**다.
현재 스크립트·시험이 판의 자유 모서리(y=±20) 전체에 `uy=0`을, 양 z면에 `uz=0`을 걸고 있다. Coreform 정의는 **왼쪽 면 pinned**뿐이다.
같은 코드로 extension 없이 왼쪽 면만 고정하면 오차가 다음처럼 바뀐다.

| 면내 h | 기존 BC | 오차 | 왼쪽 면만 고정 | 오차 |
|---|---|---|---|---|
| 10 | 3.299e-6 | −23.6% | 3.814e-6 | −11.7% |
| 5 | 3.482e-6 | −19.4% | 4.586e-6 | +6.2% |
| 2.5 | 3.330e-6 | −22.9% | 4.502e-6 | +4.3% |

- 기준값: 4.3175e-6 (순단면 Kt 2.15875 × σ_nom 2e-6)
- Coreform immersed p=2 결과(`roark-17p1-7a-qoi.csv`): h=6 −24.6%, h=3 −5.7%, h=1.5 +0.1%, h=0.375 +0.8%

따라서 "G1 V02 23% → Phase 2 extension이 해결"이라는 조건부 종료 근거는 **무효**다. 현재 Phase 2 시험은 틀린 BC 위에서 extension 효과를 측정하고 있다.

## 1. 필수 작업 (G2 보고 전)

### T1. V02 경계조건 수정
- 대상: `scripts/v02_roark.py`, `tests/**/test_p2stab_v02_fcm.py`, `tests/**/test_p2mesh_v02.py`, 그리고 V02 BC를 만드는 공용 코드가 있으면 그것까지.
- 자유 모서리 `y=±20`의 `uy=0` 삭제, 양 z면의 `uz=0` 삭제.
- 다음 중 하나로 구성하고, 보고서에 어느 쪽을 썼는지 적는다.
  - (a) Coreform과 같게: 왼쪽 면 pinned(ux=uy=uz=0), 오른쪽 면 인장 하중
  - (b) 1/4 대칭 모델: x·y 대칭면에서 법선 변위 0, 강체 z 구속은 한 점만
- BC 정의를 V02 전용 헬퍼 1곳으로 모아서 스크립트와 시험이 같은 함수를 쓰게 한다. 복사본을 두지 않는다.

### T2. extension 없이 재기준선 측정
- Phase 1 경로(extension 끔)로 면내 h = 10, 5, 3, 2.5, 1.5, 0.75 (가능한 범위까지) 측정.
- 표: h, 자유도 수, peak, 오차, 소형/fallback 셀 수, 구멍 정점 주변 fallback 셀 수.

### T3. extension 효과 분리
- 같은 격자·같은 BC에서 extension on/off를 비교한다.
- θ(소형 셀 임계값)와 α를 바꿔 가며(sweep) 비교한다.
- G2 보고서에는 다음 두 값을 **따로** 적는다.
  - BC 수정 효과 = (기존 BC, ext off) → (수정 BC, ext off)
  - extension 효과 = (수정 BC, ext off) → (수정 BC, ext on)
- extension이 개선을 못 보이면 그대로 보고한다. 수치를 맞추려고 θ·α를 V02에만 맞춰 조정하지 않는다.

### T4. 잘못 고정된 단언 교체
- `test_p2mesh_v02.py`의 `0.5×G1 < peak < 1.5×G1` 단언은 틀린 G1 값을 정답으로 쓰고 있다. 삭제한다.
- 기준값 4.3175e-6에 대한 상대오차 허용치로 바꾼다. 허용치는 T2·T3 실측 후 정하고, 근거(격자, Coreform 같은 h의 오차)를 주석으로 단다.
- `V02-ext.md` 등 문서에 있는 "~23% 재현 기대" 문구를 정정한다.

### T5. 수렴 기준 (G2 판정용)
- 4단계 이상 격자에서 최세 격자 오차 ≤1%, 단조 수열 구간에서 관측 수렴률 ≥1.8.
- 미달이면 미달로 보고한다. 기준을 낮추지 않는다.
- 측정은 체적 probe와 경계(구멍 CAD 면) probe를 따로 보고한다.

### T6. BC 수정 후 재확인
- penalty β, p=3 결과를 수정된 BC로 다시 확인한다(기존 결론이 틀린 BC에 의존했는지 점검).

## 2. Phase 1 파일 결함 (Phase 2가 이 파일들을 수정 중이므로 이쪽에서 처리)

| # | 위치 | 결함 | 수정 | 시험 |
|---|---|---|---|---|
| F1 | `solver/trimmed.py:228-237` `boundary_peak_stress` | D·ε(ε는 공학 전단변형률)로 이미 실제 전단응력을 얻었는데 다시 /2 함 → 전단 절반 | `/=2` 6줄 삭제 | 순수 전단 장에서 τ = G·γ 확인 |
| F2 | `solver/trimmed.py:141-150` `peak_max_principal` | Voigt 순서를 [xx,yy,zz,yz,xz,xy]로 가정. Kratos는 [xx,yy,zz,xy,yz,xz] | 순서 매핑 수정, 상수 1곳에 정의 | 알려진 텐서의 최대 주응력 일치 |
| F3 | `pipeline.py:171-180`, `db.py:165` | 수렴 확인 없음, 결과 DB에 `converged=1, n_iter=1` 고정 기록(꾸며낸 값) | 솔버의 실제 수렴 여부·반복 수 기록, 미수렴이면 경고 + DB 반영 | 미수렴 강제 사례에서 `converged=0` 기록 |
| F4 | `solver/validate.py:166-170` | `QUAD_FALLBACK_LOADED`가 `face:` 엔티티만 검사, 브리지는 `cell:`로 보고 → 경고가 절대 안 뜸 | 하중 면과 fallback 셀의 인접 관계로 검사 | 하중 면에 fallback 셀이 있는 사례에서 경고 발생 |
| F5 | `pipeline.py:134-157`, `trimmed.py:81` | 고정 ID 오프셋(+10000, +40000, 50000) → extension 구속이 약 4만 개를 넘으면 ID 충돌 | 단조 증가 ID 할당기 하나로 통일 | 대량 ID 할당 시 중복 없음 |
| F6 | `trimmed.py:179-184` | 구멍 면을 V02 전용 경험식으로 찾음 | 측정할 면 ID를 인자로 받음 | — |

F1·F2는 V02 정점(전단 0)의 오차는 설명하지 못한다. 다만 다른 검증 사례의 응력 값을 틀리게 만드므로 G2 전에 고친다.

## 3. 규칙 (위반 금지)

- 환경: micromamba `whtiga`. 패키지는 `uv pip`으로만 설치. conda 명령, pip Kratos wheel 설치 금지.
- `E:\kratos-build`, kratos-src 수정·재빌드 금지. `wht_topopt`·`aieng` import 금지.
- Kratos 구속 상수 c는 항상 0 (K-규칙). 구속 사슬 금지.
- `skip`, `xfail`, 허용치 완화로 실패를 숨기지 않는다. 실패는 실패로 보고한다.
- commit·push 하지 않는다.
- 결과 보고·문서는 한국어로 작성한다.

## 4. 완료 보고 형식

`WHT_IGA/docs/phase2/G2_report.md`(또는 기존 G2 보고서)에 다음을 포함한다.
1. 사용한 BC(a/b)와 BC 헬퍼 위치
2. T2 표(ext off 재기준선)
3. T3 표(ext on/off, θ·α sweep). BC 효과와 extension 효과를 분리한 수치
4. T5 수렴 판정(수렴률 계산식 포함), 체적/경계 probe 별도
5. F1–F6 각각의 수정 위치(파일:줄)와 추가한 시험 이름
6. 전체 시험 실행 결과(통과/실패 수, 실행 명령)
7. 해결 못 한 항목과 이유
