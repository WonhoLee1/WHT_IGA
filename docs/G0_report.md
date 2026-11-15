# G0 게이트 보고서 (Phase 0 완료 판정)

- 날짜: 2026-11-14 / 브랜치: WHTIGA / Kratos 10.4.3 (`E:\kratos-build\install`)
- 판정: **PASS** — p=3 교차비교 보완 완료 (아래 G2). **사용자 승인됨 (2026-11-14)**.
- 경로 제안: R1 (근거 `docs/adr/ADR-001-route.md`)

## 게이트 기준 대조 (v2 Phase 0 + 요청서 §9 DoD)

| # | 기준 | 결과 |
|---|---|---|
| G1 | S01–S05 결론 확정 | PASS — 전부 예 (S02 4.4e-16, S03 27점/가중치합, S04 K 1e-16·패치 3.1e-15, S05 4.2e-17·2.7s) |
| G2 | untrimmed 큐브 p=1,2,3 Kratos vs 참조 조립기 변위차 ≤1e-10 | **PASS** — p=1,2,3 `test_s04_sm_elements.py` parametrize로 확인 (본 작업 실행): K 상대오차 p=3 3.1e-16/4.2e-16, 패치 최대오차 3.1e-15 |
| G3 | `check_env.py` 종료 코드 0 | PASS — S01이 동일 검사를 pytest로 수행, 13/13 (HasMKL=True) |
| G4 | `pytest WHT_IGA/tests` 결과 첨부 | PASS — kratos_spikes 61/61 (통합 담당 재실행), `not kratos` 103/103 (재실행). 실패 0 |
| G5 | TODO 스텁·skip·xfail 없음 | PASS — 신규 45개 파일 grep 무적중 (skip/xfail/TODO/FIXME/NotImplemented) |
| G6 | G0 보고서 + ADR-001 제출 → 사용자 승인 대기, 승인 전 Phase 1 금지 | 본 문서 + ADR-001로 충족. **사용자 승인됨 (2026-11-14)** — ADR-001 상태와 일치 (REVIEW #10) |

## S01–S18 요약 (상세: `docs/spikes/SNN.md`, `TRACK_A/B.md`)

| spike | 결과 | 핵심 수치 |
|---|---|---|
| S01 env | 예 | 13/13, HasMKL=True |
| S02 항등사상 | 예 | 4.4e-16, detJ 1.1e-15, 회전 동일 |
| S03 셀 QP 기하 | 예 | 기하 1개/셀, 가중치합=셀체적, SFV행합 1 |
| S04 SM 요소 | 예 | K 8.7e-17~2.3e-16, 패치 3.1e-15, TL/UL 초기화 OK |
| S05 구속식 | 예 | 4.2e-17, 1e5건 2.7s, DofSet 포함 |
| S06 고유치 | 부분 | Kratos QP 질량=0 → solver 실패. scipy T-투영은 hex 대비 1e-9, 구속 전후 1e-8 |
| S07 표면 QP | 예 | 면적합 1.0 (1e-13). 점당 기하 1개 (v2 정정) |
| S08 하중 | 예 | 반력합 정확, follower 1.56e-6, +p=외법선 |
| S09 penalty | 예 | 3.4e-12. 표면 규칙 ≥2p 필수 (2차 규칙 고착오차 발견·해결) |
| S10 nodal/spring | 예 | 2.0 / 4/3 (<1e-9) |
| S11 열 | 부분 | 전도 50.0, film 수렴 66.659 (해 66.667), radiation 2.8350. Laplacian조건 불가, BDF 질량에 rho·c_p 없음 |
| S12 열팽창 | 예 | 모서리 1e-3, 반력 2e-16 |
| S13 혼합 | 아니오 (예상됨) | MixedVolumetric 초기화 크래시 (결정적). MixedStrain OK → C++ MixedUP0 Phase 5 |
| S14 접촉 | 부분 | 조건 등록 OK, 스프링-tie 힘전달 -2.0·점프 5e-11. 완전 마찰풀이 보류 |
| S15 테이블 | 부분 | 보간 OK이나 법칙이 E(T) 미반영 → Properties 직접 갱신 failover 확인 |
| S16 스케일 | 예 | 97336 요소 52.5s (<60s), +2.92GB (<8GB) |
| S17 출력 | 부분 | QP JSON 6자리 절삭, VTK는 Brep 입력 요구 → 자체 writer 유지 |
| S18 Nitsche | 아니오 (정보성) | 3D segfault, 멤브레인 전용 → C++ 트랙 |

## Phase 1 진입 체크리스트

1. [x] S04 p=3 교차비교 — 완료. 케이스가 S04 파일 parametrize에 포함됨 (본 작업 실행: K 3.1e-16/4.2e-16, 패치 p=2 3.1e-15)
2. [ ] 빌더 규약 반영: volume 수명 보관, 표면 규칙 ≥2p, 점당 표면 기하, Dirichlet c=0 규칙, ROTATION Fix, thermal 기본 settings
3. [ ] REVIEW 후속 3건: F1 운용안 확정 (Python consistent M 경로 — C++ 불필요),
   F2 상수 안전장치 (`add_linear_tie` 가드 + Newton 검증 시험),
   F3 Accessor spike (공간 변화 온도, Phase 4 전)
4. [x] 사용자 승인됨 (2026-11-14, 본 보고서 + ADR-001)

## 트랙 C (S19–S24) 요약 (REVIEW #9)

상세 `docs/spikes/S19–S24.md`, 시험 `tests/kratos_spikes/test_s19_*`–`test_s24_*`, 경로 `docs/adr/ADR-001-route.md` "트랙 C" 절.

| spike | 결과 요약 |
|---|---|
| S23 비적합 매퍼 | 본 작업 실행: `nearest_element` 분할단위성 0.0·선형 재현 1.1e-16, 결합 patch 간극 0.085987 (대역 단언), utility 48 구속·선형 오차 7.942e-02 (대역 단언) → 매퍼 경로 채택 |
| S20 TableAccessor | 본 작업 실행 1건: E(T) clamped bar 해석해 대비 6.0e-9 (<1e-8). 나머지 S20 시험은 본 작업 미실행 |
| S19/S21/S22/S24 | 본 작업 미실행 — 시험 파일·문서 참조. R1 변경 없음 (S22 접촉·S19 일부 법칙은 C++ 트랙으로) |

## v2 수정 제안 (본문 미수정, X8 준수)

ADR-001 §"채택하는 failover와 v2 정정 제안" 12항 참조.
공유 모듈 변경 제안 1건: `geometry.make_nurbs_volume`의 volume 보관
(TRACK_B §공유 모듈 변경 제안 #1, 대안으로 빌더 규약 문서화 가능).

## 산출물

- `docs/adr/ADR-001-route.md` (R1 제안)
- `docs/spikes/S01–S18.md` + `TRACK_A/B.md` (18+2)
- 시험: `tests/kratos_spikes/test_s*.py` 18개 + 신규 모듈 4개
  (`elements.py`, `constraints.py`, `surface.py`, `thermal.py`)
- 에이전트 commit 없음 (HEAD 8bff77a 유지 확인).
- 추가 검증 기록: `docs/spikes/TRACK_A_REVIEW.md` (독립 재실행 27/27 in 99.3s,
  F1 근본원인 확정 + lumped/consistent 실측, F2–F3 후속 항목). 본 보고서에 반영됨.
