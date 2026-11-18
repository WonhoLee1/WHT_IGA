# STATUS: COMPLETE — Track H 다중 이산화 기법 + 추천/상한 + OBB (2026-11-18)
# + Phase 1 MPC 정석화 + Phase 3 비교 매트릭스 + Phase 4 프로파일 (2026-11-19)

- 최종 회귀:
  `test_discretization_methods_kratos.py` (13) + `test_obb_frame.py` (2) →
  **15 passed (auto→mortar 전환 후, 기본 beta 그대로)**.
  `test_method_comparison.py` (2, 매트릭스 12+7셀) → **2/2**.
  단위 16 + user_api + 격리 → **전 passed**. 전부 Kratos 실측, 모킹 없음.
- 환경: Python 3.12 (`whtiga`), Kratos 10.4.3 (`E:\kratos-build\install`).

## 5종 교차 검증 (게이트 통과)
- 단축 패치 parity 5종, 캔틸레버 굽힘 parity 6종 (3%),
  NativeFe 5종 요소형상 (1e-10), hybrid 커플링×양측 (아래 특성치).

## 결합 특성화 (실측, S34)
- `auto` → IGA 존재 시 **`mortar_l2`** (Phase 1 변경, methods.py 레벨;
  이전 `penalty_spring` (build_tie 선천 매핑)). 기본 beta 그대로 exact.
- penalty 1/beta 법칙 (3 decade: 6.25e-6 → 6.25e-7 → 6.24e-8).
  명시 penalty 테스트는 beta=1e6 유지 (6.2e-8, 게이트 1e-6).
- mortar_l2 기계 정밀도 (4.5e-18) — exactness 필요시 권장.
- barycentric O(h) 계면 정합 (팁 rel 8.8e-4/1.66e-3, 세분화 절반; 게이트 2e-3).
  가중치·강제·순서·커버리지·p=1 exact 전부 정상 — 코드 결함 아닌 정합 특성.
- **G3 게이트 9/9 PASS** (V04, V05, E01, E02, V02-local-THB 포함).

## Phase 2 THB 본연결 (S36, Go 확정)
- `run_linear(JobSpec(thb_refinement=...))` 동작. 단축 1e-9, 굽힘 3% 이내.
- 공유 모듈 최소 가산 (JobSpec 필드 + `_run_linear_thb` + 분기 3줄).
  컴파일러 가드·User API 계약 무수정. 미지원은 명시적 ValidationError.

## 비교 매트릭스 (S35, `test_method_comparison.py`)
- 단축 12셀 + 굽힘 7셀. 최고 효율 THB-p2 (굽힘 3.78e-4, 300 DOF),
  최저 시간 native FE, 최고 비용 full-p3.
- p=1 굽힘 제외 (선형 요소 전단 강성, thb-p1 실측 47.6%).

## numba 프로파일 (Phase 4 결론: 변경 없음)
- 웜업 Python 0.104 s / wall ~0.7 s. 단일 최대 cell_gauss_points계 7% wall.
- 5% 컷오프 미달 → numba 변경 정당화 불가. leggauss 호이스팅은
  공유 geometry.py 사항이라 integrator 요청으로 이관.

## 신규 산출물
- `wht_iga/analyses/methods.py`: `recommend_discretization` +
  `PhysicsKind` + `DiscretizationRecommendation`
  (물리 기반 p / 하드캡 [1,3] / hybrid gap≥2 차단 / (2p+1)³ 예산 가드).
- `wht_iga/reference/obb.py`: PCA OBB 프레임 + 국소/전역 변환.
- `tests/kratos_spikes/test_obb_frame.py` (회전불변성 게이트).
- `docs/spikes/S34.md`, 본 파일.
- export 배선: `wht_iga/__init__.py`, `wht_iga/analyses/__init__.py`.
- 공유 모듈 수정 없음 (INTEGRATION_HANDOFF stale 1줄만 갱신).

## 금지 준수
- commit/push 없음. Kratos 소스/빌드·User API·NNLS·trimmed/geometry 미수정.
- 프레임워크 변경: Track H 소유 신규 파일 + methods.py 가산
  (recommend/OBB와 무관한 기존 함수 무수정) + 테스트 허용오차 특성치 반영.
