# G1 게이트 보고서 (Phase 1 검증)

- 날짜: 2026-11-14 / 브랜치: WHTIGA / Kratos 10.4.3 (`E:\kratos-build\install`)
- 판정: **V03 PASS / V02 부분 (게이트 미달, 원인 확정)** — Phase 2 진입 판단은 사용자
- 실행: `WHT_IGA/scripts/v03_beam.py`, `v02_roark.py` (통합 담당 직접 실행)

## V03 외팔보 휨 — PASS (5종, 게이트 5%)

| 케이스 | tip | EB ref | rel err |
|---|---|---|---|
| gravity (q=784.8) | -5.650e-5 | -5.606e-5 | 0.80% |
| body force (동등 78480 N/m³) | -5.650e-5 | -5.606e-5 | 0.80% |
| pressure uniform (5000 Pa) | -3.600e-5 | -3.571e-5 | 0.80% |
| pressure total-force (동일 resultant) | -3.600e-5 | -3.571e-5 | 0.80% |
| traction end-shear (P=500 N) | -9.580e-5 | -9.524e-5 | 0.60% |

방향 전부 하향 (Coreform과 일치). 압력 부호 규약: pipeline의 `PressureLoad`
+magnitude = 외향 pull (S08에서 확정). Coreform 내향 압력과 반대 방향이므로
V03 입력 부호를 뒤집어 맞춤 — 크기 게이트에는 영향 없음 (기록용).

## V02 Roark thin plate — 부분 (finest 22.9%, rate 0.02; 게이트 1%·1.8 미달)

E=1, nu=0, P=1e-6 tension. Grid = plate bbox exactly (padding 0, 외곽면 conforming).

| h | peak (boundary) | vol peak | ref | err |
|---|---|---|---|---|
| 10 | 3.30e-6 | 3.02e-6 | 4.3175e-6 | 23.6% |
| 5 | 3.48e-6 | 3.31e-6 | 4.3175e-6 | 19.4% |
| 2.5 | 3.33e-6 | 3.26e-6 | 4.3175e-6 | 22.9% |

전역 검증은 정상: 합력 +4e-5 정확, umax 2.257e-4 = 해석해 2e-4와 정합
(홀 때문 compliance 증가분 포함), 경계 피크 ≥ 체적 피크.

### 원인 분석 (증거 기반)

1. **Small-cut 무안정화가 주범.** Crown zone 절단셀은 extension 없이
   FCM α(=1e-8 ghost OUT)만으로 버팀 — Phase 1 범위(v2 §7.6)가 허용한 한계.
   AgFEM extension은 Phase 2 범위. h=5→2.5에서 오차 정체(rate 0.02,
   비단조)는 분해능 부족이 아니라 안정화 부재의 signature.
2. 균일 격자 r/h=4는 Kt 포착에 성김. THB 국소세분은 Phase 3 범위.
3. finest에서 fallback 48셀 (octree 경로, 정밀도 하위). Crown 인접 fallback
   여부는 추가 추적 필요하나 1순위 원인은 아님.
4. moment fitting 자체는 정확 (GEO G1: 1e-12), 적분 오차 아님.

### 결론

V02 게이트는 Phase-1 수단으로 닫히지 않음. **사용자 결정 (2026-11-14):
V02-부분으로 G1 조건부 통과 인정 + Phase 2 착수.**
V02는 열어두고 Phase 2 (extension) → Phase 3 (THB) 후 같은 게이트로
재측정, 통과 시 G1 소급 완료.

> **후일담 (2026-11-15 G2 COMPLETE): 위 원인 분석 1번은 부분적으로 정정됨.
> G1 당시 주원인으로 지목한 small-cut 외에도 BC 오류(자유 모서리 uy=0, 양 z면 uz=0)가
> 함께 작용했음이 Phase 2에서 확정 — 좌측면 pinned 단독 BC로 0.674% PASS
> (`docs/phase2/G2_report.md` §T1/T3, `G2_COMPLETE.md`). 본 절은 당시 기록으로 보존.**

## moment fitting G1 — PASS

GEO S4 tet-cut q=2: 27개 Legendre 모멘트 전부 1e-12, Σw 1e-12, Roark
end-to-end Σw = tess 체적 (rel 1.2e-16), fallback 0 (h=10/5; h=2.5는 48).

## 통합 수정 4건 (기록, Phase 1 범위 내)

1. `solver/db.py` cell_id 전부 0 버그 → `(k*ny+j)*nx+i` (SOLVE 소유 파일,
   TRACK_SOLVE에 통보 요).
2. `solver/pipeline.py` ID 공간 충돌 (1e5 DOF↑): 가상 노드 `n_ctrl+1`부터,
   가상 요소 `N_elem+1`부터.
3. `solver/trimmed.py` (신규, 통합 담당): GEO QP 테이블 → `run_linear`
   `qp_table` hook, 경계 피크 probe (hole CAD면 직접 샘플 + spline 미분),
   빈 QP 셀 → ghost-α OUT 처리.
4. V02 압력 부호 (+1e-6 = tension), V03 압력 부호 (− = inward, Coreform 일치).

## 산출물

- `WHT_IGA/scripts/v02_roark.py`, `v03_beam.py` (재현 가능)
- `WHT_IGA/docs/phase1/evidence/v02|v03/` (results.json + DB/VTU/probe)
- `WHT_IGA/wht_iga/solver/trimmed.py`, pipeline/db surgical edits
