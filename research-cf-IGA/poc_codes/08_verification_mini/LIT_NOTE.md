# 08 — 미니 검증: Kirsch + 외팔보 (LIT_NOTE)

## 문헌 상관
- Coreform 검증매뉴얼 수집본:
  `user-guide/verification-manual/problems/solid_mechanics/linear_elastic_stress/`
  `infinite-plate-with-hole/infinite-plate-with-hole.md` (Kirsch, Kt=3.0),
  `beam-bending/beam-bending.md` + `generated-downloads/.../beam-bending/*tip-deflection.csv`
  (하중별 끝단처짐 vs Euler-Bernoulli/ABAQUS 기준값).
- master-plan 테스트 V-02/V-04의 축소 선행판.
- P05 적분 정확도가 Kt 복원율을 좌우 (수렴률 O(h^{p+1})).

## 오픈소스 상관
- 수집 CSV 자체가 검증 데이터 버스 역할 (전처리 없이 바로 대조).

## 핵심 내용 (본 PoC)
1. Kirsch 해석해: σθθ=σ(1−2cos2θ) → 정점 Kt=3.0, 구멍 근방 감쇠 확인.
2. 수집 traction-tip-deflection.csv 로드 → 미세화별 Euler-Bernoulli 오차
   단조감소 + 최종 <5% + 수렴율 추산.

## 본 PoC 범위/한계
- 해석해·데이터 대조만 수행 (실제 Immersed 해석은 본개발 V-02/V-04).
