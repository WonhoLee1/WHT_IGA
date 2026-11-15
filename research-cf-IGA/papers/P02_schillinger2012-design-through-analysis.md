# P02 — An Isogeometric design-through-analysis methodology based on adaptive hierarchical refinement of NURBS, immersed boundary methods, and T-spline CAD surfaces

- 저자: Schillinger, Dedè, Scott, Evans, Borden, Rank, Hughes
- 연도/게재: 2012 / CMAME 249–252, pp.116–150
- DOI: 10.1016/j.cma.2012.03.017
- 출처: ResearchGate Pub.256918134
- 확보상태: 메타데이터만 — 유료 저널, 기관 접근 필요

## 핵심 수치 기여

계층적 B-스플라인 + FCM 결합, T-스플라인 산업 CAD(프로펠러·휠) 직접 해석, Design-through-Analysis 최초 정립

## 핵심 구현 내용 추출

- 계층적 세분화: 트리 구조로 경계 피처 셀만 국소 세분, 선형독립성 유지.
- FCM 지시함수: 물리영역 내 α=1, 허구영역 α=ε(≈1e-8~1e-10)로 약재료 처리 후 표준 적분.
- T-스플라인 표면을 배경셀과 교차 판정하여 내/외/절단 분류.

## Coreform IGA 연관성

Coreform 계열(Scott→Coreform)의 원류 — 트리밍된 CAD 직접 해석 철학의 출발점.

## 재현 체크리스트

- [ ] 허구영역 α 감도 스터디
- [ ] 세분화 레벨별 변형에너지 수렴 곡선