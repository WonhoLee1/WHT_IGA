# P23 — Numerical Analysis of a Cut Finite Element Approach for Fully Eulerian Fluid-Structure Interaction with Fixed Interface

- 저자: Frei, Knoke, Steinbach, Wenske, Wick
- 연도/게재: 2026 / arXiv:2603.25279 [math.NA]
- arXiv: 2603.25279
- 출처: arXiv
- 확보상태: FULL-TEXT 확보 — pdf/2603_25279.pdf (31pp, ~0.9MB, arXiv 검증: 2026-03-26)

## 핵심 수치 기여

유체+구조 단일 오일러 배경격자 CutFEM 완전오일러 FSI: 체적절단비 적응 가중 Ghost Penalty로 계면접촉·연성 강건화

## 핵심 구현 내용 추출

- 가중 페널티: 절단비 η→0 셀로 갈수록 가중 상향.
- 고정계면 가정 하 오차수렴 분석 → 이동계면 확장 전제.
- 유체·고체 서로 다른 구성식을 동일 격자 위 Heaviside 분기로 처리.

## Coreform IGA 연관성

FSI 확장 시 오일러 기술 경로.

## 재현 체크리스트

- [ ] 고정계면 수렴율 재현