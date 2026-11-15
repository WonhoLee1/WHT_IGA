# P24 — PCIC: Cylindrical Volume Moment Calculation and Interface Reconstruction for Sub-Grid Modeling in VOF Methods

- 저자: Cahaly, Wasquel, Zou, Desjardins, Evrard
- 연도/게재: 2026 / JCP / arXiv:2506.18811
- arXiv: 2506.18811
- 출처: arXiv comp-phys
- 확보상태: FULL-TEXT 확보 — pdf/2506_18811.pdf (33pp, ~13.1MB, arXiv 검증: 2025-06-23)

## 핵심 수치 기여

VOF 아격자 계면재구성용 원통형 체적모멘트 고정밀 산출 PCIC

## 핵심 구현 내용 추출

- 원통 좌표 체적모멘트의 해석식+수치적분 하이브리드.
- 아격자 계면 다항식 재구성 → 곡률·법선 정밀도 향상.

## Coreform IGA 연관성

다상유동 계면모멘트 계산의 기하 커널(연계 참고).

## 재현 체크리스트

- [ ] 구/원통 모멘트 해석해 대조