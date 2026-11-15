# P06 — Quadrature-free immersed isogeometric analysis

- 저자: Antolin, Hirschler
- 연도/게재: 2022 / Eng. with Computers 38(5):4475–4499 / arXiv:2107.09024
- DOI: 10.1007/s00366-022-01644-3
- arXiv: 2107.09024
- 출처: arXiv + ResearchGate 360184654 + EPFL Infoscience
- 확보상태: FULL-TEXT 확보 — pdf/2107_09024.pdf (24pp, ~12.2MB, arXiv 검증: 2021-07-19)

## 핵심 수치 기여

절단셀 적분점 제로: Galerkin 텐서를 일관 다항식 공간으로 변환 → 발산정리 2회 적용 → 3D 체적=2D B-Rep 면적분=1D 경계 선적분, 닫힌해 기계정밀도

## 핵심 구현 내용 추출

- 피적분 다항식 차수 p에서 체적항을 경계 선적분으로 축소하는 변환 테이블 사전 구축.
- B-Rep 트림 곡면의 파라미터 경계 루프를 따라 1D Gauss 적분만 수행.
- 구현 조건: 커널이 다항식이어야 하며, 비다항식 재료·하중은 투영 후 적용.

## Coreform IGA 연관성

Coreform 무적분 경로의 이론 원형 — 적분점 폭발이 문제될 때 대안.

## 재현 체크리스트

- [ ] 구/원기둥 체적·질량행렬 기계정밀도 테스트
- [ ] 다항식 차수별 변환 정확도