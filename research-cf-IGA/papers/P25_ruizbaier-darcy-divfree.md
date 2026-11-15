# P25 — Divergence-free unfitted finite element discretisations for the Darcy problem

- 저자: Ruiz-Baier, de Prenter, Larson et al.
- 연도/게재: 2024 / Monash/Authors Monograph
- 출처: 저자 저장소
- 확보상태: 메타데이터만 — 모노그래프 공개본 확인 필요

## 핵심 수치 기여

Darcy 문제 발산보존 비적합 이산화: H(div) 적합공간+AgFEM 확장으로 무진동 체적보존

## 핵심 구현 내용 추출

- Raviart-Thomas류 H(div) 공간을 절단셀에 확장.
- 압력-유속 inf-sup은 AgFEM 구속과 병행.

## Coreform IGA 연관성

다공성·침투류 체적보존 처방.

## 재현 체크리스트

- [ ] 발산오차 기계정밀도 테스트