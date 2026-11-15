# P26 — Isogeometric Topology Optimization Based on Topological Derivatives

- 저자: Teixeira, Krenn, Gangl, Marussig
- 연도/게재: 2025 / Struct. Multidisc. Optim. / arXiv:2509.09236
- arXiv: 2509.09236
- 출처: arXiv
- 확보상태: FULL-TEXT 확보 — pdf/2509_09236.pdf (19pp, ~1.5MB, arXiv 검증: 2025-09-11)

## 핵심 수치 기여

위상도함수+IGA 배경격자 재격자 없는 위상최적화: Algoim 체적비 산출 + 절단강성 완충전이

## 핵심 구현 내용 추출

- 위상도함수 필드로 구멍 생성 위치 결정 → 레벨셋 업데이트.
- 절단셀 강성은 체적비 η의 완충함수 E(η)로 전이(0/1 불연속 완화).
- Algoim으로 η 고정밀 산출이 수렴 핵심.

## Coreform IGA 연관성

Coreform 결과→경량설계 연계 시 방법론.

## 재현 체크리스트

- [ ] 캔틸레버 컴플라이언스 벤치마크