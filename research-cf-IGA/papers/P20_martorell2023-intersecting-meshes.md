# P20 — Intersecting surface meshes and Cartesian background meshes for embedded simulations

- 저자: Martorell, Badia, Verdugo
- 연도/게재: 2022/2023 / UPC UPCommons Report
- 출처: UPCommons
- 확보상태: 메타데이터만 — 기관 리포트 공개본 확인 필요

## 핵심 수치 기여

STL 표면메시×데카르트 배경격자 교차의 폴리토프 셀 분할·2단계 적분메시 체계화

## 핵심 구현 내용 추출

- 삼각형-셀 교차 → 폴리토프 분해 → 테트라 분할 후 표준 적분.
- 2단계: 1차 위상분류 메시 + 2차 고정밀 적분메시.

## Coreform IGA 연관성

STL 직접 입력 경로의 기하 커널.

## 재현 체크리스트

- [ ] 복잡 STL 교차 견고성 테스트