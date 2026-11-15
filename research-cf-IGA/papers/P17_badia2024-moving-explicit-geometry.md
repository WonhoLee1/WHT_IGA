# P17 — Space-time unfitted finite elements on moving explicit geometry representations

- 저자: Badia, Martorell, Verdugo
- 연도/게재: 2024 / JCP 511, 113127 / arXiv:2401.12649
- DOI: 10.1016/j.jcp.2024.113127
- arXiv: 2401.12649
- 출처: arXiv comp-phys
- 확보상태: FULL-TEXT 확보 — pdf/2401_12649.pdf (19pp, ~7.1MB, arXiv 검증: 2024-01-23)

## 핵심 수치 기여

명시적 다면체 B-Rep 이동계면에 시공간 비적합 FEM + 고차 안정화 파이프라인

## 핵심 구현 내용 추출

- B-Rep↔배경격자 교차의 폴리토프 분할 → 2단계 적분메시.
- 이동 시 위상변화(합병/분리) 이벤트의 슬랩 경계 처리.

## Coreform IGA 연관성

STL/B-Rep 직접 입력 동적 문제의 기하 커널.

## 재현 체크리스트

- [ ] 이동계면 질량보존 오차
- [ ] 위상변화 케이스