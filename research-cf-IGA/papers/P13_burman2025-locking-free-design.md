# P13 — On the design of locking free ghost penalty stabilization and the relation to CutFEM with discrete extension

- 저자: Burman, Hansbo, Larson
- 연도/게재: 2025/2026 / Numer. Math. 158:249–275 / arXiv:2205.01340
- DOI: 10.1007/s00211-025-01502-6
- arXiv: 2205.01340
- 출처: Springer OA + arXiv
- 확보상태: FULL-TEXT 확보 — pdf/2205_01340.pdf (31pp, ~4.6MB, arXiv 검증: 2022-05-03)

## 핵심 수치 기여

이산확장 CutFEM 공간–페널티 커널 관계 일반화, 임의 고차에서도 잠김 없는 함수형(DOF기반) Ghost Penalty 설계 완성

## 핵심 구현 내용 추출

- 셀 내부 DOF 기반 벌점(함수형)으로 설계 자유도 확보.
- 차수 p별 스케일링 h^{2k-1} 유지, 상수 γ의 p 의존성 표기.
- 이산확장 연산자와의 동치성으로 안정성 증명 재사용.

## Coreform IGA 연관성

고차(p≥3) 안정화 설계 지침.

## 재현 체크리스트

- [ ] p=1..4 조건수-오차 매트릭스