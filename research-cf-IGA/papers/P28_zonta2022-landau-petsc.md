# P28 — The Landau Collision Integral in the Particle Basis in the PETSc Library

- 저자: Zonta, Pusztay, Hirvijoki, Adams, Knepley
- 연도/게재: 2022/2023 / Phys. Plasmas 29(12):123906
- DOI: 10.1063/5.0105182
- 출처: ResearchGate 371786455
- 확보상태: 메타데이터만 — 유료 저널

## 핵심 수치 기여

Landau 충돌적분자의 보존형(FE/입자) 이산화 + Algoim 다변수 다항식 암묵영역 적분을 PETSc 백엔드 통합

## 핵심 구현 내용 추출

- 에너지·운동량·양성 보존 이산화.
- Algoim 적분 커널을 PETSc 분산 조립에 연결.

## Coreform IGA 연관성

Algoim이 플라즈마까지 쓰인다는 성숙도 증거.

## 재현 체크리스트

- [ ] 보존량 장시간 추적