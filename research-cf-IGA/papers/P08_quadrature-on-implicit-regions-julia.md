# P08 — QuadratureOnImplicitRegions (Julia package)

- 저자: Julia 전산역학 커뮤니티
- 연도/게재: 2023 / Julia Packages Registry
- 출처: Julia Packages
- 확보상태: 메타데이터만 + 코드 공개 — Pkg.add로 확보 가능

## 핵심 수치 기여

Julia JIT 환경에서 레벨셋/암시적 영역 고차 적분룰 자동 생성, 비적합 해석 라이브러리 연동 인터페이스

## 핵심 구현 내용 추출

- Saye류 축소법을 Julia로 래핑, 다중 디스패치로 기하 유형별 분기.
- GridapEmbedded 등과의 파이프라인 예제 포함.

## Coreform IGA 연관성

빠른 프로토타이핑용 — 아이디어 검증 후 C++/Coreform 경로로 이식.

## 재현 체크리스트

- [ ] Julia 재현 스크립트 버전 고정
- [ ] 적분룰 차수 감도