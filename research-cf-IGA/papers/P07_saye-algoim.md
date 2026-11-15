# P07 — Algoim: Algorithms for implicitly defined geometry, level set methods, and Voronoi implicit interface methods

- 저자: Robert I. Saye (LBNL/UC Berkeley)
- 연도/게재: 2015/2022 / SISC 2015; JCP 2022 + C++ 오픈소스 (algoim.github.io)
- DOI: 10.1137/140966290; 10.1016/j.jcp.2021.110720
- 출처: Algoim 공식 포털
- 확보상태: 메타데이터만 + 코드 공개 — 깃허브 클론으로 확보 가능

## 핵심 수치 기여

초직사각형 내 레벨셋 절단영역을 다가지 높이함수 그래프로 재구성 → 재귀 차원축소 → 1D 적분 텐서곱, 고차 가중치 기계정밀도; deal.II·PETSc 백엔드 채택

## 핵심 구현 내용 추출

- 셀 내 레벨셋 부호로 위상 분류 → 높이함수 방향 자동 선택.
- 특이 방향은 재귀 세분 후 축소 반복.
- 출력: 점 좌표+가중치 테이블 → 기존 FEM 조립 루프에 그대로 주입.

## Coreform IGA 연관성

절단셀 적분의 사실상 표준 라이브러리 — Coreform 메시 적분과 상호 검증용.

## 재현 체크리스트

- [ ] 레벨셋 구/토러스 적분 오차
- [ ] deal.II 튜토리얼 스텝 재현