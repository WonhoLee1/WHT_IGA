# P01 — Immersed IGA 심층 연구 보고서 (Bézier Extraction, FCM/CutFEM/IBM 비교, Algoim/무적분/Ghost/AgFEM, 오픈소스 생태계)

- 저자: Gemini Notebook Research Synthesis
- 연도/게재: 2026 / Technical Report (프로젝트 지식베이스)
- 출처: 프로젝트 지식베이스
- 확보상태: 메타데이터만 — 원문은 프로젝트 내부 합성 보고서, 별도 PDF 없음

## 핵심 수치 기여

기하-해석 분리, Bézier Extraction으로 상용 솔버 재활용, FCM/CutFEM/IBM 비교, Algoim·무적분·Ghost·AgFEM·오픈소스(ngsxfem, GridapEmbedded, deal.II) 집대성

## 핵심 구현 내용 추출

- Bézier Extraction: NURBS 패치 = Bernstein 기저 × 추출행렬 C, 요소 강성은 표준 Bernstein 루틴으로 조립 후 C로 변환 (솔버 수정 불필요).
- 워크플로우: B-Rep/복셀 입력 → 배경 직교격자 침강 → 절단셀 분류 → 적분(Algoim/옥트리/무적분) → 안정화(Ghost/AgFEM) → 표준 조립/해석.
- 적응 세분화는 경계 교차셀 중심 트리 구조로 국소 해상도 확보.

## Coreform IGA 연관성

Coreform IGA 전체 아키텍처의 개념 지도 — 아래 P02~P30 각 모듈의 상위 설계도.

## 재현 체크리스트

- [ ] 추출행렬 C가 단위분할을 유지하는지 단위 테스트
- [ ] 동일 모델을 경계적합 격자 결과와 변형에너지 비교