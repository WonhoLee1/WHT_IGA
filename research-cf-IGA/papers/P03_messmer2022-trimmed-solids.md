# P03 — Efficient CAD-integrated isogeometric analysis of trimmed solids

- 저자: Meßmer, Teschemacher, Leidinger, Wüchner, Bletzinger
- 연도/게재: 2022 / CMAME 400, 115584
- DOI: 10.1016/j.cma.2022.115584
- 출처: TUM mediaTUM 1767547
- 확보상태: 메타데이터만 — 유료 저널

## 핵심 수치 기여

B-Rep 솔리드 → 트리밍 3변량 B-스플라인 배경격자 매핑 파이프라인, QuESo 전처리, 선형 FEM 대비 12배 적은 자유도로 동등 정확도

## 핵심 구현 내용 추출

- QuESo(Kratos) 전처리: B-Rep 면과 배경셀 교차 → 절단셀 적분점/가중치 생성.
- 트리밍된 체적은 경계 B-Rep 서피스 적분으로 체적항을 변환하는 발산정리 경로와 병행 가능.
- 자유도당 효율 지표(변형에너지 오차/자유도)로 검증.

## Coreform IGA 연관성

B-Rep 직접 입력 파이프라인의 실용 구현 레퍼런스.

## 재현 체크리스트

- [ ] 동일 B-Rep에서 자유도-오차 곡선 재현
- [ ] QuESo 버전·파라미터 기록