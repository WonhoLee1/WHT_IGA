# P04 — Scan-based immersed isogeometric analysis (PhD thesis)

- 저자: Sai Chandana Divi (Verhoosel/Auricchio/Reali/van Brummelen)
- 연도/게재: 2021/2022 / TU/e & Univ. Pavia PhD thesis, ISBN 978-90-386-5469-0
- 출처: Univ. Pavia IRIS
- 확보상태: 메타데이터만 — 학위논문 공개본 기관 저장소 확인 필요

## 핵심 수치 기여

복셀 스캔→절단셀 Immersed IGA, Strang 제1보조정리 기반 오차추정 적응 옥트리, 오일러지표 위상오류 방지 THB 세분화

## 핵심 구현 내용 추출

- 스캔 복셀 레벨셋을 배경 THB-스플라인 격자에 침강.
- 셀별 오차지표 기반 옥트리 적분 깊이 적응 조절.
- 필터링으로 생기는 고립 위상 성분은 오일러지표 검사로 제거 후 해석.

## Coreform IGA 연관성

이미지 기반 해석(의료·적층제조 CT)의 이론 근거.

## 재현 체크리스트

- [ ] 복셀 해상도별 수렴 스터디
- [ ] 위상검사 on/off 비교