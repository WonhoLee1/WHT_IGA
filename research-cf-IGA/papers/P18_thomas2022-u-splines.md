# P18 — U-splines: splines over unstructured meshes

- 저자: Thomas, Engvall, Schmidt, Tew, Scott
- 연도/게재: 2022 / CMAME 403, 115747 / Coreform·BYU
- DOI: 10.1016/j.cma.2022.115747
- 출처: Coreform 기술자료
- 확보상태: 메타데이터만 — 유료 저널, Coreform 라이브러리 요약본 확인 권장

## 핵심 수치 기여

비구조격자 U-스플라인: h·p·연속성 국소 가변, T-접합·혼합셀 지원, 국소 영공간 해법으로 분할단위·양성·국소선형독립 증명

## 핵심 구현 내용 추출

- 셀별 요구 연속성(r,θ) → 베지어 계수 구속식 → 국소 영공간 기저 추출.
- 추출행렬은 희소 저장, 요소 루프는 Bernstein 루틴 재활용(P01 연계).
- 산업 B-Rep → 해석적합 U-스플라인 변환이 Coreform 전처리 핵심.

## Coreform IGA 연관성

Coreform IGA 기저 기술 그 자체 — 최우선 정독 대상.

## 재현 체크리스트

- [ ] 단위분할·양성 수치 검증
- [ ] 특이점 패치테스트