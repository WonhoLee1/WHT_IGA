# P05 — Comparative study of different quadrature methods for cut elements

- 저자: Loibl, Teixeira, Toprak, Shishkina, Miao, Kiendl, Marussig, Kummer
- 연도/게재: 2026 / Arch. Comput. Methods Eng. / arXiv:2602.18130
- arXiv: 2602.18130
- 출처: arXiv math.NA
- 확보상태: FULL-TEXT 확보 — pdf/2602_18130.pdf (52pp, ~9.3MB, arXiv 검증: 2026-02-20)

## 핵심 수치 기여

Algoim·BoSSS·FCMLab·QUGaR·Gridap/STLCutters·mlhp·ngsxfem·Nutils 망라 2D/3D 벤치마크: 정확도·적분점수·효율 비교, 모멘트피팅/고차매핑/높이함수 축소 정량 비교

## 핵심 구현 내용 추출

- 벤치마크 축: 기하 난이도(곡률·박판·특이점) × 요구 정확도 × 적분점 수.
- 모멘트피팅: 단항식 모멘트 일치 가중치 해법, 점 수는 적지만 조건수 관리 필요.
- 높이함수 축소(Saye): 재귀 차원축소 + 1D Gauss-Legendre/Tanh-Sinh 텐서곱, 기계정밀도 경로.
- 실무 선택표: 곡면 B-Rep이면 축소법, 다면체 근사면이면 옥트리/피팅.

## Coreform IGA 연관성

적분 백엔드 선택의 근거표 — Coreform 메시 생성 옵션 해석에 직접 활용.

## 재현 체크리스트

- [ ] 2D 원/3D 구 벤치마크 재현
- [ ] 적분점수-오차 파레토 곡선