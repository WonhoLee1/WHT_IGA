# P12 — Linking ghost penalty and aggregated unfitted methods

- 저자: Badia, Neiva, Verdugo
- 연도/게재: 2022 / CMAME 388, 114232 / arXiv:2106.13728
- DOI: 10.1016/j.cma.2021.114232
- arXiv: 2106.13728
- 출처: arXiv math.NA
- 확보상태: FULL-TEXT 확보 — pdf/2106_13728.pdf (23pp, ~9.5MB, arXiv 검증: 2021-06-19)

## 핵심 수치 기여

Ghost↔AgFEM 수학적 연결 최초 규명: 페이스기반 페널티의 무한대 극한 잠김 규명 → Locking-Free Ghost Penalty 정립

## 핵심 구현 내용 추출

- 집합(aggregate) 구속: 불량 절단 DOF를 내부셀 DOF의 외삽으로 제거(강한 안정화).
- 새 페널티는 해-집합확장 거리 제어 → γ→∞에서도 잠김 없음.
- 구현: 셀 집합 그래프 구축 → 구속행렬 적용 → 축소 시스템 해석.

## Coreform IGA 연관성

γ 튜닝 부담을 없애는 실무 해법.

## 재현 체크리스트

- [ ] γ→∞ 극한 수렴 테스트
- [ ] 집합 패턴 시각화