# P21 — Isogeometric and immersogeometric analysis of incompressible flow problems (PhD)

- 저자: Tuong Hoang (Verhoosel/van Brummelen/Reali)
- 연도/게재: 2018 / TU/e & IUSS PhD, ISBN 978-90-386-4449-3
- 출처: TU/e Pure
- 확보상태: 메타데이터만 — 학위논문 공개본 확인 필요

## 핵심 수치 기여

비압축성 NS/Stokes Immersogeometric: 속도-압력 동일 고차 B-스플라인 + 스켈레톤 안정화 + 속도 Ghost Penalty, 다공성 유동 실증

## 핵심 구현 내용 추출

- 스켈레톤 안정화: 요소 골격 위 고차도함수 도약 penalize → inf-sup 우회.
- 등차수(Equal-order) 속도-압력 공간 사용 가능.
- 다공성: 고체 골격은 FCM 지시함수로 처리.

## Coreform IGA 연관성

유동·다공성 확장 시 안정화 처방.

## 재현 체크리스트

- [ ] 원주 유동 벤치마크
- [ ] 압력 진동 여부