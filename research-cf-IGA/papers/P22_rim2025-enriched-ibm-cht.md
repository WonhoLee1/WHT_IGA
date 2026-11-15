# P22 — Enriched immersed boundary method for convective conjugate heat transfer with level set-based surface reconstruction

- 저자: Rim, Zhao, Yan
- 연도/게재: 2025 / AIMS ACSE 3:1–23
- DOI: 10.3934/acse.2025002
- 출처: AIMS Press
- 확보상태: 메타데이터만 — OA 저널이나 이번 배치 PDF 확보 실패(500), 재시도 필요

## 핵심 수치 기여

열유속 양립성 문제 해결: 절단셀 온도 자유도를 유체/고체에 복제 부여(EIBM) + 레벨셋 표면재구성

## 핵심 구현 내용 추출

- 절단셀 Heaviside 풍부화: T = H·T_f + (1-H)·T_s 형태 자유도 복제.
- 계면에 Nitsche/페널티로 온도·유속 연속성 약부과.
- 레벨셋 재구성 표면에서 열전달계수 평가.

## Coreform IGA 연관성

공액열전달(CHT) 문제의 절단셀 처방.

## 재현 체크리스트

- [ ] 평판 대류 벤치마크 Nusselt 비교