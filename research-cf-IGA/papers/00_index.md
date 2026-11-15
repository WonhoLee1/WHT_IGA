# Papers — 문헌 확보 현황 및 핵심 구현 인덱스

출처: ../문헌종합.md (번호 P01–P30은 해당 문서의 1–30번과 일치)

## 확보 요약

- FULL-TEXT 9편: P05, P06(=P29, P30), P12, P13, P16, P17, P23, P24, P26 → `pdf/` (arXiv API 검증済)
- 메타데이터+공개코드: P07(Algoim), P08(Julia)
- 메타데이터만(유료·기관저장소): P02, P03, P04, P09, P10, P11, P14, P15(2차), P18, P19, P20, P21, P22(OA 재시도), P25, P27, P28
- P29/P30은 P06과 동일 원문이므로 별도 파일 없이 P06 PDF로 커버

## 파일 목록

- P01 `immersed-iga-deep-report` — 2026 — 메타데이터만 — 원문은 프로젝트 내부 합성 보고서, 별도 PDF 없음
- P02 `schillinger2012-design-through-analysis` — 2012 — 메타데이터만 — 유료 저널, 기관 접근 필요
- P03 `messmer2022-trimmed-solids` — 2022 — 메타데이터만 — 유료 저널
- P04 `divi2021-scan-based` — 2021/2022 — 메타데이터만 — 학위논문 공개본 기관 저장소 확인 필요
- P05 `loibl2026-quadrature-comparison` — 2026 — FULL-TEXT 확보 — pdf/2602_18130.pdf (52pp, ~9.3MB, arXiv 검증: 2
- P06 `antolin2022-quadrature-free` — 2022 — FULL-TEXT 확보 — pdf/2107_09024.pdf (24pp, ~12.2MB, arXiv 검증: 
- P07 `saye-algoim` — 2015/2022 — 메타데이터만 + 코드 공개 — 깃허브 클론으로 확보 가능
- P08 `quadrature-on-implicit-regions-julia` — 2023 — 메타데이터만 + 코드 공개 — Pkg.add로 확보 가능
- P09 `thesis-curved-geometry` — 미상 — 미확보 — 서지정보 불완전(저자·연도 미상), ProQuest 접근 필요
- P10 `boundary-conformal-finite-elasticity` — 2023 — 메타데이터만 — PMC ID로 전문 확보 가능 (이번 배치 미수집)
- P11 `burman2010-ghost-penalty` — 2010 — 메타데이터만 — Numdam OA로 확보 가능 (이번 배치 미수집)
- P12 `badia2022-linking-ghost-agfem` — 2022 — FULL-TEXT 확보 — pdf/2106_13728.pdf (23pp, ~9.5MB, arXiv 검증: 2
- P13 `burman2025-locking-free-design` — 2025/2026 — FULL-TEXT 확보 — pdf/2205_01340.pdf (31pp, ~4.6MB, arXiv 검증: 2
- P14 `gross2020-cutfem-preconditioner` — 2020 — 메타데이터만 — 기관 리포트 공개본 확인 필요
- P15 `ghost-penalty-survey` — 2026 — 2차 종합자료 — 원문 대체 불가, 입문용으로만 사용
- P16 `badia2022-space-time-agfem` — 2022 — FULL-TEXT 확보 — pdf/2206_03626.pdf (23pp, ~1.2MB, arXiv 검증: 2
- P17 `badia2024-moving-explicit-geometry` — 2024 — FULL-TEXT 확보 — pdf/2401_12649.pdf (19pp, ~7.1MB, arXiv 검증: 2
- P18 `thomas2022-u-splines` — 2022 — 메타데이터만 — 유료 저널, Coreform 라이브러리 요약본 확인 권장
- P19 `thomas2024-extraordinary-spectrum` — 2024 — 메타데이터만 — 유료 저널
- P20 `martorell2023-intersecting-meshes` — 2022/2023 — 메타데이터만 — 기관 리포트 공개본 확인 필요
- P21 `hoang2018-immersogeometric-flow` — 2018 — 메타데이터만 — 학위논문 공개본 확인 필요
- P22 `rim2025-enriched-ibm-cht` — 2025 — 메타데이터만 — OA 저널이나 이번 배치 PDF 확보 실패(500), 재시도 필요
- P23 `frei2026-eulerian-fsi` — 2026 — FULL-TEXT 확보 — pdf/2603_25279.pdf (31pp, ~0.9MB, arXiv 검증: 2
- P24 `cahaly2025-pcic-vof` — 2026 — FULL-TEXT 확보 — pdf/2506_18811.pdf (33pp, ~13.1MB, arXiv 검증: 
- P25 `ruizbaier-darcy-divfree` — 2024 — 메타데이터만 — 모노그래프 공개본 확인 필요
- P26 `teixeira2025-iga-topopt` — 2025 — FULL-TEXT 확보 — pdf/2509_09236.pdf (19pp, ~1.5MB, arXiv 검증: 2
- P27 `bath2025-qmc-helmholtz` — 2025 — 메타데이터만 — 기관 모노그래프 확인 필요
- P28 `zonta2022-landau-petsc` — 2022/2023 — 메타데이터만 — 유료 저널
- P29 `antolin2021-qfree-preprint` — 2021 — P06과 동일 원문 — pdf/2107_09024.pdf로 확보済, 별도 파일 미생성
- P30 `antolin2022-qfree-infoscience` — 2022 — P06과 동일 원문 — pdf/2107_09024.pdf로 확보済, 별도 파일 미생성

## 구현 스택별 읽기 순서 (상세는 implementation_synthesis.md)
- 기하입력: P02 → P03 → P20 → P04
- 적분: P07 → P05 → P06
- 안정화: P11 → P12 → P13 → P14(+P16 시공간)
- 기저: P18 → P19
- 다중물리: P21 → P22 → P23 → P25 → P24
- 최적화·UQ: P26 → P27 → P28
- 시공간 이동계면: P16 → P17