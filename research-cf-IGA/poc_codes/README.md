# poc_codes — Immersed IGA 핵심 알고리즘 PoC 스위트

master-plan_v1.md의 오픈소스 스택과 papers/ 문헌을 코드로 연결하는 최소 실행형 근거.
전제: `numpy/scipy/matplotlib/h5py`만 사용 (splinepy·Kratos 미설치 환경에서도 전부 통과).

| 폴더 | 알고리즘 | 문헌 | 오픈소스 대응 |
|---|---|---|---|
| 01_bezier_extraction | Bézier Extraction, K^e=C^T K_B C, 패치테스트 | P01, P18 (Borden 2011) | splinepy C^e 산출 대체 |
| 02_cut_quadrature | 모멘트피팅 vs 높이함수 축소 적분 | P05, P07 (Saye) | Algoim/QUGaR 적분 백엔드 |
| 03_ghost_agfem | 소형셀 조건수, Ghost Penalty, AgFEM 집합구속 | P11, P12, P13 | ngsxfem, GridapEmbedded.jl |
| 04_nitsche_bc | Nitsche 약형식 Dirichlet, γ 민감도 | P11, P20 계열 | ngsxfem Nitsche |
| 05_thb_refinement | Knot insertion, 계층기저 분할단위 | P18, P04 | G+Smo gsHSplines, thbsplines |
| 06_data_bus_mdpa | mesh_data.sqlite BLOB 버스, mdpa 직렬화 | master-plan §2.1 | rusqlite/SQLite, Kratos mdpa |
| 07_probe_post | 스플라인 역매핑, hdf5/hvascii 후처리 | master-plan Track A/B/C | splinepy 역매핑, h5py/VTKHDF |
| 08_verification_mini | Kirsch Kt=3, 보 처짐 (수집 CSV 대조) | Coreform 검증매뉴얼 | user-guide 수집 데이터 |

실행: 각 폴더에서 `python test_*.py` (assert 실패 시 PoC 깨짐).
공통 한계: 전부 1D/2D 축소 모델. 3D B-Rep·Kratos 연동은 본개발 단계.

## 실행 결과 (2026-11-14, numpy/scipy만, 전수 PASS)

| PoC | 핵심 수치 |
|---|---|
| 01 | 재구성 7.8e-16, congruence 2.6e-10, Greville 재현 1.1e-16 |
| 02 | 피팅 4차 정확·nnz 15, 높이함수 면적 2.9e-05 (접선스트립은 Tanh-Sinh 필요) |
| 03 | η 0.5→1e-6: 미안정 45→3.5e6, Ghost ~1e2, AgFEM ~29 |
| 04 | γ<100 비강압성, γ=2000 회복, γ=20000 조건수 단조악화 (γ≳1/η 법칙) |
| 05 | 삽입 불변 1.1e-16, naive HB 0.33 vs truncated 4.4e-16 |
| 06 | sqlite BLOB·mdpa 라운드트립 일치 |
| 07 | 역매핑 5.6e-16, hdf5/hvascii 일치 |
| 08 | Kt=3.0 정확, 수집CSV 수렴율 1.82·최종오차 0.43% |
