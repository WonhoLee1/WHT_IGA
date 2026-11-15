# TRACK A — 체적·요소·구속 요약

- 범위: S02, S03, S04, S05, S06, S13, S15, S16 (트랙 A 소유 파일만 생성·수정)
- 브랜치: WHTIGA. 환경: Python 3.12.16, Kratos 10.4.3 (`E:\kratos-build\install`, 96bad120), 2026-11-14
- 실행: 각 spike 파일만 지정 실행 (`pytest WHT_IGA/tests/kratos_spikes/test_sNN_*.py`). 전체 suite 미실행(트랙 B 병행)

## 결과표

| spike | 결과 | 핵심 수치 (허용 기준) |
|---|---|---|
| S02 항등 사상 | 예 | GlobalCoordinates 4.44e-16, detJ 1.11e-15, 회전 3.33e-16/1.11e-15 (1e-14) |
| S03 셀 QP 기하 | 예 | 셀당 기하 1개·27점, 가중치합 0.125=셀체적, SFV행합 1, 지지 27노드 일치. 면 위 점은 왼쪽 셀 귀속 |
| S04 SM 요소 | 예 | K 상대오차 8.7e-17/1.1e-16/2.3e-16 (1e-12), 패치 3.1e-15 (1e-10), σ=(1,2,3,0,0,0), von Mises=√3, TL/UL 초기화 OK |
| S05 구속식 | 예 | TᵀKT 대비 4.2e-17 (1e-10), 복원 2.5e-18, 1e5 생성 2.7초 (30초), 조립+풀이 32초, master DofSet 포함 |
| S06 고유치 | 부분 | Kratos solver 실패(QP 질량행렬 0). scipy T-투영은 hex-Kratos와 ~1e-9, 구속 전후 1e-8 일치 |
| S13 혼합 요소 | 아니오 | MixedVolumetric 초기화 크래시 rc=3221225477 (결정적). MixedStrain은 초기화 OK |
| S15 테이블 | 부분 | 부착·보간 OK(1.2@20°C). 법칙이 E(T) 미반영. Properties 직접 갱신 failover ×2 스케일 확인 |
| S16 스케일 | 예 | 97336 요소 52.5초 (60초), +2.92GB (8GB). p=3 893개/초 |

## 경로 결정(R1/R2/R3)에 미치는 영향

- S03·S04·S05 전부 통과 → ADR-001 결정 트리상 **R1 = (a)+(c) (셀 단위 적분점 기하 + SM 요소 + 구속식 T)** 확정 제안. R2 불필요, R3 불필요
- S06 부분: 정적 R1에 영향 없음. Phase 4 고유치/동역학은 "Kratos K 추출 + 자체 M + scipy eigsh" failover로 진행 (v2 명시 대응). C++ 승격 사유 아님
- S13 아니오: v2 예상대로 → Phase 5 C++ MixedUP0 승격 1순위 유지. 정적 선형 경로 영향 없음
- S15 부분: Phase 4 열 관련은 Properties 직접 갱신 방식으로 진행
- S16 예: 1e5 규모 Python 생성 가능하나 52.5초로 목표 근접 → Phase 6 C++ 생성기 후보 유지

## Kratos 동작 특이점 (v2/요청서 대비 신규 발견)

1. **Dirichlet 상수 이중 계산**: 고정 master를 관계식에서 제외하고 c에 미리 접으면 slave = Ru + 2c. 고정 master는 관계식에 그대로 두고 c=0 (Kratos가 자체 접음). `constraints.py`에 규칙 기록, S05-A 4.2e-17로 검증
2. **BlockBuilder는 slave를 시스템에 잔류**시킴 (108×108, slave행 대각화). TᵀKT 소거형이 아님. 해는 정확하므로 R1 수학(TᵀKT 동등)은 scipy로 별도 검증됨
3. **QP 기하 consistent 질량행렬 = 0** (`CalculateMassMatrix`), lumped는 정상. EigensolverDynamicScheme·Bossak 동역학에 영향
4. **E(T) 테이블 미반영**: LinearElastic3DLaw·ThermalElasticIsotropic3D 모두 기본값 사용
5. **생성 순서**: `CreateNewMasterSlaveConstraint`는 master 리스트가 먼저. 0-크기 관계 행렬은 크래시
6. **추출 K 주의**: Dirichlet 포함 `GetSystemMatrix`는 고정행이 1로 치환됨. K 비교는 무구속 상태에서

## 공유 모듈 변경 제안

- 없음. 공유 모듈(`geometry.py`, `_bootstrap.py`, `grid/*`, `reference/*`, `conftest.py`, `pyproject.toml`, `uv.lock`) 수정 없이 완료
- 신규 모듈(트랙 A 소유): `wht_iga/kratos_backend/elements.py` (셀 요소 부착, 선형 풀이, CSR 추출, 순열), `wht_iga/kratos_backend/constraints.py` (tie 헬퍼, Dirichlet 규칙)
- 통합 담당 참고: `elements.system_matrix_csr` + `equation_permutation`은 G0 무구속 큐브 비교에 재사용 가능

## 생성·수정 파일

- `WHT_IGA/wht_iga/kratos_backend/elements.py` (신규)
- `WHT_IGA/wht_iga/kratos_backend/constraints.py` (신규)
- `WHT_IGA/tests/kratos_spikes/test_s02_identity.py`, `test_s03_cell_quadrature.py`, `test_s04_sm_elements.py`, `test_s05_constraints.py`, `test_s06_eigen.py`, `test_s13_mixed.py`, `test_s15_tables.py`, `test_s16_scaling.py` (신규)
- `WHT_IGA/docs/spikes/S02.md`, `S03.md`, `S04.md`, `S05.md`, `S06.md`, `S13.md`, `S15.md`, `S16.md` (신규)
- 금지 파일(ADR-001, G0 보고서, 공유 모듈, 트랙 B 파일) 미작성·미수정, commit·push 없음
