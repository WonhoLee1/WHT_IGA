# TRACK B — boundary / loads / physics spikes (S01, S07–S12, S14, S17, S18)

- 브랜치: WHTIGA. 실행: `C:\Users\GOODMAN\AppData\Local\micromamba\micromamba.exe run -n whtiga python -m pytest WHT_IGA/tests/kratos_spikes/test_s01_*.py test_s07_*.py ...` (자기 파일만 지정 실행).
- 일시: 2026-11-14. 총 34 passed (S01:13, S07:3, S08:3, S09:2, S10:2, S11:5, S12:1, S14:2, S17:2, S18:1).

## 결과표

| spike | 질문 | 결과 | 근거 숫자 |
|---|---|---|---|
| S01 | env + HasMKL | 예 | 13/13, HasMKL=True, py3.12.16 |
| S07 | 표면 QP 기하 detJ·법선 | 예 | 면적합 1.0 (1e-13), 법선 권선 추종. 단 적분점당 기하 1개 (v2 가정과 다름) |
| S08 | PRESSURE/SURFACE_LOAD/FOLLOWER | 예 | 반력합 -2.0 / +3.0 정확; follower 상대차 1.56e-6; +p = 외향법선 |
| S09 | penalty BC | 예 | 오차 3.4e-12 (=1/β). 단 표면규칙 ≥2p차 필수 (2차규칙 고착오차 4.5e-5 발견·해결) |
| S10 | nodal/spring 요소 | 예 | 2.0 / 4/3 정확 (<1e-9) |
| S11 | 열 (Laplacian/조건/BDF/film/radiation) | 부분 | 전도 50.0 정확; film n=1 오차 0.19%→n=2 0.012% + 에너지보존 1e-6; radiation 2.835 일치; Laplacian조건 불가; BDF 질량행렬에 rho·c_p 없음 |
| S12 | 열팽창 | 예 | 모서리 1e-3 정확, 반력 2e-16 |
| S14 | mortar 접촉 | 부분 | 양 조건 등록·생성; 스프링-tie 힘전달 정확(-2.0, 점프 5e-11); 완전 마찰풀이 보류 |
| S17 | 출력 프로세스 | 부분 | quadrature JSON 좌표 일치(6자리); VTK는 Brep 입력 요구 |
| S18 | Nitsche 3D | 아니오 | 멤브레인 전용; 3D 호출 segfault |

## 경로 결정 (R1/R2/R3)에 미치는 영향

- **R1 유지** (주 경로 = 셀 QP 기하 + SM 요소 + 구속식 T): S07–S10, S12가 경계·하중·물리 부착을 전부 확인. S09의 규칙차수 제약(≥2p)과 S07의 점당 기하 1개는 빌더 구현 규칙으로 반영하면 됨.
- **R2/R3 전환 사유 없음**: S05(Track A) 실패 시에만 R2, S04(Track A) 실패 시에만 R3 — 모두 Track A 관할. Track B에서 R1을 깨는 실패 없음.
- **v2 수정 제안 3건** (통합 담당 판단용, v2 본문은 수정하지 않음):
  1. §5.4 "표면 적분점 기하" → 적분점당 기하 1개씩 생성됨 (조건도 점당 부착).
  2. §5.4 penalty/§7.4 표면 규칙: penalty 피적분식 차수(2p)에 맞는 Dunavant 규칙 강제 (부족 시 β-무관 고착오차).
  3. §5.6 열: `SupportLaplacianCondition` 삭제 (사용 불가) → 열 Dirichlet는 강한 부과; 과도 열은 참조 조립기 담당 (LaplacianElement 질량행렬에 rho·c_p 없음).

## 전 트랙 공통 중요 발견 (치명 함정 4건)

1. **Python 객체 수명 (S08)**: `NurbsVolumeGeometry` 파이썬 참조를 놓으면 표면조건 조립 RHS가 조용히 0이 됨 (0-잔차 오수렴, 반력 NaN, 예외 없음). 빌더는 volume + wrapper를 해석 끝까지 보관해야 함. Track A의 표면 래퍼 사용처도 동일 주의.
2. **미등록 변수 segfault (S11)**: `ConvectionDiffusionSettings`에 gradient/surface-source 변수 미설정 시 `ThermalFace`·`LaplacianElement` 후처리에서 access violation. `thermal.py`가 기본 포함.
3. **미사용 회전 자유도 (S10/S14)**: 스프링 요소 사용 시 체적 노드의 ROTATION은 Fix 필수 (singular).
4. **조건별 기하 요구 상이**: Load/Penalty(역학)는 점당 QP에서 동작, Laplacian/Nitsche/ThermalFace는 QP에서 segfault/예외 → carrier 일반 삼각형 + MPC(BC-C) 사용.

## 공유 모듈 변경 제안 (바꾸지 않고 기록만)

1. `wht_iga/kratos_backend/geometry.py` `make_nurbs_volume`: 모델파트에 volume 객체 보관 (예: `model_part._wht_iga_volumes` 리스트에 추가) — 함정 1의 근본 방지. 또는 빌더 규약으로 문서화.
2. `wht_iga/grid/*`, `wht_iga/reference/*`, `tests/conftest.py`, `pyproject.toml`, `uv.lock`: 변경 필요 없음.
3. `surface.py` 신규 규칙 `triangle_quadrature(4)`는 Phase 1 빌더가 그대로 사용 가능.

## 산출 파일

- 시험: `WHT_IGA/tests/kratos_spikes/test_s01_*.py`, `test_s07_*.py`, `test_s08_*.py`, `test_s09_*.py`, `test_s10_*.py`, `test_s11_*.py`, `test_s12_*.py`, `test_s14_*.py`, `test_s17_*.py`, `test_s18_*.py`
- 신규 모듈: `WHT_IGA/wht_iga/kratos_backend/surface.py` (표면 QP + Dunavant 1/2/4차 + 수명 경고), `thermal.py` (열 변수/DOF/속성/settings 기본값)
- 문서: `WHT_IGA/docs/spikes/S01/S07/S08/S09/S10/S11/S12/S14/S17/S18.md`, 본서
- 금지 준수: 공유 모듈·Track A 파일·ADR·G0 미작성, commit/push 없음
