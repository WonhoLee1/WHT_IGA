# master-plan v2 Addendum 1 — Coreform IGA User-Guide 전수 대조 결과와 초과 달성 계획

- 문서 버전: v2-add_1 (2026-11-14)
- 모문서: `master-plan_v2.md` (수정하지 않음, 본 addendum이 우선 적용되는 항목은 §4에 명시)
- 근거 1: `user-guide_combined.md` / `user-guide_combined.html` (82개 md 원천 결합, `build_combined_guide.py` 생성, `generated-downloads/`·`site_libs/`·`.omc/` 제외)
- 근거 2: `master-plan_v2.md` 전문 (§0–§14, 부록 A/B)
- 근거 3: Kratos 소스 실측 (`D:\PythonCodeStudy\kratos-src\applications\ConstitutiveLawsApplication\custom_constitutive\`)
- 상태: **계획 보완안. 본 addendum의 G-1–G-8, S-항목은 master-plan v2 Phase 계획에 편입 검토 대상이다.**

---

## 0. 대조 방법

1. 결합 문서 82개 섹션을 기능군으로 분류: 설치/구성(6) · 이론(2) · 사용자 가이드 핵심(7) · Python API(11) · 릴리스 노트(5) · 응답 스펙트럼(4) · 스크립팅 가이드(5) · 검증 매뉴얼(12 문제 + 인덱스) · 예제(5) · 키워드 지원 가이드(49 케이스) · CIFA 역량(JS 렌더, 내용 없음).
2. 각 Coreform 기능을 master-plan v2 §1.3 매트릭스·§5·§7·§10과 1:1 대조. 판정은 `동등` / `부분` / `미포함(의도적)` / `미포함(갭)` 4단계.
3. `미포함(갭)`은 구현 가능성(Yes/No)과 함께 §2에, 설계상 이미 우위인 항목은 §3에 분리 기록. 의도적 제외는 §5에서 재확인한다.

---

## 1. 커버리지 판정 매트릭스

### 1.1 격자·세분화·소형 셀 (동등)

| Coreform 기능 (출처) | master-plan v2 대응 | 판정 |
|---|---|---|
| 직교 배경 격자, p=1–3, C0–C^{p-1}, element size/interval count/축별 지정, edge/node_centered, padding, Frame(원점·기저·범위, 회전 포함) (`python-api/coreform/mesh`, `release-notes-2026-3 §2`) | §1.3 격자 5행, §7.2 (knot 물리좌표, Greville, Bézier 추출, Frame 회전 detJ=1) | 동등 |
| tessellation_size, 삼각형 분할 절단 (`python-api/coreform/mesh`) | §7.1 (OCP BRepMesh, face_id 보존), §7.4 | 동등 |
| 계층 세분화 + 2:1 균형 + CadEntities/Sphere/Box region + offset p·h/2^L (2026-11: Offset 공란 시 degree·level 기반 기본값) | §7.2 THB + §5.2.3 구속식, E02/V02-local 게이트 | 동등 |
| small_cell_volume_ratio 기본 0.2·범위 [0, 0.5) (릴리스 노트·mesh 문서) | §7.5 (기본 0.2, 기저 지지 기준 extended B-spline MPC) | 동등 |
| U-spline | §1.4 범위 밖(특허) — 의도적 제외 | 미포함(의도적) |
| Body-fit mesh (정렬 hex) | §1.3 Phase 3 Alpha (정렬 블록 한정) | 동등(범위 일치) |
| Selective integration (nearly-incompressible, 2026-3 §2.2) | §7.8 reduced/SRI + B-bar/F-bar 옵션 | 동등 |

### 1.2 경계·하중·결합 (동등, 1건 부분)

| Coreform 기능 (출처) | 대응 | 판정 |
|---|---|---|
| penalty 필수 BC (`bc_penalty_value`, 검증 1e10) + 변위/핀/대칭/사용자 함수 변위 | §5.4 BC-A/B/C, §7.6, Phase 1 게이트 | 동등 |
| 압력(P)·총하중 압력·트랙션(TRVEC)·중력/체적력(GRAV/BX·BY·BZ)·RP 모멘트+kinematic coupling (`abaqus_keyword_support_guide` *DLOAD 10건·*DSLOAD 7건 Supported 행) | §5.4 표 (PRESSURE/FOLLOWER/SURFACE_LOAD/총하중 p=F/Σw/VOLUME_ACCELERATION/RP·RBE3) | 동등 |
| amplitude(시간 함수)·다단계 step·하중 변경/리셋 (2026-11 §3, 키워드 가이드) | §5.6 step 루프, §7.9 (Python 시간함수 조건값 갱신, Phase 3) | 동등 |
| multi-cell 불리언 병합(공유 스플라인 장)·surface tie·kinematic/distributing coupling·IGA+FE 혼합 (`release-notes-2026-7 §1.3`, 키워드 가이드) | §7.7, §5.4 BC-C, §1.3 결합 4행 | 동등 |
| follower 트랙션·resultant 트랙션 | §5.4 `PRESSURE_FOLLOWER_LOAD`, §7.7 | 동등 |
| film·radiation (경계 삼각망 위 적분, amplitude 없이 heat step) (키워드 가이드 *SFILM·*SRADIATE Supported 행) | §7.11 (`ThermalFace3D3N`, Newton) | 동등 |
| Surface-to-surface 접촉, penalty Coulomb μ, ALL EXTERIOR inclusions (키워드 가이드 *CONTACT INCLUSIONS Supported 2건; 자기접촉 제외) | §1.3 접촉행 Phase 5 Beta + S14 spike | 동등 |

### 1.3 해석 종류 (동등)

| Coreform 기능 (출처) | 대응 | 판정 |
|---|---|---|
| 선형 정적·다단계, NLGEOM (키워드 가이드 NLGEOM Recognized; 2026-가을 릴리스 목표 nonlinear statics) | §5.6 (선형/Newton/TL·UL), Phase 1/5 | 동등 |
| 고유진동수 Lanczos·AMS, 질량정규화 고유벡터 (response-spectrum 문서) | §5.6 `spectra_sym_g_eigs_shift`/FEAST + 질량정규화 확인, V09 | 동등 |
| 암시적 동역학 HHT (V10 vibrating-beam; Coreform 절차 단축어) | Kratos Bossak/Newmark/BDF로 대체 (사용자 결정 D5, §7.9, R11). 물리량 기준 비교로 전환 | 동등(방법 대체, §5.6) |
| 응답 스펙트럼 ABS/SRSS/NRL/TENP/CQC/GRP/DSC + 방향 ALGEBRAIC/SRSS/R40/R30 + log-log 보간 (`response-spectrum/` 4문서) | §7.10 전부 (조합식·보간·S_D=S_V/ω=S_A/ω² 포함) | 동등 |
| 랜덤 응답 PSD BASE/FORCE/DB + Segalman RMS von Mises (`random-response.md`) | §7.10 (H_i(ω), 사다리꼴 적분, 입력단위 환산) | 동등 |
| 정상/과도 열전달 (전도·film·radiation·열용량, BDF) + 열변형 (NAFEMS T1) (verification-manual heat_transfer 3건 + thermal_strain 1건) | §7.11, Phase 4, V11/V12 | 동등 |
| Explicit 동역학·Abaqus/Aqua | 양쪽 모두 미지원 | 미포함(의도적, §1.4) |
| 좌굴 고유치 | 양쪽 모두 없음 (검증 매뉴얼에 buckling 스터디 없음) | 해당 없음(동등하게 없음) |
| 에너지 이력 (kinetic/stored/plastic, 2026-11 §4 Alpha) | §4.3 `energy` 표 (kinetic/strain/plastic per frame) | 동등 |

### 1.4 재료 (1건 부분, §2 G-4로 이관)

| Coreform (B-bar/F-bar Alpha: "selected Abaqus elastic, hyperelastic, plastic") | 대응 | 판정 |
|---|---|---|
| 선형 탄성 + 온도 의존 테이블 | §1.3 + `TableAccessor` (S15) | 동등 |
| Neo-Hookean + J2 다선형 소성(온도 의존) (E03/E04/E05 게이트) | §1.3 Phase 5 (HyperElastic NeoHookean, SmallStrainJ2+테이블) | 동등 |
| 그 밖의 초탄성/소성 바리에이션 (Mooney-Rivlin·Yeoh·Ogden-3D·Marlow·Drucker 등 Abaqus 재료군) | master-plan 미기재. Kratos 실측: Hencky·Kirchhoff·Simo-Taylor Neo-Hookean·Ogden-**1D**·J2 유한변형 팩토리·damage·viscoplasticity 보유, Mooney-Rivlin/Yeoh/Marlow-3D 없음 | 부분 → **G-4** |

### 1.5 출력·후처리·API·검증 방식 (동등 또는 우위, §3으로 이관)

| Coreform 기능 | 대응 | 판정 |
|---|---|---|
| U·S·PEEQ·E·NE·LE·NT11 + PRESS_IGA + ALL/PRESELECT 확장 (키워드 가이드 *OUTPUT 27 Supported 행) | §1.3 출력행 + `CalculateOnIntegrationPoints` (§5.6) | 동등 |
| point/line/multi-point probe + probe JSON | §5.4/§7.12 가중치 0 적분점 probe + `<job>_probe_data.json` | 동등 |
| conforming 시각화 proxy 격자 (untrimmed hex/trimmed tet, 2026-10 §2.1) + ParaView 플러그인 직접 읽기 | §7.12 (IN-hex/CUT-tet 분할 + VTKHDF/VTU) | 동등 |
| results/quadrature SQLite (`<job>.sql`, `_iga_results.sql`, `_iga_quadrature.sql`) | §4.2–4.4 동일 스키마 (`wiga.sqlite`, results/quad SQLite) | 동등 |
| portable spec + validate + journaling/replay + CAE Job Monitor/백그라운드 메시 (2026-7 §1, 2026-11 §1) | `Job.to_json/from_json/validate→ValidationError`, runner subprocess (§3.3) | 동등(기능), **G-2**(진행 콜백) |
| 검증 매뉴얼 12건 + 예제 5건 + CSV 기준값 + 문서 자동 생성 | §10 매핑표(V01–V12 + E01–E05) + §10.6 자동 리포트 | 동등 이상 → **S-7** |
| 일관 단위계 문서 (`theory/consistent-unit-systems`) | meta `units`만 존재, 단위 정합 검사 규칙 없음 | 부분 → **G-6** |

---

## 2. 갭 — 구현 가능하나 master-plan v2에 없는 항목 (G-1–G-8)

### G-1. Abaqus `.inp` 읽기 전용 서브셋 (선택, Phase 6+)
- **내용**: Coreform의 keyword-aware translator와 반대 방향이 아니라, 벤치마크 이관을 위한 **읽기 전용** 임포터. 적용 범위: *DLOAD(GRAV/P/TRVEC)·*DSLOAD(P/TRVEC)·*SFILM/*SRADIATE(amplitude 없음)·*SOLID SECTION(MATERIAL 지정)·*STEP(NLGEOM/UNSYMM 플래그 인식)·*CONTACT INCLUSIONS(ALL EXTERIOR만). 지원 외는 `validate()` issue + 거부 (Coreform의 "reject with diagnostic"와 동일 정책).
- **근거**: `abaqus_keyword_support_guide.md` Supported 27행이 곧 스펙이다. master-plan §1.4는 `.inp` 번역기를 범위 밖으로 두었으나, **읽기**는 UEL/odb 의존이 없어 사용자 결정 1·4와 충돌하지 않는다.
- **수용 기준**: V03 6하중·V04·E03 입력 재현 시험 (hand-written JSON 대비 diff空).

### G-2. Job 진행 콜백/취소 (Phase 3)
- **내용**: `Job.submit()`에 `on_progress(stage, frac, eta)` + `cancel()` (Coreform Job Monitor 대응). 전처리 셀 블록 단위, Kratos step은 increment 콜백.
- **수용 기준**: E03 규모에서 진행률 단조 증가 + cancel 시 부분 산출물 정리 단위시험.

### G-3. Coreform 스펙트럼 case-file 호환 리더 (Phase 4)
- **내용**: `response-spectrum/case-files` 형식의 스펙트럼 테이블 리더 → §7.10 numpy 조합기에 그대로 입력. Coreform case 파일 재사용 + 동일한 조합 규칙 이름(ABS/SRSS/NRL/TENP/CQC/GRP/DSC, R40/R30) 사용.
- **수용 기준**: V09 랜덤/스펙트럼 값을 Coreform case 파일 입력으로도 재현.

### G-4. 재료 바리에이션 래더 (Phase 5, C++ 승격과 연계)
- **내용**: Kratos 실측 보유분을 순서대로 활성화. ① Kirchhoff/Hencky/Simo-Taylor Neo-Hookean (헤더 보유 확인), ② 유한변형 J2 isotropic/kinematic 팩토리, ③ small-strain damage·viscoplasticity. Mooney-Rivlin/Yeoh/Ogden-3D는 Kratos에 없으므로 **신규 구성법칙 = C++ 승격 사유 추가** (§5.7에 "G-4 재료" 항 추가).
- **수용 기준**: 각 재료별 patch 시험 1건 + E04급 Cook 막 1% 이내 (Neo-Hookean 이외 1종 이상).

### G-5. THIN_FEATURE 자동 세분 루프 (Phase 3)
- **내용**: 2026-11 "highly local refinement near geometric features is more robust" 대응. `validate()`의 THIN_FEATURE issue → 영향 영역 BoxRegion 자동 생성 → L+1 세분 → 재해석 (최대 2회, DOF 상한内). 수동 대응(현재 §7.5 권고)을 자동 루프로 승격.
- **수용 기준**: 얇은 형상 케이스에서 수동 대비 동일 정확도, 개입 0회.

### G-6. 일관 단위계 가드 (Phase 1)
- **내용**: `theory/consistent-unit-systems` 대응. `ModelSpec`에 `unit_system` (MMKS/SI 등) + 재료·하중·기하 단위 정합 검사 (`validate()` issue). 계산은 무차원 일관 단위로만 수행.
- **수용 기준**: 단위 불일치 모델 3종이 `validate()`에서 전부 차단.

### G-7. Abaqus 규약 출력 이름 매핑표 (Phase 3)
- **내용**: 변형률 Abaqus 규약 추종(2026-11 §4) 대응. 후처리 변수명 매핑표 (E↔LE, S분량, PRESS_IGA) + VTKHDF 필드 별칭. 벤치마크 숫자 비교 시 이름 차이로 인한 오판 방지.
- **수용 기준**: V02/V06 후처리 값이 이름 매핑 후 CSV와 직접 비교 통과.

### G-8. 성능 비교 하네스 (Phase 6)
- **내용**: Coreform 공개 수치와 동일 조건 비교표 자동 생성 (glass-to-metal seal: mesh 318.8 s + solve 691.0 s @16CPU; 2026-11 대형 어셈블리 메모리/시간 절반). `docs/verification/<ID>.md`에 `perf_vs_coreform` 섹션 추가.
- **수용 기준**: E03 1회 측정 + Phase 6 게이트(1M DOF 30분/64GB)와 함께 리포트.

---

## 3. 설계상 이미 우위인 항목 (S-1–S-10, 벤치마크 초과 달성 포인트)

- **S-1. 솔버 라이선스 독립**: Coreform은 Abaqus/Standard+HF5+CAE 플러그인 종속. WHT_IGA는 Kratos 소스빌드 + Python만으로 동일 물리 재현 (§0.4-1). 벤치마크 초과 조건 1순위.
- **S-2. Oracle 교차검증**: Numba 참조 조립기와 1e-10 일치 (§5.5, §10.4). Coreform에는 공개된 독립 교차검증 경로가 없다.
- **S-3. 조건수 게이트**: V-COND (η∈[1e-8,1] 평탄 + h 기울기 −2±0.3). Coreform 공개 기준에 없는 정량 안정성 게이트.
- **S-4. Dirty-CAD 내성**: 3축 투표 + winding fallback + octree fallback + `QUAD_FALLBACK` issue (§7.3–7.4). Coreform 2026-11도 "복잡 CAD 강건화"를 과제로 병기.
- **S-5. 개방형 산출물**: VTKHDF/VTU + SQLite + probe JSON (§4). Coreform은 ODB + 전용 ParaView 플러그인 필요.
- **S-6. FCM α 기준선 + SBM 비교 트랙** (§7.5, §10.5). 단일 주 경로에 대한 대조군을 상시 유지.
- **S-7. 문서-증거 동기**: §10.6 케이스별 자동 리포트 + CSV 기준값. Coreform "verification becomes documentation"와 동등 이상 (oracle diff 포함).
- **S-8. 시간적분 선택권**: Bossak 기본 + Newmark/BDF2 (§7.9). Coreform=HHT 단일. V10에서 물리량 기준 동등성으로 입증 계획.
- **S-9. Nitsche 선택권 + β sweep 확정 절차** (§7.6, ADR-002). Coreform 기본 penalty 일변도 대비 검증 폭이 넓다.
- **S-10. C++ 승격의 측정 기반 규칙** (§2.1-3, §5.7): 50% 조립 오버헤드·30% 핫스팟 같은 수치 게이트. Coreform 비공개 NativeTarget과 달리 재현 가능.

---

## 4. master-plan v2에 대한 수정 제안 (본 addendum 편입분)

1. §5.7 C++ 승격 기준에 "G-4 재료(신규 구성법칙 필요 시)" 추가.
2. §10.6 리포트 템플릿에 `perf_vs_coreform`(G-8)·`unit_check`(G-6)·`mapping_alias`(G-7) 섹션 추가.
3. Phase 3 산출물에 G-2(진행 콜백)·G-5(자동 세분 루프)·G-7(매핑표) 추가. Phase 1에 G-6, Phase 4에 G-3, Phase 5에 G-4, Phase 6+ 선택에 G-1 추가.
4. §14 열린 결정에 D11(`.inp` 읽기 서브셋 채택 여부, 기본안: Phase 6+ 선택) 추가.

---

## 5. 의도적 제외 재확인 (벤치마크 대상이 해도 안 하는 것 + 우리가 안 하는 것)

- U-spline (특허, §1.4) · Abaqus `.inp` 쓰기/UEL·`.odb` 출력 (사용자 결정 1) · Explicit 동역학 (양쪽 미지원) · 자기접촉 (Coreform 제외) · 셸/보 IGA (Coreform 3D 솔리드 전용, §1.4) · missing-mass/Gupta/Lindley-Yow (양쪽 범위 밖, §7.10) · Algoim 정확형상 적분 (Phase 6 선택, §7.4) · HHT 직접 구현 (D5 확정, Bossak 상위호환).
- 키워드 가이드의 "Not supported for IGA" 16행 중 IGA 영역 관련 핵심 (U994–U999 세트, A/RF 출력, RE/PE 변수)은 우리 출력 스키마(§4.3: U·S·PEEQ·E·NE·LE·NT11·PRESS_IGA·에너지)에도 동일하게 두지 않는다. 즉 Coreform이 IGA에 안 주는 변수는 우리도 주지 않는다.
