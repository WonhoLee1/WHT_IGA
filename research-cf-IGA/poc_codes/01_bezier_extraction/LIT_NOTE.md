# 01 — Bézier Extraction (LIT_NOTE)

## 문헌 상관
- P01 (Immersed IGA 심층보고): R(ξ)=C^e·B(ξ) 분해로 상용 솔버 조립 루틴 재활용.
- P18 (U-splines): 셀별 추출행렬 희소 저장, 요소 루프는 Bernstein 루틴 재활용.
- Borden et al. 2011 (원전): isogeometric → Bézier 요소 변환 정식화.
- papers/: `../papers/P01_immersed-iga-deep-report.md`, `../papers/P18_thomas2022-u-splines.md`

## 오픈소스 상관
- **splinepy** (MIT, https://github.com/tataratat/splinepy): B-스플라인/NURBS 볼륨
  배경격자 + 초기 C^e 산출 + 역매핑. 본 PoC는 splinepy 미설치 환경을 가정하고
  Cox–de Boor + 최소자승 추출로 동일 수학을 재현 (인터페이스 스케치):
  ```python
  import splinepy  # 버전별 API 확인 필요
  b = splinepy.BSpline(degrees=[2], knot_vectors=[[0,0,0,0.5,1,1,1]],
                       control_points=[[0.],[0.33],[0.66],[1.]])
  # 요소별 C^e = b.extract_bezier(...) 계열 API로 획득 (버전 확인 필요)
  ```
- **Kratos IgaApplication**: C^e 주입 시 표준 Bernstein 요소 강성 루틴 그대로 사용.
- master-plan §4.1: K^e=(C^e)^T·K_B·C^e, AgFEM 단계에서 C̃^e로 교체.

## 핵심 알고리즘
1. 전역 스플라인 기저 N_i를 요소별로 Bernstein 기저 B_j와 행렬 C^e로 분해.
2. 요소 강성은 Bernstein 공간에서 표준 적분 후 C^e로 congruence 변환.
3. Greville 점 보간으로 선형장 재현 (패치테스트).

## 본 PoC 범위/한계
- 1D p=2, 2요소 모델. C^e를 최소자승으로 복원해 splinepy 대체 검증.
- 확인 항목: 재구성 오차, 분할단위 보존, K congruence 일치, Greville 선형재현.
- 한계: 다변수·유리NURBS·THB 추출은 본개발 (05 폴더와 결합 예정).
