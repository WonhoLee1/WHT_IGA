# 02 — Cut-cell 고차 적분 (LIT_NOTE)

## 문헌 상관
- P07 (Saye, Algoim): 레벨셋 절단영역을 다가지 높이함수 그래프로 재구성 →
  재귀 차원축소 → 1D Gauss-Legendre/Tanh-Sinh 텐서곱, 기계정밀도.
- P05 (Loibl 2026 비교리뷰): Algoim·BoSSS·FCMLib·QUGaR·Gridap·mlhp·ngsxfem·Nutils의
  정확도/적분점수/효율 비교. 곡면 B-Rep=축소법, 다면체=옥트리·피팅이 실무 선택.
- P06 (Antolin 무적분): 발산정리 2회 적용으로 체적→면→선적분 축소 (다항식 커널 한정).
- papers/: `../papers/P07_saye-algoim.md`, `../papers/P05_loibl2026-quadrature-comparison.md`

## 오픈소스 상관
- **Algoim** (BSD/MIT, https://algoim.github.io): C++17 헤더온리.
  셀 위상분류 → 높이함수 방향선택 → 재귀 축소 → 점/가중치 테이블을 FEM 조립에 주입.
- **QUGaR** (MIT): Algoim 내장 비적합격자 전용 적분·재매개변수화 파이썬 인터페이스.
- **QuESo** (TUM/Kratos): STL/B-Rep ↔ 배경셀 옥트리·모멘트피팅 적분점 산출.
- master-plan §4.2: ∫_K f = Σ w_k f(x_k), 1D q점당 2q차 수렴, O(h^{p+1}) 회복 목표.

## 핵심 알고리즘 (본 PoC)
1. 모멘트피팅: 후보점(셀 내부 조밀격자) × 단항식 행렬 A, 기준모멘트 b → nnls로
   희소 가중치 w. 점 수는 적지만 조건수 관리 필요.
2. 높이함수 축소(2D Saye식): x-스트립별 GL점 → y구간 해석적분(단항식 원시함수) →
   GL 가중합. 경계가 매끄러우면 고정밀.

## 본 PoC 범위/한계
- 원형 레벨셋 단일셀, 단항식 적분 검증. 3D 재귀축소·Tanh-Sinh 특이적응은 본개발.
- 확인 항목: 면적·고차모멘트 오차, 피팅 가중치 개수, 스트립 수렴.
