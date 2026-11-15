# 05 — 계층 B-스플라인 세분화 (LIT_NOTE)

## 문헌 상관
- P18 (U-splines): h·p·연속성 국소 가변, T-접합/혼합셀, 국소 영공간 해법으로
  분할단위·양성·국소선형독립 증명.
- P04 (Divi 스캔논문): THB-스플라인 적응세분화 + 오일러지표 위상검사.
- P02 (Schillinger 2012): 계층 B-스플라인 + FCM Design-through-Analysis.
- papers/: `../papers/P18_thomas2022-u-splines.md`, `../papers/P04_divi2021-scan-based.md`

## 오픈소스 상관
- **G+Smo** (MPL-2.0): gsHSplines 계층 스플라인 C++ 레퍼런스.
- **thbsplines**: 절단 THB 옥트리 관리, 활성기저 식별.
- **splinepy**: knot insertion 기반 세분화·베지에 추출.
  본 PoC는 Boehm 삽입을 numpy로 직접 구현해 기하 불변성 증명.

## 핵심 알고리즘 (본 PoC)
1. Boehm knot insertion → 제어점 갱신, 곡선 형상 불변.
2. 계층기저 선택: 세분영역 밖 coarse + 영역 안 fine → 분할단위 유지,
   전역 fine 대비 기저수 절감.

## 본 PoC 범위/한계
- 1D p=2. THB 절단(truncation) 연산자·2D T-접합은 본개발.
