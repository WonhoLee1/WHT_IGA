# 03 — 소형 절단셀 안정화: Ghost Penalty & AgFEM (LIT_NOTE)

## 문헌 상관
- P11 (Burman 2010, 원조): 공유면 법선 고차도함수 도약 페널티 →
  κ ≤ C·h^{-2p}, 강압성 증명.
- P12 (Badia 2022): 페이스기반 페널티의 γ→∞ 잠김 규명 → Locking-Free Ghost Penalty.
- P13 (Burman 2025): 이산확장 CutFEM–페널티 커널 관계 일반화, 임의고차 함수형 설계.
- P14 (Groß/Reusken): 부분공간분해 최적 전처리기, 격자·절단위치 독립 수렴.
- papers/: `../papers/P11_burman2010-ghost-penalty.md` 등 P11–P14.

## 오픈소스 상관
- **ngsxfem** (LGPL, NGSolve 애드온): 고차 레벨셋 매핑(Θ_h), Nitsche 약형식,
  Ghost Penalty 구현 레퍼런스. 면도약 행렬 조립 패턴을 본 PoC가 1D로 축소 재현.
- **GridapEmbedded.jl** (MIT): AgFEM/Extended B-spline 집합구속
  (불량 DOF → 건전 root DOF 외삽) 레퍼런스. 본 PoC의 T행렬 소거와 동형.
- master-plan §4.3: η_e<1e-2 셀 병합, C̃^e=C^e·T_e, κ 상한 보장.

## 핵심 알고리즘 (본 PoC)
1. 물리영역만 정확 적분한 1D Poisson 강성에서 η→0 시 κ 폭발 재현.
2. Ghost: 면도약 γ·h·⟦u'⟧² 부가 → κ 바운딩.
3. AgFEM: 외부 노드를 내부 노선형외삽으로 구속 소거 → 축소계 κ 바운딩 + 정확도 유지.

## 본 PoC 범위/한계
- 1D 선형요소. 2D/3D 면·능선 집합그래프(Rust aggregation_graph)는 본개발.
- 확인 항목: η 스위프 κ 표, 안정화 후 κ 상한, AgFEM L2 오차.
