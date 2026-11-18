# 기법 선택 가이드 (사용자용)

## 결정 트리
1. 기존 FE 벌크 재활용? → **iga_fe_hybrid** (`mortar_l2`).
2. 접촉·굽힘·응력집중? → **full IGA p=2** (THB 국부 세분화 우선).
3. 열·질량·벌크? → **p=1** (full IGA 또는 native FE).
4. 메모리 초과? → `recommend_discretization(..., budget_bytes=...)` 자동 하향.
5. exactness 필요 (patch 수준)? → **mortar_l2** (barycentric 금지).

## 프리셋
```python
from wht_iga.analyses import presets as P
problem, option, rec = P.contact_iga(E=2.0e5)   # full IGA p=2
problem, option, rec = P.bending_thb(E=2.0e11)  # THB 적응
problem, option, rec = P.bulk_linear(1.0e5)     # native FE
```

## 실측 근거 (S35, OMP=1)
- 굽힘 최고 효율: THB-p2 (3.78e-4, 300 DOF).
- 최저 시간: native FE (hex20 굽힘 0.06 s).
- 최고 비용: full-p3 (peak 4e8 B, solve 1.36 s).
- p=1은 굽힘 금지 (전단 강성, thb-p1 실측 47.6%).
- barycentric은 O(h) 정합 (게이트 2e-3), exact 필요시 mortar.
