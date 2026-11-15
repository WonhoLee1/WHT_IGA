# 04 — Nitsche 약형식 Dirichlet (LIT_NOTE)

## 문헌 상관
- master-plan §4.4: 대칭 Nitsche 정식화 + γ/h 페널티항.
  γ 미달→강압성 상실(발산), 과다→조건수 악화. 고차 연속성(C^{p-1}) 덕에
  경계 트랙션 진동이 표준 FEM 대비 억제.
- P11/P20 계열: Nitsche는 비적합법의 표준 약형식 경계 부과 수단.
- P06 무적분법의 경계 선적분과 결합 시 Nitsche 표면항도 동일 백엔드 재활용.

## 오픈소스 상관
- **ngsxfem**: Nitsche 약형식 + Ghost Penalty 결합 구현 레퍼런스.
  본 PoC는 1D 대칭 Nitsche의 γ 민감도를 그대로 축소 재현.

## 핵심 알고리즘 (본 PoC)
a_N(u,v) = a(u,v) −u'(L)v(L) −u(L)v'(L) +(γ/h)u(L)v(L), 우단 절단점 x=L 부과.
- γ 스위프 → 최소고유값 부호 + L2 오차 + κ 동시 추적.
- 실무 처방(master-plan §5.6): γ_e=C·trace(K^e) 요소별 자동산정.

## 본 PoC 범위/한계
- 1D 선형, 고정 η=1e-3. 3D 곡면 Nitsche 표면적분(Algoim 경계점)은 본개발.
