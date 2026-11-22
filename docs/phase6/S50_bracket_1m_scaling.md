# S50 1M DOF 게이트 실측: Generator_bracket 스케일링 (2026-11-20)

## 질문
Phase 6 게이트 "1M DOF 3D 브래킷 end-to-end < 30분(16코어), 메모리 < 64 GB"를 이 PC에서 실측할 수 있는가? 기존 판단(G6_report, perf_vs_coreform)은 소규모 측정의 외삽이었다(AMGCL 1M: 약 22.9분, 10.71 GB).

## 조건
- 형상 `wht_iga/resources/Generator_bracket.step`(단일 솔리드, 얇은 벽), 테셀레이션 `0.003·diag`(11,860 삼각형).
- PCA OBB 배경 격자 `(11,7,4)·k`, p=2, q=2, extension θ=0.2, `SmallDisplacementElement3D4N` 선형 탄성, OBB x 양 끝 강제변위(S45와 같은 설정, 시험 helper 재사용).
- 풀이 `AMGCL_CG_ILU0`(tol 1e-10, S45 `_make_amgcl_or_pardiso`). k ≤ 2는 `PardisoLU`와 비교.
- PC: 32 GB RAM, 8코어/16스레드(게이트 사양 64 GB/16코어와 다름). 측정 중 가용 메모리 약 17 GB.
- 스크립트: `scripts/bracket_1m_scaling.py <k> <solver> [min_avail_gb]` (cwd `WHT_IGA`).

## 결과 (실측)

| k | 격자 | 활성 셀 / 절단 / fallback | 활성 DOF | 트림+분류+extension [s] | Kratos 모델 [s] | Build+풀이 [s] | 합계 [s] | 최대 working set [GB] |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 11×7×4 | 141 / 141 / 0 | 1,398 | 1.4 | 0.2 | 0.4 | 5.1 | 0.49 |
| 2 | 22×14×8 | 614 / 600 / 3 | 4,626 | 2.1 | 0.5 | 1.1 | 6.7 | 0.66 |
| 4 | 44×28×16 | 3,414 / 2,641 / 18 | 19,614 | 6.9 | 2.3 | 2.6 | 15.5 | 1.07 |
| 8 | 88×56×32 | 21,621 / 10,900 / 33 | 100,194 | 25.9 | 15.6 | 15.4 | 60.1 | 3.89 |
| 10 | 110×70×40 | 39,908 / 17,067 / 31 | 174,417 | 53.1 | 35.8 | 31.3 | 124.0 | 6.65 |
| **12** | 132×84×48 | 66,516 / 24,716 / 38 | **278,112** | 80.0 | 70.5 | 70.0 | **223.9** | **10.56** |

- 정확성: 모든 k에서 해 유한, 반력 합 균형 `|ΣR|/Σ|R|` ≤ 2.3e-16. AMGCL vs Pardiso 상대 차이 3.1e-9(k=1), 6.3e-9(k=2). 체적 오차 0.032 %(B-Rep 대비).
- k=12 Kratos 모델 단계 내역: NURBS 볼륨 2.4 s, 요소 생성 17.1 s, extension tie 43.3 s(가장 큼).
- 얇은 벽 형상이라 활성 DOF는 k³가 아니라 약 k^2.5로 늘어난다(k 8→12).

## 1M DOF 외삽 (실측 3점 k=8,10,12에 log-log 맞춤, DOF ≥ 5e4)

| 단계 | 지수 | 1M DOF 추정 |
|---|---:|---:|
| 트림+분류+extension | 1.11 | 345 s |
| Kratos 모델(요소+tie) | 1.48 | 470 s |
| Build+AMGCL 풀이 | 1.48 | 447 s |
| **합계** | 1.29 | **약 1,170 s (19.5분, 8코어)** |
| **최대 메모리** | 0.98 | **약 37 GB** |

## 결론
- **1M DOF end-to-end는 이 PC에서 실행 불가**: 최대 메모리가 DOF에 거의 비례(지수 0.98)해 1M에서 약 37 GB로, 가용 약 17 GB를 넘는다. 이 PC에서 안전한 상한은 약 30만 DOF(k=12 실측 27.8만, 10.6 GB).
- 게이트 사양(64 GB, 16코어) 기준으로는 시간(약 19.5분 < 30분, 8코어 측정이라 보수적)과 메모리(약 37 GB < 64 GB) 모두 **통과로 추정**된다. 단, 실행 확인이 아니라 외삽이다(1M 실행은 64 GB 장비 필요).
- 기존 외삽(G6_report / perf_vs_coreform: AMGCL 1M 10.71 GB)은 메모리를 **약 3.5배 과소평가**했다. 실측 최대 working set은 풀이 단계(Kratos Build + AMGCL 설정 + ILU0)에서 생기며, 모델 생성 후 RSS(k=12 2.97 GB)의 약 3.5배다.
- 시간 상위 항목(1M 추정): Kratos 모델 생성(특히 extension tie, k=12에서 43 s, 지수 > 1) ≈ Build+풀이 > 트림. Phase 6 C++/Numba 최적화 우선순위는 tie 생성과 요소 생성 루프, 그리고 풀이 단계 메모리(블록 행렬 AMGCL, 행렬 사본 축소).

## 재현
```powershell
cd D:\PythonCodeStudy\CAD-CAE-COPILOT\WHT_IGA
D:\mamba\envs\whtiga\python.exe -u scripts\bracket_1m_scaling.py 12 AMGCL_CG_ILU0 6
```
(k=12 약 4분, 최대 10.6 GB. 실행 전 가용 메모리 확인. 마지막 인자는 각 풀이 전 최소 가용 메모리 GB, 미달 시 중단.)
