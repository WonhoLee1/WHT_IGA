# Phase 0 트랙 C 실측 요약 (S19–S24, 2026-11-14)

- 실행 환경: Python 3.12.16 (`D:\mamba\envs\whtiga\python.exe`), Kratos 10.4.3 (`E:\kratos-build\install`)
- 전체 테스트: **40/40 통과** (`test_s19` 16, `test_s20` 6, `test_s21` 3, `test_s22` 4, `test_s23` 5, `test_s24` 6)

## 1. 결과 표

| Spike | 주제 | 판정 | 핵심 실측 근거 | 아키텍처 결정·영향 |
|---|---|---|---|---|
| **S19** | `ConstitutiveModelsApplication` 초탄성 (Mooney-Rivlin, Yeoh, Ogden) | **실패 (대안 `HyperElastic3DLaw` 통과)** | `MooneyRivlin`/`Yeoh` Python 미노출; 바인딩된 `LargeStrain3DLaw`는 `QuadraturePointGeometry::Volume()` 미오버라이드로 `RuntimeError`. 대안 `CL.HyperElastic3DLaw` + `UpdatedLagrangianElement3D4N`은 `6.66e-16`~`3.33e-15` 일치 (`TotalLagrangian`은 `DETJ0` 결함 관측) | Neo-Hookean은 `UpdatedLagrangian` + `CL.HyperElastic3DLaw` 사용. Mooney-Rivlin/Ogden 등 다항 초탄성은 C++(`WhtIgaApplication` 또는 `QuadraturePointGeometry::Volume()`) 승격 필요 |
| **S20** | `Accessor` 온도 의존 물성 $E(T), \alpha(T)$ | **부분 성공 (`E(T)` 통과, `α(T)` 대안 통과)** | `ReadMaterialsUtility` JSON으로 `TableAccessor` 부착 시 `CL.ThermalElasticIsotropic3D`에서 $E(T)$ 상대오차 `6.013e-09` (`< 1e-8`). 단 `SM.LinearElastic3DLaw`와 `α(T)`는 `TableAccessor` 무시 → 셀별 `Properties` 복제로 $O(h^2)$ 수렴 (`1.049e-05` @ `n_x=32`) | $E(T)$ 단독은 `TableAccessor`, $\alpha(T)$ 동시 의존 시 셀별 `Properties` 복제 채택 |
| **S21** | `LinearMasterSlaveConstraint` ($c=0$) + 비선형 Newton | **통과 (`1e-15` 일치)** | `TotalLagrangian` 대변형 외팔보/블록(고정·비영 Dirichlet master 포함)에서 참조 $T^T K_t T$ Newton 대비 `1.44e-15`~`2.12e-16` 일치. 단 **`MoveMeshFlag=True`** 및 **`AndCriteria(ResidualCriteria, DisplacementCriteria)`** 필수 | 비선형 전략 표준에 `MoveMeshFlag=True` + `AndCriteria` 강제 |
| **S22** | `ContactSM` `mesh_tying_process` / `penalty_contact_process` + 표면층 | **고수준 프로세스 실패 (표면층 스프링 tie 통과)** | `InterfacePreprocessCondition`이 `QuadraturePointGeometry::GenerateFaces()` 호출로 `RuntimeError` 실패(별도 표면 파트도 부모 요소 부재로 실패). 대안 `Triangle3D3` 표면층 + MPC + `SpringDamperElement3D` tie는 단일 블록 대비 `5.220e-11` 일치 | Phase 5 접촉은 고수준 `ContactSM` 프로세스 대신 표면층 + 직접 `PairedCondition` 생성 또는 C++/Python 활성집합 경로 채택 |
| **S23** | `MappingApplication` `nearest_element` & `AssignNeighboursUtility` | **부분 성공 (`nearest_element` 운동학적 패치 통과)** | `nearest_element` mapper: 단위분할 `0.0`, 좌표/일정/운동학적 선형 패치 `5.55e-17`~`4.44e-16` (`< 1e-12`). 자유 계면 응력 패치는 비변분 오차 `0.086`. `AssignMasterSlaveConstraintsToNeighboursUtility`는 `7.94e-02`로 미달 | 비정합 계면 가중치 생성은 `AssignNeighboursUtility` 대신 `nearest_element` 또는 mortar $L^2$ projection(Python Numba) 사용 |
| **S24** | Explicit 동역학 + HRZ lumped $M$ + $\tilde{\mathbf{C}}^e$ $\Delta t_{\text{crit}}$ 복원 | **Kratos MPC 비호환 확인 / Python Explicit + $T$ 축약 통과** | `MechanicalExplicitStrategy`는 `PreComputeExplicitConstraintMassAndInertia`가 비변분 $T^T m_{\text{slave}}$를 더해 음수 확장 가중치 시 질량 `2.78e-17` 소멸. Python Explicit + HRZ lumped $M$ + $T^T f^{\text{int}}$는 1D 파동 전파 $L^2$ 오차 `1.16e-02`, $\eta=10^{-6}$ 절단 셀 $\Delta t_{\text{crit}}$ **714,131배 복원** | 침지 IGA Explicit는 Python Explicit 루프(또는 $\tilde{\mathbf{C}}^e$ 직접 반영 C++ 요소) + HRZ/Optimal Scaled lumped $M$을 주 경로로 확정 |
