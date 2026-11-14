# Kratos 소스빌드 재현 절차

WHT_IGA는 공식 pip wheel이 아니라 이 PC의 Kratos 소스빌드를 쓴다 (공식 10.4.x wheel에는 IgaApplication이 없다).
이 문서는 그 빌드를 다시 만들 수 있도록 확인된 설정을 기록한다. **이 문서를 작성하면서 재빌드는 하지 않았다.**

## 1. 현재 빌드 (2026-11-14 확인)

| 항목 | 값 |
|---|---|
| 소스 | `D:\PythonCodeStudy\kratos-src`, 태그 v10.4.3, 커밋 `96bad120` |
| 빌드 트리 | `E:\kratos-build\Release` (Ninja) |
| 설치 | `E:\kratos-build\install` (`KratosCore.dll` 2026-09-28, `KratosIgaApplication.pyd` 2026-09-27) |
| 빌드 형식 | `CMAKE_BUILD_TYPE=Release`, 플래그 `/O2 /Ob2 /DNDEBUG` |
| 컴파일러 | MSVC 14.41 (`cl`), VS 2022 Build Tools `vcvars64.bat` |
| 빌드 도구 | `D:\PythonCodeStudy\kratos-deps\buildenv\Scripts\ninja.exe` (CMake 4.4도 같은 폴더) |
| Python | `C:\Users\GOODMAN\miniconda3\envs\vdmc\python.exe` (3.12.10) → **런타임 Python도 3.12.x여야 한다** |
| Boost | `D:\PythonCodeStudy\kratos-deps\boost_1_86_0` |
| MKL | `USE_EIGEN_MKL=ON`, `MKL_RT_LIB=D:/PythonCodeStudy/kratos-deps/condadeps/Library/lib/mkl_rt.lib` (런타임 MKL 2026.0) |
| SuiteSparse | `USE_EIGEN_SUITESPARSE=ON` (`SUITESPARSE_LIBRARY` 캐시값은 NOTFOUND로 남아 있음) |
| MMG | `INCLUDE_MMG=ON`, `MMG_ROOT=D:/PythonCodeStudy/kratos-deps/mmg-5.8-novtk` |
| 기타 | `USE_TRIANGLE_NONFREE_TPL=ON`, `KRATOS_BUILD_TESTING=OFF`, `KRATOS_GENERATE_PYTHON_STUBS=OFF` |
| 런타임 DLL | `E:\kratos-build\install\libs`, `D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin` |
| 앱 | 47개 (아래 목록) |

빌드된 앱(설치본 `libs\Kratos*Application.pyd` 기준): CableNet, Chimera, CompressiblePotentialFlow, ConstitutiveLaws, ConstitutiveModels, ContactMechanics, ContactStructuralMechanics, ConvectionDiffusion, CoSimulation, Dam, DelaunayMeshing, DEM, DemStructuresCoupling, DropletDynamics, FemToDem, FluidDynamics, FluidDynamicsBiomedical, FluidDynamicsHydraulics, FreeSurface, FSI, GeoMechanics, HDF5, **Iga**, LinearSolvers, Mapping, Med, Meshing, MeshMoving, MPM, Optimization, PFEM2, Pfem, PfemFluidDynamics, PfemSolidMechanics, Poromechanics, RANS, Rom, ShallowWater, ShapeOptimization, SolidMechanics, Statistics, **StructuralMechanics**, SwimmingDEM, SystemIdentification, ThermalDEM, TopologyOptimization, WindEngineering (+ 각 이름 뒤 `Application`).

## 2. 소스의 로컬 수정

- `git status --short` 기준 295개 항목 (추적 파일 수정 88개: 448줄 추가, 328줄 삭제 + 미추적 파일 다수).
- 주로 `ConstitutiveModelsApplication`(166), `SolidMechanicsApplication`(62), 그 외 MSVC 호환 수정.
- **IgaApplication, `kratos/geometries`, `kratos/integration`은 수정 없음.**
- 상세: `kratos_build/local_changes.txt`(`git diff --stat`), `kratos_build/local_status.txt`(`git status --short`).
- 재빌드 전에 이 수정들을 보존해야 한다(브랜치 또는 patch 파일로 저장 권장).

## 3. 빌드 스크립트

- 원래 스크립트 `D:\PythonCodeStudy\kratos-deps\build_kratos_full.bat`는 **1바이트(내용 소실)**.
- 복원 초안: `kratos_build/build_kratos_full.bat.draft` (CMakeCache와 설치된 앱 목록에서 재구성, Kratos `scripts/ninja_configure.bat` 형식).
- 초안은 실행하지 않았다. `KRATOS_APPLICATIONS`는 CMake 캐시에 남지 않는 값이라 설치본의 `.pyd` 목록에서 추정했다.
- 원본 캐시: `kratos_build/CMakeCache.txt`.

## 4. 재빌드 절차 (필요할 때만: 예) C++ 커스텀 앱 추가)

1. 로컬 수정 보존: `git -C D:\PythonCodeStudy\kratos-src stash` 대신 `git diff > local.patch` + 미추적 파일 백업 (stash는 쓰지 말 것).
2. 초안을 검토해 `build_kratos_full.bat`으로 복사하고, 커스텀 앱 경로를 `add_app`에 추가.
3. 새 cmd 창(conda/micromamba 미활성)에서 실행. 증분 빌드면 기존 `E:\kratos-build\Release`를 그대로 쓴다.
4. 병렬도는 **`-j 6` 이하** (32 GB RAM에서 그 이상은 `cl.exe` 충돌 기록).
5. 확인: `micromamba run -n whtiga python WHT_IGA/scripts/check_env.py`.

## 5. 주의

- `whtiga` env에 pip Kratos wheel을 설치하지 않는다 (DLL·버전 충돌).
- `D:\PythonCodeStudy\kratos-playground`의 venv들(10.0/10.3/10.4.3/10.4.4 wheel)과 섞지 않는다. 그 안의 `Kratos` 클론은 이 빌드의 소스가 아니다.
