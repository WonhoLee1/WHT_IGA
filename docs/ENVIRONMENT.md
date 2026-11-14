# WHT_IGA 개발 환경 (단일 원본)

다른 문서와 이 문서가 다르면 이 문서가 맞다. 바꿀 때는 이 문서부터 고친다.

## 1. 요약

| 항목 | 값 |
|---|---|
| 환경 관리자 | **micromamba** 2.10 (`%LOCALAPPDATA%\micromamba\micromamba.exe`), root `D:\mamba` |
| env 이름 | `whtiga` (`D:\mamba\envs\whtiga`) |
| Python | 3.12.x (Kratos 빌드 ABI와 일치해야 함) |
| Python 패키지 설치 | **uv pip만** |
| Kratos | 소스빌드 `E:\kratos-build\install` (10.4.3 + IgaApplication). 상세: [`KRATOS_BUILD.md`](KRATOS_BUILD.md) |
| 금지 | `conda install/create`, uv 없이 `pip install`, `whtiga`에 pip Kratos wheel 설치 |

## 2. 실행

```powershell
micromamba activate whtiga          # 대화형 셸
micromamba run -n whtiga python ... # 단발 실행 (에이전트는 이 형태 권장)
D:\mamba\envs\whtiga\python.exe ... # 직접 경로도 가능
```

셸 프로파일(PowerShell 7/5, cmd AutoRun, Git Bash `~/.bash_profile`)은 2026-11-14에 conda 훅에서 micromamba 훅으로 바뀌었다. 백업: `C:\Users\GOODMAN\env-backup-2026-11-14\`.

기존 conda env(vdmc, aieng311, calculix-env 등)도 `micromamba activate <이름>`으로 쓸 수 있다 (`~/.mambarc`의 `envs_dirs`). miniconda 설치본과 PATH 항목은 유지한다 — aieng 백엔드가 CalculiX를 `conda run -n calculix-env ccx`로 찾기 때문이다.

## 3. 처음부터 다시 만들기

```powershell
# 1) env
micromamba create -y -n whtiga -c conda-forge "python=3.12" mkl hdf5
# 2) 패키지 (저장소 루트에서)
uv pip install --python D:\mamba\envs\whtiga\python.exe -e "WHT_IGA[cad,viz,solve,test]"
# 3) 확인
micromamba run -n whtiga python WHT_IGA/scripts/check_env.py
```

- 버전 고정: `WHT_IGA/uv.lock` (`uv lock --python D:\mamba\envs\whtiga\python.exe`로 갱신).
- conda 쪽 spec: `WHT_IGA/env/whtiga.yml`.

## 4. Kratos 연결 방식

- `wht_iga/kratos_backend/_bootstrap.py`의 `load_kratos()`만 Kratos를 import한다.
  - `sys.path` 앞에 `E:\kratos-build\install`, `...\libs`를 넣는다.
  - DLL 디렉토리 `E:\kratos-build\install\libs`, `D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin`을 `os.add_dll_directory`로 등록한다.
  - import 후 `KratosMultiphysics.__file__`이 설치 경로 아래인지 검사하고, 아니면 `RuntimeError`.
- 경로 바꾸기: 환경변수 `WHT_IGA_KRATOS_PATH`, `WHT_IGA_KRATOS_DLL_PATH`(`;` 구분).
- `.pth` 파일은 쓰지 않는다. Kratos가 필요 없는 코드(전처리, 후처리, 참조 조립기)는 Kratos 경로 없이 동작해야 하기 때문이다.

## 5. 확인된 동작 (2026-11-14)

`check_env.py` 13/13 통과:
- Python 3.12.16, numpy 2.5.4, scipy 1.18.1, numba 0.68.1, OCP 8.0.1.1, h5py 3.16.0, meshio 5.3.5, pyvista 0.49.2, pypardiso 0.4.7
- 위 패키지를 먼저 불러온 뒤 Kratos 로드 성공 (DLL 충돌 없음)
- Kratos 앱: Iga, StructuralMechanics, LinearSolvers, ConstitutiveLaws, ConvectionDiffusion, ContactStructuralMechanics
- `LinearSolversApplication.HasMKL()` = True (PARDISO 사용 가능)

## 6. 자주 나는 문제

| 증상 | 원인 | 해결 |
|---|---|---|
| `wrong Kratos loaded: ...` | pip Kratos wheel이 env에 설치됨 | `uv pip uninstall --python ... KratosMultiphysics ...` |
| `Kratos build not found` | 설치 경로 변경 | `WHT_IGA_KRATOS_PATH` 설정 |
| Kratos import 시 `DLL load failed` | DLL 디렉토리 누락, Python 버전 불일치 | `WHT_IGA_KRATOS_DLL_PATH` 확인, Python 3.12인지 확인 |
| cmd에서 `micromamba` 인식 안 됨 | 부모 프로세스가 `CONDA_SHLVL`을 물려줌 (예: 예전 conda 셸에서 띄운 경우) | 새 창을 열거나 `set CONDA_SHLVL=` 후 다시 |
| `micromamba activate` 안 됨 | 셸 훅 미적용 | 새 창. 또는 `micromamba run -n whtiga ...` 사용 |
