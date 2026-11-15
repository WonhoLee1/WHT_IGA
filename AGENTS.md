# WHT_IGA — Agent Rules

Kratos-based immersed IGA solver (WHTIGA work stream). These rules apply to every agent working under `WHT_IGA/`.

## Must
1. Run Python only via `micromamba run -n whtiga python ...` (or `D:\mamba\envs\whtiga\python.exe`). Python is 3.12.x.
2. Install Python packages only with uv: `uv pip install --python D:\mamba\envs\whtiga\python.exe ...`.
3. Import Kratos only through `wht_iga.kratos_backend.load_kratos()`; it pins `E:\kratos-build\install` and rejects any other Kratos.
4. Compare every numerical test against a reference value (analytic solution, reference assembler, Coreform CSV). No "it ran" tests.
5. Report numbers you actually ran. Record spikes in `docs/spikes/SNN.md`.
6. Stop and ask the user when the plan (`research-cf-IGA/master-plan_v2.1.md`) disagrees with what Kratos actually does.

## Must not
1. `conda install` / `conda create`, or `pip install` without uv.
2. Install a pip Kratos wheel into `whtiga`.
3. Import `wht_topopt` or `aieng` (copy-and-adapt with a source comment is fine). `tests/test_isolation.py` enforces this.
4. Modify `D:\PythonCodeStudy\kratos-src` or anything under `E:\kratos-build`, or rebuild Kratos.
5. Hide failures with `skip`/`xfail`, or report TODO stubs as done.
6. Commit or push unless the user asks.

## Read
- Environment (single source): [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md)
- Kratos build: [`docs/KRATOS_BUILD.md`](docs/KRATOS_BUILD.md)
- Plan: [`research-cf-IGA/master-plan_v2.1.md`](research-cf-IGA/master-plan_v2.1.md)
- Current task scope, Kratos pitfalls K1–K14, acceptance criteria: [`research-cf-IGA/implementation-request_v1.md`](research-cf-IGA/implementation-request_v1.md)
