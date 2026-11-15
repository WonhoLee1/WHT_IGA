# V02-ext — Roark plate peak with extension ON

> **최종 상태 (2026-11-15 G2 COMPLETE): PASS.** 아래 본문은 STAB 대기 시점(PENDING) 기록으로 보존.
> 최종 근거: `docs/phase1/evidence/v02_cf3_off` (ext OFF: finest h=0.375 peak_A=4.3466e-6,
> rel_err_A=0.674%, JSON 직접 확인) 및 `v02_cf3_0.3` (θ=0.3 ext ON: 0.678%),
> 분포 rate 2.34 — `G2_COMPLETE.md` §게이트 판정.

- Status (STAB 대기 시점): **PENDING integrator rerun** (operators landed; extension hook in
  `run_linear` ready; BC fixed to left-face pinned — old per-component BCs
  are void, see G2).
- Gate (G2): V02 p = 2 finest-grid boundary peak relative error <= 1%
  vs 4.3175e-6 (`roark-17p1-7a-accepted-peak.csv`), observed rate >= 1.8.
- Phase-1 baseline (G1, OLD BCs): finest 22.9%, rate 0.02 — root cause was
  the wrong BC set (free-edge uy/uz constraints), NOT small cells.
  See `WHT_IGA/docs/G1_report.md`.

## Input summary

- Same as `WHT_IGA/scripts/v02_roark.py`: plate bbox grid, p = 2,
  h = {10, 5, 2.5}, E = 1, nu = 0, P = 1e-6 tension, GEO cut QP tables
  via `trimmed.build_trimmed_table` (q = 4), boundary peak via
  `trimmed.boundary_peak_stress` on the hole CAD face.
- Extension: Track STAB T-P2.1/T-P2.2/T-P2.3 operators applied to the
  crown-zone small-cut band (support-volume-ratio criterion, theta = 0.2).

## Why pending (update: STAB landed, rerun unblocked but not run by MESH)

- Track STAB has since landed `wht_iga/stabilize/` and `wht_iga/operators/`
  (read-only API check: `compose_thb`, `extension_matrix`, `apply_TKT`).
  Not executed by Track MESH (STAB files are read-only for this track).

## Mesh-file input path (T-P2.6, separate from extension)

- The O-grid mesh-file input rerun is Track MESH scope and is recorded in
  `WHT_IGA/docs/phase2/P26.md` (tests `test_p2mesh_v02.py::
  test_v02_mesh_input_peak[_fine]`): 64-gon peak 1.89e-06, 128-gon 1.95e-06.
  Reference these against the Coreform accepted value (4.3175e-6), NOT
  against any intermediate CAD-path number. The mesh path is consistent
  machinery with faceted-geometry error; the Roark gate needs extension +
  curved-geometry handling.

## Re-run (for the integrator, once STAB tests pass)

```powershell
C:\Users\GOODMAN\AppData\Local\micromamba\micromamba.exe run -n whtiga python WHT_IGA/scripts/v02_roark.py
```

(with extension enabled; plus the mesh-input variant in
`tests/kratos_spikes/test_p2mesh_v02.py`).

## Result

- NOT RUN by Track MESH (extension path). No numbers claimed.
