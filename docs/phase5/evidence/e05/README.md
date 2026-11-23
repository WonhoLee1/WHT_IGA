# E05 door-lock assembly — result fields and renders (2026-11-24)

Model: Coreform `Door_Lock_Assembly.sat` (hinge, pin with strong core and weak ends, link), three
trimmed p=2 B-spline bodies, NLGEOM + J2 hardening, three friction contact pairs (mu = 0.3); hinge
mounting holes driven to u_z = 5 mm, link pads fixed. The gate test
`tests/kratos_spikes/test_p5_e05_sat_door_lock.py` was run unchanged; `e05_render.py` only wraps the
body operator and contact step to copy each part's spline displacement.

Run (`E05_run.log`): RF3_link(5.0 mm) = -46.763 kN vs Coreform -46.946 kN (0.39 %), same as the gate
record; step 8 needed three bisections (4.350 -> 4.431 -> 4.512 -> 4.675 -> 5.000 mm).

| File | Content |
|---|---|
| `E05_parts_undeformed_iso.png`, `E05_parts_deformed_uz5.0mm_iso.png` | parts (true-scale deformation) |
| `E05_final_uz5.0mm_U_{x,y,z,iso}.png` | |U| on the deformed SAT surfaces at u_z = 5 mm (true scale) |
| `E05_mid_step_U_{x,y,z,iso}.png` | same at u_z = 3.667 mm (6th of 11 captured solves) |
| `E05_final_uz5.0mm.vtp`, `E05_mid_step.vtp` | SAT tessellation with `U`, `|U| [mm]`, `part` (0 hinge, 1 pin, 2 link) |

Stress / PEEQ fields were not captured for E05 (only displacement); the run takes about 35 min.
