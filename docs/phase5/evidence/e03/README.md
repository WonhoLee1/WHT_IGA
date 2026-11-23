# E03 glass-to-metal seal — result fields and renders (2026-11-24)

Model: real STEP connector (304L header, S-8061 glass, 11 Alloy 52 pins), p=2 trimmed B-spline on the
coarse development background (9 x 6 x 6 cells, ~3 mm), RBE3 annular mount at x = -25.5 mm, NLGEOM,
Coreform accepted-increment schedule (45 increments), Kratos `tol_rel = 1e-6` (see `docs/verification/E03.md`).

Step ends: `1_cool_to_room_20C` (450 -> 20 C), `2_cool_to_minus55C`, `3_heat_to_200C`.

| File pattern | Content |
|---|---|
| `E03_parts_iso.png`, `E03_parts_cutaway_iso.png` | geometry by part; cutaway shows the glass around the pins |
| `E03_<step>_U_{x,y,z,iso}.png` | displacement magnitude on the deformed CAD surface (deformation scaled, factor in the title) |
| `E03_<step>_vonMises_{x,y,z,iso}.png` | von Mises stress, all parts |
| `E03_<step>_glass_s1_{iso,x}.png` | glass maximum principal stress (tension screening) |
| `E03_<step>_PEEQ_{iso,x}.png` | equivalent plastic strain in header and pins (glass translucent) |
| `E03_<step>_section_vonMises_{y,iso}.png` | section y = 0 |
| `E03_<step>_fields.vtp` | CAD surface with `U`, `|U| [mm]`, `von Mises [MPa]`, `max principal [MPa]`, `PEEQ`, `part` (open in ParaView) |
| `E03_<step>_qp.vtp` | Kratos integration points (reference positions) with the same fields |
| `E03_summary.json` | maxima per step |

How surface and section values are made: displacement is the spline field evaluated at the CAD
vertices. Stress and PEEQ live at Kratos integration points; a surface vertex (or section grid point
inside a solid) takes the value of the nearest integration point of the same part. On this coarse
model that gives a cell-sized patchy look; a few small blobs right of the glass in the section are
inside-test errors near the glass holes.

Check: maxima equal the gate-test snapshots (20 C: pins von Mises 252.04 MPa, glass max principal
100.85 MPa, PEEQ 0.00652).

Scripts: `e03_render2.py` (solve + capture + render), `e03_rerender.py` (re-render parts/section/PEEQ
from the VTP files), `e03_render.py` (first displacement-only pass). Paths inside point to the session
scratchpad; edit them before re-running.
