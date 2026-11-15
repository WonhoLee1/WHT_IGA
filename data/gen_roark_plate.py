"""Roark finite-width plate with central hole (Coreform V02, v2 §10.1).

Dimensions (``roark-17p1-7a``): L=200, D=40, r=10, t=1, tension P=1e-6.
Plate centred at origin, hole axis = z. Reproducible: rerunning this
script regenerates ``roark_plate_hole.step`` byte-equivalent geometry
(STEP headers may carry timestamps; geometry compared by tessellation).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

L, D, R, T = 200.0, 40.0, 10.0, 1.0


def build_plate(length: float = L, width: float = D,
                radius: float = R, thick: float = T):  # type: ignore[no-untyped-def]
    """Box minus through-cylinder, healed + sewn."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2

    from wht_iga.geometry import cad

    box = BRepPrimAPI_MakeBox(
        gp_Pnt(-length / 2, -width / 2, -thick / 2),
        length, width, thick).Shape()
    cyl = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0.0, 0.0, -thick), gp_Dir(0, 0, 1)),
        radius, 2 * thick).Shape()
    cut = cad.cut_shape(box, cyl)
    return cad.heal(cut)


def main() -> None:
    """Generate the Roark STEP next to this script."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "roark_plate_hole.step"))
    args = ap.parse_args()
    from wht_iga.geometry import cad

    shape = build_plate()
    cad.write_step(shape, args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
