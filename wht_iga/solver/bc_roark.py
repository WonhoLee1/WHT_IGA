# -*- coding: utf-8 -*-
"""V02 Roark-coupon boundary conditions: shared single source (T1).

Option (a) of the fix order: left face pinned (strong Fix, all components),
tension via uniform pressure on the right (+x) face. No uy/uz constraints on
free edges or z-faces: with nu = 0 they are exactly zero anyway, and the
former per-component penalty set spuriously stiffened the transverse
direction (G1 V02 -23% root cause).

Used by: WHT_IGA/scripts/v02_roark.py, tests/kratos_spikes/test_p2mesh_v02.py,
tests/kratos_spikes/test_p2stab_v02_fcm.py. Do NOT copy these BCs elsewhere.
"""

from __future__ import annotations


def roark_left_pinned(pressure: float = 1e-6) -> dict:
    """JobSpec fragments for the Coreform V02 load case.

    ``+pressure`` pulls outward (+x traction); the pinned face carries the
    reaction. Returns {"fixed_faces": ..., "pressures": ...}.
    """
    from wht_iga.solver.validate import PressureLoad

    return {
        "fixed_faces": [(0, 0)],
        "pressures": [PressureLoad(face=(0, 1), magnitude=float(pressure))],
    }
