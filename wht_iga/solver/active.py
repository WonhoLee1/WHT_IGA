"""Active/inactive basis classification for immersed grids (plan v2.1 §7.5).

A basis function is *inactive* when its support misses Omega, i.e. every cell
of its support is OUT. ``out_mode="drop"`` creates no element in OUT cells
and fixes the inactive basis to zero instead of filling OUT cells with the
FCM ghost stiffness alpha*E. Ghost fill keeps those DOFs in the system: on the
V06 sphere octant (~20 % of the box) it carried 72 % of the DOFs and 58 % of
the QPs at level 2 p=2, and the factorization outgrew RAM (2026-11-15).
It also caps the condition number near 1/alpha, which hides what V-COND
is meant to measure.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def inactive_ctrl_mask(box, out_cells: Iterable[tuple[int, int, int]]) -> np.ndarray:
    """Boolean mask (n_ctrl,) of control points whose support holds only OUT cells.

    Cell (i, j, k) of a C^(p-1) uniform box is supported by control points
    i..i+p per axis; control points are numbered x fastest
    (I = a + ncx * (b + ncy * c), as ``BoxSpline.control_points``).
    """
    p = int(box.degree)
    nx, ny, nz = (int(v) for v in box.n_cells)
    ncx, ncy, ncz = (int(v) for v in box.n_cp)
    if (ncx, ncy, ncz) != (nx + p, ny + p, nz + p):
        raise ValueError(f"inactive_ctrl_mask needs C^(p-1) knots, got n_cp={box.n_cp}")
    act_cell = np.ones((nx, ny, nz), dtype=bool)
    for c in out_cells:
        act_cell[tuple(int(v) for v in c)] = False
    act = np.zeros((ncx, ncy, ncz), dtype=bool)
    idx = np.argwhere(act_cell)
    for di in range(p + 1):
        for dj in range(p + 1):
            for dk in range(p + 1):
                act[idx[:, 0] + di, idx[:, 1] + dj, idx[:, 2] + dk] = True
    return ~act.transpose(2, 1, 0).ravel()


def fix_inactive(KM, mp, box, out_cells) -> int:
    """Fix all displacement DOFs of inactive control-point nodes to 0. Returns count.

    Control point I is Kratos node I + 1 (``make_nurbs_volume`` numbering).
    """
    inactive = np.flatnonzero(inactive_ctrl_mask(box, out_cells))
    comps = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for I in inactive:
        node = mp.GetNode(int(I) + 1)
        for var in comps:
            node.Fix(var)
            node.SetSolutionStepValue(var, 0.0)
    return int(inactive.size)
