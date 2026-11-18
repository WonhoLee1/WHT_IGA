# -*- coding: utf-8 -*-
"""Cook membrane tapered spline (Track 3, E04 precursor).

Classic Cook vertices (mm): A(0,0), B(48,44), C(48,60), D(0,44)
(left height 44, right height 16, upper-edge rise 16 — e04.json / Coreform doc).
Bilinear map from the unit square::

    x = 48 * xi
    y = 44 * xi + 44 * eta - 28 * xi * eta
    z = t * zeta

exact for any degree p >= 1 (map is degree <= 1 per direction).

Implementation: subclass :class:`BoxSpline` over the unit cube and override
only ``control_points`` (tapered positions). Knots, cells, ordering and all
parametric tables stay identical, so ``make_nurbs_volume``,
``attach_cell_elements`` and ``cell_gauss_points`` work unchanged.
"""

from __future__ import annotations

import numpy as np

from wht_iga.reference.box import BoxSpline

L_MM = 48.0
HL_MM = 44.0
HR_MM = 16.0
T_MM = 1.91516
U_Y_REF_MM = 6.98946


class CookSpline(BoxSpline):
    """Unit-cube BoxSpline with Cook-tapered control points."""

    def control_points(self) -> np.ndarray:
        X = super().control_points()
        xi, eta, zeta = X[:, 0], X[:, 1], X[:, 2]
        x = L_MM * xi
        y = HL_MM * xi + HL_MM * eta - (HL_MM - HR_MM) * xi * eta
        return np.column_stack([x, y, T_MM * zeta])


def right_face_traction_forces(ys: np.ndarray, resultant: float) -> np.ndarray:
    """Lumped nodal forces for a uniform upward traction (unit thickness scaling).

    ``ys``: physical y of the loaded-face nodes (any order). Tributary length
    from sorted neighbour midpoints; edge nodes get half intervals. Returns
    forces in the input order summing exactly to ``resultant``.
    """
    ys = np.asarray(ys, dtype=float)
    order = np.argsort(ys)
    s = ys[order]
    mid = 0.5 * (s[:-1] + s[1:])
    lo = np.concatenate([[s[0]], mid])
    hi = np.concatenate([mid, [s[-1]]])
    trib = hi - lo
    f = resultant * trib / trib.sum()
    out = np.empty_like(f)
    out[order] = f
    return out
