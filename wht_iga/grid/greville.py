"""Greville abscissae.

With control points placed at the Greville points of a clamped knot vector the
B-spline map reproduces the identity x(u) = u (linear precision), so parameter
coordinates equal physical coordinates and det J = 1 (master-plan v2 §5.2.1).
"""

from __future__ import annotations

import numpy as np


def greville(knots: np.ndarray, degree: int) -> np.ndarray:
    knots = np.asarray(knots, dtype=float)
    p = int(degree)
    n_cp = len(knots) - p - 1
    if p == 0:
        return 0.5 * (knots[:-1] + knots[1:])
    return np.array([knots[i + 1 : i + p + 1].mean() for i in range(n_cp)])
