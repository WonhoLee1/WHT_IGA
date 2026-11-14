"""Knot vectors for the rectilinear background grid.

Conventions
-----------
* Standard ("open/clamped") knot vectors repeat each end knot ``p + 1`` times.
  Length = n_cp + p + 1.
* Kratos ``NurbsVolumeGeometry`` uses a reduced convention that repeats each end
  knot only ``p`` times (length = n_cp + p - 1); see
  ``IgaApplication/tests/test_nurbs_volume_element.py`` where one span of degree
  p is ``[0]*p + [1]*p``.
* Knots are physical coordinates along the (possibly rotated) grid axis, which
  together with Greville control points gives x(xi) = xi (identity map, v2 §5.2.1).
"""

from __future__ import annotations

import numpy as np


def uniform_clamped(lo: float, hi: float, n_cells: int, degree: int, continuity: int | None = None) -> np.ndarray:
    """Standard clamped knot vector with ``n_cells`` equal spans on [lo, hi].

    ``continuity`` k (default p-1) sets interior knot multiplicity to p - k.
    """
    p = int(degree)
    if p < 1:
        raise ValueError("degree must be >= 1")
    if n_cells < 1:
        raise ValueError("n_cells must be >= 1")
    k = p - 1 if continuity is None else int(continuity)
    if not 0 <= k <= p - 1:
        raise ValueError(f"continuity must be in [0, {p - 1}] for degree {p}")
    breaks = np.linspace(lo, hi, n_cells + 1)
    mult = p - k
    interior = np.repeat(breaks[1:-1], mult)
    return np.concatenate([np.full(p + 1, lo), interior, np.full(p + 1, hi)])


def n_control_points(knots: np.ndarray, degree: int) -> int:
    return len(knots) - degree - 1


def to_kratos(knots: np.ndarray, degree: int) -> np.ndarray:
    """Standard clamped knots -> Kratos reduced convention (drop one end knot each side)."""
    knots = np.asarray(knots, dtype=float)
    p = int(degree)
    if not (np.all(knots[: p + 1] == knots[0]) and np.all(knots[-(p + 1):] == knots[-1])):
        raise ValueError("knot vector is not clamped with multiplicity p+1")
    return knots[1:-1].copy()


def from_kratos(knots: np.ndarray, degree: int) -> np.ndarray:
    """Kratos reduced convention -> standard clamped knots."""
    knots = np.asarray(knots, dtype=float)
    return np.concatenate([[knots[0]], knots, [knots[-1]]])


def breakpoints(knots: np.ndarray) -> np.ndarray:
    """Unique knot values (cell boundaries)."""
    return np.unique(np.asarray(knots, dtype=float))
