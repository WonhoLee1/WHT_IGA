"""Tensor-product B-spline discretization of an axis-aligned box (untrimmed).

Control points sit at Greville abscissae, so the geometry map is the identity and
parameter coordinates are physical coordinates (master-plan v2 §5.2.1).

Ordering
--------
* control point (a, b, c) -> index ``a + n_u * (b + n_v * c)``
* displacement DOF of control point I, component d -> ``3 * I + d``
* cell (i, j, k) -> ``i + nx * (j + ny * k)``
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from wht_iga.grid.bspline import basis_1d
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import n_control_points, uniform_clamped


@dataclass
class BoxSpline:
    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    n_cells: tuple[int, int, int]
    degree: int
    continuity: int | None = None
    knots: list[np.ndarray] = field(init=False)
    n_cp: tuple[int, int, int] = field(init=False)

    def __post_init__(self) -> None:
        self.knots = [
            uniform_clamped(self.lo[d], self.hi[d], self.n_cells[d], self.degree, self.continuity) for d in range(3)
        ]
        self.n_cp = tuple(n_control_points(k, self.degree) for k in self.knots)  # type: ignore[assignment]

    @property
    def n_ctrl(self) -> int:
        return int(np.prod(self.n_cp))

    @property
    def n_dofs(self) -> int:
        return 3 * self.n_ctrl

    def greville_axes(self) -> list[np.ndarray]:
        return [greville(k, self.degree) for k in self.knots]

    def control_points(self) -> np.ndarray:
        gu, gv, gw = self.greville_axes()
        W, V, U = np.meshgrid(gw, gv, gu, indexing="ij")
        return np.column_stack([U.ravel(), V.ravel(), W.ravel()])

    def breakpoints(self, axis: int) -> np.ndarray:
        return np.unique(self.knots[axis])

    def axis_tables(self, axis: int, gauss_x: np.ndarray):
        """Per-cell Gauss data along one axis.

        Returns (first_cp[n_cells], N[n_cells, n_gp, p+1], dN[...], x[n_cells, n_gp], jac[n_cells]).
        ``gauss_x`` are reference points on [-1, 1].
        """
        p = self.degree
        br = self.breakpoints(axis)
        nc = len(br) - 1
        ngp = len(gauss_x)
        first = np.zeros(nc, dtype=np.int64)
        N = np.zeros((nc, ngp, p + 1))
        dN = np.zeros((nc, ngp, p + 1))
        X = np.zeros((nc, ngp))
        jac = np.zeros(nc)
        for c in range(nc):
            a, b = br[c], br[c + 1]
            jac[c] = 0.5 * (b - a)
            for g, xi in enumerate(gauss_x):
                u = a + (xi + 1.0) * jac[c]
                span, ders = basis_1d(self.knots[axis], p, u, 1)
                first[c] = span - p
                N[c, g] = ders[0]
                dN[c, g] = ders[1]
                X[c, g] = u
        return first, N, dN, X, jac

    def evaluate(self, coeffs: np.ndarray, points: np.ndarray) -> np.ndarray:
        """Evaluate a vector field (n_ctrl, 3 coefficients) at physical points (n, 3)."""
        p = self.degree
        out = np.zeros((len(points), coeffs.shape[1]))
        nu, nv, _ = self.n_cp
        for n, x in enumerate(points):
            spans, vals = [], []
            for d in range(3):
                s, ders = basis_1d(self.knots[d], p, x[d], 0)
                spans.append(s - p)
                vals.append(ders[0])
            for c in range(p + 1):
                for b in range(p + 1):
                    for a in range(p + 1):
                        I = (spans[0] + a) + nu * ((spans[1] + b) + nv * (spans[2] + c))
                        out[n] += vals[0][a] * vals[1][b] * vals[2][c] * coeffs[I]
        return out

    def boundary_control_points(self, axis: int, side: int) -> np.ndarray:
        """Indices of control points on face ``axis`` = lo (side 0) or hi (side 1)."""
        nu, nv, nw = self.n_cp
        a, b, c = np.meshgrid(np.arange(nu), np.arange(nv), np.arange(nw), indexing="ij")
        idx = [a, b, c][axis]
        target = 0 if side == 0 else self.n_cp[axis] - 1
        mask = idx == target
        return (a[mask] + nu * (b[mask] + nv * c[mask])).astype(np.int64)
