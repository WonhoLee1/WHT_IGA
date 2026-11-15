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

from numba import njit

from wht_iga.grid.bspline import basis_1d, ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import n_control_points, uniform_clamped


@njit(cache=True)
def _eval_box_field_njit(
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
    p: int,
    nu: int,
    nv: int,
    nw: int,
    coeffs: np.ndarray,
    points: np.ndarray,
) -> np.ndarray:
    """Batch-evaluate a trivariate B-spline vector field at ``points`` (M, 3)."""
    m = points.shape[0]
    ncomp = coeffs.shape[1]
    out = np.zeros((m, ncomp))
    for n in range(m):
        px = points[n, 0]
        py = points[n, 1]
        pz = points[n, 2]
        sx = find_span(nu, p, px, kx)
        sy = find_span(nv, p, py, ky)
        sz = find_span(nw, p, pz, kz)
        Nx = ders_basis_funs(sx, px, p, 0, kx)[0]
        Ny = ders_basis_funs(sy, py, p, 0, ky)[0]
        Nz = ders_basis_funs(sz, pz, p, 0, kz)[0]
        ix0 = sx - p
        iy0 = sy - p
        iz0 = sz - p
        for c in range(p + 1):
            wz = Nz[c]
            zbase = (iz0 + c) * nv
            for b in range(p + 1):
                wyz = Ny[b] * wz
                yzbase = nu * (iy0 + b + zbase)
                for a in range(p + 1):
                    w = Nx[a] * wyz
                    I = ix0 + a + yzbase
                    for d in range(ncomp):
                        out[n, d] += w * coeffs[I, d]
    return out


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
        pts = np.ascontiguousarray(points, dtype=np.float64)
        cf = np.ascontiguousarray(coeffs, dtype=np.float64)
        nu, nv, nw = (int(v) for v in self.n_cp)
        return _eval_box_field_njit(
            np.ascontiguousarray(self.knots[0], dtype=np.float64),
            np.ascontiguousarray(self.knots[1], dtype=np.float64),
            np.ascontiguousarray(self.knots[2], dtype=np.float64),
            int(self.degree), nu, nv, nw, cf, pts,
        )

    def boundary_control_points(self, axis: int, side: int) -> np.ndarray:
        """Indices of control points on face ``axis`` = lo (side 0) or hi (side 1)."""
        nu, nv, nw = self.n_cp
        a, b, c = np.meshgrid(np.arange(nu), np.arange(nv), np.arange(nw), indexing="ij")
        idx = [a, b, c][axis]
        target = 0 if side == 0 else self.n_cp[axis] - 1
        mask = idx == target
        return (a[mask] + nu * (b[mask] + nv * c[mask])).astype(np.int64)
