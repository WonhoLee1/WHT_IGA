"""THB-spline hierarchical refinement engine (Giannelli, Jüttler, Speleers 2012).

Provides:
* 1D dyadic knot refinement and exact Oslo Algorithm 1 two-scale subdivision
  matrices ``S_1d`` satisfying ``N^L(u) = S_1d^T N^{L+1}(u)``.
* 3D tensor-product sparse subdivision ``S_3d = kron(S_w, kron(S_v, S_u))``
  matching the ``i + nx * (j + ny * k)`` linear indexing convention.
* Hierarchical cell mesh with bounding-box, sphere, and explicit cell marking,
  plus BFS 2:1 balance enforcement across all face/edge/vertex-adjacent leaves.
* Active HB/THB basis selection and Giannelli (2012) truncation operator
  ``C_THB`` (sparse CSR matrix) and direct THB basis/gradient evaluation.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import scipy.sparse as sp

from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.knots import breakpoints, n_control_points, uniform_clamped


def midpoint_refine_knots(knots: np.ndarray) -> np.ndarray:
    """Insert the midpoint of every non-zero knot span into ``knots``.

    Parameters
    ----------
    knots:
        1D non-decreasing knot vector (typically clamped).

    Returns
    -------
    np.ndarray
        Sorted refined knot vector with one new knot at each span midpoint.
    """
    u = np.asarray(knots, dtype=float)
    breaks = breakpoints(u)
    if len(breaks) < 2:
        return u.copy()
    midpoints = 0.5 * (breaks[:-1] + breaks[1:])
    return np.sort(np.concatenate([u, midpoints]))


def subdivision_matrix_1d(
    coarse_knots: np.ndarray,
    fine_knots: np.ndarray,
    degree: int,
) -> sp.csr_matrix:
    """Compute the exact 1D knot-insertion subdivision matrix via Oslo Algorithm 1.

    Parameters
    ----------
    coarse_knots:
        Coarse knot vector ``U^L`` of length ``n_coarse + degree + 1``.
    fine_knots:
        Refined knot vector ``U^{L+1}`` containing ``U^L`` as a subsequence,
        of length ``n_fine + degree + 1``.
    degree:
        Polynomial degree ``p >= 1``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Matrix ``S`` of shape ``(n_fine, n_coarse)`` such that
        ``N^L(u) = S.T @ N^{L+1}(u)`` and ``c^{L+1} = S @ c^L`` identically.
    """
    u = np.asarray(coarse_knots, dtype=float)
    tau = np.asarray(fine_knots, dtype=float)
    p = int(degree)
    if p < 1:
        raise ValueError("degree must be >= 1")

    n_coarse = n_control_points(u, p)
    n_fine = n_control_points(tau, p)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    for j in range(n_fine):
        t_j = tau[j]
        # Locate mu such that u[mu] <= t_j < u[mu + 1], clamped to [p, n_coarse - 1].
        if t_j >= u[n_coarse]:
            mu = n_coarse - 1
        elif t_j <= u[p]:
            mu = p
        else:
            mu = int(np.searchsorted(u, t_j, side="right")) - 1
            mu = max(p, min(mu, n_coarse - 1))

        b = np.array([1.0], dtype=float)
        for k in range(1, p + 1):
            b_new = np.zeros(k + 1, dtype=float)
            tau_jk = tau[j + k]
            for r in range(k + 1):
                i = mu - k + r
                val = 0.0
                if r >= 1:
                    denom1 = u[i + k] - u[i]
                    if denom1 > 0.0:
                        val += ((tau_jk - u[i]) / denom1) * b[r - 1]
                if r <= k - 1:
                    denom2 = u[i + k + 1] - u[i + 1]
                    if denom2 > 0.0:
                        val += ((u[i + k + 1] - tau_jk) / denom2) * b[r]
                b_new[r] = val
            b = b_new

        for r in range(p + 1):
            val = float(b[r])
            if val != 0.0:
                rows.append(j)
                cols.append(mu - p + r)
                data.append(val)

    return sp.csr_matrix((data, (rows, cols)), shape=(n_fine, n_coarse), dtype=float)


def subdivision_matrix_3d(
    su: sp.csr_matrix,
    sv: sp.csr_matrix,
    sw: sp.csr_matrix,
) -> sp.csr_matrix:
    """Build 3D tensor-product subdivision matrix matching ``i + nu * (j + nv * k)``.

    Parameters
    ----------
    su, sv, sw:
        1D subdivision matrices along ``u``, ``v``, ``w``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse matrix ``kron(sw, kron(sv, su))``.
    """
    return sp.kron(sw, sp.kron(sv, su, format="csr"), format="csr")


@dataclass(frozen=True)
class LeafCell:
    """An active leaf cell in the hierarchical mesh.

    Parameters
    ----------
    level:
        Hierarchy level ``L >= 0``.
    ix, iy, iz:
        Integer cell indices at level ``L``.
    bounds_u, bounds_v, bounds_w:
        Parametric intervals ``(lo, hi)`` along each axis.
    """

    level: int
    ix: int
    iy: int
    iz: int
    bounds_u: tuple[float, float]
    bounds_v: tuple[float, float]
    bounds_w: tuple[float, float]


def _normalize_degrees(degree: int | tuple[int, int, int]) -> tuple[int, int, int]:
    if isinstance(degree, (tuple, list)):
        if len(degree) != 3:
            raise ValueError("degree tuple must have length 3")
        return (int(degree[0]), int(degree[1]), int(degree[2]))
    p = int(degree)
    return (p, p, p)


def _prefix_sum_3d(mask: np.ndarray) -> np.ndarray:
    """Compute 3D summed-area table with zero-padded boundary."""
    nx, ny, nz = mask.shape
    sat = np.zeros((nx + 1, ny + 1, nz + 1), dtype=np.int64)
    sat[1:, 1:, 1:] = np.cumsum(
        np.cumsum(np.cumsum(mask.astype(np.int64), axis=0), axis=1),
        axis=2,
    )
    return sat


def _box_all_true(
    sat: np.ndarray,
    x0: np.ndarray,
    x1: np.ndarray,
    y0: np.ndarray,
    y1: np.ndarray,
    z0: np.ndarray,
    z1: np.ndarray,
) -> np.ndarray:
    """Vectorized check whether every cell in ``[x0:x1, y0:y1, z0:z1]`` is True."""
    ix0 = x0[:, None, None]
    ix1 = x1[:, None, None]
    iy0 = y0[None, :, None]
    iy1 = y1[None, :, None]
    iz0 = z0[None, None, :]
    iz1 = z1[None, None, :]

    sums = (
        sat[ix1, iy1, iz1]
        - sat[ix0, iy1, iz1]
        - sat[ix1, iy0, iz1]
        - sat[ix1, iy1, iz0]
        + sat[ix0, iy0, iz1]
        + sat[ix0, iy1, iz0]
        + sat[ix1, iy0, iz0]
        - sat[ix0, iy0, iz0]
    )
    expected = (ix1 - ix0) * (iy1 - iy0) * (iz1 - iz0)
    return (sums == expected) & (expected > 0)


class THBHierarchy:
    """Truncated Hierarchical B-spline (THB-spline) refinement engine in 3D.

    Parameters
    ----------
    knots_u, knots_v, knots_w:
        Level-0 clamped knot vectors along ``u``, ``v``, ``w``.
    degree:
        Polynomial degree ``p`` (scalar or ``(pu, pv, pw)``).
    """

    def __init__(
        self,
        knots_u: np.ndarray,
        knots_v: np.ndarray,
        knots_w: np.ndarray,
        degree: int | tuple[int, int, int],
    ) -> None:
        self.pu, self.pv, self.pw = _normalize_degrees(degree)
        ku0 = np.asarray(knots_u, dtype=float).copy()
        kv0 = np.asarray(knots_v, dtype=float).copy()
        kw0 = np.asarray(knots_w, dtype=float).copy()

        self._knots_u: list[np.ndarray] = [ku0]
        self._knots_v: list[np.ndarray] = [kv0]
        self._knots_w: list[np.ndarray] = [kw0]

        self._breaks_u: list[np.ndarray] = [breakpoints(ku0)]
        self._breaks_v: list[np.ndarray] = [breakpoints(kv0)]
        self._breaks_w: list[np.ndarray] = [breakpoints(kw0)]

        self._s1d_u: list[sp.csr_matrix] = []
        self._s1d_v: list[sp.csr_matrix] = []
        self._s1d_w: list[sp.csr_matrix] = []
        self._s3d: list[sp.csr_matrix] = []

        nx0 = len(self._breaks_u[0]) - 1
        ny0 = len(self._breaks_v[0]) - 1
        nz0 = len(self._breaks_w[0]) - 1
        if nx0 < 1 or ny0 < 1 or nz0 < 1:
            raise ValueError("Each axis must have at least 1 non-zero knot span")

        self._omega_mask: list[np.ndarray] = [np.ones((nx0, ny0, nz0), dtype=bool)]
        self._refined: list[np.ndarray] = [np.zeros((nx0, ny0, nz0), dtype=bool)]

    @classmethod
    def from_uniform(
        cls,
        bounds_u: tuple[float, float],
        bounds_v: tuple[float, float],
        bounds_w: tuple[float, float],
        n_cells: int | tuple[int, int, int],
        degree: int | tuple[int, int, int],
    ) -> THBHierarchy:
        """Construct a hierarchy from uniform clamped level-0 knot vectors.

        Parameters
        ----------
        bounds_u, bounds_v, bounds_w:
            Parametric domain ``[lo, hi]`` along each axis.
        n_cells:
            Number of level-0 cells per axis (scalar or ``(nx, ny, nz)``).
        degree:
            Polynomial degree per axis (scalar or ``(pu, pv, pw)``).
        """
        if isinstance(n_cells, (tuple, list)):
            nx, ny, nz = int(n_cells[0]), int(n_cells[1]), int(n_cells[2])
        else:
            nx = ny = nz = int(n_cells)
        pu, pv, pw = _normalize_degrees(degree)
        ku = uniform_clamped(bounds_u[0], bounds_u[1], nx, pu)
        kv = uniform_clamped(bounds_v[0], bounds_v[1], ny, pv)
        kw = uniform_clamped(bounds_w[0], bounds_w[1], nz, pw)
        return cls(ku, kv, kw, (pu, pv, pw))

    @property
    def max_level(self) -> int:
        """Finest allocated level index ``L_max >= 0``."""
        return len(self._omega_mask) - 1

    @property
    def degrees(self) -> tuple[int, int, int]:
        """Polynomial degrees ``(pu, pv, pw)``."""
        return (self.pu, self.pv, self.pw)

    def cell_shape(self, level: int) -> tuple[int, int, int]:
        """Number of cells ``(nx, ny, nz)`` at ``level``."""
        self.ensure_level(level)
        return self._omega_mask[level].shape

    def dof_shape(self, level: int) -> tuple[int, int, int]:
        """Number of B-spline DOFs ``(nu, nv, nw)`` at ``level``."""
        self.ensure_level(level)
        return (
            n_control_points(self._knots_u[level], self.pu),
            n_control_points(self._knots_v[level], self.pv),
            n_control_points(self._knots_w[level], self.pw),
        )

    @property
    def n_fine_dofs(self) -> int:
        """Total number of B-spline DOFs at ``max_level``."""
        nu, nv, nw = self.dof_shape(self.max_level)
        return nu * nv * nw

    @property
    def n_thb_dofs(self) -> int:
        """Total number of active THB-spline DOFs across all levels."""
        return sum(len(self.active_dofs_at_level(lvl)) for lvl in range(self.max_level + 1))

    def ensure_level(self, level: int) -> None:
        """Allocate knot vectors, subdivision matrices, and masks up to ``level``."""
        if level < 0:
            raise ValueError("level must be >= 0")
        while self.max_level < level:
            prev = self.max_level
            ku_next = midpoint_refine_knots(self._knots_u[prev])
            kv_next = midpoint_refine_knots(self._knots_v[prev])
            kw_next = midpoint_refine_knots(self._knots_w[prev])

            su = subdivision_matrix_1d(self._knots_u[prev], ku_next, self.pu)
            sv = subdivision_matrix_1d(self._knots_v[prev], kv_next, self.pv)
            sw = subdivision_matrix_1d(self._knots_w[prev], kw_next, self.pw)
            s3d = subdivision_matrix_3d(su, sv, sw)

            self._knots_u.append(ku_next)
            self._knots_v.append(kv_next)
            self._knots_w.append(kw_next)
            self._breaks_u.append(breakpoints(ku_next))
            self._breaks_v.append(breakpoints(kv_next))
            self._breaks_w.append(breakpoints(kw_next))
            self._s1d_u.append(su)
            self._s1d_v.append(sv)
            self._s1d_w.append(sw)
            self._s3d.append(s3d)

            nx_next = len(self._breaks_u[-1]) - 1
            ny_next = len(self._breaks_v[-1]) - 1
            nz_next = len(self._breaks_w[-1]) - 1
            self._omega_mask.append(np.zeros((nx_next, ny_next, nz_next), dtype=bool))
            self._refined.append(np.zeros((nx_next, ny_next, nz_next), dtype=bool))

    def _sync_omega_from_refined(self) -> None:
        """Recompute ``omega_mask[L+1]`` from ``refined[L]`` for all levels."""
        for lvl in range(self.max_level):
            ref = self._refined[lvl] & self._omega_mask[lvl]
            self._refined[lvl] = ref
            child_mask = np.repeat(np.repeat(np.repeat(ref, 2, axis=0), 2, axis=1), 2, axis=2)
            self._omega_mask[lvl + 1] = child_mask
        # Finest level cannot be refined further.
        self._refined[self.max_level][:] = False

    def mark_box(
        self,
        lo: tuple[float, float, float],
        hi: tuple[float, float, float],
        target_level: int,
    ) -> int:
        """Mark all cells overlapping ``[lo, hi]`` for refinement up to ``target_level``.

        Parameters
        ----------
        lo, hi:
            Axis-aligned bounding box ``(u, v, w)`` bounds.
        target_level:
            Desired refinement level ``>= 0``.

        Returns
        -------
        int
            Number of newly marked cell refinements across levels ``0 .. target_level - 1``.
        """
        if target_level <= 0:
            return 0
        self.ensure_level(target_level)
        lo_arr = np.asarray(lo, dtype=float)
        hi_arr = np.asarray(hi, dtype=float)
        eps = 1e-13
        added = 0

        for lvl in range(target_level):
            bu = self._breaks_u[lvl]
            bv = self._breaks_v[lvl]
            bw = self._breaks_w[lvl]

            ux_lo, ux_hi = bu[:-1], bu[1:]
            vy_lo, vy_hi = bv[:-1], bv[1:]
            wz_lo, wz_hi = bw[:-1], bw[1:]

            if hi_arr[0] > lo_arr[0] + eps:
                mx = (ux_hi > lo_arr[0] + eps) & (ux_lo < hi_arr[0] - eps)
            else:
                mx = (ux_lo <= lo_arr[0] + eps) & (lo_arr[0] - eps <= ux_hi)

            if hi_arr[1] > lo_arr[1] + eps:
                my = (vy_hi > lo_arr[1] + eps) & (vy_lo < hi_arr[1] - eps)
            else:
                my = (vy_lo <= lo_arr[1] + eps) & (lo_arr[1] - eps <= vy_hi)

            if hi_arr[2] > lo_arr[2] + eps:
                mz = (wz_hi > lo_arr[2] + eps) & (wz_lo < hi_arr[2] - eps)
            else:
                mz = (wz_lo <= lo_arr[2] + eps) & (lo_arr[2] - eps <= wz_hi)

            box_mask = mx[:, None, None] & my[None, :, None] & mz[None, None, :]
            new_marks = box_mask & self._omega_mask[lvl] & (~self._refined[lvl])
            added += int(np.count_nonzero(new_marks))
            self._refined[lvl] |= new_marks
            self._sync_omega_from_refined()

        return added

    def mark_sphere(
        self,
        center: tuple[float, float, float],
        radius: float,
        target_level: int,
    ) -> int:
        """Mark all cells intersecting the closed ball ``||x - center|| <= radius``.

        Parameters
        ----------
        center:
            Sphere center ``(u, v, w)``.
        radius:
            Sphere radius ``>= 0``.
        target_level:
            Desired refinement level ``>= 0``.

        Returns
        -------
        int
            Number of newly marked cell refinements.
        """
        if target_level <= 0 or radius < 0.0:
            return 0
        self.ensure_level(target_level)
        cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
        r2 = float(radius) ** 2
        added = 0

        for lvl in range(target_level):
            bu = self._breaks_u[lvl]
            bv = self._breaks_v[lvl]
            bw = self._breaks_w[lvl]

            dx = np.maximum(0.0, np.maximum(bu[:-1] - cx, cx - bu[1:]))
            dy = np.maximum(0.0, np.maximum(bv[:-1] - cy, cy - bv[1:]))
            dz = np.maximum(0.0, np.maximum(bw[:-1] - cz, cz - bw[1:]))

            dist2 = (dx**2)[:, None, None] + (dy**2)[None, :, None] + (dz**2)[None, None, :]
            sphere_mask = dist2 <= r2 + 1e-14
            new_marks = sphere_mask & self._omega_mask[lvl] & (~self._refined[lvl])
            added += int(np.count_nonzero(new_marks))
            self._refined[lvl] |= new_marks
            self._sync_omega_from_refined()

        return added

    def mark_cells(
        self,
        cells: Sequence[tuple[int, int, int, int]],
        target_level: int | None = None,
    ) -> int:
        """Mark explicit ``(level, ix, iy, iz)`` cells for refinement.

        Parameters
        ----------
        cells:
            Sequence of ``(level, ix, iy, iz)`` tuples.
        target_level:
            If given, refine descendants of each cell down to ``target_level``.
            If ``None``, refine each given cell by one level (``level + 1``).

        Returns
        -------
        int
            Number of newly marked cell refinements.
        """
        if not cells:
            return 0
        max_req = max((lvl + 1 if target_level is None else target_level) for lvl, _, _, _ in cells)
        self.ensure_level(max_req)
        added = 0

        for lvl, ix, iy, iz in cells:
            t_lvl = lvl + 1 if target_level is None else target_level
            if t_lvl <= lvl:
                continue
            if not self._omega_mask[lvl][ix, iy, iz]:
                raise ValueError(f"Cell ({lvl}, {ix}, {iy}, {iz}) is not in Omega^{lvl}")
            # Mark cell at lvl, then its descendants up to t_lvl - 1.
            x0, x1 = ix, ix + 1
            y0, y1 = iy, iy + 1
            z0, z1 = iz, iz + 1
            for cur in range(lvl, t_lvl):
                sub = self._refined[cur][x0:x1, y0:y1, z0:z1]
                added += int(np.count_nonzero(~sub))
                sub[:] = True
                self._sync_omega_from_refined()
                x0, x1 = 2 * x0, 2 * x1
                y0, y1 = 2 * y0, 2 * y1
                z0, z1 = 2 * z0, 2 * z1

        return added

    def enforce_2_to_1_balance(self) -> int:
        """Enforce 2:1 level balance across all face/edge/vertex-adjacent leaves via BFS.

        A cell ``c`` at level ``L >= 1`` that is refined requires every 26-neighbor
        of ``c`` at level ``L`` to belong to ``Omega^L``, i.e. that neighbor's
        parent at level ``L - 1`` must be refined.

        Returns
        -------
        int
            Number of additional cells refined to restore 2:1 balance.
        """
        if self.max_level <= 1:
            return 0

        queue: deque[tuple[int, int, int, int]] = deque()
        for lvl in range(1, self.max_level):
            xs, ys, zs = np.nonzero(self._refined[lvl])
            for ix, iy, iz in zip(xs.tolist(), ys.tolist(), zs.tolist()):
                queue.append((lvl, ix, iy, iz))

        added = 0
        while queue:
            lvl, ix, iy, iz = queue.popleft()
            if lvl < 1:
                continue
            nx_l, ny_l, nz_l = self._omega_mask[lvl].shape
            parent_lvl = lvl - 1
            for dx in (-1, 0, 1):
                nx = ix + dx
                if not 0 <= nx < nx_l:
                    continue
                px = nx // 2
                for dy in (-1, 0, 1):
                    ny = iy + dy
                    if not 0 <= ny < ny_l:
                        continue
                    py = ny // 2
                    for dz in (-1, 0, 1):
                        if dx == 0 and dy == 0 and dz == 0:
                            continue
                        nz = iz + dz
                        if not 0 <= nz < nz_l:
                            continue
                        pz = nz // 2
                        if not self._refined[parent_lvl][px, py, pz]:
                            self._refined[parent_lvl][px, py, pz] = True
                            added += 1
                            if parent_lvl >= 1:
                                queue.append((parent_lvl, px, py, pz))

        if added > 0:
            self._sync_omega_from_refined()
        return added

    def leaf_cells(self) -> list[LeafCell]:
        """Return all active leaf cells ordered by level and ``i + nx * (j + ny * k)``."""
        leaves: list[LeafCell] = []
        for lvl in range(self.max_level + 1):
            leaf_mask = self._omega_mask[lvl] & (~self._refined[lvl])
            nx, ny, nz = leaf_mask.shape
            bu = self._breaks_u[lvl]
            bv = self._breaks_v[lvl]
            bw = self._breaks_w[lvl]
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        if leaf_mask[ix, iy, iz]:
                            leaves.append(
                                LeafCell(
                                    level=lvl,
                                    ix=ix,
                                    iy=iy,
                                    iz=iz,
                                    bounds_u=(float(bu[ix]), float(bu[ix + 1])),
                                    bounds_v=(float(bv[iy]), float(bv[iy + 1])),
                                    bounds_w=(float(bw[iz]), float(bw[iz + 1])),
                                )
                            )
        return leaves

    def active_dofs_at_level(self, level: int) -> np.ndarray:
        """Return sorted 1D array of linear DOF indices active at ``level``.

        Criterion: ``supp(B_i^L) <= Omega^L`` and ``supp(B_i^L) not <= Omega^{L+1}``.
        """
        self.ensure_level(level)
        nu, nv, nw = self.dof_shape(level)
        ku = self._knots_u[level]
        kv = self._knots_v[level]
        kw = self._knots_w[level]
        bu = self._breaks_u[level]
        bv = self._breaks_v[level]
        bw = self._breaks_w[level]

        x0 = np.searchsorted(bu, ku[:nu])
        x1 = np.searchsorted(bu, ku[self.pu + 1 : nu + self.pu + 1])
        y0 = np.searchsorted(bv, kv[:nv])
        y1 = np.searchsorted(bv, kv[self.pv + 1 : nv + self.pv + 1])
        z0 = np.searchsorted(bw, kw[:nw])
        z1 = np.searchsorted(bw, kw[self.pw + 1 : nw + self.pw + 1])

        sat_omega = _prefix_sum_3d(self._omega_mask[level])
        in_omega = _box_all_true(sat_omega, x0, x1, y0, y1, z0, z1)

        if level < self.max_level:
            sat_ref = _prefix_sum_3d(self._refined[level])
            in_next = _box_all_true(sat_ref, x0, x1, y0, y1, z0, z1)
            active_3d = in_omega & (~in_next)
        else:
            active_3d = in_omega

        # Linear index: i + nu * (j + nv * k), i.e. Fortran-order flattening.
        return np.flatnonzero(active_3d.ravel(order="F")).astype(np.int64)

    def active_dofs(self) -> list[tuple[int, int]]:
        """Return ordered list of ``(level, linear_dof_index)`` for all THB DOFs."""
        out: list[tuple[int, int]] = []
        for lvl in range(self.max_level + 1):
            for idx in self.active_dofs_at_level(lvl).tolist():
                out.append((lvl, idx))
        return out

    def build_c_thb_matrix(self) -> sp.csr_matrix:
        """Assemble exact Giannelli (2012) truncation/prolongation matrix ``C_THB``.

        Returns
        -------
        scipy.sparse.csr_matrix
            Matrix of shape ``(n_fine_dofs, n_thb_dofs)`` mapping global THB
            coefficients to ``max_level`` B-spline coefficients, such that
            ``N_THB(xi) = C_THB.T @ N_fine(xi)``.
        """
        l_max = self.max_level
        n0 = int(np.prod(self.dof_shape(0)))
        act0 = self.active_dofs_at_level(0)
        c_mat = sp.csr_matrix(
            (np.ones(len(act0), dtype=float), (act0, np.arange(len(act0), dtype=np.int64))),
            shape=(n0, len(act0)),
        )

        for lvl in range(1, l_max + 1):
            nl = int(np.prod(self.dof_shape(lvl)))
            act_l = self.active_dofs_at_level(lvl)

            # Truncate active DOFs at level `lvl` from coarse contributions.
            mask_diag = np.ones(nl, dtype=float)
            mask_diag[act_l] = 0.0
            m_l = sp.diags(mask_diag, offsets=0, shape=(nl, nl), format="csr")

            p_l = sp.csr_matrix(
                (np.ones(len(act_l), dtype=float), (act_l, np.arange(len(act_l), dtype=np.int64))),
                shape=(nl, len(act_l)),
            )

            c_mat = m_l @ (self._s3d[lvl - 1] @ c_mat)
            c_mat = sp.hstack([c_mat, p_l], format="csr")

        c_mat.eliminate_zeros()
        return c_mat

    def _eval_fine_basis_sparse(
        self,
        points_uvw: np.ndarray,
        n_ders: int = 0,
    ) -> list[sp.csr_matrix]:
        """Evaluate ``max_level`` B-spline basis (and derivatives) at ``points_uvw``."""
        pts = np.atleast_2d(np.asarray(points_uvw, dtype=float))
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError("points_uvw must have shape (n_pts, 3)")

        l_max = self.max_level
        ku = self._knots_u[l_max]
        kv = self._knots_v[l_max]
        kw = self._knots_w[l_max]
        nu, nv, nw = self.dof_shape(l_max)
        n_fine = nu * nv * nw
        n_pts = pts.shape[0]

        n_comp = 1 if n_ders == 0 else 4
        rows_list: list[list[int]] = [[] for _ in range(n_comp)]
        cols_list: list[list[int]] = [[] for _ in range(n_comp)]
        vals_list: list[list[float]] = [[] for _ in range(n_comp)]

        for p_idx in range(n_pts):
            u_val, v_val, w_val = float(pts[p_idx, 0]), float(pts[p_idx, 1]), float(pts[p_idx, 2])
            su = find_span(nu, self.pu, u_val, ku)
            sv = find_span(nv, self.pv, v_val, kv)
            sw = find_span(nw, self.pw, w_val, kw)

            du = ders_basis_funs(su, u_val, self.pu, n_ders, ku)
            dv = ders_basis_funs(sv, v_val, self.pv, n_ders, kv)
            dw = ders_basis_funs(sw, w_val, self.pw, n_ders, kw)

            iu = np.arange(su - self.pu, su + 1, dtype=np.int64)
            iv = np.arange(sv - self.pv, sv + 1, dtype=np.int64)
            iw = np.arange(sw - self.pw, sw + 1, dtype=np.int64)

            col_idx = (
                iu[:, None, None]
                + nu * (iv[None, :, None] + nv * iw[None, None, :])
            ).ravel()

            b0 = (du[0, :, None, None] * dv[0, None, :, None] * dw[0, None, None, :]).ravel()
            rows_list[0].extend([p_idx] * len(col_idx))
            cols_list[0].extend(col_idx.tolist())
            vals_list[0].extend(b0.tolist())

            if n_ders >= 1:
                bx = (du[1, :, None, None] * dv[0, None, :, None] * dw[0, None, None, :]).ravel()
                by = (du[0, :, None, None] * dv[1, None, :, None] * dw[0, None, None, :]).ravel()
                bz = (du[0, :, None, None] * dv[0, None, :, None] * dw[1, None, None, :]).ravel()
                for c_idx, arr in ((1, bx), (2, by), (3, bz)):
                    rows_list[c_idx].extend([p_idx] * len(col_idx))
                    cols_list[c_idx].extend(col_idx.tolist())
                    vals_list[c_idx].extend(arr.tolist())

        return [
            sp.csr_matrix((vals_list[c], (rows_list[c], cols_list[c])), shape=(n_pts, n_fine))
            for c in range(n_comp)
        ]

    def evaluate_thb_basis(self, points_uvw: np.ndarray) -> np.ndarray:
        """Evaluate all THB basis functions at ``points_uvw``.

        Parameters
        ----------
        points_uvw:
            Parametric coordinates of shape ``(n_pts, 3)`` or ``(3,)``.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_pts, n_thb_dofs)`` containing ``N_THB(u, v, w)``.
        """
        c_thb = self.build_c_thb_matrix()
        b_fine = self._eval_fine_basis_sparse(points_uvw, n_ders=0)[0]
        return (b_fine @ c_thb).toarray()

    def evaluate_thb_basis_and_grad(self, points_uvw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate THB basis functions and parametric gradients at ``points_uvw``.

        Parameters
        ----------
        points_uvw:
            Parametric coordinates of shape ``(n_pts, 3)`` or ``(3,)``.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            ``(values, grads)`` with shapes ``(n_pts, n_thb_dofs)`` and
            ``(n_pts, n_thb_dofs, 3)``.
        """
        c_thb = self.build_c_thb_matrix()
        mats = self._eval_fine_basis_sparse(points_uvw, n_ders=1)
        vals = (mats[0] @ c_thb).toarray()
        du = (mats[1] @ c_thb).toarray()
        dv = (mats[2] @ c_thb).toarray()
        dw = (mats[3] @ c_thb).toarray()
        grads = np.stack([du, dv, dw], axis=-1)
        return vals, grads

    def project_coarse_coefficients(self, coarse_coeffs: np.ndarray) -> np.ndarray:
        """Map level-0 B-spline coefficients to exact THB coefficients.

        Parameters
        ----------
        coarse_coeffs:
            1D array of length ``prod(dof_shape(0))``.

        Returns
        -------
        np.ndarray
            1D array of length ``n_thb_dofs`` such that
            ``C_THB @ thb_coeffs == fine_coeffs`` for any field in ``V^0``.
        """
        c = np.asarray(coarse_coeffs, dtype=float)
        parts: list[np.ndarray] = [c[self.active_dofs_at_level(0)]]
        for lvl in range(1, self.max_level + 1):
            c = self._s3d[lvl - 1] @ c
            parts.append(c[self.active_dofs_at_level(lvl)])
        return np.concatenate(parts)
