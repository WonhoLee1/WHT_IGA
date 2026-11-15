"""Consistent mass matrix assembly and shift-invert eigenmode solver for 3D B-spline IGA.

Provides Kratos-independent assembly of the 3D consistent mass matrix
:math:`M = \\int_\\Omega \\rho N^T N \\, dV` from either general quadrature data
or axis-aligned B-spline tensor-product grids, plus shift-invert generalized
eigenproblem solving with mass-orthonormalization and modal participation factors.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.knots import breakpoints, n_control_points

DensityInput = float | Callable[[np.ndarray], np.ndarray | float]


@dataclass(frozen=True)
class EigenmodeResult:
    """Result container for generalized structural eigenmode analysis.

    Attributes:
        frequencies_hz: Natural frequencies in Hz, shape ``(n_modes,)``.
        omega: Circular natural frequencies in rad/s, shape ``(n_modes,)``.
        eigenvalues: Generalized eigenvalues ``omega**2``, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes expanded to full DOF space,
            shape ``(n_dofs, n_modes)``, satisfying ``modes.T @ M @ modes = I``.
        participation_factors: Modal participation factors ``Phi.T @ M @ r_d``
            for rigid-body translations in ``(x, y, z)``, shape ``(n_modes, 3)``.
        effective_modal_masses: Effective modal masses ``Gamma**2`` per mode and
            direction, shape ``(n_modes, 3)``.
        total_mass: Total translational mass ``r_d.T @ M @ r_d`` in ``(x, y, z)``,
            shape ``(3,)``.
        free_dofs: Sorted array of free DOF indices used in the solve.
    """

    frequencies_hz: np.ndarray
    omega: np.ndarray
    eigenvalues: np.ndarray
    modes: np.ndarray
    participation_factors: np.ndarray
    effective_modal_masses: np.ndarray
    total_mass: np.ndarray
    free_dofs: np.ndarray


def _evaluate_density(rho: DensityInput, points: np.ndarray) -> np.ndarray:
    """Evaluate scalar or spatially varying mass density at quadrature points.

    Args:
        rho: Constant scalar density or callable taking ``(n_pts, 3)`` coordinates.
        points: Physical quadrature coordinates of shape ``(n_pts, 3)``.

    Returns:
        Density array of shape ``(n_pts,)``.
    """
    if callable(rho):
        vals = np.asarray(rho(points), dtype=float)
        if vals.ndim == 0:
            return np.full(points.shape[0], float(vals), dtype=float)
        return vals.reshape(-1)
    return np.full(points.shape[0], float(rho), dtype=float)


def assemble_consistent_mass(
    cell_points: np.ndarray | Sequence[np.ndarray],
    weights_detj: np.ndarray | Sequence[np.ndarray],
    basis_values: np.ndarray | Sequence[np.ndarray],
    dof_indices: np.ndarray | Sequence[np.ndarray],
    rho: DensityInput,
    n_dofs: int | None = None,
) -> sp.csr_matrix:
    """Assemble 3D B-spline consistent mass matrix from quadrature data.

    Computes :math:`M = \\int_\\Omega \\rho(\\mathbf{x}) N^T N \\, dV` with
    3 displacement DOFs ``(u_x, u_y, u_z)`` per control point.

    Args:
        cell_points: Quadrature point coordinates, either ``(n_cells, n_gp, 3)``,
            ``(n_gp_total, 3)``, or a sequence of per-cell ``(n_gp_c, 3)`` arrays.
        weights_detj: Quadrature weights multiplied by ``|det J|``, matching the
            cell/point structure of ``cell_points``.
        basis_values: Scalar B-spline basis values ``N_a(x_q)`` of shape
            ``(n_cells, n_gp, n_loc)``, ``(n_gp_total, n_loc)``, or sequence of
            ``(n_gp_c, n_loc_c)`` arrays.
        dof_indices: Either control-point indices ``(..., n_loc)`` or 3-DOF
            indices ``(..., 3 * n_loc)`` ordered as ``[3*I, 3*I+1, 3*I+2, ...]``.
        rho: Constant mass density or callable ``rho(points) -> density``.
        n_dofs: Total number of DOFs in the global matrix. If ``None``, inferred
            as ``max_dof + 1``.

    Returns:
        Symmetric positive semi-definite sparse consistent mass matrix in CSR
        format of shape ``(n_dofs, n_dofs)``.
    """
    if isinstance(weights_detj, np.ndarray) and weights_detj.ndim == 1:
        pts_arr = np.asarray(cell_points, dtype=float)
        w_arr = np.asarray(weights_detj, dtype=float)
        n_arr = np.asarray(basis_values, dtype=float)
        d_arr = np.asarray(dof_indices, dtype=np.int64)
        if pts_arr.ndim == 2:
            pts_seq: Sequence[np.ndarray] = [pts_arr[i : i + 1] for i in range(len(w_arr))]
            w_seq: Sequence[np.ndarray] = [w_arr[i : i + 1] for i in range(len(w_arr))]
            n_seq: Sequence[np.ndarray] = [n_arr[i : i + 1] for i in range(len(w_arr))]
            d_seq: Sequence[np.ndarray] = [d_arr[i] for i in range(len(w_arr))]
        else:
            pts_seq = list(pts_arr)
            w_seq = list(w_arr)
            n_seq = list(n_arr)
            d_seq = list(d_arr)
    else:
        pts_seq = [np.asarray(p, dtype=float) for p in cell_points]
        w_seq = [np.asarray(w, dtype=float) for w in weights_detj]
        n_seq = [np.asarray(n, dtype=float) for n in basis_values]
        d_seq = [np.asarray(d, dtype=np.int64) for d in dof_indices]

    rows_list: list[np.ndarray] = []
    cols_list: list[np.ndarray] = []
    vals_list: list[np.ndarray] = []
    max_dof = -1

    for pts_c, w_c, n_c, d_c in zip(pts_seq, w_seq, n_seq, d_seq, strict=True):
        if pts_c.size == 0 or w_c.size == 0:
            continue
        pts_c = pts_c.reshape(-1, 3)
        w_c = w_c.reshape(-1)
        n_loc = n_c.shape[-1]
        n_c = n_c.reshape(-1, n_loc)
        d_c = d_c.reshape(-1)

        rho_c = _evaluate_density(rho, pts_c)
        eff_w = rho_c * w_c
        # Scalar cell mass matrix: M_scalar[a, b] = sum_q eff_w[q] * N_a(q) * N_b(q)
        m_scalar = (n_c * eff_w[:, None]).T @ n_c

        if d_c.size == 3 * n_loc:
            dofs_3d = d_c
        elif d_c.size == n_loc:
            dofs_3d = np.empty(3 * n_loc, dtype=np.int64)
            dofs_3d[0::3] = 3 * d_c
            dofs_3d[1::3] = 3 * d_c + 1
            dofs_3d[2::3] = 3 * d_c + 2
        else:
            raise ValueError(
                f"dof_indices length ({d_c.size}) must equal n_loc ({n_loc}) or 3*n_loc ({3 * n_loc})."
            )

        max_dof = max(max_dof, int(np.max(dofs_3d)))

        # Build 3-component block diagonal Kronecker structure in local DOFs
        rr, cc = np.nonzero(np.ones((n_loc, n_loc), dtype=bool))
        m_flat = m_scalar[rr, cc]
        for comp in range(3):
            rows_list.append(dofs_3d[3 * rr + comp])
            cols_list.append(dofs_3d[3 * cc + comp])
            vals_list.append(m_flat)

    total_dofs = max_dof + 1 if n_dofs is None else int(n_dofs)
    if not rows_list:
        return sp.csr_matrix((total_dofs, total_dofs), dtype=float)

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    return sp.coo_matrix((vals, (rows, cols)), shape=(total_dofs, total_dofs), dtype=float).tocsr()


def _axis_quadrature_tables(
    knots: np.ndarray,
    degree: int,
    bound: tuple[float, float],
    n_gp: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build 1D Gauss-Legendre basis tables across all non-empty knot spans.

    Args:
        knots: Clamped 1D knot vector.
        degree: B-spline polynomial degree ``p``.
        bound: Physical domain interval ``(lo, hi)`` along this axis.
        n_gp: Number of Gauss-Legendre points per knot span.

    Returns:
        Tuple ``(first_cp, basis, phys_pts, weights_detj)`` where:
            - ``first_cp`` has shape ``(n_cells,)``
            - ``basis`` has shape ``(n_cells, n_gp, degree + 1)``
            - ``phys_pts`` has shape ``(n_cells, n_gp)``
            - ``weights_detj`` has shape ``(n_cells, n_gp)``
    """
    knots_arr = np.asarray(knots, dtype=float)
    p = int(degree)
    n_cp = n_control_points(knots_arr, p)
    u_breaks = breakpoints(knots_arr)
    u_lo, u_hi = float(knots_arr[0]), float(knots_arr[-1])
    phys_lo, phys_hi = float(bound[0]), float(bound[1])
    if u_hi <= u_lo:
        raise ValueError("Knot vector must have positive parametric length.")
    if phys_hi <= phys_lo:
        raise ValueError("Physical domain bound must satisfy hi > lo.")
    scale = (phys_hi - phys_lo) / (u_hi - u_lo)

    xg, wg = np.polynomial.legendre.leggauss(n_gp)
    n_cells = len(u_breaks) - 1
    first_cp = np.empty(n_cells, dtype=np.int64)
    basis = np.empty((n_cells, n_gp, p + 1), dtype=float)
    phys_pts = np.empty((n_cells, n_gp), dtype=float)
    weights_detj = np.empty((n_cells, n_gp), dtype=float)

    for c in range(n_cells):
        ua, ub = float(u_breaks[c]), float(u_breaks[c + 1])
        jac_param = 0.5 * (ub - ua)
        u_mid = 0.5 * (ua + ub)
        span = find_span(n_cp, p, u_mid, knots_arr)
        first_cp[c] = span - p
        for g in range(n_gp):
            u = u_mid + jac_param * float(xg[g])
            basis[c, g, :] = ders_basis_funs(span, u, p, 0, knots_arr)[0]
            phys_pts[c, g] = phys_lo + (u - u_lo) * scale
            weights_detj[c, g] = float(wg[g]) * jac_param * scale

    return first_cp, basis, phys_pts, weights_detj


@njit(cache=True, parallel=True)
def _assemble_const_mass_njit(rho_val: float,
                              fu: np.ndarray, fv: np.ndarray, fw: np.ndarray,
                              mu_cells: np.ndarray, mv_cells: np.ndarray, mw_cells: np.ndarray,
                              nu_cp: int, nv_cp: int,
                              pu: int, pv: int, pw: int,
                              rows: np.ndarray, cols: np.ndarray, vals: np.ndarray) -> None:
    """Parallel COO triplet fill for constant-density tensor-product B-spline mass."""
    n_cu, n_cv, n_cw = fu.shape[0], fv.shape[0], fw.shape[0]
    du, dv, dw = pu + 1, pv + 1, pw + 1
    n_loc = du * dv * dw
    stride = 3 * n_loc * n_loc
    n_cells = n_cu * n_cv * n_cw

    for cid in prange(n_cells):
        kw = cid // (n_cu * n_cv)
        rem = cid - kw * (n_cu * n_cv)
        jv = rem // n_cu
        iu = rem - jv * n_cu

        u0, v0, w0 = fu[iu], fv[jv], fw[kw]
        base = cid * stride

        cps = np.empty(n_loc, dtype=np.int64)
        idx = 0
        for c in range(dw):
            wc = nu_cp * nv_cp * (w0 + c)
            for b in range(dv):
                vb = nu_cp * (v0 + b) + wc
                for a in range(du):
                    cps[idx] = (u0 + a) + vb
                    idx += 1

        nn = n_loc * n_loc
        for r_loc in range(n_loc):
            cr = r_loc // (du * dv)
            rem_r = r_loc - cr * (du * dv)
            br = rem_r // du
            ar = rem_r - br * du
            cp_r = 3 * cps[r_loc]

            for c_loc in range(n_loc):
                cc = c_loc // (du * dv)
                rem_c = c_loc - cc * (du * dv)
                bc = rem_c // du
                ac = rem_c - bc * du
                val = rho_val * mu_cells[iu, ar, ac] * mv_cells[jv, br, bc] * mw_cells[kw, cr, cc]
                cp_c = 3 * cps[c_loc]
                k_loc = r_loc * n_loc + c_loc
                for comp in range(3):
                    p_idx = base + comp * nn + k_loc
                    rows[p_idx] = cp_r + comp
                    cols[p_idx] = cp_c + comp
                    vals[p_idx] = val


def assemble_consistent_mass_from_bspline(
    knots_uvw: tuple[np.ndarray, np.ndarray, np.ndarray],
    degrees: tuple[int, int, int] | int,
    rho: DensityInput,
    bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None = None,
    n_gp: int | tuple[int, int, int] | None = None,
) -> sp.csr_matrix:
    """Assemble 3D consistent mass matrix on an axis-aligned B-spline domain.

    Uses exact tensor-product Gauss-Legendre quadrature over every knot span
    cell in ``(u, v, w)``. Control point linear indexing follows
    ``I = i_u + n_u * (i_v + n_v * i_w)`` with 3 DOFs ``(3*I, 3*I+1, 3*I+2)``.

    Args:
        knots_uvw: Tuple ``(knots_u, knots_v, knots_w)`` of clamped knot vectors.
        degrees: Polynomial degrees ``(p_u, p_v, p_w)`` or a single integer.
        rho: Constant scalar density or callable ``rho(points) -> density``.
        bounds: Physical box bounds ``((x0, x1), (y0, y1), (z0, z1))``. If
            ``None``, uses the knot vector endpoints ``(knots[0], knots[-1])``.
        n_gp: Number of Gauss-Legendre points per axis per knot span. Defaults
            to ``p_d + 1`` per axis (exact for degree-``2*p_d`` mass integrand).

    Returns:
        Sparse CSR consistent mass matrix of shape ``(3 * n_cp_total, 3 * n_cp_total)``.
    """
    deg_tuple = (int(degrees), int(degrees), int(degrees)) if isinstance(degrees, int) else tuple(int(d) for d in degrees)
    if len(deg_tuple) != 3 or len(knots_uvw) != 3:
        raise ValueError("knots_uvw and degrees must each have 3 components.")

    if bounds is None:
        bounds_tuple = tuple((float(k[0]), float(k[-1])) for k in knots_uvw)
    else:
        bounds_tuple = tuple((float(b[0]), float(b[1])) for b in bounds)

    if n_gp is None:
        ngp_tuple = tuple(p + 1 for p in deg_tuple)
    elif isinstance(n_gp, int):
        ngp_tuple = (int(n_gp), int(n_gp), int(n_gp))
    else:
        ngp_tuple = tuple(int(g) for g in n_gp)

    nu_cp = n_control_points(np.asarray(knots_uvw[0], dtype=float), deg_tuple[0])
    nv_cp = n_control_points(np.asarray(knots_uvw[1], dtype=float), deg_tuple[1])
    nw_cp = n_control_points(np.asarray(knots_uvw[2], dtype=float), deg_tuple[2])
    total_dofs = 3 * nu_cp * nv_cp * nw_cp

    fu, nu_tab, xu, wu = _axis_quadrature_tables(knots_uvw[0], deg_tuple[0], bounds_tuple[0], ngp_tuple[0])
    fv, nv_tab, xv, wv = _axis_quadrature_tables(knots_uvw[1], deg_tuple[1], bounds_tuple[1], ngp_tuple[1])
    fw, nw_tab, xw, ww = _axis_quadrature_tables(knots_uvw[2], deg_tuple[2], bounds_tuple[2], ngp_tuple[2])

    pu, pv, pw = deg_tuple
    n_loc = (pu + 1) * (pv + 1) * (pw + 1)
    n_cells = len(fu) * len(fv) * len(fw)

    # Fast path for constant scalar density using tensor-product separability
    if not callable(rho):
        rho_val = float(rho)
        rows = np.empty(n_cells * 3 * n_loc * n_loc, dtype=np.int64)
        cols = np.empty_like(rows)
        vals = np.empty_like(rows, dtype=float)

        # Precompute 1D cell mass matrices: (n_cells_d, p_d+1, p_d+1)
        mu_cells = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", nu_tab, wu, nu_tab), dtype=np.float64)
        mv_cells = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", nv_tab, wv, nv_tab), dtype=np.float64)
        mw_cells = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", nw_tab, ww, nw_tab), dtype=np.float64)

        _assemble_const_mass_njit(
            rho_val,
            np.ascontiguousarray(fu, dtype=np.int64),
            np.ascontiguousarray(fv, dtype=np.int64),
            np.ascontiguousarray(fw, dtype=np.int64),
            mu_cells, mv_cells, mw_cells,
            int(nu_cp), int(nv_cp), int(pu), int(pv), int(pw),
            rows, cols, vals,
        )
        return sp.coo_matrix((vals, (rows, cols)), shape=(total_dofs, total_dofs), dtype=float).tocsr()

    # General path for spatially varying density rho(x, y, z)
    cell_points_list: list[np.ndarray] = []
    weights_list: list[np.ndarray] = []
    basis_list: list[np.ndarray] = []
    dofs_list: list[np.ndarray] = []

    a_idx = np.arange(pu + 1)
    b_idx = np.arange(pv + 1)
    c_idx = np.arange(pw + 1)

    for kw in range(len(fw)):
        for jv in range(len(fv)):
            for iu in range(len(fu)):
                gx, gy, gz = np.meshgrid(xu[iu], xv[jv], xw[kw], indexing="ij")
                pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
                w3 = (wu[iu][:, None, None] * wv[jv][None, :, None] * ww[kw][None, None, :]).ravel()
                n3 = np.einsum(
                    "ga,hb,ic->ghicba",
                    nu_tab[iu],
                    nv_tab[jv],
                    nw_tab[kw],
                ).reshape(-1, n_loc)

                cp_u = fu[iu] + a_idx
                cp_v = fv[jv] + b_idx
                cp_w = fw[kw] + c_idx
                cp_grid = (
                    cp_u[None, None, :]
                    + nu_cp * (cp_v[None, :, None] + nv_cp * cp_w[:, None, None])
                ).reshape(-1)

                cell_points_list.append(pts)
                weights_list.append(w3)
                basis_list.append(n3)
                dofs_list.append(cp_grid)

    return assemble_consistent_mass(
        cell_points=cell_points_list,
        weights_detj=weights_list,
        basis_values=basis_list,
        dof_indices=dofs_list,
        rho=rho,
        n_dofs=total_dofs,
    )


def solve_eigenmodes(
    K: sp.spmatrix | np.ndarray,
    M: sp.spmatrix | np.ndarray,
    n_modes: int = 20,
    sigma: float = 0.0,
    free_dofs: np.ndarray | Sequence[int] | None = None,
) -> EigenmodeResult:
    """Solve the generalized structural eigenproblem ``K phi = omega^2 M phi``.

    Uses shift-invert Lanczos (``scipy.sparse.linalg.eigsh``) around ``sigma``
    when ``n_modes < n_free - 1``, falling back to dense ``scipy.linalg.eigh``
    for very small systems where ``n_modes >= n_free - 1``. Mode shapes are
    mass-orthonormalized so that ``Phi.T @ M @ Phi = I``.

    Args:
        K: Global stiffness matrix of shape ``(n_dofs, n_dofs)``.
        M: Global consistent mass matrix of shape ``(n_dofs, n_dofs)``.
        n_modes: Number of lowest eigenmodes to compute.
        sigma: Shift frequency-squared (rad/s)^2 for shift-invert ``eigsh``.
        free_dofs: Optional indices of unconstrained DOFs. If ``None``, all DOFs
            are treated as free.

    Returns:
        ``EigenmodeResult`` containing frequencies in Hz, circular frequencies,
        eigenvalues, mass-orthonormal modes, modal participation factors, and
        effective modal masses.
    """
    K_csr = sp.csr_matrix(K, dtype=float)
    M_csr = sp.csr_matrix(M, dtype=float)
    n_dofs = K_csr.shape[0]
    if K_csr.shape != (n_dofs, n_dofs) or M_csr.shape != (n_dofs, n_dofs):
        raise ValueError("K and M must be square matrices of identical shape.")

    if free_dofs is None:
        free = np.arange(n_dofs, dtype=np.int64)
    else:
        free = np.unique(np.asarray(free_dofs, dtype=np.int64))

    n_free = len(free)
    if n_free == 0:
        raise ValueError("free_dofs must contain at least one DOF.")
    if n_modes < 1:
        raise ValueError("n_modes must be >= 1.")

    k_req = min(int(n_modes), n_free)
    K_ff = K_csr[free, :][:, free]
    M_ff = M_csr[free, :][:, free]
    # Enforce exact symmetry against floating-point assembly rounding
    K_ff = 0.5 * (K_ff + K_ff.T)
    M_ff = 0.5 * (M_ff + M_ff.T)

    if k_req >= n_free - 1:
        eigvals, eigvecs = la.eigh(K_ff.toarray(), M_ff.toarray())
        eigvals = eigvals[:k_req]
        eigvecs = eigvecs[:, :k_req]
    else:
        eigvals, eigvecs = spla.eigsh(K_ff, k=k_req, M=M_ff, sigma=float(sigma), which="LM")

    order = np.argsort(eigvals)
    eigvals = np.asarray(eigvals[order], dtype=float)
    eigvecs = np.asarray(eigvecs[:, order], dtype=float)

    # Clamp tiny numerical round-off negative zero-modes (e.g. -1e-13)
    eigvals = np.where((eigvals < 0.0) & (eigvals > -1e-9), 0.0, eigvals)
    if np.any(eigvals < 0.0):
        raise ValueError(f"Negative eigenvalue encountered: min={np.min(eigvals):.3e}")

    # Expand to full DOF vector and enforce exact M-orthonormality
    modes = np.zeros((n_dofs, k_req), dtype=float)
    modes[free, :] = eigvecs

    for j in range(k_req):
        modal_mass = float(modes[:, j] @ (M_csr @ modes[:, j]))
        if modal_mass <= 0.0:
            raise ValueError(f"Mode {j} has non-positive modal mass ({modal_mass:.3e}).")
        modes[:, j] /= np.sqrt(modal_mass)

    omega = np.sqrt(eigvals)
    frequencies_hz = omega / (2.0 * np.pi)

    # Modal participation factors for rigid-body unit translations in x, y, z
    participation = np.zeros((k_req, 3), dtype=float)
    total_mass = np.zeros(3, dtype=float)
    for d in range(3):
        r_d = np.zeros(n_dofs, dtype=float)
        r_d[d::3] = 1.0
        mr_d = M_csr @ r_d
        total_mass[d] = float(r_d @ mr_d)
        participation[:, d] = modes.T @ mr_d

    effective_modal_masses = participation**2

    return EigenmodeResult(
        frequencies_hz=frequencies_hz,
        omega=omega,
        eigenvalues=eigvals,
        modes=modes,
        participation_factors=participation,
        effective_modal_masses=effective_modal_masses,
        total_mass=total_mass,
        free_dofs=free,
    )


def align_degenerate_modes(
    eig: EigenmodeResult,
    direction: np.ndarray | Sequence[float],
    rtol: float = 1e-6,
) -> EigenmodeResult:
    """Rotate every cluster of degenerate modes so one mode carries the whole participation.

    An eigensolver returns an arbitrary orthonormal basis of a degenerate eigenspace
    (a square cross-section gives pairs of equal bending frequencies).  Sum rules that
    drop the cross terms between modes (``ABS``, ``SRSS``, ``NRL``, ...) are not
    invariant under such a rotation, so their result would depend on the solver's random
    starting vector.  Rotating inside each cluster so that the participation factor in
    ``direction`` ends up in the first mode of the cluster (and is zero in the others)
    gives a unique, excitation-aligned basis.  The rotation is orthogonal, so modes stay
    mass-orthonormal and every rotation-invariant result (CQC with equal damping, random
    response, effective modal mass per cluster) is unchanged.

    Args:
        eig: Eigenmodes.
        direction: Excitation direction (need not be normalised), shape ``(3,)``.
        rtol: Modes whose neighbouring frequencies differ by at most ``rtol * f`` form a cluster.

    Returns:
        New :class:`EigenmodeResult` with rotated ``modes``, ``participation_factors`` and
        ``effective_modal_masses``; all other fields are copied.
    """
    d = np.asarray(direction, dtype=float).reshape(3)
    norm = float(np.linalg.norm(d))
    if norm == 0.0:
        raise ValueError("direction must be non-zero")
    d = d / norm
    f = np.asarray(eig.frequencies_hz, dtype=float)
    modes = np.array(eig.modes, dtype=float, copy=True)
    gam = np.array(eig.participation_factors, dtype=float, copy=True)
    k = 0
    while k < f.size:
        end = k + 1
        while end < f.size and abs(f[end] - f[end - 1]) <= rtol * max(abs(f[end - 1]), 1e-300):
            end += 1
        if end - k > 1:
            g = gam[k:end] @ d
            gn = float(np.linalg.norm(g))
            if gn > 0.0:
                m = end - k
                basis = np.eye(m)
                basis[:, 0] = g / gn
                q, r = np.linalg.qr(basis)
                q = q * np.sign(r[0, 0])  # first column must be +g/|g|
                modes[:, k:end] = modes[:, k:end] @ q
                gam[k:end] = q.T @ gam[k:end]
        k = end
    return EigenmodeResult(
        frequencies_hz=eig.frequencies_hz,
        omega=eig.omega,
        eigenvalues=eig.eigenvalues,
        modes=modes,
        participation_factors=gam,
        effective_modal_masses=gam**2,
        total_mass=eig.total_mass,
        free_dofs=eig.free_dofs,
    )

