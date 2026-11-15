"""Heat transfer and thermal-strain analysis on box B-splines (P4-C, V11/V12).

What is done where
------------------
* The conduction matrix ``K_th`` comes from Kratos ``LaplacianElement`` (see
  ``wht_iga.kratos_backend.thermal.kratos_box_conductivity``); :func:`assemble_conduction`
  is an independent python assembly used as its oracle.
* Kratos' ``LaplacianElement::CalculateMassMatrix`` carries ``k N^T N`` and no ``rho c``
  (spike S11), and ``SupportLaplacianCondition`` / quadrature-point ``ThermalFace`` are
  unusable on immersed geometry.  The capacity ``C_th = int rho c N^T N dV``, the film
  and radiation face terms and the time integration are therefore assembled here.
* Thermal strain of an immersed hot spot (V11): both regions share ``E`` and ``nu``, so the
  stiffness is that of the uniform plate and the eigenstrain enters only as the load
  ``f = int B^T D eps* dV`` over the hot-spot quadrature (:func:`eigenstrain_load`).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.grid.bspline import basis_1d, ders_basis_funs
from wht_iga.reference.assemble_elastic import isotropic_D
from wht_iga.reference.box import BoxSpline

STEFAN_BOLTZMANN = 5.67e-8  # W / (m^2 K^4), value of the NAFEMS T2 data and of Kratos' ThermalFace


@dataclass(frozen=True)
class FaceRule:
    """Quadrature of one box face.

    Attributes:
        basis: Sparse ``(n_q, n_ctrl)`` matrix of basis values at the face points.
        weights: Quadrature weights (area), shape ``(n_q,)``.
        points: Physical coordinates ``(n_q, 3)``.
    """

    basis: sp.csr_matrix
    weights: np.ndarray
    points: np.ndarray


def _cp_index(box: BoxSpline, i: int, j: int, k: int) -> int:
    nu, nv, _nw = box.n_cp
    return i + nu * (j + nv * k)


@njit(cache=True, parallel=True)
def _assemble_conduction_njit(cond: float,
                              fu: np.ndarray, fv: np.ndarray, fw: np.ndarray,
                              mx: np.ndarray, gx_m: np.ndarray,
                              my: np.ndarray, gy_m: np.ndarray,
                              mz: np.ndarray, gz_m: np.ndarray,
                              nu_cp: int, nv_cp: int, d: int,
                              rows: np.ndarray, cols: np.ndarray, vals: np.ndarray) -> None:
    """Parallel COO triplet fill for isotropic 3D scalar conduction K_th."""
    nx, ny, nz = fu.shape[0], fv.shape[0], fw.shape[0]
    n_loc = d * d * d
    stride = n_loc * n_loc
    n_cells = nx * ny * nz

    for cid in prange(n_cells):
        k = cid // (nx * ny)
        rem = cid - k * (nx * ny)
        j = rem // nx
        i = rem - j * nx

        u0, v0, w0 = fu[i], fv[j], fw[k]
        base = cid * stride

        cps = np.empty(n_loc, dtype=np.int64)
        idx = 0
        for c in range(d):
            wc = nu_cp * nv_cp * (w0 + c)
            for b in range(d):
                vb = nu_cp * (v0 + b) + wc
                for a in range(d):
                    cps[idx] = (u0 + a) + vb
                    idx += 1

        for r_loc in range(n_loc):
            cr = r_loc // (d * d)
            rem_r = r_loc - cr * (d * d)
            br = rem_r // d
            ar = rem_r - br * d
            cp_r = cps[r_loc]

            for c_loc in range(n_loc):
                cc = c_loc // (d * d)
                rem_c = c_loc - cc * (d * d)
                bc = rem_c // d
                ac = rem_c - bc * d
                val = cond * (
                    gx_m[i, ar, ac] * my[j, br, bc] * mz[k, cr, cc]
                    + mx[i, ar, ac] * gy_m[j, br, bc] * mz[k, cr, cc]
                    + mx[i, ar, ac] * my[j, br, bc] * gz_m[k, cr, cc]
                )
                p_idx = base + r_loc * n_loc + c_loc
                rows[p_idx] = cp_r
                cols[p_idx] = cps[c_loc]
                vals[p_idx] = val


def assemble_conduction(box: BoxSpline, conductivity: float, n_gauss: int | None = None) -> sp.csr_matrix:
    """Independent assembly of ``K_th = int k grad(N)^T grad(N) dV`` (one DOF per control point).

    Args:
        box: Box B-spline.
        conductivity: Isotropic conductivity ``k``.
        n_gauss: Gauss points per axis and cell (default ``p + 1``).

    Returns:
        CSR matrix ``(n_ctrl, n_ctrl)``.
    """
    p = int(box.degree)
    d = p + 1
    ng = d if n_gauss is None else int(n_gauss)
    gx, gw = np.polynomial.legendre.leggauss(ng)
    (fu, Nu, dNu, _, ju), (fv, Nv, dNv, _, jv), (fw, Nw, dNw, _, jw) = [box.axis_tables(ax, gx) for ax in range(3)]
    wu = gw[None, :] * ju[:, None]
    wv = gw[None, :] * jv[:, None]
    ww = gw[None, :] * jw[:, None]
    mx = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", Nu, wu, Nu), dtype=np.float64)
    gx_m = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", dNu, wu, dNu), dtype=np.float64)
    my = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", Nv, wv, Nv), dtype=np.float64)
    gy_m = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", dNv, wv, dNv), dtype=np.float64)
    mz = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", Nw, ww, Nw), dtype=np.float64)
    gz_m = np.ascontiguousarray(np.einsum("cga,cg,cgb->cab", dNw, ww, dNw), dtype=np.float64)

    nx, ny, nz = (int(v) for v in box.n_cells)
    n_loc = d * d * d
    nnz_trip = nx * ny * nz * n_loc * n_loc
    rows = np.empty(nnz_trip, dtype=np.int64)
    cols = np.empty(nnz_trip, dtype=np.int64)
    vals = np.empty(nnz_trip, dtype=np.float64)
    _assemble_conduction_njit(
        float(conductivity),
        np.ascontiguousarray(fu, dtype=np.int64),
        np.ascontiguousarray(fv, dtype=np.int64),
        np.ascontiguousarray(fw, dtype=np.int64),
        mx, gx_m, my, gy_m, mz, gz_m,
        int(box.n_cp[0]), int(box.n_cp[1]), d,
        rows, cols, vals,
    )
    n = box.n_ctrl
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def assemble_capacity(box: BoxSpline, rho_c: float) -> sp.csr_matrix:
    """Consistent capacity ``C_th = int rho c N^T N dV`` (scalar field, one DOF per control point).

    Args:
        box: Box B-spline.
        rho_c: Volumetric heat capacity ``rho * c``.

    Returns:
        CSR matrix ``(n_ctrl, n_ctrl)``.
    """
    from wht_iga.dynamics import assemble_consistent_mass_from_bspline

    m3 = assemble_consistent_mass_from_bspline(box.knots, degrees=box.degree, rho=float(rho_c))
    return sp.csr_matrix(m3[0::3, 0::3])


def face_rule(box: BoxSpline, axis: int, side: int, n_gauss: int | None = None) -> FaceRule:
    """Gauss rule of the box face ``axis``/``side`` with the basis values at its points.

    Args:
        box: Box B-spline.
        axis: Face normal axis ``0, 1, 2``.
        side: ``0`` low face, ``1`` high face.
        n_gauss: Gauss points per tangential axis and cell (default ``p + 1``).

    Returns:
        :class:`FaceRule`.
    """
    p = box.degree
    ng = p + 1 if n_gauss is None else int(n_gauss)
    gx, gw = np.polynomial.legendre.leggauss(ng)
    tangential = [d for d in range(3) if d != axis]
    t1, t2 = tangential
    x_n = float(box.hi[axis] if side else box.lo[axis])
    span, ders = basis_1d(box.knots[axis], p, x_n, 0)
    first_n, Nn = int(span) - p, np.asarray(ders[0])
    f1, N1, _d1, X1, j1 = box.axis_tables(t1, gx)
    f2, N2, _d2, X2, j2 = box.axis_tables(t2, gx)
    a = np.arange(p + 1)
    rows, cols, vals, w_all, pts = [], [], [], [], []
    q = 0
    for c2 in range(len(f2)):
        for c1 in range(len(f1)):
            for g2 in range(ng):
                for g1 in range(ng):
                    idx = {axis: first_n + a, t1: f1[c1] + a, t2: f2[c2] + a}
                    val = {axis: Nn, t1: N1[c1, g1], t2: N2[c2, g2]}
                    cps = (idx[0][None, None, :] + box.n_cp[0] * (idx[1][None, :, None] + box.n_cp[1] * idx[2][:, None, None])).reshape(-1)
                    bv = (val[0][None, None, :] * val[1][None, :, None] * val[2][:, None, None]).reshape(-1)
                    rows.append(np.full(cps.size, q))
                    cols.append(cps)
                    vals.append(bv)
                    w_all.append(gw[g1] * j1[c1] * gw[g2] * j2[c2])
                    pt = np.empty(3)
                    pt[axis], pt[t1], pt[t2] = x_n, X1[c1, g1], X2[c2, g2]
                    pts.append(pt)
                    q += 1
    basis = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(q, box.n_ctrl)).tocsr()
    return FaceRule(basis=basis, weights=np.asarray(w_all), points=np.asarray(pts))


def solve_with_dirichlet(A: sp.spmatrix, b: np.ndarray, fixed: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Solve ``A x = b`` with ``x[fixed] = values`` by elimination.

    Args:
        A: Square matrix.
        b: Right-hand side.
        fixed: Indices of prescribed DOFs.
        values: Prescribed values (same length as ``fixed``).

    Returns:
        Full solution vector.
    """
    n = A.shape[0]
    x = np.zeros(n)
    fixed = np.asarray(fixed, dtype=np.int64)
    x[fixed] = values
    free = np.setdiff1d(np.arange(n), fixed)
    A = sp.csr_matrix(A)
    x[free] = spla.spsolve(sp.csc_matrix(A[free][:, free]), b[free] - A[free][:, fixed] @ x[fixed])
    return x


def steady_thermal(
    K: sp.spmatrix,
    n_ctrl: int,
    films: Sequence[tuple[FaceRule, float, float]] = (),
    radiations: Sequence[tuple[FaceRule, float, float]] = (),
    fixed: np.ndarray = np.zeros(0, dtype=np.int64),
    fixed_values: np.ndarray = np.zeros(0),
    source: np.ndarray | None = None,
    t_guess: np.ndarray | None = None,
    tol: float = 1e-10,
    max_iter: int = 30,
) -> tuple[np.ndarray, int]:
    """Steady conduction with film (``h``) and radiation (``eps``) face terms by Newton iteration.

    Residual ``R(T) = K T + sum h int N^T (T - T_inf) dA + sum eps sigma int N^T (T^4 - T_inf^4) dA - f``.
    Temperatures are absolute (Kelvin) whenever a radiation term is present.

    Args:
        K: Conduction matrix ``(n_ctrl, n_ctrl)``.
        n_ctrl: Number of control points.
        films: ``(face rule, h, T_inf)`` triples.
        radiations: ``(face rule, emissivity, T_inf)`` triples.
        fixed: Prescribed control points.
        fixed_values: Prescribed temperatures.
        source: Optional nodal heat source vector.
        t_guess: Initial iterate (default: mean of the prescribed values).
        tol: Relative update tolerance.
        max_iter: Newton iteration cap.

    Returns:
        ``(temperatures, iterations)``.
    """
    fixed = np.asarray(fixed, dtype=np.int64)
    fixed_values = np.asarray(fixed_values, dtype=float)
    f = np.zeros(n_ctrl) if source is None else np.asarray(source, dtype=float)
    t = np.full(n_ctrl, float(fixed_values.mean()) if fixed_values.size else 0.0) if t_guess is None else np.array(t_guess, dtype=float)
    t[fixed] = fixed_values
    free = np.setdiff1d(np.arange(n_ctrl), fixed)
    for it in range(1, max_iter + 1):
        res = K @ t - f
        jac = sp.csr_matrix(K, copy=True)
        for rule, h, t_inf in films:
            tq = rule.basis @ t
            res += rule.basis.T @ (rule.weights * h * (tq - t_inf))
            jac = jac + h * (rule.basis.T @ sp.diags(rule.weights) @ rule.basis)
        for rule, eps, t_inf in radiations:
            tq = rule.basis @ t
            res += rule.basis.T @ (rule.weights * eps * STEFAN_BOLTZMANN * (tq**4 - t_inf**4))
            jac = jac + rule.basis.T @ sp.diags(rule.weights * 4.0 * eps * STEFAN_BOLTZMANN * tq**3) @ rule.basis
        dt = np.zeros(n_ctrl)
        dt[free] = spla.spsolve(sp.csc_matrix(jac[free][:, free]), -res[free])
        t += dt
        if np.abs(dt[free]).max() <= tol * max(1.0, np.abs(t[free]).max()):
            return t, it
    raise RuntimeError(f"steady thermal Newton did not converge in {max_iter} iterations")


def transient_thermal(
    C: sp.spmatrix,
    K: sp.spmatrix,
    t0: np.ndarray,
    dt: float,
    n_steps: int,
    fixed: np.ndarray,
    fixed_value_fn: Callable[[float], np.ndarray],
    scheme: str = "CN",
    source_fn: Callable[[float], np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate ``C dT/dt + K T = f(t)`` with time-dependent Dirichlet values.

    ``scheme``: ``"BE"`` (backward Euler), ``"CN"`` (Crank-Nicolson, theta = 1/2) or
    ``"BDF2"`` (the first step is Crank-Nicolson).  Prescribed temperatures enter the
    discrete equations exactly (no differentiation of the boundary data).

    Args:
        C: Capacity matrix.
        K: Conduction matrix.
        t0: Initial temperatures, consistent with the Dirichlet values at ``t = 0``.
        dt: Constant time step.
        n_steps: Number of steps.
        fixed: Prescribed control points.
        fixed_value_fn: ``t -> values`` of the prescribed temperatures.
        scheme: ``"BE"``, ``"CN"`` or ``"BDF2"``.
        source_fn: Optional ``t -> nodal source vector``.

    Returns:
        ``(times (n_steps + 1,), temperatures (n_steps + 1, n_ctrl))``.
    """
    scheme = scheme.upper()
    if scheme not in ("BE", "CN", "BDF2"):
        raise ValueError(f"unknown scheme {scheme!r}")
    C, K = sp.csr_matrix(C), sp.csr_matrix(K)
    n = K.shape[0]
    fixed = np.asarray(fixed, dtype=np.int64)
    free = np.setdiff1d(np.arange(n), fixed)
    src = (lambda t: np.zeros(n)) if source_fn is None else source_fn  # noqa: E731
    out = np.zeros((n_steps + 1, n))
    out[0] = t0
    out[0, fixed] = fixed_value_fn(0.0)
    times = np.arange(n_steps + 1) * dt

    def lu_of(theta: float | None, bdf: bool):
        a = (1.5 / dt) * C + K if bdf else C / dt + float(theta) * K
        a = sp.csr_matrix(a)
        return a, spla.splu(sp.csc_matrix(a[free][:, free]))

    a_one, lu_one = lu_of(0.5 if scheme in ("CN", "BDF2") else 1.0, False)
    a_bdf, lu_bdf = lu_of(None, True) if scheme == "BDF2" else (None, None)
    for step in range(1, n_steps + 1):
        t_new, t_old = times[step], times[step - 1]
        prev = out[step - 1]
        new_fixed = np.asarray(fixed_value_fn(t_new), dtype=float)
        if scheme == "BDF2" and step >= 2:
            rhs = (2.0 / dt) * (C @ prev) - (0.5 / dt) * (C @ out[step - 2]) + src(t_new)
            a, lu = a_bdf, lu_bdf
        else:
            theta = 1.0 if scheme == "BE" else 0.5
            rhs = (C / dt - (1.0 - theta) * K) @ prev + theta * src(t_new) + (1.0 - theta) * src(t_old)
            a, lu = a_one, lu_one
        x = np.zeros(n)
        x[fixed] = new_fixed
        x[free] = lu.solve(rhs[free] - a[free][:, fixed] @ new_fixed)
        out[step] = x
    return times, out


def eigenstrain_load(
    box: BoxSpline,
    D: np.ndarray,
    eigenstrain: np.ndarray,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray],
    full_cells: Sequence[tuple[int, int, int]] = (),
) -> np.ndarray:
    """Load ``f = int_region B^T D eps* dV`` of a region with prescribed eigenstrain.

    The region is the union of ``full_cells`` (Gauss rule ``p + 1``) and the cut cells of
    ``cell_tables`` (moment-fitted points ``(x, y, z, w)``).

    Args:
        box: Box B-spline (3 displacement DOFs per control point).
        D: ``(6, 6)`` elasticity matrix, Voigt ``[xx, yy, zz, xy, yz, xz]`` with engineering shears.
        eigenstrain: Eigenstrain ``(6,)`` in the same Voigt order (engineering shears).
        cell_tables: ``{cell: (nq, 4)}`` for cut cells.
        full_cells: Cells entirely inside the region.

    Returns:
        Force vector ``(box.n_dofs,)`` in reference DOF order.
    """
    p = box.degree
    ng = p + 1
    gx, gw = np.polynomial.legendre.leggauss(ng)
    br = [box.breakpoints(d) for d in range(3)]
    sigma_star = np.asarray(D, dtype=float) @ np.asarray(eigenstrain, dtype=float)
    f = np.zeros(box.n_dofs)
    a = np.arange(p + 1)

    def add(cell: tuple[int, int, int], table: np.ndarray) -> None:
        spans = []
        for d in range(3):
            mid = 0.5 * (br[d][cell[d]] + br[d][cell[d] + 1])
            span, _ = basis_1d(box.knots[d], p, mid, 0)
            spans.append(int(span))
        cps = ((spans[0] - p + a)[None, None, :] + box.n_cp[0] * ((spans[1] - p + a)[None, :, None] + box.n_cp[1] * (spans[2] - p + a)[:, None, None])).reshape(-1)
        for x, y, z, w in table:
            grads = []
            vals = []
            for d, xd in enumerate((x, y, z)):
                dd = ders_basis_funs(spans[d], float(xd), p, 1, np.ascontiguousarray(box.knots[d], dtype=float))
                vals.append(dd[0])
                grads.append(dd[1])
            n_x, n_y, n_z = vals
            g = [
                np.einsum("a,b,c->cba", grads[0], n_y, n_z).reshape(-1),
                np.einsum("a,b,c->cba", n_x, grads[1], n_z).reshape(-1),
                np.einsum("a,b,c->cba", n_x, n_y, grads[2]).reshape(-1),
            ]
            for local, cp in enumerate(cps):
                gx_, gy_, gz_ = g[0][local], g[1][local], g[2][local]
                # B^T sigma* with Voigt [xx, yy, zz, xy, yz, xz], engineering shear strains
                f[3 * cp + 0] += w * (gx_ * sigma_star[0] + gy_ * sigma_star[3] + gz_ * sigma_star[5])
                f[3 * cp + 1] += w * (gy_ * sigma_star[1] + gx_ * sigma_star[3] + gz_ * sigma_star[4])
                f[3 * cp + 2] += w * (gz_ * sigma_star[2] + gy_ * sigma_star[4] + gx_ * sigma_star[5])

    for cell in full_cells:
        pts = []
        for d in range(3):
            lo, hi = br[d][cell[d]], br[d][cell[d] + 1]
            pts.append((lo + 0.5 * (hi - lo) * (gx + 1.0), 0.5 * (hi - lo) * gw))
        table = np.array(
            [(px, py, pz, wx * wy * wz) for (px, wx) in zip(*pts[0]) for (py, wy) in zip(*pts[1]) for (pz, wz) in zip(*pts[2])]
        )
        add(tuple(cell), table)
    for cell, table in cell_tables.items():
        add(tuple(cell), np.asarray(table, dtype=float))
    return f


def plane_stress_as_plane_strain(E: float, nu: float, eigenstrain_in_plane: float) -> tuple[float, float, float]:
    """3D plane-strain (``u_z = 0``) properties that reproduce a plane-stress problem in-plane.

    A plane-strain solid ``(E_p, nu_p)`` with isotropic 3D eigenstrain ``e_p`` has in-plane
    behaviour ``eps = (1 + nu_p) e_p`` + the plane-stress law of ``E_s = E_p / (1 - nu_p**2)``,
    ``nu_s = nu_p / (1 - nu_p)``.  Inverting for given ``(E_s, nu_s, e_s)``:
    ``nu_p = nu_s / (1 + nu_s)``, ``E_p = E_s (1 - nu_p**2)``, ``e_p = e_s / (1 + nu_p)``.

    Args:
        E: Plane-stress Young's modulus.
        nu: Plane-stress Poisson's ratio.
        eigenstrain_in_plane: Plane-stress isotropic in-plane eigenstrain.

    Returns:
        ``(E_p, nu_p, e_p)`` for the 3D plane-strain solid.
    """
    nu_p = nu / (1.0 + nu)
    return E * (1.0 - nu_p**2), nu_p, eigenstrain_in_plane / (1.0 + nu_p)


def quarter_disc_prism_soup(
    radius: float,
    z_lo: float,
    z_hi: float,
    n_seg: int = 360,
    area_preserving: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Closed outward-wound triangle soup of the quarter cylinder ``x, y >= 0``, ``r <= radius``.

    The solid lies on the high side of the grid's lower planes ``x = 0`` and ``y = 0`` and its
    walls on those planes are coplanar with the grid faces, exactly like the V06 octant.  This
    is required by ``wht_iga.geometry.classify``: its ray-parity scan starts every ray at
    ``grid minimum - 1e-6`` and assumes that start point is *outside* the solid (a cylinder
    centred on the origin made every interior vertex come out OUT).  Choose ``z_lo`` at the
    grid's lower plane and ``z_hi`` above the grid's upper plane.

    The polygon radius is enlarged so that the polygon has the area of the true sector.  The
    cap triangulation fans from an interior point, which keeps cap edges off the grid lines.

    Args:
        radius: Circle radius.
        z_lo: Lower cap height (the grid's lower plane).
        z_hi: Upper cap height (above the grid's upper plane).
        n_seg: Number of arc segments.
        area_preserving: Scale the arc radius to preserve the sector's area.

    Returns:
        ``(nodes (n, 3), triangles (m, 3))`` with outward windings.
    """
    n = int(n_seg)
    ang = 0.5 * np.pi / n
    r_poly = radius * (ang / np.sin(ang)) ** 0.5 if area_preserving else radius
    th = ang * np.arange(n + 1)
    poly = np.vstack([[0.0, 0.0], np.column_stack([r_poly * np.cos(th), r_poly * np.sin(th)])])  # CCW
    m = poly.shape[0]
    fan = (0.2137 * radius, 0.1789 * radius)  # interior point, off every grid line
    nodes = np.vstack([
        np.column_stack([poly, np.full(m, z_lo)]),
        np.column_stack([poly, np.full(m, z_hi)]),
        [[fan[0], fan[1], z_lo], [fan[0], fan[1], z_hi]],
    ])
    bot_c, top_c = 2 * m, 2 * m + 1
    tris = []
    for i in range(m):
        j = (i + 1) % m
        tris.append((i, j, m + j))          # wall, outward for a CCW polygon
        tris.append((i, m + j, m + i))
        tris.append((top_c, m + i, m + j))  # top cap, +z
        tris.append((bot_c, j, i))          # bottom cap, -z
    return nodes, np.asarray(tris, dtype=np.int64)


@dataclass(frozen=True)
class ThermoelasticTdepResult:
    """Solution of a spatially varying temperature-dependent thermo-elastic problem.

    Attributes:
        u: Displacement vector in reference DOF order, shape ``(box.n_dofs,)``.
        K: Assembled tangent/secant stiffness matrix ``int B^T D(T(x)) B dV``.
        f: Total right-hand-side load vector ``(box.n_dofs,)``.
    """

    u: np.ndarray
    K: sp.csr_matrix
    f: np.ndarray


def _resolve_scalar_material(fn_or_val: float | Callable[[float], float]) -> Callable[[float], float]:
    """Return a scalar ``T -> value`` evaluator from a constant or callable."""
    if callable(fn_or_val):
        return fn_or_val
    val = float(fn_or_val)
    return lambda _t, v=val: v


def _cell_span_and_cps(box: BoxSpline, cell: tuple[int, int, int]) -> tuple[list[int], np.ndarray]:
    """Return ``(spans, global_cp_indices)`` for one grid cell."""
    p = box.degree
    br = [box.breakpoints(d) for d in range(3)]
    a = np.arange(p + 1)
    spans: list[int] = []
    for d in range(3):
        mid = 0.5 * (br[d][cell[d]] + br[d][cell[d] + 1])
        span, _ = basis_1d(box.knots[d], p, mid, 0)
        spans.append(int(span))
    cps = (
        (spans[0] - p + a)[None, None, :]
        + box.n_cp[0] * ((spans[1] - p + a)[None, :, None] + box.n_cp[1] * (spans[2] - p + a)[:, None, None])
    ).reshape(-1)
    return spans, cps


def _eval_basis_and_B_batch(
    box: BoxSpline, spans: list[int], pts: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate trivariate B-spline basis ``N (nq, n_loc)`` and ``B (nq, 6, 3*n_loc)``."""
    p = box.degree
    nq = pts.shape[0]
    axis_vals = []
    axis_grads = []
    for d in range(3):
        knots_d = np.ascontiguousarray(box.knots[d], dtype=float)
        uniq, inv = np.unique(pts[:, d], return_inverse=True)
        u_vals = np.empty((uniq.size, p + 1), dtype=float)
        u_grads = np.empty((uniq.size, p + 1), dtype=float)
        for u_idx, xd in enumerate(uniq):
            dd = ders_basis_funs(spans[d], float(xd), p, 1, knots_d)
            u_vals[u_idx] = dd[0]
            u_grads[u_idx] = dd[1]
        axis_vals.append(u_vals[inv])
        axis_grads.append(u_grads[inv])

    nx, ny, nz = axis_vals
    N = np.einsum("qa,qb,qc->qcba", nx, ny, nz).reshape(nq, -1)
    gx = np.einsum("qa,qb,qc->qcba", axis_grads[0], ny, nz).reshape(nq, -1)
    gy = np.einsum("qa,qb,qc->qcba", nx, axis_grads[1], nz).reshape(nq, -1)
    gz = np.einsum("qa,qb,qc->qcba", nx, ny, axis_grads[2]).reshape(nq, -1)
    n_loc = N.shape[1]
    B = np.zeros((nq, 6, 3 * n_loc), dtype=float)
    c0 = 3 * np.arange(n_loc)
    B[:, 0, c0] = gx
    B[:, 1, c0 + 1] = gy
    B[:, 2, c0 + 2] = gz
    B[:, 3, c0] = gy
    B[:, 3, c0 + 1] = gx
    B[:, 4, c0 + 1] = gz
    B[:, 4, c0 + 2] = gy
    B[:, 5, c0] = gz
    B[:, 5, c0 + 2] = gx
    return N, B


def _isotropic_D_batch(E: np.ndarray, nu: np.ndarray) -> np.ndarray:
    """Vectorized ``(nq, 6, 6)`` isotropic elasticity matrices in Voigt ``[xx, yy, zz, xy, yz, xz]``."""
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    D = np.zeros((E.size, 6, 6), dtype=float)
    D[:, :3, :3] = lam[:, None, None]
    D[:, 0, 0] = D[:, 1, 1] = D[:, 2, 2] = lam + 2.0 * mu
    D[:, 3, 3] = D[:, 4, 4] = D[:, 5, 5] = mu
    return D


def _iter_box_quadrature(
    box: BoxSpline,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray] | None,
    full_cells: Sequence[tuple[int, int, int]] | None,
    n_gauss: int | None,
):
    """Yield ``(cell, table (nq, 4))`` across full cells and optional cut-cell tables."""
    p = box.degree
    ng = p + 1 if n_gauss is None else int(n_gauss)
    gx, gw = np.polynomial.legendre.leggauss(ng)
    br = [box.breakpoints(d) for d in range(3)]
    if cell_tables is None and full_cells is None:
        nx, ny, nz = box.n_cells
        full_iter: Sequence[tuple[int, int, int]] = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    else:
        full_iter = full_cells or ()
    for cell in full_iter:
        pts = []
        for d in range(3):
            lo, hi = br[d][cell[d]], br[d][cell[d] + 1]
            pts.append((lo + 0.5 * (hi - lo) * (gx + 1.0), 0.5 * (hi - lo) * gw))
        table = np.array(
            [(px, py, pz, wx * wy * wz) for (px, wx) in zip(*pts[0]) for (py, wy) in zip(*pts[1]) for (pz, wz) in zip(*pts[2])],
            dtype=float,
        )
        yield tuple(cell), table
    if cell_tables is not None:
        for cell, table in cell_tables.items():
            yield tuple(cell), np.asarray(table, dtype=float)


def _eval_temperature_fn(
    box: BoxSpline,
    temperature: np.ndarray | Callable[[float, float, float], float],
) -> Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    """Build batch point-temperature evaluator ``(N, cps, pts) -> T_q (nq,)``."""
    if callable(temperature):
        def at_quad(_N: np.ndarray, _cps: np.ndarray, pts: np.ndarray) -> np.ndarray:
            return np.array([float(temperature(float(p[0]), float(p[1]), float(p[2]))) for p in pts], dtype=float)

        return at_quad

    t_cp = np.asarray(temperature, dtype=float).reshape(box.n_ctrl)

    def at_quad(N: np.ndarray, cps: np.ndarray, _pts: np.ndarray) -> np.ndarray:
        return N @ t_cp[cps]

    return at_quad


def _assemble_tdep_system(
    box: BoxSpline,
    temperature: np.ndarray | Callable[[float, float, float], float],
    E_fn: float | Callable[[float], float],
    nu_fn: float | Callable[[float], float] = 0.0,
    alpha_fn: float | Callable[[float], float] = 0.0,
    T_ref: float = 0.0,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray] | None = None,
    full_cells: Sequence[tuple[int, int, int]] | None = None,
    n_gauss: int | None = None,
    assemble_K: bool = True,
    assemble_f: bool = True,
) -> tuple[sp.csr_matrix | None, np.ndarray | None]:
    """Assemble spatially varying temperature-dependent stiffness ``K`` and/or thermal load ``f``."""
    eval_E = _resolve_scalar_material(E_fn)
    eval_nu = _resolve_scalar_material(nu_fn)
    eval_alpha = _resolve_scalar_material(alpha_fn)
    eval_T_quad = _eval_temperature_fn(box, temperature)

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    f = np.zeros(box.n_dofs, dtype=float) if assemble_f else None

    for cell, table in _iter_box_quadrature(box, cell_tables, full_cells, n_gauss):
        if table.shape[0] == 0:
            continue
        spans, cps = _cell_span_and_cps(box, cell)
        dofs = (3 * cps[:, None] + np.arange(3)[None, :]).reshape(-1)
        n_edof = dofs.size
        pts = table[:, :3]
        wts = table[:, 3]
        N, B = _eval_basis_and_B_batch(box, spans, pts)
        Tq = eval_T_quad(N, cps, pts)
        Eq = np.array([float(eval_E(float(t))) for t in Tq], dtype=float)
        nuq = np.array([float(eval_nu(float(t))) for t in Tq], dtype=float)
        Dq = _isotropic_D_batch(Eq, nuq)

        if assemble_K:
            wDB = wts[:, None, None] * np.einsum("qrs,qsj->qrj", Dq, B)
            Ke = B.reshape(-1, n_edof).T @ wDB.reshape(-1, n_edof)
            rows.append(np.repeat(dofs, n_edof))
            cols.append(np.tile(dofs, n_edof))
            vals.append(Ke.reshape(-1))
        if f is not None:
            alphaq = np.array([float(eval_alpha(float(t))) for t in Tq], dtype=float)
            eps_th = alphaq * (Tq - float(T_ref))
            eps_star = np.zeros((Tq.size, 6), dtype=float)
            eps_star[:, :3] = eps_th[:, None]
            sig_star = np.einsum("qrs,qs->qr", Dq, eps_star)
            f[dofs] += np.einsum("q,qri,qr->i", wts, B, sig_star)

    K = (
        sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(box.n_dofs, box.n_dofs)).tocsr()
        if assemble_K and vals
        else (sp.csr_matrix((box.n_dofs, box.n_dofs), dtype=float) if assemble_K else None)
    )
    return K, f


def assemble_stiffness_tdep(
    box: BoxSpline,
    temperature: np.ndarray | Callable[[float, float, float], float],
    E_fn: float | Callable[[float], float],
    nu_fn: float | Callable[[float], float] = 0.0,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray] | None = None,
    full_cells: Sequence[tuple[int, int, int]] | None = None,
    n_gauss: int | None = None,
) -> sp.csr_matrix:
    """Assemble ``K = int B^T D(T(x)) B dV`` for spatially varying temperature-dependent ``E(T), nu(T)``.

    Args:
        box: Box B-spline (3 displacement DOFs per control point).
        temperature: Either control-point temperature vector ``(box.n_ctrl,)`` or callable ``T(x, y, z)``.
        E_fn: Constant Young's modulus or callable ``T -> E(T)``.
        nu_fn: Constant Poisson's ratio or callable ``T -> nu(T)``.
        cell_tables: Optional ``{cell: (nq, 4)}`` cut-cell quadrature tables.
        full_cells: Optional sequence of full cells (defaults to all box cells when ``cell_tables`` is ``None``).
        n_gauss: Gauss points per axis on full cells (defaults to ``box.degree + 1``).

    Returns:
        Sparse CSR stiffness matrix ``(box.n_dofs, box.n_dofs)``.
    """
    K, _ = _assemble_tdep_system(
        box,
        temperature,
        E_fn,
        nu_fn,
        cell_tables=cell_tables,
        full_cells=full_cells,
        n_gauss=n_gauss,
        assemble_K=True,
        assemble_f=False,
    )
    assert K is not None
    return K


def temperature_dependent_eigenstrain_load(
    box: BoxSpline,
    temperature: np.ndarray | Callable[[float, float, float], float],
    E_fn: float | Callable[[float], float],
    nu_fn: float | Callable[[float], float],
    alpha_fn: float | Callable[[float], float],
    T_ref: float = 0.0,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray] | None = None,
    full_cells: Sequence[tuple[int, int, int]] | None = None,
    n_gauss: int | None = None,
) -> np.ndarray:
    """Load ``f = int B^T D(T(x)) eps*(x) dV`` for spatially varying ``eps*(x) = alpha(T)(T - T_ref) I``.

    Args:
        box: Box B-spline (3 displacement DOFs per control point).
        temperature: Either control-point temperature vector ``(box.n_ctrl,)`` or callable ``T(x, y, z)``.
        E_fn: Constant Young's modulus or callable ``T -> E(T)``.
        nu_fn: Constant Poisson's ratio or callable ``T -> nu(T)``.
        alpha_fn: Constant thermal expansion coefficient or callable ``T -> alpha(T)``.
        T_ref: Stress-free reference temperature.
        cell_tables: Optional ``{cell: (nq, 4)}`` cut-cell quadrature tables.
        full_cells: Optional sequence of full cells (defaults to all box cells when ``cell_tables`` is ``None``).
        n_gauss: Gauss points per axis on full cells (defaults to ``box.degree + 1``).

    Returns:
        Force vector ``(box.n_dofs,)`` in reference DOF order.
    """
    _, f = _assemble_tdep_system(
        box,
        temperature,
        E_fn,
        nu_fn,
        alpha_fn=alpha_fn,
        T_ref=T_ref,
        cell_tables=cell_tables,
        full_cells=full_cells,
        n_gauss=n_gauss,
        assemble_K=False,
        assemble_f=True,
    )
    assert f is not None
    return f


def recover_thermoelastic_stress_tdep(
    box: BoxSpline,
    u: np.ndarray,
    temperature: np.ndarray | Callable[[float, float, float], float],
    E_fn: float | Callable[[float], float],
    nu_fn: float | Callable[[float], float],
    alpha_fn: float | Callable[[float], float],
    points: np.ndarray,
    T_ref: float = 0.0,
) -> np.ndarray:
    """Recover Cauchy stress ``sigma(x) = D(T(x)) (B(x) u - eps*(x))`` at physical points.

    Args:
        box: Box B-spline.
        u: Displacement vector in reference DOF order, shape ``(box.n_dofs,)``.
        temperature: Either control-point temperature vector ``(box.n_ctrl,)`` or callable ``T(x, y, z)``.
        E_fn: Constant Young's modulus or callable ``T -> E(T)``.
        nu_fn: Constant Poisson's ratio or callable ``T -> nu(T)``.
        alpha_fn: Constant thermal expansion coefficient or callable ``T -> alpha(T)``.
        points: Physical evaluation points, shape ``(n_pts, 3)``.
        T_ref: Stress-free reference temperature.

    Returns:
        Stress array ``(n_pts, 6)`` in Kratos Voigt order ``[xx, yy, zz, xy, yz, xz]``.
    """
    p = box.degree
    u_vec = np.asarray(u, dtype=float).reshape(-1)
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    eval_E = _resolve_scalar_material(E_fn)
    eval_nu = _resolve_scalar_material(nu_fn)
    eval_alpha = _resolve_scalar_material(alpha_fn)
    eval_T_quad = _eval_temperature_fn(box, temperature)
    nu_cp, nv_cp, _nw = box.n_cp
    out = np.zeros((pts.shape[0], 6), dtype=float)
    a = np.arange(p + 1)

    for n, x in enumerate(pts):
        coords = np.array([[min(max(float(x[d]), float(box.lo[d])), float(box.hi[d])) for d in range(3)]], dtype=float)
        spans = [int(basis_1d(box.knots[d], p, float(coords[0, d]), 0)[0]) for d in range(3)]
        cps = (
            (spans[0] - p + a)[None, None, :]
            + nu_cp * ((spans[1] - p + a)[None, :, None] + nv_cp * (spans[2] - p + a)[:, None, None])
        ).reshape(-1)
        N, B = _eval_basis_and_B_batch(box, spans, coords)
        dofs = (3 * cps[:, None] + np.arange(3)[None, :]).reshape(-1)
        eps_total = B[0] @ u_vec[dofs]
        Tq = float(eval_T_quad(N, cps, coords)[0])
        eps_th = float(eval_alpha(Tq)) * (Tq - float(T_ref))
        eps_star = np.array([eps_th, eps_th, eps_th, 0.0, 0.0, 0.0], dtype=float)
        Dq = isotropic_D(float(eval_E(Tq)), float(eval_nu(Tq)))
        out[n] = Dq @ (eps_total - eps_star)
    return out


def solve_thermoelastic_tdep(
    box: BoxSpline,
    temperature: np.ndarray | Callable[[float, float, float], float],
    E_fn: float | Callable[[float], float],
    nu_fn: float | Callable[[float], float],
    alpha_fn: float | Callable[[float], float],
    fixed: Sequence[int] | np.ndarray,
    fixed_values: Sequence[float] | np.ndarray | None = None,
    T_ref: float = 0.0,
    external_load: np.ndarray | None = None,
    cell_tables: Mapping[tuple[int, int, int], np.ndarray] | None = None,
    full_cells: Sequence[tuple[int, int, int]] | None = None,
    n_gauss: int | None = None,
) -> ThermoelasticTdepResult:
    """Solve linear thermo-elasticity with spatially varying temperature-dependent ``E(T), nu(T), alpha(T)``.

    Args:
        box: Box B-spline.
        temperature: Either control-point temperature vector ``(box.n_ctrl,)`` or callable ``T(x, y, z)``.
        E_fn: Constant Young's modulus or callable ``T -> E(T)``.
        nu_fn: Constant Poisson's ratio or callable ``T -> nu(T)``.
        alpha_fn: Constant thermal expansion coefficient or callable ``T -> alpha(T)``.
        fixed: Reference DOF indices with prescribed Dirichlet displacements.
        fixed_values: Prescribed displacement values (defaults to zeros).
        T_ref: Stress-free reference temperature.
        external_load: Optional mechanical nodal load vector ``(box.n_dofs,)`` added to thermal load.
        cell_tables: Optional ``{cell: (nq, 4)}`` cut-cell quadrature tables.
        full_cells: Optional sequence of full cells (defaults to all box cells when ``cell_tables`` is ``None``).
        n_gauss: Gauss points per axis on full cells (defaults to ``box.degree + 1``).

    Returns:
        :class:`ThermoelasticTdepResult`.
    """
    K, f_th = _assemble_tdep_system(
        box,
        temperature,
        E_fn,
        nu_fn,
        alpha_fn=alpha_fn,
        T_ref=T_ref,
        cell_tables=cell_tables,
        full_cells=full_cells,
        n_gauss=n_gauss,
        assemble_K=True,
        assemble_f=True,
    )
    assert K is not None and f_th is not None
    f_total = f_th if external_load is None else f_th + np.asarray(external_load, dtype=float).reshape(box.n_dofs)
    fixed_arr = np.asarray(fixed, dtype=np.int64)
    vals_arr = np.zeros(fixed_arr.size, dtype=float) if fixed_values is None else np.asarray(fixed_values, dtype=float)
    u = solve_with_dirichlet(K, f_total, fixed_arr, vals_arr)
    return ThermoelasticTdepResult(u=u, K=K, f=f_total)
