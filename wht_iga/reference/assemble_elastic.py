"""Reference (oracle) assembler for small-strain isotropic linear elasticity.

Independent of Kratos; used to cross-check the Kratos path (master-plan v2 §5.5)
and as the R3 fallback. Untrimmed box only in Phase 0.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from numba import njit

from wht_iga.reference.box import BoxSpline


def isotropic_D(E: float, nu: float) -> np.ndarray:
    """Voigt order [xx, yy, zz, xy, yz, xz] with engineering shear strains."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


@njit(cache=True)
def _assemble(fu, Nu, dNu, ju, fv, Nv, dNv, jv, fw, Nw, dNw, jw, wg, nu_cp, nv_cp, p, D):
    ncu, ngp, _ = Nu.shape
    ncv = Nv.shape[0]
    ncw = Nw.shape[0]
    nloc = (p + 1) ** 3
    ndof_e = 3 * nloc
    n_cells = ncu * ncv * ncw
    rows = np.empty(n_cells * ndof_e * ndof_e, dtype=np.int64)
    cols = np.empty_like(rows)
    vals = np.empty(n_cells * ndof_e * ndof_e)
    dofs = np.empty(ndof_e, dtype=np.int64)
    grad = np.empty((nloc, 3))
    B = np.zeros((6, ndof_e))
    pos = 0
    for k in range(ncw):
        for j in range(ncv):
            for i in range(ncu):
                Ke = np.zeros((ndof_e, ndof_e))
                for c in range(p + 1):
                    for b in range(p + 1):
                        for a in range(p + 1):
                            loc = a + (p + 1) * (b + (p + 1) * c)
                            I = (fu[i] + a) + nu_cp * ((fv[j] + b) + nv_cp * (fw[k] + c))
                            for d in range(3):
                                dofs[3 * loc + d] = 3 * I + d
                for gz in range(ngp):
                    for gy in range(ngp):
                        for gx in range(ngp):
                            wdet = wg[gx] * wg[gy] * wg[gz] * ju[i] * jv[j] * jw[k]
                            for c in range(p + 1):
                                for b in range(p + 1):
                                    for a in range(p + 1):
                                        loc = a + (p + 1) * (b + (p + 1) * c)
                                        grad[loc, 0] = dNu[i, gx, a] * Nv[j, gy, b] * Nw[k, gz, c]
                                        grad[loc, 1] = Nu[i, gx, a] * dNv[j, gy, b] * Nw[k, gz, c]
                                        grad[loc, 2] = Nu[i, gx, a] * Nv[j, gy, b] * dNw[k, gz, c]
                            B[:, :] = 0.0
                            for loc in range(nloc):
                                gx_, gy_, gz_ = grad[loc, 0], grad[loc, 1], grad[loc, 2]
                                c0 = 3 * loc
                                B[0, c0] = gx_
                                B[1, c0 + 1] = gy_
                                B[2, c0 + 2] = gz_
                                B[3, c0] = gy_
                                B[3, c0 + 1] = gx_
                                B[4, c0 + 1] = gz_
                                B[4, c0 + 2] = gy_
                                B[5, c0] = gz_
                                B[5, c0 + 2] = gx_
                            Ke += wdet * (B.T @ (D @ B))
                for r in range(ndof_e):
                    for s in range(ndof_e):
                        rows[pos] = dofs[r]
                        cols[pos] = dofs[s]
                        vals[pos] = Ke[r, s]
                        pos += 1
    return rows, cols, vals


def gauss_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(n)


def assemble_stiffness(box: BoxSpline, E: float, nu: float, n_gauss: int | None = None) -> sp.csr_matrix:
    ng = box.degree + 1 if n_gauss is None else n_gauss
    xg, wg = gauss_legendre(ng)
    tabs = [box.axis_tables(d, xg) for d in range(3)]
    (fu, Nu, dNu, _, ju), (fv, Nv, dNv, _, jv), (fw, Nw, dNw, _, jw) = tabs
    rows, cols, vals = _assemble(
        fu, Nu, dNu, ju, fv, Nv, dNv, jv, fw, Nw, dNw, jw, wg, box.n_cp[0], box.n_cp[1], box.degree, isotropic_D(E, nu)
    )
    n = box.n_dofs
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def assemble_face_traction(box: BoxSpline, axis: int, side: int, traction, n_gauss: int | None = None) -> np.ndarray:
    """Consistent load vector for a constant traction vector on a box face."""
    p = box.degree
    ng = p + 1 if n_gauss is None else n_gauss
    xg, wg = gauss_legendre(ng)
    others = [d for d in range(3) if d != axis]
    t1 = box.axis_tables(others[0], xg)
    t2 = box.axis_tables(others[1], xg)
    fixed_cp = 0 if side == 0 else box.n_cp[axis] - 1
    nu, nv, _ = box.n_cp
    f = np.zeros(box.n_dofs)
    traction = np.asarray(traction, dtype=float)
    for c1 in range(len(t1[4])):
        for c2 in range(len(t2[4])):
            for g1 in range(ng):
                for g2 in range(ng):
                    w = wg[g1] * wg[g2] * t1[4][c1] * t2[4][c2]
                    for a1 in range(p + 1):
                        for a2 in range(p + 1):
                            idx = [0, 0, 0]
                            idx[axis] = fixed_cp
                            idx[others[0]] = t1[0][c1] + a1
                            idx[others[1]] = t2[0][c2] + a2
                            I = idx[0] + nu * (idx[1] + nv * idx[2])
                            f[3 * I : 3 * I + 3] += w * t1[1][c1, g1, a1] * t2[1][c2, g2, a2] * traction
    return f
