# -*- coding: utf-8 -*-
"""Generalised mean dilatation: cell-wise discontinuous polynomial pressure
(P0 / P1 / P2-disc), Total-Lagrangian Neo-Hookean, Numba + CSR (S49).

Three-field Hu-Washizu (u, theta, p) with theta and p in the same cell-wise
discontinuous space ``Q_c = span{phi_k}``::

    Pi(u, theta, p) = int W_dev(F) dV + int U(theta) dV + int p (J - theta) dV
    W_dev = mu/2 (tr C - 3),   U(J) = lam/2 (ln J)^2 - mu ln J   (as fbar.py)

Stationarity in p and theta is local to each cell and has a closed form
(static condensation, no global pressure DOFs)::

    theta = Pi_Q J,     p = Pi_Q U'(theta)
    E(u)  = sum_q w_q mu/2 (tr C_q - 3) + sum_q w_q U(theta_q)

with ``Pi_Q`` the cell L2 projection (cell quadrature, reference volume).
``Q = constants`` gives ``theta = Jbar``: exactly ``wht_iga.analyses.fbar``
(Simo-Taylor-Pister mean dilatation). ``Q = P1`` is the Q2/P1disc-type
generalised mean dilatation.

Basis. ``Psi`` (nq, m) per cell, orthonormal in the cell quadrature
(``sum_q w_q Psi_qk Psi_ql = delta_kl``), so ``Pi_Q f = Psi (Psi^T W f)``.
Built once from complete monomials of total degree <= k, centred at the cell
centroid and scaled per axis (conditioning), in

* ``coords="physical"`` (default): reference coordinates X ("unmapped"). P_k
  then contains every polynomial of degree <= k on any cell shape, so the
  projection error is O(h^(k+1)) on curved / non-affine cells too.
* ``coords="parametric"``: knot-span coordinates xi ("mapped"). Optimal only
  on asymptotically affine meshes (Arnold, Boffi & Falk 2002; Boffi &
  Gastaldi 2002 for Q2-P1 Stokes); kept for comparison. On the smoothly
  mapped V07 plate both agree (p=2 C^0 P1, h=0.125: sigma_xx error
  1.017e-3 physical, 1.001e-3 parametric).

Residual and consistent tangent (H_qI = F_q^-T G_qI)::

    R_I = sum_q w_q [mu F_q G_qI + phat_q J_q H_qI],  phat = Pi_Q U'(theta)
    K   = K_dev + A^T C A + sum_q w_q phat_q d(J_q H_qI)/du_J
    A_k(Ia) = sum_q w_q Psi_qk J_q H_qIa,   C_kl = sum_q w_q Psi_qk Psi_ql U''(theta_q)

K_dev and the d(J H)/du block are those of ``fbar._cell_residual_tangent``
with the cell constant pbar replaced by the QP value phat_q; the K is
symmetric (hyperelastic). Cauchy stress: ``sigma_q = mu/J_q B_q + phat_q I``.

Inf-sup (S49, V07 Kirsch plate, discrete beta_h of the pair at h = 0.5,
1/3, 0.25): stable, beta_h h-independent, for p=2 C^1 / P0 (= fbar),
p=2 C^0 / P1 (Q2/P1disc-like) and p=3 C^0 / P2 (Q3/P2disc-like); weakly
unstable for p=3 C^1 / P1 (beta_h 5.7e-3 -> 4.8e-3, sigma error 2x the
nu=0.3 value); locking for p=2 C^1 / P1, p=3 C^2 / P1 (beta_h ~ 2e-5, about
2 in-plane displacement DOFs per cell against 3 in-plane P1 modes),
p=2 C^0 / P2 and p=3 C^1 / P2 (beta_h decaying ~ h^1.4 .. h^1.6 from h = 0.5
to 1/3). The space is not chosen automatically: callers pick a stable pair.

Geometry, quadrature, cell -> control-point map and the CSR pattern are
those of ``fbar.build_model`` (imported, not copied). No Kratos dependency.
"""

from __future__ import annotations

import itertools

import numpy as np
from numba import njit

from wht_iga.analyses import fbar
from wht_iga.analyses.fbar import _det3, _inv3

PRESSURE_SPACES = ("P0", "P1", "P2")


def pressure_exponents(space) -> tuple[tuple[int, int, int], ...]:
    """Monomial exponents (i, j, k) of a named space or an explicit list.

    ``"P0"`` (1 mode), ``"P1"`` (4: 1, x, y, z), ``"P2"`` (10: total degree
    <= 2). An explicit sequence of 3-tuples is used as given (must contain
    the constant (0, 0, 0)).
    """
    if isinstance(space, str):
        deg = {"P0": 0, "P1": 1, "P2": 2}.get(space.upper())
        if deg is None:
            raise ValueError(f"unknown pressure space {space!r}; use {PRESSURE_SPACES} or exponents")
        exps = [e for e in itertools.product(range(deg + 1), repeat=3) if sum(e) <= deg]
        exps.sort(key=lambda e: (sum(e), [-v for v in e]))
        return tuple(tuple(int(v) for v in e) for e in exps)
    exps = tuple(tuple(int(v) for v in e) for e in space)
    if (0, 0, 0) not in exps or any(len(e) != 3 or min(e) < 0 for e in exps):
        raise ValueError("explicit pressure exponents must be 3-tuples >= 0 incl. (0, 0, 0)")
    if len(set(exps)) != len(exps):
        raise ValueError("duplicate pressure exponents")
    return exps


def _qp_points(box, p: int, ng: int, conn_ref: np.ndarray):
    """Reference (X) and parametric (xi) QP coordinates, fbar cell order.

    Same loops and ``fbar._axis_cells`` tables as ``fbar.build_model``; the
    rebuilt connectivity is checked against ``conn_ref`` so the QP order
    cannot drift from the model's G/W.
    """
    xg, _ = np.polynomial.legendre.leggauss(ng)
    knots = [np.ascontiguousarray(k, dtype=np.float64) for k in box.knots]
    nu_, nv_, _ = (int(v) for v in box.n_cp)
    cp = np.asarray(box.control_points(), dtype=float)
    tx, ty, tz = (fbar._axis_cells(k, p, xg) for k in knots)
    nq, ncn = ng ** 3, (p + 1) ** 3
    ax = np.arange(p + 1)
    X_all = np.empty((len(conn_ref), nq, 3))
    xi_all = np.empty((len(conn_ref), nq, 3))
    c = 0
    for zl, zh, fz, Nw, _ in tz:
        for yl, yh, fy, Nv, _ in ty:
            for xl, xh, fx, Nu, _ in tx:
                idx = ((fx + ax)[None, None, :]
                       + nu_ * ((fy + ax)[None, :, None] + nv_ * (fz + ax)[:, None, None])).reshape(ncn)
                if not np.array_equal(idx, conn_ref[c]):
                    raise RuntimeError("cell order differs from fbar.build_model")
                N3 = np.einsum("zk,yj,xi->zyxkji", Nw, Nv, Nu).reshape(nq, ncn)
                X_all[c] = N3 @ cp[idx]
                px = xl + 0.5 * (xg + 1.0) * (xh - xl)
                py = yl + 0.5 * (xg + 1.0) * (yh - yl)
                pz = zl + 0.5 * (xg + 1.0) * (zh - zl)
                Z, Y, Xp = np.meshgrid(pz, py, px, indexing="ij")  # q = gx + ng*(gy + ng*gz)
                xi_all[c] = np.column_stack([Xp.ravel(), Y.ravel(), Z.ravel()])
                c += 1
    if c != len(conn_ref):
        raise RuntimeError("cell count differs from fbar.build_model")
    return X_all, xi_all


def orthonormal_basis(Y: np.ndarray, W: np.ndarray, exps) -> np.ndarray:
    """Per-cell basis Psi (nc, nq, m), ``Psi^T diag(W) Psi = I`` per cell.

    Y: (nc, nq, 3) coordinates the monomials are taken in, W: (nc, nq)
    weights. Monomials of ``(Y - centroid) / spread`` (per axis), then
    Cholesky orthonormalisation of the cell Gram matrix.
    """
    Wsum = W.sum(axis=1)
    Yc = np.einsum("cq,cqd->cd", W, Y) / Wsum[:, None]
    D = Y - Yc[:, None, :]
    s = np.sqrt(np.einsum("cq,cqd->cd", W, D * D) / Wsum[:, None])
    if np.any(s <= 0.0):
        raise ValueError("degenerate cell: zero spread of quadrature points along an axis")
    Dn = D / s[:, None, :]
    Phi = np.ones(Y.shape[:2] + (len(exps),))
    for k, e in enumerate(exps):
        for d in range(3):
            if e[d]:
                Phi[:, :, k] *= Dn[:, :, d] ** e[d]
    M = np.einsum("cqk,cq,cql->ckl", Phi, W, Phi)
    ev = np.linalg.eigvalsh(M)
    if np.any(ev[:, 0] <= 1e-12 * ev[:, -1]):
        raise ValueError("pressure space not unisolvent on the cell quadrature")
    L = np.linalg.cholesky(M)
    Psi = np.linalg.solve(L, np.transpose(Phi, (0, 2, 1)))  # (nc, m, nq) = L^-1 Phi^T
    return np.ascontiguousarray(np.transpose(Psi, (0, 2, 1)))


@njit(cache=True)
def _kinematics(G, u_cell):
    nq, nen, _ = G.shape
    F = np.zeros((nq, 3, 3), dtype=np.float64)
    detJ = np.zeros(nq, dtype=np.float64)
    H = np.zeros((nq, nen, 3), dtype=np.float64)
    for q in range(nq):
        for a in range(3):
            F[q, a, a] = 1.0
            for I in range(nen):
                for m in range(3):
                    F[q, a, m] += u_cell[I, a] * G[q, I, m]
        detJ[q] = _det3(F[q])
        FinvT = _inv3(F[q]).T
        for I in range(nen):
            for a in range(3):
                h = 0.0
                for m in range(3):
                    h += FinvT[a, m] * G[q, I, m]
                H[q, I, a] = h
    return F, detJ, H


@njit(cache=True)
def _project(Psi, W, f):
    """Pi_Q f at the QPs: Psi (Psi^T W f)."""
    nq, m = Psi.shape
    coef = np.zeros(m, dtype=np.float64)
    for q in range(nq):
        for k in range(m):
            coef[k] += W[q] * Psi[q, k] * f[q]
    out = np.zeros(nq, dtype=np.float64)
    for q in range(nq):
        for k in range(m):
            out[q] += Psi[q, k] * coef[k]
    return out


@njit(cache=True)
def _volumetric(Psi, W, detJ, mu, lam):
    """theta = Pi J, phat = Pi U'(theta), U''(theta); ok=False if theta <= 0."""
    nq = W.shape[0]
    theta = _project(Psi, W, detJ)
    Up = np.zeros(nq, dtype=np.float64)
    Upp = np.zeros(nq, dtype=np.float64)
    ok = True
    for q in range(nq):
        th = theta[q]
        if not th > 0.0:
            ok = False
            th = 1.0
        lt = np.log(th)
        Up[q] = (lam * lt - mu) / th
        Upp[q] = (lam - lam * lt + mu) / (th * th)
    phat = _project(Psi, W, Up)
    return theta, phat, Upp, ok


@njit(cache=True)
def _cell_residual_tangent(G, W, Psi, u_cell, mu, lam):
    """Residual (3*nen,) and tangent (3*nen, 3*nen) for one cell.

    G: (nq, nen, 3) dN/dX, W: (nq,) reference weights, Psi: (nq, m)
    orthonormal pressure basis, u_cell: (nen, 3). Non-finite R when some
    theta_q <= 0 (inadmissible state; the Newton guard then shortens).
    """
    nq, nen, _ = G.shape
    m = Psi.shape[1]
    nd = 3 * nen
    R = np.zeros(nd, dtype=np.float64)
    K = np.zeros((nd, nd), dtype=np.float64)
    F, detJ, H = _kinematics(G, u_cell)
    theta, phat, Upp, ok = _volumetric(Psi, W, detJ, mu, lam)
    if not ok:
        R[:] = np.nan
        return R, K
    # residual
    for I in range(nen):
        for a in range(3):
            val = 0.0
            for q in range(nq):
                fa = 0.0
                for mm in range(3):
                    fa += F[q, a, mm] * G[q, I, mm]
                val += W[q] * (mu * fa + phat[q] * detJ[q] * H[q, I, a])
            R[3 * I + a] = val
    # deviatoric tangent (as fbar): mu * sum_q w G_I.G_J on the diagonal blocks
    for I in range(nen):
        for J in range(I, nen):
            s = 0.0
            for q in range(nq):
                g = 0.0
                for mm in range(3):
                    g += G[q, I, mm] * G[q, J, mm]
                s += W[q] * mu * g
            for a in range(3):
                K[3 * I + a, 3 * J + a] += s
                if J != I:
                    K[3 * J + a, 3 * I + a] += s
    # A^T C A
    A = np.zeros((m, nd), dtype=np.float64)
    for q in range(nq):
        wj = W[q] * detJ[q]
        for k in range(m):
            c = wj * Psi[q, k]
            for I in range(nen):
                for a in range(3):
                    A[k, 3 * I + a] += c * H[q, I, a]
    C = np.zeros((m, m), dtype=np.float64)
    for q in range(nq):
        for k in range(m):
            for ll in range(m):
                C[k, ll] += W[q] * Psi[q, k] * Psi[q, ll] * Upp[q]
    CA = np.zeros((m, nd), dtype=np.float64)
    for k in range(m):
        for ll in range(m):
            ckl = C[k, ll]
            for s_ in range(nd):
                CA[k, s_] += ckl * A[ll, s_]
    for r in range(nd):
        for k in range(m):
            ar = A[k, r]
            if ar != 0.0:
                for s_ in range(nd):
                    K[r, s_] += ar * CA[k, s_]
    # phat * d(J F^-T G_I)/du_J: antisymmetric 3x3 blocks, zero diagonal
    # blocks (see fbar._cell_residual_tangent), pbar -> phat_q.
    for I in range(nen):
        r0, r1, r2 = 3 * I, 3 * I + 1, 3 * I + 2
        for J in range(I + 1, nen):
            c0, c1, c2 = 3 * J, 3 * J + 1, 3 * J + 2
            for q in range(nq):
                cq = phat[q] * W[q] * detJ[q]
                hi0, hi1, hi2 = H[q, I, 0], H[q, I, 1], H[q, I, 2]
                hj0, hj1, hj2 = H[q, J, 0], H[q, J, 1], H[q, J, 2]
                v01 = cq * (hi0 * hj1 - hj0 * hi1)
                v02 = cq * (hi0 * hj2 - hj0 * hi2)
                v12 = cq * (hi1 * hj2 - hj1 * hi2)
                K[r0, c1] += v01
                K[r1, c0] -= v01
                K[r0, c2] += v02
                K[r2, c0] -= v02
                K[r1, c2] += v12
                K[r2, c1] -= v12
                K[c1, r0] += v01
                K[c0, r1] -= v01
                K[c2, r0] += v02
                K[c0, r2] -= v02
                K[c2, r1] += v12
                K[c1, r2] -= v12
    return R, K


@njit(cache=True)
def _assemble_csr(G_all, W_all, Psi_all, conn, u, mu, lam, indptr, node_pos, nnz):
    """Global R and CSR data of K on the ``fbar._csr_pattern`` layout."""
    nc, nen = conn.shape
    R = np.zeros(u.shape[0], dtype=np.float64)
    data = np.zeros(nnz, dtype=np.float64)
    uc = np.empty((nen, 3), dtype=np.float64)
    for c in range(nc):
        for I in range(nen):
            n = conn[c, I]
            for a in range(3):
                uc[I, a] = u[3 * n + a]
        Rc, Kc = _cell_residual_tangent(G_all[c], W_all[c], Psi_all[c], uc, mu, lam)
        for I in range(nen):
            n = conn[c, I]
            for a in range(3):
                R[3 * n + a] += Rc[3 * I + a]
                rs = indptr[3 * n + a]
                for J in range(nen):
                    p0 = rs + 3 * node_pos[c, I, J]
                    for b in range(3):
                        data[p0 + b] += Kc[3 * I + a, 3 * J + b]
    return R, data


@njit(cache=True)
def _cell_stress(G, W, Psi, u_cell, mu, lam):
    nq = W.shape[0]
    F, detJ, _ = _kinematics(G, u_cell)
    theta, phat, _, ok = _volumetric(Psi, W, detJ, mu, lam)
    S = np.zeros((nq, 6), dtype=np.float64)
    if not ok:
        S[:, :] = np.nan
        return S, theta, phat
    B = np.zeros((3, 3), dtype=np.float64)
    for q in range(nq):
        for a in range(3):
            for b in range(3):
                B[a, b] = F[q, a, 0] * F[q, b, 0] + F[q, a, 1] * F[q, b, 1] + F[q, a, 2] * F[q, b, 2]
        cj = mu / detJ[q]
        S[q, 0] = cj * B[0, 0] + phat[q]
        S[q, 1] = cj * B[1, 1] + phat[q]
        S[q, 2] = cj * B[2, 2] + phat[q]
        S[q, 3] = cj * B[0, 1]
        S[q, 4] = cj * B[1, 2]
        S[q, 5] = cj * B[0, 2]
    return S, theta, phat


def build_model(box, degree: int, E: float, nu: float, pressure="P1",
                coords: str = "physical", n_gauss: int | None = None) -> dict:
    """``fbar.build_model`` plus the per-cell orthonormal pressure basis.

    ``pressure``: ``"P0"``, ``"P1"``, ``"P2"`` or explicit exponents (see
    :func:`pressure_exponents`); ``coords``: ``"physical"`` (reference X,
    default) or ``"parametric"`` (knot-span xi).
    """
    p = int(degree)
    ng = p + 1 if n_gauss is None else int(n_gauss)
    model = fbar.build_model(box, p, E, nu, n_gauss=ng)
    exps = pressure_exponents(pressure)
    X_all, xi_all = _qp_points(box, p, ng, model["conn"])
    kind = str(coords).lower()
    if kind not in ("physical", "parametric"):
        raise ValueError(f"coords must be 'physical' or 'parametric', not {coords!r}")
    Y = X_all if kind == "physical" else xi_all
    model.update({
        "Psi_all": orthonormal_basis(Y, model["W_all"], exps),
        "X_all": X_all, "pressure": pressure if isinstance(pressure, str) else "custom",
        "pressure_exponents": exps, "n_pressure": len(exps), "coords": kind,
    })
    return model


def pressure_projection(model: dict, f_q: np.ndarray) -> np.ndarray:
    """Cell-wise L2 projection Pi_Q of a QP field f_q (nc, nq)."""
    Psi, W = model["Psi_all"], model["W_all"]
    coef = np.einsum("cqk,cq,cq->ck", Psi, W, np.asarray(f_q, dtype=float))
    return np.einsum("cqk,ck->cq", Psi, coef)


def residual_tangent(model: dict, u: np.ndarray):
    """Dense global R, K (oracle for the sparse path; small models only)."""
    ndof = int(model["n_dof"])
    R = np.zeros(ndof)
    K = np.zeros((ndof, ndof))
    uc = np.asarray(u, dtype=float).reshape(-1, 3)
    mu, lam = float(model["mu"]), float(model["lam"])
    for c, nodes in enumerate(model["conn"]):
        Rc, Kc = _cell_residual_tangent(model["G_all"][c], model["W_all"][c], model["Psi_all"][c],
                                        np.ascontiguousarray(uc[nodes]), mu, lam)
        dofs = np.repeat(3 * nodes, 3) + np.tile([0, 1, 2], len(nodes))
        R[dofs] += Rc
        K[np.ix_(dofs, dofs)] += Kc
    return R, K


def residual_tangent_sparse(model: dict, u: np.ndarray):
    """Global R (ndarray) and K (scipy CSR, sorted indices, fbar pattern)."""
    import scipy.sparse as sp

    ndof = int(model["n_dof"])
    s = model["csr"]
    R, data = _assemble_csr(model["G_all"], model["W_all"], model["Psi_all"], model["conn"],
                            np.ascontiguousarray(u, dtype=np.float64).reshape(ndof),
                            float(model["mu"]), float(model["lam"]),
                            s["indptr"], s["node_pos"], s["nnz"])
    K = sp.csr_matrix((data, s["indices"], s["indptr"]), shape=(ndof, ndof))
    K.has_sorted_indices = True
    return R, K


def cauchy_stress(model: dict, u: np.ndarray, with_pressure: bool = False):
    """Per-QP Cauchy stress, list over cells of (nq, 6) Voigt xx,yy,zz,xy,yz,xz.

    ``sigma_q = mu/J_q B_q + phat_q I`` with ``phat = Pi_Q U'(Pi_Q J)``, the
    condensed Hu-Washizu pressure. ``with_pressure`` also returns the
    (nc, nq) arrays theta and phat.
    """
    uc = np.asarray(u, dtype=float).reshape(-1, 3)
    mu, lam = float(model["mu"]), float(model["lam"])
    out, th_all, ph_all = [], [], []
    for c, nodes in enumerate(model["conn"]):
        S, th, ph = _cell_stress(model["G_all"][c], model["W_all"][c], model["Psi_all"][c],
                                 np.ascontiguousarray(uc[nodes]), mu, lam)
        out.append(S)
        th_all.append(th)
        ph_all.append(ph)
    if with_pressure:
        return out, np.stack(th_all), np.stack(ph_all)
    return out


def solve_linearized(model: dict, f_ext: np.ndarray, fixed: np.ndarray,
                     prescribed: dict | None = None, linear_solver: str = "auto") -> np.ndarray:
    """Small-strain limit: one solve with the consistent tangent at u = 0.

    ``K(0) u = f`` with Dirichlet lifting (the reference state is stress
    free, R(0) = 0). This is the linear-elastic problem of the same
    discretisation, the bilinear form being::

        mu (grad u, grad v) + mu (grad u^T, grad v) - mu (div u, div v)
        + (lam + mu) (Pi div u, Pi div v)

    (= 2 mu eps:eps + lam div div when Pi = I). V07 is specified as linear
    elastic; the Neo-Hookean Newton solution differs from it by O(|grad u|).
    """
    ndof = int(model["n_dof"])
    prescribed = {int(k): float(v) for k, v in (prescribed or {}).items()}
    free = np.ones(ndof, dtype=bool)
    free[np.asarray(fixed, dtype=np.int64)] = False
    pk = np.array(sorted(prescribed), dtype=np.int64)
    free[pk] = False
    u = np.zeros(ndof)
    u[pk] = [prescribed[k] for k in pk]
    kind = str(linear_solver).lower()
    if kind == "dense":
        _, K = residual_tangent(model, np.zeros(ndof))
        rhs = np.asarray(f_ext, dtype=float)[free] - K[free] @ u
        u[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)
        return u
    if kind not in ("auto", "pardiso", "superlu"):
        raise ValueError(f"unknown linear_solver {linear_solver!r}")
    _, K = residual_tangent_sparse(model, np.zeros(ndof))
    fi = np.flatnonzero(free)
    Kf = K[fi]
    rhs = np.asarray(f_ext, dtype=float)[fi] - Kf @ u
    A = Kf[:, fi]
    del Kf, K
    if kind != "superlu":
        try:
            import pypardiso
        except ImportError:
            if kind == "pardiso":
                raise
            pypardiso = None
        if pypardiso is not None:
            ps = pypardiso.PyPardisoSolver()
            try:
                u[fi] = pypardiso.spsolve(A, rhs, solver=ps)
            finally:
                ps.free_memory(everything=True)
            return u
    import scipy.sparse.linalg as spla
    u[fi] = spla.spsolve(A.tocsc(), rhs)
    return u


def cauchy_stress_linearized(model: dict, u: np.ndarray) -> np.ndarray:
    """Small-strain stress of :func:`solve_linearized`, (nc, nq, 6) Voigt.

    Linearisation of ``mu/J B + phat I`` at u = 0::

        sigma = mu (grad u + grad u^T) - mu div u I + (lam + mu) Pi(div u) I
    """
    uc = np.asarray(u, dtype=float).reshape(-1, 3)
    mu, lam = float(model["mu"]), float(model["lam"])
    Gu = np.einsum("cqnm,cna->cqam", model["G_all"], uc[model["conn"]])  # du_a/dX_m
    div = np.einsum("cqaa->cq", Gu)
    pdiv = pressure_projection(model, div)
    sym = mu * (Gu + np.transpose(Gu, (0, 1, 3, 2)))
    vol = -mu * div + (lam + mu) * pdiv
    S = np.empty(Gu.shape[:2] + (6,))
    S[..., 0] = sym[..., 0, 0] + vol
    S[..., 1] = sym[..., 1, 1] + vol
    S[..., 2] = sym[..., 2, 2] + vol
    S[..., 3] = sym[..., 0, 1]
    S[..., 4] = sym[..., 1, 2]
    S[..., 5] = sym[..., 0, 2]
    return S


def solve_static(model: dict, f_ext: np.ndarray, fixed: np.ndarray,
                 n_steps: int = 1, rtol: float = 0.0, atol: float = 0.0,
                 max_it: int = 25, prescribed: dict | None = None,
                 dtol: float = 1e-10, min_frac: float = 1e-4,
                 history: list | None = None,
                 linear_solver: str = "auto",
                 predictor: bool = True) -> np.ndarray:
    """Adaptive load-stepped Newton, the loop of ``fbar.solve_static`` plus a
    consistent linear predictor.

    Adapted from ``wht_iga.analyses.fbar.solve_static`` (2026-11-19): that
    function binds fbar's own assemblers, so the loop is repeated here with
    this module's ``residual_tangent(_sparse)``. Convergence per increment:
    ``|R| <= atol + rtol * max(|R_start|, |f|)`` or ``|du| <= dtol * |u|``
    (max norms); guard-only line search (halve while R is non-finite, i.e.
    some theta <= 0); failed increments are halved, easy ones grow back.

    Two deliberate differences from fbar (S49, measured 2026-11-20):

    * ``predictor=True``: each increment starts from
      ``u + K(u)^-1 (f - R(u) - K du_D)`` (tangent of the last converged
      state, prescribed increment du_D as Dirichlet lifting) instead of
      moving only the prescribed DOFs. From the latter start (boundary moved,
      interior not) the first full-step iterate of an affine patch at
      nu = 0.49999, J = 1.02 diverges (|du|/|u| ~ 12) and P1 / P2 needed
      16 / ~150 cut-back increments; with the predictor Newton reaches the
      round-off floor in 3 iterations. From the stress-free state the
      predictor is exactly :func:`solve_linearized`.
    * ``rtol = atol = 0`` (increment test only): fbar's ``rtol=1e-6``
      relative to the huge first residual of a prescribed displacement
      stops early (V07 Kirsch p=2 C^1 P0, h=0.125: sigma_xx error 1.379e-2
      instead of the converged 1.251e-2).

    ``linear_solver``: ``"auto"`` / ``"pardiso"`` / ``"superlu"`` (CSR) or
    ``"dense"`` (oracle).
    """
    ndof = int(model["n_dof"])
    u = np.zeros(ndof)
    fixed = np.asarray(fixed, dtype=np.int64)
    prescribed = {int(k): float(v) for k, v in (prescribed or {}).items()}
    free = np.ones(ndof, dtype=bool)
    free[fixed] = False
    for k in prescribed:
        if k < 0 or k >= ndof:
            raise ValueError(f"prescribed dof {k} out of range")
        free[k] = False
    f = np.asarray(f_ext, dtype=float)
    kind = str(linear_solver).lower()
    pardiso = None
    if kind == "dense":
        assemble = residual_tangent

        def lin_solve(K, b):
            return np.linalg.solve(K[np.ix_(free, free)], b)
    elif kind in ("auto", "pardiso", "superlu"):
        assemble = residual_tangent_sparse
        free_idx = np.flatnonzero(free)
        if kind != "superlu":
            try:
                import pypardiso
                pardiso = pypardiso.PyPardisoSolver()
            except ImportError:
                if kind == "pardiso":
                    raise

        def lin_solve(K, b):
            A = K[free_idx][:, free_idx]
            if pardiso is not None:
                return np.asarray(pypardiso.spsolve(A, b, solver=pardiso), dtype=np.float64)
            import scipy.sparse.linalg as spla
            return np.asarray(spla.spsolve(A.tocsc(), b), dtype=np.float64)
    else:
        raise ValueError(f"unknown linear_solver {linear_solver!r}")
    d_nom = 1.0 / max(1, int(n_steps))
    lam_f, dlam = 0.0, d_nom
    pk = np.array(sorted(prescribed), dtype=np.int64)
    pv = np.array([prescribed[k] for k in pk], dtype=np.float64)

    def _increment(lam_new: float) -> tuple[bool, int]:
        target = f * lam_new
        if predictor:
            # consistent linear predictor: tangent of the last converged
            # state, prescribed increment as Dirichlet lifting
            R0, K0 = assemble(model, u)
            d = np.zeros(ndof)
            d[pk] = pv * lam_new - u[pk]
            Kd = K0 @ d
            rhs = -(R0[free] - target[free]) - Kd[free]
            u_try = u.copy()
            u_try[pk] = pv * lam_new
            u_try[free] += lin_solve(K0, rhs)
            del R0, K0, Kd, rhs  # keep one tangent alive (peak memory)
            R, K = assemble(model, u_try)
            if np.all(np.isfinite(R[free])):
                u[:] = u_try
            else:  # inadmissible predictor -> plain start
                u[pk] = pv * lam_new
                R, K = assemble(model, u)
        else:
            u[pk] = pv * lam_new
            R, K = assemble(model, u)
        rr = R[free] - target[free]
        if not np.all(np.isfinite(rr)):
            return False, 0
        gate = atol + rtol * max(float(np.max(np.abs(rr))),
                                 float(np.max(np.abs(target[free]))), 1e-300)
        for it in range(1, max_it + 1):
            if float(np.max(np.abs(rr))) < gate:
                return True, it - 1
            du = lin_solve(K, -rr)
            step = 1.0
            while True:
                u_try = u.copy()
                u_try[free] += step * du
                R, K = assemble(model, u_try)
                if np.all(np.isfinite(R[free])):
                    break
                step *= 0.5
                if step < 1.0 / 1024.0:
                    return False, it
            u[:] = u_try
            rr = R[free] - target[free]
            if float(np.max(np.abs(step * du))) <= dtol * max(float(np.max(np.abs(u[free]))), 1e-300):
                return True, it
        return float(np.max(np.abs(rr))) < gate, max_it

    try:
        while lam_f < 1.0 - 1e-14:
            lam_new = min(1.0, lam_f + dlam)
            u_save = u.copy()
            ok, its = _increment(lam_new)
            if ok:
                lam_f = lam_new
                if history is not None:
                    history.append((lam_f, its))
                if its <= 6:
                    dlam = min(d_nom, 1.5 * dlam)
            else:
                u[:] = u_save
                dlam *= 0.5
                if dlam < min_frac:
                    raise RuntimeError(f"mixed_p1 Newton failed: increment below {min_frac} at load factor {lam_f:.6g}")
    finally:
        if pardiso is not None:
            pardiso.free_memory(everything=True)
    return u
