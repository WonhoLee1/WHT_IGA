# -*- coding: utf-8 -*-
"""Mean-dilatation (F-bar style) Total-Lagrangian Neo-Hookean, Numba (Track H).

Orthodox Simo-Taylor-Pister mean dilatation, NOT a C++ element:
deviatoric energy integrated normally, volumetric energy sees only the
cell-averaged ``Jbar``::

    W_cell = sum_q w_q * mu/2 * (tr(B_q) - 3) + W_cell * U(Jbar)
    U(J) = lam/2 * (ln J)^2 - mu * ln J        (same law as CL.HyperElastic3DLaw)

Uniform deformation => Jbar = J pointwise => identical to stock NH
(unit-cube test must match the S19 analytic to machine precision).

Residual per node I (3-vector) and consistent tangent (dense, small models)::

    R_I = sum_q w_q * mu * F_q @ G_qI + U'(Jbar) * V_I
    V_I = sum_r w_r * J_r * (FinvT_r @ G_rI)
    K = K_dev + U''(Jbar) * dJbar/du (x) V + U'(Jbar) * dV/du

with ``dJbar/du_J = (1/Wcell) sum_s w_s J_s (FinvT_s G_sJ)`` and the
``d[J FinvT G]/du`` terms in explicit 3x3 algebra (verified against
finite differences in ``tests/unit/test_fbar_mean_dilatation.py``).

Scope: untrimmed single-patch boxes (E04 precursor). No Kratos dependency.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from wht_iga.grid.bspline import ders_basis_funs, find_span


@njit(cache=True)
def _det3(a: np.ndarray) -> float:
    return (a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1])
            - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0])
            + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]))


@njit(cache=True)
def _inv3(a: np.ndarray) -> np.ndarray:
    det = _det3(a)
    out = np.empty((3, 3), dtype=np.float64)
    out[0, 0] = (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) / det
    out[0, 1] = (a[0, 2] * a[2, 1] - a[0, 1] * a[2, 2]) / det
    out[0, 2] = (a[0, 1] * a[1, 2] - a[0, 2] * a[1, 1]) / det
    out[1, 0] = (a[1, 2] * a[2, 0] - a[1, 0] * a[2, 2]) / det
    out[1, 1] = (a[0, 0] * a[2, 2] - a[0, 2] * a[2, 0]) / det
    out[1, 2] = (a[0, 2] * a[1, 0] - a[0, 0] * a[1, 2]) / det
    out[2, 0] = (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]) / det
    out[2, 1] = (a[0, 1] * a[2, 0] - a[0, 0] * a[2, 1]) / det
    out[2, 2] = (a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]) / det
    return out


@njit(cache=True)
def _cell_residual_tangent(G: np.ndarray, W: np.ndarray, u_cell: np.ndarray,
                           mu: float, lam: float):
    """Residual (3*nen,) and tangent (3*nen, 3*nen) for one cell.

    G: (nq, nen, 3) reference gradients dN/dX. W: (nq,) reference weights.
    u_cell: (nen, 3) nodal displacements.
    """
    nq, nen, _ = G.shape
    nd = 3 * nen
    R = np.zeros(nd, dtype=np.float64)
    K = np.zeros((nd, nd), dtype=np.float64)
    F = np.zeros((nq, 3, 3), dtype=np.float64)
    detJ = np.zeros(nq, dtype=np.float64)
    for q in range(nq):
        for a in range(3):
            F[q, a, a] = 1.0
            for I in range(nen):
                for m in range(3):
                    F[q, a, m] += u_cell[I, a] * G[q, I, m]
        detJ[q] = _det3(F[q])
    Wcell = 0.0
    for q in range(nq):
        Wcell += W[q]
    Jbar = 0.0
    for q in range(nq):
        Jbar += W[q] * detJ[q]
    Jbar /= Wcell
    pbar = (lam * np.log(Jbar) - mu) / Jbar
    ppbar = (lam - lam * np.log(Jbar) + mu) / (Jbar * Jbar)
    # V_I and dJbar/du per node
    V = np.zeros((nen, 3), dtype=np.float64)
    dJbar = np.zeros((nen, 3), dtype=np.float64)
    FinvT = np.zeros((nq, 3, 3), dtype=np.float64)
    for q in range(nq):
        FinvT[q] = _inv3(F[q]).T.copy()
    for I in range(nen):
        for q in range(nq):
            wJ = W[q] * detJ[q]
            for a in range(3):
                h = 0.0
                for m in range(3):
                    h += FinvT[q, a, m] * G[q, I, m]
                V[I, a] += wJ * h
                dJbar[I, a] += wJ * h
        for a in range(3):
            dJbar[I, a] /= Wcell
    # residual: deviatoric + volumetric
    for I in range(nen):
        for a in range(3):
            val = pbar * V[I, a]
            for q in range(nq):
                fa = 0.0
                for m in range(3):
                    fa += F[q, a, m] * G[q, I, m]
                val += W[q] * mu * fa
            R[3 * I + a] = val
    # tangent part 1: deviatoric (constant)
    for I in range(nen):
        for J in range(nen):
            s = 0.0
            for q in range(nq):
                g = 0.0
                for m in range(3):
                    g += G[q, I, m] * G[q, J, m]
                s += W[q] * mu * g
            for a in range(3):
                K[3 * I + a, 3 * J + a] += s
    # tangent part 2: U'' * dJbar (x) V  +  pbar * dV/du
    for I in range(nen):
        for J in range(nen):
            for a in range(3):
                for b in range(3):
                    K[3 * I + a, 3 * J + b] += ppbar * V[I, a] * dJbar[J, b]
    for I in range(nen):
        for J in range(nen):
            for q in range(nq):
                wJ = W[q] * detJ[q]
                # H_I, H_J at this qp (F^{-T} G); dV term is
                # J*(outer(H_I,H_J) - outer(H_J,H_I)) -- verified vs FD
                HI = np.zeros(3, dtype=np.float64)
                HJ = np.zeros(3, dtype=np.float64)
                for a in range(3):
                    for m in range(3):
                        HI[a] += FinvT[q, a, m] * G[q, I, m]
                        HJ[a] += FinvT[q, a, m] * G[q, J, m]
                for a in range(3):
                    for b in range(3):
                        K[3 * I + a, 3 * J + b] += pbar * wJ * (HI[a] * HJ[b] - HJ[a] * HI[b])
    return R, K, Jbar


def _lagrange_1d(xi: float, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Values + derivatives of Lagrange basis through *nodes* at *xi*."""
    n = len(nodes)
    N = np.ones(n)
    dN = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            N[i] *= (xi - nodes[j]) / (nodes[i] - nodes[j])
        for k in range(n):
            if k == i:
                continue
            term = 1.0 / (nodes[i] - nodes[k])
            for j in range(n):
                if j == i or j == k:
                    continue
                term *= (xi - nodes[j]) / (nodes[i] - nodes[j])
            dN[i] += term
    return N, dN


def build_model(box, degree: int, E: float, nu: float,
                n_gauss: int | None = None) -> dict:
    """Precompute per-cell QP data. *box*: BoxSpline-like (pure Python ok)."""
    import numpy as np

    p = int(degree)
    ng = p + 1 if n_gauss is None else int(n_gauss)
    xg, wg = np.polynomial.legendre.leggauss(ng)
    ku, kv, kw = (np.asarray(k, dtype=float) for k in box.knots)
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    cp = np.asarray(box.control_points(), dtype=float)
    n_ctrl = cp.shape[0]
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    br = [np.unique(k) for k in (ku, kv, kw)]
    cells: list = []
    for ck in range(box.n_cells[2]):
        for cj in range(box.n_cells[1]):
            for ci in range(box.n_cells[0]):
                los = [br[0][ci], br[1][cj], br[2][ck]]
                his = [br[0][ci + 1], br[1][cj + 1], br[2][ck + 1]]
                # NOTE: uniform clamped knots assumed (BoxSpline/CookSpline).
                # cell (ci,cj,ck) spans knot intervals [ci+p, ci+p+1] etc.,
                # i.e. control points [ci:ci+p+1] per axis.
                first = [ci, cj, ck]
                sx, sy, sz = ci + p, cj + p, ck + p
                idx = []
                for ck2 in range(p + 1):
                    for cj2 in range(p + 1):
                        for ci2 in range(p + 1):
                            I = ((first[0] + ci2)
                                 + nu_ * ((first[1] + cj2) + nv_ * (first[2] + ck2)))
                            idx.append(int(I))
                ncn = len(idx)
                nq = ng ** 3
                G = np.zeros((nq, ncn, 3))
                W = np.zeros(nq)
                q = 0
                for gz, wz in zip(xg, wg):
                    uz = los[2] + 0.5 * (gz + 1.0) * (his[2] - los[2])
                    dz = ders_basis_funs(sz, uz, p, 1, kw)
                    for gy, wy in zip(xg, wg):
                        uy = los[1] + 0.5 * (gy + 1.0) * (his[1] - los[1])
                        dy = ders_basis_funs(sy, uy, p, 1, kv)
                        for gx, wx in zip(xg, wg):
                            ux = los[0] + 0.5 * (gx + 1.0) * (his[0] - los[0])
                            dx = ders_basis_funs(sx, ux, p, 1, ku)
                            Nu, dNu = np.asarray(dx[0]), np.asarray(dx[1])
                            Nv, dNv = np.asarray(dy[0]), np.asarray(dy[1])
                            Nw, dNw = np.asarray(dz[0]), np.asarray(dz[1])
                            # geometry Jacobian dX/dxi (cols = dX/dxi,dX/deta,dX/dzeta)
                            Jg = np.zeros((3, 3))
                            N3 = np.zeros((ncn,))
                            d3 = np.zeros((ncn, 3))
                            ii = 0
                            for ak in range(p + 1):
                                for aj in range(p + 1):
                                    for ai in range(p + 1):
                                        X = cp[idx[ii]]
                                        Jg[:, 0] += X * (dNu[ai] * Nv[aj] * Nw[ak])
                                        Jg[:, 1] += X * (Nu[ai] * dNv[aj] * Nw[ak])
                                        Jg[:, 2] += X * (Nu[ai] * Nv[aj] * dNw[ak])
                                        N3[ii] = Nu[ai] * Nv[aj] * Nw[ak]
                                        d3[ii, 0] = dNu[ai] * Nv[aj] * Nw[ak]
                                        d3[ii, 1] = Nu[ai] * dNv[aj] * Nw[ak]
                                        d3[ii, 2] = Nu[ai] * Nv[aj] * dNw[ak]
                                        ii += 1
                            detJ = float(np.linalg.det(Jg))
                            invJg = np.linalg.inv(Jg)
                            G[q] = (invJg.T @ d3.T).T
                            W[q] = float(wx * wy * wz * 0.125
                                         * (his[0] - los[0]) * (his[1] - los[1]) * (his[2] - los[2])
                                         * abs(detJ))
                            q += 1
                cells.append({"nodes": np.asarray(idx, dtype=np.int64), "G": G, "W": W})
    return {"cells": cells, "n_ctrl": n_ctrl, "n_dof": 3 * n_ctrl,
            "mu": float(mu), "lam": float(lam), "degree": p}


def residual_tangent(model: dict, u: np.ndarray):
    """Assemble global R, K (dense)."""
    ndof = int(model["n_dof"])
    R = np.zeros(ndof)
    K = np.zeros((ndof, ndof))
    uc = np.asarray(u, dtype=float).reshape(-1, 3)
    for cell in model["cells"]:
        nodes = cell["nodes"]
        Rc, Kc, _ = _cell_residual_tangent(
            np.ascontiguousarray(cell["G"]), np.ascontiguousarray(cell["W"]),
            np.ascontiguousarray(uc[nodes]), float(model["mu"]), float(model["lam"]))
        dofs = np.repeat(3 * nodes, 3) + np.tile([0, 1, 2], len(nodes))
        R[dofs] += Rc
        K[np.ix_(dofs, dofs)] += Kc
    return R, K


def cauchy_stress(model: dict, u: np.ndarray):
    """Per-QP Cauchy stress (list over cells of (nq, 6) Voigt xx,yy,zz,xy,yz,xz)."""
    uc = np.asarray(u, dtype=float).reshape(-1, 3)
    mu, lam = float(model["mu"]), float(model["lam"])
    out = []
    for cell in model["cells"]:
        G, W, nodes = cell["G"], cell["W"], cell["nodes"]
        nq = W.shape[0]
        S = np.zeros((nq, 6))
        F = np.zeros((nq, 3, 3))
        for q in range(nq):
            F[q] = np.eye(3) + uc[nodes].T @ G[q]
        J = np.linalg.det(F)
        Jbar = float(np.sum(W * J) / np.sum(W))
        pbar = (lam * np.log(Jbar) - mu) / Jbar
        for q in range(nq):
            B = F[q] @ F[q].T
            # NB: deviatoric part is (mu/J)*B -- the -mu/J*I lives inside
            # pbar = U'(Jbar). Splitting (B-I) here double-subtracts it.
            # (Verified: uniform field then matches S19 analytic exactly.)
            sig = mu / J[q] * B + pbar * np.eye(3)
            S[q] = [sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]]
        out.append(S)
    return out


def solve_static(model: dict, f_ext: np.ndarray, fixed: np.ndarray,
                 n_steps: int = 1, rtol: float = 1e-6, atol: float = 1e-9,
                 max_it: int = 25, prescribed: dict | None = None,
                 dtol: float = 1e-10, min_frac: float = 1e-4,
                 history: list | None = None) -> np.ndarray:
    """Adaptive load-stepped Newton with Dirichlet fixing. Returns full u.

    ``prescribed`` maps DOF -> final value (ramped with the load factor);
    plain ``fixed`` DOFs stay at 0. The nominal increment is ``1/n_steps``.

    Convergence per increment: ``|R| <= atol + rtol * max(|R0|, |f|)`` (max
    norm) OR the displacement-increment test ``|du| <= dtol * |u|`` (max norm,
    after at least one update). The second test is needed because at
    nu ~ 0.5 the assembled residual stalls at a round-off floor (~1e-7 to
    1e-6 at 1e6-scale bulk stiffness) that can sit above the residual gate on
    fine meshes (measured 2026-11-19), the same reason Kratos E04 runs use
    ``DisplacementCriteria``.

    On failure (no convergence in ``max_it`` or no admissible step) the
    increment is restored and halved; after an easy increment it grows back
    toward ``1/n_steps``. Measured 2026-11-19 on Cook nu=0.49999: whether a
    fixed-step Newton converges depends on the step count (4x4 p=1 converges
    with 5 steps, diverges with 20), so a fixed schedule is not robust.

    Line search is guard-only: the step is shortened only when the trial
    state is inadmissible (non-finite residual, i.e. Jbar <= 0). A monotone
    residual test is wrong here: the first iterate at nu ~ 0.5 raises the
    residual ~1e4x and then converges quadratically; the old max-norm test
    stalled at step 1/64.

    ``history`` (optional list) receives ``(load_factor, iterations)`` for
    each converged increment.
    """
    ndof = int(model["n_dof"])
    u = np.zeros(ndof)
    fixed = np.asarray(fixed, dtype=np.int64)
    prescribed = {int(k): float(v) for k, v in (prescribed or {}).items()}
    free = np.ones(ndof, dtype=bool)
    free[fixed] = False
    for k in list(prescribed.keys()):
        if k < 0 or k >= ndof:
            raise ValueError(f"prescribed dof {k} out of range")
        free[k] = False
    f = np.asarray(f_ext, dtype=float)
    d_nom = 1.0 / max(1, int(n_steps))
    lam, dlam = 0.0, d_nom

    def _increment(lam_new: float) -> tuple[bool, int]:
        target = f * lam_new
        for k, v in prescribed.items():
            u[k] = v * lam_new
        R, K = residual_tangent(model, u)
        rr = R[free] - target[free]
        if not np.all(np.isfinite(rr)):
            return False, 0
        gate = atol + rtol * max(float(np.max(np.abs(rr))),
                                 float(np.max(np.abs(target[free]))), 1e-300)
        for it in range(1, max_it + 1):
            if float(np.max(np.abs(rr))) < gate:
                return True, it - 1
            du = np.linalg.solve(K[np.ix_(free, free)], -rr)
            step = 1.0
            while True:
                u_try = u.copy()
                u_try[free] += step * du
                R, K = residual_tangent(model, u_try)
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

    while lam < 1.0 - 1e-14:
        lam_new = min(1.0, lam + dlam)
        u_save = u.copy()
        ok, its = _increment(lam_new)
        if ok:
            lam = lam_new
            if history is not None:
                history.append((lam, its))
            if its <= 6:
                dlam = min(d_nom, 1.5 * dlam)
        else:
            u[:] = u_save
            dlam *= 0.5
            if dlam < min_frac:
                raise RuntimeError(f"fbar Newton failed: increment below {min_frac} at load factor {lam:.6g}")
    return u
