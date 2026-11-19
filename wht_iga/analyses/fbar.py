# -*- coding: utf-8 -*-
"""Mean-dilatation (F-bar style) Total-Lagrangian Neo-Hookean, Numba (Track H).

Orthodox Simo-Taylor-Pister mean dilatation, NOT a C++ element:
deviatoric energy integrated normally, volumetric energy sees only the
cell-averaged ``Jbar``::

    W_cell = sum_q w_q * mu/2 * (tr(B_q) - 3) + W_cell * U(Jbar)
    U(J) = lam/2 * (ln J)^2 - mu * ln J        (same law as CL.HyperElastic3DLaw)

Uniform deformation => Jbar = J pointwise => identical to stock NH
(unit-cube test must match the S19 analytic to machine precision).

Residual per node I (3-vector) and consistent tangent::

    R_I = sum_q w_q * mu * F_q @ G_qI + U'(Jbar) * V_I
    V_I = sum_r w_r * J_r * (FinvT_r @ G_rI)
    K = K_dev + U''(Jbar) * dJbar/du (x) V + U'(Jbar) * dV/du

with ``dJbar/du_J = (1/Wcell) sum_s w_s J_s (FinvT_s G_sJ)`` and the
``d[J FinvT G]/du`` terms in explicit 3x3 algebra (verified against
finite differences in ``tests/unit/test_fbar_mean_dilatation.py``).

Assembly: ``residual_tangent_sparse`` scatters the cell blocks into a CSR
pattern precomputed in ``build_model`` (node-pair scatter map, one Numba
pass over all cells); ``solve_static`` factorises it with Pardiso / SuperLU.
The dense ``residual_tangent`` is kept as the test oracle (S40).

Scope: untrimmed single-patch boxes (E04 precursor), clamped knots with any
interior multiplicity 1..p per knot, i.e. C^(p-1)..C^0, uniform or mixed
(Coreform E04 uses C^(p-2)). No Kratos dependency.
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
    # H[q, I] = F^{-T} G_qI, reused by V, dJbar and the dV/du term
    H = np.zeros((nq, nen, 3), dtype=np.float64)
    for q in range(nq):
        FinvT = _inv3(F[q]).T
        for I in range(nen):
            for a in range(3):
                h = 0.0
                for m in range(3):
                    h += FinvT[a, m] * G[q, I, m]
                H[q, I, a] = h
    V = np.zeros((nen, 3), dtype=np.float64)
    dJbar = np.zeros((nen, 3), dtype=np.float64)
    for I in range(nen):
        for q in range(nq):
            wJ = W[q] * detJ[q]
            for a in range(3):
                V[I, a] += wJ * H[q, I, a]
                dJbar[I, a] += wJ * H[q, I, a]
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
    # tangent part 1: deviatoric (constant). s_IJ == s_JI bit for bit
    # (products commute exactly, same summation order) -> J >= I only.
    for I in range(nen):
        for J in range(I, nen):
            s = 0.0
            for q in range(nq):
                g = 0.0
                for m in range(3):
                    g += G[q, I, m] * G[q, J, m]
                s += W[q] * mu * g
            for a in range(3):
                K[3 * I + a, 3 * J + a] += s
                if J != I:
                    K[3 * J + a, 3 * I + a] += s
    # tangent part 2: U'' * dJbar (x) V  +  pbar * dV/du
    for I in range(nen):
        for J in range(nen):
            for a in range(3):
                for b in range(3):
                    K[3 * I + a, 3 * J + b] += ppbar * V[I, a] * dJbar[J, b]
    # dV term: J*(outer(H_I,H_J) - outer(H_J,H_I)) -- verified vs FD.
    # Each 3x3 block is antisymmetric with an exactly zero diagonal, block
    # (J,I) is the transpose of block (I,J), and block (I,I) is exactly zero
    # (IEEE: x*y - z*w == -(z*w - x*y), products commute). So only J > I and
    # (a,b) = (0,1), (0,2), (1,2) are formed; the result is bit-identical to
    # the full I, J, a, b loop.
    for I in range(nen):
        r0, r1, r2 = 3 * I, 3 * I + 1, 3 * I + 2
        for J in range(I + 1, nen):
            c0, c1, c2 = 3 * J, 3 * J + 1, 3 * J + 2
            for q in range(nq):
                wJ = W[q] * detJ[q]
                cq = pbar * wJ
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


@njit(cache=True)
def _assemble_csr(G_all: np.ndarray, W_all: np.ndarray, conn: np.ndarray, u: np.ndarray,
                  mu: float, lam: float, indptr: np.ndarray, node_pos: np.ndarray, nnz: int):
    """Global R and CSR data of K (pattern from ``_csr_pattern``).

    Entry (3n+a, 3m+b) of a coupled node pair sits at
    ``indptr[3n+a] + 3*node_pos[c, I, J] + b``. Cells are added in the same
    order as the dense oracle ``residual_tangent``.
    """
    nc, nen = conn.shape
    R = np.zeros(u.shape[0], dtype=np.float64)
    data = np.zeros(nnz, dtype=np.float64)
    uc = np.empty((nen, 3), dtype=np.float64)
    for c in range(nc):
        for I in range(nen):
            n = conn[c, I]
            for a in range(3):
                uc[I, a] = u[3 * n + a]
        Rc, Kc, _ = _cell_residual_tangent(G_all[c], W_all[c], uc, mu, lam)
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


def _csr_pattern(conn: np.ndarray, n_ctrl: int) -> dict:
    """CSR pattern of the 3-DOF-per-node tangent and the cell scatter map.

    Built per node pair (nc * nen^2 entries, 9x less than per DOF pair): nodes
    n, m are coupled iff they share a cell, and then all 3x3 DOF pairs are
    stored. Row 3n+a lists ``3m+b`` for the sorted neighbours m of n, so
    entry (3n+a, 3m+b) is at ``indptr[3n+a] + 3k + b`` with k the rank of m
    among the neighbours of n; ``node_pos[c, I, J]`` holds k for cell c.
    """
    nc, nen = conn.shape
    key = (np.repeat(conn, nen, axis=1).astype(np.int64) * n_ctrl
           + np.tile(conn, (1, nen)).astype(np.int64)).ravel()
    uniq, inv = np.unique(key, return_inverse=True)
    rows, cols = uniq // n_ctrl, uniq % n_ctrl
    deg = np.bincount(rows, minlength=n_ctrl).astype(np.int64)
    nstart = np.zeros(n_ctrl + 1, dtype=np.int64)
    np.cumsum(deg, out=nstart[1:])
    rank = np.arange(uniq.size, dtype=np.int64) - nstart[rows]
    node_pos = rank[inv.ravel()].reshape(nc, nen, nen)
    row_node = np.repeat(np.arange(n_ctrl, dtype=np.int64), 3)
    row_len = 3 * deg[row_node]
    indptr = np.zeros(3 * n_ctrl + 1, dtype=np.int64)
    np.cumsum(row_len, out=indptr[1:])
    nnz = int(indptr[-1])
    nb_cols = (3 * cols[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
    gather = np.repeat(3 * nstart[row_node] - indptr[:-1], row_len) + np.arange(nnz, dtype=np.int64)
    return {"indptr": indptr, "indices": nb_cols[gather],
            "node_pos": np.ascontiguousarray(node_pos), "nnz": nnz}


def _axis_cells(knots: np.ndarray, p: int, xg: np.ndarray) -> list:
    """Per-cell 1D data along one axis for any clamped knot vector.

    Returns ``[(lo, hi, first_cp, N(ng, p+1), dN(ng, p+1)), ...]``. The span
    comes from ``find_span`` at the cell midpoint, so interior knot
    multiplicity m (continuity C^(p-m)) is handled: first_cp = span - p.
    """
    n_cp = len(knots) - p - 1
    br = np.unique(knots)
    out = []
    for c in range(len(br) - 1):
        lo, hi = float(br[c]), float(br[c + 1])
        span = int(find_span(n_cp, p, 0.5 * (lo + hi), knots))
        if knots[span] != lo or knots[span + 1] != hi:
            raise ValueError(f"cell {c}: span {span} does not cover [{lo}, {hi}]")
        ng = len(xg)
        N = np.empty((ng, p + 1))
        dN = np.empty((ng, p + 1))
        for g in range(ng):
            d = ders_basis_funs(span, lo + 0.5 * (xg[g] + 1.0) * (hi - lo), p, 1, knots)
            N[g] = d[0]
            dN[g] = d[1]
        out.append((lo, hi, span - p, N, dN))
    return out


def build_model(box, degree: int, E: float, nu: float,
                n_gauss: int | None = None) -> dict:
    """Precompute per-cell QP data and the CSR pattern.

    *box*: BoxSpline-like (``knots``, ``n_cp``, ``control_points()``), clamped
    knots with any interior multiplicity 1..p per knot (C^(p-1)..C^0, uniform
    or mixed); ``_axis_cells`` maps each cell to its span.
    """
    p = int(degree)
    ng = p + 1 if n_gauss is None else int(n_gauss)
    xg, wg = np.polynomial.legendre.leggauss(ng)
    knots = [np.ascontiguousarray(k, dtype=np.float64) for k in box.knots]
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    for k, n in zip(knots, (nu_, nv_, nw_)):
        if len(k) - p - 1 != n:
            raise ValueError("box.n_cp does not match its knots for this degree")
    cp = np.asarray(box.control_points(), dtype=float)
    n_ctrl = cp.shape[0]
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    tx, ty, tz = (_axis_cells(k, p, xg) for k in knots)
    nq, ncn = ng ** 3, (p + 1) ** 3
    wq = np.einsum("z,y,x->zyx", wg, wg, wg).reshape(nq)
    ax = np.arange(p + 1)
    cells: list = []
    for zl, zh, fz, Nw, dNw in tz:
        for yl, yh, fy, Nv, dNv in ty:
            for xl, xh, fx, Nu, dNu in tx:
                # local node ii = ai + (p+1)*(aj + (p+1)*ak); q = gx + ng*(gy + ng*gz)
                idx = ((fx + ax)[None, None, :]
                       + nu_ * ((fy + ax)[None, :, None] + nv_ * (fz + ax)[:, None, None])).reshape(ncn)
                d3 = np.empty((nq, ncn, 3))
                d3[:, :, 0] = np.einsum("zk,yj,xi->zyxkji", Nw, Nv, dNu).reshape(nq, ncn)
                d3[:, :, 1] = np.einsum("zk,yj,xi->zyxkji", Nw, dNv, Nu).reshape(nq, ncn)
                d3[:, :, 2] = np.einsum("zk,yj,xi->zyxkji", dNw, Nv, Nu).reshape(nq, ncn)
                Jg = np.einsum("qnc,ni->qic", d3, cp[idx])  # dX_i/dxi_c
                G = np.einsum("qnc,qcd->qnd", d3, np.linalg.inv(Jg))
                W = wq * (0.125 * (xh - xl) * (yh - yl) * (zh - zl)) * np.abs(np.linalg.det(Jg))
                cells.append({"nodes": idx.astype(np.int64), "G": G, "W": W})
    conn = np.ascontiguousarray(np.stack([c["nodes"] for c in cells]))
    G_all = np.ascontiguousarray(np.stack([c["G"] for c in cells]))
    W_all = np.ascontiguousarray(np.stack([c["W"] for c in cells]))
    for i, cell in enumerate(cells):  # views: one copy of the QP data, not two
        cell["G"], cell["W"] = G_all[i], W_all[i]
    return {"cells": cells, "n_ctrl": n_ctrl, "n_dof": 3 * n_ctrl,
            "mu": float(mu), "lam": float(lam), "degree": p,
            "G_all": G_all, "W_all": W_all, "conn": conn,
            "csr": _csr_pattern(conn, n_ctrl)}


def residual_tangent(model: dict, u: np.ndarray):
    """Assemble global R, K (dense). Oracle for the sparse path; small models only."""
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


def residual_tangent_sparse(model: dict, u: np.ndarray):
    """Assemble global R (ndarray) and K (scipy CSR, sorted indices)."""
    import scipy.sparse as sp

    ndof = int(model["n_dof"])
    s = model["csr"]
    R, data = _assemble_csr(model["G_all"], model["W_all"], model["conn"],
                            np.ascontiguousarray(u, dtype=np.float64).reshape(ndof),
                            float(model["mu"]), float(model["lam"]),
                            s["indptr"], s["node_pos"], s["nnz"])
    K = sp.csr_matrix((data, s["indices"], s["indptr"]), shape=(ndof, ndof))
    K.has_sorted_indices = True
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
                 history: list | None = None,
                 linear_solver: str = "auto") -> np.ndarray:
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

    ``linear_solver``: ``"pardiso"`` (pypardiso, own handle, freed on exit),
    ``"superlu"`` (scipy), ``"auto"`` (pardiso if importable, else superlu;
    both are direct solvers) on the CSR tangent of ``residual_tangent_sparse``,
    or ``"dense"`` (dense oracle ``residual_tangent`` + ``np.linalg.solve``,
    the pre-2026-11-19 path; small models only).
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
    lam, dlam = 0.0, d_nom

    def _increment(lam_new: float) -> tuple[bool, int]:
        target = f * lam_new
        for k, v in prescribed.items():
            u[k] = v * lam_new
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
    finally:
        if pardiso is not None:
            pardiso.free_memory(everything=True)
    return u
