"""T-P2.3 acceptance: TᵀKT vs reference assembler, c=0 + chain rules (Kratos-free).

Oracle: ``wht_iga.reference.assemble_elastic`` (T-P0.5 assembler) with a
hand-rolled Lagrange-multiplier (KKT) elimination coded independently in the
test — two different algebraic paths to the same constrained solution.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.operators import compose as C
from wht_iga.operators import numbering as N
from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline
from wht_iga.stabilize import extension as E
from wht_iga.stabilize import support as S


def _coef1d(knots, p, f):
    """B-spline coefficients of f by Greville collocation (independent oracle)."""
    from wht_iga.grid.bspline import basis_1d
    from wht_iga.grid.greville import greville

    g = greville(knots, p)
    M = np.zeros((len(g), len(g)))
    for r, x in enumerate(g):
        span, N_ = basis_1d(knots, p, float(x), 0)
        M[r, span - p:span + 1] = np.atleast_2d(N_)[0]
    return np.linalg.solve(M, f(g))


def _v01_reduced_case():
    """Unit cube, p=2, 2³ cells; slab i==0 cut-like (eta 0.05); affine patch."""
    from wht_iga.grid.knots import uniform_clamped

    p = 2
    box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2, 2, 2), p)
    knots = [uniform_clamped(0.0, 1.0, 2, p) for _ in range(3)]
    cv = np.full(8, 0.125)
    eta = np.ones(8)
    eta[[0, 2, 4, 6]] = 0.05  # cells (0, j, k): order i + 2*(j + 2k)
    cls = S.classify_basis(box.n_cp, knots, eta * cv, cv, theta=0.2)
    op = E.build_extension(box.n_cp, cls.labels, p=p, knots=knots)
    assert len(op.slave_ids) > 0
    # slaves are exactly the a == 0 face (their support is the cut slab)
    nu = box.n_cp[0]
    assert set(op.slave_ids.tolist()) == {nu * (b + nu * c)
                                           for b in range(4)
                                           for c in range(4)}
    K = assemble_stiffness(box, 1.0, 0.0)
    # Manufactured field: degree-p polynomial in PHYSICAL coordinates
    # (1, x, y, z, xy, z^2). Coefficients: Greville points for linear terms,
    # tensor products for xy, Greville collocation for z^2 (independent of
    # the extension code). Index-space fields are NOT polynomials near a
    # clamped boundary (V01-ext defect, 2026-11-14).
    nu = box.n_cp[0]
    X = box.control_points()
    cz2 = _coef1d(knots[2], p, lambda z: z ** 2)
    c_idx = np.arange(box.n_ctrl) // (nu * nu)
    rng = np.random.default_rng(21)
    P = rng.uniform(-1, 1, size=(3, 6))
    mono = np.column_stack([np.ones(len(X)), X[:, 0], X[:, 1],
                            X[:, 2], X[:, 0] * X[:, 1], cz2[c_idx]])
    u_exact = (mono @ P.T).ravel()
    f = K @ u_exact
    # fixed face a == 3 (x = 1): all masters (no fixed slaves allowed)
    fixed_cps = [3 + nu * (b + nu * c) for b in range(4) for c in range(4)]
    assert not set(fixed_cps) & set(op.slave_ids.tolist())
    fixed = np.sort(np.concatenate(
        [np.array(fixed_cps) * 3 + d for d in range(3)]))
    return box, op, K, f, u_exact, fixed


def test_tkt_agreement_vs_kkt_1e10():
    box, op, K, f, u_exact, fixed = _v01_reduced_case()
    Td = C.expand_dof_operator(C.compose_thb(op.E), n_comp=3)
    assert Td.shape == (box.n_dofs, len(op.master_ids) * 3)
    assert isinstance(Td, sp.csr_matrix)
    master_dofs = (op.master_ids[:, None] * 3 + np.arange(3)).ravel()
    fixed_set = set(fixed.tolist())
    var = np.array([c for c, md in enumerate(master_dofs)
                    if int(md) not in fixed_set], dtype=np.int64)
    con = np.array([c for c, md in enumerate(master_dofs)
                    if int(md) in fixed_set], dtype=np.int64)
    free = np.setdiff1d(np.arange(box.n_dofs), fixed)
    Kff = K[free][:, free]
    rhs_f = f[free] - K[free][:, fixed] @ u_exact[fixed]
    Tfm_full = Td[free, :]
    Tfm = Tfm_full[:, var]
    cvec = (Tfm_full[:, con] @ u_exact[master_dofs[con]]
            if len(con) else np.zeros(Tfm_full.shape[0]))
    Kred = C.apply_TKT(Kff, Tfm)
    fred = C.apply_Tt_vec(rhs_f - Kff @ cvec, Tfm)
    um = C.solve_reduced(Kred, fred)
    u_a = np.zeros(box.n_dofs)
    u_a[fixed] = u_exact[fixed]
    u_a[free] = Tfm @ um + cvec
    # Path B: textbook KKT on free DOFs (independent elimination).
    free_set = set(free.tolist())
    slave_dofs = np.array([int(3 * g + d) for g in op.slave_ids
                           for d in range(3) if int(3 * g + d) in free_set])
    assert len(slave_dofs) == 3 * len(op.slave_ids)  # no fixed slaves
    free_list = free.tolist()
    pos = {d: k for k, d in enumerate(free_list)}
    n_f, n_c = len(free_list), len(slave_dofs)
    Cmat = sp.lil_matrix((n_c, n_f))
    rhs_c = np.zeros(n_c)
    Tm = Td[free, :].tocsr()
    for r, sd in enumerate(slave_dofs):
        Cmat[r, pos[int(sd)]] = 1.0
        row = Tm.getrow(pos[int(sd)])
        for cc, w in zip(row.indices, row.data):
            md = int(master_dofs[cc])
            if md in pos:
                Cmat[r, pos[md]] -= w
            else:
                rhs_c[r] += w * u_exact[md]
    KKT = sp.bmat([[Kff, Cmat.T], [Cmat, None]], format="csc")
    rhs = np.concatenate([rhs_f, rhs_c])
    sol = spla.spsolve(KKT, rhs)
    u_b = np.zeros(box.n_dofs)
    u_b[fixed] = u_exact[fixed]
    u_b[free] = sol[:n_f]
    np.testing.assert_allclose(u_a, u_b, rtol=0, atol=1e-10)
    np.testing.assert_allclose(u_a, u_exact, rtol=0, atol=1e-8)


def test_compose_identity_and_thb_path():
    _box, op, K, _f, _u, _fixed = _v01_reduced_case()
    T = C.compose_thb(op.E, None)  # Phase-2 path: C_THB = I
    assert sp.issparse(T) and isinstance(T, sp.csr_matrix)
    assert T.nnz == op.E.nnz
    np.testing.assert_allclose((T - op.E).toarray(), 0.0, atol=0)
    T2 = C.compose_thb(op.E, 2.0 * sp.identity(op.E.shape[0], format="csr"))
    np.testing.assert_allclose((T2 - 2 * op.E).toarray(), 0.0, atol=0)
    with pytest.raises(ValueError):
        C.compose_thb(op.E, sp.identity(3, format="csr"))
    with pytest.raises(TypeError):
        C.apply_TKT(np.eye(3), op.E)
    with pytest.raises(TypeError):
        C.compose_thb(np.eye(op.E.shape[0]))


def test_extension_matrix_assembles():
    _box, op, _K, _f, _u, _fixed = _v01_reduced_case()
    n_full = int(op.E.shape[0])
    sids = op.slave_ids
    Srows = op.E[sids, :]
    E2 = C.extension_matrix(n_full, op.master_ids, sids, Srows)
    assert E2.nnz == op.E.nnz
    np.testing.assert_allclose((E2 - op.E).toarray(), 0.0, atol=0)
    with pytest.raises(TypeError):
        C.extension_matrix(n_full, op.master_ids, sids, np.eye(2))


def test_c_zero_refusal():
    N.require_zero_constant(None)
    N.require_zero_constant(np.zeros(4))
    for bad in (np.array([0.0, 1e-9]), np.array([1.0]),
                np.array([[0.0], [2e-12]])):
        with pytest.raises(ValueError):
            N.require_zero_constant(bad)
    with pytest.raises(ValueError):
        N.require_zero_constant(np.array([0.0, -3e-11]), tol=1e-12)


def test_chain_rejection():
    N.validate_no_chains(np.array([1, 2]), [np.array([3, 4])])
    with pytest.raises(ValueError):
        N.validate_no_chains(np.array([1, 2]), [np.array([2, 5])])
    with pytest.raises(ValueError):
        N.build_numbering(6, np.array([4, 1]), np.array([1, 2]))
    num = N.build_numbering(6, np.array([4, 1]), np.array([5, 0]))
    assert num.n_master == 2 and list(num.master_ids) == [1, 4]


def test_scalar_specs_and_batch_timing():
    _box, op, _K, _f, _u, _fixed = _v01_reduced_case()
    num = N.build_numbering(int(op.E.shape[0]), op.master_ids, op.slave_ids)
    specs = N.make_scalar_specs(num, op.E, op.slave_ids, n_comp=3)
    assert len(specs) == 3 * len(op.slave_ids)
    for spc in specs:
        assert spc.constant is not None and float(spc.constant[0]) == 0.0
        assert len(spc.master_dofs) == len(spc.row) > 0
    with pytest.raises(ValueError):
        N.make_scalar_specs(num, op.E, op.slave_ids,
                            constant=np.array([0.0, 1e-6]))
    # 1e5-spec batch generation timing (acceptance: <= 30 s).
    rng = np.random.default_rng(7)
    n_s, n_m = 33334, 36000
    n_full = n_s + n_m  # disjoint slave/master id ranges (no chains)
    sids = np.arange(n_s)
    mids = np.arange(n_s, n_s + n_m)
    mid = np.array([rng.choice(mids, size=8, replace=False)
                    for _ in range(n_s)])
    w = rng.uniform(-1, 1, size=(n_s, 8))
    w /= np.abs(w).sum(axis=1, keepdims=True)
    rows = np.repeat(np.arange(n_s), 8)
    big = sp.csr_matrix((w.ravel(), (rows, mid.ravel() - n_s)),
                        shape=(n_s, n_m))
    Tbig = sp.vstack(
        [big,
         sp.coo_matrix(([], ([], [])), shape=(n_m, n_m))],
        format="csr")
    # NOTE: Tbig rows beyond n_s are zero (inactive-style); specs use big rows.
    num2 = N.build_numbering(n_full, mids, sids)
    t0 = time.perf_counter()
    specs2 = N.make_scalar_specs(num2, Tbig, sids, n_comp=3)
    dt = time.perf_counter() - t0
    assert len(specs2) == 3 * n_s
    print(f"1e5 scalar-spec generation: {dt:.2f} s (acceptance <= 30 s)")
    assert dt < 30.0
