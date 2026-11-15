"""P2-STAB Kratos spike K-A: TᵀKT Kratos-vs-reference + 1e5 creation timing.

Kratos-side analog of S05 (Track A): the STAB extension operator supplies a
real T (cut-slab V01-reduced case); Kratos assembles K with
``LinearMasterSlaveConstraint`` ties and the test compares against the
TᵀKT-projected reference matrix (acceptance ≤ 1e-10). Creation timing is
reported against the S05 2.7 s baseline (acceptance: 1e5 ≤ 30 s).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.operators import compose as C
from wht_iga.operators import numbering as N
from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline
from wht_iga.stabilize import extension as E
from wht_iga.stabilize import support as S

pytestmark = pytest.mark.kratos


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


def _model(KM, name: str):
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def test_kratos_tkt_matches_reference_1e10():
    """Tied-model Kratos solution vs T-reduced scipy solution (S05 pattern).

    With ties present, GetSystemMatrix is the constrained BlockBuilder system
    (slave rows retained/diagonalized — Track A finding #2), NOT raw K, so a
    raw-K comparison on the tied model is invalid. The acceptance runs in two
    halves: (a) twin model WITHOUT ties: extracted K == reference K
    (<= 1e-10); (b) tied model: Kratos solution == S05-style T-reduced scipy
    solution (<= 1e-10) with slave recovery (<= 1e-12).
    """
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    from wht_iga.kratos_backend.elements import solve_linear
    Emod, nu = 2.5, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(2, 2, 2), degree=2)
    knots = [np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]) for _ in range(3)]
    cv = np.full(8, 0.125)
    eta = np.ones(8)
    eta[[0, 2, 4, 6]] = 0.05
    cls = S.classify_basis(box.n_cp, knots, eta * cv, cv, theta=0.2)
    op = E.build_extension(box.n_cp, cls.labels, p=2, knots=knots)
    assert len(op.slave_ids) == 16
    Td = C.expand_dof_operator(C.compose_thb(op.E), n_comp=3).tocsr()
    master_dofs = (op.master_ids[:, None] * 3 + np.arange(3)).ravel()
    # Manufactured physical-coordinate quadratic (1, x, y, z, xy, z^2);
    # extension reproduces it exactly on the true knots (V01-ext fix).
    nu_ = box.n_cp[0]
    rng = np.random.default_rng(21)
    Pm = rng.uniform(-1, 1, size=(3, 6))
    X = box.control_points()
    cz2 = _coef1d(knots[2], 2, lambda z: z ** 2)
    c_idx = np.arange(box.n_ctrl) // (nu_ * nu_)
    mono = np.column_stack([np.ones(len(X)), X[:, 0], X[:, 1],
                            X[:, 2], X[:, 0] * X[:, 1], cz2[c_idx]])
    u_manu = (mono @ Pm.T).ravel()
    fixed_cps = [3 + nu_ * (b + nu_ * c) for b in range(4)
                 for c in range(4)]  # x=1 face: all masters
    prescribed = {int(3 * I + d): float(u_manu[3 * I + d])
                  for I in fixed_cps for d in range(3)}

    def _build(name: str, ties: bool):
        # NOTE: the NurbsVolume keepalive MUST stay referenced (Phase-1
        # lifetime rule) — dropping it dangling-frees the QP geometries
        # (observed: element cond=inf at Build).
        _, mp = _model(KM, name)
        vol, _ = make_nurbs_volume(KM, mp, box)
        cl = SM.LinearElastic3DLaw()
        props = mp.GetProperties()[1]
        props.SetValue(KM.YOUNG_MODULUS, Emod)
        props.SetValue(KM.POISSON_RATIO, nu)
        props.SetValue(KM.CONSTITUTIVE_LAW, cl)
        attach_cell_elements(KM, mp, vol, box, props)
        dof_of = {}
        for node in mp.Nodes:
            I = node.Id - 1
            for var, d in zip((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                               KM.DISPLACEMENT_Z), range(3)):
                dof_of[3 * I + d] = node.GetDof(var)
        for d, v in prescribed.items():  # strong Dirichlet (S05 pattern)
            node = mp.GetNode(d // 3 + 1)
            var = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                   KM.DISPLACEMENT_Z)[d % 3]
            node.Fix(var)
            node.SetSolutionStepValue(var, v)
        if ties:
            cid = 1
            for sd in sorted(int(3 * g + d) for g in op.slave_ids
                             for d in range(3)):
                r = Td.getrow(sd)
                mdofs = [int(master_dofs[c]) for c in r.indices]
                N.require_zero_constant(np.zeros(1))
                add_linear_tie(KM, mp, cid, [dof_of[sd]],
                               [dof_of[m] for m in mdofs],
                               r.data.reshape(1, -1), np.zeros(1))
                cid += 1
        return mp, vol

    # Half (a): twin model without ties — raw K agreement on free DOFs
    # (Dirichlet-fixed rows are folded to 1 by the builder — Track A #6).
    mp0, _keep0 = _build("P2STAB_A0", ties=False)
    Kk = system_matrix_csr(solve_linear(KM, mp0))
    perm = equation_permutation(KM, mp0, box)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    Kref = assemble_stiffness(box, Emod, nu)
    free = np.setdiff1d(np.arange(box.n_dofs),
                        np.array(sorted(prescribed)))
    np.testing.assert_allclose(Kk[inv][:, inv][free][:, free].toarray(),
                               Kref[free][:, free].toarray(),
                               rtol=0, atol=1e-10)
    # Half (b): tied model — solution agreement (S05-A pattern).
    mp1, _keep1 = _build("P2STAB_A1", ties=True)
    solve_linear(KM, mp1)
    uk = np.zeros(box.n_dofs)
    for node in mp1.Nodes:
        I = node.Id - 1
        uk[3 * I] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
        uk[3 * I + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
        uk[3 * I + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)
    master_list = master_dofs.tolist()
    fixed_set = set(prescribed)
    free_masters = [d for d in master_list if d not in fixed_set]
    fcol = {d: c for c, d in enumerate(free_masters)}
    T = sp.lil_matrix((box.n_dofs, len(free_masters)))
    rhs_const = np.zeros(box.n_dofs)
    slave_set = {int(3 * g + d) for g in op.slave_ids for d in range(3)}
    for d in range(box.n_dofs):
        if d in fixed_set:
            rhs_const[d] = prescribed[d]
        elif d in slave_set:
            r = Td.getrow(d)
            for cc, w in zip(r.indices, r.data):
                md = int(master_dofs[cc])
                if md in fixed_set:
                    rhs_const[d] += w * prescribed[md]
                else:
                    T[d, fcol[md]] = w
        else:
            T[d, fcol[d]] = 1.0
    T = T.tocsr()
    Kred = (T.T @ Kref @ T).tocsr()
    fred = -T.T @ (Kref @ rhs_const)
    from scipy.sparse.linalg import spsolve
    um = spsolve(Kred.tocsc(), fred)
    uref = T @ um + rhs_const
    err = float(np.max(np.abs(uk - uref)))
    rec = max(abs(uk[sd] - float(np.asarray(Td.getrow(sd)
                                              @ uk[master_dofs]).ravel()[0]))
              for sd in sorted(slave_set))
    print(f"K-A: tied-solution err={err:.3e}, slave recovery={rec:.3e}")
    assert err < 1e-10, err
    assert rec < 1e-12, rec


def test_kratos_1e5_creation_timing():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(32.0, 32.0, 32.0),
                    n_cells=(32, 32, 32), degree=1)  # 107811 DOFs >= 1e5
    _, mp = _model(KM, "P2STAB_B")
    vol, _ = make_nurbs_volume(KM, mp, box)
    assert mp.NumberOfNodes() == box.n_ctrl
    dof_of = {}
    for node in mp.Nodes:
        I = node.Id - 1
        for var, d in zip((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                            KM.DISPLACEMENT_Z), range(3)):
            dof_of[3 * I + d] = node.GetDof(var)
    rng = np.random.default_rng(7)
    n_spec = 100000
    slaves = rng.choice(box.n_dofs, size=n_spec, replace=False)
    t0 = time.perf_counter()
    for cid, sd in enumerate(slaves, start=1):
        md = int(rng.integers(0, box.n_dofs))
        if md == int(sd):
            md = (md + 1) % box.n_dofs
        add_linear_tie(KM, mp, cid, [dof_of[int(sd)]], [dof_of[md]],
                       np.array([[1.0]]), np.zeros(1))
    dt = time.perf_counter() - t0
    print(f"K-B: 1e5 LinearMasterSlaveConstraint creation: {dt:.2f} s "
          f"(S05 baseline 2.7 s, acceptance <= 30 s)")
    assert dt < 30.0
