"""S05: LinearMasterSlaveConstraint — accuracy vs T^T K T, 1e5 creation, DofSet.

Reference: scipy-assembled reduced system (same K as S04 + same T), exact
slave recovery, analytic tie constants.
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
from wht_iga.reference import assemble_elastic as ref
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

VARS3 = None  # set per test from KM


def _model(KM, name: str):
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def _props(KM, SM, mp, E=2.5, nu=0.3):
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    return props


def _fix(KM, mp, ids, values):
    for node in mp.Nodes:
        if node.Id in ids:
            for var, d in zip(
                (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z), range(3)
            ):
                node.Fix(var)
                node.SetSolutionStepValue(var, float(values[node.Id][d]))


def test_constrained_solution_matches_reduced_system():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    E, nu = 2.5, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(3.0, 2.0, 2.0), n_cells=(3, 2, 2), degree=1)
    assert box.n_dofs == 108
    _, mp = _model(KM, "S05A")
    vol, _ = make_nurbs_volume(KM, mp, box)
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp, E, nu))

    rng = np.random.default_rng(11)
    nu_, nv, _ = box.n_cp
    interior = [I for I in range(box.n_ctrl) if 0 < I % nu_ < nu_ - 1]
    slave_ctrl = rng.choice(interior, size=7, replace=False)
    slave_dofs = sorted([3 * int(I) + d for I in slave_ctrl for d in range(3)])[:20]
    slave_set = set(slave_dofs)
    fixed_ctrl = [I for I in range(box.n_ctrl) if I % nu_ == 0]
    fixed_dofs = sorted([3 * int(I) + d for I in fixed_ctrl for d in range(3)])
    driven_ctrl = [I for I in range(box.n_ctrl) if I % nu_ == nu_ - 1]
    driven_dofs = sorted([3 * int(I) + d for I in driven_ctrl for d in range(3)])
    fixed_set = set(fixed_dofs)
    free_pool = [d for d in range(box.n_dofs) if d not in slave_set and d not in fixed_set]

    # T: 20 slaves x 27 masters, masters drawn from non-slave DOFs (may be fixed).
    n_s, n_m = len(slave_dofs), 27
    Tm = np.zeros((n_s, n_m))
    master_dofs: list[int] = []
    per_row: list[tuple[np.ndarray, np.ndarray]] = []
    for r in range(n_s):
        masters = rng.choice(free_pool + fixed_dofs, size=n_m, replace=False)
        w = rng.uniform(-1.0, 1.0, size=n_m)
        w /= np.sum(np.abs(w)) + 0.5
        Tm[r] = w
        per_row.append((np.array([int(v) for v in masters]), w))
        master_dofs.extend(int(v) for v in masters)
    master_dofs = sorted(set(master_dofs))
    col_of = {d: c for c, d in enumerate(master_dofs)}
    Tall = np.zeros((n_s, len(master_dofs)))
    for r in range(n_s):
        for d, wi in zip(*per_row[r]):
            Tall[r, col_of[int(d)]] = wi

    g = np.zeros(box.n_dofs)
    g[driven_dofs[0::3]] = 0.1
    g[driven_dofs[1::3]] = -0.05
    g[driven_dofs[2::3]] = 0.02
    prescribed = dict(zip(fixed_dofs, g[fixed_dofs]))
    prescribed.update(zip(driven_dofs, g[driven_dofs]))

    # Kratos side: fix prescribed DOFs, add ties (constant absorbs fixed masters).
    dof_of = {}
    for node in mp.Nodes:
        I = node.Id - 1
        for var, d in zip((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z), range(3)):
            dof_of[3 * I + d] = node.GetDof(var)
    fixed_nodes = {1 + d // 3 for d in prescribed}
    node_vals = {nid: [0.0, 0.0, 0.0] for nid in fixed_nodes}
    for d, v in prescribed.items():
        node_vals[1 + d // 3][d % 3] = v
    _fix(KM, mp, fixed_nodes, node_vals)
    cid = 1
    for r, sd in enumerate(slave_dofs):
        # Dirichlet rule: keep prescribed masters in the relation, c = 0.
        # Kratos folds prescribed values itself; pre-folding double-counts.
        row_masters = [d for d in master_dofs if Tall[r, col_of[d]] != 0.0]
        add_linear_tie(
            KM, mp, cid,
            [dof_of[sd]],
            [dof_of[d] for d in row_masters],
            Tall[r, [col_of[d] for d in row_masters]].reshape(1, -1),
            np.zeros(1),
        )
        cid += 1
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    strat = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.Solve()
    uk = np.zeros(box.n_dofs)
    for node in mp.Nodes:
        I = node.Id - 1
        uk[3 * I] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
        uk[3 * I + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
        uk[3 * I + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)

    # scipy side: T maps master DOFs (+ prescribed) to full vector.
    K = ref.assemble_stiffness(box, E, nu).tocsr()
    master_list = master_dofs
    mcol = {d: c for c, d in enumerate(master_list)}
    free_masters = [d for d in master_list if d not in prescribed]
    fcol = {d: c for c, d in enumerate(free_masters)}
    T = sp.lil_matrix((box.n_dofs, len(free_masters)))
    rhs_const = np.zeros(box.n_dofs)
    for d in range(box.n_dofs):
        if d in prescribed:
            rhs_const[d] = prescribed[d]
        elif d in slave_set:
            r = slave_dofs.index(d)
            for md in master_list:
                w = Tall[r, mcol[md]]
                if w != 0.0:
                    if md in prescribed:
                        rhs_const[d] += w * prescribed[md]
                    else:
                        T[d, fcol[md]] = w
        else:
            T[d, fcol[d]] = 1.0
    T = T.tocsr()
    Kred = (T.T @ K @ T).tocsr()
    fred = -T.T @ (K @ rhs_const)
    from scipy.sparse.linalg import spsolve

    um = spsolve(Kred.tocsc(), fred)
    uref = T @ um + rhs_const
    err = float(np.max(np.abs(uk - uref)))
    rec = max(abs(uk[sd] - (Tall[r] @ np.array([uk[d] for d in master_list]))) for r, sd in enumerate(slave_dofs))
    print("S05-A max err:", err, "slave recovery:", rec)
    assert err < 1e-10, err
    # slave recovery: constrained DOFs follow their rows exactly
    assert rec < 1e-12, rec


def test_master_only_nodes_in_dofset():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    _, mp = _model(KM, "S05C")
    lone = mp.CreateNewNode(1000, 5.0, 5.0, 5.0)  # before AddDof so it joins (K14)
    vol, _ = make_nurbs_volume(KM, mp, box)
    from wht_iga.kratos_backend.geometry import cell_gauss_points

    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, 0, 0, 0))
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qps[0], props)
    mp.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1,
        [lone.GetDof(KM.DISPLACEMENT_X), lone.GetDof(KM.DISPLACEMENT_Y), lone.GetDof(KM.DISPLACEMENT_Z)],
        [mp.GetNode(8).GetDof(KM.DISPLACEMENT_X), mp.GetNode(8).GetDof(KM.DISPLACEMENT_Y), mp.GetNode(8).GetDof(KM.DISPLACEMENT_Z)],
        KM.Matrix(3, 3), KM.Vector(3),
    )
    mp.SetBufferSize(1)
    strat = KM.ResidualBasedLinearStrategy(
        mp, KM.ResidualBasedIncrementalUpdateStaticScheme(),
        KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver()),
        False, False, False, False,
    )
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.Solve()
    eqs = [lone.GetDof(v).EquationId for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)]
    assert all(e >= 0 for e in eqs), eqs
    assert eqs == [24, 25, 26], eqs  # 24 element DOFs + 3 master DOFs


@pytest.mark.slow
def test_create_1e5_constraints_and_solve():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    import KratosMultiphysics.LinearSolversApplication as LA

    _, mp = _model(KM, "S05B")
    masters = [mp.CreateNewNode(1 + i, float(i), 0.0, 0.0) for i in range(27)]
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    stiff = KM.Array3()
    stiff[0], stiff[1], stiff[2] = 1.0, 1.0, 1.0
    props = mp.GetProperties()[1]
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    for i, mn in enumerate(masters):
        el = mp.CreateNewElement("NodalConcentratedElement3D1N", 1 + i, [mn.Id], props)
        el.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, stiff)
    master_dofs = []
    for mn in masters:
        master_dofs += [mn.GetDof(KM.DISPLACEMENT_X), mn.GetDof(KM.DISPLACEMENT_Y), mn.GetDof(KM.DISPLACEMENT_Z)]
    rng = np.random.default_rng(5)
    rel = KM.Matrix(3, 81)
    for i in range(3):
        for j in range(81):
            rel[i, j] = 0.0
    W = rng.uniform(-1.0, 1.0, size=(3, 81))
    W /= np.sum(np.abs(W), axis=1, keepdims=True) + 1.0
    for i in range(3):
        for j in range(81):
            rel[i, j] = float(W[i, j])
    consts = rng.uniform(-0.5, 0.5, size=(100_000, 3))
    t0 = time.perf_counter()
    for n in range(100_000):
        sn = mp.CreateNewNode(1000 + n, float(n) * 1e-4, 1.0, 1.0)
        sn.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
        sn.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
        sn.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)
        c = KM.Vector(3)
        c[0], c[1], c[2] = (float(consts[n, 0]), float(consts[n, 1]), float(consts[n, 2]))
        mp.CreateNewMasterSlaveConstraint(
            "LinearMasterSlaveConstraint", 1 + n, master_dofs,
            [sn.GetDof(KM.DISPLACEMENT_X), sn.GetDof(KM.DISPLACEMENT_Y), sn.GetDof(KM.DISPLACEMENT_Z)],
            rel, c,
        )
    t_create = time.perf_counter() - t0
    assert t_create < 30.0, t_create
    print(f"1e5 constraints created in {t_create:.1f} s")
    mp.SetBufferSize(1)
    t0 = time.perf_counter()
    strat = KM.ResidualBasedLinearStrategy(
        mp, KM.ResidualBasedIncrementalUpdateStaticScheme(),
        KM.ResidualBasedBlockBuilderAndSolver(LA.PardisoLUSolver()),
        False, False, False, False,
    )
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.Solve()
    t_solve = time.perf_counter() - t0
    print(f"assemble+solve in {t_solve:.1f} s")
    worst = 0.0
    for n in range(0, 100_000, 997):
        sn = mp.GetNode(1000 + n)
        got = np.array([sn.GetSolutionStepValue(KM.DISPLACEMENT_X),
                        sn.GetSolutionStepValue(KM.DISPLACEMENT_Y),
                        sn.GetSolutionStepValue(KM.DISPLACEMENT_Z)])
        worst = max(worst, float(np.max(np.abs(got - consts[n]))))
    assert worst < 1e-9, worst
