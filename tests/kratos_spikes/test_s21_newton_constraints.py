"""S21: LinearMasterSlaveConstraint (c=0 rule) + ResidualBasedNewtonRaphsonStrategy + TotalLagrangian.

Question: Do slave DOF recovery, deformed coordinate update (MoveMeshFlag=True),
and quadratic Newton convergence work properly when LinearMasterSlaveConstraint (c=0)
is combined with TotalLagrangianElement3D4N?

Reference: Pure-Python Total Lagrangian (St. Venant-Kirchhoff) assembler with
explicit T^T K_t T Newton-Raphson reduction (target tolerance: 1e-8).
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.assemble_elastic import gauss_legendre, isotropic_D
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _assemble_tl_st_venant(
    box: BoxSpline, E: float, nu: float, u: np.ndarray, n_gauss: int | None = None
) -> tuple[np.ndarray, sp.csr_matrix]:
    """Pure-Python Total Lagrangian (St. Venant-Kirchhoff) internal force & tangent stiffness."""
    ng = box.degree + 1 if n_gauss is None else n_gauss
    xg, wg = gauss_legendre(ng)
    tabs = [box.axis_tables(d, xg) for d in range(3)]
    (fu, Nu, dNu, _, ju), (fv, Nv, dNv, _, jv), (fw, Nw, dNw, _, jw) = tabs
    p = box.degree
    nloc = (p + 1) ** 3
    ndof_e = 3 * nloc
    nu_cp, nv_cp, _ = box.n_cp
    D = isotropic_D(E, nu)

    f_int = np.zeros(box.n_dofs)
    ncu, ncv, ncw = box.n_cells
    n_cells = ncu * ncv * ncw
    rows = np.empty(n_cells * ndof_e * ndof_e, dtype=np.int64)
    cols = np.empty_like(rows)
    vals = np.empty(n_cells * ndof_e * ndof_e)
    pos = 0

    eye3 = np.eye(3)
    for k in range(ncw):
        for j in range(ncv):
            for i in range(ncu):
                dofs = np.empty(ndof_e, dtype=np.int64)
                for c in range(p + 1):
                    for b in range(p + 1):
                        for a in range(p + 1):
                            loc = a + (p + 1) * (b + (p + 1) * c)
                            I = (fu[i] + a) + nu_cp * ((fv[j] + b) + nv_cp * (fw[k] + c))
                            dofs[3 * loc : 3 * loc + 3] = [3 * I, 3 * I + 1, 3 * I + 2]
                ue = u[dofs].reshape(nloc, 3)
                fe = np.zeros(ndof_e)
                Ke = np.zeros((ndof_e, ndof_e))

                for gz in range(ng):
                    for gy in range(ng):
                        for gx in range(ng):
                            wdet = wg[gx] * wg[gy] * wg[gz] * ju[i] * jv[j] * jw[k]
                            grad = np.empty((nloc, 3))
                            for c in range(p + 1):
                                for b in range(p + 1):
                                    for a in range(p + 1):
                                        loc = a + (p + 1) * (b + (p + 1) * c)
                                        grad[loc, 0] = dNu[i, gx, a] * Nv[j, gy, b] * Nw[k, gz, c]
                                        grad[loc, 1] = Nu[i, gx, a] * dNv[j, gy, b] * Nw[k, gz, c]
                                        grad[loc, 2] = Nu[i, gx, a] * Nv[j, gy, b] * dNw[k, gz, c]

                            # Deformation gradient F = I + sum_a u_a x grad_a
                            F = eye3 + ue.T @ grad  # (3, 3), F[k, i] = dx_k / dX_i
                            C_cauchy = F.T @ F
                            E_gl = 0.5 * (C_cauchy - eye3)
                            E_voigt = np.array(
                                [
                                    E_gl[0, 0],
                                    E_gl[1, 1],
                                    E_gl[2, 2],
                                    2.0 * E_gl[0, 1],
                                    2.0 * E_gl[1, 2],
                                    2.0 * E_gl[0, 2],
                                ]
                            )
                            S_voigt = D @ E_voigt
                            S_mat = np.array(
                                [
                                    [S_voigt[0], S_voigt[3], S_voigt[5]],
                                    [S_voigt[3], S_voigt[1], S_voigt[4]],
                                    [S_voigt[5], S_voigt[4], S_voigt[2]],
                                ]
                            )

                            BNL = np.zeros((6, ndof_e))
                            for loc in range(nloc):
                                gx_, gy_, gz_ = grad[loc, 0], grad[loc, 1], grad[loc, 2]
                                c0 = 3 * loc
                                BNL[0, c0 : c0 + 3] = F[:, 0] * gx_
                                BNL[1, c0 : c0 + 3] = F[:, 1] * gy_
                                BNL[2, c0 : c0 + 3] = F[:, 2] * gz_
                                BNL[3, c0 : c0 + 3] = F[:, 0] * gy_ + F[:, 1] * gx_
                                BNL[4, c0 : c0 + 3] = F[:, 1] * gz_ + F[:, 2] * gy_
                                BNL[5, c0 : c0 + 3] = F[:, 0] * gz_ + F[:, 2] * gx_

                            fe += wdet * (BNL.T @ S_voigt)
                            Ke += wdet * (BNL.T @ (D @ BNL))

                            # Geometric stiffness: (grad_a^T S grad_b) I_3
                            gSg = grad @ S_mat @ grad.T  # (nloc, nloc)
                            for d in range(3):
                                Ke[d::3, d::3] += wdet * gSg

                f_int[dofs] += fe
                rr, cc = np.meshgrid(dofs, dofs, indexing="ij")
                n2 = ndof_e * ndof_e
                rows[pos : pos + n2] = rr.ravel()
                cols[pos : pos + n2] = cc.ravel()
                vals[pos : pos + n2] = Ke.ravel()
                pos += n2

    K_t = sp.coo_matrix((vals, (rows, cols)), shape=(box.n_dofs, box.n_dofs)).tocsr()
    return f_int, K_t


def _solve_reference_newton(
    box: BoxSpline,
    E: float,
    nu: float,
    f_ext: np.ndarray,
    T_map: sp.csr_matrix,
    u_init: np.ndarray,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> tuple[np.ndarray, int]:
    """Reference T^T K_t T Newton-Raphson solver for u = u_init + T_map @ q."""
    q = np.zeros(T_map.shape[1])
    u = u_init.copy()
    for it in range(1, max_iter + 1):
        f_int, K_t = _assemble_tl_st_venant(box, E, nu, u)
        r_red = T_map.T @ (f_ext - f_int)
        if np.linalg.norm(r_red) < tol:
            return u, it - 1
        K_red = (T_map.T @ K_t @ T_map).tocsc()
        dq = spla.spsolve(K_red, r_red)
        q += dq
        u = u_init + T_map @ q
    raise RuntimeError(f"Reference Newton did not converge in {max_iter} iterations; |r|={np.linalg.norm(r_red)}")


def _solve_kratos_tl(
    box: BoxSpline,
    E: float,
    nu: float,
    prescribed: dict[int, float],
    point_loads: dict[int, tuple[float, float, float]],
    constraints: list[tuple[int, list[int], list[float]]],
    move_mesh: bool = True,
    criteria: str = "and",
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Build and solve Kratos TotalLagrangian model with LinearMasterSlaveConstraint.

    ``move_mesh`` / ``criteria`` ("and" | "residual") exist for the negative tests that
    back the two S21 Kratos notes (MoveMeshFlag=True, AndCriteria required).
    """
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    CL = mods["ConstitutiveLawsApplication"]

    model = KM.Model()
    mp = model.CreateModelPart("S21")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, CL.KirchhoffSaintVenant3DLaw())
    attach_cell_elements(KM, mp, vol, box, props, element_name="TotalLagrangianElement3D4N")

    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    pl_vars3 = (SM.POINT_LOAD_X, SM.POINT_LOAD_Y, SM.POINT_LOAD_Z)
    dof_of: dict[int, object] = {}
    for node in mp.Nodes:
        I = node.Id - 1
        for d, var in enumerate(vars3):
            dof = 3 * I + d
            dof_of[dof] = node.GetDof(var)
            if dof in prescribed:
                node.Fix(var)
                node.SetSolutionStepValue(var, float(prescribed[dof]))
        if I in point_loads:
            fx, fy, fz = point_loads[I]
            node.SetSolutionStepValue(pl_vars3[0], float(fx))
            node.SetSolutionStepValue(pl_vars3[1], float(fy))
            node.SetSolutionStepValue(pl_vars3[2], float(fz))
            mp.CreateNewCondition("PointLoadCondition3D1N", mp.NumberOfConditions() + 1, [node.Id], props)

    for cid, (sd, mdofs, weights) in enumerate(constraints, start=1):
        add_linear_tie(
            KM,
            mp,
            cid,
            [dof_of[sd]],
            [dof_of[m] for m in mdofs],
            np.asarray(weights, dtype=float).reshape(1, -1),
            np.zeros(1),
        )

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    if criteria == "and":
        crit = KM.AndCriteria(
            KM.ResidualCriteria(1e-12, 1e-13),
            KM.DisplacementCriteria(1e-12, 1e-13),
        )
    elif criteria == "residual":
        crit = KM.ResidualCriteria(1e-12, 1e-13)
    else:
        raise ValueError(f"unknown criteria {criteria!r}")
    strat = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 30, True, False, bool(move_mesh))
    strat.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strat.Solve()

    u = np.zeros(box.n_dofs)
    dx_pos = np.zeros(box.n_dofs)
    for node in mp.Nodes:
        I = node.Id - 1
        u[3 * I : 3 * I + 3] = [
            node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            node.GetSolutionStepValue(KM.DISPLACEMENT_Z),
        ]
        dx_pos[3 * I : 3 * I + 3] = [
            node.X - node.X0,
            node.Y - node.Y0,
            node.Z - node.Z0,
        ]
    return u, dx_pos, bool(strat.IsConverged())


def test_s21_total_lagrangian_matches_reference_no_constraints() -> None:
    """Baseline: unconstrained TotalLagrangian cantilever matches Python T^T K_t T Newton."""
    E, nu = 100.0, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(3.0, 1.0, 1.0), n_cells=(3, 1, 1), degree=1)
    X = box.control_points()
    prescribed = {}
    point_loads = {}
    f_ext = np.zeros(box.n_dofs)
    for I in range(box.n_ctrl):
        if abs(X[I, 0]) < 1e-12:
            for d in range(3):
                prescribed[3 * I + d] = 0.0
        elif abs(X[I, 0] - 3.0) < 1e-12:
            point_loads[I] = (0.0, -0.2, 0.0)
            f_ext[3 * I + 1] = -0.2

    free_dofs = [d for d in range(box.n_dofs) if d not in prescribed]
    T_map = sp.coo_matrix(
        (np.ones(len(free_dofs)), (free_dofs, np.arange(len(free_dofs)))),
        shape=(box.n_dofs, len(free_dofs)),
    ).tocsr()
    u_init = np.zeros(box.n_dofs)
    u_ref, ref_iters = _solve_reference_newton(box, E, nu, f_ext, T_map, u_init)

    u_k, dx_pos, converged = _solve_kratos_tl(box, E, nu, prescribed, point_loads, [])
    assert converged
    err = float(np.max(np.abs(u_k - u_ref)))
    pos_err = float(np.max(np.abs(dx_pos - u_k)))
    print(f"S21 baseline: ref_iters={ref_iters}, max|u|={np.max(np.abs(u_k)):.6f}, max|u_k-u_ref|={err:.3e}, pos_err={pos_err:.3e}")
    assert err < 1e-8
    assert pos_err < 1e-12


def test_s21_newton_with_linear_constraints_and_dirichlet_masters() -> None:
    """Constrained TotalLagrangian cantilever (free + fixed masters, c=0 rule) matches T^T K_t T Newton."""
    E, nu = 100.0, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(3.0, 1.0, 1.0), n_cells=(3, 2, 1), degree=1)
    X = box.control_points()

    prescribed: dict[int, float] = {}
    point_loads: dict[int, tuple[float, float, float]] = {}
    f_ext = np.zeros(box.n_dofs)
    for I in range(box.n_ctrl):
        if abs(X[I, 0]) < 1e-12:
            for d in range(3):
                prescribed[3 * I + d] = 0.0
        elif abs(X[I, 0] - 3.0) < 1e-12:
            point_loads[I] = (0.0, -0.15, 0.05)
            f_ext[3 * I + 1] = -0.15
            f_ext[3 * I + 2] = 0.05

    # Choose interior x=1, y=0.5 nodes as slaves; tie them to both fixed (x=0) and free (x=2) masters
    slave_nodes = [I for I in range(box.n_ctrl) if abs(X[I, 0] - 1.0) < 1e-12 and abs(X[I, 1] - 0.5) < 1e-12]
    assert len(slave_nodes) == 2
    slave_dofs = [3 * I + d for I in slave_nodes for d in range(3)]
    slave_set = set(slave_dofs)

    fixed_pool = sorted(prescribed.keys())
    free_pool = [d for d in range(box.n_dofs) if d not in prescribed and d not in slave_set]

    rng = np.random.default_rng(21)
    constraints: list[tuple[int, list[int], list[float]]] = []
    for sd in slave_dofs:
        comp = sd % 3
        m_fixed = [d for d in fixed_pool if d % 3 == comp][:2]
        m_free = list(rng.choice([d for d in free_pool if d % 3 == comp], size=3, replace=False))
        mdofs = m_fixed + [int(v) for v in m_free]
        w = rng.uniform(0.1, 0.5, size=len(mdofs))
        w /= np.sum(w)
        constraints.append((sd, mdofs, w.tolist()))

    # Build T_map: full u = u_init + T_map @ q (where q are independent free DOFs)
    n_free = len(free_pool)
    free_col = {d: idx for idx, d in enumerate(free_pool)}
    T_dense = np.zeros((box.n_dofs, n_free))
    u_init = np.zeros(box.n_dofs)
    for d, val in prescribed.items():
        u_init[d] = val
    for d, col in free_col.items():
        T_dense[d, col] = 1.0
    for sd, mdofs, w in constraints:
        for m, wi in zip(mdofs, w):
            if m in prescribed:
                u_init[sd] += wi * prescribed[m]
            else:
                T_dense[sd, free_col[m]] += wi
    T_map = sp.csr_matrix(T_dense)

    u_ref, ref_iters = _solve_reference_newton(box, E, nu, f_ext, T_map, u_init)
    u_k, dx_pos, converged = _solve_kratos_tl(box, E, nu, prescribed, point_loads, constraints)

    assert converged
    tie_err = max(
        abs(u_k[sd] - sum(wi * u_k[m] for m, wi in zip(mdofs, w)))
        for sd, mdofs, w in constraints
    )
    sol_err = float(np.max(np.abs(u_k - u_ref)))
    pos_err = float(np.max(np.abs(dx_pos - u_k)))
    print(
        f"S21 constrained: ref_iters={ref_iters}, max|u|={np.max(np.abs(u_k)):.6f}, "
        f"tie_err={tie_err:.3e}, pos_err={pos_err:.3e}, max|u_k-u_ref|={sol_err:.3e}"
    )
    assert tie_err < 1e-13
    assert pos_err < 1e-13
    assert sol_err < 1e-8


def _nonzero_dirichlet_case():
    """Inputs of the non-zero Dirichlet master case (no external force, motion by BCs only)."""
    E, nu = 100.0, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(2, 1, 1), degree=1)
    X = box.control_points()

    prescribed: dict[int, float] = {}
    for I in range(box.n_ctrl):
        if abs(X[I, 0]) < 1e-12:
            for d in range(3):
                prescribed[3 * I + d] = 0.0
        elif abs(X[I, 0] - 2.0) < 1e-12:
            prescribed[3 * I + 0] = 0.12
            prescribed[3 * I + 1] = 0.25
            prescribed[3 * I + 2] = -0.08

    # Tie x=1 node (I where X[I,0]==1, y==0, z==0) to a non-zero Dirichlet master at x=2 and a free node at x=1,y=1,z=0
    s_idx = next(I for I in range(box.n_ctrl) if abs(X[I, 0] - 1.0) < 1e-12 and abs(X[I, 1]) < 1e-12 and abs(X[I, 2]) < 1e-12)
    m_dir = next(I for I in range(box.n_ctrl) if abs(X[I, 0] - 2.0) < 1e-12 and abs(X[I, 1]) < 1e-12 and abs(X[I, 2]) < 1e-12)
    m_free = next(I for I in range(box.n_ctrl) if abs(X[I, 0] - 1.0) < 1e-12 and abs(X[I, 1] - 1.0) < 1e-12 and abs(X[I, 2]) < 1e-12)

    constraints = []
    slave_set = set()
    for d in range(3):
        sd = 3 * s_idx + d
        slave_set.add(sd)
        constraints.append((sd, [3 * m_dir + d, 3 * m_free + d], [0.4, 0.6]))

    free_pool = [d for d in range(box.n_dofs) if d not in prescribed and d not in slave_set]
    free_col = {d: idx for idx, d in enumerate(free_pool)}
    T_dense = np.zeros((box.n_dofs, len(free_pool)))
    u_init = np.zeros(box.n_dofs)
    for d, val in prescribed.items():
        u_init[d] = val
    for d, col in free_col.items():
        T_dense[d, col] = 1.0
    for sd, mdofs, w in constraints:
        for m, wi in zip(mdofs, w):
            if m in prescribed:
                u_init[sd] += wi * prescribed[m]
            else:
                T_dense[sd, free_col[m]] += wi
    T_map = sp.csr_matrix(T_dense)

    f_ext = np.zeros(box.n_dofs)
    u_ref, ref_iters = _solve_reference_newton(box, E, nu, f_ext, T_map, u_init)
    return box, E, nu, prescribed, constraints, u_ref, ref_iters


def test_s21_nonzero_dirichlet_master_c0_rule() -> None:
    """Non-zero prescribed master DOFs with c=0 must not double-count in nonlinear Newton."""
    box, E, nu, prescribed, constraints, u_ref, ref_iters = _nonzero_dirichlet_case()
    u_k, dx_pos, converged = _solve_kratos_tl(box, E, nu, prescribed, {}, constraints)

    assert converged
    tie_err = max(
        abs(u_k[sd] - sum(wi * u_k[m] for m, wi in zip(mdofs, w)))
        for sd, mdofs, w in constraints
    )
    sol_err = float(np.max(np.abs(u_k - u_ref)))
    pos_err = float(np.max(np.abs(dx_pos - u_k)))
    print(
        f"S21 nonzero Dirichlet master: ref_iters={ref_iters}, max|u|={np.max(np.abs(u_k)):.6f}, "
        f"tie_err={tie_err:.3e}, pos_err={pos_err:.3e}, max|u_k-u_ref|={sol_err:.3e}"
    )
    assert tie_err < 1e-13
    assert pos_err < 1e-13
    assert sol_err < 1e-8


def test_s21_negative_residual_criteria_alone_false_convergence() -> None:
    """S21 note 2: with motion by Dirichlet only, ResidualCriteria alone stops early.

    Recorded value (S21.md): max|u_k - u_ref| = 1.646e-02 after one iteration.
    The AndCriteria run of the same case reaches < 1e-8 (test above).
    """
    box, E, nu, prescribed, constraints, u_ref, _ = _nonzero_dirichlet_case()
    u_k, _dx, converged = _solve_kratos_tl(box, E, nu, prescribed, {}, constraints, criteria="residual")
    err = float(np.max(np.abs(u_k - u_ref)))
    print(f"S21 negative (ResidualCriteria only): converged={converged}, max|u_k-u_ref|={err:.3e}")
    assert converged  # Kratos reports success ...
    assert err > 1e-3  # ... on a wrong answer


def test_s21_negative_move_mesh_false() -> None:
    """S21 note 1: TotalLagrangian builds F from current node coordinates (J) and the
    initial ones (J0); with MoveMeshFlag=False the coordinates never move, so F = I and
    the Newton answer is wrong. Same case as the AndCriteria test above."""
    box, E, nu, prescribed, constraints, u_ref, _ = _nonzero_dirichlet_case()
    u_k, dx_pos, converged = _solve_kratos_tl(box, E, nu, prescribed, {}, constraints, move_mesh=False)
    err = float(np.max(np.abs(u_k - u_ref)))
    moved = float(np.max(np.abs(dx_pos)))
    print(f"S21 negative (MoveMeshFlag=False): converged={converged}, max|u_k-u_ref|={err:.3e}, max|X-X0|={moved:.3e}")
    assert moved == 0.0
    assert (not converged) or err > 1e-3
