"""S24: Explicit dynamics with Kratos QP elements, HRZ lumped mass, and T-extension.

Investigates:
1. Kratos MechanicalExplicitStrategy + ExplicitCentralDifferencesScheme behavior
   on QuadraturePointGeometry SM elements (NODAL_MASS, FORCE_RESIDUAL, MPC).
2. Python explicit central-difference loop (Kratos-extracted K + HRZ lumped mass
   + T^T f_int reduction) on 1D axial step-loaded wave propagation vs analytic.
3. Cut-cell (eta = 1e-4, 1e-6) critical time step dt_crit collapse without
   extension vs restoration to O(h/c) with polynomial extension T.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import eigh

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    solve_linear,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference import assemble_elastic as ref
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _explicit_model(KM, SM, name: str = "S24"):
    """Create a ModelPart configured for Kratos MechanicalExplicitStrategy."""
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.VELOCITY)
    mp.AddNodalSolutionStepVariable(KM.ACCELERATION)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.MIDDLE_VELOCITY)
    mp.AddNodalSolutionStepVariable(KM.NODAL_MASS)
    mp.AddNodalSolutionStepVariable(KM.FORCE_RESIDUAL)
    mp.AddNodalSolutionStepVariable(KM.RESIDUAL_VECTOR)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.COMPUTE_LUMPED_MASS_MATRIX, True)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_ALPHA, 0.0)
    mp.ProcessInfo.SetValue(SM.RAYLEIGH_BETA, 0.0)
    return model, mp


def _static_model(KM, name: str = "S24K"):
    """Create a ModelPart for Kratos K-matrix extraction."""
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def _props(KM, SM, mp, E: float = 100.0, nu: float = 0.0, rho: float = 1.0):
    """Create linear elastic 3D material properties."""
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.DENSITY, rho)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    return props


def _lagrange_weights(x_eval: float, x_nodes: np.ndarray) -> np.ndarray:
    """Compute 1D Lagrange interpolation/extrapolation weights at x_eval."""
    n = len(x_nodes)
    w = np.ones(n, dtype=float)
    for j in range(n):
        for k in range(n):
            if j != k:
                w[j] *= (x_eval - x_nodes[k]) / (x_nodes[j] - x_nodes[k])
    return w


def _build_x_extension_T(box: BoxSpline) -> tuple[np.ndarray, int]:
    """Build polynomial extension T eliminating the last x-layer control point.

    Maps u_full = T @ u_master where the last control point in x (supported only
    on the trimmed last cell) is extrapolated from the p+1 nearest interior
    Greville abscissae.
    """
    p = box.degree
    nu_cp, nv_cp, nw_cp = box.n_cp
    kv = box.knots[0]
    gx = np.array([np.mean(kv[i + 1 : i + 1 + p]) for i in range(nu_cp)])
    n_good_x = nu_cp - 1
    Tx = np.zeros((nu_cp, n_good_x), dtype=float)
    Tx[:n_good_x, :n_good_x] = np.eye(n_good_x)
    src = np.arange(n_good_x - (p + 1), n_good_x)
    Tx[nu_cp - 1, src] = _lagrange_weights(gx[nu_cp - 1], gx[src])
    T_scalar = np.kron(np.eye(nw_cp * nv_cp), Tx)
    T = np.kron(T_scalar, np.eye(3))
    return T, n_good_x


def _trimmed_cell_gauss_points(box: BoxSpline, i: int, j: int, k: int, eta_x: float) -> list[list[float]]:
    """Build [[x, y, z, w], ...] for a cell trimmed to [0, eta_x] in x."""
    pts = cell_gauss_points(box, i, j, k)
    if eta_x == 1.0:
        return pts
    u0 = float(box.breakpoints(0)[i])
    return [[u0 + eta_x * (x - u0), y, z, eta_x * w] for x, y, z, w in pts]


def _assemble_trimmed_consistent_mass(box: BoxSpline, rho: float, eta_last_x: float = 1.0) -> np.ndarray:
    """Assemble consistent mass matrix M_c on a box with optional x-trimming on the last cell."""
    p = box.degree
    ng = p + 1
    xg, wg = np.polynomial.legendre.leggauss(ng)
    tabs = [box.axis_tables(d, xg) for d in range(3)]
    (fu, Nu, _, _, ju), (fv, Nv, _, _, jv), (fw, Nw, _, _, jw) = tabs
    xi_cut = -1.0 + eta_last_x * (xg + 1.0)
    Nu_cut = box.axis_tables(0, xi_cut)[1]
    nx, ny, nz = box.n_cells
    nu_cp, nv_cp, _ = box.n_cp
    M = np.zeros((box.n_dofs, box.n_dofs), dtype=float)
    nloc = (p + 1) ** 3

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                is_cut = (i == nx - 1) and (eta_last_x < 1.0)
                Nu_i = Nu_cut[i] if is_cut else Nu[i]
                ju_i = ju[i] * eta_last_x if is_cut else ju[i]
                dofs = []
                for c_ in range(p + 1):
                    for b in range(p + 1):
                        for a in range(p + 1):
                            I = (fu[i] + a) + nu_cp * ((fv[j] + b) + nv_cp * (fw[k] + c_))
                            dofs.extend([3 * I, 3 * I + 1, 3 * I + 2])
                Me = np.zeros((3 * nloc, 3 * nloc), dtype=float)
                for gz in range(ng):
                    for gy in range(ng):
                        for gx in range(ng):
                            wdet = wg[gx] * wg[gy] * wg[gz] * ju_i * jv[j] * jw[k]
                            Nvec = np.zeros(nloc, dtype=float)
                            for c_ in range(p + 1):
                                for b in range(p + 1):
                                    for a in range(p + 1):
                                        loc = a + (p + 1) * (b + (p + 1) * c_)
                                        Nvec[loc] = Nu_i[gx, a] * Nv[j, gy, b] * Nw[k, gz, c_]
                            Nmat = np.zeros((3, 3 * nloc), dtype=float)
                            for d in range(3):
                                Nmat[d, d::3] = Nvec
                            Me += wdet * rho * (Nmat.T @ Nmat)
                M[np.ix_(dofs, dofs)] += Me
    return M


def _hrz_lump_diagonal(M_consistent: np.ndarray, total_mass: float) -> np.ndarray:
    """Apply HRZ (Hinton-Rock-Zienkiewicz) diagonal scaling per spatial component."""
    diag = np.diag(M_consistent).copy()
    for d in range(3):
        s = float(np.sum(diag[d::3]))
        assert s > 0.0, f"HRZ diagonal sum must be positive, got {s}"
        diag[d::3] *= total_mass / s
    return diag


@pytest.mark.parametrize("degree", [1, 2])
def test_kratos_explicit_qp_nodal_mass_and_force_residual(degree: int):
    """Verify how Kratos MechanicalExplicitStrategy handles NODAL_MASS and FORCE_RESIDUAL on QP SM elements."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    model, mp = _explicit_model(KM, SM, f"S24_QP_p{degree}")
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(2, 1, 1), degree=degree)
    vol, _ = make_nurbs_volume(KM, mp, box)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)

    E, nu, rho = 100.0, 0.0, 2.0
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp, E=E, nu=nu, rho=rho))

    # Prescribe x=2 face displacement u_x = 0.1 and fix it so initial FORCE_RESIDUAL is non-zero
    for node in mp.Nodes:
        if abs(node.X0 - 2.0) < 1e-12:
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.1)
            node.Fix(KM.DISPLACEMENT_X)

    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1e-8)
    strat.Solve()

    # 1. NODAL_MASS is stored as NonHistorical (GetValue), NOT Historical (GetSolutionStepValue)
    hist_masses = [node.GetSolutionStepValue(KM.NODAL_MASS) for node in mp.Nodes]
    nonhist_masses = [node.GetValue(KM.NODAL_MASS) for node in mp.Nodes]
    assert all(m == 0.0 for m in hist_masses)
    expected_total_mass = rho * 2.0 * 1.0 * 1.0  # rho * V = 4.0
    assert abs(sum(nonhist_masses) - expected_total_mass) < 1e-12
    assert min(nonhist_masses) > 0.0

    # 2. FORCE_RESIDUAL matches -K @ u_init (at tiny dt=1e-8 before significant motion)
    f_res = np.zeros(box.n_dofs)
    for node in mp.Nodes:
        I = node.Id - 1
        fr = node.GetSolutionStepValue(KM.FORCE_RESIDUAL)
        f_res[3 * I : 3 * I + 3] = [fr[0], fr[1], fr[2]]
    u0 = np.zeros(box.n_dofs)
    for node in mp.Nodes:
        I = node.Id - 1
        u0[3 * I] = 0.1 if abs(node.X0 - 2.0) < 1e-12 else 0.0
    K_ref = ref.assemble_stiffness(box, E, nu).toarray()
    np.testing.assert_allclose(f_res, -(K_ref @ u0), atol=1e-4)
    print(
        f"S24 Q1 (p={degree}): sum(NODAL_MASS)={sum(nonhist_masses):.6f}, "
        f"min_mass={min(nonhist_masses):.6f}, max|f_res+Ku|={np.max(np.abs(f_res + K_ref @ u0)):.3e}"
    )


def test_kratos_explicit_mpc_mass_and_residual_behavior():
    """Verify Kratos MechanicalExplicitStrategy MPC behavior: T^T f_res vs non-variational T^T m_slave."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    # Part A: positive tie u_slave = 2.0 * u_master
    model, mp = _explicit_model(KM, SM, "S24_MPC_Pos")
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0), n_cells=(2, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp, E=100.0, nu=0.0, rho=2.0))

    nodes = {n.Id: n for n in mp.Nodes}
    master_ids = []
    slave_ids = []
    for j in range(2):
        for k in range(2):
            m_id = 1 + 1 + 3 * (j + 2 * k)
            s_id = 1 + 2 + 3 * (j + 2 * k)
            master_ids.append(m_id)
            slave_ids.append(s_id)
            add_linear_tie(
                KM,
                mp,
                s_id,
                [nodes[s_id].GetDof(KM.DISPLACEMENT_X)],
                [nodes[m_id].GetDof(KM.DISPLACEMENT_X)],
                np.array([[2.0]]),
            )
            nodes[m_id].SetSolutionStepValue(KM.DISPLACEMENT_X, 0.1)
            nodes[s_id].SetSolutionStepValue(KM.DISPLACEMENT_X, 0.2)

    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1e-8)
    strat.Solve()

    # Check 1: ApplyConstraints restores u_slave = 2.0 * u_master
    for m_id, s_id in zip(master_ids, slave_ids):
        um = nodes[m_id].GetSolutionStepValue(KM.DISPLACEMENT_X)
        us = nodes[s_id].GetSolutionStepValue(KM.DISPLACEMENT_X)
        assert abs(us - 2.0 * um) < 1e-12

    # Check 2: FORCE_RESIDUAL on master includes T^T f_slave = 0 + 2 * (-2.5) = -5.0
    f_res_m = nodes[master_ids[0]].GetSolutionStepValue(KM.FORCE_RESIDUAL)[0]
    f_res_s = nodes[slave_ids[0]].GetSolutionStepValue(KM.FORCE_RESIDUAL)[0]
    assert abs(f_res_s - (-2.5)) < 1e-4
    assert abs(f_res_m - (-5.0)) < 1e-4

    # Check 3: NODAL_MASS on master is m_m + T * m_s = 0.5 + 2.0 * 0.25 = 1.0
    # instead of variational T^T M T = 0.5 + 2.0^2 * 0.25 = 1.5, and slave mass is NOT zeroed!
    m_master = nodes[master_ids[0]].GetValue(KM.NODAL_MASS)
    m_slave = nodes[slave_ids[0]].GetValue(KM.NODAL_MASS)
    total_mass_kratos = sum(n.GetValue(KM.NODAL_MASS) for n in mp.Nodes)
    assert abs(m_master - 1.0) < 1e-12
    assert abs(m_slave - 0.25) < 1e-12
    assert abs(total_mass_kratos - 6.0) < 1e-12  # != true mass 4.0

    # Part B: Negative extension weight (T = [-1, 2] as in p=1 B-spline extension)
    model2, mp2 = _explicit_model(KM, SM, "S24_MPC_Neg")
    vol2, _ = make_nurbs_volume(KM, mp2, box)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp2)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp2)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp2)
    attach_cell_elements(KM, mp2, vol2, box, _props(KM, SM, mp2, E=100.0, nu=0.0, rho=2.0))
    nodes2 = {n.Id: n for n in mp2.Nodes}
    for j in range(2):
        for k in range(2):
            m0_id = 1 + 0 + 3 * (j + 2 * k)
            m1_id = 1 + 1 + 3 * (j + 2 * k)
            s_id = 1 + 2 + 3 * (j + 2 * k)
            add_linear_tie(
                KM,
                mp2,
                s_id,
                [nodes2[s_id].GetDof(KM.DISPLACEMENT_X)],
                [nodes2[m0_id].GetDof(KM.DISPLACEMENT_X), nodes2[m1_id].GetDof(KM.DISPLACEMENT_X)],
                np.array([[-1.0, 2.0]]),
            )
    scheme2 = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat2 = SM.MechanicalExplicitStrategy(mp2, scheme2, False, False, False)
    strat2.SetEchoLevel(0)
    strat2.Initialize()
    mp2.CloneTimeStep(1e-8)
    strat2.Solve()
    # Node i=0 has base mass 0.25, receives (-1.0) * 0.25 = 0.0!
    m_i0_neg = nodes2[1].GetValue(KM.NODAL_MASS)
    assert abs(m_i0_neg - 0.0) < 1e-12, f"Expected 0.0 mass due to T=-1, got {m_i0_neg}"
    print(
        f"S24 Q1 MPC: m_master={m_master:.3f} (variational=1.500), "
        f"total_mass={total_mass_kratos:.3f} (true=4.0), m_i0_with_neg_T={m_i0_neg:.3e}"
    )


def test_python_explicit_1d_wave_propagation_hrz_and_mpc():
    """Verify Python explicit central-difference loop with Kratos K, HRZ mass, and T^T f_int."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    L, E, nu, rho = 10.0, 100.0, 0.0, 1.0
    c = float(np.sqrt(E / rho))  # 10.0
    nx = 40
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L, 1.0, 1.0), n_cells=(nx, 1, 1), degree=2)

    # Extract K from Kratos QP elements
    _, mp = _static_model(KM, "S24_Wave")
    vol, _ = make_nurbs_volume(KM, mp, box)
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp, E=E, nu=nu, rho=rho))
    strat = solve_linear(KM, mp)
    A = system_matrix_csr(strat)
    perm = equation_permutation(KM, mp, box)
    inv_perm = np.argsort(perm)
    K_full = A[np.ix_(inv_perm, inv_perm)].toarray()

    # Run 1: Untrimmed bar with HRZ lumped mass
    M_c = _assemble_trimmed_consistent_mass(box, rho=rho, eta_last_x=1.0)
    m_diag = _hrz_lump_diagonal(M_c, total_mass=rho * L * 1.0 * 1.0)

    sigma0 = 1.0
    f_ext = ref.assemble_face_traction(box, axis=0, side=1, traction=[sigma0, 0.0, 0.0])

    nu_cp = box.n_cp[0]
    fixed = np.zeros(box.n_dofs, dtype=bool)
    for I in range(box.n_ctrl):
        if I % nu_cp == 0:
            fixed[3 * I] = True
        fixed[3 * I + 1] = True
        fixed[3 * I + 2] = True
    free = ~fixed

    dt = 0.005  # h = 0.25, h/c = 0.025 => CFL = 0.2
    t_end = 2.0
    n_steps = int(round(t_end / dt))

    def _run_cd(K_mat, M_diag_vec, f_ext_vec, free_mask, tip_dofs, mid_dofs):
        u = np.zeros(len(f_ext_vec), dtype=float)
        m_inv = np.zeros_like(M_diag_vec)
        m_inv[free_mask] = 1.0 / M_diag_vec[free_mask]
        a = m_inv * (f_ext_vec - K_mat @ u)
        v_half = -0.5 * dt * a
        t_arr = np.empty(n_steps, dtype=float)
        u_tip_hist = np.empty(n_steps, dtype=float)
        v_mid_hist = np.empty(n_steps, dtype=float)
        for s in range(n_steps):
            v_half[free_mask] += dt * a[free_mask]
            u[free_mask] += dt * v_half[free_mask]
            t_arr[s] = (s + 1) * dt
            u_tip_hist[s] = float(np.mean(u[tip_dofs]))
            v_mid_hist[s] = float(np.mean(v_half[mid_dofs]))
            a = m_inv * (f_ext_vec - K_mat @ u)
        return t_arr, u_tip_hist, v_mid_hist

    tip_ids = [3 * I for I in range(box.n_ctrl) if I % nu_cp == nu_cp - 1]
    mid_ids = [3 * I for I in range(box.n_ctrl) if I % nu_cp == nu_cp // 2]
    t_hist, u_tip, v_mid = _run_cd(K_full, m_diag, f_ext, free, tip_ids, mid_ids)

    u_exact_tip = (sigma0 * c / E) * t_hist
    rel_l2 = float(np.linalg.norm(u_tip - u_exact_tip) / np.linalg.norm(u_exact_tip))
    assert rel_l2 < 0.02, f"Tip L2 relative error too large: {rel_l2}"

    v0 = sigma0 * c / E  # 0.1
    idx_arr = int(np.argmax(v_mid > 0.5 * v0))
    t_arr_meas = float(t_hist[idx_arr])
    assert abs(t_arr_meas - 0.5) < 0.03, f"Midpoint wavefront arrival {t_arr_meas} vs 0.5"

    # Run 2: Trimmed last cell (eta = 1e-4) with T-extension and T^T f_int reduction, SAME dt = 0.005
    eta = 1e-4
    box_t = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L, 1.0, 1.0), n_cells=(nx, 1, 1), degree=2)
    _, mp_t = _static_model(KM, "S24_WaveTrim")
    vol_t, _ = make_nurbs_volume(KM, mp_t, box_t)
    props_t = _props(KM, SM, mp_t, E=E, nu=nu, rho=rho)
    eid = 1
    for k in range(1):
        for j in range(1):
            for i in range(nx):
                qps = KM.GeometriesVector()
                pts = _trimmed_cell_gauss_points(box_t, i, j, k, eta if i == nx - 1 else 1.0)
                vol_t.CreateQuadraturePointGeometries(qps, 2, pts)
                mp_t.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props_t)
                eid += 1
    strat_t = solve_linear(KM, mp_t)
    A_t = system_matrix_csr(strat_t)
    perm_t = equation_permutation(KM, mp_t, box_t)
    inv_perm_t = np.argsort(perm_t)
    K_trim = A_t[np.ix_(inv_perm_t, inv_perm_t)].toarray()

    M_trim = _assemble_trimmed_consistent_mass(box_t, rho=rho, eta_last_x=eta)
    T, n_good_x = _build_x_extension_T(box_t)
    K_red = T.T @ K_trim @ T
    M_red_c = T.T @ M_trim @ T
    h_cell = L / nx
    L_trim = (nx - 1 + eta) * h_cell
    m_red_diag = _hrz_lump_diagonal(M_red_c, total_mass=rho * L_trim * 1.0 * 1.0)

    # Consistent traction & displacement evaluation at the physical trimmed tip x = L_trim
    # Map L_trim into last cell local coord xi_tip in [-1, 1]
    xi_tip = np.array([-1.0 + 2.0 * eta])
    fu_t, Nu_t, _, _, _ = box_t.axis_tables(0, xi_tip)
    # 1D basis weights at x = L_trim
    w_x_tip = np.zeros(box_t.n_cp[0], dtype=float)
    w_x_tip[fu_t[nx - 1] : fu_t[nx - 1] + box_t.degree + 1] = Nu_t[nx - 1, 0, :]
    # Distribute sigma0 * Area across y,z control points using unit-integral partition of unity
    # (for constant traction, f_ext_yz sums to 1.0 over yz)
    f_unit_yz = ref.assemble_face_traction(box_t, axis=0, side=1, traction=[1.0, 0.0, 0.0])
    # Extract yz weights from side=1 face (where only last x-layer is nonzero)
    yz_weights = f_unit_yz[3 * (box_t.n_cp[0] - 1) :: 3 * box_t.n_cp[0]]
    f_ext_trim = np.zeros(box_t.n_dofs, dtype=float)
    for k_cp in range(box_t.n_cp[2]):
        for j_cp in range(box_t.n_cp[1]):
            yz_idx = j_cp + box_t.n_cp[1] * k_cp
            for i_cp in range(box_t.n_cp[0]):
                I = i_cp + box_t.n_cp[0] * yz_idx
                f_ext_trim[3 * I] = sigma0 * w_x_tip[i_cp] * yz_weights[yz_idx]
    f_ext_red = T.T @ f_ext_trim

    # Physical tip displacement evaluator: u_x(L_trim) = sum_i w_x_tip[i] * mean_yz((T @ u_m)[i])
    tip_eval_vec = np.zeros(box_t.n_dofs, dtype=float)
    n_yz = box_t.n_cp[1] * box_t.n_cp[2]
    for k_cp in range(box_t.n_cp[2]):
        for j_cp in range(box_t.n_cp[1]):
            yz_idx = j_cp + box_t.n_cp[1] * k_cp
            for i_cp in range(box_t.n_cp[0]):
                I = i_cp + box_t.n_cp[0] * yz_idx
                tip_eval_vec[3 * I] = w_x_tip[i_cp] / n_yz
    tip_eval_red = T.T @ tip_eval_vec

    n_master_ctrl = T.shape[1] // 3
    free_red = np.ones(T.shape[1], dtype=bool)
    for I in range(n_master_ctrl):
        if I % n_good_x == 0:
            free_red[3 * I] = False
        free_red[3 * I + 1] = False
        free_red[3 * I + 2] = False

    # Run CD up to 2 * L_trim / c (before reflection returns to trimmed tip)
    n_steps_trim = int(round((2.0 * L_trim / c) / dt))
    u_m = np.zeros(T.shape[1], dtype=float)
    m_inv_r = np.zeros(T.shape[1], dtype=float)
    m_inv_r[free_red] = 1.0 / m_red_diag[free_red]
    a_m = m_inv_r * (f_ext_red - K_red @ u_m)
    v_half_m = -0.5 * dt * a_m
    t_hist2 = np.empty(n_steps_trim, dtype=float)
    u_tip2 = np.empty(n_steps_trim, dtype=float)
    for s in range(n_steps_trim):
        v_half_m[free_red] += dt * a_m[free_red]
        u_m[free_red] += dt * v_half_m[free_red]
        t_hist2[s] = (s + 1) * dt
        u_tip2[s] = float(tip_eval_red @ u_m)
        a_m = m_inv_r * (f_ext_red - K_red @ u_m)

    u_exact_tip2 = (sigma0 * c / E) * t_hist2
    rel_l2_trim = float(np.linalg.norm(u_tip2 - u_exact_tip2) / np.linalg.norm(u_exact_tip2))
    assert rel_l2_trim < 0.02, f"Trimmed+T tip L2 error too large: {rel_l2_trim}"
    print(
        f"S24 Q2: untrimmed rel_L2={rel_l2:.4e}, t_arr={t_arr_meas:.4f}s (exact 0.5000s); "
        f"trimmed(eta=1e-4)+T rel_L2={rel_l2_trim:.4e} at same dt={dt}"
    )


@pytest.mark.parametrize("degree", [1, 2])
def test_cut_cell_dt_crit_collapse_and_extension_restoration(degree: int):
    """Compare cut-cell w_max and dt_crit with vs without extension T using Kratos-extracted K."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    E, nu, rho = 100.0, 0.0, 1.0
    c = float(np.sqrt(E / rho))  # 10.0
    nx = 4
    h = 1.0
    dt_ref_hc = h / c  # 0.1

    results = {}
    for eta in (1.0, 1e-2, 1e-4, 1e-6):
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(float(nx), 1.0, 1.0), n_cells=(nx, 1, 1), degree=degree)
        _, mp = _static_model(KM, f"S24_Cut_p{degree}_{eta:.0e}")
        vol, _ = make_nurbs_volume(KM, mp, box)
        props = _props(KM, SM, mp, E=E, nu=nu, rho=rho)
        eid = 1
        for k in range(1):
            for j in range(1):
                for i in range(nx):
                    qps = KM.GeometriesVector()
                    pts = _trimmed_cell_gauss_points(box, i, j, k, eta if i == nx - 1 else 1.0)
                    vol.CreateQuadraturePointGeometries(qps, 2, pts)
                    mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                    eid += 1
        strat = solve_linear(KM, mp)
        A = system_matrix_csr(strat)
        perm = equation_permutation(KM, mp, box)
        inv_perm = np.argsort(perm)
        K_kratos = A[np.ix_(inv_perm, inv_perm)].toarray()

        M_c = _assemble_trimmed_consistent_mass(box, rho=rho, eta_last_x=eta)
        total_vol = (nx - 1) + eta

        # 1. Without extension (T = I)
        m_hrz_no = _hrz_lump_diagonal(M_c, total_mass=rho * total_vol)
        nu_cp = box.n_cp[0]
        free_no = [3 * I + d for I in range(box.n_ctrl) if I % nu_cp != 0 for d in range(3)]
        vals_no = eigh(
            K_kratos[np.ix_(free_no, free_no)],
            np.diag(m_hrz_no[free_no]),
            eigvals_only=True,
        )
        wmax_no = float(np.sqrt(np.max(vals_no)))
        dt_no = 2.0 / wmax_no

        # 2. With polynomial extension T
        T, n_good_x = _build_x_extension_T(box)
        K_ext = T.T @ K_kratos @ T
        M_ext_c = T.T @ M_c @ T
        m_hrz_ext = _hrz_lump_diagonal(M_ext_c, total_mass=rho * total_vol)
        n_master_ctrl = T.shape[1] // 3
        free_ext = [3 * I + d for I in range(n_master_ctrl) if I % n_good_x != 0 for d in range(3)]
        vals_ext = eigh(
            K_ext[np.ix_(free_ext, free_ext)],
            np.diag(m_hrz_ext[free_ext]),
            eigvals_only=True,
        )
        wmax_ext = float(np.sqrt(np.max(vals_ext)))
        dt_ext = 2.0 / wmax_ext
        results[eta] = (wmax_no, dt_no, wmax_ext, dt_ext)
        print(
            f"S24 Q3 p={degree} eta={eta:.0e}: "
            f"no-ext dt_crit={dt_no:.4e} (wmax={wmax_no:.4e}) | "
            f"with-ext dt_crit={dt_ext:.4e} (wmax={wmax_ext:.4e})"
        )

    # Assertions: severe collapse without T, full restoration to O(h/c) with T
    _, dt_no_1e4, _, dt_ext_1e4 = results[1e-4]
    _, dt_no_1e6, _, dt_ext_1e6 = results[1e-6]
    assert dt_no_1e4 < 2e-4 * dt_ref_hc, f"Expected O(eta) collapse at 1e-4, got {dt_no_1e4}"
    assert dt_ext_1e4 > 0.5 * dt_ref_hc, f"Expected O(h/c) restoration at 1e-4, got {dt_ext_1e4}"
    assert dt_ext_1e6 > 0.5 * dt_ref_hc, f"Expected O(h/c) restoration at 1e-6, got {dt_ext_1e6}"
    assert dt_ext_1e4 / dt_no_1e4 > 5000.0
