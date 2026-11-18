# -*- coding: utf-8 -*-
"""V08 Cases R1 & R2: deformable-deformable surface Gauss SSN contact.

Tests:
1. Deformable-deformable flat p=2 blocks patch test (1D two-spring + penalty series,
   cut-face reaction equilibrium) & unilateral separation.
2. Deformable-deformable Coulomb friction patch test (mu = 0.3 -> |F_T|/F_N = 0.300000
   and cut-face horizontal reaction equilibrium).
3. V08 Case R1 (sym="half", trimmed hollow hemisphere indenter + deformable block,
   h0=0.5 and h1=0.25, t=0->0.1): n_fallback=0, equilibrium, comparison vs Hertz
   (481,579.3 lbf), Sneddon (471,927.7 lbf), and R0 rigid indenter (S31).
4. V08 Case R2 (sym="half", mu=0.3, press t=0->0.1 then slide t=0.1->0.40):
   F_T/F_N plateau, cut-face equilibrium, and relative_rms over 0.15 <= t <= 0.40.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import (
    DeformableBoxFaceMaster,
    sample_surface_triangles_gauss,
    solve_deformable_ssn_contact,
)
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import equation_permutation
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

pytestmark = [pytest.mark.kratos]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"


def _assemble_kratos_partitioned_stiffness(
    KM,
    SM,
    box: BoxSpline,
    E: float,
    nu: float,
    *,
    zero_fn=None,
    presc_fn=None,
    trim_soup: tuple[np.ndarray, np.ndarray] | None = None,
    q: int = 4,
    ghost_alpha: float = 1e-8,
):
    """Assemble unconstrained Kratos stiffness for a (trimmed) BoxSpline & partition DOFs.

    Uses ``builder.BuildLHS`` without fixing Kratos DOFs so ``K_ff``, ``K_fp``,
    ``K_pf``, and ``K_pp`` retain the exact coupling needed for prescribed
    Dirichlet loads ``f_ext = -(K_fp @ u_presc)`` and cut-face reaction recovery
    ``R_presc = K_pf @ u_free + K_pp @ u_presc + f_c_s_full[presc_mask]``.
    Inactive control points (``inactive_ctrl_mask``) are placed in ``zero_mask``.
    """
    import KratosMultiphysics.scipy_conversion_tools as sct

    if trim_soup is not None:
        nodes_soup, tris_soup = trim_soup
        qp_table, out_cells, _issues, _rec, diag = build_trimmed_table(
            nodes_soup, tris_soup, box, q=int(q), cache_dir="off"
        )
        inact_ctrl = inactive_ctrl_mask(box, out_cells)
    else:
        qp_table = {}
        out_cells = []
        diag = {"n_cut": 0, "n_fallback": 0, "n_relaxed": 0, "n_qp_cut": 0}
        inact_ctrl = np.zeros(box.n_ctrl, dtype=bool)

    out_set = {tuple(c) for c in out_cells}

    model = KM.Model()
    mp = model.CreateModelPart("SSNPart")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    vol, _ = make_nurbs_volume(KM, mp, box)
    vol_props = mp.CreateNewProperties(1)
    vol_props.SetValue(KM.YOUNG_MODULUS, float(E))
    vol_props.SetValue(KM.POISSON_RATIO, float(nu))
    vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    ghost_props = None
    if out_set and ghost_alpha > 0.0:
        ghost_props = mp.CreateNewProperties(2)
        ghost_props.SetValue(KM.YOUNG_MODULUS, float(E) * float(ghost_alpha))
        ghost_props.SetValue(KM.POISSON_RATIO, float(nu))
        ghost_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    nx, ny, nz = (int(v) for v in box.n_cells)
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (i, j, k) in out_set:
                    if ghost_alpha <= 0.0:
                        continue
                    pts = cell_gauss_points(box, i, j, k)
                    props = ghost_props
                elif (i, j, k) in qp_table:
                    pts = np.asarray(qp_table[(i, j, k)], dtype=float).tolist()
                    props = vol_props
                else:
                    pts = cell_gauss_points(box, i, j, k)
                    props = vol_props
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                assert len(qps) == 1
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                eid += 1

    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    strat = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.InitializeSolutionStep()
    A = strat.GetSystemMatrix()
    builder.BuildLHS(scheme, mp, A)
    Kkr = sct.to_csr(A)
    perm = equation_permutation(KM, mp, box)
    inv_perm = np.argsort(perm)
    K_full = Kkr[inv_perm, :][:, inv_perm].tocsr()

    Xcp = box.control_points()
    zero_mask = np.zeros(box.n_dofs, dtype=bool)
    presc_mask = np.zeros(box.n_dofs, dtype=bool)

    for I in np.flatnonzero(inact_ctrl):
        zero_mask[3 * I : 3 * I + 3] = True
    if zero_fn is not None:
        for I in range(box.n_ctrl):
            for d in zero_fn(Xcp[I]):
                zero_mask[3 * I + d] = True
    if presc_fn is not None:
        for I in range(box.n_ctrl):
            if inact_ctrl[I]:
                continue
            for d in presc_fn(Xcp[I]):
                if not zero_mask[3 * I + d]:
                    presc_mask[3 * I + d] = True

    free_mask = ~(zero_mask | presc_mask)
    return {
        "K_ff": K_full[free_mask, :][:, free_mask].tocsr(),
        "K_fp": K_full[free_mask, :][:, presc_mask].tocsr(),
        "K_pf": K_full[presc_mask, :][:, free_mask].tocsr(),
        "K_pp": K_full[presc_mask, :][:, presc_mask].tocsr(),
        "free_mask": free_mask,
        "presc_mask": presc_mask,
        "zero_mask": zero_mask,
        "inact_ctrl": inact_ctrl,
        "diag": diag,
    }


def _box_face_triangles(box: BoxSpline, axis: int = 1, side: int = 0):
    """Build triangulated mesh ``(nodes, tris)`` on a BoxSpline boundary face."""
    others = [d for d in range(3) if d != axis]
    b0 = box.breakpoints(others[0])
    b1 = box.breakpoints(others[1])
    x_fixed = float(box.lo[axis] if side == 0 else box.hi[axis])
    n0, n1 = len(b0), len(b1)
    nodes = np.empty((n0 * n1, 3), dtype=np.float64)
    for j, v1 in enumerate(b1):
        for i, v0 in enumerate(b0):
            idx = i + n0 * j
            nodes[idx, axis] = x_fixed
            nodes[idx, others[0]] = float(v0)
            nodes[idx, others[1]] = float(v1)
    tris = []
    for j in range(n1 - 1):
        for i in range(n0 - 1):
            a = i + n0 * j
            b = (i + 1) + n0 * j
            c = (i + 1) + n0 * (j + 1)
            d = i + n0 * (j + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    return nodes, np.asarray(tris, dtype=np.int64)


def _prescribed_vector(presc_mask: np.ndarray, ux: float, uy: float, uz: float = 0.0):
    """Build prescribed-DOF vector matching ``presc_mask`` component indices."""
    idx = np.flatnonzero(presc_mask)
    u_p = np.zeros(idx.size, dtype=np.float64)
    comp = idx % 3
    u_p[comp == 0] = float(ux)
    u_p[comp == 1] = float(uy)
    u_p[comp == 2] = float(uz)
    return u_p


def _cut_face_reaction(sys_part: dict, u_free: np.ndarray, u_presc: np.ndarray, f_c_full: np.ndarray) -> np.ndarray:
    """Recover total ``(Fx, Fy, Fz)`` reaction on the prescribed cut face."""
    pm = sys_part["presc_mask"]
    R_p = sys_part["K_pf"] @ u_free + sys_part["K_pp"] @ u_presc + f_c_full[pm]
    idx = np.flatnonzero(pm) % 3
    return np.array(
        [float(R_p[idx == 0].sum()), float(R_p[idx == 1].sum()), float(R_p[idx == 2].sum())],
        dtype=np.float64,
    )


def test_deformable_flat_blocks_patch_and_separation():
    """Test 1: two p=2 cubes in series-spring contact & unilateral separation."""
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    E1, E2, nu = 1000.0, 2000.0, 0.0
    delta, eps_n = 0.05, 1e8

    box_m = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    box_s = BoxSpline(lo=(0.0, 1.0, 0.0), hi=(1.0, 2.0, 1.0), n_cells=(2, 2, 2), degree=2)

    sys_m = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_m, E1, nu,
        zero_fn=lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
    )
    sys_s = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_s, E2, nu,
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - 2.0) < 1e-12 else (),
    )

    tri_nodes, tris = _box_face_triangles(box_s, axis=1, side=0)
    X_s, w_q, B_s = sample_surface_triangles_gauss(tri_nodes, tris, box_s, order=2)
    assert float(w_q.sum()) == pytest.approx(1.0, rel=1e-14)

    master_face = DeformableBoxFaceMaster(box_m, axis=1, side=1)

    # 1a. Compression delta = 0.05
    u_presc = _prescribed_vector(sys_s["presc_mask"], 0.0, -delta, 0.0)
    u_fixed_s = np.zeros(box_s.n_dofs, dtype=np.float64)
    u_fixed_s[sys_s["presc_mask"]] = u_presc
    f_ext_s = -(sys_s["K_fp"] @ u_presc)

    res = solve_deformable_ssn_contact(
        sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, eps_n,
        free_m=sys_m["free_mask"],
        free_s=sys_s["free_mask"], u_fixed_s=u_fixed_s, f_ext_s=f_ext_s,
        use_pardiso=False,
    )

    k1, k2, k_pen = E1 * 1.0 / 1.0, E2 * 1.0 / 1.0, eps_n * 1.0
    F_series = delta / (1.0 / k1 + 1.0 / k2 + 1.0 / k_pen)
    F_N = float(-res["force_total"][1])
    R_cut = _cut_face_reaction(sys_s, res["u_free_s"], u_presc, res["f_c_s_full"])

    assert res["n_active"] == len(w_q)
    assert F_N == pytest.approx(F_series, rel=1e-9)
    assert R_cut[1] == pytest.approx(float(res["force_total"][1]), rel=1e-10)
    assert abs(R_cut[0]) < 1e-8 and abs(R_cut[2]) < 1e-8

    # Check exact interface displacements on master top & slave bottom
    u_m_top_y = (res["B_m"] @ res["u_ctrl_m"])[:, 1]
    u_s_bot_y = (B_s @ res["u_ctrl_s"])[:, 1]
    np.testing.assert_allclose(u_m_top_y, -F_series / k1, atol=1e-8)
    np.testing.assert_allclose(u_s_bot_y, -delta + F_series / k2, atol=1e-8)

    # 1b. Pull upward +0.02 -> clean unilateral separation
    u_presc_up = _prescribed_vector(sys_s["presc_mask"], 0.0, 0.02, 0.0)
    u_fixed_up = np.zeros(box_s.n_dofs, dtype=np.float64)
    u_fixed_up[sys_s["presc_mask"]] = u_presc_up
    res_up = solve_deformable_ssn_contact(
        sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, eps_n,
        free_m=sys_m["free_mask"],
        free_s=sys_s["free_mask"], u_fixed_s=u_fixed_up,
        f_ext_s=-(sys_s["K_fp"] @ u_presc_up),
        u_init_m=res["u_free_m"], u_init_s=res["u_free_s"],
        use_pardiso=False,
    )
    assert res_up["n_active"] == 0
    assert float(np.linalg.norm(res_up["force_total"])) == 0.0
    assert float(np.linalg.norm(res_up["u_free_m"])) < 1e-13
    np.testing.assert_allclose(res_up["u_ctrl_s"][:, 1], 0.02, atol=1e-12)
    print(f"\n[Test 1] F_N={F_N:.6f}, F_series={F_series:.6f}, R_cut_y={R_cut[1]:.6f}, sep_n_act={res_up['n_active']}")


def test_deformable_coulomb_friction_patch_mu():
    """Test 2: deformable-deformable Coulomb friction patch test (mu = 0.3)."""
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    mu, eps_n, eps_t = 0.3, 1e6, 1e6
    delta = 0.01

    # 2a. Comparable soft deformable slabs (E1=1000, E2=2000): exact local Coulomb slip & equilibrium
    box_m_soft = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 0.2, 1.0), n_cells=(2, 1, 2), degree=2)
    box_s_soft = BoxSpline(lo=(0.0, 0.2, 0.0), hi=(1.0, 0.4, 1.0), n_cells=(2, 1, 2), degree=2)
    sys_m_soft = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_m_soft, 1000.0, 0.0,
        zero_fn=lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
    )
    sys_s_soft = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_s_soft, 2000.0, 0.0,
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - 0.4) < 1e-12 else (),
    )
    tn_s, tr_s = _box_face_triangles(box_s_soft, axis=1, side=0)
    Xs_s, wq_s, Bs_s = sample_surface_triangles_gauss(tn_s, tr_s, box_s_soft, order=2)
    mf_s = DeformableBoxFaceMaster(box_m_soft, axis=1, side=1)

    up = _prescribed_vector(sys_s_soft["presc_mask"], 0.0, -delta, 0.0)
    uf = np.zeros(box_s_soft.n_dofs, dtype=np.float64)
    uf[sys_s_soft["presc_mask"]] = up
    res = solve_deformable_ssn_contact(
        sys_m_soft["K_ff"], sys_s_soft["K_ff"], Xs_s, wq_s, Bs_s, mf_s, eps_n,
        free_m=sys_m_soft["free_mask"], free_s=sys_s_soft["free_mask"],
        u_fixed_s=uf, f_ext_s=-(sys_s_soft["K_fp"] @ up),
        mu=mu, eps_t=eps_t, use_pardiso=False,
    )
    um, us, slip_state = res["u_free_m"], res["u_free_s"], res["slip_state"]
    for step in range(1, 21):
        up = _prescribed_vector(sys_s_soft["presc_mask"], 0.002 * step, -delta, 0.0)
        uf = np.zeros(box_s_soft.n_dofs, dtype=np.float64)
        uf[sys_s_soft["presc_mask"]] = up
        res = solve_deformable_ssn_contact(
            sys_m_soft["K_ff"], sys_s_soft["K_ff"], Xs_s, wq_s, Bs_s, mf_s, eps_n,
            free_m=sys_m_soft["free_mask"], free_s=sys_s_soft["free_mask"],
            u_fixed_s=uf, f_ext_s=-(sys_s_soft["K_fp"] @ up),
            u_init_m=um, u_init_s=us, mu=mu, eps_t=eps_t,
            slip_state=slip_state, use_pardiso=False,
        )
        um, us, slip_state = res["u_free_m"], res["u_free_s"], res["slip_state"]
    assert res["n_slip"] == res["n_active"] and res["n_active"] > 0
    act = res["active"]
    t_n_vec = (wq_s[act] * eps_n * res["gap"][act])[:, None] * res["normal"][act]
    t_t_vec = res["traction_res"][act] - t_n_vec
    local_ratios = np.linalg.norm(t_t_vec, axis=1) / np.linalg.norm(t_n_vec, axis=1)
    np.testing.assert_allclose(local_ratios, mu, rtol=1e-10)
    R_cut_soft = _cut_face_reaction(sys_s_soft, us, up, res["f_c_s_full"])
    assert R_cut_soft[0] == pytest.approx(float(res["force_total"][0]), rel=1e-9)

    # 2b. Flat-master deformable patch test (E1=1e9, E2=2000): global |F_T|/F_N == 0.300000
    box_m = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    box_s = BoxSpline(lo=(0.0, 1.0, 0.0), hi=(1.0, 1.25, 1.0), n_cells=(2, 1, 2), degree=2)
    sys_m = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_m, 1.0e9, 0.0,
        zero_fn=lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
    )
    sys_s = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_s, 2000.0, 0.0,
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - 1.25) < 1e-12 else (),
    )
    tri_nodes, tris = _box_face_triangles(box_s, axis=1, side=0)
    X_s, w_q, B_s = sample_surface_triangles_gauss(tri_nodes, tris, box_s, order=2)
    master_face = DeformableBoxFaceMaster(box_m, axis=1, side=1)

    u_presc = _prescribed_vector(sys_s["presc_mask"], 0.0, -0.05, 0.0)
    u_fixed_s = np.zeros(box_s.n_dofs, dtype=np.float64)
    u_fixed_s[sys_s["presc_mask"]] = u_presc
    res = solve_deformable_ssn_contact(
        sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, 1e6,
        free_m=sys_m["free_mask"], free_s=sys_s["free_mask"],
        u_fixed_s=u_fixed_s, f_ext_s=-(sys_s["K_fp"] @ u_presc),
        mu=mu, eps_t=1e5, use_pardiso=False,
    )
    assert res["n_active"] == len(w_q)
    um, us, slip_state = res["u_free_m"], res["u_free_s"], res["slip_state"]

    ratios = []
    slipped_any = False
    for step in range(1, 26):
        ux = 0.002 * step
        u_presc = _prescribed_vector(sys_s["presc_mask"], ux, -0.05, 0.0)
        u_fixed_s = np.zeros(box_s.n_dofs, dtype=np.float64)
        u_fixed_s[sys_s["presc_mask"]] = u_presc
        res = solve_deformable_ssn_contact(
            sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, 1e6,
            free_m=sys_m["free_mask"], free_s=sys_s["free_mask"],
            u_fixed_s=u_fixed_s, f_ext_s=-(sys_s["K_fp"] @ u_presc),
            u_init_m=um, u_init_s=us, mu=mu, eps_t=1e5,
            slip_state=slip_state, use_pardiso=False,
        )
        um, us, slip_state = res["u_free_m"], res["u_free_s"], res["slip_state"]
        slipped_any |= (res["n_slip"] > 0)
        fn = float(-res["force_total"][1])
        ft = float(np.hypot(res["force_total"][0], res["force_total"][2]))
        ratios.append(ft / fn)

    R_cut = _cut_face_reaction(sys_s, us, u_presc, res["f_c_s_full"])
    print(f"\n[Test 2] stick ratio[0]={ratios[0]:.6f}, plateau ratio[-1]={ratios[-1]:.8f}, "
          f"n_slip={res['n_slip']}/{len(w_q)}, R_cut_x={R_cut[0]:.6f}, F_T={res['force_total'][0]:.6f}")
    assert slipped_any
    assert ratios[0] < mu
    assert res["n_slip"] == len(w_q)
    assert all(r <= mu * (1.0 + 1e-5) for r in ratios)
    assert ratios[-1] == pytest.approx(mu, rel=1e-5)
    assert R_cut[0] == pytest.approx(float(res["force_total"][0]), rel=1e-10)
    assert R_cut[1] == pytest.approx(float(res["force_total"][1]), rel=1e-10)


def _setup_v08_half(KM, SM, inp: H.V08Inputs, h: float, n_soup: int = 24):
    """Build trimmed deformable indenter + deformable block systems for V08 sym='half'."""
    lo_b, hi_b, nc_b = H.v08_block_grid(h, sym="half")
    box_blk = BoxSpline(lo_b, hi_b, nc_b, degree=2)

    def _blk_zero(x):
        if abs(x[1] - lo_b[1]) < 1e-12:
            return (0, 1, 2)
        if abs(x[2] - lo_b[2]) < 1e-12:
            return (2,)
        return ()

    sys_m = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_blk, inp.E_block, inp.nu_block, zero_fn=_blk_zero
    )

    nodes_ind, tris_ind, fids_ind = H.hollow_hemisphere_soup(
        inp.R, inp.Ri, inp.center, n=n_soup, sym="half"
    )
    lo_i, hi_i, nc_i = H.v08_indenter_grid(h, inp.center, inp.R, sym="half")
    box_ind = BoxSpline(lo_i, hi_i, nc_i, degree=2)

    sys_s = _assemble_kratos_partitioned_stiffness(
        KM, SM, box_ind, inp.E_indenter, inp.nu_indenter,
        zero_fn=lambda x: (2,) if abs(x[2] - lo_i[2]) < 1e-12 else (),
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - hi_i[1]) < 1e-12 else (),
        trim_soup=(nodes_ind, tris_ind),
        q=4,
        ghost_alpha=1e-8,
    )
    assert sys_s["diag"]["n_fallback"] == 0

    tris_outer = tris_ind[fids_ind == H.FACE_OUTER]
    X_s, w_q, B_s = sample_surface_triangles_gauss(nodes_ind, tris_outer, box_ind, order=2)
    master_face = DeformableBoxFaceMaster(box_blk, axis=1, side=1)
    return box_blk, box_ind, sys_m, sys_s, X_s, w_q, B_s, master_face


def test_v08_case_r1_deformable_half():
    """Test 3: V08 Case R1 (deformable trimmed hollow hemisphere on deformable block, sym='half')."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    F_hertz = H.v08_reference(inp)["normal_force_F_N_lbf"]
    F_sneddon = H.v08_model_form_gaps(inp)["exact_F_N_lbf"]

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    t_steps = np.linspace(0.03, 0.1, 8)
    print(f"\n[Test 3] V08 Case R1 Deformable-Deformable (Hertz={F_hertz:,.1f} lbf, Sneddon={F_sneddon:,.1f} lbf):")

    results_by_h = {}
    for h in (0.5, 0.25):
        t0 = time.perf_counter()
        _bb, box_ind, sys_m, sys_s, X_s, w_q, B_s, master_face = _setup_v08_half(KM, SM, inp, h)
        t_setup = time.perf_counter() - t0

        um = np.zeros(int(sys_m["free_mask"].sum()), dtype=np.float64)
        us = np.zeros(int(sys_s["free_mask"].sum()), dtype=np.float64)
        eps_n = 1e9
        t_s = time.perf_counter()
        curve = []
        for t in t_steps:
            ux, uy = H.v08_cut_face_displacement(float(t))
            u_presc = _prescribed_vector(sys_s["presc_mask"], float(ux), float(uy), 0.0)
            u_fixed_s = np.zeros(box_ind.n_dofs, dtype=np.float64)
            u_fixed_s[sys_s["presc_mask"]] = u_presc
            res = solve_deformable_ssn_contact(
                sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, eps_n,
                free_m=sys_m["free_mask"],
                free_s=sys_s["free_mask"], u_fixed_s=u_fixed_s,
                f_ext_s=-(sys_s["K_fp"] @ u_presc),
                u_init_m=um, u_init_s=us,
                rtol=1e-9, use_pardiso=True,
            )
            um, us = res["u_free_m"], res["u_free_s"]
            fn_iface = 2.0 * float(-res["force_total"][1])
            R_cut = _cut_face_reaction(sys_s, us, u_presc, res["f_c_s_full"])
            fn_cut = 2.0 * float(-R_cut[1])
            curve.append((float(H.v08_delta(float(t), inp.g0)), res["iterations"], res["n_active"], fn_iface, fn_cut))
        dt_solve = time.perf_counter() - t_s

        fn_iface_final = curve[-1][3]
        fn_cut_final = curve[-1][4]
        eq_rel = abs(fn_cut_final - fn_iface_final) / fn_iface_final
        err_h = (fn_iface_final / F_hertz - 1.0) * 100.0
        err_s = (fn_iface_final / F_sneddon - 1.0) * 100.0
        results_by_h[h] = (fn_iface_final, fn_cut_final, err_h, err_s, curve)

        print(
            f"  h={h:.2f} (dofs={int(sys_m['free_mask'].sum())}+{int(sys_s['free_mask'].sum())}, "
            f"n_fb={sys_s['diag']['n_fallback']}, setup={t_setup:.2f}s, solve={dt_solve:.2f}s): "
            f"it={curve[-1][1]}, n_act={curve[-1][2]}, "
            f"F_N_iface={fn_iface_final:,.1f} lbf, F_N_cut={fn_cut_final:,.1f} lbf (eq_rel={eq_rel:.2e}), "
            f"vs Hertz {err_h:+.2f}%, vs Sneddon {err_s:+.2f}%"
        )
        assert eq_rel < 1e-4
        assert max(r[1] for r in curve) <= 12

    # NOT a G5 pass (S41): at h=0.25, n=24, eps_n=1e9 the soft-penalty penetration
    # (about 5 % of delta), the coarse soup and the coarse mesh cancel. With
    # eps_n -> inf and n -> inf the small-strain result at h=0.125 is about +5 %.
    # Regression pin of this configuration only (measured 2026-11-18/19):
    assert results_by_h[0.25][0] == pytest.approx(486248.6, rel=1e-5)


def test_v08_case_r2_sliding_friction_half():
    """Test 4: V08 Case R2 (deformable-deformable sliding friction mu=0.3, sym='half').

    Caveat (S41): runs at eps_n=1e9 and soup n=24 like R1, whose h=0.25 pass is an
    error cancellation. In the slip window F_T = mu F_N, so the RMS gates below are
    expected to inherit the same penalty dependence (not yet measured for R2).
    """
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    print("\n[Test 4] V08 Case R2 Sliding Friction (mu=0.3, t=0->0.40, window 0.15<=t<=0.40):")
    configs = (
        (0.5, 1e9, 1e8, np.linspace(0.03, 0.1, 8), np.linspace(0.112, 0.40, 25), False),
        (0.25, 1e9, 1e8, np.linspace(0.03, 0.1, 6), np.linspace(0.115, 0.40, 20), True),
    )
    for h, eps_n, eps_t, t_indent, t_slide, use_p in configs:
        t_all = np.concatenate([t_indent, t_slide])
        t0 = time.perf_counter()
        _bb, box_ind, sys_m, sys_s, X_s, w_q, B_s, master_face = _setup_v08_half(KM, SM, inp, h)
        t_setup = time.perf_counter() - t0

        um = np.zeros(int(sys_m["free_mask"].sum()), dtype=np.float64)
        us = np.zeros(int(sys_s["free_mask"].sum()), dtype=np.float64)
        slip_state = None
        t_s = time.perf_counter()

        rows = []
        for t in t_all:
            ux, uy = H.v08_cut_face_displacement(float(t))
            u_presc = _prescribed_vector(sys_s["presc_mask"], float(ux), float(uy), 0.0)
            u_fixed_s = np.zeros(box_ind.n_dofs, dtype=np.float64)
            u_fixed_s[sys_s["presc_mask"]] = u_presc
            res = solve_deformable_ssn_contact(
                sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, eps_n,
                free_m=sys_m["free_mask"],
                free_s=sys_s["free_mask"], u_fixed_s=u_fixed_s,
                f_ext_s=-(sys_s["K_fp"] @ u_presc),
                u_init_m=um, u_init_s=us,
                mu=inp.mu, eps_t=eps_t, slip_state=slip_state,
                rtol=1e-8, use_pardiso=use_p,
            )
            um, us, slip_state = res["u_free_m"], res["u_free_s"], res["slip_state"]
            fn_iface = 2.0 * float(-res["force_total"][1])
            ft_iface = 2.0 * float(res["force_total"][0])
            R_cut = _cut_face_reaction(sys_s, us, u_presc, res["f_c_s_full"])
            fn_cut = 2.0 * float(-R_cut[1])
            ft_cut = 2.0 * float(R_cut[0])
            rows.append((float(t), res["iterations"], res["n_active"], res["n_slip"], fn_iface, ft_iface, fn_cut, ft_cut))
        dt_solve = time.perf_counter() - t_s

        t_arr = np.array([r[0] for r in rows])
        fn_arr = np.array([r[4] for r in rows])
        ft_arr = np.array([r[5] for r in rows])
        fn_cut_arr = np.array([r[6] for r in rows])
        ft_cut_arr = np.array([r[7] for r in rows])

        win = (t_arr >= 0.15) & (t_arr <= 0.40)
        ref = H.v08_reference_history(t_arr[win], inp)
        rms_fn = H.relative_rms(fn_arr[win], ref["F_N"])
        rms_ft = H.relative_rms(ft_arr[win], ref["F_T"])
        ratios_win = ft_arr[win] / fn_arr[win]
        eq_fn = float(np.max(np.abs(fn_cut_arr[win] - fn_arr[win]) / fn_arr[win]))
        eq_ft = float(np.max(np.abs(ft_cut_arr[win] - ft_arr[win]) / ft_arr[win]))

        print(
            f"  h={h:.2f} (setup={t_setup:.2f}s, solve={dt_solve:.2f}s, max_it={max(r[1] for r in rows)}): "
            f"F_N_win={fn_arr[win].mean():,.1f} lbf, F_T_win={ft_arr[win].mean():,.1f} lbf, "
            f"F_T/F_N in [{ratios_win.min():.4f}, {ratios_win.max():.4f}] (mean={ratios_win.mean():.4f}), "
            f"eq_err=(N:{eq_fn:.2e}, T:{eq_ft:.2e}), "
            f"rel_RMS=(F_N:{rms_fn*100:.2f}%, F_T:{rms_ft*100:.2f}%)"
        )
        assert eq_fn < 1e-4
        assert eq_ft < 1e-4
        assert abs(float(ratios_win.mean()) - inp.mu) < 0.015
        if h == 0.5:
            assert abs(float(ratios_win[-1]) - inp.mu) < 0.035
        else:
            assert abs(float(ratios_win[-1]) - inp.mu) < 0.015
            assert rms_fn < 0.045
            assert rms_ft < 0.045


def test_v08_r1_penalty_converged_small_strain_misses_gate():
    """S41: with a converged penalty the small-strain R1 at h=0.25 misses the 4 % gate.

    h=0.25, soup n=48: eps_n 1e10 -> 1e11 changes F_N by 0.5 % (penalty converged
    to < 1 %), and F_N is +6.6 % above Hertz. Pins the honest status of the
    small-strain SSN path; the S32 pass at eps_n=1e9 is penetration-softened.
    Measured 2026-11-19: 510,891.8 (1e10), 513,553.7 (1e11) lbf.
    """
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))
    F_hertz = H.v08_reference(inp)["normal_force_F_N_lbf"]
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    _bb, box_ind, sys_m, sys_s, X_s, w_q, B_s, master_face = _setup_v08_half(KM, SM, inp, 0.25, n_soup=48)
    fn = {}
    for eps_n in (1e10, 1e11):
        um = np.zeros(int(sys_m["free_mask"].sum()), dtype=np.float64)
        us = np.zeros(int(sys_s["free_mask"].sum()), dtype=np.float64)
        for t in np.linspace(0.03, 0.1, 8):
            ux, uy = H.v08_cut_face_displacement(float(t))
            u_presc = _prescribed_vector(sys_s["presc_mask"], float(ux), float(uy), 0.0)
            u_fixed_s = np.zeros(box_ind.n_dofs, dtype=np.float64)
            u_fixed_s[sys_s["presc_mask"]] = u_presc
            res = solve_deformable_ssn_contact(
                sys_m["K_ff"], sys_s["K_ff"], X_s, w_q, B_s, master_face, eps_n,
                free_m=sys_m["free_mask"], free_s=sys_s["free_mask"], u_fixed_s=u_fixed_s,
                f_ext_s=-(sys_s["K_fp"] @ u_presc), u_init_m=um, u_init_s=us,
                rtol=1e-9, use_pardiso=True,
            )
            um, us = res["u_free_m"], res["u_free_s"]
        fn[eps_n] = 2.0 * float(-res["force_total"][1])
    assert abs(fn[1e11] / fn[1e10] - 1.0) < 0.01  # penalty converged
    err = fn[1e11] / F_hertz - 1.0
    assert 0.06 < err < 0.075, err  # small-strain, h=0.25: gate (4 %) NOT met
