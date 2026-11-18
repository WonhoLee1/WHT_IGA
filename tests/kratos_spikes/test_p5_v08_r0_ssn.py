# -*- coding: utf-8 -*-
"""Priority b-1: Gauss-weighted Python SSN contact solver & V08 Case R0 verification.

Tests:
1. Unit/analytic flat rigid punch on a p=2 C^1 block via Kratos stiffness K.
2. Unilateral separation (sphere press -> lift -> zero active GPs & zero force).
3. V08 Case R0 (1/4 p=2 C^1 block, rigid sphere R=0.5, E' = (1-nu^2)E*, delta=0->0.1):
   penalty convergence (eps_n = 1e8, 1e9, 1e10), multi-increment curve & delta^1.5 scaling,
   h0/h1/h2 exact comparisons vs Hertz (481,579.26 lbf) and Sneddon (471,927.70 lbf).
4. Coulomb friction patch test (mu = 0.3): normal press + tangential slide -> |F_T|/F_N = 0.3000.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import (
    RigidPlaneObstacle,
    RigidSphereObstacle,
    bspline_basis_matrix,
    expand_vector_basis,
    sample_box_face_gauss,
    solve_ssn_contact,
)
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    solve_linear,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"


def _assemble_kratos_stiffness(
    KM,
    SM,
    box: BoxSpline,
    E: float,
    nu: float,
    fix_fn,
    *,
    use_pardiso: bool = True,
):
    """Build Kratos IGA model part, extract free-DOF stiffness Kf (in BoxSpline order)."""
    model = KM.Model()
    mp = model.CreateModelPart("SSN")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    attach_cell_elements(KM, mp, vol, box, props, element_name="SmallDisplacementElement3D4N")

    Xcp = box.control_points()
    fixed = np.zeros(box.n_dofs, dtype=bool)
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for I in range(box.n_ctrl):
        nd = mp.GetNode(I + 1)
        for d in fix_fn(Xcp[I]):
            fixed[3 * I + d] = True
            nd.Fix(comp_vars[d])

    strat = solve_linear(KM, mp, use_pardiso=use_pardiso)
    Kkr = system_matrix_csr(strat).tocsr()
    perm = equation_permutation(KM, mp, box)
    free = ~fixed
    free_eq = perm[np.flatnonzero(free)]
    Kf = Kkr[free_eq, :][:, free_eq].tocsr()
    return Kf, free


def _v08_quarter_fix(lo):
    def _fix(x):
        dofs = []
        if abs(x[1] - lo[1]) < 1e-12:
            return (0, 1, 2)
        if abs(x[0] - lo[0]) < 1e-12:
            dofs.append(0)
        if abs(x[2] - lo[2]) < 1e-12:
            dofs.append(2)
        return dofs

    return _fix


def test_flat_rigid_punch_p2_analytic():
    """Test 1: uniform indentation on p=2 C1 block matches analytic uniaxial stress & displacement."""
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    E, nu, delta = 1000.0, 0.3, 0.05
    H_blk, A_blk = 1.0, 1.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)

    def _fix_uniaxial(x):
        dofs = []
        if abs(x[1]) < 1e-12:
            dofs.append(1)
        if abs(x[0]) < 1e-12:
            dofs.append(0)
        if abs(x[2]) < 1e-12:
            dofs.append(2)
        return dofs

    Kf, free = _assemble_kratos_stiffness(KM, SM, box, E, nu, _fix_uniaxial, use_pardiso=False)
    X_q, w_q, B_scalar = sample_box_face_gauss(box, axis=1, side=1, n_gp=3)

    # Verify partition of unity, vector basis shape, and exact face quadrature area.
    row_sums = np.asarray(B_scalar.sum(axis=1)).ravel()
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-14)
    assert expand_vector_basis(B_scalar, 3).shape == (3 * len(w_q), box.n_dofs)
    assert float(w_q.sum()) == pytest.approx(A_blk, rel=1e-14)

    eps_n = 1e8
    punch = RigidPlaneObstacle(point=(0.0, 1.0 - delta, 0.0), normal=(0.0, -1.0, 0.0))
    res = solve_ssn_contact(Kf, X_q, w_q, B_scalar, punch, eps_n, free_mask=free, use_pardiso=False)

    F_N = float(res["force_total"][1])
    k_blk = E * A_blk / H_blk
    k_pen = eps_n * A_blk
    F_series = delta / (1.0 / k_blk + 1.0 / k_pen)
    F_rigid = E * A_blk * delta / H_blk

    assert res["n_active"] == len(w_q)
    assert F_N == pytest.approx(F_series, rel=1e-11)
    assert abs(F_N - F_rigid) / F_rigid < 2.0e-5  # O(E / (eps_n * H))

    # Surface displacement evaluated through B-spline basis matches exact 3D Poisson field.
    u_surf = B_scalar @ res["u_ctrl"]
    eps_yy = -F_series / (E * A_blk)
    np.testing.assert_allclose(u_surf[:, 1], eps_yy * H_blk, atol=1e-11)
    np.testing.assert_allclose(u_surf[:, 0], -nu * eps_yy * X_q[:, 0], atol=1e-11)
    np.testing.assert_allclose(u_surf[:, 2], -nu * eps_yy * X_q[:, 2], atol=1e-11)
    print(f"\n[Test 1] F_N={F_N:.6f}, F_series={F_series:.6f}, F_rigid={F_rigid:.6f}, it={res['iterations']}")


def test_unilateral_sphere_separation():
    """Test 2: unilateral contact engages on indentation and drops cleanly to 0 on lift."""
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(3, 3, 3), degree=2)

    def _fix(x):
        dofs = []
        if abs(x[1]) < 1e-12:
            return (0, 1, 2)
        if abs(x[0]) < 1e-12:
            dofs.append(0)
        if abs(x[2]) < 1e-12:
            dofs.append(2)
        return dofs

    Kf, free = _assemble_kratos_stiffness(KM, SM, box, 1000.0, 0.3, _fix, use_pardiso=False)
    X_q, w_q, B_scalar = sample_box_face_gauss(box, axis=1, side=1, n_gp=3, bounds=((0.0, 0.6), (0.0, 0.6)), n_sub=2)

    R = 0.5
    sphere_in = RigidSphereObstacle(center=(0.0, 1.0 + R - 0.05, 0.0), radius=R)
    res_in = solve_ssn_contact(Kf, X_q, w_q, B_scalar, sphere_in, eps_n=1e5, free_mask=free, use_pardiso=False)
    assert res_in["n_active"] > 0
    assert res_in["force_total"][1] > 0.0

    # Lift sphere above flat undeformed surface (gap = +0.05 > 0)
    sphere_out = RigidSphereObstacle(center=(0.0, 1.0 + R + 0.05, 0.0), radius=R)
    res_out = solve_ssn_contact(
        Kf, X_q, w_q, B_scalar, sphere_out, eps_n=1e5,
        free_mask=free, u_init=res_in["u_free"], use_pardiso=False,
    )
    assert res_out["n_active"] == 0
    assert res_out["force_total"][1] == 0.0
    assert float(np.linalg.norm(res_out["u_free"])) < 1e-13
    print(f"\n[Test 2] press n_act={res_in['n_active']}, F_N={res_in['force_total'][1]:.4f}; "
          f"lift n_act={res_out['n_active']}, F_N={res_out['force_total'][1]:.1f}")


def test_v08_case_r0_ssn():
    """Test 3: V08 Case R0 (quarter p=2 C1 block, rigid sphere R=0.5, delta=0->0.1)."""
    with open(_V08_JSON, encoding="utf-8") as f:
        inp = H.V08Inputs.from_json(json.load(f))

    Eb = H.rigid_equivalent_block_modulus(inp.E_star, inp.nu_block)
    F_hertz = H.v08_reference(inp)["normal_force_F_N_lbf"]
    F_sneddon = H.v08_model_form_gaps(inp)["exact_F_N_lbf"]

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    # 3a. Penalty convergence at h0 = 0.5 and h1 = 0.25 (multi-step warm-started)
    deltas = np.linspace(0.0125, 0.1, 8)
    print("\n[Test 3] V08 Case R0 Measurements (Hertz F_N=481,579.3 lbf, Sneddon F_N=471,927.7 lbf):")
    h1_curve = None
    for h in (0.5, 0.25):
        t0 = time.perf_counter()
        lo, hi, nc = H.v08_block_grid(h, sym="quarter")
        box = BoxSpline(lo, hi, nc, degree=2)
        Kf, free = _assemble_kratos_stiffness(KM, SM, box, Eb, inp.nu_block, _v08_quarter_fix(lo), use_pardiso=True)
        X_q, w_q, B_scalar = sample_box_face_gauss(
            box, axis=1, side=1, n_gp=3, bounds=((-1.5, -1.0), (0.0, 0.5)), n_sub=4
        )
        t_setup = time.perf_counter() - t0

        fn_by_eps = []
        for eps_n in (1e8, 1e9, 1e10):
            t_s = time.perf_counter()
            u = np.zeros(int(free.sum()))
            step_rows = []
            for d in deltas:
                sph = RigidSphereObstacle(center=(-1.5, 1.0 + inp.R - float(d), 0.0), radius=inp.R)
                res = solve_ssn_contact(
                    Kf, X_q, w_q, B_scalar, sph, eps_n=eps_n,
                    free_mask=free, u_init=u, rtol=1e-9, use_pardiso=True,
                )
                u = res["u_free"]
                step_rows.append((float(d), res["iterations"], res["n_active"], 4.0 * float(res["force_total"][1])))
            dt_solve = time.perf_counter() - t_s
            fn_final = step_rows[-1][3]
            fn_by_eps.append(fn_final)
            err_h = (fn_final / F_hertz - 1.0) * 100.0
            err_s = (fn_final / F_sneddon - 1.0) * 100.0
            print(f"  h={h:.3f} (dofs={int(free.sum())}, setup={t_setup:.2f}s) eps={eps_n:.0e} ({dt_solve:.2f}s): "
                  f"it={step_rows[-1][1]}, n_act={step_rows[-1][2]}, F_N={fn_final:,.1f} lbf "
                  f"(vs Hertz {err_h:+.2f}%, vs Sneddon {err_s:+.2f}%)")
            if h == 0.25 and eps_n == 1e9:
                h1_curve = step_rows

        # Monotonic convergence with increasing penalty stiffness
        assert fn_by_eps[0] < fn_by_eps[1] < fn_by_eps[2]

    # 3b. Verify h1 = 0.25 curve shape (~ delta^1.5) and Newton iteration efficiency
    assert h1_curve is not None
    d_arr = np.array([r[0] for r in h1_curve])
    fn_arr = np.array([r[3] for r in h1_curve])
    slope, _ = np.polyfit(np.log(d_arr), np.log(fn_arr), 1)
    max_it = max(r[1] for r in h1_curve)
    print(f"  h1=0.25 (eps=1e9) log-log slope F_N vs delta = {slope:.3f} (Hertz=1.500), max Newton it/step = {max_it}")
    assert 1.30 <= slope <= 1.60
    assert max_it <= 12

    # 3c. Single-point h2 = 0.125 check at delta=0.1 (eps=3e8 & eps=1e9; fast <10s, <1GB RSS)
    t0 = time.perf_counter()
    lo2, hi2, nc2 = H.v08_block_grid(0.125, sym="quarter")
    box2 = BoxSpline(lo2, hi2, nc2, degree=2)
    Kf2, free2 = _assemble_kratos_stiffness(KM, SM, box2, Eb, inp.nu_block, _v08_quarter_fix(lo2), use_pardiso=True)
    X_q2, w_q2, B_scalar2 = sample_box_face_gauss(
        box2, axis=1, side=1, n_gp=3, bounds=((-1.5, -1.0), (0.0, 0.5)), n_sub=2
    )
    sph2 = RigidSphereObstacle(center=(-1.5, 1.0 + inp.R - 0.1, 0.0), radius=inp.R)
    for eps_n in (3e8, 1e9):
        res2 = solve_ssn_contact(
            Kf2, X_q2, w_q2, B_scalar2, sph2, eps_n=eps_n,
            free_mask=free2, rtol=1e-9, use_pardiso=True,
        )
        fn2 = 4.0 * float(res2["force_total"][1])
        print(f"  h=0.125 (dofs={int(free2.sum())}) eps={eps_n:.0e}: it={res2['iterations']}, "
              f"n_act={res2['n_active']}, F_N={fn2:,.1f} lbf "
              f"(vs Hertz {(fn2/F_hertz-1)*100:+.2f}%, vs Sneddon {(fn2/F_sneddon-1)*100:+.2f}%)")
        assert res2["iterations"] <= 12
    print(f"  h2=0.125 total elapsed = {time.perf_counter() - t0:.2f}s")


def test_coulomb_friction_patch_mu():
    """Test 4: Coulomb friction patch test (mu = 0.3) via solve_ssn_contact."""
    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    Kf, free = _assemble_kratos_stiffness(
        KM, SM, box, 1000.0, 0.3,
        lambda x: (0, 1, 2) if abs(x[1]) < 1e-12 else (),
        use_pardiso=False,
    )
    X_q, w_q, B_scalar = sample_box_face_gauss(box, axis=1, side=1, n_gp=3)

    mu = 0.3
    eps_n, eps_t = 1e6, 1e5
    punch = RigidPlaneObstacle(point=(0.0, 1.0 - 0.05, 0.0), normal=(0.0, -1.0, 0.0))

    # Step 1: pure normal indentation
    res = solve_ssn_contact(
        Kf, X_q, w_q, B_scalar, punch, eps_n,
        free_mask=free, mu=mu, eps_t=eps_t, use_pardiso=False,
    )
    u = res["u_free"]
    slip_state = res["slip_state"]
    assert res["n_active"] == len(w_q)
    assert res["force_total"][1] > 0.0

    # Steps 2..51: incremental tangential slide of punch along +x (dx_total = 0.20)
    ratios = []
    slipped_any = False
    for _ in range(50):
        res = solve_ssn_contact(
            Kf, X_q, w_q, B_scalar, punch, eps_n,
            free_mask=free, u_init=u, mu=mu, eps_t=eps_t,
            slip_state=slip_state, obstacle_disp=(0.004, 0.0, 0.0),
            use_pardiso=False,
        )
        u = res["u_free"]
        slip_state = res["slip_state"]
        slipped_any |= (res["n_slip"] > 0)
        fn = float(res["force_total"][1])
        ft = float(np.hypot(res["force_total"][0], res["force_total"][2]))
        ratios.append(ft / fn)

    print(f"\n[Test 4] stick ratio[0]={ratios[0]:.6f}, plateau ratio[-1]={ratios[-1]:.8f}, n_slip={res['n_slip']}/{len(w_q)}")
    assert slipped_any
    assert ratios[0] < mu
    assert res["n_slip"] == len(w_q)
    assert all(r <= mu * (1.0 + 1e-12) for r in ratios)
    assert ratios[-1] == pytest.approx(mu, rel=1e-7)
