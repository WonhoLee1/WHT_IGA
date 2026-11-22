# -*- coding: utf-8 -*-
"""Unit tests for F-bar mean-dilatation pipeline and User API wiring (Kratos-free)."""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError

import wht_iga.spec as s
from wht_iga.analyses import fbar as FB
from wht_iga.analyses.methods import (
    BoxMethodProblem,
    FullIgaOption,
    parse_discretization_option,
    recommend_discretization,
    solve_fbar_box,
    solve_with_discretization_method,
)
from wht_iga.backends import ReferenceNumbaBackend
from wht_iga.plan import compile_plan
from wht_iga.reference.assemble_elastic import assemble_face_traction
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.cook import CookSpline, HR_MM, L_MM, T_MM
from wht_iga.runner import AnalysisSession


def test_recommend_discretization_surfaces_fbar_for_near_incompressible():
    """nu >= 0.49 recommends locking_mitigation='fbar' with reason and warning."""
    r_default = recommend_discretization("general")
    assert r_default.locking_mitigation == "none"

    r_comp = recommend_discretization("bulk", nu=0.3)
    assert r_comp.locking_mitigation == "none"

    r_edge = recommend_discretization("bending", nu=0.49)
    assert r_edge.locking_mitigation == "fbar"
    assert any("F-bar" in reason for reason in r_edge.reasons)
    assert any("locking_mitigation='fbar'" in w for w in r_edge.warnings)

    r_near = recommend_discretization("general", nu=0.49999)
    assert r_near.locking_mitigation == "fbar"

    with pytest.raises(ValueError, match="nu"):
        recommend_discretization("general", nu=0.5)
    with pytest.raises(ValueError, match="nu"):
        recommend_discretization("general", nu=-1.0)


def test_full_iga_option_locking_mitigation_validation():
    """FullIgaOption validates locking_mitigation and round-trips JSON."""
    opt = FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(3, 2, 2), locking_mitigation="fbar")
    parsed = parse_discretization_option(opt.model_dump_json())
    assert parsed == opt

    # Since 2026-11-20 immersed F-bar routes to fbar_trimmed; only configurations outside its scope fail.
    FullIgaOption(kind="full_iga_immersed", locking_mitigation="fbar")
    with pytest.raises(ValidationError, match="out_mode"):
        FullIgaOption(kind="full_iga_immersed", locking_mitigation="fbar", out_mode="ghost")


@pytest.mark.parametrize("mitigation", ["fbar", "auto"])
def test_methods_pipeline_fbar_box_matches_fbar_solve_static(mitigation):
    """solve_with_discretization_method routes nu=0.49999 BoxSpline to fbar.solve_static."""
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(4.0, 1.0, 1.0),
        E=240.565,
        nu=0.49999,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(0.0, 0.5, 0.0),
    )
    opt = FullIgaOption(
        kind="full_iga_body_fit",
        degree=2,
        continuity=0,
        n_cells=(3, 2, 2),
        locking_mitigation=mitigation,
    )
    sol = solve_with_discretization_method(problem, opt, use_pardiso=False)
    assert sol.converged
    assert sol.native["locking_mitigation"] == "fbar"

    box = BoxSpline(lo=problem.lo, hi=problem.hi, n_cells=opt.n_cells, degree=2, continuity=0)
    model = FB.build_model(box, 2, E=problem.E, nu=problem.nu)
    f_ref = assemble_face_traction(box, 0, 1, np.array(problem.traction))
    fixed_cp = box.boundary_control_points(0, 0)
    fixed_ref = np.array([3 * int(i) + d for i in fixed_cp for d in range(3)], dtype=np.int64)
    u_ref = FB.solve_static(model, f_ref, fixed_ref, n_steps=1, linear_solver="superlu")

    assert np.allclose(sol.native["u"], u_ref, rtol=1e-12, atol=1e-12)
    tip_ref = box.evaluate(u_ref.reshape(-1, 3), problem.tip_point().reshape(1, 3))[0]
    assert np.allclose(sol.tip_displacement, tip_ref, rtol=1e-12, atol=1e-12)


def test_solve_fbar_box_cook_spline_c_p_minus_2_matches_fbar_direct():
    """solve_fbar_box solves a C^(p-2) CookSpline at nu=0.49999 and matches fbar.solve_static."""
    p = 2
    n_cells = (4, 4, 1)
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=p, continuity=0)
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(1.0, 1.0, 1.0),
        E=240.565,
        nu=0.49999,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(0.0, 6.25, 0.0),
    )
    sol = solve_fbar_box(
        problem,
        n_steps=4,
        use_pardiso=False,
        box=box,
        fixed_axes=(2,),
    )
    assert sol.converged
    assert sol.native["locking_mitigation"] == "fbar"

    # Direct fbar.solve_static oracle using exact consistent Cook right-face traction
    model = FB.build_model(box, p, E=problem.E, nu=problem.nu)
    cp = np.asarray(box.control_points())
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    kv, kw = (np.asarray(k, dtype=float) for k in box.knots[1:])
    a = np.array([(kv[j + p + 1] - kv[j]) / (p + 1) for j in range(nv_)])
    c = np.array([(kw[k + p + 1] - kw[k]) / (p + 1) for k in range(nw_)])
    f_ref = np.zeros(model["n_dof"])
    for k in range(nw_):
        for j in range(nv_):
            ctrl_id = (nu_ - 1) + nu_ * (j + nv_ * k)
            assert abs(cp[ctrl_id, 0] - L_MM) < 1e-9
            f_ref[3 * ctrl_id + 1] = 6.25 * HR_MM * T_MM * a[j] * c[k]

    fixed = set()
    for ctrl_id in range(cp.shape[0]):
        fixed.add(3 * ctrl_id + 2)
        if abs(cp[ctrl_id, 0]) < 1e-9:
            fixed.add(3 * ctrl_id)
            fixed.add(3 * ctrl_id + 1)
    u_ref = FB.solve_static(
        model,
        f_ref,
        np.array(sorted(fixed), dtype=np.int64),
        n_steps=4,
        linear_solver="superlu",
    )

    assert np.allclose(sol.native["f_ext"], f_ref, rtol=1e-12, atol=1e-12)
    assert np.allclose(sol.native["u"], u_ref, rtol=1e-12, atol=1e-12)
    stresses = FB.cauchy_stress(sol.native["model"], sol.native["u"])
    assert len(stresses) == n_cells[0] * n_cells[1] * n_cells[2]
    assert all(np.all(np.isfinite(s_cell)) for s_cell in stresses)


@pytest.mark.parametrize("mitigation", ["fbar", "auto"])
def test_user_api_study_spec_fbar_wiring(tmp_path, mitigation):
    """StudySpec -> compile_plan -> ReferenceNumbaBackend routes F-bar cleanly."""
    study = s.StudySpec(
        name=f"fbar_bar_{mitigation}",
        units="mm-N-tonne-s",
        parts=[s.PartSpec.from_box((0.0, 0.0, 0.0), (4.0, 1.0, 1.0), "rubber")],
        materials={"rubber": s.LinearElasticSpec(E=240.565, nu=0.49999)},
        mesh=s.BodyFitSpec(degree=2, continuity=1, n_cells=(3, 2, 2)),
        steps=[
            s.StaticStepSpec(
                name="bend",
                bcs=[s.DisplacementBCSpec(set="xmin", ux=0.0, uy=0.0, uz=0.0)],
                loads=[s.TractionLoadSpec(set="xmax", vector=(0.0, 0.5, 0.0))],
                locking_mitigation=mitigation,
            )
        ],
    )
    plan = compile_plan(study, output_dir=tmp_path / study.name)
    with AnalysisSession(plan, backend=ReferenceNumbaBackend()) as ses:
        res = ses.run()

    step_res = res.step("bend")
    assert step_res.converged
    assert "fbar.solve_static" in step_res.solver_log

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(4.0, 1.0, 1.0), n_cells=(3, 2, 2), degree=2, continuity=1)
    model = FB.build_model(box, 2, E=240.565, nu=0.49999)
    f_ref = assemble_face_traction(box, 0, 1, np.array([0.0, 0.5, 0.0]))
    fixed_cp = box.boundary_control_points(0, 0)
    fixed_ref = np.array([3 * int(i) + d for i in fixed_cp for d in range(3)], dtype=np.int64)
    u_ref = FB.solve_static(model, f_ref, fixed_ref, linear_solver="auto")
    assert np.allclose(step_res.U, u_ref, rtol=1e-12, atol=1e-12)
