"""User API (plan v2.1 §8, D18): spec / plan / ports / runner / compat, Kratos-free.

Numerical checks compare against closed-form solutions that the B-spline
spaces reproduce exactly (uniaxial bar: linear field; gravity column with
nu = 0: quadratic field at p = 2), against analytic DOF/nnz counts, and
against finite differences (gradient kernel). The Kratos path is covered by
``tests/kratos_spikes/test_user_api_kratos.py``.
"""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pytest
from pydantic import ValidationError

import wht_iga.compat.coreform as cf
import wht_iga.spec as s
from wht_iga.backends import GridMeshEngine, KratosPythonBackend, ReferenceNumbaBackend
from wht_iga.plan import PlanError, compile_plan
from wht_iga.ports import AnalysisCancelled, CancelToken, MeshEnginePort, SolverBackendPort
from wht_iga.runner import AnalysisSession, SessionError

E_MOD, NU = 2.0e5, 0.3  # mm-N-tonne-s steel
L = 10.0
T_X = 2.0  # traction (MPa)
ROLLERS = (s.DisplacementBCSpec(set="xmin", ux=0.0), s.DisplacementBCSpec(set="ymin", uy=0.0),
           s.DisplacementBCSpec(set="zmin", uz=0.0))


def bar_study(*, degree=2, n_cells=(8, 2, 2), loads=None, bcs=ROLLERS, probes=(), nu=NU, density=0.0,
              mesh=None, steps=None, units="mm-N-tonne-s", E=E_MOD, name="bar") -> s.StudySpec:
    """Box bar [0, L] x [0, 1] x [0, 1] with roller supports on the three lo faces."""
    if loads is None:
        loads = (s.TractionLoadSpec(set="xmax", vector=(T_X, 0.0, 0.0)),)
    return s.StudySpec(
        name=name, units=units, parts=[s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel")],
        materials={"steel": s.LinearElasticSpec(E=E, nu=nu, density=density)},
        mesh=mesh or s.BodyFitSpec(degree=degree, n_cells=n_cells),
        steps=steps or [s.StaticStepSpec(name="pull", bcs=bcs, loads=loads)], probes=probes)


LINE = s.LineProbeSpec(name="axis", start=(0.0, 0.3, 0.7), end=(L, 0.3, 0.7), n_points=5,
                       variables=("U", "S", "E", "MISES"))


def run_reference(study, tmp_path, **session_kw):
    """Compile + run on the oracle backend."""
    plan = compile_plan(study, output_dir=tmp_path / study.name)
    with AnalysisSession(plan, backend=ReferenceNumbaBackend(), **session_kw) as ses:
        return plan, ses.run()


# ============================================================================ Layer 1: spec
def test_spec_frozen_and_extra_forbidden():
    st = bar_study()
    with pytest.raises(ValidationError):
        st.name = "other"  # frozen
    with pytest.raises(ValidationError):
        s.ImmersedGridSpec(degree=2, element_size=1.0, elemnt_size=2.0)  # typo -> extra forbidden
    assert isinstance(st.parts, tuple) and isinstance(st.steps[0].bcs, tuple)


def test_with_updates_copies_and_revalidates():
    st = bar_study()
    fine = st.with_updates(mesh=st.mesh.with_updates(n_cells=(16, 4, 4)))
    assert fine.mesh.n_cells == (16, 4, 4) and st.mesh.n_cells == (8, 2, 2)
    with pytest.raises(ValidationError):
        st.mesh.with_updates(continuity=2)  # continuity must stay < degree
    with pytest.raises(ValueError, match="unknown field"):
        st.with_updates(nonsense=1)
    with pytest.raises(ValidationError):  # cross references re-checked on update
        st.with_updates(materials={"alu": s.LinearElasticSpec(E=7e4, nu=0.33)})


def test_json_roundtrip_and_schema():
    st = bar_study(probes=(LINE,), mesh=s.ImmersedGridSpec(degree=2, element_size=1.0))
    back = s.StudySpec.from_json(st.to_json())
    assert back == st and isinstance(back.mesh, s.ImmersedGridSpec)
    schema = s.StudySpec.json_schema()
    assert schema == s.StudySpec.model_json_schema()
    mesh_ref = schema["properties"]["mesh"]
    assert set(mesh_ref["discriminator"]["mapping"]) == {"immersed", "body_fit"}
    assert {"name", "units", "parts", "materials", "mesh", "steps", "probes"} <= set(schema["properties"])


def test_cross_reference_validation():
    with pytest.raises(ValidationError, match="unknown material"):
        s.StudySpec(name="x", parts=[s.PartSpec.from_box((0, 0, 0), (1, 1, 1), "nope")],
                    materials={"steel": s.LinearElasticSpec(E=2e5, nu=0.3)},
                    mesh=s.BodyFitSpec(n_cells=(1, 1, 1)),
                    steps=[s.StaticStepSpec(name="a", bcs=ROLLERS)])
    with pytest.raises(ValidationError, match="unknown set"):
        bar_study(bcs=(s.DisplacementBCSpec(set="hole", ux=0.0),))
    with pytest.raises(ValidationError, match="duplicate step"):
        bar_study(steps=[s.StaticStepSpec(name="a", bcs=ROLLERS), s.StaticStepSpec(name="a", bcs=ROLLERS)])
    with pytest.raises(ValidationError, match="reserved"):
        bar_study().with_updates(sets={"xmin": s.GridFaceSetSpec(axis=0, side=0)})


def test_unit_system_plausibility():
    with pytest.raises(ValidationError, match="check the unit system"):
        bar_study(E=2.0e11)  # SI steel modulus typed into mm-N-tonne-s
    with pytest.raises(ValidationError, match="check the unit system"):
        bar_study(density=7850.0)  # SI density typed into mm-N-tonne-s
    bar_study(density=7.85e-9)  # mm-N-tonne-s steel
    bar_study(units="m-N-kg-s", E=2.0e11, density=7850.0)  # SI steel
    with pytest.raises(ValidationError, match="check the unit system"):
        bar_study(units="m-N-kg-s", E=2.0e11, density=7.85e-9)


def test_physical_range_validation():
    with pytest.raises(ValidationError):
        s.LinearElasticSpec(E=2e5, nu=0.5)
    with pytest.raises(ValidationError):
        s.LinearElasticSpec(E=0.0, nu=0.3)
    with pytest.raises(ValidationError):
        s.ImmersedGridSpec(element_size=1.0, small_cell_volume_ratio=0.5)
    with pytest.raises(ValidationError):
        s.ImmersedGridSpec(element_size=1.0, padding=(-1, 0, 0))
    with pytest.raises(ValidationError, match="exactly one"):
        s.ImmersedGridSpec(element_size=1.0, n_cells=(2, 2, 2))
    with pytest.raises(ValidationError, match="constrains no component"):
        s.DisplacementBCSpec(set="xmin")
    with pytest.raises(ValidationError, match="must differ"):
        s.LineProbeSpec(name="l", start=(0, 0, 0), end=(0, 0, 0))
    with pytest.raises(ValidationError, match="unknown probe variable"):
        s.PointProbeSpec(name="p", location=(0, 0, 0), variables=("SIGMA",))
    with pytest.raises(ValidationError):
        s.ImmersedGridSpec(element_size=float("nan"))
    with pytest.raises(ValidationError, match="exactly one of sets"):
        s.RefinementSpec(name="r", level=1)


def test_probe_variables_match_resultview_aliases():
    from wht_iga.query.results import FIELD_ALIASES

    assert set(s.PROBE_VARIABLE_CANONICAL) == set(FIELD_ALIASES)


# ============================================================================ Layer 2: plan
def test_plan_grid_and_dof_counts(tmp_path):
    rep = compile_plan(bar_study(), output_dir=tmp_path).report
    assert rep.is_valid and not rep.grid.trimmed and rep.grid.n_cells == (8, 2, 2)
    assert rep.total_dofs == 3 * (8 + 2) * (2 + 2) * (2 + 2)  # C^{p-1}: n + p control points per axis
    c0 = compile_plan(bar_study(mesh=s.BodyFitSpec(degree=2, continuity=0, n_cells=(8, 2, 2))),
                      output_dir=tmp_path).report
    assert c0.total_dofs == 3 * (2 * 8 + 1) * (2 * 2 + 1) * (2 * 2 + 1)  # C^0: p*n + 1 per axis
    bf = compile_plan(bar_study(mesh=s.BodyFitSpec(degree=1, element_size=3.0)), output_dir=tmp_path).report
    assert bf.grid.n_cells == (3, 1, 1) and bf.grid.hi == (L, 1.0, 1.0)  # round(10/3) = 3
    ex = compile_plan(bar_study(mesh=s.ImmersedGridSpec(degree=2, element_size=0.5)), output_dir=tmp_path).report
    assert not ex.grid.trimmed and "IMMERSED_EXACT_FIT" in ex.codes() and ex.grid.n_cells == (20, 2, 2)


def test_plan_padding_off_domain_and_coverage(tmp_path):
    pad = compile_plan(bar_study(mesh=s.ImmersedGridSpec(degree=1, element_size=1.0, padding=(1, 0, 0))),
                       output_dir=tmp_path).report
    g = pad.grid
    assert g.trimmed and g.n_cells == (12, 1, 1) and g.lo[0] == -1.0 and g.hi[0] == L + 1.0
    off = [d for d in pad.diagnostics if d.code == "SET_OFF_DOMAIN"]
    assert len(off) == 2 and all("lies in padding" in d.message for d in off)  # xmin BC, xmax traction
    # zmax lies on the part but the part covers only y in [0, 1] of the y in [0, 2] grid face
    cov = compile_plan(bar_study(
        mesh=s.ImmersedGridSpec(degree=1, n_cells=(10, 2, 1), bounds=((0, 0, 0), (L, 2, 1))),
        loads=(s.TractionLoadSpec(set="zmax", vector=(0, 0, 1.0)),)), output_dir=tmp_path).report
    assert "LOAD_FACE_NOT_COVERED" in cov.codes() and not cov.is_valid
    with pytest.raises(PlanError, match="LOAD_FACE_NOT_COVERED"):
        cov.raise_if_invalid()


def test_plan_unsupported_combinations(tmp_path):
    ref = bar_study(mesh=s.ImmersedGridSpec(degree=2, element_size=1.0, refinements=[
        s.RefinementSpec(name="r", level=1, region=s.SphereRegionSpec(center=(5, 0.5, 0.5), radius=1.0))]))
    assert "UNSUPPORTED_REFINEMENT" in compile_plan(ref, output_dir=tmp_path).report.codes()
    nt = bar_study(probes=(s.PointProbeSpec(name="t", location=(1, 0.5, 0.5), variables=("NT11",)),))
    assert "PROBE_VARIABLE_UNSUPPORTED" in compile_plan(nt, output_dir=tmp_path).report.codes()
    fol = bar_study(steps=[s.StaticStepSpec(name="p", bcs=ROLLERS, linear_solver="pardiso_ldlt", loads=[
        s.PressureLoadSpec(set="xmax", magnitude=1.0, follower=True)])])
    assert "FOLLOWER_NEEDS_LU" in compile_plan(fol, output_dir=tmp_path).report.codes()
    stp = bar_study().with_updates(parts=(s.PartSpec.from_step(tmp_path / "missing.step", "steel"),))
    rep = compile_plan(stp, output_dir=tmp_path).report
    assert rep.codes() == ("BODYFIT_REQUIRES_BOX",) and not rep.is_valid
    imm = stp.with_updates(mesh=s.ImmersedGridSpec(element_size=1.0))
    assert "PART_FILE_MISSING" in compile_plan(imm, output_dir=tmp_path).report.codes()
    two = bar_study().with_updates(parts=(s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel", name="a"),
                                          s.PartSpec.from_box((0, 0, 1), (L, 1, 2), "steel", name="b")))
    assert "UNSUPPORTED_MULTI_PART" in compile_plan(two, output_dir=tmp_path).report.codes()


def test_plan_memory_and_solver_prediction(tmp_path):
    from wht_iga.solver.memory import est_amgcl_bytes, est_pardiso_bytes

    rep = compile_plan(bar_study(), output_dir=tmp_path, avail_bytes=10 ** 12).report
    n_dof = 3 * 10 * 4 * 4
    nnz = n_dof * 3 * (5 * 4 * 4)  # overlap: min(2p+1, n+p) per axis = 5, 4, 4
    assert (rep.est_active_dofs, rep.est_nnz) == (n_dof, nnz)
    st = rep.steps[0]
    assert st.kratos_solver == "pardiso_ldlt" and not st.switched
    assert st.est_pardiso_bytes == est_pardiso_bytes(n_dof, nnz, 2) == st.est_peak_bytes
    tiny = compile_plan(bar_study(), output_dir=tmp_path, avail_bytes=1000).report
    assert tiny.steps[0].switched and tiny.steps[0].kratos_solver == "amgcl_cg"
    assert tiny.steps[0].est_peak_bytes == est_amgcl_bytes(n_dof, nnz) and "AUTO_SOLVER_SWITCH" in tiny.codes()
    expl = bar_study(steps=[s.StaticStepSpec(name="p", bcs=ROLLERS, linear_solver="pardiso_ldlt")])
    assert "MEMORY_OVER_BUDGET" in compile_plan(expl, output_dir=tmp_path, avail_bytes=1000).report.codes()
    fol = bar_study(steps=[s.StaticStepSpec(name="p", bcs=ROLLERS, loads=[
        s.PressureLoadSpec(set="xmax", magnitude=1.0, follower=True)])])
    assert compile_plan(fol, output_dir=tmp_path, avail_bytes=10 ** 12).report.steps[0].kratos_solver == "pardiso_lu"


def test_plan_job_spec_mapping(tmp_path):
    from wht_iga.solver.validate import validate

    st = bar_study(density=7.85e-9, mesh=s.BodyFitSpec(degree=2, n_cells=(4, 1, 1), bc_penalty=1e9), steps=[
        s.StaticStepSpec(name="a", bcs=ROLLERS, loads=[s.GravityLoadSpec(acceleration=(0, 0, -9810.0)),
                                                       s.PressureLoadSpec(set="xmax", magnitude=-1.5)]),
        s.StaticStepSpec(name="b", bcs=ROLLERS, linear_solver="pardiso_lu")])
    plan = compile_plan(st, output_dir=tmp_path)
    js = plan.job_spec("a")
    assert js.name == "bar_a" and plan.artifacts.step("a").mesh_db.name == "bar_a.wiga.sqlite"
    assert [b.face for b in js.bcs] == [(0, 0), (1, 0), (2, 0)]
    assert js.bcs[1].components == (False, True, False)
    assert js.pressures[0].face == (0, 1) and js.pressures[0].magnitude == -1.5
    assert js.gravity == (0.0, 0.0, -9810.0) and js.beta == 1e9 and js.density == 7.85e-9
    assert plan.job_spec("b").linear_solver == "pardiso_lu"
    assert not [i for i in validate(js) if i.severity == "error"]
    single = compile_plan(bar_study(), output_dir=tmp_path)
    assert single.job_spec("pull").name == "bar" and single.artifacts.step("pull").vtkhdf.name == "bar.vtkhdf"


# ============================================================================ Layers 3-4: reference backend + session
@pytest.mark.parametrize("degree", [1, 2])
def test_reference_uniaxial_bar_matches_closed_form(tmp_path, degree):
    _, res = run_reference(bar_study(degree=degree, probes=(LINE,)), tmp_path)
    pts = LINE.points()
    pr = res.probe("axis")
    u_exact = np.column_stack([T_X * pts[:, 0] / E_MOD, -NU * T_X * pts[:, 1] / E_MOD,
                               -NU * T_X * pts[:, 2] / E_MOD])
    np.testing.assert_allclose(pr["U"], u_exact, rtol=0, atol=1e-12 * T_X * L / E_MOD)
    np.testing.assert_allclose(pr["S"], np.tile([T_X, 0, 0, 0, 0, 0], (5, 1)), rtol=0, atol=1e-10 * T_X)
    np.testing.assert_allclose(pr["E"], np.tile([T_X / E_MOD, -NU * T_X / E_MOD, -NU * T_X / E_MOD, 0, 0, 0],
                                                (5, 1)), rtol=0, atol=1e-10 * T_X / E_MOD)
    np.testing.assert_allclose(pr["MISES"], np.full(5, T_X), rtol=1e-10)
    step = res.step()
    np.testing.assert_allclose(step.total_force, [T_X * 1.0, 0, 0], rtol=1e-12, atol=1e-12)
    assert step.reaction_residual < 1e-10 * T_X


def test_reference_gravity_column_quadratic_exact(tmp_path):
    rho, g = 7850.0, 9.81  # SI, nu = 0: u_x = rho g / E (x^2 / 2 - L x) is exactly quadratic
    E_si = 2.0e11
    st = bar_study(units="m-N-kg-s", E=E_si, nu=0.0, density=rho, degree=2, n_cells=(3, 1, 1),
                   loads=(s.GravityLoadSpec(acceleration=(-g, 0.0, 0.0)),), probes=(LINE,))
    _, res = run_reference(st, tmp_path)
    x = LINE.points()[:, 0]
    pr = res.probe("axis")
    np.testing.assert_allclose(pr["U"][:, 0], rho * g / E_si * (x ** 2 / 2 - L * x), rtol=1e-9,
                               atol=1e-12 * rho * g * L ** 2 / E_si)
    np.testing.assert_allclose(pr["S"][:, 0], rho * g * (x - L), rtol=1e-9, atol=1e-9 * rho * g * L)
    assert res.step().reaction_residual < 1e-9 * rho * g * L


def test_reference_pressure_sign_and_prescribed_displacement(tmp_path):
    p = 3.0  # repo convention: traction = +p * n_out (S08) -> tension on xmax
    _, res = run_reference(bar_study(loads=(s.PressureLoadSpec(set="xmax", magnitude=p),),
                                     probes=(s.PointProbeSpec(name="tip", location=(L, 0, 0)),)), tmp_path / "p")
    np.testing.assert_allclose(res.probe("tip")["U"][0, 0], p * L / E_MOD, rtol=1e-10)
    d = 0.01
    bcs = ROLLERS + (s.DisplacementBCSpec(set="xmax", ux=d),)
    _, res2 = run_reference(bar_study(bcs=bcs, loads=(), probes=(LINE,)), tmp_path / "d")
    pr = res2.probe("axis")
    np.testing.assert_allclose(pr["S"][:, 0], E_MOD * d / L, rtol=1e-9)
    np.testing.assert_allclose(pr["U"][:, 1], -NU * d / L * LINE.points()[:, 1], rtol=1e-9, atol=1e-15)


def test_solver_log_on_off(tmp_path, monkeypatch):
    monkeypatch.delenv("WHT_SOLVER_LOG", raising=False)
    plan, res = run_reference(bar_study(name="logon"), tmp_path, enable_solver_log=True)
    log = plan.artifacts.step("pull").solver_log
    text = log.read_text(encoding="utf-8")
    for k in range(1, 7):
        assert f"[{k}/6]" in text
    assert "HWT-IGA SOLVER" in text and "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in text
    assert res.step().solver_log == text
    plan2, res2 = run_reference(bar_study(name="logoff"), tmp_path, enable_solver_log=False)
    assert not plan2.artifacts.step("pull").solver_log.exists() and res2.step().solver_log == ""


def test_progress_monotone_and_cancel_cleans_partial_artifacts(tmp_path):
    events = []
    plan, _ = run_reference(bar_study(name="prog"), tmp_path, on_progress=events.append)
    fr = [e.fraction for e in events]
    assert fr == sorted(fr) and fr[0] == 0.0 and fr[-1] == 1.0
    assert {"session", "mesh", "solve:pull", "post:pull", "done"} <= {e.stage for e in events}

    plan = compile_plan(bar_study(name="cancel"), output_dir=tmp_path / "c")
    out = plan.artifacts.output_dir
    out.mkdir(parents=True)
    stale = plan.artifacts.step("pull").probe_json
    stale.write_text("old", encoding="utf-8")  # pre-existing artifact, never touched -> kept
    (out / "notes.txt").write_text("keep", encoding="utf-8")  # unrelated file -> kept
    token = CancelToken()

    def cb(ev):
        if ev.stage.startswith("solve:"):
            token.cancel("test cancel")

    with pytest.raises(AnalysisCancelled, match="test cancel"):
        with AnalysisSession(plan, backend=ReferenceNumbaBackend(), on_progress=cb, cancel_token=token) as ses:
            ses.run()
    assert not plan.artifacts.spec_json.exists()
    assert not plan.artifacts.step("pull").solver_log.exists()  # written during the solve -> removed
    assert stale.read_text(encoding="utf-8") == "old" and (out / "notes.txt").exists()
    assert set(ses.removed_artifacts) == {plan.artifacts.spec_json, plan.artifacts.step("pull").solver_log}


def test_session_guards(tmp_path):
    plan = compile_plan(bar_study(), output_dir=tmp_path)
    with pytest.raises(SessionError, match="context manager"):
        AnalysisSession(plan, backend=ReferenceNumbaBackend()).run()
    trimmed = compile_plan(_immersed_gravity_study(), output_dir=tmp_path)
    assert trimmed.report.is_valid and trimmed.report.grid.trimmed
    with pytest.raises(SessionError, match="trimmed"):
        with AnalysisSession(trimmed, backend=ReferenceNumbaBackend()):
            pass
    bad = compile_plan(bar_study(bcs=(s.DisplacementBCSpec(set="xmin", ux=0.0),),
                                 mesh=s.ImmersedGridSpec(degree=1, element_size=1.0, padding=(1, 1, 1))),
                       output_dir=tmp_path)
    with pytest.raises(PlanError):
        with AnalysisSession(bad, backend=ReferenceNumbaBackend()):
            pass


def test_field_query_through_resultview(tmp_path):
    _, res = run_reference(bar_study(), tmp_path)
    u = res.field("U")
    np.testing.assert_array_equal(u, res.field("DISPLACEMENT"))
    np.testing.assert_allclose(u[:, 0].max(), T_X * L / E_MOD, rtol=1e-10)
    np.testing.assert_allclose(res.field("S")[:, 0], T_X, rtol=1e-9)


def test_gradient_kernel_matches_finite_differences():
    from wht_iga.reference.box import BoxSpline
    from wht_iga.runner.post import displacement_gradient

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.5), n_cells=(3, 2, 2), degree=2, continuity=1)
    rng = np.random.default_rng(7)
    U = rng.standard_normal(box.n_dofs)
    pts = rng.uniform([0.1, 0.1, 0.1], [1.9, 0.9, 1.4], size=(6, 3))
    G = displacement_gradient(box, U, pts)
    h = 1e-6
    for d in range(3):
        e = np.zeros(3)
        e[d] = h
        fd = (box.evaluate(U.reshape(-1, 3), pts + e) - box.evaluate(U.reshape(-1, 3), pts - e)) / (2 * h)
        np.testing.assert_allclose(G[:, :, d], fd, rtol=1e-6, atol=1e-6)


def test_ports_protocol_conformance():
    assert isinstance(ReferenceNumbaBackend(), SolverBackendPort)
    assert isinstance(KratosPythonBackend(threads=2), SolverBackendPort)
    assert isinstance(GridMeshEngine(), MeshEnginePort)
    with pytest.raises(TypeError, match="factories"):
        KratosPythonBackend(custom_processes=[object()])


def _immersed_gravity_study():
    """Box part inside a taller grid (OUT cells at y in [1, 1.5]); body load only."""
    return s.StudySpec(
        name="imm", units="m-N-kg-s", parts=[s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel")],
        materials={"steel": s.LinearElasticSpec(E=2.0e11, nu=0.0, density=7850.0)},
        mesh=s.ImmersedGridSpec(degree=2, n_cells=(10, 3, 2), bounds=((0, 0, 0), (L, 1.5, 1.0))),
        steps=[s.StaticStepSpec(name="g", bcs=ROLLERS, loads=[s.GravityLoadSpec(acceleration=(-9.81, 0, 0))])])


def test_mesh_engine_trimmed_box_volume(tmp_path):
    plan = compile_plan(_immersed_gravity_study(), output_dir=tmp_path)
    mesh = GridMeshEngine(workers=1, cache_dir="off").build(plan)
    assert mesh.trimmed and len(mesh.out_cells) == 10 * 1 * 2  # cells with y in [1, 1.5]
    cell_vol = 1.0 * 0.5 * 0.5
    vol = sum(float(r["eta"]) * cell_vol for r in mesh.cell_records)
    np.testing.assert_allclose(vol, L * 1.0 * 1.0, rtol=1e-9)


# ============================================================================ Layer 5 / package surface
def test_lazy_public_api_does_not_load_kratos_or_pydantic_eagerly():
    code = (
        "import sys, wht_iga\n"
        "assert 'pydantic' not in sys.modules, 'pydantic imported eagerly'\n"
        "import wht_iga as wi\n"
        "assert wi.StudySpec is __import__('wht_iga.spec', fromlist=['x']).StudySpec\n"
        "b = wi.KratosPythonBackend(); p = wi.compile_plan; v = wi.ResultView; c = wi.AnalysisSession\n"
        "import wht_iga.compat.coreform\n"
        "bad = [m for m in sys.modules if m.split('.')[0] in ('KratosMultiphysics', 'aieng', 'wht_topopt', 'gmsh')]\n"
        "print(bad); sys.exit(1 if bad else 0)\n")
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr


# ============================================================================ compat: Coreform facade
def test_coreform_mesh_facade_maps_to_core():
    m = cf.ImmersedRectilinearMesh(degree=2, continuity=1, bc_penalty_value=1e6, tessellation_size=0.1,
                                   sizing=cf.ElementSizeSizing(element_size=0.5), padding=(1, 1, 1),
                                   small_cell_volume_ratio=0.25)
    core = m.to_core()
    assert isinstance(core, s.ImmersedGridSpec)
    assert (core.degree, core.continuity, core.bc_penalty, core.element_size, core.padding,
            core.small_cell_volume_ratio, core.tessellation_size) == (2, 1, 1e6, 0.5, (1, 1, 1), 0.25, 0.1)
    assert cf.immersed_rectilinear(degree=1, sizing=cf.ElementIntervalCountSizing(
        element_interval_count=4)).to_core().n_cells == (4, 4, 4)
    by_axis = cf.body_fit(degree=2, continuity=1, sizing=cf.ElementIntervalCountByAxisSizing(
        element_interval_count_by_axis=(8, 2, 2))).to_core()
    assert isinstance(by_axis, s.BodyFitSpec) and by_axis.n_cells == (8, 2, 2)
    assert cf.ImmersedRectilinearMesh.from_dict(m.to_dict()) == m == cf.ImmersedRectilinearMesh.from_json(m.to_json())
    assert m.copy_with(degree=3, continuity=2).degree == 3 and m.degree == 2 and m.validate() == []
    assert m.to_dict()["meshType"] == "immersed_rectilinear" and m.to_dict()["schemaVersion"] == 1
    for bad in ({"continuity": 1, "degree": 1}, {"small_cell_volume_ratio": 0.5},
                {"frame": "local"}, {"element_layout": ("node_centered",) * 3}):
        with pytest.raises(ValidationError):
            cf.ImmersedRectilinearMesh(sizing=cf.ElementSizeSizing(element_size=1.0), **bad)
    with pytest.raises(ValueError, match="sizing is required"):
        cf.ImmersedRectilinearMesh(degree=1).to_core()
    r = cf.Refinement(name="edge", level=2, region=cf.SphereRegion(center=(0, 0, 0), radius=1.0)).to_core()
    assert isinstance(r.region, s.SphereRegionSpec) and r.level == 2
    assert isinstance(cf.BoxRegion(lo=(0, 0, 0), hi=(1, 1, 1)).to_core(), s.BoxRegionSpec)
    assert cf.CadEntitiesRegion(entities=("@face:f001",)).to_core().entities == ("@face:f001",)


def test_coreform_probe_facade_maps_to_core():
    lp = cf.LineProbe(name="l", start_location=(0, 0, 0), stop_location=(4, 0, 0), num_eval_points=5,
                      variables=("U",)).to_core()
    np.testing.assert_allclose(lp.points()[:, 0], [0, 1, 2, 3, 4])
    mp = cf.MultiPointProbe(name="m", locations=[(1, 2, 3), (4, 5, 6)], variables=("S",)).to_core()
    assert mp.points().shape == (2, 3) and mp.variables == ("S",)
    assert cf.PointProbe(name="p", location=(1, 1, 1)).to_core().variables == ()  # Coreform default
    with pytest.raises(ValidationError):
        cf.PointProbe(name="p", location=(1, 1, 1), region="instances")


def _job(tmp_path, **kw):
    return cf.Job("cfbar", parts=[s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel")],
                  materials={"steel": s.LinearElasticSpec(E=E_MOD, nu=NU)},
                  mesh=cf.BodyFitMesh(degree=2, continuity=1, sizing=cf.ElementIntervalCountByAxisSizing(
                      element_interval_count_by_axis=(8, 2, 2))),
                  steps=[s.StaticStepSpec(name="pull", bcs=ROLLERS,
                                          loads=[s.TractionLoadSpec(set="xmax", vector=(T_X, 0, 0))])],
                  probes=[cf.PointProbe(name="tip", location=(L, 1.0, 1.0), variables=("U", "S"))],
                  output_dir=tmp_path / "cf", backend=ReferenceNumbaBackend(), **kw)


def test_coreform_job_lifecycle(tmp_path):
    job = _job(tmp_path)
    assert job.status is cf.JobStatus.CREATED
    path = job.write_input()
    assert s.StudySpec.from_json(path.read_text(encoding="utf-8")) == job.spec
    with pytest.raises(RuntimeError, match="no results"):
        _ = job.results
    job.submit()
    assert job.wait_for_completion(timeout=300) is cf.JobStatus.COMPLETED, job.error
    pr = job.results.probe("tip")
    np.testing.assert_allclose(pr["U"][0], [T_X * L / E_MOD, -NU * T_X / E_MOD, -NU * T_X / E_MOD], rtol=1e-10)
    np.testing.assert_allclose(pr["S"][0, 0], T_X, rtol=1e-10)
    data = json.loads(job.plan.artifacts.step("pull").probe_json.read_text(encoding="utf-8"))
    np.testing.assert_allclose(data["tip"]["U"][0][0], T_X * L / E_MOD, rtol=1e-10)
    finer = job.copy_with(mesh=cf.BodyFitMesh(degree=1, sizing=cf.ElementSizeSizing(element_size=2.5)))
    assert finer.spec.mesh.element_size == 2.5 and job.spec.mesh.n_cells == (8, 2, 2)
    with pytest.raises(RuntimeError, match="CREATED"):
        job.submit()


def test_coreform_job_kill_aborts_and_cleans(tmp_path):
    holder = {}

    def cb(ev):
        if ev.stage.startswith("solve:"):
            holder["job"].kill("stop now")

    job = _job(tmp_path, on_progress=cb)
    holder["job"] = job
    job.submit(wait=True)
    assert job.status is cf.JobStatus.ABORTED and isinstance(job.error, AnalysisCancelled)
    assert not job.plan.artifacts.existing()  # everything this run wrote was removed
    with pytest.raises(RuntimeError, match="no results"):
        _ = job.results
    idle = _job(tmp_path / "idle")
    idle.kill()
    assert idle.status is cf.JobStatus.ABORTED
