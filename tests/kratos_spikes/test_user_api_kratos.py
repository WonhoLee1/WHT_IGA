"""User API on the Kratos backend (small smoke cases, seconds each).

* body-fit bar: KratosPythonBackend vs closed form and vs ReferenceNumbaBackend
  (penalty vs strong Dirichlet), escape hatch objects and hook/process calls;
* immersed box (OUT cells, out_mode="drop") under gravity vs the closed-form
  quadratic solution, plus the artifacts the pipeline and session write.
"""

from __future__ import annotations

import numpy as np
import pytest

import wht_iga.spec as s
from wht_iga.backends import GridMeshEngine, KratosPythonBackend, ReferenceNumbaBackend
from wht_iga.plan import compile_plan
from wht_iga.runner import AnalysisSession

pytestmark = [pytest.mark.kratos]

E_MOD, NU, L, T_X = 2.0e5, 0.3, 10.0, 2.0
ROLLERS = (s.DisplacementBCSpec(set="xmin", ux=0.0), s.DisplacementBCSpec(set="ymin", uy=0.0),
           s.DisplacementBCSpec(set="zmin", uz=0.0))
LINE = s.LineProbeSpec(name="axis", start=(0.0, 0.5, 0.5), end=(L, 0.5, 0.5), n_points=5,
                       variables=("U", "S"))


def _bar():
    return s.StudySpec(
        name="kbar", parts=[s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel")],
        materials={"steel": s.LinearElasticSpec(E=E_MOD, nu=NU)},
        mesh=s.BodyFitSpec(degree=2, n_cells=(8, 2, 2)),
        steps=[s.StaticStepSpec(name="pull", bcs=ROLLERS,
                                loads=[s.TractionLoadSpec(set="xmax", vector=(T_X, 0, 0))])],
        probes=[LINE])


def test_kratos_bar_matches_closed_form_reference_and_exposes_native(tmp_path):
    from wht_iga.kratos_backend import load_kratos

    KM = load_kratos()["KratosMultiphysics"]
    calls: list[str] = []

    class Recorder(KM.Process):
        """Records the lifecycle calls the backend makes."""

        def __init__(self, model, params):
            super().__init__()
            self.tag = params["tag"].GetString()

        def ExecuteInitialize(self):
            calls.append(f"{self.tag}:init")

        def ExecuteInitializeSolutionStep(self):
            calls.append(f"{self.tag}:step_init")

        def ExecuteFinalizeSolutionStep(self):
            calls.append(f"{self.tag}:step_final")

    seen = {}
    backend = KratosPythonBackend(
        threads=2, custom_processes=[lambda model, prm: Recorder(model, prm)],
        kratos_parameters={"tag": "rec"},
        on_before_solve=lambda mp: seen.setdefault("before", mp.Name),
        on_step_end=lambda mp, k: seen.setdefault("end", (mp.Name, k)))
    plan = compile_plan(_bar(), output_dir=tmp_path / "k")
    with AnalysisSession(plan, backend=backend) as ses:
        res = ses.run()
        assert ses.backend.root_model_part.Name == "LinearJob"
        assert ses.backend.kratos_model is not None and ses.backend.strategy is not None
    assert calls == ["rec:init", "rec:step_init", "rec:step_final"]
    assert seen == {"before": "LinearJob", "end": ("LinearJob", 1)}
    x = LINE.points()[:, 0]
    pk = res.probe("axis")
    np.testing.assert_allclose(pk["U"][:, 0], T_X * x / E_MOD, rtol=0, atol=1e-6 * T_X * L / E_MOD)
    np.testing.assert_allclose(pk["S"][:, 0], T_X, rtol=1e-5)
    plan_r = compile_plan(_bar(), output_dir=tmp_path / "r")
    with AnalysisSession(plan_r, backend=ReferenceNumbaBackend()) as ses_r:
        pr = ses_r.run().probe("axis")
    np.testing.assert_allclose(pk["U"], pr["U"], rtol=0, atol=1e-6 * T_X * L / E_MOD)
    step = res.step()
    assert step.converged and step.backend == "kratos_python"
    # surface-QP weight sums of the face triangulation: ~3e-11 relative (measured)
    np.testing.assert_allclose(step.total_force, [T_X, 0, 0], rtol=1e-9, atol=1e-9)
    arts = plan.artifacts.step("pull")
    for p in (arts.mesh_db, arts.results_db, arts.vtu, arts.vtkhdf, arts.probe_json, arts.solver_log):
        assert p.exists(), p
    assert "[6/6]" in arts.solver_log.read_text(encoding="utf-8")


def test_kratos_immersed_box_gravity_closed_form(tmp_path):
    rho, g, E_si = 7850.0, 9.81, 2.0e11
    study = s.StudySpec(
        name="kimm", units="m-N-kg-s", parts=[s.PartSpec.from_box((0, 0, 0), (L, 1, 1), "steel")],
        materials={"steel": s.LinearElasticSpec(E=E_si, nu=0.0, density=rho)},
        mesh=s.ImmersedGridSpec(degree=2, n_cells=(10, 3, 2), bounds=((0, 0, 0), (L, 1.5, 1.0))),
        steps=[s.StaticStepSpec(name="g", bcs=ROLLERS, loads=[s.GravityLoadSpec(acceleration=(-g, 0, 0))])],
        probes=[s.LineProbeSpec(name="axis", start=(0, 0.5, 0.5), end=(L, 0.5, 0.5), n_points=5,
                                variables=("U", "S"))])
    plan = compile_plan(study, output_dir=tmp_path)
    assert plan.report.grid.trimmed
    with AnalysisSession(plan, backend=KratosPythonBackend(),
                         mesh_engine=GridMeshEngine(workers=1, cache_dir="off")) as ses:
        mesh = ses.build_mesh()
        res = ses.run()
    assert len(mesh.out_cells) == 20
    x = np.linspace(0, L, 5)
    pr = res.probe("axis")
    scale = rho * g * L ** 2 / E_si
    np.testing.assert_allclose(pr["U"][:, 0], rho * g / E_si * (x ** 2 / 2 - L * x), rtol=0, atol=1e-6 * scale)
    np.testing.assert_allclose(pr["S"][:, 0], rho * g * (x - L), rtol=0, atol=1e-5 * rho * g * L)
