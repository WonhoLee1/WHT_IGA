# -*- coding: utf-8 -*-
"""Probe S / MISES / E of an F-bar (Neo-Hookean mean-dilatation) solve.

Before 2026-11-19 the session post-processed every solve with linear
``D(E, nu):eps`` -- for an F-bar solve at nu = 0.49999 that gave S_xx = 13.52
where the Cauchy stress is about 0.5 (found by the aieng bridge verifier).

Oracles (independent of ``runner.post``):
* homogeneous uniaxial Neo-Hookean: F = diag(a, b, b) from P11 = t, P22 = 0
  solved by bisection here; Cauchy sigma_xx = t / b**2, others 0;
  Green-Lagrange E = diag(a**2 - 1, b**2 - 1, b**2 - 1) / 2.
* non-uniform field: ``fbar.cauchy_stress`` (per-QP loop, separate code) at the
  QP coordinates of every cell -- checks the point -> cell -> Jbar mapping.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import wht_iga.spec as s
from wht_iga.analyses import fbar as FB
from wht_iga.backends import ReferenceNumbaBackend
from wht_iga.plan import compile_plan
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.assemble_elastic import assemble_face_traction
from wht_iga.runner import AnalysisSession
from wht_iga.runner.post import evaluate_point_fields, fbar_point_stress_strain

E_RUB, NU_RUB, T_RUB = 240.565, 0.49999, 0.5
# corner, interior, and points on interior cell faces (x = 1 of a 2-cell x split, y/z faces of 2x2)
PTS = [[2.0, 1.0, 1.0], [1.0, 0.5, 0.25], [0.3, 0.9, 0.6], [1.0, 0.5, 0.5], [0.0, 0.0, 0.0]]


def _nh_uniaxial_ab(E: float, nu: float, t: float) -> tuple[float, float]:
    """Stretches (a, b) of homogeneous uniaxial NH: P11 = t, P22 = P33 = 0 (nested bisection)."""
    mu, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))

    def a_of_b(b: float) -> float:  # P11 - P22 * (b/a)... use P22 = 0 => g = -mu b^2; P11 = mu a + g / a = t
        g = -mu * b * b
        # mu a^2 - t a + g = 0 -> positive root
        return (t + math.sqrt(t * t - 4.0 * mu * g)) / (2.0 * mu)

    def h(b: float) -> float:  # P22 = 0  <=>  lam ln(a b^2) - mu + mu b^2 = 0
        return lam * math.log(a_of_b(b) * b * b) - mu + mu * b * b

    lo, hi = 0.5, 1.0
    assert h(lo) < 0.0 < h(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    b = 0.5 * (lo + hi)
    return a_of_b(b), b


def _uniaxial_study(locking: str) -> s.StudySpec:
    return s.StudySpec(
        name=f"uni_{locking}",
        units="mm-N-tonne-s",
        parts=[s.PartSpec.from_box((0.0, 0.0, 0.0), (2.0, 1.0, 1.0), "rubber")],
        materials={"rubber": s.LinearElasticSpec(E=E_RUB, nu=NU_RUB)},
        mesh=s.BodyFitSpec(degree=2, continuity=1, n_cells=(2, 2, 2)),
        steps=[s.StaticStepSpec(
            name="pull",
            bcs=[s.DisplacementBCSpec(set="xmin", ux=0.0), s.DisplacementBCSpec(set="ymin", uy=0.0),
                 s.DisplacementBCSpec(set="zmin", uz=0.0)],
            loads=[s.TractionLoadSpec(set="xmax", vector=(T_RUB, 0.0, 0.0))],
            locking_mitigation=locking,
        )],
        probes=[s.MultiPointProbeSpec(name="pts", locations=PTS, variables=("U", "S", "MISES", "E"))],
    )


def test_session_fbar_probe_is_cauchy_stress_of_homogeneous_uniaxial(tmp_path):
    a, b = _nh_uniaxial_ab(E_RUB, NU_RUB, T_RUB)
    plan = compile_plan(_uniaxial_study("fbar"), output_dir=tmp_path / "fbar")
    with AnalysisSession(plan, backend=ReferenceNumbaBackend()) as ses:
        res = ses.run()
    r = res.step("pull")
    assert r.stress_measure == "fbar_neo_hookean_cauchy"
    assert r.n_iter > 1  # real Newton iteration count, not the old constant 1
    vals = r.probes["pts"]
    x = np.asarray(PTS)
    u_ref = np.column_stack([(a - 1) * x[:, 0], (b - 1) * x[:, 1], (b - 1) * x[:, 2]])
    assert np.max(np.abs(vals["U"] - u_ref)) < 1e-9 * np.max(np.abs(u_ref))
    # Stress gate: the Newton stops at |du| <= 1e-10 |u| (fbar.solve_static dtol); lam/t ~ 8e6 turns
    # a strain error ~1e-10 * |u|/L ~ 2e-13 into ~1.6e-6 t. Gate 1e-5 t (6x margin), far below the
    # 27x error of the old linear post-processing (S_xx 13.52 vs 0.5).
    S = np.asarray(vals["S"])
    s_ref = np.zeros_like(S)
    s_ref[:, 0] = T_RUB / b ** 2
    assert np.max(np.abs(S - s_ref)) < 1e-5 * T_RUB
    assert np.max(np.abs(np.asarray(vals["MISES"]) - T_RUB / b ** 2)) < 1e-5 * T_RUB
    E_ref = np.zeros_like(S)
    E_ref[:, 0], E_ref[:, 1], E_ref[:, 2] = 0.5 * (a * a - 1), 0.5 * (b * b - 1), 0.5 * (b * b - 1)
    assert np.max(np.abs(np.asarray(vals["E"]) - E_ref)) < 1e-9 * abs(E_ref[0, 0])


def test_session_linear_probe_unchanged(tmp_path):
    """No F-bar: S = D(E, nu) eps of the small-strain solution, stress_measure linear."""
    nu = 0.3
    study = _uniaxial_study("none")
    study.materials["rubber"] = s.LinearElasticSpec(E=E_RUB, nu=nu)
    plan = compile_plan(study, output_dir=tmp_path / "lin")
    with AnalysisSession(plan, backend=ReferenceNumbaBackend()) as ses:
        r = ses.run().step("pull")
    assert r.stress_measure == "linear_elastic"
    S = np.asarray(r.probes["pts"]["S"])
    s_ref = np.zeros_like(S)
    s_ref[:, 0] = T_RUB
    assert np.max(np.abs(S - s_ref)) < 1e-9 * T_RUB


@pytest.mark.parametrize("cont", [1, 0])
def test_point_stress_matches_fbar_qp_cauchy_in_every_cell(cont):
    """Non-uniform cantilever: point evaluation at each cell's QPs == fbar.cauchy_stress per cell."""
    n_cells, p = (3, 2, 2), 2
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(4.0, 1.0, 1.0), n_cells=n_cells, degree=p, continuity=cont)
    model = FB.build_model(box, p, E=E_RUB, nu=0.49)
    f = assemble_face_traction(box, 0, 1, np.array([0.0, 0.5, 0.2]))  # wiring-test load level; large-rotation Newton robustness is not tested here
    fixed = np.array([3 * int(i) + d for i in box.boundary_control_points(0, 0) for d in range(3)], dtype=np.int64)
    u = FB.solve_static(model, f, fixed, linear_solver="auto")
    ref = FB.cauchy_stress(model, u)
    xg, _ = np.polynomial.legendre.leggauss(p + 1)
    br = [np.unique(k) for k in box.knots]
    worst, spread = 0.0, []
    c = 0
    for iz in range(n_cells[2]):
        for iy in range(n_cells[1]):
            for ix in range(n_cells[0]):
                lo = np.array([br[0][ix], br[1][iy], br[2][iz]])
                hi = np.array([br[0][ix + 1], br[1][iy + 1], br[2][iz + 1]])
                g = np.array([[xg[a], xg[b_], xg[c_]] for c_ in range(p + 1) for b_ in range(p + 1)
                              for a in range(p + 1)])  # q = a + ng*(b + ng*c), fbar order
                pts = lo + 0.5 * (g + 1.0) * (hi - lo)
                sig, _e = fbar_point_stress_strain(box, u, pts, model)
                worst = max(worst, float(np.max(np.abs(sig - ref[c]))))
                spread.append(float(np.mean(ref[c][:, 0])))
                c += 1
    scale = max(float(np.max(np.abs(r))) for r in ref)
    assert worst < 1e-10 * scale
    assert np.ptp(spread) > 1e-2 * scale  # cells really differ: a wrong cell -> Jbar map would show

    # Interior cell face x = 4/3: the pressure (cell Jbar) jumps; documented side = +x cell (find_span).
    xf = float(br[0][1])
    on, plus, minus = (fbar_point_stress_strain(box, u, np.array([[x, 0.4, 0.3]]), model)[0][0]
                       for x in (xf, xf + 1e-12, xf - 1e-12))
    jump = float(np.max(np.abs(plus - minus)))
    assert jump > 1e-3 * scale
    assert float(np.max(np.abs(on - plus))) < 1e-6 * jump


def test_evaluate_point_fields_without_model_is_linear():
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=2, continuity=1)
    U = np.zeros(box.n_dofs)
    U[0::3] = 1e-3 * box.control_points()[:, 0]
    out = evaluate_point_fields(box, U, [[0.5, 0.5, 0.5]], 100.0, 0.25, ("S",))
    lam, mu = 100.0 * 0.25 / (1.25 * 0.5), 100.0 / 2.5
    assert out["S"][0, 0] == pytest.approx((lam + 2 * mu) * 1e-3, rel=1e-12)
