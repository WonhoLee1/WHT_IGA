# -*- coding: utf-8 -*-
"""aieng-ui <-> WHT_IGA bridge (master plan v2.1 section 14, D8): end to end, Kratos-free.

The bridge client ``aieng-ui/backend/app/iga_bridge.py`` is stdlib-only and is
loaded here by file path (the aieng-ui ``app`` package needs the aieng311 env,
and ``wht_iga`` itself never imports aieng). Every solve below crosses the real
process boundary: the bridge writes a JSON request, spawns this whtiga
interpreter on ``app/iga_runner.py`` and parses the JSON result DTO.

Independent references and tolerances:

* uniaxial tension with roller BCs, linear elasticity: ``u = (t x/E, -nu t y/E,
  -nu t z/E)``, ``sigma = (t, 0, 0, 0, 0, 0)``. The field is linear, hence in
  the B-spline space, so the discrete solution is exact up to direct-solver
  round-off (cond ~ lam/mu <= 1e5 at nu = 0.49999): tolerance 1e-9 relative.
* the same problem with F-bar: the compressible Neo-Hookean law of
  ``wht_iga.analyses.fbar`` (W = mu/2 (I1 - 3) - mu ln J + lam/2 (ln J)^2, dead
  nominal traction) has the homogeneous solution F = diag(a, b, b) with
  P11 = t, P22 = 0 (solved by a 2x2 Newton below). Tolerance 1e-6 relative,
  the residual gate rtol = 1e-6 of ``fbar.solve_static``.
* in-process wht_iga calls on the same inputs in this interpreter: tolerance
  1e-10 relative (same code; only MKL/Pardiso threading round-off differs, and
  JSON floats round-trip exactly).
"""

from __future__ import annotations

import ast
import dataclasses
import importlib.util
import json
import math
import os
import sys
import textwrap
import time
from pathlib import Path

import numpy as np
import psutil
import pytest

_REPO = Path(__file__).resolve().parents[3]
_APP = _REPO / "aieng-ui" / "backend" / "app"


def _load_bridge():
    path = _APP / "iga_bridge.py"
    spec = importlib.util.spec_from_file_location("aieng_iga_bridge_under_test", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses resolve the module through sys.modules
    spec.loader.exec_module(module)
    return module


B = _load_bridge()
CFG = B.resolve_runner_config(python=sys.executable)  # this interpreter IS the whtiga env

E_LIN, NU_LIN, T_LIN = 210000.0, 0.3, 100.0
E_RUB, NU_RUB, T_RUB = 240.565, 0.49999, 0.5
PROBE_PTS = [[2.0, 1.0, 1.0], [1.0, 0.5, 0.25], [0.3, 0.9, 0.6]]


def _uniaxial_study(name: str, E: float, nu: float, t: float, locking: str = "none") -> dict:
    """Box [0,2]x[0,1]x[0,1], rollers on xmin/ymin/zmin, traction t e_x on xmax (plain JSON)."""
    return {
        "name": name,
        "parts": [{"name": "bar", "material": "m", "geometry": {"kind": "box", "lo": [0, 0, 0], "hi": [2, 1, 1]}}],
        "materials": {"m": {"kind": "linear_elastic", "E": E, "nu": nu}},
        "mesh": {"kind": "body_fit", "degree": 2, "continuity": 1, "n_cells": [2, 1, 1]},
        "steps": [{"name": "pull", "locking_mitigation": locking,
                   "bcs": [{"set": "xmin", "ux": 0.0}, {"set": "ymin", "uy": 0.0}, {"set": "zmin", "uz": 0.0}],
                   "loads": [{"kind": "traction", "set": "xmax", "vector": [t, 0.0, 0.0]}]}],
        "probes": [{"kind": "multi_point", "name": "pts", "locations": PROBE_PTS, "variables": ["U", "S", "MISES"]}],
    }


def _linear_uniaxial_u(E: float, nu: float, t: float) -> np.ndarray:
    x = np.asarray(PROBE_PTS)
    return np.column_stack([t * x[:, 0] / E, -nu * t * x[:, 1] / E, -nu * t * x[:, 2] / E])


def _neo_hookean_uniaxial_u(E: float, nu: float, t: float) -> np.ndarray:
    """Homogeneous F = diag(a, b, b): P11 = mu a + (lam ln J - mu)/a = t, P22 = mu b + (lam ln J - mu)/b = 0."""
    mu, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    a, b = 1.0 + t / E, 1.0 - nu * t / E
    for _ in range(50):
        J = a * b * b
        g = lam * math.log(J) - mu
        r = np.array([mu * a + g / a - t, mu * b + g / b])
        # dJ/da = b^2, dJ/db = 2ab -> d(ln J)/da = 1/a, d(ln J)/db = 2/b
        jac = np.array([[mu - g / a**2 + lam / a**2, 2 * lam / (a * b)],
                        [lam / (a * b), mu - g / b**2 + 2 * lam / b**2]])
        step = np.linalg.solve(jac, -r)
        a, b = a + step[0], b + step[1]
        if np.max(np.abs(step)) < 1e-16:
            break
    J = a * b * b
    g = lam * math.log(J) - mu
    # round-off floor of lam*ln(J) is ~ lam * 1e-16 ~ 4e-10 at nu = 0.49999, i.e. ~1e-9 t
    assert abs(mu * a + g / a - t) < 1e-8 * t and abs(mu * b + g / b) < 1e-8 * t
    x = np.asarray(PROBE_PTS)
    return np.column_stack([(a - 1) * x[:, 0], (b - 1) * x[:, 1], (b - 1) * x[:, 2]])


def _rel(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


# ---------------------------------------------------------------------------
# Structure of the boundary
# ---------------------------------------------------------------------------

def _imports(tree: ast.AST, *, top_level_only: bool) -> set[str]:
    nodes = tree.body if top_level_only else list(ast.walk(tree))
    roots: set[str] = set()
    for node in nodes:
        if isinstance(node, ast.Import):
            roots |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            roots.add(node.module.split(".")[0])
    return roots


def test_process_boundary_is_kept_in_the_source():
    """The bridge is stdlib-only, the runner imports wht_iga only inside ops, nobody imports aieng(-ui)."""
    stdlib = set(sys.stdlib_module_names) | {"__future__"}
    bridge = ast.parse((_APP / "iga_bridge.py").read_text(encoding="utf-8"))
    runner = ast.parse((_APP / "iga_runner.py").read_text(encoding="utf-8"))
    assert _imports(bridge, top_level_only=False) <= stdlib, _imports(bridge, top_level_only=False) - stdlib
    assert _imports(runner, top_level_only=True) <= stdlib, _imports(runner, top_level_only=True) - stdlib
    runner_all = _imports(runner, top_level_only=False)
    assert "wht_iga" in runner_all
    assert not runner_all & {"app", "aieng", "fastapi", "KratosMultiphysics"}
    for tree in (bridge, runner):  # relative imports would pull in the aieng-ui app package
        assert not [n for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.level > 0]

    reg = ast.parse((_APP / "runtime_registry" / "iga.py").read_text(encoding="utf-8"))
    gates = next(ast.literal_eval(n.value) for n in reg.body
                 if isinstance(n, ast.AnnAssign) and getattr(n.target, "id", "") == "_GATES")
    assert tuple(gates) == B.TOOL_NAMES == tuple(B.TOOL_INPUT_SCHEMAS) == tuple(B.TOOL_DESCRIPTIONS)
    assert gates["iga.run_study"] == (True, False) and gates["iga.run_method"] == (True, False)  # solver: approval
    assert gates["iga.runner_status"] == (False, True) and gates["iga.recommend_discretization"] == (False, True)


# ---------------------------------------------------------------------------
# ping / recommend
# ---------------------------------------------------------------------------

def test_runner_status_reports_the_whtiga_runtime_without_kratos():
    import wht_iga
    from wht_iga.analyses.methods import list_discretization_methods

    res = B.runner_status(config=CFG)
    prov = res["provenance"]
    assert res["status"] == "success" and res["ok"] is True and res["bridge"]["exit_code"] == 0
    assert prov["pid"] != os.getpid()  # really another process
    assert os.path.normcase(prov["python"]) == os.path.normcase(sys.executable)
    assert prov["wht_iga_version"] == wht_iga.__version__
    assert Path(prov["wht_iga_path"]) == Path(wht_iga.__file__).resolve().parent
    assert prov["kratos_loaded"] is False and res["capabilities"]["kratos"] == "blocked"
    assert [m["id"] for m in res["methods"]] == [m.id for m in list_discretization_methods()]


@pytest.mark.parametrize("physics, degree, nu, expected", [
    # rule 2 (cap at 3) + rule 5 (nu >= 0.49 -> F-bar), documented in recommend_discretization
    ("bending", 5, 0.49999, {"degree": 3, "capped": True, "locking_mitigation": "fbar"}),
    # rule 1 (bulk -> p=1), compressible -> no mitigation
    ("bulk", None, 0.3, {"degree": 1, "capped": False, "locking_mitigation": "none"}),
])
def test_recommendation_crosses_the_boundary_unchanged(physics, degree, nu, expected):
    from wht_iga.analyses.methods import recommend_discretization

    res = B.recommend_discretization(physics, degree, nu=nu, config=CFG)
    rec = res["recommendation"]
    assert {k: rec[k] for k in expected} == expected
    direct = dataclasses.asdict(recommend_discretization(physics, degree, nu=nu))
    assert rec == {k: list(v) if isinstance(v, tuple) else v for k, v in direct.items()}
    assert res["option_patch"] == {"degree": expected["degree"], "locking_mitigation": expected["locking_mitigation"]}
    assert res["credibility"] == "heuristic"


# ---------------------------------------------------------------------------
# StudySpec -> ReferenceNumbaBackend
# ---------------------------------------------------------------------------

def test_run_study_tool_linear_uniaxial_matches_the_analytic_solution(tmp_path):
    """MCP adapter path (handle_tool): JSON StudySpec in, analytic field out, artifacts under data_root."""
    res = B.handle_tool("iga.run_study", {"study": _uniaxial_study("uni_lin", E_LIN, NU_LIN, T_LIN),
                                          "run_id": "lin01"}, data_root=tmp_path, config=CFG)
    assert res["status"] == "success", res
    step = res["steps"][0]
    assert res["credibility"] == "computed" and step["converged"] and step["equilibrium_check"]["passed"]
    assert step["formulation"] == "linear_elastic" and step["solver"] == "reference.solve_with_dirichlet"
    assert step["locking_mitigation"] == {"in_spec": "none", "override": None, "requested": "none", "effective": "none"}
    assert Path(res["output_dir"]) == (tmp_path / "iga_runs" / "lin01").resolve()

    vals = step["probes"]["pts"]["values"]
    assert _rel(vals["U"], _linear_uniaxial_u(E_LIN, NU_LIN, T_LIN)) < 1e-9
    sig = np.asarray(vals["S"])
    assert np.max(np.abs(sig - np.array([T_LIN, 0, 0, 0, 0, 0]))) < 1e-9 * T_LIN
    assert np.max(np.abs(np.asarray(vals["MISES"]) - T_LIN)) < 1e-9 * T_LIN
    assert step["probes"]["pts"]["invalid_variables"] == {}
    assert _rel(step["applied_surface_force"], [T_LIN * 1.0 * 1.0, 0, 0]) < 1e-12  # t * area(xmax)
    assert step["n_dofs"] == 3 * (2 + 2) * (1 + 2) * (1 + 2)  # C^1 quadratic: n_cells + p per axis

    on_disk = json.loads(Path(step["artifacts"]["probe_json"]).read_text(encoding="utf-8"))
    assert np.array_equal(np.asarray(on_disk["pts"]["U"]), np.asarray(vals["U"]))  # DTO == artifact


@pytest.mark.parametrize("override, effective", [(None, "none"), ("auto", "fbar"), ("fbar", "fbar")])
def test_run_study_locking_override_reaches_the_backend(tmp_path, override, effective):
    """nu = 0.49999, spec says 'none': only the bridge override can switch the backend to F-bar."""
    lin_ref = _linear_uniaxial_u(E_RUB, NU_RUB, T_RUB)
    nh_ref = _neo_hookean_uniaxial_u(E_RUB, NU_RUB, T_RUB)
    assert _rel(lin_ref, nh_ref) > 1e-4  # the two references are distinguishable at the 1e-6 gate

    study = _uniaxial_study("uni_rub", E_RUB, NU_RUB, T_RUB, locking="none")
    res = B.run_study(study, output_dir=tmp_path / "bridge", locking_mitigation=override, config=CFG)
    step = res["steps"][0]
    assert step["locking_mitigation"]["effective"] == effective
    assert step["locking_mitigation"]["requested"] == (override or "none")
    assert res["credibility"] == "computed" and step["equilibrium_check"]["passed"]
    U = np.asarray(step["probes"]["pts"]["values"]["U"])
    if effective == "none":
        assert step["formulation"] == "linear_elastic"
        assert _rel(U, lin_ref) < 1e-9
        assert step["probes"]["pts"]["invalid_variables"] == {}
        return
    assert step["formulation"] == "fbar_mean_dilatation_neo_hookean"
    assert _rel(U, nh_ref) < 1e-6
    assert set(step["probes"]["pts"]["invalid_variables"]) == {"S", "MISES"}  # linear post-processing of F-bar
    assert any("not the F-bar Cauchy stress" in w for w in res["warnings"])

    import wht_iga.spec as s
    from wht_iga.backends import ReferenceNumbaBackend
    from wht_iga.plan import compile_plan
    from wht_iga.runner import AnalysisSession

    study["steps"][0]["locking_mitigation"] = override
    plan = compile_plan(s.StudySpec.model_validate(study), output_dir=tmp_path / "direct")
    with AnalysisSession(plan, backend=ReferenceNumbaBackend(), write_vtkhdf=False) as session:
        direct = session.run().step("pull")
    assert _rel(U, direct.probes["pts"]["U"]) < 1e-10
    assert _rel(step["max_abs_cp_displacement_component"], np.max(np.abs(direct.U))) < 1e-10


# ---------------------------------------------------------------------------
# BoxMethodProblem: recommend -> solve_with_discretization_method
# ---------------------------------------------------------------------------

def test_recommend_then_run_method_matches_in_process_fbar_solve():
    from wht_iga.analyses.methods import (
        BoxMethodProblem,
        parse_discretization_option,
        solve_with_discretization_method,
    )

    problem = {"lo": [0.0, 0.0, 0.0], "hi": [4.0, 1.0, 1.0], "E": E_RUB, "nu": NU_RUB,
               "fixed_face": [0, 0], "load_face": [0, 1], "traction": [0.0, 0.5, 0.0]}
    n_cells = [3, 2, 2]
    rec = B.recommend_discretization("bending", nu=NU_RUB, n_cells=n_cells, config=CFG)
    option = {"kind": "full_iga_body_fit", "n_cells": n_cells, **rec["option_patch"]}
    assert option["locking_mitigation"] == "fbar" and option["degree"] == 2
    probe = [[4.0, 0.5, 0.5], [2.0, 0.25, 0.75], [4.0, 1.0, 0.0]]
    res = B.run_method(problem, option, probe_points=probe, config=CFG)

    p = option["degree"]
    n_cp = [n + p for n in n_cells]
    assert res["total_dofs"] == 3 * n_cp[0] * n_cp[1] * n_cp[2]
    assert res["active_dofs"] == res["total_dofs"] - 3 * n_cp[1] * n_cp[2]  # clamped x = 0 face
    assert res["locking_mitigation"]["effective"] == "fbar" and res["method"] == "full_iga_body_fit"
    assert res["equilibrium_check"]["evaluated"] and res["equilibrium_check"]["passed"]
    assert res["credibility"] == "computed"

    direct = solve_with_discretization_method(
        BoxMethodProblem(lo=tuple(problem["lo"]), hi=tuple(problem["hi"]), E=E_RUB, nu=NU_RUB,
                         fixed_face=(0, 0), load_face=(0, 1), traction=(0.0, 0.5, 0.0)),
        parse_discretization_option(option), use_pardiso=False)
    assert _rel(res["tip_displacement"], direct.tip_displacement) < 1e-10
    assert _rel(res["probes"]["U"], direct.displacement_at(np.asarray(probe))) < 1e-10
    assert res["tip_point"] == [4.0, 0.5, 0.5]  # centre of the loaded x = 4 face


def test_kratos_route_is_refused_not_started():
    """nu = 0.3 without F-bar routes to the Kratos pipeline: explicit kratos_required, Kratos never loaded."""
    problem = {"lo": [0, 0, 0], "hi": [4, 1, 1], "E": 200.0, "nu": 0.3, "traction": [0, 0.5, 0]}
    with pytest.raises(B.IgaBridgeError) as err:
        B.run_method(problem, {"kind": "full_iga_body_fit", "degree": 2, "n_cells": [3, 2, 2]}, config=CFG)
    assert err.value.code == "kratos_required"
    assert err.value.details["provenance"]["kratos_loaded"] is False


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

def test_bad_inputs_are_refused_before_spawning(tmp_path):
    study = _uniaxial_study("x", E_LIN, NU_LIN, T_LIN)
    cases = [
        (lambda: B.run_study(study, output_dir=tmp_path, locking_mitigation="bogus", config=CFG), "invalid_request"),
        (lambda: B.run_study(study, output_dir=tmp_path, backend="kratos_python", config=CFG), "invalid_request"),
        (lambda: B.run_study("{not json", output_dir=tmp_path, config=CFG), "invalid_request"),
        (lambda: B.run_method({"lo": [0, 0], "hi": [1, 1, 1], "E": 1, "nu": 0.3}, {"kind": "x"}, config=CFG),
         "invalid_request"),
        (lambda: B.recommend_discretization("bulk", True, config=CFG), "invalid_request"),
    ]
    for call, code in cases:
        with pytest.raises(B.IgaBridgeError) as err:
            call()
        assert err.value.code == code
    assert B.handle_tool("iga.run_study", {"study": study, "extra": 1}, data_root=tmp_path,
                         config=CFG)["code"] == "invalid_request"
    assert B.handle_tool("iga.run_study", {"study": study, "run_id": "../escape"}, data_root=tmp_path,
                         config=CFG)["code"] == "invalid_request"
    taken = tmp_path / "iga_runs" / "taken"
    taken.mkdir(parents=True)
    (taken / "keep.txt").write_text("previous run", encoding="utf-8")
    res = B.handle_tool("iga.run_study", {"study": study, "run_id": "taken"}, data_root=tmp_path, config=CFG)
    assert res["status"] == "error" and res["code"] == "run_id_exists"
    assert (taken / "keep.txt").read_text(encoding="utf-8") == "previous run"
    assert B.handle_tool("iga.nope", {}, data_root=tmp_path, config=CFG)["code"] == "unknown_tool"
    missing = B.resolve_runner_config(env={B.ENV_PYTHON: str(tmp_path / "no_python.exe")})
    res = B.handle_tool("iga.runner_status", {}, data_root=tmp_path, config=missing)
    assert res["code"] == "runner_python_not_found"
    assert res["details"]["python_source"] == f"env:{B.ENV_PYTHON}"


def _study_with(**changes) -> dict:
    study = _uniaxial_study("bad", E_LIN, NU_LIN, T_LIN)
    study.update(changes)
    return study


@pytest.mark.parametrize("study, code, needle", [
    (_study_with(materials={"steel": {"kind": "linear_elastic", "E": 1.0, "nu": 0.3}}),
     "invalid_request", "unknown material"),
    (_study_with(parts=[{"name": n, "material": "m", "geometry": {"kind": "box", "lo": [0, 0, 0], "hi": [2, 1, 1]}}
                        for n in ("a", "b")]),
     "plan_invalid", "UNSUPPORTED_MULTI_PART"),
    # a VALID trimmed plan (grid wider than the part in y/z, gravity load): the plan passes, the
    # Kratos-free reference backend refuses trimmed grids at session entry
    (_study_with(mesh={"kind": "immersed", "degree": 2, "n_cells": [4, 4, 4],
                       "bounds": [[0.0, -0.5, -0.5], [2.0, 1.5, 1.5]]},
                 materials={"m": {"kind": "linear_elastic", "E": E_LIN, "nu": NU_LIN, "density": 7.85e-9}},
                 steps=[{"name": "sag", "bcs": [{"set": "xmin", "ux": 0.0, "uy": 0.0, "uz": 0.0}],
                         "loads": [{"kind": "gravity", "acceleration": [0.0, 0.0, -9810.0]}]}],
                 probes=[]),
     "backend_unsupported", "cannot run trimmed grids"),
])
def test_runner_reports_wht_iga_rejections_with_codes(tmp_path, study, code, needle):
    with pytest.raises(B.IgaBridgeError) as err:
        B.run_study(study, output_dir=tmp_path, config=CFG)
    assert err.value.code == code
    assert needle in err.value.message
    assert err.value.details["exit_code"] == 2  # handled failure, result file written


def _fake_runner(tmp_path: Path, body: str) -> object:
    """A stand-in runner process: exercises process-level failure handling only (no physics)."""
    script = tmp_path / "fake_runner.py"
    script.write_text(textwrap.dedent("""
        import json, sys, time, subprocess
        from pathlib import Path
        args = sys.argv[1:]
        result = Path(args[args.index("--result") + 1])
    """) + textwrap.dedent(body), encoding="utf-8")
    return B.resolve_runner_config(python=sys.executable, runner=script)


@pytest.mark.parametrize("body, code", [
    ("sys.stderr.write('boom: runner crashed'); sys.exit(3)\n", "runner_failed"),
    ("result.write_text('not json'); sys.exit(0)\n", "runner_protocol_error"),
    ("result.write_text(json.dumps({'schema': 'other', 'op': 'ping', 'ok': True})); sys.exit(0)\n",
     "runner_protocol_error"),
    ("result.write_text(json.dumps({'schema': 'aieng.iga_bridge.result/v1', 'op': 'ping', 'ok': True}));"
     " sys.exit(1)\n", "runner_protocol_error"),
])
def test_runner_process_failures_are_explicit(tmp_path, body, code):
    with pytest.raises(B.IgaBridgeError) as err:
        B.runner_status(config=_fake_runner(tmp_path, body))
    assert err.value.code == code
    if code == "runner_failed":
        assert err.value.details["exit_code"] == 3 and "boom: runner crashed" in err.value.details["stderr_tail"]


def test_timeout_kills_the_runner_process_tree(tmp_path):
    pid_file = tmp_path / "pids.json"
    cfg = _fake_runner(tmp_path, f"""
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                                 stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        Path({str(pid_file)!r}).write_text(json.dumps([__import__("os").getpid(), child.pid]))
        time.sleep(120)
    """)
    t0 = time.perf_counter()
    pids: list[int] = []
    try:
        with pytest.raises(B.IgaBridgeError) as err:
            B.runner_status(config=cfg, timeout_s=4.0)
        elapsed = time.perf_counter() - t0
        pids = json.loads(pid_file.read_text(encoding="utf-8"))
        assert err.value.code == "runner_timeout" and err.value.details["timeout_s"] == 4.0
        assert err.value.details["tree_killed"] is True
        assert elapsed < 60.0
        deadline = time.time() + 15.0
        while time.time() < deadline and any(psutil.pid_exists(p) for p in pids):
            time.sleep(0.2)
        alive = [p for p in pids if psutil.pid_exists(p)]
        assert alive == [], f"runner or its child survived the timeout: {alive}"
    finally:
        for p in pids or (json.loads(pid_file.read_text()) if pid_file.exists() else []):
            if psutil.pid_exists(p):
                psutil.Process(p).kill()
