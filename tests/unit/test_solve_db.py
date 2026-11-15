"""Unit tests for solver.db roundtrip (Kratos-free, bit-exact BLOBs)."""

from __future__ import annotations

import numpy as np

from wht_iga.solver import db as meshdb


def test_mesh_db_roundtrip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rng = np.random.default_rng(11)
    xw = rng.random((4, 8, 4))
    xwn = rng.random((2, 6, 7))
    path = str(tmp_path / "job.wiga.sqlite")
    spec = {"name": "job", "degree": 2}
    meshdb.write_mesh_db(path, spec, meta_extra={"kratos": "test"},
                         degree=2, continuity=1, n_cells=(2, 1, 1),
                         cells=[{"i": 0, "j": 0, "k": 0, "status": 1, "eta": 1.0, "n_qp": 27}],
                         qp_volumes={0: xw}, qp_surfaces={(2, 1): xwn},
                         probes=[{"name": "pk", "kind": "point", "xyz": (1.0, 0.0, 0.0), "variables": ["U"]}],
                         issues=[{"severity": "warning", "code": "W", "message": "m", "entity": "e"}])
    back = meshdb.read_mesh_db(path)
    assert back["meta"]["schema_version"] == meshdb.MESH_SCHEMA_VERSION
    assert back["meta"]["kratos"] == "test"
    assert back["grid"]["degree"] == 2
    np.testing.assert_array_equal(back["qp_volumes"][0][1], xw)
    np.testing.assert_array_equal(back["qp_surfaces"][0][1], xwn)
    assert back["probes"][0]["name"] == "pk"
    assert back["issues"][0]["code"] == "W"
    assert json_spec_roundtrip(back["meta"]) == spec


def json_spec_roundtrip(meta: dict) -> dict:  # type: ignore[no-untyped-def]
    import json

    return json.loads(meta["spec_json"])


def test_results_db_roundtrip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rng = np.random.default_rng(12)
    U = rng.random(81)
    S = rng.random((8, 6))
    path = str(tmp_path / "job_iga_results.sqlite")
    meshdb.write_results_db(path, {"name": "job"}, cp_fields={(0, "U"): U}, qp_fields={(0, "S"): S})
    back = meshdb.read_results_db(path)
    assert back["meta"]["schema_version"] == meshdb.RESULTS_SCHEMA_VERSION
    np.testing.assert_array_equal(back["cp_fields"][(0, "U")], U)
    np.testing.assert_array_equal(back["qp_fields"][(0, "S")], S)
    # default frames mean "not recorded", never "converged" (F3)
    assert back["frames"][0]["converged"] == 0
    assert back["frames"][0]["n_iter"] == 0


def test_results_db_explicit_frames_roundtrip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = str(tmp_path / "job2_iga_results.sqlite")
    meshdb.write_results_db(
        path, {"name": "job2"},
        frames=[{"frame_id": 0, "step_id": 0, "time": 1.0,
                 "increment": 3, "converged": 1, "n_iter": 3}])
    back = meshdb.read_results_db(path)
    assert back["frames"][0]["converged"] == 1
    assert back["frames"][0]["n_iter"] == 3
    assert back["frames"][0]["increment"] == 3
