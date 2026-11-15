"""run_linear output artifacts: mesh/results DB + VTU + probe JSON (Track SOLVE).

G1 evidence: DB roundtrip green through the real pipeline; VTU verified by
meshio re-read (ParaView not available in this environment -- recorded in
TRACK_SOLVE.md).
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec

pytestmark = pytest.mark.kratos


def test_run_linear_writes_artifacts(tmp_path) -> None:  # type: ignore[no-untyped-def]
    from wht_iga.solver import db as meshdb
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.vtu import read_vtu

    out = str(tmp_path)
    spec = JobSpec(name="art", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(1, 1, 1), degree=1, E=1.0, nu=0.0, output_dir=out,
                   bcs=[DirichletBC(face=(a, s), field=AffineField(
                       A=((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
                       b=(0.0, 0.0, 0.0))) for a in range(3) for s in (0, 1)],
                   probes=[{"name": "pk", "kind": "point", "xyz": (1.0, 1.0, 1.0), "variables": ["U"]}])
    res = run_linear(spec)
    assert res.mesh_db and res.results_db and res.vtu_path and res.probe_json
    for p in (res.mesh_db, res.results_db, res.vtu_path, res.probe_json):
        assert os.path.exists(p), p
    mesh = meshdb.read_mesh_db(res.mesh_db)
    assert mesh["meta"]["schema_version"] == meshdb.MESH_SCHEMA_VERSION
    assert mesh["grid"]["degree"] == 1
    back = meshdb.read_results_db(res.results_db)
    np.testing.assert_array_equal(back["cp_fields"][(0, "U")], res.U)
    np.testing.assert_array_equal(back["qp_fields"][(0, "S")], np.concatenate(res.cauchy_stress, axis=0))
    vtu = read_vtu(res.vtu_path)
    assert vtu["cells"][0][0] == "hexahedron"
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.vtu import hex_proxy_grid

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    pts, _ = hex_proxy_grid(np.zeros(3), np.ones(3), (1, 1, 1))
    np.testing.assert_array_equal(vtu["points"], pts)
    np.testing.assert_allclose(vtu["point_data"]["U"], box.evaluate(res.U.reshape(-1, 3), pts),
                               rtol=0.0, atol=1e-15)
    probe = json.load(open(res.probe_json, encoding="utf-8"))
    np.testing.assert_allclose(probe["pk"], [1.0, 2.0, 3.0], rtol=1e-4)
    print(f"\nP1 artifacts: mesh={os.path.basename(res.mesh_db or '')} "
          f"results + vtu + probe.json OK; probe pk={probe['pk']}")
