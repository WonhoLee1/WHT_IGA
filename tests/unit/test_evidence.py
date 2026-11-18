# -*- coding: utf-8 -*-
"""Evidence export roundtrips (Kratos-free) + one Kratos end-to-end."""

from __future__ import annotations

import os

import numpy as np
import pytest


def _bar_mesh():
    xs = np.linspace(0, 4, 3)
    ys = np.linspace(0, 1, 2)
    zs = np.linspace(0, 1, 2)
    pts = np.array([[x, y, z] for z in zs for y in ys for x in xs], dtype=float)

    def nid(i, j, k):
        return i + 3 * (j + 2 * k)

    conn = np.array([
        [nid(0, 0, 0), nid(1, 0, 0), nid(1, 1, 0), nid(0, 1, 0),
         nid(0, 0, 1), nid(1, 0, 1), nid(1, 1, 1), nid(0, 1, 1)],
        [nid(1, 0, 0), nid(2, 0, 0), nid(2, 1, 0), nid(1, 1, 0),
         nid(1, 0, 1), nid(2, 0, 1), nid(2, 1, 1), nid(1, 1, 1)],
    ])
    u = np.column_stack([0.01 * pts[:, 0], np.zeros(12), np.zeros(12)])
    return pts, conn, {"DISPLACEMENT": u}


def test_evidence_all_formats_roundtrip(tmp_path):
    import meshio

    from wht_iga.post.evidence import exterior_triangles, write_evidence_set
    from wht_iga.post.vtk_writer import read_vtkhdf

    pts, conn, fields = _bar_mesh()
    assert exterior_triangles(conn).shape == (20, 3)  # 2 cells: 10 quads x2
    out = write_evidence_set(str(tmp_path), "bar", pts, conn, fields)
    assert set(out) == {"vtkhdf", "vtu", "msh", "fem", "stl", "png"}
    for p in out.values():
        assert os.path.isfile(p) and os.path.getsize(p) > 0, p
    for key in ("vtu", "msh", "fem", "stl"):
        m = meshio.read(out[key])
        assert m.points.shape == pts.shape, key
    assert meshio.read(out["stl"]).cells_dict["triangle"].shape == (20, 3)
    back = read_vtkhdf(out["vtkhdf"])
    assert back is not None


def test_evidence_rejects_bad_input(tmp_path):
    import numpy as np

    from wht_iga.post.evidence import write_evidence_set

    pts, conn, fields = _bar_mesh()
    with pytest.raises(ValueError):
        write_evidence_set(str(tmp_path), "empty", np.zeros((0, 3)), conn, fields)
    with pytest.raises(ValueError):
        write_evidence_set(str(tmp_path), "badidx", pts, conn + 100, fields)
    with pytest.raises(ValueError):
        write_evidence_set(str(tmp_path), "ragged", pts, conn,
                           {"DISPLACEMENT": np.zeros((5, 3))})


@pytest.mark.kratos
def test_evidence_native_fe_end_to_end(tmp_path):
    import numpy as np

    from wht_iga.analyses.methods import NativeFeOption, solve_with_discretization_method
    from wht_iga.analyses.presets import bulk_linear
    from wht_iga.post.evidence import write_evidence_set

    problem, _, _ = bulk_linear(1.0e5, 0.0, span=4.0, traction=10.0)
    sol = solve_with_discretization_method(
        problem, NativeFeOption(element_type="hex8", n_cells=(2, 1, 1)))
    assert sol.converged
    mesh, u = sol.native["mesh"], sol.native["u_nodal"]
    out = write_evidence_set(str(tmp_path), "bulk", np.asarray(mesh.nodes),
                             np.asarray(mesh.blocks[0].conn),
                             {"DISPLACEMENT": np.asarray(u).reshape(-1, 3)})
    assert all(os.path.getsize(p) > 0 for p in out.values())
