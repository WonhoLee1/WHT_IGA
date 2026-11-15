"""Quadrature DB tests: write → read roundtrip lossless (v2 §4.4)."""

from __future__ import annotations

import numpy as np

from wht_iga.geometry import quad_db, quadrature

V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
              [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
T = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
              [0, 1, 5], [0, 5, 4], [2, 3, 7], [2, 7, 6],
              [0, 4, 7], [0, 7, 3], [1, 2, 6], [1, 6, 5]])


def test_write_read_roundtrip(tmp_path) -> None:
    """Two cells (IN + tet-cut): meta, counts, ΣdV, U/S blobs roundtrip."""
    r_in = quadrature.cell_quadrature(V, T, np.array([0.2, 0.2, 0.2]),
                                      np.array([0.5, 0.5, 0.5]), q=1)
    s = 0.9
    TV = np.array([[0, 0, 0], [s, 0, 0], [0, s, 0], [0, 0, s]])
    TT = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    r_cut = quadrature.cell_quadrature(TV, TT, np.full(3, -0.1),
                                       np.full(3, 1.1), q=1)
    n_in, n_cut = len(r_in.w), len(r_cut.w)
    assert n_in > 0 and n_cut > 0
    u_in = np.tile([1.0, 2.0, 3.0], (n_in, 1))
    s_in = np.tile(np.arange(6.0), (n_in, 1))
    cells = [{"cell_id": 0, "material": 0, "x": r_in.x, "w": r_in.w,
              "U": u_in, "S": s_in},
             {"cell_id": 5, "material": 1, "x": r_cut.x, "w": r_cut.w}]
    db = str(tmp_path / "t_iga_quadrature.sqlite")
    n = quad_db.write_quad_db(db, cells, meta={"job": "t", "p": "2"})
    assert n == n_in + n_cut
    back = quad_db.read_quad_db(db)
    assert back["meta"]["schema_version"] == "1"
    assert back["meta"]["job"] == "t"
    assert len(back["qp"]) == n
    tot, cnt = quad_db.qp_totals(db)
    # SQL SUM adds rows sequentially, numpy .sum() pairwise: equal to roundoff,
    # not bit-for-bit (rtol=0 passed only by summation-order luck).
    np.testing.assert_allclose(tot, float(r_in.w.sum() + r_cut.w.sum()),
                               rtol=1e-14)
    assert cnt == n
    # blob payloads decode exactly
    got_u = np.frombuffer(back["qp_val"][0][2], dtype=np.float64)
    got_s = np.frombuffer(back["qp_val"][0][3], dtype=np.float64)
    np.testing.assert_allclose(got_u, [1.0, 2.0, 3.0], rtol=0)
    np.testing.assert_allclose(got_s, np.arange(6.0), rtol=0)
