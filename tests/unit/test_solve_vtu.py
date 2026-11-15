"""Unit tests for solver.vtu proxy writer (Kratos-free, meshio re-read)."""

from __future__ import annotations

import numpy as np

from wht_iga.solver.vtu import hex_proxy_grid, read_vtu, write_vtu


def test_hex_proxy_grid() -> None:
    pts, conn = hex_proxy_grid(np.zeros(3), np.ones(3), (2, 1, 1))
    assert pts.shape == (12, 3) and conn.shape == (2, 8)
    np.testing.assert_allclose(pts.min(axis=0), [0, 0, 0])
    np.testing.assert_allclose(pts.max(axis=0), [1, 1, 1])


def test_vtu_write_reread(tmp_path) -> None:  # type: ignore[no-untyped-def]
    pts, conn = hex_proxy_grid(np.zeros(3), np.ones(3), (2, 1, 1))
    U = pts[:, 0] + 2 * pts[:, 1] + 3 * pts[:, 2]
    path = write_vtu(str(tmp_path / "m.vtu"), pts, conn, {"U": U})
    back = read_vtu(path)
    np.testing.assert_array_equal(back["points"], pts)
    assert back["cells"][0][0] == "hexahedron"
    np.testing.assert_array_equal(back["cells"][0][1], conn)
    np.testing.assert_array_equal(back["point_data"]["U"], U)
