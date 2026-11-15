"""Gmsh .msh reader: physical groups, solid mapping, lower-d exclusion."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.read_msh import read_msh

MSH = """$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
1
3 1 "solid"
$EndPhysicalNames
$Nodes
8
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0 0 1
6 1 0 1
7 1 1 1
8 0 1 1
$EndNodes
$Elements
3
1 5 2 1 1 2 3 4 5 6 7 8
2 5 2 1 1 2 3 4 5 6 7 8
3 2 2 0 1 1 2 3
$EndElements
"""


def test_read_msh_reference(tmp_path) -> None:
    path = tmp_path / "cube.msh"
    path.write_text(MSH, encoding="utf-8")
    model, issues = read_msh(str(path), "m")
    assert model.nodes.shape == (8, 3)
    np.testing.assert_allclose(model.nodes[6], [1.0, 1.0, 1.0], rtol=0, atol=0)
    # 2 hexahedra kept; the triangle (physical 0) excluded
    assert model.n_elem == 2
    assert [b.elem_type for b in model.blocks] == ["hex8"]
    np.testing.assert_array_equal(model.blocks[0].conn[0], [0, 1, 2, 3, 4, 5, 6, 7])
    assert model.parts == {1: "solid"}
    codes = [i["code"] for i in issues]
    assert "NON_SOLID_EXCLUDED" in codes


def test_read_msh_roundtrip_lossless(tmp_path) -> None:
    path = tmp_path / "cube.msh"
    path.write_text(MSH, encoding="utf-8")
    model, _ = read_msh(str(path), "m")
    back = model.roundtrip()
    np.testing.assert_array_equal(back.nodes, model.nodes)
    assert [b.conn.tolist() for b in back.blocks] == [b.conn.tolist() for b in model.blocks]
    assert back.parts == model.parts


def test_read_msh_unit_system_required(tmp_path) -> None:
    path = tmp_path / "cube.msh"
    path.write_text(MSH, encoding="utf-8")
    with pytest.raises(ValueError):
        read_msh(str(path), "")
