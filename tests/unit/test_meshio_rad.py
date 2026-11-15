"""Radioss .rad reader: exact mesh/sets against an inline reference deck."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.read_rad import read_rad

DECK = """/NODE
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
5 0.0 0.0 1.0
6 1.0 0.0 1.0
7 1.0 1.0 1.0
8 0.0 1.0 1.0
9 2.0 0.0 0.0
/BRICK
101 1 2 3 4 5 6 7 8 7
/TETRA4
102 2 9 3 6 7
/SHELL
201 1 2 6 5
/PART
7 BulkSteel
/MAT/PLAS_JOHNS/1
/GRNOD/PART/1
loaded nodes
2 3
/GRBRIC
1 101
/END
"""


def test_read_rad_reference(tmp_path) -> None:
    path = tmp_path / "plate_0000.rad"
    path.write_text(DECK, encoding="utf-8")
    model, issues = read_rad(str(path), "mm")
    assert model.nodes.shape == (9, 3)
    np.testing.assert_allclose(model.nodes[1], [1.0, 0.0, 0.0], rtol=0, atol=0)
    assert model.n_elem == 2  # /SHELL excluded
    kinds = sorted(b.elem_type for b in model.blocks)
    assert kinds == ["hex8", "tet4"]
    hex_block = next(b for b in model.blocks if b.elem_type == "hex8")
    np.testing.assert_array_equal(hex_block.conn[0], [0, 1, 2, 3, 4, 5, 6, 7])
    assert model.parts[hex_block.part_id] == "BulkSteel"
    assert model.node_sets["GRNOD_0"] == [1, 2]
    assert model.elem_sets["GRBRIC_1"] == [0]
    codes = [i["code"] for i in issues]
    assert "NON_SOLID_EXCLUDED" in codes
    assert "UNSUPPORTED_DECK_KEYWORD" in codes


def test_read_rad_roundtrip_lossless(tmp_path) -> None:
    path = tmp_path / "plate_0000.rad"
    path.write_text(DECK, encoding="utf-8")
    model, _ = read_rad(str(path), "mm")
    back = model.roundtrip()
    np.testing.assert_array_equal(back.nodes, model.nodes)
    assert [b.conn.tolist() for b in back.blocks] == [b.conn.tolist() for b in model.blocks]
    assert back.node_sets == model.node_sets and back.elem_sets == model.elem_sets


def test_read_rad_unit_system_required(tmp_path) -> None:
    path = tmp_path / "plate_0000.rad"
    path.write_text(DECK, encoding="utf-8")
    with pytest.raises(ValueError):
        read_rad(str(path), "")
