"""Abaqus .inp reader: exact mesh/sets against an inline reference deck."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.read_inp import read_inp

DECK = """** plate deck: mesh + sets only
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9, 2., 0., 0.
*ELEMENT, TYPE=C3D8, ELSET=BULK
101, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D4, ELSET=FILL
102, 2, 9, 3, 6
*ELEMENT, TYPE=S4R, ELSET=SKIN
201, 1, 2, 6, 5
*NSET, NSET=LOADED
2, 3
*ELSET, ELSET=EALL, GENERATE
101, 102, 1
*SOLID SECTION, ELSET=BULK, MATERIAL=STEEL
*MATERIAL, NAME=STEEL
*ELASTIC
2.1e5, 0.3
*STEP
*STATIC
*END STEP
"""


def test_read_inp_reference(tmp_path) -> None:
    path = tmp_path / "plate.inp"
    path.write_text(DECK, encoding="utf-8")
    model, issues = read_inp(str(path), "mm")
    assert model.unit_system == "mm"
    np.testing.assert_allclose(
        model.nodes[1], [1.0, 0.0, 0.0], rtol=0, atol=1e-15
    )
    assert model.nodes.shape == (9, 3)
    assert model.n_elem == 2  # S4R shell excluded
    kinds = sorted(b.elem_type for b in model.blocks)
    assert kinds == ["hex8", "tet4"]
    hex_block = next(b for b in model.blocks if b.elem_type == "hex8")
    np.testing.assert_array_equal(hex_block.conn[0], [0, 1, 2, 3, 4, 5, 6, 7])
    tet_block = next(b for b in model.blocks if b.elem_type == "tet4")
    np.testing.assert_array_equal(tet_block.conn[0], [1, 8, 2, 5])
    # part mapping from SOLID SECTION material
    assert model.parts[hex_block.part_id] == "STEEL"
    # sets
    assert model.node_sets["LOADED"] == [1, 2]
    assert model.elem_sets["EALL"] == [0, 1]
    assert model.elem_sets["BULK"] == [0]
    assert model.elem_sets["FILL"] == [1]
    # scope issues: shell excluded + deck keywords ignored
    codes = [i["code"] for i in issues]
    assert "NON_SOLID_EXCLUDED" in codes
    assert "UNSUPPORTED_DECK_KEYWORD" in codes
    deck_msg = next(i["message"] for i in issues if i["code"] == "UNSUPPORTED_DECK_KEYWORD")
    assert "MATERIAL" in deck_msg and "STEP" in deck_msg


def test_read_inp_roundtrip_lossless(tmp_path) -> None:
    path = tmp_path / "plate.inp"
    path.write_text(DECK, encoding="utf-8")
    model, _ = read_inp(str(path), "mm")
    back = model.roundtrip()
    np.testing.assert_array_equal(back.nodes, model.nodes)
    assert [b.conn.tolist() for b in back.blocks] == [b.conn.tolist() for b in model.blocks]
    assert back.node_sets == model.node_sets and back.elem_sets == model.elem_sets


def test_read_inp_wrapped_hex20(tmp_path) -> None:
    """C3D20 connectivity wrapped across lines parses exactly."""
    nodes = "\n".join(f"{i}, {float(i)}, 0., 0." for i in range(1, 21))
    deck = ("*NODE\n" + nodes + "\n*ELEMENT, TYPE=C3D20, ELSET=BIG\n"
            "501, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11\n"
            "12, 13, 14, 15, 16, 17, 18, 19, 20\n")
    path = tmp_path / "big.inp"
    path.write_text(deck, encoding="utf-8")
    model, _ = read_inp(str(path), "m")
    assert model.n_elem == 1
    assert model.blocks[0].elem_type == "hex20"
    np.testing.assert_array_equal(model.blocks[0].conn[0], np.arange(20))
    assert model.elem_sets["BIG"] == [0]
    path = tmp_path / "plate.inp"
    path.write_text(DECK, encoding="utf-8")
    with pytest.raises(ValueError):
        read_inp(str(path), "")
