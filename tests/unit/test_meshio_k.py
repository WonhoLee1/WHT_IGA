"""LS-DYNA .k reader: exact mesh/sets + *INCLUDE against inline decks."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.read_k import read_k

MAIN = """*KEYWORD
*TITLE
plate main
*PART
$ pid secid mid
1 1 1
BulkSteel
*NODE
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
5 0.0 0.0 1.0
6 1.0 0.0 1.0
7 1.0 1.0 1.0
8 0.0 1.0 1.0
*ELEMENT_SOLID
101 1 1 2 3 4 5 6 7 8
*INCLUDE
extra.k
*ELEMENT_SHELL
201 1 1 2 6 5
*SET_NODE_LIST
2 3
*SET_SOLID_LIST
101
*SET_SEGMENT
1 2 6 5
*MAT_ELASTIC
*END
"""

EXTRA = """*KEYWORD
*NODE
9 2.0 0.0 0.0
*ELEMENT_SOLID
102 1 2 9 3 6 0 0 0 0
*END
"""


def _write(tmp_path) -> str:
    (tmp_path / "extra.k").write_text(EXTRA, encoding="utf-8")
    main = tmp_path / "plate.k"
    main.write_text(MAIN, encoding="utf-8")
    return str(main)


def test_read_k_reference(tmp_path) -> None:
    model, issues = read_k(_write(tmp_path), "mm")
    assert model.nodes.shape == (9, 3)
    np.testing.assert_allclose(model.nodes[1], [1.0, 0.0, 0.0], rtol=0, atol=0)
    # hex8 from main + tet4 from *INCLUDE; shell excluded
    assert model.n_elem == 2
    kinds = sorted(b.elem_type for b in model.blocks)
    assert kinds == ["hex8", "tet4"]
    hex_block = next(b for b in model.blocks if b.elem_type == "hex8")
    np.testing.assert_array_equal(hex_block.conn[0], [0, 1, 2, 3, 4, 5, 6, 7])
    assert model.parts[hex_block.part_id] == "BulkSteel"
    tet_block = next(b for b in model.blocks if b.elem_type == "tet4")
    np.testing.assert_array_equal(tet_block.conn[0], [1, 8, 2, 5])
    assert model.node_sets["SET_0"] == [1, 2]
    assert model.elem_sets["SET_1"] == [0]
    # segment (1,2,6,5) zero-based = y=0 face of the hex -> (elem 0, local face 4)
    assert model.face_sets["SEGMENT"] == [(0, 4)]
    codes = [i["code"] for i in issues]
    assert "NON_SOLID_EXCLUDED" in codes
    assert "UNSUPPORTED_DECK_KEYWORD" in codes
    assert "INCLUDE_FOLLOWED" in codes


def test_read_k_roundtrip_lossless(tmp_path) -> None:
    model, _ = read_k(_write(tmp_path), "mm")
    back = model.roundtrip()
    np.testing.assert_array_equal(back.nodes, model.nodes)
    assert [b.conn.tolist() for b in back.blocks] == [b.conn.tolist() for b in model.blocks]
    assert back.face_sets == model.face_sets


def test_read_k_unit_system_required(tmp_path) -> None:
    with pytest.raises(ValueError):
        read_k(_write(tmp_path), "")
