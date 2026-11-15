"""Nastran .bdf reader: pyNastran parse, CORD2R frame, PSOLID parts."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.read_nastran import read_nastran

DECK = """$$HMNAME COMP \\1\\"bulk"
CORD2R  10      0       0.      0.      0.      0.      1.      0.
        0.      0.      1.
GRID    1       0       0.      0.      0.
GRID    2       0       1.      0.      0.
GRID    3       0       1.      1.      0.
GRID    4       0       0.      1.      0.
GRID    5       0       0.      0.      1.
GRID    6       0       1.      0.      1.
GRID    7       0       1.      1.      1.
GRID    8       0       0.      1.      1.
GRID    9       10      1.      0.      0.
CHEXA   101     1       1       2       3       4       5       6
        7       8
CTETRA  102     1       1       2       4       5
CQUAD4  201     1       1       2       6       5
PSOLID  1       1
MAT1    1       2.1E5           .3      7.85-9
ENDDATA
"""


def test_read_nastran_reference(tmp_path) -> None:
    path = tmp_path / "plate.bdf"
    path.write_text(DECK, encoding="utf-8")
    model, issues = read_nastran(str(path), "mm")
    assert model.nodes.shape == (9, 3)
    np.testing.assert_allclose(model.nodes[8], [0.0, 0.0, 1.0], rtol=0, atol=1e-12)
    # CORD2R 10 maps local x -> global z (z-dir (0,1,0), xz-pt (0,0,1))
    assert model.n_elem == 2  # CQUAD4 excluded
    kinds = sorted(b.elem_type for b in model.blocks)
    assert kinds == ["hex8", "tet4"]
    hex_block = next(b for b in model.blocks if b.elem_type == "hex8")
    np.testing.assert_array_equal(hex_block.conn[0], [0, 1, 2, 3, 4, 5, 6, 7])
    assert model.parts[hex_block.part_id] == "bulk"
    codes = [i["code"] for i in issues]
    assert "NON_SOLID_EXCLUDED" in codes


def test_read_nastran_roundtrip_lossless(tmp_path) -> None:
    path = tmp_path / "plate.bdf"
    path.write_text(DECK, encoding="utf-8")
    model, _ = read_nastran(str(path), "mm")
    back = model.roundtrip()
    np.testing.assert_allclose(back.nodes, model.nodes, rtol=0, atol=0)
    assert [b.conn.tolist() for b in back.blocks] == [b.conn.tolist() for b in model.blocks]


def test_read_nastran_unit_system_required(tmp_path) -> None:
    path = tmp_path / "plate.bdf"
    path.write_text(DECK, encoding="utf-8")
    with pytest.raises(ValueError):
        read_nastran(str(path), "")
