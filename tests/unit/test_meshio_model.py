"""MeshModel: normalization tables, JSON roundtrip, boundary extraction."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.meshio_ext.model import (
    CANONICAL_COUNTS,
    HEX8_FACES,
    NON_SOLID_TOKENS,
    NORMALIZE,
    TET4_FACES,
    WEDGE6_FACES,
    ElemBlock,
    MeshModel,
    apply_order,
    outer_faces,
    require_unit_system,
)

REF_CORNERS: dict[str, np.ndarray] = {
    "hex8": np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    ),
    "tet4": np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    "wedge6": np.array(
        [[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]
    ),
}
REF_FACES = {"hex8": HEX8_FACES, "tet4": TET4_FACES, "wedge6": WEDGE6_FACES}


def _tri_normal(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    n = np.cross(b - a, c - a)
    return n / float(np.linalg.norm(n))


def test_normalize_tokens() -> None:
    assert NORMALIZE["C3D8"] == "hex8"
    assert NORMALIZE["C3D20R"] == "hex20"
    assert NORMALIZE["C3D10"] == "tet10"
    assert NORMALIZE["C3D6"] == "wedge6"
    assert NORMALIZE["CHEXA20"] == "hex20"
    assert NORMALIZE["CTETRA"] == "tet4"
    assert NORMALIZE["CPENTA"] == "wedge6"
    assert NORMALIZE["BRICK"] == "hex8"
    assert NORMALIZE["HEXAHEDRON"] == "hex8"
    # shells / beams / springs are NOT normalizable (K27: issue + exclude)
    assert NORMALIZE.get("S4R") is None
    assert NON_SOLID_TOKENS["S4R"] == "shell"
    assert NON_SOLID_TOKENS["CBAR"] == "beam"
    assert NON_SOLID_TOKENS["CELAS1"] == "spring"
    for canon, count in CANONICAL_COUNTS.items():
        assert count == {"hex8": 8, "tet4": 4, "tet10": 10, "hex20": 20, "wedge6": 6}[canon]


def test_apply_order_identity_on_reference() -> None:
    for canon in ("hex8", "tet4", "tet10", "hex20", "wedge6"):
        conn = list(range(CANONICAL_COUNTS[canon]))
        assert apply_order(conn, "inp", canon) == conn
    with pytest.raises(ValueError):
        apply_order([0, 1, 2], "inp", "hex8")


def test_face_tables_outward_on_reference() -> None:
    for canon, faces in REF_FACES.items():
        xyz = REF_CORNERS[canon]
        center = xyz.mean(axis=0)
        for face in faces:
            pts = xyz[list(face)]
            n = _tri_normal(pts[0], pts[1], pts[2])
            outward = pts.mean(axis=0) - center
            assert float(n @ outward) > 0.0, (canon, face, n)


def test_require_unit_system() -> None:
    assert require_unit_system("mm") == "mm"
    with pytest.raises(ValueError):
        require_unit_system("")
    with pytest.raises(ValueError):
        require_unit_system(None)  # type: ignore[arg-type]


def _demo_model() -> MeshModel:
    nodes = np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 0, 1],
            [2, 1, 1],
        ]
    )
    hex_conn = np.array([[0, 1, 2, 3, 4, 5, 6, 7], [1, 8, 9, 2, 5, 10, 11, 6]])
    tet_conn = np.array([[1, 2, 3, 0]])
    return MeshModel(
        nodes=nodes,
        blocks=[ElemBlock("hex8", hex_conn, 1), ElemBlock("tet4", tet_conn, 2)],
        parts={1: "bulk", 2: "fill"},
        node_sets={"loaded": [1, 2]},
        elem_sets={"all_hex": [0, 1]},
        face_sets={"wall": [(0, 5)]},
        frames={"global": [0.0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]},
        unit_system="mm",
    )


def test_json_roundtrip_lossless() -> None:
    model = _demo_model()
    back = model.roundtrip()
    np.testing.assert_array_equal(back.nodes, model.nodes)
    assert back.unit_system == "mm"
    assert back.parts == model.parts
    assert back.node_sets == model.node_sets
    assert back.elem_sets == model.elem_sets
    assert back.face_sets == model.face_sets
    assert back.frames == model.frames
    assert len(back.blocks) == 2
    for a, b in zip(back.blocks, model.blocks):
        assert a.elem_type == b.elem_type and a.part_id == b.part_id
        np.testing.assert_array_equal(a.conn, b.conn)
    assert back.n_elem == 3


def test_outer_faces_two_hexes() -> None:
    nodes = np.array(
        [
            [0.0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 0, 1],
            [2, 1, 1],
        ]
    )
    conn = np.array([[0, 1, 2, 3, 4, 5, 6, 7], [1, 8, 9, 2, 5, 10, 11, 6]])
    model = MeshModel(nodes=nodes, blocks=[ElemBlock("hex8", conn, 1)], unit_system="m")
    tnodes, tris, notes = outer_faces(model)
    assert notes == []
    # 2x1x1 brick: surface 2*(2+2+1) = 10, shared face cancels: 10 quads -> 20 tris
    assert tris.shape == (20, 3)
    area = 0.0
    center = np.array([1.0, 0.5, 0.5])
    for t in tris:
        a, b, c = tnodes[t]
        n = np.cross(b - a, c - a)
        area += float(np.linalg.norm(n)) / 2.0
        assert float(n @ ((a + b + c) / 3.0 - center)) > 0.0
    assert area == pytest.approx(10.0, rel=1e-12)
