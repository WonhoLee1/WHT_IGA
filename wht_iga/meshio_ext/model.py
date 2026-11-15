"""Internal mesh model + element normalization (v2.1 section 5.10.2).

Kratos-free. Canonical solid types: ``hex8``, ``tet4``, ``tet10``,
``hex20``, ``wedge6``. Every reader maps its native node order onto the
canonical order through the explicit ``NODE_ORDER`` tables below
(identity where the conventions agree -- verified numerically in
``tests/unit/test_meshio_model.py`` by checking positive corner volumes
and edge-midpoint positions on reference elements).

Canonical corner coordinates (unit reference elements):

* hex8: 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
         4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
* tet4: 0:(0,0,0) 1:(1,0,0) 2:(0,1,0) 3:(0,0,1)
* tet10: corners as tet4, mids 4:(0-1) 5:(1-2) 6:(2-0) 7:(0-3) 8:(1-3) 9:(2-3)
* hex20: corners as hex8, mids 8:(0-1) 9:(1-2) 10:(2-3) 11:(3-0)
         12:(4-5) 13:(5-6) 14:(6-7) 15:(7-4) 16:(0-4) 17:(1-5) 18:(2-6) 19:(3-7)
* wedge6: 0:(0,0,0) 1:(1,0,0) 2:(0,1,0) 3:(0,0,1) 4:(1,0,1) 5:(0,1,1)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

Severity = Literal["info", "warning", "error"]

#: canonical element type -> node count
CANONICAL_COUNTS: dict[str, int] = {
    "hex8": 8,
    "tet4": 4,
    "tet10": 10,
    "hex20": 20,
    "wedge6": 6,
}

#: native element token (upper-cased, format-specific) -> canonical type.
#: Tokens not listed here are unsupported (caller must issue + exclude).
NORMALIZE: dict[str, str] = {
    # Abaqus .inp C3D family
    "C3D8": "hex8", "C3D8R": "hex8", "C3D8H": "hex8", "C3D8I": "hex8",
    "C3D4": "tet4", "C3D4H": "tet4",
    "C3D10": "tet10", "C3D10M": "tet10",
    "C3D20": "hex20", "C3D20R": "hex20", "C3D20H": "hex20",
    "C3D6": "wedge6",
    # Nastran CHEXA / CTETRA / CPENTA
    "CHEXA": "hex8", "CHEXA8": "hex8", "CHEXA20": "hex20",
    "CTETRA": "tet4", "CTETRA4": "tet4", "CTETRA10": "tet10",
    "CPENTA": "wedge6", "CPENTA6": "wedge6",
    # LS-DYNA *ELEMENT_SOLID (by node count, resolved by reader)
    "SOLID8": "hex8", "SOLID4": "tet4", "SOLID10": "tet10",
    "SOLID20": "hex20", "SOLID6": "wedge6",
    # Radioss starter
    "BRICK": "hex8", "BRICK8": "hex8", "BRICK20": "hex20",
    "TETRA4": "tet4", "TETRA10": "tet10", "PENTA6": "wedge6",
    # meshio / Gmsh cell types
    "HEXAHEDRON": "hex8", "HEXAHEDRON20": "hex20",
    "TETRA": "tet4", "TETRA10": "tet10", "WEDGE": "wedge6",
    # explicit canonical names (pass-through)
    "HEX8": "hex8", "TET4": "tet4", "TET10": "tet10",
    "HEX20": "hex20", "WEDGE6": "wedge6",
}

#: shell/beam/spring tokens: read (counted) then issued + excluded (K27).
NON_SOLID_TOKENS: dict[str, str] = {
    "S3": "shell", "S4": "shell", "S4R": "shell", "STRI3": "shell",
    "CQUAD4": "shell", "CTRIA3": "shell", "CQUAD8": "shell", "CTRIA6": "shell",
    "SHELL": "shell", "QUAD": "shell", "TRIA": "shell",
    "B31": "beam", "BEAM": "beam", "CBAR": "beam", "CBEAM": "beam",
    "CROD": "beam", "TRUSS": "beam",
    "SPRINGA": "spring", "SPRING1": "spring", "SPRING2": "spring",
    "CELAS1": "spring", "CELAS2": "spring", "CONM2": "mass",
}

_IDENTITY_8 = [0, 1, 2, 3, 4, 5, 6, 7]
_IDENTITY_4 = [0, 1, 2, 3]
_IDENTITY_10 = list(range(10))
_IDENTITY_20 = list(range(20))
_IDENTITY_6 = list(range(6))

#: (format, canonical) -> permutation raw -> canonical.
#: canonical[i] = raw[perm[i]]. All identity today (conventions agree);
#: the table exists so a future dialect difference is a one-line fix.
NODE_ORDER: dict[tuple[str, str], list[int]] = {
    ("inp", "hex8"): list(_IDENTITY_8), ("inp", "tet4"): list(_IDENTITY_4),
    ("inp", "tet10"): list(_IDENTITY_10), ("inp", "hex20"): list(_IDENTITY_20),
    ("inp", "wedge6"): list(_IDENTITY_6),
    ("nastran", "hex8"): list(_IDENTITY_8), ("nastran", "tet4"): list(_IDENTITY_4),
    ("nastran", "tet10"): list(_IDENTITY_10), ("nastran", "hex20"): list(_IDENTITY_20),
    ("nastran", "wedge6"): list(_IDENTITY_6),
    ("k", "hex8"): list(_IDENTITY_8), ("k", "tet4"): list(_IDENTITY_4),
    ("k", "tet10"): list(_IDENTITY_10), ("k", "hex20"): list(_IDENTITY_20),
    ("k", "wedge6"): list(_IDENTITY_6),
    ("rad", "hex8"): list(_IDENTITY_8), ("rad", "tet4"): list(_IDENTITY_4),
    ("rad", "tet10"): list(_IDENTITY_10), ("rad", "hex20"): list(_IDENTITY_20),
    ("rad", "wedge6"): list(_IDENTITY_6),
    ("msh", "hex8"): list(_IDENTITY_8), ("msh", "tet4"): list(_IDENTITY_4),
    ("msh", "tet10"): list(_IDENTITY_10), ("msh", "hex20"): list(_IDENTITY_20),
    ("msh", "wedge6"): list(_IDENTITY_6),
}


def make_issue(severity: Severity, code: str, message: str, entity: str = "") -> dict[str, str]:
    """Issue record (same shape as Track GEO issues)."""
    return {"severity": severity, "code": code, "message": message, "entity": entity}


def require_unit_system(unit_system: str | None) -> str:
    """Enforce the mandatory unit_system argument (K26: never guess units)."""
    if not unit_system or not str(unit_system).strip():
        raise ValueError("unit_system is required (e.g. 'mm', 'm', 'tonne_mm_s'); refusing to guess units")
    return str(unit_system)


def apply_order(conn: list[int], fmt: str, canonical: str) -> list[int]:
    """Map raw connectivity onto canonical node order."""
    perm = NODE_ORDER[(fmt, canonical)]
    if len(conn) != len(perm):
        raise ValueError(f"{fmt}/{canonical}: connectivity length {len(conn)} != {len(perm)}")
    return [conn[p] for p in perm]


@dataclass
class ElemBlock:
    """One homogeneous element block."""

    elem_type: str  # canonical: hex8 | tet4 | tet10 | hex20 | wedge6
    conn: np.ndarray  # (n_elem, n_nodes) 0-based node indices
    part_id: int = 1

    def __post_init__(self) -> None:
        if self.elem_type not in CANONICAL_COUNTS:
            raise ValueError(f"non-canonical elem_type {self.elem_type!r}")
        arr = np.asarray(self.conn, dtype=np.int64)
        if arr.ndim != 2 or arr.shape[1] != CANONICAL_COUNTS[self.elem_type]:
            raise ValueError(f"{self.elem_type} conn shape {arr.shape}")
        self.conn = arr


@dataclass
class MeshModel:
    """Internal mesh model (v2.1 section 5.10.2). Mesh + sets only."""

    nodes: np.ndarray  # (n_nodes, 3) float
    blocks: list[ElemBlock] = field(default_factory=list)
    parts: dict[int, str] = field(default_factory=dict)  # part_id -> name
    node_sets: dict[str, list[int]] = field(default_factory=dict)
    elem_sets: dict[str, list[int]] = field(default_factory=dict)
    face_sets: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    frames: dict[str, list[float]] = field(default_factory=dict)
    unit_system: str = ""

    def __post_init__(self) -> None:
        self.nodes = np.asarray(self.nodes, dtype=float).reshape(-1, 3)

    @property
    def n_elem(self) -> int:
        """Total element count across blocks."""
        return int(sum(b.conn.shape[0] for b in self.blocks))

    def to_dict(self) -> dict[str, Any]:
        """Lossless JSON-safe serialization (roundtrip acceptance)."""
        return {
            "unit_system": self.unit_system,
            "nodes": np.asarray(self.nodes, dtype=float).tolist(),
            "blocks": [
                {"elem_type": b.elem_type, "conn": np.asarray(b.conn, dtype=int).tolist(),
                 "part_id": int(b.part_id)}
                for b in self.blocks
            ],
            "parts": {str(k): v for k, v in self.parts.items()},
            "node_sets": {k: [int(i) for i in v] for k, v in self.node_sets.items()},
            "elem_sets": {k: [int(i) for i in v] for k, v in self.elem_sets.items()},
            "face_sets": {k: [[int(e), int(f)] for e, f in v] for k, v in self.face_sets.items()},
            "frames": {k: [float(x) for x in v] for k, v in self.frames.items()},
        }

    @staticmethod
    def from_dict(data: dict[str, Any]) -> MeshModel:
        """Inverse of to_dict (exact)."""
        return MeshModel(
            nodes=np.asarray(data["nodes"], dtype=float),
            blocks=[ElemBlock(elem_type=b["elem_type"],
                              conn=np.asarray(b["conn"], dtype=np.int64),
                              part_id=int(b.get("part_id", 1))) for b in data.get("blocks", [])],
            parts={int(k): str(v) for k, v in data.get("parts", {}).items()},
            node_sets={k: [int(i) for i in v] for k, v in data.get("node_sets", {}).items()},
            elem_sets={k: [int(i) for i in v] for k, v in data.get("elem_sets", {}).items()},
            face_sets={k: [(int(e), int(f)) for e, f in v] for k, v in data.get("face_sets", {}).items()},
            frames={k: [float(x) for x in v] for k, v in data.get("frames", {}).items()},
            unit_system=str(data.get("unit_system", "")),
        )

    def roundtrip(self) -> MeshModel:
        """Serialize -> deserialize (acceptance: lossless)."""
        return MeshModel.from_dict(self.to_dict())


#: canonical outward faces (unit reference elements above).
#: quads listed as 4-tuples; boundary extraction triangulates them.
HEX8_FACES: list[tuple[int, ...]] = [
    (0, 3, 2, 1),  # z = 0
    (4, 5, 6, 7),  # z = 1
    (0, 4, 7, 3),  # x = 0
    (1, 2, 6, 5),  # x = 1
    (0, 1, 5, 4),  # y = 0
    (2, 3, 7, 6),  # y = 1
]
TET4_FACES: list[tuple[int, ...]] = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]
WEDGE6_FACES: list[tuple[int, ...]] = [
    (0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (0, 3, 5, 2), (1, 2, 5, 4),
]
_FACE_TABLE: dict[str, list[tuple[int, ...]]] = {
    "hex8": HEX8_FACES, "tet4": TET4_FACES, "wedge6": WEDGE6_FACES,
}


def outer_faces(model: MeshModel) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Boundary triangulation of hex8/tet4/wedge6 blocks (use-case 3).

    Returns (nodes, tris, notes) with outward winding. tet10/hex20 blocks
    fall back to their corner nodes as linear elements (note recorded).
    Interior faces shared by two cells cancel exactly.
    """
    notes: list[str] = []
    counts: dict[tuple[int, ...], int] = {}
    keep: dict[tuple[int, ...], tuple[int, ...]] = {}
    for b in model.blocks:
        if b.elem_type in ("tet10", "hex20"):
            corner = {"tet10": [0, 1, 2, 3], "hex20": [0, 1, 2, 3, 4, 5, 6, 7]}[b.elem_type]
            lin = "tet4" if b.elem_type == "tet10" else "hex8"
            sub = b.conn[:, corner]
            notes.append(f"{b.elem_type} linearized to {lin} on boundary")
        elif b.elem_type in _FACE_TABLE:
            sub = b.conn
            lin = b.elem_type
        else:  # pragma: no cover - constructor already guards
            raise ValueError(b.elem_type)
        for elem in sub:
            for face in _FACE_TABLE[lin]:
                nodes_f = tuple(int(elem[i]) for i in face)
                key = tuple(sorted(nodes_f))
                counts[key] = counts.get(key, 0) + 1
                keep[key] = nodes_f
    tris: list[tuple[int, int, int]] = []
    for key, n in counts.items():
        if n == 1:
            face = keep[key]
            if len(face) == 3:
                tris.append((face[0], face[1], face[2]))
            else:
                tris.append((face[0], face[1], face[2]))
                tris.append((face[0], face[2], face[3]))
    return (np.asarray(model.nodes, dtype=float),
            np.asarray(tris, dtype=np.int64), notes)
