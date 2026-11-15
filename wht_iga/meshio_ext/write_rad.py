"""Radioss starter ``.rad`` mesh writer preserving parts, groups, and surfaces.

Emits ``#RADIOSS STARTER``, ``/NODE``, ``/BRICK``, ``/TETRA4``, ``/TETRA10``,
``/PENTA6``, ``/PART``, ``/GRNOD/NODE/<id>``, ``/GRBRIC/BRIC/<id>``,
``/SURF/SEG/<id>``, and ``/END`` aligned with
:func:`wht_iga.meshio_ext.read_rad.read_rad`.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from wht_iga.meshio_ext.model import HEX8_FACES, TET4_FACES, WEDGE6_FACES, MeshModel
from wht_iga.meshio_ext.write_inp import _merge_extra_sets

_CANON_TO_RAD: dict[str, str] = {
    "hex8": "BRICK",
    "hex20": "BRICK",
    "tet4": "TETRA4",
    "tet10": "TETRA10",
    "wedge6": "PENTA6",
}

_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}


def _wrap_space_ids(ids: Iterable[int], per_line: int = 10) -> list[str]:
    items = [str(int(i)) for i in ids]
    return [" ".join(items[k : k + per_line]) for k in range(0, len(items), per_line)]


def write_rad(
    model: MeshModel,
    path: str | Path,
    extra_sets: Any = None,
) -> Path:
    """Write a :class:`MeshModel` to a Radioss starter ``.rad`` ASCII file.

    Parameters
    ----------
    model:
        Mesh model to serialize (0-based node/element indices become 1-based).
    path:
        Output ``.rad`` file path.
    extra_sets:
        Optional ``SetRegistry`` or iterable of ``SelectionSet`` items.

    Returns
    -------
    Path
        Written output file path.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    node_sets, elem_sets = _merge_extra_sets(model.node_sets, model.elem_sets, extra_sets)

    elem_records: list[tuple[int, int, str, list[int]]] = []
    part_ids_seen: list[int] = []
    cursor = 0
    for block in model.blocks:
        if block.elem_type not in _CANON_TO_RAD:
            raise ValueError(f"Unsupported element type for Radioss .rad: {block.elem_type!r}")
        if block.part_id not in part_ids_seen:
            part_ids_seen.append(block.part_id)
        for row in block.conn:
            cursor += 1
            elem_records.append(
                (cursor, int(block.part_id), block.elem_type, [int(n) + 1 for n in row])
            )

    lines: list[str] = [
        "#RADIOSS STARTER",
        f"$ unit_system={model.unit_system}",
        "/NODE",
    ]
    for idx, xyz in enumerate(model.nodes):
        lines.append(f"{idx + 1} {xyz[0]:.17g} {xyz[1]:.17g} {xyz[2]:.17g}")

    # Group contiguous elements of the same Radioss card
    current_card = ""
    for eid, pid, canon, conn1 in elem_records:
        card = _CANON_TO_RAD[canon]
        if card != current_card:
            lines.append(f"/{card}")
            current_card = card
        lines.append(f"{eid} " + " ".join(str(n) for n in conn1) + f" {pid}")

    all_pids = sorted(set(part_ids_seen) | set(model.parts.keys()))
    for pid in all_pids:
        pname = model.parts.get(pid, f"PART-{pid}")
        lines.append("/PART")
        lines.append(f"{pid} {pname}")

    gid = 0
    for sname, nids in node_sets.items():
        gid += 1
        lines.append(f"/GRNOD/NODE/{gid}")
        lines.append(sname)
        lines.extend(_wrap_space_ids(int(n) + 1 for n in nids))

    for sname, eids in elem_sets.items():
        gid += 1
        lines.append(f"/GRBRIC/BRIC/{gid}")
        lines.append(sname)
        lines.extend(_wrap_space_ids(int(e) + 1 for e in eids))

    sid = 0
    seg_id = 0
    for sname, pairs in model.face_sets.items():
        sid += 1
        lines.append(f"/SURF/SEG/{sid}")
        lines.append(sname)
        for ge, lf in pairs:
            if 0 <= ge < len(elem_records):
                _eid, _pid, canon, conn1 = elem_records[ge]
                faces_tbl = _FACE_DEF_MAP.get(canon, ())
                if 0 <= lf < len(faces_tbl):
                    seg_id += 1
                    fnodes = [conn1[m] for m in faces_tbl[lf]]
                    lines.append(f"{seg_id} " + " ".join(str(n) for n in fnodes))

    lines.append("/END")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


__all__ = ["write_rad"]
