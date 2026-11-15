"""Nastran ``.bdf`` / OptiStruct ``.fem`` mesh writer preserving parts and sets.

Emits high-precision ``GRID*`` (16-char large field) nodes, ``CHEXA``/``CTETRA``/
``CPENTA`` solid elements, ``PSOLID`` + ``MAT1``, HyperMesh ``$$HMNAME COMP`` /
``$$HMNAME SETS`` metadata, ``SET1`` / ``SET3`` (for ``.bdf``) or ``SET``
(for ``.fem``), and ``SURF`` (``ELFACE``) / ``BSURF`` cards matching
:func:`wht_iga.meshio_ext.read_nastran.read_nastran`.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal

from wht_iga.meshio_ext.model import HEX8_FACES, TET4_FACES, WEDGE6_FACES, MeshModel
from wht_iga.meshio_ext.write_inp import _merge_extra_sets

_CANON_TO_NASTRAN: dict[str, str] = {
    "hex8": "CHEXA",
    "hex20": "CHEXA",
    "tet4": "CTETRA",
    "tet10": "CTETRA",
    "wedge6": "CPENTA",
}

_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}


def _compress_thru(ids: Iterable[int]) -> list[str]:
    """Compress sorted positive integers into Nastran ``A THRU B`` tokens."""
    nums = sorted({int(i) for i in ids})
    if not nums:
        return []
    tokens: list[str] = []
    start = nums[0]
    prev = nums[0]
    for x in nums[1:]:
        if x == prev + 1:
            prev = x
            continue
        if prev - start >= 2:
            tokens.extend([str(start), "THRU", str(prev)])
        elif prev > start:
            tokens.extend([str(start), str(prev)])
        else:
            tokens.append(str(start))
        start = prev = x
    if prev - start >= 2:
        tokens.extend([str(start), "THRU", str(prev)])
    elif prev > start:
        tokens.extend([str(start), str(prev)])
    else:
        tokens.append(str(start))
    return tokens


def _format_elem_card_csv(card: str, eid: int, pid: int, nids_1based: list[int]) -> list[str]:
    """Format CHEXA/CTETRA/CPENTA with CSV continuation lines (<= 8 nodes on line 1)."""
    if len(nids_1based) <= 6:
        body = ",".join(str(n) for n in nids_1based)
        return [f"{card},{eid},{pid},{body}"]
    lines = [f"{card},{eid},{pid}," + ",".join(str(n) for n in nids_1based[:6]) + ","]
    rem = nids_1based[6:]
    for k in range(0, len(rem), 8):
        chunk = rem[k : k + 8]
        suffix = "," if (k + 8 < len(rem)) else ""
        lines.append("," + ",".join(str(n) for n in chunk) + suffix)
    return lines


def _format_grid_large(nid: int, x: float, y: float, z: float) -> list[str]:
    """Format exact high-precision ``GRID*`` 16-character large-field Nastran card."""
    try:
        from pyNastran.bdf.field_writer_16 import print_card_16

        return [print_card_16(["GRID", nid, 0, float(x), float(y), float(z)]).rstrip("\n")]
    except ImportError:
        return [f"GRID,{nid},0,{float(x):.16g},{float(y):.16g},{float(z):.16g}"]


def _format_set_tokens_csv(head: str, tokens: list[str], cont_prefix: str = "+") -> list[str]:
    if not tokens:
        return [head]
    lines: list[str] = []
    first = True
    i = 0
    while i < len(tokens):
        # Never split an "A THRU B" triplet across continuation lines.
        take = min(7, len(tokens) - i)
        if i + take < len(tokens) and tokens[i + take] == "THRU" and take > 1:
            take -= 1
        elif tokens[i + take - 1] == "THRU" and take > 2:
            take -= 2
        chunk = tokens[i : i + take]
        i += take
        if first:
            lines.append(f"{head}," + ",".join(chunk))
            first = False
        else:
            lines.append(f"{cont_prefix}," + ",".join(chunk))
    return lines


def write_nastran(
    model: MeshModel,
    path: str | Path,
    format_dialect: Literal["auto", "nastran", "optistruct"] = "auto",
    extra_sets: Any = None,
) -> Path:
    """Write a :class:`MeshModel` to Nastran ``.bdf`` or OptiStruct ``.fem``.

    Parameters
    ----------
    model:
        Mesh model to serialize (0-based node/element indices become 1-based).
    path:
        Output ``.bdf`` or ``.fem`` path.
    format_dialect:
        ``"auto"`` chooses ``"optistruct"`` for ``.fem`` files and ``"nastran"``
        for ``.bdf`` / other extensions.
    extra_sets:
        Optional ``SetRegistry`` or iterable of ``SelectionSet`` objects.

    Returns
    -------
    Path
        Written output file path.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if format_dialect == "auto":
        dialect = "optistruct" if out.suffix.lower() == ".fem" else "nastran"
    else:
        dialect = format_dialect

    node_sets, elem_sets = _merge_extra_sets(model.node_sets, model.elem_sets, extra_sets)

    # Flatten blocks in order and track 1-based eid + connectivity.
    elem_records: list[tuple[int, int, str, list[int]]] = []  # (eid, pid, canon, 1based_conn)
    part_ids_ordered: list[int] = []
    cursor = 0
    for block in model.blocks:
        if block.elem_type not in _CANON_TO_NASTRAN:
            raise ValueError(f"Unsupported element type for Nastran: {block.elem_type!r}")
        if block.part_id not in part_ids_ordered:
            part_ids_ordered.append(block.part_id)
        for row in block.conn:
            cursor += 1
            elem_records.append(
                (cursor, int(block.part_id), block.elem_type, [int(n) + 1 for n in row])
            )

    lines: list[str] = [
        "$ Generated by WHT_IGA meshio_ext.write_nastran",
        f"$ unit_system={model.unit_system}",
    ]

    for pid in sorted(part_ids_ordered):
        pname = model.parts.get(pid, f"PSOLID-{pid}")
        lines.append(f'$$HMNAME COMP \\{pid}\\"{pname}"')

    set_sid = 0
    nset_sids: list[tuple[int, str, list[int]]] = []
    for sname, nids in node_sets.items():
        set_sid += 1
        lines.append(f'$$HMNAME SETS \\{set_sid}\\1\\"{sname}"')
        nset_sids.append((set_sid, sname, [int(n) + 1 for n in nids]))

    eset_sids: list[tuple[int, str, list[int]]] = []
    for sname, eids in elem_sets.items():
        set_sid += 1
        lines.append(f'$$HMNAME SETS \\{set_sid}\\2\\"{sname}"')
        eset_sids.append((set_sid, sname, [int(e) + 1 for e in eids]))

    fset_sids: list[tuple[int, str, list[tuple[int, int]]]] = []
    for sname, pairs in model.face_sets.items():
        set_sid += 1
        lines.append(f'$$HMNAME SETS \\{set_sid}\\3\\"{sname}"')
        fset_sids.append((set_sid, sname, [(int(e), int(lf)) for e, lf in pairs]))

    # Nodes (GRID*)
    for idx, xyz in enumerate(model.nodes):
        lines.extend(_format_grid_large(idx + 1, xyz[0], xyz[1], xyz[2]))

    # Elements
    for eid, pid, canon, conn1 in elem_records:
        card = _CANON_TO_NASTRAN[canon]
        lines.extend(_format_elem_card_csv(card, eid, pid, conn1))

    # Properties & placeholder MAT1 so pyNastran cross-reference succeeds
    for pid in sorted(part_ids_ordered):
        lines.append(f"PSOLID,{pid},{pid}")
        lines.append(f"MAT1,{pid},2.1E5,,0.3,7.85E-9")

    # Node & Element sets
    for sid, _sname, nids1 in nset_sids:
        if not nids1:
            lines.append(f"SET,{sid},GRID,LIST")
            continue
        if dialect == "optistruct":
            toks = [str(n) for n in sorted(set(nids1))]
            lines.extend(_format_set_tokens_csv(f"SET,{sid},GRID,LIST", toks))
        else:
            toks = _compress_thru(nids1)
            lines.extend(_format_set_tokens_csv(f"SET1,{sid}", toks))

    for sid, _sname, eids1 in eset_sids:
        if not eids1:
            lines.append(f"SET,{sid},ELEM,LIST")
            continue
        if dialect == "optistruct":
            toks = [str(e) for e in sorted(set(eids1))]
            lines.extend(_format_set_tokens_csv(f"SET,{sid},ELEM,LIST", toks))
        else:
            toks = _compress_thru(eids1)
            lines.extend(_format_set_tokens_csv(f"SET3,{sid},ELEM", toks))

    # Element face sets (OptiStruct SURF / ELFACE carrying exact element + face nodes)
    for sid, _sname, pairs in fset_sids:
        lines.append(f"SURF,{sid},ELFACE")
        for ge, lf in pairs:
            if 0 <= ge < len(elem_records):
                eid, _pid, canon, conn1 = elem_records[ge]
                faces_tbl = _FACE_DEF_MAP.get(canon, ())
                if 0 <= lf < len(faces_tbl):
                    fnodes = [conn1[m] for m in faces_tbl[lf]]
                    lines.append(f"+,{eid}," + ",".join(str(n) for n in fnodes))

    lines.append("ENDDATA")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


__all__ = ["write_nastran"]
