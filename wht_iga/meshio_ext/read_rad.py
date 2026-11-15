"""Radioss starter ``.rad`` mesh reader (lightweight parser, mesh + sets only).

Documented subset (free-format fields; fixed 10-char columns also split
cleanly, so both layouts parse):

* ``/NODE`` -- ``nid x y z``
* ``/BRICK`` -- ``bid n1..n8 pid`` (8-node hex; 20-node ``bid`` lines with
  20 nodes + pid map to hex20)
* ``/TETRA4`` -- ``tid n1 n2 n3 n4 pid``; ``/TETRA10``; ``/PENTA6``
* ``/PART`` -- ``pid`` on the first data line (title kept as part name)
* ``/GRNOD/PART/...`` -- node ids -> node_sets; ``/GRBRIC`` -- brick ids
  -> elem_sets (by position index)

Card names are matched case-insensitively on the part before any ``/``
option (``/GRNOD/PART/1`` -> ``GRNOD``). Anything else (materials,
properties, boundary conditions, ...) is aggregated into one
``UNSUPPORTED_DECK_KEYWORD`` info issue (K25). Shell/beam/spring element
blocks are counted, issued and excluded (K27).
"""

from __future__ import annotations

import os

import numpy as np

from wht_iga.meshio_ext.model import (
    HEX8_FACES,
    TET4_FACES,
    WEDGE6_FACES,
    ElemBlock,
    MeshModel,
    apply_order,
    make_issue,
    require_unit_system,
)

_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}


def _cyclic_match(seq_a: tuple[int, ...], seq_b: tuple[int, ...]) -> bool:
    """Return True if ``seq_a`` is a cyclic rotation of ``seq_b`` (preserving winding)."""
    n = len(seq_a)
    if n != len(seq_b) or n == 0:
        return False
    return any(seq_a == seq_b[k:] + seq_b[:k] for k in range(n))

_SOLID_CARDS = {
    "BRICK": ("rad", ("hex8", "hex20")),
    "TETRA4": ("rad", ("tet4",)),
    "TETRA10": ("rad", ("tet10",)),
    "PENTA6": ("rad", ("wedge6",)),
}
_NON_SOLID_CARDS = ("SHELL", "SH3N", "BEAM", "TRUSS", "SPRING", "MASS")


def _nums(line: str) -> list[float]:
    return [float(t) for t in line.replace(",", " ").split()]


def read_rad(path: str, unit_system: str) -> tuple[MeshModel, list[dict]]:
    """Read a Radioss starter ``.rad`` mesh file. Returns (model, issues)."""
    require_unit_system(unit_system)
    issues: list[dict] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        content = [ln.rstrip("\n") for ln in f]

    nodes: dict[int, list[float]] = {}
    order: list[int] = []
    elems: list[tuple[str, int, list[int]]] = []  # (canon, pid, conn)
    eorder: list[int] = []
    parts: dict[int, str] = {}
    nsets: dict[str, list[int]] = {}
    esets: dict[str, list[int]] = {}
    surf_segs: dict[str, list[list[int]]] = {}
    skipped: dict[str, int] = {}
    nonsolid: dict[str, int] = {}
    # pending_group: (kind: "node"|"elem"|"surf", name, items, expect_title)
    pending_group: list | None = None

    card = ""

    def card_parts(line: str) -> list[str]:
        head = line.strip()[1:]
        return [p.strip().upper() for p in head.split("/") if p.strip()]

    def flush_group() -> None:
        nonlocal pending_group
        if pending_group is None:
            return
        gkind, gname, items, _ = pending_group
        if gkind == "node":
            nsets.setdefault(gname, []).extend(items)
        elif gkind == "elem":
            esets.setdefault(gname, []).extend(items)
        elif gkind == "surf":
            surf_segs.setdefault(gname, []).extend(items)
        pending_group = None

    for ln in content:
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("$"):
            continue
        if s.startswith("/"):
            flush_group()
            parts_tok = card_parts(s)
            card = parts_tok[0] if parts_tok else ""
            subcard = parts_tok[1] if len(parts_tok) > 1 else ""
            if card == "END":
                card = ""
            elif card in ("NODE", "BRICK", "TETRA4", "TETRA10", "PENTA6", "PART"):
                pass
            elif card == "GRNOD":
                if subcard == "NODE":
                    pending_group = ["node", "", [], True]
                else:
                    pending_group = [
                        "node",
                        f"{card}_{len(nsets) + len(esets)}",
                        [],
                        False,
                    ]
            elif card == "GRBRIC":
                if subcard == "BRIC":
                    pending_group = ["elem", "", [], True]
                else:
                    pending_group = [
                        "elem",
                        f"{card}_{len(nsets) + len(esets)}",
                        [],
                        False,
                    ]
            elif card == "SURF":
                pending_group = ["surf", f"SURF_{len(surf_segs)}", [], subcard == "SEG"]
            elif card.startswith(_NON_SOLID_CARDS):
                nonsolid[card] = nonsolid.get(card, 0)
            else:
                skipped[card] = skipped.get(card, 0) + 1
                card = ""
            continue
        if card == "NODE":
            try:
                v = _nums(s)
            except ValueError:
                continue
            if len(v) >= 4:
                nid = int(v[0])
                if nid not in nodes:
                    order.append(nid)
                nodes[nid] = [v[1], v[2], v[3]]
        elif card in _SOLID_CARDS:
            fmt, _ = _SOLID_CARDS[card]
            try:
                v = _nums(s)
            except ValueError:
                continue
            ids = [int(x) for x in v]
            if len(ids) < 5:
                continue
            eid = ids[0]
            if card == "BRICK":
                if len(ids) >= 22:  # bid + 20 nodes + pid
                    canon, conn, pid = "hex20", ids[1:21], ids[21]
                else:  # bid + 8 nodes [+ pid]
                    canon = "hex8"
                    conn = ids[1:9]
                    pid = ids[9] if len(ids) > 9 else 1
            elif card == "TETRA10":
                canon, conn = "tet10", ids[1:11]
                pid = ids[11] if len(ids) > 11 else 1
            elif card == "TETRA4":
                canon, conn = "tet4", ids[1:5]
                pid = ids[5] if len(ids) > 5 else 1
            else:  # PENTA6
                canon, conn = "wedge6", ids[1:7]
                pid = ids[7] if len(ids) > 7 else 1
            elems.append((canon, pid, apply_order(conn, fmt, canon)))
            eorder.append(eid)
        elif card == "PART":
            v = s.split()
            try:
                pid = int(float(v[0]))
            except (ValueError, IndexError):
                continue
            parts[pid] = " ".join(v[1:]) or f"PART-{pid}"
            card = ""
        elif pending_group is not None:
            if pending_group[3]:
                pending_group[1] = s
                pending_group[3] = False
                continue
            gkind = pending_group[0]
            if gkind == "surf":
                try:
                    segs = [int(float(t)) for t in s.replace(",", " ").split()]
                except ValueError:
                    continue
                if len(segs) >= 4:
                    # Segment format: seg_id n1 n2 n3 [n4] (strip seg_id when >= 4 tokens)
                    seg_nodes = segs[1:] if len(segs) >= 4 else segs
                    pending_group[2].append(seg_nodes)
            else:
                try:
                    pending_group[2].extend(
                        int(float(t)) for t in s.replace(",", " ").split()
                    )
                except ValueError:
                    continue
        elif card and card in nonsolid:
            nonsolid[card] += 1

    flush_group()

    if nonsolid:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(nonsolid.items()))
        issues.append(make_issue("warning", "NON_SOLID_EXCLUDED",
                                 f"non-solid elements excluded: {summary}",
                                 entity=f"rad:{os.path.basename(path)}"))
    if skipped:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(skipped.items()))
        issues.append(make_issue("info", "UNSUPPORTED_DECK_KEYWORD",
                                 f"analysis-deck keywords ignored (mesh+sets scope): {summary}",
                                 entity=f"rad:{os.path.basename(path)}"))

    idx = {nid: k for k, nid in enumerate(order)}
    nodes_arr = np.array([nodes[nid] for nid in order], dtype=float)
    part_ids = sorted({p for _, p, _ in elems} | set(parts))
    pmap = {p: i + 1 for i, p in enumerate(part_ids)} if part_ids else {}
    parts_out = {pmap[p]: parts.get(p, f"PART-{p}") for p in part_ids} or {1: "PART-1"}
    by_type_indices: dict[tuple[str, int], list[int]] = {}
    for k, (canon, pid, _conn) in enumerate(elems):
        by_type_indices.setdefault((canon, pmap.get(pid, 1)), []).append(k)

    blocks: list[ElemBlock] = []
    raw_k_to_model_idx: dict[int, int] = {}
    elem_conn0_list: list[tuple[str, list[int]]] = [("", [])] * len(elems)
    cursor = 0
    for (canon, p), raw_indices in sorted(by_type_indices.items()):
        rows: list[list[int]] = []
        for k in raw_indices:
            conn0 = [idx[v] for v in elems[k][2]]
            rows.append(conn0)
            raw_k_to_model_idx[k] = cursor
            elem_conn0_list[cursor] = (canon, conn0)
            cursor += 1
        blocks.append(ElemBlock(elem_type=canon, conn=np.asarray(rows, dtype=np.int64), part_id=p))

    node_sets = {k: sorted(idx[v] for v in vs if v in idx) for k, vs in nsets.items()}
    elem_pos = {eid: raw_k_to_model_idx[k] for k, eid in enumerate(eorder)}
    elem_sets = {k: sorted(elem_pos[v] for v in vs if v in elem_pos) for k, vs in esets.items()}

    # Build lookup for /SURF/SEG segment -> (model_elem_idx, local_face)
    face_sets: dict[str, list[tuple[int, int]]] = {}
    if surf_segs:
        oriented_lookup: dict[tuple[int, ...], tuple[int, int]] = {}
        unordered_lookup: dict[tuple[int, ...], list[tuple[int, int]]] = {}
        for ge, (canon, conn0) in enumerate(elem_conn0_list):
            for lf, f_locals in enumerate(_FACE_DEF_MAP.get(canon, ())):
                ftup = tuple(conn0[m] for m in f_locals)
                oriented_lookup[ftup] = (ge, lf)
                unordered_lookup.setdefault(tuple(sorted(ftup)), []).append((ge, lf))

        for sname, seg_list in surf_segs.items():
            pairs: list[tuple[int, int]] = []
            for raw_seg in seg_list:
                seg0 = tuple(idx[v] for v in raw_seg if v in idx)
                if len(seg0) < 3:
                    continue
                # Try exact cyclic winding match first (disambiguates shared interior faces).
                matched: tuple[int, int] | None = None
                cands = unordered_lookup.get(tuple(sorted(seg0)), [])
                for ge, lf in cands:
                    canon, conn0 = elem_conn0_list[ge]
                    ftup = tuple(conn0[m] for m in _FACE_DEF_MAP[canon][lf])
                    if _cyclic_match(seg0, ftup):
                        matched = (ge, lf)
                        break
                if matched is None and cands:
                    matched = cands[0]
                if matched is not None:
                    pairs.append(matched)
            face_sets[sname] = pairs

    return (MeshModel(nodes=nodes_arr, blocks=blocks, parts=parts_out,
                      node_sets=node_sets, elem_sets=elem_sets,
                      face_sets=face_sets,
                      unit_system=unit_system), issues)
