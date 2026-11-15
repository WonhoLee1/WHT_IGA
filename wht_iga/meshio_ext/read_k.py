"""LS-DYNA ``.k`` / ``.key`` mesh reader (lightweight parser, mesh + sets only).

Supported: ``*NODE``, ``*ELEMENT_SOLID`` (solid type by trailing-field
count: 8 -> hex8, 4 -> tet4, 6 -> wedge6, 10 -> tet10, 20 -> hex20),
``*PART`` (``pid`` -> part mapping), ``*SET_NODE_LIST``,
``*SET_SOLID_LIST`` / ``*SET_SOLID``, ``*SET_SEGMENT`` (mapped to
``face_sets`` as ``(elem, local_face)`` via exact node-set match on hex8
blocks), ``*INCLUDE`` (followed relative to the including file, cycle
guarded). ``$`` comments skipped.

Unsupported cards (materials, loads, contacts, ...) are aggregated into
one ``UNSUPPORTED_DECK_KEYWORD`` info issue (K25). Shell/beam/spring
elements are counted, issued and excluded (K27).
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

_SOLID_BY_COUNT = {8: "hex8", 4: "tet4", 6: "wedge6", 10: "tet10", 20: "hex20"}
_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}

_NON_SOLID_CARDS = ("*ELEMENT_SHELL", "*ELEMENT_BEAM", "*ELEMENT_DISCRETE",
                    "*ELEMENT_MASS", "*ELEMENT_SEATBELT", "*CONSTRAINED")


def _ints(line: str) -> list[int]:
    toks = line.replace(",", " ").split()
    out: list[int] = []
    for t in toks:
        out.append(int(t))
    return out


def _cyclic_match(seq_a: tuple[int, ...], seq_b: tuple[int, ...]) -> bool:
    """Return True if ``seq_a`` is a cyclic rotation of ``seq_b``."""
    n = len(seq_a)
    if n != len(seq_b) or n == 0:
        return False
    return any(seq_a == seq_b[k:] + seq_b[:k] for k in range(n))


def _floats(line: str) -> list[float]:
    return [float(t) for t in line.replace(",", " ").split()]


def _is_numeric_line(line: str) -> bool:
    toks = line.replace(",", " ").split()
    if not toks:
        return False
    try:
        [float(t) for t in toks]
        return True
    except ValueError:
        return False


def read_k(path: str, unit_system: str) -> tuple[MeshModel, list[dict]]:
    """Read an LS-DYNA ``.k`` mesh file. Returns (model, issues)."""
    require_unit_system(unit_system)
    issues: list[dict] = []
    seen: set[str] = set()
    acc: dict = {"nodes": {}, "order": [], "elems": [], "eorder": [],
                 "parts": {}, "nsets": {}, "esets": {}, "segments": {},
                 "skipped": {}, "nonsolid": {}, "files": []}

    def ingest(p: str) -> None:
        key = os.path.abspath(p)
        if key in seen:
            issues.append(make_issue("warning", "INCLUDE_CYCLE",
                                     f"*INCLUDE cycle guarded at {p}",
                                     entity=f"k:{os.path.basename(p)}"))
            return
        seen.add(key)
        acc["files"].append(p)
        with open(p, encoding="utf-8", errors="replace") as f:
            content = [ln.rstrip("\n") for ln in f]
        base = os.path.dirname(key)
        card: str | None = None
        pending_set: list | None = None  # [kind, name, list_or_segs, expect_title]
        pending_part_title: str | None = None

        def flush_pending() -> None:
            nonlocal pending_set
            if pending_set is None:
                return
            kind, sname, payload, _ = pending_set
            if kind == "node":
                acc["nsets"].setdefault(sname, []).extend(payload)
            elif kind == "elem":
                acc["esets"].setdefault(sname, []).extend(payload)
            elif kind == "seg":
                acc["segments"].setdefault(sname, []).extend(payload)
            pending_set = None

        for ln in content:
            s = ln.strip()
            if not s or s.startswith("$"):
                continue
            if s.startswith("*"):
                flush_pending()
                pending_part_title = None
                head = s.split()[0].upper()
                if head == "*INCLUDE":
                    card = "INCLUDE_PATH"
                    continue
                if head == "*END":
                    card = None
                    continue
                card = head
                if head == "*KEYWORD":
                    acc["skipped"]["*KEYWORD"] = acc["skipped"].get("*KEYWORD", 0) + 1
                    card = None
                elif head.startswith("*TITLE"):
                    card = "TITLE"
                elif head.startswith("*SET_NODE_LIST"):
                    if head.endswith("_TITLE"):
                        pending_set = ["node", "", [], True]
                    else:
                        pending_set = ["node", f"SET_{len(acc['nsets']) + len(acc['esets'])}", [], False]
                    card = "SET_IDS"
                elif head.startswith("*SET_SOLID"):
                    if head.endswith("_TITLE"):
                        pending_set = ["elem", "", [], True]
                    else:
                        pending_set = ["elem", f"SET_{len(acc['nsets']) + len(acc['esets'])}", [], False]
                    card = "SET_IDS"
                elif head.startswith("*SET_SEGMENT"):
                    if head.endswith("_TITLE"):
                        pending_set = ["seg", "", [], True]
                    else:
                        pending_set = ["seg", "SEGMENT", [], False]
                    card = "SEGMENT"
                elif head.startswith("*PART"):
                    card = "PART_PID"
                elif head.startswith("*NODE"):
                    card = "NODE"
                elif head.startswith("*ELEMENT_SOLID"):
                    card = "ESOLID"
                elif head.startswith(_NON_SOLID_CARDS):
                    card = "NONSOLID:" + head
                    acc["nonsolid"][card] = acc["nonsolid"].get(card, 0)
                elif head.startswith(("*MAT", "*SECTION", "*CONTROL", "*BOUNDARY",
                                      "*LOAD", "*CONTACT", "*DATABASE", "*DEFINE",
                                      "*INITIAL", "*HOURGLASS", "*DAMPING")):
                    acc["skipped"][head.split()[0]] = acc["skipped"].get(head.split()[0], 0) + 1
                    card = None
                else:
                    acc["skipped"][head] = acc["skipped"].get(head, 0) + 1
                    card = None
                continue
            # data line
            if card == "TITLE":
                card = None
            elif card == "PART_PID":
                if _is_numeric_line(s):
                    vals = _ints(s)
                    if vals:
                        pid = vals[0]
                        if pending_part_title is not None:
                            acc["parts"][pid] = pending_part_title
                            pending_part_title = None
                            card = None
                        else:
                            acc["parts"][pid] = f"PART-{pid}"
                            acc["_last_part"] = pid
                            card = "PART_TITLE"
                else:
                    pending_part_title = s.strip()
            elif card == "PART_TITLE":
                pid = acc.pop("_last_part", None)
                if pid is not None and not _is_numeric_line(s):
                    acc["parts"][pid] = s.strip()
                card = None
            elif card == "NODE":
                try:
                    v = _floats(s)
                except ValueError:
                    continue
                if len(v) >= 4:
                    nid = int(v[0])
                    if nid not in acc["nodes"]:
                        acc["order"].append(nid)
                    acc["nodes"][nid] = [v[1], v[2], v[3]]
            elif card == "ESOLID":
                try:
                    vals = _ints(s)
                except ValueError:
                    continue
                if len(vals) >= 3:
                    eid, pid = vals[0], vals[1]
                    conn = [c for c in vals[2:] if c != 0]
                    canon = _SOLID_BY_COUNT.get(len(conn))
                    if canon is None:
                        acc["skipped"][f"ESOLID_N{len(conn)}"] = \
                            acc["skipped"].get(f"ESOLID_N{len(conn)}", 0) + 1
                        issues.append(make_issue(
                            "warning", "NON_SOLID_EXCLUDED",
                            f"*ELEMENT_SOLID eid={eid} with {len(conn)} nodes excluded",
                            entity=f"k:{os.path.basename(p)}"))
                    else:
                        acc["elems"].append((canon, pid, apply_order(conn, "k", canon)))
                        acc["eorder"].append(eid)
            elif card == "SET_IDS" and pending_set is not None:
                if pending_set[3]:
                    pending_set[1] = s.strip()
                    pending_set[3] = False
                    continue
                try:
                    vals = _ints(s)
                except ValueError:
                    continue
                pending_set[2].extend(v for v in vals if v != 0)
            elif card == "SEGMENT" and pending_set is not None:
                if pending_set[3]:
                    pending_set[1] = s.strip()
                    pending_set[3] = False
                    continue
                try:
                    vals = _ints(s)
                except ValueError:
                    continue
                if len(vals) >= 3:
                    # Deduplicate repeated 4th vertex on triangular faces (n1 n2 n3 n3)
                    seg = [v for v in vals[:4] if v != 0]
                    if len(seg) == 4 and seg[3] == seg[2]:
                        seg = seg[:3]
                    if len(seg) >= 3:
                        pending_set[2].append(seg)
            elif card is not None and card.startswith("NONSOLID"):
                acc["nonsolid"][card] = acc["nonsolid"].get(card, 0) + 1
            elif card == "INCLUDE_PATH":
                inc = s.strip()
                ingest(os.path.join(base, inc) if not os.path.isabs(inc) else inc)
                card = None

        flush_pending()

    # handle *INCLUDE path lines: re-scan hook (card set above needs the marker)
    ingest(path)

    if acc["nonsolid"]:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(acc["nonsolid"].items()))
        issues.append(make_issue("warning", "NON_SOLID_EXCLUDED",
                                 f"non-solid elements excluded: {summary}",
                                 entity=f"k:{os.path.basename(path)}"))
    if acc["skipped"]:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(acc["skipped"].items()))
        issues.append(make_issue("info", "UNSUPPORTED_DECK_KEYWORD",
                                 f"analysis-deck keywords ignored (mesh+sets scope): {summary}",
                                 entity=f"k:{os.path.basename(path)}"))
    if len(acc["files"]) > 1:
        issues.append(make_issue("info", "INCLUDE_FOLLOWED",
                                 f"*INCLUDE files read: {acc['files']}",
                                 entity=f"k:{os.path.basename(path)}"))

    idx = {nid: k for k, nid in enumerate(acc["order"])}
    nodes = np.array([acc["nodes"][nid] for nid in acc["order"]], dtype=float)
    part_ids = sorted(acc["parts"])
    pid_map = {p: i + 1 for i, p in enumerate(part_ids)} if part_ids else {}
    parts = {pid_map[p]: acc["parts"][p].strip() or f"PART-{p}" for p in part_ids} or {1: "PART-1"}
    by_type_indices: dict[tuple[str, int], list[int]] = {}
    for k, (canon, pid, _conn) in enumerate(acc["elems"]):
        by_type_indices.setdefault((canon, pid_map.get(pid, 1)), []).append(k)

    blocks: list[ElemBlock] = []
    raw_k_to_model_idx: dict[int, int] = {}
    elem_conn0_list: list[tuple[str, list[int]]] = [("", [])] * len(acc["elems"])
    cursor = 0
    for (canon, p), raw_indices in sorted(by_type_indices.items()):
        rows: list[list[int]] = []
        for k in raw_indices:
            conn0 = [idx[v] for v in acc["elems"][k][2]]
            rows.append(conn0)
            raw_k_to_model_idx[k] = cursor
            elem_conn0_list[cursor] = (canon, conn0)
            cursor += 1
        blocks.append(ElemBlock(elem_type=canon, conn=np.asarray(rows, dtype=np.int64), part_id=p))

    node_sets = {k: sorted(idx[v] for v in vs if v in idx) for k, vs in acc["nsets"].items()}
    elem_pos = {eid: raw_k_to_model_idx[k] for k, eid in enumerate(acc["eorder"])}
    elem_sets = {k: sorted(elem_pos[v] for v in vs if v in elem_pos) for k, vs in acc["esets"].items()}

    # *SET_SEGMENT[_TITLE] -> face_sets via exact cyclic + unordered match across all solid blocks.
    face_sets: dict[str, list[tuple[int, int]]] = {}
    if acc["segments"]:
        unordered_lookup: dict[tuple[int, ...], list[tuple[int, int]]] = {}
        for ge, (canon, conn0) in enumerate(elem_conn0_list):
            for lf, f_locals in enumerate(_FACE_DEF_MAP.get(canon, ())):
                ftup = tuple(conn0[m] for m in f_locals)
                unordered_lookup.setdefault(tuple(sorted(ftup)), []).append((ge, lf))

        for sname, seg_list in acc["segments"].items():
            matched: list[tuple[int, int]] = []
            unmatched = 0
            for raw_seg in seg_list:
                seg0 = tuple(idx[v] for v in raw_seg if v in idx)
                cands = unordered_lookup.get(tuple(sorted(seg0)), [])
                found: tuple[int, int] | None = None
                for ge, lf in cands:
                    canon, conn0 = elem_conn0_list[ge]
                    ftup = tuple(conn0[m] for m in _FACE_DEF_MAP[canon][lf])
                    if _cyclic_match(seg0, ftup):
                        found = (ge, lf)
                        break
                if found is None and cands:
                    found = cands[0]
                if found is not None:
                    matched.append(found)
                else:
                    unmatched += 1
            face_sets[sname] = matched
            if unmatched:
                issues.append(make_issue(
                    "warning", "SEGMENT_UNMATCHED",
                    f"{unmatched}/{len(seg_list)} *SET_SEGMENT entries in {sname!r} matched no solid face",
                    entity=f"k:{os.path.basename(path)}"))

    return (MeshModel(nodes=nodes, blocks=blocks, parts=parts,
                      node_sets=node_sets, elem_sets=elem_sets, face_sets=face_sets,
                      unit_system=unit_system), issues)
