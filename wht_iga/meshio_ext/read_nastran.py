"""Nastran ``.bdf`` / OptiStruct ``.fem`` mesh reader via pyNastran.

pyNastran (BSD-3, ``mesh`` extra) parses cards; this module maps
CHEXA (8/20), CTETRA (4/10) and CPENTA (6) onto the canonical types and
applies the global frame (GRID + CORD2R, K26). Part mapping: ``PSOLID``
``pid`` -> part; OptiStruct ``$$HMNAME COMP`` comments
(``$$HMNAME COMP \\\\<pid>\\\\"<name>"``) supply component names when present.

Shell (CQUAD4/CTRIA3/...), beam (CBAR/CBEAM/CROD), spring (CELAS*/CBUSH)
and mass (CONM2) elements are counted, issued and excluded (K27).
``unit_system`` is mandatory (K26); deck cards (loads, SPC, MAT, ...)
are aggregated into one info issue (K25).
"""

from __future__ import annotations

import io
import os
import re

import numpy as np

from wht_iga.meshio_ext.model import (
    HEX8_FACES,
    TET4_FACES,
    WEDGE6_FACES,
    NON_SOLID_TOKENS,
    ElemBlock,
    MeshModel,
    apply_order,
    make_issue,
    require_unit_system,
)

_HMNAME_RE = re.compile(r"\$\$HMNAME\s+COMP\s+\\(\d+)\\?\"([^\"]*)\"?", re.IGNORECASE)
_HMSET_RE = re.compile(
    r"\$\$HMNAME\s+SETS\s+\\(\d+)\\(\d+)\\?\"([^\"]*)\"?", re.IGNORECASE
)

_SOLID_TYPES = {"CHEXA", "CTETRA", "CPENTA"}
_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}


def _split_nastran_fields(line: str) -> list[str]:
    """Split a bulk-data card line (free-field CSV or 8-char fixed field)."""
    s = line.split("$", 1)[0].rstrip("\r\n")
    if "," in s:
        return [tok.strip() for tok in s.split(",")]
    return [s[k : k + 8].strip() for k in range(0, len(s), 8)]


def _expand_nastran_int_tokens(tokens: list[str]) -> list[int]:
    """Expand Nastran integer list tokens including ``A THRU B`` ranges."""
    out: list[int] = []
    i = 0
    while i < len(tokens):
        tok = tokens[i].strip()
        if not tok:
            i += 1
            continue
        if (
            i + 2 < len(tokens)
            and tokens[i + 1].strip().upper() == "THRU"
            and tok.lstrip("-").isdigit()
            and tokens[i + 2].strip().lstrip("-").isdigit()
        ):
            a, b = int(tok), int(tokens[i + 2].strip())
            step = 1 if b >= a else -1
            out.extend(range(a, b + step, step))
            i += 3
        elif tok.lstrip("-").isdigit():
            out.append(int(tok))
            i += 1
        else:
            i += 1
    return out


def _extract_optistruct_sets_and_filter(
    raw_text: str,
) -> tuple[
    str,
    dict[int, tuple[int, str]],
    list[tuple[int, str, list[int]]],
    list[tuple[int, list[tuple[int, list[int]]]]],
]:
    """Extract OptiStruct ``$$HMNAME SETS``, ``SET``, and ``SURF`` cards; strip them for pyNastran."""
    set_meta: dict[int, tuple[int, str]] = {}
    for m in _HMSET_RE.finditer(raw_text):
        set_meta[int(m.group(1))] = (int(m.group(2)), m.group(3).strip())

    opti_sets: list[tuple[int, str, list[int]]] = []
    opti_surfs: list[tuple[int, list[tuple[int, list[int]]]]] = []
    filtered_lines: list[str] = []

    raw_lines = raw_text.splitlines()
    n = len(raw_lines)
    i = 0
    while i < n:
        ln = raw_lines[i]
        stripped = ln.strip()
        if not stripped or stripped.startswith("$"):
            filtered_lines.append(ln)
            i += 1
            continue

        # Collect continuation lines belonging to current card.
        card_lines = [ln]
        j = i + 1
        while j < n:
            nxt = raw_lines[j]
            nxt_s = nxt.strip()
            if not nxt_s or nxt_s.startswith("$"):
                break
            fields = _split_nastran_fields(nxt)
            head = fields[0].upper() if fields else ""
            if nxt[0] in (" ", "\t", "+", "*") or head == "":
                card_lines.append(nxt)
                j += 1
            else:
                break

        first_upper = stripped.upper()
        # Handle "SET <sid> = ..." case-control / OptiStruct syntax or CSV/fixed "SET,sid,..."
        if first_upper.startswith("SET ") and "=" in first_upper:
            lhs, rhs = stripped.split("=", 1)
            sid_tok = lhs.strip().split()[1] if len(lhs.strip().split()) > 1 else "0"
            sid = int(sid_tok) if sid_tok.lstrip("-").isdigit() else 0
            toks = rhs.replace(",", " ").split()
            for cl in card_lines[1:]:
                toks.extend(cl.lstrip("+*, \t").replace(",", " ").split())
            opti_sets.append((sid, "UNKNOWN", _expand_nastran_int_tokens(toks)))
            i = j
            continue

        fields0 = _split_nastran_fields(ln)
        card_name = fields0[0].upper() if fields0 else ""
        if card_name == "SET":
            sid = int(fields0[1]) if len(fields0) > 1 and fields0[1].isdigit() else 0
            stype = fields0[2].upper() if len(fields0) > 2 else "GRID"
            toks = [t for t in fields0[4:] if t] if len(fields0) > 4 else []
            for cl in card_lines[1:]:
                cf = _split_nastran_fields(cl)
                toks.extend(t for t in cf[1:] if t)
            opti_sets.append((sid, stype, _expand_nastran_int_tokens(toks)))
            i = j
            continue

        if card_name == "SURF":
            sid = int(fields0[1]) if len(fields0) > 1 and fields0[1].isdigit() else 0
            entries: list[tuple[int, list[int]]] = []
            for cl in card_lines[1:]:
                cf = _split_nastran_fields(cl)
                nums = [int(t) for t in cf[1:] if t and t.lstrip("-").isdigit()]
                if len(nums) >= 2:
                    entries.append((nums[0], nums[1:]))
            opti_surfs.append((sid, entries))
            i = j
            continue

        filtered_lines.extend(card_lines)
        i = j

    return "\n".join(filtered_lines) + "\n", set_meta, opti_sets, opti_surfs


def _node_count(elem: object) -> int:
    ids = list(elem.node_ids)
    return len(ids)


def read_nastran(path: str, unit_system: str) -> tuple[MeshModel, list[dict]]:
    """Read a Nastran ``.bdf`` / OptiStruct ``.fem`` mesh. Returns (model, issues)."""
    require_unit_system(unit_system)
    issues: list[dict] = []
    try:
        from pyNastran.bdf.bdf import BDF
    except ImportError as exc:  # pragma: no cover - env provides pyNastran
        raise ImportError("pyNastran is required for .bdf/.fem (mesh extra)") from exc

    with open(path, encoding="utf-8", errors="replace") as f:
        raw_text = f.read()

    filtered_bdf, set_meta, opti_sets, opti_surfs = _extract_optistruct_sets_and_filter(raw_text)

    model = BDF(debug=False)
    # pure bulk-data mesh decks (no exec/case control) parse with punch=True
    model.read_bdf(io.StringIO(filtered_bdf), punch=True)

    base = os.path.basename(path)
    nodes: dict[int, list[float]] = {}
    order: list[int] = []
    for nid in sorted(model.nodes):
        node = model.nodes[nid]
        try:
            pos = node.get_position()
        except Exception:
            pos = node.xyz
        pos = np.asarray(pos, dtype=float).ravel()
        if nid not in nodes:
            order.append(nid)
        nodes[nid] = [float(pos[0]), float(pos[1]), float(pos[2])]

    comp_names: dict[int, str] = {}
    for m in _HMNAME_RE.finditer(raw_text):
        comp_names[int(m.group(1))] = m.group(2).strip()

    elems: list[tuple[int, str, int, list[int]]] = []
    nonsolid: dict[str, int] = {}
    for eid in sorted(model.elements):
        elem = model.elements[eid]
        etype = elem.type.upper()
        if etype not in _SOLID_TYPES:
            nonsolid[etype] = nonsolid.get(etype, 0) + 1
            continue
        nids = list(elem.node_ids)
        count = len(nids)
        if etype == "CHEXA":
            canon = "hex8" if count == 8 else "hex20" if count == 20 else ""
        elif etype == "CTETRA":
            canon = "tet4" if count == 4 else "tet10" if count == 10 else ""
        else:
            canon = "wedge6" if count == 6 else ""
        if not canon:
            nonsolid[f"{etype}N{count}"] = nonsolid.get(f"{etype}N{count}", 0) + 1
            issues.append(make_issue(
                "warning", "NON_SOLID_EXCLUDED",
                f"{etype} eid={eid} with {count} nodes excluded",
                entity=f"bdf:{base}"))
            continue
        pid = int(getattr(elem, "Pid", lambda: 1)())
        elems.append((int(eid), canon, pid, apply_order([int(v) for v in nids], "nastran", canon)))

    if nonsolid:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(nonsolid.items()))
        cat = next((NON_SOLID_TOKENS.get(t.split('N')[0], "") for t in nonsolid), "")
        issues.append(make_issue("warning", "NON_SOLID_EXCLUDED",
                                 f"non-solid elements excluded ({cat or 'mixed'}): {summary}",
                                 entity=f"bdf:{base}"))

    pid_of: dict[int, int] = {}
    for pid, prop in model.properties.items():
        if prop.type.upper() == "PSOLID":
            pid_of[int(pid)] = int(pid)
    part_names: dict[int, str] = {}
    for _eid, _canon, pid, _c in elems:
        if pid not in part_names:
            part_names[pid] = comp_names.get(pid, f"PSOLID-{pid}")
    prop_pids = sorted(part_names)
    pmap = {p: i + 1 for i, p in enumerate(prop_pids)} if prop_pids else {}
    parts = {pmap[p]: part_names[p] for p in prop_pids} or {1: "PART-1"}
    if pid_of.keys() - set(part_names):
        issues.append(make_issue(
            "info", "UNREFERENCED_PSOLID",
            f"PSOLID pids with no solid elements: {sorted(pid_of.keys() - set(part_names))}",
            entity=f"bdf:{base}"))

    idx = {nid: k for k, nid in enumerate(order)}
    nodes_arr = np.array([nodes[nid] for nid in order], dtype=float)
    by_type_indices: dict[tuple[str, int], list[int]] = {}
    for k, (_eid, canon, pid, _conn) in enumerate(elems):
        by_type_indices.setdefault((canon, pmap.get(pid, 1)), []).append(k)

    blocks: list[ElemBlock] = []
    elem_pos: dict[int, int] = {}
    elem_canon: dict[int, str] = {}
    elem_conn0: dict[int, list[int]] = {}
    cursor = 0
    for (canon, p), raw_indices in sorted(by_type_indices.items()):
        rows: list[list[int]] = []
        for k in raw_indices:
            eid, _, _pid, conn = elems[k]
            conn0 = [idx[v] for v in conn]
            rows.append(conn0)
            elem_pos[eid] = cursor
            elem_canon[eid] = canon
            elem_conn0[eid] = conn0
            cursor += 1
        blocks.append(ElemBlock(elem_type=canon, conn=np.asarray(rows, dtype=np.int64), part_id=p))

    node_sets: dict[str, list[int]] = {}
    elem_sets: dict[str, list[int]] = {}
    face_sets: dict[str, list[tuple[int, int]]] = {}

    # 1) Nastran SET1 / SET3 parsed by pyNastran
    for sid, s_obj in sorted(model.sets.items()):
        stype = s_obj.type.upper()
        desc = (getattr(s_obj, "desc", "") or "").upper()
        hm_kind, hm_name = set_meta.get(int(sid), (0, ""))
        raw_ids = [int(v) for v in getattr(s_obj, "ids", [])]
        if stype == "SET1" or (stype == "SET3" and desc in ("GRID", "NODE", "NODES")) or hm_kind == 1:
            sname = hm_name or f"SET_{sid}"
            node_sets[sname] = [idx[v] for v in raw_ids if v in idx]
        elif stype == "SET3" and desc in ("ELEM", "ELEMENT", "ELEMENTS") or hm_kind == 2:
            sname = hm_name or f"ELSET_{sid}"
            elem_sets[sname] = sorted(elem_pos[v] for v in raw_ids if v in elem_pos)

    # 2) OptiStruct SET cards pre-parsed above
    for sid, stype, raw_ids in opti_sets:
        hm_kind, hm_name = set_meta.get(sid, (0, ""))
        if stype in ("ELEM", "ELEMENT", "ELEMENTS") or hm_kind == 2:
            sname = hm_name or f"ELSET_{sid}"
            elem_sets[sname] = sorted(elem_pos[v] for v in raw_ids if v in elem_pos)
        else:
            sname = hm_name or f"SET_{sid}"
            node_sets[sname] = [idx[v] for v in raw_ids if v in idx]

    # 3) OptiStruct SURF cards pre-parsed above (matched by face nodes or 1-based face index)
    for sid, entries in opti_surfs:
        _hm_kind, hm_name = set_meta.get(sid, (3, ""))
        sname = hm_name or f"SURF_{sid}"
        pairs: list[tuple[int, int]] = []
        for eid, tail in entries:
            if eid not in elem_pos:
                continue
            ge = elem_pos[eid]
            canon = elem_canon[eid]
            faces_tbl = _FACE_DEF_MAP.get(canon, ())
            if len(tail) == 1:
                lf = tail[0] - 1
                if 0 <= lf < len(faces_tbl):
                    pairs.append((ge, lf))
            else:
                conn0 = elem_conn0[eid]
                target_nids = {idx[v] for v in tail if v in idx}
                for lf, f_locals in enumerate(faces_tbl):
                    face_nids = {conn0[m] for m in f_locals}
                    if target_nids == face_nids:
                        pairs.append((ge, lf))
                        break
        face_sets[sname] = pairs

    # 4) Empty sets declared only via $$HMNAME SETS metadata
    for sid, (hm_kind, hm_name) in sorted(set_meta.items()):
        if not hm_name:
            continue
        if hm_kind == 1 and hm_name not in node_sets:
            node_sets[hm_name] = []
        elif hm_kind == 2 and hm_name not in elem_sets:
            elem_sets[hm_name] = []
        elif hm_kind == 3 and hm_name not in face_sets:
            face_sets[hm_name] = []

    deck_cards = sorted({c for c in model.card_count if c not in
                         {"GRID", "CHEXA", "CTETRA", "CPENTA", "PSOLID",
                          "CORD2R", "CORD2C", "CORD2S", "CORD1R", "CORD1C",
                          "CORD1S", "GRDSET", "PARAM", "MDLPRM",
                          "SET1", "SET3", "BSURF"}})
    if deck_cards:
        issues.append(make_issue(
            "info", "UNSUPPORTED_DECK_KEYWORD",
            f"analysis-deck cards present but ignored (mesh+sets scope): {deck_cards}",
            entity=f"bdf:{base}"))
    return (MeshModel(nodes=nodes_arr, blocks=blocks, parts=parts,
                      node_sets=node_sets, elem_sets=elem_sets,
                      face_sets=face_sets,
                      unit_system=unit_system), issues)
