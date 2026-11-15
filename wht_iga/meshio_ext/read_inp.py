"""Abaqus ``.inp`` mesh reader (lightweight parser, mesh + sets only).

Supported: ``*NODE``, ``*ELEMENT`` (solid C3D family; ``TYPE`` mapped
through :data:`NORMALIZE`), ``*NSET``/``*ELSET`` (explicit lists and
``GENERATE``), ``*PART``/``*END PART`` + ``*INSTANCE`` (identity
transform; any placement token is issued, never applied), ``*SOLID
SECTION, ELSET=`` for part mapping.

Everything else (``*MATERIAL``, ``*STEP``, ``*BOUNDARY``, loads, ...) is
collected into one aggregated ``UNSUPPORTED_DECK_KEYWORD`` info issue and
otherwise ignored (K25). Shell/beam/spring element types are counted,
issued and excluded (K27).
"""

from __future__ import annotations

import os

import numpy as np

from wht_iga.meshio_ext.model import (
    CANONICAL_COUNTS,
    NON_SOLID_TOKENS,
    NORMALIZE,
    ElemBlock,
    MeshModel,
    apply_order,
    make_issue,
    require_unit_system,
)

# Abaqus surface label S<k> (1-based) -> canonical local_face index matching
# HEX8_FACES / TET4_FACES / WEDGE6_FACES in model.py.
_INP_SURFACE_TO_LF: dict[str, dict[int, int]] = {
    "hex8": {1: 0, 2: 1, 3: 4, 4: 3, 5: 5, 6: 2},
    "hex20": {1: 0, 2: 1, 3: 4, 4: 3, 5: 5, 6: 2},
    "tet4": {1: 0, 2: 1, 3: 2, 4: 3},
    "tet10": {1: 0, 2: 1, 3: 2, 4: 3},
    "wedge6": {1: 0, 2: 1, 3: 2, 4: 4, 5: 3},
}


def _parse_floats(line: str) -> list[float]:
    return [float(t) for t in line.replace(",", " ").split()]


def _parse_id_list(lines: list[str]) -> list[int]:
    out: list[int] = []
    for line in lines:
        for tok in line.replace(",", " ").split():
            out.append(int(tok))
    return out


def _parse_generate(lines: list[str]) -> list[int]:
    out: list[int] = []
    for line in lines:
        parts = [p.strip() for p in line.split(",") if p.strip()]
        if len(parts) >= 2:
            a, b = int(parts[0]), int(parts[1])
            inc = int(parts[2]) if len(parts) > 2 else 1
            out.extend(range(a, b + (1 if inc > 0 else -1), inc))
        else:
            out.extend(int(p) for p in parts)
    return out


def _split_keyword(line: str) -> tuple[str, dict[str, str]]:
    body = line.strip()[1:]
    parts = [p.strip() for p in body.split(",")]
    name = parts[0].upper()
    opts: dict[str, str] = {}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            opts[k.strip().upper()] = v.strip()
        else:
            opts[p.upper()] = "1"
    return name, opts


def read_inp(path: str, unit_system: str) -> tuple[MeshModel, list[dict]]:
    """Read an Abaqus ``.inp`` mesh file. Returns (model, issues)."""
    require_unit_system(unit_system)
    issues: list[dict] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        raw = [ln.rstrip("\n") for ln in f]

    # Strip comment lines ("** ...").
    lines = [ln for ln in raw if not ln.lstrip().startswith("**")]

    node_xyz: dict[int, list[float]] = {}
    node_order: list[int] = []
    raw_elems: list[tuple[str, int, list[int]]] = []  # (TYPE, id, conn)
    elem_id_order: list[int] = []
    nsets: dict[str, list[int]] = {}
    elsets: dict[str, list[int]] = {}
    raw_surfaces: dict[str, list[tuple[str, int]]] = {}
    solid_section_elset: dict[str, str] = {}
    part_blocks: list[tuple[str, int, int]] = []  # (name, first_elem_idx, count)
    skipped_deck: dict[str, int] = {}
    nonsolid: dict[str, int] = {}

    cur_part: str | None = None
    i = 0
    n = len(lines)

    def data_until_next_kw(j: int) -> tuple[list[str], int]:
        out: list[str] = []
        k = j
        while k < n and not lines[k].lstrip().startswith("*"):
            if lines[k].strip():
                out.append(lines[k])
            k += 1
        return out, k

    while i < n:
        line = lines[i].strip()
        if not line or not line.startswith("*"):
            skipped_deck["FREE_TEXT"] = skipped_deck.get("FREE_TEXT", 0) + 1
            i += 1
            continue
        name, opts = _split_keyword(line)
        data, i = data_until_next_kw(i + 1)
        if name == "NODE":
            for dl in data:
                vals = _parse_floats(dl)
                nid = int(vals[0])
                if nid not in node_xyz:
                    node_order.append(nid)
                node_xyz[nid] = [float(vals[1]), float(vals[2]), float(vals[3])]
        elif name == "ELEMENT":
            etype = opts.get("TYPE", "").upper()
            elset_opt = opts.get("ELSET", "").strip()
            canon = NORMALIZE.get(etype)
            noncat = NON_SOLID_TOKENS.get(etype)
            if canon is None:
                if noncat is not None:
                    nonsolid[etype] = nonsolid.get(etype, 0) + len(data)
                    issues.append(make_issue(
                        "warning", "NON_SOLID_EXCLUDED",
                        f"{etype} ({noncat}) elements excluded: 3D solids only",
                        entity=f"inp:{os.path.basename(path)}"))
                else:
                    skipped_deck[f"ELEMENT_TYPE_{etype}"] = len(data)
                continue
            base = len(raw_elems)
            want = CANONICAL_COUNTS[canon]
            # connectivity may wrap across lines: accumulate id tokens first
            toks: list[str] = []
            for dl in data:
                toks.extend(dl.replace(",", " ").split())
            k = 0
            n_elem = 0
            while k < len(toks):
                if k + 1 + want > len(toks):
                    issues.append(make_issue(
                        "warning", "ELEMENT_TRUNCATED",
                        f"{etype} trailing {len(toks) - k} tokens ignored (need {want + 1})",
                        entity=f"inp:{os.path.basename(path)}"))
                    break
                eid = int(toks[k])
                conn = [int(v) for v in toks[k + 1:k + 1 + want]]
                k += 1 + want
                raw_elems.append((canon, eid, apply_order(conn, "inp", canon)))
                elem_id_order.append(eid)
                if elset_opt:
                    elsets.setdefault(elset_opt, []).append(eid)
                n_elem += 1
            if cur_part is not None:
                part_blocks.append((cur_part, base, n_elem))
        elif name in ("NSET", "ELSET"):
            sname = opts.get("NSET", opts.get("ELSET", "")).strip()
            if not sname:
                issues.append(make_issue("warning", "UNNAMED_SET",
                                         f"{name} without a name ignored",
                                         entity=f"inp:{os.path.basename(path)}"))
                continue
            ids = _parse_generate(data) if "GENERATE" in opts else _parse_id_list(data)
            target = nsets if name == "NSET" else elsets
            target[sname] = target.get(sname, []) + ids
        elif name == "PART":
            cur_part = opts.get("NAME", "PART").strip() or "PART"
            in_part_elem_base = len(raw_elems)
        elif name == "END PART":
            cur_part = None
        elif name == "INSTANCE":
            pname = opts.get("PART", "").strip()
            extra = {k: v for k, v in opts.items() if k not in ("NAME", "PART")}
            if extra:
                issues.append(make_issue(
                    "warning", "INSTANCE_TRANSFORM_IGNORED",
                    f"*INSTANCE placement {extra} applied as identity (not supported)",
                    entity=f"inp:{os.path.basename(path)}"))
            if pname:
                issues.append(make_issue(
                    "info", "INSTANCE_REFERENCED",
                    f"*INSTANCE of *PART {pname} merged with identity transform",
                    entity=f"inp:{os.path.basename(path)}"))
        elif name == "SOLID SECTION":
            elset = opts.get("ELSET", "").strip()
            mat = opts.get("MATERIAL", "").strip()
            if elset:
                solid_section_elset[elset.upper()] = mat or elset
        elif name == "SURFACE":
            surf_type = opts.get("TYPE", "ELEMENT").strip().upper()
            sname = opts.get("NAME", "").strip()
            if not sname:
                issues.append(make_issue(
                    "warning", "UNNAMED_SET",
                    "SURFACE without NAME ignored",
                    entity=f"inp:{os.path.basename(path)}"))
                continue
            if surf_type != "ELEMENT":
                skipped_deck[f"SURFACE_TYPE_{surf_type}"] = (
                    skipped_deck.get(f"SURFACE_TYPE_{surf_type}", 0) + len(data)
                )
                continue
            entries = raw_surfaces.setdefault(sname, [])
            for dl in data:
                toks = [t.strip() for t in dl.split(",") if t.strip()]
                if len(toks) < 2:
                    continue
                target_tok, face_tok = toks[0], toks[1].upper()
                if not face_tok.startswith("S") or not face_tok[1:].isdigit():
                    continue
                entries.append((target_tok, int(face_tok[1:])))
        else:
            skipped_deck[name] = skipped_deck.get(name, 0) + len(data)

    if skipped_deck:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(skipped_deck.items()))
        issues.append(make_issue("info", "UNSUPPORTED_DECK_KEYWORD",
                                 f"analysis-deck keywords ignored (mesh+sets scope): {summary}",
                                 entity=f"inp:{os.path.basename(path)}"))

    idx = {nid: k for k, nid in enumerate(node_order)}
    nodes = np.array([node_xyz[nid] for nid in node_order], dtype=float)

    # Part mapping: PART blocks first (deduplicated by part name), else SOLID SECTION ELSET, else single part.
    part_of_elem = [1] * len(raw_elems)
    parts: dict[int, str] = {1: "PART-1"}
    elset_members_upper: dict[str, set[int]] = {e.upper(): set(v) for e, v in elsets.items()}
    if part_blocks:
        parts = {}
        part_name_to_id: dict[str, int] = {}
        for pname, b0, cnt in part_blocks:
            if pname not in part_name_to_id:
                pid = len(part_name_to_id) + 1
                part_name_to_id[pname] = pid
                parts[pid] = pname
            pid = part_name_to_id[pname]
            for k in range(b0, b0 + cnt):
                part_of_elem[k] = pid
    elif solid_section_elset:
        parts = {}
        for pi, (_elset, mat) in enumerate(solid_section_elset.items()):
            parts[pi + 1] = mat
        for k, (_canon, eid, _c) in enumerate(raw_elems):
            for pi, elset in enumerate(solid_section_elset):
                if eid in elset_members_upper.get(elset.upper(), set()):
                    part_of_elem[k] = pi + 1

    by_type_indices: dict[tuple[str, int], list[int]] = {}
    for k, (canon, _eid, _conn) in enumerate(raw_elems):
        by_type_indices.setdefault((canon, part_of_elem[k]), []).append(k)

    blocks: list[ElemBlock] = []
    raw_k_to_model_idx: dict[int, int] = {}
    model_elem_cursor = 0
    for (canon, pid), raw_indices in sorted(by_type_indices.items()):
        rows = [[idx[v] for v in raw_elems[k][2]] for k in raw_indices]
        blocks.append(ElemBlock(elem_type=canon, conn=np.asarray(rows, dtype=np.int64), part_id=pid))
        for k in raw_indices:
            raw_k_to_model_idx[k] = model_elem_cursor
            model_elem_cursor += 1

    node_sets = {k: [idx[v] for v in vs if v in idx] for k, vs in nsets.items()}
    elem_pos = {eid: raw_k_to_model_idx[k] for k, eid in enumerate(elem_id_order)}
    elem_canon = {eid: raw_elems[k][0] for k, eid in enumerate(elem_id_order)}
    elem_sets = {k: sorted(elem_pos[v] for v in vs if v in elem_pos) for k, vs in elsets.items()}

    face_sets: dict[str, list[tuple[int, int]]] = {}
    for sname, entries in raw_surfaces.items():
        pairs: list[tuple[int, int]] = []
        for target_tok, snum in entries:
            if target_tok.lstrip("-").isdigit():
                eids = [int(target_tok)]
            else:
                eids = elsets.get(target_tok)
                if eids is None:
                    eids = list(elset_members_upper.get(target_tok.upper(), ()))
            for eid in eids:
                if eid not in elem_pos:
                    continue
                canon = elem_canon[eid]
                lf_map = _INP_SURFACE_TO_LF.get(canon, {})
                if snum in lf_map:
                    pairs.append((elem_pos[eid], lf_map[snum]))
        face_sets[sname] = pairs

    return (MeshModel(nodes=nodes, blocks=blocks, parts=parts,
                      node_sets=node_sets, elem_sets=elem_sets,
                      face_sets=face_sets,
                      unit_system=unit_system), issues)
