"""Gmsh ``.msh`` mesh reader via meshio (mesh + sets only).

meshio (MIT, ``viz`` env extra) parses the file; this module maps
``hexahedron``/``hexahedron20``/``tetra``/``tetra10``/``wedge`` onto the
canonical types, takes Gmsh physical groups as parts, meshio cell sets
as element sets and point sets as node sets. Lower-dimensional cells
(quad/triangle/line/vertex) are counted, issued and excluded (K27, they
are boundary/post cells, not 3D solids).
"""

from __future__ import annotations

import json
import os

import numpy as np

from wht_iga.meshio_ext.model import (
    HEX8_FACES,
    NORMALIZE,
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
    """Return True if ``seq_a`` is a cyclic rotation of ``seq_b``."""
    n = len(seq_a)
    if n != len(seq_b) or n == 0:
        return False
    return any(seq_a == seq_b[k:] + seq_b[:k] for k in range(n))


def _read_wht_comments_metadata(path: str) -> dict:
    """Parse optional WHT_IGA JSON metadata embedded inside a Gmsh ``$Comments`` block."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            in_comments = False
            buf: list[str] = []
            for raw in f:
                s = raw.strip()
                if s == "$Comments":
                    in_comments = True
                    buf.clear()
                    continue
                if s == "$EndComments":
                    if in_comments and buf:
                        try:
                            payload = json.loads("\n".join(buf))
                            if isinstance(payload, dict) and payload.get("wht_iga_v") == 1:
                                return payload
                        except json.JSONDecodeError:
                            pass
                    in_comments = False
                    continue
                if in_comments:
                    buf.append(raw)
    except OSError:
        return {}
    return {}


def read_msh(path: str, unit_system: str) -> tuple[MeshModel, list[dict]]:
    """Read a Gmsh ``.msh`` file via meshio. Returns (model, issues)."""
    require_unit_system(unit_system)
    issues: list[dict] = []
    try:
        import meshio
    except ImportError as exc:  # pragma: no cover - env provides meshio
        raise ImportError("meshio is required for .msh (viz extra)") from exc

    mesh = meshio.read(path)
    base = os.path.basename(path)
    nodes = np.asarray(mesh.points, dtype=float)
    if nodes.shape[1] == 2:
        nodes = np.column_stack([nodes, np.zeros(nodes.shape[0])])

    meta = _read_wht_comments_metadata(path)

    # Physical groups: mesh.field_data maps name -> (id, dim)
    phys_name_by_dim_id: dict[tuple[int, int], str] = {}
    name_of_id: dict[int, str] = {}
    if mesh.field_data:
        for pname, pval in mesh.field_data.items():
            pid = int(pval[0])
            pdim = int(pval[1]) if len(pval) > 1 else 3
            phys_name_by_dim_id[(pdim, pid)] = str(pname)
            name_of_id[pid] = str(pname)

    raw_phys_blocks = mesh.cell_data.get("gmsh:physical", [])

    nonsolid: dict[str, int] = {}
    blocks: list[ElemBlock] = []
    kept_cb_indices: list[tuple[int, str, np.ndarray]] = []
    surf_2d_raw: dict[str, list[tuple[int, ...]]] = {}

    for cb_idx, cb in enumerate(mesh.cells):
        canon = NORMALIZE.get(cb.type.upper())
        tags_arr: np.ndarray | None = None
        if cb_idx < len(raw_phys_blocks):
            tags_arr = np.atleast_1d(np.asarray(raw_phys_blocks[cb_idx], dtype=np.int64).ravel())
        elif cb.type in mesh.cell_data_dict.get("gmsh:physical", {}):
            tags_arr = np.atleast_1d(
                np.asarray(mesh.cell_data_dict["gmsh:physical"][cb.type], dtype=np.int64).ravel()
            )

        if canon is None:
            # Capture named dim=2 physical surface groups as candidate face_sets,
            # while still counting ungrouped / other lower-dim cells as NON_SOLID_EXCLUDED.
            unmapped_count = 0
            if cb.type in ("triangle", "quad", "triangle6", "quad8") and tags_arr is not None:
                n_corner = 3 if cb.type.startswith("triangle") else 4
                for row_idx, row in enumerate(cb.data):
                    ptag = int(tags_arr[row_idx]) if row_idx < len(tags_arr) else 0
                    sname = phys_name_by_dim_id.get((2, ptag)) if ptag > 0 else None
                    if sname is not None:
                        surf_2d_raw.setdefault(sname, []).append(
                            tuple(int(v) for v in row[:n_corner])
                        )
                    else:
                        unmapped_count += 1
            else:
                unmapped_count = len(cb.data)
            if unmapped_count:
                nonsolid[cb.type] = nonsolid.get(cb.type, 0) + unmapped_count
            continue

        conn = np.asarray([[int(v) for v in row] for row in cb.data], dtype=np.int64)
        conn = np.array([apply_order(list(row), "msh", canon) for row in conn], dtype=np.int64)

        if tags_arr is not None and len(tags_arr) == len(conn) and len(conn) > 0:
            # Split contiguous runs by physical tag so global element ordering is preserved.
            run_start = 0
            cur_tag = int(tags_arr[0])
            for r_i in range(1, len(conn) + 1):
                if r_i == len(conn) or int(tags_arr[r_i]) != cur_tag:
                    blocks.append(
                        ElemBlock(
                            elem_type=canon,
                            conn=conn[run_start:r_i],
                            part_id=cur_tag if cur_tag else 0,
                        )
                    )
                    kept_cb_indices.append((cb_idx, cb.type, np.arange(run_start, r_i)))
                    if r_i < len(conn):
                        run_start = r_i
                        cur_tag = int(tags_arr[r_i])
        else:
            phys = int(tags_arr[0]) if (tags_arr is not None and tags_arr.size) else 0
            blocks.append(ElemBlock(elem_type=canon, conn=conn, part_id=phys))
            kept_cb_indices.append((cb_idx, cb.type, np.arange(len(conn))))

    if nonsolid:
        summary = ", ".join(f"{k}({v})" for k, v in sorted(nonsolid.items()))
        issues.append(make_issue("warning", "NON_SOLID_EXCLUDED",
                                 f"lower-dimensional cells excluded: {summary}",
                                 entity=f"msh:{base}"))

    # Parts from physical groups (3D), else single part.
    pids = sorted({b.part_id for b in blocks if b.part_id})
    if pids:
        pmap = {p: i + 1 for i, p in enumerate(pids)}
        parts = {
            pmap[p]: phys_name_by_dim_id.get((3, p), name_of_id.get(p, f"PHYSICAL-{p}"))
            for p in pids
        }
        for b in blocks:
            b.part_id = pmap.get(b.part_id, 1)
    else:
        default_name = phys_name_by_dim_id.get((3, 1), name_of_id.get(1, "PART-1"))
        parts = {1: default_name if len(blocks) <= 1 else "PART-1"}
        for b in blocks:
            b.part_id = 1

    elem_sets: dict[str, list[int]] = {}
    try:
        sets_dict = dict(mesh.cell_sets_dict)
    except Exception:
        sets_dict = {}
    if sets_dict:
        offsets: list[int] = []
        total = 0
        for b in blocks:
            offsets.append(total)
            total += b.conn.shape[0]
        for name, per_type in sets_dict.items():
            collected: list[int] = []
            items = per_type.items() if isinstance(per_type, dict) else [("", per_type)]
            for ctype, idxs in items:
                arr = np.atleast_1d(np.asarray(idxs)).ravel()
                for bi, (_cb_i, mtype, local_slice) in enumerate(kept_cb_indices):
                    if ctype and str(ctype) != mtype:
                        continue
                    local_pos = {int(orig): pos for pos, orig in enumerate(local_slice)}
                    collected.extend(
                        offsets[bi] + local_pos[int(i)]
                        for i in arr
                        if int(i) in local_pos
                    )
            if collected:
                elem_sets[str(name)] = sorted(set(collected))

    node_sets: dict[str, list[int]] = {}
    for k, v in getattr(mesh, "point_sets", {}).items():
        node_sets[str(k)] = sorted(int(i) for i in np.atleast_1d(np.asarray(v).ravel()))

    # Resolve 2D physical surface cells onto (elem_idx, local_face)
    face_sets: dict[str, list[tuple[int, int]]] = {}
    if surf_2d_raw:
        unordered_lookup: dict[tuple[int, ...], list[tuple[int, int]]] = {}
        elem_info: list[tuple[str, np.ndarray]] = []
        ge = 0
        for b in blocks:
            faces_tbl = _FACE_DEF_MAP.get(b.elem_type, ())
            for row in b.conn:
                elem_info.append((b.elem_type, row))
                for lf, f_locals in enumerate(faces_tbl):
                    ftup = tuple(int(row[m]) for m in f_locals)
                    unordered_lookup.setdefault(tuple(sorted(ftup)), []).append((ge, lf))
                ge += 1
        for sname, raw_faces in surf_2d_raw.items():
            pairs: list[tuple[int, int]] = []
            for ftup in raw_faces:
                cands = unordered_lookup.get(tuple(sorted(ftup)), [])
                matched: tuple[int, int] | None = None
                for cand_ge, cand_lf in cands:
                    etype, row = elem_info[cand_ge]
                    canon_face = tuple(int(row[m]) for m in _FACE_DEF_MAP[etype][cand_lf])
                    if _cyclic_match(ftup, canon_face):
                        matched = (cand_ge, cand_lf)
                        break
                if matched is None and cands:
                    matched = cands[0]
                if matched is not None:
                    pairs.append(matched)
            face_sets[sname] = pairs

    # Merge lossless set metadata from embedded $Comments if present
    if meta:
        for sname, nids in meta.get("node_sets", {}).items():
            node_sets[str(sname)] = [int(i) for i in nids]
        for sname, eids in meta.get("elem_sets", {}).items():
            elem_sets[str(sname)] = [int(i) for i in eids]
        for sname, fpairs in meta.get("face_sets", {}).items():
            face_sets[str(sname)] = [(int(e), int(lf)) for e, lf in fpairs]

    return (MeshModel(nodes=nodes, blocks=blocks, parts=parts,
                      node_sets=node_sets, elem_sets=elem_sets,
                      face_sets=face_sets,
                      unit_system=unit_system), issues)
