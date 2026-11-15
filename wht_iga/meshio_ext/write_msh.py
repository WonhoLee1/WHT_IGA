"""Gmsh ``.msh`` (MSH 2.2 ASCII) mesh writer via meshio preserving parts and sets.

Emits 3D solid cell blocks with ``$PhysicalNames`` (dim=3), 2D boundary face
blocks for ``face_sets`` with ``$PhysicalNames`` (dim=2), and a standard Gmsh
``$Comments`` metadata section preserving ``node_sets``, ``elem_sets``, and
empty ``face_sets`` keys losslessly for :func:`wht_iga.meshio_ext.read_msh.read_msh`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from wht_iga.meshio_ext.model import HEX8_FACES, TET4_FACES, WEDGE6_FACES, MeshModel
from wht_iga.meshio_ext.write_inp import _merge_extra_sets

_CANON_TO_MESHIO: dict[str, str] = {
    "hex8": "hexahedron",
    "hex20": "hexahedron20",
    "tet4": "tetra",
    "tet10": "tetra10",
    "wedge6": "wedge",
}

_FACE_DEF_MAP = {
    "hex8": HEX8_FACES,
    "hex20": HEX8_FACES,
    "tet4": TET4_FACES,
    "tet10": TET4_FACES,
    "wedge6": WEDGE6_FACES,
}


def write_msh(
    model: MeshModel,
    path: str | Path,
    extra_sets: Any = None,
) -> Path:
    """Write a :class:`MeshModel` to a Gmsh MSH 2.2 ASCII ``.msh`` file.

    Parameters
    ----------
    model:
        Mesh model to serialize.
    path:
        Output ``.msh`` file path.
    extra_sets:
        Optional ``SetRegistry`` or iterable of ``SelectionSet`` objects.

    Returns
    -------
    Path
        Written output file path.
    """
    try:
        import meshio
    except ImportError as exc:  # pragma: no cover - env provides meshio
        raise ImportError("meshio is required for .msh (viz extra)") from exc

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    node_sets, elem_sets = _merge_extra_sets(model.node_sets, model.elem_sets, extra_sets)

    cells: list[tuple[str, np.ndarray]] = []
    phys_tags: list[np.ndarray] = []
    geom_tags: list[np.ndarray] = []
    field_data: dict[str, np.ndarray] = {}

    elem_info: list[tuple[str, np.ndarray]] = []
    used_pids: set[int] = set()
    for block in model.blocks:
        mtype = _CANON_TO_MESHIO.get(block.elem_type)
        if mtype is None:
            raise ValueError(f"Unsupported element type for Gmsh .msh: {block.elem_type!r}")
        pid = int(block.part_id) if block.part_id else 1
        used_pids.add(pid)
        n_b = int(block.conn.shape[0])
        cells.append((mtype, np.asarray(block.conn, dtype=np.int64)))
        phys_tags.append(np.full(n_b, pid, dtype=np.int64))
        geom_tags.append(np.full(n_b, pid, dtype=np.int64))
        for row in block.conn:
            elem_info.append((block.elem_type, row))

    for pid in sorted(used_pids):
        pname = model.parts.get(pid, f"PART-{pid}")
        field_data[pname] = np.array([pid, 3], dtype=np.int64)

    next_surf_pid = max(used_pids, default=0) + 1
    tri_rows: list[list[int]] = []
    tri_tags: list[int] = []
    quad_rows: list[list[int]] = []
    quad_tags: list[int] = []

    for sname, pairs in model.face_sets.items():
        # Ensure field_data name does not collide with a 3D part name.
        fd_key = sname if sname not in field_data else f"SURF_{sname}"
        spid = next_surf_pid
        next_surf_pid += 1
        field_data[fd_key] = np.array([spid, 2], dtype=np.int64)
        for ge, lf in pairs:
            if not (0 <= ge < len(elem_info)):
                continue
            etype, row = elem_info[ge]
            faces_tbl = _FACE_DEF_MAP.get(etype, ())
            if not (0 <= lf < len(faces_tbl)):
                continue
            fnodes = [int(row[m]) for m in faces_tbl[lf]]
            if len(fnodes) == 3:
                tri_rows.append(fnodes)
                tri_tags.append(spid)
            elif len(fnodes) == 4:
                quad_rows.append(fnodes)
                quad_tags.append(spid)

    if tri_rows:
        cells.append(("triangle", np.asarray(tri_rows, dtype=np.int64)))
        arr = np.asarray(tri_tags, dtype=np.int64)
        phys_tags.append(arr)
        geom_tags.append(arr.copy())
    if quad_rows:
        cells.append(("quad", np.asarray(quad_rows, dtype=np.int64)))
        arr = np.asarray(quad_tags, dtype=np.int64)
        phys_tags.append(arr)
        geom_tags.append(arr.copy())

    mesh = meshio.Mesh(
        points=np.asarray(model.nodes, dtype=np.float64),
        cells=cells,
        field_data=field_data,
        cell_data={
            "gmsh:physical": phys_tags,
            "gmsh:geometrical": geom_tags,
        },
    )
    meshio.write(str(out), mesh, file_format="gmsh22", binary=False)

    # Preserve full node_sets / elem_sets / empty face_sets inside Gmsh $Comments.
    meta = {
        "wht_iga_v": 1,
        "node_sets": {str(k): [int(v) for v in vs] for k, vs in node_sets.items()},
        "elem_sets": {str(k): [int(v) for v in vs] for k, vs in elem_sets.items()},
        "face_sets": {
            str(k): [[int(ge), int(lf)] for ge, lf in pairs]
            for k, pairs in model.face_sets.items()
        },
    }
    base_txt = out.read_text(encoding="utf-8")
    out.write_text(
        base_txt.rstrip() + "\n$Comments\n" + json.dumps(meta) + "\n$EndComments\n",
        encoding="utf-8",
    )
    return out


__all__ = ["write_msh"]
