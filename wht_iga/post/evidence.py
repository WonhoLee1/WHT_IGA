# -*- coding: utf-8 -*-
"""Multi-format evidence export (Track H, user-requested).

Writes one solved mesh to every checkable format:
``.vtkhdf`` (in-house), ``.vtu`` / ``.msh`` / ``.fem`` (meshio; ``.fem`` is
Nastran bulk data, directly readable by OptiStruct), ``.stl`` (exterior
triangulation), plus a ``.png`` deformed-shape contour.
``hvascii`` (HyperView ASCII) has no writer in-repo or in meshio and is
explicitly out of scope (no format is ever fabricated).

All functions are Kratos-free; inputs are plain numpy arrays.
"""

from __future__ import annotations

import os as _os

import numpy as np

_VTK_HEXAHEDRON = 12

_HEX_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def exterior_triangles(conn_hex8: np.ndarray) -> np.ndarray:
    """Boundary quad faces of a hex8 mesh as triangles ``(T, 3)``.

    Counts face usage (sorted 4-tuples); keeps faces used once; splits each
    quad ``(a, b, c, d)`` into ``(a, b, c)``, ``(a, c, d)``.
    """
    conn = np.asarray(conn_hex8, dtype=np.int64).reshape(-1, 8)
    counts: dict[tuple[int, int, int, int], int] = {}
    order: list[tuple[int, int, int, int]] = []
    for el in conn:
        for f in _HEX_FACES:
            key = tuple(sorted(int(el[v]) for v in f))
            if key not in counts:
                counts[key] = 0
                order.append(tuple(int(el[v]) for v in f))
            counts[key] += 1
    tris: list[list[int]] = []
    for a, b, c, d in order:
        if counts[tuple(sorted((a, b, c, d)))] == 1:
            tris.append([a, b, c])
            tris.append([a, c, d])
    return np.asarray(tris, dtype=np.int64).reshape(-1, 3)


def write_evidence_set(
    out_dir: str,
    stem: str,
    points: np.ndarray,
    conn_hex8: np.ndarray,
    point_data: dict[str, np.ndarray] | None = None,
    *,
    png_scale: float = 1.0,
    png_field: str = "DISPLACEMENT",
) -> dict[str, str]:
    """Write ``.vtkhdf`` + ``.vtu`` + ``.msh`` + ``.fem`` + ``.stl`` + ``.png``.

    Returns ``{format: path}``. Raises on empty mesh or ragged fields (no mocks).
    """
    import meshio

    from wht_iga.post.proxy_mesh import ProxyMesh
    from wht_iga.post.vtk_writer import write_vtkhdf

    pts = np.asarray(points, dtype=np.float64)
    conn = np.asarray(conn_hex8, dtype=np.int64).reshape(-1, 8)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError(f"points must be (N>=1, 3), got {pts.shape}")
    if conn.shape[0] == 0:
        raise ValueError("conn_hex8 must contain at least one element")
    if int(conn.max()) >= pts.shape[0] or int(conn.min()) < 0:
        raise ValueError("conn_hex8 indices out of range")
    fields = {k: np.asarray(v, dtype=np.float64) for k, v in (point_data or {}).items()}
    for k, v in fields.items():
        if v.shape[0] != pts.shape[0]:
            raise ValueError(f"field {k!r} has {v.shape[0]} rows, need {pts.shape[0]}")

    _os.makedirs(out_dir, exist_ok=True)
    paths: dict[str, str] = {}
    mesh = meshio.Mesh(pts, [("hexahedron", conn)])
    for key, ext, fmt in (("vtu", "vtu", None), ("msh", "msh", "gmsh"), ("fem", "fem", "nastran")):
        path = _os.path.join(out_dir, f"{stem}.{ext}")
        if fmt is None:
            meshio.write(path, mesh)
        else:
            meshio.write(path, mesh, file_format=fmt)
        paths[key] = path

    tris = exterior_triangles(conn)
    tri_mesh = meshio.Mesh(pts, [("triangle", tris)])
    stl_path = _os.path.join(out_dir, f"{stem}.stl")
    meshio.write(stl_path, tri_mesh)
    paths["stl"] = stl_path

    proxy = ProxyMesh(
        points=pts,
        connectivity=conn.ravel(),
        cell_types=np.full(conn.shape[0], _VTK_HEXAHEDRON, dtype=np.uint8),
        point_data={k: np.ascontiguousarray(v) for k, v in fields.items()},
        block_name=stem,
    )
    paths["vtkhdf"] = write_vtkhdf(_os.path.join(out_dir, f"{stem}.vtkhdf"), proxy)
    paths["png"] = _contour_png(_os.path.join(out_dir, f"{stem}.png"), pts, tris,
                                fields.get(png_field), png_scale, stem)
    return paths


def export_bulk_bar_evidence(out_dir: str) -> dict[str, str]:
    """Solve the canonical bulk bar (native FE) and write the full evidence set.

    Requires Kratos (real solve, no mocks). Returns ``{format: path}``.
    """
    from wht_iga.analyses.methods_solve import NativeFeOption, solve_with_discretization_method
    from wht_iga.analyses.presets import bulk_linear

    problem, _, _ = bulk_linear(1.0e5, 0.0, span=4.0, traction=10.0)
    sol = solve_with_discretization_method(
        problem, NativeFeOption(element_type="hex8", n_cells=(2, 1, 1)))
    assert sol.converged
    mesh, u = sol.native["mesh"], sol.native["u_nodal"]
    return write_evidence_set(
        out_dir, "bulk_bar_ux",
        np.asarray(mesh.nodes), np.asarray(mesh.blocks[0].conn),
        {"DISPLACEMENT": np.asarray(u).reshape(-1, 3)})


def _contour_png(path: str, points: np.ndarray, tris: np.ndarray,
                 disp: np.ndarray | None, scale: float, title: str) -> str:
    """Deformed-shape surface contour. True scale iff ``scale == 1`` (labeled)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if disp is None:
        disp = np.zeros_like(points)
    disp = np.asarray(disp, dtype=float).reshape(points.shape[0], -1)
    mag = np.linalg.norm(disp, axis=1) if disp.shape[1] > 1 else np.abs(disp[:, 0])
    deformed = points + float(scale) * (disp if disp.shape[1] == 3
                                        else np.zeros_like(points))
    face_mag = mag[np.asarray(tris, dtype=np.int64)].mean(axis=1)
    fig = plt.figure(figsize=(7, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("jet")
    vmin, vmax = float(face_mag.min()), float(face_mag.max())
    norm = plt.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1e-30)
    coll = Poly3DCollection(deformed[tris], alpha=0.9, edgecolor="k", linewidths=0.15,
                            facecolors=cmap(norm(face_mag)), shade=False)
    ax.add_collection3d(coll)
    import matplotlib.cm as _cm
    fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, label="|u|")
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"{title} |u| (scale={scale:g}, max={mag.max():.3g})", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
