"""Minimal VTU visualization writer: proxy mesh + fields (Track SOLVE).

Phase-1 proxy (v2 §7.12): IN cells as linear hexahedra, values sampled by
evaluating the spline field at mesh nodes. ParaView-openable VTU via meshio.
"""

from __future__ import annotations

import numpy as np


def hex_proxy_grid(lo: np.ndarray, hi: np.ndarray, n_cells: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Nodal grid + hexahedron connectivity for an axis-aligned box mesh."""
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    axes = [np.linspace(lo[d], hi[d], n_cells[d] + 1) for d in range(3)]
    W, V, U = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    points = np.column_stack([U.ravel(), V.ravel(), W.ravel()])
    nx, ny, nz = (n_cells[d] + 1 for d in range(3))
    conn = []
    for k in range(n_cells[2]):
        for j in range(n_cells[1]):
            for i in range(n_cells[0]):
                n0 = i + nx * (j + ny * k)
                conn.append([n0, n0 + 1, n0 + 1 + nx, n0 + nx,
                             n0 + nx * ny, n0 + 1 + nx * ny, n0 + 1 + nx + nx * ny, n0 + nx + nx * ny])
    return points, np.asarray(conn, dtype=np.int64)


def write_vtu(path: str, points: np.ndarray, connectivity: np.ndarray,
              point_data: dict[str, np.ndarray] | None = None) -> str:
    """Write a hexahedron VTU (meshio). Returns path."""
    import meshio

    points = np.asarray(points, dtype=float)
    cells = [("hexahedron", np.asarray(connectivity, dtype=np.int64))]
    mesh = meshio.Mesh(points, cells,
                       point_data={k: np.asarray(v, dtype=float) for k, v in (point_data or {}).items()})
    meshio.write(path, mesh, file_format="vtu")
    return path


def read_vtu(path: str) -> dict:
    """Re-read a VTU (ParaView substitute check when ParaView is unavailable)."""
    import meshio

    mesh = meshio.read(path)
    out = {"points": np.asarray(mesh.points), "cells": [], "point_data": {k: np.asarray(v) for k, v in mesh.point_data.items()}}
    for block in mesh.cells:
        out["cells"].append((block.type, np.asarray(block.data)))
    return out
