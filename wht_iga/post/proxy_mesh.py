"""Visualization proxy mesh generation and B-spline field evaluation.

Generates ParaView-compatible visualization proxy meshes from immersed IGA
background grid cells:
- ``IN`` cells -> subdivided linear hexahedra (``VTK_HEXAHEDRON = 12``)
- ``CUT`` cells -> 6-tetrahedra Kuhn decomposition (``VTK_TETRA = 10``),
  optionally filtered by an inside predicate, clipped by a signed level-set
  function, or built from pre-supplied cut tetrahedra.
- Optional tensor-product B-spline evaluator so proxy point fields match
  exact control-point re-evaluation to machine precision.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from wht_iga.grid.bspline import ders_basis_funs, find_span

VTK_TETRA: int = 10
VTK_HEXAHEDRON: int = 12

CELL_OUT: int = 0
CELL_IN: int = 1
CELL_CUT: int = 2

# 6-tetrahedra Kuhn decomposition of a VTK_HEXAHEDRON [v0..v7], all sharing diagonal 0-6
# and all having positive Jacobian determinant for standard right-handed hex ordering.
_KUHN_TETS: tuple[tuple[int, int, int, int], ...] = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)

_NODES_PER_TYPE: dict[int, int] = {
    VTK_TETRA: 4,
    VTK_HEXAHEDRON: 8,
}


@dataclass
class ProxyMesh:
    """Container for mixed-topology visualization proxy meshes.

    Parameters
    ----------
    points:
        Nodal coordinates of shape ``(n_points, 3)``.
    connectivity:
        Flat 1D integer connectivity array of vertex indices, or 2D array
        (automatically flattened in ``__post_init__``).
    cell_types:
        1D ``uint8`` array of VTK cell types (e.g., ``VTK_HEXAHEDRON=12``,
        ``VTK_TETRA=10``) of length ``n_cells``.
    point_data:
        Mapping from field name to array of shape ``(n_points, ...)``.
    cell_data:
        Mapping from field name to array of shape ``(n_cells, ...)``.
    block_name:
        Logical block identifier (e.g., ``"iga_proxy"``, ``"in_cells"``).
    offsets:
        Optional 1D ``int64`` offsets array of length ``n_cells + 1``. Computed
        automatically from ``cell_types`` when omitted.
    """

    points: np.ndarray
    connectivity: np.ndarray
    cell_types: np.ndarray
    point_data: dict[str, np.ndarray] = field(default_factory=dict)
    cell_data: dict[str, np.ndarray] = field(default_factory=dict)
    block_name: str = "iga_proxy"
    offsets: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=np.float64)
        if self.points.ndim != 2 or self.points.shape[1] != 3:
            raise ValueError(f"points must have shape (N, 3), got {self.points.shape}")
        self.cell_types = np.asarray(self.cell_types, dtype=np.uint8)
        conn = np.asarray(self.connectivity, dtype=np.int64)
        if conn.ndim == 2:
            conn = conn.ravel()
        self.connectivity = conn

        n_cells = len(self.cell_types)
        if self.offsets is None:
            counts = np.empty(n_cells, dtype=np.int64)
            for idx, ctype in enumerate(self.cell_types):
                n_nodes = _NODES_PER_TYPE.get(int(ctype))
                if n_nodes is None:
                    raise ValueError(f"Unsupported VTK cell type: {int(ctype)}")
                counts[idx] = n_nodes
            offsets = np.zeros(n_cells + 1, dtype=np.int64)
            np.cumsum(counts, out=offsets[1:])
            self.offsets = offsets
        else:
            self.offsets = np.asarray(self.offsets, dtype=np.int64)
            if len(self.offsets) != n_cells + 1:
                raise ValueError(
                    f"offsets length {len(self.offsets)} must equal n_cells + 1 ({n_cells + 1})"
                )

        if len(self.connectivity) != int(self.offsets[-1]):
            raise ValueError(
                f"connectivity length {len(self.connectivity)} != offsets[-1] ({int(self.offsets[-1])})"
            )

        self.point_data = {k: np.asarray(v, dtype=np.float64) for k, v in self.point_data.items()}
        self.cell_data = {k: np.asarray(v, dtype=np.float64) for k, v in self.cell_data.items()}

    @property
    def n_points(self) -> int:
        """Number of mesh vertices."""
        return int(self.points.shape[0])

    @property
    def n_cells(self) -> int:
        """Number of mesh cells."""
        return int(self.cell_types.shape[0])

    def cell_nodes(self, cell_idx: int) -> np.ndarray:
        """Return vertex indices for cell ``cell_idx``."""
        assert self.offsets is not None
        return self.connectivity[self.offsets[cell_idx] : self.offsets[cell_idx + 1]]


class BSplineFieldEvaluator:
    """Exact 3D tensor-product B-spline field evaluator.

    Parameters
    ----------
    knots:
        Tuple of 3 clamped knot vectors ``(knots_u, knots_v, knots_w)`` in
        physical coordinates (matching ``wht_iga.grid.knots``).
    degrees:
        Tuple of polynomial degrees ``(p_u, p_v, p_w)``.
    """

    def __init__(
        self,
        knots: tuple[np.ndarray, np.ndarray, np.ndarray],
        degrees: tuple[int, int, int],
    ) -> None:
        self.knots = tuple(np.asarray(k, dtype=np.float64) for k in knots)
        self.degrees = tuple(int(p) for p in degrees)
        if len(self.knots) != 3 or len(self.degrees) != 3:
            raise ValueError("knots and degrees must each have 3 components (u, v, w)")
        self.n_cp = tuple(len(self.knots[d]) - self.degrees[d] - 1 for d in range(3))
        for d in range(3):
            if self.n_cp[d] <= self.degrees[d]:
                raise ValueError(f"Invalid knot vector along axis {d}: n_cp={self.n_cp[d]}, p={self.degrees[d]}")

    def _reshape_cp(self, cp_values: np.ndarray) -> tuple[np.ndarray, tuple[int, ...]]:
        arr = np.asarray(cp_values, dtype=np.float64)
        nx, ny, nz = self.n_cp
        n_total = nx * ny * nz
        if arr.shape[:3] == (nx, ny, nz):
            tail_shape = arr.shape[3:]
            return arr, tail_shape
        if arr.shape[0] == n_total:
            tail_shape = arr.shape[1:]
            return arr.reshape((nx, ny, nz) + tail_shape), tail_shape
        raise ValueError(
            f"cp_values shape {arr.shape} incompatible with control grid {self.n_cp} ({n_total} CPs)"
        )

    def evaluate(self, points: np.ndarray, cp_values: np.ndarray) -> np.ndarray:
        """Evaluate a scalar or vector/tensor B-spline field at ``points``.

        Parameters
        ----------
        points:
            Physical coordinates of shape ``(n_pts, 3)``.
        cp_values:
            Control-point values of shape ``(nx, ny, nz, ...)`` or ``(nx*ny*nz, ...)``.

        Returns
        -------
        np.ndarray
            Evaluated field values of shape ``(n_pts, ...)``.
        """
        pts = np.asarray(points, dtype=np.float64)
        if pts.ndim == 1:
            pts = pts[None, :]
        cp_grid, tail_shape = self._reshape_cp(cp_values)
        n_pts = pts.shape[0]
        out = np.zeros((n_pts,) + tail_shape, dtype=np.float64)

        ku, kv, kw = self.knots
        pu, pv, pw = self.degrees
        nu, nv, nw = self.n_cp

        for idx in range(n_pts):
            x, y, z = float(pts[idx, 0]), float(pts[idx, 1]), float(pts[idx, 2])
            su = find_span(nu, pu, x, ku)
            sv = find_span(nv, pv, y, kv)
            sw = find_span(nw, pw, z, kw)
            Nu = ders_basis_funs(su, x, pu, 0, ku)[0]
            Nv = ders_basis_funs(sv, y, pv, 0, kv)[0]
            Nw = ders_basis_funs(sw, z, pw, 0, kw)[0]

            local_cp = cp_grid[su - pu : su + 1, sv - pv : sv + 1, sw - pw : sw + 1]
            # Tensor contraction over local (pu+1, pv+1, pw+1) support
            weights = Nu[:, None, None] * Nv[None, :, None] * Nw[None, None, :]
            if tail_shape:
                out[idx] = np.tensordot(weights, local_cp, axes=([0, 1, 2], [0, 1, 2]))
            else:
                out[idx] = float(np.sum(weights * local_cp))
        return out

    def evaluate_gradient(self, points: np.ndarray, cp_values: np.ndarray) -> np.ndarray:
        """Evaluate spatial gradients ``d(field)/d(x,y,z)`` at ``points``.

        Parameters
        ----------
        points:
            Physical coordinates of shape ``(n_pts, 3)``.
        cp_values:
            Control-point values of shape ``(nx, ny, nz, ...)`` or ``(nx*ny*nz, ...)``.

        Returns
        -------
        np.ndarray
            Gradient array of shape ``(n_pts, ..., 3)`` where the last axis
            holds ``(d/dx, d/dy, d/dz)``.
        """
        pts = np.asarray(points, dtype=np.float64)
        if pts.ndim == 1:
            pts = pts[None, :]
        cp_grid, tail_shape = self._reshape_cp(cp_values)
        n_pts = pts.shape[0]
        out = np.zeros((n_pts,) + tail_shape + (3,), dtype=np.float64)

        ku, kv, kw = self.knots
        pu, pv, pw = self.degrees
        nu, nv, nw = self.n_cp

        for idx in range(n_pts):
            x, y, z = float(pts[idx, 0]), float(pts[idx, 1]), float(pts[idx, 2])
            su = find_span(nu, pu, x, ku)
            sv = find_span(nv, pv, y, kv)
            sw = find_span(nw, pw, z, kw)
            du = ders_basis_funs(su, x, pu, 1, ku)
            dv = ders_basis_funs(sv, y, pv, 1, kv)
            dw = ders_basis_funs(sw, z, pw, 1, kw)

            local_cp = cp_grid[su - pu : su + 1, sv - pv : sv + 1, sw - pw : sw + 1]
            wx = du[1][:, None, None] * dv[0][None, :, None] * dw[0][None, None, :]
            wy = du[0][:, None, None] * dv[1][None, :, None] * dw[0][None, None, :]
            wz = du[0][:, None, None] * dv[0][None, :, None] * dw[1][None, None, :]

            if tail_shape:
                out[idx, ..., 0] = np.tensordot(wx, local_cp, axes=([0, 1, 2], [0, 1, 2]))
                out[idx, ..., 1] = np.tensordot(wy, local_cp, axes=([0, 1, 2], [0, 1, 2]))
                out[idx, ..., 2] = np.tensordot(wz, local_cp, axes=([0, 1, 2], [0, 1, 2]))
            else:
                out[idx, 0] = float(np.sum(wx * local_cp))
                out[idx, 1] = float(np.sum(wy * local_cp))
                out[idx, 2] = float(np.sum(wz * local_cp))
        return out


def _normalize_sub(n_sub: int | tuple[int, int, int]) -> tuple[int, int, int]:
    if isinstance(n_sub, int):
        if n_sub < 1:
            raise ValueError(f"n_sub must be >= 1, got {n_sub}")
        return (n_sub, n_sub, n_sub)
    if len(n_sub) != 3 or any(v < 1 for v in n_sub):
        raise ValueError(f"n_sub must be 3 positive ints, got {n_sub}")
    return (int(n_sub[0]), int(n_sub[1]), int(n_sub[2]))


def _normalize_bounds(cell_bounds: np.ndarray) -> np.ndarray:
    b = np.asarray(cell_bounds, dtype=np.float64)
    if b.size == 0:
        return b.reshape((0, 2, 3))
    if b.ndim == 2 and b.shape[1] == 6:
        b = b.reshape((-1, 2, 3))
    if b.ndim != 3 or b.shape[1:] != (2, 3):
        raise ValueError(f"cell_bounds must have shape (N, 2, 3) or (N, 6), got {b.shape}")
    return b


def decompose_hex_to_6_tets(hex_nodes: np.ndarray) -> np.ndarray:
    """Decompose hexahedron vertex index tuples ``(M, 8)`` into ``(6*M, 4)`` tets."""
    h = np.asarray(hex_nodes, dtype=np.int64)
    if h.ndim == 1:
        h = h[None, :]
    if h.shape[1] != 8:
        raise ValueError(f"hex_nodes must have shape (M, 8), got {h.shape}")
    tets = np.empty((h.shape[0] * 6, 4), dtype=np.int64)
    for k, (a, b, c, d) in enumerate(_KUHN_TETS):
        tets[k::6, 0] = h[:, a]
        tets[k::6, 1] = h[:, b]
        tets[k::6, 2] = h[:, c]
        tets[k::6, 3] = h[:, d]
    return tets


def _subdivide_box(lo: np.ndarray, hi: np.ndarray, sub: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(points, hex_conn)`` for a single axis-aligned box."""
    nsx, nsy, nsz = sub
    xs = np.linspace(lo[0], hi[0], nsx + 1)
    ys = np.linspace(lo[1], hi[1], nsy + 1)
    zs = np.linspace(lo[2], hi[2], nsz + 1)
    Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    nx, ny = nsx + 1, nsy + 1
    conn = np.empty((nsx * nsy * nsz, 8), dtype=np.int64)
    idx = 0
    for k in range(nsz):
        for j in range(nsy):
            for i in range(nsx):
                n0 = i + nx * (j + ny * k)
                n1 = n0 + 1
                n2 = n0 + 1 + nx
                n3 = n0 + nx
                n4 = n0 + nx * ny
                n5 = n4 + 1
                n6 = n4 + 1 + nx
                n7 = n4 + nx
                conn[idx] = (n0, n1, n2, n3, n4, n5, n6, n7)
                idx += 1
    return pts, conn


def _clip_tet_negative_halfspace(
    pts_list: list[np.ndarray],
    tet_coords: np.ndarray,
    phi_vals: np.ndarray,
    tol: float = 1e-12,
) -> list[tuple[int, int, int, int]]:
    """Clip a single tetrahedron to ``phi <= 0`` via exact linear edge interpolation.

    Returns a list of 4-tuples indexing into ``pts_list``.
    """
    inside = phi_vals <= tol
    n_in = int(np.sum(inside))
    if n_in == 0:
        return []

    def _add_pt(p: np.ndarray) -> int:
        idx = len(pts_list)
        pts_list.append(p)
        return idx

    if n_in == 4:
        base = len(pts_list)
        for v in tet_coords:
            pts_list.append(v)
        return [(base, base + 1, base + 2, base + 3)]

    def _interp(ia: int, ib: int) -> np.ndarray:
        pa, pb = tet_coords[ia], tet_coords[ib]
        fa, fb = float(phi_vals[ia]), float(phi_vals[ib])
        denom = fa - fb
        t = 0.5 if abs(denom) < 1e-15 else fa / denom
        t = float(np.clip(t, 0.0, 1.0))
        return pa + t * (pb - pa)

    in_idx = [i for i in range(4) if inside[i]]
    out_idx = [i for i in range(4) if not inside[i]]

    if n_in == 1:
        i0 = in_idx[0]
        j0, j1, j2 = out_idx
        v0 = _add_pt(tet_coords[i0])
        q0 = _add_pt(_interp(i0, j0))
        q1 = _add_pt(_interp(i0, j1))
        q2 = _add_pt(_interp(i0, j2))
        return [(v0, q0, q1, q2)]

    if n_in == 2:
        i0, i1 = in_idx
        j0, j1 = out_idx
        v0 = _add_pt(tet_coords[i0])
        v1 = _add_pt(tet_coords[i1])
        q00 = _add_pt(_interp(i0, j0))
        q01 = _add_pt(_interp(i0, j1))
        q10 = _add_pt(_interp(i1, j0))
        q11 = _add_pt(_interp(i1, j1))
        # Triangular prism (v0, q00, q01) -> (v1, q10, q11) split into 3 tets
        return [
            (v0, q00, q01, v1),
            (q00, q01, v1, q10),
            (q01, v1, q10, q11),
        ]

    # n_in == 3: 1 vertex outside -> truncated pyramid (prism) split into 3 tets
    j0 = out_idx[0]
    i0, i1, i2 = in_idx
    v0 = _add_pt(tet_coords[i0])
    v1 = _add_pt(tet_coords[i1])
    v2 = _add_pt(tet_coords[i2])
    q0 = _add_pt(_interp(i0, j0))
    q1 = _add_pt(_interp(i1, j0))
    q2 = _add_pt(_interp(i2, j0))
    return [
        (v0, v1, v2, q0),
        (v1, v2, q0, q1),
        (v2, q0, q1, q2),
    ]


def generate_hex_proxy_for_cells(
    cell_bounds: np.ndarray,
    n_sub: int | tuple[int, int, int] = 2,
    block_name: str = "in_cells",
    deduplicate_tolerance: float | None = 1e-12,
) -> ProxyMesh:
    """Generate a subdivided ``VTK_HEXAHEDRON`` proxy mesh for ``IN`` cells.

    Parameters
    ----------
    cell_bounds:
        Array of shape ``(n_cells, 2, 3)`` or ``(n_cells, 6)`` giving
        ``[lo, hi]`` per cell.
    n_sub:
        Subdivision count per cell axis (int or 3-tuple).
    block_name:
        Block label stored on the returned :class:`ProxyMesh`.
    deduplicate_tolerance:
        Optional coordinate rounding tolerance to merge shared vertices across
        adjacent cells. Pass ``None`` to skip deduplication.
    """
    bounds = _normalize_bounds(cell_bounds)
    sub = _normalize_sub(n_sub)
    if bounds.shape[0] == 0:
        return ProxyMesh(
            points=np.zeros((0, 3), dtype=np.float64),
            connectivity=np.zeros(0, dtype=np.int64),
            cell_types=np.zeros(0, dtype=np.uint8),
            block_name=block_name,
        )

    all_pts: list[np.ndarray] = []
    all_conn: list[np.ndarray] = []
    offset = 0
    for c_idx in range(bounds.shape[0]):
        pts, conn = _subdivide_box(bounds[c_idx, 0], bounds[c_idx, 1], sub)
        all_pts.append(pts)
        all_conn.append(conn + offset)
        offset += pts.shape[0]

    points = np.vstack(all_pts)
    conn_2d = np.vstack(all_conn)

    if deduplicate_tolerance is not None and points.shape[0] > 0:
        scale = 1.0 / float(deduplicate_tolerance)
        keys = np.round(points * scale).astype(np.int64)
        _, first_idx, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
        order = np.argsort(first_idx)
        remap = np.empty_like(order)
        remap[order] = np.arange(len(order), dtype=np.int64)
        points = points[first_idx[order]]
        conn_2d = remap[inverse[conn_2d]]

    n_hex = conn_2d.shape[0]
    return ProxyMesh(
        points=points,
        connectivity=conn_2d.ravel(),
        cell_types=np.full(n_hex, VTK_HEXAHEDRON, dtype=np.uint8),
        block_name=block_name,
    )


def generate_cut_tet_proxy(
    cell_bounds: np.ndarray | None = None,
    n_sub: int | tuple[int, int, int] = 2,
    inside_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    phi_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    cut_points: np.ndarray | None = None,
    cut_tets: np.ndarray | None = None,
    block_name: str = "cut_cells",
) -> ProxyMesh:
    """Generate a ``VTK_TETRA`` proxy mesh for ``CUT`` cells.

    Supports either pre-supplied cut tetrahedra (``cut_points`` + ``cut_tets``)
    or automatic 6-tet Kuhn decomposition of subdivided ``cell_bounds``,
    filtered by ``inside_fn(pts) -> bool`` or clipped by ``phi_fn(pts) <= 0``.
    """
    if cut_points is not None and cut_tets is not None:
        pts = np.asarray(cut_points, dtype=np.float64)
        tets = np.asarray(cut_tets, dtype=np.int64)
        if tets.ndim != 2 or tets.shape[1] != 4:
            raise ValueError(f"cut_tets must have shape (M, 4), got {tets.shape}")
        return ProxyMesh(
            points=pts,
            connectivity=tets.ravel(),
            cell_types=np.full(tets.shape[0], VTK_TETRA, dtype=np.uint8),
            block_name=block_name,
        )

    if cell_bounds is None:
        raise ValueError("Provide either (cut_points, cut_tets) or cell_bounds")

    bounds = _normalize_bounds(cell_bounds)
    sub = _normalize_sub(n_sub)
    if bounds.shape[0] == 0:
        return ProxyMesh(
            points=np.zeros((0, 3), dtype=np.float64),
            connectivity=np.zeros(0, dtype=np.int64),
            cell_types=np.zeros(0, dtype=np.uint8),
            block_name=block_name,
        )

    all_pts: list[np.ndarray] = []
    all_tets: list[np.ndarray] = []
    offset = 0
    for c_idx in range(bounds.shape[0]):
        pts, hex_conn = _subdivide_box(bounds[c_idx, 0], bounds[c_idx, 1], sub)
        tets = decompose_hex_to_6_tets(hex_conn)
        all_pts.append(pts)
        all_tets.append(tets + offset)
        offset += pts.shape[0]

    points = np.vstack(all_pts)
    tets_all = np.vstack(all_tets)

    if phi_fn is not None:
        phi = np.asarray(phi_fn(points), dtype=np.float64).ravel()
        clipped_pts: list[np.ndarray] = []
        clipped_tets: list[tuple[int, int, int, int]] = []
        for t in tets_all:
            clipped_tets.extend(
                _clip_tet_negative_halfspace(clipped_pts, points[t], phi[t])
            )
        if not clipped_tets:
            return ProxyMesh(
                points=np.zeros((0, 3), dtype=np.float64),
                connectivity=np.zeros(0, dtype=np.int64),
                cell_types=np.zeros(0, dtype=np.uint8),
                block_name=block_name,
            )
        points = np.vstack(clipped_pts)
        tets_all = np.asarray(clipped_tets, dtype=np.int64)
    elif inside_fn is not None:
        centroids = points[tets_all].mean(axis=1)
        mask = np.asarray(inside_fn(centroids), dtype=bool).ravel()
        tets_all = tets_all[mask]
        if tets_all.shape[0] == 0:
            return ProxyMesh(
                points=np.zeros((0, 3), dtype=np.float64),
                connectivity=np.zeros(0, dtype=np.int64),
                cell_types=np.zeros(0, dtype=np.uint8),
                block_name=block_name,
            )
        used, inv = np.unique(tets_all.ravel(), return_inverse=True)
        points = points[used]
        tets_all = inv.reshape((-1, 4))

    return ProxyMesh(
        points=points,
        connectivity=tets_all.ravel(),
        cell_types=np.full(tets_all.shape[0], VTK_TETRA, dtype=np.uint8),
        block_name=block_name,
    )


def merge_proxy_meshes(meshes: list[ProxyMesh], block_name: str = "iga_proxy") -> ProxyMesh:
    """Merge multiple :class:`ProxyMesh` blocks into one combined mesh.

    Stores integer ``BlockId`` in ``cell_data`` when more than one block is
    provided.
    """
    non_empty = [m for m in meshes if m.n_cells > 0]
    if not non_empty:
        return ProxyMesh(
            points=np.zeros((0, 3), dtype=np.float64),
            connectivity=np.zeros(0, dtype=np.int64),
            cell_types=np.zeros(0, dtype=np.uint8),
            block_name=block_name,
        )

    pts_list: list[np.ndarray] = []
    conn_list: list[np.ndarray] = []
    types_list: list[np.ndarray] = []
    block_ids: list[np.ndarray] = []
    pt_offset = 0

    shared_pt_keys = set.intersection(*(set(m.point_data.keys()) for m in non_empty))
    shared_cell_keys = set.intersection(*(set(m.cell_data.keys()) for m in non_empty))

    for bid, m in enumerate(non_empty):
        pts_list.append(m.points)
        conn_list.append(m.connectivity + pt_offset)
        types_list.append(m.cell_types)
        block_ids.append(np.full(m.n_cells, float(bid), dtype=np.float64))
        pt_offset += m.n_points

    point_data = {
        k: np.concatenate([m.point_data[k] for m in non_empty], axis=0)
        for k in sorted(shared_pt_keys)
    }
    cell_data = {
        k: np.concatenate([m.cell_data[k] for m in non_empty], axis=0)
        for k in sorted(shared_cell_keys)
    }
    if len(non_empty) > 1 and "BlockId" not in cell_data:
        cell_data["BlockId"] = np.concatenate(block_ids, axis=0)

    return ProxyMesh(
        points=np.vstack(pts_list),
        connectivity=np.concatenate(conn_list, axis=0),
        cell_types=np.concatenate(types_list, axis=0),
        point_data=point_data,
        cell_data=cell_data,
        block_name=block_name,
    )


def build_proxy_mesh(
    in_bounds: np.ndarray | None = None,
    cut_bounds: np.ndarray | None = None,
    n_sub_in: int | tuple[int, int, int] = 2,
    n_sub_cut: int | tuple[int, int, int] = 2,
    inside_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    phi_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    cut_points: np.ndarray | None = None,
    cut_tets: np.ndarray | None = None,
    evaluator: BSplineFieldEvaluator | None = None,
    cp_fields: dict[str, np.ndarray] | None = None,
    block_name: str = "iga_proxy",
) -> ProxyMesh:
    """Build a combined hex (IN) + tet (CUT) proxy mesh and evaluate B-spline fields.

    Parameters
    ----------
    in_bounds:
        Bounds ``(N_in, 2, 3)`` or ``(N_in, 6)`` of ``IN`` cells.
    cut_bounds:
        Bounds ``(N_cut, 2, 3)`` or ``(N_cut, 6)`` of ``CUT`` cells.
    n_sub_in, n_sub_cut:
        Subdivision counts for ``IN`` and ``CUT`` cells.
    inside_fn, phi_fn:
        Optional implicit geometry predicate or signed level-set for ``CUT`` cells.
    cut_points, cut_tets:
        Optional pre-computed cut tetrahedra overriding ``cut_bounds`` subdivision.
    evaluator:
        Optional :class:`BSplineFieldEvaluator` to sample ``cp_fields`` onto
        the proxy mesh vertices.
    cp_fields:
        Mapping of field name to control-point array evaluated at proxy points.
    block_name:
        Name of the resulting merged proxy block.
    """
    blocks: list[ProxyMesh] = []
    if in_bounds is not None and len(in_bounds) > 0:
        blocks.append(generate_hex_proxy_for_cells(in_bounds, n_sub=n_sub_in, block_name="in_cells"))
    if (cut_points is not None and cut_tets is not None) or (cut_bounds is not None and len(cut_bounds) > 0):
        blocks.append(
            generate_cut_tet_proxy(
                cell_bounds=cut_bounds,
                n_sub=n_sub_cut,
                inside_fn=inside_fn,
                phi_fn=phi_fn,
                cut_points=cut_points,
                cut_tets=cut_tets,
                block_name="cut_cells",
            )
        )

    mesh = merge_proxy_meshes(blocks, block_name=block_name)
    if evaluator is not None and cp_fields and mesh.n_points > 0:
        for name, cp_arr in cp_fields.items():
            mesh.point_data[name] = evaluator.evaluate(mesh.points, cp_arr)
    return mesh
