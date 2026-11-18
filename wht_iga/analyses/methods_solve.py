"""Discretization solvers and unified facade (split from methods.py)."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable, Literal

import numpy as np

from wht_iga.analyses.methods_options import (
    DiscretizationOption,
    FullIgaOption,
    IgaFeHybridOption,
    MethodId,
    NativeFeOption,
    ThbAdaptiveIgaOption,
)


@dataclass(frozen=True)
class BoxMethodProblem:
    """Axis-aligned 3D linear elasticity benchmark problem for method comparison.

    Parameters
    ----------
    lo:
        Lower box corner ``(x0, y0, z0)``.
    hi:
        Upper box corner ``(x1, y1, z1)``.
    E:
        Young's modulus (``> 0``).
    nu:
        Poisson's ratio in ``(-1, 0.5)``.
    fixed_face:
        Dirichlet clamped face ``(axis, side)`` where ``side`` is ``0`` (lo) or ``1`` (hi).
    load_face:
        Loaded face ``(axis, side)`` where ``traction`` is applied.
    traction:
        Uniform surface traction vector ``(tx, ty, tz)`` over ``load_face``.
    tess_soup:
        Optional ``(nodes, tris)`` surface triangle soup used by ``"full_iga_immersed"``.
    """

    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    E: float
    nu: float
    fixed_face: tuple[int, int] = (0, 0)
    load_face: tuple[int, int] = (0, 1)
    traction: tuple[float, float, float] = (0.0, 0.0, 0.0)
    tess_soup: tuple[np.ndarray, np.ndarray] | None = None

    def __post_init__(self) -> None:
        lo_t = tuple(float(v) for v in self.lo)
        hi_t = tuple(float(v) for v in self.hi)
        tr_t = tuple(float(v) for v in self.traction)
        if len(lo_t) != 3 or len(hi_t) != 3 or len(tr_t) != 3:
            raise ValueError("lo, hi, and traction must each have length 3")
        if any(not np.isfinite(v) for v in (*lo_t, *hi_t, *tr_t)):
            raise ValueError("lo, hi, and traction values must be finite")
        if any(h <= l for l, h in zip(lo_t, hi_t)):
            raise ValueError(f"hi must exceed lo on every axis, got lo={lo_t}, hi={hi_t}")
        if not (np.isfinite(self.E) and float(self.E) > 0.0):
            raise ValueError(f"E must be finite and > 0, got {self.E}")
        if not (np.isfinite(self.nu) and -1.0 < float(self.nu) < 0.5):
            raise ValueError(f"nu must be in (-1, 0.5), got {self.nu}")
        for name, face in (("fixed_face", self.fixed_face), ("load_face", self.load_face)):
            if len(face) != 2 or face[0] not in (0, 1, 2) or face[1] not in (0, 1):
                raise ValueError(f"{name} must be (axis in 0..2, side in 0..1), got {face}")
        object.__setattr__(self, "lo", lo_t)
        object.__setattr__(self, "hi", hi_t)
        object.__setattr__(self, "traction", tr_t)
        object.__setattr__(self, "E", float(self.E))
        object.__setattr__(self, "nu", float(self.nu))
        object.__setattr__(self, "fixed_face", (int(self.fixed_face[0]), int(self.fixed_face[1])))
        object.__setattr__(self, "load_face", (int(self.load_face[0]), int(self.load_face[1])))
        if self.tess_soup is not None:
            nodes, tris = self.tess_soup
            n_arr = np.asarray(nodes, dtype=np.float64)
            t_arr = np.asarray(tris, dtype=np.int64)
            if n_arr.ndim != 2 or n_arr.shape[1] != 3:
                raise ValueError(f"tess_soup nodes must have shape (N, 3), got {n_arr.shape}")
            if t_arr.ndim != 2 or t_arr.shape[1] != 3:
                raise ValueError(f"tess_soup tris must have shape (M, 3), got {t_arr.shape}")

    def tip_point(self) -> np.ndarray:
        """Return the geometric center of ``load_face``."""
        pt = [0.5 * (self.lo[d] + self.hi[d]) for d in range(3)]
        axis, side = self.load_face
        pt[axis] = self.lo[axis] if side == 0 else self.hi[axis]
        return np.asarray(pt, dtype=np.float64)


@dataclass
class MethodSolution:
    """Unified solution object returned by :func:`solve_with_discretization_method`.

    Attributes
    ----------
    method:
        Discretization method identifier used for the solve.
    total_dofs:
        Total displacement DOFs in the assembled discretization.
    active_dofs:
        Active/reduced displacement DOFs solved.
    tip_displacement:
        Evaluated displacement vector ``(3,)`` at the center of ``problem.load_face``.
    max_displacement:
        Maximum nodal/control-point displacement magnitude.
    converged:
        Whether the linear solver converged.
    timings:
        Wall-clock timing breakdown in seconds.
    native:
        Raw native solver outputs and mesh/model handles.
    """

    method: MethodId
    total_dofs: int
    active_dofs: int
    tip_displacement: np.ndarray
    max_displacement: float
    converged: bool
    timings: dict[str, Any] = field(default_factory=dict)
    native: dict[str, Any] = field(default_factory=dict)
    _eval_fn: Callable[[np.ndarray], np.ndarray] | None = field(default=None, repr=False)

    def displacement_at(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the solved displacement field at physical ``points`` ``(M, 3)`` or ``(3,)``."""
        pts = np.asarray(points, dtype=np.float64)
        single = pts.ndim == 1
        if single:
            pts = pts.reshape(1, 3)
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError(f"points must have shape (M, 3) or (3,), got {pts.shape}")
        if self._eval_fn is None:
            raise RuntimeError("MethodSolution has no displacement evaluator attached")
        vals = np.asarray(self._eval_fn(pts), dtype=np.float64)
        return vals[0] if single else vals


# ---------------------------------------------------------------------------
# Structured FE Mesh Generation, Consistent Loads, and Exact Interpolation
# ---------------------------------------------------------------------------

_HEX_CORNER_SIGNS = np.array(
    [
        [-1.0, -1.0, -1.0],
        [+1.0, -1.0, -1.0],
        [+1.0, +1.0, -1.0],
        [-1.0, +1.0, -1.0],
        [-1.0, -1.0, +1.0],
        [+1.0, -1.0, +1.0],
        [+1.0, +1.0, +1.0],
        [-1.0, +1.0, +1.0],
    ],
    dtype=np.float64,
)

# Edge midpoints in MeshModel canonical hex20 order (indices 8..19):
# (axis_zero, sign_a, sign_b) relative to the two non-zero axes in (xi, eta, zeta).
_HEX20_MID_SPECS: tuple[tuple[int, float, float], ...] = (
    # Bottom (zeta = -1): 8:(0-1), 9:(1-2), 10:(2-3), 11:(3-0)
    (0, -1.0, -1.0),
    (1, +1.0, -1.0),
    (0, +1.0, -1.0),
    (1, -1.0, -1.0),
    # Top (zeta = +1): 12:(4-5), 13:(5-6), 14:(6-7), 15:(7-4)
    (0, -1.0, +1.0),
    (1, +1.0, +1.0),
    (0, +1.0, +1.0),
    (1, -1.0, +1.0),
    # Vertical (zeta = 0): 16:(0-4), 17:(1-5), 18:(2-6), 19:(3-7)
    (2, -1.0, -1.0),
    (2, +1.0, -1.0),
    (2, +1.0, +1.0),
    (2, -1.0, +1.0),
)

_FE_FACE_TABLES: dict[str, tuple[tuple[str, tuple[int, ...]], ...]] = {
    "hex8": (
        ("quad4", (0, 3, 2, 1)),
        ("quad4", (4, 5, 6, 7)),
        ("quad4", (0, 1, 5, 4)),
        ("quad4", (1, 2, 6, 5)),
        ("quad4", (2, 3, 7, 6)),
        ("quad4", (3, 0, 4, 7)),
    ),
    "hex20": (
        ("quad8", (0, 3, 2, 1, 11, 10, 9, 8)),
        ("quad8", (4, 5, 6, 7, 12, 13, 14, 15)),
        ("quad8", (0, 1, 5, 4, 8, 17, 12, 16)),
        ("quad8", (1, 2, 6, 5, 9, 18, 13, 17)),
        ("quad8", (2, 3, 7, 6, 10, 19, 14, 18)),
        ("quad8", (3, 0, 4, 7, 11, 16, 15, 19)),
    ),
    "tet4": (
        ("tri3", (0, 2, 1)),
        ("tri3", (0, 1, 3)),
        ("tri3", (1, 2, 3)),
        ("tri3", (2, 0, 3)),
    ),
    "tet10": (
        ("tri6", (0, 2, 1, 6, 5, 4)),
        ("tri6", (0, 1, 3, 4, 8, 7)),
        ("tri6", (1, 2, 3, 5, 9, 8)),
        ("tri6", (2, 0, 3, 6, 7, 9)),
    ),
    "wedge6": (
        ("tri3", (0, 2, 1)),
        ("tri3", (3, 4, 5)),
        ("quad4", (0, 1, 4, 3)),
        ("quad4", (1, 2, 5, 4)),
        ("quad4", (2, 0, 3, 5)),
    ),
}


def _structured_box_mesh(
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    n_cells: tuple[int, int, int],
    element_type: FeElementType,
    part_id: int = 1,
) -> Any:
    """Build a canonical structured :class:`MeshModel` for any of the 5 solid element types."""
    from wht_iga.meshio_ext.model import ElemBlock, MeshModel

    nx, ny, nz = (int(v) for v in n_cells)
    x0, y0, z0 = (float(v) for v in lo)
    x1, y1, z1 = (float(v) for v in hi)

    if element_type == "hex20":
        xs = np.linspace(x0, x1, 2 * nx + 1)
        ys = np.linspace(y0, y1, 2 * ny + 1)
        zs = np.linspace(z0, z1, 2 * nz + 1)
        coord_to_id: dict[tuple[int, int, int], int] = {}
        nodes_list: list[list[float]] = []

        def _pt(i: int, j: int, k: int) -> int:
            key = (i, j, k)
            if key not in coord_to_id:
                coord_to_id[key] = len(nodes_list)
                nodes_list.append([float(xs[i]), float(ys[j]), float(zs[k])])
            return coord_to_id[key]

        conn20 = []
        for kz in range(nz):
            for jy in range(ny):
                for ix in range(nx):
                    cx, cy, cz = 2 * ix, 2 * jy, 2 * kz
                    conn20.append([
                        _pt(cx, cy, cz), _pt(cx + 2, cy, cz), _pt(cx + 2, cy + 2, cz), _pt(cx, cy + 2, cz),
                        _pt(cx, cy, cz + 2), _pt(cx + 2, cy, cz + 2), _pt(cx + 2, cy + 2, cz + 2), _pt(cx, cy + 2, cz + 2),
                        _pt(cx + 1, cy, cz), _pt(cx + 2, cy + 1, cz), _pt(cx + 1, cy + 2, cz), _pt(cx, cy + 1, cz),
                        _pt(cx + 1, cy, cz + 2), _pt(cx + 2, cy + 1, cz + 2), _pt(cx + 1, cy + 2, cz + 2), _pt(cx, cy + 1, cz + 2),
                        _pt(cx, cy, cz + 1), _pt(cx + 2, cy, cz + 1), _pt(cx + 2, cy + 2, cz + 1), _pt(cx, cy + 2, cz + 1),
                    ])
        return MeshModel(
            nodes=np.asarray(nodes_list, dtype=np.float64),
            blocks=[ElemBlock("hex20", np.asarray(conn20, dtype=np.int64), part_id)],
        )

    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    zs = np.linspace(z0, z1, nz + 1)
    nodes = np.array([[x, y, z] for z in zs for y in ys for x in xs], dtype=np.float64)

    def _nid(i: int, j: int, k: int) -> int:
        return i + (nx + 1) * (j + (ny + 1) * k)

    hex_conn: list[list[int]] = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hex_conn.append([
                    _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1), _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                ])

    if element_type == "hex8":
        return MeshModel(nodes=nodes, blocks=[ElemBlock("hex8", np.asarray(hex_conn, dtype=np.int64), part_id)])

    if element_type == "wedge6":
        wconn: list[list[int]] = []
        for c in hex_conn:
            wconn.append([c[0], c[1], c[2], c[4], c[5], c[6]])
            wconn.append([c[0], c[2], c[3], c[4], c[6], c[7]])
        return MeshModel(nodes=nodes, blocks=[ElemBlock("wedge6", np.asarray(wconn, dtype=np.int64), part_id)])

    if element_type in ("tet4", "tet10"):
        tets4: list[list[int]] = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    v000 = _nid(i, j, k)
                    v100 = _nid(i + 1, j, k)
                    v010 = _nid(i, j + 1, k)
                    v110 = _nid(i + 1, j + 1, k)
                    v001 = _nid(i, j, k + 1)
                    v101 = _nid(i + 1, j, k + 1)
                    v011 = _nid(i, j + 1, k + 1)
                    v111 = _nid(i + 1, j + 1, k + 1)
                    raw_tets = [
                        [v000, v100, v110, v111],
                        [v000, v110, v010, v111],
                        [v000, v010, v011, v111],
                        [v000, v011, v001, v111],
                        [v000, v001, v101, v111],
                        [v000, v101, v100, v111],
                    ]
                    for rt in raw_tets:
                        jac = np.column_stack([
                            nodes[rt[1]] - nodes[rt[0]],
                            nodes[rt[2]] - nodes[rt[0]],
                            nodes[rt[3]] - nodes[rt[0]],
                        ])
                        if np.linalg.det(jac) < 0.0:
                            rt = [rt[0], rt[2], rt[1], rt[3]]
                        tets4.append(rt)
        if element_type == "tet4":
            return MeshModel(nodes=nodes, blocks=[ElemBlock("tet4", np.asarray(tets4, dtype=np.int64), part_id)])

        nodes_list = nodes.tolist()
        edge_mid: dict[tuple[int, int], int] = {}

        def _mid(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edge_mid:
                edge_mid[key] = len(nodes_list)
                nodes_list.append((0.5 * (nodes[a] + nodes[b])).tolist())
            return edge_mid[key]

        tets10 = [
            [a, b, c, d, _mid(a, b), _mid(b, c), _mid(c, a), _mid(a, d), _mid(b, d), _mid(c, d)]
            for a, b, c, d in tets4
        ]
        return MeshModel(
            nodes=np.asarray(nodes_list, dtype=np.float64),
            blocks=[ElemBlock("tet10", np.asarray(tets10, dtype=np.int64), part_id)],
        )

    raise ValueError(f"Unsupported FE element_type {element_type!r}")


def _hex8_shape(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Evaluate 8-node trilinear hexahedron shape functions on ``[-1, 1]^3``."""
    return (
        0.125
        * (1.0 + _HEX_CORNER_SIGNS[:, 0] * xi)
        * (1.0 + _HEX_CORNER_SIGNS[:, 1] * eta)
        * (1.0 + _HEX_CORNER_SIGNS[:, 2] * zeta)
    )


def _hex20_shape(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Evaluate 20-node serendipity hexahedron shape functions in canonical order."""
    N = np.empty(20, dtype=np.float64)
    sx = _HEX_CORNER_SIGNS[:, 0]
    sy = _HEX_CORNER_SIGNS[:, 1]
    sz = _HEX_CORNER_SIGNS[:, 2]
    N[:8] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta) * (
        sx * xi + sy * eta + sz * zeta - 2.0
    )
    coords = (xi, eta, zeta)
    for idx, (zero_ax, s_a, s_b) in enumerate(_HEX20_MID_SPECS):
        other = [a for a in (0, 1, 2) if a != zero_ax]
        N[8 + idx] = (
            0.25
            * (1.0 - coords[zero_ax] ** 2)
            * (1.0 + s_a * coords[other[0]])
            * (1.0 + s_b * coords[other[1]])
        )
    return N


def _eval_structured_fe_displacement(
    mesh: Any,
    element_type: FeElementType,
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    n_cells: tuple[int, int, int],
    u_nodal: np.ndarray,
    points: np.ndarray,
) -> np.ndarray:
    """Evaluate exact isoparametric FE displacement field at arbitrary physical ``points``."""
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    lo_arr = np.asarray(lo, dtype=np.float64)
    hi_arr = np.asarray(hi, dtype=np.float64)
    n_arr = np.asarray(n_cells, dtype=np.int64)
    h = (hi_arr - lo_arr) / n_arr
    conn = mesh.blocks[0].conn
    coords = mesh.nodes
    out = np.empty((len(pts), 3), dtype=np.float64)

    for pi, p in enumerate(pts):
        rel = (p - lo_arr) / h
        ix = int(np.clip(np.floor(rel[0]), 0, n_arr[0] - 1))
        iy = int(np.clip(np.floor(rel[1]), 0, n_arr[1] - 1))
        iz = int(np.clip(np.floor(rel[2]), 0, n_arr[2] - 1))
        cell_idx = ix + int(n_arr[0]) * (iy + int(n_arr[1]) * iz)

        cell_lo = lo_arr + np.array([ix, iy, iz], dtype=np.float64) * h
        uvw = np.clip((p - cell_lo) / h, 0.0, 1.0)
        xi, eta, zeta = 2.0 * uvw - 1.0

        if element_type == "hex8":
            N = _hex8_shape(float(xi), float(eta), float(zeta))
            out[pi] = N @ u_nodal[conn[cell_idx]]
        elif element_type == "hex20":
            N = _hex20_shape(float(xi), float(eta), float(zeta))
            out[pi] = N @ u_nodal[conn[cell_idx]]
        elif element_type == "wedge6":
            best_min = -float("inf")
            best_N = np.zeros(6, dtype=np.float64)
            best_row = conn[2 * cell_idx]
            for sub in range(2):
                row = conn[2 * cell_idx + sub]
                tri = coords[row[:3], :2]
                mat = np.column_stack([tri[1] - tri[0], tri[2] - tri[0]])
                l1, l2 = np.linalg.solve(mat, p[:2] - tri[0])
                l0 = 1.0 - l1 - l2
                m_val = min(float(l0), float(l1), float(l2))
                if m_val > best_min:
                    best_min = m_val
                    w = float(uvw[2])
                    best_N = np.array(
                        [l0 * (1.0 - w), l1 * (1.0 - w), l2 * (1.0 - w), l0 * w, l1 * w, l2 * w],
                        dtype=np.float64,
                    )
                    best_row = row
            out[pi] = best_N @ u_nodal[best_row]
        else:  # tet4 / tet10
            best_min = -float("inf")
            best_b = np.zeros(4, dtype=np.float64)
            best_row = conn[6 * cell_idx]
            for sub in range(6):
                row = conn[6 * cell_idx + sub]
                t_pts = coords[row[:4]]
                jac = np.column_stack([t_pts[1] - t_pts[0], t_pts[2] - t_pts[0], t_pts[3] - t_pts[0]])
                l1, l2, l3 = np.linalg.solve(jac, p - t_pts[0])
                l0 = 1.0 - l1 - l2 - l3
                m_val = min(float(l0), float(l1), float(l2), float(l3))
                if m_val > best_min:
                    best_min = m_val
                    best_b = np.array([l0, l1, l2, l3], dtype=np.float64)
                    best_row = row
            l0, l1, l2, l3 = best_b
            if element_type == "tet4":
                N_tet = best_b
            else:
                N_tet = np.array(
                    [
                        l0 * (2.0 * l0 - 1.0),
                        l1 * (2.0 * l1 - 1.0),
                        l2 * (2.0 * l2 - 1.0),
                        l3 * (2.0 * l3 - 1.0),
                        4.0 * l0 * l1,
                        4.0 * l1 * l2,
                        4.0 * l2 * l0,
                        4.0 * l0 * l3,
                        4.0 * l1 * l3,
                        4.0 * l2 * l3,
                    ],
                    dtype=np.float64,
                )
            out[pi] = N_tet @ u_nodal[best_row]

    return out


def _integrate_fe_face_traction(
    mesh: Any,
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    axis: int,
    side: int,
    traction: tuple[float, float, float],
    mesh_to_kratos: dict[int, int],
    tol: float = 1e-10,
) -> dict[int, np.ndarray]:
    """Compute exact consistent nodal traction loads on a planar box boundary face."""
    from wht_iga.kratos_backend.interactions import _fe_face_shape_and_deriv, _unit_tri_gauss

    t_vec = np.asarray(traction, dtype=np.float64)
    if float(np.linalg.norm(t_vec)) <= 1e-18:
        return {}

    target = float(lo[axis] if side == 0 else hi[axis])
    coords = mesh.nodes
    loads: dict[int, np.ndarray] = {}

    gx2, gw2 = np.polynomial.legendre.leggauss(2)
    gx3, gw3 = np.polynomial.legendre.leggauss(3)
    tri2_pts, tri2_wts = _unit_tri_gauss(2)
    tri4_pts, tri4_wts = _unit_tri_gauss(4)

    for block in mesh.blocks:
        ftab = _FE_FACE_TABLES[block.elem_type]
        for row in block.conn:
            for fkind, lidx in ftab:
                n_corner = 3 if fkind.startswith("tri") else 4
                fnodes = [int(row[k]) for k in lidx]
                if not all(abs(float(coords[v, axis]) - target) <= tol for v in fnodes[:n_corner]):
                    continue
                fpts = coords[fnodes]
                if fkind in ("quad4", "quad8"):
                    gx, gw = (gx2, gw2) if fkind == "quad4" else (gx3, gw3)
                    for eta, we in zip(gx, gw):
                        for xi, wx in zip(gx, gw):
                            N, dxi, deta = _fe_face_shape_and_deriv(fkind, float(xi), float(eta))
                            dA = float(np.linalg.norm(np.cross(dxi @ fpts, deta @ fpts)))
                            w_q = float(wx * we) * dA
                            for a, m_idx in enumerate(fnodes):
                                knid = int(mesh_to_kratos[m_idx])
                                loads[knid] = loads.get(knid, np.zeros(3)) + float(N[a]) * w_q * t_vec
                else:
                    pts_ref, wts_ref = (tri2_pts, tri2_wts) if fkind == "tri3" else (tri4_pts, tri4_wts)
                    for bary, wt in zip(pts_ref, wts_ref):
                        xi, eta = float(bary[1]), float(bary[2])
                        N, dxi, deta = _fe_face_shape_and_deriv(fkind, xi, eta)
                        dA2 = float(np.linalg.norm(np.cross(dxi @ fpts, deta @ fpts)))
                        w_q = float(wt) * dA2
                        for a, m_idx in enumerate(fnodes):
                            knid = int(mesh_to_kratos[m_idx])
                            loads[knid] = loads.get(knid, np.zeros(3)) + float(N[a]) * w_q * t_vec

    return loads


def _integrate_iga_face_traction(
    box: Any,
    axis: int,
    side: int,
    traction: tuple[float, float, float],
) -> dict[int, np.ndarray]:
    """Compute consistent control-point traction loads on a ``BoxSpline`` face."""
    from wht_iga.grid.bspline import basis_1d

    t_vec = np.asarray(traction, dtype=np.float64)
    if float(np.linalg.norm(t_vec)) <= 1e-18:
        return {}

    in_axes = [a for a in (0, 1, 2) if a != axis]
    a0, a1 = in_axes[0], in_axes[1]
    bps0 = box.breakpoints(a0)
    bps1 = box.breakpoints(a1)
    p = int(box.degree)
    gx, gw = np.polynomial.legendre.leggauss(p + 1)
    nu, nv, _ = (int(v) for v in box.n_cp)
    knots_list = box.knots
    c_fixed = float(box.lo[axis] if side == 0 else box.hi[axis])

    s_fix, n_fix_arr = basis_1d(knots_list[axis], p, c_fixed, 0)
    i_fix0 = s_fix - p
    n_fix = n_fix_arr[0]

    loads: dict[int, np.ndarray] = {}
    for j_cell in range(len(bps1) - 1):
        v0, v1 = float(bps1[j_cell]), float(bps1[j_cell + 1])
        pts1 = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gx
        wts1 = 0.5 * (v1 - v0) * gw
        for i_cell in range(len(bps0) - 1):
            u0, u1 = float(bps0[i_cell]), float(bps0[i_cell + 1])
            pts0 = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gx
            wts0 = 0.5 * (u1 - u0) * gw
            for q1, w1 in zip(pts1, wts1):
                s1, n1_arr = basis_1d(knots_list[a1], p, float(q1), 0)
                i1_0 = s1 - p
                n1 = n1_arr[0]
                for q0, w0 in zip(pts0, wts0):
                    s0, n0_arr = basis_1d(knots_list[a0], p, float(q0), 0)
                    i0_0 = s0 - p
                    n0 = n0_arr[0]
                    w_area = float(w0 * w1)
                    for lf in range(p + 1):
                        if abs(float(n_fix[lf])) < 1e-16:
                            continue
                        idx_fix = i_fix0 + lf
                        for l1 in range(p + 1):
                            idx1 = i1_0 + l1
                            for l0 in range(p + 1):
                                idx0 = i0_0 + l0
                                weight = float(n_fix[lf] * n1[l1] * n0[l0]) * w_area
                                if abs(weight) <= 1e-18:
                                    continue
                                coords = [0, 0, 0]
                                coords[axis] = idx_fix
                                coords[a0] = idx0
                                coords[a1] = idx1
                                flat = coords[0] + nu * (coords[1] + nv * coords[2])
                                loads[flat] = loads.get(flat, np.zeros(3)) + weight * t_vec
    return loads


def _attach_point_loads(
    SM: Any,
    built: Any,
    loads_by_kratos_node: dict[int, np.ndarray],
    props: Any,
    scope: str = "Loads.point",
) -> None:
    """Attach nodal ``POINT_LOAD`` + ``PointLoadCondition3D1N`` entries via ``built.id_table``."""
    active = [(nid, fv) for nid, fv in loads_by_kratos_node.items() if float(np.linalg.norm(fv)) > 1e-18]
    if not active:
        return
    cid_entry = built.id_table.allocate_conditions(scope, len(active), kind="PointLoadCondition3D1N")
    cid = cid_entry.start_id
    for nid, fvec in active:
        nd = built.structure.GetNode(int(nid))
        cur = np.array(nd.GetSolutionStepValue(SM.POINT_LOAD, 0), dtype=np.float64)
        nd.SetSolutionStepValue(SM.POINT_LOAD, 0, [float(x) for x in (cur + fvec)])
        built.structure.CreateNewCondition("PointLoadCondition3D1N", cid, [int(nid)], props)
        cid += 1


def _fix_node_all_disp(KM: Any, node: Any) -> None:
    """Fix all three translational displacement DOFs on a Kratos node."""
    node.Fix(KM.DISPLACEMENT_X)
    node.Fix(KM.DISPLACEMENT_Y)
    node.Fix(KM.DISPLACEMENT_Z)


# ---------------------------------------------------------------------------
# Per-Method Solver Routing
# ---------------------------------------------------------------------------

def _solve_full_iga(
    problem: BoxMethodProblem,
    option: FullIgaOption,
    use_pardiso: bool,
) -> MethodSolution:
    """Route ``full_iga_body_fit`` and ``full_iga_immersed`` to ``wht_iga.solver.pipeline.run_linear``."""
    from wht_iga.backends.mesh_engine import box_triangle_soup
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.trimmed import build_trimmed_table
    from wht_iga.solver.validate import DirichletBC, JobSpec, TractionLoad

    t_start = time.perf_counter()
    box = BoxSpline(
        lo=problem.lo,
        hi=problem.hi,
        n_cells=option.n_cells,
        degree=option.degree,
        continuity=option.resolved_continuity,
    )

    qp_table = None
    out_cells: list[tuple[int, int, int]] = []
    geo_issues: list[dict] | None = None
    cell_records: list[dict] | None = None

    if option.kind == "full_iga_immersed":
        if problem.tess_soup is not None:
            tess_nodes, tess_tris = problem.tess_soup
        else:
            tess_nodes, tess_tris = box_triangle_soup(problem.lo, problem.hi)
        qp_table, out_cells, geo_issues, cell_records, _ = build_trimmed_table(
            tess_nodes,
            tess_tris,
            box,
            q=option.resolved_quadrature_degree,
        )

    tractions = [TractionLoad(face=problem.load_face, vector=problem.traction)]
    job_spec = JobSpec(
        name=option.kind,
        lo=problem.lo,
        hi=problem.hi,
        n_cells=option.n_cells,
        degree=option.degree,
        continuity=option.resolved_continuity,
        E=problem.E,
        nu=problem.nu,
        beta=option.bc_penalty,
        bcs=[DirichletBC(face=problem.fixed_face, components=(True, True, True), value=(0.0, 0.0, 0.0))],
        tractions=tractions,
        out_cells=out_cells,
        ghost_alpha=option.ghost_alpha,
        out_mode=option.out_mode,
        linear_solver="pardiso_ldlt" if use_pardiso else "auto",
        enable_solver_log=False,
    )
    res = run_linear(
        job_spec,
        geo_issues=geo_issues,
        qp_table=qp_table,
        cell_records=cell_records,
    )
    coeffs = res.U.reshape(-1, 3)

    def _eval(pts: np.ndarray) -> np.ndarray:
        return box.evaluate(coeffs, pts)

    tip_disp = _eval(problem.tip_point().reshape(1, 3))[0]
    total_wall = time.perf_counter() - t_start
    timings = dict(res.timings)
    timings["total_s"] = total_wall

    return MethodSolution(
        method=option.kind,
        total_dofs=int(res.total_dofs or box.n_dofs),
        active_dofs=int(res.active_dofs or box.n_dofs),
        tip_displacement=tip_disp,
        max_displacement=float(np.max(np.abs(res.U))) if res.U.size else 0.0,
        converged=bool(res.converged),
        timings=timings,
        native={"linear_result": res, "box": box},
        _eval_fn=_eval,
    )


def _solve_thb_adaptive(
    problem: BoxMethodProblem,
    option: ThbAdaptiveIgaOption,
    use_pardiso: bool,
) -> MethodSolution:
    """Route ``thb_adaptive_iga`` to ``THBHierarchy`` + ``AssemblyBuilder`` + ``assemble_and_solve_thb_reduced``."""
    from wht_iga.grid.hierarchy import THBHierarchy
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import AssemblyBuilder, AssemblySpec, InstanceSpec, PartSpec
    from wht_iga.solver.penalty import attach_penalty_bcs
    from wht_iga.solver.validate import DirichletBC

    t0 = time.perf_counter()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    IGA = mods["IgaApplication"]

    thb = THBHierarchy.from_uniform(
        (problem.lo[0], problem.hi[0]),
        (problem.lo[1], problem.hi[1]),
        (problem.lo[2], problem.hi[2]),
        option.coarse_n_cells,
        degree=option.degree,
    )

    if option.refine_box is not None:
        thb.mark_box(option.refine_box[0], option.refine_box[1], target_level=option.refinement_level)
    if option.refine_sphere is not None:
        thb.mark_sphere(option.refine_sphere[0], option.refine_sphere[1], target_level=option.refinement_level)
    if option.refine_box is None and option.refine_sphere is None:
        # Default adaptive zone: 40% adjacent to fixed_face where gradients/bending concentrate.
        r_lo = list(problem.lo)
        r_hi = list(problem.hi)
        fix_ax, fix_side = problem.fixed_face
        span = r_hi[fix_ax] - r_lo[fix_ax]
        if fix_side == 0:
            r_hi[fix_ax] = r_lo[fix_ax] + 0.4 * span
        else:
            r_lo[fix_ax] = r_hi[fix_ax] - 0.4 * span
        thb.mark_box(tuple(r_lo), tuple(r_hi), target_level=option.refinement_level)

    if option.balance == "two_to_one":
        thb.enforce_2_to_1_balance()

    part = PartSpec(
        name="domain",
        kind="thb",
        thb=thb,
        young_modulus=problem.E,
        poisson_ratio=problem.nu,
    )
    spec = AssemblySpec(
        job_name="thb_adaptive_iga",
        parts={"domain": part},
        instances=[InstanceSpec("domain", "1")],
        use_pardiso=use_pardiso,
    )
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(spec)
    rt = built.get_instance("domain_1")
    assert rt.box is not None
    box_fine = rt.box

    keep_bc = None
    if option.penalty_beta is not None:
        bc = DirichletBC(face=problem.fixed_face, components=(True, True, True), value=(0.0, 0.0, 0.0))
        keep_bc = attach_penalty_bcs(
            KM,
            IGA,
            built.structure,
            rt.volume,
            box_fine,
            [bc],
            float(option.penalty_beta),
            first_cond_id=built.id_table.next_id("conditions") + 100,
            first_node_id=built.id_table.next_id("nodes") + 100,
            first_prop_id=built.id_table.next_id("properties") + 10,
        )
    else:
        for idx in box_fine.boundary_control_points(*problem.fixed_face):
            _fix_node_all_disp(KM, built.structure.GetNode(rt.node_ids[int(idx)]))

    f_local = _integrate_iga_face_traction(box_fine, *problem.load_face, problem.traction)
    f_kratos = {rt.node_ids[flat]: vec for flat, vec in f_local.items()}
    _attach_point_loads(SM, built, f_kratos, rt.properties_map[0], scope="Loads.thb")

    t1 = time.perf_counter()
    res = builder.assemble_and_solve_thb_reduced(built, "domain_1", use_pardiso=use_pardiso)
    t2 = time.perf_counter()

    u_hier = np.asarray(res["u_hier"], dtype=np.float64).reshape(-1, 3)

    def _eval(pts: np.ndarray) -> np.ndarray:
        return box_fine.evaluate(u_hier, pts)

    tip_disp = _eval(problem.tip_point().reshape(1, 3))[0]
    n_thb = int(thb.n_thb_dofs)

    return MethodSolution(
        method="thb_adaptive_iga",
        total_dofs=3 * n_thb,
        active_dofs=int(res["u_red"].shape[0]),
        tip_displacement=tip_disp,
        max_displacement=float(np.max(np.abs(u_hier))) if u_hier.size else 0.0,
        converged=True,
        timings={"setup_s": t1 - t0, "solve_s": t2 - t1, "total_s": t2 - t0},
        native={"built": built, "thb": thb, "box_fine": box_fine, "thb_result": res, "keep_bc": keep_bc},
        _eval_fn=_eval,
    )


def _solve_native_fe(
    problem: BoxMethodProblem,
    option: NativeFeOption,
    use_pardiso: bool,
) -> MethodSolution:
    """Route ``native_fe`` to structured ``MeshModel`` + ``AssemblyBuilder`` + ``solve_linear``."""
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import AssemblyBuilder, AssemblySpec, InstanceSpec, PartSpec
    from wht_iga.kratos_backend.elements import solve_linear

    t0 = time.perf_counter()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    mesh = _structured_box_mesh(problem.lo, problem.hi, option.n_cells, option.element_type)
    part = PartSpec(
        name="fe_body",
        kind="native_fe",
        mesh=mesh,
        young_modulus=problem.E,
        poisson_ratio=problem.nu,
    )
    spec = AssemblySpec(
        job_name="native_fe",
        parts={"fe_body": part},
        instances=[InstanceSpec("fe_body", "1")],
        use_pardiso=use_pardiso,
    )
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(spec)
    rt = built.get_instance("fe_body_1")
    assert rt.fe_handle is not None

    fix_ax, fix_side = problem.fixed_face
    fix_coord = problem.lo[fix_ax] if fix_side == 0 else problem.hi[fix_ax]
    for nid in rt.node_ids:
        nd = built.structure.GetNode(nid)
        pos = (nd.X0, nd.Y0, nd.Z0)
        if abs(pos[fix_ax] - fix_coord) <= 1e-10:
            _fix_node_all_disp(KM, nd)

    fe_loads = _integrate_fe_face_traction(
        mesh,
        problem.lo,
        problem.hi,
        problem.load_face[0],
        problem.load_face[1],
        problem.traction,
        rt.fe_handle.mesh_to_kratos_node,
    )
    _attach_point_loads(SM, built, fe_loads, rt.properties_map[0], scope="Loads.fe")

    t1 = time.perf_counter()
    strategy = solve_linear(KM, built.structure, use_pardiso=use_pardiso, echo=0)
    t2 = time.perf_counter()

    u_nodal = np.array(
        [
            list(built.structure.GetNode(rt.fe_handle.mesh_to_kratos_node[m_idx]).GetSolutionStepValue(KM.DISPLACEMENT, 0))
            for m_idx in range(len(mesh.nodes))
        ],
        dtype=np.float64,
    )

    def _eval(pts: np.ndarray) -> np.ndarray:
        return _eval_structured_fe_displacement(
            mesh,
            option.element_type,
            problem.lo,
            problem.hi,
            option.n_cells,
            u_nodal,
            pts,
        )

    tip_disp = _eval(problem.tip_point().reshape(1, 3))[0]
    total_dofs = 3 * len(rt.node_ids)

    return MethodSolution(
        method="native_fe",
        total_dofs=total_dofs,
        active_dofs=total_dofs,
        tip_displacement=tip_disp,
        max_displacement=float(np.max(np.abs(u_nodal))) if u_nodal.size else 0.0,
        converged=True,
        timings={"setup_s": t1 - t0, "solve_s": t2 - t1, "total_s": t2 - t0},
        native={"built": built, "mesh": mesh, "strategy": strategy, "u_nodal": u_nodal},
        _eval_fn=_eval,
    )


def _subdomain_touches_face(
    face: tuple[int, int],
    split_axis: int,
    side_on_split: Literal["lo", "hi"],
) -> bool:
    """Return whether a subdomain on ``side_on_split`` of ``split_axis`` touches box ``face``."""
    ax, side = face
    if ax != split_axis:
        return True
    return (side == 0 and side_on_split == "lo") or (side == 1 and side_on_split == "hi")


def _apply_iga_barycentric_tie(
    KM: Any,
    built: Any,
    surf_iga: Any,
    surf_fe: Any,
    scope: str = "Couplings.hyb_barycentric",
) -> list[Any]:
    """Attach exact spline-evaluation barycentric MPC ties from IGA master to FE slave nodes."""
    from wht_iga.kratos_backend.interactions import LinearRelation, apply_relations_to_model_part

    assert surf_iga.eval_basis_fn is not None and surf_iga.dof_node_ids is not None
    assert surf_fe.node_ids is not None
    Nm = surf_iga.eval_basis_fn(surf_fe.coords).tocsr()
    comp_names = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    rels: list[LinearRelation] = []
    for si, snid in enumerate(surf_fe.node_ids):
        s_node = built.structure.GetNode(int(snid))
        row_start, row_end = Nm.indptr[si], Nm.indptr[si + 1]
        cols = Nm.indices[row_start:row_end]
        vals = Nm.data[row_start:row_end]
        keep = np.abs(vals) > 1e-14
        cols = cols[keep]
        vals = vals[keep]
        if len(cols) == 0:
            continue
        for cname, cvar in zip(comp_names, comp_vars):
            if s_node.IsFixed(cvar):
                continue
            terms = {(int(surf_iga.dof_node_ids[int(c)]), cname): float(w) for c, w in zip(cols, vals)}
            rels.append(
                LinearRelation(
                    slave_dof=(int(snid), cname),
                    terms=terms,
                    constant=0.0,
                    tag="hyb_barycentric",
                )
            )
    return apply_relations_to_model_part(KM, built.structure, rels, id_table=built.id_table, scope=scope)


def _solve_iga_fe_hybrid(
    problem: BoxMethodProblem,
    option: IgaFeHybridOption,
    use_pardiso: bool,
) -> MethodSolution:
    """Route ``iga_fe_hybrid`` to ``AssemblyBuilder`` + ``SurfaceMesh`` + ``solve_hybrid_iga_fe``."""
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import AssemblyBuilder, AssemblySpec, InstanceSpec, PartSpec
    from wht_iga.kratos_backend.interactions import SurfaceMesh, TieSpec
    from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe
    from wht_iga.reference.box import BoxSpline

    t0 = time.perf_counter()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    ax = int(option.split_axis)
    x_split = float(problem.lo[ax] + option.split_fraction * (problem.hi[ax] - problem.lo[ax]))
    lo_sub = list(problem.lo)
    hi_sub = list(problem.hi)

    lo_box_lo = tuple(lo_sub)
    lo_box_hi_l = list(hi_sub)
    lo_box_hi_l[ax] = x_split
    lo_box_hi = tuple(lo_box_hi_l)

    hi_box_lo_l = list(lo_sub)
    hi_box_lo_l[ax] = x_split
    hi_box_lo = tuple(hi_box_lo_l)
    hi_box_hi = tuple(hi_sub)

    if option.fe_side == "lo":
        fe_lo, fe_hi = lo_box_lo, lo_box_hi
        iga_lo, iga_hi = hi_box_lo, hi_box_hi
        iga_iface_side: Literal["lo", "hi"] = "lo"
        iga_side_tag: Literal["lo", "hi"] = "hi"
    else:
        iga_lo, iga_hi = lo_box_lo, lo_box_hi
        fe_lo, fe_hi = hi_box_lo, hi_box_hi
        iga_iface_side = "hi"
        iga_side_tag = "lo"

    fe_n_cells = option.fe_n_cells or (2, 2, 2)
    iga_n_cells = option.iga_n_cells or (2, 2, 2)

    mesh_fe = _structured_box_mesh(fe_lo, fe_hi, fe_n_cells, option.fe_element_type)
    box_iga = BoxSpline(lo=iga_lo, hi=iga_hi, n_cells=iga_n_cells, degree=option.iga_degree)

    spec = AssemblySpec(
        job_name="iga_fe_hybrid",
        parts={
            "fe_part": PartSpec(
                name="fe_part",
                kind="native_fe",
                mesh=mesh_fe,
                young_modulus=problem.E,
                poisson_ratio=problem.nu,
            ),
            "iga_part": PartSpec(
                name="iga_part",
                kind="iga",
                box=box_iga,
                young_modulus=problem.E,
                poisson_ratio=problem.nu,
            ),
        },
        instances=[InstanceSpec("fe_part", "1"), InstanceSpec("iga_part", "1")],
        use_pardiso=use_pardiso,
    )
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(spec)
    rt_fe = built.get_instance("fe_part_1")
    rt_iga = built.get_instance("iga_part_1")
    assert rt_fe.fe_handle is not None

    # Dirichlet BCs on whichever subdomain(s) touch problem.fixed_face
    fix_ax, fix_side = problem.fixed_face
    fix_coord = problem.lo[fix_ax] if fix_side == 0 else problem.hi[fix_ax]
    if _subdomain_touches_face(problem.fixed_face, ax, option.fe_side):
        for nid in rt_fe.node_ids:
            nd = built.structure.GetNode(nid)
            pos = (nd.X0, nd.Y0, nd.Z0)
            if abs(pos[fix_ax] - fix_coord) <= 1e-10:
                _fix_node_all_disp(KM, nd)
    if _subdomain_touches_face(problem.fixed_face, ax, iga_side_tag):
        for idx in box_iga.boundary_control_points(fix_ax, fix_side):
            _fix_node_all_disp(KM, built.structure.GetNode(rt_iga.node_ids[int(idx)]))

    # Traction loads on whichever subdomain(s) touch problem.load_face
    load_ax, load_side = problem.load_face
    if _subdomain_touches_face(problem.load_face, ax, option.fe_side):
        fe_loads = _integrate_fe_face_traction(
            mesh_fe,
            fe_lo,
            fe_hi,
            load_ax,
            load_side,
            problem.traction,
            rt_fe.fe_handle.mesh_to_kratos_node,
        )
        _attach_point_loads(SM, built, fe_loads, rt_fe.properties_map[0], scope="Loads.hyb_fe")
    if _subdomain_touches_face(problem.load_face, ax, iga_side_tag):
        iga_local_loads = _integrate_iga_face_traction(box_iga, load_ax, load_side, problem.traction)
        iga_kratos_loads = {rt_iga.node_ids[flat]: vec for flat, vec in iga_local_loads.items()}
        _attach_point_loads(SM, built, iga_kratos_loads, rt_iga.properties_map[0], scope="Loads.hyb_iga")

    # Interface Tie on x[ax] == x_split
    in_axes = [a for a in (0, 1, 2) if a != ax]
    n_div = (
        max(4, 2 * max(int(fe_n_cells[in_axes[0]]), int(iga_n_cells[in_axes[0]]))),
        max(4, 2 * max(int(fe_n_cells[in_axes[1]]), int(iga_n_cells[in_axes[1]]))),
    )
    surf_fe = SurfaceMesh.from_fe(
        "fe_iface",
        mesh_fe,
        rt_fe.fe_handle.mesh_to_kratos_node,
        lambda c: abs(float(c[ax]) - x_split) <= 1e-10,
    )
    surf_iga = SurfaceMesh.from_iga_face(
        "iga_iface",
        box_iga,
        axis=ax,  # type: ignore[arg-type]
        side=iga_iface_side,
        dof_node_ids=rt_iga.node_ids,
        n_div=n_div,
    )

    t1 = time.perf_counter()
    ties: list[TieSpec] | None = None
    # "auto" resolves to the Kratos-native exact path (mortar_l2): penalty is
    # only approximate (1/beta) and barycentric is O(h)-consistent at the
    # interface, while mortar reproduces linears to machine precision on both
    # matching and non-matching interfaces (S34). mpc_matching stays available
    # as an explicit option for coincident FE-FE interfaces.
    coupling = "mortar_l2" if option.coupling_kind == "auto" else option.coupling_kind
    if coupling == "barycentric":
        _apply_iga_barycentric_tie(KM, built, surf_iga, surf_fe)
        res_h = solve_hybrid_iga_fe(KM, SM, built, ties=None, use_pardiso=use_pardiso)
    else:
        if coupling == "mortar_l2":
            # Provide master Kratos control-point node IDs matching eval_basis_fn columns.
            surf_iga.node_ids = surf_iga.dof_node_ids
        tie_spec = TieSpec(
            name="hyb_iface_tie",
            master=surf_iga,
            slave=surf_fe,
            method=coupling,  # type: ignore[arg-type]
            E_max=problem.E,
            penalty_beta=option.penalty_beta,
        )
        ties = [tie_spec]
        res_h = solve_hybrid_iga_fe(KM, SM, built, ties=ties, use_pardiso=use_pardiso)
    t2 = time.perf_counter()

    u_fe = np.array(
        [
            list(built.structure.GetNode(rt_fe.fe_handle.mesh_to_kratos_node[m_idx]).GetSolutionStepValue(KM.DISPLACEMENT, 0))
            for m_idx in range(len(mesh_fe.nodes))
        ],
        dtype=np.float64,
    )
    u_iga = np.array(
        [list(built.structure.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0)) for nid in rt_iga.node_ids],
        dtype=np.float64,
    )

    def _eval(pts: np.ndarray) -> np.ndarray:
        pts_arr = np.atleast_2d(np.asarray(pts, dtype=np.float64))
        out = np.empty((len(pts_arr), 3), dtype=np.float64)
        coord = pts_arr[:, ax]
        fe_mask = (coord <= x_split + 1e-12) if option.fe_side == "lo" else (coord >= x_split - 1e-12)
        if np.any(fe_mask):
            out[fe_mask] = _eval_structured_fe_displacement(
                mesh_fe,
                option.fe_element_type,
                fe_lo,
                fe_hi,
                fe_n_cells,
                u_fe,
                pts_arr[fe_mask],
            )
        iga_mask = ~fe_mask
        if np.any(iga_mask):
            out[iga_mask] = box_iga.evaluate(u_iga, pts_arr[iga_mask])
        return out

    tip_disp = _eval(problem.tip_point().reshape(1, 3))[0]
    total_dofs = 3 * (len(rt_fe.node_ids) + len(rt_iga.node_ids))

    return MethodSolution(
        method="iga_fe_hybrid",
        total_dofs=total_dofs,
        active_dofs=total_dofs,
        tip_displacement=tip_disp,
        max_displacement=float(res_h.max_displacement),
        converged=True,
        timings={"setup_s": t1 - t0, "solve_s": t2 - t1, "total_s": t2 - t0},
        native={
            "built": built,
            "hybrid_result": res_h,
            "mesh_fe": mesh_fe,
            "box_iga": box_iga,
            "x_split": x_split,
        },
        _eval_fn=_eval,
    )


def solve_with_discretization_method(
    problem: BoxMethodProblem,
    option: DiscretizationOption,
    use_pardiso: bool = False,
) -> MethodSolution:
    """Solve a :class:`BoxMethodProblem` using the selected :data:`DiscretizationOption`.

    Parameters
    ----------
    problem:
        Box geometry, material properties, Dirichlet face, and traction load specification.
    option:
        Validated discretization method option (any of the 5 supported methods).
    use_pardiso:
        Prefer MKL Pardiso solver when available.

    Returns
    -------
    MethodSolution
        Unified displacement evaluator, tip displacement, DOF counts, timings, and native handles.
    """
    if isinstance(option, FullIgaOption):
        return _solve_full_iga(problem, option, use_pardiso=use_pardiso)
    if isinstance(option, ThbAdaptiveIgaOption):
        return _solve_thb_adaptive(problem, option, use_pardiso=use_pardiso)
    if isinstance(option, NativeFeOption):
        return _solve_native_fe(problem, option, use_pardiso=use_pardiso)
    if isinstance(option, IgaFeHybridOption):
        return _solve_iga_fe_hybrid(problem, option, use_pardiso=use_pardiso)
    raise TypeError(f"Unsupported discretization option type: {type(option).__name__}")

def compare_method_solutions(
    ref: MethodSolution,
    other: MethodSolution,
    points: np.ndarray,
) -> dict[str, float]:
    """Compare two :class:`MethodSolution` fields over probe ``points``.

    Returns tip/max/total displacement differences plus the relative L2
    displacement difference over ``points``. Pure function (no Kratos).
    """
    pts = np.asarray(points, dtype=np.float64)
    u_ref = np.asarray(ref.displacement_at(pts), dtype=np.float64).reshape(-1, 3)
    u_o = np.asarray(other.displacement_at(pts), dtype=np.float64).reshape(-1, 3)
    diff = u_o - u_ref
    denom = float(np.sqrt(np.sum(u_ref**2)))
    return {
        "tip_abs": float(np.linalg.norm(other.tip_displacement - ref.tip_displacement)),
        "tip_rel": float(np.linalg.norm(other.tip_displacement - ref.tip_displacement)
                         / max(1e-300, float(np.linalg.norm(ref.tip_displacement)))),
        "max_abs": float(abs(other.max_displacement - ref.max_displacement)),
        "l2_abs": float(np.sqrt(np.sum(diff**2))),
        "l2_rel": float(np.sqrt(np.sum(diff**2)) / denom) if denom > 0.0 else 0.0,
    }
