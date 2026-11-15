"""Body-fit structured all-hex to NURBS volume converter (T-P3.4 / v2.1 §5.10.4).

Reconstructs 3D structured ``(i, j, k)`` grids from conforming ``hex8``
:class:`~wht_iga.meshio_ext.model.MeshModel` meshes and converts each block
into a clamped B-Spline / NURBS volume patch (:class:`BodyFitSplinePatch`).

Unstructured all-hex meshes (T-junctions / hanging nodes, irregular interior
vertex valence != 8, irregular interior edge valence != 4, pinch-outs, holes,
or mixed element types) are detected explicitly and reported with
``ValidationIssue(code="UNSTRUCTURED_HEX_UNSUPPORTED", ...)`` so callers can
route them to the immersed IGA pipeline instead.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline

from wht_iga.grid.knots import to_kratos, uniform_clamped
from wht_iga.meshio_ext.model import HEX8_FACES, MeshModel
from wht_iga.solver.validate import ValidationIssue

# Canonical hex8 local node coordinates in (u, v, w) in {0, 1}^3.
# Matches HEX8_FACES winding in model.py.
_CORNER_COORDS: tuple[tuple[int, int, int], ...] = (
    (0, 0, 0),  # 0
    (1, 0, 0),  # 1
    (1, 1, 0),  # 2
    (0, 1, 0),  # 3
    (0, 0, 1),  # 4
    (1, 0, 1),  # 5
    (1, 1, 1),  # 6
    (0, 1, 1),  # 7
)
_COORD_TO_LOCAL: dict[tuple[int, int, int], int] = {
    coord: idx for idx, coord in enumerate(_CORNER_COORDS)
}
_HEX8_EDGES: tuple[tuple[int, int], ...] = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
)


class BodyFitValidationError(ValueError):
    """Raised when a mesh cannot be represented as body-fit structured NURBS patches."""

    def __init__(self, issues: list[ValidationIssue]) -> None:
        self.issues = list(issues)
        summary = "; ".join(f"[{iss.code}] {iss.message}" for iss in self.issues)
        super().__init__(summary)


@dataclass
class PatchInterface:
    """Conforming shared boundary interface between two structured NURBS patches.

    Parameters
    ----------
    patch_a, patch_b:
        Zero-based patch indices in :attr:`BodyFitConversionResult.patches`.
    shared_nodes:
        Sorted global 0-based mesh node indices shared across the interface.
    """

    patch_a: int
    patch_b: int
    shared_nodes: list[int] = field(default_factory=list)


@dataclass
class BodyFitSplinePatch:
    """Single 3D structured NURBS / B-Spline volume patch extracted from a hex8 block.

    Control points follow Fortran linear ordering ``I = a + nu * (b + nv * c)``
    to match :class:`~wht_iga.reference.box.BoxSpline` and Kratos
    ``NurbsVolumeGeometry``.
    """

    patch_id: int
    part_id: int
    part_name: str
    degree: tuple[int, int, int]
    knots: tuple[np.ndarray, np.ndarray, np.ndarray]
    control_points: np.ndarray  # shape (nu, nv, nw, 3)
    weights: np.ndarray  # shape (nu, nv, nw)
    node_grid: np.ndarray  # shape (nx+1, ny+1, nz+1) int64 node indices
    elem_grid: np.ndarray  # shape (nx, ny, nz) int64 element indices

    @property
    def shape(self) -> tuple[int, int, int]:
        """Number of control points ``(nu, nv, nw)`` per parameter direction."""
        return (
            int(self.control_points.shape[0]),
            int(self.control_points.shape[1]),
            int(self.control_points.shape[2]),
        )

    @property
    def n_cells(self) -> tuple[int, int, int]:
        """Number of hex8 cells ``(nx, ny, nz)`` per parameter direction."""
        return (
            int(self.elem_grid.shape[0]),
            int(self.elem_grid.shape[1]),
            int(self.elem_grid.shape[2]),
        )

    def control_points_flat(self) -> np.ndarray:
        """Return control points of shape ``(nu * nv * nw, 3)`` in Fortran order."""
        return self.control_points.reshape(-1, 3, order="F").copy()

    def weights_flat(self) -> np.ndarray:
        """Return NURBS weights of shape ``(nu * nv * nw,)`` in Fortran order."""
        return self.weights.reshape(-1, order="F").copy()

    def kratos_knots(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Convert clamped knot vectors to Kratos reduced-end convention."""
        return (
            to_kratos(self.knots[0], self.degree[0]),
            to_kratos(self.knots[1], self.degree[1]),
            to_kratos(self.knots[2], self.degree[2]),
        )

    def evaluate(
        self, u: np.ndarray | float, v: np.ndarray | float, w: np.ndarray | float
    ) -> np.ndarray:
        """Evaluate the NURBS volume map ``S(u, v, w)`` at parameter coordinates.

        Parameters
        ----------
        u, v, w:
            Scalar or 1D array parameter coordinates in ``[0, 1]``. When arrays
            are passed, they must have matching shape ``(M,)`` and evaluate
            pointwise into ``(M, 3)``. Scalars return shape ``(3,)``.
        """
        scalar_mode = np.ndim(u) == 0 and np.ndim(v) == 0 and np.ndim(w) == 0
        uu = np.atleast_1d(np.asarray(u, dtype=np.float64))
        vv = np.atleast_1d(np.asarray(v, dtype=np.float64))
        ww = np.atleast_1d(np.asarray(w, dtype=np.float64))
        if not (uu.shape == vv.shape == ww.shape and uu.ndim == 1):
            raise ValueError("u, v, w must all be scalars or equal-length 1D arrays")

        bu = BSpline.design_matrix(uu, self.knots[0], self.degree[0]).toarray()
        bv = BSpline.design_matrix(vv, self.knots[1], self.degree[1]).toarray()
        bw = BSpline.design_matrix(ww, self.knots[2], self.degree[2]).toarray()

        wp = self.weights[..., None] * self.control_points
        num = np.einsum("ma,mb,mc,abcd->md", bu, bv, bw, wp, optimize=True)
        den = np.einsum("ma,mb,mc,abc->m", bu, bv, bw, self.weights, optimize=True)
        pts = num / den[:, None]
        return pts[0] if scalar_mode else pts


@dataclass
class BodyFitConversionResult:
    """Result returned by :func:`mesh_to_bodyfit_patches`.

    Parameters
    ----------
    patches:
        Extracted :class:`BodyFitSplinePatch` objects (empty if invalid).
    interfaces:
        Conforming shared-node interfaces between adjacent patches.
    issues:
        Validation diagnostics (``UNSTRUCTURED_HEX_UNSUPPORTED`` when invalid).
    """

    patches: list[BodyFitSplinePatch] = field(default_factory=list)
    interfaces: list[PatchInterface] = field(default_factory=list)
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """``True`` when conversion succeeded with zero error-severity issues."""
        return not any(iss.severity == "error" for iss in self.issues)


def _issue(code: str, message: str, entity: str = "mesh") -> ValidationIssue:
    return ValidationIssue(
        severity="error",
        code=code,
        message=(
            f"{message} Use the immersed IGA pipeline instead of body-fit conversion."
        ),
        entity=entity,
    )


def _hex_corner_jacobian_sign(nodes: np.ndarray, cell: np.ndarray) -> float:
    """Return signed triple product at canonical local corner 0 of a hex8."""
    p0 = nodes[int(cell[0])]
    e1 = nodes[int(cell[1])] - p0
    e2 = nodes[int(cell[3])] - p0
    e3 = nodes[int(cell[4])] - p0
    return float(np.dot(e1, np.cross(e2, e3)))


def _canonical_quad(q: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(sorted(q))  # type: ignore[return-value]


def _check_t_junctions(
    nodes: np.ndarray,
    conn: np.ndarray,
    boundary_faces: list[tuple[int, int, int, int]],
    elem_of_node: list[set[int]],
    tol_rel: float = 1e-6,
) -> ValidationIssue | None:
    """Detect hanging/T-junction nodes lying strictly inside boundary edges or faces."""
    if not boundary_faces:
        return None

    coords_arr = np.asarray(nodes, dtype=np.float64)
    # Build face bounding boxes and corner sets.
    face_verts = np.array(boundary_faces, dtype=np.int64)
    face_pts = coords_arr[face_verts]  # (F, 4, 3)
    f_lo = face_pts.min(axis=1)
    f_hi = face_pts.max(axis=1)
    f_diag = np.linalg.norm(f_hi - f_lo, axis=1)
    max_diag = float(np.max(f_diag)) if len(f_diag) else 1.0
    pad = np.maximum(f_diag * tol_rel * 100.0, 1e-9 * max(max_diag, 1.0))

    # Candidate query: for each node appearing in boundary faces (or any node),
    # check if it sits inside the interior of another element's boundary face or edge.
    n_nodes = len(coords_arr)
    for nid in range(n_nodes):
        if not elem_of_node[nid]:
            continue
        pt = coords_arr[nid]
        mask = (
            (pt[0] >= f_lo[:, 0] - pad)
            & (pt[0] <= f_hi[:, 0] + pad)
            & (pt[1] >= f_lo[:, 1] - pad)
            & (pt[1] <= f_hi[:, 1] + pad)
            & (pt[2] >= f_lo[:, 2] - pad)
            & (pt[2] <= f_hi[:, 2] + pad)
        )
        cand_indices = np.nonzero(mask)[0]
        for f_idx in cand_indices:
            quad = boundary_faces[int(f_idx)]
            if nid in quad:
                continue
            q_pts = face_pts[int(f_idx)]  # (4, 3)
            tol = max(float(f_diag[int(f_idx)]) * tol_rel, 1e-10)

            # 1) Check strictly inside one of the 4 quad edges.
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
                va, vb = q_pts[a], q_pts[b]
                seg = vb - va
                seg_len2 = float(np.dot(seg, seg))
                if seg_len2 <= 0.0:
                    continue
                t = float(np.dot(pt - va, seg)) / seg_len2
                if 1e-4 < t < 1.0 - 1e-4:
                    dist = float(np.linalg.norm(pt - (va + t * seg)))
                    if dist <= tol:
                        return _issue(
                            "UNSTRUCTURED_HEX_UNSUPPORTED",
                            f"T-junction / hanging node {nid} lies on edge "
                            f"({quad[a]}, {quad[b]}) of hex face {quad}.",
                            entity=f"node:{nid}",
                        )

            # 2) Check strictly inside quad interior via projection onto quad plane + 2 triangles.
            for tri in ((0, 1, 2), (0, 2, 3)):
                p0, p1, p2 = q_pts[tri[0]], q_pts[tri[1]], q_pts[tri[2]]
                v0 = p1 - p0
                v1 = p2 - p0
                normal = np.cross(v0, v1)
                n_norm = float(np.linalg.norm(normal))
                if n_norm <= 0.0:
                    continue
                unit_n = normal / n_norm
                plane_dist = abs(float(np.dot(pt - p0, unit_n)))
                if plane_dist > tol:
                    continue
                proj = pt - float(np.dot(pt - p0, unit_n)) * unit_n
                # Barycentric coordinates in triangle (p0, p1, p2)
                d00 = float(np.dot(v0, v0))
                d01 = float(np.dot(v0, v1))
                d11 = float(np.dot(v1, v1))
                v2 = proj - p0
                d20 = float(np.dot(v2, v0))
                d21 = float(np.dot(v2, v1))
                denom = d00 * d11 - d01 * d01
                if abs(denom) <= 1e-24:
                    continue
                beta = (d11 * d20 - d01 * d21) / denom
                gamma = (d00 * d21 - d01 * d20) / denom
                alpha = 1.0 - beta - gamma
                if alpha > -1e-4 and beta > -1e-4 and gamma > -1e-4:
                    # Ensure it isn't actually one of the vertices (already guarded by nid in quad
                    # plus check non-coincident vertex distance).
                    if min(np.linalg.norm(q_pts - pt, axis=1)) > tol:
                        return _issue(
                            "UNSTRUCTURED_HEX_UNSUPPORTED",
                            f"T-junction / hanging node {nid} lies inside quad face {quad}.",
                            entity=f"node:{nid}",
                        )
    return None


def _propagate_component_grid(
    comp_elems: list[int],
    conn: np.ndarray,
    nodes: np.ndarray,
    face_map: dict[tuple[int, int, int, int], list[tuple[int, int]]],
    comp_label: str,
) -> tuple[np.ndarray, np.ndarray] | ValidationIssue:
    """Flood-fill a face-connected hex8 component into structured ``(node_grid, elem_grid)``."""
    comp_set = set(comp_elems)
    seed = min(comp_elems)

    # Align seed corner assignment so local (i, j, k) is right-handed in XYZ.
    seed_conn = conn[seed]
    flip_k = _hex_corner_jacobian_sign(nodes, seed_conn) < 0.0

    elem_node_coords: dict[int, dict[int, tuple[int, int, int]]] = {}
    elem_cell_origin: dict[int, tuple[int, int, int]] = {}

    seed_map: dict[int, tuple[int, int, int]] = {}
    for local_idx, (cu, cv, cw) in enumerate(_CORNER_COORDS):
        nid = int(seed_conn[local_idx])
        k_coord = 1 - cw if flip_k else cw
        if nid in seed_map:
            return _issue(
                "UNSTRUCTURED_HEX_UNSUPPORTED",
                f"Degenerate hex8 element {seed} has repeated node {nid}.",
                entity=comp_label,
            )
        seed_map[nid] = (cu, cv, k_coord)

    elem_node_coords[seed] = seed_map
    elem_cell_origin[seed] = (0, 0, 0)

    queue: deque[int] = deque([seed])
    while queue:
        cur = queue.popleft()
        cur_map = elem_node_coords[cur]
        cur_org = elem_cell_origin[cur]
        cur_conn = conn[cur]

        for lf_idx, lf_nodes in enumerate(HEX8_FACES):
            face_nids = tuple(int(cur_conn[m]) for m in lf_nodes)
            f_key = _canonical_quad(face_nids)
            owners = [e for e, _ in face_map[f_key] if e in comp_set]
            if len(owners) > 2:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Non-manifold face {f_key} shared by {len(owners)} hex elements.",
                    entity=comp_label,
                )
            if len(owners) != 2:
                continue

            nb = owners[1] if owners[0] == cur else owners[0]
            nb_conn = conn[nb]
            # Find which local face of `nb` matches `f_key`.
            nb_matched_locals: tuple[int, int, int, int] | None = None
            for nb_lf in HEX8_FACES:
                if _canonical_quad(tuple(int(nb_conn[m]) for m in nb_lf)) == f_key:  # type: ignore[arg-type]
                    nb_matched_locals = nb_lf  # type: ignore[assignment]
                    break
            if nb_matched_locals is None:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Neighbor element {nb} does not share quad face {f_key} cleanly.",
                    entity=comp_label,
                )

            # Determine constant axis (axis, shared_val, nb_step) of the shared face.
            face_coords = [cur_map[nid] for nid in face_nids]
            const_axes = [
                ax
                for ax in range(3)
                if len({fc[ax] for fc in face_coords}) == 1
            ]
            if len(const_axes) != 1:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Irregular grid coordinates on shared face {f_key}.",
                    entity=comp_label,
                )
            ax = const_axes[0]
            shared_plane_val = face_coords[0][ax]
            if shared_plane_val == cur_org[ax] + 1:
                step_sign = 1
                nb_org_ax = cur_org[ax] + 1
            elif shared_plane_val == cur_org[ax]:
                step_sign = -1
                nb_org_ax = cur_org[ax] - 1
            else:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Face {f_key} coordinates inconsistent with cell origin {cur_org}.",
                    entity=comp_label,
                )

            # Reconstruct the remaining 4 vertices of `nb` opposite its matched face:
            # Each corner `m` in `nb` has a unique cube-opposite corner across the
            # local axis perpendicular to `nb_matched_locals`.
            nb_face_set = set(nb_matched_locals)
            nb_face_local_coords = [_CORNER_COORDS[m] for m in nb_matched_locals]
            nb_const_local_axes = [
                lax
                for lax in range(3)
                if len({lc[lax] for lc in nb_face_local_coords}) == 1
            ]
            if len(nb_const_local_axes) != 1:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Pinch-out / invalid face topology in element {nb}.",
                    entity=comp_label,
                )
            nb_perp_ax = nb_const_local_axes[0]

            nb_map: dict[int, tuple[int, int, int]] = {}
            for m_face in nb_matched_locals:
                nid_face = int(nb_conn[m_face])
                if nid_face not in cur_map:
                    return _issue(
                        "UNSTRUCTURED_HEX_UNSUPPORTED",
                        f"Mismatched face vertices between elements {cur} and {nb}.",
                        entity=comp_label,
                    )
                gcoord = cur_map[nid_face]
                nb_map[nid_face] = gcoord

                lc = list(_CORNER_COORDS[m_face])
                lc[nb_perp_ax] = 1 - lc[nb_perp_ax]
                m_opp = _COORD_TO_LOCAL[tuple(lc)]  # type: ignore[arg-type]
                nid_opp = int(nb_conn[m_opp])
                opp_coord = list(gcoord)
                opp_coord[ax] += step_sign
                nb_map[nid_opp] = tuple(opp_coord)  # type: ignore[assignment]

            if len(nb_map) != 8 or len(set(nb_map.values())) != 8:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Pinch-out or degenerate node mapping on element {nb}.",
                    entity=comp_label,
                )

            nb_org = tuple(
                nb_org_ax if a == ax else min(c[a] for c in nb_map.values())
                for a in range(3)
            )
            if nb in elem_node_coords:
                # Topological loop closure check (catches irregular valence / O-grid / torus).
                if elem_node_coords[nb] != nb_map or elem_cell_origin[nb] != nb_org:
                    return _issue(
                        "UNSTRUCTURED_HEX_UNSUPPORTED",
                        f"Irregular hex valence / topological singularity around element {nb}.",
                        entity=comp_label,
                    )
            else:
                elem_node_coords[nb] = nb_map
                elem_cell_origin[nb] = nb_org  # type: ignore[assignment]
                queue.append(nb)

    if len(elem_cell_origin) != len(comp_elems):
        return _issue(
            "UNSTRUCTURED_HEX_UNSUPPORTED",
            f"{comp_label} flood-fill visited {len(elem_cell_origin)} of {len(comp_elems)} elements.",
            entity=comp_label,
        )

    origins = np.array(list(elem_cell_origin.values()), dtype=np.int64)
    lo = origins.min(axis=0)
    hi = origins.max(axis=0)
    nx, ny, nz = int(hi[0] - lo[0] + 1), int(hi[1] - lo[1] + 1), int(hi[2] - lo[2] + 1)

    if len(comp_elems) != nx * ny * nz:
        return _issue(
            "UNSTRUCTURED_HEX_UNSUPPORTED",
            f"Non-rectangular or incomplete hex topology in {comp_label}: "
            f"{len(comp_elems)} elements cannot fill a ({nx}, {ny}, {nz}) structured block.",
            entity=comp_label,
        )

    elem_grid = np.full((nx, ny, nz), -1, dtype=np.int64)
    node_grid = np.full((nx + 1, ny + 1, nz + 1), -1, dtype=np.int64)

    for e in comp_elems:
        ox, oy, oz = elem_cell_origin[e]
        ci, cj, ck = int(ox - lo[0]), int(oy - lo[1]), int(oz - lo[2])
        if elem_grid[ci, cj, ck] != -1:
            return _issue(
                "UNSTRUCTURED_HEX_UNSUPPORTED",
                f"Overlapping cells at structured grid index ({ci}, {cj}, {ck}) in {comp_label}.",
                entity=comp_label,
            )
        elem_grid[ci, cj, ck] = e

        for nid, (gx, gy, gz) in elem_node_coords[e].items():
            ni, nj, nk = int(gx - lo[0]), int(gy - lo[1]), int(gz - lo[2])
            prev = int(node_grid[ni, nj, nk])
            if prev == -1:
                node_grid[ni, nj, nk] = nid
            elif prev != nid:
                return _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Conflicting node IDs ({prev} vs {nid}) at lattice vertex ({ni}, {nj}, {nk}) "
                    f"in {comp_label}.",
                    entity=comp_label,
                )

    if np.any(elem_grid < 0) or np.any(node_grid < 0):
        return _issue(
            "UNSTRUCTURED_HEX_UNSUPPORTED",
            f"Holes or missing entries in structured grid for {comp_label}.",
            entity=comp_label,
        )

    expected_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    if len(np.unique(node_grid)) != expected_nodes:
        return _issue(
            "UNSTRUCTURED_HEX_UNSUPPORTED",
            f"Pinch-out / collapsed vertices in {comp_label}: expected {expected_nodes} unique "
            f"nodes, found {len(np.unique(node_grid))}.",
            entity=comp_label,
        )

    return node_grid, elem_grid


def _elevate_1d_piecewise_linear(pts: np.ndarray, degree: int) -> np.ndarray:
    """Degree-elevate 1D piecewise-linear control points along axis 0 from 1 -> ``degree``.

    ``pts`` has shape ``(n_cells + 1, ...)``; output has shape
    ``(n_cells * degree + 1, ...)`` matching ``uniform_clamped(0.0, 1.0, n_cells, degree, continuity=0)``.
    """
    n_cells = pts.shape[0] - 1
    if degree == 1:
        return pts.copy()
    rest_shape = pts.shape[1:]
    out = np.empty((n_cells * degree + 1, *rest_shape), dtype=np.float64)
    out[0] = pts[0]
    for e in range(n_cells):
        p0 = pts[e]
        p1 = pts[e + 1]
        for m in range(1, degree + 1):
            alpha = m / float(degree)
            out[e * degree + m] = (1.0 - alpha) * p0 + alpha * p1
    return out


def _fit_axis_spline(
    cp: np.ndarray, n_cells: int, degree: int, continuity: Literal["c0", "smooth"]
) -> tuple[np.ndarray, np.ndarray]:
    """Elevate or smoothly interpolate control net along axis 0.

    Returns ``(knots_1d, new_cp)``.
    """
    if degree == 1:
        return uniform_clamped(0.0, 1.0, n_cells, degree=1), cp.copy()

    if continuity == "c0" or n_cells < degree:
        knots = uniform_clamped(0.0, 1.0, n_cells, degree=degree, continuity=0)
        new_cp = _elevate_1d_piecewise_linear(cp, degree)
        return knots, new_cp

    # Smooth C^{p-1} B-spline interpolation using scipy.interpolate.make_interp_spline.
    t_params = np.linspace(0.0, 1.0, n_cells + 1)
    flat = cp.reshape(n_cells + 1, -1)
    spl = make_interp_spline(t_params, flat, k=degree)
    knots = np.asarray(spl.t, dtype=np.float64)
    new_cp = np.asarray(spl.c, dtype=np.float64).reshape(-1, *cp.shape[1:])
    return knots, new_cp


def _node_grid_to_nurbs(
    nodes: np.ndarray,
    node_grid: np.ndarray,
    degrees: tuple[int, int, int],
    continuity: Literal["c0", "smooth"],
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray, np.ndarray]:
    """Build clamped knot vectors, control points, and unit weights from ``node_grid``."""
    nx, ny, nz = (node_grid.shape[0] - 1, node_grid.shape[1] - 1, node_grid.shape[2] - 1)
    cp = nodes[node_grid].astype(np.float64)  # (nx+1, ny+1, nz+1, 3)

    ku, cp = _fit_axis_spline(cp, nx, degrees[0], continuity)
    cp = np.swapaxes(cp, 0, 1)
    kv, cp = _fit_axis_spline(cp, ny, degrees[1], continuity)
    cp = np.swapaxes(cp, 0, 1)
    cp = np.swapaxes(cp, 0, 2)
    kw, cp = _fit_axis_spline(cp, nz, degrees[2], continuity)
    cp = np.swapaxes(cp, 0, 2)

    weights = np.ones(cp.shape[:3], dtype=np.float64)
    return (ku, kv, kw), cp, weights


def mesh_to_bodyfit_patches(
    model: MeshModel,
    degree: int | tuple[int, int, int] = 1,
    continuity: Literal["c0", "smooth"] = "c0",
    strict: bool = True,
) -> BodyFitConversionResult:
    """Convert a structured all-hex8 :class:`MeshModel` into NURBS volume patches.

    Parameters
    ----------
    model:
        Input :class:`MeshModel`. Must contain only ``hex8`` blocks.
    degree:
        Polynomial degree per axis (scalar ``p`` or ``(pu, pv, pw)``, each >= 1).
    continuity:
        ``"c0"`` performs exact trilinear degree elevation via ``continuity=0``
        clamped knots (reproducing the hex8 mesh identically for all ``nx >= 1``).
        ``"smooth"`` uses ``scipy.interpolate.make_interp_spline`` when
        ``n_cells >= degree`` (falling back to ``"c0"`` when ``n_cells < degree``).
    strict:
        When ``True``, raise :class:`BodyFitValidationError` if any
        error-severity validation issue occurs.

    Returns
    -------
    BodyFitConversionResult
        Extracted patches, conforming interfaces, and validation issues.
    """
    if isinstance(degree, int):
        degrees = (int(degree), int(degree), int(degree))
    else:
        degrees = (int(degree[0]), int(degree[1]), int(degree[2]))
    if any(p < 1 for p in degrees):
        raise ValueError(f"All spline degrees must be >= 1, got {degrees}")
    if continuity not in ("c0", "smooth"):
        raise ValueError(f"continuity must be 'c0' or 'smooth', got {continuity!r}")

    issues: list[ValidationIssue] = []

    if not model.blocks:
        issues.append(
            _issue("UNSTRUCTURED_HEX_UNSUPPORTED", "MeshModel contains zero element blocks.")
        )
        if strict:
            raise BodyFitValidationError(issues)
        return BodyFitConversionResult(issues=issues)

    non_hex = [b for b in model.blocks if b.elem_type != "hex8"]
    if non_hex:
        bad_types = sorted({b.elem_type for b in non_hex})
        issues.append(
            _issue(
                "UNSTRUCTURED_HEX_UNSUPPORTED",
                f"Mixed or non-hex8 element types found ({bad_types}); body-fit conversion "
                "requires pure hex8 structured grids.",
            )
        )
        if strict:
            raise BodyFitValidationError(issues)
        return BodyFitConversionResult(issues=issues)

    # Flatten hex8 connectivity while keeping block order / part_id ownership.
    conn_list: list[np.ndarray] = []
    part_id_per_elem: list[int] = []
    part_elem_groups: dict[int, list[int]] = {}
    offset = 0
    for block in model.blocks:
        n_b = block.conn.shape[0]
        if n_b == 0:
            continue
        conn_list.append(block.conn.astype(np.int64, copy=False))
        for local_e in range(n_b):
            ge = offset + local_e
            part_id_per_elem.append(block.part_id)
            part_elem_groups.setdefault(block.part_id, []).append(ge)
        offset += n_b

    if not conn_list:
        issues.append(
            _issue("UNSTRUCTURED_HEX_UNSUPPORTED", "All hex8 element blocks are empty.")
        )
        if strict:
            raise BodyFitValidationError(issues)
        return BodyFitConversionResult(issues=issues)

    conn = np.vstack(conn_list)
    n_elems = int(conn.shape[0])
    nodes = np.asarray(model.nodes, dtype=np.float64)

    # Global node -> incident elements and global face adjacency.
    elem_of_node: list[set[int]] = [set() for _ in range(len(nodes))]
    face_map: dict[tuple[int, int, int, int], list[tuple[int, int]]] = {}
    for e in range(n_elems):
        cell = conn[e]
        if len({int(v) for v in cell}) != 8:
            issues.append(
                _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Element {e} has repeated node indices (pinch-out/degenerate).",
                    entity=f"elem:{e}",
                )
            )
            if strict:
                raise BodyFitValidationError(issues)
            return BodyFitConversionResult(issues=issues)
        for v in cell:
            elem_of_node[int(v)].add(e)
        for lf_idx, lf_nodes in enumerate(HEX8_FACES):
            fkey = _canonical_quad(tuple(int(cell[m]) for m in lf_nodes))  # type: ignore[arg-type]
            face_map.setdefault(fkey, []).append((e, lf_idx))

    # Non-manifold faces (> 2 elements sharing a face).
    for fkey, owners in face_map.items():
        if len(owners) > 2:
            issues.append(
                _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Non-manifold quad face {fkey} shared by {len(owners)} hex elements.",
                )
            )
            if strict:
                raise BodyFitValidationError(issues)
            return BodyFitConversionResult(issues=issues)

    # Boundary faces for T-junction detection + interior node valence check.
    boundary_faces: list[tuple[int, int, int, int]] = []
    boundary_nodes: set[int] = set()
    edge_incidents: dict[tuple[int, int], set[int]] = {}
    boundary_edges: set[tuple[int, int]] = set()

    for fkey, owners in face_map.items():
        if len(owners) == 1:
            e0, lf0 = owners[0]
            ordered = tuple(int(conn[e0, m]) for m in HEX8_FACES[lf0])
            boundary_faces.append(ordered)  # type: ignore[arg-type]
            boundary_nodes.update(ordered)
            for idx in range(4):
                ea, eb = ordered[idx], ordered[(idx + 1) % 4]
                boundary_edges.add((min(ea, eb), max(ea, eb)))

    for e in range(n_elems):
        cell = conn[e]
        for la, lb in _HEX8_EDGES:
            va, vb = int(cell[la]), int(cell[lb])
            ek = (min(va, vb), max(va, vb))
            edge_incidents.setdefault(ek, set()).add(e)

    # Check interior node valence == 8 and interior edge valence == 4 globally.
    for nid, inc in enumerate(elem_of_node):
        if inc and nid not in boundary_nodes and len(inc) != 8:
            issues.append(
                _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Irregular interior vertex valence at node {nid}: "
                    f"expected 8 incident hexes, found {len(inc)}.",
                    entity=f"node:{nid}",
                )
            )
            if strict:
                raise BodyFitValidationError(issues)
            return BodyFitConversionResult(issues=issues)

    for ek, inc in edge_incidents.items():
        if ek not in boundary_edges and len(inc) != 4:
            issues.append(
                _issue(
                    "UNSTRUCTURED_HEX_UNSUPPORTED",
                    f"Irregular interior edge valence at edge {ek}: "
                    f"expected 4 incident hexes, found {len(inc)}.",
                    entity=f"edge:{ek[0]}-{ek[1]}",
                )
            )
            if strict:
                raise BodyFitValidationError(issues)
            return BodyFitConversionResult(issues=issues)

    # T-junction / hanging node detection across boundary faces.
    tj_issue = _check_t_junctions(nodes, conn, boundary_faces, elem_of_node)
    if tj_issue is not None:
        issues.append(tj_issue)
        if strict:
            raise BodyFitValidationError(issues)
        return BodyFitConversionResult(issues=issues)

    # Decompose each part into face-connected components and reconstruct (i, j, k) blocks.
    patches: list[BodyFitSplinePatch] = []
    for part_id in sorted(part_elem_groups.keys()):
        p_elems = set(part_elem_groups[part_id])
        part_name = model.parts.get(part_id, f"part_{part_id}")
        visited: set[int] = set()
        components: list[list[int]] = []

        for e_start in sorted(p_elems):
            if e_start in visited:
                continue
            comp: list[int] = []
            q: deque[int] = deque([e_start])
            visited.add(e_start)
            while q:
                cur = q.popleft()
                comp.append(cur)
                cell = conn[cur]
                for lf_nodes in HEX8_FACES:
                    fkey = _canonical_quad(tuple(int(cell[m]) for m in lf_nodes))  # type: ignore[arg-type]
                    for nb, _ in face_map[fkey]:
                        if nb in p_elems and nb not in visited:
                            visited.add(nb)
                            q.append(nb)
            components.append(comp)

        # Within a single part, if two face-disconnected components share nodes,
        # they meet only along a corner/edge singularity rather than cleanly.
        if len(components) > 1:
            comp_nodes = [
                {int(v) for e in comp for v in conn[e]} for comp in components
            ]
            for a_idx in range(len(components)):
                for b_idx in range(a_idx + 1, len(components)):
                    shared = comp_nodes[a_idx] & comp_nodes[b_idx]
                    if shared:
                        issues.append(
                            _issue(
                                "UNSTRUCTURED_HEX_UNSUPPORTED",
                                f"Part {part_name!r} has disconnected hex components sharing "
                                f"nodes {sorted(shared)[:8]} without conforming faces.",
                                entity=part_name,
                            )
                        )
                        if strict:
                            raise BodyFitValidationError(issues)
                        return BodyFitConversionResult(issues=issues)

        for comp_idx, comp in enumerate(components):
            comp_label = (
                part_name if len(components) == 1 else f"{part_name}#block{comp_idx}"
            )
            res = _propagate_component_grid(comp, conn, nodes, face_map, comp_label)
            if isinstance(res, ValidationIssue):
                issues.append(res)
                if strict:
                    raise BodyFitValidationError(issues)
                return BodyFitConversionResult(issues=issues)

            node_grid, elem_grid = res
            knots, cp, weights = _node_grid_to_nurbs(nodes, node_grid, degrees, continuity)
            patches.append(
                BodyFitSplinePatch(
                    patch_id=len(patches),
                    part_id=part_id,
                    part_name=part_name,
                    degree=degrees,
                    knots=knots,
                    control_points=cp,
                    weights=weights,
                    node_grid=node_grid,
                    elem_grid=elem_grid,
                )
            )

    # Record shared-node conforming interfaces between patches.
    patch_node_sets = [set(int(v) for v in p.node_grid.ravel()) for p in patches]
    interfaces: list[PatchInterface] = []
    for i in range(len(patches)):
        for j in range(i + 1, len(patches)):
            shared = sorted(patch_node_sets[i] & patch_node_sets[j])
            if shared:
                interfaces.append(
                    PatchInterface(patch_a=i, patch_b=j, shared_nodes=shared)
                )

    return BodyFitConversionResult(patches=patches, interfaces=interfaces, issues=issues)


def build_kratos_nurbs_volumes(
    KM: Any,
    model_part: Any,
    patches: list[BodyFitSplinePatch],
    first_node_id: int = 1,
    add_dofs: bool = True,
) -> list[tuple[Any, Any]]:
    """Build Kratos ``NurbsVolumeGeometry`` objects from :class:`BodyFitSplinePatch` list.

    Note: caller passes ``KM`` (``KratosMultiphysics``) so ``meshio_ext`` remains
    100% Kratos-free at import time.
    """
    built: list[tuple[Any, Any]] = []
    next_id = int(first_node_id)
    for patch in patches:
        pts = patch.control_points_flat()
        nodes_vec = KM.NodesVector()
        for idx, (x, y, z) in enumerate(pts):
            node = model_part.CreateNewNode(next_id + idx, float(x), float(y), float(z))
            nodes_vec.append(node)
        next_id += len(pts)

        if add_dofs:
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        ku, kv, kw = patch.kratos_knots()
        vol = KM.NurbsVolumeGeometry(
            nodes_vec,
            patch.degree[0],
            patch.degree[1],
            patch.degree[2],
            KM.Vector(ku.tolist()),
            KM.Vector(kv.tolist()),
            KM.Vector(kw.tolist()),
        )
        model_part.AddGeometry(vol)
        built.append((vol, nodes_vec))
    return built


__all__ = [
    "BodyFitConversionResult",
    "BodyFitSplinePatch",
    "BodyFitValidationError",
    "PatchInterface",
    "ValidationIssue",
    "build_kratos_nurbs_volumes",
    "mesh_to_bodyfit_patches",
]
