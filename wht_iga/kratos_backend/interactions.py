"""Coupling interactions, Tie, RBE2, RBE3, and constraint-graph validation.

Implements Track P3-C (`T-P3.2`) for multi-part Kratos assemblies:

* **Constraint Graph Validator & Compiler** (:func:`validate_constraint_graph`,
  :func:`compile_constraint_graph`): enforces S04 (``constant == 0.0``),
  double-slave detection, cycle detection, and optional automatic chain
  flattening (``u_s = A u_m, u_m = B u_r => u_s = (AB) u_r``).
* **Tie Constraints** (:class:`TieSpec`, :func:`build_tie`, :func:`apply_tie`):
  1:1 node-matching ``LinearMasterSlaveConstraint``, barycentric closest-point
  projection, exact Sutherland-Hodgman clipped surface Gauss ``L^2`` mortar
  projection (variational consistency for free-interface stress patches), and
  surface quadrature penalty spring ``K_tie`` for IGA-involved ties.
* **RBE2** (:class:`RBE2Spec`, :func:`build_rbe2`, :func:`apply_rbe2`): rigid
  kinematic coupling ``u_i = u_ref + theta_ref x (x_i - x_ref)``.
* **RBE3** (:class:`RBE3Spec`, :func:`build_rbe3`, :func:`apply_rbe3`):
  non-stiffening weighted least-squares distributing coupling and exact
  force/moment virtual-work distribution.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

import numpy as np
import scipy.sparse as sp

from wht_iga.meshio_ext.model import MeshModel, outer_faces
from wht_iga.selection.contact_check import _closest_point_on_triangle


DofKey = tuple[int, str]


@dataclass(frozen=True)
class LinearRelation:
    """Scalar linear master-slave constraint equation.

    Represents ``slave_dof = sum(coeff * master_dof) + constant``.
    Kratos S04 requires ``constant == 0.0``.
    """

    slave_dof: DofKey
    terms: dict[DofKey, float]
    constant: float = 0.0
    tag: str = ""


@dataclass
class ConstraintGraphReport:
    """Validation diagnostics produced by :func:`validate_constraint_graph`."""

    valid: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    slave_count: int = 0
    master_count: int = 0
    double_slaves: list[DofKey] = field(default_factory=list)
    nonzero_constants: list[DofKey] = field(default_factory=list)
    chain_edges: list[tuple[DofKey, DofKey]] = field(default_factory=list)
    cycles: list[list[DofKey]] = field(default_factory=list)


def validate_constraint_graph(
    relations: list[LinearRelation],
    allow_chains: bool = False,
    tol: float = 1e-12,
    raise_on_error: bool = False,
) -> ConstraintGraphReport:
    """Inspect a list of linear constraints for S04, double-slaves, chains, and cycles.

    Parameters
    ----------
    relations:
        Candidate scalar master-slave constraint equations.
    allow_chains:
        If ``False``, any ``slave -> master`` dependency chain is flagged as an
        error; if ``True``, chains are recorded as warnings so the compiler may
        flatten them.
    tol:
        Numerical tolerance for ``constant == 0.0`` and non-zero coefficient checks.
    raise_on_error:
        If ``True``, raise ``ValueError`` when any validation error is present.

    Returns
    -------
    ConstraintGraphReport
        Detailed validity and graph structure summary.
    """
    errors: list[str] = []
    warnings: list[str] = []
    double_slaves: list[DofKey] = []
    nonzero_constants: list[DofKey] = []
    chain_edges: list[tuple[DofKey, DofKey]] = []
    cycles: list[list[DofKey]] = []

    slave_map: dict[DofKey, LinearRelation] = {}
    all_masters: set[DofKey] = set()

    for rel in relations:
        if abs(float(rel.constant)) > tol:
            nonzero_constants.append(rel.slave_dof)
            errors.append(
                f"S04 violation on slave {rel.slave_dof}: constant={rel.constant:.3e} != 0"
            )
        if rel.slave_dof in slave_map:
            if rel.slave_dof not in double_slaves:
                double_slaves.append(rel.slave_dof)
            errors.append(f"Double-slave conflict on DOF {rel.slave_dof}")
        else:
            slave_map[rel.slave_dof] = rel
        for mdof, coeff in rel.terms.items():
            if abs(float(coeff)) > 1e-15:
                all_masters.add(mdof)

    # Detect cycles via 3-color DFS over slave -> active master edges.
    state: dict[DofKey, int] = {}  # 0=unvisited, 1=visiting, 2=done
    stack: list[DofKey] = []

    def _dfs(node: DofKey) -> None:
        state[node] = 1
        stack.append(node)
        rel = slave_map.get(node)
        if rel is not None:
            for mdof, coeff in rel.terms.items():
                if abs(float(coeff)) <= 1e-15:
                    continue
                if mdof in slave_map:
                    chain_edges.append((node, mdof))
                    st = state.get(mdof, 0)
                    if st == 1:
                        idx = stack.index(mdof)
                        cycles.append(stack[idx:] + [mdof])
                    elif st == 0:
                        _dfs(mdof)
        stack.pop()
        state[node] = 2

    for sdof in slave_map:
        if state.get(sdof, 0) == 0:
            _dfs(sdof)

    for cyc in cycles:
        errors.append(f"Cyclic constraint dependency detected: {cyc}")

    if chain_edges and not allow_chains:
        errors.append(
            f"Detected {len(chain_edges)} chained slave->master edge(s) (first: {chain_edges[0]})"
        )
    elif chain_edges:
        warnings.append(
            f"Detected {len(chain_edges)} chained slave->master edge(s) pending flattening"
        )

    pure_masters = all_masters - set(slave_map.keys())
    report = ConstraintGraphReport(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        slave_count=len(slave_map),
        master_count=len(pure_masters),
        double_slaves=double_slaves,
        nonzero_constants=nonzero_constants,
        chain_edges=chain_edges,
        cycles=cycles,
    )
    if raise_on_error and not report.valid:
        raise ValueError("; ".join(report.errors))
    return report


def compile_constraint_graph(
    relations: list[LinearRelation],
    chain_policy: Literal["reject", "flatten"] = "reject",
    prune_tol: float = 1e-14,
    s04_tol: float = 1e-12,
) -> list[LinearRelation]:
    """Validate and optionally flatten chained master-slave linear relations.

    Parameters
    ----------
    relations:
        Input linear constraint equations.
    chain_policy:
        ``"reject"`` forbids ``slave -> master`` chains; ``"flatten"`` performs
        exact algebraic substitution ``u_s = A u_m, u_m = B u_r => u_s = (AB) u_r``.
    prune_tol:
        Threshold below which cancelled master coefficients are dropped.
    s04_tol:
        Zero-constant verification tolerance.

    Returns
    -------
    list[LinearRelation]
        Chain-free, cycle-free, S04-compliant constraint relations.
    """
    report = validate_constraint_graph(
        relations,
        allow_chains=(chain_policy == "flatten"),
        tol=s04_tol,
        raise_on_error=True,
    )
    if chain_policy == "reject" or not report.chain_edges:
        return [
            LinearRelation(
                slave_dof=rel.slave_dof,
                terms={
                    k: float(v)
                    for k, v in rel.terms.items()
                    if abs(float(v)) > prune_tol
                },
                constant=0.0,
                tag=rel.tag,
            )
            for rel in relations
        ]

    slave_map: dict[DofKey, LinearRelation] = {r.slave_dof: r for r in relations}
    resolved: dict[DofKey, dict[DofKey, float]] = {}

    def _resolve(dof: DofKey) -> dict[DofKey, float]:
        if dof in resolved:
            return resolved[dof]
        rel = slave_map.get(dof)
        if rel is None:
            return {dof: 1.0}
        accum: dict[DofKey, float] = {}
        for mdof, coeff in rel.terms.items():
            c = float(coeff)
            if abs(c) <= prune_tol:
                continue
            sub = _resolve(mdof) if mdof in slave_map else {mdof: 1.0}
            for root_dof, root_w in sub.items():
                accum[root_dof] = accum.get(root_dof, 0.0) + c * root_w
        cleaned = {
            k: float(v) for k, v in accum.items() if abs(float(v)) > prune_tol
        }
        resolved[dof] = cleaned
        return cleaned

    compiled: list[LinearRelation] = []
    for rel in relations:
        terms = _resolve(rel.slave_dof)
        compiled.append(
            LinearRelation(
                slave_dof=rel.slave_dof,
                terms=terms,
                constant=0.0,
                tag=rel.tag,
            )
        )

    validate_constraint_graph(compiled, allow_chains=False, tol=s04_tol, raise_on_error=True)
    return compiled


def apply_relations_to_model_part(
    KM: Any,
    model_part: Any,
    relations: list[LinearRelation],
    id_table: Any | None = None,
    scope: str = "Couplings",
    start_constraint_id: int = 1,
) -> list[Any]:
    """Attach compiled :class:`LinearRelation`s as Kratos ``LinearMasterSlaveConstraint``s."""
    from wht_iga.kratos_backend.constraints import add_linear_tie

    compiled = compile_constraint_graph(relations, chain_policy="reject")
    if not compiled:
        return []
    if id_table is not None:
        entry = id_table.allocate_constraints(scope, len(compiled))
        cid = entry.start_id
    else:
        cid = int(start_constraint_id)

    created: list[Any] = []
    for rel in compiled:
        s_node = model_part.GetNode(int(rel.slave_dof[0]))
        s_var = getattr(KM, rel.slave_dof[1])
        s_dofs = [s_node.GetDof(s_var)]

        sorted_masters = sorted(rel.terms.items(), key=lambda kv: (kv[0][0], kv[0][1]))
        m_dofs = []
        weights = np.empty((1, len(sorted_masters)), dtype=np.float64)
        for j, (mdof, w) in enumerate(sorted_masters):
            m_node = model_part.GetNode(int(mdof[0]))
            m_var = getattr(KM, mdof[1])
            m_dofs.append(m_node.GetDof(m_var))
            weights[0, j] = float(w)
        msc = add_linear_tie(KM, model_part, cid, s_dofs, m_dofs, weights, constant=None)
        created.append(msc)
        cid += 1
    return created


def _unit_tri_gauss(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return all-positive quadrature points ``(Q, 3)`` barycentric & weights on unit tri."""
    # Area of reference unit triangle is 0.5; all weights are strictly positive and sum to 0.5.
    if order <= 1:
        bary = np.array([[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]], dtype=np.float64)
        wts = np.array([0.5], dtype=np.float64)
        return bary, wts
    if order == 2:
        pts = np.array(
            [
                [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
            ],
            dtype=np.float64,
        )
        wts = np.full(3, 1.0 / 6.0, dtype=np.float64)
        return pts, wts
    # Positive-weight Duffy tensor-product Gauss quadrature for order >= 3
    # (avoids the negative center weight of 4-pt degree-3 Dunavant rules so D^T W D is PSD).
    n = max(3, int(order))
    xi, wi = np.polynomial.legendre.leggauss(n)
    eta, wj = np.polynomial.legendre.leggauss(n)
    s_pts = 0.5 * (1.0 + xi)
    s_wts = 0.5 * wi
    t_pts = 0.5 * (1.0 + eta)
    t_wts = 0.5 * wj
    bary_list: list[list[float]] = []
    wt_list: list[float] = []
    for s, ws in zip(s_pts, s_wts):
        for t, wt in zip(t_pts, t_wts):
            l1 = 1.0 - s
            l2 = s * (1.0 - t)
            l3 = s * t
            bary_list.append([float(l1), float(l2), float(l3)])
            wt_list.append(float(ws * wt * s))
    return np.asarray(bary_list, dtype=np.float64), np.asarray(wt_list, dtype=np.float64)


def _clip_polygon_2d(subject: np.ndarray, clip_tri: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Sutherland-Hodgman 2D convex polygon clipping against CCW triangle ``clip_tri``."""
    area_clip = float(np.cross(clip_tri[1] - clip_tri[0], clip_tri[2] - clip_tri[0]))
    tri = clip_tri if area_clip >= 0.0 else clip_tri[[0, 2, 1]]

    poly = list(subject)
    for e_idx in range(3):
        a = tri[e_idx]
        b = tri[(e_idx + 1) % 3]
        edge = b - a
        if not poly:
            return np.empty((0, 2), dtype=np.float64)
        out: list[np.ndarray] = []

        def _inside(p: np.ndarray) -> float:
            return float(edge[0] * (p[1] - a[1]) - edge[1] * (p[0] - a[0]))

        prev = poly[-1]
        prev_d = _inside(prev)
        for cur in poly:
            cur_d = _inside(cur)
            if cur_d >= -tol:
                if prev_d < -tol:
                    denom = prev_d - cur_d
                    alpha = prev_d / denom if abs(denom) > 1e-15 else 0.0
                    out.append(prev + alpha * (cur - prev))
                out.append(cur)
            elif prev_d >= -tol:
                denom = prev_d - cur_d
                alpha = prev_d / denom if abs(denom) > 1e-15 else 0.0
                out.append(prev + alpha * (cur - prev))
            prev = cur
            prev_d = cur_d
        poly = out

    if len(poly) < 3:
        return np.empty((0, 2), dtype=np.float64)
    return np.asarray(poly, dtype=np.float64)


def _triangle_barycentric_in_plane(tri_pts: np.ndarray, q_pts: np.ndarray) -> np.ndarray:
    """Exact affine barycentric coordinates ``(Q, 3)`` on 3D planar triangle ``tri_pts (3,3)``."""
    v0 = tri_pts[1] - tri_pts[0]
    v1 = tri_pts[2] - tri_pts[0]
    d00 = float(np.dot(v0, v0))
    d01 = float(np.dot(v0, v1))
    d11 = float(np.dot(v1, v1))
    denom = d00 * d11 - d01 * d01
    if abs(denom) < 1e-24:
        return np.full((len(q_pts), 3), 1.0 / 3.0, dtype=np.float64)
    inv_det = 1.0 / denom
    v2 = q_pts - tri_pts[0][None, :]
    d20 = v2 @ v0
    d21 = v2 @ v1
    b1 = (d11 * d20 - d01 * d21) * inv_det
    b2 = (d00 * d21 - d01 * d20) * inv_det
    b0 = 1.0 - b1 - b2
    return np.column_stack([b0, b1, b2])


def _fe_face_shape_and_deriv(fkind: str, xi: float, eta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate FE surface shape functions ``N`` and derivatives ``dN/dxi, dN/deta``."""
    if fkind == "quad4":
        xm, xp = 1.0 - xi, 1.0 + xi
        em, ep = 1.0 - eta, 1.0 + eta
        N = 0.25 * np.array([xm * em, xp * em, xp * ep, xm * ep], dtype=np.float64)
        dxi = 0.25 * np.array([-em, em, ep, -ep], dtype=np.float64)
        deta = 0.25 * np.array([-xm, -xp, xp, xm], dtype=np.float64)
        return N, dxi, deta
    if fkind == "quad8":
        xm, xp = 1.0 - xi, 1.0 + xi
        em, ep = 1.0 - eta, 1.0 + eta
        N = np.array(
            [
                0.25 * xm * em * (-xi - eta - 1.0),
                0.25 * xp * em * (xi - eta - 1.0),
                0.25 * xp * ep * (xi + eta - 1.0),
                0.25 * xm * ep * (-xi + eta - 1.0),
                0.5 * (1.0 - xi * xi) * em,
                0.5 * xp * (1.0 - eta * eta),
                0.5 * (1.0 - xi * xi) * ep,
                0.5 * xm * (1.0 - eta * eta),
            ],
            dtype=np.float64,
        )
        dxi = np.array(
            [
                0.25 * em * (2.0 * xi + eta),
                0.25 * em * (2.0 * xi - eta),
                0.25 * ep * (2.0 * xi + eta),
                0.25 * ep * (2.0 * xi - eta),
                -xi * em,
                0.5 * (1.0 - eta * eta),
                -xi * ep,
                -0.5 * (1.0 - eta * eta),
            ],
            dtype=np.float64,
        )
        deta = np.array(
            [
                0.25 * xm * (xi + 2.0 * eta),
                0.25 * xp * (2.0 * eta - xi),
                0.25 * xp * (xi + 2.0 * eta),
                0.25 * xm * (2.0 * eta - xi),
                -0.5 * (1.0 - xi * xi),
                -eta * xp,
                0.5 * (1.0 - xi * xi),
                -eta * xm,
            ],
            dtype=np.float64,
        )
        return N, dxi, deta
    if fkind == "tri6":
        l1 = xi
        l2 = eta
        l0 = 1.0 - l1 - l2
        N = np.array(
            [
                l0 * (2.0 * l0 - 1.0),
                l1 * (2.0 * l1 - 1.0),
                l2 * (2.0 * l2 - 1.0),
                4.0 * l0 * l1,
                4.0 * l1 * l2,
                4.0 * l2 * l0,
            ],
            dtype=np.float64,
        )
        dxi = np.array(
            [
                -(4.0 * l0 - 1.0),
                4.0 * l1 - 1.0,
                0.0,
                4.0 * (l0 - l1),
                4.0 * l2,
                -4.0 * l2,
            ],
            dtype=np.float64,
        )
        deta = np.array(
            [
                -(4.0 * l0 - 1.0),
                0.0,
                4.0 * l2 - 1.0,
                -4.0 * l1,
                4.0 * l1,
                4.0 * (l0 - l2),
            ],
            dtype=np.float64,
        )
        return N, dxi, deta
    # tri3
    l0 = 1.0 - xi - eta
    N = np.array([l0, xi, eta], dtype=np.float64)
    dxi = np.array([-1.0, 1.0, 0.0], dtype=np.float64)
    deta = np.array([-1.0, 0.0, 1.0], dtype=np.float64)
    return N, dxi, deta


def _project_point_onto_fe_face(
    fkind: str,
    face_pts: np.ndarray,
    q: np.ndarray,
) -> tuple[np.ndarray, float, bool]:
    """Invert 3D physical point ``q`` onto parametric ``(xi, eta)`` of an FE face.

    Returns ``(N_vals, distance, inside_domain)``.
    """
    if fkind == "tri3":
        cp = _closest_point_on_triangle(q, face_pts[0], face_pts[1], face_pts[2])
        b = _triangle_barycentric_in_plane(face_pts, cp[None, :])[0]
        dist = float(np.linalg.norm(q - cp))
        inside = bool(np.all(b >= -1e-7) and np.all(b <= 1.0 + 1e-7))
        return b, dist, inside

    if fkind == "tri6":
        b0 = _triangle_barycentric_in_plane(face_pts[:3], q[None, :])[0]
        xi, eta = float(b0[1]), float(b0[2])
    else:
        center = face_pts[:4].mean(axis=0)
        v1 = 0.25 * (-face_pts[0] + face_pts[1] + face_pts[2] - face_pts[3])
        v2 = 0.25 * (-face_pts[0] - face_pts[1] + face_pts[2] + face_pts[3])
        n11 = max(float(np.dot(v1, v1)), 1e-24)
        n22 = max(float(np.dot(v2, v2)), 1e-24)
        xi = float(np.clip(np.dot(q - center, v1) / n11, -1.2, 1.2))
        eta = float(np.clip(np.dot(q - center, v2) / n22, -1.2, 1.2))

    for _ in range(8):
        N, dxi, deta = _fe_face_shape_and_deriv(fkind, xi, eta)
        x_cur = N @ face_pts
        r = q - x_cur
        J0 = dxi @ face_pts
        J1 = deta @ face_pts
        G = np.array(
            [
                [float(np.dot(J0, J0)), float(np.dot(J0, J1))],
                [float(np.dot(J0, J1)), float(np.dot(J1, J1))],
            ],
            dtype=np.float64,
        )
        rhs = np.array([float(np.dot(J0, r)), float(np.dot(J1, r))], dtype=np.float64)
        det = G[0, 0] * G[1, 1] - G[0, 1] * G[1, 0]
        if abs(det) < 1e-24:
            break
        d_xi = (G[1, 1] * rhs[0] - G[0, 1] * rhs[1]) / det
        d_eta = (-G[1, 0] * rhs[0] + G[0, 0] * rhs[1]) / det
        xi += float(d_xi)
        eta += float(d_eta)
        if abs(d_xi) + abs(d_eta) < 1e-13:
            break

    if fkind == "tri6":
        l0 = 1.0 - xi - eta
        inside = bool(l0 >= -1e-7 and xi >= -1e-7 and eta >= -1e-7)
        if not inside:
            xi_c = max(0.0, min(1.0, xi))
            eta_c = max(0.0, min(1.0 - xi_c, eta))
            xi, eta = xi_c, eta_c
    else:
        inside = bool(abs(xi) <= 1.0 + 1e-7 and abs(eta) <= 1.0 + 1e-7)
        xi = float(np.clip(xi, -1.0, 1.0))
        eta = float(np.clip(eta, -1.0, 1.0))

    N, _, _ = _fe_face_shape_and_deriv(fkind, xi, eta)
    dist = float(np.linalg.norm(q - N @ face_pts))
    return N, dist, inside


def _build_fe_surface_evaluator(
    mesh: MeshModel,
    coords_all: np.ndarray,
    selected_tri_sets: set[frozenset[int]],
    remap: dict[int, int],
    n_surf_verts: int,
) -> Callable[[np.ndarray], sp.csr_matrix] | None:
    """Build exact isoparametric FE surface basis evaluator (Q4/Q8/P1/P2)."""
    block_faces: dict[str, list[tuple[str, tuple[int, ...]]]] = {
        "hex8": [
            ("quad4", (0, 3, 2, 1)),
            ("quad4", (4, 5, 6, 7)),
            ("quad4", (0, 1, 5, 4)),
            ("quad4", (1, 2, 6, 5)),
            ("quad4", (2, 3, 7, 6)),
            ("quad4", (3, 0, 4, 7)),
        ],
        "hex20": [
            ("quad8", (0, 3, 2, 1, 11, 10, 9, 8)),
            ("quad8", (4, 5, 6, 7, 12, 13, 14, 15)),
            ("quad8", (0, 1, 5, 4, 8, 17, 12, 16)),
            ("quad8", (1, 2, 6, 5, 9, 18, 13, 17)),
            ("quad8", (2, 3, 7, 6, 10, 19, 14, 18)),
            ("quad8", (3, 0, 4, 7, 11, 16, 15, 19)),
        ],
        "tet4": [
            ("tri3", (0, 2, 1)),
            ("tri3", (0, 1, 3)),
            ("tri3", (1, 2, 3)),
            ("tri3", (2, 0, 3)),
        ],
        "tet10": [
            ("tri6", (0, 2, 1, 6, 5, 4)),
            ("tri6", (0, 1, 3, 4, 8, 7)),
            ("tri6", (1, 2, 3, 5, 9, 8)),
            ("tri6", (2, 0, 3, 6, 7, 9)),
        ],
        "wedge6": [
            ("tri3", (0, 2, 1)),
            ("tri3", (3, 4, 5)),
            ("quad4", (0, 1, 4, 3)),
            ("quad4", (1, 2, 5, 4)),
            ("quad4", (2, 0, 3, 5)),
        ],
    }

    face_counts: dict[tuple[int, ...], int] = {}
    face_records: list[tuple[tuple[int, ...], str, tuple[int, ...]]] = []
    for block in mesh.blocks:
        templates = block_faces.get(block.elem_type)
        if templates is None:
            continue
        for row in block.conn:
            for fkind, lidx in templates:
                g_nodes = tuple(int(row[k]) for k in lidx)
                n_corner = 3 if fkind.startswith("tri") else 4
                ckey = tuple(sorted(g_nodes[:n_corner]))
                face_counts[ckey] = face_counts.get(ckey, 0) + 1
                face_records.append((ckey, fkind, g_nodes))

    selected_faces: list[tuple[str, list[int], np.ndarray, np.ndarray, np.ndarray]] = []
    for ckey, fkind, fv in face_records:
        if face_counts.get(ckey, 0) != 1:
            continue
        if fkind == "tri3":
            matched = frozenset(fv) in selected_tri_sets
        elif fkind == "quad4":
            matched = any(
                frozenset(tri) in selected_tri_sets
                for tri in (
                    (fv[0], fv[1], fv[2]),
                    (fv[0], fv[2], fv[3]),
                    (fv[1], fv[2], fv[3]),
                    (fv[0], fv[1], fv[3]),
                )
            )
        elif fkind == "quad8":
            matched = any(
                frozenset((fv[k], fv[4 + k], fv[4 + ((k - 1) % 4)])) in selected_tri_sets
                for k in range(4)
            )
        else:  # tri6
            matched = any(
                frozenset((fv[k], fv[3 + k], fv[3 + ((k - 1) % 3)])) in selected_tri_sets
                for k in range(3)
            )
        if not matched:
            continue
        if any(v not in remap for v in fv):
            continue
        pts = np.asarray(coords_all[list(fv)], dtype=np.float64)
        lcols = [remap[v] for v in fv]
        diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
        pad = max(0.05 * diag, 1e-6)
        selected_faces.append((fkind, lcols, pts, pts.min(axis=0) - pad, pts.max(axis=0) + pad))

    if not selected_faces:
        return None

    bmins = np.vstack([f[3] for f in selected_faces])
    bmaxs = np.vstack([f[4] for f in selected_faces])
    centers = np.vstack([f[2].mean(axis=0) for f in selected_faces])

    def _eval(query_pts: np.ndarray) -> sp.csr_matrix:
        q_arr = np.asarray(query_pts, dtype=np.float64)
        n_q = len(q_arr)
        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []
        for qi, q in enumerate(q_arr):
            cand = np.where(np.all(q >= bmins, axis=1) & np.all(q <= bmaxs, axis=1))[0]
            if len(cand) == 0:
                cd = np.linalg.norm(centers - q[None, :], axis=1)
                cand = np.argsort(cd)[: min(8, len(selected_faces))]
            best_d = float("inf")
            best_n: np.ndarray | None = None
            best_cols: list[int] = []
            for fi in cand:
                fkind, lcols, pts, _, _ = selected_faces[int(fi)]
                N_f, dist, inside = _project_point_onto_fe_face(fkind, pts, q)
                if dist < best_d:
                    best_d = dist
                    best_n = N_f
                    best_cols = lcols
                    if inside and dist < 1e-11:
                        break
            if best_d > 1e-6 and len(cand) < len(selected_faces):
                cd = np.linalg.norm(centers - q[None, :], axis=1)
                for fi in np.argsort(cd)[: min(16, len(selected_faces))]:
                    fkind, lcols, pts, _, _ = selected_faces[int(fi)]
                    N_f, dist, _ = _project_point_onto_fe_face(fkind, pts, q)
                    if dist < best_d:
                        best_d = dist
                        best_n = N_f
                        best_cols = lcols
            assert best_n is not None
            for c_idx, w in zip(best_cols, best_n):
                if abs(float(w)) > 1e-15:
                    rows.append(qi)
                    cols.append(c_idx)
                    vals.append(float(w))
        return sp.csr_matrix((vals, (rows, cols)), shape=(n_q, n_surf_verts))

    return _eval


@dataclass
class SurfaceMesh:
    """Surface discretization descriptor supporting FE meshes and IGA/THB faces.

    Parameters
    ----------
    name:
        Surface identifier.
    kind:
        ``"fe"``, ``"iga"``, or ``"thb"``.
    coords:
        ``(N_vert, 3)`` surface triangulation vertex physical coordinates.
    triangles:
        ``(M_tri, 3)`` 0-based vertex indices into ``coords``.
    node_ids:
        ``(N_vert,)`` Kratos node IDs for ``kind == "fe"`` vertices.
    dof_node_ids:
        ``(N_dof,)`` Kratos DOF node IDs (equals ``node_ids`` for FE, or control-point
        / THB master node IDs for IGA/THB).
    eval_basis_fn:
        Optional callable ``fn(global_pts (Q, 3)) -> sp.csr_matrix (Q, N_dof)``
        evaluating basis weights in ``dof_node_ids`` ordering.
    """

    name: str
    kind: Literal["fe", "iga", "thb"]
    coords: np.ndarray
    triangles: np.ndarray
    node_ids: np.ndarray | None = None
    dof_node_ids: np.ndarray | None = None
    eval_basis_fn: Callable[[np.ndarray], sp.csr_matrix] | None = None

    @classmethod
    def from_fe(
        cls,
        name: str,
        mesh: MeshModel,
        mesh_to_kratos_node: dict[int, int],
        face_selector: Callable[[np.ndarray], bool] | None = None,
        coords_override: np.ndarray | None = None,
    ) -> "SurfaceMesh":
        """Extract boundary surface triangles from a :class:`MeshModel`.

        Parameters
        ----------
        name:
            Surface name.
        mesh:
            Source volumetric :class:`MeshModel`.
        mesh_to_kratos_node:
            Mapping from 0-based ``mesh.nodes`` index to Kratos node ID.
        face_selector:
            Optional predicate on face centroid ``fn(centroid) -> bool``.
        coords_override:
            Optional transformed node coordinates matching ``mesh.nodes``.
        """
        coords_all = (
            np.asarray(coords_override, dtype=np.float64)
            if coords_override is not None
            else np.asarray(mesh.nodes, dtype=np.float64)
        )
        has_quad_elems = any(b.elem_type in ("hex20", "tet10") for b in mesh.blocks)
        selected_tris: list[tuple[int, int, int]] = []
        used_verts: set[int] = set()

        if not has_quad_elems:
            _, tri_indices, _ = outer_faces(mesh)
            for t in tri_indices:
                t_verts = (int(t[0]), int(t[1]), int(t[2]))
                centroid = coords_all[list(t_verts)].mean(axis=0)
                if face_selector is not None and not face_selector(centroid):
                    continue
                selected_tris.append(t_verts)
                used_verts.update(t_verts)
        else:
            # Preserve quadratic midside nodes on exterior boundary faces.
            _HEX20_FACES = (
                (0, 3, 2, 1, 11, 10, 9, 8),
                (4, 5, 6, 7, 12, 13, 14, 15),
                (0, 1, 5, 4, 8, 17, 12, 16),
                (2, 3, 7, 6, 10, 19, 14, 18),
                (0, 4, 7, 3, 16, 15, 19, 11),
                (1, 2, 6, 5, 9, 18, 13, 17),
            )
            _TET10_FACES = (
                (0, 2, 1, 6, 5, 4),
                (0, 1, 3, 4, 8, 7),
                (1, 2, 3, 5, 9, 8),
                (2, 0, 3, 6, 7, 9),
            )
            _LIN_TABLES = {
                "hex8": ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (0, 4, 7, 3), (1, 2, 6, 5)),
                "tet4": ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
                "wedge6": ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
            }
            counts: dict[tuple[int, ...], int] = {}
            keep: dict[tuple[int, ...], tuple[str, tuple[int, ...]]] = {}
            for blk in mesh.blocks:
                if blk.elem_type == "hex20":
                    ftab, n_corner, fkind = _HEX20_FACES, 4, "quad8"
                elif blk.elem_type == "tet10":
                    ftab, n_corner, fkind = _TET10_FACES, 3, "tri6"
                else:
                    ftab = _LIN_TABLES[blk.elem_type]
                    n_corner, fkind = None, "lin"
                for elem in blk.conn:
                    for fdef in ftab:
                        fnodes = tuple(int(elem[idx]) for idx in fdef)
                        key = tuple(sorted(fnodes[:n_corner] if n_corner else fnodes))
                        counts[key] = counts.get(key, 0) + 1
                        keep[key] = (fkind, fnodes)
            for key, cnt in counts.items():
                if cnt != 1:
                    continue
                fkind, fv = keep[key]
                centroid = coords_all[list(fv)].mean(axis=0)
                if face_selector is not None and not face_selector(centroid):
                    continue
                if fkind == "quad8":
                    # Sub-triangulate 8-node face preserving all corner + midside nodes
                    c0, c1, c2, c3, m01, m12, m23, m30 = fv
                    tris = [
                        (c0, m01, m30),
                        (m01, c1, m12),
                        (m23, m12, c2),
                        (m30, m23, c3),
                        (m01, m12, m30),
                        (m12, m23, m30),
                    ]
                elif fkind == "tri6":
                    c0, c1, c2, m01, m12, m20 = fv
                    tris = [
                        (c0, m01, m20),
                        (m01, c1, m12),
                        (m20, m12, c2),
                        (m01, m12, m20),
                    ]
                elif len(fv) == 3:
                    tris = [(fv[0], fv[1], fv[2])]
                else:
                    tris = [(fv[0], fv[1], fv[2]), (fv[0], fv[2], fv[3])]
                for t in tris:
                    selected_tris.append(t)
                    used_verts.update(t)

        if not selected_tris:
            raise ValueError(f"Surface '{name}' matched 0 boundary faces")

        vert_order = sorted(used_verts)
        remap = {v: idx for idx, v in enumerate(vert_order)}
        surf_coords = coords_all[vert_order]
        surf_node_ids = np.array(
            [int(mesh_to_kratos_node[v]) for v in vert_order], dtype=np.int64
        )
        tri_arr = np.array(
            [[remap[t[0]], remap[t[1]], remap[t[2]]] for t in selected_tris],
            dtype=np.int64,
        )
        selected_tri_sets = {
            frozenset((int(t[0]), int(t[1]), int(t[2]))) for t in selected_tris
        }
        eval_fn = _build_fe_surface_evaluator(
            mesh, coords_all, selected_tri_sets, remap, len(vert_order)
        )
        return cls(
            name=name,
            kind="fe",
            coords=surf_coords,
            triangles=tri_arr,
            node_ids=surf_node_ids,
            dof_node_ids=surf_node_ids,
            eval_basis_fn=eval_fn,
        )

    @classmethod
    def from_iga_face(
        cls,
        name: str,
        box: Any,
        axis: Literal[0, 1, 2],
        side: Literal["lo", "hi"],
        dof_node_ids: np.ndarray | list[int],
        n_div: tuple[int, int] | None = None,
        transform_points: Callable[[np.ndarray], np.ndarray] | None = None,
        inv_transform_points: Callable[[np.ndarray], np.ndarray] | None = None,
        T_basis: sp.spmatrix | None = None,
    ) -> "SurfaceMesh":
        """Build a triangulated axis face on a :class:`BoxSpline` with exact spline evaluation.

        Parameters
        ----------
        name:
            Surface name.
        box:
            Part-local :class:`BoxSpline`.
        axis:
            Normal axis index ``0``, ``1``, or ``2``.
        side:
            ``"lo"`` (``lo[axis]``) or ``"hi"`` (``hi[axis]``).
        dof_node_ids:
            Active Kratos DOF node IDs (either ``box.n_dofs`` control points or
            ``T_basis.shape[1]`` THB/extension reduced DOFs).
        n_div:
            Surface quadrature subdivision counts along the two in-plane axes.
            Defaults to twice the spline cell count per axis.
        transform_points:
            Local-to-global point mapping for transformed instances.
        inv_transform_points:
            Global-to-local point mapping for evaluating spline basis at global points.
        T_basis:
            Optional prolongation operator ``T`` (``(box.n_dofs, n_red)``) when
            DOFs are in reduced THB/extended coordinates.
        """
        in_axes = [a for a in (0, 1, 2) if a != axis]
        a0, a1 = in_axes[0], in_axes[1]
        if n_div is None:
            n0 = max(2, int(box.n_cells[a0]) * 2)
            n1 = max(2, int(box.n_cells[a1]) * 2)
        else:
            n0, n1 = int(n_div[0]), int(n_div[1])

        u_vals = np.linspace(float(box.lo[a0]), float(box.hi[a0]), n0 + 1)
        v_vals = np.linspace(float(box.lo[a1]), float(box.hi[a1]), n1 + 1)
        c_val = float(box.lo[axis] if side == "lo" else box.hi[axis])

        local_pts = np.empty(((n0 + 1) * (n1 + 1), 3), dtype=np.float64)
        for j, x1 in enumerate(v_vals):
            for i, x0 in enumerate(u_vals):
                idx = i + (n0 + 1) * j
                local_pts[idx, axis] = c_val
                local_pts[idx, a0] = x0
                local_pts[idx, a1] = x1

        global_pts = transform_points(local_pts) if transform_points is not None else local_pts

        tris: list[tuple[int, int, int]] = []
        for j in range(n1):
            for i in range(n0):
                p00 = i + (n0 + 1) * j
                p10 = (i + 1) + (n0 + 1) * j
                p01 = i + (n0 + 1) * (j + 1)
                p11 = (i + 1) + (n0 + 1) * (j + 1)
                tris.append((p00, p10, p11))
                tris.append((p00, p11, p01))

        T_csr = T_basis.tocsr() if T_basis is not None else None
        dof_ids_arr = np.asarray(dof_node_ids, dtype=np.int64)

        def _eval(q_global: np.ndarray) -> sp.csr_matrix:
            from wht_iga.grid.bspline import basis_1d

            q_local = (
                inv_transform_points(q_global)
                if inv_transform_points is not None
                else q_global
            )
            q_clamped = np.clip(
                q_local,
                np.asarray(box.lo, dtype=np.float64),
                np.asarray(box.hi, dtype=np.float64),
            )
            ku, kv, kw = box.knots
            p = int(box.degree)
            nu, nv, _ = (int(v) for v in box.n_cp)
            rows: list[int] = []
            cols: list[int] = []
            vals: list[float] = []
            for qi, (xu, yv, zw) in enumerate(q_clamped):
                su, nu_vals = basis_1d(ku, p, float(xu), 0)
                sv, nv_vals = basis_1d(kv, p, float(yv), 0)
                sw, nw_vals = basis_1d(kw, p, float(zw), 0)
                bu = nu_vals[0]
                bv = nv_vals[0]
                bw = nw_vals[0]
                iu0 = su - p
                iv0 = sv - p
                iw0 = sw - p
                for c_w in range(p + 1):
                    kw_idx = iw0 + c_w
                    ww = float(bw[c_w])
                    for c_v in range(p + 1):
                        kv_idx = iv0 + c_v
                        wv = ww * float(bv[c_v])
                        for c_u in range(p + 1):
                            ku_idx = iu0 + c_u
                            w_val = wv * float(bu[c_u])
                            if abs(w_val) > 1e-16:
                                rows.append(qi)
                                cols.append(int(ku_idx + nu * (kv_idx + nv * kw_idx)))
                                vals.append(w_val)
            N_fine = sp.csr_matrix(
                (vals, (rows, cols)), shape=(len(q_clamped), int(box.n_ctrl))
            )
            if T_csr is not None:
                return (N_fine @ T_csr).tocsr()
            return N_fine

        return cls(
            name=name,
            kind="thb" if T_csr is not None else "iga",
            coords=global_pts,
            triangles=np.asarray(tris, dtype=np.int64),
            node_ids=None,
            dof_node_ids=dof_ids_arr,
            eval_basis_fn=_eval,
        )

    def characteristic_edge_length(self) -> float:
        """Compute median edge length across all triangles."""
        t = self.triangles
        p0 = self.coords[t[:, 0]]
        p1 = self.coords[t[:, 1]]
        p2 = self.coords[t[:, 2]]
        edges = np.concatenate(
            [
                np.linalg.norm(p1 - p0, axis=1),
                np.linalg.norm(p2 - p1, axis=1),
                np.linalg.norm(p0 - p2, axis=1),
            ]
        )
        return float(np.median(edges))


@dataclass
class TieSpec:
    """Interface Tie coupling specification between two surfaces."""

    name: str
    master: SurfaceMesh
    slave: SurfaceMesh
    method: Literal["auto", "mpc_matching", "mortar_l2", "penalty_spring"] = "auto"
    tol: float = 1e-6
    penalty_beta: float = 1e4
    E_max: float = 1e5
    h_char: float | None = None
    gauss_order: int = 3


@dataclass
class TieResult:
    """Compiled Tie coupling output (MPC relations and/or sparse penalty matrix)."""

    spec: TieSpec
    resolved_method: str
    relations: list[LinearRelation] = field(default_factory=list)
    penalty_matrix: sp.csr_matrix | None = None
    penalty_node_ids: np.ndarray | None = None
    alpha_tie: float = 0.0
    max_gap: float = 0.0


def _project_points_to_surface(
    surface: SurfaceMesh,
    query_pts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find closest triangle index, barycentric weights ``(Q, 3)``, and distances."""
    tris = surface.triangles
    coords = surface.coords
    n_q = len(query_pts)
    tri_idx = np.zeros(n_q, dtype=np.int64)
    bary = np.zeros((n_q, 3), dtype=np.float64)
    dists = np.zeros(n_q, dtype=np.float64)

    tri_centroids = coords[tris].mean(axis=1)
    for qi, q in enumerate(query_pts):
        # Order candidate triangles by centroid distance for fast pruning while checking close set.
        cdists = np.linalg.norm(tri_centroids - q[None, :], axis=1)
        cand_order = np.argsort(cdists)[: min(24, len(tris))]
        best_d = float("inf")
        best_t = 0
        best_cp = coords[tris[0, 0]]
        for ti in cand_order:
            t = tris[ti]
            a, b, c = coords[t[0]], coords[t[1]], coords[t[2]]
            cp = _closest_point_on_triangle(q, a, b, c)
            d = float(np.linalg.norm(q - cp))
            if d < best_d:
                best_d = d
                best_t = int(ti)
                best_cp = cp
                if d < 1e-13:
                    break
        # Fallback full scan if nearest centroid heuristic had gap > tol.
        if best_d > 1e-6:
            for ti in range(len(tris)):
                t = tris[ti]
                a, b, c = coords[t[0]], coords[t[1]], coords[t[2]]
                cp = _closest_point_on_triangle(q, a, b, c)
                d = float(np.linalg.norm(q - cp))
                if d < best_d:
                    best_d = d
                    best_t = int(ti)
                    best_cp = cp
        best_w = _triangle_barycentric_in_plane(
            coords[tris[best_t]], best_cp[None, :]
        )[0]
        tri_idx[qi] = best_t
        bary[qi] = best_w
        dists[qi] = best_d
    return tri_idx, bary, dists


def _evaluate_surface_basis(
    surface: SurfaceMesh,
    query_pts: np.ndarray,
) -> sp.csr_matrix:
    """Evaluate basis matrix ``N(query_pts)`` of shape ``(Q, N_dof)``."""
    if surface.eval_basis_fn is not None:
        return surface.eval_basis_fn(query_pts).tocsr()
    if surface.node_ids is None:
        raise ValueError(f"Surface '{surface.name}' lacks node_ids and eval_basis_fn")
    tri_idx, bary, _ = _project_points_to_surface(surface, query_pts)
    n_q = len(query_pts)
    n_dof = len(surface.node_ids)
    rows = np.repeat(np.arange(n_q, dtype=np.int64), 3)
    cols = surface.triangles[tri_idx].ravel()
    vals = bary.ravel()
    return sp.csr_matrix((vals, (rows, cols)), shape=(n_q, n_dof))


def _integrate_surface_quadrature(
    surface: SurfaceMesh,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate physical surface Gauss points ``(Q, 3)`` and area weights ``(Q,)``."""
    bary_ref, w_ref = _unit_tri_gauss(order)
    tris = surface.triangles
    coords = surface.coords
    pts_list: list[np.ndarray] = []
    wts_list: list[np.ndarray] = []
    for t in tris:
        v0, v1, v2 = coords[t[0]], coords[t[1]], coords[t[2]]
        area2 = float(np.linalg.norm(np.cross(v1 - v0, v2 - v0)))
        pts_list.append(bary_ref @ np.vstack([v0, v1, v2]))
        wts_list.append(w_ref * area2)
    return np.vstack(pts_list), np.concatenate(wts_list)


def _clipped_interface_quadrature(
    slave: SurfaceMesh,
    master: SurfaceMesh,
    order: int,
    normal_tol: float = 1e-4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute exact sub-triangulated overlap quadrature between two triangulated surfaces.

    Projects nearly coplanar master triangles onto each slave triangle's plane,
    performs Sutherland-Hodgman polygon clipping, and evaluates Gauss rules on
    every clipped sub-triangle so piecewise-linear basis products ``N_s^T N_m``
    are integrated to machine precision.

    Falls back to slave triangle Gauss quadrature + closest-point master projection
    for any slave triangle not covered by coplanar polygon overlap (e.g., curved).
    """
    bary_ref, w_ref = _unit_tri_gauss(order)
    s_tris = slave.triangles
    m_tris = master.triangles
    s_coords = slave.coords
    m_coords = master.coords

    # Precompute master triangle bounding boxes and normals.
    m_v = m_coords[m_tris]  # (Mm, 3, 3)
    m_bmin = m_v.min(axis=1)
    m_bmax = m_v.max(axis=1)
    m_normals = np.cross(m_v[:, 1] - m_v[:, 0], m_v[:, 2] - m_v[:, 0])
    m_area2 = np.linalg.norm(m_normals, axis=1)
    m_normals = np.where(
        m_area2[:, None] > 1e-24, m_normals / np.maximum(m_area2[:, None], 1e-24), 0.0
    )

    q_pts_all: list[np.ndarray] = []
    q_wts_all: list[np.ndarray] = []
    s_tidx_all: list[int] = []
    s_bary_all: list[np.ndarray] = []
    m_tidx_all: list[int] = []
    m_bary_all: list[np.ndarray] = []

    char_h = max(slave.characteristic_edge_length(), master.characteristic_edge_length(), 1e-12)
    plane_dist_tol = max(normal_tol * char_h, 1e-6)

    for si, st in enumerate(s_tris):
        sv = s_coords[st]  # (3, 3)
        e1 = sv[1] - sv[0]
        e2 = sv[2] - sv[0]
        cr = np.cross(e1, e2)
        s_a2 = float(np.linalg.norm(cr))
        if s_a2 <= 1e-24:
            continue
        sn = cr / s_a2
        u_ax = e1 / np.linalg.norm(e1)
        v_ax = np.cross(sn, u_ax)

        sv_2d = np.column_stack([(sv - sv[0]) @ u_ax, (sv - sv[0]) @ v_ax])
        s_bmin = sv.min(axis=0) - plane_dist_tol
        s_bmax = sv.max(axis=0) + plane_dist_tol

        cand_mask = np.all(m_bmax >= s_bmin, axis=1) & np.all(m_bmin <= s_bmax, axis=1)
        covered_area = 0.0
        for mi in np.where(cand_mask)[0]:
            if m_area2[mi] <= 1e-24:
                continue
            if abs(float(np.dot(sn, m_normals[mi]))) < 0.999:
                continue
            mv = m_v[mi]
            plane_d = np.max(np.abs((mv - sv[0]) @ sn))
            if plane_d > plane_dist_tol:
                continue
            mv_2d = np.column_stack([(mv - sv[0]) @ u_ax, (mv - sv[0]) @ v_ax])
            poly_2d = _clip_polygon_2d(mv_2d, sv_2d)
            if len(poly_2d) < 3:
                continue
            centroid_2d = poly_2d.mean(axis=0)
            n_poly = len(poly_2d)
            for k in range(n_poly):
                p0 = centroid_2d
                p1 = poly_2d[k]
                p2 = poly_2d[(k + 1) % n_poly]
                sub_a2 = abs(float(np.cross(p1 - p0, p2 - p0)))
                if sub_a2 <= 1e-22:
                    continue
                covered_area += 0.5 * sub_a2
                p0_3 = sv[0] + p0[0] * u_ax + p0[1] * v_ax
                p1_3 = sv[0] + p1[0] * u_ax + p1[1] * v_ax
                p2_3 = sv[0] + p2[0] * u_ax + p2[1] * v_ax
                sub_verts = np.vstack([p0_3, p1_3, p2_3])
                q_sub = bary_ref @ sub_verts
                w_sub = w_ref * sub_a2
                sb_sub = _triangle_barycentric_in_plane(sv, q_sub)
                mb_sub = _triangle_barycentric_in_plane(mv, q_sub)
                q_pts_all.append(q_sub)
                q_wts_all.append(w_sub)
                s_tidx_all.extend([si] * len(q_sub))
                s_bary_all.append(sb_sub)
                m_tidx_all.extend([int(mi)] * len(q_sub))
                m_bary_all.append(mb_sub)

        if covered_area < 0.5 * 0.5 * s_a2:
            # Non-coplanar fallback: sample slave triangle directly.
            q_sub = bary_ref @ sv
            w_sub = w_ref * s_a2
            mi_arr, mb_sub, _ = _project_points_to_surface(master, q_sub)
            q_pts_all.append(q_sub)
            q_wts_all.append(w_sub)
            s_tidx_all.extend([si] * len(q_sub))
            s_bary_all.append(bary_ref.copy())
            m_tidx_all.extend(mi_arr.tolist())
            m_bary_all.append(mb_sub)

    if not q_pts_all:
        raise RuntimeError("No interface overlap quadrature points generated between Tie surfaces")

    n_q = sum(len(x) for x in q_pts_all)
    n_s = len(slave.coords)
    n_m = len(master.coords)
    s_bary = np.vstack(s_bary_all)
    m_bary = np.vstack(m_bary_all)
    s_idx = np.asarray(s_tidx_all, dtype=np.int64)
    m_idx = np.asarray(m_tidx_all, dtype=np.int64)

    q_rows = np.repeat(np.arange(n_q, dtype=np.int64), 3)
    Ns = sp.csr_matrix((s_bary.ravel(), (q_rows, s_tris[s_idx].ravel())), shape=(n_q, n_s))
    Nm = sp.csr_matrix((m_bary.ravel(), (q_rows, m_tris[m_idx].ravel())), shape=(n_q, n_m))
    return np.vstack(q_pts_all), np.concatenate(q_wts_all), Ns, Nm


def build_tie(spec: TieSpec) -> TieResult:
    """Compile a :class:`TieSpec` into exact MPC relations or penalty stiffness matrix.

    Parameters
    ----------
    spec:
        Tie constraint descriptor between master and slave surfaces.

    Returns
    -------
    TieResult
        Contains compiled :class:`LinearRelation`s (for FE-FE matching or mortar)
        and/or exact symmetric PSD penalty matrix ``penalty_matrix`` (for IGA).
    """
    method = spec.method
    has_iga = spec.master.kind in ("iga", "thb") or spec.slave.kind in ("iga", "thb")

    if method == "auto":
        if has_iga:
            method = "penalty_spring"
        else:
            assert spec.slave.node_ids is not None and spec.master.node_ids is not None
            dists = np.linalg.norm(
                spec.slave.coords[:, None, :] - spec.master.coords[None, :, :],
                axis=2,
            )
            if np.all(dists.min(axis=1) <= spec.tol):
                method = "mpc_matching"
            else:
                method = "mortar_l2"

    comp_names = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")

    if method == "mpc_matching":
        assert spec.slave.node_ids is not None and spec.master.node_ids is not None
        dists = np.linalg.norm(
            spec.slave.coords[:, None, :] - spec.master.coords[None, :, :],
            axis=2,
        )
        min_idx = dists.argmin(axis=1)
        min_d = dists[np.arange(len(dists)), min_idx]
        max_gap = float(min_d.max()) if len(min_d) else 0.0
        rels: list[LinearRelation] = []
        if max_gap <= spec.tol:
            for si, mi in enumerate(min_idx):
                snid = int(spec.slave.node_ids[si])
                mnid = int(spec.master.node_ids[mi])
                if snid == mnid:
                    continue
                for cname in comp_names:
                    rels.append(
                        LinearRelation(
                            slave_dof=(snid, cname),
                            terms={(mnid, cname): 1.0},
                            constant=0.0,
                            tag=f"tie_1to1:{spec.name}",
                        )
                    )
        else:
            tri_idx, bary, proj_d = _project_points_to_surface(spec.master, spec.slave.coords)
            max_gap = float(proj_d.max()) if len(proj_d) else 0.0
            for si in range(len(spec.slave.node_ids)):
                snid = int(spec.slave.node_ids[si])
                mt = spec.master.triangles[tri_idx[si]]
                bw = bary[si]
                for cname in comp_names:
                    terms: dict[DofKey, float] = {}
                    for lv, w in zip(mt, bw):
                        if abs(float(w)) > 1e-14:
                            mnid = int(spec.master.node_ids[int(lv)])
                            terms[(mnid, cname)] = terms.get((mnid, cname), 0.0) + float(w)
                    rels.append(
                        LinearRelation(
                            slave_dof=(snid, cname),
                            terms=terms,
                            constant=0.0,
                            tag=f"tie_bary:{spec.name}",
                        )
                    )
        return TieResult(
            spec=spec,
            resolved_method="mpc_matching",
            relations=rels,
            max_gap=max_gap,
        )

    if method == "mortar_l2":
        if spec.slave.node_ids is None or spec.master.node_ids is None:
            raise ValueError("mortar_l2 requires node_ids on both FE surfaces")
        q_pts, q_wts, Ns_tri, Nm_tri = _clipped_interface_quadrature(
            spec.slave, spec.master, order=max(spec.gauss_order, 4)
        )
        Ns = (
            spec.slave.eval_basis_fn(q_pts).tocsr()
            if spec.slave.eval_basis_fn is not None
            else Ns_tri
        )
        Nm = (
            spec.master.eval_basis_fn(q_pts).tocsr()
            if spec.master.eval_basis_fn is not None
            else Nm_tri
        )
        W = sp.diags(q_wts, format="csr")
        M_ss = (Ns.T @ W @ Ns).toarray()
        C_sm = (Ns.T @ W @ Nm).toarray()
        # Handle near-zero diagonal safely if an isolated vertex had tiny numerical area
        diag_m = np.diag(M_ss)
        max_diag = float(np.max(diag_m)) if len(diag_m) else 1.0
        for idx_d in range(len(diag_m)):
            if diag_m[idx_d] < 1e-14 * max_diag:
                M_ss[idx_d, idx_d] = 1.0
        P = np.linalg.solve(M_ss, C_sm)
        _, _, gaps = _project_points_to_surface(spec.master, spec.slave.coords)
        max_gap = float(gaps.max()) if len(gaps) else 0.0

        rels = []
        for si, snid in enumerate(spec.slave.node_ids):
            row = P[si]
            nz = np.where(np.abs(row) > 1e-13)[0]
            if len(nz) == 0:
                continue
            weights = row[nz]
            w_sum = float(np.sum(weights))
            if abs(w_sum - 1.0) < 1e-6 and abs(w_sum) > 1e-15:
                weights = weights / w_sum
            for cname in comp_names:
                terms = {
                    (int(spec.master.node_ids[int(mi)]), cname): float(w)
                    for mi, w in zip(nz, weights)
                }
                rels.append(
                    LinearRelation(
                        slave_dof=(int(snid), cname),
                        terms=terms,
                        constant=0.0,
                        tag=f"tie_mortar_l2:{spec.name}",
                    )
                )
        return TieResult(
            spec=spec,
            resolved_method="mortar_l2",
            relations=rels,
            max_gap=max_gap,
        )

    if method == "penalty_spring":
        q_pts, q_wts = _integrate_surface_quadrature(spec.slave, order=max(spec.gauss_order, 4))
        Ns = _evaluate_surface_basis(spec.slave, q_pts)
        Nm = _evaluate_surface_basis(spec.master, q_pts)
        assert spec.slave.dof_node_ids is not None and spec.master.dof_node_ids is not None

        h = spec.h_char
        if h is None:
            h = min(
                spec.slave.characteristic_edge_length(),
                spec.master.characteristic_edge_length(),
            )
        h = max(float(h), 1e-12)
        alpha = float(spec.penalty_beta) * float(spec.E_max) / h

        # Build union node ID table so K_tie applies even if DOF sets overlap or differ.
        all_nids = np.unique(
            np.concatenate([spec.slave.dof_node_ids, spec.master.dof_node_ids])
        )
        nid_to_pos = {int(nid): idx for idx, nid in enumerate(all_nids)}
        s_map = [nid_to_pos[int(n)] for n in spec.slave.dof_node_ids]
        m_map = [nid_to_pos[int(n)] for n in spec.master.dof_node_ids]

        # Embed Ns, Nm into common column space (Q, len(all_nids))
        Ns_coo = Ns.tocoo()
        Ns_emb = sp.csr_matrix(
            (Ns_coo.data, (Ns_coo.row, np.asarray(s_map, dtype=np.int64)[Ns_coo.col])),
            shape=(len(q_pts), len(all_nids)),
        )
        Nm_coo = Nm.tocoo()
        Nm_emb = sp.csr_matrix(
            (Nm_coo.data, (Nm_coo.row, np.asarray(m_map, dtype=np.int64)[Nm_coo.col])),
            shape=(len(q_pts), len(all_nids)),
        )
        D = (Ns_emb - Nm_emb).tocsr()
        W_alpha = sp.diags(alpha * q_wts, format="csr")
        K_scalar = (D.T @ W_alpha @ D).tocsr()
        K_scalar = 0.5 * (K_scalar + K_scalar.T)
        K_tie_3d = sp.kron(K_scalar, sp.eye(3, format="csr"), format="csr")

        _, _, gaps = _project_points_to_surface(spec.master, q_pts)
        max_gap = float(gaps.max()) if len(gaps) else 0.0
        return TieResult(
            spec=spec,
            resolved_method="penalty_spring",
            relations=[],
            penalty_matrix=K_tie_3d,
            penalty_node_ids=all_nids,
            alpha_tie=alpha,
            max_gap=max_gap,
        )

    raise ValueError(f"Unsupported Tie method '{spec.method}'")


def apply_tie(
    KM: Any,
    model_part: Any,
    spec: TieSpec,
    id_table: Any | None = None,
    couplings_smp: Any | None = None,
    start_constraint_id: int = 1,
    skip_fixed_slaves: bool = True,
) -> TieResult:
    """Build and attach a :class:`TieSpec` to ``model_part``.

    For ``mpc_matching`` and ``mortar_l2``, inserts ``LinearMasterSlaveConstraint``s
    into ``model_part`` (filtering already-fixed Dirichlet slave DOFs if
    ``skip_fixed_slaves=True`` so Kratos BlockBuilder never over-constrains).
    For ``penalty_spring``, returns the exact sparse penalty matrix ``K_tie`` in
    :class:`TieResult` for solver assembly.
    """
    res = build_tie(spec)
    if res.relations:
        rels = res.relations
        if skip_fixed_slaves:
            filtered: list[LinearRelation] = []
            for rel in rels:
                s_node = model_part.GetNode(int(rel.slave_dof[0]))
                s_var = getattr(KM, rel.slave_dof[1])
                if not s_node.IsFixed(s_var):
                    filtered.append(rel)
            rels = filtered
            res.relations = rels
        target_mp = model_part
        scope = f"Couplings.{spec.name}"
        if couplings_smp is not None:
            if not couplings_smp.HasSubModelPart(spec.name):
                couplings_smp.CreateSubModelPart(spec.name)
        apply_relations_to_model_part(
            KM,
            target_mp,
            rels,
            id_table=id_table,
            scope=scope,
            start_constraint_id=start_constraint_id,
        )
    return res


@dataclass
class RBE2Spec:
    """Rigid kinematic coupling (RBE2) from one 6-DOF reference node to slave points."""

    name: str
    ref_node_id: int
    ref_coord: tuple[float, float, float]
    slave_node_ids: np.ndarray | list[int]
    slave_coords: np.ndarray


@dataclass
class RBE2Result:
    """Compiled RBE2 kinematic matrix ``C_rbe2`` and constraint relations."""

    spec: RBE2Spec
    operator: np.ndarray
    relations: list[LinearRelation]


def build_rbe2(spec: RBE2Spec, prune_tol: float = 1e-14) -> RBE2Result:
    """Construct RBE2 rigid kinematic relation ``u_i = u_ref + theta_ref x (x_i - x_ref)``.

    Parameters
    ----------
    spec:
        RBE2 descriptor with reference point and slave nodes.
    prune_tol:
        Threshold below which zero lever-arm entries are omitted from sparse MPC terms.

    Returns
    -------
    RBE2Result
        Dense ``(3*N_slave, 6)`` matrix and S04-compliant :class:`LinearRelation` list.
    """
    x_ref = np.asarray(spec.ref_coord, dtype=np.float64)
    coords = np.asarray(spec.slave_coords, dtype=np.float64)
    n_s = len(spec.slave_node_ids)
    if coords.shape != (n_s, 3):
        raise ValueError(f"slave_coords shape {coords.shape} != ({n_s}, 3)")

    C = np.zeros((3 * n_s, 6), dtype=np.float64)
    master_keys: list[DofKey] = [
        (int(spec.ref_node_id), "DISPLACEMENT_X"),
        (int(spec.ref_node_id), "DISPLACEMENT_Y"),
        (int(spec.ref_node_id), "DISPLACEMENT_Z"),
        (int(spec.ref_node_id), "ROTATION_X"),
        (int(spec.ref_node_id), "ROTATION_Y"),
        (int(spec.ref_node_id), "ROTATION_Z"),
    ]
    disp_names = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")

    rels: list[LinearRelation] = []
    for i, snid in enumerate(spec.slave_node_ids):
        rx, ry, rz = coords[i] - x_ref
        # theta x r = [tz*(-ry) + ty*rz, tx*rz - tz*rx, ty*rx - tx*ry]
        block = np.array(
            [
                [1.0, 0.0, 0.0, 0.0, rz, -ry],
                [0.0, 1.0, 0.0, -rz, 0.0, rx],
                [0.0, 0.0, 1.0, ry, -rx, 0.0],
            ],
            dtype=np.float64,
        )
        C[3 * i : 3 * i + 3, :] = block
        for r_idx, cname in enumerate(disp_names):
            terms = {
                master_keys[c_idx]: float(block[r_idx, c_idx])
                for c_idx in range(6)
                if abs(float(block[r_idx, c_idx])) > prune_tol
            }
            rels.append(
                LinearRelation(
                    slave_dof=(int(snid), cname),
                    terms=terms,
                    constant=0.0,
                    tag=f"rbe2:{spec.name}",
                )
            )
    return RBE2Result(spec=spec, operator=C, relations=rels)


def _ensure_ref_node(
    KM: Any,
    model_part: Any,
    target_smp: Any,
    ref_node_id: int,
    ref_coord: tuple[float, float, float],
    with_rotation: bool = True,
) -> Any:
    """Fetch existing reference node or create it with 6 DOFs."""
    if model_part.HasNode(int(ref_node_id)):
        node = model_part.GetNode(int(ref_node_id))
    else:
        node = target_smp.CreateNewNode(
            int(ref_node_id),
            float(ref_coord[0]),
            float(ref_coord[1]),
            float(ref_coord[2]),
        )
    for dv, rv in (
        (KM.DISPLACEMENT_X, KM.REACTION_X),
        (KM.DISPLACEMENT_Y, KM.REACTION_Y),
        (KM.DISPLACEMENT_Z, KM.REACTION_Z),
    ):
        if not node.HasDofFor(dv):
            node.AddDof(dv, rv)
    if with_rotation:
        for dv, rv in (
            (KM.ROTATION_X, KM.REACTION_MOMENT_X),
            (KM.ROTATION_Y, KM.REACTION_MOMENT_Y),
            (KM.ROTATION_Z, KM.REACTION_MOMENT_Z),
        ):
            if not node.HasDofFor(dv):
                node.AddDof(dv, rv)
    return node


def apply_rbe2(
    KM: Any,
    model_part: Any,
    spec: RBE2Spec,
    id_table: Any | None = None,
    couplings_smp: Any | None = None,
    start_constraint_id: int = 1,
    skip_fixed_slaves: bool = True,
) -> RBE2Result:
    """Attach RBE2 rigid kinematic constraints into ``model_part``."""
    scope = f"Couplings.{spec.name}"
    target_smp = model_part
    if couplings_smp is not None:
        target_smp = (
            couplings_smp.GetSubModelPart(spec.name)
            if couplings_smp.HasSubModelPart(spec.name)
            else couplings_smp.CreateSubModelPart(spec.name)
        )
    if not model_part.HasNode(int(spec.ref_node_id)) and id_table is not None:
        # Track externally allocated ref node only if requested id matches next_id.
        if int(spec.ref_node_id) == id_table.next_id("nodes"):
            id_table.allocate_nodes(scope, 1, provenance="reference_point")
    _ensure_ref_node(KM, model_part, target_smp, spec.ref_node_id, spec.ref_coord, with_rotation=True)

    res = build_rbe2(spec)
    rels = res.relations
    if skip_fixed_slaves:
        rels = [
            rel
            for rel in rels
            if not model_part.GetNode(int(rel.slave_dof[0])).IsFixed(getattr(KM, rel.slave_dof[1]))
        ]
        res.relations = rels
    apply_relations_to_model_part(
        KM,
        model_part,
        rels,
        id_table=id_table,
        scope=scope,
        start_constraint_id=start_constraint_id,
    )
    return res


@dataclass
class RBE3Spec:
    """Non-stiffening distributing coupling (RBE3) descriptor."""

    name: str
    ref_node_id: int
    ref_coord: tuple[float, float, float]
    slave_coords: np.ndarray
    slave_node_ids: np.ndarray | list[int] | None = None
    weights: np.ndarray | None = None
    weight_mode: Literal["uniform", "area", "inverse_distance"] = "uniform"
    slave_basis_mat: sp.spmatrix | np.ndarray | None = None


@dataclass
class RBE3Result:
    """Compiled RBE3 weighted least-squares operator ``W_rbe3`` and MPCs."""

    spec: RBE3Spec
    operator: np.ndarray
    weights: np.ndarray
    centroid: np.ndarray
    inertia_tensor: np.ndarray
    dof_node_ids: np.ndarray
    relations: list[LinearRelation]


def _skew(v: np.ndarray) -> np.ndarray:
    """Return 3x3 skew-symmetric cross-product matrix ``[v]_x`` so ``[v]_x w = v x w``."""
    return np.array(
        [
            [0.0, -float(v[2]), float(v[1])],
            [float(v[2]), 0.0, -float(v[0])],
            [-float(v[1]), float(v[0]), 0.0],
        ],
        dtype=np.float64,
    )


def build_rbe3(spec: RBE3Spec, prune_tol: float = 1e-14) -> RBE3Result:
    """Build non-stiffening RBE3 least-squares operator ``q_ref = W_rbe3 u_surf``.

    Satisfies exact rigid-body reproduction (zero strain and exact ``(u_0, omega)``
    recovery for any rigid motion) and exact virtual-work static force/moment
    equivalence ``f = W_rbe3^T [F; M]`` with zero inter-node stiffening.
    """
    pts = np.asarray(spec.slave_coords, dtype=np.float64)
    n_pts = pts.shape[0]
    if pts.ndim != 2 or pts.shape[1] != 3 or n_pts < 3:
        raise ValueError(f"RBE3 requires at least 3 3D slave points, got shape {pts.shape}")

    x_ref = np.asarray(spec.ref_coord, dtype=np.float64)
    if spec.weights is not None:
        w = np.asarray(spec.weights, dtype=np.float64)
    elif spec.weight_mode == "uniform":
        w = np.ones(n_pts, dtype=np.float64)
    elif spec.weight_mode == "inverse_distance":
        d = np.linalg.norm(pts - x_ref[None, :], axis=1)
        scale = max(float(np.median(d)), 1e-12)
        w = 1.0 / np.maximum(d, 1e-3 * scale)
    elif spec.weight_mode == "area":
        raise ValueError("weight_mode='area' requires explicit tributary area weights array")
    else:
        raise ValueError(f"Unknown RBE3 weight_mode '{spec.weight_mode}'")

    if w.shape != (n_pts,) or np.any(w <= 0.0):
        raise ValueError("RBE3 weights must be strictly positive 1D array matching slave_coords")

    w_norm = w / float(np.sum(w))
    x_bar = (w_norm[:, None] * pts).sum(axis=0)
    d_vecs = pts - x_bar[None, :]
    c_vec = x_ref - x_bar

    J = np.zeros((3, 3), dtype=np.float64)
    for wi, di in zip(w, d_vecs):
        J += float(wi) * (float(np.dot(di, di)) * np.eye(3) - np.outer(di, di))
    eig_min = float(np.min(np.linalg.eigvalsh(J)))
    eig_max = float(np.max(np.linalg.eigvalsh(J)))
    if eig_max <= 0.0 or eig_min <= 1e-12 * eig_max:
        raise ValueError("RBE3 slave points are collinear/degenerate (singular inertia tensor)")
    J_inv = np.linalg.inv(J)
    skew_c = _skew(c_vec)

    W_pts = np.zeros((6, 3 * n_pts), dtype=np.float64)
    for i in range(n_pts):
        Ri = float(w[i]) * (J_inv @ _skew(d_vecs[i]))
        Ti = float(w_norm[i]) * np.eye(3) - (skew_c @ Ri)
        W_pts[0:3, 3 * i : 3 * i + 3] = Ti
        W_pts[3:6, 3 * i : 3 * i + 3] = Ri

    if spec.slave_basis_mat is not None:
        N_mat = (
            spec.slave_basis_mat.toarray()
            if sp.issparse(spec.slave_basis_mat)
            else np.asarray(spec.slave_basis_mat, dtype=np.float64)
        )
        if N_mat.shape[0] != n_pts:
            raise ValueError(
                f"slave_basis_mat rows {N_mat.shape[0]} != slave_coords count {n_pts}"
            )
        W_op = W_pts @ np.kron(N_mat, np.eye(3))
        n_dofs = N_mat.shape[1]
    else:
        W_op = W_pts
        n_dofs = n_pts

    if spec.slave_node_ids is None:
        dof_ids = np.arange(1, n_dofs + 1, dtype=np.int64)
    else:
        dof_ids = np.asarray(spec.slave_node_ids, dtype=np.int64)
        if len(dof_ids) != n_dofs:
            raise ValueError(
                f"slave_node_ids length {len(dof_ids)} != DOF node count {n_dofs}"
            )

    ref_dofs: list[DofKey] = [
        (int(spec.ref_node_id), "DISPLACEMENT_X"),
        (int(spec.ref_node_id), "DISPLACEMENT_Y"),
        (int(spec.ref_node_id), "DISPLACEMENT_Z"),
        (int(spec.ref_node_id), "ROTATION_X"),
        (int(spec.ref_node_id), "ROTATION_Y"),
        (int(spec.ref_node_id), "ROTATION_Z"),
    ]
    comp_names = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")

    rels: list[LinearRelation] = []
    for r_idx, sdof in enumerate(ref_dofs):
        terms: dict[DofKey, float] = {}
        for j, nid in enumerate(dof_ids):
            for c_idx, cname in enumerate(comp_names):
                val = float(W_op[r_idx, 3 * j + c_idx])
                if abs(val) > prune_tol:
                    terms[(int(nid), cname)] = terms.get((int(nid), cname), 0.0) + val
        if terms:
            rels.append(
                LinearRelation(
                    slave_dof=sdof,
                    terms=terms,
                    constant=0.0,
                    tag=f"rbe3:{spec.name}",
                )
            )

    return RBE3Result(
        spec=spec,
        operator=W_op,
        weights=w,
        centroid=x_bar,
        inertia_tensor=J,
        dof_node_ids=dof_ids,
        relations=rels,
    )


def distribute_rbe3_loads(
    rbe3_result: RBE3Result,
    force: tuple[float, float, float] = (0.0, 0.0, 0.0),
    moment: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Compute virtual-work statically equivalent nodal loads ``f = W_rbe3^T [F; M]``."""
    gen_load = np.array(
        [
            float(force[0]),
            float(force[1]),
            float(force[2]),
            float(moment[0]),
            float(moment[1]),
            float(moment[2]),
        ],
        dtype=np.float64,
    )
    f_vec = rbe3_result.operator.T @ gen_load
    return f_vec.reshape(-1, 3)


def apply_rbe3(
    KM: Any,
    SM: Any,
    model_part: Any,
    spec: RBE3Spec,
    force: tuple[float, float, float] = (0.0, 0.0, 0.0),
    moment: tuple[float, float, float] = (0.0, 0.0, 0.0),
    mode: Literal["mpc", "load_distribution"] = "mpc",
    id_table: Any | None = None,
    couplings_smp: Any | None = None,
    start_constraint_id: int = 1,
    start_condition_id: int = 1,
) -> RBE3Result:
    """Apply RBE3 distributing coupling and optional reference force/moment to ``model_part``.

    In ``"mpc"`` mode, the reference node is created with 6 DOFs and attached as
    the **slave** of surface nodes via ``LinearMasterSlaveConstraint``s (so no
    artificial stiffness is added between surface nodes), and ``PointLoadCondition3D1N``
    / ``PointMomentCondition3D1N`` are created at the reference node.
    In ``"load_distribution"`` mode, statically equivalent nodal forces
    :func:`distribute_rbe3_loads` are added directly as ``POINT_LOAD`` +
    ``PointLoadCondition3D1N`` on the surface nodes.
    """
    scope = f"Couplings.{spec.name}"
    target_smp = model_part
    if couplings_smp is not None:
        target_smp = (
            couplings_smp.GetSubModelPart(spec.name)
            if couplings_smp.HasSubModelPart(spec.name)
            else couplings_smp.CreateSubModelPart(spec.name)
        )

    res = build_rbe3(spec)
    f_vec = np.asarray(force, dtype=np.float64)
    m_vec = np.asarray(moment, dtype=np.float64)

    if mode == "mpc":
        if not model_part.HasNode(int(spec.ref_node_id)) and id_table is not None:
            if int(spec.ref_node_id) == id_table.next_id("nodes"):
                id_table.allocate_nodes(scope, 1, provenance="reference_point")
        ref_node = _ensure_ref_node(
            KM, model_part, target_smp, spec.ref_node_id, spec.ref_coord, with_rotation=True
        )
        apply_relations_to_model_part(
            KM,
            model_part,
            res.relations,
            id_table=id_table,
            scope=scope,
            start_constraint_id=start_constraint_id,
        )
        n_conds = int(np.any(f_vec != 0.0)) + int(np.any(m_vec != 0.0))
        if n_conds > 0:
            cid = (
                id_table.allocate_conditions(scope, n_conds, kind="PointLoadMoment").start_id
                if id_table is not None
                else start_condition_id
            )
            props = next(iter(model_part.Properties), None)
            if props is None:
                props = model_part.CreateNewProperties(
                    id_table.allocate_properties(scope, 1).start_id
                    if id_table is not None
                    else 99999
                )
            if np.any(f_vec != 0.0):
                ref_node.SetSolutionStepValue(SM.POINT_LOAD, 0, [float(x) for x in f_vec])
                target_smp.CreateNewCondition("PointLoadCondition3D1N", cid, [int(ref_node.Id)], props)
                cid += 1
            if np.any(m_vec != 0.0):
                ref_node.SetSolutionStepValue(SM.POINT_MOMENT, 0, [float(x) for x in m_vec])
                target_smp.CreateNewCondition("PointMomentCondition3D1N", cid, [int(ref_node.Id)], props)
        return res

    if mode == "load_distribution":
        nodal_f = distribute_rbe3_loads(res, force=force, moment=moment)
        nz_rows = [
            idx for idx in range(len(res.dof_node_ids)) if np.linalg.norm(nodal_f[idx]) > 1e-18
        ]
        if nz_rows:
            cid = (
                id_table.allocate_conditions(scope, len(nz_rows), kind="PointLoadCondition3D1N").start_id
                if id_table is not None
                else start_condition_id
            )
            props = next(iter(model_part.Properties), None)
            if props is None:
                props = model_part.CreateNewProperties(
                    id_table.allocate_properties(scope, 1).start_id
                    if id_table is not None
                    else 99999
                )
            for idx in nz_rows:
                nid = int(res.dof_node_ids[idx])
                nd = model_part.GetNode(nid)
                cur = np.array(nd.GetSolutionStepValue(SM.POINT_LOAD, 0), dtype=np.float64)
                upd = cur + nodal_f[idx]
                nd.SetSolutionStepValue(SM.POINT_LOAD, 0, [float(x) for x in upd])
                target_smp.CreateNewCondition("PointLoadCondition3D1N", cid, [nid], props)
                cid += 1
        return res

    raise ValueError(f"Unsupported RBE3 mode '{mode}'")


__all__ = [
    "ConstraintGraphReport",
    "DofKey",
    "LinearRelation",
    "RBE2Result",
    "RBE2Spec",
    "RBE3Result",
    "RBE3Spec",
    "SurfaceMesh",
    "TieResult",
    "TieSpec",
    "apply_rbe2",
    "apply_rbe3",
    "apply_relations_to_model_part",
    "apply_tie",
    "build_rbe2",
    "build_rbe3",
    "build_tie",
    "compile_constraint_graph",
    "distribute_rbe3_loads",
    "validate_constraint_graph",
]
