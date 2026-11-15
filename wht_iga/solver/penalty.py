"""Penalty Dirichlet BCs: BC-A (full 3-comp) + BC-B (per-component) (Track SOLVE).

BC-A: IgaApp ``SupportPenaltyCondition`` = integral beta (u - g).v dA, penalizes all
3 components simultaneously (v2 section 5.4, pitfall K6; S09 proven).
BC-B: per-component/per-direction (symmetry planes, v2 section 5.4): one virtual
node per surface QP tied by ``LinearMasterSlaveConstraint`` to the spline
field, loaded by ``NodalConcentratedElement3D1N`` (k = beta.dA) +
``PointLoadCondition3D1N`` (F = k.g) (S10 proven element pair; avoids spring
ROTATION dofs per ADR-001 failover #9).

Surface rule order >= 2p enforced (S09: order-2 rule left a beta-independent
4.5e-5 floor; degree-4 Dunavant required for p = 2).

Face triangles are emitted PER FACE CELL (one pair per knot-span rectangle):
B-spline products N_i.N_j are only piecewise smooth, so whole-face triangles
make the Dunavant rule inexact across knot lines (found 2026-11-14: 28%
penalty-matrix error on 2x2 faces, exact on single-cell faces; S09 1-cell
result reproduced to 12 digits with the same builder).
"""

from __future__ import annotations

import numpy as np

from wht_iga.kratos_backend.surface import make_surface_quadrature, triangle_quadrature
from wht_iga.solver.validate import DirichletBC


#: Penalty scale constant beta = DEFAULT_C_PEN . E . p^2 / h.
#: Fixed by ADR-002 from the patch-cube beta sweep (exact 1/beta, C ~ 4-8):
#: c = 3e5 holds V01 max-relative error <= 1e-5 with >= 3x margin (p = 1 worst).
DEFAULT_C_PEN = 3.0e5


def default_beta(E: float, h: float, p: int, c_pen: float = DEFAULT_C_PEN) -> float:
    """Penalty default beta = c_pen . E . p^2 / h (v2 section 5.4; c_pen tuned in ADR-002)."""
    return c_pen * E * p * p / h


def eval_g(bc: DirichletBC, x: np.ndarray) -> np.ndarray:
    """Evaluate the prescribed value g at physical point x (3,)."""
    x = np.asarray(x, dtype=float)
    if bc.field is not None:
        return np.asarray(bc.field.A, dtype=float) @ x + np.asarray(bc.field.b, dtype=float)
    return np.asarray(bc.value, dtype=float)


def rule_for_degree(p: int) -> int:
    """Dunavant rule order: per-face-cell points >= (p+1)^2 (rank rule).

    The 1-cell p=3 case with order 4 (12 pts < 16) showed a beta-independent
    1.8e-4 floor; order 6 (24 pts >= 16) recovers exact 1/beta scaling
    (S09 correction 2026-11-14). Degree >= 2p remains an accuracy
    recommendation, not the floor mechanism.
    """
    if p <= 1:
        return 2
    if p <= 2:
        return 4
    if p == 3:
        return 6
    raise ValueError(f"no Dunavant rule available for p = {p} (triangle_quadrature max: 6)")


def face_subdiv(box: object, axis: int) -> list[np.ndarray]:
    """Breakpoint coordinates along the two face axes (per-cell triangulation)."""
    other = [d for d in range(3) if d != axis]
    return [np.asarray(box.breakpoints(o), dtype=float) for o in other]


def box_face_triangles(
    lo: np.ndarray,
    hi: np.ndarray,
    axis: int,
    side: int,
    subdiv: list[np.ndarray] | None = None,
) -> list[np.ndarray]:
    """Outward-wound triangle pair per face cell (see module docstring for why).

    ``subdiv = [coords_ax0, coords_ax1]`` gives one pair per knot-span
    rectangle; ``None`` gives a single pair (single-cell faces only).
    Local copy (Track GEO owns geometry/); winding verified against the S08
    convention: +PRESSURE pushes along the outward face normal.
    """
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    other = [d for d in range(3) if d != axis]
    bounds = [np.array([lo[o], hi[o]]) for o in other] if subdiv is None else [np.asarray(s) for s in subdiv]
    tris: list[np.ndarray] = []
    for c0 in range(len(bounds[0]) - 1):
        for c1 in range(len(bounds[1]) - 1):
            quad = []
            for s1 in (bounds[0][c0], bounds[0][c0 + 1]):
                for s2 in (bounds[1][c1], bounds[1][c1 + 1]):
                    v = np.array([lo[0], lo[1], lo[2]])
                    v[axis] = hi[axis] if side == 1 else lo[axis]
                    v[other[0]] = s1
                    v[other[1]] = s2
                    quad.append(v)
            cand = [np.array([quad[0], quad[1], quad[2]]), np.array([quad[1], quad[3], quad[2]])]
            n = np.cross(cand[0][1] - cand[0][0], cand[0][2] - cand[0][0])
            outward = np.zeros(3)
            outward[axis] = 1.0 if side == 1 else -1.0
            if float(n @ outward) < 0.0:
                cand = [t[::-1] for t in cand]
            tris.extend(cand)
    return tris


def _qp_physical_xy(triangles: list[np.ndarray], rule_order: int) -> tuple[np.ndarray, np.ndarray]:
    """Physical coords + dA weights of surface QPs in make_surface_quadrature order."""
    from wht_iga.kratos_backend.surface import triangle_area

    ref_pts, ref_w = triangle_quadrature(rule_order)
    xs, dAs = [], []
    for tri in triangles:
        area = triangle_area(np.asarray(tri))
        for (xi, eta), w in zip(ref_pts, ref_w):
            xs.append((1.0 - xi - eta) * tri[0] + xi * tri[1] + eta * tri[2])
            dAs.append(2.0 * area * float(w))
    return np.array(xs), np.array(dAs)


def attach_penalty_bcs(
    KM: object,
    IGA: object,
    model_part: object,
    volume: object,
    box: object,
    bcs: list[DirichletBC],
    beta: float,
    first_cond_id: int = 100,
    first_node_id: int = 10_000,
    first_prop_id: int = 100,
) -> dict:
    """Attach BC-A conditions for full-3-component BCs. Returns keepalive refs.

    Raises ValueError if a BC selects a strict subset of components (use
    ``attach_component_bcs`` -- K6: SupportPenaltyCondition is all-or-nothing).
    """
    keep: dict = {"wrappers": [], "qp_geoms": [], "n_cond": 0}
    cid = first_cond_id
    nid = first_node_id
    pid = first_prop_id
    for bc in bcs:
        if not all(bc.components):
            raise ValueError(f"BC-A requires all 3 components (K6); face {bc.face} -> use attach_component_bcs")
        order = rule_for_degree(box.degree)
        tris = box_face_triangles(np.asarray(box.lo), np.asarray(box.hi), *bc.face,
                                  subdiv=face_subdiv(box, bc.face[0]))
        qp_geoms, wrappers, _, _ = make_surface_quadrature(
            KM, model_part, volume, tris, rule_order=order, first_node_id=nid)
        nid += 3 * len(tris)
        xs, _ = _qp_physical_xy(tris, order)
        assert len(qp_geoms) == len(xs), "surface QP count must match local prediction"
        props = model_part.CreateNewProperties(pid)
        pid += 1
        props.SetValue(IGA.PENALTY_FACTOR, float(beta))
        for qp, x in zip(qp_geoms, xs):
            cond = model_part.CreateNewCondition("SupportPenaltyCondition", cid, qp, props)
            g = KM.Vector(3)
            gv = eval_g(bc, x)
            g[0], g[1], g[2] = float(gv[0]), float(gv[1]), float(gv[2])
            cond.SetValue(KM.DISPLACEMENT, g)
            cid += 1
        keep["wrappers"].extend(wrappers)
        keep["qp_geoms"].extend(qp_geoms)
        keep["n_cond"] += len(qp_geoms)
    return keep


def attach_component_bcs(
    KM: object,
    SM: object,
    model_part: object,
    volume: object,
    box: object,
    bcs: list[DirichletBC],
    beta: float,
    first_cond_id: int = 100,
    first_node_id: int = 10_000,
    first_elem_id: int = 1_000,
    first_constr_id: int = 50_000,
    first_prop_id: int = 100,
) -> dict:
    """Attach BC-B per-component penalty via virtual nodes (v2 section 5.4).

    Per surface QP: virtual node V with slave relation V_d = sum N_i(x_qp) u_i,d
    over the QP support nodes (relation row = ShapeFunctionsValues, checked
    for partition of unity), spring k_d = beta.dA on selected components d via
    ``NodalConcentratedElement3D1N`` and load F_d = k_d.g_d via
    ``PointLoadCondition3D1N``. Minimizer: k(V - g)^2 per component.
    """
    from wht_iga.kratos_backend.constraints import add_linear_tie

    keep: dict = {"wrappers": [], "qp_geoms": [], "virtual_nodes": [], "n_cond": 0}
    cid = first_cond_id
    nid = first_node_id
    eid = first_elem_id
    tid = first_constr_id
    pid = first_prop_id
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    reac_vars = (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z)
    load_vars = (SM.POINT_LOAD_X, SM.POINT_LOAD_Y, SM.POINT_LOAD_Z)
    for bc in bcs:
        order = rule_for_degree(box.degree)
        tris = box_face_triangles(np.asarray(box.lo), np.asarray(box.hi), *bc.face,
                                  subdiv=face_subdiv(box, bc.face[0]))
        qp_geoms, wrappers, _, _ = make_surface_quadrature(
            KM, model_part, volume, tris, rule_order=order, first_node_id=nid)
        nid += 3 * len(tris)
        xs, dAs = _qp_physical_xy(tris, order)
        assert len(qp_geoms) == len(xs)
        for qp, x, dA in zip(qp_geoms, xs, dAs):
            row = np.array(qp.ShapeFunctionsValues(), dtype=float).ravel()
            assert abs(float(row.sum()) - 1.0) < 1e-9, "partition of unity at surface QP"
            support = [nd for nd in qp]
            assert len(support) == len(row), "SFV columns must match support nodes"
            vnode = model_part.CreateNewNode(nid, float(x[0]), float(x[1]), float(x[2]))
            nid += 1
            for d in range(3):
                vnode.AddDof(comp_vars[d], reac_vars[d])  # K14, targeted (no VariableUtils re-sweep)
            slave = [vnode.GetDof(comp_vars[d]) for d in range(3) if bc.components[d]]
            master = [nd.GetDof(comp_vars[d]) for nd in support for d in range(3) if bc.components[d]]
            n_sel = sum(1 for c in bc.components if c)
            rel = np.zeros((n_sel, n_sel * len(support)))
            for k in range(n_sel):
                rel[k, k::n_sel] = row
            add_linear_tie(KM, model_part, tid, slave, master, rel)
            tid += 1
            stiff = [0.0, 0.0, 0.0]
            gv = eval_g(bc, x)
            for d in range(3):
                if bc.components[d]:
                    stiff[d] = float(beta) * float(dA)
            elem = model_part.CreateNewElement("NodalConcentratedElement3D1N", eid, [vnode.Id],
                                               model_part.CreateNewProperties(pid))
            pid += 1
            elem.SetValue(SM.NODAL_DISPLACEMENT_STIFFNESS, stiff)
            eid += 1
            model_part.CreateNewCondition("PointLoadCondition3D1N", cid, [vnode.Id],
                                          model_part.CreateNewProperties(pid))
            pid += 1
            cid += 1
            for d in range(3):
                if bc.components[d]:
                    vnode.SetSolutionStepValue(load_vars[d], stiff[d] * float(gv[d]))
            keep["virtual_nodes"].append(vnode)
        keep["wrappers"].extend(wrappers)
        keep["qp_geoms"].extend(qp_geoms)
        keep["n_cond"] += len(qp_geoms)
    return keep


def analyze_sweep(betas: np.ndarray, errors: np.ndarray) -> dict:
    """Plateau analysis of a beta sweep: O(1/beta) trend -> discretization floor.

    Returns plateau start index (first point within 2x of the minimum error),
    the fitted pre-plateau slope d log(err)/d log(beta), and the recommended
    default (10x above plateau start, rounded to a power of ten).
    """
    betas = np.asarray(betas, dtype=float)
    errors = np.asarray(errors, dtype=float)
    emin = float(errors.min())
    plateau = np.where(errors <= 2.0 * emin)[0]
    start = int(plateau[0]) if len(plateau) else int(np.argmin(errors))
    if start >= 2:
        slope = float(np.polyfit(np.log(betas[:start]), np.log(errors[:start]), 1)[0])
    else:
        slope = float("nan")
    rec = 10.0 ** np.ceil(np.log10(betas[start] * 10.0))
    return {"plateau_start": start, "plateau_beta": float(betas[start]),
            "pre_plateau_slope": slope, "recommended": float(rec), "min_error": emin}
