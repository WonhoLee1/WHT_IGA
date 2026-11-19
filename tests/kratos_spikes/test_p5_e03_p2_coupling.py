# -*- coding: utf-8 -*-
"""P5 / E03: Degree p=2 trimmed 3-region connector with annular RBE3 mounting coupling.

Extends the S33 coarse p=1 E03 proof-of-concept to:
  1. Degree ``p = 2`` C^1 trimmed B-spline partition on a compact ``(9, 6, 6)``
     background ``BoxSpline`` (``q = 2``, ``n_fallback == 0``, volume error < 1%).
  2. Multi-material union extended B-spline stabilization (``theta = 0.2``) so
     small-support boundary basis functions are smoothly extrapolated with zero
     ``THIN_FEATURE`` fallbacks.
  3. Coreform-equivalent annular distributing coupling (``DISTRIBUTED_COUPLING``
     with fixed reference point at ``(-25.5, 0, 0)``) on the physical header
     mounting annulus ``x = -25.50 mm``, enforced via ``add_linear_tie`` after
     flattening the constraint graph (``c = 0``, chain-free).
  4. Zero-stress uniform free thermal expansion verification (``< 1e-6 MPa``)
     under both minimal 3-2-1 restraint and annular RBE3 coupling.
  5. Full 3-step thermoplastic cycle (450 -> 20 -> -55 -> 200 C) with
     ``TotalLagrangianElement3D4N``, ``finite_strain=True``, and
     ``TP.e03_materials`` (isotropic hardening + temperature-dependent tables)
     using ``PardisoLUSolver``, verifying 0 cutbacks, hardening, monotonic PEEQ,
     glass clamping/reheat tension, and peak RSS < 1.5 GB.
"""

from __future__ import annotations

import gc
import os
import time
from pathlib import Path
from typing import Any

import numpy as np
import psutil
import pytest
import scipy.linalg

from wht_iga.analyses import nonlinear as NL
from wht_iga.analyses import thermal_plasticity as TP
from wht_iga.geometry.cad import brep_area, brep_volume, shape_bbox, tessellate, validate_shape
from wht_iga.geometry.occ_io import _count_free_edges, load_cad_file
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.interactions import (
    LinearRelation,
    RBE3Spec,
    apply_relations_to_model_part,
    build_rbe3,
    compile_constraint_graph,
)
from wht_iga.kratos_backend.solvers import build_kratos_linear_solver
from wht_iga.kratos_backend.surface import triangle_quadrature
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import fix_inactive, inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table
from wht_iga.stabilize import extension as EXT
from wht_iga.stabilize import support as SUP

pytestmark = [pytest.mark.kratos]

_E03_ROOT = (
    Path(__file__).resolve().parents[2]
    / "research-cf-IGA"
    / "user-guide"
    / "example-problems"
    / "glass-to-metal-seal"
)
_STP_PATH = _E03_ROOT / "downloads" / "glass_to_metal_seal_connector.stp"
_DATA_DIR = _E03_ROOT / "data"
_COMP_NAMES = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")
_GIB = float(1024 ** 3)


def _rss_stats_mb() -> tuple[float, float]:
    """Return ``(current_rss_mb, peak_working_set_mb)`` for the current process."""
    mi = psutil.Process(os.getpid()).memory_info()
    rss_mb = float(mi.rss) / (1024.0 ** 2)
    peak_mb = float(getattr(mi, "peak_wset", mi.rss)) / (1024.0 ** 2)
    return rss_mb, peak_mb


def _merge_tessellations(solids: list[Any], deflection: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Tessellate disjoint solids individually and concatenate into one triangle soup."""
    nodes_list: list[np.ndarray] = []
    tris_list: list[np.ndarray] = []
    offset = 0
    for solid in solids:
        tess = tessellate(solid, deflection=deflection)
        nodes_list.append(tess.nodes)
        tris_list.append(tess.tris + offset)
        offset += len(tess.nodes)
    return np.vstack(nodes_list), np.vstack(tris_list)


def _classify_e03_solids(solids: list[Any]) -> dict[str, list[int]]:
    """Identify header, glass seal, and 11 pins geometrically (independent of STEP order)."""
    vols = [brep_volume(s) for s in solids]
    header = [i for i, v in enumerate(vols) if v > 1000.0]
    glass = [i for i, v in enumerate(vols) if 100.0 < v < 500.0]
    pins = [i for i, v in enumerate(vols) if 10.0 < v < 30.0]
    if len(header) != 1 or len(glass) != 1 or len(pins) != 11:
        raise AssertionError(
            f"unexpected solid classification: header={header}, glass={glass}, pins={pins}"
        )
    return {"header": header, "glass": glass, "pins": pins}


def _header_mounting_face_area(header_solid: Any, x_mount: float = -25.50) -> float:
    """Exact B-Rep surface area of the planar mounting annulus at ``x = x_mount``."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    area_sum = 0.0
    exp = TopExp_Explorer(header_solid, TopAbs_FACE)
    while exp.More():
        face = TopoDS.Face(exp.Current())
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        com = props.CentreOfMass()
        if abs(float(com.X()) - x_mount) < 1e-3:
            area_sum += float(props.Mass())
        exp.Next()
    return area_sum


def _extract_mounting_surface_quadrature(
    h_nodes: np.ndarray,
    h_tris: np.ndarray,
    box: BoxSpline,
    x_mount: float = -25.50,
    tri_order: int = 2,
) -> dict[str, Any]:
    """Build surface quadrature points and p=2 B-spline basis matrix on ``x = x_mount``."""
    tri_pts = h_nodes[h_tris]
    mask = np.all(np.abs(tri_pts[:, :, 0] - x_mount) < 1e-3, axis=1)
    mount_tris = tri_pts[mask]
    if len(mount_tris) == 0:
        raise AssertionError(f"no header triangles found on x = {x_mount}")

    ref_pts, ref_w = triangle_quadrature(tri_order)
    q_coords: list[np.ndarray] = []
    q_weights: list[float] = []
    for tri in mount_tris:
        e1 = tri[1] - tri[0]
        e2 = tri[2] - tri[0]
        twice_area = float(np.linalg.norm(np.cross(e1, e2)))
        for (xi, eta), w in zip(ref_pts, ref_w):
            q_coords.append(tri[0] + xi * e1 + eta * e2)
            q_weights.append(float(w) * twice_area)

    q_arr = np.asarray(q_coords, dtype=float)
    w_arr = np.asarray(q_weights, dtype=float)
    centroid = (w_arr[:, None] * q_arr).sum(axis=0) / float(w_arr.sum())

    nu, nv, nw = (int(v) for v in box.n_cp)
    p = int(box.degree)
    surf_cp_set: set[int] = set()
    basis_rows: list[dict[int, float]] = []
    for px, py, pz in q_arr:
        sx = find_span(nu, p, float(px), box.knots[0])
        sy = find_span(nv, p, float(py), box.knots[1])
        sz = find_span(nw, p, float(pz), box.knots[2])
        nx_v = ders_basis_funs(sx, float(px), p, 0, box.knots[0])[0]
        ny_v = ders_basis_funs(sy, float(py), p, 0, box.knots[1])[0]
        nz_v = ders_basis_funs(sz, float(pz), p, 0, box.knots[2])[0]
        row_dict: dict[int, float] = {}
        for c in range(p + 1):
            for b in range(p + 1):
                for a in range(p + 1):
                    idx = (sx - p + a) + nu * ((sy - p + b) + nv * (sz - p + c))
                    val = float(nx_v[a] * ny_v[b] * nz_v[c])
                    if val > 1e-15:
                        row_dict[idx] = val
                        surf_cp_set.add(idx)
        basis_rows.append(row_dict)

    return {
        "n_tris": int(len(mount_tris)),
        "q_coords": q_arr,
        "q_weights": w_arr,
        "area": float(w_arr.sum()),
        "centroid": centroid,
        "surf_cp_set": surf_cp_set,
        "basis_rows": basis_rows,
    }


def _principal_and_mises(cauchy_voigt: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(von_mises, max_principal, min_principal)`` for ``(N, 6)`` Kratos Voigt stress."""
    sxx, syy, szz, sxy, syz, sxz = cauchy_voigt.T
    vm = np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy ** 2 + syz ** 2 + sxz ** 2)
    )
    mat3 = np.zeros((len(cauchy_voigt), 3, 3), dtype=float)
    mat3[:, 0, 0], mat3[:, 1, 1], mat3[:, 2, 2] = sxx, syy, szz
    mat3[:, 0, 1] = mat3[:, 1, 0] = sxy
    mat3[:, 1, 2] = mat3[:, 2, 1] = syz
    mat3[:, 0, 2] = mat3[:, 2, 0] = sxz
    eigs = np.linalg.eigvalsh(mat3)
    return vm, eigs[:, 2], eigs[:, 0]


def _make_pardiso_newton_strategy(
    KM: Any,
    LA: Any,
    mp: Any,
    max_iters: int = 30,
    tol_rel: float = 1e-8,
    tol_abs: float = 1e-10,
    move_mesh: bool = True,
) -> Any:
    """Newton-Raphson strategy using MKL ``PardisoLUSolver`` for unsymmetric p=2+MPC systems."""
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = build_kratos_linear_solver(KM, "pardiso_lu", LA=LA)
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    rc = KM.ResidualCriteria(tol_rel, tol_abs)
    rc.SetEchoLevel(0)
    dc = KM.DisplacementCriteria(tol_rel, tol_abs)
    dc.SetEchoLevel(0)
    conv = KM.AndCriteria(rc, dc)
    conv.SetEchoLevel(0)
    strat = KM.ResidualBasedNewtonRaphsonStrategy(
        mp, scheme, conv, builder, int(max_iters), True, False, bool(move_mesh)
    )
    strat.SetEchoLevel(0)
    return strat


def _build_extension_relations(op: EXT.ExtensionOperator) -> list[LinearRelation]:
    """Convert sparse extension operator ``u_slave = E u_master`` into ``LinearRelation``s."""
    e_csr = op.E.tocsr()
    rels: list[LinearRelation] = []
    for s_id in op.slave_ids:
        row = e_csr.getrow(int(s_id))
        for cname in _COMP_NAMES:
            terms = {
                (int(op.master_ids[m_col]) + 1, cname): float(w)
                for m_col, w in zip(row.indices, row.data)
                if abs(float(w)) > 1e-15
            }
            rels.append(
                LinearRelation(
                    slave_dof=(int(s_id) + 1, cname),
                    terms=terms,
                    constant=0.0,
                    tag="ext_p2",
                )
            )
    return rels


def _build_annular_rbe3(
    bundle: dict[str, Any],
    inactive_mask: np.ndarray,
    op: EXT.ExtensionOperator,
) -> tuple[Any, list[int], list[LinearRelation], float]:
    """Construct exact continuum RBE3 distributing constraint ``W_op u_surf = 0``.

    Reduces ``W_op`` onto active master DOFs via ``E_surf``, selects 6 well-conditioned
    pivot DOFs on the ``i = 0`` mounting layer via column-pivoted QR, and emits
    6 ``LinearRelation``s ``u_piv = -(W_S^{-1} W_M) u_M`` (``c = 0``).
    """
    box: BoxSpline = bundle["box"]
    mount = bundle["mount"]
    nu = int(box.n_cp[0])

    surf_cps = sorted(I for I in mount["surf_cp_set"] if not bool(inactive_mask[I]))
    cp_col = {I: idx for idx, I in enumerate(surf_cps)}
    n_mat = np.zeros((len(mount["q_coords"]), len(surf_cps)), dtype=float)
    for r_i, row_dict in enumerate(mount["basis_rows"]):
        for idx, val in row_dict.items():
            if idx in cp_col:
                n_mat[r_i, cp_col[idx]] = val

    x_bar = tuple(float(v) for v in mount["centroid"])
    rbe3 = build_rbe3(
        RBE3Spec(
            name="E03_mount_annulus",
            ref_node_id=0,
            ref_coord=x_bar,
            slave_coords=mount["q_coords"],
            weights=mount["q_weights"],
            slave_basis_mat=n_mat,
            slave_node_ids=[I + 1 for I in surf_cps],
        )
    )

    e_surf = op.E.tocsr()[surf_cps, :].toarray()
    w_master = rbe3.operator @ np.kron(e_surf, np.eye(3))
    i0_cols = [c for c in range(w_master.shape[1]) if (int(op.master_ids[c // 3]) % nu) == 0]
    if len(i0_cols) < 6:
        raise RuntimeError(f"need >= 6 i=0 master DOFs for RBE3 pivots, got {len(i0_cols)}")
    _, _, piv = scipy.linalg.qr(w_master[:, i0_cols], pivoting=True)
    piv_cols = [i0_cols[int(idx)] for idx in piv[:6]]
    piv_set = set(piv_cols)
    master_cols = [c for c in range(w_master.shape[1]) if c not in piv_set]
    w_s = w_master[:, piv_cols]
    cond_ws = float(np.linalg.cond(w_s))
    r_mat = -np.linalg.solve(w_s, w_master[:, master_cols])

    rbe3_rels: list[LinearRelation] = []
    for row_i, sc in enumerate(piv_cols):
        sdof = (int(op.master_ids[sc // 3]) + 1, _COMP_NAMES[sc % 3])
        terms: dict[tuple[int, str], float] = {}
        for col_j, mc in enumerate(master_cols):
            val = float(r_mat[row_i, col_j])
            if abs(val) > 1e-14:
                terms[(int(op.master_ids[mc // 3]) + 1, _COMP_NAMES[mc % 3])] = val
        rbe3_rels.append(
            LinearRelation(
                slave_dof=sdof,
                terms=terms,
                constant=0.0,
                tag="mount_rbe3",
            )
        )
    return rbe3, surf_cps, rbe3_rels, cond_ws


@pytest.fixture(scope="module")
def e03_p2_bundle() -> dict[str, Any]:
    """Load E03 13-solid STEP once per module and build p=2 trimmed partition + extension."""
    t0 = time.perf_counter()
    cad = load_cad_file(_STP_PATH, heal=False)
    groups = _classify_e03_solids(cad.solids)

    solid_info = []
    for idx, solid in enumerate(cad.solids):
        tess = tessellate(solid, deflection=0.05)
        rep = validate_shape(solid, tess, n_samples=10)
        solid_info.append({
            "index": idx,
            "volume": brep_volume(solid),
            "area": brep_area(solid),
            "bbox": shape_bbox(solid),
            "free_edges": _count_free_edges(solid),
            "watertight": rep.watertight,
            "tess": tess,
        })

    h_tess = solid_info[groups["header"][0]]["tess"]
    g_tess = solid_info[groups["glass"][0]]["tess"]
    p_nodes, p_tris = _merge_tessellations(
        [cad.solids[i] for i in groups["pins"]], deflection=0.05
    )

    lo, hi = cad.bbox
    pad = np.array([0.25, 0.05, 0.05], dtype=float)
    box = BoxSpline(lo=tuple(lo - pad), hi=tuple(hi + pad), n_cells=(9, 6, 6), degree=2)
    q_order = 2

    nx, ny, nz = (int(v) for v in box.n_cells)
    n_cells_tot = nx * ny * nz
    cell_vol = np.full(
        n_cells_tot,
        float(np.prod([(box.hi[d] - box.lo[d]) / box.n_cells[d] for d in range(3)])),
    )

    region_tables: dict[str, dict[str, Any]] = {}
    mat_etas: list[np.ndarray] = []
    for name, nds, trs, idxs in [
        ("header", h_tess.nodes, h_tess.tris, groups["header"]),
        ("glass", g_tess.nodes, g_tess.tris, groups["glass"]),
        ("pins", p_nodes, p_tris, groups["pins"]),
    ]:
        qp_table, out_cells, issues, recs, diag = build_trimmed_table(
            nds, trs, box, q=q_order, workers=1, cache_dir="off"
        )
        eta = np.zeros(n_cells_tot, dtype=float)
        for rec in recs:
            eta[rec["i"] + nx * (rec["j"] + ny * rec["k"])] = float(rec["eta"])
        mat_etas.append(eta)
        region_tables[name] = {
            "qp_table": qp_table,
            "out_cells": out_cells,
            "issues": issues,
            "recs": recs,
            "diag": diag,
            "brep_volume": sum(solid_info[i]["volume"] for i in idxs),
        }

    cls = SUP.classify_basis(
        box.n_cp,
        box.knots,
        np.zeros(n_cells_tot, dtype=float),
        cell_vol,
        theta=0.2,
        material_etas=mat_etas,
    )
    ext_op = EXT.build_extension(box.n_cp, cls.labels, p=2, knots=box.knots)
    mount = _extract_mounting_surface_quadrature(
        h_tess.nodes, h_tess.tris, box, x_mount=-25.50, tri_order=2
    )
    mount["brep_area"] = _header_mounting_face_area(cad.solids[groups["header"][0]], x_mount=-25.50)

    return {
        "cad": cad,
        "groups": groups,
        "solid_info": solid_info,
        "box": box,
        "q_order": q_order,
        "regions": region_tables,
        "cls": cls,
        "ext_op": ext_op,
        "mount": mount,
        "prep_seconds": time.perf_counter() - t0,
    }


def _build_and_run_e03_p2(
    bundle: dict[str, Any],
    mats: list[TP.ThermoMaterial],
    increments: list[int],
    *,
    bc_mode: str,
    tol_rel: float = 1e-8,
    tol_abs: float = 1e-10,
) -> dict[str, Any]:
    """Assemble and run the p=2 3-region E03 thermal model under ``321`` or ``rbe3`` BCs."""
    t0 = time.perf_counter()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    LA = mods["LinearSolversApplication"]
    box: BoxSpline = bundle["box"]
    cls: SUP.SupportClassification = bundle["cls"]
    ext_op: EXT.ExtensionOperator = bundle["ext_op"]
    nx, ny, nz = (int(v) for v in box.n_cells)

    region_order = ("header", "glass", "pins")
    reg_pairs = [
        (bundle["regions"][name]["qp_table"], bundle["regions"][name]["out_cells"])
        for name in region_order
    ]

    model = KM.Model()
    mp = model.CreateModelPart("E03_P2_Coupled")
    for var in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(var)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    volume, _nodes = make_nurbs_volume(KM, mp, box)
    elements, elem_mat_idx, active_union = TP.attach_region_elements(
        KM,
        mp,
        volume,
        box,
        reg_pairs,
        element_name="TotalLagrangianElement3D4N",
        n_gauss=bundle["q_order"] + 1,
    )

    all_out = [
        (i, j, k)
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
        if (i, j, k) not in active_union
    ]
    inact = inactive_ctrl_mask(box, all_out)
    n_inact = fix_inactive(KM, mp, box, all_out)

    rels = _build_extension_relations(ext_op)
    rbe3_obj = None
    surf_cps: list[int] = []
    cond_ws = 0.0
    p012: tuple[int, int, int] | None = None

    if bc_mode == "321":
        p012 = TP.fix_isostatic_321(
            KM, mp, box, inactive_mask=(cls.labels != SUP.MASTER), axis=0, side=0
        )
        compiled = compile_constraint_graph(rels, chain_policy="reject")
        apply_relations_to_model_part(KM, mp, compiled)
    elif bc_mode == "rbe3":
        rbe3_obj, surf_cps, rbe3_rels, cond_ws = _build_annular_rbe3(bundle, inact, ext_op)
        rels.extend(rbe3_rels)
        compiled = compile_constraint_graph(rels, chain_policy="flatten")
        apply_relations_to_model_part(KM, mp, compiled)
    else:
        raise ValueError(f"unsupported bc_mode={bc_mode!r}")

    mp.CloneTimeStep(0.0)
    states = [450.0, 20.0, -55.0, 200.0]
    groups = TP.CellThermalGroups(
        KM,
        mods,
        mp,
        elements,
        mats,
        elem_mat_idx,
        states,
        law_factory=TP.default_law_factory(mods, finite_strain=True),
    )

    snapshots: list[dict[str, Any]] = []
    peak_rss_mb = _rss_stats_mb()[1]

    def observer(level: float) -> None:
        nonlocal peak_rss_mb
        if abs(level - round(level)) > 1e-9:
            return
        peeq_by_id = groups.equivalent_plastic_strain()
        entry: dict[str, Any] = {"level": float(round(level))}
        for m_idx, name in enumerate(region_order):
            m_els = [el for el, mi in zip(elements, elem_mat_idx) if mi == m_idx]
            sigs = np.vstack([
                np.array([
                    list(v)
                    for v in el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)
                ])
                for el in m_els
            ])
            vm, s1, s3 = _principal_and_mises(sigs)
            diss_max = 0.0
            if mats[m_idx].plastic:
                diss_max = max(
                    float(v)
                    for el in m_els
                    for v in el.CalculateOnIntegrationPoints(KM.PLASTIC_DISSIPATION, mp.ProcessInfo)
                )
            peeq_max = 0.0
            if mats[m_idx].hardening is not None:
                peeq_max = max(float(peeq_by_id[int(el.Id)].max()) for el in m_els)
            entry[name] = {
                "vm_max": float(vm.max()),
                "s1_max": float(s1.max()),
                "s3_min": float(s3.min()),
                "abs_max": float(np.abs(sigs).max()),
                "diss_max": diss_max,
                "peeq_max": peeq_max,
            }

        active_ids = np.flatnonzero(~inact)
        u_act = np.array(
            [list(mp.GetNode(int(I) + 1).GetSolutionStepValue(KM.DISPLACEMENT)) for I in active_ids],
            dtype=float,
        )
        entry["u_ctrl_max"] = float(np.linalg.norm(u_act, axis=1).max())

        if rbe3_obj is not None:
            u_surf = np.array(
                [list(mp.GetNode(int(I) + 1).GetSolutionStepValue(KM.DISPLACEMENT)) for I in surf_cps],
                dtype=float,
            ).ravel()
            entry["rbe3_residual"] = float(np.linalg.norm(rbe3_obj.operator @ u_surf))
            u_surf_mat = u_surf.reshape(-1, 3)
            entry["mount_ctrl_max_yz"] = float(np.linalg.norm(u_surf_mat[:, 1:], axis=1).max())

        peak_rss_mb = max(peak_rss_mb, _rss_stats_mb()[1])
        snapshots.append(entry)

    out = TP.run_thermal_steps(
        mp,
        KM,
        lambda: _make_pardiso_newton_strategy(
            KM, LA, mp, max_iters=30, tol_rel=tol_rel, tol_abs=tol_abs
        ),
        groups,
        increments=increments,
        on_converged=observer,
    )
    rss_now_mb, peak_wset_mb = _rss_stats_mb()
    peak_rss_mb = max(peak_rss_mb, peak_wset_mb)
    wall_s = time.perf_counter() - t0

    return {
        "out": out,
        "snapshots": snapshots,
        "n_elements": len(elements),
        "n_active_cells": len(active_union),
        "n_inactive_ctrl": int(n_inact),
        "n_master_ctrl": int(cls.n_master),
        "n_extended_ctrl": int(cls.n_extended),
        "n_constraints": len(compiled),
        "cond_ws": cond_ws,
        "p012": p012,
        "wall_seconds": wall_s,
        "rss_mb": rss_now_mb,
        "peak_rss_mb": peak_rss_mb,
    }


# --------------------------------------------------------------------------- tests

def test_e03_p2_trimmed_partition_and_extension_support(e03_p2_bundle: dict[str, Any]) -> None:
    """Verify p=2 (9, 6, 6) 3-region partition (0 fallbacks, <1% vol err), extension, & mount face."""
    cad = e03_p2_bundle["cad"]
    info = e03_p2_bundle["solid_info"]
    box: BoxSpline = e03_p2_bundle["box"]
    cls: SUP.SupportClassification = e03_p2_bundle["cls"]
    ext_op: EXT.ExtensionOperator = e03_p2_bundle["ext_op"]
    mount = e03_p2_bundle["mount"]
    nx, ny, nz = (int(v) for v in box.n_cells)
    n_gauss = e03_p2_bundle["q_order"] + 1

    assert len(cad.solids) == 13
    assert cad.healing_report.watertight is True
    assert all(s["free_edges"] == 0 and s["watertight"] for s in info)
    assert box.degree == 2 and box.n_cells == (9, 6, 6)

    active_union: set[tuple[int, int, int]] = set()
    total_qp_vol = 0.0
    total_brep_vol = 0.0
    for name, reg in e03_p2_bundle["regions"].items():
        diag = reg["diag"]
        assert diag["n_fallback"] == 0, f"{name} had fallback cells: {diag}"
        qp_table = reg["qp_table"]
        out_set = set(reg["out_cells"])
        vol_qp = 0.0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    cell = (i, j, k)
                    if cell in out_set:
                        continue
                    active_union.add(cell)
                    pts = (
                        qp_table[cell]
                        if cell in qp_table
                        else np.asarray(cell_gauss_points(box, i, j, k, n_gauss=n_gauss), dtype=float)
                    )
                    assert np.all(pts[:, 3] > 0.0)
                    vol_qp += float(pts[:, 3].sum())
        bvol = float(reg["brep_volume"])
        rel_err = abs(vol_qp - bvol) / bvol
        assert rel_err < 0.01, f"{name} volume error {rel_err:.4e} >= 1%"
        total_qp_vol += vol_qp
        total_brep_vol += bvol

    tot_rel_err = abs(total_qp_vol - total_brep_vol) / total_brep_vol
    assert tot_rel_err < 0.01

    assert ext_op.stats["n_constant_fallback"] == 0
    assert ext_op.stats["n_degree_reduced"] == 0
    assert ext_op.stats["n_thin_feature"] == 0
    row_sums = np.asarray(ext_op.E[ext_op.slave_ids].sum(axis=1)).ravel()
    assert row_sums == pytest.approx(1.0, abs=1e-12)

    # Annular mounting surface checks at x = -25.50 mm
    area_rel = abs(mount["area"] - mount["brep_area"]) / mount["brep_area"]
    assert area_rel < 0.005
    assert mount["centroid"][0] == pytest.approx(-25.50, abs=1e-6)
    assert np.hypot(mount["centroid"][1], mount["centroid"][2]) < 0.01

    _, peak_mb = _rss_stats_mb()
    assert peak_mb / 1024.0 < 1.5
    print(
        f"\n[p=2 partition] n_cells={box.n_cells}, n_cp={box.n_cp} ({box.n_dofs} DOFs), "
        f"active_cells={len(active_union)}/{nx*ny*nz}, vol_err={tot_rel_err:.4e}, "
        f"master/ext/inact={cls.n_master}/{cls.n_extended}/{cls.n_inactive}, "
        f"mount_area={mount['area']:.4f} (brep={mount['brep_area']:.4f}), "
        f"prep={e03_p2_bundle['prep_seconds']:.2f}s, peak_rss={peak_mb:.1f} MB"
    )


def test_e03_p2_free_thermal_expansion_zero_stress(e03_p2_bundle: dict[str, Any]) -> None:
    """Uniform SS304L across all 3 regions must remain stress-free (< 1e-6 MPa) for 3-2-1 and RBE3."""
    mats_dict = TP.e03_materials(_DATA_DIR)
    uniform_mats = [mats_dict["SS304L_header"]] * 3

    for bc_mode in ("321", "rbe3"):
        res = _build_and_run_e03_p2(
            e03_p2_bundle,
            uniform_mats,
            increments=[1, 1, 1],
            bc_mode=bc_mode,
            tol_rel=1e-10,
            tol_abs=1e-12,
        )
        assert res["out"]["steps"] == 3
        assert res["out"]["n_solves"] == 3, f"{bc_mode} free expansion had cutbacks"
        worst_stress = 0.0
        for snap in res["snapshots"]:
            for rname in ("header", "glass", "pins"):
                worst_stress = max(worst_stress, snap[rname]["abs_max"])
                assert snap[rname]["abs_max"] < 1e-6
                assert snap[rname]["peeq_max"] == pytest.approx(0.0, abs=1e-15)
            if bc_mode == "rbe3":
                assert snap["rbe3_residual"] < 1e-11
                assert snap["mount_ctrl_max_yz"] > 1e-3

        assert res["peak_rss_mb"] / 1024.0 < 1.5
        print(
            f"\n[p=2 free expansion {bc_mode}] n_solves={res['out']['n_solves']}, "
            f"worst_stress={worst_stress:.3e} MPa, wall={res['wall_seconds']:.2f}s, "
            f"peak_rss={res['peak_rss_mb']:.1f} MB"
        )
        del res
        gc.collect()


def test_e03_p2_annular_coupling_3step_hardening_nlgeom(e03_p2_bundle: dict[str, Any]) -> None:
    """Full 450 -> 20 -> -55 -> 200 C p=2 E03 with RBE3 annular coupling, hardening, and NLGEOM."""
    mats_dict = TP.e03_materials(_DATA_DIR)
    increments = [5, 5, 5]
    res = _build_and_run_e03_p2(
        e03_p2_bundle,
        [mats_dict["SS304L_header"], mats_dict["S8061_glass"], mats_dict["Alloy52_pins"]],
        increments=increments,
        bc_mode="rbe3",
        tol_rel=1e-8,
        tol_abs=1e-10,
    )

    assert all(n <= 2 * 15 for n in increments)
    assert res["out"]["steps"] == sum(increments)
    assert res["out"]["n_solves"] == sum(increments), "unexpected Newton cutback in p=2 RBE3 cycle"
    assert res["out"]["level"] == pytest.approx(3.0)
    assert res["peak_rss_mb"] / 1024.0 < 1.5

    s_room, s_cold, s_hot = res["snapshots"]
    assert [s["level"] for s in res["snapshots"]] == [1.0, 2.0, 3.0]

    # Annular distributing constraint W_op u_surf == 0 satisfied while allowing local radial contraction
    for snap in res["snapshots"]:
        assert snap["rbe3_residual"] < 1e-10
        assert snap["mount_ctrl_max_yz"] > 1e-3

    # Isotropic hardening raises vm_max strictly above initial yield at all 3 endpoints
    assert s_room["header"]["vm_max"] > 162.0
    assert s_cold["header"]["vm_max"] > 194.4
    assert s_hot["header"]["vm_max"] > 120.2
    for snap in (s_room, s_cold, s_hot):
        assert snap["pins"]["vm_max"] > 234.4

    # Monotonically non-decreasing PEEQ_max across the 3-step history
    assert s_room["header"]["peeq_max"] > 0.0
    assert s_cold["header"]["peeq_max"] >= s_room["header"]["peeq_max"]
    assert s_hot["header"]["peeq_max"] >= s_cold["header"]["peeq_max"]
    assert s_room["pins"]["peeq_max"] > 0.0
    assert s_cold["pins"]["peeq_max"] > s_room["pins"]["peeq_max"]
    assert s_hot["pins"]["peeq_max"] >= s_cold["pins"]["peeq_max"]

    # Glass thermoelastic compressive clamping and reheat tensile residual stress rise
    assert s_room["glass"]["s3_min"] < -300.0
    assert s_cold["glass"]["s3_min"] < s_room["glass"]["s3_min"]
    assert s_room["glass"]["s1_max"] > 0.0
    assert s_cold["glass"]["s1_max"] > s_room["glass"]["s1_max"]
    assert s_hot["glass"]["s1_max"] > s_cold["glass"]["s1_max"]

    print(
        f"\n[p=2 RBE3 3-step NLGEOM+hardening] n_elements={res['n_elements']}, "
        f"active_cells={res['n_active_cells']}, master/ext/inact="
        f"{res['n_master_ctrl']}/{res['n_extended_ctrl']}/{res['n_inactive_ctrl']}, "
        f"n_constraints={res['n_constraints']}, cond(Ws)={res['cond_ws']:.2f}, "
        f"steps={res['out']['steps']}, n_solves={res['out']['n_solves']}, "
        f"wall={res['wall_seconds']:.2f}s, rss={res['rss_mb']:.1f} MB, "
        f"peak_rss={res['peak_rss_mb']:.1f} MB"
    )
    for snap in res["snapshots"]:
        print(f"  level={snap['level']}: {snap}")
