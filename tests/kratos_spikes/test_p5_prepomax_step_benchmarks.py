# -*- coding: utf-8 -*-
"""P5 / S44: PrePoMax STEP Trimmed IGA & Multi-Body Contact Benchmarks.

Benchmarks 4 PrePoMax v2.1.4 STEP models:
1. ``Dvizna_plosca_podmodel.STEP``: Coarse trimmed BoxSpline (q=2, 0 fallbacks,
   <2% volume error) + Kratos ``SmallDisplacementElement3D4N`` +
   ``fix_inactive`` + ``add_linear_tie`` linear elastic solve.
2. ``Bearing_bracket.STEP``: Axis-aligned vs PCA OBB (``compute_pca_obb``)
   active-cell fraction and quadrature volume comparison + compact Kratos
   ``SmallDisplacementElement3D4N`` linear elastic solve.
3. ``Plates_bolt.STEP``: 5-solid watertight STEP partition on a shared coarse
   background ``BoxSpline`` via ``build_trimmed_table`` (q=2, 0 fallbacks,
   per-solid and total volume error < 2%).
4. ``Pin_on_a_plate.STEP``: 2-solid trimmed quadrature + Kratos
   ``SmallDisplacementElement3D4N`` stiffness assembly (``KratosBodyOperator``)
   + ``solve_deformable_ssn_contact`` with Coulomb friction (mu=0.3),
   verifying ``n_active > 0``, ``F_N > 0``, and interface vs cut-face reaction
   equilibrium < 1e-6.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import (
    DeformableBoxFaceMaster,
    sample_surface_triangles_gauss,
    solve_deformable_ssn_contact,
)
from wht_iga.analyses.contact_nlgeom import KratosBodyOperator
from wht_iga.analyses.methods import compute_pca_obb
from wht_iga.geometry.cad import brep_volume, shape_bbox, tessellate, validate_shape
from wht_iga.geometry.occ_io import _count_free_edges, load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import fix_inactive, inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

_REPO_RESOURCES = Path(__file__).resolve().parents[2] / "wht_iga" / "resources"
_EXTERNAL_PREPOMAX = Path(r"D:\SOFTWARE\PrePoMax v2.1.4 dev\Models")
_PREPOMAX_MODELS = _REPO_RESOURCES if _REPO_RESOURCES.is_dir() else _EXTERNAL_PREPOMAX

pytestmark = pytest.mark.skipif(
    not _PREPOMAX_MODELS.is_dir(),
    reason=f"PrePoMax models directory not found: {_PREPOMAX_MODELS}",
)


def _quadrature_volume(
    box: BoxSpline,
    qp_table: dict,
    out_cells: list[tuple[int, int, int]],
    q: int,
) -> float:
    """Compute total quadrature volume from trimmed + untrimmed active cells."""
    nx, ny, nz = (int(v) for v in box.n_cells)
    out_set = {tuple(c) for c in out_cells}
    n_gauss = int(q) + 1
    vol = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = (i, j, k)
                if cell in out_set:
                    continue
                if cell in qp_table:
                    pts = qp_table[cell]
                else:
                    pts = np.asarray(cell_gauss_points(box, i, j, k, n_gauss=n_gauss), dtype=float)
                vol += float(np.sum(pts[:, 3]))
    return vol


def _solve_kratos_trimmed_elastic(
    box: BoxSpline,
    qp_table: dict,
    out_cells: list[tuple[int, int, int]],
    *,
    E: float = 210_000.0,
    nu: float = 0.3,
    q: int = 2,
    disp_x: float = 0.1,
    use_linear_tie: bool = False,
) -> dict:
    """Assemble and solve a compact Kratos SmallDisplacementElement3D4N model."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    model = KM.Model()
    mp = model.CreateModelPart("PrePoMaxTrimmed")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    volume, _nodes = make_nurbs_volume(KM, mp, box)
    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    nx, ny, nz = (int(v) for v in box.n_cells)
    out_set = {tuple(c) for c in out_cells}
    n_gauss = int(q) + 1
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = (i, j, k)
                if cell in out_set:
                    continue
                if cell in qp_table:
                    pts = np.asarray(qp_table[cell], dtype=float).tolist()
                else:
                    pts = cell_gauss_points(box, i, j, k, n_gauss=n_gauss)
                qps = KM.GeometriesVector()
                vol_qps = qps
                volume.CreateQuadraturePointGeometries(vol_qps, 2, pts)
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, vol_qps[0], props)
                eid += 1

    for v, r in (
        (KM.DISPLACEMENT_X, KM.REACTION_X),
        (KM.DISPLACEMENT_Y, KM.REACTION_Y),
        (KM.DISPLACEMENT_Z, KM.REACTION_Z),
    ):
        KM.VariableUtils().AddDof(v, r, mp)

    fix_inactive(KM, mp, box, out_cells)
    inact = inactive_ctrl_mask(box, out_cells)
    cp = box.control_points()
    x_min = float(box.lo[0])
    x_max = float(box.hi[0])
    tol = 1e-9 * max(1.0, x_max - x_min)

    active_ctrl = [i for i in range(box.n_ctrl) if not inact[i]]
    left_ids = [i for i in active_ctrl if abs(float(cp[i, 0]) - x_min) < tol]
    right_ids = [i for i in active_ctrl if abs(float(cp[i, 0]) - x_max) < tol]
    assert len(left_ids) > 0 and len(right_ids) > 0

    for idx in left_ids:
        node = mp.GetNode(idx + 1)
        for comp in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            node.Fix(comp)
            node.SetSolutionStepValue(comp, 0.0)

    if use_linear_tie and len(right_ids) > 1:
        master_node = mp.GetNode(right_ids[0] + 1)
        master_node.Fix(KM.DISPLACEMENT_X)
        master_node.SetSolutionStepValue(KM.DISPLACEMENT_X, float(disp_x))
        m_dof = [master_node.GetDof(KM.DISPLACEMENT_X)]
        rel = np.ones((1, 1), dtype=float)
        for cid, idx in enumerate(right_ids[1:], start=1):
            s_node = mp.GetNode(idx + 1)
            add_linear_tie(KM, mp, cid, [s_node.GetDof(KM.DISPLACEMENT_X)], m_dof, rel)
    else:
        for idx in right_ids:
            node = mp.GetNode(idx + 1)
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, float(disp_x))

    builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    strat = KM.ResidualBasedLinearStrategy(mp, scheme, builder, True, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.Solve()

    ctrl_nodes = [mp.GetNode(i + 1) for i in range(box.n_ctrl)]
    U = np.array(
        [
            [
                float(n.GetSolutionStepValue(KM.DISPLACEMENT_X)),
                float(n.GetSolutionStepValue(KM.DISPLACEMENT_Y)),
                float(n.GetSolutionStepValue(KM.DISPLACEMENT_Z)),
            ]
            for n in ctrl_nodes
        ],
        dtype=float,
    )
    R = np.array(
        [
            [
                float(n.GetSolutionStepValue(KM.REACTION_X)),
                float(n.GetSolutionStepValue(KM.REACTION_Y)),
                float(n.GetSolutionStepValue(KM.REACTION_Z)),
            ]
            for n in ctrl_nodes
        ],
        dtype=float,
    )
    return {
        "U": U,
        "R": R,
        "n_elements": eid - 1,
        "n_left": len(left_ids),
        "n_right": len(right_ids),
    }


def test_prepomax_dvizna_plosca_podmodel_immersed_iga():
    """Trimmed IGA + Kratos SmallDisplacementElement3D4N + add_linear_tie on Dvizna_plosca_podmodel.STEP."""
    t0 = time.perf_counter()
    cad = load_cad_file(_PREPOMAX_MODELS / "Dvizna_plosca_podmodel.STEP", heal=False)
    assert len(cad.solids) == 1
    solid = cad.solids[0]
    assert _count_free_edges(solid) == 0

    tess = tessellate(solid, deflection=0.05)
    rep = validate_shape(solid, tess, n_samples=5)
    assert rep.watertight

    lo, hi = cad.bbox
    box = BoxSpline(lo=tuple(lo - 0.5), hi=tuple(hi + 0.5), n_cells=(8, 4, 2), degree=1)
    q = 2
    qp_table, out_cells, _issues, _recs, diag = build_trimmed_table(
        tess.nodes, tess.tris, box, q=q, workers=1, cache_dir="off"
    )
    assert diag["n_fallback"] == 0

    brep_vol = brep_volume(solid)
    quad_vol = _quadrature_volume(box, qp_table, out_cells, q=q)
    vol_err = abs(quad_vol - brep_vol) / brep_vol
    assert vol_err < 0.02

    res = _solve_kratos_trimmed_elastic(
        box, qp_table, out_cells, E=210_000.0, nu=0.3, q=q, disp_x=0.1, use_linear_tie=True
    )
    U, R = res["U"], res["R"]
    assert np.all(np.isfinite(U))
    assert float(np.max(np.abs(U))) == pytest.approx(0.1, rel=1e-2) or float(np.max(np.abs(U))) >= 0.1 - 1e-9
    R_sum = R.sum(axis=0)
    assert float(np.linalg.norm(R_sum)) < 1e-6

    dt = time.perf_counter() - t0
    print(
        f"\n[Dvizna_plosca_podmodel] time={dt:.2f}s, n_fallback={diag['n_fallback']}, "
        f"vol_err={vol_err*100:.4f}%, n_elems={res['n_elements']}, "
        f"max_u={np.max(np.abs(U)):.4f}, |sum(R)|={np.linalg.norm(R_sum):.3e}"
    )


def test_prepomax_bearing_bracket_obb_immersed():
    """Axis-aligned vs PCA OBB active-cell fraction & Kratos elastic solve on Bearing_bracket.STEP."""
    t0 = time.perf_counter()
    cad = load_cad_file(_PREPOMAX_MODELS / "Bearing_bracket.STEP", heal=False)
    assert len(cad.solids) == 1
    solid = cad.solids[0]
    assert _count_free_edges(solid) == 0

    tess = tessellate(solid, deflection=0.1)
    rep = validate_shape(solid, tess, n_samples=5)
    assert rep.watertight
    brep_vol = brep_volume(solid)

    n_cells = (6, 8, 3)
    q = 2
    lo_aa, hi_aa = cad.bbox
    box_aa = BoxSpline(lo=tuple(lo_aa - 0.5), hi=tuple(hi_aa + 0.5), n_cells=n_cells, degree=1)
    qp_aa, out_aa, _, _, diag_aa = build_trimmed_table(
        tess.nodes, tess.tris, box_aa, q=q, workers=1, cache_dir="off"
    )
    assert diag_aa["n_fallback"] == 0
    vol_aa = _quadrature_volume(box_aa, qp_aa, out_aa, q=q)
    err_aa = abs(vol_aa - brep_vol) / brep_vol
    assert err_aa < 0.02

    frame, lo_obb, hi_obb = compute_pca_obb(tess.nodes, padding_frac=0.02)
    box_obb = BoxSpline(lo=lo_obb, hi=hi_obb, n_cells=n_cells, degree=1)
    pts_loc = frame.to_local_points(tess.nodes)
    qp_obb, out_obb, _, _, diag_obb = build_trimmed_table(
        pts_loc, tess.tris, box_obb, q=q, workers=1, cache_dir="off"
    )
    assert diag_obb["n_fallback"] == 0
    vol_obb = _quadrature_volume(box_obb, qp_obb, out_obb, q=q)
    err_obb = abs(vol_obb - brep_vol) / brep_vol
    assert err_obb < 0.02

    total_cells = int(np.prod(n_cells))
    active_frac_aa = 1.0 - len(out_aa) / total_cells
    active_frac_obb = 1.0 - len(out_obb) / total_cells
    assert active_frac_obb >= active_frac_aa - 1e-12

    res = _solve_kratos_trimmed_elastic(
        box_obb, qp_obb, out_obb, E=210_000.0, nu=0.3, q=q, disp_x=0.05, use_linear_tie=False
    )
    U, R = res["U"], res["R"]
    assert np.all(np.isfinite(U))
    assert float(np.max(np.abs(U))) >= 0.05 - 1e-9
    R_sum = R.sum(axis=0)
    assert float(np.linalg.norm(R_sum)) < 1e-6

    dt = time.perf_counter() - t0
    print(
        f"\n[Bearing_bracket] time={dt:.2f}s, AA active={active_frac_aa:.3f} (err={err_aa*100:.4f}%), "
        f"OBB active={active_frac_obb:.3f} (err={err_obb*100:.4f}%), "
        f"max_u={np.max(np.abs(U)):.4f}, |sum(R)|={np.linalg.norm(R_sum):.3e}"
    )


def test_prepomax_plates_bolt_5solid_partition():
    """5-solid watertight STEP partition on a shared coarse BoxSpline (Plates_bolt.STEP)."""
    t0 = time.perf_counter()
    cad = load_cad_file(_PREPOMAX_MODELS / "Plates_bolt.STEP", heal=False)
    assert len(cad.solids) == 5

    lo, hi = cad.bbox
    box = BoxSpline(lo=tuple(lo - 0.5), hi=tuple(hi + 0.5), n_cells=(8, 10, 6), degree=1)
    q = 2

    tot_brep = 0.0
    tot_quad = 0.0
    per_solid_errs = []
    for idx, solid in enumerate(cad.solids):
        assert _count_free_edges(solid) == 0
        tess = tessellate(solid, deflection=0.05)
        rep = validate_shape(solid, tess, n_samples=5)
        assert rep.watertight

        qp_table, out_cells, _, _, diag = build_trimmed_table(
            tess.nodes, tess.tris, box, q=q, workers=1, cache_dir="off"
        )
        assert diag["n_fallback"] == 0, f"solid[{idx}] had fallback cells"

        bvol = brep_volume(solid)
        qvol = _quadrature_volume(box, qp_table, out_cells, q=q)
        err = abs(qvol - bvol) / bvol
        assert err < 0.02, f"solid[{idx}] volume error {err:.4f} >= 2%"
        per_solid_errs.append(err)
        tot_brep += bvol
        tot_quad += qvol

    tot_err = abs(tot_quad - tot_brep) / tot_brep
    assert tot_err < 0.02

    dt = time.perf_counter() - t0
    errs_str = ", ".join(f"{e*100:.4f}%" for e in per_solid_errs)
    print(
        f"\n[Plates_bolt] time={dt:.2f}s, tot_err={tot_err*100:.4f}%, per_solid=[{errs_str}]"
    )


def _extract_pin_bottom_triangles(nodes: np.ndarray, tris: np.ndarray, y_max: float = 1.5):
    """Extract downward-facing triangles on the lower curved contact zone of the pin."""
    v0 = nodes[tris[:, 0]]
    v1 = nodes[tris[:, 1]]
    v2 = nodes[tris[:, 2]]
    cy = (v0[:, 1] + v1[:, 1] + v2[:, 1]) / 3.0
    cross = np.cross(v1 - v0, v2 - v0)
    norms = np.maximum(np.linalg.norm(cross, axis=1), 1e-12)
    ny = cross[:, 1] / norms
    keep = (cy < y_max) & (ny < -0.1)
    sub_tris = tris[keep]
    used, inv = np.unique(sub_tris.ravel(), return_inverse=True)
    return nodes[used], inv.reshape(-1, 3)


def test_prepomax_pin_on_plate_step_contact():
    """2-solid trimmed Kratos stiffness + deformable SSN contact with mu=0.3 on Pin_on_a_plate.STEP."""
    t0 = time.perf_counter()
    cad = load_cad_file(_PREPOMAX_MODELS / "Pin_on_a_plate.STEP", heal=False)
    assert len(cad.solids) == 2

    # Identify pin (y in [0, 5]) and plate (y in [-5, 0]) geometrically
    bboxes = [shape_bbox(s) for s in cad.solids]
    pin_idx = 0 if bboxes[0][1][1] > bboxes[1][1][1] else 1
    plate_idx = 1 - pin_idx
    pin_solid = cad.solids[pin_idx]
    plate_solid = cad.solids[plate_idx]

    t_pin = tessellate(pin_solid, deflection=0.05)
    t_plate = tessellate(plate_solid, deflection=0.05)

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    box_plate = BoxSpline(lo=(-15.0, -5.0, -5.0), hi=(15.0, 0.0, 0.0), n_cells=(6, 2, 2), degree=2)
    box_pin = BoxSpline(lo=(-10.0, 0.0, -5.0), hi=(0.0, 5.0, 0.0), n_cells=(3, 2, 2), degree=2)

    E_plate, E_pin, nu = 210_000.0, 210_000.0, 0.3
    op_plate = KratosBodyOperator(
        KM, SM, box_plate, E_plate, nu,
        element="SmallDisplacementElement3D4N",
        zero_fn=lambda x: (0, 1, 2) if abs(x[1] + 5.0) < 1e-9 else (),
        trim_soup=(t_plate.nodes, t_plate.tris),
        q=2,
        ghost_alpha=1e-8,
    )
    op_pin = KratosBodyOperator(
        KM, SM, box_pin, E_pin, nu,
        element="SmallDisplacementElement3D4N",
        presc_fn=lambda x: (0, 1, 2) if abs(x[1] - 5.0) < 1e-9 else (),
        trim_soup=(t_pin.nodes, t_pin.tris),
        q=2,
        ghost_alpha=1e-8,
    )
    assert op_plate.diag["n_fallback"] == 0
    assert op_pin.diag["n_fallback"] == 0

    K_full_m, _ = op_plate.eval(np.zeros(box_plate.n_dofs))
    K_full_s, _ = op_pin.eval(np.zeros(box_pin.n_dofs))

    fm = op_plate.free_mask
    fs = op_pin.free_mask
    pm = op_pin.presc_mask

    Kf_m = K_full_m[fm, :][:, fm].tocsr()
    Kf_s = K_full_s[fs, :][:, fs].tocsr()
    K_pf = K_full_s[pm, :][:, fs].tocsr()
    K_pp = K_full_s[pm, :][:, pm].tocsr()
    K_fp = K_full_s[fs, :][:, pm].tocsr()

    pin_bot_nodes, pin_bot_tris = _extract_pin_bottom_triangles(t_pin.nodes, t_pin.tris)
    assert len(pin_bot_tris) > 0
    X_s, w_q, B_s = sample_surface_triangles_gauss(pin_bot_nodes, pin_bot_tris, box_pin, order=2)
    master_face = DeformableBoxFaceMaster(box_plate, axis=1, side=1)

    # Normal press (-0.1 in y) + tangential slide (+0.03 in x) with Coulomb friction mu=0.3
    idx_p = np.flatnonzero(pm)
    comp_p = idx_p % 3
    u_presc = np.zeros(idx_p.size, dtype=float)
    u_presc[comp_p == 0] = 0.03
    u_presc[comp_p == 1] = -0.10

    u_fixed_s = np.zeros(box_pin.n_dofs, dtype=float)
    u_fixed_s[pm] = u_presc
    f_ext_s = -(K_fp @ u_presc)

    eps_n = 1e6
    res = solve_deformable_ssn_contact(
        Kf_m, Kf_s, X_s, w_q, B_s, master_face, eps_n,
        free_m=fm,
        free_s=fs,
        u_fixed_s=u_fixed_s,
        f_ext_s=f_ext_s,
        mu=0.3,
        eps_t=1e6,
        rtol=1e-12,
        use_pardiso=False,
    )

    assert res["n_active"] > 0
    F_N = float(-res["force_total"][1])
    assert F_N > 0.0

    R_p = K_pf @ res["u_free_s"] + K_pp @ u_presc + res["f_c_s_full"][pm]
    R_cut = np.array(
        [float(R_p[comp_p == 0].sum()), float(R_p[comp_p == 1].sum()), float(R_p[comp_p == 2].sum())],
        dtype=float,
    )
    abs_eq_err = float(np.linalg.norm(R_cut - res["force_total"]))
    rel_eq_err = abs_eq_err / max(1.0, float(np.linalg.norm(res["force_total"])))
    assert rel_eq_err < 1e-6

    dt = time.perf_counter() - t0
    print(
        f"\n[Pin_on_a_plate] time={dt:.2f}s, n_active={res['n_active']}/{len(w_q)}, "
        f"F_N={F_N:.4f}, force_total={res['force_total']}, abs_eq_err={abs_eq_err:.3e}, rel_eq_err={rel_eq_err:.3e}"
    )
