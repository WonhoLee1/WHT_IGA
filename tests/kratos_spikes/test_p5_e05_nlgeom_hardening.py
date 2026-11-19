# -*- coding: utf-8 -*-
"""P5 / E05: 3D Curved Surface Contact + NLGEOM + J2 Isotropic Hardening.

Tests covering :mod:`wht_iga.analyses.contact_surface`:
1. ``test_curved_triangle_surface_master_patch``: patch / isoparametric consistency
   (``sum(B_m) == 1``, ``B_m @ (X0 + U_m) == x_m`` to ``< 1e-12``), smooth
   angle-weighted normals, closest-point gap ``g_n``, and open-boundary exclusion
   on deformed curved 3D triangulated surfaces.
2. ``test_prepomax_contact_step_curved_surface_nlgeom``: 2-solid curved B-Rep
   contact on ``wht_iga/resources/Contact.STEP`` using trimmed
   :class:`NonlinearBodyOperator` (``HyperElastic3DLaw``, ``NLGEOM``) and
   :func:`solve_multibody_nlgeom_contact_step` with ``mu = 0.3``.
3. ``test_e05_parametric_nlgeom_j2_hardening``: E05 5-solid watertight B-Rep
   assembly (``hinge``, ``link``, 3-region partitioned ``pin``) under
   ``NonlinearBodyOperator`` + 2-point J2 isotropic hardening (``e05.json``
   ``weak_steel`` / ``strong_steel``) and 3-pair frictional (``mu = 0.3``)
   curved-surface contact up to yielding (``u_z = 1.0 mm``), compared against
   a purely hyperelastic NLGEOM baseline.
"""

from __future__ import annotations

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import psutil
import pytest

from tests.kratos_spikes.test_p5_e05_door_lock_precursor import build_door_lock_brep_parts
from wht_iga.analyses.contact import bspline_basis_matrix, sample_surface_triangles_gauss
from wht_iga.analyses.contact_surface import (
    DeformableTriangleSurfaceMaster,
    NonlinearBodyOperator,
    filter_triangles,
    release_solver_memory,
    solve_multibody_nlgeom_contact_step,
    subdivide_long_edges,
)
from wht_iga.analyses.thermal_plasticity import HardeningTable
from wht_iga.geometry.cad import shape_bbox, tessellate, validate_shape
from wht_iga.geometry.occ_io import _count_free_edges, load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_WHT_ROOT = Path(__file__).resolve().parents[2]
_CONTACT_STEP = _WHT_ROOT / "wht_iga" / "resources" / "Contact.STEP"
_E05_JSON = _WHT_ROOT / "tests" / "verification" / "data" / "e05.json"


def _rss_mb() -> float:
    """Return current process resident set size in MiB."""
    return psutil.Process(os.getpid()).memory_info().rss / (1024.0 * 1024.0)


def _sum_dirichlet_rf(body: NonlinearBodyOperator, r_full: np.ndarray) -> np.ndarray:
    """Sum reaction force vector ``(RF1, RF2, RF3)`` over Dirichlet DOFs."""
    mask = body.zero_mask | body.presc_mask
    dofs = np.arange(body.n_dofs, dtype=np.int64) % 3
    return np.array(
        [float(r_full[mask & (dofs == d)].sum()) for d in range(3)],
        dtype=np.float64,
    )


def _make_curved_grid_triangles(
    u_vals: np.ndarray,
    v_vals: np.ndarray,
    xyz_fn,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a structured triangle mesh from a parametric surface ``xyz_fn(u, v)``."""
    nu_pts = len(u_vals)
    nv_pts = len(v_vals)
    nodes = np.empty((nu_pts * nv_pts, 3), dtype=np.float64)
    for iu, u in enumerate(u_vals):
        for jv, v in enumerate(v_vals):
            nodes[iu * nv_pts + jv] = xyz_fn(float(u), float(v))

    tris = []
    for iu in range(nu_pts - 1):
        for jv in range(nv_pts - 1):
            a = iu * nv_pts + jv
            b = (iu + 1) * nv_pts + jv
            c = (iu + 1) * nv_pts + (jv + 1)
            d = iu * nv_pts + (jv + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    return nodes, np.asarray(tris, dtype=np.int64)


# --------------------------------------------------------------------------------------
# Test 1: Curved DeformableTriangleSurfaceMaster patch & closest-point projection
# --------------------------------------------------------------------------------------


def test_curved_triangle_surface_master_patch() -> None:
    """Verify DeformableTriangleSurfaceMaster B-spline consistency, smooth normals, and gaps."""
    box = BoxSpline(lo=(-8.0, -5.0, -5.0), hi=(8.0, 5.0, 5.0), n_cells=(4, 4, 4), degree=2)
    r_cyl = 4.0

    # Cylindrical half-shell r = 4.0 around X-axis, outward normal in +YZ
    x_vals = np.linspace(-7.0, 7.0, 5)
    th_vals = np.linspace(-0.4 * np.pi, 0.4 * np.pi, 25)
    nodes, tris = _make_curved_grid_triangles(
        x_vals,
        th_vals,
        lambda x, th: (x, r_cyl * np.sin(th), r_cyl * np.cos(th)),
    )
    # Subdivide long straight axial edges so B-spline deformation varies along X
    nodes, tris = subdivide_long_edges(nodes, tris, max_edge=2.0)

    master = DeformableTriangleSurfaceMaster(
        box,
        nodes,
        tris,
        outward_flip=False,
        smooth_normals=True,
        exclude_boundary=True,
    )

    # Smooth non-affine deformation field evaluated at control points
    cp = box.control_points()
    U_m = np.column_stack(
        [
            0.04 * cp[:, 1] - 0.02 * cp[:, 2],
            0.06 * np.sin(0.15 * cp[:, 0]) + 0.03 * cp[:, 2],
            0.12 - 0.04 * (cp[:, 0] / 8.0) ** 2 + 0.02 * cp[:, 1],
        ]
    )

    # Verify reference B-spline coordinate identity on undeformed surface
    X0_eval = master.B_verts @ cp
    assert np.max(np.abs(X0_eval - master.tri_nodes)) < 1e-12

    # Build slave query points near deformed surface (penetrating, on-surface, gapped)
    V_def = master.tri_nodes + master.B_verts @ U_m
    tri_centers = V_def[master.tris].mean(axis=1)
    v0 = V_def[master.tris[:, 0]]
    v1 = V_def[master.tris[:, 1]]
    v2 = V_def[master.tris[:, 2]]
    fn_def = np.cross(v1 - v0, v2 - v0)
    fn_def /= np.linalg.norm(fn_def, axis=1, keepdims=True)

    # Pick interior triangles away from patch boundary
    inner_ids = np.arange(0, len(master.tris), max(1, len(master.tris) // 24))
    inner_ids = [
        k
        for k in inner_ids
        if abs(tri_centers[k, 0]) < 5.0 and abs(tri_centers[k, 1]) < 3.0
    ]
    assert len(inner_ids) >= 8

    offsets = np.array([-0.15, 0.0, +0.20], dtype=np.float64)
    x_queries = []
    expected_gn = []
    for k in inner_ids:
        for d_off in offsets:
            x_queries.append(tri_centers[k] + d_off * fn_def[k])
            expected_gn.append(d_off)
    x_s = np.asarray(x_queries, dtype=np.float64)
    expected_gn_arr = np.asarray(expected_gn, dtype=np.float64)

    g_n, n_m, x_m, ab_m, B_m = master.evaluate(x_s, U_m)

    # 1. Partition of unity & exact isoparametric B-spline interpolation
    row_sums = np.asarray(B_m.sum(axis=1)).ravel()
    assert np.max(np.abs(row_sums - 1.0)) < 1e-12
    x_m_from_B = B_m @ (cp + U_m)
    assert np.max(np.abs(x_m_from_B - x_m)) < 1e-12

    # 2. Barycentric coordinates in [0, 1] and g_n == dot(x_s - x_m, n_m)
    b0 = 1.0 - ab_m[:, 0] - ab_m[:, 1]
    assert np.all(b0 >= -1e-12) and np.all(ab_m >= -1e-12)
    assert np.all(b0 <= 1.0 + 1e-12) and np.all(ab_m <= 1.0 + 1e-12)
    gn_recomputed = np.sum((x_s - x_m) * n_m, axis=1)
    assert np.max(np.abs(g_n - gn_recomputed)) < 1e-12

    # 3. Signed gap signs & accuracy on penetrating / on-surface / gapped points
    pen_mask = expected_gn_arr < -0.05
    gap_mask = expected_gn_arr > 0.05
    assert np.all(g_n[pen_mask] < 0.0)
    assert np.all(g_n[gap_mask] > 0.0)
    assert np.max(np.abs(g_n - expected_gn_arr)) < 5e-3

    # 4. Open-boundary exclusion check: overhanging point outside x > 7.0 stays inactive (g_n > 0)
    x_overhang = np.array([[7.8, 0.0, r_cyl - 0.3]], dtype=np.float64) + (
        bspline_basis_matrix(box, np.array([[7.0, 0.0, r_cyl]])) @ U_m
    )
    g_over, _, _, _, _ = master.evaluate(x_overhang, U_m)
    assert g_over[0] > 0.0

    # 5. Pure rigid-body shift on curved dome verifies exact analytical gap preservation
    box_dome = BoxSpline(lo=(-4.0, -4.0, -1.0), hi=(4.0, 4.0, 3.0), n_cells=(4, 4, 2), degree=2)
    uv = np.linspace(-3.0, 3.0, 21)
    d_nodes, d_tris = _make_curved_grid_triangles(
        uv, uv, lambda x, y: (x, y, 0.08 * (x * x + y * y))
    )
    master_dome = DeformableTriangleSurfaceMaster(box_dome, d_nodes, d_tris, smooth_normals=True)
    shift = np.array([0.35, -0.25, 0.40], dtype=np.float64)
    U_rigid = np.broadcast_to(shift, (box_dome.n_ctrl, 3)).copy()
    x_dome_q = np.array([[0.0, 0.0, -0.10], [0.0, 0.0, 0.25]], dtype=np.float64) + shift
    g_dome, n_dome, x_dome_m, _, B_dome = master_dome.evaluate(x_dome_q, U_rigid)
    assert np.max(np.abs(B_dome @ (box_dome.control_points() + U_rigid) - x_dome_m)) < 1e-12
    assert g_dome[0] == pytest.approx(-0.10, abs=1e-12)
    assert g_dome[1] == pytest.approx(+0.25, abs=1e-3)
    assert np.allclose(n_dome[0], [0.0, 0.0, 1.0], atol=1e-12)
    assert np.allclose(n_dome[1], [0.0, 0.0, 1.0], atol=2e-3)

    print(
        f"\n[Curved Triangle Master Patch] M={len(x_s)}, "
        f"row_sum_err={np.max(np.abs(row_sums - 1.0)):.2e}, "
        f"xm_interp_err={np.max(np.abs(x_m_from_B - x_m)):.2e}, "
        f"gap_max_err={np.max(np.abs(g_n - expected_gn_arr)):.2e}, "
        f"overhang_gn={g_over[0]:+.4f}"
    )


# --------------------------------------------------------------------------------------
# Test 2: Contact.STEP 3D curved-surface NLGEOM contact (HyperElastic3DLaw, mu=0.3)
# --------------------------------------------------------------------------------------


def test_prepomax_contact_step_curved_surface_nlgeom() -> None:
    """Run 3-step NLGEOM frictional curved-surface contact on Contact.STEP."""
    if not _CONTACT_STEP.is_file():
        pytest.skip(f"Contact.STEP not found: {_CONTACT_STEP}")

    t0 = time.perf_counter()
    cad = load_cad_file(_CONTACT_STEP, heal=False)
    assert len(cad.solids) == 2

    s_hook, s_cyl = cad.solids[0], cad.solids[1]
    assert _count_free_edges(s_hook) == 0 and _count_free_edges(s_cyl) == 0
    t_hook = tessellate(s_hook, deflection=0.05)
    t_cyl = tessellate(s_cyl, deflection=0.05)
    assert validate_shape(s_hook, t_hook, n_samples=5).watertight
    assert validate_shape(s_cyl, t_cyl, n_samples=5).watertight

    bb0 = shape_bbox(s_hook)
    bb1 = shape_bbox(s_cyl)
    box_hook = BoxSpline(lo=tuple(bb0[0] - 0.25), hi=tuple(bb0[1] + 0.25), n_cells=(5, 6, 3), degree=2)
    box_cyl = BoxSpline(lo=tuple(bb1[0] - 0.25), hi=tuple(bb1[1] + 0.25), n_cells=(5, 5, 4), degree=2)

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    E, nu = 210_000.0, 0.3

    op_hook = None
    op_cyl = None
    try:
        op_hook = NonlinearBodyOperator(
            KM,
            mods,
            box_hook,
            [{"trim_soup": (t_hook.nodes, t_hook.tris), "E": E, "nu": nu}],
            presc_fn=lambda x: (0, 1, 2) if x[0] <= -11.5 else (),
            q=2,
            ghost_alpha=0.0,
        )
        op_cyl = NonlinearBodyOperator(
            KM,
            mods,
            box_cyl,
            [{"trim_soup": (t_cyl.nodes, t_cyl.tris), "E": E, "nu": nu}],
            zero_fn=lambda x: (0, 1, 2) if x[0] >= 16.5 else (),
            q=2,
            ghost_alpha=0.0,
        )
        assert op_hook.diag["n_fallback"] == 0
        assert op_cyl.diag["n_fallback"] == 0

        # Master: outer curved cylindrical wall of solid[1] (r8 = 10 centered at (8, 0))
        v1 = t_cyl.nodes[t_cyl.tris]
        c1 = v1.mean(axis=1)
        n1 = np.cross(v1[:, 1] - v1[:, 0], v1[:, 2] - v1[:, 0])
        n1 /= np.linalg.norm(n1, axis=1, keepdims=True)
        r8_1 = np.sqrt((c1[:, 0] - 8.0) ** 2 + c1[:, 1] ** 2)
        m_nodes, m_tris = filter_triangles(
            t_cyl.nodes, t_cyl.tris, (r8_1 > 9.5) & (np.abs(n1[:, 2]) < 0.1)
        )
        m_nodes, m_tris = subdivide_long_edges(m_nodes, m_tris, max_edge=3.0)
        master_cyl = DeformableTriangleSurfaceMaster(
            box_cyl, m_nodes, m_tris, smooth_normals=True, exclude_boundary=True
        )

        # Slave: inner curved hook arms of solid[0] facing solid[1]
        v0 = t_hook.nodes[t_hook.tris]
        c0 = v0.mean(axis=1)
        n0 = np.cross(v0[:, 1] - v0[:, 0], v0[:, 2] - v0[:, 0])
        n0 /= np.linalg.norm(n0, axis=1, keepdims=True)
        r8_0 = np.sqrt((c0[:, 0] - 8.0) ** 2 + c0[:, 1] ** 2)
        s_nodes, s_tris = filter_triangles(
            t_hook.nodes, t_hook.tris, (r8_0 < 12.0) & (c0[:, 0] > 0.0) & (np.abs(n0[:, 2]) < 0.1)
        )
        s_nodes, s_tris = subdivide_long_edges(s_nodes, s_tris, max_edge=3.0)
        X_s, w_q, B_s = sample_surface_triangles_gauss(s_nodes, s_tris, box_hook, order=2)
        assert len(X_s) > 0

        pairs = [
            {
                "name": "cyl_to_hook",
                "m_idx": 1,
                "s_idx": 0,
                "master_face": master_cyl,
                "X_s": X_s,
                "w_q": w_q,
                "B_s": B_s,
            }
        ]

        u_hook = np.zeros(op_hook.n_dofs, dtype=np.float64)
        u_cyl = np.zeros(op_cyl.n_dofs, dtype=np.float64)
        slip_states = None
        pm = np.flatnonzero(op_hook.presc_mask)

        print("\n[Contact.STEP Curved-Surface NLGEOM Contact, mu=0.3]")
        last_stats = None
        for step, ux in enumerate((0.7, 1.0, 1.3), start=1):
            u_hook[pm] = 0.0
            u_hook[pm[pm % 3 == 0]] = ux
            res = solve_multibody_nlgeom_contact_step(
                [op_hook, op_cyl],
                [u_hook, u_cyl],
                pairs,
                eps_n=5e5,
                step_time=float(step),
                mu=0.3,
                eps_t=5e5,
                slip_states=slip_states,
                outer_tol=1e-8,
                rtol=1e-8,
                use_pardiso=True,
            )
            u_hook, u_cyl = res["u_full_list"]
            slip_states = res["slip_states"]

            ps = res["pair_stats"]["cyl_to_hook"]
            last_stats = ps
            rf_hook = _sum_dirichlet_rf(op_hook, res["f_int_list"][0] + res["f_c_full_list"][0])
            rf_cyl = _sum_dirichlet_rf(op_cyl, res["f_int_list"][1] + res["f_c_full_list"][1])
            eq_rel = float(np.linalg.norm(rf_hook + rf_cyl) / max(np.linalg.norm(rf_hook), 1.0))

            print(
                f"  step {step} (ux={ux:.1f}mm): outer={res['n_outer']}, "
                f"n_act={ps['n_active']:3d}, F_N={ps['normal_pressure_integral']:7.2f} N, "
                f"|Ft|={np.linalg.norm(ps['F_t']):6.2f} N, "
                f"RF_x=({rf_hook[0]:+7.2f}, {rf_cyl[0]:+7.2f}) N, eq_rel={eq_rel:.2e}"
            )
            assert ps["n_active"] > 0
            assert ps["normal_pressure_integral"] > 0.0
            assert float(np.linalg.norm(ps["F_t"])) > 0.0
            assert eq_rel < 1e-6

        assert last_stats is not None
        assert last_stats["normal_pressure_integral"] > 400.0
        elapsed = time.perf_counter() - t0
        rss = _rss_mb()
        print(f"  elapsed={elapsed:.2f}s, RSS={rss:.1f} MB")
        assert rss < 1536.0
    finally:
        if op_hook is not None:
            op_hook.close()
        if op_cyl is not None:
            op_cyl.close()
        release_solver_memory()
        gc.collect()


# --------------------------------------------------------------------------------------
# Test 3: E05 5-solid B-Rep + NLGEOM + 2-point J2 hardening + 3-pair curved contact
# --------------------------------------------------------------------------------------


def _extract_e05_contact_surfaces(parts: dict, box_h: BoxSpline, box_p: BoxSpline, box_l: BoxSpline):
    """Build triangulated master surfaces and slave Gauss rules for E05's 3 contact pairs."""
    tess = {k: tessellate(s, deflection=0.2) for k, s in parts.items()}

    # Pin outer cylindrical surface (x in [-16, 16], |nx| < 0.5)
    vp = tess["pin_cyl"].nodes[tess["pin_cyl"].tris]
    np_ = np.cross(vp[:, 1] - vp[:, 0], vp[:, 2] - vp[:, 0])
    np_ /= np.linalg.norm(np_, axis=1, keepdims=True)
    pm_nodes, pm_tris = filter_triangles(
        tess["pin_cyl"].nodes, tess["pin_cyl"].tris, np.abs(np_[:, 0]) < 0.5
    )
    pm_nodes, pm_tris = subdivide_long_edges(pm_nodes, pm_tris, max_edge=4.0)
    master_pin = DeformableTriangleSurfaceMaster(
        box_p, pm_nodes, pm_tris, smooth_normals=True, exclude_boundary=True, max_penetration=2.0
    )

    # Hinge lower pin bore (slave for pin_to_hinge_contact)
    vh = tess["hinge"].nodes[tess["hinge"].tris]
    ch = vh.mean(axis=1)
    nh = np.cross(vh[:, 1] - vh[:, 0], vh[:, 2] - vh[:, 0])
    nh /= np.linalg.norm(nh, axis=1, keepdims=True)
    rh = np.sqrt(ch[:, 1] ** 2 + ch[:, 2] ** 2)
    sh_nodes, sh_tris = filter_triangles(
        tess["hinge"].nodes,
        tess["hinge"].tris,
        (rh < 4.3) & (ch[:, 2] < -0.2) & (np.abs(nh[:, 0]) < 0.5),
    )
    sh_nodes, sh_tris = subdivide_long_edges(sh_nodes, sh_tris, max_edge=4.0)
    X_ph, w_ph, B_ph = sample_surface_triangles_gauss(sh_nodes, sh_tris, box_h, order=2)

    # Link upper pin bore (slave for pin_to_link_contact)
    vl = tess["link"].nodes[tess["link"].tris]
    cl = vl.mean(axis=1)
    nl = np.cross(vl[:, 1] - vl[:, 0], vl[:, 2] - vl[:, 0])
    nl /= np.linalg.norm(nl, axis=1, keepdims=True)
    rl = np.sqrt(cl[:, 1] ** 2 + cl[:, 2] ** 2)
    sl_nodes, sl_tris = filter_triangles(
        tess["link"].nodes,
        tess["link"].tris,
        (rl < 4.3) & (cl[:, 2] > 0.2) & (np.abs(nl[:, 0]) < 0.5),
    )
    sl_nodes, sl_tris = subdivide_long_edges(sl_nodes, sl_tris, max_edge=4.0)
    X_pl, w_pl, B_pl = sample_surface_triangles_gauss(sl_nodes, sl_tris, box_l, order=2)

    # Hinge ear inner walls (master) vs link saddle side walls (slave) for hinge_to_link_contact
    hm_mask = (np.abs(ch[:, 0]) > 10.5) & (np.abs(ch[:, 0]) < 11.5) & (ch[:, 1] > -10.0) & (np.abs(nh[:, 0]) > 0.5)
    hm_nodes, hm_tris = filter_triangles(tess["hinge"].nodes, tess["hinge"].tris, hm_mask)
    hm_nodes, hm_tris = subdivide_long_edges(hm_nodes, hm_tris, max_edge=4.0)
    master_hinge = DeformableTriangleSurfaceMaster(
        box_h, hm_nodes, hm_tris, smooth_normals=True, exclude_boundary=True, max_distance=2.0
    )

    hl_mask = (np.abs(cl[:, 0]) > 9.5) & (cl[:, 1] < 8.5) & (np.abs(nl[:, 0]) > 0.5)
    hl_nodes, hl_tris = filter_triangles(tess["link"].nodes, tess["link"].tris, hl_mask)
    hl_nodes, hl_tris = subdivide_long_edges(hl_nodes, hl_tris, max_edge=4.0)
    X_hl, w_hl, B_hl = sample_surface_triangles_gauss(hl_nodes, hl_tris, box_l, order=2)

    pairs = [
        {
            "name": "pin_to_hinge_contact",
            "m_idx": 1,
            "s_idx": 0,
            "master_face": master_pin,
            "X_s": X_ph,
            "w_q": w_ph,
            "B_s": B_ph,
        },
        {
            "name": "pin_to_link_contact",
            "m_idx": 1,
            "s_idx": 2,
            "master_face": master_pin,
            "X_s": X_pl,
            "w_q": w_pl,
            "B_s": B_pl,
        },
        {
            "name": "hinge_to_link_contact",
            "m_idx": 0,
            "s_idx": 2,
            "master_face": master_hinge,
            "X_s": X_hl,
            "w_q": w_hl,
            "B_s": B_hl,
        },
    ]
    return tess, pairs


def _run_e05_multistep(
    KM,
    mods,
    parts: dict,
    tess: dict,
    box_h: BoxSpline,
    box_p: BoxSpline,
    box_l: BoxSpline,
    pairs: list[dict],
    *,
    E: float,
    nu: float,
    weak_h: HardeningTable | None,
    strong_h: HardeningTable | None,
    uz_steps: tuple[float, ...],
    mu: float = 0.3,
) -> dict:
    """Execute multi-step coupled 3-part NLGEOM contact on E05 with or without J2 hardening."""
    op_h = None
    op_p = None
    op_l = None
    try:
        op_h = NonlinearBodyOperator(
            KM,
            mods,
            box_h,
            [{"trim_soup": (tess["hinge"].nodes, tess["hinge"].tris), "E": E, "nu": nu, "hardening": weak_h}],
            presc_fn=lambda x: (0, 1, 2) if x[2] >= 33.0 else (),
            q=2,
            ghost_alpha=0.0,
        )
        # Tiny X/Y/Z rigid-mode stabilization on pin (Z reg 1e-6 keeps Z balance < 1e-6)
        pin_reg = np.zeros(box_p.n_dofs, dtype=np.float64)
        pin_reg[0::3] = 1e-2
        pin_reg[1::3] = 1e-2
        pin_reg[2::3] = 1e-6
        op_p = NonlinearBodyOperator(
            KM,
            mods,
            box_p,
            [
                {"trim_soup": (tess["pin_cyl"].nodes, tess["pin_cyl"].tris), "E": E, "nu": nu, "hardening": strong_h},
                {"trim_soup": (tess["pin_cap_l"].nodes, tess["pin_cap_l"].tris), "E": E, "nu": nu, "hardening": weak_h},
                {"trim_soup": (tess["pin_cap_r"].nodes, tess["pin_cap_r"].tris), "E": E, "nu": nu, "hardening": weak_h},
            ],
            q=2,
            ghost_alpha=0.0,
            stiffness_reg=pin_reg,
        )
        op_l = NonlinearBodyOperator(
            KM,
            mods,
            box_l,
            [{"trim_soup": (tess["link"].nodes, tess["link"].tris), "E": E, "nu": nu, "hardening": weak_h}],
            zero_fn=lambda x: (0, 1, 2) if x[1] >= 22.0 else (),
            q=2,
            ghost_alpha=0.0,
        )
        assert op_h.diag["n_fallback"] == 0
        assert op_p.diag["n_fallback"] == 0
        assert op_l.diag["n_fallback"] == 0

        u_h = np.zeros(op_h.n_dofs, dtype=np.float64)
        u_p = np.zeros(op_p.n_dofs, dtype=np.float64)
        u_l = np.zeros(op_l.n_dofs, dtype=np.float64)
        slip_states = None
        pm = np.flatnonzero(op_h.presc_mask)

        history = []
        for step, uz in enumerate(uz_steps, start=1):
            u_h[pm] = 0.0
            u_h[pm[pm % 3 == 2]] = uz
            res = solve_multibody_nlgeom_contact_step(
                [op_h, op_p, op_l],
                [u_h, u_p, u_l],
                pairs,
                eps_n=2e5,
                step_time=float(step),
                mu=mu,
                eps_t=4e4,
                slip_states=slip_states,
                outer_tol=1e-8,
                rtol=1e-8,
                use_pardiso=True,
            )
            u_h, u_p, u_l = res["u_full_list"]
            slip_states = res["slip_states"]

            rf_h = _sum_dirichlet_rf(op_h, res["f_int_list"][0] + res["f_c_full_list"][0])
            rf_l = _sum_dirichlet_rf(op_l, res["f_int_list"][2] + res["f_c_full_list"][2])
            bal_z = abs(rf_h[2] + rf_l[2]) / max(abs(rf_l[2]), 1.0)

            history.append(
                {
                    "step": step,
                    "uz": float(uz),
                    "n_outer": int(res["n_outer"]),
                    "outer_rel": float(res["outer_rel"]),
                    "peeq_max_list": [float(v) for v in res["peeq_max_list"]],
                    "rf_hinge_z": float(rf_h[2]),
                    "rf_link_z": float(rf_l[2]),
                    "bal_z": float(bal_z),
                    "n_active_ph": int(res["pair_stats"]["pin_to_hinge_contact"]["n_active"]),
                    "n_active_pl": int(res["pair_stats"]["pin_to_link_contact"]["n_active"]),
                    "n_active_hl": int(res["pair_stats"]["hinge_to_link_contact"]["n_active"]),
                }
            )
        return {"history": history, "final_peeq_max": max(history[-1]["peeq_max_list"])}
    finally:
        if op_h is not None:
            op_h.close()
        if op_p is not None:
            op_p.close()
        if op_l is not None:
            op_l.close()
        release_solver_memory()
        gc.collect()


def test_e05_parametric_nlgeom_j2_hardening() -> None:
    """E05 5-solid B-Rep + NLGEOM + 2-point J2 hardening + 3-pair mu=0.3 contact."""
    t0 = time.perf_counter()
    spec = json.loads(_E05_JSON.read_text(encoding="utf-8"))
    E = float(spec["materials"]["weak_steel"]["E_MPa"])
    nu = float(spec["materials"]["weak_steel"]["nu"])
    assert E == 200_000.0 and nu == 0.29

    # 2-point linearly extrapolated isotropic hardening from e05.json
    weak_h = HardeningTable((0.0, 10.0), ((1500.0, 150_000.0),))
    strong_h = HardeningTable((0.0, 10.0), ((2000.0, 200_000.0),))

    parts = build_door_lock_brep_parts()
    bb_h = shape_bbox(parts["hinge"])
    box_h = BoxSpline(lo=tuple(bb_h[0] - 0.5), hi=tuple(bb_h[1] + 0.5), n_cells=(6, 4, 5), degree=2)
    p_lo = np.min([shape_bbox(parts[k])[0] for k in ("pin_cyl", "pin_cap_l", "pin_cap_r")], axis=0)
    p_hi = np.max([shape_bbox(parts[k])[1] for k in ("pin_cyl", "pin_cap_l", "pin_cap_r")], axis=0)
    box_p = BoxSpline(lo=tuple(p_lo - 0.25), hi=tuple(p_hi + 0.25), n_cells=(8, 3, 3), degree=2)
    bb_l = shape_bbox(parts["link"])
    box_l = BoxSpline(lo=tuple(bb_l[0] - 0.5), hi=tuple(bb_l[1] + 0.5), n_cells=(4, 5, 3), degree=2)

    tess, pairs = _extract_e05_contact_surfaces(parts, box_h, box_p, box_l)
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]

    uz_steps = (0.25, 0.50, 0.75, 1.00)
    res_j2 = _run_e05_multistep(
        KM, mods, parts, tess, box_h, box_p, box_l, pairs,
        E=E, nu=nu, weak_h=weak_h, strong_h=strong_h, uz_steps=uz_steps, mu=0.3,
    )
    res_el = _run_e05_multistep(
        KM, mods, parts, tess, box_h, box_p, box_l, pairs,
        E=E, nu=nu, weak_h=None, strong_h=None, uz_steps=uz_steps, mu=0.3,
    )

    print("\n[E05 Parametric 5-Solid NLGEOM + J2 Hardening vs HyperElastic, mu=0.3]")
    prev_peeq = 0.0
    for r_j2, r_el in zip(res_j2["history"], res_el["history"]):
        pmax = max(r_j2["peeq_max_list"])
        print(
            f"  step {r_j2['step']} uz={r_j2['uz']:.2f}mm: outer={r_j2['n_outer']}, "
            f"peeq_max={pmax:.6f}, RF3_j2={r_j2['rf_link_z']/1e3:+7.4f} kN, "
            f"RF3_el={r_el['rf_link_z']/1e3:+7.4f} kN, bal_z={r_j2['bal_z']:.2e}, "
            f"act(ph/pl/hl)={r_j2['n_active_ph']}/{r_j2['n_active_pl']}/{r_j2['n_active_hl']}"
        )
        assert r_j2["outer_rel"] <= 1e-8
        assert r_el["outer_rel"] <= 1e-8
        assert r_j2["rf_hinge_z"] > 0.0
        assert r_j2["rf_link_z"] < 0.0
        assert r_j2["bal_z"] < 1e-6
        assert r_el["bal_z"] < 1e-6
        assert r_j2["n_active_ph"] > 0 and r_j2["n_active_pl"] > 0
        assert pmax >= prev_peeq - 1e-15
        prev_peeq = pmax

    final_j2 = res_j2["history"][-1]
    final_el = res_el["history"][-1]
    assert res_j2["final_peeq_max"] > 0.0
    assert abs(final_j2["rf_link_z"]) < abs(final_el["rf_link_z"])

    elapsed = time.perf_counter() - t0
    rss = _rss_mb()
    print(
        f"  ratio |RF3_j2|/|RF3_el|={abs(final_j2['rf_link_z'])/abs(final_el['rf_link_z']):.4f}, "
        f"elapsed={elapsed:.2f}s, RSS={rss:.1f} MB"
    )
    assert rss < 1536.0
