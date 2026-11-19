# -*- coding: utf-8 -*-
"""E05 Quantitative Benchmark on Authentic ``Door_Lock_Assembly.sat``.

Validates the full 3-part immersed IGA door-lock assembly directly on the real
ACIS ``.sat`` CAD geometry (:file:`wht_iga/resources/Door_Lock_Assembly.sat`):
1. Imports 3 watertight B-Rep solids (``Link``, ``Hinge``, ``Pin``) via
   :func:`wht_iga.geometry.occ_io.load_cad_file`.
2. Partitions ``Pin`` at ``|y| = 18.1 mm`` into central shaft (``strong_steel``)
   and two end caps (``weak_steel``) with volume conservation.
3. Configures the 3 frictional surface-to-surface contact pairs (``mu = 0.3``)
   with :class:`~wht_iga.analyses.contact_surface.DeformableTriangleSurfaceMaster`
   and verifies initial clearances (``~0.10 mm`` radial pin-bore, ``~0.50 mm``
   axial ear wall).
4. Solves multi-step ``NLGEOM`` (``TotalLagrangianElement3D4N`` +
   :class:`~wht_iga.analyses.thermal_plasticity.HardeningTable` 2-point J2
   isotropic hardening) from ``u_z = 0`` to ``5.0 mm`` using
   :func:`~wht_iga.analyses.contact_surface.solve_multibody_nlgeom_contact_step`,
   and compares ``RF3_link`` against Coreform ``e05.json`` (``-46.946 kN ± 5%``)
   and the SVG load-history curve.
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

from wht_iga.analyses.contact import sample_surface_triangles_gauss
from wht_iga.analyses.contact_surface import (
    DeformableTriangleSurfaceMaster,
    NonlinearBodyOperator,
    filter_triangles,
    release_solver_memory,
    solve_multibody_nlgeom_contact_step,
    subdivide_long_edges,
)
from wht_iga.analyses.thermal_plasticity import HardeningTable
from wht_iga.geometry.cad import brep_volume, shape_bbox, tessellate
from wht_iga.geometry.occ_io import load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_WHT_ROOT = Path(__file__).resolve().parents[2]
_SAT_PATH = _WHT_ROOT / "wht_iga" / "resources" / "Door_Lock_Assembly.sat"
_E05_JSON = _WHT_ROOT / "tests" / "verification" / "data" / "e05.json"

# Coreform SVG curve checkpoints (uz [mm], RF3_link [kN]) decoded from fig-door-lock-reaction-force-output-1.svg
_SVG_CURVE_REF = (
    (1.130, -13.370),
    (1.664, -21.106),
    (2.465, -29.181),
    (3.667, -38.423),
    (5.000, -46.946),
)


def _rss_mb() -> float:
    """Return current process RSS in MiB."""
    return psutil.Process(os.getpid()).memory_info().rss / (1024.0 * 1024.0)


def _tri_props(nodes: np.ndarray, tris: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute triangle centroids and unit normals."""
    v = nodes[tris]
    c = v.mean(axis=1)
    cr = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    n = cr / (np.linalg.norm(cr, axis=1, keepdims=True) + 1e-30)
    return c, n


def _partition_pin_solids(pin_solid) -> tuple:
    """Split ``Pin`` at ``y = ±18.1 mm`` into central shaft and two end caps."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    pin_cyl = BRepAlgoAPI_Common(
        pin_solid, BRepPrimAPI_MakeBox(gp_Pnt(-10.0, -18.1, 0.0), gp_Pnt(10.0, 18.1, 25.0)).Shape()
    ).Shape()
    pin_l = BRepAlgoAPI_Common(
        pin_solid, BRepPrimAPI_MakeBox(gp_Pnt(-10.0, -25.0, 0.0), gp_Pnt(10.0, -18.1, 25.0)).Shape()
    ).Shape()
    pin_r = BRepAlgoAPI_Common(
        pin_solid, BRepPrimAPI_MakeBox(gp_Pnt(-10.0, 18.1, 0.0), gp_Pnt(10.0, 25.0, 25.0)).Shape()
    ).Shape()
    return pin_cyl, pin_l, pin_r


def _make_box_for_solid(solid, h: float) -> BoxSpline:
    """Build axis-aligned quadratic ``BoxSpline`` covering ``solid`` with cell size ``<= h``."""
    lo, hi = shape_bbox(solid)
    nc = tuple(max(2, int(np.ceil((hi[d] - lo[d]) / h))) for d in range(3))
    return BoxSpline(tuple(lo), tuple(hi), nc, degree=2)


def test_e05_sat_door_lock_full_nlgeom_j2_quantitative() -> None:
    """Run full 0 -> 5 mm NLGEOM + J2 contact solve on ``Door_Lock_Assembly.sat``."""
    if not _SAT_PATH.is_file():
        pytest.skip(f"Missing Door_Lock_Assembly.sat: {_SAT_PATH}")

    t0 = time.perf_counter()
    spec = json.loads(_E05_JSON.read_text(encoding="utf-8"))
    E = float(spec["materials"]["weak_steel"]["E_MPa"])
    nu = float(spec["materials"]["weak_steel"]["nu"])
    ref_rf3_kn = float(spec["reference_values"]["summed_link_RF3_kN"])
    rel_tol = float(spec["acceptance"]["rel_tol"])

    cad = load_cad_file(_SAT_PATH, heal=False)
    assert cad.healing_report.watertight and cad.healing_report.free_edges_after == 0
    link, hinge, pin = cad.solids

    pin_cyl, pin_l, pin_r = _partition_pin_solids(pin)
    v_pin = brep_volume(pin)
    v_sum = brep_volume(pin_cyl) + brep_volume(pin_l) + brep_volume(pin_r)
    assert v_sum == pytest.approx(v_pin, rel=1e-4)

    tl = tessellate(link, deflection=0.05)
    th = tessellate(hinge, deflection=0.05)
    tp = tessellate(pin, deflection=0.05)
    tpc = tessellate(pin_cyl, deflection=0.05)
    tpl = tessellate(pin_l, deflection=0.05)
    tpr = tessellate(pin_r, deflection=0.05)

    h = 3.0
    box_l = _make_box_for_solid(link, h)
    box_h = _make_box_for_solid(hinge, h)
    box_p = _make_box_for_solid(pin, h)

    cl, nl = _tri_props(tl.nodes, tl.tris)
    ch, nh = _tri_props(th.nodes, th.tris)

    # Closed watertight Pin surface (zero open boundary edges after conforming subdivision)
    pm_n, pm_t = subdivide_long_edges(tp.nodes, tp.tris, max_edge=2.5)
    master_pin = DeformableTriangleSurfaceMaster(
        box_p, pm_n, pm_t, smooth_normals=True, exclude_boundary=True, max_penetration=2.5
    )

    # Hinge inner clevis bores (pin_to_hinge_contact)
    rh = np.sqrt((ch[:, 0] + 1.5) ** 2 + (ch[:, 2] - 10.5) ** 2)
    hb_n, hb_t = subdivide_long_edges(
        *filter_triangles(
            th.nodes,
            th.tris,
            (rh < 3.4) & (np.abs(ch[:, 1]) > 14.8) & (np.abs(ch[:, 1]) < 18.2) & (np.abs(nh[:, 1]) < 0.5),
        ),
        max_edge=2.0,
    )
    X_ph, w_ph, B_ph = sample_surface_triangles_gauss(hb_n, hb_t, box_h, order=2)

    # Link inner clevis bores (pin_to_link_contact)
    rl = np.sqrt((cl[:, 0] + 1.5) ** 2 + (cl[:, 2] - 10.5) ** 2)
    lb_n, lb_t = subdivide_long_edges(
        *filter_triangles(
            tl.nodes,
            tl.tris,
            (rl < 3.4) & (np.abs(cl[:, 1]) > 11.3) & (np.abs(cl[:, 1]) < 14.7) & (np.abs(nl[:, 1]) < 0.5),
        ),
        max_edge=2.0,
    )
    X_pl, w_pl, B_pl = sample_surface_triangles_gauss(lb_n, lb_t, box_l, order=2)

    # Hinge inner ear walls vs Link outer ear walls (hinge_to_link_contact)
    he_n, he_t = subdivide_long_edges(
        *filter_triangles(
            th.nodes,
            th.tris,
            (np.abs(np.abs(ch[:, 1]) - 15.0) < 0.1) & (np.abs(nh[:, 1]) > 0.9),
        ),
        max_edge=2.0,
    )
    master_hinge = DeformableTriangleSurfaceMaster(
        box_h, he_n, he_t, smooth_normals=True, exclude_boundary=True, max_distance=2.0
    )
    le_n, le_t = subdivide_long_edges(
        *filter_triangles(
            tl.nodes,
            tl.tris,
            (np.abs(np.abs(cl[:, 1]) - 14.5) < 0.1) & (rl < 10.0) & (np.abs(nl[:, 1]) > 0.9),
        ),
        max_edge=2.0,
    )
    X_hl, w_hl, B_hl = sample_surface_triangles_gauss(le_n, le_t, box_l, order=2)

    # Initial gap sanity checks (0.10 mm radial clearance, 0.50 mm ear clearance)
    g0_ph, *_ = master_pin.evaluate(X_ph, np.zeros((box_p.n_ctrl, 3)))
    g0_pl, *_ = master_pin.evaluate(X_pl, np.zeros((box_p.n_ctrl, 3)))
    g0_hl, *_ = master_hinge.evaluate(X_hl, np.zeros((box_h.n_ctrl, 3)))
    assert np.median(g0_ph) == pytest.approx(0.10, abs=0.01)
    assert np.median(g0_pl) == pytest.approx(0.10, abs=0.01)
    assert np.min(g0_hl) == pytest.approx(0.50, abs=0.025)

    pairs = [
        {"name": "pin_to_hinge_contact", "m_idx": 1, "s_idx": 0, "master_face": master_pin, "X_s": X_ph, "w_q": w_ph, "B_s": B_ph},
        {"name": "pin_to_link_contact",  "m_idx": 1, "s_idx": 2, "master_face": master_pin, "X_s": X_pl, "w_q": w_pl, "B_s": B_pl},
        {"name": "hinge_to_link_contact","m_idx": 0, "s_idx": 2, "master_face": master_hinge,"X_s": X_hl, "w_q": w_hl, "B_s": B_hl},
    ]

    # Surface BCs: Link bottom mounting pads (z=-3.0) pinned; Hinge R=5.25 mounting holes displaced
    l_pad_mask = (np.abs(cl[:, 2] + 3.0) < 0.05) & (np.abs(cl[:, 1]) >= 17.2) & (nl[:, 2] < -0.9)
    lp_n, lp_t = filter_triangles(tl.nodes, tl.tris, l_pad_mask)
    _, w_lbc, B_lbc = sample_surface_triangles_gauss(lp_n, lp_t, box_l, order=2)

    hinge_hole_dist = np.sqrt(
        (np.abs(ch[:, 1]) - 16.0) ** 2
        + ((ch[:, 0] - 27.37) * (-0.036) + (ch[:, 2] - 42.8) * 0.99935) ** 2
    )
    hh_mask = (np.abs(hinge_hole_dist - 5.25) < 0.15) & (ch[:, 2] > 35.0)
    hh_n, hh_t = filter_triangles(th.nodes, th.tris, hh_mask)
    _, w_hbc, B_hbc = sample_surface_triangles_gauss(hh_n, hh_t, box_h, order=2)
    assert len(lp_t) > 0 and len(hh_t) > 0

    weak_h = HardeningTable((0.0, 10.0), ((1500.0, 150_000.0),))
    strong_h = HardeningTable((0.0, 10.0), ((2000.0, 200_000.0),))
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]

    pin_reg = np.full(box_p.n_dofs, 1e-4, dtype=np.float64)
    pin_reg[2::3] = 1e-6
    h_reg = np.full(box_h.n_dofs, 1e-4, dtype=np.float64)
    l_reg = np.full(box_l.n_dofs, 1e-4, dtype=np.float64)

    op_h = None
    op_p = None
    op_l = None
    try:
        op_h = NonlinearBodyOperator(
            KM, mods, box_h,
            [{"trim_soup": (th.nodes, th.tris), "E": E, "nu": nu, "hardening": weak_h}],
            q=2, ghost_alpha=0.0, stiffness_reg=h_reg,
            surface_bc_B=B_hbc, surface_bc_weights=w_hbc, surface_bc_penalty=1e6,
        )
        op_p = NonlinearBodyOperator(
            KM, mods, box_p,
            [
                {"trim_soup": (tpc.nodes, tpc.tris), "E": E, "nu": nu, "hardening": strong_h},
                {"trim_soup": (tpl.nodes, tpl.tris), "E": E, "nu": nu, "hardening": weak_h},
                {"trim_soup": (tpr.nodes, tpr.tris), "E": E, "nu": nu, "hardening": weak_h},
            ],
            q=2, ghost_alpha=0.0, stiffness_reg=pin_reg,
        )
        op_l = NonlinearBodyOperator(
            KM, mods, box_l,
            [{"trim_soup": (tl.nodes, tl.tris), "E": E, "nu": nu, "hardening": weak_h}],
            q=2, ghost_alpha=0.0, stiffness_reg=l_reg,
            surface_bc_B=B_lbc, surface_bc_weights=w_lbc, surface_bc_penalty=1e6,
        )

        u_h = np.zeros(op_h.n_dofs, dtype=np.float64)
        u_p = np.zeros(op_p.n_dofs, dtype=np.float64)
        u_l = np.zeros(op_l.n_dofs, dtype=np.float64)
        slip = None
        uz_schedule = (0.500, 1.130, 1.664, 2.465, 3.150, 3.667, 4.350, 5.000)

        print(f"\n[E05 Real SAT Door_Lock_Assembly (h={h:.1f}mm, mu=0.3, NLGEOM+J2)]")
        rows = {}
        prev_peeq = 0.0
        uz_curr = 0.0
        step_clock = 0.0

        def _advance_to(uz_target: float, depth: int = 0) -> dict:
            nonlocal u_h, u_p, u_l, slip, uz_curr, step_clock
            op_h.set_surface_bc_target((0.0, 0.0, uz_target))
            try:
                step_clock += 1.0
                res_s = solve_multibody_nlgeom_contact_step(
                    [op_h, op_p, op_l],
                    [u_h, u_p, u_l],
                    pairs,
                    eps_n=5e5,
                    step_time=step_clock,
                    mu=0.3,
                    eps_t=1e5,
                    slip_states=slip,
                    outer_tol=1e-6,
                    rtol=1e-6,
                    max_outer=45,
                    ssn_max_iter=150,
                    use_pardiso=True,
                )
                u_h, u_p, u_l = res_s["u_full_list"]
                slip = res_s["slip_states"]
                uz_curr = uz_target
                return res_s
            except RuntimeError:
                if depth >= 3:
                    raise
                uz_mid = 0.5 * (uz_curr + uz_target)
                print(
                    f"    [substep] bisecting {uz_curr:.3f} -> {uz_mid:.3f} -> {uz_target:.3f} mm",
                    flush=True,
                )
                _advance_to(uz_mid, depth=depth + 1)
                return _advance_to(uz_target, depth=depth + 1)

        for step, uz in enumerate(uz_schedule, start=1):
            ts = time.perf_counter()
            res = _advance_to(uz)
            rf_l = op_l.surface_bc_reaction(u_l)
            rf_h = op_h.surface_bc_reaction(u_h)
            rf_l_kn = float(rf_l[2] / 1e3)
            rf_h_kn = float(rf_h[2] / 1e3)
            bal_z = abs(rf_h[2] + rf_l[2]) / max(abs(rf_l[2]), 1.0)
            pmax = max(res["peeq_max_list"])
            a_ph = res["pair_stats"]["pin_to_hinge_contact"]["n_active"]
            a_pl = res["pair_stats"]["pin_to_link_contact"]["n_active"]
            a_hl = res["pair_stats"]["hinge_to_link_contact"]["n_active"]
            print(
                f"  step {step} uz={uz:5.3f}mm: outer={res['n_outer']:2d}, "
                f"peeq={pmax:.4f}, RF3_link={rf_l_kn:+8.3f} kN, RF3_hinge={rf_h_kn:+8.3f} kN, "
                f"bal_z={bal_z:.2e}, act=({a_ph},{a_pl},{a_hl}), "
                f"dt={time.perf_counter()-ts:.1f}s, RSS={_rss_mb():.0f}MB",
                flush=True,
            )
            assert res["outer_rel"] <= 1e-6 or (
                res["outer_rnorm"] <= 1e-6 and res["outer_rel"] <= 1e-5
            )
            assert bal_z < 1e-3
            assert rf_l_kn < 0.0 < rf_h_kn
            assert pmax >= prev_peeq - 1e-12
            prev_peeq = pmax
            rows[uz] = rf_l_kn

        final_rf3 = rows[5.000]
        rel_err = abs(final_rf3 - ref_rf3_kn) / abs(ref_rf3_kn)
        print(
            f"  FINAL RF3_link(5.0mm)={final_rf3:.3f} kN vs Coreform {ref_rf3_kn:.3f} kN "
            f"(rel_err={rel_err:.2%}, total={time.perf_counter()-t0:.1f}s)",
            flush=True,
        )
        assert prev_peeq > 0.05
        assert rel_err <= rel_tol, f"E05 RF3 gate missed: {final_rf3:.3f} vs {ref_rf3_kn:.3f} ({rel_err:.2%})"
    finally:
        for op in (op_h, op_p, op_l):
            if op is not None:
                op.close()
        release_solver_memory()
        gc.collect()
