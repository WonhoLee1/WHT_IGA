# -*- coding: utf-8 -*-
"""P5 / S33: E03 13-solid STEP inspection & reduced 3-step thermoplastic model.

Inspects ``glass_to_metal_seal_connector.stp`` (13 watertight B-Rep solids:
1 S-8061 glass seal + 11 Alloy 52 pins + 1 304L VAR header), partitions a coarse
background ``BoxSpline`` across the three material regions with ``build_trimmed_table``
(``q=2``, ``cache_dir="off"``), verifies quadrature volume conservation with zero
fallbacks, checks zero-stress free thermal expansion under minimal 3-2-1 kinematic
restraint, and runs all 3 thermal steps (450 -> 20 -> -55 -> 200 C) via
``thermal_plasticity.run_thermal_steps``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses import nonlinear as NL
from wht_iga.analyses import thermal_plasticity as TP
from wht_iga.geometry.cad import brep_area, brep_volume, shape_bbox, tessellate, validate_shape
from wht_iga.geometry.occ_io import _count_free_edges, load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.active import fix_inactive, inactive_ctrl_mask
from wht_iga.solver.trimmed import build_trimmed_table

_E03_ROOT = (
    Path(__file__).resolve().parents[2]
    / "research-cf-IGA"
    / "user-guide"
    / "example-problems"
    / "glass-to-metal-seal"
)
_STP_PATH = _E03_ROOT / "downloads" / "glass_to_metal_seal_connector.stp"
_DATA_DIR = _E03_ROOT / "data"


def _solid_centroid(solid) -> np.ndarray:
    """Exact B-Rep center of mass via ``BRepGProp.VolumeProperties_s``."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props)
    c = props.CentreOfMass()
    return np.array([float(c.X()), float(c.Y()), float(c.Z())], dtype=float)


def _merge_tessellations(solids, deflection: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
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


def _classify_e03_solids(solids: list) -> dict[str, list[int]]:
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


@pytest.fixture(scope="module")
def e03_step_bundle():
    """Load E03 STEP once per module and build the coarse 3-region trimmed partition."""
    cad = load_cad_file(_STP_PATH, heal=False)
    groups = _classify_e03_solids(cad.solids)

    solid_info = []
    for idx, solid in enumerate(cad.solids):
        tess = tessellate(solid, deflection=0.05)
        rep = validate_shape(solid, tess, n_samples=10)
        lo_v, hi_v = tess.nodes.min(axis=0), tess.nodes.max(axis=0)
        solid_info.append({
            "index": idx,
            "volume": brep_volume(solid),
            "area": brep_area(solid),
            "centroid": _solid_centroid(solid),
            "bbox": shape_bbox(solid),
            "vertex_bbox": (lo_v, hi_v),
            "free_edges": _count_free_edges(solid),
            "watertight": rep.watertight,
            "tess": tess,
        })

    h_tess = solid_info[groups["header"][0]]["tess"]
    g_tess = solid_info[groups["glass"][0]]["tess"]
    p_nodes, p_tris = _merge_tessellations([cad.solids[i] for i in groups["pins"]], deflection=0.05)

    lo, hi = cad.bbox
    pad = 0.05
    box = BoxSpline(lo=tuple(lo - pad), hi=tuple(hi + pad), n_cells=(8, 6, 6), degree=1)
    q_order = 2

    region_tables = {}
    for name, nds, trs, idxs in [
        ("header", h_tess.nodes, h_tess.tris, groups["header"]),
        ("glass", g_tess.nodes, g_tess.tris, groups["glass"]),
        ("pins", p_nodes, p_tris, groups["pins"]),
    ]:
        qp_table, out_cells, issues, recs, diag = build_trimmed_table(
            nds, trs, box, q=q_order, workers=1, cache_dir="off"
        )
        region_tables[name] = {
            "qp_table": qp_table,
            "out_cells": out_cells,
            "issues": issues,
            "recs": recs,
            "diag": diag,
            "brep_volume": sum(solid_info[i]["volume"] for i in idxs),
        }

    return {
        "cad": cad,
        "groups": groups,
        "solid_info": solid_info,
        "box": box,
        "q_order": q_order,
        "regions": region_tables,
    }


def _build_and_run_e03(bundle, mats: list[TP.ThermoMaterial], increments: list[int]):
    """Build the 3-region bonded Kratos model on the coarse BoxSpline and run thermal steps."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]
    box: BoxSpline = bundle["box"]
    nx, ny, nz = box.n_cells

    region_order = ("header", "glass", "pins")
    reg_pairs = [
        (bundle["regions"][name]["qp_table"], bundle["regions"][name]["out_cells"])
        for name in region_order
    ]

    model = KM.Model()
    mp = model.CreateModelPart("E03_Reduced")
    for var in (KM.DISPLACEMENT, KM.REACTION):
        mp.AddNodalSolutionStepVariable(var)
    NL.enforce_static_settings(mp)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    volume, _nodes = make_nurbs_volume(KM, mp, box)
    elements, elem_mat_idx, active_union = TP.attach_region_elements(
        KM, mp, volume, box, reg_pairs, element_name="SmallDisplacementElement3D4N",
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
    fix_inactive(KM, mp, box, all_out)
    p012 = TP.fix_isostatic_321(KM, mp, box, inactive_mask=inact, axis=0, side=0)

    mp.CloneTimeStep(0.0)
    states = [450.0, 20.0, -55.0, 200.0]
    groups = TP.CellThermalGroups(
        KM, mods, mp, elements, mats, elem_mat_idx, states,
        law_factory=TP.default_law_factory(mods, finite_strain=False),
    )

    snapshots: list[dict] = []

    def observer(level: float) -> None:
        if abs(level - round(level)) > 1e-9:
            return
        entry: dict[str, object] = {"level": float(round(level))}
        for m_idx, name in enumerate(region_order):
            m_els = [el for el, mi in zip(elements, elem_mat_idx) if mi == m_idx]
            sigs = np.vstack([
                np.array([list(v) for v in el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)])
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
            entry[name] = {
                "vm_max": float(vm.max()),
                "s1_max": float(s1.max()),
                "s3_min": float(s3.min()),
                "abs_max": float(np.abs(sigs).max()),
                "diss_max": diss_max,
            }
        snapshots.append(entry)

    out = TP.run_thermal_steps(
        mp,
        KM,
        lambda: NL.make_newton_strategy(KM, mp, max_iters=30),
        groups,
        increments=increments,
        on_converged=observer,
    )
    return {
        "out": out,
        "groups": groups,
        "snapshots": snapshots,
        "n_elements": len(elements),
        "n_active_cells": len(active_union),
        "n_inactive_ctrl": int(inact.sum()),
        "p012": p012,
        "CL": CL,
    }


# --------------------------------------------------------------------------- tests

def test_e03_step_13_solids_inspection(e03_step_bundle):
    """Inspect glass_to_metal_seal_connector.stp: 13 watertight solids matching e03.json."""
    cad = e03_step_bundle["cad"]
    groups = e03_step_bundle["groups"]
    info = e03_step_bundle["solid_info"]

    assert len(cad.solids) == 13
    assert cad.healing_report.watertight is True
    assert cad.healing_report.free_edges_after == 0
    assert all(s["free_edges"] == 0 and s["watertight"] for s in info)

    # Zero overlap between the 13 solids: sum of solid volumes == compound B-Rep volume.
    sum_vol = sum(s["volume"] for s in info)
    assert sum_vol == pytest.approx(cad.volume, rel=1e-4)
    assert cad.volume == pytest.approx(2590.49, rel=1e-3)

    # 1 header + 1 glass + 11 pins
    h = info[groups["header"][0]]
    g = info[groups["glass"][0]]
    pins = [info[i] for i in groups["pins"]]

    assert h["volume"] == pytest.approx(2226.00, rel=1e-3)
    assert g["volume"] == pytest.approx(156.46, rel=1e-3)
    assert sum(p["volume"] for p in pins) == pytest.approx(208.02, rel=1e-3)

    # MMTS mm dimensions: header mounting face at x = -25.50 mm, glass x in [-13.05, -9.95] mm
    assert h["vertex_bbox"][0][0] == pytest.approx(-25.50, abs=0.05)
    assert h["vertex_bbox"][1][0] == pytest.approx(0.00, abs=0.05)
    assert g["vertex_bbox"][0][0] == pytest.approx(-13.05, abs=0.05)
    assert g["vertex_bbox"][1][0] == pytest.approx(-9.95, abs=0.05)

    # 11 pins split cleanly into 8 outer-ring (x <= 0.45) + 3 inner-ring (x <= 2.45), z=0 symmetric
    outer_pins = [p for p in pins if p["vertex_bbox"][1][0] < 1.0]
    inner_pins = [p for p in pins if p["vertex_bbox"][1][0] > 2.0]
    assert len(outer_pins) == 8
    assert len(inner_pins) == 3

    outer_r = [float(np.hypot(p["centroid"][1], p["centroid"][2])) for p in outer_pins]
    assert np.allclose(outer_r, outer_r[0], rtol=1e-2)
    assert outer_r[0] == pytest.approx(3.241, abs=0.02)

    # Mirror symmetry across z = 0 plane: sum of pin z-centroids == 0
    assert abs(sum(p["centroid"][2] for p in pins)) < 1e-3


def test_e03_multibody_trimmed_quadrature_partition(e03_step_bundle):
    """Coarse 3-region trimmed quadrature recovers B-Rep solid volumes with 0 fallbacks."""
    box: BoxSpline = e03_step_bundle["box"]
    nx, ny, nz = box.n_cells
    n_gauss = e03_step_bundle["q_order"] + 1

    active_union: set[tuple[int, int, int]] = set()
    total_qp_vol = 0.0
    total_brep_vol = 0.0

    for name, reg in e03_step_bundle["regions"].items():
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
                    if cell in qp_table:
                        pts = qp_table[cell]
                    else:
                        pts = np.asarray(cell_gauss_points(box, i, j, k, n_gauss=n_gauss), dtype=float)
                    assert np.all(pts[:, 3] > 0.0)
                    assert np.all(pts[:, :3] >= np.asarray(box.lo) - 1e-12)
                    assert np.all(pts[:, :3] <= np.asarray(box.hi) + 1e-12)
                    vol_qp += float(pts[:, 3].sum())

        bvol = reg["brep_volume"]
        rel_err = abs(vol_qp - bvol) / bvol
        assert rel_err < 0.01, f"{name} volume mismatch: {vol_qp:.4f} vs {bvol:.4f}"
        total_qp_vol += vol_qp
        total_brep_vol += bvol

    assert 100 <= len(active_union) <= 600
    assert abs(total_qp_vol - total_brep_vol) / total_brep_vol < 0.01


def test_e03_reduced_3step_thermal_cycle(e03_step_bundle):
    """Run 450 -> 20 -> -55 -> 200 C: zero free-expansion artifact + physical J2/glass response."""
    mats_dict = TP.e03_materials_initial_yield(_DATA_DIR)

    # 1. Free-thermal-expansion artifact check: uniform material across all 3 regions
    #    with minimal 3-2-1 restraint must remain stress-free (< 1e-8 MPa) across all 3 steps.
    free_res = _build_and_run_e03(
        e03_step_bundle,
        [mats_dict["SS304L_header"], mats_dict["SS304L_header"], mats_dict["SS304L_header"]],
        increments=[1, 1, 1],
    )
    assert free_res["out"]["steps"] == 3
    assert free_res["out"]["n_solves"] == 3
    for snap in free_res["snapshots"]:
        for name in ("header", "glass", "pins"):
            assert snap[name]["abs_max"] < 1e-8
            assert snap[name]["diss_max"] == pytest.approx(0.0, abs=1e-25)

    # 2. Full 3-material E03 thermal cycle (450 -> 20 -> -55 -> 200 C)
    res = _build_and_run_e03(
        e03_step_bundle,
        [mats_dict["SS304L_header"], mats_dict["S8061_glass"], mats_dict["Alloy52_pins"]],
        increments=[3, 3, 3],
    )
    assert res["out"]["steps"] == 9
    assert res["out"]["n_solves"] == 9, "unexpected Newton cutback in coarse E03 cycle"
    assert res["out"]["level"] == pytest.approx(3.0)
    assert res["groups"].n_groups == 3

    s_room, s_cold, s_hot = res["snapshots"]
    assert [s["level"] for s in res["snapshots"]] == [1.0, 2.0, 3.0]

    # Metal yielding at initial yield stress (304L: 162 -> 194.4 -> 120.2 MPa; Alloy 52: 234.4 MPa)
    assert s_room["header"]["vm_max"] == pytest.approx(162.0, rel=1e-3)
    assert s_cold["header"]["vm_max"] == pytest.approx(194.4, rel=1e-3)
    assert s_hot["header"]["vm_max"] == pytest.approx(120.2, rel=1e-3)
    for snap in (s_room, s_cold, s_hot):
        assert snap["pins"]["vm_max"] == pytest.approx(234.4, rel=1e-3)

    # Monotonically increasing plastic dissipation in both header and pins
    assert s_room["header"]["diss_max"] > 0.0
    assert s_cold["header"]["diss_max"] > s_room["header"]["diss_max"]
    assert s_hot["header"]["diss_max"] >= s_cold["header"]["diss_max"] - 1e-30
    assert s_room["pins"]["diss_max"] > 0.0
    assert s_cold["pins"]["diss_max"] > s_room["pins"]["diss_max"]
    assert s_hot["pins"]["diss_max"] > s_cold["pins"]["diss_max"]

    # Glass thermoelastic response: compressive clamping during cooling + residual tension on reheat
    assert s_room["glass"]["s3_min"] < -300.0
    assert s_cold["glass"]["s3_min"] < s_room["glass"]["s3_min"]
    assert s_room["glass"]["s1_max"] > 0.0
    assert s_cold["glass"]["s1_max"] > s_room["glass"]["s1_max"]
    assert s_hot["glass"]["s1_max"] > s_cold["glass"]["s1_max"]
