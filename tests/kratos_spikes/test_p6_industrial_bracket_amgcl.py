# -*- coding: utf-8 -*-
"""Phase 6 S45: Industrial 3D CAD Trimmed IGA + Extension + AMGCL Scaling.

Exercises real CAD-to-trimmed-IGA-to-AMGCL end-to-end on medium industrial
STEP models in ``wht_iga/resources/`` (``p=2``, ``q=2``, ``theta=0.2`` extension):

1. Inspects ``Generator_bracket.step``, ``Flange.STEP``, and ``Support_plate.STEP``
   (B-Rep healing/solid counts, watertightness, exact B-Rep volume/area, and
   tessellation quality).
2. Compares axis-aligned vs PCA OBB (``compute_pca_obb``) on
   ``Generator_bracket.step`` (active-cell ratio, ``n_fallback == 0``,
   quadrature volume error ``< 2%``, ``n_master / n_extended / n_inactive``),
   plus trimmed ``p=2`` + extension Kratos solves comparing ``PardisoLU``,
   ``AMGCL_CG_ILU0``, and ``AMGCL_GMRES_ILU0``.
3. Benchmarks ``Flange.STEP`` (solid 0) trimmed ``p=2`` + ``theta=0.2`` extension
   linear elastic Kratos solves at ~5.3k and ~8.5k active DOFs comparing
   ``PardisoLU``, ``AMGCL_CG_ILU0``, and ``AMGCL_GMRES_ILU0``
   (``rel_l2_diff < 1e-7``, stage timings, Peak RSS ``< 1.5 GB``, and global
   reaction equilibrium ``|sum(R)| < 1e-6``).
"""

from __future__ import annotations

import gc
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import psutil
import pytest

from wht_iga.analyses.methods import ObbFrame, compute_pca_obb
from wht_iga.geometry import cad
from wht_iga.geometry.occ_io import load_cad_file
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.solver.active import fix_inactive, inactive_ctrl_mask
from wht_iga.stabilize import extension as EXT
from wht_iga.stabilize import support as SUP

_RESOURCES = Path(__file__).resolve().parents[2] / "wht_iga" / "resources"
_PROC = psutil.Process()

pytestmark = [
    pytest.mark.kratos,
    pytest.mark.skipif(not _RESOURCES.is_dir(), reason=f"Resources dir missing: {_RESOURCES}"),
]

RSS_LIMIT_MB = 1500.0
DEGREE = 2
QUAD_DEG = 2
THETA = 0.2
E_MOD = 210_000.0
NU_POISSON = 0.3
SOLVER_KEYS = ("PardisoLU", "AMGCL_CG_ILU0", "AMGCL_GMRES_ILU0")


def _rss_mb() -> float:
    """Return current process resident set size (RSS) in decimal MB (10^6 bytes)."""
    return float(_PROC.memory_info().rss) / 1e6


def _cleanup(*objs: Any) -> float:
    """Drop references, run ``gc.collect()``, and return post-GC RSS in MB."""
    del objs
    gc.collect()
    return _rss_mb()


def _make_amgcl_or_pardiso(KM: Any, LA: Any, solver_key: str, tol: float = 1e-10, max_iter: int = 5000) -> Any:
    """Instantiate a Kratos ``PardisoLUSolver`` or ``AMGCLSolver`` by key.

    Parameters
    ----------
    KM:
        Imported ``KratosMultiphysics`` module.
    LA:
        Imported ``KratosMultiphysics.LinearSolversApplication`` module.
    solver_key:
        One of ``"PardisoLU"``, ``"AMGCL_CG_ILU0"``, or ``"AMGCL_GMRES_ILU0"``.
    tol:
        Relative residual tolerance for iterative AMGCL solvers.
    max_iter:
        Maximum iteration cap for iterative AMGCL solvers.
    """
    if solver_key == "PardisoLU":
        return LA.PardisoLUSolver()
    if solver_key == "AMGCL_CG_ILU0":
        params = {
            "solver_type": "amgcl",
            "krylov_type": "cg",
            "smoother_type": "ilu0",
            "coarsening_type": "aggregation",
            "max_iteration": int(max_iter),
            "tolerance": float(tol),
            "verbosity": 0,
            "use_block_matrices_if_possible": False,
        }
        return KM.AMGCLSolver(KM.Parameters(json.dumps(params)))
    if solver_key == "AMGCL_GMRES_ILU0":
        params = {
            "solver_type": "amgcl",
            "krylov_type": "gmres",
            "smoother_type": "ilu0",
            "coarsening_type": "aggregation",
            "max_iteration": int(max_iter),
            "tolerance": float(tol),
            "verbosity": 0,
            "gmres_krylov_space_dimension": 100,
            "use_block_matrices_if_possible": False,
        }
        return KM.AMGCLSolver(KM.Parameters(json.dumps(params)))
    raise ValueError(f"Unsupported solver_key {solver_key!r}")


def _build_trimmed_and_classify(
    nodes: np.ndarray,
    tris: np.ndarray,
    box: BoxSpline,
    brep_vol: float,
    *,
    q: int = QUAD_DEG,
    theta: float = THETA,
    workers: int = 8,
) -> dict[str, Any]:
    """Run ``build_trimmed_table``, basis classification, and extension build.

    Parameters
    ----------
    nodes:
        Tessellation vertex coordinates ``(N, 3)``.
    tris:
        Tessellation triangle connectivity ``(M, 3)``.
    box:
        Background ``BoxSpline`` grid.
    brep_vol:
        Reference exact B-Rep solid volume for quadrature error check.
    q:
        Moment quadrature degree.
    theta:
        Basis-support small-cell threshold for extension classification.
    workers:
        Thread-pool worker count passed to ``build_trimmed_table``.
    """
    t0 = time.perf_counter()
    qp_table, out_cells, issues, cell_records, diag = T.build_trimmed_table(
        nodes, tris, box, q=int(q), workers=int(workers), cache_dir="off"
    )
    t_trim = time.perf_counter() - t0

    nx, ny, nz = (int(v) for v in box.n_cells)
    n_total_cells = nx * ny * nz
    n_active_cells = n_total_cells - len(out_cells)
    h = (np.asarray(box.hi, float) - np.asarray(box.lo, float)) / np.asarray([nx, ny, nz])
    cell_vol = float(np.prod(h))

    eta = np.zeros(n_total_cells, dtype=np.float64)
    for rec in cell_records:
        eta[int(rec["i"]) + nx * (int(rec["j"]) + ny * int(rec["k"]))] = float(rec.get("eta", 1.0))

    quad_vol = float(eta.sum() * cell_vol)
    vol_err = abs(quad_vol - brep_vol) / brep_vol
    knots = [np.asarray(k, dtype=float) for k in box.knots]

    t1 = time.perf_counter()
    cls = SUP.classify_basis(
        tuple(int(v) for v in box.n_cp),
        knots,
        eta * cell_vol,
        np.full(n_total_cells, cell_vol),
        theta=float(theta),
    )
    t_cls = time.perf_counter() - t1

    t2 = time.perf_counter()
    ext_op = EXT.build_extension(
        tuple(int(v) for v in box.n_cp),
        cls.labels,
        int(box.degree),
        knots=knots,
    )
    t_ext = time.perf_counter() - t2

    n_master = int(cls.n_master)
    n_extended = int(cls.n_extended)
    n_inactive = int(cls.n_inactive)
    active_dofs = 3 * (n_master + n_extended)

    return {
        "box": box,
        "qp_table": qp_table,
        "out_cells": out_cells,
        "issues": issues,
        "cell_records": cell_records,
        "diag": diag,
        "n_total_cells": n_total_cells,
        "n_active_cells": n_active_cells,
        "active_cell_ratio": float(n_active_cells / n_total_cells),
        "n_cut": int(diag["n_cut"]),
        "n_fallback": int(diag["n_fallback"]),
        "n_relaxed": int(diag["n_relaxed"]),
        "n_qp_cut": int(diag["n_qp_cut"]),
        "brep_vol": float(brep_vol),
        "quad_vol": quad_vol,
        "vol_err": float(vol_err),
        "n_master": n_master,
        "n_extended": n_extended,
        "n_inactive": n_inactive,
        "active_dofs": active_dofs,
        "n_thin_feature": int(ext_op.stats.get("n_thin_feature", 0)),
        "build_trimmed_table_s": float(t_trim),
        "classify_basis_s": float(t_cls),
        "build_extension_s": float(t_ext),
    }


def _populate_kratos_model(
    KM: Any,
    SM: Any,
    model_name: str,
    box: BoxSpline,
    qp_table: dict[tuple[int, int, int], np.ndarray],
    out_cells: list[tuple[int, int, int]],
    cell_records: list[dict[str, Any]],
    *,
    axis_idx: int,
    disp_prescribed: float = 0.05,
    q: int = QUAD_DEG,
    theta: float = THETA,
) -> dict[str, Any]:
    """Build a fresh Kratos ``ModelPart`` with trimmed elements, BCs, and extension ties."""
    model = KM.Model()
    mp = model.CreateModelPart(model_name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    t0 = time.perf_counter()
    volume, nodes_vec = make_nurbs_volume(KM, mp, box)
    t_vol = time.perf_counter() - t0

    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E_MOD))
    props.SetValue(KM.POISSON_RATIO, float(NU_POISSON))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    nx, ny, nz = (int(v) for v in box.n_cells)
    out_set = {tuple(c) for c in out_cells}
    n_gauss = int(q) + 1
    eid = 1
    elements: list[Any] = []
    t1 = time.perf_counter()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = (i, j, k)
                if cell in out_set:
                    continue
                pts = (
                    qp_table[cell].tolist()
                    if cell in qp_table
                    else cell_gauss_points(box, i, j, k, n_gauss=n_gauss)
                )
                qps = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qps, 2, pts)
                elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))
                eid += 1
    t_elem = time.perf_counter() - t1

    fix_inactive(KM, mp, box, out_cells)
    inact = inactive_ctrl_mask(box, out_cells)
    cp = box.control_points()
    act_ids = [i for i in range(box.n_ctrl) if not inact[i]]
    coords = cp[act_ids, axis_idx]
    c_min = float(coords.min())
    c_max = float(coords.max())
    tol_band = 0.05 * max(c_max - c_min, 1e-12)

    d_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    load_var = d_vars[axis_idx]
    n_fixed = 0
    n_loaded = 0
    for idx in act_ids:
        node = mp.GetNode(idx + 1)
        val = float(cp[idx, axis_idx])
        if val <= c_min + tol_band:
            for comp in d_vars:
                node.Fix(comp)
                node.SetSolutionStepValue(comp, 0.0)
            n_fixed += 1
        elif val >= c_max - tol_band:
            node.Fix(load_var)
            node.SetSolutionStepValue(load_var, float(disp_prescribed))
            n_loaded += 1
    assert n_fixed > 0 and n_loaded > 0

    t2 = time.perf_counter()
    tie_info = T.apply_extension_ties(KM, mp, box, cell_records, theta=float(theta))
    t_ties = time.perf_counter() - t2

    mp.SetBufferSize(2)
    mp.CloneTimeStep(1.0)
    return {
        "model": model,
        "mp": mp,
        "volume": volume,
        "nodes_vec": nodes_vec,
        "props": props,
        "elements": elements,
        "tie_info": tie_info,
        "make_vol_s": float(t_vol),
        "attach_elem_s": float(t_elem),
        "apply_ties_s": float(t_ties),
        "n_elements": len(elements),
        "n_ties": int(tie_info["n_ties"]),
        "n_fixed_nodes": n_fixed,
        "n_loaded_nodes": n_loaded,
    }


def _benchmark_kratos_solvers(
    label: str,
    prep: dict[str, Any],
    *,
    axis_idx: int,
    obb_frame: ObbFrame | None = None,
) -> dict[str, Any]:
    """Compare ``PardisoLU``, ``AMGCL_CG_ILU0``, and ``AMGCL_GMRES_ILU0`` on fresh ModelParts."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    LA = mods["LinearSolversApplication"]
    import KratosMultiphysics.scipy_conversion_tools as sct

    box: BoxSpline = prep["box"]
    qp_table = prep["qp_table"]
    out_cells = prep["out_cells"]
    cell_records = prep["cell_records"]

    peak_rss = _rss_mb()
    u_ref: np.ndarray | None = None
    assembly_s = 0.0
    n_eq = 0
    nnz = 0
    solver_rows: list[dict[str, Any]] = []
    stage_timings: dict[str, float] = {}

    for s_idx, solver_key in enumerate(SOLVER_KEYS):
        gc.collect()
        rss_before = _rss_mb()
        built = _populate_kratos_model(
            KM,
            SM,
            f"{label}_{solver_key}",
            box,
            qp_table,
            out_cells,
            cell_records,
            axis_idx=axis_idx,
        )
        mp = built["mp"]
        if s_idx == 0:
            stage_timings = {
                "make_vol_s": built["make_vol_s"],
                "attach_elem_s": built["attach_elem_s"],
                "apply_ties_s": built["apply_ties_s"],
            }
            scheme_g = KM.ResidualBasedIncrementalUpdateStaticScheme()
            solver_g = _make_amgcl_or_pardiso(KM, LA, solver_key)
            builder_g = KM.ResidualBasedBlockBuilderAndSolver(solver_g)
            strat_g = KM.ResidualBasedLinearStrategy(mp, scheme_g, builder_g, False, False, False, False)
            strat_g.SetEchoLevel(0)
            strat_g.Initialize()
            A_tmp = KM.CompressedMatrix()
            Dx_tmp = KM.Vector()
            b_tmp = KM.Vector()
            builder_g.SetUpDofSet(scheme_g, mp)
            builder_g.SetUpSystem(mp)
            builder_g.ResizeAndInitializeVectors(scheme_g, A_tmp, Dx_tmp, b_tmp, mp)
            scheme_g.InitializeSolutionStep(mp, A_tmp, Dx_tmp, b_tmp)
            builder_g.InitializeSolutionStep(mp, A_tmp, Dx_tmp, b_tmp)
            scheme_g.InitializeNonLinIteration(mp, A_tmp, Dx_tmp, b_tmp)
            tb0 = time.perf_counter()
            builder_g.Build(scheme_g, mp, A_tmp, b_tmp)
            if mp.NumberOfMasterSlaveConstraints() > 0:
                builder_g.ApplyConstraints(scheme_g, mp, A_tmp, b_tmp)
            builder_g.ApplyDirichletConditions(scheme_g, mp, A_tmp, Dx_tmp, b_tmp)
            assembly_s = time.perf_counter() - tb0
            k_csr = sct.to_csr(A_tmp)
            n_eq = int(k_csr.shape[0])
            nnz = int(k_csr.nnz)
            del k_csr, A_tmp, Dx_tmp, b_tmp, strat_g, builder_g, solver_g, scheme_g
            gc.collect()

        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        solver = _make_amgcl_or_pardiso(KM, LA, solver_key)
        builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
        strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, True, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()

        t0 = time.perf_counter()
        strategy.Solve()
        strategy_s = time.perf_counter() - t0
        rss_after = _rss_mb()
        peak_rss = max(peak_rss, rss_after)

        d_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
        r_vars = (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z)
        U_loc = np.array(
            [[float(mp.GetNode(i + 1).GetSolutionStepValue(c)) for c in d_vars] for i in range(box.n_ctrl)],
            dtype=np.float64,
        )
        R_loc = np.array(
            [[float(mp.GetNode(i + 1).GetSolutionStepValue(c)) for c in r_vars] for i in range(box.n_ctrl)],
            dtype=np.float64,
        )
        if obb_frame is not None:
            U = obb_frame.to_global_vectors(U_loc)
            R = obb_frame.to_global_vectors(R_loc)
        else:
            U = U_loc
            R = R_loc

        assert np.all(np.isfinite(U))
        assert np.all(np.isfinite(R))
        rxn_sum_norm = float(np.linalg.norm(R.sum(axis=0)))
        if u_ref is None:
            u_ref = U.copy()
        rel_diff = float(np.linalg.norm(U - u_ref) / max(float(np.linalg.norm(u_ref)), 1e-30))
        max_u = float(np.max(np.abs(U)))

        solver_rows.append(
            {
                "solver": solver_key,
                "strategy_s": float(strategy_s),
                "assembly_s": float(assembly_s),
                "solve_net_s": float(max(0.0, strategy_s - assembly_s)),
                "rel_l2_diff": rel_diff,
                "rxn_sum_norm": rxn_sum_norm,
                "max_u": max_u,
                "rss_before_mb": float(rss_before),
                "rss_after_mb": float(rss_after),
                "rss_delta_mb": float(rss_after - rss_before),
            }
        )

        del U, R, U_loc, R_loc, strategy, builder, solver, scheme, built, mp
        _cleanup()

    rss_end = _cleanup(u_ref)
    return {
        "label": label,
        "n_eq": n_eq,
        "nnz": nnz,
        "assembly_s": float(assembly_s),
        "stage_timings": stage_timings,
        "solvers": solver_rows,
        "peak_rss_mb": float(peak_rss),
        "rss_end_mb": float(rss_end),
    }


def test_industrial_cad_inspection() -> None:
    """Inspect Generator_bracket, Flange, and Support_plate STEP models."""
    t0 = time.perf_counter()
    expected = {
        "Generator_bracket.step": {"min_solids": 1, "min_faces": 200, "deflection_frac": 0.003},
        "Flange.STEP": {"min_solids": 12, "min_faces": 150, "deflection_frac": 0.003},
        "Support_plate.STEP": {"min_solids": 1, "min_faces": 70, "deflection_frac": 0.001},
    }
    stats: list[dict[str, Any]] = []
    for fname, cfg in expected.items():
        geom = load_cad_file(_RESOURCES / fname, heal=True)
        hr = geom.healing_report
        assert hr.solids_after >= cfg["min_solids"]
        assert hr.faces_after >= cfg["min_faces"]
        assert hr.watertight is True
        assert hr.free_edges_after == 0

        lo, hi = cad.shape_bbox(geom.shape)
        span = hi - lo
        diag = float(np.linalg.norm(span))
        bvol = cad.brep_volume(geom.shape)
        barea = cad.brep_area(geom.shape)
        assert bvol > 0.0 and barea > 0.0

        tess = cad.tessellate(geom.shape, deflection=cfg["deflection_frac"] * diag)
        svol = cad.signed_volume(tess)
        be = cad.boundary_edge_count(tess)
        assert be == 0
        assert svol > 0.0
        tess_vol_err = abs(svol - bvol) / bvol
        assert tess_vol_err < 0.01

        stats.append(
            {
                "file": fname,
                "solids": hr.solids_after,
                "faces": hr.faces_after,
                "span": [round(float(x), 3) for x in span],
                "brep_vol": bvol,
                "brep_area": barea,
                "tess_nodes": len(tess.nodes),
                "tess_tris": len(tess.tris),
                "tess_vol_err": tess_vol_err,
            }
        )

    dt = time.perf_counter() - t0
    assert _rss_mb() < RSS_LIMIT_MB
    for s in stats:
        print(
            f"\n[CAD] {s['file']}: solids={s['solids']}, faces={s['faces']}, "
            f"span={s['span']}, brep_vol={s['brep_vol']:.2f}, "
            f"tris={s['tess_tris']}, tess_vol_err={s['tess_vol_err']*100:.4f}%"
        )
    print(f"[CAD inspection] wall={dt:.2f}s, rss={_rss_mb():.1f}MB")


def test_generator_bracket_aa_vs_pca_obb() -> None:
    """Compare AA vs PCA OBB on Generator_bracket.step and run p=2 + extension AMGCL solves."""
    t0 = time.perf_counter()
    geom = load_cad_file(_RESOURCES / "Generator_bracket.step", heal=False)
    assert len(geom.solids) == 1
    bvol = cad.brep_volume(geom.shape)
    lo_aa, hi_aa = cad.shape_bbox(geom.shape)
    span_aa = hi_aa - lo_aa
    diag = float(np.linalg.norm(span_aa))
    t_tess0 = time.perf_counter()
    tess = cad.tessellate(geom.shape, deflection=0.003 * diag)
    t_tess = time.perf_counter() - t_tess0
    assert cad.boundary_edge_count(tess) == 0

    # Axis-aligned background BoxSpline (p=2, q=2, theta=0.2, ~3.1k active DOFs)
    pad_aa = 0.03 * span_aa
    box_aa = BoxSpline(
        lo=tuple((lo_aa - pad_aa).tolist()),
        hi=tuple((hi_aa + pad_aa).tolist()),
        n_cells=(8, 12, 14),
        degree=DEGREE,
    )
    prep_aa = _build_trimmed_and_classify(tess.nodes, tess.tris, box_aa, bvol)
    assert prep_aa["n_fallback"] == 0
    assert prep_aa["vol_err"] < 0.02
    assert 3000 <= prep_aa["active_dofs"] <= 15000
    assert (
        prep_aa["n_master"] + prep_aa["n_extended"] + prep_aa["n_inactive"]
        == box_aa.n_ctrl
    )

    # PCA OBB background BoxSpline (compute_pca_obb, p=2, q=2, theta=0.2)
    frame, lo_obb, hi_obb = compute_pca_obb(tess.nodes, padding_frac=0.05)
    pts_local = frame.to_local_points(tess.nodes)
    box_obb = BoxSpline(lo=lo_obb, hi=hi_obb, n_cells=(11, 7, 4), degree=DEGREE)
    prep_obb = _build_trimmed_and_classify(pts_local, tess.tris, box_obb, bvol)
    assert prep_obb["n_fallback"] == 0
    assert prep_obb["vol_err"] < 0.02
    assert (
        prep_obb["n_master"] + prep_obb["n_extended"] + prep_obb["n_inactive"]
        == box_obb.n_ctrl
    )
    # PCA OBB packs the inclined bracket with substantially higher active-cell ratio
    assert prep_obb["active_cell_ratio"] > 1.5 * prep_aa["active_cell_ratio"]

    bench_aa = _benchmark_kratos_solvers("GB_AA", prep_aa, axis_idx=2)
    assert bench_aa["peak_rss_mb"] < RSS_LIMIT_MB
    for row in bench_aa["solvers"]:
        assert row["rel_l2_diff"] < 1e-7, row
        assert row["rxn_sum_norm"] < 1e-6, row

    bench_obb = _benchmark_kratos_solvers("GB_OBB", prep_obb, axis_idx=0, obb_frame=frame)
    assert bench_obb["peak_rss_mb"] < RSS_LIMIT_MB
    for row in bench_obb["solvers"]:
        assert row["rel_l2_diff"] < 1e-7, row
        assert row["rxn_sum_norm"] < 1e-6, row

    dt = time.perf_counter() - t0
    print(
        f"\n[GB AA vs OBB] tess={t_tess:.2f}s (tris={len(tess.tris)}) | "
        f"AA (8,12,14): act={prep_aa['n_active_cells']}/{prep_aa['n_total_cells']} "
        f"({prep_aa['active_cell_ratio']:.4f}), fb={prep_aa['n_fallback']}, "
        f"vol_err={prep_aa['vol_err']*100:.4f}%, M/E/I={prep_aa['n_master']}/{prep_aa['n_extended']}/{prep_aa['n_inactive']}, "
        f"act_dofs={prep_aa['active_dofs']}, trim={prep_aa['build_trimmed_table_s']:.2f}s"
    )
    print(
        f"[GB OBB (11,7,4)]: act={prep_obb['n_active_cells']}/{prep_obb['n_total_cells']} "
        f"({prep_obb['active_cell_ratio']:.4f}), fb={prep_obb['n_fallback']}, "
        f"vol_err={prep_obb['vol_err']*100:.4f}%, M/E/I={prep_obb['n_master']}/{prep_obb['n_extended']}/{prep_obb['n_inactive']}, "
        f"act_dofs={prep_obb['active_dofs']}, trim={prep_obb['build_trimmed_table_s']:.2f}s"
    )
    for b in (bench_aa, bench_obb):
        print(f"  [{b['label']}] n_eq={b['n_eq']}, nnz={b['nnz']}, asm={b['assembly_s']:.3f}s, peak_rss={b['peak_rss_mb']:.1f}MB")
        for r in b["solvers"]:
            print(
                f"    {r['solver']}: strat={r['strategy_s']:.3f}s (net={r['solve_net_s']:.3f}s), "
                f"diff={r['rel_l2_diff']:.3e}, |sum(R)|={r['rxn_sum_norm']:.3e}, "
                f"rss={r['rss_before_mb']:.1f}->{r['rss_after_mb']:.1f}MB (delta={r['rss_delta_mb']:+.1f}MB)"
            )
    print(f"[GB total wall]={dt:.2f}s")
    _cleanup(prep_aa, prep_obb, bench_aa, bench_obb, tess, geom)


def test_flange_trimmed_extension_amgcl_scaling() -> None:
    """Flange.STEP solid 0 trimmed p=2 + extension AMGCL scaling at ~5.3k and ~8.5k active DOFs."""
    t0 = time.perf_counter()
    geom = load_cad_file(_RESOURCES / "Flange.STEP", heal=False)
    solid = geom.solids[0]
    bvol = cad.brep_volume(solid)
    lo, hi = cad.shape_bbox(solid)
    span = hi - lo
    diag = float(np.linalg.norm(span))
    t_tess0 = time.perf_counter()
    tess = cad.tessellate(solid, deflection=0.003 * diag)
    t_tess = time.perf_counter() - t_tess0
    assert cad.boundary_edge_count(tess) == 0

    configs = [
        ("Flange_5k", (22, 10, 11), 0.05, 0),
        ("Flange_8k5", (27, 12, 13), 0.04, 2),
    ]
    for label, grid, pad_frac, max_fb in configs:
        pad = pad_frac * span
        box = BoxSpline(
            lo=tuple((lo - pad).tolist()),
            hi=tuple((hi + pad).tolist()),
            n_cells=grid,
            degree=DEGREE,
        )
        prep = _build_trimmed_and_classify(tess.nodes, tess.tris, box, bvol)
        assert prep["n_fallback"] <= max_fb
        assert prep["vol_err"] < 0.02
        assert 3000 <= prep["active_dofs"] <= 15000
        assert (
            prep["n_master"] + prep["n_extended"] + prep["n_inactive"]
            == box.n_ctrl
        )

        bench = _benchmark_kratos_solvers(label, prep, axis_idx=1)
        assert bench["peak_rss_mb"] < RSS_LIMIT_MB
        for row in bench["solvers"]:
            assert row["rel_l2_diff"] < 1e-7, row
            assert row["rxn_sum_norm"] < 1e-6, row

        print(
            f"\n[{label} {grid}] tess={t_tess:.2f}s (tris={len(tess.tris)}), "
            f"act={prep['n_active_cells']}/{prep['n_total_cells']} ({prep['active_cell_ratio']:.4f}), "
            f"cut={prep['n_cut']}, fb={prep['n_fallback']}, vol_err={prep['vol_err']*100:.4f}%, "
            f"M/E/I={prep['n_master']}/{prep['n_extended']}/{prep['n_inactive']}, act_dofs={prep['active_dofs']}"
        )
        print(
            f"  stages: trim={prep['build_trimmed_table_s']:.3f}s, "
            f"cls={prep['classify_basis_s']:.4f}s, ext={prep['build_extension_s']:.4f}s, "
            f"vol={bench['stage_timings']['make_vol_s']:.3f}s, "
            f"elem={bench['stage_timings']['attach_elem_s']:.3f}s, "
            f"ties={bench['stage_timings']['apply_ties_s']:.3f}s, "
            f"asm={bench['assembly_s']:.3f}s"
        )
        print(f"  system: n_eq={bench['n_eq']}, nnz={bench['nnz']}, peak_rss={bench['peak_rss_mb']:.1f}MB")
        for r in bench["solvers"]:
            print(
                f"    {r['solver']}: strat={r['strategy_s']:.3f}s (net={r['solve_net_s']:.3f}s), "
                f"diff={r['rel_l2_diff']:.3e}, |sum(R)|={r['rxn_sum_norm']:.3e}, "
                f"rss={r['rss_before_mb']:.1f}->{r['rss_after_mb']:.1f}MB (delta={r['rss_delta_mb']:+.1f}MB)"
            )
        _cleanup(prep, bench, box)

    dt = time.perf_counter() - t0
    print(f"[Flange total wall]={dt:.2f}s")
    _cleanup(tess, solid, geom)
