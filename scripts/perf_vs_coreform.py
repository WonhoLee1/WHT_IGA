# -*- coding: utf-8 -*-
"""Phase 6 G-8 Performance Profiling & Iterative Solver Harness vs Coreform.

Measures per-stage wall time and RSS memory across compact-to-medium 3D B-spline
models (p=1, p=2), benchmarks Kratos direct vs iterative solvers (including
IgaApplication AdditiveSchwarzPreconditioner and AMGCL variants), profiles
``wht_iga.analyses.fbar`` (``residual_tangent_sparse`` and ``solve_static`` with
``pypardiso`` vs ``superlu``) at E04/V07 sizing, and compares/extrapolates against
Coreform's ``glass-to-metal-seal`` baselines and the Phase 6 1M-DOF gate
(< 30 min, < 64 GB).

Peak RSS is kept strictly < 1.0 GB via compact/medium scaling points (1k-15k
DOFs) and explicit cleanup + ``gc.collect()`` between runs.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import psutil
import scipy.sparse as sp
import scipy.sparse.linalg as spla

_REPO_WHT_IGA = Path(__file__).resolve().parents[1]
if str(_REPO_WHT_IGA) not in sys.path:
    sys.path.insert(0, str(_REPO_WHT_IGA))

from wht_iga.analyses import fbar as FB
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import n_control_points, uniform_clamped
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements, system_matrix_csr
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.solvers import solve_linear_csr
from wht_iga.linalg.pardiso import PardisoFactor
from wht_iga.reference.box import BoxSpline
from wht_iga.reference.cook import CookSpline, HR_MM, T_MM
from wht_iga.solver import trimmed as T
from wht_iga.solver.active import fix_inactive
from wht_iga.solver.memory import est_amgcl_bytes, est_pardiso_bytes, probe_system_memory
from wht_iga.stabilize import extension as EXT
from wht_iga.stabilize import support as SUP


_PROC = psutil.Process()


def get_rss_mb() -> float:
    """Return current process resident set size (RSS) in decimal MB (10^6 bytes)."""
    return float(_PROC.memory_info().rss) / 1e6


def get_rss_mib() -> float:
    """Return current process resident set size (RSS) in MiB (1024^2 bytes)."""
    return float(_PROC.memory_info().rss) / (1024.0 ** 2)


def _cleanup(*objs: Any) -> float:
    """Drop references, run ``gc.collect()``, and return post-GC RSS in MB."""
    del objs
    gc.collect()
    return get_rss_mb()


def load_coreform_baselines(repo_root: Path | None = None) -> dict[str, Any]:
    """Load Coreform glass-to-metal-seal performance and problem-size CSV baselines.

    Parameters
    ----------
    repo_root:
        Optional path to ``WHT_IGA/``. Defaults to the repository root inferred
        from this script's location.

    Returns
    -------
    dict[str, Any]
        Parsed performance timings, problem sizes, and derived throughput metrics.
    """
    root = Path(repo_root) if repo_root is not None else _REPO_WHT_IGA
    base = (
        root
        / "research-cf-IGA"
        / "user-guide"
        / "generated-downloads"
        / "example-problems"
        / "glass-to-metal-seal"
    )
    perf_csv = base / "glass-to-metal-seal-performance.csv"
    size_csv = base / "glass-to-metal-seal-problem-size.csv"

    phases: dict[str, dict[str, Any]] = {}
    with open(perf_csv, "r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            phase_name = row["phase"].strip()
            phases[phase_name] = {
                "resources": row["resources"].strip(),
                "wall_s": float(row["wall_s"]),
            }

    sizes: dict[str, int] = {}
    with open(size_csv, "r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            sizes[row["metric"].strip()] = int(float(row["value"]))

    mesh_s = float(phases["IGA mesh"]["wall_s"])
    solve_s = float(phases["Abaqus solve"]["wall_s"])
    total_s = float(sum(v["wall_s"] for v in phases.values()))
    bg_elems = int(sizes["Background elements"])
    active_elems = int(sizes["Active elements"])
    cut_elems = int(sizes["Cut elements"])
    abaqus_eq = int(sizes["Abaqus equations"])

    return {
        "phases": phases,
        "sizes": sizes,
        "cae_script_s": float(phases["CAE script"]["wall_s"]),
        "iga_mesh_s": mesh_s,
        "interop_s": float(phases["Interop"]["wall_s"]),
        "abaqus_solve_s": solve_s,
        "total_s": total_s,
        "background_elements": bg_elems,
        "active_elements": active_elems,
        "cut_elements": cut_elems,
        "abaqus_elements": int(sizes["Abaqus elements"]),
        "abaqus_nodes": int(sizes["Abaqus nodes"]),
        "abaqus_variables": int(sizes["Abaqus variables"]),
        "abaqus_equations": abaqus_eq,
        "mesh_rate_bg_elem_per_s": bg_elems / mesh_s,
        "mesh_rate_active_elem_per_s": active_elems / mesh_s,
        "solve_rate_eq_per_s": abaqus_eq / solve_s,
        "phase6_gate": {
            "target_dofs": 1_000_000,
            "max_wall_s": 1800.0,
            "max_ram_gb": 64.0,
            "reference_cores": 16,
        },
        "s16_reference": {
            "degree": 2,
            "n_elements": 97_336,
            "n_nodes": 110_592,
            "creation_s": 52.5,
            "rate_elem_per_s": 97_336 / 52.5,
            "rss_delta_gb": 2.92,
        },
    }


def _sphere_octant_tessellation(
    r_in: float = 0.45,
    r_out: float = 0.95,
    n_seg: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a closed, outward-wound triangulated spherical shell octant in [0, 1]^3.

    Parameters
    ----------
    r_in:
        Inner radius of the spherical shell.
    r_out:
        Outer radius of the spherical shell (kept < 1.0 so outer cut cells appear).
    n_seg:
        Angular/radial subdivision count per patch edge.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        ``(nodes, tris)`` arrays suitable for ``build_trimmed_table``.
    """
    th = np.linspace(0.0, 0.5 * np.pi, n_seg + 1)
    ph = np.linspace(0.0, 0.5 * np.pi, n_seg + 1)

    def _sphere_patch(radius: float, outward_sign: float) -> tuple[np.ndarray, np.ndarray]:
        idx: dict[tuple[int, int], int] = {}
        xyz: list[list[float]] = []
        tri: list[tuple[int, int, int]] = []

        def vid(i: int, j: int) -> int:
            key = (0, 0) if i == 0 else (i, j)
            if key not in idx:
                st, ct = float(np.sin(th[i])), float(np.cos(th[i]))
                idx[key] = len(xyz)
                xyz.append([radius * st * float(np.cos(ph[j])), radius * st * float(np.sin(ph[j])), radius * ct])
            return idx[key]

        for i in range(n_seg):
            for j in range(n_seg):
                v00, v10 = vid(i, j), vid(i + 1, j)
                v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
                cands = [(v00, v10, v11)] if (i == 0 and v00 == v01) else [(v00, v10, v11), (v00, v11, v01)]
                for a, b, c in cands:
                    pa, pb, pc = np.array(xyz[a]), np.array(xyz[b]), np.array(xyz[c])
                    nrm = np.cross(pb - pa, pc - pa)
                    cent = (pa + pb + pc) / 3.0
                    if float(np.dot(nrm, outward_sign * cent)) < 0.0:
                        tri.append((a, c, b))
                    else:
                        tri.append((a, b, c))
        return np.asarray(xyz, dtype=np.float64), np.asarray(tri, dtype=np.int64)

    def _plane_patch(axis: int) -> tuple[np.ndarray, np.ndarray]:
        rr = np.linspace(r_in, r_out, n_seg + 1)
        aa = np.linspace(0.0, 0.5 * np.pi, n_seg + 1)
        idx: dict[tuple[int, int], int] = {}
        xyz: list[list[float]] = []
        tri: list[tuple[int, int, int]] = []
        u = (axis + 1) % 3
        v = (axis + 2) % 3
        want = np.zeros(3, dtype=np.float64)
        want[axis] = -1.0

        def vid(i: int, j: int) -> int:
            if (i, j) not in idx:
                pt = [0.0, 0.0, 0.0]
                pt[u] = float(rr[i] * np.cos(aa[j]))
                pt[v] = float(rr[i] * np.sin(aa[j]))
                idx[(i, j)] = len(xyz)
                xyz.append(pt)
            return idx[(i, j)]

        for i in range(n_seg):
            for j in range(n_seg):
                v00, v10 = vid(i, j), vid(i + 1, j)
                v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
                for a, b, c in ((v00, v10, v11), (v00, v11, v01)):
                    pa, pb, pc = np.array(xyz[a]), np.array(xyz[b]), np.array(xyz[c])
                    if float(np.dot(np.cross(pb - pa, pc - pa), want)) < 0.0:
                        tri.append((a, c, b))
                    else:
                        tri.append((a, b, c))
        return np.asarray(xyz, dtype=np.float64), np.asarray(tri, dtype=np.int64)

    patches = [
        _sphere_patch(r_in, -1.0),
        _sphere_patch(r_out, 1.0),
        _plane_patch(0),
        _plane_patch(1),
        _plane_patch(2),
    ]
    all_xyz = []
    all_tris = []
    offset = 0
    for px, pt in patches:
        all_xyz.append(px)
        all_tris.append(pt + offset)
        offset += len(px)
    nodes = np.vstack(all_xyz)
    tris = np.vstack(all_tris)

    # Weld coinciding vertices along patch seams for watertight classification.
    rounded = np.round(nodes, decimals=11)
    _, uniq_idx, inv = np.unique(rounded, axis=0, return_index=True, return_inverse=True)
    order = np.argsort(uniq_idx)
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    welded_nodes = nodes[uniq_idx[order]]
    welded_tris = remap[inv[tris]]
    return np.ascontiguousarray(welded_nodes), np.ascontiguousarray(welded_tris, dtype=np.int64)


def _apply_cube_bcs(KM: Any, mp: Any, box: BoxSpline, u_prescribed: float = 0.01) -> None:
    """Fix ``x=0`` face completely and prescribe ``u_x = u_prescribed`` on ``x=1``."""
    for idx in box.boundary_control_points(0, 0):
        nd = mp.GetNode(int(1 + idx))
        nd.Fix(KM.DISPLACEMENT_X)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
    for idx in box.boundary_control_points(0, 1):
        nd = mp.GetNode(int(1 + idx))
        if not nd.IsFixed(KM.DISPLACEMENT_X):
            nd.Fix(KM.DISPLACEMENT_X)
            nd.SetSolutionStepValue(KM.DISPLACEMENT_X, 0, float(u_prescribed))


def _make_kratos_solver(KM: Any, LA: Any, IGA: Any, solver_key: str, tol: float, max_iter: int) -> Any:
    """Instantiate a Kratos LinearSolver by benchmark key."""
    if solver_key == "PardisoLU":
        return LA.PardisoLUSolver()
    if solver_key == "PardisoLDLT":
        return LA.PardisoLDLTSolver()
    if solver_key == "SkylineLU":
        return KM.SkylineLUFactorizationSolver()
    if solver_key == "SparseLU":
        return LA.SparseLUSolver()
    if solver_key == "AMGCL_CG_ILU0":
        return KM.AMGCLSolver(
            KM.Parameters(
                json.dumps(
                    {
                        "solver_type": "amgcl",
                        "krylov_type": "cg",
                        "smoother_type": "ilu0",
                        "coarsening_type": "aggregation",
                        "max_iteration": int(max_iter),
                        "tolerance": float(tol),
                        "verbosity": 0,
                        "use_block_matrices_if_possible": False,
                    }
                )
            )
        )
    if solver_key == "AMGCL_CG_SPAI0_Block":
        return KM.AMGCLSolver(
            KM.Parameters(
                json.dumps(
                    {
                        "solver_type": "amgcl",
                        "krylov_type": "cg",
                        "smoother_type": "spai0",
                        "coarsening_type": "smoothed_aggregation",
                        "max_iteration": int(max_iter),
                        "tolerance": float(tol),
                        "verbosity": 0,
                        "use_block_matrices_if_possible": True,
                    }
                )
            )
        )
    if solver_key == "AMGCL_GMRES_ILU0":
        return KM.AMGCLSolver(
            KM.Parameters(
                json.dumps(
                    {
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
                )
            )
        )
    if solver_key == "CG_Diag":
        return KM.CGSolver(float(tol), int(max_iter), KM.DiagonalPreconditioner())
    if solver_key == "CG_ILU0":
        return KM.CGSolver(float(tol), int(max_iter), KM.ILU0Preconditioner())
    if solver_key == "CG_AdditiveSchwarz":
        return KM.CGSolver(float(tol), int(max_iter), IGA.AdditiveSchwarzPreconditioner())
    if solver_key == "Scaling_CG_Diag":
        return KM.ScalingSolver(KM.CGSolver(float(tol), int(max_iter), KM.DiagonalPreconditioner()), True)
    if solver_key == "Scaling_AMGCL_CG":
        inner = _make_kratos_solver(KM, LA, IGA, "AMGCL_CG_ILU0", tol, max_iter)
        return KM.ScalingSolver(inner, True)
    raise ValueError(f"Unknown solver_key {solver_key!r}")


def _extract_u_vector(KM: Any, mp: Any, n_ctrl: int) -> np.ndarray:
    """Extract control-point displacement vector ``U`` (length ``3 * n_ctrl``)."""
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    u = np.zeros(3 * n_ctrl, dtype=np.float64)
    for node in mp.Nodes:
        idx = int(node.Id) - 1
        if 0 <= idx < n_ctrl:
            for d, var in enumerate(vars3):
                u[3 * idx + d] = float(node.GetSolutionStepValue(var))
    return u


def _measure_kratos_assembly_and_solve(
    KM: Any,
    LA: Any,
    mp: Any,
    box: BoxSpline,
    use_pardiso: bool = True,
) -> dict[str, Any]:
    """Measure granular Kratos setup, BuildRHS, Build+Dirichlet, to_csr, and Strategy Solve."""
    import KratosMultiphysics.scipy_conversion_tools as sct

    mp.SetBufferSize(2)
    mp.CloneTimeStep(1.0)

    # Granular assembly breakdown using a dedicated builder/strategy instance.
    scheme_g = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver_g = LA.PardisoLUSolver() if use_pardiso else KM.SkylineLUFactorizationSolver()
    builder_g = KM.ResidualBasedBlockBuilderAndSolver(solver_g)
    strategy_g = KM.ResidualBasedLinearStrategy(mp, scheme_g, builder_g, False, False, False, False)
    strategy_g.SetEchoLevel(0)
    strategy_g.Initialize()

    A = KM.CompressedMatrix()
    Dx = KM.Vector()
    b = KM.Vector()

    t0 = time.perf_counter()
    builder_g.SetUpDofSet(scheme_g, mp)
    builder_g.SetUpSystem(mp)
    builder_g.ResizeAndInitializeVectors(scheme_g, A, Dx, b, mp)
    scheme_g.InitializeSolutionStep(mp, A, Dx, b)
    builder_g.InitializeSolutionStep(mp, A, Dx, b)
    scheme_g.InitializeNonLinIteration(mp, A, Dx, b)
    t_setup = time.perf_counter() - t0

    b_rhs = KM.Vector(len(b), 0.0)
    t0 = time.perf_counter()
    builder_g.BuildRHS(scheme_g, mp, b_rhs)
    t_residual = time.perf_counter() - t0
    del b_rhs

    t0 = time.perf_counter()
    builder_g.Build(scheme_g, mp, A, b)
    if mp.NumberOfMasterSlaveConstraints() > 0:
        builder_g.ApplyConstraints(scheme_g, mp, A, b)
    builder_g.ApplyDirichletConditions(scheme_g, mp, A, Dx, b)
    t_build = time.perf_counter() - t0

    t0 = time.perf_counter()
    K_csr = sct.to_csr(A)
    t_csr = time.perf_counter() - t0
    n_eq = int(K_csr.shape[0])
    nnz = int(K_csr.nnz)

    del K_csr, A, Dx, b, strategy_g, builder_g, solver_g, scheme_g
    gc.collect()

    # Fresh strategy for clean end-to-end Solve() timing.
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = LA.PardisoLUSolver() if use_pardiso else KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
    strategy.SetEchoLevel(0)
    strategy.Initialize()

    rss_before_solve = get_rss_mb()
    t0 = time.perf_counter()
    strategy.Solve()
    t_strategy = time.perf_counter() - t0
    rss_after_solve = get_rss_mb()
    t_solve_net = max(0.0, t_strategy - t_build)

    del strategy, builder, solver, scheme
    gc.collect()

    return {
        "n_eq": n_eq,
        "nnz": nnz,
        "kratos_setup_s": t_setup,
        "kratos_residual_s": t_residual,
        "kratos_build_s": t_build,
        "kratos_system_matrix_csr_s": t_csr,
        "kratos_strategy_solve_s": t_strategy,
        "kratos_solve_net_s": t_solve_net,
        "rss_before_solve_mb": rss_before_solve,
        "rss_after_solve_mb": rss_after_solve,
    }


def profile_untrimmed_stages(
    n_cells: tuple[int, int, int],
    degree: int,
    use_pardiso: bool = True,
) -> dict[str, Any]:
    """Profile untrimmed 3D B-spline stages (S16 element rate + Kratos assembly/solve).

    Parameters
    ----------
    n_cells:
        Cell counts ``(nx, ny, nz)`` on the unit cube.
    degree:
        B-spline polynomial degree (e.g. 1 or 2).
    use_pardiso:
        Whether to use ``PardisoLUSolver`` (True) or ``SkylineLUFactorizationSolver``.
    """
    gc.collect()
    rss_start = get_rss_mb()
    peak_rss = rss_start

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    LA = mods["LinearSolversApplication"]
    peak_rss = max(peak_rss, get_rss_mb())

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=int(degree))
    model = KM.Model()
    mp = model.CreateModelPart("PerfUntrimmed")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    r0 = get_rss_mb()
    t0 = time.perf_counter()
    volume, _nodes = make_nurbs_volume(KM, mp, box)
    t_vol = time.perf_counter() - t0
    r1 = get_rss_mb()
    peak_rss = max(peak_rss, r1)

    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, 1000.0)
    props.SetValue(KM.POISSON_RATIO, 0.3)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    t0 = time.perf_counter()
    elements = attach_cell_elements(KM, mp, volume, box, props)
    t_elem = time.perf_counter() - t0
    r2 = get_rss_mb()
    peak_rss = max(peak_rss, r2)

    _apply_cube_bcs(KM, mp, box, u_prescribed=0.01)
    asm = _measure_kratos_assembly_and_solve(KM, LA, mp, box, use_pardiso=use_pardiso)
    peak_rss = max(peak_rss, asm["rss_before_solve_mb"], asm["rss_after_solve_mb"])

    n_el = len(elements)
    del elements, _nodes, volume, props, mp, model
    rss_end = _cleanup()

    return {
        "mode": "untrimmed",
        "degree": int(degree),
        "n_cells": list(n_cells),
        "n_elements": int(n_el),
        "n_ctrl": int(box.n_ctrl),
        "n_dofs": int(box.n_dofs),
        "make_nurbs_volume_s": float(t_vol),
        "make_nurbs_volume_rss_delta_mb": float(r1 - r0),
        "attach_cell_elements_s": float(t_elem),
        "attach_cell_elements_rss_delta_mb": float(r2 - r1),
        "element_rate_eps": float(n_el / max(t_elem, 1e-12)),
        **asm,
        "rss_start_mb": float(rss_start),
        "rss_end_mb": float(rss_end),
        "peak_rss_mb": float(peak_rss),
    }


def profile_trimmed_stages(
    n_cells: tuple[int, int, int],
    degree: int,
    n_geo: int = 12,
    q: int | None = None,
    theta: float = 0.2,
    workers: int | None = None,
) -> dict[str, Any]:
    """Profile full trimmed 3D B-spline pipeline stages on a spherical-shell octant."""
    gc.collect()
    rss_start = get_rss_mb()
    peak_rss = rss_start

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    LA = mods["LinearSolversApplication"]
    peak_rss = max(peak_rss, get_rss_mb())

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=int(degree))
    q_eff = min(2 * int(degree), 4) if q is None else int(q)

    t0 = time.perf_counter()
    nodes, tris = _sphere_octant_tessellation(r_in=0.45, r_out=0.95, n_seg=int(n_geo))
    t_tess = time.perf_counter() - t0

    r0 = get_rss_mb()
    t0 = time.perf_counter()
    qp_table, out_cells, _issues, cell_records, diag = T.build_trimmed_table(
        nodes, tris, box, q=q_eff, workers=workers, cache_dir="off"
    )
    t_trim = time.perf_counter() - t0
    r1 = get_rss_mb()
    peak_rss = max(peak_rss, r1)

    nx, ny, nz = (int(v) for v in box.n_cells)
    ncell = nx * ny * nz
    h = (np.asarray(box.hi, float) - np.asarray(box.lo, float)) / np.asarray([nx, ny, nz])
    cv = float(np.prod(h))
    eta = np.zeros(ncell, dtype=np.float64)
    for c in cell_records:
        eta[int(c["i"]) + nx * (int(c["j"]) + ny * int(c["k"]))] = float(c.get("eta", 1.0))
    knots = [np.asarray(k, dtype=float) for k in box.knots]

    t0 = time.perf_counter()
    cls = SUP.classify_basis(tuple(int(v) for v in box.n_cp), knots, eta * cv, np.full(ncell, cv), theta=float(theta))
    t_cls = time.perf_counter() - t0

    t0 = time.perf_counter()
    ext_op = EXT.build_extension(tuple(int(v) for v in box.n_cp), cls.labels, int(degree), knots=knots)
    t_ext = time.perf_counter() - t0
    r2 = get_rss_mb()
    peak_rss = max(peak_rss, r2)

    n_master = int(cls.n_master)
    n_extended = int(cls.n_extended)
    n_inactive_cls = int(cls.n_inactive)
    n_thin = sum(1 for iss in ext_op.issues if iss.get("code") == EXT.THIN_FEATURE)

    model = KM.Model()
    mp = model.CreateModelPart("PerfTrimmed")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    t0 = time.perf_counter()
    volume, _nodes = make_nurbs_volume(KM, mp, box)
    t_vol = time.perf_counter() - t0
    peak_rss = max(peak_rss, get_rss_mb())

    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, 1000.0)
    props.SetValue(KM.POISSON_RATIO, 0.3)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    out_set = set(out_cells)
    elements = []
    eid = 1
    t0 = time.perf_counter()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (i, j, k) in out_set:
                    continue
                pts = qp_table[(i, j, k)].tolist() if (i, j, k) in qp_table else cell_gauss_points(box, i, j, k)
                qps = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qps, 2, pts)
                elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))
                eid += 1
    t_elem = time.perf_counter() - t0
    r3 = get_rss_mb()
    peak_rss = max(peak_rss, r3)

    n_inact = fix_inactive(KM, mp, box, out_cells)
    _apply_cube_bcs(KM, mp, box, u_prescribed=0.01)

    t0 = time.perf_counter()
    tie_info = T.apply_extension_ties(KM, mp, box, cell_records, theta=float(theta))
    t_ties = time.perf_counter() - t0
    n_ties = int(tie_info["n_ties"])
    peak_rss = max(peak_rss, get_rss_mb())

    asm = _measure_kratos_assembly_and_solve(KM, LA, mp, box, use_pardiso=True)
    peak_rss = max(peak_rss, asm["rss_before_solve_mb"], asm["rss_after_solve_mb"])

    n_active_cells = len(elements)
    n_active_dofs = int(3 * (box.n_ctrl - n_inact))

    del elements, tie_info, ext_op, cls, _nodes, volume, props, mp, model
    rss_end = _cleanup()

    return {
        "mode": "trimmed",
        "degree": int(degree),
        "q": int(q_eff),
        "n_cells": list(n_cells),
        "n_total_cells": int(ncell),
        "n_active_cells": int(n_active_cells),
        "n_cut_cells": int(diag["n_cut"]),
        "n_out_cells": int(len(out_cells)),
        "n_qp_cut": int(diag["n_qp_cut"]),
        "n_fallback": int(diag["n_fallback"]),
        "n_relaxed": int(diag["n_relaxed"]),
        "n_ctrl": int(box.n_ctrl),
        "n_total_dofs": int(box.n_dofs),
        "n_active_dofs": int(n_active_dofs),
        "n_master": n_master,
        "n_extended": n_extended,
        "n_inactive": n_inactive_cls,
        "n_ties": n_ties,
        "n_thin_feature": int(n_thin),
        "tessellate_s": float(t_tess),
        "build_trimmed_table_s": float(t_trim),
        "build_trimmed_table_rss_delta_mb": float(r1 - r0),
        "classify_basis_s": float(t_cls),
        "build_extension_s": float(t_ext),
        "build_extension_total_s": float(t_cls + t_ext),
        "build_extension_rss_delta_mb": float(r2 - r1),
        "make_nurbs_volume_s": float(t_vol),
        "attach_trimmed_elements_s": float(t_elem),
        "element_rate_eps": float(n_active_cells / max(t_elem, 1e-12)),
        "apply_extension_ties_s": float(t_ties),
        **asm,
        "rss_start_mb": float(rss_start),
        "rss_end_mb": float(rss_end),
        "peak_rss_mb": float(peak_rss),
    }


def benchmark_linear_solvers(
    n_cells: tuple[int, int, int],
    degree: int,
    trimmed: bool = False,
    use_extension: bool = True,
    solvers: list[str] | tuple[str, ...] | None = None,
    tol: float = 1e-10,
    max_iter: int = 2000,
    n_geo: int = 12,
    theta: float = 0.2,
    workers: int | None = None,
    include_csr_isolation: bool = True,
) -> dict[str, Any]:
    """Compare Kratos linear solvers (direct vs AMGCL vs CG/AdditiveSchwarz) on identical 3D B-spline elastic systems.

    Always runs ``PardisoLU`` first as the reference displacement solution ``U_ref``
    and computes ``rel_l2_diff_vs_pardiso`` and ``max_abs_diff_vs_pardiso`` on a
    fresh ``ModelPart`` per solver so initial guesses and preconditioner state never leak.
    """
    gc.collect()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    IGA = mods["IgaApplication"]
    LA = mods["LinearSolversApplication"]
    import KratosMultiphysics.scipy_conversion_tools as sct

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=int(degree))
    nx, ny, nz = (int(v) for v in box.n_cells)

    qp_table: dict[tuple[int, int, int], np.ndarray] | None = None
    out_cells: list[tuple[int, int, int]] = []
    cell_records: list[dict[str, Any]] | None = None
    if trimmed:
        q_eff = min(2 * int(degree), 4)
        nodes, tris = _sphere_octant_tessellation(r_in=0.45, r_out=0.95, n_seg=int(n_geo))
        qp_table, out_cells, _issues, cell_records, _diag = T.build_trimmed_table(
            nodes, tris, box, q=q_eff, workers=workers, cache_dir="off"
        )
    out_set = set(out_cells)

    if solvers is None:
        base_solvers = [
            "PardisoLU",
            "PardisoLDLT",
            "AMGCL_CG_ILU0",
            "AMGCL_CG_SPAI0_Block",
            "AMGCL_GMRES_ILU0",
            "CG_Diag",
            "CG_ILU0",
            "CG_AdditiveSchwarz",
            "Scaling_CG_Diag",
            "Scaling_AMGCL_CG",
        ]
        if box.n_dofs <= 8500:
            base_solvers.insert(2, "SkylineLU")
        solver_list = base_solvers
    else:
        solver_list = list(solvers)
        if "PardisoLU" not in solver_list:
            solver_list = ["PardisoLU"] + solver_list

    # Measure deterministic assembly time once on the exact same model setup.
    u_ref: np.ndarray | None = None
    k_csr_ref: sp.csr_matrix | None = None
    b_np_ref: np.ndarray | None = None
    assembly_ref_s = 0.0
    n_eq = 0
    nnz = 0
    peak_rss = get_rss_mb()

    results: list[dict[str, Any]] = []
    for s_idx, solver_key in enumerate(solver_list):
        gc.collect()
        rss_before = get_rss_mb()
        model = KM.Model()
        mp = model.CreateModelPart(f"Bench_{s_idx}")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        volume, _nodes = make_nurbs_volume(KM, mp, box)
        props = mp.CreateNewProperties(1)
        props.SetValue(KM.YOUNG_MODULUS, 1000.0)
        props.SetValue(KM.POISSON_RATIO, 0.3)
        props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

        eid = 1
        elements = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if trimmed and (i, j, k) in out_set:
                        continue
                    if trimmed and qp_table is not None and (i, j, k) in qp_table:
                        pts = qp_table[(i, j, k)].tolist()
                    else:
                        pts = cell_gauss_points(box, i, j, k)
                    qps = KM.GeometriesVector()
                    volume.CreateQuadraturePointGeometries(qps, 2, pts)
                    elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))
                    eid += 1

        if trimmed:
            fix_inactive(KM, mp, box, out_cells)
        _apply_cube_bcs(KM, mp, box, u_prescribed=0.01)
        tie_keep = None
        if trimmed and use_extension and cell_records is not None:
            tie_keep = T.apply_extension_ties(KM, mp, box, cell_records, theta=float(theta))

        mp.SetBufferSize(2)
        mp.CloneTimeStep(1.0)

        # On first (PardisoLU) run, measure standalone Build+Constraints+Dirichlet assembly time
        # using a separate builder_g/strategy_g so the main strategy runs on clean state.
        if s_idx == 0:
            scheme_g = KM.ResidualBasedIncrementalUpdateStaticScheme()
            solver_g = _make_kratos_solver(KM, LA, IGA, solver_key, tol=tol, max_iter=max_iter)
            builder_g = KM.ResidualBasedBlockBuilderAndSolver(solver_g)
            strategy_g = KM.ResidualBasedLinearStrategy(mp, scheme_g, builder_g, False, False, False, False)
            strategy_g.SetEchoLevel(0)
            strategy_g.Initialize()
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
            assembly_ref_s = time.perf_counter() - tb0
            k_csr_ref = sct.to_csr(A_tmp)
            b_np_ref = np.array([b_tmp[i] for i in range(len(b_tmp))], dtype=np.float64)
            n_eq = int(k_csr_ref.shape[0])
            nnz = int(k_csr_ref.nnz)
            del A_tmp, Dx_tmp, b_tmp, strategy_g, builder_g, solver_g, scheme_g
            gc.collect()

        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        solver = _make_kratos_solver(KM, LA, IGA, solver_key, tol=tol, max_iter=max_iter)
        builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
        strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()

        t0 = time.perf_counter()
        err_msg: str | None = None
        try:
            strategy.Solve()
            dt_strategy = time.perf_counter() - t0
            rss_after = get_rss_mb()
            peak_rss = max(peak_rss, rss_after)
            u_sol = _extract_u_vector(KM, mp, box.n_ctrl)
            try:
                it_num = int(solver.GetIterationsNumber())
            except Exception:
                it_num = 0
            try:
                kratos_conv = bool(strategy.IsConverged())
            except Exception:
                kratos_conv = True
        except Exception as exc:
            dt_strategy = time.perf_counter() - t0
            rss_after = get_rss_mb()
            peak_rss = max(peak_rss, rss_after)
            u_sol = np.zeros(box.n_dofs, dtype=np.float64)
            it_num = -1
            kratos_conv = False
            err_msg = str(exc)

        if s_idx == 0:
            u_ref = u_sol.copy()

        assert u_ref is not None
        ref_norm = max(float(np.linalg.norm(u_ref)), 1e-30)
        rel_l2 = float(np.linalg.norm(u_sol - u_ref) / ref_norm)
        max_abs = float(np.max(np.abs(u_sol - u_ref)))
        finite_ok = bool(np.all(np.isfinite(u_sol)))
        # Iterative solvers like CG_ILU0 / CG_AdditiveSchwarz may return non-converged values when diverging.
        solution_accurate = bool(finite_ok and rel_l2 < 1e-5 and err_msg is None)

        results.append(
            {
                "solver": solver_key,
                "kratos_strategy_s": float(dt_strategy),
                "assembly_ref_s": float(assembly_ref_s),
                "solve_net_s": float(max(0.0, dt_strategy - assembly_ref_s)),
                "iterations": int(it_num),
                "iterations_reported": solver_key.startswith(("CG_", "Scaling_CG")),
                "kratos_is_converged": bool(kratos_conv),
                "solution_accurate": solution_accurate,
                "rel_l2_diff_vs_pardiso": rel_l2,
                "max_abs_diff_vs_pardiso": max_abs,
                "rss_before_mb": float(rss_before),
                "rss_after_mb": float(rss_after),
                "rss_delta_mb": float(rss_after - rss_before),
                "error": err_msg,
            }
        )

        del u_sol, strategy, builder, solver, scheme, tie_keep, elements, _nodes, volume, props, mp, model
        _cleanup()

    csr_results: list[dict[str, Any]] = []
    if include_csr_isolation and k_csr_ref is not None and b_np_ref is not None:
        csr_specs = ["pardiso_lu", "pardiso_ldlt", "amgcl", "cg", "splu"]
        u_csr_ref = None
        for spec in csr_specs:
            gc.collect()
            r0 = get_rss_mb()
            res = solve_linear_csr(
                k_csr_ref,
                b_np_ref,
                solver_spec=spec,
                degree=int(degree),
                tol=float(tol),
                max_iter=int(max_iter),
            )
            r1 = get_rss_mb()
            peak_rss = max(peak_rss, r1)
            if u_csr_ref is None:
                u_csr_ref = res.solution.copy()
            rel_vs_lu = float(
                np.linalg.norm(res.solution - u_csr_ref) / max(float(np.linalg.norm(u_csr_ref)), 1e-30)
            )
            csr_results.append(
                {
                    "solver": f"csr:{spec}",
                    "wall_time_s": float(res.wall_time_s),
                    "iterations": int(res.iterations),
                    "rel_residual": float(res.rel_residual),
                    "converged": bool(res.converged),
                    "rel_l2_diff_vs_pardiso_lu": rel_vs_lu,
                    "rss_after_mb": float(r1),
                    "rss_delta_mb": float(r1 - r0),
                }
            )

    del k_csr_ref, b_np_ref, u_ref
    rss_end = _cleanup()

    return {
        "trimmed": bool(trimmed),
        "use_extension": bool(use_extension) if trimmed else False,
        "degree": int(degree),
        "n_cells": list(n_cells),
        "n_dofs_total": int(box.n_dofs),
        "n_eq": int(n_eq),
        "nnz": int(nnz),
        "assembly_ref_s": float(assembly_ref_s),
        "solvers": results,
        "csr_isolated": csr_results,
        "peak_rss_mb": float(peak_rss),
        "rss_end_mb": float(rss_end),
    }


def _cook_e04_case(n: int, degree: int = 2, nu: float = 0.49999) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    """Construct E04 Cook's membrane ``fbar`` model, consistent right-face load, and fixed DOFs."""
    box = CookSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, 1), degree=degree, continuity=None)
    model = FB.build_model(box, degree, E=240.565, nu=nu)
    cp = np.asarray(box.control_points(), dtype=np.float64)
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    kv, kw = (np.asarray(k, dtype=np.float64) for k in box.knots[1:])
    a = np.array([(kv[j + degree + 1] - kv[j]) / (degree + 1) for j in range(nv_)])
    c = np.array([(kw[k + degree + 1] - kw[k]) / (degree + 1) for k in range(nw_)])
    f = np.zeros(model["n_dof"], dtype=np.float64)
    fixed: list[int] = []
    for idx, (x, _y, _z) in enumerate(cp):
        fixed.append(3 * idx + 2)
        if abs(x) < 1e-9:
            fixed.extend([3 * idx, 3 * idx + 1])
    for k in range(nw_):
        for j in range(nv_):
            idx = (nu_ - 1) + nu_ * (j + nv_ * k)
            f[3 * idx + 1] = 6.25 * HR_MM * T_MM * a[j] * c[k]
    return model, f, np.asarray(fixed, dtype=np.int64)


def _collocate_1d(knots: np.ndarray, p: int, params: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Solve 1D B-spline interpolation ``sum_b N_b(params[m]) * c_b = values[m]``."""
    n_cp = n_control_points(knots, p)
    mat = np.zeros((n_cp, n_cp), dtype=np.float64)
    for m, xi in enumerate(params):
        span = int(find_span(n_cp, p, float(xi), knots))
        N = ders_basis_funs(span, float(xi), p, 0, knots)[0]
        mat[m, span - p : span + 1] = N
    return np.linalg.solve(mat, values)


class _MappedKirschQuarterPlate:
    """Body-fitted p=2 quarter-plate Kirsch geometry (V07 sizing)."""

    def __init__(self, n_u: int, n_v_half: int, n_w: int, degree: int = 2) -> None:
        self.n_u = int(n_u)
        self.n_v_half = int(n_v_half)
        self.n_w = int(n_w)
        self.degree = int(degree)
        a_val, b_val, t_thick = 1.0, 4.0, 0.25
        p = self.degree
        ku = uniform_clamped(0.0, 1.0, self.n_u, p, continuity=p - 1)
        br_l = np.linspace(0.0, 0.5, self.n_v_half + 1)[1:-1]
        br_r = np.linspace(0.5, 1.0, self.n_v_half + 1)[1:-1]
        kv = np.concatenate([np.full(p + 1, 0.0), br_l, np.full(p, 0.5), br_r, np.full(p + 1, 1.0)])
        kw = uniform_clamped(0.0, t_thick, self.n_w, p, continuity=p - 1)
        self.knots = [ku, kv, kw]
        self.n_cp = tuple(n_control_points(k, p) for k in self.knots)
        self.n_cells = (self.n_u, 2 * self.n_v_half, self.n_w)

        gu = greville(ku, p)
        gv = greville(kv, p)
        gw = greville(kw, p)
        th = 0.5 * np.pi * gv
        xy_in = np.column_stack([a_val * np.cos(th), a_val * np.sin(th)])
        xy_out = np.zeros((len(gv), 2), dtype=np.float64)
        for m, vi in enumerate(gv):
            xy_out[m] = [b_val, b_val * (2.0 * vi)] if vi <= 0.5 else [b_val * (2.0 * (1.0 - vi)), b_val]
        cp_in = _collocate_1d(kv, p, gv, xy_in)
        cp_out = _collocate_1d(kv, p, gv, xy_out)

        nu, nv, nw = self.n_cp
        cp = np.zeros((nu * nv * nw, 3), dtype=np.float64)
        for k_idx, z in enumerate(gw):
            for j_idx in range(nv):
                for i_idx, u_val in enumerate(gu):
                    xy = (1.0 - u_val) * cp_in[j_idx] + u_val * cp_out[j_idx]
                    idx = i_idx + nu * (j_idx + nv * k_idx)
                    cp[idx] = [xy[0], xy[1], z]
        self._cp = cp

    def control_points(self) -> np.ndarray:
        return self._cp


def _kirsch_exact_disp(x: np.ndarray, y: np.ndarray, a: float = 1.0, E: float = 100_000.0, nu: float = 0.49999, T_load: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Exact plane-strain Kirsch analytical displacement."""
    mu = E / (2.0 * (1.0 + nu))
    kappa = 3.0 - 4.0 * nu
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    q = a / r
    q3 = q ** 3
    fac = T_load * a / (8.0 * mu)
    ux = fac * ((r / a) * (kappa + 1.0) * np.cos(th) + 2.0 * q * ((kappa + 1.0) * np.cos(th) + np.cos(3.0 * th)) - 2.0 * q3 * np.cos(3.0 * th))
    uy = fac * ((r / a) * (kappa - 3.0) * np.sin(th) + 2.0 * q * ((1.0 - kappa) * np.sin(th) + np.sin(3.0 * th)) - 2.0 * q3 * np.sin(3.0 * th))
    return ux, uy


def _v07_kirsch_case(n_u: int, n_v_half: int, n_w: int) -> tuple[dict[str, Any], np.ndarray, dict[int, float]]:
    """Construct V07 Kirsch ``fbar`` model, fixed DOFs, and prescribed outer-boundary DOFs."""
    E_val, nu_val = 100_000.0, 0.49999
    plate = _MappedKirschQuarterPlate(n_u=n_u, n_v_half=n_v_half, n_w=n_w, degree=2)
    model = FB.build_model(plate, degree=2, E=E_val, nu=nu_val)
    nu, nv, nw = plate.n_cp
    cp = plate.control_points()
    fixed_set: set[int] = set()
    prescribed: dict[int, float] = {}
    for idx in range(nu * nv * nw):
        fixed_set.add(3 * idx + 2)
    for k_idx in range(nw):
        for i_idx in range(nu):
            idx_sym_y = i_idx + nu * (0 + nv * k_idx)
            idx_sym_x = i_idx + nu * ((nv - 1) + nv * k_idx)
            fixed_set.add(3 * idx_sym_y + 1)
            fixed_set.add(3 * idx_sym_x + 0)
    gv = greville(plate.knots[1], 2)
    gw = greville(plate.knots[2], 2)
    xy_out_pts = np.zeros((len(gv), 2), dtype=np.float64)
    for m, vi in enumerate(gv):
        xy_out_pts[m] = [4.0, 4.0 * (2.0 * vi)] if vi <= 0.5 else [4.0 * (2.0 * (1.0 - vi)), 4.0]
    ux_g, uy_g = _kirsch_exact_disp(xy_out_pts[:, 0], xy_out_pts[:, 1], E=E_val, nu=nu_val)
    cx = _collocate_1d(plate.knots[1], 2, gv, ux_g)
    cy = _collocate_1d(plate.knots[1], 2, gv, uy_g)
    for k_idx in range(len(gw)):
        for j_idx in range(nv):
            idx = (nu - 1) + nu * (j_idx + nv * k_idx)
            if (3 * idx + 0) not in fixed_set:
                prescribed[3 * idx + 0] = float(cx[j_idx])
            if (3 * idx + 1) not in fixed_set:
                prescribed[3 * idx + 1] = float(cy[j_idx])
    del cp
    return model, np.asarray(sorted(fixed_set), dtype=np.int64), prescribed


def profile_fbar_sizing(
    cook_sizes: tuple[int, ...] = (4, 8, 16),
    v07_levels: tuple[tuple[str, int, int, int], ...] = (
        ("h0=0.5", 6, 8, 2),
        ("h1=0.25", 12, 16, 2),
        ("h2=0.125", 24, 32, 2),
    ),
    run_superlu_max_dofs: int = 6500,
    cook_n_steps: int = 10,
) -> dict[str, Any]:
    """Profile ``fbar.residual_tangent_sparse`` and ``fbar.solve_static`` (pypardiso vs superlu) on E04/V07 sizing."""
    gc.collect()
    peak_rss = get_rss_mb()

    # Warm up Numba JIT on tiny 1x1x1 box so benchmark timings exclude one-time compile latency.
    warm_box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=2)
    warm_m = FB.build_model(warm_box, 2, E=100.0, nu=0.3)
    _ = FB.residual_tangent_sparse(warm_m, np.zeros(warm_m["n_dof"]))
    del warm_m, warm_box
    gc.collect()

    cook_rows: list[dict[str, Any]] = []
    for n in cook_sizes:
        gc.collect()
        t0 = time.perf_counter()
        model, f, fixed = _cook_e04_case(int(n), degree=2, nu=0.49999)
        t_build = time.perf_counter() - t0
        ndof = int(model["n_dof"])
        free = np.ones(ndof, dtype=bool)
        free[fixed] = False
        free_idx = np.flatnonzero(free)

        t0 = time.perf_counter()
        R0, K0 = FB.residual_tangent_sparse(model, np.zeros(ndof))
        t_asm = time.perf_counter() - t0
        nnz = int(K0.nnz)
        A_free = K0[free_idx][:, free_idx].tocsr()
        b_free = f[free_idx]

        # Single CSR linear solve timing (PardisoFactor auto vs SuperLU).
        t0 = time.perf_counter()
        pf = PardisoFactor(A_free, kind="auto")
        x_par = pf.solve(b_free)
        pf.free()
        t_lin_pardiso = time.perf_counter() - t0

        t_lin_superlu: float | None = None
        lin_diff: float | None = None
        if ndof <= run_superlu_max_dofs:
            t0 = time.perf_counter()
            x_slu = spla.spsolve(A_free.tocsc(), b_free)
            t_lin_superlu = time.perf_counter() - t0
            lin_diff = float(np.linalg.norm(x_par - x_slu) / max(float(np.linalg.norm(x_par)), 1e-30))
            del x_slu
        del R0, K0, A_free, x_par

        # Full solve_static comparison (pypardiso vs superlu).
        hist_par: list[tuple[float, int]] = []
        t0 = time.perf_counter()
        u_par = FB.solve_static(model, f, fixed, n_steps=int(cook_n_steps), linear_solver="pardiso", history=hist_par)
        t_stat_pardiso = time.perf_counter() - t0
        peak_rss = max(peak_rss, get_rss_mb())

        t_stat_superlu: float | None = None
        stat_diff: float | None = None
        if ndof <= run_superlu_max_dofs:
            hist_slu: list[tuple[float, int]] = []
            t0 = time.perf_counter()
            u_slu = FB.solve_static(model, f, fixed, n_steps=int(cook_n_steps), linear_solver="superlu", history=hist_slu)
            t_stat_superlu = time.perf_counter() - t0
            stat_diff = float(np.linalg.norm(u_par - u_slu) / max(float(np.linalg.norm(u_par)), 1e-30))
            peak_rss = max(peak_rss, get_rss_mb())
            del u_slu

        cook_rows.append(
            {
                "case": f"E04_Cook_{n}x{n}x1",
                "n_cells": [int(n), int(n), 1],
                "n_dofs": ndof,
                "n_free_dofs": int(free_idx.size),
                "nnz": nnz,
                "build_model_s": float(t_build),
                "residual_tangent_sparse_s": float(t_asm),
                "single_solve_pardiso_s": float(t_lin_pardiso),
                "single_solve_superlu_s": float(t_lin_superlu) if t_lin_superlu is not None else None,
                "single_solve_rel_diff": lin_diff,
                "solve_static_pardiso_s": float(t_stat_pardiso),
                "solve_static_superlu_s": float(t_stat_superlu) if t_stat_superlu is not None else None,
                "solve_static_rel_diff": stat_diff,
                "n_increments": int(len(hist_par)),
                "total_newton_iters": int(sum(it for _, it in hist_par)),
                "rss_mb": get_rss_mb(),
            }
        )
        del u_par, model, f, fixed
        _cleanup()

    v07_rows: list[dict[str, Any]] = []
    for label, nu_c, nv_h, nw_c in v07_levels:
        gc.collect()
        t0 = time.perf_counter()
        model, fixed, prescribed = _v07_kirsch_case(nu_c, nv_h, nw_c)
        t_build = time.perf_counter() - t0
        ndof = int(model["n_dof"])
        free = np.ones(ndof, dtype=bool)
        free[fixed] = False
        if prescribed:
            free[np.fromiter(prescribed.keys(), dtype=np.int64)] = False
        free_idx = np.flatnonzero(free)

        t0 = time.perf_counter()
        R0, K0 = FB.residual_tangent_sparse(model, np.zeros(ndof))
        t_asm = time.perf_counter() - t0
        nnz = int(K0.nnz)
        A_free = K0[free_idx][:, free_idx].tocsr()
        rhs_dummy = np.ones(free_idx.size, dtype=np.float64)

        t0 = time.perf_counter()
        pf = PardisoFactor(A_free, kind="auto")
        x_par = pf.solve(rhs_dummy)
        pf.free()
        t_lin_pardiso = time.perf_counter() - t0

        t_lin_superlu: float | None = None
        lin_diff: float | None = None
        if ndof <= run_superlu_max_dofs:
            t0 = time.perf_counter()
            x_slu = spla.spsolve(A_free.tocsc(), rhs_dummy)
            t_lin_superlu = time.perf_counter() - t0
            lin_diff = float(np.linalg.norm(x_par - x_slu) / max(float(np.linalg.norm(x_par)), 1e-30))
            del x_slu
        del R0, K0, A_free, x_par

        f_ext = np.zeros(ndof, dtype=np.float64)
        hist_par_v: list[tuple[float, int]] = []
        t0 = time.perf_counter()
        u_par = FB.solve_static(
            model, f_ext, fixed, n_steps=1, prescribed=prescribed, linear_solver="pardiso", history=hist_par_v
        )
        t_stat_pardiso = time.perf_counter() - t0
        peak_rss = max(peak_rss, get_rss_mb())

        t_stat_superlu: float | None = None
        stat_diff: float | None = None
        if ndof <= run_superlu_max_dofs:
            hist_slu_v: list[tuple[float, int]] = []
            t0 = time.perf_counter()
            u_slu = FB.solve_static(
                model, f_ext, fixed, n_steps=1, prescribed=prescribed, linear_solver="superlu", history=hist_slu_v
            )
            t_stat_superlu = time.perf_counter() - t0
            stat_diff = float(np.linalg.norm(u_par - u_slu) / max(float(np.linalg.norm(u_par)), 1e-30))
            peak_rss = max(peak_rss, get_rss_mb())
            del u_slu

        v07_rows.append(
            {
                "case": f"V07_Kirsch_{label}",
                "n_cells": [int(nu_c), int(2 * nv_h), int(nw_c)],
                "n_dofs": ndof,
                "n_free_dofs": int(free_idx.size),
                "nnz": nnz,
                "build_model_s": float(t_build),
                "residual_tangent_sparse_s": float(t_asm),
                "single_solve_pardiso_s": float(t_lin_pardiso),
                "single_solve_superlu_s": float(t_lin_superlu) if t_lin_superlu is not None else None,
                "single_solve_rel_diff": lin_diff,
                "solve_static_pardiso_s": float(t_stat_pardiso),
                "solve_static_superlu_s": float(t_stat_superlu) if t_stat_superlu is not None else None,
                "solve_static_rel_diff": stat_diff,
                "n_increments": int(len(hist_par_v)),
                "total_newton_iters": int(sum(it for _, it in hist_par_v)),
                "rss_mb": get_rss_mb(),
            }
        )
        del u_par, model, f_ext, fixed, prescribed
        _cleanup()

    return {
        "e04_cook": cook_rows,
        "v07_kirsch": v07_rows,
        "peak_rss_mb": float(peak_rss),
    }


def _fit_power_law(xs: list[float], ys: list[float]) -> tuple[float, float]:
    """Fit ``y = c * x^alpha`` via OLS in log-log space. Returns ``(c, alpha)``."""
    pairs = [(float(x), float(y)) for x, y in zip(xs, ys) if x > 0.0 and y > 0.0]
    if len(pairs) < 2:
        if len(pairs) == 1:
            return pairs[0][1] / pairs[0][0], 1.0
        return 0.0, 1.0
    lx = np.log([p[0] for p in pairs])
    ly = np.log([p[1] for p in pairs])
    alpha, intercept = np.polyfit(lx, ly, 1)
    return float(math.exp(intercept)), float(alpha)


def _eval_power_law(c: float, alpha: float, x: float) -> float:
    return float(c * (float(x) ** alpha))


def extrapolate_to_1m_and_coreform(
    untrimmed_p2: list[dict[str, Any]],
    trimmed_p2: list[dict[str, Any]],
    solver_cmp_p2: dict[str, Any],
    coreform: dict[str, Any],
    amgcl_scaling_p2: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Extrapolate measured p=2 3D B-spline timings and RAM to Coreform GTM-seal and Phase 6 1M-DOF gate."""
    # Fit untrimmed element creation rate & Kratos assembly/solve scaling vs DOFs.
    dofs_u = [float(r["n_dofs"]) for r in untrimmed_p2]
    elems_u = [float(r["n_elements"]) for r in untrimmed_p2]
    t_vol_u = [float(r["make_nurbs_volume_s"]) for r in untrimmed_p2]
    t_el_u = [float(r["attach_cell_elements_s"]) for r in untrimmed_p2]
    t_bld_u = [float(r["kratos_build_s"]) for r in untrimmed_p2]
    t_par_u = [float(r["kratos_solve_net_s"]) for r in untrimmed_p2]

    c_vol, a_vol = _fit_power_law(dofs_u, t_vol_u)
    c_el, a_el = _fit_power_law(elems_u, t_el_u)
    c_bld, a_bld = _fit_power_law(dofs_u, t_bld_u)
    c_par, a_par = _fit_power_law(dofs_u, t_par_u)

    # Trimmed stages: prefer multi-threaded points (>= 400 active cells) when >= 2 exist.
    parallel_trim = [r for r in trimmed_p2 if int(r.get("n_active_cells", 0)) >= 400]
    trim_fit_rows = parallel_trim if len(parallel_trim) >= 2 else trimmed_p2
    cut_t = [float(r["n_cut_cells"]) for r in trim_fit_rows]
    t_trim = [float(r["build_trimmed_table_s"]) for r in trim_fit_rows]
    dofs_t = [float(r["n_total_dofs"]) for r in trimmed_p2]
    t_ext = [float(r["build_extension_total_s"]) for r in trimmed_p2]
    t_ties = [float(r["apply_extension_ties_s"]) for r in trimmed_p2]

    c_trim, a_trim = _fit_power_law(cut_t, t_trim)
    c_ext, a_ext = _fit_power_law(dofs_t, t_ext)
    c_ties, a_ties = _fit_power_law(dofs_t, t_ties)

    # Measure mean element creation rate on largest p=2 untrimmed point.
    last_u = untrimmed_p2[-1]
    elem_rate_p2 = float(last_u["element_rate_eps"])

    # AMGCL scaling curve (if multi-point scaling supplied, fit empirically; else O(N^1.15) from solver_cmp_p2).
    sol_map = {r["solver"]: r for r in solver_cmp_p2.get("solvers", [])}
    par_s = float(sol_map.get("PardisoLU", {}).get("solve_net_s", 1.0))
    amg_s = float(sol_map.get("AMGCL_CG_ILU0", {}).get("solve_net_s", par_s))
    amg_ratio = amg_s / max(par_s, 1e-9)

    has_amg_fit = bool(amgcl_scaling_p2 and len(amgcl_scaling_p2) >= 2)
    if has_amg_fit:
        assert amgcl_scaling_p2 is not None
        c_amg, a_amg = _fit_power_law(
            [float(r["n_dofs"]) for r in amgcl_scaling_p2],
            [float(r["solve_s"]) for r in amgcl_scaling_p2],
        )
    else:
        n_ref_dof = max(float(solver_cmp_p2["n_dofs_total"]), 1.0)
        a_amg = 1.15
        c_amg = amg_s / (n_ref_dof ** a_amg)

    # Extrapolation targets: Coreform GTM-seal (~178,278 eq, 32,098 active, 3,420 cut) & 1M DOF (~310k cells p=2).
    targets = {}
    for label, n_dof, n_bg_el, n_act_el, n_cut_el in [
        (
            "Coreform_GTM_Seal_178k",
            int(coreform["abaqus_equations"]),
            int(coreform["background_elements"]),
            int(coreform["active_elements"]),
            int(coreform["cut_elements"]),
        ),
        (
            "Phase6_Gate_1M_DOF",
            1_000_000,
            310_000,
            310_000,
            15_000,
        ),
    ]:
        # Average 3D p=2 nonzeros per DOF observed on largest p=2 cube.
        avg_nnz_per_dof = float(last_u["nnz"]) / max(float(last_u["n_dofs"]), 1.0)
        est_nnz = int(round(n_dof * avg_nnz_per_dof))
        pardiso_gb = est_pardiso_bytes(n_dof, est_nnz, degree=2) / (1024.0 ** 3)
        amgcl_gb = est_amgcl_bytes(n_dof, est_nnz) / (1024.0 ** 3)

        vol_s = _eval_power_law(c_vol, a_vol, n_dof)
        el_s = _eval_power_law(c_el, a_el, n_act_el)
        trim_s = _eval_power_law(c_trim, a_trim, n_cut_el)
        ext_s = _eval_power_law(c_ext, a_ext, n_dof)
        ties_s = _eval_power_law(c_ties, a_ties, n_dof)
        bld_s = _eval_power_law(c_bld, a_bld, n_dof)
        solve_pardiso_s = _eval_power_law(c_par, a_par, n_dof)
        solve_amgcl_s = _eval_power_law(c_amg, a_amg, n_dof)
        total_with_amgcl_s = vol_s + el_s + trim_s + ext_s + ties_s + bld_s + solve_amgcl_s
        total_with_pardiso_s = vol_s + el_s + trim_s + ext_s + ties_s + bld_s + solve_pardiso_s

        targets[label] = {
            "n_dof": n_dof,
            "n_bg_cells": n_bg_el,
            "n_active_cells": n_act_el,
            "n_cut_cells": n_cut_el,
            "est_nnz": est_nnz,
            "est_pardiso_gb": float(pardiso_gb),
            "est_amgcl_gb": float(amgcl_gb),
            "extrapolated_stages_s": {
                "make_nurbs_volume_s": vol_s,
                "attach_cell_elements_s": el_s,
                "build_trimmed_table_s": trim_s,
                "build_extension_s": ext_s,
                "apply_extension_ties_s": ties_s,
                "kratos_build_s": bld_s,
                "solve_pardiso_s": solve_pardiso_s,
                "solve_amgcl_est_s": solve_amgcl_s,
            },
            "total_with_pardiso_s": float(total_with_pardiso_s),
            "total_with_amgcl_s": float(total_with_amgcl_s),
            "amgcl_vs_pardiso_ratio_at_measured": float(amg_ratio),
            "meets_1800s_with_amgcl": bool(total_with_amgcl_s < 1800.0),
            "meets_64gb_with_amgcl": bool(amgcl_gb < 64.0),
            "meets_64gb_with_pardiso": bool(pardiso_gb < 64.0),
        }

    # Rank hotspots at 1M-DOF projection.
    proj_1m = targets["Phase6_Gate_1M_DOF"]["extrapolated_stages_s"]
    tot_1m = targets["Phase6_Gate_1M_DOF"]["total_with_amgcl_s"]
    hotspots = sorted(
        [
            {"stage": k, "extrapolated_1m_s": float(v), "share_pct": 100.0 * float(v) / max(tot_1m, 1e-9)}
            for k, v in proj_1m.items()
            if k != "solve_pardiso_s"
        ],
        key=lambda item: item["extrapolated_1m_s"],
        reverse=True,
    )

    return {
        "fits_p2": {
            "make_nurbs_volume": {"c": c_vol, "alpha": a_vol},
            "attach_cell_elements": {"c": c_el, "alpha": a_el, "measured_rate_eps": elem_rate_p2},
            "build_trimmed_table": {"c": c_trim, "alpha": a_trim, "parallel_points_used": len(parallel_trim) >= 2},
            "build_extension": {"c": c_ext, "alpha": a_ext},
            "apply_extension_ties": {"c": c_ties, "alpha": a_ties},
            "kratos_build": {"c": c_bld, "alpha": a_bld},
            "kratos_pardiso_solve": {"c": c_par, "alpha": a_par},
            "amgcl_cg_ilu0_solve": {"c": c_amg, "alpha": a_amg, "empirical_fit": has_amg_fit},
        },
        "projections": targets,
        "hotspot_ranking_1m": hotspots,
    }


def _measure_amgcl_p2_scaling(nc_list: list[tuple[int, int, int]]) -> list[dict[str, Any]]:
    """Measure pure ``solve_linear_csr(..., solver_spec='amgcl')`` & Pardiso scaling across p=2 sizes."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    LA = mods["LinearSolversApplication"]
    import KratosMultiphysics.scipy_conversion_tools as sct

    rows: list[dict[str, Any]] = []
    for nc in nc_list:
        gc.collect()
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=nc, degree=2)
        model = KM.Model()
        mp = model.CreateModelPart("AmgScale")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        volume, _nodes = make_nurbs_volume(KM, mp, box)
        props = mp.CreateNewProperties(1)
        props.SetValue(KM.YOUNG_MODULUS, 1000.0)
        props.SetValue(KM.POISSON_RATIO, 0.3)
        props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
        els = attach_cell_elements(KM, mp, volume, box, props)
        _apply_cube_bcs(KM, mp, box, u_prescribed=0.01)
        mp.SetBufferSize(2)
        mp.CloneTimeStep(1.0)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        solver = LA.PardisoLUSolver()
        builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
        strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        A = KM.CompressedMatrix()
        Dx = KM.Vector()
        b = KM.Vector()
        builder.SetUpDofSet(scheme, mp)
        builder.SetUpSystem(mp)
        builder.ResizeAndInitializeVectors(scheme, A, Dx, b, mp)
        scheme.InitializeSolutionStep(mp, A, Dx, b)
        builder.InitializeSolutionStep(mp, A, Dx, b)
        scheme.InitializeNonLinIteration(mp, A, Dx, b)
        builder.Build(scheme, mp, A, b)
        builder.ApplyDirichletConditions(scheme, mp, A, Dx, b)
        K_csr = sct.to_csr(A)
        b_np = np.array([b[i] for i in range(len(b))], dtype=np.float64)
        del A, Dx, b, strategy, builder, solver, scheme, els, _nodes, volume, props, mp, model
        gc.collect()

        res_amg = solve_linear_csr(K_csr, b_np, solver_spec="amgcl", degree=2, tol=1e-10, max_iter=2000)
        res_par = solve_linear_csr(K_csr, b_np, solver_spec="pardiso_ldlt", degree=2, tol=1e-10)
        rows.append(
            {
                "n_cells": list(nc),
                "n_dofs": int(box.n_dofs),
                "nnz": int(K_csr.nnz),
                "solve_s": float(res_amg.wall_time_s),
                "amgcl_rel_residual": float(res_amg.rel_residual),
                "pardiso_ldlt_s": float(res_par.wall_time_s),
                "rss_mb": get_rss_mb(),
            }
        )
        del K_csr, b_np
        _cleanup()
    return rows


def run_full_benchmark(quick: bool = False) -> dict[str, Any]:
    """Run the complete G-8 benchmark suite under strict < 1.0 GB RSS limits."""
    mem_info = probe_system_memory()
    coreform = load_coreform_baselines()

    if quick:
        p1_untrimmed = [(6, 6, 6)]
        p2_untrimmed = [(6, 6, 6)]
        p1_trimmed = [(6, 6, 6)]
        p2_trimmed = [(6, 6, 6)]
        solver_untrimmed_nc = (6, 6, 6)
        solver_trimmed_nc = (6, 6, 6)
        amg_scale_ncs = [(6, 6, 6)]
        fbar_cook = (4,)
        fbar_v07 = (("h0=0.5", 6, 8, 2),)
    else:
        p1_untrimmed = [(6, 6, 6), (10, 10, 10), (14, 14, 14), (16, 16, 16)]
        p2_untrimmed = [(6, 6, 6), (9, 9, 9), (12, 12, 12), (13, 13, 13)]
        p1_trimmed = [(6, 6, 6), (10, 10, 10), (14, 14, 14)]
        p2_trimmed = [(6, 6, 6), (9, 9, 9), (12, 12, 12), (13, 13, 13)]
        solver_untrimmed_nc = (9, 9, 9)
        solver_trimmed_nc = (8, 8, 8)
        amg_scale_ncs = [(6, 6, 6), (9, 9, 9), (12, 12, 12), (13, 13, 13)]
        fbar_cook = (4, 8, 16)
        fbar_v07 = (("h0=0.5", 6, 8, 2), ("h1=0.25", 12, 16, 2), ("h2=0.125", 24, 32, 2))

    untrimmed_p1_rows = [profile_untrimmed_stages(nc, degree=1) for nc in p1_untrimmed]
    untrimmed_p2_rows = [profile_untrimmed_stages(nc, degree=2) for nc in p2_untrimmed]
    trimmed_p1_rows = [profile_trimmed_stages(nc, degree=1, n_geo=12, workers=None) for nc in p1_trimmed]
    trimmed_p2_rows = [profile_trimmed_stages(nc, degree=2, n_geo=12, workers=None) for nc in p2_trimmed]

    solver_untrimmed_p2 = benchmark_linear_solvers(solver_untrimmed_nc, degree=2, trimmed=False)
    solver_trimmed_p2_ext = benchmark_linear_solvers(
        solver_trimmed_nc, degree=2, trimmed=True, use_extension=True
    )
    solver_trimmed_p2_noext = benchmark_linear_solvers(
        solver_trimmed_nc, degree=2, trimmed=True, use_extension=False
    )
    amgcl_scaling_p2 = _measure_amgcl_p2_scaling(amg_scale_ncs)

    fbar_rows = profile_fbar_sizing(cook_sizes=fbar_cook, v07_levels=fbar_v07)
    extrap = extrapolate_to_1m_and_coreform(
        untrimmed_p2_rows, trimmed_p2_rows, solver_untrimmed_p2, coreform, amgcl_scaling_p2=amgcl_scaling_p2
    )

    all_peaks = (
        [r["peak_rss_mb"] for r in untrimmed_p1_rows]
        + [r["peak_rss_mb"] for r in untrimmed_p2_rows]
        + [r["peak_rss_mb"] for r in trimmed_p1_rows]
        + [r["peak_rss_mb"] for r in trimmed_p2_rows]
        + [solver_untrimmed_p2["peak_rss_mb"], solver_trimmed_p2_ext["peak_rss_mb"], solver_trimmed_p2_noext["peak_rss_mb"]]
        + [r["rss_mb"] for r in amgcl_scaling_p2]
        + [fbar_rows["peak_rss_mb"]]
    )

    return {
        "host": {
            "total_ram_gb": float(mem_info.total_gb),
            "available_ram_gb": float(mem_info.available_gb),
            "cpu_count": int(os.cpu_count() or 1),
            "overall_peak_rss_mb": float(max(all_peaks)),
            "overall_peak_rss_gb": float(max(all_peaks) / 1000.0),
        },
        "coreform_baselines": coreform,
        "untrimmed_p1": untrimmed_p1_rows,
        "untrimmed_p2": untrimmed_p2_rows,
        "trimmed_p1": trimmed_p1_rows,
        "trimmed_p2": trimmed_p2_rows,
        "solver_comparison_untrimmed_p2": solver_untrimmed_p2,
        "solver_comparison_trimmed_p2_ext": solver_trimmed_p2_ext,
        "solver_comparison_trimmed_p2_noext": solver_trimmed_p2_noext,
        "amgcl_scaling_p2": amgcl_scaling_p2,
        "fbar_sizing": fbar_rows,
        "extrapolation": extrap,
    }


def format_markdown_summary(report: dict[str, Any]) -> str:
    """Render concise Markdown tables from a ``run_full_benchmark`` report dict."""
    lines: list[str] = ["# WHT_IGA Phase 6 G-8 Benchmark Summary", ""]
    lines.append("## Untrimmed Element Rate & Kratos Stages")
    lines.append("| p | Grid | Elements | DOFs | make_vol (s) | attach_elem (s) | Rate (elem/s) | Build (s) | Pardiso (s) | Peak RSS (MB) |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in report.get("untrimmed_p1", []) + report.get("untrimmed_p2", []):
        lines.append(
            f"| {row['degree']} | {row['n_cells']} | {row['n_elements']} | {row['n_dofs']} | "
            f"{row['make_nurbs_volume_s']:.4f} | {row['attach_cell_elements_s']:.4f} | "
            f"{row['element_rate_eps']:.1f} | {row['kratos_build_s']:.4f} | "
            f"{row['kratos_solve_net_s']:.4f} | {row['peak_rss_mb']:.1f} |"
        )
    lines.append("")
    lines.append("## Solver Comparison (Untrimmed p=2)")
    lines.append("| Solver | Strategy (s) | Solve Net (s) | Iters | Rel L2 Diff | Accurate | RSS After (MB) |")
    lines.append("|---|---:|---:|---:|---:|---|---:|")
    for row in report.get("solver_comparison_untrimmed_p2", {}).get("solvers", []):
        lines.append(
            f"| {row['solver']} | {row['kratos_strategy_s']:.4f} | {row['solve_net_s']:.4f} | "
            f"{row['iterations']} | {row['rel_l2_diff_vs_pardiso']:.3e} | "
            f"{row['solution_accurate']} | {row['rss_after_mb']:.1f} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for G-8 performance profiling."""
    parser = argparse.ArgumentParser(description="WHT_IGA Phase 6 G-8 Performance Profiling vs Coreform")
    parser.add_argument("--quick", action="store_true", help="Run compact smoke sizing only")
    parser.add_argument("--json-out", type=str, default=None, help="Optional JSON output path (UTF-8)")
    args = parser.parse_args(argv)

    report = run_full_benchmark(quick=args.quick)
    text = json.dumps(report, indent=2)
    if args.json_out:
        out_path = Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
        print(f"Wrote G-8 benchmark JSON ({len(text)} bytes) to {out_path}")
    else:
        print(text)
    print(
        f"[G-8] Peak RSS = {report['host']['overall_peak_rss_mb']:.1f} MB "
        f"({report['host']['overall_peak_rss_gb']:.3f} GB)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
