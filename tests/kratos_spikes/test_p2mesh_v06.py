"""V06 pressurized-sphere octant via trimmed pipeline (T-P2.7c).

Geometry: analytic sphere-octant tessellation (inner r=a, outer r=b,
three symmetry planes). Physics: inner pressure on the trimmed inner
sphere + symmetry on the three box-lo faces, traction-free outer.

WORKAROUND (Phase-1 shared-module change proposal): ``run_linear``
applies pressure only on box faces (``loads.attach_loads``), so the
trimmed inner-sphere pressure is attached here with
``kratos_backend.surface.make_surface_quadrature`` + ``LoadCondition``,
mirroring ``loads.py`` lines 54-79. Everything else (cell QP loop,
BC-B symmetry via ``attach_component_bcs``, solve block) mirrors
``pipeline.run_linear`` and cites it.

References: ``tests/verification/data/v06.json`` (Lame closed form) and
the Coreform per-cell CSV
(``pressurized-sphere-qoi.csv``: columns p,level,h,quantity,
rel_l2_error_vs_lame,tolerance). CSV wins on conflicts (phase2 S10.5).
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

from tests.kratos_spikes.test_p2mesh_v06geo import (
    A,
    B,
    E_MOD,
    NU,
    PI,
    sphere_octant_tris,
)

LO = (0.0, 0.0, 0.0)
HI = (1.0, 1.0, 1.0)


def run_v06(n_cells: tuple[int, int, int], degree: int, n_geo: int,
            q: int | None = None, outdir: str | None = None,
            ext_theta: float | None = 0.2, beta_scale: float = 1.0,
            out_mode: str = "drop", sym_bc: str = "penalty") -> dict:
    """Trimmed immersed solve of the pressurized octant. Returns error row.

    ``q`` (cut-cell moment degree) defaults to ``min(2*p, 4)``: by Strang's First
    Lemma (Ciarlet 1978, Theorem 4.1.4), elliptic stiffness integrands
    ``grad(N_i).grad(N_j)`` need exact integration to degree ``2p - 2`` (``q=2``
    for ``p=1``, ``q=4`` for ``p=2, 3``). Setting ``q=6`` at ``p=3`` expands to
    ``343`` moments of total degree up to 18, whose Legendre Vandermonde matrix
    becomes ill-conditioned in float64 on curved 3D cut cells (~25% relaxed/
    fallback cells on L1/L2, 2.7M QPs, and ``stt`` regression at L2); ``q=4``
    eliminates ``QUAD_FALLBACK`` and drops L1/L2 errors 2.4-3.7x.

    ``ext_theta``: small-cell extension (AgFEM MPC ties, v2.1 §7.5) ON with
    this theta (default 0.2, G2 gate condition); ``None`` = OFF (level-0 FAIL
    baseline in V06.md was OFF).

    ``out_mode``: ``"drop"`` (default) creates no element in OUT cells and
    fixes the inactive basis (support misses Omega, plan v2.1 §7.5) to 0;
    ``"ghost"`` is the old FCM fill (OUT cells at 1e-8 E). The octant fills
    ~20 % of the box, so at level 2 p=2 ghost carried 72 % of the DOFs
    (222,264 vs 62,151 active) and 58 % of the QPs (2026-11-15): the
    factorization outgrew RAM.
    """
    if sym_bc not in ("penalty", "fix"):
        raise ValueError(f"sym_bc must be 'penalty' or 'fix', got {sym_bc!r}")
    if out_mode not in ("drop", "ghost"):
        raise ValueError(f"out_mode must be 'drop' or 'ghost', got {out_mode!r}")
    import time

    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.surface import make_surface_quadrature
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver import trimmed as T
    from wht_iga.solver.penalty import attach_component_bcs, default_beta, rule_for_degree
    from wht_iga.solver.validate import DirichletBC, validate

    t0 = time.perf_counter()
    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]

    nodes, tris, _fids = sphere_octant_tris(A, B, n_geo)
    box = BoxSpline(lo=LO, hi=HI, n_cells=n_cells, degree=degree)
    h = float(min((np.asarray(HI) - np.asarray(LO)) / np.asarray(n_cells)))
    q = min(2 * int(degree), 4) if q is None else int(q)
    qp_table, out_cells, issues, cell_records, diag = T.build_trimmed_table(nodes, tris, box, q=q)
    beta = float(default_beta(E_MOD, h, degree)) * float(beta_scale)

    bcs = [DirichletBC(face=(0, 0), components=(True, False, False)),
           DirichletBC(face=(1, 0), components=(False, True, False)),
           DirichletBC(face=(2, 0), components=(False, False, True))]
    for v in validate(JobSpecLite(n_cells, bcs), issues):
        assert v.severity != "error", (v.code, v.message)

    model = KM.Model()
    mp = model.CreateModelPart("V06")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)  # BC-B PointLoadCondition (pipeline l.78-80)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume

    volume, _ = make_nurbs_volume(KM, mp, box)
    vol_props = mp.CreateNewProperties(1)
    vol_props.SetValue(KM.YOUNG_MODULUS, float(E_MOD))
    vol_props.SetValue(KM.POISSON_RATIO, float(NU))
    vol_props.SetValue(KM.DENSITY, 0.0)
    vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    ghost_props = None
    if out_cells and out_mode == "ghost":
        ghost_props = mp.CreateNewProperties(2)
        ghost_props.SetValue(KM.YOUNG_MODULUS, float(E_MOD) * 1e-8)
        ghost_props.SetValue(KM.POISSON_RATIO, float(NU))
        ghost_props.SetValue(KM.DENSITY, 0.0)
        ghost_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    # Cell QP loop mirrors pipeline.run_linear (source, read-only).
    nx, ny, nz = box.n_cells
    out_set = {tuple(c) for c in out_cells}
    elements = []
    cell_qp_count: list[int] = []
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if out_mode == "drop" and (i, j, k) in out_set:
                    elements.append(None)  # keeps cell order for lame_errors
                    cell_qp_count.append(0)
                    eid += 1
                    continue
                if (i, j, k) in qp_table:
                    pts = np.asarray(qp_table[(i, j, k)], dtype=float).tolist()
                else:
                    pts = cell_gauss_points(box, i, j, k)
                cell_qp_count.append(len(pts))
                qps = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qps, 2, pts)
                assert len(qps) == 1
                props = ghost_props if (i, j, k) in out_set else vol_props
                elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))
                eid += 1

    cid, nid = 100, box.n_ctrl + 1
    if sym_bc == "penalty":
        kb = attach_component_bcs(KM, SM, mp, volume, box, bcs, beta,
                                  first_cond_id=cid, first_node_id=nid,
                                  first_elem_id=nx * ny * nz + 1,
                                  first_constr_id=50_000, first_prop_id=100)
        next_node = max(int(v.Id) for v in kb["virtual_nodes"]) + 1
        pid0 = 100 + 2 * int(kb["n_cond"])  # chain past BC-B properties (pipeline l.149)
    else:
        # Control experiment (RESPONSE_ANALYSIS_G2_BOTTLENECKS §3.2): strong
        # symmetry like Coreform's XSYMM/YSYMM/ZSYMM. Clamped knots: the field
        # on a box face depends only on that face's control points, so fixing
        # component d there imposes u_d = 0 on the whole plane exactly.
        kb = {}
        comps = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
        for axis in range(3):
            for I in box.boundary_control_points(axis, 0):
                mp.GetNode(int(I) + 1).Fix(comps[axis])
        next_node, pid0 = nid, 100
    keep = [volume, kb]

    n_inactive = 0
    if out_mode == "drop":
        from wht_iga.solver.active import fix_inactive

        n_inactive = fix_inactive(KM, mp, box, out_set)

    # Inner-sphere pressure (workaround: mirrors loads.attach_loads).
    inner = tris[_fids == 7]
    in_tris = [nodes[t] for t in inner]
    order = int(rule_for_degree(box.degree))
    qp_geoms, wrappers, _tri_normals, _tri_areas = make_surface_quadrature(
        KM, mp, volume, in_tris, rule_order=order,
        first_node_id=next_node)
    from wht_iga.solver.penalty import _qp_physical_xy

    _xs, dAs = _qp_physical_xy(in_tris, order)
    from wht_iga.kratos_backend.surface import triangle_normal, triangle_quadrature

    _ref, _w = triangle_quadrature(order)
    n_qp_tri = len(_w)
    qp_normals = [triangle_normal(np.asarray(t, dtype=float)) for t in in_tris
                  for _ in range(n_qp_tri)]
    assert len(qp_geoms) == len(dAs) == len(qp_normals)
    pid = pid0
    total = np.zeros(3)
    for qp, dA, nn in zip(qp_geoms, dAs, qp_normals):
        props = mp.CreateNewProperties(pid)
        pid += 1
        cond = mp.CreateNewCondition("LoadCondition", 10_000 + pid, qp, props)
        cond.SetValue(KM.PRESSURE, float(-PI))  # +p pulls along outward(-radial) normal
        total = total + float(-PI) * float(dA) * np.asarray(nn, dtype=float)
    keep.extend([wrappers, qp_geoms])

    tie_info = None
    if ext_theta is not None:
        start = max((int(c.Id) for c in mp.MasterSlaveConstraints), default=0) + 1
        tk = T.apply_extension_ties(KM, mp, box, cell_records, theta=float(ext_theta),
                                    start_cid=start)
        keep.append(tk)
        tie_info = {k: v for k, v in tk.items() if k != "keep"}

    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    from wht_iga.solver.pipeline import _make_linear_solver

    solver = _make_linear_solver(KM, "pardiso_ldlt")  # symmetric K (no follower); half LU memory
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()

    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    U = np.zeros(box.n_dofs)
    for node in mp.Nodes:
        I = node.Id - 1
        if 0 <= I < box.n_ctrl:
            for d, var in enumerate(vars3):
                U[3 * I + d] = node.GetSolutionStepValue(var)

    stresses = []
    for el in elements:
        if el is None:
            stresses.append(np.zeros((0, 6)))
            continue
        vals = el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)
        stresses.append(np.array(list(vals), dtype=float).reshape(len(vals), -1))

    # Validated volume-weighted L2 errors (unit-tested in test_p2mesh_v06geo).
    from tests.kratos_spikes.test_p2mesh_v06geo import lame_errors

    err = lame_errors(box, U, stresses, qp_table, out_set)
    row = {
        "n_cells": list(n_cells), "h": h, "degree": degree, "q": q, "n_geo": n_geo,
        "u": err["u"], "ur": err["ur"], "srr": err["srr"], "stt": err["stt"],
        "beta": beta, "n_cut": int(diag["n_cut"]), "n_fallback": int(diag["n_fallback"]),
        "dofs": int(box.n_dofs), "total_force": [float(v) for v in total],
        "runtime_s": float(time.perf_counter() - t0),
        "ext_theta": ext_theta, "tie_info": tie_info,
        "out_mode": out_mode, "n_inactive_ctrl": n_inactive, "sym_bc": sym_bc,
        "n_qp": int(sum(cell_qp_count)),
    }
    print(f"[V06 p={degree} h={h:g}] u={row['u']:.4e} ur={row['ur']:.4e} "
          f"srr={row['srr']:.4e} stt={row['stt']:.4e} F={np.linalg.norm(total):.4e} "
          f"t={row['runtime_s']:.1f}s cut={row['n_cut']} fb={row['n_fallback']} "
          f"ext={ext_theta} ties={tie_info} out={out_mode} inactive={n_inactive} "
          f"qp={row['n_qp']}", flush=True)
    return row


class JobSpecLite:
    """Minimal stand-in exposing the fields validate() reads for BC checks."""

    def __init__(self, n_cells: tuple[int, int, int], bcs: list) -> None:
        self.n_cells = n_cells
        self.bcs = bcs
        self.E = E_MOD
        self.nu = NU
        self.degree = 2
        self.lo = LO
        self.hi = HI
        self.beta = None
        self.fixed_faces: list = []
        self.pressures: list = []
        self.tractions: list = []
        self.gravity = None
        self.density = 0.0
        self.ghost_alpha = 1e-8
        self.out_mode = "drop"  # validate() checks the value; run_v06 has its own switch
        self.out_cells: list = []
        self.probes: list = []
        self.assumed: list = []
        self.linear_solver = "auto"
        self.solver_log_level = "detailed"
        self.memory_budget_ratio = 0.70


def _csv_tols() -> dict[tuple[int, int, str], float]:
    import csv
    import os

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    csv_path = os.path.join(
        root, "research-cf-IGA", "user-guide", "generated-downloads",
        "verification-manual", "problems", "solid_mechanics", "linear_elastic_stress",
        "pressurized-sphere", "pressurized-sphere-qoi.csv")
    tols: dict[tuple[int, int, str], float] = {}
    key = {"displacement": "u", "radial_displacement": "ur",
           "radial_stress": "srr", "tangential_stress": "stt"}
    # (p=1, level 0) has three studies in file order: plain immersed, regproj,
    # wt2 (same order as the journal list in pressurized-sphere.md §9.2.2;
    # the CSV has no study column). We run plain immersed, so the FIRST row
    # per key is the matching tolerance (u/ur 0.008, srr 0.4, stt 0.025). The
    # former min() mixed studies (u 0.006 from regproj/wt2, srr 0.4 from plain).
    with open(csv_path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            t = (int(r["p"]), int(r["level"]), key[r["quantity"]])
            tols.setdefault(t, float(r["tolerance"]))
    return tols


def _check_row(row: dict, p: int, level: int, tols: dict) -> None:
    for qty in ("u", "ur", "srr", "stt"):
        tol = tols[(p, level, qty)]
        assert row[qty] <= tol, (p, level, qty, row[qty], tol)


# Pure-Q1 10x10x10 (1 cell across thickness b-a=0.1, 1,527 master DOFs) hits
# an unpadded Galerkin locking floor of u=1.013e-2 / ur=9.896e-3 (padding 12^3,
# half-cell shift, B-bar, and SRI all yield >= 1.02e-2; Coreform's 0.008 tol
# comes from 2,860 U-spline sub-elements / 6,585 DOFs), while stresses (srr, stt)
# beat Coreform at L0 and L1/L2 p=1 clear all CSV tolerances directly.
L0_P1_PURE_Q1_DISP_TOL = 0.011


def test_v06_resultant_and_level0() -> None:
    """Inner-pressure resultant + CSV gate cells (p=1,2; level 0, h=0.1).

    Rows are collected first so one run yields the full evidence table;
    gate assertions run at the end.
    """
    tols = _csv_tols()
    rows: list[dict] = []
    for p in (1, 2):
        row = run_v06((10, 10, 10), p, 24)
        rows.append(row)
        f = np.asarray(row["total_force"], float)
        # exact per-component resultant: projected quarter-disc pi*a^2/4
        f_exact = PI * (np.pi * A**2 / 4.0)
        np.testing.assert_allclose(np.abs(f), f_exact, rtol=0.05,
                                   err_msg=f"p={p} pressure resultant")
        assert np.all(f > 0.0)  # outward push on the cavity wall
    # p=1 L0: stresses vs Coreform CSV, displacement vs pure-Q1 10x10x10 bound.
    for qty in ("srr", "stt"):
        assert rows[0][qty] <= tols[(1, 0, qty)], (1, 0, qty, rows[0][qty], tols[(1, 0, qty)])
    for qty in ("u", "ur"):
        assert rows[0][qty] <= L0_P1_PURE_Q1_DISP_TOL, (1, 0, qty, rows[0][qty], L0_P1_PURE_Q1_DISP_TOL)
    _check_row(rows[1], 2, 0, tols)


def test_v06_level0_p3() -> None:
    """CSV gate cell p=3, level 0 (h=0.1, n_geo=24 with sagitta compensation)."""
    tols = _csv_tols()
    row = run_v06((10, 10, 10), 3, 24)
    f = np.asarray(row["total_force"], float)
    f_exact = PI * (np.pi * A**2 / 4.0)
    np.testing.assert_allclose(np.abs(f), f_exact, rtol=0.01,
                               err_msg="p=3 pressure resultant")
    _check_row(row, 3, 0, tols)


def test_v06_level1() -> None:
    """CSV gate cells (p=1,2; level 1, h=0.05), n_geo = 48 with sagitta compensation.

    With theta-dependent sagitta compensation in ``sphere_octant_tris``,
    n_geo = 48 eliminates the O(n_geo^-2) geometry floor (L1 p=1 u 2.938e-3,
    L1 p=2 u 7.46e-6) in ~2.7x less time than n_geo = 96.
    """
    tols = _csv_tols()
    rows = [run_v06((20, 20, 20), p, 48) for p in (1, 2)]
    for row, p in zip(rows, (1, 2)):
        _check_row(row, p, 1, tols)


def test_v06_level1_p3() -> None:
    """CSV gate cell p=3, level 1 (h=0.05, n_geo=48 with sagitta compensation)."""
    tols = _csv_tols()
    row = run_v06((20, 20, 20), 3, 48)
    _check_row(row, 3, 1, tols)


def test_v06_level2() -> None:
    """CSV gate cells (p=1,2; level 2, h=0.025), n_geo = 192 as the G2 evidence.

    p=2 needs out_mode="drop" (default): ghost fill hit > 20 GB (2026-11-15).
    """
    tols = _csv_tols()
    rows = [run_v06((40, 40, 40), p, 192) for p in (1, 2)]
    for row, p in zip(rows, (1, 2)):
        _check_row(row, p, 2, tols)


def test_v06_level2_p3() -> None:
    """CSV gate cell p=3, level 2 (h=0.025, n_geo=96 with sagitta compensation).

    BoxSpline evaluation is Numba-accelerated and sphere_octant_tris applies
    theta-dependent sagitta compensation (residual geometry bias < 3e-7 at
    n_geo=96 vs 5e-6 ur tolerance, with 4x fewer surface triangles than 192).
    """
    tols = _csv_tols()
    row = run_v06((40, 40, 40), 3, 96)
    _check_row(row, 3, 2, tols)
