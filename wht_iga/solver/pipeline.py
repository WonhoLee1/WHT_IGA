"""Phase-1 linear static pipeline: grid → Kratos model → solve → probe → DBs.

Entry ``run_linear(spec)`` is callable by the integrator for V01/V02/V03
acceptance. Flow (v2 §5.2–§5.4, R1): BoxSpline grid → ``make_nurbs_volume``
(keepalive held on the result) → per-cell QP elements (IN cells full
stiffness, OUT cells FCM ghost α·E, v2 §7.6) → Dirichlet (BC-A/BC-B) → loads
→ Newton/linear solve → probe → mesh/results DB + VTU + probe JSON.

Kratos is imported only inside ``run_linear`` via ``load_kratos``; module
import stays Kratos-free.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from wht_iga.reference.box import BoxSpline
from wht_iga.solver.loads import apply_gravity, attach_loads
from wht_iga.solver.memory import get_available_ram_bytes, select_linear_solver
from wht_iga.solver.penalty import (
    DEFAULT_C_PEN,
    attach_component_bcs,
    attach_penalty_bcs,
    box_face_triangles,
    default_beta,
    face_subdiv,
)
from wht_iga.solver.solver_log import SolverLogger
from wht_iga.solver.validate import JobSpec, ValidationError, ValidationIssue, validate



@dataclass
class LinearResult:
    U: np.ndarray  # (n_dofs,) control-point displacements
    beta_used: float
    total_force: np.ndarray
    face_averages: dict[tuple[int, int], np.ndarray] = field(default_factory=dict)
    probe_values: dict[str, np.ndarray] = field(default_factory=dict)
    cauchy_stress: list[np.ndarray] = field(default_factory=list)
    issues: list = field(default_factory=list)
    mesh_db: str | None = None
    results_db: str | None = None
    vtu_path: str | None = None
    probe_json: str | None = None
    keepalive: tuple = ()  # volume/wrappers/geoms/model: drop only after post-process
    strategy: object | None = None  # Kratos strategy (exposed for V-COND κ extraction)
    tie_info: dict | None = None  # extension hook summary (n_ties/n_slave/n_master/theta)
    converged: bool = True  # Kratos IsConverged() after Solve
    n_iter: int = 0  # Kratos NL_ITERATION_NUMBER after Solve (-1 if unreadable)
    timings: dict = field(default_factory=dict)  # wall seconds: setup / solve / post
    solver_log: str = ""  # 6-stage solver log text (empty if enable_solver_log=False)
    solver_log_file: str | None = None
    active_dofs: int = 0
    total_dofs: int = 0
    solver_decision: dict[str, Any] | None = None


def _cell_size(box: BoxSpline) -> float:
    br = [box.breakpoints(d) for d in range(3)]
    return float(min(float(np.diff(b).min()) for b in br))


def _next_id(container, start: int) -> int:
    """Monotonic ID allocator: one counter per Kratos container (F5).

    Returns max(existing ids) + 1 so BC/load/extension blocks can never
    collide no matter how many objects each creates (replaces the old
    fixed +10000/+40000 offsets).
    """
    try:
        return max(int(o.Id) for o in container) + 1
    except ValueError:
        return start


_LINEAR_SOLVERS = (
    "pardiso",
    "pardiso_lu",
    "pardiso_ldlt",
    "pardiso_llt",
    "skyline_lu",
    "amgcl",
    "amgcl_cg",
    "cg",
    "splu",
)
# AMGCL + CG (Kratos core): no factor fill, memory ~ matrix + AMG hierarchy. For grids
# whose LDLT factor does not fit RAM (V02 level 4: ~18 GB with LDLT, 2026-11-15).
# Symmetric positive definite K only; tolerance is the relative residual.
_AMGCL_SETTINGS = """{
    "solver_type": "amgcl", "krylov_type": "cg", "smoother_type": "ilu0",
    "coarsening_type": "aggregation", "max_iteration": 5000, "tolerance": 1e-10,
    "verbosity": 0, "use_block_matrices_if_possible": false
}"""
_PARDISO_CLASS = {
    "pardiso": "PardisoLDLTSolver",
    "pardiso_lu": "PardisoLUSolver",
    "pardiso_ldlt": "PardisoLDLTSolver",
    "pardiso_llt": "PardisoLLTSolver",
}


def _make_linear_solver(
    KM,
    name: str,
    num_threads_mkl: int | str | None = None,
    amgcl_params: dict[str, Any] | None = None,
):
    """Kratos linear solver by name. Unknown names raise (no silent fallback).

    Pardiso variants wrap Eigen's PardisoSupport (iparm fixed in C++: METIS
    ordering, 2 refinement steps, scaling/matching only for LU). The
    controllable knobs are the factorization type and ``num_threads_mkl``
    (int, "consistent" = Kratos thread count, None = MKL default).
    Kratos' block builder zeroes fixed rows AND columns, so K stays symmetric:
    "pardiso_ldlt" / "pardiso_llt" are valid for linear elasticity with penalty
    and MPC ties (TᵀKT); follower loads (unsymmetric) need "pardiso_lu".
    """
    if name in _PARDISO_CLASS:
        import KratosMultiphysics.LinearSolversApplication as LA

        if not LA.HasMKL():
            raise RuntimeError(f"linear_solver={name!r} needs MKL (LinearSolversApplication.HasMKL() is False)")
        cls = getattr(LA, _PARDISO_CLASS[name])
        if num_threads_mkl is None:
            return cls()
        prm = KM.Parameters("{}")
        if isinstance(num_threads_mkl, str):
            prm.AddEmptyValue("num_threads_mkl").SetString(num_threads_mkl)
        else:
            prm.AddEmptyValue("num_threads_mkl").SetInt(int(num_threads_mkl))
        return cls(prm)
    if name in ("skyline_lu", "amgcl", "amgcl_cg", "cg", "splu"):
        from wht_iga.kratos_backend.solvers import build_kratos_linear_solver

        merged = json.loads(_AMGCL_SETTINGS)
        if amgcl_params:
            merged.update(amgcl_params)
        return build_kratos_linear_solver(
            KM,
            name,
            num_threads_mkl=num_threads_mkl,
            amgcl_params=merged,
        )
    raise ValueError(f"unknown linear_solver {name!r}; use one of {_LINEAR_SOLVERS}")


def _run_linear_thb(spec: JobSpec) -> LinearResult:
    """THB-adaptive linear static solve (Track H Phase 2, precursor scope).

    Mirrors ``_solve_thb_adaptive`` (methods.py) but driven by ``JobSpec``:
    strong ``fixed_faces`` clamping (or penalty ``bcs``), uniform ``tractions``
    on the fine-grid faces, reduced solve via ``assemble_and_solve_thb_reduced``.
    Pressures, gravity, OUT cells and named-set regions are out of scope and
    raise :class:`ValidationError` instead of silently falling back.
    """
    import time as _time

    from wht_iga.grid.hierarchy import THBHierarchy
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import (
        AssemblyBuilder,
        AssemblySpec,
        InstanceSpec,
        PartSpec,
    )
    from wht_iga.solver.validate import DirichletBC, ValidationIssue

    _t0 = _time.perf_counter()

    def _err(code: str, message: str) -> ValidationIssue:
        return ValidationIssue(severity="error", code=code, message=message,
                               entity="thb_refinement")

    cfg = dict(spec.thb_refinement or {})
    if cfg.get("sets"):
        raise ValidationError([_err("THB_SETS_UNSUPPORTED",
                                    f"thb: named-set regions unsupported: {list(cfg['sets'])}")])
    level = cfg.get("level", 1)
    balance = cfg.get("balance", "two_to_one")
    boxes = [tuple(map(tuple, b)) for b in cfg.get("boxes", [])]
    spheres = []
    for s in cfg.get("spheres", []):
        if not (isinstance(s, (list, tuple)) and len(s) == 2):
            issues.append(_err("THB_REGION", f"thb sphere must be ((cx,cy,cz), r), got {s!r}"))
            continue
        spheres.append((tuple(float(v) for v in s[0]), float(s[1])))
    issues = []
    if not isinstance(level, int) or not (1 <= level <= 5):
        issues.append(_err("THB_LEVEL", f"thb level must be int in [1, 5], got {level!r}"))
    if balance not in ("two_to_one", "none"):
        issues.append(_err("THB_BALANCE", f"thb balance must be 'two_to_one'|'none', got {balance!r}"))
    if not boxes and not spheres:
        issues.append(_err("THB_REGION", "thb requires at least one box or sphere region"))
    if spec.out_cells:
        issues.append(_err("THB_OUT_CELLS", "thb path does not support out_cells (uniform FCM path only)"))
    if spec.pressures:
        issues.append(_err("THB_LOADS", "thb path supports tractions only (no pressures)"))
    if spec.gravity is not None:
        issues.append(_err("THB_LOADS", "thb path does not support gravity"))
    if any(not all(b.components) for b in spec.bcs):
        issues.append(_err("THB_BCS", "thb path supports full-component bcs only (no component BCs)"))
    if issues:
        raise ValidationError(issues)

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]

    lo = tuple(float(v) for v in spec.lo)
    hi = tuple(float(v) for v in spec.hi)
    thb = THBHierarchy.from_uniform(
        (lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2]),
        tuple(int(v) for v in spec.n_cells), degree=int(spec.degree))
    for r_lo, r_hi in boxes:
        thb.mark_box(tuple(float(v) for v in r_lo), tuple(float(v) for v in r_hi),
                     target_level=int(level))
    for entry in spheres:
        (cx, cy, cz), rad = entry
        thb.mark_sphere((float(cx), float(cy), float(cz)), float(rad),
                        target_level=int(level))
    if balance == "two_to_one":
        thb.enforce_2_to_1_balance()

    part = PartSpec(name="domain", kind="thb", thb=thb,
                    young_modulus=float(spec.E), poisson_ratio=float(spec.nu))
    spec_asm = AssemblySpec(job_name="thb_pipeline",
                            parts={"domain": part},
                            instances=[InstanceSpec("domain", "1")],
                            use_pardiso=False)
    builder = AssemblyBuilder(KM=KM, SM=SM)
    built = builder.build(spec_asm)
    rt = built.get_instance("domain_1")
    assert rt.box is not None
    box_fine = rt.box

    from wht_iga.analyses.methods_solve import _attach_point_loads, _fix_node_all_disp, _integrate_iga_face_traction
    from wht_iga.solver.penalty import attach_penalty_bcs, default_beta

    for ax, side in spec.fixed_faces:
        for idx in box_fine.boundary_control_points(int(ax), int(side)):
            _fix_node_all_disp(KM, built.structure.GetNode(int(rt.node_ids[int(idx)])))
    if spec.bcs:
        h0 = min(float(b) - float(a) for a, b in zip(lo, hi)) / max(1, min(int(v) for v in spec.n_cells))
        beta = float(spec.beta) if spec.beta is not None else default_beta(float(spec.E), h0, int(spec.degree), 10.0)
        attach_penalty_bcs(KM, IGA, built.structure, rt.volume, box_fine,
                           [DirichletBC(face=tuple(b.face), components=tuple(b.components),
                                        value=tuple(b.value), field=b.field) for b in spec.bcs],
                           beta,
                           first_cond_id=built.id_table.next_id("conditions") + 100,
                           first_node_id=built.id_table.next_id("nodes") + 100,
                           first_prop_id=built.id_table.next_id("properties") + 10)
    total_f = np.zeros(3, dtype=np.float64)
    for tr in spec.tractions:
        ax, side = int(tr.face[0]), int(tr.face[1])
        f_local = _integrate_iga_face_traction(box_fine, ax, side, tuple(float(v) for v in tr.vector))
        f_kratos = {int(rt.node_ids[int(flat)]): np.asarray(vec, dtype=float)
                    for flat, vec in f_local.items()}
        for vec in f_kratos.values():
            total_f = total_f + vec
        _attach_point_loads(SM, built, f_kratos, rt.properties_map[0], scope="Loads.thb")

    _t1 = _time.perf_counter()
    res = builder.assemble_and_solve_thb_reduced(built, "domain_1", use_pardiso=False)
    _t2 = _time.perf_counter()

    u_hier = np.asarray(res["u_hier"], dtype=np.float64).reshape(-1)
    n_thb = int(thb.n_thb_dofs)
    return LinearResult(
        U=u_hier,
        beta_used=float(spec.beta) if spec.beta is not None else 0.0,
        total_force=total_f,
        keepalive=(thb, box_fine, rt.node_ids, built, res),
        converged=True,
        timings={"setup_s": _t1 - _t0, "solve_s": _t2 - _t1, "total_s": _t2 - _t0},
        active_dofs=3 * n_thb,
        total_dofs=3 * n_thb,
    )


def run_linear(spec: JobSpec, geo_issues: list[dict] | None = None,
             qp_table: dict[tuple[int, int, int], np.ndarray] | None = None,
             cell_records: list[dict] | None = None,
             extension: dict | None = None,
             hooks: Any = None) -> LinearResult:
    """Linear static solve (see module docstring).

    ``hooks`` (optional, user API escape hatch, plan v2.1 §8.4): an object with
    ``before_solve(model, model_part, strategy)`` and/or
    ``after_solve(model, model_part, strategy)`` methods, called immediately
    before and after ``strategy.Solve()``. ``None`` (default) changes nothing.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume

    import time as _time

    _t0 = _time.perf_counter()
    issues = validate(spec, geo_issues)
    if any(i.severity == "error" for i in issues):
        raise ValidationError(issues)

    if spec.thb_refinement is not None:
        return _run_linear_thb(spec)

    mods = load_kratos()
    KM, SM, IGA = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"], mods["IgaApplication"]

    box = BoxSpline(lo=spec.lo, hi=spec.hi, n_cells=spec.n_cells,
                    degree=spec.degree, continuity=spec.continuity)
    h = _cell_size(box)
    beta = float(spec.beta) if spec.beta is not None else default_beta(spec.E, h, spec.degree, DEFAULT_C_PEN)

    model = KM.Model()
    mp = model.CreateModelPart("LinearJob")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    need_point_load = any(not all(b.components) for b in spec.bcs)
    if need_point_load:
        mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    volume, _nodes = make_nurbs_volume(KM, mp, box)
    vol_props = mp.CreateNewProperties(1)
    vol_props.SetValue(KM.YOUNG_MODULUS, float(spec.E))
    vol_props.SetValue(KM.POISSON_RATIO, float(spec.nu))
    vol_props.SetValue(KM.DENSITY, float(spec.density))
    vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    if spec.gravity is not None and spec.density > 0.0:
        apply_gravity(vol_props, KM, spec.density, spec.gravity)
    ghost_props = None
    drop_out = spec.out_mode == "drop"
    if spec.out_cells and not drop_out:
        ghost_props = mp.CreateNewProperties(2)
        ghost_props.SetValue(KM.YOUNG_MODULUS, float(spec.E) * float(spec.ghost_alpha))
        ghost_props.SetValue(KM.POISSON_RATIO, float(spec.nu))
        ghost_props.SetValue(KM.DENSITY, float(spec.density))
        ghost_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    # Per-cell QP elements (K2 one-call-per-cell pattern). A local loop is used
    # instead of elements.attach_cell_elements so OUT cells can carry the FCM
    # ghost-alpha properties (uniform case is identical to the helper).
    nx, ny, nz = box.n_cells
    out_set = {tuple(c) for c in spec.out_cells}
    from wht_iga.progress import Ticker, log

    elements = []
    eid = 1
    log(f"run_linear[{spec.name}]: creating {nx * ny * nz} cell elements")
    tick = Ticker("elements", nx * ny * nz, every_frac=0.25)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                tick.advance()
                if drop_out and (i, j, k) in out_set:
                    elements.append(None)  # keeps the cell order of cauchy_stress
                    eid += 1
                    continue
                if qp_table is not None and (i, j, k) in qp_table:
                    pts = np.asarray(qp_table[(i, j, k)], dtype=float).tolist()
                else:
                    pts = cell_gauss_points(box, i, j, k)
                qps = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qps, 2, pts)
                assert len(qps) == 1
                props = ghost_props if (i, j, k) in out_set else vol_props
                elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))
                eid += 1

    n_inactive = 0
    if drop_out:
        from wht_iga.solver.active import fix_inactive

        n_inactive = fix_inactive(KM, mp, box, spec.out_cells)
        log(f"run_linear[{spec.name}]: out_mode=drop, {len(out_set)} OUT cells skipped, "
            f"{n_inactive} inactive control points fixed")

    for axis, side in spec.fixed_faces:
        for I in box.boundary_control_points(axis, side):
            node = mp.GetNode(int(1 + I))
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)

    keep_parts: list = [volume]
    bcs_full = [b for b in spec.bcs if all(b.components)]
    bcs_part = [b for b in spec.bcs if not all(b.components)]

    cid = _next_id(mp.Conditions, 100)
    nid = _next_id(mp.Nodes, box.n_ctrl + 1)
    pid = _next_id(mp.Properties, 100)
    eid2 = _next_id(mp.Elements, nx * ny * nz + 1)
    tid = _next_id(mp.MasterSlaveConstraints, 1)
    if bcs_full:
        ka = attach_penalty_bcs(KM, IGA, mp, volume, box, bcs_full, beta,
                                first_cond_id=cid, first_node_id=nid, first_prop_id=pid)
        keep_parts.append(ka)
        cid, nid, pid = (_next_id(mp.Conditions, cid), _next_id(mp.Nodes, nid),
                         _next_id(mp.Properties, pid))
    if bcs_part:
        kb = attach_component_bcs(KM, SM, mp, volume, box, bcs_part, beta,
                                  first_cond_id=cid, first_node_id=nid,
                                  first_elem_id=eid2, first_constr_id=tid, first_prop_id=pid)
        keep_parts.append(kb)
        cid, nid, pid = (_next_id(mp.Conditions, cid), _next_id(mp.Nodes, nid),
                         _next_id(mp.Properties, pid))
        eid2, tid = _next_id(mp.Elements, eid2), _next_id(mp.MasterSlaveConstraints, tid)
    load_keep = attach_loads(KM, SM, IGA, mp, volume, box, spec.pressures, spec.tractions,
                             first_cond_id=cid, first_node_id=nid, first_prop_id=pid)
    keep_parts.append(load_keep)
    cid, nid, pid = (_next_id(mp.Conditions, cid), _next_id(mp.Nodes, nid),
                     _next_id(mp.Properties, pid))

    tie_info = None
    if extension is not None:
        from wht_iga.solver import trimmed as _T
        tie_keep = _T.apply_extension_ties(
            KM, mp, box, cell_records,
            theta=float(extension.get("theta", 0.2)),
            start_cid=_next_id(mp.MasterSlaveConstraints, 1))
        keep_parts.append(tie_keep)
        tie_info = {k: v for k, v in tie_keep.items() if k != "keep"}

    move_mesh = any(p.follower for p in spec.pressures)
    n_total_cells = nx * ny * nz
    n_out_count = len(out_set)
    if cell_records is not None:
        n_cut_count = sum(
            1 for r in cell_records
            if r.get("status") == 2 or 0.0 < float(r.get("eta", 1.0)) < 1.0
        )
        n_in_count = max(len(cell_records) - n_out_count - n_cut_count, 0)
    else:
        n_cut_count = 0
        n_in_count = max(n_total_cells - n_out_count, 0)
    active_ctrl = max(box.n_ctrl - n_inactive, 0)
    total_dofs = int(box.n_dofs)
    active_dofs = int(3 * active_ctrl)
    overlap_cp = (
        min(2 * box.degree + 1, nx + box.degree)
        * min(2 * box.degree + 1, ny + box.degree)
        * min(2 * box.degree + 1, nz + box.degree)
    )
    est_nnz = int(active_dofs * 3 * overlap_cp)
    n_qp_total = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if drop_out and (i, j, k) in out_set:
                    continue
                if qp_table is not None and (i, j, k) in qp_table:
                    n_qp_total += int(len(qp_table[(i, j, k)]))
                else:
                    n_qp_total += int((box.degree + 1) ** 3)

    decision = select_linear_solver(
        spec.linear_solver,
        active_dofs,
        est_nnz,
        degree=spec.degree,
        memory_budget_ratio=spec.memory_budget_ratio,
        pypardiso_available=True,
    )
    if spec.linear_solver == "auto":
        if decision.switched and not move_mesh:
            solver_name = "amgcl_cg"
        else:
            solver_name = "pardiso_lu" if move_mesh else "pardiso_ldlt"
    else:
        solver_name = spec.linear_solver

    if move_mesh and solver_name in ("pardiso", "pardiso_ldlt", "pardiso_llt", "amgcl", "amgcl_cg", "cg"):
        raise ValueError(f"follower pressure gives an unsymmetric tangent; "
                         f"linear_solver={solver_name!r} needs a symmetric matrix (use 'pardiso_lu')")

    log_file_path = spec.solver_log_file
    if log_file_path is None and spec.output_dir:
        log_file_path = os.path.join(spec.output_dir, f"{spec.name}_solver.log")
    slogger = SolverLogger(
        enable_solver_log=spec.enable_solver_log,
        solver_log_file=log_file_path,
        solver_log_stdout=spec.solver_log_stdout,
        solver_log_level=spec.solver_log_level,
    )
    solver_params_used = dict(spec.amgcl_params) if spec.amgcl_params else None

    with slogger:
        try:
            kratos_ver = str(KM.KratosGlobals.Kernel.Version())
        except Exception:
            kratos_ver = ""
        try:
            from wht_iga.solver.memory import get_total_ram_bytes

            total_ram_gb = get_total_ram_bytes() / (1024.0 ** 3)
        except Exception:
            total_ram_gb = float(decision.total_bytes) / (1024.0 ** 3)
        avail_ram_gb = get_available_ram_bytes() / (1024.0 ** 3)

        slogger.log_stage_1_environment(
            job_name=spec.name,
            total_ram_gb=total_ram_gb,
            avail_ram_gb=avail_ram_gb,
            mkl_threads=spec.num_threads_mkl,
            kratos_version=kratos_ver,
        )
        slogger.log_stage_2_model(
            n_cells=spec.n_cells,
            degree=spec.degree,
            continuity=spec.continuity,
            n_ctrl_pts=box.n_ctrl,
            n_total_dofs=total_dofs,
            n_active_dofs=active_dofs,
            n_in_cells=n_in_count,
            n_cut_cells=n_cut_count,
            n_out_cells=n_out_count,
            out_mode=spec.out_mode,
            n_quad_pts=n_qp_total,
            nnz=est_nnz,
        )
        slogger.log_stage_3_constraints(
            n_dirichlet_bcs=len(spec.bcs),
            n_fixed_faces=len(spec.fixed_faces),
            n_inactive_fixed=3 * n_inactive,
            n_ties=int(tie_info.get("n_ties", 0)) if tie_info else 0,
            n_slave_dofs=int(tie_info.get("n_slave", 0)) if tie_info else 0,
            n_master_dofs=int(tie_info.get("n_master", 0)) if tie_info else 0,
            beta=beta,
            c_pen=DEFAULT_C_PEN,
            theta_thb=float(tie_info["theta"]) if tie_info and "theta" in tie_info else None,
        )
        slogger.log_stage_4_solver_selection(
            decision=decision,
            selected_solver=solver_name,
            solver_params=solver_params_used,
        )
        slogger.log_stage_5_header()

        mp.SetBufferSize(2)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        solver = _make_linear_solver(KM, solver_name, spec.num_threads_mkl, amgcl_params=spec.amgcl_params)
        builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
        crit = KM.ResidualCriteria(1e-9, 1e-12)
        strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, move_mesh)
        strategy.SetEchoLevel(0)
        mp.CloneTimeStep(1.0)
        _t1 = _time.perf_counter()
        log(f"run_linear[{spec.name}]: solve start ({solver_name}, {len(mp.Nodes)} nodes, "
            f"{len(mp.MasterSlaveConstraints)} MPC)")
        if hooks is not None and hasattr(hooks, "before_solve"):
            hooks.before_solve(model, mp, strategy)
        strategy.Solve()
        if hooks is not None and hasattr(hooks, "after_solve"):
            hooks.after_solve(model, mp, strategy)
        _t2 = _time.perf_counter()
        log(f"run_linear[{spec.name}]: solve done in {_t2 - _t1:.1f}s")
        try:
            n_iter = int(mp.ProcessInfo[KM.NL_ITERATION_NUMBER])
        except Exception:
            n_iter = -1  # iteration count unavailable from this Kratos build
        try:
            converged = bool(strategy.IsConverged())
        except Exception:
            converged = False
        if not converged:
            slogger.log_stage_5_row(
                step=1, inc=1, att=1, time=1.0, dt=1.0,
                lin_it=max(n_iter, 0), rel_res=float("nan"), energy_norm=0.0, wall_s=_t2 - _t1,
            )
            slogger.log_stage_6_summary(
                converged=False,
                max_disp=0.0,
                total_force=np.asarray(load_keep["total_force"]),
                setup_s=_t1 - _t0,
                solve_s=_t2 - _t1,
                post_s=0.0,
                warnings=[f"{i.code}: {i.message}" for i in issues],
                error_message=f"solver did not report convergence (iterations={n_iter})",
            )
            issues.append(ValidationIssue(
                "error", "SOLVER_NOT_CONVERGED",
                f"solver did not report convergence (iterations={n_iter})",
                "solver"))
            import warnings
            warnings.warn("Kratos solve did not report convergence; "
                          "run_linear raises instead of writing results", RuntimeWarning)
            err = ValidationError(issues)
            # The assembled system is still valid for matrix-only studies (V-COND
            # measures kappa of K, including deliberately singular unstabilised
            # cases that cannot converge). Results are NOT written; the caller
            # must opt in by catching this error and reading err.strategy.
            err.strategy = strategy
            err.keepalive = (volume, keep_parts, model, mp, elements)
            raise err

        vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
        U = np.zeros(box.n_dofs)
        for node in mp.Nodes:
            I = node.Id - 1
            if 0 <= I < box.n_ctrl:
                for d, var in enumerate(vars3):
                    U[3 * I + d] = node.GetSolutionStepValue(var)

        face_avgs = {}
        load_faces = {tuple(p.face) for p in spec.pressures} | {tuple(t.face) for t in spec.tractions}
        for axis, side in {tuple(b.face) for b in spec.bcs} | {tuple(f) for f in spec.fixed_faces} | load_faces:
            ids = box.boundary_control_points(axis, side)
            face_avgs[(axis, side)] = U.reshape(-1, 3)[ids].mean(axis=0)

        coeffs = U.reshape(-1, 3)
        probe_values = {}
        for p in spec.probes:
            xyz = np.asarray(p["xyz"], dtype=float).reshape(1, 3)
            probe_values[p.get("name", "p")] = box.evaluate(coeffs, xyz)[0]

        stresses = []
        log(f"run_linear[{spec.name}]: stress extraction ({len(elements)} elements)")
        for el in elements:
            if el is None:  # OUT cell under out_mode="drop"
                stresses.append(np.zeros((0, 6)))
                continue
            vals = el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)
            # one buffer conversion per element (was a float() per entry: 62% of
            # the post-solve time at 59k elements; identical values)
            stresses.append(np.array(list(vals), dtype=float).reshape(len(vals), -1))

        _t3 = _time.perf_counter()
        rxn_sum = np.zeros(3, dtype=np.float64)
        for node in mp.Nodes:
            if 1 <= node.Id <= box.n_ctrl:
                rxn = node.GetSolutionStepValue(KM.REACTION)
                rxn_sum += np.array([float(rxn[0]), float(rxn[1]), float(rxn[2])], dtype=np.float64)
        total_f = np.asarray(load_keep["total_force"], dtype=np.float64)
        rxn_balance = float(np.linalg.norm(rxn_sum + total_f)) if not spec.bcs else None
        max_disp = float(np.max(np.abs(U))) if U.size else 0.0
        lin_it_log = int(solver.GetIterationsNumber()) if solver_name == "cg" else max(n_iter, 1)
        energy_est = float(0.5 * abs(float(np.dot(U.reshape(-1, 3).mean(axis=0), total_f)))) if U.size else 0.0
        slogger.log_stage_5_row(
            step=1,
            inc=1,
            att=1,
            time=1.0,
            dt=1.0,
            lin_it=lin_it_log,
            rel_res=0.0 if converged else float("nan"),
            energy_norm=energy_est,
            wall_s=_t2 - _t1,
        )
        slogger.log_stage_6_summary(
            converged=converged,
            max_disp=max_disp,
            total_force=total_f,
            reaction_residual=rxn_balance,
            setup_s=_t1 - _t0,
            solve_s=_t2 - _t1,
            post_s=_t3 - _t2,
            warnings=[f"{i.code}: {i.message}" for i in issues] if issues else None,
        )

    solver_log_text = slogger.get_text()
    resolved_log_file = str(slogger.solver_log_file) if (slogger.enabled and slogger.solver_log_file is not None) else None

    res = LinearResult(U=U, beta_used=beta, total_force=np.asarray(load_keep["total_force"]),
                       face_averages=face_avgs, probe_values=probe_values,
                       cauchy_stress=stresses, issues=issues,
                       keepalive=(volume, keep_parts, model, mp, elements),
                       strategy=strategy, tie_info=tie_info,
                       converged=converged, n_iter=n_iter,
                       timings={"setup_s": _t1 - _t0, "solve_s": _t2 - _t1, "post_s": _t3 - _t2,
                                "linear_solver": solver_name, "active_dofs": active_dofs,
                                "total_dofs": total_dofs},
                       solver_log=solver_log_text,
                       solver_log_file=resolved_log_file,
                       active_dofs=active_dofs,
                       total_dofs=total_dofs,
                       solver_decision={
                           "requested_solver": decision.requested_solver,
                           "selected_solver": decision.selected_solver,
                           "kratos_solver": solver_name,
                           "est_pardiso_bytes": decision.est_pardiso_bytes,
                           "est_amgcl_bytes": decision.est_amgcl_bytes,
                           "avail_bytes": decision.avail_bytes,
                           "switched": decision.switched,
                           "switched_reason": decision.switched_reason,
                       })
    if spec.output_dir:
        os.makedirs(spec.output_dir, exist_ok=True)
        log(f"run_linear[{spec.name}]: writing outputs")
        res.mesh_db, res.results_db, res.vtu_path, res.probe_json = _write_outputs(spec, box, res, beta, cell_records)
    return res


def _write_outputs(spec: JobSpec, box: BoxSpline, res: LinearResult, beta: float,
                   cell_records: list[dict] | None = None) -> tuple[str, str, str, str]:
    from wht_iga.solver import db as meshdb
    from wht_iga.solver.vtu import hex_proxy_grid, write_vtu

    name = spec.name
    d = spec.output_dir or "."
    nx, ny, nz = box.n_cells
    out = {tuple(c) for c in spec.out_cells}
    if cell_records is not None:
        cells = cell_records
    else:
        cells = [{"i": i, "j": j, "k": k, "status": 0 if (i, j, k) in out else 1,
                  "eta": 0.0 if (i, j, k) in out else 1.0,
                  "n_qp": (box.degree + 1) ** 3}
             for k in range(nz) for j in range(ny) for i in range(nx)]
    mesh_path = os.path.join(d, f"{name}.wiga.sqlite")
    meshdb.write_mesh_db(mesh_path, spec.to_json(),
                         meta_extra={"beta_used": repr(beta), "kratos": "E:/kratos-build/install"},
                         degree=box.degree, continuity=box.continuity,
                         lo=spec.lo, hi=spec.hi, n_cells=spec.n_cells,
                         cells=cells,
                         probes=spec.probes,
                         issues=[{"severity": i.severity, "code": i.code, "message": i.message, "entity": i.entity} for i in res.issues])
    results_path = os.path.join(d, f"{name}_iga_results.sqlite")
    s_all = np.concatenate(res.cauchy_stress, axis=0) if res.cauchy_stress else np.zeros((0, 6))
    meshdb.write_results_db(results_path, spec.to_json(),
                            frames=[{"frame_id": 0, "step_id": 0, "time": 1.0,
                                     "increment": int(res.n_iter),
                                     "converged": 1 if res.converged else 0,
                                     "n_iter": int(res.n_iter)}],
                            cp_fields={(0, "U"): res.U},
                            qp_fields={(0, "S"): s_all})
    pts, conn = hex_proxy_grid(np.asarray(spec.lo), np.asarray(spec.hi), spec.n_cells)
    gu, gv, gw = [np.linspace(spec.lo[dd], spec.hi[dd], spec.n_cells[dd] + 1) for dd in range(3)]
    W, V, Ugrid = np.meshgrid(gw, gv, gu, indexing="ij")
    grid_pts = np.column_stack([Ugrid.ravel(), V.ravel(), W.ravel()])
    field = box.evaluate(res.U.reshape(-1, 3), grid_pts)
    vtu_path = os.path.join(d, f"{name}.vtu")
    write_vtu(vtu_path, pts, conn, {"U": field})
    probe_path = os.path.join(d, f"{name}_probe_data.json")
    with open(probe_path, "w", encoding="utf-8") as f:
        json.dump({k: [float(v) for v in np.atleast_1d(val)] for k, val in res.probe_values.items()}, f)
    return mesh_path, results_path, vtu_path, probe_path
