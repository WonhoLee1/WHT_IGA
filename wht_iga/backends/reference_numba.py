"""``ReferenceNumbaBackend``: Kratos-free oracle for untrimmed box grids.

Assembles the exact B-spline stiffness with the Numba reference assembler
(``wht_iga.reference.assemble_elastic``), consistent face loads and body force,
imposes Dirichlet data STRONGLY on the face control points (exact for constant
values on clamped faces) and solves with ``wht_iga.reference.solve``. It is the
cross-verification partner of ``KratosPythonBackend`` (which uses penalty BCs),
so agreement is expected up to the penalty consistency error.

Scope: untrimmed grids (no CUT/OUT cells), no follower pressure.
"""

from __future__ import annotations

import time

import numpy as np

from wht_iga.ports.contracts import BackendCapabilities, BackendSolution, StaticJob


class BackendUnsupported(RuntimeError):
    """The job needs a feature this backend does not implement."""


def _axis_integrals(box, axis: int) -> np.ndarray:
    """∫ N_a(x) dx for every 1D basis function along ``axis``."""
    xg, wg = np.polynomial.legendre.leggauss(box.degree + 1)
    first, N, _dN, _X, jac = box.axis_tables(axis, xg)
    out = np.zeros(box.n_cp[axis])
    for c in range(len(jac)):
        for a in range(box.degree + 1):
            out[first[c] + a] += jac[c] * float(np.dot(wg, N[c, :, a]))
    return out


def body_force_vector(box, accel, density: float) -> np.ndarray:
    """Consistent body-force vector ``f_I = density * accel * ∫ N_I dV``."""
    ix, iy, iz = (_axis_integrals(box, d) for d in range(3))
    vol_int = np.einsum("k,j,i->kji", iz, iy, ix).ravel()  # I = i + nu*(j + nv*k)
    f = np.zeros(box.n_dofs)
    for d in range(3):
        f[d::3] = density * float(accel[d]) * vol_int
    return f


class ReferenceNumbaBackend:
    """Oracle backend (strong Dirichlet, exact Gauss, sparse direct solve)."""

    name = "reference_numba"

    def capabilities(self) -> BackendCapabilities:
        """Untrimmed boxes with BCs, pressure/traction and gravity."""
        return BackendCapabilities(name=self.name, trimmed=False, partial_component_bcs=True,
                                   follower_pressure=False, gravity=True, writes_solver_artifacts=False)

    def solve_static(self, job: StaticJob) -> BackendSolution:
        """Assemble and solve one static step on the untrimmed grid."""
        from wht_iga.reference.assemble_elastic import assemble_face_traction, assemble_stiffness
        from wht_iga.reference.solve import solve_with_dirichlet
        from wht_iga.solver.memory import get_available_ram_bytes, get_total_ram_bytes
        from wht_iga.solver.solver_log import SolverLogger

        js = job.job_spec
        mesh = job.mesh
        if mesh.trimmed or mesh.out_cells:
            raise BackendUnsupported("ReferenceNumbaBackend handles untrimmed grids only "
                                     "(use KratosPythonBackend for immersed parts)")
        if any(p.follower for p in js.pressures):
            raise BackendUnsupported("ReferenceNumbaBackend has no follower-pressure (geometric) stiffness")
        t0 = time.perf_counter()
        box = mesh.box()
        f_surf = np.zeros(box.n_dofs)
        for t in js.tractions:
            f_surf += assemble_face_traction(box, t.face[0], t.face[1], np.asarray(t.vector, float))
        for p in js.pressures:
            n_out = np.zeros(3)
            n_out[p.face[0]] = 1.0 if p.face[1] == 1 else -1.0
            f_surf += assemble_face_traction(box, p.face[0], p.face[1], p.magnitude * n_out)
        f = f_surf.copy()
        if js.gravity is not None and js.density > 0.0:
            f += body_force_vector(box, js.gravity, js.density)
        fixed: dict[int, float] = {}
        for bc in js.bcs:
            for I in box.boundary_control_points(*bc.face):
                for d in range(3):
                    if bc.components[d]:
                        fixed[3 * int(I) + d] = float(bc.value[d])
        for axis, side in js.fixed_faces:
            for I in box.boundary_control_points(axis, side):
                for d in range(3):
                    fixed[3 * int(I) + d] = 0.0
        if not fixed:
            raise BackendUnsupported("no Dirichlet data: rigid-body singular")
        idx = np.fromiter(fixed.keys(), dtype=np.int64)
        vals = np.fromiter(fixed.values(), dtype=float)

        mitigation = str(getattr(js, "locking_mitigation", "none")).lower()
        use_fbar = mitigation == "fbar" or (mitigation == "auto" and float(js.nu) >= 0.49)
        if use_fbar:
            from wht_iga.analyses import fbar

            fbar_model = fbar.build_model(box, box.degree, E=js.E, nu=js.nu)
            zero_idx = np.array([k for k, v in fixed.items() if v == 0.0], dtype=np.int64)
            prescribed = {k: v for k, v in fixed.items() if v != 0.0}
            t1 = time.perf_counter()
            U = fbar.solve_static(
                fbar_model,
                f,
                zero_idx,
                prescribed=prescribed or None,
                linear_solver="auto",
            )
            t2 = time.perf_counter()
            R_int, K_tan = fbar.residual_tangent_sparse(fbar_model, U)
            R = R_int - f
            nnz_val = int(K_tan.nnz)
            solver_label = "fbar.solve_static"
        else:
            K = assemble_stiffness(box, js.E, js.nu)
            t1 = time.perf_counter()
            U = solve_with_dirichlet(K.tocsr(), f, idx, vals)
            t2 = time.perf_counter()
            R = K @ U - f  # reactions live on the fixed DOFs; free residual ~ 0
            nnz_val = int(K.nnz)
            solver_label = "reference.solve_with_dirichlet"
        free = np.setdiff1d(np.arange(box.n_dofs), idx)
        rel_res = float(np.linalg.norm(R[free]) / max(np.linalg.norm(f[free]), 1e-300)) if free.size else 0.0
        rxn = np.array([R[idx[idx % 3 == d]].sum() for d in range(3)])
        f_tot = np.array([f[d::3].sum() for d in range(3)])
        rxn_res = float(np.linalg.norm(rxn + f_tot))
        total_force = np.array([f_surf[d::3].sum() for d in range(3)])

        logger = SolverLogger(enable_solver_log=js.enable_solver_log, solver_log_file=js.solver_log_file,
                              solver_log_stdout=js.solver_log_stdout, solver_log_level=js.solver_log_level)
        with logger:
            logger.log_stage_1_environment(job_name=js.name, total_ram_gb=get_total_ram_bytes() / 1024.0 ** 3,
                                           avail_ram_gb=get_available_ram_bytes() / 1024.0 ** 3,
                                           kratos_version="(reference backend, Kratos not used)")
            nx, ny, nz = box.n_cells
            logger.log_stage_2_model(n_cells=box.n_cells, degree=box.degree, continuity=box.continuity,
                                     n_ctrl_pts=box.n_ctrl, n_total_dofs=box.n_dofs,
                                     n_active_dofs=box.n_dofs - len(idx), n_in_cells=nx * ny * nz,
                                     n_cut_cells=0, n_out_cells=0, out_mode="none",
                                     n_quad_pts=nx * ny * nz * (box.degree + 1) ** 3, nnz=nnz_val)
            logger.log_stage_3_constraints(n_dirichlet_bcs=len(js.bcs), n_fixed_faces=len(js.fixed_faces),
                                           extra_notes=[f"strong Dirichlet on {len(idx)} face DOFs"])
            logger.log_stage_4_solver_selection(requested_solver="reference",
                                                selected_solver=solver_label)
            logger.log_stage_5_row(lin_it=1, rel_res=rel_res, wall_s=t2 - t1)
            logger.log_stage_6_summary(converged=True, max_disp=float(np.max(np.abs(U))),
                                       total_force=total_force, reaction_residual=rxn_res,
                                       setup_s=t1 - t0, solve_s=t2 - t1, post_s=time.perf_counter() - t2)
        log_file = str(logger.solver_log_file) if (logger.enabled and logger.solver_log_file) else None
        return BackendSolution(
            backend=self.name, U=U, converged=True, n_iter=1, total_force=total_force,
            timings={"setup_s": t1 - t0, "solve_s": t2 - t1, "linear_solver": "reference",
                     "rel_residual": rel_res},
            solver_log=logger.get_text(), solver_log_file=log_file,
            artifacts={"solver_log": log_file} if log_file else {}, reaction_residual=rxn_res)

    def close(self) -> None:
        """Nothing native to release."""
