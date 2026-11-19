# -*- coding: utf-8 -*-
"""Geometrically nonlinear (NLGEOM) bulk around the surface-Gauss SSN contact solver.

Why (S41): V08 R1 with a small-strain bulk converges to about +5 % over Hertz at
h=0.125 once the penalty and the surface soup are converged; Coreform runs the
case with ``nlgeom=YES``. With a Total-Lagrangian Neo-Hookean bulk the same
model gives +0.58 % (h=0.125, soup n=48, eps_n=1e11).

Scheme (no change to :mod:`wht_iga.analyses.contact`): an outer Newton on the
bulk wraps :func:`solve_deformable_ssn_contact`. At outer iterate ``u_k`` the
Kratos TL elements give ``K_T(u_k)`` and ``f_int(u_k)`` (mesh moved to
``X0 + u_k``); the inner contact solve gets ``Kf = K_T_ff`` and
``f_ext' = K_T_ff u_k - f_int_free(u_k)``, i.e. it solves
``f_int(u_k) + K_T (u - u_k) + f_c(u) = 0`` with the contact fully nonlinear.
The outer loop stops when ``|u_{k+1} - u_k| <= outer_tol |u|``. Its fixed point
is the exact NLGEOM equilibrium ``f_int(u) + f_c(u) = 0``.

Material: use a hyperelastic law (``CL.HyperElastic3DLaw``, compressible
Neo-Hookean). ``SM.LinearElastic3DLaw`` under TL is St. Venant-Kirchhoff, which
softens strongly in compression and gave -8.2 % on V08 (measured 2026-11-19);
Abaqus ``*Elastic`` with ``nlgeom`` is hypoelastic and does not soften that way.
Ghost (OUT) cells use ``ghost_element`` (default ``SmallDisplacementElement3D4N``)
with the linear law at ``ghost_alpha * E``. Inactive control points are fixed at
zero, so a ghost cell between them and active points that slide far is strongly
stretched: as TL + linear law (= St. Venant-Kirchhoff) its force grows like
(u/h)**3 and leaked 1,891 lbf into the zero DOFs on V08 R2 h=0.25 at t=0.70
(cut-face tangential equilibrium error 1.3e-2). The small-strain ghost element
leaks 5.1 lbf (eq error 3.5e-5) with the contact history unchanged to 0.013 %
(measured 2026-11-19, scratchpad v08/ghost_r2.py). Neo-Hookean is never used on
ghost cells (NaN on inverted cells).

Friction: the slip history (``slip_state``) is the state at the START of the
load step and is passed unchanged to every outer iteration; the step-end state
returned by the last inner solve becomes the next step's start state.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from wht_iga.analyses.contact import solve_deformable_ssn_contact
from wht_iga.kratos_backend.elements import check_qp_element_name, equation_permutation
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume


class KratosBodyOperator:
    """Internal force and tangent of one spline body from Kratos elements.

    ``eval(u_full) -> (K_full csr, f_int_full)`` in reference DOF order
    ``3*I + d`` (I = control point). DOF classes: ``zero`` (homogeneous
    Dirichlet incl. inactive control points), ``presc`` (prescribed), ``free``.
    """

    def __init__(
        self,
        KM: Any,
        SM: Any,
        box: Any,
        E: float,
        nu: float,
        *,
        law_factory: Callable[[], Any] | None = None,
        element: str = "TotalLagrangianElement3D4N",
        zero_fn: Callable[[np.ndarray], tuple[int, ...]] | None = None,
        presc_fn: Callable[[np.ndarray], tuple[int, ...]] | None = None,
        trim_soup: tuple[np.ndarray, np.ndarray] | None = None,
        q: int = 4,
        ghost_alpha: float = 1e-8,
        ghost_element: str = "SmallDisplacementElement3D4N",
    ) -> None:
        from wht_iga.solver.active import inactive_ctrl_mask
        from wht_iga.solver.trimmed import build_trimmed_table

        check_qp_element_name(element)
        check_qp_element_name(ghost_element)
        self.KM = KM
        self.box = box
        if trim_soup is not None:
            qp_table, out_cells, _issues, _rec, self.diag = build_trimmed_table(
                trim_soup[0], trim_soup[1], box, q=int(q), cache_dir="off")
            inact = inactive_ctrl_mask(box, out_cells)
        else:
            qp_table, out_cells = {}, []
            self.diag = {"n_cut": 0, "n_fallback": 0}
            inact = np.zeros(box.n_ctrl, dtype=bool)
        out_set = {tuple(c) for c in out_cells}

        self.model = KM.Model()
        mp = self.model.CreateModelPart("NLGEOMBody")
        self.mp = mp
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        vol, _ = make_nurbs_volume(KM, mp, box)
        law_factory = law_factory or SM.LinearElastic3DLaw
        pv = mp.CreateNewProperties(1)
        pv.SetValue(KM.YOUNG_MODULUS, float(E))
        pv.SetValue(KM.POISSON_RATIO, float(nu))
        pv.SetValue(KM.CONSTITUTIVE_LAW, law_factory())
        pg = mp.CreateNewProperties(2)
        pg.SetValue(KM.YOUNG_MODULUS, float(E) * float(ghost_alpha))
        pg.SetValue(KM.POISSON_RATIO, float(nu))
        pg.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
        nx, ny, nz = (int(v) for v in box.n_cells)
        eid = 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    elem_name = element
                    if (i, j, k) in out_set:
                        if ghost_alpha <= 0.0:
                            continue
                        pts, props, elem_name = cell_gauss_points(box, i, j, k), pg, ghost_element
                    elif (i, j, k) in qp_table:
                        pts, props = np.asarray(qp_table[(i, j, k)], dtype=float).tolist(), pv
                    else:
                        pts, props = cell_gauss_points(box, i, j, k), pv
                    qps = KM.GeometriesVector()
                    vol.CreateQuadraturePointGeometries(qps, 2, pts)
                    mp.CreateNewElement(elem_name, eid, qps[0], props)
                    eid += 1
        for v, r in ((KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                     (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
            KM.VariableUtils().AddDof(v, r, mp)
        self._scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        self._builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
        self._strat = KM.ResidualBasedLinearStrategy(mp, self._scheme, self._builder, False, False, False, False)
        self._strat.SetEchoLevel(0)
        self._strat.Initialize()
        mp.CloneTimeStep(1.0)
        self._strat.InitializeSolutionStep()
        self._A = self._strat.GetSystemMatrix()
        self._b = self._strat.GetSystemVector()
        self._space = KM.UblasSparseSpace()
        self._inv = np.argsort(equation_permutation(KM, mp, box))
        self._ctrl_nodes = sorted((n for n in mp.Nodes if 0 <= n.Id - 1 < box.n_ctrl), key=lambda n: n.Id)
        self._X0 = np.array([[n.X0, n.Y0, n.Z0] for n in self._ctrl_nodes])

        Xcp = box.control_points()
        zero = np.zeros(box.n_dofs, dtype=bool)
        presc = np.zeros(box.n_dofs, dtype=bool)
        for I in np.flatnonzero(inact):
            zero[3 * I:3 * I + 3] = True
        if zero_fn is not None:
            for I in range(box.n_ctrl):
                for d in zero_fn(Xcp[I]):
                    zero[3 * I + d] = True
        if presc_fn is not None:
            for I in range(box.n_ctrl):
                if inact[I]:
                    continue
                for d in presc_fn(Xcp[I]):
                    if not zero[3 * I + d]:
                        presc[3 * I + d] = True
        self.zero_mask, self.presc_mask = zero, presc
        self.free_mask = ~(zero | presc)

    def eval(self, u_full: np.ndarray):
        """Return ``(K_T csr, f_int)`` at ``u_full`` (moves the Kratos mesh to ``X0 + u``)."""
        import KratosMultiphysics.scipy_conversion_tools as sct

        KM = self.KM
        U = np.asarray(u_full, dtype=float).reshape(-1, 3)
        for n, x0, uu in zip(self._ctrl_nodes, self._X0, U):
            n.SetSolutionStepValue(KM.DISPLACEMENT, KM.Array3([float(uu[0]), float(uu[1]), float(uu[2])]))
            n.X, n.Y, n.Z = x0[0] + uu[0], x0[1] + uu[1], x0[2] + uu[2]
        self._space.SetToZeroMatrix(self._A)
        self._space.SetToZeroVector(self._b)
        self._builder.Build(self._scheme, self.mp, self._A, self._b)
        K = sct.to_csr(self._A)[self._inv, :][:, self._inv].tocsr()
        f_int = -np.array(self._b, copy=True)[self._inv]
        return K, f_int


def solve_nlgeom_contact_step(
    body_m: KratosBodyOperator,
    body_s: KratosBodyOperator,
    X_s: np.ndarray,
    w_q: np.ndarray,
    B_s,
    master_face,
    eps_n: float,
    um_full: np.ndarray,
    us_full: np.ndarray,
    *,
    mu: float = 0.0,
    eps_t: float | None = None,
    slip_state: dict | None = None,
    outer_tol: float = 1e-8,
    max_outer: int = 30,
    rtol: float = 1e-9,
    use_pardiso: bool = True,
) -> dict:
    """One load step of NLGEOM bulk + SSN contact (outer Newton, see module doc).

    ``um_full`` / ``us_full`` are full DOF vectors; prescribed entries of
    ``us_full`` must already hold this step's values. Returns the last inner
    result dict extended with ``um_full``, ``us_full``, ``n_outer``, ``outer_rel``
    and ``f_int_s`` (slave internal force at the converged state, for cut-face
    reactions ``R = f_int_s[presc] + f_c_s_full[presc]``).
    """
    um_full = np.asarray(um_full, dtype=float).copy()
    us_full = np.asarray(us_full, dtype=float).copy()
    fm_mask, fs_mask = body_m.free_mask, body_s.free_mask
    rel = np.inf
    res = None
    for n_outer in range(1, int(max_outer) + 1):
        Km, fm = body_m.eval(um_full)
        Ks, fs = body_s.eval(us_full)
        Kmf = Km[fm_mask][:, fm_mask].tocsr()
        Ksf = Ks[fs_mask][:, fs_mask].tocsr()
        um, us = um_full[fm_mask].copy(), us_full[fs_mask].copy()
        res = solve_deformable_ssn_contact(
            Kmf, Ksf, X_s, w_q, B_s, master_face, eps_n,
            free_m=fm_mask, free_s=fs_mask,
            u_fixed_m=um_full.copy(), u_fixed_s=us_full.copy(),
            f_ext_m=Kmf @ um - fm[fm_mask], f_ext_s=Ksf @ us - fs[fs_mask],
            u_init_m=um, u_init_s=us, mu=mu, eps_t=eps_t, slip_state=slip_state,
            rtol=rtol, use_pardiso=use_pardiso,
        )
        du = np.concatenate([res["u_free_m"] - um, res["u_free_s"] - us])
        um_full[fm_mask] = res["u_free_m"]
        us_full[fs_mask] = res["u_free_s"]
        rel = float(np.linalg.norm(du) / max(np.linalg.norm(np.concatenate([um_full, us_full])), 1e-300))
        if rel <= outer_tol:
            break
    else:
        raise RuntimeError(f"NLGEOM outer Newton did not converge in {max_outer} iterations (rel={rel:.2e})")
    _, fs = body_s.eval(us_full)
    out = dict(res)
    out.update(um_full=um_full, us_full=us_full, n_outer=n_outer, outer_rel=rel, f_int_s=fs)
    return out
