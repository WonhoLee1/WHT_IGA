"""Kratos stiffness + Python consistent mass + eigenmodes for box B-spline models (P4-A).

Kratos' ``QuadraturePointGeometry`` elements return an identically zero consistent
mass matrix (spike S06), so ``EigensolverStrategy`` cannot be used.  Only the lumped
mass vector is valid.  The route here is therefore:

* ``K``  from Kratos (route R1: one QP-geometry element per cell, block builder);
* ``M``  assembled in Python from the same quadrature (:func:`assemble_box_mass`);
* modes from ``K phi = w**2 M phi`` (``wht_iga.dynamics.solve_eigenmodes``).

DOF ordering is the reference ordering of :class:`~wht_iga.reference.box.BoxSpline`
(``3 * I + d``, control point ``I = i + nu * (j + nv * k)``); the Kratos matrix is
permuted into it with :func:`~wht_iga.kratos_backend.elements.equation_permutation`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp

from wht_iga.dynamics.consistent_mass import (
    EigenmodeResult,
    assemble_consistent_mass,
    assemble_consistent_mass_from_bspline,
    solve_eigenmodes,
)
from wht_iga.grid.bspline import basis_1d, ders_basis_funs, find_span
from wht_iga.reference.assemble_elastic import isotropic_D
from wht_iga.reference.box import BoxSpline


@dataclass(frozen=True)
class KratosModalModel:
    """Stiffness, mass and eigenmodes of one clamped box B-spline model.

    Attributes:
        K: Stiffness in reference DOF order (Kratos-assembled), ``(n_dofs, n_dofs)``.
            Rows/columns of clamped DOFs are zeroed by the block builder.
        M: Consistent mass in reference DOF order, ``(n_dofs, n_dofs)``.
        free_dofs: Reference indices of the unconstrained DOFs.
        eig: Eigenmodes (:class:`~wht_iga.dynamics.EigenmodeResult`).
    """

    K: sp.csr_matrix
    M: sp.csr_matrix
    free_dofs: np.ndarray
    eig: EigenmodeResult


def face_dofs(box: BoxSpline, axis: int, side: int, components: Sequence[int] = (0, 1, 2)) -> np.ndarray:
    """Reference DOF indices ``3 I + d`` of the control points on one box face.

    Args:
        box: Box B-spline.
        axis: Face normal axis ``0, 1, 2``.
        side: ``0`` for the low face, ``1`` for the high face.
        components: Displacement components to include.

    Returns:
        Sorted DOF indices.
    """
    cps = np.asarray(box.boundary_control_points(axis, side), dtype=np.int64)
    return np.sort(np.concatenate([3 * cps + int(d) for d in components]))


def _build_kratos_box(box: BoxSpline, E: float, nu: float, density: float, dynamic: bool = False):
    """Create a Kratos model part with one QP element per cell and DISPLACEMENT DOFs.

    Args:
        box: Box B-spline.
        E: Young's modulus.
        nu: Poisson's ratio.
        density: Mass density.
        dynamic: Also register VELOCITY, ACCELERATION and POINT_LOAD and use a buffer of
            two steps.  Solution-step variables must exist *before* the nodes are
            created, which is why this is a constructor flag and not a later call.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.elements import attach_cell_elements
    from wht_iga.kratos_backend.geometry import make_nurbs_volume

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("P4Box")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    if dynamic:
        mp.AddNodalSolutionStepVariable(KM.VELOCITY)
        mp.AddNodalSolutionStepVariable(KM.ACCELERATION)
        mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
        mp.SetBufferSize(2)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    volume, _nodes = make_nurbs_volume(KM, mp, box)
    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.DENSITY, float(density))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    elements = attach_cell_elements(KM, mp, volume, box, props)
    utils = KM.VariableUtils()
    utils.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    utils.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    utils.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
    return KM, mp, volume, elements


def kratos_box_stiffness(
    box: BoxSpline,
    E: float,
    nu: float,
    fixed_dofs: Sequence[int] | np.ndarray = (),
    density: float = 1.0,
) -> sp.csr_matrix:
    """Stiffness of a box B-spline assembled by Kratos, returned in reference DOF order.

    The block builder keeps clamped DOFs in the matrix with zeroed rows and columns, so
    only the free-free block is meaningful; pass ``free_dofs = complement(fixed_dofs)``
    to :func:`~wht_iga.dynamics.solve_eigenmodes`.

    Args:
        box: Box B-spline (untrimmed; every cell is a full Gauss cell).
        E: Young's modulus.
        nu: Poisson's ratio.
        fixed_dofs: Reference DOF indices ``3 I + d`` to clamp (needed so Kratos' linear
            solve of the zero-load problem is non-singular).
        density: Mass density stored on the properties (unused by K).

    Returns:
        CSR matrix of shape ``(box.n_dofs, box.n_dofs)``.
    """
    from wht_iga.kratos_backend.elements import equation_permutation, solve_linear, system_matrix_csr

    KM, mp, _volume, _elements = _build_kratos_box(box, E, nu, density)
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for dof in np.asarray(fixed_dofs, dtype=np.int64).reshape(-1):
        mp.GetNode(int(dof) // 3 + 1).Fix(comp_vars[int(dof) % 3])
    strategy = solve_linear(KM, mp)
    perm = equation_permutation(KM, mp, box)
    k_kr = system_matrix_csr(strategy)
    return sp.csr_matrix(k_kr[perm][:, perm])


def kratos_total_mass(box: BoxSpline, density: float, E: float = 1.0, nu: float = 0.0) -> float:
    """Total mass from Kratos' (valid) lumped mass vectors: ``sum(lumped) / 3``.

    Args:
        box: Box B-spline.
        density: Mass density.
        E: Young's modulus (irrelevant to mass).
        nu: Poisson's ratio (irrelevant to mass).

    Returns:
        ``density * volume`` as Kratos computes it.
    """
    KM, mp, _volume, elements = _build_kratos_box(box, E, nu, density)
    total = 0.0
    for el in elements:
        el.Initialize(mp.ProcessInfo)
        lump = KM.Vector(0)
        el.CalculateLumpedMassVector(lump, mp.ProcessInfo)
        total += float(np.asarray(lump).sum())
    return total / 3.0


def _cell_basis(box: BoxSpline, cell: tuple[int, int, int], pts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Scalar basis values at ``pts`` and control-point indices of one cell.

    Local ordering is ``a + (p+1) * (b + (p+1) * c)`` for the control point
    ``(i0 + a, j0 + b, k0 + c)`` -- the same order used for the returned indices.

    Args:
        box: Box B-spline.
        cell: Cell indices ``(i, j, k)`` (selects the knot span of every point).
        pts: Points inside (or on the boundary of) the cell, shape ``(n, 3)``.

    Returns:
        ``(N, cps)`` with ``N`` of shape ``(n, (p+1)**3)`` and ``cps`` of shape ``((p+1)**3,)``.
    """
    p = box.degree
    nu, nv, nw = box.n_cp
    br = [box.breakpoints(d) for d in range(3)]
    n_cps = (nu, nv, nw)
    spans = []
    for d in range(3):
        mid = 0.5 * (br[d][cell[d]] + br[d][cell[d] + 1])
        spans.append(int(find_span(n_cps[d], p, mid, np.ascontiguousarray(box.knots[d], dtype=float))))
    vals = np.empty((pts.shape[0], 3, p + 1))
    for n, x in enumerate(pts):
        for d in range(3):
            vals[n, d] = ders_basis_funs(spans[d], float(x[d]), p, 0, np.ascontiguousarray(box.knots[d], dtype=float))[0]
    n_loc = (p + 1) ** 3
    basis = np.einsum("na,nb,nc->ncba", vals[:, 0], vals[:, 1], vals[:, 2]).reshape(pts.shape[0], n_loc)
    idx = np.arange(p + 1)
    cps = (spans[0] - p + idx)[None, None, :] + nu * ((spans[1] - p + idx)[None, :, None] + nv * (spans[2] - p + idx)[:, None, None])
    return basis, cps.reshape(-1).astype(np.int64)


def _cell_gauss(box: BoxSpline, cell: tuple[int, int, int], n_gauss: int) -> np.ndarray:
    """Tensor Gauss points and weights ``(n, 4)`` of one cell (weight includes the cell Jacobian)."""
    xg, wg = np.polynomial.legendre.leggauss(n_gauss)
    br = [box.breakpoints(d) for d in range(3)]
    pts, wts = [], []
    for d in range(3):
        a, b = br[d][cell[d]], br[d][cell[d] + 1]
        pts.append(a + 0.5 * (b - a) * (xg + 1.0))
        wts.append(0.5 * (b - a) * wg)
    X, Y, Z = np.meshgrid(*pts, indexing="ij")
    W = np.einsum("i,j,k->ijk", *wts)
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel(), W.ravel()])


def assemble_box_mass(
    box: BoxSpline,
    rho: float,
    qp_table: Mapping[tuple[int, int, int], np.ndarray] | None = None,
    out_cells: Sequence[tuple[int, int, int]] | None = None,
    ghost_alpha: float | None = None,
    n_gauss: int | None = None,
) -> sp.csr_matrix:
    """Consistent mass ``M = int rho N^T N dV`` of a (possibly trimmed) box B-spline.

    * Cells in ``qp_table`` (cut cells) use their moment-fitted points ``(x, y, z, w)``.
    * Cells in ``out_cells`` are skipped (``ghost_alpha=None``, "drop" mode) or carry the
      ghost density ``ghost_alpha * rho`` on a full Gauss rule (FCM fill).
    * All other cells are full cells with ``n_gauss`` Gauss points per axis
      (default ``p + 1``, exact for the mass integrand).

    Without ``qp_table``/``out_cells`` this reduces to the verified
    :func:`~wht_iga.dynamics.assemble_consistent_mass_from_bspline`.

    Args:
        box: Box B-spline.
        rho: Mass density of the material.
        qp_table: ``{(i, j, k): array (nq, 4)}`` of cut-cell quadrature points and weights.
        out_cells: Cells entirely outside the body.
        ghost_alpha: Ghost density factor for OUT cells, or ``None`` to drop them.
        n_gauss: Gauss points per axis for full cells.

    Returns:
        CSR matrix ``(box.n_dofs, box.n_dofs)``.
    """
    if not qp_table and not out_cells:
        return assemble_consistent_mass_from_bspline(box.knots, degrees=box.degree, rho=rho, n_gp=n_gauss)
    ng = box.degree + 1 if n_gauss is None else int(n_gauss)
    cut = {tuple(int(v) for v in k): np.asarray(t, dtype=float) for k, t in (qp_table or {}).items()}
    outs = {tuple(int(v) for v in c) for c in (out_cells or [])}
    nx, ny, nz = box.n_cells
    cell_pts, cell_w, cell_n, cell_d, cell_rho = [], [], [], [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = (i, j, k)
                if cell in cut:
                    table, scale = cut[cell], 1.0
                elif cell in outs:
                    if ghost_alpha is None:
                        continue
                    table, scale = _cell_gauss(box, cell, ng), float(ghost_alpha)
                else:
                    table, scale = _cell_gauss(box, cell, ng), 1.0
                if table.shape[0] == 0:
                    continue
                basis, cps = _cell_basis(box, cell, table[:, :3])
                cell_pts.append(table[:, :3])
                cell_w.append(table[:, 3] * scale)
                cell_n.append(basis)
                cell_d.append(cps)
    return assemble_consistent_mass(cell_pts, cell_w, cell_n, cell_d, rho, n_dofs=box.n_dofs)


def reduce_by_operator(A: sp.spmatrix | np.ndarray, T: sp.spmatrix | np.ndarray) -> sp.csr_matrix:
    """Algebraic reduction ``T^T A T`` (extension / THB operator), symmetrised.

    Args:
        A: Square matrix in the fine (full) DOF space.
        T: Operator ``(n_fine, n_reduced)`` mapping reduced DOFs to fine DOFs.

    Returns:
        ``(n_reduced, n_reduced)`` CSR matrix.
    """
    Tc = sp.csr_matrix(T)
    out = sp.csr_matrix(Tc.T @ sp.csr_matrix(A) @ Tc)
    return sp.csr_matrix(0.5 * (out + out.T))


def solve_kratos_eigenmodes(
    box: BoxSpline,
    E: float,
    nu: float,
    density: float,
    clamped_dofs: Sequence[int] | np.ndarray,
    n_modes: int = 20,
    sigma: float = 0.0,
    mass: sp.spmatrix | None = None,
) -> KratosModalModel:
    """Eigenmodes of a clamped box B-spline: Kratos ``K`` + Python consistent ``M``.

    Args:
        box: Box B-spline.
        E: Young's modulus.
        nu: Poisson's ratio.
        density: Mass density.
        clamped_dofs: Reference DOF indices fixed to zero.
        n_modes: Number of modes.
        sigma: Shift-invert shift (rad/s)**2.
        mass: Optional precomputed consistent mass (e.g. a trimmed one).

    Returns:
        :class:`KratosModalModel`.
    """
    clamped = np.unique(np.asarray(clamped_dofs, dtype=np.int64))
    K = kratos_box_stiffness(box, E, nu, clamped, density)
    M = sp.csr_matrix(mass) if mass is not None else assemble_box_mass(box, density)
    free = np.setdiff1d(np.arange(box.n_dofs, dtype=np.int64), clamped)
    eig = solve_eigenmodes(K, M, n_modes=n_modes, sigma=sigma, free_dofs=free)
    return KratosModalModel(K=K, M=M, free_dofs=free, eig=eig)


def kratos_recover_modal_stresses(
    box: BoxSpline,
    E: float,
    nu: float,
    modes: np.ndarray,
    points: np.ndarray,
    *,
    verify_zero_contribution: bool = True,
) -> np.ndarray:
    """Stress mode shapes ``psi_k`` via Kratos zero-weight virtual probe elements.

    Creates one zero-weight (``w = 0.0``) ``SmallDisplacementElement3D4N`` probe
    element per evaluation point on the box ``NurbsVolumeGeometry`` so each point
    resolves its exact knot span and adds zero stiffness/load to the system.
    Modal Cauchy stresses are read back through
    ``el.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, ProcessInfo)``.

    Args:
        box: Box B-spline.
        E: Young's modulus.
        nu: Poisson's ratio.
        modes: Mode shapes in reference DOF order, shape ``(n_dofs, n_modes)``.
        points: Physical points, shape ``(n_pts, 3)``, inside the box.
        verify_zero_contribution: If ``True``, check ``CalculateLocalSystem`` returns
            exact zeros for every probe element.

    Returns:
        Array ``(n_pts, n_modes, 6)`` in Kratos Voigt order ``[xx, yy, zz, xy, yz, xz]``.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.geometry import make_nurbs_volume

    modes_arr = np.asarray(modes, dtype=float)
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    n_modes = modes_arr.shape[1]
    if pts.shape[0] == 0 or n_modes == 0:
        return np.zeros((pts.shape[0], n_modes, 6), dtype=float)

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("ModalStressProbe")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    volume, _nodes = make_nurbs_volume(KM, mp, box, add_dofs=False)
    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    probe_elements = []
    for idx, x in enumerate(pts):
        clamped = [min(max(float(x[d]), float(box.lo[d])), float(box.hi[d])) for d in range(3)]
        qps = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(qps, 2, [[clamped[0], clamped[1], clamped[2], 0.0]])
        if len(qps) != 1:
            raise RuntimeError(f"Expected 1 probe QuadraturePointGeometry, got {len(qps)}")
        el = mp.CreateNewElement("SmallDisplacementElement3D4N", 900001 + idx, qps[0], props)
        el.Initialize(mp.ProcessInfo)
        if verify_zero_contribution:
            lhs = KM.Matrix(0, 0)
            rhs = KM.Vector(0)
            el.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            lhs_max = float(np.max(np.abs(np.asarray(lhs)))) if lhs.Size1() and lhs.Size2() else 0.0
            rhs_max = float(np.max(np.abs(np.asarray(rhs)))) if len(rhs) else 0.0
            if lhs_max != 0.0 or rhs_max != 0.0:
                raise AssertionError(f"Zero-weight probe element produced non-zero LHS/RHS ({lhs_max}, {rhs_max})")
        probe_elements.append(el)

    nodes_sorted = sorted(mp.Nodes, key=lambda n: n.Id)
    out = np.zeros((pts.shape[0], n_modes, 6), dtype=float)
    for k in range(n_modes):
        phi = modes_arr[:, k].reshape(box.n_ctrl, 3)
        for nd in nodes_sorted:
            nd.SetSolutionStepValue(KM.DISPLACEMENT, 0, phi[nd.Id - 1].tolist())
        for idx, el in enumerate(probe_elements):
            sig = np.asarray(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo), dtype=float)
            out[idx, k] = sig[0]
    return out


def recover_modal_stresses(
    box: BoxSpline,
    E: float,
    nu: float,
    modes: np.ndarray,
    points: np.ndarray,
    *,
    backend: str = "analytic",
    verify_zero_contribution: bool = True,
) -> np.ndarray:
    """Stress mode shapes ``psi_k = D B(x) phi_k`` at physical points.

    With ``backend="analytic"`` (default), evaluates the symmetric gradient of the
    B-spline displacement field directly (right-limit values on knot lines) and
    multiplies by ``D`` in Kratos Voigt order ``[xx, yy, zz, xy, yz, xz]``.
    With ``backend="kratos"``, delegates to :func:`kratos_recover_modal_stresses`
    using zero-weight Kratos probe elements.

    Args:
        box: Box B-spline.
        E: Young's modulus.
        nu: Poisson's ratio.
        modes: Mode shapes in reference DOF order, shape ``(n_dofs, n_modes)``.
        points: Physical points, shape ``(n_pts, 3)``, inside the box.
        backend: ``"analytic"`` (default) or ``"kratos"``.
        verify_zero_contribution: Passed to :func:`kratos_recover_modal_stresses`
            when ``backend="kratos"``.

    Returns:
        Array ``(n_pts, n_modes, 6)``.
    """
    if backend == "kratos":
        return kratos_recover_modal_stresses(
            box, E, nu, modes, points, verify_zero_contribution=verify_zero_contribution
        )
    if backend != "analytic":
        raise ValueError(f"unknown modal stress recovery backend {backend!r}")
    p = box.degree
    modes = np.asarray(modes, dtype=float)
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    n_modes = modes.shape[1]
    coef = modes.reshape(box.n_ctrl, 3, n_modes)
    D = isotropic_D(E, nu)
    nu_cp, nv_cp, _nw = box.n_cp
    out = np.zeros((pts.shape[0], n_modes, 6))
    for n, x in enumerate(pts):
        spans, ders = [], []
        for d in range(3):
            xd = min(max(float(x[d]), float(box.lo[d])), float(box.hi[d]))
            span, dd = basis_1d(box.knots[d], p, xd, 1)
            spans.append(int(span) - p)
            ders.append(np.asarray(dd))  # (2, p+1): values, first derivatives
        # per-axis (value, derivative) tensors over the local (c, b, a) grid
        val = [ders[d][0] for d in range(3)]
        der = [ders[d][1] for d in range(3)]
        a_idx, b_idx, c_idx = np.arange(p + 1), np.arange(p + 1), np.arange(p + 1)
        cps = ((spans[0] + a_idx)[None, None, :] + nu_cp * ((spans[1] + b_idx)[None, :, None] + nv_cp * (spans[2] + c_idx)[:, None, None])).reshape(-1)
        grads = []
        for j in range(3):
            fx = der[0] if j == 0 else val[0]
            fy = der[1] if j == 1 else val[1]
            fz = der[2] if j == 2 else val[2]
            grads.append(np.einsum("a,b,c->cba", fx, fy, fz).reshape(-1))
        g = np.stack([np.einsum("l,lim->im", dn, coef[cps]) for dn in grads], axis=1)  # (3 comps, 3 dirs, n_modes)
        eps = np.stack(
            [g[0, 0], g[1, 1], g[2, 2], g[0, 1] + g[1, 0], g[1, 2] + g[2, 1], g[0, 2] + g[2, 0]], axis=0
        )  # (6, n_modes)
        out[n] = (D @ eps).T
    return out
