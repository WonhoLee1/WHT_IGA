"""Master plan v2.1 §10.4: Kratos path vs reference assembler on reduced V cases.

Plan criterion: "for reduced versions of every V case, the Kratos path and the reference
assembler differ by <= 1e-10". Status per V case and all measured numbers:
``docs/verification/cross_validation_10_4.md``.

What is compared for the ``run_linear`` (pipeline) cases:
* matrix: the matrix Kratos factorized (strategy system matrix, reference DOF order) against the
  reference matrix on the DOFs of the solved system, max|dK| / max|K| <= 1e-12 (TOL_K);
* load vector: Kratos' own assembled right-hand side (``_kratos_rhs``: T^T f at u = 0) against the
  reference load vector, <= 1e-10;
* solution: Kratos U against the accurate reference solve ``_ref_solve``, <= 1e-10, plus Cauchy
  stress at every Kratos QP and clamp reactions where present;
* normalized backward error of the Kratos U inside the reference system (correctly rounded
  residual) <= 1e-10. This is a consistency check only: it barely sees load-vector errors (a
  1 + 1e-9 scaled load leaves it below 1e-10 when ||K|| ||u|| >> ||f||), so the identity of the
  discrete system rests on the matrix and load comparisons, never on the backward error alone.

Reference side (no Kratos anywhere in it):
* ``wht_iga.reference.assemble_elastic`` (box stiffness, face traction), ``body_force_vector``,
  ``wht_iga.analyses.thermal.assemble_conduction``, ``kratos_backend.mass.assemble_box_mass``
  (Python consistent mass), scipy eigen / Bossak solves.
* The reference assembler has no trimmed / quadrature-table / penalty path, so this file adds
  three small numpy oracles: ``_qp_stiffness`` (B^T D B over given [x, y, z, w] points, checked
  against ``assemble_stiffness`` to 1e-13 in ``test_oracle_matches_reference_assembler``),
  ``_penalty_exact`` (beta * int N N dA by tensor Gauss, exact) and ``_penalty_same_rule``
  (the same integrand on the pipeline's own face triangles and Dunavant rule).
* ``_ref_solve``: symmetric Jacobi scaling D K D + dense Cholesky, then iterative refinement with
  correctly rounded residuals (``_residual``: Dekker TwoProduct + ``math.fsum``), repeated until
  the correction stops shrinking (see ``_ref_solve`` for why the residual alone is not the stop
  signal). It fails unless the last correction is <= REF_ACC = 1e-13
  (relative), i.e. the reference error sits three orders below the criterion. Both parts are
  checked against exact arithmetic in ``test_exact_residual_and_refined_solve``. An unscaled
  scipy ``spsolve`` is not accurate enough as a reference here: on the FCM ghost alpha = 1e-8
  systems (kappa 3e9..5e9) and the p = 3 trimmed system (kappa 5e8) it is 6e-8 / 1.1e-10 away from
  the refined solution, while the Kratos U is <= 1e-12 away.

Tolerances: the plan's 1e-10 (relative to the max of the reference field) everywhere. Two
exceptions are recorded with their measured values; neither loosens a 1e-10 check:
* BC-A / BC-B penalty faces: the pipeline integrates beta N_i N_j (total degree 4p on a
  triangle) with a Dunavant rule of order 2/4/6 (exact for degree 2p only). Against the
  exact-integration reference the matrix differs by 3e-3..7e-3 and the displacement by
  1e-9..2e-8 (O(1/beta) suppressed); printed as ``OPEN`` and reported as an open issue. The 1e-10
  checks are made against the same-rule reference.
* V07 nu = 0.49999 (kappa 2.4e9): matrix and load vector agree to round-off, but the Kratos U is
  3.5e-9 away from the refined reference (a 50-digit decimal elimination of the same system
  agrees with the refined reference to 6e-25). That gap is the error of Kratos' own
  double-precision linear solve; it is bounded by Kratos' residual, ||u_k - u||_2 <=
  ||K u_k - f||_2 / lambda_min(K), which is asserted. The forward 1e-10 is asserted for the
  nu = 0.3 twin.

The Dunavant order-4 table (``triangle_quadrature(4)``) had a transcription typo in one weight
(0.0549758718227661 for 0.054975871827661; weight sum 1.5e-11 low), fixed in
``kratos_backend/surface.py`` on 2026-11-20; the penalty-total checks below now hold to 1e-13.
"""

from __future__ import annotations

import functools
import math
import time
from fractions import Fraction

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.backends.reference_numba import body_force_vector
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.kratos_backend.surface import triangle_quadrature
from wht_iga.reference.assemble_elastic import assemble_face_traction, assemble_stiffness, isotropic_D
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.penalty import box_face_triangles, face_subdiv, rule_for_degree
from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec, PressureLoad, TractionLoad

pytestmark = pytest.mark.kratos

TOL = 1e-10  # master plan v2.1 §10.4
TOL_K = 1e-12  # assembled matrices: element-by-element sums of identical products (round-off only)
REF_ACC = 1e-3 * TOL  # the reference solve's own error must sit three orders below the criterion


def _rel(a: np.ndarray, ref: np.ndarray) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(ref)).max() / np.abs(np.asarray(ref)).max())


def _mrel(A, ref) -> float:
    return float(abs(A - ref).max() / abs(ref).max())


def _report(tag: str, **vals) -> None:
    body = " ".join(f"{k}={v:.3e}" if isinstance(v, float) else f"{k}={v}" for k, v in vals.items())
    print(f"\nXVAL {tag}: {body}")


# --------------------------------------------------------------------------------------------
# Accurate reference solve (no Kratos)
# --------------------------------------------------------------------------------------------
_SPLIT = 134217729.0  # 2**27 + 1, Dekker splitting constant for binary64


def _two_prod(a: np.ndarray, b: np.ndarray):
    """Dekker TwoProduct: p + e == a * b exactly (binary64, round to nearest, no over/underflow)."""
    p = a * b
    t = _SPLIT * a
    ah = t - (t - a)
    al = a - ah
    t = _SPLIT * b
    bh = t - (t - b)
    bl = b - bh
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def _residual(K, u, f) -> np.ndarray:
    """f - K u with every entry correctly rounded (exact products, exactly rounded fsum per row)."""
    K = sp.csr_matrix(K)
    p, e = _two_prod(K.data, np.asarray(u, dtype=float)[K.indices])
    f = np.asarray(f, dtype=float)
    r = np.empty(K.shape[0])
    for i in range(K.shape[0]):
        s, t = K.indptr[i], K.indptr[i + 1]
        r[i] = math.fsum([float(f[i]), *(-p[s:t]).tolist(), *(-e[s:t]).tolist()])
    return r


def _ref_solve(K, f, max_steps: int = 10):
    """Reference solution of the SPD system K u = f, accurate far beyond eps * kappa.

    D K D (D = diag(K)^-1/2) is factorized by dense Cholesky; the solution is refined with
    correctly rounded residuals until the correction stops shrinking (it no longer halves) or
    reaches eps. The residual alone is not a usable stop signal: for an ill-conditioned K the
    first Cholesky solution already has a residual at the round-off floor (~eps |K| |u|), which
    the refinement does not lower although it improves u by orders (pascal(12): first correction
    4.7e-7 at an unchanged residual). Returns u and {"steps": refinement steps, "dx": relative
    size of the last correction (error left in u is below it), "res": residual max-norm history}.
    Raises unless dx <= REF_ACC, so a reference that did not converge cannot pass silently.
    """
    K = sp.csr_matrix(K)
    f = np.asarray(f, dtype=float)
    d = 1.0 / np.sqrt(K.diagonal())
    cf = sla.cho_factor((sp.diags(d) @ K @ sp.diags(d)).toarray())
    u = d * sla.cho_solve(cf, d * f)
    r = _residual(K, u, f)
    hist = [float(np.abs(r).max())]
    steps, dx, prev = 0, float("inf"), float("inf")
    while steps < max_steps:
        du = d * sla.cho_solve(cf, d * r)
        dx = float(np.abs(du).max() / np.abs(u).max())
        u = u + du
        r = _residual(K, u, f)
        hist.append(float(np.abs(r).max()))
        steps += 1
        if dx <= 2.0 ** -52 or dx > 0.5 * prev:
            break
        prev = dx
    info = {"steps": steps, "dx": dx, "res": hist}
    assert dx <= REF_ACC, info
    return u, info


def _solve_fixed(K, f, fixed: np.ndarray):
    u = np.zeros(K.shape[0])
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    u[free], info = _ref_solve(sp.csr_matrix(K)[free][:, free], f[free])
    return u, info


def _backward_error(K, u, f) -> float:
    """Normalized backward error ||K u - f|| / (||K|| ||u|| + ||f||) (inf norms, rounded residual)."""
    r = _residual(K, u, f)
    nK = float(abs(K).sum(axis=1).max())
    return float(np.abs(r).max() / (nK * np.abs(u).max() + np.abs(f).max()))


# --------------------------------------------------------------------------------------------
# Kratos-free oracles
# --------------------------------------------------------------------------------------------
def _gauss_points(box: BoxSpline, cell: tuple[int, int, int]) -> np.ndarray:
    """(p+1)^3 tensor Gauss points [x, y, z, w] of one cell, x fastest (own implementation)."""
    xg, wg = np.polynomial.legendre.leggauss(box.degree + 1)
    br = [box.breakpoints(d) for d in range(3)]
    ax = [0.5 * (br[d][cell[d]] + br[d][cell[d] + 1]) + 0.5 * (br[d][cell[d] + 1] - br[d][cell[d]]) * xg
          for d in range(3)]
    aw = [0.5 * (br[d][cell[d] + 1] - br[d][cell[d]]) * wg for d in range(3)]
    Z, Y, X = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")
    W = np.einsum("k,j,i->kji", aw[2], aw[1], aw[0])
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel(), W.ravel()])


def _cell_basis(box: BoxSpline, cell, pts: np.ndarray):
    """Basis values, gradients (n, nloc[, 3]) and control points of one cell at ``pts``."""
    p = box.degree
    br = [box.breakpoints(d) for d in range(3)]
    spans, val, der = [], [], []
    for d in range(3):
        kn = np.ascontiguousarray(box.knots[d], dtype=float)
        s = int(find_span(box.n_cp[d], p, 0.5 * (br[d][cell[d]] + br[d][cell[d] + 1]), kn))
        dd = np.array([ders_basis_funs(s, float(x), p, 1, kn) for x in pts[:, d]])
        spans.append(s)
        val.append(dd[:, 0, :])
        der.append(dd[:, 1, :])
    n, nl = pts.shape[0], (p + 1) ** 3
    N = np.einsum("na,nb,nc->ncba", val[0], val[1], val[2]).reshape(n, nl)
    G = np.stack([np.einsum("na,nb,nc->ncba", der[0], val[1], val[2]).reshape(n, nl),
                  np.einsum("na,nb,nc->ncba", val[0], der[1], val[2]).reshape(n, nl),
                  np.einsum("na,nb,nc->ncba", val[0], val[1], der[2]).reshape(n, nl)], axis=2)
    i = np.arange(p + 1)
    nu, nv = box.n_cp[0], box.n_cp[1]
    cps = ((spans[0] - p + i)[None, None, :]
           + nu * ((spans[1] - p + i)[None, :, None] + nv * (spans[2] - p + i)[:, None, None])).ravel()
    return N, G, cps


def _B(G: np.ndarray) -> np.ndarray:
    """Voigt [xx, yy, zz, xy, yz, xz] strain-displacement matrices, engineering shear."""
    n, nl, _ = G.shape
    B = np.zeros((n, 6, 3 * nl))
    gx, gy, gz = G[:, :, 0], G[:, :, 1], G[:, :, 2]
    B[:, 0, 0::3], B[:, 1, 1::3], B[:, 2, 2::3] = gx, gy, gz
    B[:, 3, 0::3], B[:, 3, 1::3] = gy, gx
    B[:, 4, 1::3], B[:, 4, 2::3] = gz, gy
    B[:, 5, 0::3], B[:, 5, 2::3] = gz, gx
    return B


def _qp_stiffness(box: BoxSpline, rules: dict, E: float, nu: float) -> sp.csr_matrix:
    """K = sum_cells sum_qp w * scale * B^T D B for rules {cell: (pts (n, 4), scale)}."""
    D = isotropic_D(E, nu)
    rows, cols, vals = [], [], []
    for cell, (pts, scale) in rules.items():
        pts = np.asarray(pts, dtype=float)
        _, G, cps = _cell_basis(box, cell, pts[:, :3])
        B = _B(G)
        Ke = np.einsum("n,nik,ij,njl->kl", pts[:, 3] * scale, B, D, B)
        dofs = (3 * cps[:, None] + np.arange(3)).ravel()
        rows.append(np.repeat(dofs, len(dofs)))
        cols.append(np.tile(dofs, len(dofs)))
        vals.append(Ke.ravel())
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(box.n_dofs, box.n_dofs)).tocsr()


def _stress(box: BoxSpline, U: np.ndarray, cell, pts: np.ndarray, E: float, nu: float) -> np.ndarray:
    """Cauchy stress D B u at ``pts`` of ``cell`` (Kratos Voigt order)."""
    _, G, cps = _cell_basis(box, cell, np.asarray(pts, dtype=float)[:, :3])
    ue = U.reshape(-1, 3)[cps].ravel()
    return np.einsum("ij,njk,k->ni", isotropic_D(E, nu), _B(G), ue)


def _face_cells(box: BoxSpline, axis: int, side: int):
    br = [box.breakpoints(d) for d in range(3)]
    o = [d for d in range(3) if d != axis]
    fixed = 0 if side == 0 else len(br[axis]) - 2
    for c0 in range(len(br[o[0]]) - 1):
        for c1 in range(len(br[o[1]]) - 1):
            c = [0, 0, 0]
            c[axis], c[o[0]], c[o[1]] = fixed, c0, c1
            yield tuple(c)


def _penalty_from_points(box, beta, comps, g_fn, pts_by_cell) -> tuple[sp.csr_matrix, np.ndarray]:
    rows, cols, vals = [], [], []
    f = np.zeros(box.n_dofs)
    for cell, (X, dA) in pts_by_cell:
        N, _, cps = _cell_basis(box, cell, X)
        for x, w, Nq in zip(X, dA, N):
            g = g_fn(x)
            for d in range(3):
                if comps[d]:
                    dofs = 3 * cps + d
                    rows.append(np.repeat(dofs, len(cps)))
                    cols.append(np.tile(dofs, len(cps)))
                    vals.append((beta * w * np.outer(Nq, Nq)).ravel())
                    f[dofs] += beta * w * g[d] * Nq
    P = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(box.n_dofs, box.n_dofs)).tocsr()
    return P, f


def _penalty_exact(box, face, beta, comps, g_fn):
    """beta int_face N N dA (x) e_d e_d^T and beta int g N dA, tensor Gauss p+1 (exact)."""
    axis, side = face
    xg, wg = np.polynomial.legendre.leggauss(box.degree + 1)
    o = [d for d in range(3) if d != axis]
    br = [box.breakpoints(d) for d in range(3)]
    xf = box.lo[axis] if side == 0 else box.hi[axis]
    out = []
    for cell in _face_cells(box, axis, side):
        a = [0.5 * (br[d][cell[d]] + br[d][cell[d] + 1]) + 0.5 * (br[d][cell[d] + 1] - br[d][cell[d]]) * xg for d in o]
        w = [0.5 * (br[d][cell[d] + 1] - br[d][cell[d]]) * wg for d in o]
        A0, A1 = np.meshgrid(a[0], a[1], indexing="ij")
        X = np.zeros((A0.size, 3))
        X[:, axis], X[:, o[0]], X[:, o[1]] = xf, A0.ravel(), A1.ravel()
        out.append((cell, (X, np.outer(w[0], w[1]).ravel())))
    return _penalty_from_points(box, beta, comps, g_fn, out)


def _penalty_same_rule(box, face, beta, comps, g_fn):
    """Same integrand on the pipeline's per-face-cell triangles + Dunavant rule (pure Python)."""
    axis, side = face
    tris = box_face_triangles(np.asarray(box.lo), np.asarray(box.hi), axis, side, subdiv=face_subdiv(box, axis))
    rp, rw = triangle_quadrature(rule_for_degree(box.degree))
    br = [box.breakpoints(d) for d in range(3)]
    out = []
    for tri in tris:
        tri = np.asarray(tri, dtype=float)
        area = 0.5 * float(np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0])))
        cen = tri.mean(axis=0)
        cell = tuple(min(int(np.searchsorted(br[d], cen[d], side="right")) - 1, len(br[d]) - 2) for d in range(3))
        X = np.array([(1.0 - xi - eta) * tri[0] + xi * tri[1] + eta * tri[2] for xi, eta in rp])
        out.append((cell, (X, 2.0 * area * np.asarray(rw))))
    return _penalty_from_points(box, beta, comps, g_fn, out)


def _bc_g(bc: DirichletBC):
    if bc.field is not None:
        A, b = np.asarray(bc.field.A, float), np.asarray(bc.field.b, float)
        return lambda x: A @ x + b
    v = np.asarray(bc.value, float)
    return lambda x: v


def _reference_box(spec: JobSpec):
    """Reference K (box assembler) and load vector (face tractions, pressures, gravity)."""
    box = BoxSpline(spec.lo, spec.hi, spec.n_cells, spec.degree, spec.continuity)
    K = assemble_stiffness(box, spec.E, spec.nu).tocsr()
    f = np.zeros(box.n_dofs)
    for t in spec.tractions:
        f += assemble_face_traction(box, t.face[0], t.face[1], np.asarray(t.vector, float))
    for pr in spec.pressures:  # pipeline convention (loads.py, S08): traction = +p * n_out
        n_out = np.zeros(3)
        n_out[pr.face[0]] = 1.0 if pr.face[1] == 1 else -1.0
        f += assemble_face_traction(box, pr.face[0], pr.face[1], pr.magnitude * n_out)
    if spec.gravity is not None and spec.density > 0.0:
        f += body_force_vector(box, spec.gravity, spec.density)
    return box, K, f


def _with_penalty(box, K, f, bcs, beta, rule):
    make = _penalty_same_rule if rule == "same" else _penalty_exact
    for bc in bcs:
        P, fp = make(box, tuple(bc.face), beta, tuple(bc.components), _bc_g(bc))
        K, f = K + P, f + fp
    return K.tocsr(), f


def _face_dofs(box: BoxSpline, axis: int, side: int) -> np.ndarray:
    return np.sort((3 * box.boundary_control_points(axis, side)[:, None] + np.arange(3)).ravel())


# --------------------------------------------------------------------------------------------
# Kratos-side extraction
# --------------------------------------------------------------------------------------------
def _run(spec: JobSpec, **kw):
    from wht_iga.solver.pipeline import run_linear

    t0 = time.perf_counter()
    res = run_linear(spec, **kw)
    return res, time.perf_counter() - t0


def _kratos_matrix(res, box: BoxSpline) -> sp.csr_matrix:
    """Kratos system matrix of the last solve in reference DOF order (spline DOFs only)."""
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.elements import equation_permutation, system_matrix_csr

    KM = load_kratos()["KratosMultiphysics"]
    perm = equation_permutation(KM, res.keepalive[3], box)
    return system_matrix_csr(res.strategy)[perm][:, perm].tocsr()


def _kratos_rhs(res, box: BoxSpline) -> np.ndarray:
    """Kratos' own load vector of the solved system, T^T f, in reference DOF order (spline DOFs).

    A fresh ResidualBasedBlockBuilderAndSolver assembles every element and condition of the
    solved model part at u = 0 (linear elements: RHS = external load, the -K u term vanishes) and
    applies the model part's master-slave constraints (T^T b: BC-B virtual-node loads and
    extension slaves land on their masters). Dirichlet rows are not zeroed. The displacement
    state is restored afterwards; call it after U, stresses and reactions were read.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.elements import equation_permutation

    KM = load_kratos()["KratosMultiphysics"]
    mp = res.keepalive[3]
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    perm0 = equation_permutation(KM, mp, box)
    saved = [(nd, [nd.GetSolutionStepValue(v) for v in vars3]) for nd in mp.Nodes]
    try:
        for nd, vals in saved:
            for v, x in zip(vars3, vals):
                if nd.HasDofFor(v):
                    # b(0) is the external load only if no DOF is prescribed to a non-zero value
                    assert not nd.IsFixed(v) or x == 0.0, (nd.Id, x)
                    nd.SetSolutionStepValue(v, 0.0)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        bs = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
        bs.SetUpDofSet(scheme, mp)
        bs.SetUpSystem(mp)
        space = KM.UblasSparseSpace()
        A, dx, b = space.CreateEmptyMatrixPointer(), space.CreateEmptyVectorPointer(), space.CreateEmptyVectorPointer()
        bs.ResizeAndInitializeVectors(scheme, A, dx, b, mp)
        bs.Build(scheme, mp, A, b)
        bs.ApplyConstraints(scheme, mp, A, b)
        bv = np.array([b[i] for i in range(b.Size())])
        perm = equation_permutation(KM, mp, box)
    finally:
        for nd, vals in saved:
            for v, x in zip(vars3, vals):
                if nd.HasDofFor(v):
                    nd.SetSolutionStepValue(v, x)
    assert np.array_equal(perm, perm0)  # same equation numbering as the solved system
    return bv[perm]


def _kratos_reactions(res, box: BoxSpline) -> np.ndarray:
    from wht_iga.kratos_backend import load_kratos

    KM = load_kratos()["KratosMultiphysics"]
    R = np.zeros(box.n_dofs)
    for nd in res.keepalive[3].Nodes:
        I = nd.Id - 1
        if 0 <= I < box.n_ctrl:
            r = nd.GetSolutionStepValue(KM.REACTION)
            R[3 * I:3 * I + 3] = [r[0], r[1], r[2]]
    return R


# --------------------------------------------------------------------------------------------
# 0. the test-local oracles and the reference solve (no Kratos)
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("p", [1, 2, 3])
def test_oracle_matches_reference_assembler(p: int) -> None:
    """_qp_stiffness on full Gauss rules == assemble_stiffness; penalty totals == beta * area."""
    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 1.0, 1.5), (3, 2, 2), p)
    rules = {(i, j, k): (_gauss_points(box, (i, j, k)), 1.0) for k in range(2) for j in range(2) for i in range(3)}
    K0 = assemble_stiffness(box, 210.0, 0.3)
    k_err = abs(_qp_stiffness(box, rules, 210.0, 0.3) - K0).max() / abs(K0).max()
    assert k_err <= 1e-13, k_err  # same products, different summation order
    tot = []
    for make in (_penalty_exact, _penalty_same_rule):
        P, fp = make(box, (1, 1), 7.0, (True, False, True), lambda x: np.array([x[0], 0.0, 2.0]))
        area = 2.0 * 1.5
        # sum_ij N_i N_j = 1 (partition of unity): total = beta * area per selected component;
        # g_x = x integrates to beta * area * 1.0, g_z = 2 to 2 * beta * area. Bound 1e-13: Gauss
        # and the (15-significant-digit) Dunavant tables plus round-off over <= ~10^3 terms.
        e = [abs(P.sum() / (2 * 7.0 * area) - 1.0), abs(fp[0::3].sum() / (7.0 * area * 1.0) - 1.0),
             abs(fp[2::3].sum() / (7.0 * area * 2.0) - 1.0)]
        tot.append(max(e))
        assert max(e) <= 1e-13, e
        assert np.abs(fp[1::3]).max() == 0.0
    _report(f"oracle p={p}", K_rel=float(k_err), pen_total_exact=float(tot[0]), pen_total_same=float(tot[1]))


def test_exact_residual_and_refined_solve() -> None:
    """_residual == exact rational residual rounded once; _ref_solve recovers a known exact solution."""
    rng = np.random.default_rng(7)
    A = (sp.random(40, 40, density=0.3, random_state=3, format="csr") * 1e6 + sp.eye(40)).tocsr()
    u = rng.standard_normal(40)
    f = A @ u  # residual is pure cancellation: |f - A u| ~ eps * |A| |u|
    exact = np.array([float(Fraction(f[i]) - sum(Fraction(A.data[k]) * Fraction(u[A.indices[k]])
                                                for k in range(A.indptr[i], A.indptr[i + 1])))
                      for i in range(40)])
    r = _residual(A, u, f)
    assert np.abs(exact).max() > 0.0
    assert np.array_equal(r, exact), np.abs(r - exact).max()  # both are the correctly rounded value
    # Pascal matrix: SPD, integer entries, kappa = 8.8e11; f = P u_true is exact in binary64
    P = sp.csr_matrix(sla.pascal(12).astype(float))
    u_true = np.array([(-1.0) ** i * (i + 1) for i in range(12)])
    f = P @ u_true
    u_ref, info = _ref_solve(P, f)
    err = _rel(u_ref, u_true)
    err_plain = _rel(spla.spsolve(P.tocsc(), f), u_true)
    _report("ref solve pascal(12)", err=err, err_plain_spsolve=err_plain, ref_ir=info["steps"], ref_dx=info["dx"])
    assert err <= REF_ACC, err


# --------------------------------------------------------------------------------------------
# V01 prescribed-displacement cube (E = 1, nu = 0, u = (x, 2y, 3z))
# --------------------------------------------------------------------------------------------
V01_A = ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))


@pytest.mark.parametrize("p", [1, 2, 3])
def test_v01_bca_affine_cube(p: int) -> None:
    """BC-A (SupportPenaltyCondition) with the affine field on all six faces (pipeline V01 form)."""
    bcs = [DirichletBC(face=(a, s), field=AffineField(A=V01_A, b=(0.0, 0.0, 0.0))) for a in range(3) for s in (0, 1)]
    spec = JobSpec(name=f"xval_v01a_p{p}", lo=(0, 0, 0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=p,
                   E=1.0, nu=0.0, bcs=bcs)
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    Ks, fs = _with_penalty(box, K, f, bcs, res.beta_used, "same")
    Ke, fe = _with_penalty(box, K, f, bcs, res.beta_used, "exact")
    Ak = _kratos_matrix(res, box)
    m_same, m_exact = _mrel(Ak, Ks), _mrel(Ak, Ke)
    fk = _kratos_rhs(res, box)
    f_same, f_exact = _rel(fk, fs), _rel(fk, fe)
    us, info = _ref_solve(Ks, fs)
    u_same, u_exact = _rel(res.U, us), _rel(res.U, _ref_solve(Ke, fe)[0])
    _report(f"V01 BC-A p={p}", u_same_rule=u_same, K_same_rule=m_same, f_same_rule=f_same, f_exact_rule=f_exact,
            OPEN_u_exact_rule=u_exact, OPEN_K_exact_rule=m_exact, beta=float(res.beta_used),
            ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert m_same <= TOL_K, m_same
    assert f_same <= TOL, f_same
    assert u_same <= TOL, u_same


@pytest.mark.parametrize("p", [1, 2, 3])
def test_v01_bcb_normal_components(p: int) -> None:
    """Coreform V01 form: one normal component per face (BC-B virtual node + MPC springs)."""
    bcs = []
    for a in range(3):
        for s in (0, 1):
            comps = [False, False, False]
            comps[a] = True
            val = [0.0, 0.0, 0.0]
            val[a] = float((a + 1) * s)  # u = (x, 2y, 3z) on the unit cube
            bcs.append(DirichletBC(face=(a, s), components=tuple(comps), value=tuple(val)))
    spec = JobSpec(name=f"xval_v01b_p{p}", lo=(0, 0, 0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=p,
                   E=1.0, nu=0.0, bcs=bcs)
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    Ks, fs = _with_penalty(box, K, f, bcs, res.beta_used, "same")
    Ke, fe = _with_penalty(box, K, f, bcs, res.beta_used, "exact")
    m_same = _mrel(_kratos_matrix(res, box), Ks)  # T^T K T on the spline DOFs (virtual nodes eliminated)
    f_same = _rel(_kratos_rhs(res, box), fs)  # virtual-node loads k g mapped by T^T
    us, info = _ref_solve(Ks, fs)
    u_same, u_exact = _rel(res.U, us), _rel(res.U, _ref_solve(Ke, fe)[0])
    _report(f"V01 BC-B p={p}", u_same_rule=u_same, u_exact_rule=u_exact, K_same_rule=m_same, f_same_rule=f_same,
            ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert m_same <= TOL_K, m_same
    assert f_same <= TOL, f_same
    assert u_same <= TOL, u_same
    # uniform boundary deviation: only int N_i dA enters, which both rules integrate exactly
    assert u_exact <= TOL, u_exact


# --------------------------------------------------------------------------------------------
# V03 cantilever (b = h = 0.1, L = 1, E = 200 GPa, nu = 0; coarsest h = (0.1, 0.1, 0.25))
# --------------------------------------------------------------------------------------------
def _v03_spec(p: int, load: str, clamp: str = "strong") -> JobSpec:
    kw = dict(name=f"xval_v03_{load}_p{p}", lo=(-0.05, -0.05, 0.0), hi=(0.05, 0.05, 1.0), n_cells=(1, 1, 4),
              degree=p, E=2.0e11, nu=0.0, density=8000.0)
    if clamp == "strong":
        kw["fixed_faces"] = [(2, 0)]
    else:
        kw["bcs"] = [DirichletBC(face=(2, 0))]
    if load == "traction":  # v03.json traction: 500 Pa along -y on z = L
        kw["tractions"] = [TractionLoad(face=(2, 1), vector=(0.0, -500.0, 0.0))]
    elif load == "pressure":  # v03.json pressure_uniform: 5000 Pa pushing -y on the top face
        kw["pressures"] = [PressureLoad(face=(1, 1), magnitude=-5000.0)]
    else:  # v03.json gravity: rho g = 78480 N/m^3 along -y
        kw["gravity"] = (0.0, -9.81, 0.0)
    return JobSpec(**kw)


@pytest.mark.parametrize("load", ["traction", "pressure", "gravity"])
@pytest.mark.parametrize("p", [1, 2, 3])
def test_v03_cantilever_strong_clamp(p: int, load: str) -> None:
    """Load vector, displacements, Cauchy stress at every Kratos QP and clamp reactions."""
    spec = _v03_spec(p, load)
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    fixed = _face_dofs(box, 2, 0)
    u, info = _solve_fixed(K, f, fixed)
    nx, ny, nz = box.n_cells
    cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    s_ref = np.concatenate([_stress(box, u, c, _gauss_points(box, c), spec.E, spec.nu) for c in cells])
    s_kr = np.concatenate(res.cauchy_stress)
    R_ref = -_residual(K, u, f)[fixed]  # K u - f on the clamp
    R_kr = _kratos_reactions(res, box)[fixed]
    df = _rel(_kratos_rhs(res, box), f)  # all DOFs, clamp rows included (no Dirichlet zeroing)
    du, ds, dr = _rel(res.U, u), _rel(s_kr, s_ref), _rel(R_kr, R_ref)
    _report(f"V03 {load} p={p}", u=du, stress=ds, reaction=dr, rhs=df, ref_ir=info["steps"], ref_dx=info["dx"],
            t_kratos=tk)
    assert s_kr.shape == s_ref.shape
    assert df <= TOL and du <= TOL and ds <= TOL and dr <= TOL, (df, du, ds, dr)


def test_v03_cantilever_penalty_clamp() -> None:
    """Same beam clamped by BC-A penalty (g = 0) instead of strong Fix."""
    spec = _v03_spec(2, "traction", clamp="penalty")
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    Ks, fs = _with_penalty(box, K, f, spec.bcs, res.beta_used, "same")
    Ke, fe = _with_penalty(box, K, f, spec.bcs, res.beta_used, "exact")
    m_same = _mrel(_kratos_matrix(res, box), Ks)
    f_same = _rel(_kratos_rhs(res, box), fs)
    us, info = _ref_solve(Ks, fs)
    u_same, u_exact = _rel(res.U, us), _rel(res.U, _ref_solve(Ke, fe)[0])
    _report("V03 BC-A clamp p=2", u_same_rule=u_same, K_same_rule=m_same, f_same_rule=f_same,
            OPEN_u_exact_rule=u_exact, ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert m_same <= TOL_K and f_same <= TOL and u_same <= TOL, (m_same, f_same, u_same)


def test_symmetry_rollers_bcb_bending() -> None:
    """Per-component BC-B symmetry planes (V02/V06/V07 type) under bending: non-uniform boundary field."""
    bcs = [DirichletBC(face=(0, 0), components=(True, False, False)),
           DirichletBC(face=(1, 0), components=(False, True, False)),
           DirichletBC(face=(2, 0), components=(False, False, True))]
    spec = JobSpec(name="xval_rollers", lo=(0, 0, 0), hi=(4.0, 1.0, 1.0), n_cells=(4, 2, 2), degree=2, E=2.0e5,
                   nu=0.3, bcs=bcs, tractions=[TractionLoad(face=(0, 1), vector=(1.0, 0.0, 0.5))])
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    Ks, fs = _with_penalty(box, K, f, bcs, res.beta_used, "same")
    Ke, fe = _with_penalty(box, K, f, bcs, res.beta_used, "exact")
    m_same = _mrel(_kratos_matrix(res, box), Ks)
    f_same = _rel(_kratos_rhs(res, box), fs)
    us, info = _ref_solve(Ks, fs)
    u_same, u_exact = _rel(res.U, us), _rel(res.U, _ref_solve(Ke, fe)[0])
    _report("BC-B rollers p=2", u_same_rule=u_same, K_same_rule=m_same, f_same_rule=f_same,
            OPEN_u_exact_rule=u_exact, ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert m_same <= TOL_K and f_same <= TOL and u_same <= TOL, (m_same, f_same, u_same)


# --------------------------------------------------------------------------------------------
# V07 near-incompressible (nu = 0.49999, E = 1e5, T = 10) on the 4 x 4 x 0.25 block, untrimmed
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("nu,forward", [(0.3, True), (0.49999, False)])
def test_v07_near_incompressible_block(nu: float, forward: bool) -> None:
    """Small-strain Kratos path (no F-bar) vs reference; nu = 0.3 twin isolates the conditioning.

    nu = 0.49999: the forward gap is Kratos' own double-precision solve error (see the module
    docstring); it is recorded, and the rigorous bound ||u_k - u||_2 <= ||K u_k - f||_2 / lambda_min
    (no tolerance, the reference error is 1e-13 or less) is asserted instead of 1e-10.
    """
    spec = JobSpec(name=f"xval_v07_{nu}", lo=(0, 0, 0), hi=(4.0, 4.0, 0.25), n_cells=(4, 4, 1), degree=2,
                   E=1.0e5, nu=nu, fixed_faces=[(0, 0)], tractions=[TractionLoad(face=(0, 1), vector=(10.0, 0, 0))])
    res, tk = _run(spec)
    box, K, f = _reference_box(spec)
    fixed = _face_dofs(box, 0, 0)
    free = np.setdiff1d(np.arange(box.n_dofs), fixed)
    u, info = _solve_fixed(K, f, fixed)
    Kff = K[free][:, free]
    m_err = _mrel(_kratos_matrix(res, box)[free][:, free], Kff)
    df = _rel(_kratos_rhs(res, box), f)
    eta = _backward_error(Kff, res.U[free], f[free])
    Kd = Kff.toarray()
    kappa = float(np.linalg.cond(Kd))
    lam_min = float(sla.eigvalsh(Kd, subset_by_index=[0, 0])[0])
    gap2 = float(np.linalg.norm(res.U[free] - u[free]) / np.linalg.norm(res.U[free]))
    bound2 = float(np.linalg.norm(_residual(Kff, res.U[free], f[free])) / (lam_min * np.linalg.norm(res.U[free])))
    du = _rel(res.U, u)
    u_plain = np.zeros(box.n_dofs)
    u_plain[free] = spla.spsolve(Kff.tocsc(), f[free])
    _report(f"V07 nu={nu}", u=du, u_plain_spsolve_ref=_rel(res.U, u_plain), K=m_err, rhs=df, backward=eta,
            gap2=gap2, residual_bound2=bound2, kappa=kappa, eps_kappa=2.2e-16 * kappa, ref_ir=info["steps"],
            ref_dx=info["dx"], t_kratos=tk)
    assert m_err <= TOL_K, m_err
    assert df <= TOL, df
    assert eta <= TOL, eta
    assert gap2 <= bound2, (gap2, bound2)
    if forward:
        assert du <= TOL, du


# --------------------------------------------------------------------------------------------
# V02 reduced trimmed: plate with a curved (quarter-circle) notch, cut cells + OUT + extension
# --------------------------------------------------------------------------------------------
def _notched_plate_soup(r: float, A: float, B: float, z_hi: float, n_seg: int = 48):
    """Closed outward soup of {0<=x<=A, 0<=y<=B, 0<=z<=z_hi, x^2+y^2 >= r^2}.

    Walls x = 0, y = 0, z = 0 lie on the grid's lower planes (classify's rays start at
    grid min - 1e-6, which must be outside: same layout as ``quarter_disc_prism_soup``).
    """
    th = np.linspace(0.0, 0.5 * np.pi, n_seg + 1)
    arc = np.column_stack([r * np.cos(th), r * np.sin(th)])
    poly = np.vstack([[[A, 0.0], [A, B], [0.0, B]], arc[::-1]])  # CCW around the region
    m = len(poly)
    nodes = np.vstack([np.column_stack([poly, np.zeros(m)]), np.column_stack([poly, np.full(m, z_hi)])])
    tris = []
    for i in range(m):
        j = (i + 1) % m
        tris += [(i, j, m + j), (i, m + j, m + i)]
    for k in range(m):  # caps: fan from (A, B); every region point is visible from it
        a, b = k, (k + 1) % m
        if 1 not in (a, b):
            tris += [(m + 1, m + a, m + b), (1, b, a)]
    return nodes, np.asarray(tris, dtype=np.int64)


@functools.lru_cache(maxsize=None)
def _notched_table(p: int):
    from wht_iga.solver.trimmed import build_trimmed_table

    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 2.0, 0.5), (4, 4, 1), p)
    nodes, tris = _notched_plate_soup(1.09, 2.5, 2.5, 0.75)  # r = 1.09 leaves eta = 0.004 slivers
    v = nodes[tris]
    vol = float(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2])).sum() / 6.0)
    assert abs(vol - (2.5 * 2.5 - 0.25 * np.pi * 1.09 ** 2) * 0.75) < 1e-3 * vol  # closed + outward
    qp_table, out_cells, _issues, recs, diag = build_trimmed_table(nodes, tris, box, q=min(2 * p, 4),
                                                                   workers=1, cache_dir="off")
    return box, qp_table, [tuple(c) for c in out_cells], recs, diag


def _extension_prolongation(box: BoxSpline, recs, theta: float, fixed: set[int]):
    """u = P u_free as Kratos' MPC elimination sees it: identity on every non-slave DOF,
    slave rows from the extension operator (same operator the pipeline builds)."""
    from wht_iga.operators import compose as C
    from wht_iga.stabilize import extension as E
    from wht_iga.stabilize import support as S

    nx, ny, nz = box.n_cells
    cv = float(np.prod((np.asarray(box.hi) - np.asarray(box.lo)) / np.asarray(box.n_cells)))
    eta = np.zeros(nx * ny * nz)
    for c in recs:
        eta[c["i"] + nx * (c["j"] + ny * c["k"])] = float(c["eta"])
    knots = [np.asarray(k, float) for k in box.knots]
    cls = S.classify_basis(tuple(box.n_cp), knots, eta * cv, np.full(eta.size, cv), theta=theta)
    op = E.build_extension(tuple(box.n_cp), cls.labels, box.degree, knots=knots)
    Td = C.expand_dof_operator(C.compose_thb(op.E), n_comp=3).tocsr()
    mdofs = (np.asarray(op.master_ids).reshape(-1, 1) * 3 + np.arange(3)).ravel()
    sdofs = (np.asarray(op.slave_ids).reshape(-1, 1) * 3 + np.arange(3)).ravel()
    assert len(sdofs) > 0 and not set(sdofs.tolist()) & fixed  # slivers produce slaves; none clamped
    keep = [d for d in range(box.n_dofs) if d not in set(sdofs.tolist()) and d not in fixed]
    col = {d: c for c, d in enumerate(keep)}
    P = sp.lil_matrix((box.n_dofs, len(keep)))
    for d in keep:
        P[d, col[d]] = 1.0
    for s in sdofs:
        r = Td.getrow(int(s))
        for c, w in zip(r.indices, r.data):
            if int(mdofs[c]) in col:
                P[int(s), col[int(mdofs[c])]] += w
    return P.tocsr(), len(sdofs) // 3, np.asarray(keep, dtype=np.int64)


@pytest.mark.parametrize("p,mode", [(1, "drop"), (2, "drop"), (3, "drop"), (1, "ghost"), (2, "ghost")])
def test_v02_trimmed_notched_plate(p: int, mode: str) -> None:
    """run_linear(qp_table, OUT cells, extension theta = 0.2) vs the QP-table oracle + same T.

    Ghost alpha = 1e-8 (kappa 3e9..5e9) and p = 3 (kappa 5e8) are asserted at 1e-10 like the
    rest: the refined reference solve removes the reference-side error an unscaled spsolve has.
    """
    box, qp_table, out_cells, recs, diag = _notched_table(p)
    ghost, theta, E, nu = 1e-8, 0.2, 2.0e5, 0.3
    spec = JobSpec(name=f"xval_v02t_{mode}_p{p}", lo=box.lo, hi=box.hi, n_cells=box.n_cells, degree=p, E=E, nu=nu,
                   fixed_faces=[(0, 1)], tractions=[TractionLoad(face=(1, 1), vector=(0.0, 1.0, 0.0))],
                   out_cells=out_cells, ghost_alpha=ghost, out_mode=mode)
    res, tk = _run(spec, qp_table=qp_table, cell_records=recs, extension={"theta": theta})
    nx, ny, nz = box.n_cells
    cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    rules = {}
    for c in cells:
        if c in qp_table:
            rules[c] = (np.asarray(qp_table[c]), 1.0)
        elif c in out_cells:
            if mode == "ghost":
                rules[c] = (_gauss_points(box, c), ghost)
        else:
            rules[c] = (_gauss_points(box, c), 1.0)
    K = _qp_stiffness(box, rules, E, nu)
    f = assemble_face_traction(box, 1, 1, (0.0, 1.0, 0.0))
    fixed = set(_face_dofs(box, 0, 1).tolist())
    if mode == "drop":  # control points supported only by OUT cells are fixed by the pipeline
        p_ = box.degree
        act = np.zeros(box.n_ctrl, dtype=bool)
        for (i, j, k) in (c for c in cells if c not in out_cells):
            for c_ in range(p_ + 1):
                for b_ in range(p_ + 1):
                    for a_ in range(p_ + 1):
                        act[(i + a_) + box.n_cp[0] * ((j + b_) + box.n_cp[1] * (k + c_))] = True
        fixed |= set((3 * np.flatnonzero(~act)[:, None] + np.arange(3)).ravel().tolist())
    P, n_slave, keep = _extension_prolongation(box, recs, theta, fixed)
    Kr = (P.T @ K @ P).tocsr()
    fr = P.T @ f
    x, info = _ref_solve(Kr, fr)
    u = P @ x
    m_err = _mrel(_kratos_matrix(res, box)[keep][:, keep], Kr)  # Kratos T^T K T on the reduced DOFs
    s_kr = np.concatenate([s for s in res.cauchy_stress if len(s)])
    # ghost cells carry alpha * E in Kratos, so their stress is alpha * D B u
    s_ref = np.concatenate([_stress(box, u, c, rules[c][0], E * rules[c][1], nu) for c in cells if c in rules])
    df = _rel(_kratos_rhs(res, box)[keep], fr)
    xk = res.U[keep]  # Kratos solution in the reduced coordinates (identity rows of P)
    eta = _backward_error(Kr, xk, fr)
    tie_gap = float(np.abs(P @ xk - res.U).max() / np.abs(res.U).max())
    kappa = float(np.linalg.cond(Kr.toarray()))
    du = _rel(res.U, u)
    assert s_kr.shape == s_ref.shape, (s_kr.shape, s_ref.shape)  # one Kratos QP per table row
    ds = _rel(s_kr, s_ref)
    _report(f"V02-trim {mode} p={p}", u=du, u_plain_spsolve_ref=_rel(res.U, P @ spla.spsolve(Kr.tocsc(), fr)),
            stress=ds, K=m_err, rhs=df, backward=eta, ties_gap=tie_gap, kappa=kappa, eps_kappa=2.2e-16 * kappa,
            ref_ir=info["steps"], ref_dx=info["dx"], n_cut=diag["n_cut"], n_out=len(out_cells),
            n_slave_cp=n_slave, t_kratos=tk)
    assert res.tie_info["n_slave"] == n_slave
    assert tie_gap <= 1e-12, tie_gap  # Kratos U honours u = P u_free (slave recovery)
    assert m_err <= TOL_K, m_err
    assert df <= TOL, df
    assert eta <= TOL, eta
    assert du <= TOL and ds <= TOL, (du, ds)


# --------------------------------------------------------------------------------------------
# V04 two-part Cu/W tie (AssemblyBuilder + mpc_matching) vs glued reference boxes
# --------------------------------------------------------------------------------------------
def test_v04_two_part_tie() -> None:
    """Stiffness + tie path only: the reference face-traction vector is fed to Kratos as point loads."""
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import AssemblyBuilder, AssemblySpec, InstanceSpec, PartSpec
    from wht_iga.kratos_backend.interactions import SurfaceMesh, TieSpec
    from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    E_cu, E_w, L, t_y = 1.0e5, 4.0e5, 100.0, -0.01
    b_cu = BoxSpline((-5.0, -5.0, 0.0), (5.0, 0.0, L), (2, 2, 8), 2)
    b_w = BoxSpline((-5.0, 0.0, 0.0), (5.0, 5.0, L), (2, 2, 8), 2)
    t0 = time.perf_counter()
    built = AssemblyBuilder(KM, SM).build(AssemblySpec(
        "xval_v04", parts={"cu": PartSpec("cu", "iga", box=b_cu, young_modulus=E_cu, poisson_ratio=0.0),
                           "w": PartSpec("w", "iga", box=b_w, young_modulus=E_w, poisson_ratio=0.0)},
        instances=[InstanceSpec("cu", "1"), InstanceSpec("w", "1")], use_pardiso=True))
    mp = built.structure
    icu, iw = built.instances["cu_1"], built.instances["w_1"]
    lp = mp.CreateNewProperties(99)
    cid = 900_000
    for box, inst in ((b_cu, icu), (b_w, iw)):
        for I in box.boundary_control_points(2, 0):
            nd = mp.GetNode(int(inst.node_ids[int(I)]))
            for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
                nd.Fix(var)
        fz = assemble_face_traction(box, 2, 1, (0.0, t_y, 0.0)).reshape(-1, 3)
        for I in np.flatnonzero(np.abs(fz).sum(axis=1) > 0.0):
            cond = mp.CreateNewCondition("PointLoadCondition3D1N", cid, [int(inst.node_ids[int(I)])], lp)
            cond.SetValue(SM.POINT_LOAD, [float(v) for v in fz[I]])
            cid += 1
    idx_cu, idx_w = b_cu.boundary_control_points(1, 1), b_w.boundary_control_points(1, 0)
    X_cu, X_w = b_cu.control_points()[idx_cu], b_w.control_points()[idx_w]
    assert np.abs(X_cu - X_w).max() == 0.0  # matching interface control points, same order
    tri = np.array([[0, 1, 2]], dtype=np.int64)
    sm_cu = SurfaceMesh("cu_if", "iga", X_cu, tri, node_ids=np.array([icu.node_ids[int(I)] for I in idx_cu]))
    sm_w = SurfaceMesh("w_if", "iga", X_w, tri, node_ids=np.array([iw.node_ids[int(I)] for I in idx_w]))
    solve_hybrid_iga_fe(KM, SM, built, ties=[TieSpec("xval_tie", sm_cu, sm_w, method="mpc_matching")],
                        use_pardiso=True)
    tk = time.perf_counter() - t0

    def _u(box, inst):
        return np.array([[mp.GetNode(int(n)).GetSolutionStepValue(v) for v in
                          (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)] for n in inst.node_ids]).ravel()

    # reference: global CPs = Cu CPs, then W CPs; W interface CPs merged onto the Cu ones
    n_cu = b_cu.n_ctrl
    gmap = np.empty(b_w.n_ctrl, dtype=np.int64)
    nxt = n_cu
    merged = dict(zip(idx_w.tolist(), idx_cu.tolist()))
    for I in range(b_w.n_ctrl):
        if I in merged:
            gmap[I] = merged[I]
        else:
            gmap[I], nxt = nxt, nxt + 1
    n_g = nxt
    dof_w = (3 * gmap[:, None] + np.arange(3)).ravel()
    Pw = sp.csr_matrix((np.ones(dof_w.size), (np.arange(dof_w.size), dof_w)), shape=(dof_w.size, 3 * n_g))
    Pc = sp.eye(3 * n_cu, 3 * n_g, format="csr")
    K = Pc.T @ assemble_stiffness(b_cu, E_cu, 0.0) @ Pc + Pw.T @ assemble_stiffness(b_w, E_w, 0.0) @ Pw
    f = Pc.T @ assemble_face_traction(b_cu, 2, 1, (0.0, t_y, 0.0)) + Pw.T @ assemble_face_traction(b_w, 2, 1, (0.0, t_y, 0.0))
    fixed = np.unique(np.concatenate([_face_dofs(b_cu, 2, 0), dof_w[_face_dofs(b_w, 2, 0)]]))
    ug, info = _solve_fixed(K.tocsr(), f, fixed)
    u_ref = np.concatenate([Pc @ ug, Pw @ ug])
    u_kr = np.concatenate([_u(b_cu, icu), _u(b_w, iw)])
    du = _rel(u_kr, u_ref)
    _report("V04 tie", u=du, ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert du <= TOL, du


def test_v05_cell_domains_shared_spline() -> None:
    """V05 form: one shared spline, 7 planar material slabs (CellDomainSpec), vs per-cell E oracle.

    Stiffness / cell-domain path only: the reference face-traction vector is fed to Kratos as
    point loads.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.assembly_builder import (
        AssemblyBuilder,
        AssemblySpec,
        CellDomainSpec,
        InstanceSpec,
        PartSpec,
    )
    from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    box = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 7), 2)
    E0, nu = 1.0e11, 0.3
    E_k = [E0 * (1.0 + 0.5 * k) for k in range(7)]  # distinct slab moduli (the G3 test uses one E)
    domains = [CellDomainSpec(cell_id=k, name=f"slab_{k}", young_modulus=E_k[k], poisson_ratio=nu,
                              cells=[(i, j, k) for i in range(2) for j in range(2)]) for k in range(7)]
    t0 = time.perf_counter()
    built = AssemblyBuilder(KM, SM).build(AssemblySpec(
        "xval_v05", parts={"beam": PartSpec("beam", "iga", box=box, young_modulus=E0, poisson_ratio=nu,
                                            cell_domains=domains)},
        instances=[InstanceSpec("beam", "1")], use_pardiso=True))
    mp, inst = built.structure, built.instances["beam_1"]
    for I in box.boundary_control_points(2, 0):
        nd = mp.GetNode(int(inst.node_ids[int(I)]))
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            nd.Fix(var)
    f = assemble_face_traction(box, 2, 1, (0.0, -0.01, 0.0))
    lp = mp.CreateNewProperties(99)
    for n, I in enumerate(np.flatnonzero(np.abs(f.reshape(-1, 3)).sum(axis=1) > 0.0)):
        cond = mp.CreateNewCondition("PointLoadCondition3D1N", 900_000 + n, [int(inst.node_ids[int(I)])], lp)
        cond.SetValue(SM.POINT_LOAD, [float(v) for v in f.reshape(-1, 3)[I]])
    solve_hybrid_iga_fe(KM, SM, built, use_pardiso=True)
    tk = time.perf_counter() - t0
    u_kr = np.array([[mp.GetNode(int(n)).GetSolutionStepValue(v) for v in
                      (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)] for n in inst.node_ids]).ravel()
    rules = {(i, j, k): (_gauss_points(box, (i, j, k)), E_k[k] / E0) for k in range(7) for j in range(2) for i in range(2)}
    u, info = _solve_fixed(_qp_stiffness(box, rules, E0, nu), f, _face_dofs(box, 2, 0))
    du = _rel(u_kr, u)
    _report("V05 7 slabs", u=du, ref_ir=info["steps"], ref_dx=info["dx"], t_kratos=tk)
    assert du <= TOL, du


# --------------------------------------------------------------------------------------------
# V09 modes, V10 Bossak, V-EXP explicit: Kratos K / schemes vs reference K + Python M
# --------------------------------------------------------------------------------------------
def test_v09_cantilever_modes() -> None:
    from wht_iga.dynamics import solve_eigenmodes
    from wht_iga.kratos_backend.mass import assemble_box_mass, face_dofs, solve_kratos_eigenmodes

    box = BoxSpline((-0.05, -0.05, 0.0), (0.05, 0.05, 1.0), (2, 2, 10), 2)  # V09 beam, 4x4x40 -> 2x2x10
    clamp = face_dofs(box, 2, 0)
    t0 = time.perf_counter()
    km = solve_kratos_eigenmodes(box, 2.0e11, 0.0, 8000.0, clamp, n_modes=10)
    tk = time.perf_counter() - t0
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    ref = solve_eigenmodes(assemble_stiffness(box, 2.0e11, 0.0), assemble_box_mass(box, 8000.0), n_modes=10,
                           free_dofs=free)
    fk, fr = np.asarray(km.eig.frequencies_hz), np.asarray(ref.frequencies_hz)
    df = float(np.max(np.abs(fk - fr) / fr))
    _report("V09 modes", f_rel_max=df, f1_hz=float(fr[0]), t_kratos=tk)
    assert df <= TOL, df


@pytest.mark.parametrize("alpha", [-0.05, -0.30])
def test_v10_bossak(alpha: float) -> None:
    """Kratos Bossak (Kratos K, Kratos lumped M) vs scipy Bossak (reference K, row-sum M)."""
    from wht_iga.dynamics import BossakConfig, solve_bossak_linear
    from wht_iga.kratos_backend.dynamics import run_kratos_bossak, smooth_pulse
    from wht_iga.kratos_backend.mass import assemble_box_mass, face_dofs

    E, rho, length = 1.0e6, 1.0e-6, 4.0
    box = BoxSpline((0, 0, 0), (length, 1.0, 1.0), (4, 1, 1), 2)
    clamp = face_dofs(box, 0, 0)
    f0 = assemble_face_traction(box, 0, 1, (1.0, 0.0, 0.0))
    amp = smooth_pulse(4.0 * length / np.sqrt(E / rho))
    dt, n = 4.0 * length / np.sqrt(E / rho) / 20.0, 60
    probe = np.array([3 * int(c) for c in box.boundary_control_points(0, 1)])
    t0 = time.perf_counter()
    kr = run_kratos_bossak(box, E, 0.0, rho, clamp, f0, amp, alpha, dt, n, probe)
    tk = time.perf_counter() - t0
    M = sp.diags(np.asarray(assemble_box_mass(box, rho).sum(axis=1)).ravel())
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    py = solve_bossak_linear(M, None, assemble_stiffness(box, E, 0.0), lambda t: amp(t) * f0,
                             BossakConfig(dt=dt, n_steps=n, alpha_b=alpha), free_dofs=free)
    du, dv, da = (_rel(getattr(kr, k), getattr(py, k)[:, probe]) for k in ("u", "v", "a"))
    _report(f"V10 Bossak alpha={alpha}", u=du, v=dv, a=da, t_kratos=tk)
    assert max(du, dv, da) <= TOL, (du, dv, da)


@pytest.mark.parametrize("p", [1, 2])
def test_vexp_explicit_bar(p: int) -> None:
    """Kratos central differences (S24b wave model) vs leapfrog on reference K + row-sum M."""
    from tests.kratos_spikes.test_s24b_hybrid_explicit import _E, _RHO, _SIGMA0, _apply_tip_traction, _build_wave_model
    from wht_iga.dynamics import explicit as EX
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.mass import assemble_box_mass

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    nx, n_steps = 10, 100  # CFL 0.2: t_end = 100 dt = 2 L / c (reflection included)
    box = BoxSpline((0, 0, 0), (float(nx), 1.0, 1.0), (nx, 1, 1), p)
    dt = 0.2 * 1.0 / np.sqrt(_E / _RHO)
    t0 = time.perf_counter()
    mp = _build_wave_model(KM, SM, box, 1.0, False, f"XVAL_VEXP_p{p}")["mp"]
    _apply_tip_traction(KM, SM, mp, box, nx, 1.0, False, p)
    strat = SM.MechanicalExplicitStrategy(mp, SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0),
                                          False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    nodes = [mp.GetNode(I + 1) for I in range(box.n_ctrl)]
    hist: list[list[float]] = []

    def probe(_mp):
        hist.append([nd.GetSolutionStepValue(KM.DISPLACEMENT_X) for nd in nodes])
        return hist[-1][-1]

    EX.run_explicit_history(KM, mp, strat, dt, n_steps, probe_fn=probe)
    tk = time.perf_counter() - t0
    xd = np.array([3 * I for I in range(box.n_ctrl) if I % box.n_cp[0] != 0])  # y, z fixed; x = 0 fixed
    Kxx = assemble_stiffness(box, _E, 0.0).tocsr()[xd][:, xd]
    m = np.asarray(assemble_box_mass(box, _RHO).sum(axis=1)).ravel()[xd]
    f = assemble_face_traction(box, 0, 1, (_SIGMA0, 0.0, 0.0))[xd]
    u, v = np.zeros(len(xd)), 0.5 * dt * f / m  # Kratos CD start: v_1/2 = (dt/2) a_0
    ref = np.zeros((n_steps, box.n_ctrl))
    for i in range(n_steps):
        u = u + dt * v
        ref[i, xd // 3] = u
        v = v + dt * (f - Kxx @ u) / m
    du = _rel(np.array(hist), ref)
    _report(f"V-EXP p={p}", u_history=du, t_kratos=tk)
    assert du <= TOL, du


# --------------------------------------------------------------------------------------------
# V11 region stiffness (Kratos QP elements on the hot-spot partition), V12 conduction
# --------------------------------------------------------------------------------------------
def test_v11_region_stiffness() -> None:
    """kratos_region_stiffness (V11 Nitsche solve input) vs the QP oracle, both regions."""
    from wht_iga.analyses.bimaterial import kratos_region_stiffness, quarter_disc_bimaterial_regions

    box = BoxSpline((0.0, 0.0, 0.0), (10.0, 10.0, 0.25), (11, 11, 1), 1)  # V11 44 x 44 -> 11 x 11
    reg = quarter_disc_bimaterial_regions(box, radius=1.0)
    t0 = time.perf_counter()
    errs = []
    for r in (0, 1):
        Kk = kratos_region_stiffness(box, 1.0e5, 0.3, reg.cut_tables[r], reg.full_cells[r])
        rules = {c: (_gauss_points(box, c), 1.0) for c in reg.full_cells[r]}
        rules.update({c: (np.asarray(t), 1.0) for c, t in reg.cut_tables[r].items() if len(t)})
        Kr = _qp_stiffness(box, rules, 1.0e5, 0.3)
        errs.append(float(abs(Kk - Kr).max() / abs(Kr).max()))
    _report("V11 region K", K_region0=errs[0], K_region1=errs[1], n_cut=len(reg.cut_tables[0]),
            t_kratos=time.perf_counter() - t0)
    assert max(errs) <= TOL_K, errs


@pytest.mark.parametrize("p", [1, 2])
def test_v12_conduction(p: int) -> None:
    """Kratos LaplacianElement K_th vs Python assembly, matrix and steady temperature field."""
    from wht_iga.analyses.thermal import assemble_conduction, solve_with_dirichlet
    from wht_iga.kratos_backend.thermal import kratos_box_conductivity

    box = BoxSpline((0.0, 0.0, 0.0), (0.1, 0.02, 0.01), (10, 2, 1), p)  # T3 bar, k = 35
    t0 = time.perf_counter()
    Kk = kratos_box_conductivity(box, 35.0)
    tk = time.perf_counter() - t0
    Kr = assemble_conduction(box, 35.0)
    left, top = box.boundary_control_points(0, 0), box.boundary_control_points(1, 1)
    fixed = np.unique(np.concatenate([left, top]))
    vals = np.where(np.isin(fixed, top), 100.0 * box.control_points()[fixed, 0] / 0.1, 0.0)  # 2D field
    Tk = solve_with_dirichlet(Kk, np.zeros(box.n_ctrl), fixed, vals)
    Tr = solve_with_dirichlet(Kr, np.zeros(box.n_ctrl), fixed, vals)
    dk, dT = float(abs(Kk - Kr).max() / abs(Kr).max()), _rel(Tk, Tr)
    _report(f"V12 conduction p={p}", K=dk, T=dT, t_kratos=tk)
    assert dk <= TOL_K and dT <= TOL, (dk, dT)
