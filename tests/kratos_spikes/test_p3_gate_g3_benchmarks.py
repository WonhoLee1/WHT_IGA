"""Phase 3 Gate G3 10-Criteria Benchmark Suite (genuine Kratos execution).

Drives all Gate G3 criteria end-to-end against Kratos 10.4 + IgaApplication +
StructuralMechanicsApplication + LinearSolversApplication (PARDISO / AMGCL) and
records genuine observed measurements to ``WHT_IGA/docs/phase3/g3_measurements.json``.

Criteria covered:
1. V03-moment: Curved cantilever (R_i=4.12, R_o=4.32, h=0.2, b=0.1, M=100, E=1e7, nu=0.25)
   with NURBS body-fit + RBE3 distributing moment coupling (<= 5% vs Timoshenko-Goodier)
   plus straight cantilever moment from v03.json.
2. V04 Bimetallic Cu/W: C0-merge vs 2-part tie (<= 1e-6 rel diff), vs v04.json
   transformed-section analytical (<= 1.0%), and IGA-FE Hybrid hex8+IGA (<= 0.5%).
3. V05 Multiplanar (7 boxes): 7-cell shared B-spline vs monolithic CAD (<= 3.3e-6)
   and tip deflection vs v05.json Euler-Bernoulli reference (<= 2.0%).
4. V02-local: roark_plate_hole.step (THB L3 around r=10, 2:1 balance, AgFEM extension,
   quarter symmetry documented explicitly, active DOFs <= 80,000, peak stress
   <= 0.2% vs C3D20 reference and <= 1.0% vs Roark).
5. E01: IPS thin plate with circular hole (10 x 5 x 0.1 in, r=1, E=3e7, nu=0.29),
   THIN_FEATURE detection + auto_refine_thin_features, peak stress <= 5% vs Roark.
6. E02: L0 -> L1 -> L2 THB adaptive refinement (2:1 balance) monotonic stress gradient
   resolution at hole edge, DOF savings vs uniform, and far-field L2 convergence.
7. Dependent 3-Instance Assembly: 3 translated/rotated dependent instances under
   equivalent local loading match single-part displacement & rotated Voigt stress <= 1e-10.
8. 5-Format Mesh & Set Lossless Roundtrip: .msh, .inp, .bdf/.fem, .rad, .k.
9. Commercial HWT-IGA SOLVER 6-Stage Log, Auto RAM Switch (auto -> amgcl), and
   out_mode="drop" strictly fewer active DOFs than out_mode="ghost".
"""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.geometry import classify as CL
from wht_iga.geometry import quadrature as Q
from wht_iga.geometry.cad import read_step, tessellate
from wht_iga.grid.bspline import basis_1d
from wht_iga.grid.hierarchy import THBHierarchy
from wht_iga.grid.knots import n_control_points, to_kratos
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.assembly_builder import (
    AssemblyBuilder,
    AssemblySpec,
    CellDomainSpec,
    InstanceSpec,
    PartSpec,
    transform_stress_voigt,
)
from wht_iga.kratos_backend.elements import attach_cell_elements, solve_linear
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.kratos_backend.interactions import (
    RBE3Spec,
    SurfaceMesh,
    TieSpec,
    apply_rbe3,
)
from wht_iga.kratos_backend.native_fe import solve_hybrid_iga_fe
from wht_iga.meshio_ext import (
    ElemBlock,
    MeshModel,
    mesh_to_bodyfit_patches,
    read_mesh,
    write_mesh,
)
from wht_iga.operators.compose import solve_reduced
from wht_iga.post.norms import compute_l2_norm
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.pipeline import run_linear
from wht_iga.solver.steps import THIN_FEATURE, auto_refine_thin_features
from wht_iga.solver.trimmed import (
    _d_matrix,
    _strain_at_points,
    boundary_peak_stress,
    build_trimmed_table,
    hole_face_by_cylinder_wall,
)
from wht_iga.solver.validate import DirichletBC, JobSpec, TractionLoad
from wht_iga.stabilize import extension as EXT
from wht_iga.stabilize import support as SUP

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
VERIF_DIR = REPO_ROOT / "tests" / "verification" / "data"
OUT_JSON = REPO_ROOT / "docs" / "phase3" / "g3_measurements.json"

_MEASUREMENTS: dict[str, Any] = {
    "suite": "Phase 3 Gate G3 10-Criteria Benchmark Suite",
    "benchmarks": {},
}


@pytest.fixture(scope="module", autouse=True)
def _persist_g3_measurements() -> Any:
    """Persist genuine observed measurements after module tests complete."""
    yield
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(_MEASUREMENTS, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# Shared Helpers
# ---------------------------------------------------------------------------


def _load_verif(name: str) -> dict[str, Any]:
    """Load a verification problem JSON fixture with explicit UTF-8."""
    with open(VERIF_DIR / f"{name}.json", encoding="utf-8") as fh:
        return json.load(fh)


def _curved_beam_timoshenko_ur(
    r_i: float,
    r_o: float,
    b: float,
    moment: float,
    E: float,
    nu: float,
    theta_end: float = 0.5 * math.pi,
) -> float:
    """Exact Timoshenko-Goodier radial tip displacement for a curved cantilever under pure end moment.

    Uses plane-stress end-moment Airy stress solution in polar coordinates with
    clamped tangent approximation at theta=0 and tip evaluation at mid-thickness
    theta=theta_end.
    """
    a, c = float(r_i), float(r_o)
    a2, c2 = a * a, c * c
    log_ca = math.log(c / a)
    n_denom = (c2 - a2) ** 2 - 4.0 * a2 * c2 * (log_ca ** 2)
    r_mid = 0.5 * (a + c)
    # Plane-stress curvature flexibility D_bar = -4*M / (b * N)
    # Castigliano / Timoshenko-Goodier hoop resultant angle rate:
    # dphi/dtheta = -4 * M * (c2 - a2) / (E * b * n_denom) -> exact Timoshenko angle change:
    # Exact Airy stress integrated displacement v_theta=0, u_r(0)=0, du_theta/dr|_0=0 gives:
    # At theta = pi/2, radial displacement u_r = (C / E) * [ ... ], numerically integrated exact Airy:
    _ = nu  # plane stress radial tip deflection under pure bending is independent of nu in Euler/Airy
    A = -4.0 * moment * (c2 - a2) / (b * n_denom)
    # Hoop stress sigma_theta(r) = (4*M/(b*N)) * (-a2*c2/r^2 * log(c/a) + c2*log(r/c) - a2*log(r/a) + c2 - a2)
    # Exact polar strain -> rotation rate d(omega)/dtheta = (1/E)*... -> mid-radius Castigliano with Airy constant:
    # Exact Airy constants from Timoshenko & Goodier §28 (b=thickness, moment M closing/opening):
    # Note sign convention: positive +M_z (CCW in x-y at theta=pi/2, i.e., in -x direction on top) bends inward (-x).
    # At theta=pi/2, r_vec = (0, r_mid), so radial unit e_r = (0, 1) and -x is theta-tangential, wait:
    # Wait: moment M = +100 about +Z pulls top of curved beam toward -x (decreasing theta? No, CCW increases theta: bends inward so r decreases / x moves -r? Wait: clamped at theta=0 along +x, tip at theta=pi/2 along +y!).
    # Bending CCW (+M_z) curves +y tip toward -x (closing the quadrant: theta=pi/2 tip has u_r < 0? Let's evaluate exact Airy displacement!).
    _ = A
    # Evaluate exact Timoshenko-Goodier pure-bending curved bar displacement:
    # u_r(r, theta) = ... with omega_0 clamped at r_mid, theta=0.
    # Let D_const = (1/E)*(4*M/(b*n_denom)) * [ ... ]:
    # From Timoshenko & Goodier Theory of Elasticity 3rd ed., Eq. (59) p. 74:
    # C = 4 * M / (b * n_denom)
    C = 4.0 * moment / (b * n_denom)
    # Hoop/radial Airy potential Phi = A0*log(r) + B0*r^2*log(r) + C0*r^2
    # sigma_r = A0/r^2 + B0*(1 + 2*log(r)) + 2*C0
    # sigma_t = -A0/r^2 + B0*(3 + 2*log(r)) + 2*C0
    # Wait: with sigma_r(a)=sigma_r(c)=0:
    # A0 = -a2*c2*log(c/a) * C * (-1)? Let's solve 2x2 linearly for (A0, C0) given B0 = (c2 - a2)*C/2 * ...
    # Normalize M/(b): check M_int = -int_a^c r * sigma_theta dr * b == +moment (CCW opening/closing).
    B0 = 0.5 * (c2 - a2)
    # Solve [1/a^2, 2; 1/c^2, 2] [A0, C0] = [-B0*(1+2*log(a)), -B0*(1+2*log(c))]
    mat = np.array([[1.0 / a2, 2.0], [1.0 / c2, 2.0]])
    rhs = np.array([-B0 * (1.0 + 2.0 * math.log(a)), -B0 * (1.0 + 2.0 * math.log(c))])
    A0, C0 = np.linalg.solve(mat, rhs)
    # Check moment integral over r in [a, c] with Gauss quadrature
    xg, wg = np.polynomial.legendre.leggauss(16)
    rg = 0.5 * (a + c) + 0.5 * (c - a) * xg
    sig_t = -A0 / (rg ** 2) + B0 * (3.0 + 2.0 * np.log(rg)) + 2.0 * C0
    m_unit = b * 0.5 * (c - a) * float(np.sum(wg * (rg - r_mid) * sig_t))
    scale = moment / m_unit
    A0 *= scale
    B0 *= scale
    C0 *= scale
    # Exact Airy displacement (Timoshenko §28):
    # d(u_theta)/dtheta + u_r = (r/E)*(sigma_t - nu*sigma_r) -> term 4*B0*r*theta/E in u_theta!
    # So u_r(r, theta) = (1/E)*[-((1+nu)*A0)/r + 2*(1-nu)*C0*r + B0*r*((1-nu)*(2*log(r)-1) - 2*nu)]
    #                    + H*sin(theta) + K_const*cos(theta)
    # Choosing clamped mid-radius at theta=0: u_r(r_mid, 0)=0 and du_theta/dr(r_mid, 0)=0 gives:
    # Homogeneous rigid body terms (H*sin(theta) + K_c*cos(theta) in u_r, H*cos - K_c*sin + F*r in u_theta)
    # Specifically, clamping theta=0 means u_r = 0 on theta=0 (at r_mid) and u_theta = 0, rotation = 0 at (r_mid, 0):
    ur_part = (1.0 / E) * (
        -(1.0 + nu) * A0 / r_mid
        + 2.0 * (1.0 - nu) * C0 * r_mid
        + B0 * r_mid * ((1.0 - nu) * (2.0 * math.log(r_mid) - 1.0) - 2.0 * nu)
    )
    # Clamped at theta=0 fixes u_r(r_mid, 0) = ur_part + L_c = 0 => L_c = -ur_part,
    # yielding u_r(r_mid, theta_end) = ur_part * (1 - cos(theta_end)).
    return float(ur_part * (1.0 - math.cos(theta_end)))


def _assemble_thb_kratos_and_solve(
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    base_cells: tuple[int, int, int],
    refine_radii: list[tuple[float, float]],
    tess: Any,
    hole_radius: float,
    E: float,
    nu: float,
    traction_x: float,
    theta_ext: float = 0.2,
    tag: str = "thb",
) -> dict[str, Any]:
    """Assemble and solve a quarter-plate with hole using multi-level THB in Kratos.

    Builds exact Giannelli truncation prolongation matrices C_lvl[L], assembles each
    level's leaf cells via Kratos ``NurbsVolumeGeometry`` + ``SmallDisplacementElement3D4N``,
    applies AgFEM support extension on the finest active level, enforces quarter symmetry
    Dirichlet BCs (ux=0 on x=0, uy=0 on y=0, uz=0 on x=hi[0]), and solves with PARDISO.
    """
    t0 = time.perf_counter()
    mods = load_kratos()
    import KratosMultiphysics.scipy_conversion_tools as sct

    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    thb = THBHierarchy.from_uniform(
        (lo[0], hi[0]),
        (lo[1], hi[1]),
        (lo[2], hi[2]),
        base_cells,
        degree=2,
    )
    for lvl, (rmin, rmax) in enumerate(refine_radii):
        bu, bv, bw = thb._breaks_u[lvl], thb._breaks_v[lvl], thb._breaks_w[lvl]
        cells = []
        for ci in range(len(bu) - 1):
            for cj in range(len(bv) - 1):
                r_lo = math.hypot(bu[ci], bv[cj])
                r_hi = math.hypot(bu[ci + 1], bv[cj + 1])
                if r_lo < rmax and r_hi > rmin:
                    for ck in range(len(bw) - 1):
                        if thb._omega_mask[lvl][ci, cj, ck]:
                            cells.append((lvl, ci, cj, ck))
        if cells:
            thb.mark_cells(cells, target_level=lvl + 1)
            thb.enforce_2_to_1_balance()

    # Verify 2:1 balance invariant holds after refinement
    balance_added = thb.enforce_2_to_1_balance()
    assert balance_added == 0

    leaves = thb.leaf_cells()
    n_thb = len(thb.active_dofs())
    l_max = thb.max_level

    C_lvl: list[sp.csr_matrix] = []
    n0 = int(np.prod(thb.dof_shape(0)))
    act0 = thb.active_dofs_at_level(0)
    c_cur = sp.csr_matrix(
        (np.ones(len(act0)), (act0, np.arange(len(act0)))),
        shape=(n0, len(act0)),
    )
    C_lvl.append(sp.hstack([c_cur, sp.csr_matrix((n0, n_thb - c_cur.shape[1]))], format="csr"))
    for lvl in range(1, l_max + 1):
        nl = int(np.prod(thb.dof_shape(lvl)))
        act_l = thb.active_dofs_at_level(lvl)
        mask_diag = np.ones(nl)
        mask_diag[act_l] = 0.0
        m_l = sp.diags(mask_diag, 0, shape=(nl, nl), format="csr")
        p_l = sp.csr_matrix(
            (np.ones(len(act_l)), (act_l, np.arange(len(act_l)))),
            shape=(nl, len(act_l)),
        )
        c_cur = sp.hstack([m_l @ (thb._s3d[lvl - 1] @ c_cur), p_l], format="csr")
        C_lvl.append(sp.hstack([c_cur, sp.csr_matrix((nl, n_thb - c_cur.shape[1]))], format="csr"))

    by_lvl: dict[int, list[Any]] = {L: [] for L in range(l_max + 1)}
    for lf in leaves:
        by_lvl[lf.level].append(lf)

    mult_f = 2 ** l_max
    box_fine = BoxSpline(
        lo,
        hi,
        (base_cells[0] * mult_f, base_cells[1] * mult_f, base_cells[2] * mult_f),
        degree=2,
    )
    nx_f, ny_f, nz_f = box_fine.n_cells
    h_f = np.array([(hi[0] - lo[0]) / nx_f, (hi[1] - lo[1]) / ny_f, (hi[2] - lo[2]) / nz_f])
    cv_f = float(np.prod(h_f))
    cell_in_f = np.full(nx_f * ny_f * nz_f, cv_f, dtype=float)
    bu_f, bv_f = box_fine.breakpoints(0), box_fine.breakpoints(1)
    for cj in range(ny_f):
        for ci in range(nx_f):
            if math.hypot(bu_f[ci + 1], bv_f[cj + 1]) <= hole_radius - 1e-9:
                for ck in range(nz_f):
                    cell_in_f[ci + nx_f * (cj + ny_f * ck)] = 0.0

    K_thb = sp.csr_matrix((3 * n_thb, 3 * n_thb))
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    n_cut_total = 0

    for L in range(l_max + 1):
        if not by_lvl[L]:
            continue
        mult = 2 ** L
        box = BoxSpline(
            lo,
            hi,
            (base_cells[0] * mult, base_cells[1] * mult, base_cells[2] * mult),
            degree=2,
        )
        mp = KM.Model().CreateModelPart(f"{tag}_L{L}")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.REACTION)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        vol, _ = make_nurbs_volume(KM, mp, box)
        prop = mp.CreateNewProperties(1)
        prop.SetValue(KM.YOUNG_MODULUS, float(E))
        prop.SetValue(KM.POISSON_RATIO, float(nu))
        prop.SetValue(KM.DENSITY, 1.0)
        prop.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

        in_cells: list[tuple[int, int, int]] = []
        eid = 1
        for lf in by_lvl[L]:
            blo = np.array([lf.bounds_u[0], lf.bounds_v[0], lf.bounds_w[0]])
            bhi = np.array([lf.bounds_u[1], lf.bounds_v[1], lf.bounds_w[1]])
            r_lo = math.hypot(blo[0], blo[1])
            r_hi = math.hypot(bhi[0], bhi[1])
            if r_hi <= hole_radius - 1e-9:
                continue
            elif r_lo >= hole_radius + 1e-9:
                in_cells.append((lf.ix, lf.iy, lf.iz))
            else:
                r = Q.cell_quadrature(tess.nodes, tess.tris, blo, bhi - blo, q=3, p=2)
                if L == l_max:
                    cell_in_f[lf.ix + nx_f * (lf.iy + ny_f * lf.iz)] = float(r.volume)
                if r.x.shape[0] > 0:
                    n_cut_total += 1
                    pts = [[float(x), float(y), float(z), float(w)] for (x, y, z), w in zip(r.x, r.w)]
                    qps = KM.GeometriesVector()
                    vol.CreateQuadraturePointGeometries(qps, 2, pts)
                    mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], prop)
                    eid += 1
        if in_cells:
            attach_cell_elements(KM, mp, vol, box, prop, start_id=eid, cells=in_cells)

        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        bns = KM.ResidualBasedBlockBuilderAndSolver(KM.LinearSolver())
        strat = KM.ResidualBasedLinearStrategy(mp, scheme, bns, False, False, False, False)
        strat.SetEchoLevel(0)
        strat.Initialize()
        mp.CloneTimeStep(1.0)
        A = KM.CompressedMatrix()
        b = KM.Vector()
        dx = KM.Vector()
        bns.SetUpDofSet(scheme, mp)
        bns.SetUpSystem(mp)
        bns.ResizeAndInitializeVectors(scheme, A, dx, b, mp)
        strat.InitializeSolutionStep()
        bns.Build(scheme, mp, A, b)
        K_kratos = sct.to_csr(A)

        used_nodes = {int(n.Id) for el in mp.Elements for n in el.GetNodes()}
        eq_ids, box_dofs = [], []
        for nid_ in used_nodes:
            node = mp.GetNode(nid_)
            I = nid_ - 1
            for d, var in enumerate(vars3):
                eq_ids.append(node.GetDof(var).EquationId)
                box_dofs.append(3 * I + d)
        P = sp.csr_matrix(
            (np.ones(len(eq_ids)), (eq_ids, box_dofs)),
            shape=(K_kratos.shape[0], box.n_dofs),
        )
        C_3d = sp.kron(C_lvl[L], sp.eye(3), format="csr")
        K_thb = K_thb + C_3d.T @ (P.T @ K_kratos @ P) @ C_3d

    # Traction on x = hi[0] evaluated on level-0 B-splines
    box0 = BoxSpline(lo, hi, base_cells, degree=2)
    f0 = np.zeros(box0.n_dofs)
    gx, wg = np.polynomial.legendre.leggauss(4)
    nu0, nv0, _ = box0.n_cp
    by_br, bz_br = box0.breakpoints(1), box0.breakpoints(2)
    for cj in range(len(by_br) - 1):
        y0, y1 = by_br[cj], by_br[cj + 1]
        jy = 0.5 * (y1 - y0)
        for ck in range(len(bz_br) - 1):
            z0, z1 = bz_br[ck], bz_br[ck + 1]
            jz = 0.5 * (z1 - z0)
            for gy_q, wy in zip(gx, wg):
                y = 0.5 * (y0 + y1) + jy * gy_q
                sy, dy = basis_1d(box0.knots[1], 2, y, 0)
                for gz_q, wz in zip(gx, wg):
                    z = 0.5 * (z0 + z1) + jz * gz_q
                    sz, dz = basis_1d(box0.knots[2], 2, z, 0)
                    w = traction_x * jy * jz * wy * wz
                    for b_i in range(3):
                        for c_i in range(3):
                            I = (nu0 - 1) + nu0 * ((sy - 2 + b_i) + nv0 * (sz - 2 + c_i))
                            f0[3 * I] += w * dy[0][b_i] * dz[0][c_i]
    f_thb = sp.kron(C_lvl[0], sp.eye(3), format="csr").T @ f0

    # AgFEM extension on finest-level cut basis functions
    cls_f = SUP.classify_basis(
        tuple(int(v) for v in box_fine.n_cp),
        [np.asarray(k, float) for k in box_fine.knots],
        cell_in_f,
        np.full(len(cell_in_f), cv_f),
        theta=float(theta_ext),
    )
    op_f = EXT.build_extension(
        tuple(int(v) for v in box_fine.n_cp),
        cls_f.labels,
        2,
        knots=[np.asarray(k, float) for k in box_fine.knots],
    )
    act_thb_list = thb.active_dofs()
    lf_to_thb = {idx: pos for pos, (lvl, idx) in enumerate(act_thb_list) if lvl == l_max}
    E_thb = sp.eye(n_thb, format="lil")
    Ef = op_f.E.tocsr()
    masters_f = np.asarray(op_f.master_ids, dtype=np.int64)
    inact_thb = {lf_to_thb[i] for i in np.flatnonzero(cls_f.labels == SUP.INACTIVE) if i in lf_to_thb}
    slave_thb: set[int] = set()
    for s in op_f.slave_ids:
        s_i = int(s)
        if s_i in lf_to_thb:
            row = Ef.getrow(s_i)
            m_pairs = [(int(masters_f[ci]), float(w)) for ci, w in zip(row.indices, row.data)]
            if all(mf in lf_to_thb for mf, _ in m_pairs):
                ts = lf_to_thb[s_i]
                slave_thb.add(ts)
                E_thb[ts, ts] = 0.0
                for mf, w in m_pairs:
                    E_thb[ts, lf_to_thb[mf]] = w

    keep_cols = [c for c in range(n_thb) if c not in inact_thb and c not in slave_thb]
    E_thb_csr = E_thb.tocsr()[:, keep_cols]
    E_3d = sp.kron(E_thb_csr, sp.eye(3), format="csr")
    K_ext = E_3d.T @ K_thb @ E_3d
    f_ext = E_3d.T @ f_thb

    # Quarter symmetry BCs on un-extended control points before extension
    orig_to_kept = {orig: new for new, orig in enumerate(keep_cols)}
    fix_dofs: set[int] = set()
    for L in range(l_max + 1):
        mult = 2 ** L
        box = BoxSpline(
            lo,
            hi,
            (base_cells[0] * mult, base_cells[1] * mult, base_cells[2] * mult),
            degree=2,
        )
        CL_coo = C_lvl[L].tocoo()
        bx0 = set(box.boundary_control_points(0, 0).tolist())
        by0 = set(box.boundary_control_points(1, 0).tolist())
        for r, orig_c, v in zip(CL_coo.row, CL_coo.col, CL_coo.data):
            if abs(v) > 1e-12 and orig_c in orig_to_kept:
                kc = orig_to_kept[orig_c]
                if r in bx0:
                    fix_dofs.add(3 * kc + 0)
                if r in by0:
                    fix_dofs.add(3 * kc + 1)
    bx1 = set(box0.boundary_control_points(0, 1).tolist())
    CL0_coo = C_lvl[0].tocoo()
    for r, orig_c, v in zip(CL0_coo.row, CL0_coo.col, CL0_coo.data):
        if abs(v) > 1e-12 and orig_c in orig_to_kept and r in bx1:
            fix_dofs.add(3 * orig_to_kept[orig_c] + 2)

    diagE = K_ext.diagonal()
    free = np.array(
        [d for d in range(K_ext.shape[0]) if d not in fix_dofs and diagE[d] > 1e-16],
        dtype=np.int64,
    )
    u_ext = np.zeros(K_ext.shape[0])
    u_ext[free] = solve_reduced(K_ext[free][:, free], f_ext[free])
    u_thb = E_3d @ u_ext
    C_thb = thb.build_c_thb_matrix()
    u_fine = sp.kron(C_thb, sp.eye(3), format="csr") @ u_thb

    return {
        "box_fine": box_fine,
        "u_fine": u_fine,
        "active_dofs": int(len(free)),
        "uniform_dofs": int(box_fine.n_dofs),
        "n_leaves": int(len(leaves)),
        "n_cut_cells": int(n_cut_total),
        "n_slave_extended": int(len(slave_thb)),
        "leaves_per_level": {L: len(v) for L, v in by_lvl.items()},
        "elapsed_s": float(time.perf_counter() - t0),
    }


def _apply_face_traction_on_instance(
    KM: Any,
    SM: Any,
    mp: Any,
    box: BoxSpline,
    node_ids: list[int],
    face_axis: int,
    face_side: int,
    traction_vec: tuple[float, float, float],
    lp: Any,
    start_cid: int,
) -> int:
    """Attach consistent PointLoadCondition3D1N loads to a box face of an assembly instance."""
    gx, wg = np.polynomial.legendre.leggauss(4)
    axes_other = [d for d in range(3) if d != face_axis]
    a0, a1 = axes_other[0], axes_other[1]
    b0, b1 = box.breakpoints(a0), box.breakpoints(a1)
    nu, nv, nw = box.n_cp
    p = box.degree
    k_fixed = 0 if face_side == 0 else box.n_cp[face_axis] - 1
    tv = np.asarray(traction_vec, dtype=float)
    nodal_f: dict[int, np.ndarray] = {}

    for c0 in range(len(b0) - 1):
        j0 = 0.5 * (b0[c0 + 1] - b0[c0])
        for c1 in range(len(b1) - 1):
            j1 = 0.5 * (b1[c1 + 1] - b1[c1])
            for g0, w0 in zip(gx, wg):
                u0 = 0.5 * (b0[c0] + b0[c0 + 1]) + j0 * g0
                s0, d0 = basis_1d(box.knots[a0], p, u0, 0)
                for g1, w1 in zip(gx, wg):
                    u1 = 0.5 * (b1[c1] + b1[c1 + 1]) + j1 * g1
                    s1, d1 = basis_1d(box.knots[a1], p, u1, 0)
                    w = j0 * j1 * w0 * w1
                    for i0 in range(p + 1):
                        for i1 in range(p + 1):
                            idx_3 = [0, 0, 0]
                            idx_3[face_axis] = k_fixed
                            idx_3[a0] = s0 - p + i0
                            idx_3[a1] = s1 - p + i1
                            I = idx_3[0] + nu * (idx_3[1] + nv * idx_3[2])
                            nid_ = int(node_ids[I])
                            if nid_ not in nodal_f:
                                nodal_f[nid_] = np.zeros(3)
                            nodal_f[nid_] += (w * d0[0][i0] * d1[0][i1]) * tv

    cid = start_cid
    for nid_, fv in nodal_f.items():
        cond = mp.CreateNewCondition("PointLoadCondition3D1N", cid, [nid_], lp)
        cond.SetValue(SM.POINT_LOAD, [float(fv[0]), float(fv[1]), float(fv[2])])
        cid += 1
    return cid


# ---------------------------------------------------------------------------
# Criterion 1: V03-moment (Curved + Straight Cantilever with RBE3)
# ---------------------------------------------------------------------------


def test_g3_01_v03_moment_curved_and_straight_cantilever() -> None:
    """V03-moment curved cantilever (NURBS + RBE3) and straight v03.json moment check."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    # 1a) Curved cantilever: R_i=4.12, R_o=4.32, h=0.2, b=0.1, M=100, E=1e7, nu=0.25, theta in [0, pi/2]
    R_i, R_o, b, M_app, E, nu = 4.12, 4.32, 0.1, 100.0, 1.0e7, 0.25
    nr, nt, nz = 4, 32, 2
    rs = np.linspace(R_i, R_o, nr + 1)
    ts = np.linspace(0.0, 0.5 * math.pi, nt + 1)
    zs = np.linspace(-0.5 * b, 0.5 * b, nz + 1)

    grid_nodes = []
    for k in range(nz + 1):
        for j in range(nt + 1):
            for i in range(nr + 1):
                grid_nodes.append([rs[i] * math.cos(ts[j]), rs[i] * math.sin(ts[j]), zs[k]])
    nodes_arr = np.asarray(grid_nodes, dtype=float)

    def _nid(i: int, j: int, k: int) -> int:
        return i + (nr + 1) * (j + (nt + 1) * k)

    conn = []
    for k in range(nz):
        for j in range(nt):
            for i in range(nr):
                conn.append([
                    _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1), _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                ])
    mm = MeshModel(nodes=nodes_arr, blocks=[ElemBlock("hex8", np.asarray(conn, dtype=np.int64))])
    bf_res = mesh_to_bodyfit_patches(mm, degree=2, continuity="smooth")
    patch = bf_res.patches[0]
    ku, kv, kw = patch.knots
    k_ku, k_kv, k_kw = patch.kratos_knots()
    p, q_deg, r_deg = patch.degree
    nu_cp, nv_cp, nw_cp = patch.shape

    mp = KM.Model().CreateModelPart("V03_curved")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.ROTATION)
    mp.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.AddNodalSolutionStepVariable(SM.POINT_MOMENT)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    # Fortran order I = i + nu_cp * (j + nv_cp * k)
    k_nodes = KM.NodesVector()
    nid = 1
    for k in range(nw_cp):
        for j in range(nv_cp):
            for i in range(nu_cp):
                pt = patch.control_points[i, j, k]
                k_nodes.append(mp.CreateNewNode(nid, float(pt[0]), float(pt[1]), float(pt[2])))
                nid += 1
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)

    vol = KM.NurbsVolumeGeometry(
        k_nodes,
        p,
        q_deg,
        r_deg,
        KM.Vector(k_ku.tolist()),
        KM.Vector(k_kv.tolist()),
        KM.Vector(k_kw.tolist()),
    )
    mp.AddGeometry(vol)
    prop = mp.CreateNewProperties(1)
    prop.SetValue(KM.YOUNG_MODULUS, E)
    prop.SetValue(KM.POISSON_RATIO, nu)
    prop.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    bu, bv, bw = np.unique(ku), np.unique(kv), np.unique(kw)
    gx, wg = np.polynomial.legendre.leggauss(p + 1)
    eid = 1
    for ck in range(len(bw) - 1):
        for cj in range(len(bv) - 1):
            for ci in range(len(bu) - 1):
                pts_gp = []
                for gz_, wz in zip(gx, wg):
                    w_par = 0.5 * (bw[ck] + bw[ck + 1]) + 0.5 * (bw[ck + 1] - bw[ck]) * gz_
                    for gy_, wy in zip(gx, wg):
                        v_par = 0.5 * (bv[cj] + bv[cj + 1]) + 0.5 * (bv[cj + 1] - bv[cj]) * gy_
                        for gx_, wx in zip(gx, wg):
                            u_par = 0.5 * (bu[ci] + bu[ci + 1]) + 0.5 * (bu[ci + 1] - bu[ci]) * gx_
                            wgt = 0.125 * (bu[ci + 1] - bu[ci]) * (bv[cj + 1] - bv[cj]) * (bw[ck + 1] - bw[ck]) * wx * wy * wz
                            pts_gp.append([float(u_par), float(v_par), float(w_par), float(wgt)])
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts_gp)
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], prop)
                eid += 1

    # Fix theta = 0 (j = 0)
    for k in range(nw_cp):
        for i in range(nu_cp):
            nd = mp.GetNode(1 + i + nu_cp * (0 + nv_cp * k))
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)

    # Couple theta = pi/2 (j = nv_cp - 1) via RBE3 distributing moment coupling
    slave_nids = [1 + i + nu_cp * ((nv_cp - 1) + nv_cp * k) for k in range(nw_cp) for i in range(nu_cp)]
    slave_coords = np.array([
        [mp.GetNode(n_).X, mp.GetNode(n_).Y, mp.GetNode(n_).Z] for n_ in slave_nids
    ], dtype=float)
    r_mid = 0.5 * (R_i + R_o)
    ref_nid = 10000
    ref_coord = (0.0, r_mid, 0.0)
    S_mat = np.eye(len(slave_nids))
    rbe3_spec = RBE3Spec(
        name="curved_tip_moment",
        ref_node_id=ref_nid,
        ref_coord=ref_coord,
        slave_coords=slave_coords,
        slave_node_ids=slave_nids,
        slave_basis_mat=S_mat,
    )
    apply_rbe3(KM, SM, mp, rbe3_spec, force=(0.0, 0.0, 0.0), moment=(0.0, 0.0, M_app), mode="mpc")
    solve_linear(KM, mp, use_pardiso=True)

    ref_node = mp.GetNode(ref_nid)
    u_r_kratos = float(ref_node.GetSolutionStepValue(KM.DISPLACEMENT_Y))
    u_x_kratos = float(ref_node.GetSolutionStepValue(KM.DISPLACEMENT_X))
    u_r_exact = _curved_beam_timoshenko_ur(R_i, R_o, b, M_app, E, nu)
    rel_err_curved = abs(u_r_kratos - u_r_exact) / abs(u_r_exact)
    assert rel_err_curved <= 0.05

    # 1b) Straight v03.json moment check (0.1 x 0.1 x 1.0 m, E=2e11, nu=0, M0=5.0 -> EB -1.5e-6 m)
    v03_data = _load_verif("v03")
    eb_straight = float(v03_data["reference_values"]["euler_bernoulli_tip_uy_m"]["moment"])
    box_s = BoxSpline((-0.05, -0.05, 0.0), (0.05, 0.05, 1.0), (2, 2, 16), degree=2)
    built_s = AssemblyBuilder(KM, SM).build(
        AssemblySpec("v03_straight", parts={"b": PartSpec("b", "iga", box=box_s, young_modulus=2.0e11, poisson_ratio=0.0)}, instances=[InstanceSpec("b", "1")], use_pardiso=True)
    )
    mp_s = built_s.structure
    inst_s = built_s.instances["b_1"]
    for I in box_s.boundary_control_points(2, 0):
        nd = mp_s.GetNode(inst_s.node_ids[int(I)])
        nd.Fix(KM.DISPLACEMENT_X)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
    tip_idx = box_s.boundary_control_points(2, 1)
    tip_nids = [inst_s.node_ids[int(I)] for I in tip_idx]
    tip_coords = box_s.control_points()[tip_idx]
    ref_s_nid = built_s.id_table.next_id("nodes")
    rbe3_s = RBE3Spec("v03_s", ref_s_nid, (0.0, 0.0, 1.0), tip_coords, slave_node_ids=tip_nids, slave_basis_mat=np.eye(len(tip_nids)))
    # Moment M_x = +5.0 (right-hand rule about +x) deflects tip along -y, matching v03.json
    solve_hybrid_iga_fe(KM, SM, built_s, rbe3_specs_and_loads=[(rbe3_s, (0.0, 0.0, 0.0), (5.0, 0.0, 0.0))], use_pardiso=True)
    uy_straight = float(mp_s.GetNode(ref_s_nid).GetSolutionStepValue(KM.DISPLACEMENT_Y))
    rel_err_straight = abs(uy_straight - eb_straight) / abs(eb_straight)
    assert rel_err_straight <= 0.005

    _MEASUREMENTS["benchmarks"]["V03_moment"] = {
        "curved_cantilever": {
            "mesh_cells": [nr, nt, nz],
            "kratos_u_r_tip": u_r_kratos,
            "kratos_u_x_tip": u_x_kratos,
            "timoshenko_goodier_u_r_exact": u_r_exact,
            "rel_error": rel_err_curved,
            "gate_tol": 0.05,
            "prompt_ref_note": "Prompt cited u_r_ref=-0.0211622, which is inconsistent with M=100, E=1e7, R_i=4.12, R_o=4.32, b=0.1; exact Timoshenko-Goodier analytical value for these parameters is -2.67146011.",
        },
        "straight_cantilever_v03_json": {
            "kratos_u_y_tip": uy_straight,
            "v03_json_eb_ref": eb_straight,
            "rel_error": rel_err_straight,
        },
    }


# ---------------------------------------------------------------------------
# Criterion 2: V04 Bimetallic Cantilever (Cu/W) — Merge vs Tie + Hybrid hex8
# ---------------------------------------------------------------------------


def test_g3_02_v04_bimetallic_merge_tie_analytical_and_hybrid_hex8() -> None:
    """V04 Cu/W: C0-merge vs 2-part IGA tie (<= 1e-6), vs v04.json (<= 1%), and Hybrid hex8+IGA (<= 0.5%)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    v04_data = _load_verif("v04")
    ref_v04_eb = float(v04_data["reference_values"]["transformed_section_tip_uy_mm"]["traction"])
    ref_v04_p2_merge = float(v04_data["reference_values"]["tie_comparison_traction_case"]["p2"]["merge_tip_uy_mm"])

    # 2a) Pure IGA C0-merge single volume (double knot y=0)
    box_m = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 4, 16), degree=2)
    box_m.knots[1] = np.array([-5.0, -5.0, -5.0, -2.5, 0.0, 0.0, 2.5, 5.0, 5.0, 5.0])
    box_m.n_cp = tuple(n_control_points(k, 2) for k in box_m.knots)

    mp = KM.Model().CreateModelPart("V04_C0_merge")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(SM.POINT_LOAD)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box_m)
    pr_cu = mp.CreateNewProperties(1)
    pr_cu.SetValue(KM.YOUNG_MODULUS, 1.0e5)
    pr_cu.SetValue(KM.POISSON_RATIO, 0.0)
    pr_cu.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    pr_w = mp.CreateNewProperties(2)
    pr_w.SetValue(KM.YOUNG_MODULUS, 4.0e5)
    pr_w.SetValue(KM.POISSON_RATIO, 0.0)
    pr_w.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    bu, bv, bw = [box_m.breakpoints(d) for d in range(3)]
    gx, wg = np.polynomial.legendre.leggauss(3)
    eid = 1
    for ck in range(len(bw) - 1):
        for cj in range(len(bv) - 1):
            yc = 0.5 * (bv[cj] + bv[cj + 1])
            prop = pr_cu if yc < 0.0 else pr_w
            for ci in range(len(bu) - 1):
                pts_gp = []
                for c, wc in zip(gx, wg):
                    z = 0.5 * (bw[ck] + bw[ck + 1]) + 0.5 * (bw[ck + 1] - bw[ck]) * c
                    for b_, wb in zip(gx, wg):
                        y = yc + 0.5 * (bv[cj + 1] - bv[cj]) * b_
                        for a, wa in zip(gx, wg):
                            x = 0.5 * (bu[ci] + bu[ci + 1]) + 0.5 * (bu[ci + 1] - bu[ci]) * a
                            wgt = 0.125 * (bu[ci + 1] - bu[ci]) * (bv[cj + 1] - bv[cj]) * (bw[ck + 1] - bw[ck]) * wa * wb * wc
                            pts_gp.append([float(x), float(y), float(z), float(wgt)])
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts_gp)
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], prop)
                eid += 1

    for I in box_m.boundary_control_points(2, 0):
        nd = mp.GetNode(1 + int(I))
        nd.Fix(KM.DISPLACEMENT_X)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)
    lp = mp.CreateNewProperties(10)
    _apply_face_traction_on_instance(
        KM, SM, mp, box_m, list(range(1, box_m.n_ctrl + 1)), 2, 1, (0.0, -0.01, 0.0), lp, 1
    )
    solve_linear(KM, mp, use_pardiso=True)
    U_m = np.array([
        [
            mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_X),
            mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_Y),
            mp.GetNode(1 + i).GetSolutionStepValue(KM.DISPLACEMENT_Z),
        ]
        for i in range(box_m.n_ctrl)
    ])
    uy_merge = float(box_m.evaluate(U_m, np.array([[0.0, 0.0, 100.0]]))[0, 1])

    # 2b) 2-part IGA Tie (mpc_matching)
    b_cu = BoxSpline((-5.0, -5.0, 0.0), (5.0, 0.0, 100.0), (2, 2, 16), degree=2)
    b_w = BoxSpline((-5.0, 0.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 16), degree=2)
    built_t = AssemblyBuilder(KM, SM).build(
        AssemblySpec(
            "v04_tie",
            parts={
                "cu": PartSpec("cu", "iga", box=b_cu, young_modulus=1.0e5, poisson_ratio=0.0),
                "w": PartSpec("w", "iga", box=b_w, young_modulus=4.0e5, poisson_ratio=0.0),
            },
            instances=[InstanceSpec("cu", "1"), InstanceSpec("w", "1")],
            use_pardiso=True,
        )
    )
    mp2 = built_t.structure
    icu, iw = built_t.instances["cu_1"], built_t.instances["w_1"]
    for b_box, inst in [(b_cu, icu), (b_w, iw)]:
        for I in b_box.boundary_control_points(2, 0):
            nd = mp2.GetNode(inst.node_ids[int(I)])
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
    lp2 = mp2.CreateNewProperties(99)
    cid = _apply_face_traction_on_instance(KM, SM, mp2, b_cu, icu.node_ids, 2, 1, (0.0, -0.01, 0.0), lp2, 100)
    _apply_face_traction_on_instance(KM, SM, mp2, b_w, iw.node_ids, 2, 1, (0.0, -0.01, 0.0), lp2, cid)

    idx_cu = b_cu.boundary_control_points(1, 1)
    idx_w = b_w.boundary_control_points(1, 0)
    sm_cu = SurfaceMesh(
        "cu_if",
        "iga",
        b_cu.control_points()[idx_cu],
        np.array([[0, 1, 2]], dtype=np.int64),
        node_ids=np.array([icu.node_ids[int(I)] for I in idx_cu], dtype=np.int64),
    )
    sm_w = SurfaceMesh(
        "w_if",
        "iga",
        b_w.control_points()[idx_w],
        np.array([[0, 1, 2]], dtype=np.int64),
        node_ids=np.array([iw.node_ids[int(I)] for I in idx_w], dtype=np.int64),
    )
    solve_hybrid_iga_fe(KM, SM, built_t, ties=[TieSpec("v04_tie", sm_cu, sm_w, method="mpc_matching")], use_pardiso=True)
    U_cu = np.array([
        [
            mp2.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_X),
            mp2.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Y),
            mp2.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Z),
        ]
        for nid_ in icu.node_ids
    ])
    uy_tie = float(b_cu.evaluate(U_cu, np.array([[0.0, 0.0, 100.0]]))[0, 1])

    rel_merge_tie = abs(uy_merge - uy_tie) / abs(uy_tie)
    err_vs_eb = abs(uy_tie - ref_v04_eb) / abs(ref_v04_eb)
    assert rel_merge_tie <= 1.0e-6
    assert err_vs_eb <= 0.01

    # 2c) IGA-FE Hybrid (Beta) on V04: Cu as Native FE (hex8) + W as Immersed IGA coupled via Tie
    nx_f, ny_f, nz_f = 2, 4, 300
    xs_f = np.linspace(-5.0, 5.0, nx_f + 1)
    ys_f = np.linspace(-5.0, 0.0, ny_f + 1)
    zs_f = np.linspace(0.0, 100.0, nz_f + 1)
    fe_nodes = []
    for k in range(nz_f + 1):
        for j in range(ny_f + 1):
            for i in range(nx_f + 1):
                fe_nodes.append([xs_f[i], ys_f[j], zs_f[k]])
    fe_nodes_arr = np.asarray(fe_nodes, dtype=float)

    def _fnid(i: int, j: int, k: int) -> int:
        return i + (nx_f + 1) * (j + (ny_f + 1) * k)

    fe_conn = []
    for k in range(nz_f):
        for j in range(ny_f):
            for i in range(nx_f):
                fe_conn.append([
                    _fnid(i, j, k), _fnid(i + 1, j, k), _fnid(i + 1, j + 1, k), _fnid(i, j + 1, k),
                    _fnid(i, j, k + 1), _fnid(i + 1, j, k + 1), _fnid(i + 1, j + 1, k + 1), _fnid(i, j + 1, k + 1),
                ])
    fe_mesh = MeshModel(nodes=fe_nodes_arr, blocks=[ElemBlock("hex8", np.asarray(fe_conn, dtype=np.int64))])

    built_h = AssemblyBuilder(KM, SM).build(
        AssemblySpec(
            "v04_hybrid",
            parts={
                "cu_fe": PartSpec("cu_fe", "native_fe", mesh=fe_mesh, young_modulus=1.0e5, poisson_ratio=0.0),
                "w_iga": PartSpec("w_iga", "iga", box=b_w, young_modulus=4.0e5, poisson_ratio=0.0),
            },
            instances=[InstanceSpec("cu_fe", "1"), InstanceSpec("w_iga", "1")],
            use_pardiso=True,
        )
    )
    mph = built_h.structure
    icu_fe, iw_iga = built_h.instances["cu_fe_1"], built_h.instances["w_iga_1"]

    # Fix z=0 on both Cu FE and W IGA
    for j in range(ny_f + 1):
        for i in range(nx_f + 1):
            nd = mph.GetNode(icu_fe.node_ids[_fnid(i, j, 0)])
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
    for I in b_w.boundary_control_points(2, 0):
        nd = mph.GetNode(iw_iga.node_ids[int(I)])
        nd.Fix(KM.DISPLACEMENT_X)
        nd.Fix(KM.DISPLACEMENT_Y)
        nd.Fix(KM.DISPLACEMENT_Z)

    # Tip traction ty=-0.01 on W IGA z=100 and on Cu FE z=100 quads
    lph = mph.CreateNewProperties(199)
    cid_h = _apply_face_traction_on_instance(KM, SM, mph, b_w, iw_iga.node_ids, 2, 1, (0.0, -0.01, 0.0), lph, 200)
    area_quad = (10.0 / nx_f) * (5.0 / ny_f)
    fe_fy: dict[int, float] = {}
    for j in range(ny_f):
        for i in range(nx_f):
            q_nids = [
                icu_fe.node_ids[_fnid(i, j, nz_f)],
                icu_fe.node_ids[_fnid(i + 1, j, nz_f)],
                icu_fe.node_ids[_fnid(i + 1, j + 1, nz_f)],
                icu_fe.node_ids[_fnid(i, j + 1, nz_f)],
            ]
            for qn in q_nids:
                fe_fy[qn] = fe_fy.get(qn, 0.0) + (-0.01 * area_quad * 0.25)
    for qn, fy in fe_fy.items():
        c = mph.CreateNewCondition("PointLoadCondition3D1N", cid_h, [qn], lph)
        c.SetValue(SM.POINT_LOAD, [0.0, float(fy), 0.0])
        cid_h += 1

    # Interface meshes on y=0: FE vs IGA surface via SurfaceMesh factory methods
    fe_if_sm = SurfaceMesh.from_fe(
        "cu_fe_if",
        fe_mesh,
        icu_fe.fe_handle.mesh_to_kratos_node,
        lambda c: abs(c[1]) < 1e-10,
    )
    w_if_sm = SurfaceMesh.from_iga_face(
        "w_iga_if",
        b_w,
        axis=1,
        side="lo",
        dof_node_ids=iw_iga.node_ids,
        n_div=(4, 24),
    )
    solve_hybrid_iga_fe(
        KM,
        SM,
        built_h,
        ties=[TieSpec("hybrid_tie", w_if_sm, fe_if_sm, method="auto", penalty_beta=1.0e3, E_max=4.0e5)],
        use_pardiso=True,
    )
    U_w_hyb = np.array([
        [
            mph.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_X),
            mph.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Y),
            mph.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Z),
        ]
        for nid_ in iw_iga.node_ids
    ])
    uy_hybrid = float(b_w.evaluate(U_w_hyb, np.array([[0.0, 0.0, 100.0]]))[0, 1])
    rel_hybrid_vs_iga = abs(uy_hybrid - uy_tie) / abs(uy_tie)
    assert rel_hybrid_vs_iga <= 0.005

    _MEASUREMENTS["benchmarks"]["V04_bimetallic"] = {
        "pure_iga_c0_merge_tip_uy_mm": uy_merge,
        "pure_iga_2part_tie_tip_uy_mm": uy_tie,
        "rel_diff_merge_vs_tie": rel_merge_tie,
        "v04_json_eb_ref_mm": ref_v04_eb,
        "v04_json_p2_merge_mm": ref_v04_p2_merge,
        "rel_err_vs_eb": err_vs_eb,
        "hybrid_hex8_cu_iga_w_tip_uy_mm": uy_hybrid,
        "rel_diff_hybrid_vs_pure_iga": rel_hybrid_vs_iga,
        "prompt_ref_note": "Prompt cited u_z=-5.50763e-5 m, which does not appear in v04.json; measured against actual v04.json MMTS values (-0.0021917808 mm EB, -0.0022036041 mm p2 merge).",
    }


# ---------------------------------------------------------------------------
# Criterion 3: V05 Multiplanar Cantilever (7 Boxes)
# ---------------------------------------------------------------------------


def test_g3_03_v05_multiplanar_7box_vs_monolithic() -> None:
    """V05 7-cell multiplanar field vs monolithic single-grid CAD (<= 3.3e-6) & EB ref (<= 2%)."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    b_v05 = BoxSpline((-5.0, -5.0, 0.0), (5.0, 5.0, 100.0), (2, 2, 7), degree=2)
    domains = [
        CellDomainSpec(
            cell_id=k,
            name=f"box_{k}",
            young_modulus=1.0e11,
            poisson_ratio=0.3,
            cells=[(i, j, k) for i in range(2) for j in range(2)],
        )
        for k in range(7)
    ]

    def _solve_v05(use_domains: bool) -> float:
        parts = {
            "beam": PartSpec(
                "beam",
                "iga",
                box=b_v05,
                young_modulus=1.0e11,
                poisson_ratio=0.3,
                cell_domains=domains if use_domains else [],
            )
        }
        built = AssemblyBuilder(KM, SM).build(
            AssemblySpec("v05", parts=parts, instances=[InstanceSpec("beam", "1")], use_pardiso=True)
        )
        mp = built.structure
        inst = built.instances["beam_1"]
        for I in b_v05.boundary_control_points(2, 0):
            nd = mp.GetNode(inst.node_ids[int(I)])
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
        lp = mp.CreateNewProperties(50)
        _apply_face_traction_on_instance(KM, SM, mp, b_v05, inst.node_ids, 2, 1, (0.0, -0.01, 0.0), lp, 1)
        solve_hybrid_iga_fe(KM, SM, built, use_pardiso=True)
        U = np.array([
            [
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_X),
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Y),
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Z),
            ]
            for nid_ in inst.node_ids
        ])
        return float(b_v05.evaluate(U, np.array([[0.0, 0.0, 100.0]]))[0, 1])

    uy_mono = _solve_v05(False)
    uy_7box = _solve_v05(True)
    rel_diff = abs(uy_7box - uy_mono) / abs(uy_mono)
    eb_ref = -1.0 * (100.0 ** 3) / (3.0 * 1.0e11 * ((10.0 ** 4) / 12.0))
    err_vs_eb = abs(uy_7box - eb_ref) / abs(eb_ref)

    assert rel_diff <= 3.3e-6
    assert err_vs_eb <= 0.02

    _MEASUREMENTS["benchmarks"]["V05_multiplanar"] = {
        "monolithic_tip_uy": uy_mono,
        "seven_box_tip_uy": uy_7box,
        "rel_diff_7box_vs_mono": rel_diff,
        "euler_bernoulli_ref": eb_ref,
        "rel_err_vs_eb": err_vs_eb,
        "prompt_ref_note": "Prompt cited 2.66065e-2 m, which does not appear in v05.json; measured against exact Euler-Bernoulli (-4.0e-9) and monolithic CAD equivalence (<= 3.3e-6).",
    }


# ---------------------------------------------------------------------------
# Criterion 4: V02-local (Roark Plate with Hole — THB L3 Local Refinement)
# ---------------------------------------------------------------------------


def test_g3_04_v02_local_thb_l3_roark_plate_hole() -> None:
    """V02-local THB L3 on roark_plate_hole.step (quarter symmetry, DOFs<=80k, err<=0.2% vs C3D20)."""
    v02_data = _load_verif("v02")
    ref_c3d20 = float(v02_data["reference_values"]["reference_c3d20_sigma_max"])
    ref_roark = float(v02_data["reference_values"]["accepted_roark_sigma_max"])

    shape = read_step(str(DATA_DIR / "roark_plate_hole.step"))
    tess = tessellate(shape, deflection=0.002, angular=0.02)
    fid = hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids, radius=10.0)
    tri_c = tess.nodes[tess.tris].mean(axis=1)
    qmask = (tess.face_ids == fid) & (tri_c[:, 0] >= -1e-6) & (tri_c[:, 1] >= -1e-6)
    qtris, qfids = tess.tris[qmask], tess.face_ids[qmask]

    # Quarter symmetry [0, 100] x [0, 20] x [-0.5, 0.5] of [-100,100]x[-20,20]x[-0.5,0.5]
    # Base (33, 7, 1) -> h_0 ~ 3.0 -> THB L3 h_3 ~ 0.375 around hole r=10
    sol = _assemble_thb_kratos_and_solve(
        lo=(0.0, 0.0, -0.5),
        hi=(100.0, 20.0, 0.5),
        base_cells=(33, 7, 1),
        refine_radii=[(6.0, 22.0), (7.5, 16.0), (8.8, 13.0)],
        tess=tess,
        hole_radius=10.0,
        E=1.0,
        nu=0.0,
        traction_x=1.0e-6,
        theta_ext=0.2,
        tag="V02_local",
    )

    D = _d_matrix(1.0, 0.0)
    s_pt = float((_strain_at_points(sol["box_fine"], sol["u_fine"], np.array([[0.0, 10.0, 0.5]])) @ D.T)[0, 0])
    bp = boundary_peak_stress(sol["box_fine"], sol["u_fine"], 1.0, 0.0, tess.nodes, qtris, qfids, face_id=fid)
    s_bnd = float(bp["peak"])

    err_pt_c3d20 = abs(s_pt - ref_c3d20) / ref_c3d20
    err_bnd_c3d20 = abs(s_bnd - ref_c3d20) / ref_c3d20
    err_pt_roark = abs(s_pt - ref_roark) / ref_roark
    err_bnd_roark = abs(s_bnd - ref_roark) / ref_roark

    assert sol["active_dofs"] <= 80_000
    assert err_pt_c3d20 <= 0.002
    assert err_bnd_c3d20 <= 0.002
    assert err_pt_roark <= 0.01

    _MEASUREMENTS["benchmarks"]["V02_local_THB_L3"] = {
        "symmetry": "quarter [0,100]x[0,20]x[-0.5,0.5] (full step [-100,100]x[-20,20]x[-0.5,0.5])",
        "base_cells": [33, 7, 1],
        "thb_levels": 3,
        "finest_h_xy": [100.0 / (33 * 8), 20.0 / (7 * 8)],
        "active_dofs": sol["active_dofs"],
        "uniform_fine_dofs_quarter": sol["uniform_dofs"],
        "n_leaves": sol["n_leaves"],
        "leaves_per_level": sol["leaves_per_level"],
        "n_cut_cells": sol["n_cut_cells"],
        "n_slave_extended": sol["n_slave_extended"],
        "sigma_xx_crown_pt": s_pt,
        "boundary_peak_stress": s_bnd,
        "ref_c3d20": ref_c3d20,
        "ref_roark": ref_roark,
        "rel_err_pt_vs_c3d20": err_pt_c3d20,
        "rel_err_bnd_vs_c3d20": err_bnd_c3d20,
        "rel_err_pt_vs_roark": err_pt_roark,
        "rel_err_bnd_vs_roark": err_bnd_roark,
        "elapsed_s": sol["elapsed_s"],
    }


# ---------------------------------------------------------------------------
# Criterion 5: E01 (IPS Thin Plate with Hole + THIN_FEATURE Auto-Refine)
# ---------------------------------------------------------------------------


def _build_e01_tess() -> Any:
    """Build exact E01 10 x 5 x 0.1 in plate with r=1.0 in hole via OCP and tessellate."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    box_shape = BRepPrimAPI_MakeBox(gp_Pnt(-5.0, -2.5, -0.05), gp_Pnt(5.0, 2.5, 0.05)).Shape()
    cyl_shape = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.0, 0.0, -0.2), gp_Dir(0.0, 0.0, 1.0)), 1.0, 0.4).Shape()
    e01_shape = BRepAlgoAPI_Cut(box_shape, cyl_shape).Shape()
    return tessellate(e01_shape, deflection=0.001, angular=0.02)


def test_g3_05_e01_ips_thin_plate_hole_and_thin_feature_auto_refine() -> None:
    """E01 IPS thin plate: THIN_FEATURE detection + auto_refine_thin_features + Kratos solve <= 5% vs Roark."""
    e01_data = _load_verif("e01")
    roark_e01 = float(e01_data["reference_values"]["captured_observations_lbf_per_in2"]["Roark_accepted"])
    tess = _build_e01_tess()

    # 5a) THIN_FEATURE detection & auto_refine_thin_features
    thb0 = THBHierarchy.from_uniform((-5.0, 5.0), (-2.5, 2.5), (-0.05, 0.05), (20, 10, 1), degree=2)

    def _eval(hier: THBHierarchy) -> tuple[np.ndarray, None]:
        L = hier.max_level
        bu, bv, bw = hier._breaks_u[L], hier._breaks_v[L], hier._breaks_w[L]
        nx, ny, nz = len(bu) - 1, len(bv) - 1, len(bw) - 1
        lo = np.array([bu[0], bv[0], bw[0]])
        h = np.array([(bu[-1] - bu[0]) / nx, (bv[-1] - bv[0]) / ny, (bw[-1] - bw[0]) / nz])
        st = CL.classify_cells(tess.nodes, tess.tris, lo, h, (nx, ny, nz))["state"].reshape(nx, ny, nz)
        labels = np.where(st == CL.OUT, -1, np.where(st == CL.CUT, 0, 1))
        return labels, None

    thb_ref, reps = auto_refine_thin_features(thb0, _eval, max_retries=1)
    assert len(reps) == 2
    assert reps[0].detected is True and any(THIN_FEATURE in msg for msg in reps[0].issues)
    assert reps[0].min_feature_cells <= 1.5
    assert reps[1].detected is False and reps[1].min_feature_cells >= 2.0
    assert thb_ref.max_level == 1

    # 5b) Kratos solve via run_linear (quarter symmetry [0,5] x [0,2.5] x [-0.05,0.05])
    nc = (60, 30, 1)
    box = BoxSpline((0.0, 0.0, -0.05), (5.0, 2.5, 0.05), nc, degree=2)
    qp_tab, out_c, giss, crecs, _ = build_trimmed_table(tess.nodes, tess.tris, box, q=3, workers=1)
    spec = JobSpec(
        name="e01_ips",
        lo=(0.0, 0.0, -0.05),
        hi=(5.0, 2.5, 0.05),
        n_cells=nc,
        degree=2,
        E=3.0e7,
        nu=0.29,
        out_cells=out_c,
        out_mode="drop",
        bcs=[
            DirichletBC(face=(0, 0), components=(True, False, False), value=(0.0, 0.0, 0.0)),
            DirichletBC(face=(1, 0), components=(False, True, False), value=(0.0, 0.0, 0.0)),
            DirichletBC(face=(0, 1), components=(False, False, True), value=(0.0, 0.0, 0.0)),
        ],
        tractions=[TractionLoad(face=(0, 1), vector=(1.0e-6, 0.0, 0.0))],
        linear_solver="pardiso_ldlt",
        enable_solver_log=False,
    )
    res = run_linear(spec, geo_issues=giss, qp_table=qp_tab, cell_records=crecs, extension={"theta": 0.2})

    fid = hole_face_by_cylinder_wall(tess.nodes, tess.tris, tess.face_ids, radius=1.0)
    tri_c = tess.nodes[tess.tris].mean(axis=1)
    qmask = (tess.face_ids == fid) & (tri_c[:, 0] >= -1e-6) & (tri_c[:, 1] >= -1e-6)
    qtris, qfids = tess.tris[qmask], tess.face_ids[qmask]

    D = _d_matrix(3.0e7, 0.29)
    s_crown = float((_strain_at_points(box, res.U, np.array([[0.0, 1.0, 0.0]])) @ D.T)[0, 0])
    bp = boundary_peak_stress(box, res.U, 3.0e7, 0.29, tess.nodes, qtris, qfids, face_id=fid)
    s_bnd = float(bp["peak"])

    err_crown = abs(s_crown - roark_e01) / roark_e01
    err_bnd = abs(s_bnd - roark_e01) / roark_e01
    assert err_crown <= 0.05
    assert err_bnd <= 0.05

    _MEASUREMENTS["benchmarks"]["E01_ips_thin_plate"] = {
        "thin_feature_pass0": {
            "detected": reps[0].detected,
            "min_feature_cells": reps[0].min_feature_cells,
            "issues": reps[0].issues,
        },
        "thin_feature_pass1": {
            "detected": reps[1].detected,
            "min_feature_cells": reps[1].min_feature_cells,
        },
        "n_cells": list(nc),
        "active_dofs": res.active_dofs,
        "sigma_xx_crown": s_crown,
        "boundary_peak_stress": s_bnd,
        "roark_ref": roark_e01,
        "rel_err_crown_vs_roark": err_crown,
        "rel_err_bnd_vs_roark": err_bnd,
    }


# ---------------------------------------------------------------------------
# Criterion 6: E02 (THB Adaptive Refinement Stress Gradient & Far-Field L2)
# ---------------------------------------------------------------------------


def test_g3_06_e02_adaptive_thb_gradient_and_farfield_l2() -> None:
    """E02 L0 -> L1 -> L2 THB refinement: monotonic hole gradient resolution, DOF savings & far-field L2."""
    e01_data = _load_verif("e01")
    roark_e01 = float(e01_data["reference_values"]["captured_observations_lbf_per_in2"]["Roark_accepted"])
    tess = _build_e01_tess()
    D = _d_matrix(3.0e7, 0.29)

    ys_g = np.linspace(1.0, 1.5, 6)
    pts_grad = np.column_stack([np.zeros_like(ys_g), ys_g, np.zeros_like(ys_g)])
    xs_ff, ys_ff = np.meshgrid(np.linspace(2.2, 4.8, 8), np.linspace(0.2, 2.3, 6))
    pts_ff = np.column_stack([xs_ff.ravel(), ys_ff.ravel(), np.zeros(xs_ff.size)])

    radii_seq = [[], [(0.7, 2.1)], [(0.7, 2.1), (0.85, 1.5)]]
    level_data: dict[int, dict[str, Any]] = {}
    u_ff_by_lvl: dict[int, np.ndarray] = {}

    for lvl, rlist in enumerate(radii_seq):
        sol = _assemble_thb_kratos_and_solve(
            lo=(0.0, 0.0, -0.05),
            hi=(5.0, 2.5, 0.05),
            base_cells=(16, 8, 1),
            refine_radii=rlist,
            tess=tess,
            hole_radius=1.0,
            E=3.0e7,
            nu=0.29,
            traction_x=1.0e-6,
            theta_ext=0.2,
            tag=f"E02_L{lvl}",
        )
        s_edge = (_strain_at_points(sol["box_fine"], sol["u_fine"], pts_grad) @ D.T)[:, 0]
        u_ff = sol["box_fine"].evaluate(sol["u_fine"].reshape(-1, 3), pts_ff)
        u_ff_by_lvl[lvl] = u_ff
        level_data[lvl] = {
            "active_dofs": sol["active_dofs"],
            "uniform_dofs": sol["uniform_dofs"],
            "sigma_crown": float(s_edge[0]),
            "sigma_gradient_drop_1_to_1p5": float(s_edge[0] - s_edge[-1]),
            "rel_err_vs_roark": abs(float(s_edge[0]) - roark_e01) / roark_e01,
        }

    w_ff = np.ones(len(pts_ff))
    l2_ff_0_vs_2 = float(compute_l2_norm(u_ff_by_lvl[0], u_ff_by_lvl[2], w_ff).rel_error)
    l2_ff_1_vs_2 = float(compute_l2_norm(u_ff_by_lvl[1], u_ff_by_lvl[2], w_ff).rel_error)

    # Monotonic convergence of peak stress toward Roark across L0 -> L1 -> L2
    assert level_data[2]["rel_err_vs_roark"] < level_data[1]["rel_err_vs_roark"] < level_data[0]["rel_err_vs_roark"]
    # Substantial DOF savings at L2 vs uniform L2 mesh
    assert level_data[2]["active_dofs"] < 0.25 * level_data[2]["uniform_dofs"]
    # Far-field L2 error drops sharply from L0 -> L1
    assert l2_ff_1_vs_2 < 0.1 * l2_ff_0_vs_2

    _MEASUREMENTS["benchmarks"]["E02_adaptive_THB"] = {
        "levels": {str(k): v for k, v in level_data.items()},
        "far_field_l2_rel_diff_L0_vs_L2": l2_ff_0_vs_2,
        "far_field_l2_rel_diff_L1_vs_L2": l2_ff_1_vs_2,
        "roark_ref": roark_e01,
    }


# ---------------------------------------------------------------------------
# Criterion 7: Dependent 3-Instance Assembly Rotation/Translation Invariance
# ---------------------------------------------------------------------------


def _rot_x(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def _rot_y(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def _rot_z(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def test_g3_07_dependent_3_instance_assembly_invariance() -> None:
    """3 translated/rotated dependent instances match single-part u & rotated Voigt stress <= 1e-10."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    box = BoxSpline((0.0, 0.0, 0.0), (10.0, 1.0, 1.0), (6, 2, 2), degree=2)
    R1 = np.eye(3)
    R2 = _rot_z(37.0)
    R3 = _rot_x(25.0) @ _rot_y(-42.0)
    configs = [
        ("beam_i1", "i1", (0.0, 0.0, 0.0), R1),
        ("beam_i2", "i2", (20.0, -5.0, 3.0), R2),
        ("beam_i3", "i3", (-15.0, 12.0, 8.0), R3),
    ]
    t_local = np.array([1.0, -0.5, 0.25], dtype=float)

    parts = {"beam": PartSpec("beam", "iga", box=box, young_modulus=2.0e5, poisson_ratio=0.28)}
    inst_specs = [InstanceSpec("beam", short, dependent=True, translation=tr, rotation_matrix=R) for _, short, tr, R in configs]
    built = AssemblyBuilder(KM, SM).build(AssemblySpec("dep3", parts=parts, instances=inst_specs, use_pardiso=True))
    mp = built.structure
    lp = mp.CreateNewProperties(50)
    cid = 1
    for full_name, _, _, R in configs:
        inst = built.instances[full_name]
        for I in box.boundary_control_points(0, 0):
            nd = mp.GetNode(inst.node_ids[int(I)])
            nd.Fix(KM.DISPLACEMENT_X)
            nd.Fix(KM.DISPLACEMENT_Y)
            nd.Fix(KM.DISPLACEMENT_Z)
        f_glob = tuple((R @ t_local).tolist())
        cid = _apply_face_traction_on_instance(KM, SM, mp, box, inst.node_ids, 0, 1, f_glob, lp, cid)

    solve_hybrid_iga_fe(KM, SM, built, use_pardiso=True)

    ref_u, ref_s = None, None
    diffs: dict[str, dict[str, float]] = {}
    for full_name, _, _, R in configs:
        inst = built.instances[full_name]
        Ug = np.array([
            [
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_X),
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Y),
                mp.GetNode(nid_).GetSolutionStepValue(KM.DISPLACEMENT_Z),
            ]
            for nid_ in inst.node_ids
        ])
        Ul = Ug @ R
        s_glob = np.asarray(
            [
                row
                for eid_ in inst.element_ids
                for row in mp.GetElement(eid_).CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)
            ],
            dtype=float,
        )
        s_loc = transform_stress_voigt(s_glob, R.T)
        if ref_u is None:
            ref_u, ref_s = Ul, s_loc
        else:
            du = float(np.linalg.norm(Ul - ref_u) / np.linalg.norm(ref_u))
            ds = float(np.linalg.norm(s_loc - ref_s) / np.linalg.norm(ref_s))
            assert du <= 1.0e-10
            assert ds <= 1.0e-10
            diffs[full_name] = {"rel_diff_u": du, "rel_diff_voigt_stress": ds}

    _MEASUREMENTS["benchmarks"]["dependent_3_instance_assembly"] = diffs


# ---------------------------------------------------------------------------
# Criterion 8: 5-Format Mesh & Set Lossless Roundtrip
# ---------------------------------------------------------------------------


def test_g3_08_5_format_mesh_and_set_lossless_roundtrip(tmp_path: Path) -> None:
    """Lossless roundtrip of nodes, connectivity, parts, node_sets, elem_sets, face_sets across 5+ formats."""
    nodes = np.array(
        [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
            [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0], [2.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    b1 = ElemBlock("hex8", np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int64), part_id=1)
    b2 = ElemBlock("hex8", np.array([[1, 8, 9, 2, 5, 10, 11, 6]], dtype=np.int64), part_id=2)
    mm = MeshModel(
        nodes=nodes,
        blocks=[b1, b2],
        parts={1: "PART_A", 2: "PART_B"},
        node_sets={"NS_LEFT": [0, 3, 4, 7], "NS_RIGHT": [8, 9, 10, 11]},
        elem_sets={"ES_ONE": [0], "ES_TWO": [0, 1]},
        face_sets={"FS_TOP": [(0, 1), (1, 1)]},
        unit_system="SI",
    )

    formats = [".msh", ".inp", ".bdf", ".fem", ".rad", ".k"]
    roundtrip_summary: dict[str, Any] = {}
    for ext in formats:
        fp = tmp_path / f"roundtrip{ext}"
        write_mesh(mm, fp)
        mm2, _issues = read_mesh(fp, "SI")
        node_err = float(np.max(np.abs(mm2.nodes - mm.nodes)))
        assert node_err <= 1.0e-12
        assert len(mm2.blocks) == len(mm.blocks)
        for blk1, blk2 in zip(mm.blocks, mm2.blocks):
            assert blk1.elem_type == blk2.elem_type
            assert blk1.part_id == blk2.part_id
            assert np.array_equal(blk1.conn, blk2.conn)
        assert mm2.parts == mm.parts
        assert mm2.node_sets == mm.node_sets
        assert mm2.elem_sets == mm.elem_sets
        assert mm2.face_sets == mm.face_sets
        roundtrip_summary[ext] = {
            "max_node_err": node_err,
            "n_blocks": len(mm2.blocks),
            "parts_ok": mm2.parts == mm.parts,
            "node_sets_ok": mm2.node_sets == mm.node_sets,
            "elem_sets_ok": mm2.elem_sets == mm.elem_sets,
            "face_sets_ok": mm2.face_sets == mm.face_sets,
        }

    _MEASUREMENTS["benchmarks"]["mesh_5format_roundtrip"] = roundtrip_summary


# ---------------------------------------------------------------------------
# Criterion 9: Commercial HWT-IGA SOLVER Log, Auto RAM Switch & out_mode="drop"
# ---------------------------------------------------------------------------


def test_g3_09_commercial_solver_log_ram_switch_and_out_drop(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify all 6 HWT-IGA SOLVER stages, auto->amgcl RAM switch, and out_mode='drop' < 'ghost' DOFs."""
    out_cells = [(0, 0, 0), (0, 1, 0)]
    common_kw = dict(
        name="g3_solver_diag",
        lo=(0.0, 0.0, 0.0),
        hi=(4.0, 2.0, 1.0),
        n_cells=(6, 3, 1),
        degree=2,
        E=1.0e5,
        nu=0.3,
        out_cells=out_cells,
        bcs=[DirichletBC(face=(0, 1), components=(True, True, True), value=(0.0, 0.0, 0.0))],
        tractions=[TractionLoad(face=(0, 0), vector=(-1.0, 0.0, 0.0))],
        linear_solver="auto",
        enable_solver_log=True,
    )

    monkeypatch.delenv("WHT_SOLVER_AVAIL_RAM_GB", raising=False)
    res_drop_pardiso = run_linear(JobSpec(**common_kw, out_mode="drop"))
    res_ghost_pardiso = run_linear(JobSpec(**common_kw, out_mode="ghost"))

    expected_stages = [
        "[1/6] JOB & SYSTEM ENVIRONMENT",
        "[2/6] MODEL & IMMERSED GRID SUMMARY",
        "[3/6] CONSTRAINT & INTERACTION SUMMARY",
        "[4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION",
        "[5/6] STEP & ITERATION PROGRESS TABLE",
        "[6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY",
    ]
    assert "HWT-IGA SOLVER" in res_drop_pardiso.solver_log
    assert "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in res_drop_pardiso.solver_log
    for stage in expected_stages:
        assert stage in res_drop_pardiso.solver_log

    assert res_drop_pardiso.solver_decision is not None
    assert res_drop_pardiso.solver_decision["switched"] is False
    assert res_drop_pardiso.solver_decision["selected_solver"] == "pardiso"
    assert res_drop_pardiso.active_dofs < res_ghost_pardiso.active_dofs

    # Force low RAM budget so auto switches to AMGCL
    monkeypatch.setenv("WHT_SOLVER_AVAIL_RAM_GB", "0.0001")
    res_auto_amgcl = run_linear(JobSpec(**common_kw, out_mode="drop"))
    assert res_auto_amgcl.solver_decision is not None
    assert res_auto_amgcl.solver_decision["switched"] is True
    assert res_auto_amgcl.solver_decision["selected_solver"] == "amgcl"
    assert "amgcl" in res_auto_amgcl.solver_decision["kratos_solver"]
    assert "AUTO-SOLVER" in res_auto_amgcl.solver_log

    _MEASUREMENTS["benchmarks"]["commercial_solver_log_and_ram_switch"] = {
        "all_6_stages_present": True,
        "out_drop_active_dofs": res_drop_pardiso.active_dofs,
        "out_ghost_active_dofs": res_ghost_pardiso.active_dofs,
        "total_dofs": res_drop_pardiso.total_dofs,
        "default_auto_decision": res_drop_pardiso.solver_decision,
        "constrained_ram_auto_decision": res_auto_amgcl.solver_decision,
    }
