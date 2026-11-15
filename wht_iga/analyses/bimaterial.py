"""Bimaterial Nitsche / Hansbo-Hansbo immersed IGA elasticity & eigenstrain analysis (P4 / V11).

Implements the two-field unfitted discretization of Hansbo & Hansbo (2002) with
symmetric Nitsche interface coupling, Annavarapu / volume-fraction weighted flux
averages, and per-region AgFEM basis extension (``wht_iga.stabilize``).

Volume stiffness matrices for each material region are assembled with Kratos
``SmallDisplacementElement3D4N`` quadrature-point elements so the structural
physics backend remains Kratos-authoritative.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
import math

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from wht_iga.analyses.thermal import eigenstrain_load, solve_with_dirichlet
from wht_iga.geometry.quadrature import _axis_pw, _gauss1d, _polyval, legendre_monomials
from wht_iga.grid.bspline import basis_1d, ders_basis_funs
from wht_iga.solver.trimmed import _VOIGT as _KR_VOIGT
from wht_iga.kratos_backend.mass import _cell_gauss, reduce_by_operator
from wht_iga.reference.assemble_elastic import isotropic_D
from wht_iga.reference.box import BoxSpline
from wht_iga.stabilize.extension import ExtensionOperator, build_extension
from wht_iga.stabilize.support import INACTIVE, MASTER, SupportClassification, classify_basis

# Kratos Voigt order [xx, yy, zz, xy, yz, xz] normal projection: sigma * n = N_mat @ sigma_voigt
def _traction_normal_matrix(n: np.ndarray) -> np.ndarray:
    """Return ``(3, 6)`` matrix mapping Voigt stress ``[xx, yy, zz, xy, yz, xz]`` to traction ``sigma @ n``."""
    nx, ny, nz = float(n[0]), float(n[1]), float(n[2])
    return np.array(
        [
            [nx, 0.0, 0.0, ny, 0.0, nz],
            [0.0, ny, 0.0, nx, nz, 0.0],
            [0.0, 0.0, nz, 0.0, ny, nx],
        ],
        dtype=float,
    )


@dataclass(frozen=True)
class BimaterialRegions:
    """Two-region partition and compatible volume/interface quadrature on a ``BoxSpline``.

    Attributes:
        box: Background grid B-spline.
        full_cells: ``(full_0, full_1)`` cell index lists strictly inside region 0 or 1.
        cut_tables: ``(cut_0, cut_1)`` maps ``cell -> (nq, 4)`` ``[x, y, z, weight]``.
        cell_volumes: ``(vol_0, vol_1)`` inside volumes per cell in flat ``i + nx*(j + ny*k)`` order.
        interface: Map ``cut_cell -> (pts_w (nq, 4), normals (nq, 3))`` with normals pointing ``0 -> 1``.
        cut_fractions: Map ``cut_cell -> eta_0`` in ``[0, 1]`` for region 0 inside the cell.
    """

    box: BoxSpline
    full_cells: tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]
    cut_tables: tuple[dict[tuple[int, int, int], np.ndarray], dict[tuple[int, int, int], np.ndarray]]
    cell_volumes: tuple[np.ndarray, np.ndarray]
    interface: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray]]
    cut_fractions: dict[tuple[int, int, int], float]


@dataclass(frozen=True)
class BimaterialResult:
    """Solution and diagnostics of a bimaterial Nitsche + per-region AgFEM solve."""

    u1: np.ndarray
    u2: np.ndarray
    u_reduced: np.ndarray
    A_reduced: sp.csr_matrix
    K_vol_reduced: sp.csr_matrix
    f_reduced: np.ndarray
    free_reduced: np.ndarray
    extensions: tuple[ExtensionOperator, ExtensionOperator]
    classifications: tuple[SupportClassification, SupportClassification]
    diagnostics: dict[str, float] = field(default_factory=dict)


def _moments_2d_quarter_disc_cell(
    xa: float, xb: float, ya: float, yb: float, radius: float, q: int
) -> np.ndarray:
    """Exact 2D moments ``M[a, b] = int_{C_xy cap {r<=R}} P_a(x) P_b(y) dx dy`` via divergence theorem."""
    L = legendre_monomials(q)
    Hx = _axis_pw(xa, xb - xa, L, True)
    Py = _axis_pw(ya, yb - ya, L, False)
    gq, gw = _gauss1d(max(q + 6, 8))
    M = np.zeros((q + 1, q + 1), dtype=float)
    r2 = radius * radius

    # Vertical cell edges: x = xb (n_x = +1) and x = xa (n_x = -1)
    for x_edge, sign in ((xb, 1.0), (xa, -1.0)):
        if x_edge >= radius:
            continue
        y_cap = math.sqrt(max(0.0, r2 - x_edge * x_edge))
        y_hi = min(yb, y_cap)
        if y_hi <= ya + 1e-15:
            continue
        ys = ya + (y_hi - ya) * gq
        ws = sign * (y_hi - ya) * gw
        vx = _polyval(Hx, np.full_like(ys, x_edge))
        vy = _polyval(Py, ys)
        M += np.einsum("p,ap,bp->ab", ws, vx, vy)

    # Circular arc segment inside [xa, xb] x [ya, yb]
    th_min = max(
        math.asin(min(1.0, max(0.0, ya / radius))),
        math.acos(min(1.0, max(0.0, xb / radius))),
    )
    th_max = min(
        math.acos(min(1.0, max(0.0, xa / radius))),
        math.asin(min(1.0, max(0.0, yb / radius))),
    )
    if th_max > th_min + 1e-14:
        th = th_min + (th_max - th_min) * gq
        xs = radius * np.cos(th)
        ys = radius * np.sin(th)
        # n_x * ds = cos(theta) * (R * dtheta)
        ws = (th_max - th_min) * gw * radius * np.cos(th)
        vx = _polyval(Hx, xs)
        vy = _polyval(Py, ys)
        M += np.einsum("p,ap,bp->ab", ws, vx, vy)
    return M


def _fit_2d_cell_quadrature(
    moments_2d: np.ndarray,
    xa: float,
    xb: float,
    ya: float,
    yb: float,
    q: int,
    inside_mask_fn,
    extra_pts: np.ndarray | None = None,
) -> np.ndarray:
    """Fit positive 2D moment rule ``(nq, 3)`` ``[x, y, w]`` inside ``[xa, xb] x [ya, yb]``."""
    hx, hy = xb - xa, yb - ya
    area = hx * hy
    L = legendre_monomials(q)

    orders = (max(2 * q + 2, 8), max(3 * q + 4, 14))
    best_res = float("inf")
    best_xy = np.zeros((0, 2), dtype=float)
    best_w = np.zeros(0, dtype=float)
    for ord_n in orders:
        g1, _ = _gauss1d(ord_n)
        gx, gy = np.meshgrid(xa + hx * g1, ya + hy * g1, indexing="ij")
        pts = np.column_stack([gx.reshape(-1), gy.reshape(-1)])
        if extra_pts is not None and len(extra_pts):
            pts = np.vstack([pts, extra_pts])
        pts = pts[inside_mask_fn(pts)]
        if len(pts) == 0:
            continue
        # Only columns with z=0.5 need solving; build exact 2D design matrix in cell coordinates
        # Reuse legendre_monomials for 2D orthonormal NNLS directly matching fit_moment_rule scaling
        Px = _axis_pw(xa, hx, L, False)
        Py = _axis_pw(ya, hy, L, False)
        vx = _polyval(Px, pts[:, 0])
        vy = _polyval(Py, pts[:, 1])
        A = np.einsum("ap,bp->abp", vx, vy).reshape((q + 1) ** 2, len(pts))
        rhs = moments_2d.reshape(-1)
        # Scale rows so fit tolerance is relative to box area (identical principle to fit_moment_rule)
        from scipy.optimize import nnls

        w, rnorm = nnls(A, rhs)
        keep = w > 1e-14 * area
        if keep.any():
            w = w[keep]
            pts = pts[keep]
            rnorm = float(np.linalg.norm(A[:, keep] @ w - rhs))
        else:
            w = np.zeros(0, dtype=float)
            pts = np.zeros((0, 2), dtype=float)
        if rnorm < best_res:
            best_res, best_xy, best_w = rnorm, pts, w
        if rnorm <= 1e-11 * area:
            break
    if best_res > 1e-7 * area:
        raise RuntimeError(f"2D bimaterial moment fit did not converge: residual {best_res:.3e} vs area {area:.3e}")
    return np.column_stack([best_xy, best_w])


def _extrude_2d_table(box: BoxSpline, k: int, pts_2d: np.ndarray, n_gz: int) -> np.ndarray:
    """Extrude a 2D rule ``(nq2, 3)`` across slab cell layer ``k`` in z."""
    za, zb = float(box.breakpoints(2)[k]), float(box.breakpoints(2)[k + 1])
    gz, gw = _gauss1d(n_gz)
    zs = za + (zb - za) * gz
    wz = (zb - za) * gw
    out = np.empty((len(pts_2d) * len(zs), 4), dtype=float)
    q = 0
    for z, w in zip(zs, wz, strict=True):
        out[q : q + len(pts_2d), :2] = pts_2d[:, :2]
        out[q : q + len(pts_2d), 2] = z
        out[q : q + len(pts_2d), 3] = pts_2d[:, 2] * w
        q += len(pts_2d)
    return out


def quarter_disc_bimaterial_regions(
    box: BoxSpline,
    radius: float = 1.0,
    q: int = 4,
    n_iface_gauss: int | None = None,
) -> BimaterialRegions:
    """Exact extruded quarter-disc partition ``r <= radius`` (region 0) and ``r >= radius`` (region 1).

    Computes 2D cut-cell moments analytically on the exact circular arc and extrudes along z
    using Gauss-Legendre points. Guarantees ``vol_0 + vol_1 == cell_vol`` and geometric
    compatibility between both volume rules and interface arc quadrature.
    """
    bx, by, bz = (box.breakpoints(d) for d in range(3))
    nx, ny, nz = box.n_cells
    n_gz = box.degree + 1
    ng_if = max(box.degree + 2, q + 2) if n_iface_gauss is None else int(n_iface_gauss)
    L = legendre_monomials(q)

    full_0: list[tuple[int, int, int]] = []
    full_1: list[tuple[int, int, int]] = []
    cut_0: dict[tuple[int, int, int], np.ndarray] = {}
    cut_1: dict[tuple[int, int, int], np.ndarray] = {}
    vol_0 = np.zeros(nx * ny * nz, dtype=float)
    vol_1 = np.zeros(nx * ny * nz, dtype=float)
    iface: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray]] = {}
    fracs: dict[tuple[int, int, int], float] = {}

    g_if, w_if = _gauss1d(ng_if)
    gz, gwz = _gauss1d(n_gz)

    for j in range(ny):
        ya, yb = float(by[j]), float(by[j + 1])
        for i in range(nx):
            xa, xb = float(bx[i]), float(bx[i + 1])
            area = (xb - xa) * (yb - ya)
            r2_min = xa * xa + ya * ya
            r2_max = xb * xb + yb * yb
            r2 = radius * radius

            if r2_max <= r2 + 1e-14:
                for k in range(nz):
                    c = (i, j, k)
                    cid = i + nx * (j + ny * k)
                    full_0.append(c)
                    vol_0[cid] = area * float(bz[k + 1] - bz[k])
                continue
            if r2_min >= r2 - 1e-14:
                for k in range(nz):
                    c = (i, j, k)
                    cid = i + nx * (j + ny * k)
                    full_1.append(c)
                    vol_1[cid] = area * float(bz[k + 1] - bz[k])
                continue

            m0_2d = _moments_2d_quarter_disc_cell(xa, xb, ya, yb, radius, q)
            # Full-rectangle 2D moments (orthogonal Legendre: only (0, 0) = area)
            m_rect = np.zeros((q + 1, q + 1), dtype=float)
            Px = _axis_pw(xa, xb - xa, L, False)
            Py = _axis_pw(ya, yb - ya, L, False)
            gx, gw = _gauss1d(q + 2)
            m_rect = np.einsum(
                "p,ap,r,br->ab",
                (xb - xa) * gw,
                _polyval(Px, xa + (xb - xa) * gx),
                (yb - ya) * gw,
                _polyval(Py, ya + (yb - ya) * gx),
            )
            m1_2d = m_rect - m0_2d

            th_min = max(
                math.asin(min(1.0, max(0.0, ya / radius))),
                math.acos(min(1.0, max(0.0, xb / radius))),
            )
            th_max = min(
                math.acos(min(1.0, max(0.0, xa / radius))),
                math.asin(min(1.0, max(0.0, yb / radius))),
            )
            # Add polar candidates in cut band so thin sliver corners never starve NNLS
            rs_cand = np.linspace(max(math.sqrt(r2_min), 1e-6 * radius), radius, 12)
            ths_cand = np.linspace(max(th_min, 0.0), min(th_max, 0.5 * math.pi), 12)
            rr, tt = np.meshgrid(rs_cand, ths_cand, indexing="ij")
            cand_0_extra = np.column_stack([(rr * np.cos(tt)).ravel(), (rr * np.sin(tt)).ravel()])
            in_box = (
                (cand_0_extra[:, 0] > xa + 1e-12)
                & (cand_0_extra[:, 0] < xb - 1e-12)
                & (cand_0_extra[:, 1] > ya + 1e-12)
                & (cand_0_extra[:, 1] < yb - 1e-12)
            )
            cand_0_extra = cand_0_extra[in_box]

            rs_out = np.linspace(radius, math.sqrt(r2_max), 12)
            rr2, tt2 = np.meshgrid(rs_out, ths_cand, indexing="ij")
            cand_1_extra = np.column_stack([(rr2 * np.cos(tt2)).ravel(), (rr2 * np.sin(tt2)).ravel()])
            in_box_1 = (
                (cand_1_extra[:, 0] > xa + 1e-12)
                & (cand_1_extra[:, 0] < xb - 1e-12)
                & (cand_1_extra[:, 1] > ya + 1e-12)
                & (cand_1_extra[:, 1] < yb - 1e-12)
            )
            cand_1_extra = cand_1_extra[in_box_1]

            eta0 = float(np.clip(m0_2d[0, 0] / area, 0.0, 1.0))
            pts0_2d = (
                _fit_2d_cell_quadrature(
                    m0_2d, xa, xb, ya, yb, q,
                    lambda P: (P[:, 0] ** 2 + P[:, 1] ** 2) < r2 - 1e-12,
                    extra_pts=cand_0_extra,
                )
                if m0_2d[0, 0] > 1e-13 * area
                else np.zeros((0, 3), dtype=float)
            )
            pts1_2d = (
                _fit_2d_cell_quadrature(
                    m1_2d, xa, xb, ya, yb, q,
                    lambda P: (P[:, 0] ** 2 + P[:, 1] ** 2) > r2 + 1e-12,
                    extra_pts=cand_1_extra,
                )
                if m1_2d[0, 0] > 1e-13 * area
                else np.zeros((0, 3), dtype=float)
            )

            arc_th = th_min + (th_max - th_min) * g_if
            arc_w = (th_max - th_min) * w_if * radius

            for k in range(nz):
                c = (i, j, k)
                cid = i + nx * (j + ny * k)
                za, zb = float(bz[k]), float(bz[k + 1])
                hz = zb - za
                vol_0[cid] = float(m0_2d[0, 0]) * hz
                vol_1[cid] = float(m1_2d[0, 0]) * hz
                fracs[c] = eta0
                cut_0[c] = _extrude_2d_table(box, k, pts0_2d, n_gz) if len(pts0_2d) else np.zeros((0, 4))
                cut_1[c] = _extrude_2d_table(box, k, pts1_2d, n_gz) if len(pts1_2d) else np.zeros((0, 4))
                if th_max > th_min + 1e-14:
                    zs = za + hz * gz
                    wz = hz * gwz
                    pts_w = np.empty((len(arc_th) * len(zs), 4), dtype=float)
                    nrms = np.empty((len(arc_th) * len(zs), 3), dtype=float)
                    p = 0
                    for z, wz_ in zip(zs, wz, strict=True):
                        for th_, wth in zip(arc_th, arc_w, strict=True):
                            cx, cy = math.cos(th_), math.sin(th_)
                            pts_w[p] = (radius * cx, radius * cy, z, wth * wz_)
                            nrms[p] = (cx, cy, 0.0)
                            p += 1
                    iface[c] = (pts_w, nrms)

    return BimaterialRegions(
        box=box,
        full_cells=(full_0, full_1),
        cut_tables=(cut_0, cut_1),
        cell_volumes=(vol_0, vol_1),
        interface=iface,
        cut_fractions=fracs,
    )


def planar_bimaterial_regions(
    box: BoxSpline,
    axis: int = 0,
    interface_pos: float = 0.5,
    q: int = 2,
) -> BimaterialRegions:
    """Exact planar interface partition ``x[axis] <= interface_pos`` (region 0) and ``>= interface_pos`` (region 1).

    Used for exact constant-stress / free-expansion / dissimilar-modulus patch tests and
    cut-sliver stability sweeps (``eta_0 -> 0``).
    """
    if axis not in (0, 1, 2):
        raise ValueError(f"axis must be 0, 1, or 2, got {axis}")
    br = [box.breakpoints(d) for d in range(3)]
    nx, ny, nz = box.n_cells
    n_g = max(box.degree + 1, q + 1)
    g1, w1 = _gauss1d(n_g)

    full_0: list[tuple[int, int, int]] = []
    full_1: list[tuple[int, int, int]] = []
    cut_0: dict[tuple[int, int, int], np.ndarray] = {}
    cut_1: dict[tuple[int, int, int], np.ndarray] = {}
    vol_0 = np.zeros(nx * ny * nz, dtype=float)
    vol_1 = np.zeros(nx * ny * nz, dtype=float)
    iface: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray]] = {}
    fracs: dict[tuple[int, int, int], float] = {}

    nrm_vec = np.zeros(3, dtype=float)
    nrm_vec[axis] = 1.0

    def _subbox_table(b0: tuple[float, float, float], b1: tuple[float, float, float]) -> np.ndarray:
        axes_pts = [b0[d] + (b1[d] - b0[d]) * g1 for d in range(3)]
        axes_w = [(b1[d] - b0[d]) * w1 for d in range(3)]
        xs, ys, zs = np.meshgrid(*axes_pts, indexing="ij")
        wx, wy, wz = np.meshgrid(*axes_w, indexing="ij")
        return np.column_stack([xs.ravel(), ys.ravel(), zs.ravel(), (wx * wy * wz).ravel()])

    tang = [d for d in range(3) if d != axis]
    t1, t2 = tang[0], tang[1]

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = (i, j, k)
                cid = i + nx * (j + ny * k)
                lo = [float(br[0][i]), float(br[1][j]), float(br[2][k])]
                hi = [float(br[0][i + 1]), float(br[1][j + 1]), float(br[2][k + 1])]
                cvol = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
                if hi[axis] <= interface_pos + 1e-13:
                    full_0.append(c)
                    vol_0[cid] = cvol
                elif lo[axis] >= interface_pos - 1e-13:
                    full_1.append(c)
                    vol_1[cid] = cvol
                else:
                    hi0 = list(hi)
                    hi0[axis] = interface_pos
                    lo1 = list(lo)
                    lo1[axis] = interface_pos
                    t0 = _subbox_table(tuple(lo), tuple(hi0))
                    t1_tab = _subbox_table(tuple(lo1), tuple(hi))
                    cut_0[c] = t0
                    cut_1[c] = t1_tab
                    vol_0[cid] = float(t0[:, 3].sum())
                    vol_1[cid] = float(t1_tab[:, 3].sum())
                    fracs[c] = vol_0[cid] / cvol
                    s1 = lo[t1] + (hi[t1] - lo[t1]) * g1
                    s2 = lo[t2] + (hi[t2] - lo[t2]) * g1
                    ws1 = (hi[t1] - lo[t1]) * w1
                    ws2 = (hi[t2] - lo[t2]) * w1
                    pts_w = np.empty((len(s1) * len(s2), 4), dtype=float)
                    qidx = 0
                    for a_pos, wa in zip(s1, ws1, strict=True):
                        for b_pos, wb in zip(s2, ws2, strict=True):
                            pt = [0.0, 0.0, 0.0]
                            pt[axis] = interface_pos
                            pt[t1] = float(a_pos)
                            pt[t2] = float(b_pos)
                            pts_w[qidx] = (pt[0], pt[1], pt[2], float(wa * wb))
                            qidx += 1
                    iface[c] = (pts_w, np.tile(nrm_vec, (len(pts_w), 1)))

    return BimaterialRegions(
        box=box,
        full_cells=(full_0, full_1),
        cut_tables=(cut_0, cut_1),
        cell_volumes=(vol_0, vol_1),
        interface=iface,
        cut_fractions=fracs,
    )


def kratos_region_stiffness(
    box: BoxSpline,
    E: float,
    nu: float,
    cell_tables: dict[tuple[int, int, int], np.ndarray],
    full_cells: Sequence[tuple[int, int, int]],
) -> sp.csr_matrix:
    """Assemble region elasticity stiffness using Kratos ``SmallDisplacementElement3D4N``.

    Creates QP elements ONLY on active cells (``full_cells`` with box Gauss points and
    ``cell_tables`` with cut quadrature points) and extracts Kratos element stiffnesses
    into reference DOF order ``(box.n_dofs, box.n_dofs)``.
    """
    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.geometry import make_nurbs_volume

    active_cells: list[tuple[tuple[int, int, int], np.ndarray]] = []
    for c in full_cells:
        active_cells.append((tuple(c), np.asarray(_cell_gauss(box, tuple(c), box.degree + 1), dtype=float)))
    for c, tbl in cell_tables.items():
        if len(tbl):
            active_cells.append((tuple(c), np.asarray(tbl, dtype=float)))
    if not active_cells:
        return sp.csr_matrix((box.n_dofs, box.n_dofs), dtype=float)

    mods = load_kratos()
    KM, SM = mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("BimaterialRegion")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    volume, _nodes = make_nurbs_volume(KM, mp, box)
    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.DENSITY, 1.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    lhs = KM.Matrix()
    rhs = KM.Vector()
    for eid, (_c, pts) in enumerate(active_cells, start=1):
        qps = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(qps, 2, pts.tolist())
        assert len(qps) == 1
        el = mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
        el.Initialize(mp.ProcessInfo)
        el.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        ke = np.array(lhs, dtype=float, copy=True)
        cp_ids = np.array([n.Id - 1 for n in el.GetGeometry()], dtype=np.int64)
        edofs = (3 * cp_ids[:, None] + np.arange(3, dtype=np.int64)[None, :]).reshape(-1)
        m = edofs.size
        rows.append(np.repeat(edofs, m))
        cols.append(np.tile(edofs, m))
        vals.append(ke.reshape(-1))

    K = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(box.n_dofs, box.n_dofs),
    ).tocsr()
    return sp.csr_matrix(0.5 * (K + K.T))


def _cell_shape_and_strain(
    box: BoxSpline, cell: tuple[int, int, int], pts: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate displacement matrix ``N_u (nq, 3, 3*n_loc)``, Voigt strain matrix ``B (nq, 6, 3*n_loc)`` and DOFs."""
    p = box.degree
    nu, nv, nw = box.n_cp
    br = [box.breakpoints(d) for d in range(3)]
    spans = []
    for d in range(3):
        mid = 0.5 * (br[d][cell[d]] + br[d][cell[d] + 1])
        span, _ = basis_1d(box.knots[d], p, mid, 0)
        spans.append(int(span))
    idx = np.arange(p + 1, dtype=np.int64)
    cps = (
        (spans[0] - p + idx)[None, None, :]
        + nu * ((spans[1] - p + idx)[None, :, None] + nv * (spans[2] - p + idx)[:, None, None])
    ).reshape(-1)
    n_loc = cps.size
    edofs = (3 * cps[:, None] + np.arange(3, dtype=np.int64)[None, :]).reshape(-1)

    nq = pts.shape[0]
    N_u = np.zeros((nq, 3, 3 * n_loc), dtype=float)
    B = np.zeros((nq, 6, 3 * n_loc), dtype=float)
    k0 = np.ascontiguousarray(box.knots[0], dtype=float)
    k1 = np.ascontiguousarray(box.knots[1], dtype=float)
    k2 = np.ascontiguousarray(box.knots[2], dtype=float)
    ixx, iyy, izz = _KR_VOIGT["xx"], _KR_VOIGT["yy"], _KR_VOIGT["zz"]
    ixy, iyz, ixz = _KR_VOIGT["xy"], _KR_VOIGT["yz"], _KR_VOIGT["xz"]

    for q_idx, (x, y, z) in enumerate(pts):
        d0 = ders_basis_funs(spans[0], float(x), p, 1, k0)
        d1 = ders_basis_funs(spans[1], float(y), p, 1, k1)
        d2 = ders_basis_funs(spans[2], float(z), p, 1, k2)
        val = np.einsum("a,b,c->cba", d0[0], d1[0], d2[0]).reshape(-1)
        gx = np.einsum("a,b,c->cba", d0[1], d1[0], d2[0]).reshape(-1)
        gy = np.einsum("a,b,c->cba", d0[0], d1[1], d2[0]).reshape(-1)
        gz = np.einsum("a,b,c->cba", d0[0], d1[0], d2[1]).reshape(-1)
        for comp in range(3):
            N_u[q_idx, comp, comp::3] = val
        B[q_idx, ixx, 0::3] = gx
        B[q_idx, iyy, 1::3] = gy
        B[q_idx, izz, 2::3] = gz
        B[q_idx, ixy, 0::3] = gy
        B[q_idx, ixy, 1::3] = gx
        B[q_idx, iyz, 1::3] = gz
        B[q_idx, iyz, 2::3] = gy
        B[q_idx, ixz, 0::3] = gz
        B[q_idx, ixz, 2::3] = gx
    return N_u, B, edofs


def nitsche_interface_coupling(
    regions: BimaterialRegions,
    E1: float,
    nu1: float,
    E2: float,
    nu2: float,
    C_gamma: float | None = None,
    kappa_mode: str = "volume_fraction",
    eigenstrain1: np.ndarray | None = None,
    eigenstrain2: np.ndarray | None = None,
) -> tuple[sp.csr_matrix, sp.csr_matrix, np.ndarray, float]:
    """Assemble symmetric Nitsche consistency + penalty coupling across the bimaterial interface.

    Returns ``(A_nitsche, P_penalty, f_iface, gamma)`` on the concatenated ``2 * box.n_dofs``
    space ``[u1, u2]`` where ``P_penalty`` is the SPD jump penalty contribution ``gamma J^T J``.
    """
    box = regions.box
    n_d = box.n_dofs
    p = max(box.degree, 1)
    c_gam = 10.0 * (p ** 2) if C_gamma is None else float(C_gamma)
    br = [box.breakpoints(d) for d in range(3)]
    h_in = float(min(np.min(np.diff(br[0])), np.min(np.diff(br[1]))))
    e_eff = (2.0 * E1 * E2 / (E1 + E2)) if (E1 + E2) > 0 else max(E1, E2)
    gamma = c_gam * e_eff / h_in

    D1 = isotropic_D(E1, nu1)
    D2 = isotropic_D(E2, nu2)
    sstar1 = np.zeros(6) if eigenstrain1 is None else D1 @ np.asarray(eigenstrain1, dtype=float)
    sstar2 = np.zeros(6) if eigenstrain2 is None else D2 @ np.asarray(eigenstrain2, dtype=float)

    rows_a, cols_a, vals_a = [], [], []
    rows_p, cols_p, vals_p = [], [], []
    f_iface = np.zeros(2 * n_d, dtype=float)

    for cell, (pts_w, nrms) in regions.interface.items():
        if len(pts_w) == 0:
            continue
        eta1 = float(regions.cut_fractions[cell])
        eta2 = max(0.0, 1.0 - eta1)
        if kappa_mode == "volume_fraction":
            w1 = eta1 / max(E1, 1e-30)
            w2 = eta2 / max(E2, 1e-30)
            k1 = (w1 / (w1 + w2)) if (w1 + w2) > 0.0 else 0.5
        elif kappa_mode == "half":
            k1 = 0.5
        else:
            raise ValueError(f"unknown kappa_mode {kappa_mode!r}")
        k2 = 1.0 - k1

        N_u, B, edofs = _cell_shape_and_strain(box, cell, pts_w[:, :3])
        both_dofs = np.concatenate([edofs, n_d + edofs])
        m = both_dofs.size
        A_cell = np.zeros((m, m), dtype=float)
        P_cell = np.zeros((m, m), dtype=float)
        f_cell = np.zeros(m, dtype=float)

        for q_idx in range(len(pts_w)):
            wq = float(pts_w[q_idx, 3])
            Nmat = _traction_normal_matrix(nrms[q_idx])
            Nq = N_u[q_idx]
            Bq = B[q_idx]
            J = np.hstack([Nq, -Nq])  # [[u]] = u1 - u2
            S1 = k1 * (Nmat @ D1 @ Bq)
            S2 = k2 * (Nmat @ D2 @ Bq)
            F = np.hstack([S1, S2])  # {{D eps(u) n}}
            P_q = (wq * gamma) * (J.T @ J)
            C_q = wq * (J.T @ F + F.T @ J)
            A_cell += P_q - C_q
            P_cell += P_q
            tstar = Nmat @ (k1 * sstar1 + k2 * sstar2)
            f_cell -= wq * (J.T @ tstar)

        rows_a.append(np.repeat(both_dofs, m))
        cols_a.append(np.tile(both_dofs, m))
        vals_a.append(A_cell.reshape(-1))
        rows_p.append(np.repeat(both_dofs, m))
        cols_p.append(np.tile(both_dofs, m))
        vals_p.append(P_cell.reshape(-1))
        f_iface[both_dofs] += f_cell

    if not rows_a:
        zero = sp.csr_matrix((2 * n_d, 2 * n_d), dtype=float)
        return zero, zero, f_iface, gamma

    A_nit = sp.coo_matrix(
        (np.concatenate(vals_a), (np.concatenate(rows_a), np.concatenate(cols_a))),
        shape=(2 * n_d, 2 * n_d),
    ).tocsr()
    P_pen = sp.coo_matrix(
        (np.concatenate(vals_p), (np.concatenate(rows_p), np.concatenate(cols_p))),
        shape=(2 * n_d, 2 * n_d),
    ).tocsr()
    return sp.csr_matrix(0.5 * (A_nit + A_nit.T)), sp.csr_matrix(0.5 * (P_pen + P_pen.T)), f_iface, gamma


def per_region_extension(
    box: BoxSpline, regions: BimaterialRegions, theta: float = 0.2
) -> tuple[tuple[sp.csr_matrix, sp.csr_matrix], tuple[SupportClassification, SupportClassification], tuple[ExtensionOperator, ExtensionOperator]]:
    """Build independent per-region AgFEM extension matrices ``T_0, T_1`` (material interfaces trigger extension)."""
    nx, ny, nz = box.n_cells
    bx, by, bz = (box.breakpoints(d) for d in range(3))
    hx, hy, hz = np.diff(bx), np.diff(by), np.diff(bz)
    cvol = np.einsum("i,j,k->kji", hx, hy, hz).reshape(-1)  # i + nx*(j + ny*k)

    Ts = []
    clss = []
    exts = []
    for r in (0, 1):
        cls_r = classify_basis(box.n_cp, list(box.knots), regions.cell_volumes[r], cvol, theta=float(theta))
        ext_r = build_extension(box.n_cp, cls_r.labels, box.degree, knots=list(box.knots))
        T_r = sp.kron(ext_r.E, sp.eye(3, format="csr"), format="csr")
        Ts.append(T_r)
        clss.append(cls_r)
        exts.append(ext_r)
    return (Ts[0], Ts[1]), (clss[0], clss[1]), (exts[0], exts[1])


def solve_bimaterial_thermal_stress(
    regions: BimaterialRegions,
    E1: float,
    nu1: float,
    E2: float,
    nu2: float,
    fixed_dofs_1: Sequence[int] | np.ndarray,
    fixed_dofs_2: Sequence[int] | np.ndarray,
    eigenstrain1: np.ndarray | None = None,
    eigenstrain2: np.ndarray | None = None,
    theta: float = 0.2,
    C_gamma: float | None = None,
    kappa_mode: str = "volume_fraction",
) -> BimaterialResult:
    """Solve coupled two-field bimaterial thermo-elasticity with Nitsche + per-region AgFEM."""
    box = regions.box
    n_d = box.n_dofs
    D1 = isotropic_D(E1, nu1)
    D2 = isotropic_D(E2, nu2)

    K1 = kratos_region_stiffness(box, E1, nu1, regions.cut_tables[0], regions.full_cells[0])
    K2 = kratos_region_stiffness(box, E2, nu2, regions.cut_tables[1], regions.full_cells[1])
    K_vol = sp.block_diag([K1, K2], format="csr")

    f1 = np.zeros(n_d) if eigenstrain1 is None else eigenstrain_load(box, D1, eigenstrain1, regions.cut_tables[0], regions.full_cells[0])
    f2 = np.zeros(n_d) if eigenstrain2 is None else eigenstrain_load(box, D2, eigenstrain2, regions.cut_tables[1], regions.full_cells[1])
    f_vol = np.concatenate([f1, f2])

    A_nit, _P_pen, f_iface, gamma = nitsche_interface_coupling(
        regions, E1, nu1, E2, nu2, C_gamma=C_gamma, kappa_mode=kappa_mode,
        eigenstrain1=eigenstrain1, eigenstrain2=eigenstrain2,
    )
    A_full = sp.csr_matrix(K_vol + A_nit)
    f_full = f_vol + f_iface

    (T1, T2), clss, exts = per_region_extension(box, regions, theta=theta)
    T = sp.block_diag([T1, T2], format="csr")
    A_red = reduce_by_operator(A_full, T)
    K_red = reduce_by_operator(K_vol, T)
    f_red = np.asarray(T.T @ f_full, dtype=float)

    fixed_list: list[int] = []
    offset = 0
    for ext_r, fixed_r in ((exts[0], fixed_dofs_1), (exts[1], fixed_dofs_2)):
        fixed_set = {int(v) for v in np.asarray(fixed_r, dtype=np.int64).reshape(-1)}
        for m_idx, b_id in enumerate(ext_r.master_ids):
            for comp in range(3):
                if 3 * int(b_id) + comp in fixed_set:
                    fixed_list.append(offset + 3 * m_idx + comp)
        offset += 3 * len(ext_r.master_ids)

    fixed_red = np.array(sorted(set(fixed_list)), dtype=np.int64)
    free_red = np.setdiff1d(np.arange(A_red.shape[0], dtype=np.int64), fixed_red)
    u_red = solve_with_dirichlet(A_red, f_red, fixed_red, np.zeros(fixed_red.size))
    u_full = np.asarray(T @ u_red, dtype=float)
    u1, u2 = u_full[:n_d], u_full[n_d:]

    f1_arr = np.asarray(fixed_dofs_1, dtype=np.int64).reshape(-1)
    f2_arr = np.asarray(fixed_dofs_2, dtype=np.int64).reshape(-1)
    bc_res1 = float(np.max(np.abs(u1[f1_arr]))) if f1_arr.size else 0.0
    bc_res2 = float(np.max(np.abs(u2[f2_arr]))) if f2_arr.size else 0.0
    if max(bc_res1, bc_res2) > 1e-11:
        raise RuntimeError(f"extended Dirichlet symmetry leaked: max |u_fixed| = {max(bc_res1, bc_res2):.3e}")

    return BimaterialResult(
        u1=u1,
        u2=u2,
        u_reduced=u_red,
        A_reduced=A_red,
        K_vol_reduced=K_red,
        f_reduced=f_red,
        free_reduced=free_red,
        extensions=exts,
        classifications=clss,
        diagnostics={
            "gamma": float(gamma),
            "theta": float(theta),
            "n_master_1": float(clss[0].n_master),
            "n_extended_1": float(clss[0].n_extended),
            "n_master_2": float(clss[1].n_master),
            "n_extended_2": float(clss[1].n_extended),
            "n_reduced_free": float(free_red.size),
            "bc_leak_max": max(bc_res1, bc_res2),
        },
    )


def interface_coercivity_diagnostic(
    res: BimaterialResult,
    compute_cond: bool = True,
) -> dict[str, float]:
    """Compute minimum generalized eigenvalue ``lambda_min(A_ff, K_vol_ff)`` and condition number ``kappa(A_ff)``.

    Coercivity holds iff ``lambda_min > 0`` (every reduced free mode retains positive energy after
    subtracting the Nitsche consistency terms and adding the jump penalty).
    """
    ff = res.free_reduced
    Aff = res.A_reduced[ff][:, ff].toarray()
    Kff = res.K_vol_reduced[ff][:, ff].toarray()
    evals_a = la.eigvalsh(0.5 * (Aff + Aff.T))
    evals_k = la.eigvalsh(0.5 * (Kff + Kff.T))
    # Jump modes lie in ker(K_vol) only if disconnected, so regularize tiny eps if needed; normally Kff > 0
    reg_k = 0.5 * (Kff + Kff.T) + 1e-16 * max(float(np.max(evals_k)), 1.0) * np.eye(Aff.shape[0])
    lam_gen = la.eigvalsh(0.5 * (Aff + Aff.T), reg_k)
    return {
        "lambda_min_abs": float(evals_a[0]),
        "lambda_max_abs": float(evals_a[-1]),
        "lambda_min_coercivity_ratio": float(lam_gen[0]),
        "coercive": float(evals_a[0] > 0.0 and lam_gen[0] > 0.0),
        "kappa_A": float(evals_a[-1] / evals_a[0]) if (compute_cond and evals_a[0] > 0.0) else float("inf"),
        "k_vol_min": float(evals_k[0]),
    }
