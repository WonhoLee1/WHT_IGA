# -*- coding: utf-8 -*-
"""S42: Phase 5 V07 (Kirsch infinite plate with hole, nu=0.49999) precursor.

Checks:
1. Exact Kirsch analytical reference (u_x, u_y, sigma_xx, sigma_yy, sigma_xy)
   verified against div(sigma) = 0, traction-free hole, sigma_xx(0, a) = 3*T,
   and Hooke's law from finite-difference strains.
2. Check A: Kratos built-in mixed elements on p=2 QuadraturePointGeometry
   (nu=0.49999) — subprocess-isolated so access violations do not kill pytest.
3. Check B: Body-fitted mapped quarter-plate Kirsch (0 <= x, y <= 4, r >= 1,
   t = 0.25, E = 100000, nu = 0.49999, T = 10) at h0 = 0.5 and h1 = 0.25:
   fbar (mean-dilatation, wht_iga.analyses.fbar) vs Kratos standard
   SmallDisplacementElement3D4N (volumetric locking baseline).
"""

from __future__ import annotations

import subprocess
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import pytest

from wht_iga.analyses import fbar
from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import n_control_points, to_kratos, uniform_clamped
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements, solve_linear

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

E_VAL = 100000.0
NU_VAL = 0.49999
T_VAL = 10.0
A_VAL = 1.0
B_VAL = 4.0
T_THICK = 0.25


# ---------------------------------------------------------------------------
# 1. Exact Kirsch analytical solution (plane strain, tensile-positive)
# ---------------------------------------------------------------------------

def kirsch_displacement(
    x: np.ndarray,
    y: np.ndarray,
    a: float = A_VAL,
    E: float = E_VAL,
    nu: float = NU_VAL,
    T: float = T_VAL,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact plane-strain Kirsch displacement (u_x, u_y)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mu = E / (2.0 * (1.0 + nu))
    kappa = 3.0 - 4.0 * nu
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    q = a / r
    q3 = q ** 3
    fac = T * a / (8.0 * mu)
    ux = fac * (
        (r / a) * (kappa + 1.0) * np.cos(theta)
        + 2.0 * q * ((kappa + 1.0) * np.cos(theta) + np.cos(3.0 * theta))
        - 2.0 * q3 * np.cos(3.0 * theta)
    )
    uy = fac * (
        (r / a) * (kappa - 3.0) * np.sin(theta)
        + 2.0 * q * ((1.0 - kappa) * np.sin(theta) + np.sin(3.0 * theta))
        - 2.0 * q3 * np.sin(3.0 * theta)
    )
    return ux, uy


def kirsch_stress(
    x: np.ndarray,
    y: np.ndarray,
    a: float = A_VAL,
    T: float = T_VAL,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact Kirsch in-plane Cauchy stresses (sigma_xx, sigma_yy, sigma_xy)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    q2 = (a / r) ** 2
    q4 = q2 * q2
    c2 = np.cos(2.0 * theta)
    s2 = np.sin(2.0 * theta)
    c4 = np.cos(4.0 * theta)
    s4 = np.sin(4.0 * theta)
    sxx = T * (1.0 - q2 * (1.5 * c2 + c4) + 1.5 * q4 * c4)
    syy = T * (-q2 * (0.5 * c2 - c4) - 1.5 * q4 * c4)
    sxy = T * (-q2 * (0.5 * s2 + s4) + 1.5 * q4 * s4)
    return sxx, syy, sxy


def test_kirsch_analytical_verification():
    """Verify Kirsch reference: hole traction-free, peak stress, div(sigma)=0, Hooke."""
    # (a) Traction-free hole r = a
    thetas = np.linspace(0.05, 0.5 * np.pi - 0.05, 19)
    xa, ya = A_VAL * np.cos(thetas), A_VAL * np.sin(thetas)
    sxx, syy, sxy = kirsch_stress(xa, ya)
    c, s = np.cos(thetas), np.sin(thetas)
    sigma_rr = sxx * c**2 + syy * s**2 + 2.0 * sxy * c * s
    sigma_rt = (syy - sxx) * c * s + sxy * (c**2 - s**2)
    sigma_tt = sxx * s**2 + syy * c**2 - 2.0 * sxy * c * s
    assert np.max(np.abs(sigma_rr)) < 1e-13
    assert np.max(np.abs(sigma_rt)) < 1e-13
    assert np.max(np.abs(sigma_tt - T_VAL * (1.0 - 2.0 * np.cos(2.0 * thetas)))) < 1e-13

    # (b) Peak stress at (0, a)
    sxx_top, syy_top, sxy_top = kirsch_stress(np.array([0.0]), np.array([A_VAL]))
    assert abs(float(sxx_top[0]) - 3.0 * T_VAL) < 1e-13
    assert abs(float(syy_top[0])) < 1e-13
    assert abs(float(sxy_top[0])) < 1e-13

    # (c) Equilibrium div(sigma) = 0 and Hooke's law via central FD
    pts_x = np.array([1.5, 2.2, 3.1, 1.2, 2.8])
    pts_y = np.array([0.7, 1.9, 1.1, 2.6, 3.3])
    h = 1e-6

    sxx_xp, _, sxy_xp = kirsch_stress(pts_x + h, pts_y)
    sxx_xm, _, sxy_xm = kirsch_stress(pts_x - h, pts_y)
    _, syy_yp, sxy_yp = kirsch_stress(pts_x, pts_y + h)
    _, syy_ym, sxy_ym = kirsch_stress(pts_x, pts_y - h)
    div_x = (sxx_xp - sxx_xm) / (2.0 * h) + (sxy_yp - sxy_ym) / (2.0 * h)
    div_y = (sxy_xp - sxy_xm) / (2.0 * h) + (syy_yp - syy_ym) / (2.0 * h)
    assert np.max(np.abs(div_x)) < 1e-8
    assert np.max(np.abs(div_y)) < 1e-8

    # Hooke's law check (use h = 1e-5 to balance truncation vs 1e9*roundoff)
    hu = 1e-5
    ux_xp, uy_xp = kirsch_displacement(pts_x + hu, pts_y)
    ux_xm, uy_xm = kirsch_displacement(pts_x - hu, pts_y)
    ux_yp, uy_yp = kirsch_displacement(pts_x, pts_y + hu)
    ux_ym, uy_ym = kirsch_displacement(pts_x, pts_y - hu)
    exx = (ux_xp - ux_xm) / (2.0 * hu)
    eyy = (uy_yp - uy_ym) / (2.0 * hu)
    exy = 0.5 * ((ux_yp - ux_ym) / (2.0 * hu) + (uy_xp - uy_xm) / (2.0 * hu))
    mu = E_VAL / (2.0 * (1.0 + NU_VAL))
    lam = E_VAL * NU_VAL / ((1.0 + NU_VAL) * (1.0 - 2.0 * NU_VAL))
    tr_e = exx + eyy
    sxx_fd = 2.0 * mu * exx + lam * tr_e
    syy_fd = 2.0 * mu * eyy + lam * tr_e
    sxy_fd = 2.0 * mu * exy
    sxx_ex, syy_ex, sxy_ex = kirsch_stress(pts_x, pts_y)
    assert np.max(np.abs(sxx_fd - sxx_ex)) < 1e-4
    assert np.max(np.abs(syy_fd - syy_ex)) < 1e-4
    assert np.max(np.abs(sxy_fd - sxy_ex)) < 1e-5


# ---------------------------------------------------------------------------
# 2. Check A: Kratos built-in mixed elements on p=2 QuadraturePointGeometry
# ---------------------------------------------------------------------------

_CHECK_A_VOL_CHILD = (
    "from wht_iga.kratos_backend import load_kratos;"
    "from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume;"
    "from wht_iga.reference.box import BoxSpline;"
    "mods = load_kratos();"
    "KM = mods['KratosMultiphysics']; SM = mods['StructuralMechanicsApplication'];"
    "box = BoxSpline(lo=(0.0,0.0,0.0), hi=(1.0,1.0,1.0), n_cells=(1,1,1), degree=2);"
    "model = KM.Model(); mp = model.CreateModelPart('A_vol');"
    "mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT);"
    "mp.AddNodalSolutionStepVariable(KM.REACTION);"
    "mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3);"
    "vol, _ = make_nurbs_volume(KM, mp, box);"
    "props = mp.GetProperties()[1];"
    "props.SetValue(KM.YOUNG_MODULUS, 100000.0); props.SetValue(KM.POISSON_RATIO, 0.49999);"
    "props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw());"
    "qps = KM.GeometriesVector();"
    "vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, 0, 0, 0));"
    "el = mp.CreateNewElement('SmallDisplacementMixedVolumetricStrainElement3D4N', 1, qps[0], props);"
    "print('created', flush=True);"
    "el.Initialize(mp.ProcessInfo);"
    "print('INIT-OK', flush=True)"
)

_CHECK_A_STRAIN_CHILD = """
import numpy as np
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline
mods = load_kratos()
KM = mods['KratosMultiphysics']; SM = mods['StructuralMechanicsApplication']
box = BoxSpline(lo=(0.0,0.0,0.0), hi=(1.0,1.0,0.25), n_cells=(2,2,1), degree=2)
model = KM.Model(); mp = model.CreateModelPart('A_strain')
mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
mp.AddNodalSolutionStepVariable(KM.REACTION)
mp.AddNodalSolutionStepVariable(KM.STRAIN)
mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
vol, _ = make_nurbs_volume(KM, mp, box)
props = mp.GetProperties()[1]
props.SetValue(KM.YOUNG_MODULUS, 100000.0); props.SetValue(KM.POISSON_RATIO, 0.49999)
props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
eid = 1
for k in range(1):
    for j in range(2):
        for i in range(2):
            qps = KM.GeometriesVector()
            vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, j, k))
            mp.CreateNewElement('SmallDisplacementMixedStrainElement3D4N', eid, qps[0], props)
            eid += 1
print('created_all', flush=True)
scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
solver = KM.SkylineLUFactorizationSolver()
builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
strategy.SetEchoLevel(0)
strategy.Initialize()
print('INIT-OK', flush=True)
for node in mp.Nodes:
    node.Fix(KM.DISPLACEMENT_Z)
    if abs(node.X0) < 1e-12: node.Fix(KM.DISPLACEMENT_X)
    if abs(node.Y0) < 1e-12: node.Fix(KM.DISPLACEMENT_Y)
    if abs(node.X0 - 1.0) < 1e-12:
        node.Fix(KM.DISPLACEMENT_X)
        node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0, 1e-4)
mp.CloneTimeStep(1.0)
strategy.Solve()
u = np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X), n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in mp.Nodes])
print(f'SOLVE-OK max_u={np.max(np.abs(u)):.6e}', flush=True)
"""


def test_check_a_kratos_mixed_qp_compatibility():
    """Check A: Kratos mixed elements on p=2 QuadraturePointGeometry (nu=0.49999)."""
    proc_vol = subprocess.run(
        [sys.executable, "-c", _CHECK_A_VOL_CHILD],
        capture_output=True,
        text=True,
        timeout=120,
    )
    print(f"\n[Check A] SmallDisplacementMixedVolumetricStrainElement3D4N (p=2, nu=0.49999): "
          f"returncode={proc_vol.returncode}, stdout={proc_vol.stdout.strip()!r}")
    assert "created" in proc_vol.stdout
    assert proc_vol.returncode != 0, "Expected access violation / crash on QP geometry"
    assert "INIT-OK" not in proc_vol.stdout

    proc_str = subprocess.run(
        [sys.executable, "-c", _CHECK_A_STRAIN_CHILD],
        capture_output=True,
        text=True,
        timeout=120,
    )
    print(f"[Check A] SmallDisplacementMixedStrainElement3D4N (p=2, nu=0.49999): "
          f"returncode={proc_str.returncode}, stdout={proc_str.stdout.strip()!r}, "
          f"stderr_tail={proc_str.stderr.strip()[-300:]!r}")


# ---------------------------------------------------------------------------
# 3. Mapped body-fitted quarter-plate B-spline geometry
# ---------------------------------------------------------------------------

def _collocate_1d(knots: np.ndarray, p: int, params: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Solve 1D B-spline interpolation sum_b N_b(params[m]) * c_b = values[m]."""
    n_cp = n_control_points(knots, p)
    A = np.zeros((n_cp, n_cp), dtype=np.float64)
    for m, xi in enumerate(params):
        span = int(find_span(n_cp, p, float(xi), knots))
        N = ders_basis_funs(span, float(xi), p, 0, knots)[0]
        A[m, span - p : span + 1] = N
    return np.linalg.solve(A, values)


def _outer_square_xy(v: np.ndarray, b: float = B_VAL) -> np.ndarray:
    """Piecewise-linear outer boundary: x=b for v<=0.5, y=b for v>=0.5."""
    out = np.zeros((len(v), 2), dtype=np.float64)
    for m, vi in enumerate(v):
        if vi <= 0.5:
            out[m] = [b, b * (2.0 * vi)]
        else:
            out[m] = [b * (2.0 * (1.0 - vi)), b]
    return out


def _inner_circle_xy(v: np.ndarray, a: float = A_VAL) -> np.ndarray:
    """Quarter circle r=a from theta=0 (v=0) to theta=pi/2 (v=1)."""
    theta = 0.5 * np.pi * v
    return np.column_stack([a * np.cos(theta), a * np.sin(theta)])


@dataclass
class MappedQuarterPlate:
    """Single-patch p=2 B-spline quarter plate (0<=x,y<=4, r>=1, 0<=z<=t).

    Uses a C^0 double knot at v=0.5 so the outer corner (4, 4) is exact, while
    remaining C^1 everywhere else. Compatible with ``fbar.build_model`` and
    Kratos ``NurbsVolumeGeometry``.
    """

    n_u: int
    n_v_half: int
    n_w: int
    degree: int = 2
    a: float = A_VAL
    b: float = B_VAL
    t: float = T_THICK
    knots: list[np.ndarray] = field(init=False)
    n_cp: tuple[int, int, int] = field(init=False)
    n_cells: tuple[int, int, int] = field(init=False)
    _cp: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        p = int(self.degree)
        ku = uniform_clamped(0.0, 1.0, self.n_u, p, continuity=p - 1)
        br_l = np.linspace(0.0, 0.5, self.n_v_half + 1)[1:-1]
        br_r = np.linspace(0.5, 1.0, self.n_v_half + 1)[1:-1]
        kv = np.concatenate([
            np.full(p + 1, 0.0),
            br_l,
            np.full(p, 0.5),
            br_r,
            np.full(p + 1, 1.0),
        ])
        kw = uniform_clamped(0.0, self.t, self.n_w, p, continuity=p - 1)
        self.knots = [ku, kv, kw]
        self.n_cp = tuple(n_control_points(k, p) for k in self.knots)  # type: ignore[assignment]
        self.n_cells = (self.n_u, 2 * self.n_v_half, self.n_w)

        gu = greville(ku, p)
        gv = greville(kv, p)
        gw = greville(kw, p)
        cp_in = _collocate_1d(kv, p, gv, _inner_circle_xy(gv, self.a))
        cp_out = _collocate_1d(kv, p, gv, _outer_square_xy(gv, self.b))

        nu, nv, nw = self.n_cp
        cp = np.zeros((nu * nv * nw, 3), dtype=np.float64)
        for k_idx, z in enumerate(gw):
            for j_idx in range(nv):
                for i_idx, u_val in enumerate(gu):
                    xy = (1.0 - u_val) * cp_in[j_idx] + u_val * cp_out[j_idx]
                    I = i_idx + nu * (j_idx + nv * k_idx)
                    cp[I] = [xy[0], xy[1], z]
        self._cp = cp

    def control_points(self) -> np.ndarray:
        return self._cp

    def breakpoints(self, axis: int) -> np.ndarray:
        return np.unique(self.knots[axis])


def _qp_positions_and_weights(plate: MappedQuarterPlate, n_gauss: int | None = None):
    """Return physical QP coordinates (nc, nq, 3), physical weights (nc, nq),
    and parametric cell Gauss points list."""
    p = plate.degree
    ng = p + 1 if n_gauss is None else int(n_gauss)
    xg, wg = np.polynomial.legendre.leggauss(ng)
    wq = np.einsum("z,y,x->zyx", wg, wg, wg).reshape(ng**3)
    cp = plate.control_points()
    nu, nv, nw = plate.n_cp
    tx = fbar._axis_cells(plate.knots[0], p, xg)
    ty = fbar._axis_cells(plate.knots[1], p, xg)
    tz = fbar._axis_cells(plate.knots[2], p, xg)
    ax = np.arange(p + 1)
    ncn = (p + 1) ** 3
    nq = ng**3
    X_all = []
    W_all = []
    param_cells = []
    for zl, zh, fz, Nw, dNw in tz:
        for yl, yh, fy, Nv, dNv in ty:
            for xl, xh, fx, Nu, dNu in tx:
                idx = (
                    (fx + ax)[None, None, :]
                    + nu * ((fy + ax)[None, :, None] + nv * (fz + ax)[:, None, None])
                ).reshape(ncn)
                N3 = np.einsum("zk,yj,xi->zyxkji", Nw, Nv, Nu).reshape(nq, ncn)
                d3 = np.empty((nq, ncn, 3))
                d3[:, :, 0] = np.einsum("zk,yj,xi->zyxkji", Nw, Nv, dNu).reshape(nq, ncn)
                d3[:, :, 1] = np.einsum("zk,yj,xi->zyxkji", Nw, dNv, Nu).reshape(nq, ncn)
                d3[:, :, 2] = np.einsum("zk,yj,xi->zyxkji", dNw, Nv, Nu).reshape(nq, ncn)
                Xq = N3 @ cp[idx]
                Jg = np.einsum("qnc,ni->qic", d3, cp[idx])
                Wq = wq * (0.125 * (xh - xl) * (yh - yl) * (zh - zl)) * np.abs(np.linalg.det(Jg))
                X_all.append(Xq)
                W_all.append(Wq)

                # Parametric cell Gauss points for Kratos
                pts = []
                hx, hy, hz = 0.5 * (xh - xl), 0.5 * (yh - yl), 0.5 * (zh - zl)
                for gz, wz in zip(xg, wg):
                    for gy, wy in zip(xg, wg):
                        for gx, wx in zip(xg, wg):
                            pts.append([
                                xl + (gx + 1.0) * hx,
                                yl + (gy + 1.0) * hy,
                                zl + (gz + 1.0) * hz,
                                wx * wy * wz * hx * hy * hz,
                            ])
                param_cells.append(pts)
    return np.stack(X_all), np.stack(W_all), param_cells


def _eval_u_at_qps(plate: MappedQuarterPlate, u_cp: np.ndarray, n_gauss: int | None = None) -> np.ndarray:
    """Evaluate B-spline displacement field u_h at all cell QPs -> (nc, nq, 3)."""
    p = plate.degree
    ng = p + 1 if n_gauss is None else int(n_gauss)
    xg, _ = np.polynomial.legendre.leggauss(ng)
    nu, nv, nw = plate.n_cp
    tx = fbar._axis_cells(plate.knots[0], p, xg)
    ty = fbar._axis_cells(plate.knots[1], p, xg)
    tz = fbar._axis_cells(plate.knots[2], p, xg)
    ax = np.arange(p + 1)
    ncn = (p + 1) ** 3
    nq = ng**3
    uc = np.asarray(u_cp, dtype=np.float64).reshape(-1, 3)
    out = []
    for _, _, fz, Nw, _ in tz:
        for _, _, fy, Nv, _ in ty:
            for _, _, fx, Nu, _ in tx:
                idx = (
                    (fx + ax)[None, None, :]
                    + nu * ((fy + ax)[None, :, None] + nv * (fz + ax)[:, None, None])
                ).reshape(ncn)
                N3 = np.einsum("zk,yj,xi->zyxkji", Nw, Nv, Nu).reshape(nq, ncn)
                out.append(N3 @ uc[idx])
    return np.stack(out)


def _kirsch_dirichlet_dofs(plate: MappedQuarterPlate):
    """Compute fixed DOFs (symmetry + u_z=0) and prescribed DOFs (outer u=1)."""
    nu, nv, nw = plate.n_cp
    p = plate.degree
    kv = plate.knots[1]
    gv = greville(kv, p)
    outer_xy = _outer_square_xy(gv, plate.b)
    ux_b, uy_b = kirsch_displacement(outer_xy[:, 0], outer_xy[:, 1])
    cx_b = _collocate_1d(kv, p, gv, ux_b)
    cy_b = _collocate_1d(kv, p, gv, uy_b)

    fixed_set: set[int] = set()
    prescribed: dict[int, float] = {}
    for k_idx in range(nw):
        for j_idx in range(nv):
            for i_idx in range(nu):
                I = i_idx + nu * (j_idx + nv * k_idx)
                fixed_set.add(3 * I + 2)  # u_z = 0 everywhere (plane strain)
                if j_idx == 0:
                    fixed_set.add(3 * I + 1)  # y = 0 symmetry: u_y = 0
                if j_idx == nv - 1:
                    fixed_set.add(3 * I + 0)  # x = 0 symmetry: u_x = 0
                if i_idx == nu - 1:
                    if j_idx != nv - 1:
                        prescribed[3 * I + 0] = float(cx_b[j_idx])
                    if j_idx != 0:
                        prescribed[3 * I + 1] = float(cy_b[j_idx])
    fixed = np.array(sorted(fixed_set), dtype=np.int64)
    return fixed, prescribed


def _l2_errors(X_all: np.ndarray, W_all: np.ndarray, u_q: np.ndarray, sxx_q: np.ndarray):
    """Compute relative L2 errors (e_u, e_sigma_xx) against exact Kirsch."""
    x = X_all[:, :, 0]
    y = X_all[:, :, 1]
    ux_ex, uy_ex = kirsch_displacement(x, y)
    sxx_ex, _, _ = kirsch_stress(x, y)

    du2 = (u_q[:, :, 0] - ux_ex) ** 2 + (u_q[:, :, 1] - uy_ex) ** 2 + u_q[:, :, 2] ** 2
    u2 = ux_ex**2 + uy_ex**2
    e_u = float(np.sqrt(np.sum(W_all * du2) / np.sum(W_all * u2)))

    ds2 = (sxx_q - sxx_ex) ** 2
    s2 = sxx_ex**2
    e_sxx = float(np.sqrt(np.sum(W_all * ds2) / np.sum(W_all * s2)))
    return e_u, e_sxx


def _solve_kratos_disp(plate: MappedQuarterPlate, param_cells: list, fixed: np.ndarray, prescribed: dict):
    """Solve the same body-fitted plate with Kratos SmallDisplacementElement3D4N."""
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("V07_kratos_disp")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    cp = plate.control_points()
    nodes = KM.NodesVector()
    for I, (x, y, z) in enumerate(cp):
        nodes.append(mp.CreateNewNode(I + 1, float(x), float(y), float(z)))
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)

    p = plate.degree
    kk = [KM.Vector(to_kratos(k, p).tolist()) for k in plate.knots]
    vol = KM.NurbsVolumeGeometry(nodes, p, p, p, kk[0], kk[1], kk[2])
    mp.AddGeometry(vol)

    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E_VAL)
    props.SetValue(KM.POISSON_RATIO, NU_VAL)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    elements = []
    for eid, pts in enumerate(param_cells, start=1):
        qps = KM.GeometriesVector()
        vol.CreateQuadraturePointGeometries(qps, 2, pts)
        elements.append(mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props))

    # Apply Dirichlet BCs
    disp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for dof in fixed:
        nid = int(dof // 3) + 1
        comp = int(dof % 3)
        mp.GetNode(nid).Fix(disp_vars[comp])
    for dof, val in prescribed.items():
        nid = int(dof // 3) + 1
        comp = int(dof % 3)
        node = mp.GetNode(nid)
        node.Fix(disp_vars[comp])
        node.SetSolutionStepValue(disp_vars[comp], 0, float(val))

    try:
        solve_linear(KM, mp, use_pardiso=True, echo=0)
    except Exception:
        solve_linear(KM, mp, use_pardiso=False, echo=0)

    u_cp = np.zeros(3 * len(cp), dtype=np.float64)
    for I in range(len(cp)):
        node = mp.GetNode(I + 1)
        u_cp[3 * I + 0] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
        u_cp[3 * I + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
        u_cp[3 * I + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)

    sxx_cells = []
    for el in elements:
        stresses = el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)
        sxx_cells.append([float(s[0]) for s in stresses])
    return u_cp, np.asarray(sxx_cells, dtype=np.float64)


# ---------------------------------------------------------------------------
# 4. Check B: fbar vs Kratos SmallDisplacementElement3D4N at h0=0.5, h1=0.25
# ---------------------------------------------------------------------------

def test_check_b_kirsch_fbar_vs_kratos_disp():
    """Check B: body-fitted quarter-plate Kirsch at nu=0.49999 (h0=0.5, h1=0.25)."""
    v_exact = T_THICK * (B_VAL**2 - 0.25 * np.pi * A_VAL**2)
    configs = [
        ("h0=0.5", 0.5, 6, 8, 2),
        ("h1=0.25", 0.25, 12, 16, 2),
    ]
    results_fbar = []
    results_kratos = []

    for label, h_val, nu_c, nv_half, nw_c in configs:
        plate = MappedQuarterPlate(n_u=nu_c, n_v_half=nv_half, n_w=nw_c, degree=2)
        X_all, W_all, param_cells = _qp_positions_and_weights(plate)
        vol_h = float(np.sum(W_all))
        vol_err = abs(vol_h - v_exact) / v_exact
        assert vol_err <= 0.001, f"{label}: domain volume error {vol_err:.6e} > 0.1%"

        fixed, prescribed = _kirsch_dirichlet_dofs(plate)

        # (1) fbar mean-dilatation
        model = fbar.build_model(plate, degree=2, E=E_VAL, nu=NU_VAL)
        f_ext = np.zeros(model["n_dof"], dtype=np.float64)
        u_fbar = fbar.solve_static(
            model, f_ext, fixed, n_steps=1, prescribed=prescribed, linear_solver="auto"
        )
        u_q_fbar = _eval_u_at_qps(plate, u_fbar)
        s_list_fbar = fbar.cauchy_stress(model, u_fbar)
        sxx_q_fbar = np.stack([s[:, 0] for s in s_list_fbar])
        eu_f, es_f = _l2_errors(X_all, W_all, u_q_fbar, sxx_q_fbar)
        results_fbar.append((h_val, eu_f, es_f, vol_err))

        # (2) Kratos standard displacement element
        u_kr, sxx_q_kr = _solve_kratos_disp(plate, param_cells, fixed, prescribed)
        u_q_kr = _eval_u_at_qps(plate, u_kr)
        eu_k, es_k = _l2_errors(X_all, W_all, u_q_kr, sxx_q_kr)
        results_kratos.append((h_val, eu_k, es_k, vol_err))

        print(
            f"\n[Check B {label}] n_ctrl={plate.n_cp[0]*plate.n_cp[1]*plate.n_cp[2]}, "
            f"n_cells={X_all.shape[0]}, vol_err={vol_err:.3e}\n"
            f"  fbar   : e_u={eu_f:.6e}, e_sigma_xx={es_f:.6e}\n"
            f"  kratos : e_u={eu_k:.6e}, e_sigma_xx={es_k:.6e}"
        )

    # Check convergence slopes between h0=0.5 and h1=0.25 for fbar
    slope_u = np.log(results_fbar[0][1] / results_fbar[1][1]) / np.log(2.0)
    slope_s = np.log(results_fbar[0][2] / results_fbar[1][2]) / np.log(2.0)
    print(f"\n[Check B slopes h0->h1] fbar: slope_u={slope_u:.3f}, slope_sigma_xx={slope_s:.3f}")

    # Assert fbar beats standard displacement locking and meets coarse/fine G5 gates
    assert results_fbar[0][1] < 0.01 and results_fbar[0][2] < 0.2
    assert results_fbar[1][1] < 0.005 and results_fbar[1][2] < 0.1
    assert results_fbar[1][1] < results_kratos[1][1]
    assert results_fbar[1][2] < results_kratos[1][2]


# ---------------------------------------------------------------------------
# 5. Check C: V07 3-Level Hard Gate (h0=0.5, h1=0.25, h2=0.125) with fbar
# ---------------------------------------------------------------------------

def test_check_c_kirsch_fbar_3level_hard_gate():
    """Check C: V07 3-level Master Plan v2.1 hard gate for fbar at nu=0.49999."""
    v_exact = T_THICK * (B_VAL**2 - 0.25 * np.pi * A_VAL**2)
    configs = [
        ("h0=0.5", 0.5, 6, 8, 2, 0.01, 0.2),
        ("h1=0.25", 0.25, 12, 16, 2, 0.005, 0.1),
        ("h2=0.125", 0.125, 24, 32, 2, 0.002, 0.05),
    ]
    results = []

    for label, h_val, nu_c, nv_half, nw_c, eu_tol, es_tol in configs:
        t0 = time.perf_counter()
        plate = MappedQuarterPlate(n_u=nu_c, n_v_half=nv_half, n_w=nw_c, degree=2)
        X_all, W_all, _ = _qp_positions_and_weights(plate)
        vol_h = float(np.sum(W_all))
        vol_err = abs(vol_h - v_exact) / v_exact
        assert vol_err <= 0.001, f"{label}: domain volume error {vol_err:.6e} > 0.1%"

        fixed, prescribed = _kirsch_dirichlet_dofs(plate)
        model = fbar.build_model(plate, degree=2, E=E_VAL, nu=NU_VAL)
        f_ext = np.zeros(model["n_dof"], dtype=np.float64)
        u_fbar = fbar.solve_static(
            model, f_ext, fixed, n_steps=1, prescribed=prescribed, linear_solver="auto"
        )
        u_q_fbar = _eval_u_at_qps(plate, u_fbar)
        s_list_fbar = fbar.cauchy_stress(model, u_fbar)
        sxx_q_fbar = np.stack([s[:, 0] for s in s_list_fbar])
        eu_f, es_f = _l2_errors(X_all, W_all, u_q_fbar, sxx_q_fbar)
        elapsed = time.perf_counter() - t0

        n_ctrl = plate.n_cp[0] * plate.n_cp[1] * plate.n_cp[2]
        n_cells = X_all.shape[0]
        print(
            f"\n[Check C {label}] n_ctrl={n_ctrl}, n_dof={model['n_dof']}, "
            f"n_cells={n_cells}, vol_err={vol_err:.3e}, time={elapsed:.2f}s\n"
            f"  fbar: e_u={eu_f:.6e} (tol={eu_tol}), e_sigma_xx={es_f:.6e} (tol={es_tol})"
        )
        assert eu_f <= eu_tol, f"{label}: e_u={eu_f:.6e} > {eu_tol}"
        assert es_f <= es_tol, f"{label}: e_sigma_xx={es_f:.6e} > {es_tol}"
        results.append((h_val, eu_f, es_f, vol_err, n_ctrl, model["n_dof"], n_cells, elapsed))

    slope_u_01 = float(np.log(results[0][1] / results[1][1]) / np.log(2.0))
    slope_s_01 = float(np.log(results[0][2] / results[1][2]) / np.log(2.0))
    slope_u_12 = float(np.log(results[1][1] / results[2][1]) / np.log(2.0))
    slope_s_12 = float(np.log(results[1][2] / results[2][2]) / np.log(2.0))

    hs = np.array([r[0] for r in results])
    eus = np.array([r[1] for r in results])
    ess = np.array([r[2] for r in results])
    fit_slope_u = float(np.polyfit(np.log(hs), np.log(eus), 1)[0])
    fit_slope_s = float(np.polyfit(np.log(hs), np.log(ess), 1)[0])

    print(
        f"\n[Check C slopes]\n"
        f"  h0->h1: slope_u={slope_u_01:.3f}, slope_sigma_xx={slope_s_01:.3f}\n"
        f"  h1->h2: slope_u={slope_u_12:.3f}, slope_sigma_xx={slope_s_12:.3f}\n"
        f"  3-lvl OLS fit: slope_u={fit_slope_u:.3f}, slope_sigma_xx={fit_slope_s:.3f}"
    )

    assert slope_u_01 > 0.5 and slope_s_01 > 0.5
    assert slope_u_12 > 0.5 and slope_s_12 > 0.5
