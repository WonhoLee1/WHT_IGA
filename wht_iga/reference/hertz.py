"""Hertz-Coulomb sphere-on-flat reference (V08) and the hollow-hemisphere indenter soup.

Kratos-free. Units: any consistent set; V08 is IPS (in, lbf, psi).

1. Hertz contact of a sphere on a flat, both bodies linear-elastic half-spaces with a
   paraboloidal gap ``r**2 / (2 R)`` (Johnson 1985, Contact Mechanics, sec. 4.2):
   ``1/E* = (1-nu1**2)/E1 + (1-nu2**2)/E2``, ``a = sqrt(R delta)``,
   ``F = 4/3 E* sqrt(R) delta**1.5``, ``p0 = 3F/(2 pi a**2) = 2 E* a/(pi R)``,
   ``p_mean = 2/3 p0``; Coulomb gross sliding ``F_T = mu F_N``.
2. The exact spherical gap ``R - sqrt(R**2 - r**2)`` on the same half-spaces (Segedin 1957;
   Sneddon 1965, Int. J. Eng. Sci. 3, 47): ``delta = a/2 ln((R+a)/(R-a))``,
   ``F = E*/2 [(a**2+R**2) ln((R+a)/(R-a)) - 2aR]``,
   ``p0 = E*/pi [chi2(a/R) + artanh(a/R)]``. It reduces to Hertz for ``a << R`` and is
   used to size the paraboloid model-form error at V08's ``a/R = 0.447``.
3. V08 (Coreform ``contact-indenter``): reproduction of ``v08.json`` reference values from
   its own inputs, the prescribed load history, relative RMS, grids per h level, and a
   closed, outward-wound triangle soup of the hollow lower hemisphere for the immersed
   (trimmed) pipeline.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------------------
# 1. Hertz (paraboloid) sphere on flat
# --------------------------------------------------------------------------------------


def effective_modulus(E1: float, nu1: float, E2: float = math.inf, nu2: float = 0.0) -> float:
    """Contact modulus ``E*``: ``1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2``; ``E2 = inf`` is rigid."""
    c2 = 0.0 if math.isinf(E2) else (1.0 - nu2 * nu2) / E2
    return 1.0 / ((1.0 - nu1 * nu1) / E1 + c2)


def effective_radius(R1: float, R2: float = math.inf) -> float:
    """``1/R = 1/R1 + 1/R2`` (``R2 = inf``: flat counter-face)."""
    return 1.0 / (1.0 / R1 + (0.0 if math.isinf(R2) else 1.0 / R2))


def hertz_force(delta, R: float, E_star: float):
    """Normal force ``F = 4/3 E* sqrt(R) delta^(3/2)`` (``delta <= 0`` gives 0)."""
    d = np.maximum(np.asarray(delta, float), 0.0)
    out = 4.0 / 3.0 * E_star * math.sqrt(R) * d ** 1.5
    return float(out) if out.ndim == 0 else out


def hertz_delta(F, R: float, E_star: float):
    """Approach ``delta = (9 F^2 / (16 E*^2 R))^(1/3)``."""
    f = np.asarray(F, float)
    out = (9.0 * f * f / (16.0 * E_star * E_star * R)) ** (1.0 / 3.0)
    return float(out) if out.ndim == 0 else out


def hertz_contact_radius(F, R: float, E_star: float):
    """Contact radius from force ``a = (3 F R / (4 E*))^(1/3)``."""
    f = np.asarray(F, float)
    out = (3.0 * f * R / (4.0 * E_star)) ** (1.0 / 3.0)
    return float(out) if out.ndim == 0 else out


def hertz_contact_radius_from_delta(delta, R: float):
    """Contact radius from approach ``a = sqrt(R delta)``."""
    out = np.sqrt(R * np.maximum(np.asarray(delta, float), 0.0))
    return float(out) if out.ndim == 0 else out


def hertz_max_pressure(F, R: float, E_star: float):
    """Peak pressure ``p0 = 3 F / (2 pi a^2)``."""
    a = np.asarray(hertz_contact_radius(F, R, E_star), float)
    out = 3.0 * np.asarray(F, float) / (2.0 * math.pi * a * a)
    return float(out) if out.ndim == 0 else out


def hertz_mean_pressure(F, R: float, E_star: float):
    """Mean pressure ``F / (pi a^2) = 2/3 p0``."""
    a = np.asarray(hertz_contact_radius(F, R, E_star), float)
    out = np.asarray(F, float) / (math.pi * a * a)
    return float(out) if out.ndim == 0 else out


def hertz_pressure(r, F: float, R: float, E_star: float):
    """Pressure ``p0 sqrt(1 - r^2/a^2)`` inside ``r < a``, zero outside."""
    a = hertz_contact_radius(F, R, E_star)
    p0 = hertz_max_pressure(F, R, E_star)
    rr = np.asarray(r, float)
    out = p0 * np.sqrt(np.clip(1.0 - (rr / a) ** 2, 0.0, None))
    return float(out) if out.ndim == 0 else out


def hertz_contact_stiffness(delta, R: float, E_star: float):
    """Normal contact stiffness ``dF/d delta = 2 E* a``."""
    out = 2.0 * E_star * np.asarray(hertz_contact_radius_from_delta(delta, R), float)
    return float(out) if out.ndim == 0 else out


def coulomb_sliding_force(mu: float, F_N):
    """Gross-sliding tangential force ``F_T = mu F_N``."""
    out = mu * np.asarray(F_N, float)
    return float(out) if out.ndim == 0 else out


def dundurs_beta(E1: float, nu1: float, E2: float = math.inf, nu2: float = 0.0) -> float:
    """Dundurs ``beta = 1/2 [(1-2nu1)/G1 - (1-2nu2)/G2] / [(1-nu1)/G1 + (1-nu2)/G2]``.

    ``beta = 0`` decouples normal and tangential contact (friction does not change the
    Hertz normal solution); ``|beta| <= 1/2``. ``E2 = inf`` is a rigid body 2.
    """
    G1 = E1 / (2.0 * (1.0 + nu1))
    if math.isinf(E2):
        c2a = c2b = 0.0
    else:
        G2 = E2 / (2.0 * (1.0 + nu2))
        c2a, c2b = (1.0 - 2.0 * nu2) / G2, (1.0 - nu2) / G2
    return 0.5 * ((1.0 - 2.0 * nu1) / G1 - c2a) / ((1.0 - nu1) / G1 + c2b)


def rigid_equivalent_block_modulus(E_star: float, nu_block: float) -> float:
    """Block modulus that gives the same ``E*`` against a rigid indenter: ``(1-nu^2) E*``."""
    return (1.0 - nu_block * nu_block) * E_star


# --------------------------------------------------------------------------------------
# 2. Exact spherical gap on elastic half-spaces (Segedin 1957 / Sneddon 1965)
# --------------------------------------------------------------------------------------


def _chi2(x: float) -> float:
    """Legendre chi: ``chi2(x) = (Li2(x) - Li2(-x)) / 2 = int_0^x artanh(t)/t dt``."""
    from scipy.special import spence  # spence(1 - z) = Li2(z)

    return 0.5 * (float(spence(1.0 - x)) - float(spence(1.0 + x)))


def _check_a(a: float, R: float) -> None:
    if not (0.0 < a < R):
        raise ValueError(f"exact-sphere contact radius must satisfy 0 < a < R, got a={a}, R={R}")


def sphere_exact_delta(a: float, R: float) -> float:
    """Approach for contact radius ``a``: ``delta = a/2 ln((R+a)/(R-a))``."""
    _check_a(a, R)
    return 0.5 * a * math.log((R + a) / (R - a))


def sphere_exact_force(a: float, R: float, E_star: float) -> float:
    """Normal force ``F = E*/2 [(a^2+R^2) ln((R+a)/(R-a)) - 2 a R]``."""
    _check_a(a, R)
    return 0.5 * E_star * ((a * a + R * R) * math.log((R + a) / (R - a)) - 2.0 * a * R)


def sphere_exact_p0(a: float, R: float, E_star: float) -> float:
    """Centre pressure ``p(0) = E*/pi [chi2(a/R) + artanh(a/R)]``."""
    _check_a(a, R)
    x = a / R
    return E_star / math.pi * (_chi2(x) + math.atanh(x))


def sphere_exact_contact_radius(delta: float, R: float) -> float:
    """Invert ``delta = a/2 ln((R+a)/(R-a))`` for ``a`` in ``(0, R)`` (monotone)."""
    if delta <= 0.0:
        raise ValueError(f"delta must be > 0, got {delta}")
    from scipy.optimize import brentq

    return float(brentq(lambda a: sphere_exact_delta(a, R) - delta,
                        1e-14 * R, R * (1.0 - 1e-15), xtol=1e-15 * R, rtol=1e-15))


def boussinesq_axis_displacement(F: float, E: float, nu: float, depth: float) -> float:
    """Half-space displacement on the load axis at ``depth`` under point load ``F``.

    ``w = F (1+nu)(3-2nu) / (2 pi E depth)``. Used only as an order-of-magnitude size of
    what a pinned face at ``depth`` removes from the half-space compliance.
    """
    return F * (1.0 + nu) * (3.0 - 2.0 * nu) / (2.0 * math.pi * E * depth)


# --------------------------------------------------------------------------------------
# 3. V08 contact-indenter
# --------------------------------------------------------------------------------------

#: Block extent from the Coreform Abaqus input (contact-indenter.md sec. 8.3 ``*Node``:
#: x in [-2.5, 2.5], y in [-1, 1], z in [-1.5, 1.5]; sec. 2 Table 1: 5 x 2 x 3 in).
#: ``v08.json`` stores ``block_Lx_Ly_Lz_in = [5, 2, 1]``; Lz = 1 contradicts both.
V08_BLOCK_LO = (-2.5, -1.0, -1.5)
V08_BLOCK_HI = (2.5, 1.0, 1.5)
V08_H_LEVELS = (0.5, 0.25, 0.125)
#: Prescribed indenter cut-face motion (doc sec. 4): -0.125 in down over 0 <= t <= 0.1,
#: then +4.25 in along x over 0.1 <= t <= 1 (linear amplitudes, Abaqus ``*Amplitude``).
V08_PUSH = 0.125
V08_SLIDE = 4.25
V08_T_INDENT = 0.1
V08_SLIDE_WINDOW = (0.15, 0.70)


@dataclass(frozen=True)
class V08Inputs:
    """Inputs the V08 Hertz-Coulomb reference is computed from (IPS units).

    Assumptions (doc sec. 6): both bodies linear-elastic half-spaces combined through
    ``E*``; the indenter's curvature radius is its outer radius ``R``; the block top is
    flat (``R2 = inf``); the approach is the prescribed cut-face travel minus the initial
    clearance; friction does not alter the normal solution (``dundurs_beta ~ -1e-3``).
    """

    E_block: float
    nu_block: float
    E_indenter: float
    nu_indenter: float
    R: float
    Ri: float
    g0: float
    delta: float
    mu: float
    center: tuple[float, float, float]
    block_top_y: float

    @classmethod
    def from_json(cls, data: dict) -> "V08Inputs":
        """Read the inputs of ``tests/verification/data/v08.json`` (no defaults)."""
        geo, mat, bc = data["geometry"], data["materials"], data["boundary_conditions"]
        m = re.search(r"\bmu\s*=\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", bc["contact"])
        if m is None:
            raise ValueError(f"no 'mu = <value>' in boundary_conditions.contact: {bc['contact']!r}")
        return cls(
            E_block=float(mat["block"]["E_psi"]), nu_block=float(mat["block"]["nu"]),
            E_indenter=float(mat["indenter"]["E_psi"]), nu_indenter=float(mat["indenter"]["nu"]),
            R=float(geo["indenter_outer_radius_R_in"]), Ri=float(geo["indenter_inner_radius_Ri_in"]),
            g0=float(geo["initial_clearance_g0_in"]), delta=float(data["loads"]["indentation_delta_in"]),
            mu=float(m.group(1)), center=tuple(float(v) for v in geo["indenter_initial_center"]),
            block_top_y=float(geo["block_top_surface_y"]),
        )

    @property
    def E_star(self) -> float:
        return effective_modulus(self.E_block, self.nu_block, self.E_indenter, self.nu_indenter)


def v08_reference(inp: V08Inputs) -> dict[str, float]:
    """Hertz-Coulomb values at full indentation (doc sec. 6 Table 3) from ``inp``."""
    Es = inp.E_star
    F = hertz_force(inp.delta, inp.R, Es)
    return {
        "reduced_modulus_E_star_psi": Es,
        "full_indentation_delta_in": inp.delta,
        "contact_radius_a_in": hertz_contact_radius_from_delta(inp.delta, inp.R),
        "normal_force_F_N_lbf": F,
        "gross_sliding_force_F_T_lbf": coulomb_sliding_force(inp.mu, F),
        "maximum_pressure_p0_psi": hertz_max_pressure(F, inp.R, Es),
        "mean_pressure_psi": hertz_mean_pressure(F, inp.R, Es),
    }


def v08_model_form_gaps(inp: V08Inputs) -> dict[str, float]:
    """What the Hertz reference leaves out at V08's full indentation.

    Returns ``a/R``; the exact-sphere half-space solution at the same approach and its
    relative difference to Hertz (``(exact - hertz)/hertz``); Dundurs ``beta``; and the
    Boussinesq axis displacement at the pinned block bottom (depth = block top - bottom)
    and at the indenter cut face (depth ``R``) as fractions of ``delta``. The two
    Boussinesq numbers are order-of-magnitude only: the block has free sides, and the
    indenter is a hollow hemisphere, not a half-space.
    """
    Es, R, d = inp.E_star, inp.R, inp.delta
    ref = v08_reference(inp)
    a_ex = sphere_exact_contact_radius(d, R)
    F_ex = sphere_exact_force(a_ex, R, Es)
    p0_ex = sphere_exact_p0(a_ex, R, Es)
    F_h = ref["normal_force_F_N_lbf"]
    depth_block = inp.block_top_y - V08_BLOCK_LO[1]
    return {
        "a_over_R": ref["contact_radius_a_in"] / R,
        "exact_a_in": a_ex,
        "exact_F_N_lbf": F_ex,
        "exact_p0_psi": p0_ex,
        "exact_vs_hertz_a": a_ex / ref["contact_radius_a_in"] - 1.0,
        "exact_vs_hertz_F_N": F_ex / F_h - 1.0,
        "exact_vs_hertz_p0": p0_ex / ref["maximum_pressure_p0_psi"] - 1.0,
        "dundurs_beta": dundurs_beta(inp.E_block, inp.nu_block, inp.E_indenter, inp.nu_indenter),
        "boussinesq_block_bottom_over_delta":
            boussinesq_axis_displacement(F_h, inp.E_block, inp.nu_block, depth_block) / d,
        "boussinesq_indenter_face_over_delta":
            boussinesq_axis_displacement(F_h, inp.E_indenter, inp.nu_indenter, R) / d,
    }


def v08_cut_face_displacement(t):
    """Prescribed indenter cut-face displacement ``(u_x, u_y)`` at time ``t`` in [0, 1]."""
    tt = np.asarray(t, float)
    uy = -V08_PUSH * np.clip(tt / V08_T_INDENT, 0.0, 1.0)
    ux = V08_SLIDE * np.clip((tt - V08_T_INDENT) / (1.0 - V08_T_INDENT), 0.0, 1.0)
    return ux, uy


def v08_delta(t, g0: float = 0.025):
    """Hertz approach ``delta(t) = max(-u_y(t) - g0, 0)`` (= ``max(1.25 t - 0.025, 0)``)."""
    _, uy = v08_cut_face_displacement(t)
    out = np.maximum(-uy - g0, 0.0)
    return float(out) if out.ndim == 0 else out


def v08_reference_history(t, inp: V08Inputs) -> dict[str, np.ndarray]:
    """Hertz ``F_N(t)``, ``p0(t)`` and Coulomb ``F_T(t) = mu F_N`` (gross sliding, t > 0.1)."""
    tt = np.atleast_1d(np.asarray(t, float))
    Es = inp.E_star
    d = np.asarray(v08_delta(tt, inp.g0), float)
    F = np.asarray(hertz_force(d, inp.R, Es), float)
    a = np.sqrt(inp.R * d)
    p0 = np.where(a > 0.0, 2.0 * Es * a / (math.pi * inp.R), 0.0)
    FT = np.where(tt > V08_T_INDENT, inp.mu * F, np.nan)  # reference undefined while indenting
    return {"t": tt, "delta": d, "F_N": F, "p0": p0, "F_T": FT}


def relative_rms(q, q_ref) -> float:
    """Doc sec. 7.1: ``sqrt(sum (q - q*)^2) / sqrt(sum q*^2)``."""
    q = np.asarray(q, float)
    r = np.asarray(q_ref, float)
    return float(np.sqrt(np.sum((q - r) ** 2)) / np.sqrt(np.sum(r * r)))


# --------------------------------------------------------------------------------------
# Grids per h level (p = 2, C^1 as in the Coreform IGA runs)
# --------------------------------------------------------------------------------------


def bspline_dofs(n_cells, degree: int = 2) -> int:
    """Displacement DOFs of a maximal-continuity tensor B-spline block (``3 prod(n+p)``)."""
    return 3 * int(np.prod([int(n) + degree for n in n_cells]))


def bspline_stiffness_nnz(n_cells, degree: int = 2) -> int:
    """Exact non-zeros of the 3D vector stiffness pattern (basis overlap ``|i-j| <= p``)."""
    p = degree
    out = 9
    for n in n_cells:
        m = int(n) + p
        out *= sum(min(m - 1, i + p) - max(0, i - p) + 1 for i in range(m))
    return out


#: Symmetry options: ``"full"``; ``"half"`` keeps ``z >= c_z`` (exact for the whole V08
#: history); ``"quarter"`` also keeps ``x >= c_x`` (indentation only, before sliding).
_SYMS = ("full", "half", "quarter")


def _check_sym(sym: str) -> None:
    if sym not in _SYMS:
        raise ValueError(f"sym must be one of {_SYMS}, got {sym!r}")


def v08_block_grid(h: float, sym: str = "full", x_sym: float = -1.5):
    """Body-fitted block grid ``(lo, hi, n_cells)``; symmetry planes are the lower planes."""
    _check_sym(sym)
    lo = np.array(V08_BLOCK_LO, float)
    hi = np.array(V08_BLOCK_HI, float)
    if sym != "full":
        lo[2] = 0.0
    if sym == "quarter":
        lo[0] = float(x_sym)
    n = tuple(int(round(v)) for v in (hi - lo) / h)
    if not np.allclose(np.array(n) * h, hi - lo, rtol=0.0, atol=1e-12):
        raise ValueError(f"h={h} does not divide the block {hi - lo}")
    return tuple(lo), tuple(hi), n


def v08_indenter_grid(h: float, center=(-1.5, 1.525, 0.0), R: float = 0.5,
                      sym: str = "full", margin: float | None = None):
    """Immersed indenter grid ``(lo, hi, n_cells)`` with cell size exactly ``h``.

    The top plane is the cut plane ``y = center_y`` (strong Dirichlet on the boundary
    control points moves the whole loaded face). Other sides clear the sphere by at
    least ``margin`` (default ``0.1 h``), so the contact pole is not on a grid plane.
    Symmetry planes (``sym``) are the grid's lower z (and x) planes, as for the V06 octant.
    """
    _check_sym(sym)
    c = np.asarray(center, float)
    m = 0.1 * h if margin is None else float(margin)
    n_one = math.ceil((R + m) / h - 1e-12)  # one-sided count
    n_two = math.ceil((2.0 * (R + m)) / h - 1e-12)
    nx = n_one if sym == "quarter" else n_two
    nz = n_two if sym == "full" else n_one
    ny = n_one
    lo = np.array([c[0] if sym == "quarter" else c[0] - 0.5 * nx * h,
                   c[1] - ny * h,
                   c[2] if sym != "full" else c[2] - 0.5 * nz * h])
    hi = lo + h * np.array([nx, ny, nz], float)
    hi[1] = c[1]
    return tuple(lo), tuple(hi), (nx, ny, nz)


# --------------------------------------------------------------------------------------
# Hollow lower hemisphere triangle soup
# --------------------------------------------------------------------------------------

#: face ids of :func:`hollow_hemisphere_soup`
FACE_OUTER, FACE_INNER, FACE_CUT, FACE_SYM_Z, FACE_SYM_X = 1, 2, 3, 4, 5


def _trig_table(n: int) -> tuple[np.ndarray, np.ndarray]:
    """``sin``/``cos`` of ``linspace(0, pi, n+1)`` with exact end and mid values."""
    ang = np.linspace(0.0, math.pi, n + 1)
    s, c = np.sin(ang), np.cos(ang)
    s[0] = s[-1] = 0.0
    c[0], c[-1] = 1.0, -1.0
    if n % 2 == 0:
        s[n // 2], c[n // 2] = 1.0, 0.0
    return s, c


def hollow_hemisphere_soup(R: float, Ri: float, center=(0.0, 0.0, 0.0), n: int = 24,
                           n_rad: int | None = None, sym: str = "full"
                           ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed, outward-wound triangle soup of the lower hollow hemisphere.

    Solid: ``Ri <= |x - c| <= R`` and ``y <= c_y`` (bottom half of a spherical shell cut
    through its centre by the plane ``y = c_y``, parallel to x-z; V08 geometry).
    ``sym = "half"`` keeps ``z >= c_z`` and ``"quarter"`` also ``x >= c_x``; the cut
    symmetry planes are closed by planar faces.

    Spheres use ``x = c + r (sin t cos f, -sin t sin f, cos t)`` with ``t, f`` in
    ``[0, pi]``, ``n`` steps each (``n`` even). The poles sit on the cut plane at
    ``z = c_z +- r``; the lowest point (``t = f = pi/2``, the contact point) is a regular
    vertex, not a pole fan. Vertices lie on the exact surfaces (inscribed facets) and the
    trig tables are snapped, so vertices shared by two faces have bit-identical
    coordinates and are welded by position: the soup is watertight. ``n_rad`` radial
    steps on the planar faces (default: spacing close to the outer arc step ``pi R / n``).

    Returns:
        ``(nodes (N, 3), tris (M, 3) int64, face_ids (M,) int64)`` with face ids
        ``FACE_OUTER = 1`` (contact surface), ``FACE_INNER = 2``, ``FACE_CUT = 3``
        (loaded face), ``FACE_SYM_Z = 4``, ``FACE_SYM_X = 5``; normals point out of the
        solid.
    """
    _check_sym(sym)
    if not (0.0 < Ri < R):
        raise ValueError(f"need 0 < Ri < R, got Ri={Ri}, R={R}")
    n = int(n)
    if n < 2 or n % 2:
        raise ValueError(f"n must be even and >= 2, got {n}")
    if n_rad is None:
        n_rad = max(1, math.ceil((R - Ri) / (math.pi * R / n) - 1e-9))
    n_rad = int(n_rad)
    c = np.asarray(center, float)
    s, co = _trig_table(n)
    rho = [Ri + (R - Ri) * k / n_rad for k in range(n_rad + 1)]
    rho[0], rho[-1] = Ri, R
    i_max = n // 2 if sym != "full" else n      # t range (z >= c_z for half/quarter)
    j_max = n // 2 if sym == "quarter" else n   # f range (x >= c_x for quarter)

    xyz: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}

    def vid(p: tuple[float, float, float]) -> int:
        key = (float(p[0]) + 0.0, float(p[1]) + 0.0, float(p[2]) + 0.0)  # folds -0.0
        if key not in index:
            index[key] = len(xyz)
            xyz.append(key)
        return index[key]

    tris: list[tuple[int, int, int]] = []
    fids: list[int] = []

    def patch(pt, nu: int, nv: int, flip: bool, fid: int) -> None:
        """Structured (u, v) patch; drops triangles degenerate at collapsed poles."""
        ids = [[vid(pt(u, v)) for v in range(nv + 1)] for u in range(nu + 1)]
        for u in range(nu):
            for v in range(nv):
                a, b, d, e = ids[u][v], ids[u + 1][v], ids[u + 1][v + 1], ids[u][v + 1]
                for t in ((a, b, d), (a, d, e)):
                    if len(set(t)) < 3:
                        continue
                    tris.append((t[0], t[2], t[1]) if flip else t)
                    fids.append(fid)

    def sphere(r: float):
        return lambda i, j: (c[0] + r * s[i] * co[j], c[1] - r * s[i] * s[j], c[2] + r * co[i])

    # d/dt x d/df points into the sphere: flip on the outer sphere (outward = +radial).
    patch(sphere(R), i_max, j_max, True, FACE_OUTER)
    patch(sphere(Ri), i_max, j_max, False, FACE_INNER)

    # cut plane y = c_y: rim position m runs from +z towards +x (f = 0 meridian, t = m),
    # on to -z, then back up along x < 0 (f = pi, t = 2n - m). e_rho x e_psi = +y.
    if sym == "full":
        ms = list(range(2 * n + 1))
    elif sym == "half":
        ms = list(range(3 * n // 2, 2 * n)) + list(range(0, n // 2 + 1))
    else:
        ms = list(range(0, n // 2 + 1))

    def cut(k: int, q: int):
        m = ms[q] % (2 * n)
        t, sx = (m, s[m]) if m <= n else (2 * n - m, -s[2 * n - m])
        return (c[0] + rho[k] * sx, c[1], c[2] + rho[k] * co[t])

    patch(cut, n_rad, len(ms) - 1, False, FACE_CUT)
    if sym != "full":   # z = c_z plane, outward -z: d/drho x d/df = -rho e_z
        patch(lambda k, j: (c[0] + rho[k] * co[j], c[1] - rho[k] * s[j], c[2]),
              n_rad, j_max, False, FACE_SYM_Z)
    if sym == "quarter":  # x = c_x plane, outward -x: d/drho x d/dt = +rho e_x
        patch(lambda k, i: (c[0], c[1] - rho[k] * s[i], c[2] + rho[k] * co[i]),
              n_rad, i_max, True, FACE_SYM_X)

    return (np.asarray(xyz, float), np.asarray(tris, np.int64), np.asarray(fids, np.int64))


def hollow_hemisphere_volume(R: float, Ri: float, sym: str = "full") -> float:
    """Exact solid volume ``2/3 pi (R^3 - Ri^3)`` (halved per symmetry plane)."""
    _check_sym(sym)
    return 2.0 / 3.0 * math.pi * (R ** 3 - Ri ** 3) / {"full": 1, "half": 2, "quarter": 4}[sym]


def soup_signed_volume(nodes: np.ndarray, tris: np.ndarray) -> float:
    """Divergence-theorem volume ``sum v0 . (v1 x v2) / 6`` (> 0 for outward windings)."""
    p = np.asarray(nodes, float)[np.asarray(tris, np.int64)]
    return float(np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum() / 6.0)
