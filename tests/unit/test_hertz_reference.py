"""V08 Hertz-Coulomb reference + hollow-hemisphere soup (Kratos-free).

Independent references:
- Timoshenko & Goodier, Theory of Elasticity (3rd ed.), sec. 140, nu = 0.3 coefficients
  (printed to 3-4 digits; tolerance is one unit in the last printed digit).
- Numerical Boussinesq integration of the Hertz pressure (the defining gap condition).
- Numerical Sneddon (1965) axisymmetric-punch integrals for the exact spherical gap.
- Coreform contact-indenter.md sec. 6 Table 3 and ``v08.json`` reference values.
- Divergence-theorem volume vs ``2/3 pi (R^3 - Ri^3)``; repo winding-number inside test
  and trimmed moment-fit quadrature for pipeline compatibility.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from wht_iga.reference import hertz as H

V08_JSON = Path(__file__).resolve().parents[1] / "verification" / "data" / "v08.json"
C08 = (-1.5, 1.525, 0.0)  # V08 indenter centre (doc sec. 2)
R08, RI08 = 0.5, 0.25


@pytest.fixture(scope="module")
def v08() -> dict:
    return json.loads(V08_JSON.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------------------
# Hertz formulas
# --------------------------------------------------------------------------------------


def test_effective_modulus_limits() -> None:
    """Rigid limit, symmetry, equal bodies E/(2(1-nu^2)), stiff body 2 -> rigid."""
    E, nu = 2.0e11, 0.3
    assert H.effective_modulus(E, nu) == pytest.approx(E / (1 - nu**2), rel=1e-15)
    assert H.effective_modulus(E, nu, E, nu) == pytest.approx(E / (2 * (1 - nu**2)), rel=1e-15)
    assert H.effective_modulus(1.0, 0.2, 3.0, 0.4) == pytest.approx(
        H.effective_modulus(3.0, 0.4, 1.0, 0.2), rel=1e-15)
    assert H.effective_modulus(E, nu, 1e9 * E, 0.25) == pytest.approx(E / (1 - nu**2), rel=1e-8)
    assert H.effective_radius(2.0) == 2.0
    assert H.effective_radius(2.0, 2.0) == pytest.approx(1.0, rel=1e-15)


@pytest.mark.parametrize("two_spheres, coef_a, coef_d, coef_p", [
    (False, (1.109, 1e-3), (1.231, 1e-3), (0.388, 1e-3)),   # ball on plane
    (True, (0.881, 1e-3), (1.55, 1e-2), (0.616, 1e-3)),     # two equal balls
])
def test_timoshenko_goodier_coefficients(two_spheres, coef_a, coef_d, coef_p) -> None:
    """Same steel-like material, nu = 0.3: a, approach, p0 coefficients of sec. 140."""
    E, nu, P, R = 2.1e11, 0.3, 1.0e3, 0.01
    Es = H.effective_modulus(E, nu, E, nu)
    Re = H.effective_radius(R, R) if two_spheres else H.effective_radius(R)
    a = H.hertz_contact_radius(P, Re, Es)
    d = H.hertz_delta(P, Re, Es)
    p0 = H.hertz_max_pressure(P, Re, Es)
    assert a / (P * R / E) ** (1 / 3) == pytest.approx(coef_a[0], abs=coef_a[1])
    assert d / (P**2 / (E**2 * R)) ** (1 / 3) == pytest.approx(coef_d[0], abs=coef_d[1])
    assert p0 / (P * E**2 / R**2) ** (1 / 3) == pytest.approx(coef_p[0], abs=coef_p[1])


def test_hertz_scaling_and_round_trip() -> None:
    """a ~ F^(1/3), delta ~ F^(2/3), p0 ~ F^(1/3); F(delta(F)) = F; a(F) = sqrt(R delta)."""
    R, Es, F = 0.5, 1.6e7, 1.0e5
    assert H.hertz_contact_radius(8 * F, R, Es) / H.hertz_contact_radius(F, R, Es) == pytest.approx(2.0, rel=1e-14)
    assert H.hertz_delta(8 * F, R, Es) / H.hertz_delta(F, R, Es) == pytest.approx(4.0, rel=1e-14)
    assert H.hertz_max_pressure(8 * F, R, Es) / H.hertz_max_pressure(F, R, Es) == pytest.approx(2.0, rel=1e-14)
    d = H.hertz_delta(F, R, Es)
    assert H.hertz_force(d, R, Es) == pytest.approx(F, rel=1e-14)
    assert H.hertz_contact_radius(F, R, Es) == pytest.approx(H.hertz_contact_radius_from_delta(d, R), rel=1e-14)
    assert H.hertz_mean_pressure(F, R, Es) == pytest.approx(2 / 3 * H.hertz_max_pressure(F, R, Es), rel=1e-14)
    assert H.hertz_force(-0.1, R, Es) == 0.0


def test_hertz_pressure_integrates_to_force_and_stiffness() -> None:
    """2 pi int p r dr = F (quadrature) and dF/d delta = 2 E* a (central difference)."""
    R, Es, F = 0.5, 1.6e7, 2.0e5
    a = H.hertz_contact_radius(F, R, Es)
    val, _ = quad(lambda r: 2 * math.pi * r * H.hertz_pressure(r, F, R, Es), 0.0, a, epsabs=0, epsrel=1e-12)
    assert val == pytest.approx(F, rel=1e-10)
    d, eps = 0.05, 1e-6
    dF = (H.hertz_force(d + eps, R, Es) - H.hertz_force(d - eps, R, Es)) / (2 * eps)
    assert dF == pytest.approx(H.hertz_contact_stiffness(d, R, Es), rel=1e-8)


def test_hertz_pressure_closes_the_gap_boussinesq() -> None:
    """Defining Hertz condition: u_z(r) from Boussinesq integration of p equals
    delta - r^2/(2R) inside the contact (two bodies: compliance 1/(pi E*))."""
    R, Es, F = 0.5, 1.6e7, 2.0e5
    a = H.hertz_contact_radius(F, R, Es)
    d = H.hertz_delta(F, R, Es)
    nphi = 128
    phis = (np.arange(nphi) + 0.5) * 2 * math.pi / nphi  # periodic integrand: midpoint rule

    def uz(r: float) -> float:
        tot = 0.0
        for ph in phis:
            smax = -r * math.cos(ph) + math.sqrt(a * a - (r * math.sin(ph)) ** 2)
            f = lambda s: H.hertz_pressure(math.hypot(r + s * math.cos(ph), s * math.sin(ph)), F, R, Es)
            tot += quad(f, 0.0, smax, epsabs=0, epsrel=1e-11, limit=200)[0]
        return tot * (2 * math.pi / nphi) / (math.pi * Es)

    for r in (0.0, 0.3 * a, 0.7 * a):
        assert uz(r) == pytest.approx(d - r * r / (2 * R), rel=1e-7)


def test_dundurs_beta_limits(v08) -> None:
    """beta = 0 for identical bodies and for rigid on incompressible; 1/2 for rigid on
    nu = 0 (Johnson's bound); V08 pair is nearly decoupled (|beta| < 2e-3)."""
    assert H.dundurs_beta(2e11, 0.3, 2e11, 0.3) == 0.0
    assert H.dundurs_beta(1.0, 0.5) == 0.0
    assert H.dundurs_beta(1.0, 0.0) == pytest.approx(0.5, rel=1e-15)
    inp = H.V08Inputs.from_json(v08)
    assert abs(H.dundurs_beta(inp.E_block, inp.nu_block, inp.E_indenter, inp.nu_indenter)) < 2e-3
    # rigid indenter replacement with the same E*:
    Eb = H.rigid_equivalent_block_modulus(inp.E_star, inp.nu_block)
    assert H.effective_modulus(Eb, inp.nu_block) == pytest.approx(inp.E_star, rel=1e-14)


# --------------------------------------------------------------------------------------
# Exact spherical gap (Sneddon integrals, computed here numerically)
# --------------------------------------------------------------------------------------


def _sneddon(fp, fpp, a: float, Es: float) -> tuple[float, float, float]:
    """Numerical Sneddon solution for a smooth axisymmetric gap f(r) (f' = fp, f'' = fpp).

    I(t) = int_0^{pi/2} f'(t sin s) ds; g(t) = delta - t I(t), delta = a I(a);
    F = 2 E* int_0^a g dt; p(0) = -(E*/pi) int_0^a g'(t)/t dt with
    g' = -I - t I', I'(t) = int f''(t sin s) sin s ds.
    """
    opts = dict(epsabs=0, epsrel=1e-13, limit=200)
    I = lambda t: quad(lambda s: fp(t * math.sin(s)), 0, math.pi / 2, **opts)[0]
    Ip = lambda t: quad(lambda s: fpp(t * math.sin(s)) * math.sin(s), 0, math.pi / 2, **opts)[0]
    delta = a * I(a)
    F = 2 * Es * quad(lambda t: delta - t * I(t), 0, a, **opts)[0]
    p0 = Es / math.pi * quad(lambda t: (I(t) + t * Ip(t)) / t if t > 0 else 2 * Ip(0.0), 0, a, **opts)[0]
    return delta, F, p0


def test_sneddon_integrals_reproduce_hertz_for_paraboloid() -> None:
    """The numerical machinery itself: f = r^2/(2R) gives a^2/R, 4E*a^3/(3R), 2E*a/(pi R)."""
    R, Es, a = 0.5, 1.6e7, 0.2
    d, F, p0 = _sneddon(lambda r: r / R, lambda r: 1 / R, a, Es)
    assert d == pytest.approx(a * a / R, rel=1e-11)
    assert F == pytest.approx(4 * Es * a**3 / (3 * R), rel=1e-11)
    assert p0 == pytest.approx(2 * Es * a / (math.pi * R), rel=1e-11)


@pytest.mark.parametrize("a_over_R", [0.1, 0.447, 0.8])
def test_exact_sphere_closed_forms_vs_sneddon(a_over_R) -> None:
    """Closed forms (Segedin/Sneddon) == numerical integrals with f = R - sqrt(R^2 - r^2)."""
    R, Es = 0.5, 1.6e7
    a = a_over_R * R
    d, F, p0 = _sneddon(lambda r: r / math.sqrt(R * R - r * r),
                        lambda r: R * R / (R * R - r * r) ** 1.5, a, Es)
    assert H.sphere_exact_delta(a, R) == pytest.approx(d, rel=1e-10)
    assert H.sphere_exact_force(a, R, Es) == pytest.approx(F, rel=1e-9)
    assert H.sphere_exact_p0(a, R, Es) == pytest.approx(p0, rel=1e-10)
    assert H.sphere_exact_contact_radius(d, R) == pytest.approx(a, rel=1e-12)


def test_exact_sphere_tends_to_hertz() -> None:
    """a/R = 1e-3: exact sphere and paraboloid differ by O((a/R)^2) only."""
    R, Es = 0.5, 1.6e7
    a = 1e-3 * R
    d = H.sphere_exact_delta(a, R)
    assert d == pytest.approx(a * a / R, rel=1e-5)
    assert H.sphere_exact_force(a, R, Es) == pytest.approx(H.hertz_force(a * a / R, R, Es), rel=1e-5)
    assert H.sphere_exact_p0(a, R, Es) == pytest.approx(2 * Es * a / (math.pi * R), rel=1e-5)


# --------------------------------------------------------------------------------------
# V08 reference reproduction
# --------------------------------------------------------------------------------------


def test_v08_inputs_consistent(v08) -> None:
    """Inputs: clearance = centre_y - R - top; delta = push - g0; mu parsed = F_T/F_N."""
    inp = H.V08Inputs.from_json(v08)
    assert inp.center[1] - inp.R - inp.block_top_y == pytest.approx(inp.g0, abs=1e-12)
    assert H.V08_PUSH - inp.g0 == pytest.approx(inp.delta, abs=1e-15)
    ref = v08["reference_values"]["Hertz_Coulomb_accepted_at_full_indentation"]
    assert inp.mu == pytest.approx(ref["gross_sliding_force_F_T_lbf"] / ref["normal_force_F_N_lbf"], rel=1e-6)
    # doc sec. 6: delta(t) = max(1.25 t - 0.025, 0) for t <= 0.1, held at 0.1 after
    t = np.array([0.0, 0.02, 0.05, 0.1, 0.4, 1.0])
    want = np.where(t <= 0.1, np.maximum(1.25 * t - 0.025, 0.0), 0.1)
    np.testing.assert_allclose(H.v08_delta(t, inp.g0), want, rtol=0, atol=1e-15)
    ux, uy = H.v08_cut_face_displacement(np.array([0.1, 0.55, 1.0]))
    np.testing.assert_allclose(ux, [0.0, 2.125, 4.25], atol=1e-15)
    np.testing.assert_allclose(uy, [-0.125] * 3, atol=1e-15)


def test_v08_reference_values_reproduced(v08) -> None:
    """a, F_N, F_T, p0 of v08.json from its own inputs. Tolerance 5e-6: the json prints
    6-7 significant digits (a = 0.223607 has rounding up to 2.2e-6)."""
    inp = H.V08Inputs.from_json(v08)
    got = H.v08_reference(inp)
    ref = v08["reference_values"]["Hertz_Coulomb_accepted_at_full_indentation"]
    for key in ("contact_radius_a_in", "normal_force_F_N_lbf", "gross_sliding_force_F_T_lbf",
                "maximum_pressure_p0_psi", "full_indentation_delta_in"):
        assert got[key] == pytest.approx(ref[key], rel=5e-6), key
    # E* against the source document (contact-indenter.md Table 3: 1.615266e+07 psi)
    assert got["reduced_modulus_E_star_psi"] == pytest.approx(1.615266e7, rel=1e-6)


def test_v08_json_E_star_and_block_dimensions(v08) -> None:
    """v08.json stores E* in psi (doc Table 3: 1.615266e7) and block Lx,Ly,Lz = [5, 2, 3] (Table 1)."""
    E_star = H.V08Inputs.from_json(v08).E_star
    assert E_star == pytest.approx(1.615266e7, rel=1e-6)
    for v in (v08["materials"]["reduced_modulus_E_star_psi"],
              v08["reference_values"]["Hertz_Coulomb_accepted_at_full_indentation"]["reduced_modulus_E_star_psi"]):
        assert v == pytest.approx(1.615266e7, rel=1e-6)
    assert v08["geometry"]["block_Lx_Ly_Lz_in"] == [5, 2, 3]


def test_v08_model_form_gap_exact_sphere(v08) -> None:
    """Hertz paraboloid vs exact sphere at V08's delta = 0.1 in (a/R = 0.447): the exact
    solution is recomputed here by the numerical Sneddon integrals (root on delta)."""
    from scipy.optimize import brentq

    inp = H.V08Inputs.from_json(v08)
    R, Es = inp.R, inp.E_star
    fp = lambda r: r / math.sqrt(R * R - r * r)
    fpp = lambda r: R * R / (R * R - r * r) ** 1.5
    a = brentq(lambda x: _sneddon(fp, fpp, x, Es)[0] - inp.delta, 0.1, 0.3, xtol=1e-14)
    _, F, p0 = _sneddon(fp, fpp, a, Es)
    gaps = H.v08_model_form_gaps(inp)
    ref = H.v08_reference(inp)
    assert gaps["exact_a_in"] == pytest.approx(a, rel=1e-9)
    assert gaps["exact_vs_hertz_F_N"] == pytest.approx(F / ref["normal_force_F_N_lbf"] - 1, abs=1e-9)
    assert gaps["exact_vs_hertz_p0"] == pytest.approx(p0 / ref["maximum_pressure_p0_psi"] - 1, abs=1e-9)
    assert gaps["a_over_R"] == pytest.approx(math.sqrt(0.2), rel=1e-12)


def test_v08_history_and_relative_rms(v08) -> None:
    """F_N(t) Hertz, F_T = mu F_N during sliding, doc sec. 7.1 RMS definition."""
    inp = H.V08Inputs.from_json(v08)
    t = np.linspace(0.15, 0.70, 12)
    hist = H.v08_reference_history(t, inp)
    F_ref = v08["reference_values"]["Hertz_Coulomb_accepted_at_full_indentation"]["normal_force_F_N_lbf"]
    np.testing.assert_allclose(hist["F_N"], F_ref, rtol=5e-6)
    np.testing.assert_allclose(hist["F_T"], inp.mu * hist["F_N"], rtol=1e-15)
    assert np.isnan(H.v08_reference_history([0.05], inp)["F_T"]).all()
    assert H.v08_reference_history([0.01], inp)["F_N"][0] == 0.0
    # e = sqrt(sum (q - q*)^2) / sqrt(sum q*^2): hand value for q = [1, 2], q* = [1, 1]
    assert H.relative_rms([1.0, 2.0], [1.0, 1.0]) == pytest.approx(1 / math.sqrt(2), rel=1e-15)


def test_grids_and_dof_counts() -> None:
    """Block grid divides the 5x2x3 block; DOFs vs BoxSpline; nnz vs explicit pattern."""
    import scipy.sparse as sp

    from wht_iga.reference.box import BoxSpline

    assert H.v08_block_grid(0.5) == ((-2.5, -1.0, -1.5), (2.5, 1.0, 1.5), (10, 4, 6))
    assert H.v08_block_grid(0.125, "half")[2] == (40, 16, 12)
    assert H.v08_block_grid(0.25, "quarter")[2] == (16, 8, 6)
    for n in ((10, 4, 6), (3, 2, 2), (5, 3, 3)):
        lo, hi = (0, 0, 0), tuple(float(v) for v in n)
        assert H.bspline_dofs(n) == BoxSpline(lo=lo, hi=hi, n_cells=n, degree=2).n_dofs
        band = [sp.diags([np.ones(v + 2 - abs(k)) for k in range(-2, 3)], list(range(-2, 3))).tocsr()
                for v in n]
        pat = sp.kron(sp.kron(sp.kron(band[2], band[1], format="csr"), band[0], format="csr"),
                      sp.csr_matrix(np.ones((3, 3))), format="csr")
        pat.eliminate_zeros()
        assert H.bspline_stiffness_nnz(n) == pat.nnz
    for h in H.V08_H_LEVELS:
        for sym in ("full", "half", "quarter"):
            lo, hi, n = H.v08_indenter_grid(h, C08, R08, sym=sym)
            assert hi[1] == C08[1]                                   # top plane = cut plane
            np.testing.assert_allclose(np.subtract(hi, lo), h * np.array(n), atol=1e-12)
            assert lo[1] <= C08[1] - R08 - 0.1 * h + 1e-12           # pole clears the bottom
            assert (lo[2] == C08[2]) if sym != "full" else (lo[2] <= C08[2] - R08 - 0.1 * h + 1e-12)
            assert (lo[0] == C08[0]) if sym == "quarter" else (hi[0] >= C08[0] + R08 + 0.1 * h - 1e-12)


# --------------------------------------------------------------------------------------
# Hollow hemisphere soup
# --------------------------------------------------------------------------------------


def _expected_outward(nodes, tri, c, R, Ri, sym):
    """Outward direction from the vertex geometry alone (not from face ids)."""
    p = nodes[tri]
    if np.all(p[:, 1] == c[1]):
        return np.array([0.0, 1.0, 0.0])
    if sym != "full" and np.all(p[:, 2] == c[2]):
        return np.array([0.0, 0.0, -1.0])
    if sym == "quarter" and np.all(p[:, 0] == c[0]):
        return np.array([-1.0, 0.0, 0.0])
    r = np.linalg.norm(p - c, axis=1)
    cen = p.mean(axis=0) - c
    if np.allclose(r, R, rtol=0, atol=1e-14):
        return cen / np.linalg.norm(cen)
    if np.allclose(r, Ri, rtol=0, atol=1e-14):
        return -cen / np.linalg.norm(cen)
    raise AssertionError(f"triangle {tri} lies on no expected surface")


@pytest.mark.parametrize("sym", ["full", "half", "quarter"])
def test_soup_closed_outward_on_surface(sym) -> None:
    """Every directed edge once with its reverse once (watertight, consistent); no
    degenerate triangles; vertices on the exact surfaces; every normal outward;
    closed-surface identity sum(n dA) = 0."""
    c = np.array(C08)
    nodes, tris, fids = H.hollow_hemisphere_soup(R08, RI08, C08, n=12, sym=sym)
    edges = Counter((int(a), int(b)) for t in tris for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])))
    assert all(k == 1 and edges.get((b, a)) == 1 for (a, b), k in edges.items())
    p = nodes[tris]
    nrm = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    area2 = np.linalg.norm(nrm, axis=1)
    assert area2.min() > 1e-6
    r = np.linalg.norm(nodes - c, axis=1)
    assert np.all((r >= RI08 - 1e-14) & (r <= R08 + 1e-14))
    assert np.all(nodes[:, 1] <= c[1])
    for t, nv in zip(tris, nrm):
        assert nv @ _expected_outward(nodes, t, c, R08, RI08, sym) > 0.0
    np.testing.assert_allclose(nrm.sum(axis=0), 0.0, atol=1e-13)
    assert set(np.unique(fids).tolist()) == {"full": {1, 2, 3}, "half": {1, 2, 3, 4},
                                             "quarter": {1, 2, 3, 4, 5}}[sym]


def test_soup_volume_divergence_theorem_converges() -> None:
    """Inscribed facets: volume deficit shrinks as O(n^-2); translation-invariant (closed);
    half/quarter soups are exactly 1/2, 1/4 of the full soup (same facets)."""
    Ve = H.hollow_hemisphere_volume(R08, RI08)
    assert Ve == pytest.approx(2 / 3 * math.pi * (0.125 - 0.015625), rel=1e-15)
    errs = []
    for n in (8, 16, 32):
        nodes, tris, _ = H.hollow_hemisphere_soup(R08, RI08, C08, n=n)
        V = H.soup_signed_volume(nodes, tris)
        V0 = H.soup_signed_volume(*H.hollow_hemisphere_soup(R08, RI08, (0, 0, 0), n=n)[:2])
        assert V == pytest.approx(V0, rel=1e-12)
        errs.append(abs(V - Ve) / Ve)
        for sym, frac in (("half", 0.5), ("quarter", 0.25)):
            Vs = H.soup_signed_volume(*H.hollow_hemisphere_soup(R08, RI08, C08, n=n, sym=sym)[:2])
            assert Vs == pytest.approx(frac * V, rel=1e-12)
    assert errs[0] > errs[1] > errs[2]
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(2)]
    assert all(1.9 < q < 2.1 for q in orders), orders
    assert errs[2] < 5e-3


@pytest.mark.parametrize("n", [12, 24, 48])
def test_soup_contact_region_has_no_slivers(n) -> None:
    """The lowest point (contact pole) is a vertex and outer triangles within 2a of it
    (a = 0.2236 in at full indentation) keep a minimum angle >= 25 deg at every n.
    Reason for 25 deg: a UV sphere with its pole at the contact point would have fan
    apex angles 360/(2n) deg there (7.5 deg at n = 24); this soup's fans sit on the cut
    plane rim instead (measured global minimum 7.4 deg at n = 24, near-contact 30.6 deg)."""
    nodes, tris, fids = H.hollow_hemisphere_soup(R08, RI08, C08, n=n)
    pole = np.array([C08[0], C08[1] - R08, C08[2]])
    assert np.any(np.all(nodes == pole, axis=1))
    p = nodes[tris[fids == H.FACE_OUTER]]
    near = np.linalg.norm(p.mean(axis=1) - pole, axis=1) < 2 * math.sqrt(0.05)
    assert near.sum() >= n
    e = [p[:, (k + 1) % 3] - p[:, k] for k in range(3)]
    cosang = [np.einsum("ij,ij->i", e[k], -e[(k + 2) % 3])
              / (np.linalg.norm(e[k], axis=1) * np.linalg.norm(e[(k + 2) % 3], axis=1)) for k in range(3)]
    ang = np.degrees(np.arccos(np.clip(np.stack(cosang, axis=1), -1.0, 1.0)))
    np.testing.assert_allclose(ang.sum(axis=1), 180.0, atol=1e-9)
    assert ang[near].min() >= 25.0
    assert ang.min() < 360.0 / (2 * n) + 1e-9  # the rim fans exist, away from contact


@pytest.mark.parametrize("sym", ["full", "half", "quarter"])
def test_soup_winding_inside_matches_analytic(sym) -> None:
    """Repo inside test (generalized winding number) vs exact membership on random points
    away from the surfaces by more than the facet sagitta."""
    from wht_iga.geometry import classify

    n = 16
    nodes, tris, _ = H.hollow_hemisphere_soup(R08, RI08, C08, n=n, sym=sym)
    c = np.array(C08)
    rng = np.random.default_rng(8)
    pts = c + rng.uniform(-0.6, 0.6, size=(4000, 3))
    r = np.linalg.norm(pts - c, axis=1)
    inside = (r <= R08) & (r >= RI08) & (pts[:, 1] <= c[1])
    dist = np.minimum.reduce([np.abs(r - R08), np.abs(r - RI08), np.abs(pts[:, 1] - c[1])])
    if sym != "full":
        inside &= pts[:, 2] >= c[2]
        dist = np.minimum(dist, np.abs(pts[:, 2] - c[2]))
    if sym == "quarter":
        inside &= pts[:, 0] >= c[0]
        dist = np.minimum(dist, np.abs(pts[:, 0] - c[0]))
    sag = R08 * (1 - math.cos(math.pi / n))  # > max facet-to-sphere distance
    keep = dist > sag
    got = classify.winding_inside(pts[keep], nodes, tris)
    np.testing.assert_array_equal(got, inside[keep])
    assert keep.sum() > 3000


def test_soup_trimmed_quadrature_volume_full_grid() -> None:
    """Trimmed moment-fit quadrature on the V08 h = 0.25 indenter grid integrates the soup
    volume (moment-fit residual tolerance; measured 2e-11).

    Only the FULL grid is tested here (for historical comparison; ``classify._find_cut_cells_njit``
    has since been fixed to write C order ``(i ny + j) nz + k``, and half-grid trimmed volume
    regression is covered by ``tests/unit/test_geo_classify_order.py``).
    """
    from wht_iga.reference.box import BoxSpline
    from wht_iga.solver.trimmed import build_trimmed_table

    h = 0.25
    nodes, tris, _ = H.hollow_hemisphere_soup(R08, RI08, C08, n=16)
    lo, hi, n = H.v08_indenter_grid(h, C08, R08)
    assert n[0] == n[2] and lo[0] - C08[0] == pytest.approx(lo[2] - C08[2], abs=1e-15)
    box = BoxSpline(lo=lo, hi=hi, n_cells=n, degree=2)
    qp, out_cells, _issues, _rec, diag = build_trimmed_table(nodes, tris, box, q=4, cache_dir="off")
    out = set(out_cells)
    V = 0.0
    for k in range(n[2]):
        for j in range(n[1]):
            for i in range(n[0]):
                if (i, j, k) in out:
                    continue
                V += float(np.asarray(qp[(i, j, k)])[:, 3].sum()) if (i, j, k) in qp else h**3
    assert diag["n_fallback"] == 0
    assert V == pytest.approx(H.soup_signed_volume(nodes, tris), rel=1e-8)
