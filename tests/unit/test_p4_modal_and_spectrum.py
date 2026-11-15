"""P4-B: case files, seven modal rules, four directional rules, random response.

Oracles are independent of the engine code path: closed forms (SDOF white-noise
variance ``S0 / (8 xi w0**3)``), hand-computed scalars, and a pure-Python
double loop over modes for the cross-spectral integral.  What the hand
calculations check is the *implementation of the stated formula*; the
Abaqus-specific conventions (bias formula, TENP/GRP/DSC details) are documented
assumptions in ``response_spectrum.py``, not verified against Abaqus here.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy.integrate import trapezoid

from wht_iga.dynamics.case_file import (
    CaseFileError,
    ModalDamping,
    RandomResponseCase,
    ResponseSpectrumCase,
    Spectrum,
    SpectrumType,
    case_to_dict,
    dump_case_file,
    load_case_file,
    parse_case,
)
from wht_iga.dynamics.consistent_mass import EigenmodeResult, align_degenerate_modes
from wht_iga.dynamics.response_spectrum import (
    VON_MISES_A,
    abaqus_frequency_sweep,
    combine_directions,
    combine_modal_peaks,
    random_response_from_case,
    response_spectrum_from_case,
)

# V09 Table 4 (magnitude m/s^2, frequency Hz), read from the Coreform source page.
V09_TABLE4 = [[20.0, 10.0, 0.02], [200.0, 100.0, 0.02], [200.0, 600.0, 0.02], [60.0, 2000.0, 0.02]]


def _rs_dict(**over) -> dict:
    """A valid one-excitation response-spectrum case dictionary (overridable)."""
    base = {
        "spectra": {"S": {"type": "DISPLACEMENT", "gravity": None, "rows": [[10.0, 1.0, 0.02], [10.0, 1000.0, 0.02]]}},
        "excitations": [{"spectrum": "S", "direction": [1, 0, 0], "scale": 1.0, "duration": None}],
        "modal_summation": "SRSS",
        "directional_summation": "SRSS",
        "damping": 0.02,
    }
    base.update(over)
    return base


def _random_dict(**over) -> dict:
    """A valid single-base-motion random case dictionary (overridable)."""
    base = {
        "psds": {"P": {"type": "BASE", "gravity": 9.81, "rows": [[0.01, 0.0, 10.0], [0.01, 0.0, 2000.0]]}},
        "base_motions": [{"dof": 2, "load_case": 1, "kind": "ACCELERATION"}],
        "correlations": [{"psd": "P", "type": "UNCORRELATED", "weights": {"1": [1.0, 0.0]}}],
        "lower_frequency_hz": 10.0, "upper_frequency_hz": 2000.0, "n_points": 20, "bias": 3.0,
        "frequency_scale": "LOGARITHMIC", "damping": 0.02,
    }
    base.update(over)
    return base


# --------------------------------------------------------------------------- spectrum
def test_spectrum_log_log_interpolation_and_clamp() -> None:
    """Log-log interpolation: S = 2 f between (10 Hz, 20) and (100 Hz, 200); constant beyond the ends."""
    sp = Spectrum.from_rows("T", V09_TABLE4)
    f = np.array([10.0, 31.622776601683793, 100.0, 245.0, 1095.4451150103323, 2000.0, 5.0, 9999.0])
    got = sp.tabulated(f, 0.02)
    want = [20.0, 63.245553203367585, 200.0, 200.0, 200.0 * math.sqrt(60.0 / 200.0), 60.0, 20.0, 60.0]
    np.testing.assert_allclose(got, want, rtol=1e-12)


def test_spectrum_damping_interpolation_is_logarithmic_and_clamped() -> None:
    """Between damping curves the weight is linear in xi and the magnitude interpolates logarithmically."""
    rows = [[10.0, 1.0, 0.02], [10.0, 100.0, 0.02], [2.5, 1.0, 0.06], [2.5, 100.0, 0.06]]
    sp = Spectrum.from_rows("T", rows)
    got = sp.tabulated(np.array([5.0, 5.0, 5.0, 5.0]), np.array([0.02, 0.04, 0.06, 0.5]))
    np.testing.assert_allclose(got, [10.0, math.sqrt(10.0 * 2.5), 2.5, 2.5], rtol=1e-12)
    assert sp.tabulated(np.array([5.0]), 0.0)[0] == pytest.approx(10.0)  # clamped below


def test_spectrum_type_conversions() -> None:
    """S_D = S_V / w = S_A / w**2 whichever type is supplied; G multiplies by gravity."""
    f, w = 100.0, 2.0 * math.pi * 100.0
    acc = Spectrum.from_rows("A", [[200.0, 1.0, 0.02], [200.0, 1e4, 0.02]], SpectrumType.ACCELERATION)
    assert acc.at("D", np.array([f]), 0.02)[0] == pytest.approx(200.0 / w**2, rel=1e-13)
    assert acc.at("V", np.array([f]), 0.02)[0] == pytest.approx(200.0 / w, rel=1e-13)
    vel = Spectrum.from_rows("V", [[0.3, 1.0, 0.02], [0.3, 1e4, 0.02]], SpectrumType.VELOCITY)
    assert vel.at("A", np.array([f]), 0.02)[0] == pytest.approx(0.3 * w, rel=1e-13)
    assert vel.at("D", np.array([f]), 0.02)[0] == pytest.approx(0.3 / w, rel=1e-13)
    dsp = Spectrum.from_rows("D", [[1e-3, 1.0, 0.02], [1e-3, 1e4, 0.02]], SpectrumType.DISPLACEMENT)
    assert dsp.at("A", np.array([f]), 0.02)[0] == pytest.approx(1e-3 * w * w, rel=1e-13)
    g = Spectrum.from_rows("G", [[2.0, 1.0, 0.02], [2.0, 1e4, 0.02]], SpectrumType.G, gravity=9.81)
    assert g.at("A", np.array([f]), 0.02)[0] == pytest.approx(19.62, rel=1e-13)


# --------------------------------------------------------------------------- modal rules
W3 = np.array([100.0, 108.0, 116.0])  # 8 % and 7.4 % neighbours, 16 % first-to-third
Z3 = np.full(3, 0.02)
R3 = np.array([[1.0, 2.0, 3.0]])


def test_modal_rules_abs_srss_nrl_hand_values() -> None:
    """ABS = 6, SRSS = sqrt(14), NRL = 3 + sqrt(1 + 4)."""
    assert combine_modal_peaks(R3, W3, Z3, "ABS")[0] == pytest.approx(6.0, abs=1e-14)
    assert combine_modal_peaks(R3, W3, Z3, "SRSS")[0] == pytest.approx(math.sqrt(14.0), abs=1e-14)
    assert combine_modal_peaks(R3, W3, Z3, "NRL")[0] == pytest.approx(3.0 + math.sqrt(5.0), abs=1e-14)


def test_tenp_and_grp_differ_on_chained_modes() -> None:
    """Pairs (1,2) and (2,3) are within 10 %, (1,3) is not: TENP counts both pairs, GRP groups {1,2},{3}."""
    assert combine_modal_peaks(R3, W3, Z3, "TENP")[0] == pytest.approx(math.sqrt(14.0 + 2.0 * (2.0 + 6.0)), abs=1e-13)
    assert combine_modal_peaks(R3, W3, Z3, "GRP")[0] == pytest.approx(math.sqrt((1.0 + 2.0) ** 2 + 3.0**2), abs=1e-13)


def test_tenp_grp_ignore_sign_but_cqc_does_not() -> None:
    """Two coincident modes with opposite signs: TENP/GRP add magnitudes, CQC (rho = 1) cancels."""
    w = np.array([100.0, 100.0])
    z = np.array([0.05, 0.05])
    r = np.array([[1.0, -1.0]])
    assert combine_modal_peaks(r, w, z, "TENP")[0] == pytest.approx(2.0, abs=1e-14)
    assert combine_modal_peaks(r, w, z, "GRP")[0] == pytest.approx(2.0, abs=1e-14)
    assert combine_modal_peaks(r, w, z, "CQC")[0] == pytest.approx(0.0, abs=1e-7)


def test_cqc_matches_scalar_der_kiureghian_formula() -> None:
    """Two modes: R**2 = r1**2 + r2**2 + 2 rho r1 r2 with rho from the explicit scalar formula."""
    w1, w2, z1, z2 = 100.0, 112.0, 0.03, 0.05
    r = w2 / w1
    rho = 8 * math.sqrt(z1 * z2) * (z1 + r * z2) * r**1.5 / (
        (1 - r * r) ** 2 + 4 * z1 * z2 * r * (1 + r * r) + 4 * (z1 * z1 + z2 * z2) * r * r
    )
    got = combine_modal_peaks(np.array([[2.0, -1.5]]), np.array([w1, w2]), np.array([z1, z2]), "CQC")[0]
    assert got == pytest.approx(math.sqrt(4.0 + 2.25 + 2.0 * rho * 2.0 * (-1.5)), abs=1e-13)


def test_dsc_matches_scalar_double_sum_formula() -> None:
    """DSC: R**2 = sum_ij |Ri||Rj| eps_ij with eps from the explicit scalar RG 1.92 formula."""
    w = [100.0, 104.0]
    z = [0.02, 0.04]
    td = 5.0
    wd = [w[k] * math.sqrt(1 - z[k] ** 2) for k in range(2)]
    zp = [z[k] + 2.0 / (w[k] * td) for k in range(2)]
    eps = 1.0 / (1.0 + ((wd[0] - wd[1]) / (zp[0] * w[0] + zp[1] * w[1])) ** 2)
    got = combine_modal_peaks(np.array([[3.0, -4.0]]), np.array(w), np.array(z), "DSC", duration=td)[0]
    assert got == pytest.approx(math.sqrt(9.0 + 16.0 + 2.0 * 3.0 * 4.0 * eps), abs=1e-13)
    with pytest.raises(ValueError, match="duration"):
        combine_modal_peaks(R3, W3, Z3, "DSC")


def test_modal_rules_are_order_independent_and_validate() -> None:
    """Unsorted modes give the same answer; unknown rules and non-positive frequencies fail."""
    perm = [2, 0, 1]
    for rule in ("TENP", "GRP", "NRL", "CQC"):
        a = combine_modal_peaks(R3, W3, Z3, rule)
        b = combine_modal_peaks(R3[:, perm], W3[perm], Z3[perm], rule)
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-14)
    with pytest.raises(ValueError, match="unknown"):
        combine_modal_peaks(R3, W3, Z3, "BOGUS")
    with pytest.raises(ValueError, match="positive"):
        combine_modal_peaks(R3, np.array([0.0, 1.0, 2.0]), Z3, "SRSS")


def test_directional_rules_hand_values() -> None:
    """[3,4,5]: SRSS sqrt(50); R40 = 5 + 0.4 (3 + 4); R30 = 5 + 0.3 (3 + 4)."""
    p = np.array([3.0, 4.0, 5.0])
    assert combine_directions(p, "SRSS") == pytest.approx(math.sqrt(50.0), abs=1e-14)
    assert combine_directions(p, "R40") == pytest.approx(7.8, abs=1e-14)
    assert combine_directions(p, "R30") == pytest.approx(7.1, abs=1e-14)
    # brute-force the "largest over the dominant direction" definition
    best = max(abs(p[a]) + 0.4 * sum(abs(p[b]) for b in range(3) if b != a) for a in range(3))
    assert best == pytest.approx(7.8, abs=1e-14)


# --------------------------------------------------------------------------- spectrum case
def test_response_spectrum_case_directional_orderings() -> None:
    """One mode, directions x and y: ALGEBRAIC merges first (|10 - 15| = 5), the rest combine afterwards."""
    d = _rs_dict(
        excitations=[
            {"spectrum": "S", "direction": [1, 0, 0], "scale": 1.0, "duration": None},
            {"spectrum": "S", "direction": [0, 1, 0], "scale": 1.0, "duration": None},
        ]
    )
    w = np.array([2.0 * math.pi * 50.0])
    shapes = np.array([[0.5]])
    gam = np.array([[2.0, -3.0, 0.0]])  # R_x = 2*10*0.5 = 10, R_y = -3*10*0.5 = -15
    want = {"ALGEBRAIC": 5.0, "SRSS": math.sqrt(325.0), "R40": 15.0 + 0.4 * 10.0, "R30": 15.0 + 0.3 * 10.0}
    for rule, value in want.items():
        d["directional_summation"] = rule
        res = response_spectrum_from_case(w, shapes, gam, parse_case(d))
        assert res.peaks[0] == pytest.approx(value, abs=1e-12), rule


def test_response_spectrum_case_quantities_and_scale() -> None:
    """Acceleration spectrum 200 m/s^2 at 100 Hz: U uses S/w^2, V uses S/w, A uses S; scale multiplies."""
    d = _rs_dict(
        spectra={"S": {"type": "ACCELERATION", "gravity": None, "rows": [[200.0, 1.0, 0.02], [200.0, 1e4, 0.02]]}},
        excitations=[{"spectrum": "S", "direction": [0, 1, 0], "scale": 0.5, "duration": None}],
    )
    w0 = 2.0 * math.pi * 100.0
    gam = np.array([[0.0, 4.0, 0.0]])
    shapes = np.array([[0.25]])
    case = parse_case(d)
    for q, s in (("U", 200.0 / w0**2), ("V", 200.0 / w0), ("A", 200.0), ("S", 200.0 / w0**2)):
        got = response_spectrum_from_case(np.array([w0]), shapes, gam, case, q).peaks[0]
        assert got == pytest.approx(4.0 * 0.5 * s * 0.25, rel=1e-13), q


def test_response_spectrum_stress_shape_components() -> None:
    """3D shapes (n_pts, n_modes, 6) give per-component peaks of shape (n_pts, 6)."""
    case = parse_case(_rs_dict())
    w = np.array([2.0 * math.pi * 30.0, 2.0 * math.pi * 90.0])
    psi = np.zeros((2, 2, 6))
    psi[0, 0, 0], psi[0, 1, 0] = 3.0, 4.0
    gam = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    res = response_spectrum_from_case(w, psi, gam, case, "S")
    assert res.peaks.shape == (2, 6)
    assert res.peaks[0, 0] == pytest.approx(10.0 * math.sqrt(3.0**2 + 4.0**2), rel=1e-13)  # S_D = 10 const
    assert np.all(res.peaks[1] == 0.0)


# --------------------------------------------------------------------------- case file
def test_case_file_roundtrip_and_disk(tmp_path) -> None:
    """parse -> dict -> parse is lossless; dump/load use UTF-8 JSON without BOM."""
    for raw in (_rs_dict(damping={"mode_ranges": [[1, 10, 0.05], [11, 20, 0.08]]}), _random_dict()):
        case = parse_case(raw)
        again = parse_case(json.loads(json.dumps(case_to_dict(case))))
        assert case_to_dict(again) == case_to_dict(case)
        path = tmp_path / "case.json"
        dump_case_file(case, path)
        assert not path.read_bytes().startswith(b"\xef\xbb\xbf")
        assert case_to_dict(load_case_file(path)) == case_to_dict(case)
    assert isinstance(parse_case(_rs_dict()), ResponseSpectrumCase)
    assert isinstance(parse_case(_random_dict()), RandomResponseCase)


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda d: d.pop("damping"), "missing"),
        (lambda d: d.update(extra=1), "unknown"),
        (lambda d: d["excitations"][0].update(direction=[2, 0, 0]), "unit vector"),
        (lambda d: d.update(excitations=d["excitations"] * 2), "orthogonal"),
        (lambda d: d.update(modal_summation="DSC"), "duration"),
        (lambda d: d.update(modal_summation="BOGUS"), "modal_summation"),
        (lambda d: d["spectra"]["S"]["rows"].append([5.0, 0.5, 0.02]), "strictly increasing"),
        (lambda d: d["spectra"]["S"].update(rows=[[10.0, 1.0, 0.02], [5.0, 2.0, 0.05], [5.0, 3.0, 0.02]]), "repeated"),
        (lambda d: d["spectra"]["S"].update(type="G"), "gravity"),
        (lambda d: d["excitations"][0].update(spectrum="NOPE"), "not defined"),
    ],
)
def test_response_spectrum_case_validation(mutate, message) -> None:
    """Every documented rule is enforced with a message naming the problem."""
    d = _rs_dict()
    mutate(d)
    with pytest.raises(CaseFileError, match=message):
        parse_case(d)


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda d: d.update(bias=0.5), "bias"),
        (lambda d: d.update(frequency_scale="SQRT"), "frequency_scale"),
        (lambda d: d.update(lower_frequency_hz=3000.0), "lower"),
        (lambda d: d["base_motions"][0].update(dof=5), "rotational"),
        (lambda d: d["psds"]["P"].pop("gravity"), "gravity"),
        (lambda d: d["correlations"][0].update(type="MAYBE"), "CORRELATED"),
        (lambda d: d["correlations"][0].update(weights={"7": [1.0, 0.0]}), "no base motion"),
        (lambda d: d["base_motions"].append({"dof": 1, "load_case": 1, "kind": "VELOCITY"}), "share a kind"),
    ],
)
def test_random_case_validation(mutate, message) -> None:
    """Random-case rules from the user guide are enforced."""
    d = _random_dict()
    mutate(d)
    with pytest.raises(CaseFileError, match=message):
        parse_case(d)


def test_unsupported_psd_types_fail_loudly() -> None:
    """DB tables (band-number map not published) and FORCE PSDs raise instead of guessing."""
    db = _random_dict()
    db["psds"]["P"] = {"type": "DB", "db_reference": 1e-6, "rows": [[10.0, 0.0, 1.0], [10.0, 0.0, 2.0]]}
    case = parse_case(db)
    w = np.array([2 * math.pi * 100.0])
    with pytest.raises(CaseFileError, match="band"):
        random_response_from_case(w, np.ones((1, 1)), np.ones((1, 3)), case)
    frc = _random_dict()
    frc["psds"]["P"]["type"] = "FORCE"
    with pytest.raises(CaseFileError, match="FORCE"):
        random_response_from_case(w, np.ones((1, 1)), np.ones((1, 3)), parse_case(frc))


def test_modal_damping_three_spellings() -> None:
    """Uniform, mode ranges (uncovered modes get 0) and a frequency table (linear, clamped)."""
    f = np.array([10.0, 100.0, 1000.0, 5000.0])
    np.testing.assert_allclose(ModalDamping.parse(0.05).ratios(f), 0.05)
    np.testing.assert_allclose(ModalDamping.parse({"mode_ranges": [[1, 2, 0.05], [3, 3, 0.08]]}).ratios(f), [0.05, 0.05, 0.08, 0.0])
    tab = ModalDamping.parse({"frequency_table": [[100.0, 0.02], [1000.0, 0.08]]})
    np.testing.assert_allclose(tab.ratios(f), [0.02, 0.02, 0.08, 0.08], atol=1e-14)
    np.testing.assert_allclose(tab.ratios(np.array([550.0])), [0.05], atol=1e-14)


# --------------------------------------------------------------------------- sweep
def test_sweep_breaks_at_eigenfrequencies_and_counts_points() -> None:
    """Interval ends are the in-range eigenfrequencies; a degenerate pair adds one break, not two."""
    fm = np.array([50.0, 81.6, 81.6, 492.0, 3000.0])
    f = abaqus_frequency_sweep(fm, 10.0, 2000.0, n_points=20, bias=3.0, scale="LOGARITHMIC")
    for edge in (10.0, 50.0, 81.6, 492.0, 2000.0):
        assert np.any(np.isclose(f, edge, rtol=1e-12)), edge
    assert not np.any(np.isclose(f, 3000.0))
    assert f.size == 4 * 20 - 3  # four intervals of 20 points sharing three interior ends
    assert np.all(np.diff(f) > 0.0)


def test_sweep_bias_one_is_uniform_and_bias_clusters_toward_ends() -> None:
    """bias = 1 is uniform; Abaqus bias uses |y|^(1/p) (endpoint slope 1/p), sinh uses mid/end = bias."""
    lin = abaqus_frequency_sweep(np.array([]), 10.0, 20.0, 11, 1.0, "LINEAR")
    np.testing.assert_allclose(lin, np.linspace(10.0, 20.0, 11), rtol=0, atol=1e-12)

    # Legacy sensitivity model: exact mid/end ratio = bias
    f_sinh = abaqus_frequency_sweep(np.array([]), 100.0, 200.0, 4001, 3.0, "LINEAR", model="sinh_ratio")
    d_sinh = np.diff(f_sinh)
    assert d_sinh[len(d_sinh) // 2] / d_sinh[0] == pytest.approx(3.0, rel=2e-3)
    assert np.all(d_sinh > 0.0)

    # Abaqus default model (and "power" alias): f_hat = mid + half * sign(y) * |y|^(1/p)
    f_abq = abaqus_frequency_sweep(np.array([]), 100.0, 200.0, 4001, 3.0, "LINEAR")
    f_pow = abaqus_frequency_sweep(np.array([]), 100.0, 200.0, 4001, 3.0, "LINEAR", model="power")
    np.testing.assert_allclose(f_abq, f_pow, rtol=0, atol=1e-13)
    d_abq = np.diff(f_abq)
    assert np.all(d_abq > 0.0)
    assert d_abq[0] == pytest.approx(d_abq[-1], rel=1e-12)
    # At |y|=1, d(|y|^(1/p))/dy = 1/p -> endpoint step is (1/p) of uniform spacing 100/4000
    assert d_abq[0] / (100.0 / 4000.0) == pytest.approx(1.0 / 3.0, rel=2e-3)
    assert d_abq[len(d_abq) // 2] > d_abq[0]



# --------------------------------------------------------------------------- degenerate modes
def _pair_eig(theta: float) -> EigenmodeResult:
    """Two degenerate bending modes rotated by ``theta`` inside their eigenspace, plus one distinct mode.

    Six DOFs (two nodes x 3), M = I.  Aligned modes: ``a`` bends in y, ``b`` bends in x.
    """
    c = math.sqrt(1.25)
    a = np.array([0, 1, 0, 0, 0.5, 0]) / c
    b = np.array([1, 0, 0, 0.5, 0, 0]) / c
    third = np.array([0, 0.5, 0, 0, -1, 0]) / c
    phi = np.column_stack([math.cos(theta) * a + math.sin(theta) * b, -math.sin(theta) * a + math.cos(theta) * b, third])
    r = np.zeros((6, 3))
    r[0::3, 0] = 1.0
    r[1::3, 1] = 1.0
    r[2::3, 2] = 1.0
    gam = phi.T @ r
    f = np.array([50.0, 50.0, 120.0])
    return EigenmodeResult(f, 2 * np.pi * f, (2 * np.pi * f) ** 2, phi, gam, gam**2, np.full(3, 2.0), np.arange(6))


def _y_spectrum_case() -> ResponseSpectrumCase:
    return parse_case(_rs_dict(
        spectra={"S": {"type": "DISPLACEMENT", "gravity": None, "rows": [[2.0, 1.0, 0.02], [2.0, 1000.0, 0.02]]}},
        excitations=[{"spectrum": "S", "direction": [0, 1, 0], "scale": 1.0, "duration": None}],
    ))


def test_align_degenerate_modes_makes_srss_orientation_independent() -> None:
    """SRSS on a degenerate pair depends on the eigenvector rotation; after alignment it equals the invariant sum."""
    case = _y_spectrum_case()
    results = {}
    for theta in (0.0, 0.4, 0.9, 1.3):
        eig = _pair_eig(theta)
        raw = response_spectrum_from_case(eig.omega, eig.modes[[1]], eig.participation_factors, case).peaks[0]
        al = align_degenerate_modes(eig, [0, 1, 0])
        aligned = response_spectrum_from_case(al.omega, al.modes[[1]], al.participation_factors, case).peaks[0]
        results[theta] = (raw, aligned)
        # participation: first mode of the cluster carries everything, the partner is zero
        assert abs(al.participation_factors[1, 1]) < 1e-14
        assert al.participation_factors[0, 1] == pytest.approx(1.5 / math.sqrt(1.25), rel=1e-14)  # positive by convention
    raws = [v[0] for v in results.values()]
    aligned = [v[1] for v in results.values()]
    assert max(raws) / min(raws) > 1.2  # unaligned result swings with the arbitrary orientation
    assert max(aligned) == pytest.approx(min(aligned), rel=1e-13)
    # analytic: the aligned pair contributes |Gamma_a phi_a[1]| S_D exactly (here S_D = 2) plus the third mode
    g_a, phi_a_y = 1.5 / math.sqrt(1.25), 1.0 / math.sqrt(1.25)
    third = abs(_pair_eig(0.0).participation_factors[2, 1] * _pair_eig(0.0).modes[1, 2]) * 2.0
    assert aligned[0] == pytest.approx(math.hypot(g_a * phi_a_y * 2.0, third), rel=1e-13)


def test_align_degenerate_modes_preserves_basis_and_invariant_results() -> None:
    """Rotation keeps mass-orthonormality, participation = Phi^T M r, untouched modes, and random RMS."""
    eig = _pair_eig(0.7)
    al = align_degenerate_modes(eig, [0, 1, 0])
    assert np.abs(al.modes.T @ al.modes - np.eye(3)).max() < 1e-14
    r = np.zeros((6, 3))
    r[0::3, 0] = r[1::3, 1] = r[2::3, 2] = 1.0
    np.testing.assert_allclose(al.participation_factors, al.modes.T @ r, atol=1e-14)
    np.testing.assert_array_equal(al.modes[:, 2], eig.modes[:, 2])  # distinct mode untouched
    np.testing.assert_allclose(al.effective_modal_masses, al.participation_factors**2)
    # cross terms between the rotated modes are kept by the random-response engine: invariant
    case = parse_case(_random_dict(damping=0.02))
    grid = np.geomspace(10.0, 400.0, 4001)
    a = random_response_from_case(eig.omega, eig.modes, eig.participation_factors, case, frequencies_hz=grid)
    b = random_response_from_case(al.omega, al.modes, al.participation_factors, case, frequencies_hz=grid)
    np.testing.assert_allclose(a.ru, b.ru, rtol=1e-12, atol=1e-12 * a.ru.max())  # x DOFs are exactly 0 up to 1e-20 noise
    # a cluster with zero participation in the excitation direction is left alone
    same = align_degenerate_modes(eig, [0, 0, 1])
    np.testing.assert_array_equal(same.modes[:, :2], eig.modes[:, :2])
    with pytest.raises(ValueError, match="non-zero"):
        align_degenerate_modes(eig, [0, 0, 0])


# --------------------------------------------------------------------------- random response
def test_random_sdof_white_noise_closed_form() -> None:
    """Single mode, flat S0 per Hz: sigma**2 = S0 / (8 xi w0**3) (Coreform guide check)."""
    f0, xi, s0 = 100.0, 0.02, 0.0123
    w0 = 2.0 * math.pi * f0
    d = _random_dict(
        psds={"P": {"type": "FORCE", "rows": [[s0, 0.0, 1e-3], [s0, 0.0, 1e6]]}},
        damping=xi,
    )
    d["psds"]["P"]["type"] = "BASE"
    d["psds"]["P"]["gravity"] = 1.0
    case = parse_case(d)
    grid = np.concatenate([np.geomspace(0.1, 400.0, 6000), np.linspace(f0 * (1 - 6 * xi), f0 * (1 + 6 * xi), 6001)])
    grid = np.unique(grid)
    res = random_response_from_case(
        np.array([w0]), np.array([[1.0]]), np.array([[0.0, 1.0, 0.0]]), case, frequencies_hz=grid
    )
    # tail beyond 400 Hz: int |H|^2 s0 df ~ s0 / (3 (2 pi)^4 f^3) -> relative 1e-4 here
    assert res.ru[0] ** 2 == pytest.approx(s0 / (8.0 * xi * w0**3), rel=2e-3)


def _scalar_reference_rms(w, phi, gam, zeta, power_of_f, freqs, lam):
    """Pure-Python double loop: sqrt(trapz(sum_ab phi_a phi_b G_a G_b H_a conj(H_b) P(f) (2 pi f)^(2 lam)))."""
    vals = []
    for f in freqs:
        om = 2.0 * math.pi * f
        s = 0.0
        for a in range(len(w)):
            for b in range(len(w)):
                ha = 1.0 / (w[a] ** 2 - om**2 + 2j * zeta[a] * w[a] * om)
                hb = 1.0 / (w[b] ** 2 - om**2 + 2j * zeta[b] * w[b] * om)
                s += phi[a] * phi[b] * gam[a] * gam[b] * (ha * hb.conjugate()).real
        vals.append(s * power_of_f(f) * om ** (2 * lam))
    return math.sqrt(float(trapezoid(vals, freqs)))


@pytest.mark.parametrize("kind, lam", [("ACCELERATION", 0), ("VELOCITY", 1), ("DISPLACEMENT", 2)])
def test_random_two_modes_matches_scalar_double_loop(kind, lam) -> None:
    """Cross terms, base-motion kind exponent and the g**2 factor match an independent scalar loop to 1e-12."""
    w = np.array([2 * math.pi * 80.0, 2 * math.pi * 95.0])  # close modes: cross term matters
    zeta = np.array([0.03, 0.05])
    phi = np.array([0.7, -0.4])
    gam = np.array([1.3, 0.9])
    g = 9.81
    case = parse_case(
        _random_dict(
            psds={"P": {"type": "BASE", "gravity": g, "rows": [[0.01, 0.0, 10.0], [0.04, 0.0, 200.0]]}},
            base_motions=[{"dof": 2, "load_case": 1, "kind": kind}],
            damping={"mode_ranges": [[1, 1, 0.03], [2, 2, 0.05]]},
        )
    )
    grid = np.geomspace(20.0, 400.0, 3001)
    res = random_response_from_case(
        w, phi[None, :], np.column_stack([np.zeros(2), gam, np.zeros(2)]), case, frequencies_hz=grid
    )

    def power(f: float) -> float:  # log-log 0.01 @10 Hz -> 0.04 @200 Hz, times g^2
        return g * g * math.exp(np.interp(math.log(f), [math.log(10.0), math.log(200.0)], [math.log(0.01), math.log(0.04)]))

    want = _scalar_reference_rms(w, phi, gam, zeta, power, grid, lam)
    assert res.ru[0] == pytest.approx(want, rel=1e-12)
    # velocity / acceleration carry (2 pi f)^2 and (2 pi f)^4 on top of the displacement integrand
    want_v = _scalar_reference_rms(w, phi, gam, zeta, power, grid, lam + 1)
    want_a = _scalar_reference_rms(w, phi, gam, zeta, power, grid, lam + 2)
    assert res.rv[0] == pytest.approx(want_v, rel=1e-12)
    assert res.ra[0] == pytest.approx(want_a, rel=1e-12)


def test_random_correlated_vs_uncorrelated_two_dofs() -> None:
    """CORRELATED adds participation factors (Gx + Gy)^2; UNCORRELATED adds squares Gx^2 + Gy^2."""
    w = np.array([2 * math.pi * 100.0])
    gam = np.array([[1.5, 2.0, 0.0]])
    phi = np.array([[1.0]])
    grid = np.geomspace(50.0, 200.0, 801)
    moves = [{"dof": 1, "load_case": 1, "kind": "ACCELERATION"}, {"dof": 2, "load_case": 1, "kind": "ACCELERATION"}]
    out = {}
    for ctype in ("CORRELATED", "UNCORRELATED"):
        case = parse_case(
            _random_dict(base_motions=moves, correlations=[{"psd": "P", "type": ctype, "weights": {"1": [1.0, 0.0]}}])
        )
        out[ctype] = random_response_from_case(w, phi, gam, case, frequencies_hz=grid).ru[0] ** 2
    unit = random_response_from_case(
        w, phi, np.array([[1.0, 0.0, 0.0]]), parse_case(_random_dict(base_motions=moves[:1])), frequencies_hz=grid
    ).ru[0] ** 2
    assert out["CORRELATED"] == pytest.approx(unit * (1.5 + 2.0) ** 2, rel=1e-12)
    assert out["UNCORRELATED"] == pytest.approx(unit * (1.5**2 + 2.0**2), rel=1e-12)
    weighted = parse_case(
        _random_dict(base_motions=moves[:1], correlations=[{"psd": "P", "type": "UNCORRELATED", "weights": {"1": [0.0, 3.0]}}])
    )
    got = random_response_from_case(w, phi, np.array([[1.0, 0.0, 0.0]]), weighted, frequencies_hz=grid).ru[0] ** 2
    assert got == pytest.approx(9.0 * unit, rel=1e-12)  # |w|^2 = 0^2 + 3^2


def test_segalman_rmises_hand_values() -> None:
    """Uniaxial s: RMISES = |s| sigma_q; pure shear t: sqrt(3)|t| sigma_q; two modes: V11 + V22 - V12."""
    w = np.array([2 * math.pi * 60.0, 2 * math.pi * 70.0])
    case = parse_case(_random_dict(damping=0.05, base_motions=[{"dof": 1, "load_case": 1, "kind": "ACCELERATION"}]))
    gam = np.array([[1.0, 0.0, 0.0], [0.6, 0.0, 0.0]])
    grid = np.geomspace(20.0, 200.0, 2001)
    phi = np.ones((1, 2))
    psi = np.zeros((3, 2, 6))
    psi[0, 0, 0] = 5.0  # point 0: uniaxial s = 5 in mode 1 only
    psi[1, 0, 3] = 2.0  # point 1: pure shear t = 2 in mode 1 only
    psi[2, 0, 0], psi[2, 1, 1] = 1.0, 1.0  # point 2: sxx from mode 1, syy from mode 2
    res = random_response_from_case(w, phi, gam, case, stress_modes=psi, frequencies_hz=grid)
    v = res.modal_variance
    assert res.rmises[0] == pytest.approx(5.0 * math.sqrt(v[0, 0]), rel=1e-13)
    assert res.rmises[1] == pytest.approx(math.sqrt(3.0) * 2.0 * math.sqrt(v[0, 0]), rel=1e-13)
    assert res.rmises[2] == pytest.approx(math.sqrt(v[0, 0] + v[1, 1] - v[0, 1]), rel=1e-13)  # T12 = -1/2
    assert VON_MISES_A[0, 1] == -0.5 and VON_MISES_A[3, 3] == 3.0
    with pytest.raises(ValueError, match="stress_modes"):
        random_response_from_case(w, phi, gam, case, stress_modes=np.zeros((1, 2, 5)), frequencies_hz=grid)


def test_random_sweep_converges_with_points_for_sdof() -> None:
    """SDOF sweep converges to closed form as n_points grows (band truncation [10, 2000] Hz ~ -0.26%)."""
    f0, xi, s0 = 100.0, 0.02, 0.01
    w0 = 2.0 * math.pi * f0
    closed = s0 / (8.0 * xi * w0**3)
    errs_sinh = []
    errs_abq = []
    for n in (20, 80, 320, 1280):
        case = parse_case(
            _random_dict(
                n_points=n, damping=xi,
                psds={"P": {"type": "BASE", "gravity": 1.0, "rows": [[s0, 0.0, 10.0], [s0, 0.0, 2000.0]]}},
            )
        )
        res_s = random_response_from_case(np.array([w0]), np.array([[1.0]]), np.array([[0.0, 1.0, 0.0]]),
                                          case, bias_model="sinh_ratio")
        res_a = random_response_from_case(np.array([w0]), np.array([[1.0]]), np.array([[0.0, 1.0, 0.0]]),
                                          case, bias_model="abaqus")
        errs_sinh.append(abs(res_s.ru[0] ** 2 / closed - 1.0))
        errs_abq.append(abs(res_a.ru[0] ** 2 / closed - 1.0))
    assert errs_sinh[0] > errs_sinh[1] > errs_sinh[2] > 0.0  # sinh coarse sweeps overshoot monotonically
    assert errs_sinh[-1] < 5e-3
    assert max(errs_abq[:2]) > errs_abq[-1] or errs_abq[-1] < 5e-3
    assert errs_abq[-1] < 5e-3


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
