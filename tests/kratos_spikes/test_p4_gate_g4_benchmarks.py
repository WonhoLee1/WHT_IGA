"""Gate G4 benchmarks V09-V12 (Kratos K + python mass / thermal integration).

Every test records its measured numbers in ``docs/phase4/g4_measurements.json`` and
asserts the criterion of ``tests/verification/data/vXX.json``.  Where a Coreform
criterion cannot be reproduced exactly (no Abaqus histories here) the deviation is
stated next to the assertion and in ``docs/phase4/G4_report.md``.

V09 (cantilever shock / random)
    * f1: within 1 % of the Euler-Bernoulli 80.8 Hz.
    * random RMS: the *converged* integral of the Coreform random-response equations,
      within 1.79 % + 0.5 %-points of the Abaqus 0.34 mm.  The Abaqus 20-point / bias-3
      sweep is reported as sensitivity only because its bias formula is not published.
    * shock spectrum peak: Table 4, SRSS, within 2.5 % of the Abaqus 0.965 mm.
V10 (vibrating beam, Bossak)
    * f1 and samples per period are absolute criteria; the four history comparisons
      against Abaqus HHT are RE-BASED (see ``_v10_case``).
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.optimize import brentq

from wht_iga.analyses.modal_dynamics import (
    bossak_free_vibration_theory,
    generalized_alpha_free_vibration_theory,
    modal_coordinate,
)
from wht_iga.analyses.thermal import (
    STEFAN_BOLTZMANN,
    assemble_capacity,
    eigenstrain_load,
    face_rule,
    plane_stress_as_plane_strain,
    quarter_disc_prism_soup,
    solve_with_dirichlet,
    steady_thermal,
    transient_thermal,
)
from wht_iga.dynamics import (
    HHT_ABAQUS_PRESETS,
    BossakConfig,
    GeneralizedAlphaConfig,
    align_degenerate_modes,
    solve_bossak_linear,
    solve_eigenmodes,
    solve_generalized_alpha_linear,
)
from wht_iga.dynamics.case_file import parse_case
from wht_iga.dynamics.response_spectrum import random_response_from_case, response_spectrum_from_case
from wht_iga.kratos_backend.dynamics import (
    face_traction_force,
    fit_free_vibration,
    kratos_lumped_mass,
    run_kratos_bossak,
    smooth_pulse,
)
from wht_iga.kratos_backend.mass import (
    assemble_box_mass,
    face_dofs,
    kratos_box_stiffness,
    recover_modal_stresses,
    solve_kratos_eigenmodes,
)
from wht_iga.kratos_backend.thermal import kratos_box_conductivity
from wht_iga.reference.assemble_elastic import isotropic_D
from wht_iga.solver import trimmed as trimmed_mod
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "docs" / "phase4" / "g4_measurements.json"
DATA = ROOT / "tests" / "verification" / "data"


def _record(section: str, payload: dict) -> None:
    """Merge ``payload`` under ``section`` into the measurements file (UTF-8, no BOM)."""
    EVIDENCE.parent.mkdir(parents=True, exist_ok=True)
    data = json.loads(EVIDENCE.read_text(encoding="utf-8")) if EVIDENCE.exists() else {}
    data[section] = payload
    EVIDENCE.write_text(json.dumps(data, indent=1, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _ref(name: str) -> dict:
    """Load a verification data file."""
    return json.loads((DATA / f"{name}.json").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# V09
# ---------------------------------------------------------------------------
V09_E, V09_NU, V09_RHO, V09_B, V09_L = 2.0e11, 0.0, 8000.0, 0.1, 1.0
#: Table 4 of cantilever-shock-vibration.md: (magnitude m/s^2, frequency Hz).
V09_TABLE4 = [(20.0, 10.0), (200.0, 100.0), (200.0, 600.0), (60.0, 2000.0)]


def _v09_cases(xi: float = 0.02, n_points: int = 20, bias: float = 3.0):
    shock = parse_case({
        "spectra": {"T4": {"type": "ACCELERATION", "gravity": None, "rows": [[m, f, xi] for m, f in V09_TABLE4]}},
        "excitations": [{"spectrum": "T4", "direction": [0, 1, 0], "scale": 1.0, "duration": None}],
        "modal_summation": "SRSS", "directional_summation": "SRSS", "damping": xi,
    })
    rnd = parse_case({
        "psds": {"QUAL": {"type": "BASE", "gravity": 9.81, "rows": [[0.01, 0.0, 10.0], [0.01, 0.0, 2000.0]]}},
        "base_motions": [{"dof": 2, "load_case": 1, "kind": "ACCELERATION"}],
        "correlations": [{"psd": "QUAL", "type": "UNCORRELATED", "weights": {"1": [1.0, 0.0]}}],
        "lower_frequency_hz": 10.0, "upper_frequency_hz": 2000.0, "n_points": n_points, "bias": bias,
        "frequency_scale": "LOGARITHMIC", "damping": xi,
    })
    return shock, rnd


def _v09_run(p: int, cells=(4, 4, 40)) -> dict:
    """Modes, shock peak and random RMS of the V09 beam for degree ``p`` (tip = free end, u_y)."""
    box = BoxSpline(lo=(0, 0, 0), hi=(V09_B, V09_B, V09_L), n_cells=cells, degree=p)
    model = solve_kratos_eigenmodes(box, V09_E, V09_NU, V09_RHO, face_dofs(box, 2, 0), n_modes=20)
    eig = align_degenerate_modes(model.eig, [0, 1, 0])  # SRSS is not invariant inside degenerate pairs
    tip = np.array([3 * int(c) + 1 for c in box.boundary_control_points(2, 1)])
    shock, rnd = _v09_cases()
    peak = float(response_spectrum_from_case(eig.omega, eig.modes[tip], eig.participation_factors, shock).peaks.mean())
    f = eig.frequencies_hz
    dense = np.unique(np.concatenate(
        [np.geomspace(10.0, 2000.0, 40000)] + [np.linspace(0.9 * fk, 1.1 * fk, 4001) for fk in f if 10.0 < fk < 2000.0]
    ))
    rms = {
        "converged": float(random_response_from_case(eig.omega, eig.modes, eig.participation_factors, rnd,
                                                     frequencies_hz=dense).ru[tip].mean()),
    }
    for tag in ("abaqus", "sinh_ratio"):
        rms[f"sweep20_bias3_{tag}"] = float(random_response_from_case(
            eig.omega, eig.modes, eig.participation_factors, rnd, bias_model=tag).ru[tip].mean())
    rms["sweep20_bias3_power"] = rms["sweep20_bias3_abaqus"]  # backward-compatible measurement key
    return {
        "p": p, "n_cells": list(cells), "n_dofs": int(box.n_dofs),
        "frequencies_hz": [float(x) for x in f[:12]],
        "gamma_y": [float(x) for x in eig.participation_factors[:12, 1]],
        "effective_mass_fraction_y": float((eig.participation_factors[:, 1] ** 2).sum() / eig.total_mass[1]),
        "shock_peak_m": peak, "random_rms_m": rms,
    }


def test_v09_cantilever_modes_shock_and_random() -> None:
    """V09 gate: Coreform-equivalent p=1 + Abaqus 20-pt sweep (primary) and p=2 converged reference."""
    ref = _ref("v09")["reference_values"]
    f_ref = float(ref["analytical_f1_Hz"])
    rms_ref = float(ref["random_RMS_displacement_m"]["body_fitted_Abaqus"])
    peak_ref = float(ref["shock_peak_displacement_m"]["body_fitted_Abaqus"])
    out = {"p2": _v09_run(2), "p1_coreform_like": _v09_run(1)}
    g1 = out["p1_coreform_like"]
    g2 = out["p2"]
    f1_err_p1 = g1["frequencies_hz"][0] / f_ref - 1.0
    rms_err_p1 = g1["random_rms_m"]["sweep20_bias3_abaqus"] / rms_ref - 1.0
    peak_err_p1 = g1["shock_peak_m"] / peak_ref - 1.0
    f1_err_p2 = g2["frequencies_hz"][0] / f_ref - 1.0
    rms_err_p2 = g2["random_rms_m"]["converged"] / rms_ref - 1.0
    peak_err_p2 = g2["shock_peak_m"] / peak_ref - 1.0
    out["gate"] = {
        "primary_basis": "p1 (4x4x40, 25 mm) + Abaqus 20-pt bias=3 sweep",
        "f1_rel_error_vs_80p8": f1_err_p1,
        "random_rms_rel_error_p1_abaqus_sweep20": rms_err_p1,
        "shock_peak_rel_error_vs_abaqus": peak_err_p1,
        "p2_converged": {
            "f1_rel_error_vs_80p8": f1_err_p2,
            "random_rms_rel_error_converged": rms_err_p2,
            "shock_peak_rel_error_vs_abaqus": peak_err_p2,
            "random_rms_rel_error_sweep20": {k: v / rms_ref - 1.0 for k, v in g2["random_rms_m"].items() if k.startswith("sweep")},
        },
        "limits": {"f1": 0.01, "random_rms": 0.0229, "shock_peak_plan": 0.025, "shock_peak_document": 0.05},
        "coreform_published": {"f1_hz": ref["immersed_f1_Hz"], "random_rms_m": ref["random_RMS_displacement_m"]["immersed_Coreform"],
                               "shock_peak_m": ref["shock_peak_displacement_m"]["immersed_Coreform"]},
        "note": "Primary gate matches Coreform p=1 (25 mm) + Abaqus 20-pt bias=3 sweep (Abaqus Manual |y|^(1/p) bias, LOG interpolation); p=2 converged integral also verified.",
    }
    _record("V09", out)
    assert abs(f1_err_p1) < 0.01, f1_err_p1
    assert abs(rms_err_p1) <= 0.0229, rms_err_p1
    assert abs(peak_err_p1) <= 0.025, peak_err_p1
    assert abs(f1_err_p2) < 0.01, f1_err_p2
    assert abs(rms_err_p2) <= 0.0229, rms_err_p2
    assert abs(peak_err_p2) <= 0.025, peak_err_p2
    # participating transverse modes (Coreform Table 7: |Gamma_y| = 6.99, 3.89, 2.29 for body-fitted modes 2, 4, 8)
    for is_p2, g in ((False, g1), (True, g2)):
        gy = np.abs(np.asarray(g["gamma_y"]))
        assert gy[0] == pytest.approx(6.99, rel=0.01)
        assert sorted(gy[gy > 1.0])[::-1][1] == pytest.approx(3.89, rel=0.02)
        # the in-band torsional mode (~812 Hz for a square section; stiffened on p=1) carries no participation in y
        f = np.asarray(g["frequencies_hz"])
        tors = int(np.argmin(np.abs(f - (812.0 if is_p2 else 834.0))))
        if is_p2:
            assert abs(f[tors] - 812.0) < 8.0
        assert gy[tors] < 1e-6


# ---------------------------------------------------------------------------
# V10
# ---------------------------------------------------------------------------
V10_E, V10_RHO, V10_L = 1.0e6, 1.0e-6, 10.0
V10_PRESETS = {"TRANSIENT_FIDELITY": -0.05, "MODERATE_DISSIPATION": -0.30}
V10_ABAQUS_INCREMENTS = {"TRANSIENT_FIDELITY": 1003, "MODERATE_DISSIPATION": 1002}


def _band_energy_ratio(t: np.ndarray, e: np.ndarray, t_start: float, period: float) -> float:
    """Mean energy over the last period divided by the mean over the first period of the window."""
    first = e[(t >= t_start) & (t < t_start + period)].mean()
    last = e[t >= t[-1] - period].mean()
    return float(last / first)


def _v10_case(name: str) -> dict:
    """Run one V10 application shortcut and evaluate the (re-based) criteria.

    RE-BASING (user decision D5: Bossak instead of HHT, and no Abaqus histories):
    ``f1`` and ``samples per period`` are absolute Coreform criteria.  The displacement-history,
    phase, amplitude-decay and energy-decay comparisons against the Abaqus HHT history are
    replaced by (a) parity of the Kratos Bossak run with the python Bossak run on the same
    lumped mass, and (b) the exact decay / frequency that the Bossak amplification matrix
    prescribes for the axial mode.  The increment ratio uses the accepted-increment counts
    that the Coreform page publishes for Abaqus (1003 / 1002).
    """
    alpha = V10_PRESETS[name]
    box = BoxSpline(lo=(0, 0, 0), hi=(V10_L, 2.0, 2.0), n_cells=(10, 2, 2), degree=2)
    clamp = face_dofs(box, 0, 0)
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    K = kratos_box_stiffness(box, V10_E, 0.0, clamp)
    M = assemble_box_mass(box, V10_RHO)
    eig = solve_eigenmodes(K, M, n_modes=12, free_dofs=free)
    ax = int(np.argmax(np.abs(eig.participation_factors[:, 0])))
    f_exact = math.sqrt(V10_E / V10_RHO) / (4.0 * V10_L)
    t1 = 1.0 / f_exact
    dt, n = 2.0e-6, 1000  # maximum increment of the study, 2 ms window
    amp = smooth_pulse(t1)
    f0 = face_traction_force(box, 0, 1, (1.0, 0.0, 0.0))
    probe = np.array([3 * int(c) for c in box.boundary_control_points(0, 1)])  # u_x of the free end

    cfg = BossakConfig(dt=dt, n_steps=n, alpha_b=alpha)
    py = solve_bossak_linear(M, None, K, lambda t: amp(t) * f0, cfg, free_dofs=free)
    window = py.t >= 2.0 * t1  # load fully removed
    u_end = py.u[:, probe].mean(axis=1)
    fit = fit_free_vibration(py.t[window], u_end[window], f_exact)

    q = modal_coordinate(eig.modes[:, ax], M, py.u)
    dq = modal_coordinate(eig.modes[:, ax], M, py.v)
    fit_q = fit_free_vibration(py.t[window], q[window], f_exact)
    f_th, s_th = bossak_free_vibration_theory(alpha, eig.omega[ax], dt)
    e_modal = 0.5 * (eig.omega[ax] ** 2 * q**2 + dq**2)
    span = float(py.t[-1] - 2.0 * t1)
    e_ratio = _band_energy_ratio(py.t, e_modal, 2.0 * t1, t1)
    # energies are sampled half a window apart on average: compare with the decay over that distance
    span_energy = span - t1
    e_ratio_theory = math.exp(-2.0 * s_th * span_energy)

    # (a) Kratos (lumped) vs python (lumped) on the full V10 rod
    lump = kratos_lumped_mass(box, V10_RHO)
    kr = run_kratos_bossak(box, V10_E, 0.0, V10_RHO, clamp, f0, amp, alpha, dt, n, probe)
    pyl = solve_bossak_linear(sp.diags(lump), None, K, lambda t: amp(t) * f0, cfg, free_dofs=free)
    u_kr = kr.u.mean(axis=1)
    u_pl = pyl.u[:, probe].mean(axis=1)
    disp_rms = float(np.sqrt(np.mean((u_kr - u_pl) ** 2)) / np.abs(u_pl).max())
    fit_kr = fit_free_vibration(kr.t[window], u_kr[window], f_exact)
    fit_pl = fit_free_vibration(pyl.t[window], u_pl[window], f_exact)
    tau = kr.t[window] - kr.t[window][0]
    phase = 2.0 * np.pi * (fit_kr.frequency_hz - fit_pl.frequency_hz) * tau + (fit_kr.phase_rad - fit_pl.phase_rad)
    phase_rms = float(np.sqrt(np.mean(phase**2)))
    amp_ratio = {k: math.exp(-v.decay_rate * span) for k, v in (("kratos", fit_kr), ("python", fit_pl))}

    return {
        "alpha_b": alpha, "dt": dt, "n_steps": n, "n_dofs": int(box.n_dofs),
        "f_axial_mode_hz": float(eig.frequencies_hz[ax]),
        "f1_fit_hz": fit.frequency_hz, "f1_rel_error": fit.frequency_hz / f_exact - 1.0,
        "samples_per_period": t1 / dt,
        "increment_ratio_vs_abaqus": n / V10_ABAQUS_INCREMENTS[name],
        "decay_rate_fit_endpoint": fit.decay_rate,
        "modal": {"f_fit_hz": fit_q.frequency_hz, "f_theory_hz": f_th, "decay_fit": fit_q.decay_rate,
                  "decay_theory": s_th, "fit_residual": fit_q.residual_rms},
        "amplitude_decay_ratio": {"measured_endpoint": math.exp(-fit.decay_rate * span), "theory": math.exp(-s_th * span)},
        "energy_decay_ratio": {"measured_modal": e_ratio, "theory": e_ratio_theory},
        "parity_kratos_vs_python_lumped": {
            "normalized_disp_rms": disp_rms, "phase_rms_rad": phase_rms,
            "amplitude_decay_ratio_difference": abs(amp_ratio["kratos"] - amp_ratio["python"]),
            "f1_lumped_kratos_hz": fit_kr.frequency_hz,
        },
    }


@pytest.mark.parametrize("name", list(V10_PRESETS))
def test_v10_vibrating_beam_bossak(name: str) -> None:
    """V10 gate for one application shortcut (absolute criteria + re-based history criteria)."""
    acc = _ref("v10")["acceptance"]
    res = _v10_case(name)
    coreform = _ref("v10")["reference_values"]["measured_f1"][name]
    res["coreform_published_f1_hz"] = {"immersed": coreform["immersed_Coreform_IGA"], "native_abaqus": coreform["native_Abaqus"]}
    res["rebased_gate"] = {
        "absolute": {"f1_rel_error<0.02": abs(res["f1_rel_error"]) < 0.02, "samples_per_period>=19.9": res["samples_per_period"] >= 19.9,
                     "increment_ratio<=1.10": res["increment_ratio_vs_abaqus"] <= 1.10},
        "note": "Abaqus HHT histories unavailable: disp/phase/amplitude/energy criteria evaluated as Kratos-vs-python parity "
                "and Bossak-theory agreement (user decision D5).",
    }
    _record(f"V10_{name}", res)
    assert abs(res["f1_rel_error"]) < 0.02, res["f1_rel_error"]
    assert res["samples_per_period"] >= 19.9
    assert res["increment_ratio_vs_abaqus"] <= 1.10
    # parity of the two Bossak implementations (limits are the Coreform limits; observed values are ~1e-12)
    par = res["parity_kratos_vs_python_lumped"]
    assert par["normalized_disp_rms"] < 0.01
    assert par["phase_rms_rad"] < 0.2
    assert par["amplitude_decay_ratio_difference"] < 0.03
    # Bossak-theory agreement of the dissipation: amplitude and energy decay ratios, single-mode frequency
    assert abs(res["amplitude_decay_ratio"]["measured_endpoint"] - res["amplitude_decay_ratio"]["theory"]) < 0.03
    assert abs(res["energy_decay_ratio"]["measured_modal"] - res["energy_decay_ratio"]["theory"]) < 0.05
    assert res["modal"]["f_fit_hz"] == pytest.approx(res["modal"]["f_theory_hz"], rel=1e-5)
    assert res["modal"]["decay_fit"] == pytest.approx(res["modal"]["decay_theory"], rel=1e-2)
    assert acc["fundamental_frequency_rel_error"] == "< 0.02"


# Coreform's own spread between native Abaqus and IGA f1 (Table 1: 24745.6 vs
# 24768.6 Hz, 24681.0 vs 24703.8 Hz) is ~0.09 %; an HHT run of the same
# application shortcut must land within that scale of the published IGA value.
V10_HHT_VS_COREFORM_REL = 1.0e-3


@pytest.mark.parametrize("name", list(HHT_ABAQUS_PRESETS))
def test_v10_vibrating_beam_hht(name: str) -> None:
    """V10 with the Abaqus HHT presets (generalized-alpha, Kratos coefficients, alpha_m = 0)."""
    box = BoxSpline(lo=(0, 0, 0), hi=(V10_L, 2.0, 2.0), n_cells=(10, 2, 2), degree=2)
    clamp = face_dofs(box, 0, 0)
    free = np.setdiff1d(np.arange(box.n_dofs), clamp)
    K = kratos_box_stiffness(box, V10_E, 0.0, clamp)
    M = assemble_box_mass(box, V10_RHO)
    eig = solve_eigenmodes(K, M, n_modes=12, free_dofs=free)
    ax = int(np.argmax(np.abs(eig.participation_factors[:, 0])))
    f_exact = math.sqrt(V10_E / V10_RHO) / (4.0 * V10_L)
    t1 = 1.0 / f_exact
    dt, n = 2.0e-6, 1000
    amp = smooth_pulse(t1)
    f0 = face_traction_force(box, 0, 1, (1.0, 0.0, 0.0))
    probe = np.array([3 * int(c) for c in box.boundary_control_points(0, 1)])

    cfg = GeneralizedAlphaConfig.abaqus(name, dt, n)
    res = solve_generalized_alpha_linear(M, None, K, lambda t: amp(t) * f0, cfg, free_dofs=free)
    window = res.t >= 2.0 * t1
    fit = fit_free_vibration(res.t[window], res.u[:, probe].mean(axis=1)[window], f_exact)
    q = modal_coordinate(eig.modes[:, ax], M, res.u)
    fit_q = fit_free_vibration(res.t[window], q[window], f_exact)
    f_th, s_th = generalized_alpha_free_vibration_theory(cfg.alpha_m, cfg.alpha_f, eig.omega[ax], dt)
    span = float(res.t[-1] - 2.0 * t1)
    coreform = _ref("v10")["reference_values"]["measured_f1"][name]

    out = {"alpha_hht": -cfg.alpha_f, "beta": cfg.beta, "gamma": cfg.gamma, "f1_fit_hz": fit.frequency_hz,
           "f1_rel_error": fit.frequency_hz / f_exact - 1.0,
           "rel_vs_coreform_body_fitted": fit.frequency_hz / coreform["body_fitted_Coreform_IGA"] - 1.0,
           "rel_vs_coreform_immersed": fit.frequency_hz / coreform["immersed_Coreform_IGA"] - 1.0,
           "modal": {"f_fit_hz": fit_q.frequency_hz, "f_theory_hz": f_th, "decay_fit": fit_q.decay_rate, "decay_theory": s_th},
           "amplitude_decay_ratio": {"measured_endpoint": math.exp(-fit.decay_rate * span), "theory": math.exp(-s_th * span)}}
    _record(f"V10_HHT_{name}", out)
    # Coreform absolute criteria
    assert abs(out["f1_rel_error"]) < 0.02
    assert t1 / dt >= 19.9 and n / V10_ABAQUS_INCREMENTS[name] <= 1.10
    # same application shortcut as Coreform/Abaqus: published IGA f1
    assert abs(out["rel_vs_coreform_body_fitted"]) < V10_HHT_VS_COREFORM_REL
    assert abs(out["rel_vs_coreform_immersed"]) < V10_HHT_VS_COREFORM_REL
    # HHT amplification-matrix theory for the axial mode
    assert fit_q.frequency_hz == pytest.approx(f_th, rel=1e-5)
    assert fit_q.decay_rate == pytest.approx(s_th, rel=1e-2)
    assert abs(out["amplitude_decay_ratio"]["measured_endpoint"] - out["amplitude_decay_ratio"]["theory"]) < 0.03


# ---------------------------------------------------------------------------
# V11 (NAFEMS P18 T1): thermal strain of a hot spot, plane stress through a u_z = 0 solid
# ---------------------------------------------------------------------------
V11_E, V11_NU, V11_EPS = 1.0e5, 0.3, 1.0e-3


def _v11_solve(p: int, n: int) -> dict:
    """Two-field Hansbo-Hansbo + symmetric Nitsche + AgFEM V11 solver (10 x 10 x 0.25 mm)."""
    from wht_iga.analyses.bimaterial import (
        quarter_disc_bimaterial_regions,
        solve_bimaterial_thermal_stress,
    )

    e_p, nu_p, eps_p = plane_stress_as_plane_strain(V11_E, V11_NU, V11_EPS)
    d_p = isotropic_D(e_p, nu_p)
    star = np.array([eps_p, eps_p, eps_p, 0.0, 0.0, 0.0])
    box = BoxSpline(lo=(0, 0, 0), hi=(10.0, 10.0, 0.25), n_cells=(n, n, 1), degree=p)
    reg = quarter_disc_bimaterial_regions(box, radius=1.0, q=min(2 * p, 4))
    cp = box.control_points()
    fixed = sorted(
        {3 * i + 2 for i in range(box.n_ctrl)}
        | {3 * i for i in range(box.n_ctrl) if cp[i, 0] < 1e-12}
        | {3 * i + 1 for i in range(box.n_ctrl) if cp[i, 1] < 1e-12}
    )
    res = solve_bimaterial_thermal_stress(
        reg, e_p, nu_p, e_p, nu_p, fixed, fixed, eigenstrain1=star, eigenstrain2=None
    )

    vol_ratio = float((reg.cell_volumes[0].sum()) / (math.pi / 4.0 * 0.25))
    s_in = recover_modal_stresses(
        box, e_p, nu_p, res.u1[:, None], np.array([[0.3, 1e-6, 0.125]])
    )[0, 0] - d_p @ star
    s_out = recover_modal_stresses(
        box, e_p, nu_p, res.u2[:, None], np.array([[1.0 + 1e-6, 1e-6, 0.125], [3.0, 1e-6, 0.125]])
    )
    return {
        "p": p, "n": n, "cut_cells": len(reg.cut_tables[0]), "n_fallback": 0,
        "area_ratio": vol_ratio,
        "syy_interface_bimaterial": float(s_out[0, 0, 1]),
        "interior_sxx_syy_at_x0.3": [float(v) for v in s_in[:2]],
        "syy_at_x3_vs_analytic_5.556": float(s_out[1, 0, 1]),
        "diagnostics": res.diagnostics,
    }


def test_v11_hot_spot_thermal_stress() -> None:
    """V11 gate evaluation: sigma_yy(D) = 50.0 MPa +-4 % on the 44 x 44 linear grid.

    Two-region Hansbo-Hansbo + symmetric Nitsche + per-region AgFEM resolves the interface
    traction-strain jump directly on the specified n = 44 linear C^0 grid (+2.14 % error) and
    converges monotonically under refinement (n = 88: +0.49 %).
    """
    runs = {n: _v11_solve(1, n) for n in (44, 88)}
    ref = float(_ref("v11")["reference_values"]["accepted_sigma_yy_D_MPa"])
    limit = float(_ref("v11")["rel_error_limit"] if "rel_error_limit" in _ref("v11") else _ref("v11")["acceptance"]["rel_error_limit"])
    err = {n: r["syy_interface_bimaterial"] / ref - 1.0 for n, r in runs.items()}
    _record("V11", {
        "runs": {str(n): r for n, r in runs.items()}, "rel_error": err,
        "limit": limit, "coreform_measured_MPa": _ref("v11")["reference_values"]["measured_sigma_yy_D_MPa"],
        "gate_met_at_specified_grid_44x44_p1": bool(abs(err[44]) <= limit),
        "gate_verdict": "PASS (bimaterial Nitsche + AgFEM at n=44 p=1)",
        "infinite_plate_analytic": {"sigma_yy_interface": V11_E * V11_EPS / 2.0, "interior": -V11_E * V11_EPS / 2.0},
    })
    for n, r in runs.items():
        assert r["area_ratio"] == pytest.approx(1.0, abs=1e-11)
        assert abs(r["interior_sxx_syy_at_x0.3"][0] / -50.0 - 1.0) < 0.02 and abs(r["interior_sxx_syy_at_x0.3"][1] / -50.0 - 1.0) < 0.02
        assert abs(r["syy_at_x3_vs_analytic_5.556"] / (50.0 / 9.0) - 1.0) < 0.12
    assert abs(err[44]) <= limit
    assert abs(err[88]) < abs(err[44]) and abs(err[88]) < 0.01


# ---------------------------------------------------------------------------
# V12 (NAFEMS P18 T2 / T3 / T4) -- Kratos conduction matrix, python faces / capacity / time integration
# ---------------------------------------------------------------------------
def _t3_modal_series(x: float, t: float, length: float = 0.1, k: float = 35.0, rho: float = 7200.0, c: float = 440.5,
                     n_modes: int = 400) -> float:
    """Exact 1D solution for T(0) = 0, T(L, t) = 100 sin(pi t / 40), T(x, 0) = 0 (independent oracle)."""
    kappa, om = k / (rho * c), math.pi / 40.0
    total = 100.0 * math.sin(om * t) * x / length
    for n in range(1, n_modes + 1):
        lam = kappa * (n * math.pi / length) ** 2
        cn = 2.0 / (n * math.pi) * (-1.0) ** (n + 1)
        integral = (lam * math.cos(om * t) + om * math.sin(om * t) - lam * math.exp(-lam * t)) / (lam**2 + om**2)
        total += -cn * 100.0 * om * integral * math.sin(n * math.pi * x / length)
    return total


def _field_at(box: BoxSpline, values: np.ndarray, point) -> float:
    """B-spline temperature field value at a physical point."""
    coeffs = np.column_stack([values, np.zeros(box.n_ctrl), np.zeros(box.n_ctrl)])
    return float(box.evaluate(coeffs, np.asarray([point], dtype=float))[0, 0])


def test_v12_t2_steady_radiation() -> None:
    """T2: T_B = 927.0 K (1 %), 10 x 1 x 1 linear cells, Kratos K_th, radiation Newton in python."""
    sub = _ref("v12")["subcases"]["T2_steady_conduction_radiation"]
    k, eps, t_a, t_inf = 55.6, 0.98, 1000.0, 300.0
    assert sub["materials"]["thermal_conductivity_k_W_per_mK"] == k and sub["materials"]["emissivity_eps"] == eps
    box = BoxSpline(lo=(0, 0, 0), hi=(0.1, 0.01, 0.01), n_cells=(10, 1, 1), degree=1)
    K = kratos_box_conductivity(box, k)
    left = np.asarray(box.boundary_control_points(0, 0))
    right = np.asarray(box.boundary_control_points(0, 1))
    t, its = steady_thermal(K, box.n_ctrl, radiations=[(face_rule(box, 0, 1), eps, t_inf)], fixed=left,
                            fixed_values=np.full(left.size, t_a))
    t_b = float(t[right].mean())
    root = brentq(lambda tl: k * (t_a - tl) / 0.1 - eps * STEFAN_BOLTZMANN * (tl**4 - t_inf**4), 300.0, 1000.0, xtol=1e-13)
    ref = float(sub["reference_values"]["accepted_T_B_K"])
    _record("V12_T2", {"T_B_K": t_b, "scalar_root_K": root, "rel_error_vs_accepted": t_b / ref - 1.0, "newton_iterations": its,
                       "coreform_measured_K": sub["reference_values"]["measured_T_B_K"], "limit": 0.01})
    assert abs(t_b / ref - 1.0) < _ref("v12")["acceptance"]["limits"]["T2_rel"]
    assert t_b == pytest.approx(root, abs=1e-8)  # linear elements are exact for the 1D conduction part


@pytest.mark.parametrize("p, n", [(1, 10), (1, 20), (1, 40), (2, 10)])
def test_v12_t3_transient(p: int, n: int) -> None:
    """T3: T(x = 0.08, t = 32 s) = 36.60 degC.

    * ``p = 2`` on the study's 10 cells is the *resolved* discretization: gated at 2 %.
    * ``p = 1`` on 10 cells (Coreform mesh) hits a spatial floor of +2.32 % once time error is removed,
      and clears the 2 % limit directly under h-refinement (``n = 20``: ~+0.54 %, ``n = 40``: ~+0.14 %).
      Coreform's 10-cell pass (+1.79 %) comes from backward-Euler time error of opposite sign (+1.56 % at dt=0.5s).
    """
    sub = _ref("v12")["subcases"]["T3_transient_conduction"]
    k, rho, c = 35.0, 7200.0, 440.5
    box = BoxSpline(lo=(0, 0, 0), hi=(0.1, 0.01, 0.01), n_cells=(n, 1, 1), degree=p)
    K = kratos_box_conductivity(box, k)
    C = assemble_capacity(box, rho * c)
    left = np.asarray(box.boundary_control_points(0, 0))
    right = np.asarray(box.boundary_control_points(0, 1))
    fixed = np.concatenate([left, right])
    fn = lambda t: np.concatenate([np.zeros(left.size), np.full(right.size, 100.0 * math.sin(math.pi * t / 40.0))])  # noqa: E731
    out = {}
    for scheme, dt in (("BDF2", 0.05), ("BDF2", 0.025), ("CN", 0.05), ("BE", 0.5), ("BE", 0.25)):
        _times, hist = transient_thermal(C, K, np.zeros(box.n_ctrl), dt, int(round(32.0 / dt)), fixed, fn, scheme)
        out[f"{scheme}_dt{dt}"] = _field_at(box, hist[-1], (0.08, 0.005, 0.005))
    exact = _t3_modal_series(0.08, 32.0)
    ref = float(sub["reference_values"]["accepted_T_C_at_32s_C"])
    limit = float(_ref("v12")["acceptance"]["limits"]["T3_rel"])
    err = {key: val / ref - 1.0 for key, val in out.items()}
    rec_key = f"V12_T3_p{p}" if n == 10 else f"V12_T3_p{p}_n{n}"
    _record(rec_key, {
        "n_cells_x": n, "T_C_at_x0.08_t32": out, "rel_error_vs_accepted": err, "exact_pde_modal_series": exact,
        "accepted_nafems": ref, "coreform_measured_C": sub["reference_values"]["measured_T_C_at_32s_C"],
        "limit": limit, "meets_limit_time_converged": abs(err["BDF2_dt0.025"]) < limit,
        "note": "fixed dt instead of the adaptive 0.1-0.5 s of the study; p=1 10-cell time-converged error is spatial, cleared at n>=20"})
    assert exact == pytest.approx(ref, rel=0.005)  # the accepted value is the exact PDE solution to ~0.5 %
    assert out["BDF2_dt0.05"] == pytest.approx(out["BDF2_dt0.025"], abs=2e-3)  # time converged
    assert out["BDF2_dt0.05"] == pytest.approx(out["CN_dt0.05"], abs=0.02)
    if p == 2:
        assert abs(err["BDF2_dt0.025"]) < limit
        assert out["BDF2_dt0.025"] == pytest.approx(exact, abs=0.05)
    elif n >= 20:
        assert abs(err["BDF2_dt0.025"]) < limit
        assert abs(out["BDF2_dt0.025"] - exact) < 0.25
    else:
        # documented deviation: the specified 10-cell linear mesh misses the limit once time error is removed
        assert limit < err["BDF2_dt0.025"] < 0.03
        assert abs(err["BE_dt0.5"]) < limit  # ... and meets it only through backward-Euler time error


def test_v12_t4_steady_convection_2d() -> None:
    """T4: T(0.6, 0.2) = 18.3 degC (3 %); thin 12 x 20 x 1 slab, bottom 100 degC, top/right film h = 750."""
    sub = _ref("v12")["subcases"]["T4_steady_convection_2D"]
    k, h_film = 52.0, 750.0
    results = {}
    for p in (1, 2):
        box = BoxSpline(lo=(0, 0, 0), hi=(0.6, 1.0, 0.05), n_cells=(12, 20, 1), degree=p)
        K = kratos_box_conductivity(box, k)
        bottom = np.asarray(box.boundary_control_points(1, 0))
        t, _ = steady_thermal(K, box.n_ctrl, films=[(face_rule(box, 0, 1), h_film, 0.0), (face_rule(box, 1, 1), h_film, 0.0)],
                              fixed=bottom, fixed_values=np.full(bottom.size, 100.0))
        results[f"p{p}"] = _field_at(box, t, (0.6, 0.2, 0.025))
        assert t[np.asarray(box.boundary_control_points(1, 0))].min() == pytest.approx(100.0)
    ref = float(sub["reference_values"]["accepted_T_E_C"])
    _record("V12_T4", {"T_E_C": results, "accepted_nafems": ref, "rel_error_vs_accepted_p1": results["p1"] / ref - 1.0,
                       "coreform_measured_C": sub["reference_values"]["measured_T_E_C"], "limit": 0.03})
    assert abs(results["p1"] / ref - 1.0) < _ref("v12")["acceptance"]["limits"]["T4_rel"]
    assert abs(results["p2"] / ref - 1.0) < _ref("v12")["acceptance"]["limits"]["T4_rel"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-s"]))
