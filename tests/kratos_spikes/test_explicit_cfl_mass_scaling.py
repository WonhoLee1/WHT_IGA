# -*- coding: utf-8 -*-
"""CFL and mass-scaling execution tests (Track 2 open item 3).

Small bar: 4 cells, h = 1, E = 100, rho = 1 (c = 10), y/z fixed, x fixed at
x = 0. Kratos ``MechanicalExplicitStrategy`` + ``ExplicitCentralDifferencesScheme``.

``estimate_dt_crit`` is applied to the system Kratos actually integrates:
- K from the Kratos element LHS, not from the static builder. The builder
  replaces near-zero diagonals (``CheckAndCorrectZeroDiagonalValues``). Measured
  2026-11-19 for p = 3, eta = 1e-6: the last x-layer rows are <= 1.5e-18, and
  the builder K has 19.0 on that diagonal where the element K has 3.7e-30.
- M is the Kratos NODAL_MASS read back, or the HRZ values written by
  ``overwrite_nodal_mass``. With ties, K_eff = T^T K T on the masters.

References:
- p = 1 uncut: 1D fixed-free lumped bar, w_max = (2 c / h) cos(pi / (4 n)).
- Execution: from rest, central differences give u_n = u_0 cos(n theta) per mode.
  This uses the half-step start v_{1/2} = (dt/2) a_0. Kratos starts that way:
  the V-EXP 1D reference with this start matches Kratos to round-off.
  With E_n = 1/2 u_n K u_n + 1/2 v_{n-1/2} M v_{n-1/2}, this gives
  E_n / E_0 <= 1 + w dt / 2 < 2 while w_max dt < 2. Above the limit the energy
  grows geometrically (gate: >= 1e6 within 100 steps).
- S24 claim: cut-cell dt_crit with extension >= 0.5 x uncut dt_crit.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.kratos_spikes.test_s24_explicit_dynamics import (
    _assemble_trimmed_consistent_mass,
)
from tests.kratos_spikes.test_s24b_hybrid_explicit import _build_wave_model
from wht_iga.dynamics import explicit as EX
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_NX, _H, _RHO = 4, 1.0, 1.0
_C = 10.0  # sqrt(E / rho) with the S24b model constants (E = 100, rho = 1)
_STABLE_BOUND = 2.0  # leapfrog energy bound below dt_crit (see module doc)
_BLOWUP = 1e6
_N_STABLE, _N_UNSTABLE = 500, 100
_COUNTER = [0]


def _kratos():
    mods = load_kratos()
    return mods["KratosMultiphysics"], mods["StructuralMechanicsApplication"]


def _bar(KM, SM, degree, eta, ties):
    """Explicit bar model (+ optional extension ties on the last x-layer)."""
    _COUNTER[0] += 1
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(_NX * _H, 1.0, 1.0),
                    n_cells=(_NX, 1, 1), degree=degree)
    mh = _build_wave_model(KM, SM, box, eta, ties or eta < 1.0, f"CFL{_COUNTER[0]}")
    mp = mh["mp"]
    info = EX.build_x_extension(box)
    if ties:
        nu = box.n_cp[0]
        src = [int(s) for s in info["src"]]
        w = np.array([[float(v) for v in info["masters"][0][1]]])
        cid = 100000
        for kc in range(box.n_cp[2]):
            for jc in range(box.n_cp[1]):
                row = nu * (jc + box.n_cp[1] * kc)
                add_linear_tie(KM, mp, cid,
                               [mp.GetNode(row + nu).GetDof(KM.DISPLACEMENT_X)],
                               [mp.GetNode(row + s + 1).GetDof(KM.DISPLACEMENT_X) for s in src],
                               w)
                cid += 1
    scheme = SM.ExplicitCentralDifferencesScheme(1.0, 1.0 / 3.0, 0)
    strat = SM.MechanicalExplicitStrategy(mp, scheme, False, False, False)
    strat.SetEchoLevel(0)
    return box, mh, strat, info


def _element_k_x(KM, mp, n_ctrl):
    """x-block of K assembled from the Kratos element LHS (call after Initialize)."""
    K = np.zeros((n_ctrl, n_ctrl))
    for el in mp.Elements:
        lhs = KM.Matrix()
        el.CalculateLeftHandSide(lhs, mp.ProcessInfo)
        ids = [nd.Id - 1 for nd in el.GetNodes()]
        K[np.ix_(ids, ids)] += np.array(lhs)[0::3, 0::3]
    return K


def _system(KM, SM, degree, eta, mass):
    """Effective explicit system. mass: 'kratos' | 'hrz' | 'hrz_ext' (ties + HRZ(T^T M T))."""
    ties = mass == "hrz_ext"
    box, mh, strat, info = _bar(KM, SM, degree, eta, ties)
    strat.Initialize()
    n, nu = box.n_ctrl, box.n_cp[0]
    K = _element_k_x(KM, mh["mp"], n)
    fixed = np.array([I % nu == 0 for I in range(n)])
    total = _RHO * (_NX - 1 + eta) * _H
    if mass == "kratos":
        m = EX.read_nodal_mass(KM, mh["mp"], [I + 1 for I in range(n)])
    else:
        M_c = _assemble_trimmed_consistent_mass(box, rho=_RHO, eta_last_x=eta)
        if ties:
            T = info["T"]
            m = EX.hrz_lump_diagonal(T.T @ M_c @ T, total_mass=total)[0::3]
        else:
            m = EX.hrz_lump_diagonal(M_c, total_mass=total)[0::3]
    if ties:
        Tx = info["T"][0::3, 0::3]
        lex = [I for I in range(n) if I % nu != nu - 1]
    else:
        Tx = np.eye(n)
        lex = list(range(n))
    return {"degree": degree, "eta": eta, "ties": ties, "overwrite": mass != "kratos",
            "K": K, "Keff": Tx.T @ K @ Tx, "m": np.asarray(m, dtype=float),
            "free": ~fixed[lex], "Tx": Tx, "lex": lex}


def _dt_crit(S, m=None):
    return EX.estimate_dt_crit(S["Keff"], S["m"] if m is None else m, S["free"])[0]


def _max_energy_ratio(KM, SM, S, dt, n_steps, m=None, seed=0, return_state=False):
    """Run a fresh Kratos model of system S from rest with a random u_0; max E_n / E_0."""
    m = S["m"] if m is None else np.asarray(m, dtype=float)
    box, mh, strat, _ = _bar(KM, SM, S["degree"], S["eta"], S["ties"])
    mp = mh["mp"]
    lex, free, K = S["lex"], S["free"], S["K"]
    rng = np.random.default_rng(seed)
    um = np.zeros(len(lex))
    um[free] = 1e-3 * rng.standard_normal(int(free.sum()))
    u_full = S["Tx"] @ um  # slaves consistent with the tie
    nodes = [mp.GetNode(I + 1) for I in range(box.n_ctrl)]
    for I, nd in enumerate(nodes):
        nd.SetSolutionStepValue(KM.DISPLACEMENT_X, float(u_full[I]))
    strat.Initialize()
    if S["overwrite"]:
        EX.overwrite_nodal_mass(KM, mp, {lex[i] + 1: float(m[i]) for i in range(len(lex))})
    e0 = 0.5 * float(um @ S["Keff"] @ um)
    ratios = []

    def energy(_mp):
        u = np.array([nd.GetSolutionStepValue(KM.DISPLACEMENT_X) for nd in nodes])
        v = np.array([nd.GetSolutionStepValue(SM.MIDDLE_VELOCITY_X) for nd in nodes])[lex]
        with np.errstate(over="ignore", invalid="ignore"):  # overflow = the blow-up being measured
            e = 0.5 * float(u @ K @ u) + 0.5 * float(np.sum(m[free] * v[free] ** 2))
        ratios.append(e / e0)
        return e / e0

    try:
        EX.run_explicit_history(KM, mp, strat, dt, n_steps, probe_fn=energy)
    except RuntimeError:  # non-finite energy = blow-up
        ratios.append(np.inf)
    r = float(np.max(ratios))
    if return_state:
        u_end = np.array([nd.GetSolutionStepValue(KM.DISPLACEMENT_X) for nd in nodes])
        return r, u_full, u_end
    return r


def test_dt_crit_matches_1d_lumped_bar_p1():
    """p = 1 uncut, Kratos lumped mass: dt_crit = (h/c) / cos(pi / (4 n)) (1D fixed-free).

    One cell across and nu = 0: the y/z-symmetric x-modes are the 1D lumped
    bar. The shear (antisymmetric) modes are lower.
    """
    KM, SM = _kratos()
    S = _system(KM, SM, 1, 1.0, "kratos")
    dt, wmax = EX.estimate_dt_crit(S["K"], S["m"], S["free"])
    ref = (_H / _C) / np.cos(np.pi / (4 * _NX))
    print(f"CFL p=1: dt_crit={dt:.10f} analytic={ref:.10f}")
    assert abs(dt / ref - 1.0) < 1e-9, f"dt_crit {dt} vs analytic {ref}"
    assert abs(wmax * dt - 2.0) < 1e-12


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize("case", ["uncut_kratos_mass", "uncut_hrz", "cut_ext_1e-4", "cut_ext_1e-6"])
def test_dt_crit_execution_bracket(degree: int, case: str):
    """Kratos run at 0.95 x estimate stays within the leapfrog bound; 1.05x and 1.5x blow up."""
    KM, SM = _kratos()
    mass, eta = {"uncut_kratos_mass": ("kratos", 1.0), "uncut_hrz": ("hrz", 1.0),
                 "cut_ext_1e-4": ("hrz_ext", 1e-4), "cut_ext_1e-6": ("hrz_ext", 1e-6)}[case]
    S = _system(KM, SM, degree, eta, mass)
    dtc = _dt_crit(S)
    r_stable = _max_energy_ratio(KM, SM, S, 0.95 * dtc, _N_STABLE)
    r_105 = _max_energy_ratio(KM, SM, S, 1.05 * dtc, _N_UNSTABLE)
    r_150 = _max_energy_ratio(KM, SM, S, 1.5 * dtc, _N_UNSTABLE)
    print(f"CFL bracket p={degree} {case}: dt_crit={dtc:.6f} "
          f"maxE/E0 0.95x={r_stable:.3f} 1.05x={r_105:.2e} 1.5x={r_150:.2e}")
    assert r_stable <= _STABLE_BOUND, f"unstable at 0.95 dt_crit: {r_stable}"
    assert r_105 >= _BLOWUP, f"no blow-up at 1.05 dt_crit: {r_105}"
    assert r_150 >= _BLOWUP, f"no blow-up at 1.5 dt_crit: {r_150}"


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize("eta", [1e-4, 1e-6])
def test_cut_cell_dt_crit_at_least_half_uncut(degree: int, eta: float):
    """S24 claim: extended cut-cell dt_crit >= 0.5 x uncut (same HRZ lumping), measured.

    - Uncut limit: blow-up at 1.05 x dt_U, so the true limit is < 1.05 dt_U.
    - Cut + T: stable at 0.5 x 1.05 x dt_U, which is above half the uncut limit.
    - Control without T: dt_crit collapses as O(eta) (<= 10 eta x dt_U).
    """
    KM, SM = _kratos()
    U = _system(KM, SM, degree, 1.0, "hrz")
    X = _system(KM, SM, degree, eta, "hrz_ext")
    N = _system(KM, SM, degree, eta, "hrz")
    dt_u, dt_x, dt_n = _dt_crit(U), _dt_crit(X), _dt_crit(N)
    r_u = _max_energy_ratio(KM, SM, U, 1.05 * dt_u, _N_UNSTABLE)
    r_x = _max_energy_ratio(KM, SM, X, 0.5 * 1.05 * dt_u, _N_STABLE)
    print(f"S24 claim p={degree} eta={eta:.0e}: dt_U={dt_u:.6f} dt_ext/dt_U={dt_x / dt_u:.4f} "
          f"dt_noext/dt_U={dt_n / dt_u:.2e} uncut@1.05={r_u:.2e} cut+T@0.525dt_U={r_x:.3f}")
    assert dt_x >= 0.5 * dt_u, f"extended cut dt_crit {dt_x} < 0.5 x uncut {dt_u}"
    assert dt_n <= 10.0 * eta * dt_u, f"control: no-extension dt_crit did not collapse: {dt_n}"
    assert r_u >= _BLOWUP, f"uncut did not blow up at 1.05 dt_U: {r_u}"
    assert r_x <= _STABLE_BOUND, f"cut+T unstable at 0.525 dt_U: {r_x}"


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_mass_floor_restores_dt_on_unextended_cut_cell(degree: int):
    """floor_lumped_mass on an unextended cut cell (eta = 1e-2) restores dt_crit.

    m_min = 1 % of the smallest uncut HRZ mass. Checks:
    - dt_F >= 0.5 dt_U.
    - Added mass <= 1 % of the total.
    - The floored run is stable at 0.95 dt_F and blows up at 1.05 dt_F.
    - The unfloored masses blow up at 0.95 dt_F.
    """
    KM, SM = _kratos()
    eta = 1e-2
    U = _system(KM, SM, degree, 1.0, "hrz")
    N = _system(KM, SM, degree, eta, "hrz")
    dt_u, dt_n = _dt_crit(U), _dt_crit(N)
    m_min = 0.01 * float(np.min(U["m"][U["free"]]))
    m_f = EX.floor_lumped_mass(N["m"], m_min)
    dt_f = _dt_crit(N, m_f)
    added = float(np.sum(m_f - N["m"]) / np.sum(N["m"]))
    r_f095 = _max_energy_ratio(KM, SM, N, 0.95 * dt_f, _N_STABLE, m=m_f)
    r_f105 = _max_energy_ratio(KM, SM, N, 1.05 * dt_f, _N_UNSTABLE, m=m_f)
    r_raw = _max_energy_ratio(KM, SM, N, 0.95 * dt_f, _N_UNSTABLE)
    print(f"mass floor p={degree}: dt_noext/dt_U={dt_n / dt_u:.3e} dt_F/dt_U={dt_f / dt_u:.4f} "
          f"added={added:.2e} floored 0.95x={r_f095:.3f} 1.05x={r_f105:.2e} unfloored@0.95dt_F={r_raw:.2e}")
    assert dt_n <= 10.0 * eta * dt_u, f"control: unextended dt_crit did not collapse: {dt_n}"
    assert dt_f >= 0.5 * dt_u, f"floored dt_crit {dt_f} < 0.5 x uncut {dt_u}"
    assert added <= 0.01, f"mass floor added {added:.3%} of the total mass"
    assert r_f095 <= _STABLE_BOUND, f"floored run unstable at 0.95 dt_F: {r_f095}"
    assert r_f105 >= _BLOWUP, f"floored run did not blow up at 1.05 dt_F: {r_f105}"
    assert r_raw >= _BLOWUP, f"unfloored masses did not blow up at 0.95 dt_F: {r_raw}"


def test_kratos_freezes_sub_epsilon_nodal_masses():
    """Kratos sets a = 0 when NODAL_MASS <= eps (``numerical_limit`` in the scheme).

    Unextended p = 2, eta = 1e-4: the cut-cell HRZ masses are < 2.2e-16. The run
    at 0.5 dt_U is about 4000x dt_crit, yet it stays finite because those DOFs never
    move (u_end == u_0 exactly). A bounded Kratos run is therefore no evidence that
    dt_crit was restored; this is why the controls above use eta = 1e-2.
    """
    KM, SM = _kratos()
    U = _system(KM, SM, 2, 1.0, "hrz")
    N = _system(KM, SM, 2, 1e-4, "hrz")
    eps = float(np.finfo(float).eps)
    tiny = N["free"] & (N["m"] <= eps)
    assert np.any(tiny), "premise: some cut-cell HRZ masses are <= eps"
    dt_u, dt_n = _dt_crit(U), _dt_crit(N)
    r, u0, u_end = _max_energy_ratio(KM, SM, N, 0.5 * dt_u, _N_UNSTABLE, return_state=True)
    print(f"freeze: min m={N['m'][N['free']].min():.2e} dt_noext/dt_U={dt_n / dt_u:.2e} maxE/E0={r:.3f}")
    assert 0.5 * dt_u > 1000.0 * dt_n
    assert np.isfinite(r)
    assert np.array_equal(u_end[tiny], u0[tiny]), "sub-eps-mass DOFs moved"
    assert not np.array_equal(u_end[~tiny & N["free"]], u0[~tiny & N["free"]])
