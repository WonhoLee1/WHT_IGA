# -*- coding: utf-8 -*-
"""S29b: unilateral Kratos penalty mortar contact on IGA blocks (``analyses.kratos_contact``).

Model: two p = 1 IGA blocks, unit cross-section, z in [0, 1] and [1 + g0, 2 + g0],
E = 1000, nu = 0, tied quad surface layers at the interface (lower = slave).
Bottom fixed in z, x = 0 / y = 0 faces fixed in x / y, top prescribed u_z.

References (independent of the run):
- Monolithic 1D bar: two blocks in series, k_s = E*A/(2L) = 500, force E*delta/2.
- Penalty spring of Kratos' mortar formulation (per condition e of area A_e):
  pressure_i = penalty * weighted gap_i = penalty * g * A_e/4 (integral of N_i over e),
  nodal force_i = pressure_i * A_e/4  ->  k_e = penalty * A_e**2/4, k_c = sum_e k_e
  (n x n quads on the unit face: k_c = penalty / (4 n**2)).
  Bar in series with k_c: R = d / (1/k_s + 1/k_c). Error vs E*d/2 is O(1/penalty).
- Open contact: no force, and the upper block moves rigidly (gap = g0 + u_top).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses.kratos_contact import (
    APPS,
    SurfaceLayerContact,
    add_contact_variables,
    add_tied_quad_layer,
)
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_E = 1000.0
_L = 1.0
_A = 1.0
_K_BAR = _E * _A / (2.0 * _L)  # two unit blocks in series


def _spring_force(d: float, penalty: float, n: int = 1) -> float:
    """Bar (k_s) in series with the penalty spring of n x n quad pairs (k_c = penalty*A**2/(4 n**2)).

    ``d`` is the top travel beyond first touch (> 0: closed). Valid for a uniform gap,
    which holds here (nu = 0, identical conditions).
    """
    if d <= 0.0:
        return 0.0
    k_c = penalty * _A**2 / (4.0 * n * n)
    return d / (1.0 / _K_BAR + 1.0 / k_c)


class _TwoBlocks:
    """Two stacked p=1 blocks with tied quad layers and a SurfaceLayerContact."""

    def __init__(self, penalty: float, g0: float = 0.0, n: int = 1):
        mods = load_kratos(APPS)
        KM = mods["KratosMultiphysics"]
        SM = mods["StructuralMechanicsApplication"]
        self.KM = KM
        self.model = KM.Model()
        mp = self.model.CreateModelPart("Structure")
        add_contact_variables(mp)
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        mp.SetBufferSize(2)
        props = mp.GetProperties()[1]
        props.SetValue(KM.YOUNG_MODULUS, _E)
        props.SetValue(KM.POISSON_RATIO, 0.0)
        props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
        layer_props = mp.GetProperties()[2]
        self.top_z = 2.0 + g0
        lo = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, 1), degree=1)
        hi = BoxSpline(lo=(0.0, 0.0, 1.0 + g0), hi=(1.0, 1.0, self.top_z), n_cells=(n, n, 1), degree=1)
        vol_lo, _ = make_nurbs_volume(KM, mp, lo, first_node_id=1)
        attach_cell_elements(KM, mp, vol_lo, lo, props, start_id=1)
        vol_hi, _ = make_nurbs_volume(KM, mp, hi, first_node_id=1001)
        attach_cell_elements(KM, mp, vol_hi, hi, props, start_id=1001)
        self.lo_ids = set(range(1, 1 + lo.n_ctrl))
        self.hi_ids = set(range(1001, 1001 + hi.n_ctrl))
        nu, nv, nw = lo.n_cp

        def face(first: int, c: int):
            return [[mp.GetNode(first + a + nu * (b + nv * c)) for b in range(nv)] for a in range(nu)]

        _, kid = add_tied_quad_layer(mp, "LayerLo", face(1, nw - 1), (0.0, 0.0, 1.0), layer_props,
                                     first_node_id=5001, first_condition_id=5001, first_constraint_id=1)
        add_tied_quad_layer(mp, "LayerHi", face(1001, 0), (0.0, 0.0, -1.0), layer_props,
                            first_node_id=6001, first_condition_id=6001, first_constraint_id=kid)
        for node in mp.Nodes:
            if node.Id not in self.lo_ids | self.hi_ids:
                continue
            if abs(node.X0) < 1e-12:
                node.Fix(KM.DISPLACEMENT_X)
            if abs(node.Y0) < 1e-12:
                node.Fix(KM.DISPLACEMENT_Y)
            if abs(node.Z0) < 1e-12 or abs(node.Z0 - self.top_z) < 1e-12:
                node.Fix(KM.DISPLACEMENT_Z)
        self.mp = mp
        self.contact = SurfaceLayerContact(self.model, mp, "LayerLo", "LayerHi", penalty, contact_property_id=3)
        self.time = 0.0

    def step(self, top_uz: float) -> dict:
        KM, mp = self.KM, self.mp
        top = [n for n in mp.Nodes if n.Id in self.hi_ids and abs(n.Z0 - self.top_z) < 1e-12]

        def load():
            for node in top:
                node.SetSolutionStepValue(KM.DISPLACEMENT_Z, float(top_uz))

        self.time += 1.0
        converged = self.contact.solve_step(self.time, load)
        uz = lambda name: [n.Z0 + n.GetSolutionStepValue(KM.DISPLACEMENT_Z)
                           for n in mp.GetSubModelPart(name).Nodes]
        return {
            "converged": converged,
            "R": sum(n.GetSolutionStepValue(KM.REACTION_Z) for n in mp.Nodes
                     if n.Id in self.lo_ids and abs(n.Z0) < 1e-12),
            "gap": float(np.mean(uz("LayerHi")) - np.mean(uz("LayerLo"))),
            "gap_spread": float(np.ptp(uz("LayerHi"))) + float(np.ptp(uz("LayerLo"))),
            "state": self.contact.slave_state(),
            "iterations": self.contact.iterations,
        }


def test_s29b_compression_matches_penalty_spring_and_tends_to_monolithic():
    """R(penalty) = penalty-spring model to 1e-9; |R - E d/2| / (E d/2) = 2000/penalty / (1 + ...)
    drops 10x per decade (measured 2.0e-3, 2.0e-4, 2.0e-5)."""
    d = 0.01
    ref = _E * d / 2.0
    errs = []
    for pen in (1e6, 1e7, 1e8):
        r = _TwoBlocks(pen).step(-d)
        assert r["converged"] and r["iterations"] <= 3, r
        assert all(active for active, _ in r["state"].values())
        assert r["R"] == pytest.approx(_spring_force(d, pen), rel=1e-9)
        errs.append(abs(r["R"] - ref) / ref)
    assert errs[-1] < 3e-5, errs
    for coarse, fine in zip(errs, errs[1:]):
        assert 9.0 < coarse / fine < 11.0, errs


def test_s29b_pull_apart_separates_with_zero_force():
    """Pulling the top up: no tension, rigid lift of the upper block (gap = +d),
    all slave nodes inactive, weighted gap = d * A_i with A_i = 1/4."""
    d = 0.01
    r = _TwoBlocks(1e6).step(+d)
    assert r["converged"], r
    assert abs(r["R"]) < 1e-8 * _E * d, r["R"]
    assert r["gap"] == pytest.approx(d, abs=1e-12)
    assert r["gap_spread"] < 1e-12
    for active, wg in r["state"].values():
        assert not active
        assert wg == pytest.approx(d * 0.25, rel=1e-9)


@pytest.mark.parametrize("top_uz", [-0.003, -0.0049, -0.0051, -0.01])
def test_s29b_initial_gap_closes_only_when_top_travel_exceeds_it(top_uz):
    """Upper block lifted by g0 = 0.005. |u| < g0: zero force, gap = g0 + u.
    |u| > g0: force of the closed-gap model with closure |u| - g0."""
    g0, pen = 0.005, 1e6
    r = _TwoBlocks(pen, g0=g0).step(top_uz)
    assert r["converged"], r
    closure = -top_uz - g0
    if closure < 0.0:
        assert abs(r["R"]) < 1e-8 * _E * abs(top_uz), r["R"]
        assert r["gap"] == pytest.approx(g0 + top_uz, abs=1e-12)
    else:
        assert r["R"] == pytest.approx(_spring_force(closure, pen), rel=1e-9)
        assert r["gap"] < 0.0  # penalty penetration


def test_s29b_load_path_press_release_pull_repress():
    """One strategy over a history: each step matches its own reference (no memory of
    the previous active set): press 0.005, 0.01, release to 0, pull 0.01, press 0.002, 0.01."""
    pen = 1e6
    blocks = _TwoBlocks(pen)
    for top_uz in (-0.005, -0.01, 0.0, 0.01, -0.002, -0.01):
        r = blocks.step(top_uz)
        assert r["converged"], (top_uz, r)
        ref = _spring_force(-top_uz, pen)
        if ref == 0.0:
            assert abs(r["R"]) < 1e-8 * _E * 0.01, (top_uz, r["R"])
            assert r["gap"] == pytest.approx(top_uz, abs=1e-12)
        else:
            assert r["R"] == pytest.approx(ref, rel=1e-9), top_uz


def test_s29b_two_by_two_interface_compression_and_separation():
    """2x2 cells per block (4 quads per layer, mortar pairs across shared edges).
    Compression = penalty-spring model with k_c = penalty/16 (to 1e-9) and tends to
    E d/2 as O(1/penalty) (measured 7.9e-3, 8.0e-4, 8.0e-5); pulling: zero force, gap = d."""
    d = 0.01
    ref = _E * d / 2.0
    errs = []
    for pen in (1e6, 1e7, 1e8):
        r = _TwoBlocks(pen, n=2).step(-d)
        assert r["converged"], r
        assert r["R"] == pytest.approx(_spring_force(d, pen, n=2), rel=1e-9)
        errs.append(abs(r["R"] - ref) / ref)
    assert errs[-1] < 1e-4, errs
    for coarse, fine in zip(errs, errs[1:]):
        assert 9.0 < coarse / fine < 11.0, errs
    r = _TwoBlocks(1e6, n=2).step(+d)
    assert r["converged"], r
    assert abs(r["R"]) < 1e-8 * _E * d, r["R"]
    assert r["gap"] == pytest.approx(d, abs=1e-12)
