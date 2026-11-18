# -*- coding: utf-8 -*-
"""Cross-method comparison matrix: accuracy vs cost (Track H, Phase 3).

Compares full IGA / THB / hybrid-mortar / native FE on identical problems.
Metrics per cell: tip error vs analytic, DOFs, Pardiso peak estimate,
setup/solve/total wall time. Printed as a table for S35; gates below pin the
characterized values (all real Kratos solves, no mocks).

Exclusions (evidence-based, not hidden):
- p=1 in bending: linear elements are shear-stiff (thb-p1 measured 47.6% err).
  Bending parity requires p>=2 (matches recommend_discretization).
- barycentric bending: table-only (no characterized bound yet); uniaxial
  asserts its 2e-3 characterized gate.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pytest

from wht_iga.analyses.methods import (
    BoxMethodProblem,
    FullIgaOption,
    IgaFeHybridOption,
    NativeFeOption,
    ThbAdaptiveIgaOption,
    compare_method_solutions,
    solve_with_discretization_method,
)
from wht_iga.analyses.methods_recommend import _estimate_peak_bytes

pytestmark = pytest.mark.kratos


def _uniaxial_problem() -> tuple[BoxMethodProblem, float]:
    E_mod, nu, L, t_x = 2.0e5, 0.0, 10.0, 2.0
    return (
        BoxMethodProblem(
            lo=(0.0, 0.0, 0.0), hi=(L, 1.0, 1.0), E=E_mod, nu=nu,
            fixed_face=(0, 0), load_face=(0, 1), traction=(t_x, 0.0, 0.0),
        ),
        t_x * L / E_mod,
    )


def _bending_problem() -> tuple[BoxMethodProblem, float]:
    E_mod, nu = 2.0e11, 0.0
    return (
        BoxMethodProblem(
            lo=(-0.05, -0.05, 0.0), hi=(0.05, 0.05, 1.0), E=E_mod, nu=nu,
            fixed_face=(2, 0), load_face=(2, 1), traction=(0.0, -500.0, 0.0),
        ),
        -1.0e-6,
    )


def _collect(name, sol, tip_ref, tip_comp, n_cells, degree):
    tip_v = sol.tip_displacement[int(tip_comp)]
    return {
        "method": name,
        "tip_rel": abs(float(tip_v) - tip_ref) / abs(tip_ref),
        "dofs": int(sol.total_dofs),
        "peak_B": _estimate_peak_bytes(tuple(int(v) for v in n_cells), int(degree)),
        "setup_s": float(sol.timings.get("setup_s", -1.0)),
        "solve_s": float(sol.timings.get("solve_s", -1.0)),
        "total_s": float(sol.timings.get("total_s", -1.0)),
    }


def _print_table(title, rows):
    print(f"[compare] {title}", flush=True)
    print("[compare] method | tip_rel | dofs | peak_B | setup_s | solve_s | total_s", flush=True)
    for r in rows:
        print(
            f"[compare] {r['method']} | {r['tip_rel']:.3e} | {r['dofs']} | "
            f"{r['peak_B']:.3g} | {r['setup_s']:.2f} | {r['solve_s']:.2f} | {r['total_s']:.2f}",
            flush=True,
        )


def test_compare_uniaxial_matrix():
    """Uniaxial patch: exactness + cost across methods (coarse + fine)."""
    problem, tip_ref = _uniaxial_problem()
    pts = np.array([[2.5, 0.5, 0.5], [5.0, 0.5, 0.5], [7.5, 0.5, 0.5], [10.0, 0.5, 0.5]])
    cases = [
        ("full-p1", FullIgaOption(kind="full_iga_body_fit", degree=1, n_cells=(4, 2, 2)), (4, 2, 2), 1, 1e-5),
        ("full-p2", FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(4, 2, 2)), (4, 2, 2), 2, 1e-5),
        ("full-p3", FullIgaOption(kind="full_iga_body_fit", degree=3, n_cells=(4, 2, 2)), (4, 2, 2), 3, 1e-5),
        ("full-p2-fine", FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(8, 4, 4)), (8, 4, 4), 2, 1e-5),
        ("thb-p1", ThbAdaptiveIgaOption(degree=1, coarse_n_cells=(2, 1, 1), refinement_level=1,
                                        refine_box=((0.0, 0.0, 0.0), (4.0, 1.0, 1.0))), (2, 1, 1), 1, 1e-9),
        ("thb-p2", ThbAdaptiveIgaOption(degree=2, coarse_n_cells=(2, 1, 1), refinement_level=1,
                                        refine_box=((0.0, 0.0, 0.0), (4.0, 1.0, 1.0))), (2, 1, 1), 2, 1e-9),
        ("thb-p3", ThbAdaptiveIgaOption(degree=3, coarse_n_cells=(2, 1, 1), refinement_level=1,
                                        refine_box=((0.0, 0.0, 0.0), (4.0, 1.0, 1.0))), (2, 1, 1), 3, 1e-9),
        ("hyb-mortar", IgaFeHybridOption(fe_element_type="hex20", fe_n_cells=(2, 2, 2),
                                         iga_degree=2, iga_n_cells=(2, 2, 2), split_axis=0,
                                         split_fraction=0.5, fe_side="lo",
                                         coupling_kind="mortar_l2"), (2, 2, 2), 2, 1e-6),
        ("hyb-bary", IgaFeHybridOption(fe_element_type="hex20", fe_n_cells=(2, 2, 2),
                                       iga_degree=2, iga_n_cells=(2, 2, 2), split_axis=0,
                                       split_fraction=0.5, fe_side="lo",
                                       coupling_kind="barycentric"), (2, 2, 2), 2, 2e-3),
        ("fe-hex8", NativeFeOption(element_type="hex8", n_cells=(4, 2, 2)), (4, 2, 2), 1, 1e-8),
        ("fe-hex20", NativeFeOption(element_type="hex20", n_cells=(4, 2, 2)), (4, 2, 2), 2, 1e-8),
        ("fe-hex20-fine", NativeFeOption(element_type="hex20", n_cells=(8, 4, 4)), (8, 4, 4), 2, 1e-8),
    ]
    rows = []
    ref_sol = None
    for name, opt, nc, deg, tol in cases:
        sol = solve_with_discretization_method(problem, opt)
        assert sol.converged is True
        if name == "full-p2":
            ref_sol = sol
        r = _collect(name, sol, tip_ref, 0, nc, deg)
        assert r["tip_rel"] < tol, (name, r["tip_rel"])
        if ref_sol is not None and name != "full-p2":
            cmp = compare_method_solutions(ref_sol, sol, pts)
            assert cmp["l2_rel"] < max(1e-4, tol), (name, cmp["l2_rel"])
        rows.append(r)
    _print_table("uniaxial", rows)


def test_compare_bending_matrix():
    """Cantilever bending: 3% parity + cost (p>=2 only, see module docstring)."""
    problem, tip_ref = _bending_problem()
    cases = [
        ("full-p2", FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(2, 2, 8)), (2, 2, 8), 2),
        ("full-p3", FullIgaOption(kind="full_iga_body_fit", degree=3, n_cells=(2, 2, 8)), (2, 2, 8), 3),
        ("thb-p2", ThbAdaptiveIgaOption(degree=2, coarse_n_cells=(1, 1, 4), refinement_level=1,
                                        refine_box=((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4))), (1, 1, 4), 2),
        ("thb-p3", ThbAdaptiveIgaOption(degree=3, coarse_n_cells=(1, 1, 4), refinement_level=1,
                                        refine_box=((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4))), (1, 1, 4), 3),
        ("hyb-mortar", IgaFeHybridOption(fe_element_type="hex20", fe_n_cells=(2, 2, 4),
                                         iga_degree=2, iga_n_cells=(2, 2, 4), split_axis=2,
                                         split_fraction=0.5, fe_side="lo",
                                         coupling_kind="mortar_l2"), (2, 2, 4), 2),
        ("fe-hex20", NativeFeOption(element_type="hex20", n_cells=(2, 2, 8)), (2, 2, 8), 2),
    ]
    rows = []
    ref_sol = None
    for name, opt, nc, deg in cases:
        sol = solve_with_discretization_method(problem, opt)
        assert sol.converged is True
        if name == "full-p2":
            ref_sol = sol
        r = _collect(name, sol, tip_ref, 1, nc, deg)
        assert r["tip_rel"] < 0.03, (name, r["tip_rel"])
        if ref_sol is not None and name != "full-p2":
            cmp = compare_method_solutions(ref_sol, sol, problem.tip_point().reshape(1, 3))
            assert cmp["tip_rel"] < 0.03, (name, cmp["tip_rel"])
        rows.append(r)
    # barycentric bending: table-only (no characterized bound yet)
    opt_b = IgaFeHybridOption(fe_element_type="hex20", fe_n_cells=(2, 2, 4),
                              iga_degree=2, iga_n_cells=(2, 2, 4), split_axis=2,
                              split_fraction=0.5, fe_side="lo",
                              coupling_kind="barycentric")
    sol_b = solve_with_discretization_method(problem, opt_b)
    assert sol_b.converged is True
    rows.append(_collect("hyb-bary*", sol_b, tip_ref, 1, (2, 2, 4), 2))
    _print_table("bending (* = table-only, uncharacterized)", rows)
