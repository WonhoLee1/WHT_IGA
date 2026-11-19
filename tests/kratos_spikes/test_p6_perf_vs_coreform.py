# -*- coding: utf-8 -*-
"""Kratos integration tests for Phase 6 G-8 performance profiling & solver harness."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pytest

_WHT_IGA_ROOT = Path(__file__).resolve().parents[2]
if str(_WHT_IGA_ROOT) not in sys.path:
    sys.path.insert(0, str(_WHT_IGA_ROOT))

import scripts.perf_vs_coreform as perf

pytestmark = [pytest.mark.kratos]

RSS_LIMIT_MB = 1000.0  # Strictly < 1.0 GB while Track 1 runs in parallel


def test_coreform_baselines_csv_parsed() -> None:
    """Coreform glass-to-metal-seal CSV baselines parse to exact positive metrics."""
    cf = perf.load_coreform_baselines()
    assert abs(cf["cae_script_s"] - 4.98345822840929) < 1e-9
    assert abs(cf["iga_mesh_s"] - 318.8071711361408) < 1e-9
    assert abs(cf["interop_s"] - 3.6039445847272873) < 1e-9
    assert abs(cf["abaqus_solve_s"] - 691.0) < 1e-9
    assert cf["background_elements"] == 146_461
    assert cf["active_elements"] == 32_098
    assert cf["cut_elements"] == 3_420
    assert cf["abaqus_equations"] == 178_278
    assert cf["mesh_rate_active_elem_per_s"] > 0.0
    assert cf["solve_rate_eq_per_s"] > 0.0


def test_untrimmed_stage_profiling_and_rss() -> None:
    """Untrimmed 3D p=1 and p=2 stage breakdown produces positive rates and < 1 GB RSS."""
    t0 = time.perf_counter()
    r_p1 = perf.profile_untrimmed_stages((4, 4, 4), degree=1)
    r_p2 = perf.profile_untrimmed_stages((4, 4, 4), degree=2)
    dt = time.perf_counter() - t0
    assert dt < 30.0

    for row in (r_p1, r_p2):
        assert row["n_elements"] == 64
        assert row["element_rate_eps"] > 0.0
        for key in (
            "make_nurbs_volume_s",
            "attach_cell_elements_s",
            "kratos_setup_s",
            "kratos_residual_s",
            "kratos_build_s",
            "kratos_system_matrix_csr_s",
            "kratos_strategy_solve_s",
        ):
            assert row[key] > 0.0, (row["degree"], key)
        assert row["peak_rss_mb"] < RSS_LIMIT_MB

    # Ensure JSON serialization succeeds cleanly.
    parsed = json.loads(json.dumps({"p1": r_p1, "p2": r_p2}))
    assert parsed["p2"]["n_dofs"] == 3 * (4 + 2) ** 3


def test_trimmed_stage_profiling_and_extension() -> None:
    """Trimmed 3D p=2 stage breakdown exercises cut-cell quadrature, extension, and assembly."""
    t0 = time.perf_counter()
    row = perf.profile_trimmed_stages((4, 4, 4), degree=2, n_geo=6, workers=None)
    dt = time.perf_counter() - t0
    assert dt < 30.0

    assert row["n_cut_cells"] > 0
    assert row["n_out_cells"] > 0
    assert row["n_active_cells"] > 0
    assert row["n_extended"] > 0
    assert row["n_ties"] == 3 * row["n_extended"]
    assert row["element_rate_eps"] > 0.0
    for key in (
        "tessellate_s",
        "build_trimmed_table_s",
        "classify_basis_s",
        "build_extension_s",
        "make_nurbs_volume_s",
        "attach_trimmed_elements_s",
        "apply_extension_ties_s",
        "kratos_residual_s",
        "kratos_build_s",
        "kratos_system_matrix_csr_s",
        "kratos_strategy_solve_s",
    ):
        assert row[key] > 0.0, key
    assert row["peak_rss_mb"] < RSS_LIMIT_MB


def test_linear_solvers_agree_untrimmed_and_trimmed() -> None:
    """Direct (PardisoLU/LDLT, SkylineLU) and working iterative (AMGCL, CG) solvers agree < 1e-6."""
    t0 = time.perf_counter()
    untrimmed = perf.benchmark_linear_solvers((4, 4, 4), degree=2, trimmed=False)
    trimmed_ext = perf.benchmark_linear_solvers(
        (4, 4, 4),
        degree=2,
        trimmed=True,
        use_extension=True,
        n_geo=6,
        solvers=[
            "PardisoLU",
            "PardisoLDLT",
            "SkylineLU",
            "AMGCL_CG_ILU0",
            "AMGCL_GMRES_ILU0",
            "CG_Diag",
            "CG_AdditiveSchwarz",
            "Scaling_CG_Diag",
            "Scaling_AMGCL_CG",
        ],
    )
    dt = time.perf_counter() - t0
    assert dt < 30.0
    assert untrimmed["peak_rss_mb"] < RSS_LIMIT_MB
    assert trimmed_ext["peak_rss_mb"] < RSS_LIMIT_MB

    # Untrimmed: all tested Kratos solvers + isolated CSR solvers agree within < 1e-6.
    for row in untrimmed["solvers"]:
        assert row["error"] is None, row
        assert row["solution_accurate"] is True, row
        assert row["rel_l2_diff_vs_pardiso"] < 1e-6, row
    for row in untrimmed["csr_isolated"]:
        assert row["converged"] is True, row
        assert row["rel_residual"] < 1e-6, row
        assert row["rel_l2_diff_vs_pardiso_lu"] < 1e-6, row

    # Trimmed + extension: direct, AMGCL(ilu0/aggregation), CG_Diag, Scaling solvers agree < 1e-6.
    must_agree = {
        "PardisoLU",
        "PardisoLDLT",
        "SkylineLU",
        "AMGCL_CG_ILU0",
        "AMGCL_GMRES_ILU0",
        "CG_Diag",
        "Scaling_CG_Diag",
        "Scaling_AMGCL_CG",
    }
    by_name = {r["solver"]: r for r in trimmed_ext["solvers"]}
    for name in must_agree:
        row = by_name[name]
        assert row["error"] is None, row
        assert row["solution_accurate"] is True, row
        assert row["rel_l2_diff_vs_pardiso"] < 1e-6, row

    # AdditiveSchwarz executes cleanly (its trimmed divergence is captured quantitatively).
    assert "CG_AdditiveSchwarz" in by_name
    assert by_name["CG_AdditiveSchwarz"]["error"] is None


def test_fbar_sizing_extrapolation_and_outputs() -> None:
    """fbar E04/V07 sizing, 1M-DOF extrapolation, JSON, and Markdown outputs run < 30s and < 1 GB."""
    t0 = time.perf_counter()
    coreform = perf.load_coreform_baselines()
    u_p2 = [perf.profile_untrimmed_stages((4, 4, 4), degree=2), perf.profile_untrimmed_stages((5, 5, 5), degree=2)]
    t_p2 = [perf.profile_trimmed_stages((4, 4, 4), degree=2, n_geo=6, workers=None)]
    sol_u = perf.benchmark_linear_solvers(
        (4, 4, 4), degree=2, trimmed=False, solvers=["PardisoLU", "AMGCL_CG_ILU0"], include_csr_isolation=False
    )
    fbar = perf.profile_fbar_sizing(
        cook_sizes=(4,),
        v07_levels=(("h0=0.5", 6, 8, 2),),
        cook_n_steps=10,
    )
    extrap = perf.extrapolate_to_1m_and_coreform(u_p2, t_p2, sol_u, coreform)
    dt = time.perf_counter() - t0
    assert dt < 30.0
    assert fbar["peak_rss_mb"] < RSS_LIMIT_MB

    cook0 = fbar["e04_cook"][0]
    v070 = fbar["v07_kirsch"][0]
    for row in (cook0, v070):
        assert row["build_model_s"] > 0.0
        assert row["residual_tangent_sparse_s"] > 0.0
        assert row["single_solve_pardiso_s"] > 0.0
        assert row["single_solve_superlu_s"] > 0.0
        assert row["single_solve_rel_diff"] < 1e-8
        assert row["solve_static_rel_diff"] < 1e-8

    report = {
        "coreform_baselines": coreform,
        "untrimmed_p1": [],
        "untrimmed_p2": u_p2,
        "trimmed_p2": t_p2,
        "solver_comparison_untrimmed_p2": sol_u,
        "fbar_sizing": fbar,
        "extrapolation": extrap,
    }
    md = perf.format_markdown_summary(report)
    assert "| Solver | Strategy (s) |" in md
    assert "Phase6_Gate_1M_DOF" in json.dumps(report)
    assert len(extrap["hotspot_ranking_1m"]) >= 5
