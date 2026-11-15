"""Unit tests for Phase 3 post-processing (wht_iga.post) and query (wht_iga.query)."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from wht_iga.grid.greville import greville
from wht_iga.grid.knots import uniform_clamped
from wht_iga.post import (
    VTK_HEXAHEDRON,
    VTK_TETRA,
    BSplineFieldEvaluator,
    ProxyMesh,
    build_proxy_mesh,
    compute_energy_norm,
    compute_h1_norm,
    compute_l2_norm,
    decompose_hex_to_6_tets,
    generate_cut_tet_proxy,
    generate_hex_proxy_for_cells,
    integrate_error_against_analytic,
    integrate_error_on_tets,
    read_vtkhdf,
    write_vtkhdf,
)
from wht_iga.query import ResultView, von_mises_stress
from wht_iga.solver.db import write_results_db


def _tet_signed_volumes(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    v = points[tets]
    return np.einsum("ij,ij->i", v[:, 1] - v[:, 0], np.cross(v[:, 2] - v[:, 0], v[:, 3] - v[:, 0])) / 6.0


def test_hex_and_cut_tet_proxy_generation() -> None:
    in_bounds = np.array([[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]], dtype=np.float64)
    hex_mesh = generate_hex_proxy_for_cells(in_bounds, n_sub=2)
    assert hex_mesh.n_cells == 8
    assert hex_mesh.n_points == 27
    assert np.all(hex_mesh.cell_types == VTK_HEXAHEDRON)

    # 6-tet Kuhn decomposition preserves unit volume with all positive Jacobians
    tets = decompose_hex_to_6_tets(hex_mesh.connectivity.reshape((-1, 8)))
    vols = _tet_signed_volumes(hex_mesh.points, tets)
    assert np.all(vols > 0.0)
    assert np.isclose(np.sum(vols), 1.0, atol=1e-13)

    # Cut cells clipped by plane x <= 0.5 -> exact half volume
    cut_bounds = np.array([[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]], dtype=np.float64)
    cut_mesh = generate_cut_tet_proxy(
        cell_bounds=cut_bounds,
        n_sub=2,
        phi_fn=lambda pts: pts[:, 0] - 0.5,
    )
    assert cut_mesh.n_cells > 0
    assert np.all(cut_mesh.cell_types == VTK_TETRA)
    cut_tets = cut_mesh.connectivity.reshape((-1, 4))
    cut_vol = float(np.sum(np.abs(_tet_signed_volumes(cut_mesh.points, cut_tets))))
    assert np.isclose(cut_vol, 0.5, atol=1e-12)

    # Pre-supplied cut tets path
    pts_pre = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    tets_pre = np.array([[0, 1, 2, 3]])
    pre_mesh = generate_cut_tet_proxy(cut_points=pts_pre, cut_tets=tets_pre)
    assert pre_mesh.n_cells == 1
    assert pre_mesh.n_points == 4


def test_bspline_proxy_evaluation_matches_cp_reevaluation() -> None:
    ku = uniform_clamped(0.0, 2.0, n_cells=2, degree=2)
    kv = uniform_clamped(0.0, 1.0, n_cells=2, degree=2)
    kw = uniform_clamped(0.0, 1.0, n_cells=1, degree=2)
    evaluator = BSplineFieldEvaluator((ku, kv, kw), (2, 2, 2))

    gx = greville(ku, 2)
    gy = greville(kv, 2)
    gz = greville(kw, 2)
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing="ij")

    # Quadratic displacement u = [x^2 + 2y, y^2 - z, 3xz] is exactly reproduced by p=2 B-splines
    # Verify identity map x(xi) = xi via Greville abscissae first:
    cp_xyz = np.stack([GX, GY, GZ], axis=-1)
    test_pts = np.array([[0.3, 0.4, 0.5], [1.7, 0.8, 0.2], [1.0, 0.5, 0.5]])
    assert np.allclose(evaluator.evaluate(test_pts, cp_xyz), test_pts, atol=1e-13)

    # Build linear + quadratic displacement field on control points
    cp_u = np.stack([2.0 * GX - 0.5 * GY, 3.0 * GY + GZ, -GX + 4.0 * GZ], axis=-1)
    in_bounds = np.array([[[0.0, 0.0, 0.0], [1.0, 0.5, 1.0]]], dtype=np.float64)
    cut_bounds = np.array([[[1.0, 0.0, 0.0], [2.0, 0.5, 1.0]]], dtype=np.float64)

    proxy = build_proxy_mesh(
        in_bounds=in_bounds,
        cut_bounds=cut_bounds,
        n_sub_in=2,
        n_sub_cut=2,
        inside_fn=lambda pts: pts[:, 0] <= 1.75,
        evaluator=evaluator,
        cp_fields={"U": cp_u},
    )
    assert "U" in proxy.point_data
    assert "BlockId" in proxy.cell_data

    # Re-evaluate directly from control points at all proxy vertices -> exact 1e-12 match
    u_re = evaluator.evaluate(proxy.points, cp_u)
    assert np.allclose(proxy.point_data["U"], u_re, atol=1e-12)

    # Gradient check for linear field
    grad_u = evaluator.evaluate_gradient(test_pts, cp_u)
    expected_grad = np.array([[2.0, -0.5, 0.0], [0.0, 3.0, 1.0], [-1.0, 0.0, 4.0]])
    for k in range(len(test_pts)):
        assert np.allclose(grad_u[k], expected_grad, atol=1e-12)


def test_vtkhdf_write_schema_and_control_points_roundtrip(tmp_path) -> None:
    in_bounds = np.array([[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]], dtype=np.float64)
    cut_bounds = np.array([[[1.0, 0.0, 0.0], [2.0, 1.0, 1.0]]], dtype=np.float64)

    frames: list[ProxyMesh] = []
    cp_dict: dict[int, dict[str, np.ndarray]] = {}
    for step_idx, scale in enumerate([1.0, 2.5]):
        mesh = build_proxy_mesh(in_bounds=in_bounds, cut_bounds=cut_bounds, n_sub_in=1, n_sub_cut=1)
        mesh.point_data["U"] = mesh.points * scale
        mesh.point_data["NT11"] = 20.0 + scale * mesh.points[:, 0]
        mesh.cell_data["CellTag"] = np.arange(mesh.n_cells, dtype=np.float64)
        frames.append(mesh)
        cp_dict[step_idx] = {
            "U": np.ones((8, 3), dtype=np.float64) * scale,
            "NT11": np.full(8, 20.0 * scale, dtype=np.float64),
            "state": np.array([1, 1, 2, 0], dtype=np.int64),
        }

    out_file = tmp_path / "job.vtkhdf"
    write_vtkhdf(
        out_file,
        frames,
        control_points=cp_dict,
        times=[0.0, 1.0],
        block_names=["in_cells", "cut_cells"],
    )

    # Direct h5py schema inspection
    with h5py.File(out_file, "r") as h5:
        assert "VTKHDF" in h5
        root = h5["VTKHDF"]
        assert tuple(root.attrs["Version"]) == (2, 1)
        assert root.attrs["Type"] == b"UnstructuredGrid"
        assert "Steps" in root
        assert int(root["Steps"].attrs["NSteps"]) == 2
        assert np.allclose(root["Steps"]["Values"][:], [0.0, 1.0])
        assert "U" in root["PointData"]
        assert "NT11" in root["PointData"]
        assert "CellTag" in root["CellData"]
        assert "WHT_IGA" in h5
        assert "control_points" in h5["WHT_IGA"]
        assert set(h5["WHT_IGA"]["control_points"]["0"].keys()) == {"U", "NT11", "state"}

    # Roundtrip via read_vtkhdf
    loaded = read_vtkhdf(out_file)
    assert len(loaded["frames"]) == 2
    assert loaded["block_names"] == ["in_cells", "cut_cells"]
    for f_idx, f_orig in enumerate(frames):
        f_read = loaded["frames"][f_idx]
        assert np.allclose(f_read.points, f_orig.points, atol=1e-12)
        assert np.array_equal(f_read.connectivity, f_orig.connectivity)
        assert np.array_equal(f_read.cell_types, f_orig.cell_types)
        assert np.allclose(f_read.point_data["U"], f_orig.point_data["U"], atol=1e-12)
        assert np.allclose(f_read.point_data["NT11"], f_orig.point_data["NT11"], atol=1e-12)
        assert np.allclose(
            loaded["control_points"][f_idx]["U"], cp_dict[f_idx]["U"], atol=1e-12
        )
        assert np.array_equal(
            loaded["control_points"][f_idx]["state"], cp_dict[f_idx]["state"]
        )


def test_norms_exact_polynomial_and_energy() -> None:
    bounds = np.array([[[0.0, 0.0, 0.0], [2.0, 1.0, 1.0]]], dtype=np.float64)

    # Exact polynomial u = x^2 + 2 y^2 + 3 z^2 -> n_gp=2 integrates degree 3, n_gp=3 integrates degree 5
    def u_exact(pts: np.ndarray) -> np.ndarray:
        return pts[:, 0] ** 2 + 2.0 * pts[:, 1] ** 2 + 3.0 * pts[:, 2] ** 2

    def grad_exact(pts: np.ndarray) -> np.ndarray:
        return np.column_stack([2.0 * pts[:, 0], 4.0 * pts[:, 1], 6.0 * pts[:, 2]])

    # Exact match -> zero error
    h1_zero = integrate_error_against_analytic(
        bounds, u_exact, u_exact, n_gp=3, eval_grad_fn=grad_exact, analytic_grad_fn=grad_exact
    )
    assert h1_zero.l2.abs_error < 1e-13
    assert h1_zero.h1_semi.abs_error < 1e-13
    assert h1_zero.h1.abs_error < 1e-13

    # Known constant offset e = 1.0 -> ||e||_L2 = sqrt(Vol) = sqrt(2.0)
    l2_const = integrate_error_against_analytic(
        bounds,
        eval_num_fn=lambda pts: u_exact(pts) + 1.0,
        analytic_fn=u_exact,
        n_gp=3,
    )
    assert np.isclose(l2_const.abs_error, np.sqrt(2.0), atol=1e-12)

    # Energy norm with Voigt 6-vectors: E=100, nu=0 -> sigma = 100 * eps
    w = np.array([0.5, 0.5], dtype=np.float64)
    j = np.array([2.0, 2.0], dtype=np.float64)  # total integral = 2.0
    eps_ref = np.zeros((2, 6), dtype=np.float64)
    eps_ref[:, 0] = 1e-3
    sig_ref = 100.0 * eps_ref
    eps_num = 1.1 * eps_ref
    sig_num = 1.1 * sig_ref
    en = compute_energy_norm(eps_num, sig_num, eps_ref, sig_ref, w, j)
    # (0.1 * 1e-3) * (0.1 * 0.1) * 2.0 = 2e-6 -> sqrt = sqrt(2e-6)
    expected_abs = np.sqrt(2.0 * (0.1 * 1e-3) * (0.1 * 0.1))
    assert np.isclose(en.abs_error, expected_abs, rtol=1e-12)
    assert np.isclose(en.rel_error, 0.1, rtol=1e-12)

    # Standalone compute_l2_norm & compute_h1_norm
    u_n = np.array([1.0, 2.0])
    u_r = np.array([1.0, 2.0])
    g_n = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    g_r = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    h1_res = compute_h1_norm(u_n, g_n, u_r, g_r, w, j)
    assert h1_res.h1.abs_error == 0.0
    l2_res = compute_l2_norm(u_n + 2.0, u_r, w, j)
    assert np.isclose(l2_res.abs_error, np.sqrt(2.0 * 4.0), atol=1e-12)


def test_result_view_queries_and_sqlite_integration(tmp_path) -> None:
    # Build a 2x1x1 hex mesh on [0, 2] x [0, 1] x [0, 1]
    bounds = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
            [[1.0, 0.0, 0.0], [2.0, 1.0, 1.0]],
        ],
        dtype=np.float64,
    )
    mesh = generate_hex_proxy_for_cells(bounds, n_sub=2)
    pts = mesh.points
    # Linear displacement U = [0.01*x, -0.005*y, 0.002*z]
    u_field = np.column_stack([0.01 * pts[:, 0], -0.005 * pts[:, 1], 0.002 * pts[:, 2]])
    # Bending + membrane stress: sxx(z) = 100 + 60 * (2*z - 1) -> membrane=100, bending2=+60
    s_field = np.zeros((mesh.n_points, 6), dtype=np.float64)
    s_field[:, 0] = 100.0 + 60.0 * (2.0 * pts[:, 2] - 1.0)
    # Reaction forces balancing an applied load of [100, -50, 25]
    rf_field = np.zeros((mesh.n_points, 3), dtype=np.float64)
    rf_field[0] = np.array([-100.0, 50.0, -25.0])

    mesh.point_data["DISPLACEMENT"] = u_field
    mesh.point_data["CAUCHY_STRESS_TENSOR"] = s_field
    mesh.point_data["REACTION"] = rf_field

    vtkhdf_path = tmp_path / "beam.vtkhdf"
    write_vtkhdf(
        vtkhdf_path,
        mesh,
        control_points={0: {"U": u_field[:8], "NT11": np.full(8, 25.0)}},
        times=[1.0],
        block_names=["in_cells"],
    )

    sqlite_path = tmp_path / "beam_iga_results.sqlite"
    write_results_db(
        str(sqlite_path),
        spec_json={"job": "beam"},
        meta_extra={"solver": "test"},
        frames=[{"frame_id": 0, "step_id": 0, "time": 1.0, "increment": 1, "converged": 1, "n_iter": 1}],
        cp_fields={(0, "U"): u_field[:8]},
    )

    rv = ResultView(vtkhdf_path, sqlite_path)
    assert rv.meta()["solver"] == "test"
    assert len(rv.frames()) == 1

    # Alias mapping: U <-> DISPLACEMENT, S <-> CAUCHY_STRESS_TENSOR, RF <-> REACTION
    assert np.allclose(rv.field("U"), u_field, atol=1e-12)
    assert np.allclose(rv.field("S"), s_field, atol=1e-12)
    assert np.allclose(rv.field("RF"), rf_field, atol=1e-12)
    assert np.allclose(rv.control_point_field("TEMPERATURE"), 25.0, atol=1e-12)

    # Probe arbitrary interior points -> exact linear recovery
    probe_pts = np.array([[0.37, 0.61, 0.83], [1.45, 0.25, 0.15]])
    probed_u = rv.probe(probe_pts, field_name="U")
    expected_u = np.column_stack([0.01 * probe_pts[:, 0], -0.005 * probe_pts[:, 1], 0.002 * probe_pts[:, 2]])
    assert np.allclose(probed_u, expected_u, atol=1e-12)

    # Hotspots (top von Mises at z=1 where sxx = 160 MPa)
    spots = rv.hotspots(field_name="S", top_n=5)
    assert len(spots) == 5
    assert np.isclose(spots[0]["value"], 160.0, atol=1e-12)
    assert np.isclose(spots[0]["coords"][2], 1.0, atol=1e-12)
    assert np.isclose(von_mises_stress(spots[0]["raw"]), 160.0, atol=1e-12)

    # Equilibrium check
    eq_ok = rv.check_equilibrium(applied_force=np.array([100.0, -50.0, 25.0]), tol=1e-10)
    assert eq_ok.balanced is True
    assert eq_ok.rel_error < 1e-12

    eq_bad = rv.check_equilibrium(applied_force=np.array([101.0, -50.0, 25.0]), tol=1e-10)
    assert eq_bad.balanced is False

    # Section force across x = 0.75 (area = 1.0, integral of sxx(z) over z in [0, 1] is 100.0)
    f_sec = rv.section_force(plane_point=np.array([0.75, 0.5, 0.5]), plane_normal=np.array([1.0, 0.0, 0.0]))
    assert np.allclose(f_sec, [100.0, 0.0, 0.0], atol=1e-10)

    # Stress linearization along z in [0, 1] at (x=0.5, y=0.5)
    scl = rv.stress_linearize(p1=np.array([0.5, 0.5, 0.0]), p2=np.array([0.5, 0.5, 1.0]), n_samples=41)
    assert np.isclose(scl.membrane[0], 100.0, atol=1e-12)
    assert np.isclose(scl.bending_surface1[0], -60.0, atol=1e-12)
    assert np.isclose(scl.bending_surface2[0], 60.0, atol=1e-12)
    assert np.allclose(scl.peak_surface1, 0.0, atol=1e-12)
    assert np.allclose(scl.peak_surface2, 0.0, atol=1e-12)

    # Error norm delegation
    err = rv.error_norm(
        analytic_fn=lambda p: np.column_stack([0.01 * p[:, 0], -0.005 * p[:, 1], 0.002 * p[:, 2]]),
        field_name="U",
        cell_bounds=bounds,
    )
    assert err.abs_error < 1e-12
