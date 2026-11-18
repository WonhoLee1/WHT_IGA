"""Kratos parity tests across all five selectable discretization methods."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses.methods import (
    BoxMethodProblem,
    FullIgaOption,
    IgaFeHybridOption,
    NativeFeOption,
    ThbAdaptiveIgaOption,
    solve_with_discretization_method,
)


@pytest.mark.kratos
def test_all_five_methods_uniaxial_patch_parity():
    """All 5 methods must reproduce exact uniaxial tension displacement u_x = t_x * x / E."""
    E_mod, nu, L, t_x = 2.0e5, 0.0, 10.0, 2.0
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(L, 1.0, 1.0),
        E=E_mod,
        nu=nu,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(t_x, 0.0, 0.0),
    )

    options = [
        FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(4, 2, 2)),
        FullIgaOption(kind="full_iga_immersed", degree=2, n_cells=(4, 2, 2), out_mode="drop"),
        ThbAdaptiveIgaOption(
            degree=2,
            coarse_n_cells=(2, 1, 1),
            refinement_level=1,
            refine_box=((0.0, 0.0, 0.0), (4.0, 1.0, 1.0)),
        ),
        IgaFeHybridOption(
            fe_element_type="hex20",
            fe_n_cells=(2, 2, 2),
            iga_degree=2,
            iga_n_cells=(2, 2, 2),
            split_axis=0,
            split_fraction=0.5,
            fe_side="lo",
            coupling_kind="auto",
            # auto must resolve Kratos-natively exact (mortar_l2), independent
            # of penalty_beta (left at default 1e3 to prove beta-independence).
        ),
        NativeFeOption(element_type="hex20", n_cells=(4, 2, 2)),
    ]

    pts = np.array(
        [
            [2.5, 0.5, 0.5],
            [5.0, 0.5, 0.5],
            [7.5, 0.5, 0.5],
            [10.0, 0.5, 0.5],
        ],
        dtype=np.float64,
    )
    expected_ux = t_x * pts[:, 0] / E_mod

    for opt in options:
        sol = solve_with_discretization_method(problem, opt)
        assert sol.converged is True
        assert sol.method == opt.kind
        assert sol.total_dofs > 0
        u_eval = sol.displacement_at(pts)
        np.testing.assert_allclose(u_eval[:, 0], expected_ux, rtol=1e-5, atol=1e-9)
        np.testing.assert_allclose(u_eval[:, 1:], 0.0, atol=1e-8)
        np.testing.assert_allclose(sol.tip_displacement[0], t_x * L / E_mod, rtol=1e-5)


@pytest.mark.kratos
def test_all_five_methods_cantilever_bending_parity():
    """All 5 methods must agree with V03 cantilever tip deflection (-1.0e-6 m)."""
    E_mod, nu = 2.0e11, 0.0
    problem = BoxMethodProblem(
        lo=(-0.05, -0.05, 0.0),
        hi=(0.05, 0.05, 1.0),
        E=E_mod,
        nu=nu,
        fixed_face=(2, 0),
        load_face=(2, 1),
        traction=(0.0, -500.0, 0.0),
    )
    eb_tip_uy = -1.0e-6

    options = [
        FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(2, 2, 8)),
        FullIgaOption(kind="full_iga_immersed", degree=2, n_cells=(2, 2, 8), out_mode="drop"),
        ThbAdaptiveIgaOption(
            degree=2,
            coarse_n_cells=(1, 1, 4),
            refinement_level=1,
            refine_box=((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4)),
        ),
        IgaFeHybridOption(
            fe_element_type="hex20",
            fe_n_cells=(2, 2, 4),
            iga_degree=2,
            iga_n_cells=(2, 2, 4),
            split_axis=2,
            split_fraction=0.5,
            fe_side="lo",
            coupling_kind="auto",
            penalty_beta=1e4,
        ),
        IgaFeHybridOption(
            fe_element_type="hex20",
            fe_n_cells=(2, 2, 4),
            iga_degree=2,
            iga_n_cells=(2, 2, 4),
            split_axis=2,
            split_fraction=0.5,
            fe_side="hi",
            coupling_kind="auto",
            penalty_beta=1e4,
        ),
        NativeFeOption(element_type="hex20", n_cells=(2, 2, 8)),
    ]

    tip_uy_values: list[float] = []
    for opt in options:
        sol = solve_with_discretization_method(problem, opt)
        assert sol.converged is True
        tip_uy = float(sol.tip_displacement[1])
        rel_err = abs(tip_uy - eb_tip_uy) / abs(eb_tip_uy)
        assert rel_err < 0.03, f"{opt.kind} tip_uy={tip_uy:.6e}, rel_err={rel_err:.4f}"
        tip_uy_values.append(tip_uy)

    ref_bf = tip_uy_values[0]
    for tip_val in tip_uy_values[1:]:
        assert abs(tip_val - ref_bf) / abs(ref_bf) < 0.03


@pytest.mark.kratos
@pytest.mark.parametrize("etype", ["hex8", "hex20", "tet4", "tet10", "wedge6"])
def test_native_fe_all_element_types_patch(etype: str):
    """NativeFeOption across all 5 canonical solid element types matches exact uniaxial solution."""
    E_mod, nu, L, t_x = 1.0e5, 0.0, 4.0, 10.0
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(L, 1.0, 1.0),
        E=E_mod,
        nu=nu,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(t_x, 0.0, 0.0),
    )
    sol = solve_with_discretization_method(problem, NativeFeOption(element_type=etype, n_cells=(2, 2, 2)))  # type: ignore[arg-type]
    expected_tip_ux = t_x * L / E_mod
    assert abs(float(sol.tip_displacement[0]) - expected_tip_ux) / expected_tip_ux < 1e-10


@pytest.mark.kratos
@pytest.mark.parametrize("coupling", ["penalty_spring", "mortar_l2", "barycentric"])
@pytest.mark.parametrize("fe_side", ["lo", "hi"])
def test_iga_fe_hybrid_couplings_and_sides(coupling: str, fe_side: str):
    """IgaFeHybridOption supports penalty_spring, mortar_l2, and barycentric on both sides."""
    E_mod, nu, L, t_x = 1.0e5, 0.0, 4.0, 10.0
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(L, 1.0, 1.0),
        E=E_mod,
        nu=nu,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(t_x, 0.0, 0.0),
    )
    opt = IgaFeHybridOption(
        fe_element_type="hex20",
        fe_n_cells=(2, 2, 2),
        iga_degree=2,
        iga_n_cells=(2, 2, 2),
        split_axis=0,
        split_fraction=0.5,
        fe_side=fe_side,  # type: ignore[arg-type]
        coupling_kind=coupling,  # type: ignore[arg-type]
        # Penalty accuracy scales as 1/beta (measured 6.25e-6 at 1e4 ->
        # 6.24e-8 at 1e6); beta=1e6 clears the 1e-6 gate with margin.
        # Barycentric is O(h)-consistent at the interface (nc=2 tip rel
        # 8.8e-4/1.66e-3, halving on refinement; weights and enforcement
        # verified exact, p=1/hex8 exact at 1.7e-18); mortar_l2 is the
        # exact reference. Gate barycentric at its characterized bound.
        penalty_beta=1e6 if coupling == "penalty_spring" else 1e4,
    )
    sol = solve_with_discretization_method(problem, opt)
    expected_tip_ux = t_x * L / E_mod
    tol = 2e-3 if coupling == "barycentric" else 1e-6
    assert abs(float(sol.tip_displacement[0]) - expected_tip_ux) / expected_tip_ux < tol


@pytest.mark.kratos
def test_immersed_obb_rotated_bar_invariance():
    """35-deg rotated bar solved via solve_immersed_obb matches exact global displacement."""
    from wht_iga.analyses.methods import compute_pca_obb, solve_immersed_obb
    from wht_iga.backends.mesh_engine import box_triangle_soup

    E_mod, nu, L, h, t_mag = 2.0e5, 0.0, 10.0, 1.0, 2.0
    theta = np.deg2rad(35.0)
    c, s = float(np.cos(theta)), float(np.sin(theta))
    rot_z = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    shift = np.array([2.5, -1.5, 0.5], dtype=np.float64)

    pts_local, tris = box_triangle_soup((-L / 2, -h / 2, -h / 2), (L / 2, h / 2, h / 2))
    pts_global = pts_local @ rot_z.T + shift

    # Determine which local OBB face is at -L/2 (fixed) and +L/2 (loaded) along axis 0
    frame, _, _ = compute_pca_obb(pts_global, padding_frac=0.0)
    axis0_global = np.asarray(frame.axes[0], dtype=np.float64)
    # Sign of axis0_global along rotated +x determines which end is lo vs hi
    sign = 1.0 if float(axis0_global @ rot_z[:, 0]) >= 0.0 else -1.0
    fixed_face = (0, 0) if sign > 0.0 else (0, 1)
    load_face = (0, 1) if sign > 0.0 else (0, 0)
    # Traction along rotated bar axis (+rot_z[:, 0] if sign > 0, else also +rot_z[:, 0] since load is tension)
    traction_global = tuple((t_mag * rot_z[:, 0]).tolist())

    sol = solve_immersed_obb(
        pts_global,
        tris,
        E=E_mod,
        nu=nu,
        n_cells=(4, 2, 2),
        degree=2,
        fixed_face_local=fixed_face,
        load_face_local=load_face,
        traction_global=traction_global,  # type: ignore[arg-type]
        padding_frac=0.0,
    )
    assert sol.converged is True
    assert sol.native["n_out_cells"] == 0

    expected_tip = (t_mag * L / E_mod) * rot_z[:, 0]
    np.testing.assert_allclose(sol.tip_displacement, expected_tip, rtol=1e-5, atol=1e-9)

    # Check interior points along the rotated bar centerline
    s_vals = np.array([2.5, 5.0, 7.5, 10.0])
    pts_eval = shift + np.outer(s_vals - L / 2, rot_z[:, 0])
    u_eval = sol.displacement_at(pts_eval)
    expected_u = np.outer(t_mag * s_vals / E_mod, rot_z[:, 0])
    np.testing.assert_allclose(u_eval, expected_u, rtol=1e-5, atol=1e-9)
