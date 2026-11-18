"""Kratos-free unit and validation tests for selectable discretization methods."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
from pydantic import ValidationError

from wht_iga.analyses import methods as methods_shim
from wht_iga.analyses.methods import (
    BoxMethodProblem,
    DiscretizationRecommendation,
    FullIgaOption,
    IgaFeHybridOption,
    MethodDescriptor,
    NativeFeOption,
    ThbAdaptiveIgaOption,
    compare_method_solutions,
    list_discretization_methods,
    parse_discretization_option,
    recommend_discretization,
)
from wht_iga.analyses.methods_solve import (
    MethodSolution,
    _eval_structured_fe_displacement,
    _structured_box_mesh,
)


def test_methods_shim_matches_canonical_homes():
    import wht_iga.analyses.methods_options as _o
    import wht_iga.analyses.methods_recommend as _r
    import wht_iga.analyses.methods_solve as _s
    assert methods_shim.FullIgaOption is _o.FullIgaOption
    assert methods_shim.solve_with_discretization_method is _s.solve_with_discretization_method
    assert methods_shim.recommend_discretization is _r.recommend_discretization
    assert methods_shim.compare_method_solutions is _s.compare_method_solutions
    assert methods_shim.BoxMethodProblem is _s.BoxMethodProblem
from wht_iga.kratos_backend.native_fe import ensure_positive_jacobian


def test_list_discretization_methods_metadata():
    methods = list_discretization_methods()
    assert len(methods) == 5
    assert all(isinstance(m, MethodDescriptor) for m in methods)
    ids = [m.id for m in methods]
    assert ids == [
        "full_iga_body_fit",
        "full_iga_immersed",
        "thb_adaptive_iga",
        "iga_fe_hybrid",
        "native_fe",
    ]
    by_id = {m.id: m for m in methods}
    assert by_id["full_iga_body_fit"].uses_bspline is True
    assert by_id["full_iga_immersed"].uses_bspline is True
    assert by_id["thb_adaptive_iga"].supports_adaptive_refinement is True
    assert by_id["iga_fe_hybrid"].uses_bspline is True
    assert by_id["iga_fe_hybrid"].supports_fe_bulk is True
    assert by_id["native_fe"].uses_bspline is False
    assert by_id["native_fe"].supports_fe_bulk is True
    for m in methods:
        assert m.korean_summary.strip() != ""
        assert m.engine_module.startswith("wht_iga.")


def test_option_models_frozen_and_extra_forbidden():
    opt = FullIgaOption(kind="full_iga_body_fit", degree=2, n_cells=(4, 2, 2))
    with pytest.raises(ValidationError):
        opt.degree = 3  # type: ignore[misc]
    with pytest.raises(ValidationError):
        FullIgaOption(kind="full_iga_body_fit", unknown_kw=1)  # type: ignore[call-arg]
    with pytest.raises(ValidationError):
        ThbAdaptiveIgaOption(typo_field=True)  # type: ignore[call-arg]
    with pytest.raises(ValidationError):
        NativeFeOption(element_type="hex20", bogus=42)  # type: ignore[call-arg]
    with pytest.raises(ValidationError):
        IgaFeHybridOption(fe_side="lo", extra_arg=1)  # type: ignore[call-arg]


def test_option_validation_rules():
    with pytest.raises(ValidationError):
        FullIgaOption(degree=0)
    with pytest.raises(ValidationError):
        FullIgaOption(degree=4)
    with pytest.raises(ValidationError, match="continuity"):
        FullIgaOption(degree=2, continuity=2)
    with pytest.raises(ValidationError):
        FullIgaOption(kind="full_iga_immersed", ghost_alpha=0.0)
    with pytest.raises(ValidationError):
        ThbAdaptiveIgaOption(refinement_level=0)
    with pytest.raises(ValidationError):
        ThbAdaptiveIgaOption(refinement_level=6)
    with pytest.raises(ValidationError, match="refine_box"):
        ThbAdaptiveIgaOption(refine_box=((1.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
    with pytest.raises(ValidationError):
        NativeFeOption(element_type="hex27")  # type: ignore[arg-type]
    with pytest.raises(ValidationError):
        IgaFeHybridOption(split_fraction=0.0)
    with pytest.raises(ValidationError):
        IgaFeHybridOption(split_fraction=1.0)
    with pytest.raises(ValidationError):
        IgaFeHybridOption(split_axis=3)  # type: ignore[arg-type]
    with pytest.raises(ValidationError):
        IgaFeHybridOption(fe_side="mid")  # type: ignore[arg-type]


def test_option_json_roundtrip_and_discriminated_union():
    samples = [
        FullIgaOption(kind="full_iga_body_fit", degree=3, continuity=1, n_cells=(6, 2, 2)),
        FullIgaOption(
            kind="full_iga_immersed",
            degree=2,
            n_cells=(4, 2, 2),
            out_mode="ghost",
            ghost_alpha=1e-6,
            trim_deflection=0.05,
        ),
        ThbAdaptiveIgaOption(
            degree=2,
            coarse_n_cells=(2, 2, 4),
            refinement_level=2,
            refine_box=((0.0, 0.0, 0.0), (1.0, 1.0, 0.4)),
            balance="two_to_one",
            penalty_beta=1e9,
        ),
        NativeFeOption(element_type="tet10", n_cells=(3, 2, 2)),
        IgaFeHybridOption(
            fe_element_type="hex20",
            fe_n_cells=(2, 2, 2),
            iga_degree=2,
            iga_n_cells=(2, 2, 3),
            split_axis=2,
            split_fraction=0.4,
            fe_side="hi",
            coupling_kind="mortar_l2",
            penalty_beta=5e3,
        ),
    ]
    for opt in samples:
        raw_json = opt.to_json()
        parsed_json = parse_discretization_option(raw_json)
        assert parsed_json == opt
        parsed_dict = parse_discretization_option(opt.model_dump())
        assert parsed_dict == opt

    updated = samples[0].with_updates(degree=2, continuity=1)
    assert updated.degree == 2 and samples[0].degree == 3
    assert samples[1].resolved_quadrature_degree == 4


def test_box_method_problem_validation():
    prob = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0),
        hi=(10.0, 1.0, 2.0),
        E=2.0e5,
        nu=0.3,
        fixed_face=(0, 0),
        load_face=(0, 1),
        traction=(2.0, 0.0, 0.0),
    )
    np.testing.assert_allclose(prob.tip_point(), [10.0, 0.5, 1.0])

    with pytest.raises(ValueError, match="hi must exceed lo"):
        BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(0.0, 1.0, 1.0), E=1e5, nu=0.3)
    with pytest.raises(ValueError, match="E must be"):
        BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), E=-1.0, nu=0.3)
    with pytest.raises(ValueError, match="nu must be"):
        BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), E=1e5, nu=0.5)
    with pytest.raises(ValueError, match="fixed_face"):
        BoxMethodProblem(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), E=1e5, nu=0.3, fixed_face=(3, 0))


@pytest.mark.parametrize("etype", ["hex8", "hex20", "tet4", "tet10", "wedge6"])
def test_structured_box_mesh_jacobians_and_interpolation(etype: str):
    lo = (-1.0, 0.5, 2.0)
    hi = (3.0, 2.5, 6.0)
    n_cells = (2, 2, 2)
    mesh = _structured_box_mesh(lo, hi, n_cells, etype)  # type: ignore[arg-type]
    for row in mesh.blocks[0].conn:
        _, flipped = ensure_positive_jacobian(etype, row, mesh.nodes, auto_flip=False)
        assert flipped is False

    # Linear field is reproduced exactly by all 5 element types; quadratic by hex20/tet10.
    coords = mesh.nodes
    if etype in ("hex20", "tet10"):
        u_nodal = np.column_stack([
            0.1 * coords[:, 0] ** 2 - 0.2 * coords[:, 1] * coords[:, 2] + 0.5,
            0.3 * coords[:, 1] ** 2 + 0.1 * coords[:, 0],
            -0.05 * coords[:, 2] ** 2 + 0.4 * coords[:, 0] * coords[:, 1],
        ])
        expected_fn = lambda p: np.column_stack([
            0.1 * p[:, 0] ** 2 - 0.2 * p[:, 1] * p[:, 2] + 0.5,
            0.3 * p[:, 1] ** 2 + 0.1 * p[:, 0],
            -0.05 * p[:, 2] ** 2 + 0.4 * p[:, 0] * p[:, 1],
        ])
    else:
        u_nodal = np.column_stack([
            0.2 * coords[:, 0] - 0.1 * coords[:, 1] + 0.3,
            0.4 * coords[:, 1] + 0.05 * coords[:, 2],
            -0.15 * coords[:, 2] + 0.25 * coords[:, 0],
        ])
        expected_fn = lambda p: np.column_stack([
            0.2 * p[:, 0] - 0.1 * p[:, 1] + 0.3,
            0.4 * p[:, 1] + 0.05 * p[:, 2],
            -0.15 * p[:, 2] + 0.25 * p[:, 0],
        ])

    rng = np.random.default_rng(12345)
    query = rng.uniform(low=lo, high=hi, size=(12, 3))
    query = np.vstack([query, np.array([[3.0, 1.5, 4.0], [-1.0, 0.5, 2.0]])])
    u_eval = _eval_structured_fe_displacement(
        mesh,
        etype,  # type: ignore[arg-type]
        lo,
        hi,
        n_cells,
        u_nodal,
        query,
    )
    np.testing.assert_allclose(u_eval, expected_fn(query), rtol=1e-12, atol=1e-12)


def test_lazy_exports_stay_kratos_free():
    code = (
        "import sys, wht_iga\n"
        "assert 'pydantic' not in sys.modules\n"
        "assert len(wht_iga.list_discretization_methods()) == 5\n"
        "from wht_iga.analyses import methods as m\n"
        "assert wht_iga.FullIgaOption is m.FullIgaOption\n"
        "assert wht_iga.IgaFeHybridOption is m.IgaFeHybridOption\n"
        "assert wht_iga.solve_with_discretization_method is m.solve_with_discretization_method\n"
        "bad = [k for k in sys.modules if k.split('.')[0] in ('KratosMultiphysics', 'aieng', 'wht_topopt', 'gmsh')]\n"
        "assert not bad, bad\n"
    )
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr


def test_recommend_physics_defaults():
    assert recommend_discretization("contact").degree == 2
    assert recommend_discretization("bending").degree == 2
    assert recommend_discretization("general").degree == 2
    assert recommend_discretization("thermal").degree == 1
    assert recommend_discretization("bulk").degree == 1
    assert recommend_discretization("explicit_dynamics").degree == 1
    r = recommend_discretization("contact")
    assert r.capped is False and r.warnings == ()


def test_recommend_hard_cap_and_invalid():
    r = recommend_discretization("bending", 5)
    assert r.degree == 3 and r.capped is True
    with pytest.raises(ValueError):
        recommend_discretization("bending", 0)
    with pytest.raises(ValueError):
        recommend_discretization("nope")


def test_recommend_hybrid_pairing_guard():
    r = recommend_discretization("contact", 3, fe_element_type="hex8")
    assert r.degree == 2 and r.capped is True and len(r.warnings) == 1
    r2 = recommend_discretization("contact", 2, fe_element_type="hex8")
    assert r2.degree == 2 and r2.capped is False
    r3 = recommend_discretization("contact", 3, fe_element_type="hex20")
    assert r3.degree == 3 and r3.capped is False
    r4 = recommend_discretization("bulk", 1, fe_element_type="hex20")
    assert r4.degree == 1 and any("below" in w for w in r4.warnings)


def test_recommend_budget_guard_steps_down():
    tiny = 1.0  # 1 byte: forces step-down to p=1 plus warning
    r = recommend_discretization("contact", 3, n_cells=(4, 4, 4), budget_bytes=tiny)
    assert r.degree == 1 and r.capped is True and len(r.warnings) == 1
    huge = 1e18
    r2 = recommend_discretization("contact", 3, n_cells=(4, 4, 4), budget_bytes=huge)
    assert r2.degree == 3 and r2.capped is False
    r3 = recommend_discretization("thermal", None, n_cells=(4, 4, 4), budget_bytes=huge)
    assert r3.degree == 1
    with pytest.raises(ValueError):
        recommend_discretization("bulk", 2, n_cells=(2, 2, 2), budget_bytes=0.0)


def test_compare_method_solutions_pure():
    import numpy as np

    from wht_iga.analyses.methods import MethodSolution

    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    ref = MethodSolution(
        method="native_fe", total_dofs=3, active_dofs=3,
        tip_displacement=np.array([1.0, 0.0, 0.0]), max_displacement=1.0,
        converged=True, _eval_fn=lambda p: np.column_stack([p[:, 0], p[:, 1], p[:, 2]]),
    )
    other = MethodSolution(
        method="native_fe", total_dofs=3, active_dofs=3,
        tip_displacement=np.array([1.1, 0.0, 0.0]), max_displacement=1.1,
        converged=True, _eval_fn=lambda p: np.column_stack([1.1 * p[:, 0], p[:, 1], p[:, 2]]),
    )
    out = compare_method_solutions(ref, other, pts)
    assert out["tip_abs"] == pytest.approx(0.1)
    assert out["tip_rel"] == pytest.approx(0.1)
    assert out["max_abs"] == pytest.approx(0.1)
    assert out["l2_abs"] == pytest.approx(0.1)
    assert out["l2_rel"] == pytest.approx(0.1 / 1.0)
    same = compare_method_solutions(ref, ref, pts)
    assert same["l2_abs"] == 0.0 and same["tip_rel"] == 0.0


def test_compute_pca_obb_canonical_and_transforms():
    from wht_iga.analyses.methods import ObbFrame, compute_pca_obb

    # 35-deg rotated box shifted to (3.0, -2.0, 1.5)
    theta = np.deg2rad(35.0)
    c, s = float(np.cos(theta)), float(np.sin(theta))
    rot_z = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    L, W, H = 10.0, 2.0, 1.0
    local_corners = np.array(
        [[x, y, z] for x in (-L / 2, L / 2) for y in (-W / 2, W / 2) for z in (-H / 2, H / 2)],
        dtype=np.float64,
    )
    shift = np.array([3.0, -2.0, 1.5], dtype=np.float64)
    global_corners = local_corners @ rot_z.T + shift

    frame, lo, hi = compute_pca_obb(global_corners, padding_frac=0.02)
    assert isinstance(frame, ObbFrame)
    R = frame.rotation_matrix()
    np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(frame.origin, shift, atol=1e-12)

    # Spans match (L, W, H) * (1 + 2 * 0.02)
    spans = np.array(hi) - np.array(lo)
    np.testing.assert_allclose(spans, np.array([L, W, H]) * 1.04, rtol=1e-12)

    # Canonical sign & permutation invariance
    perm = np.array([5, 2, 7, 0, 4, 1, 6, 3])
    frame_perm, lo_perm, hi_perm = compute_pca_obb(global_corners[perm], padding_frac=0.02)
    np.testing.assert_allclose(frame_perm.rotation_matrix(), R, atol=1e-12)
    np.testing.assert_allclose(lo_perm, lo, atol=1e-12)
    np.testing.assert_allclose(hi_perm, hi, atol=1e-12)

    # Point, vector, and tensor round-trip transforms
    pts_loc = frame.to_local_points(global_corners)
    np.testing.assert_allclose(frame.to_global_points(pts_loc), global_corners, atol=1e-12)
    vecs = np.array([[1.0, 2.0, -0.5], [0.0, -3.0, 4.0]])
    np.testing.assert_allclose(frame.to_global_vectors(frame.to_local_vectors(vecs)), vecs, atol=1e-12)
    sig_loc = np.diag([100.0, 20.0, -5.0])
    sig_glob = frame.to_global_tensors(sig_loc)
    np.testing.assert_allclose(sig_glob, R @ sig_loc @ R.T, atol=1e-12)


def test_compute_pca_obb_validation_errors():
    from wht_iga.analyses.methods import compute_pca_obb

    with pytest.raises(ValueError, match="shape"):
        compute_pca_obb(np.zeros((3, 3)))
    with pytest.raises(ValueError, match="padding_frac"):
        compute_pca_obb(np.eye(4, 3), padding_frac=-0.1)
    # Coplanar (rank 2) points must raise
    coplanar = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    with pytest.raises(ValueError, match="rank < 3"):
        compute_pca_obb(coplanar)
