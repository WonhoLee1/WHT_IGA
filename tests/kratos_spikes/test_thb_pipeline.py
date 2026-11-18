# -*- coding: utf-8 -*-
"""THB adaptive refinement through the solver pipeline (Track H Phase 2).

Go/No-Go: uniaxial exact + bending 3% via ``JobSpec.thb_refinement``.
Unsupported configs (named-set regions, level > 5) must raise, not silently
fall back to uniform (no-silent-substitute).
"""

from __future__ import annotations

import pytest

from wht_iga.solver.validate import JobSpec, TractionLoad


def _base_job(**kw):
    args: dict = dict(
        name="thb_pipe",
        lo=(0.0, 0.0, 0.0),
        hi=(10.0, 1.0, 1.0),
        n_cells=(4, 2, 2),
        degree=2,
        E=2.0e5,
        nu=0.0,
        fixed_faces=[(0, 0)],
        tractions=[TractionLoad(face=(0, 1), vector=(2.0, 0.0, 0.0))],
    )
    args.update(kw)
    return args


def test_thb_refinement_field_roundtrip():
    """JobSpec carries thb_refinement through to_json/from_json."""
    spec = JobSpec(
        **_base_job(thb_refinement={"level": 1, "balance": "two_to_one",
                                    "boxes": [((0.0, 0.0, 0.0), (4.0, 1.0, 1.0))]}))
    assert spec.thb_refinement is not None
    assert spec.thb_refinement["level"] == 1
    back = JobSpec.from_json(spec.to_json())
    assert back.thb_refinement is not None
    assert back.thb_refinement["level"] == 1


def test_thb_pipeline_uniaxial_exact():
    """THB pipeline reproduces uniaxial tension (machine precision)."""
    from wht_iga.solver.pipeline import run_linear

    spec = JobSpec(
        **_base_job(thb_refinement={"level": 1, "balance": "two_to_one",
                                    "boxes": [((0.0, 0.0, 0.0), (4.0, 1.0, 1.0))]}))
    res = run_linear(spec)
    assert res.converged is True
    import numpy as np

    box_fine = res.keepalive[1]
    u_tip = box_fine.evaluate(np.asarray(res.U).reshape(-1, 3), np.array([[10.0, 0.5, 0.5]]))
    assert abs(float(u_tip[0, 0]) - 2.0 * 10.0 / 2.0e5) / (2.0 * 10.0 / 2.0e5) < 1e-9, u_tip


def test_thb_pipeline_bending_gate():
    """THB pipeline cantilever bending within 3% (EB tip -1e-6 m)."""
    from wht_iga.solver.pipeline import run_linear

    spec = JobSpec(
        name="thb_bend",
        lo=(-0.05, -0.05, 0.0),
        hi=(0.05, 0.05, 1.0),
        n_cells=(1, 1, 4),
        degree=2,
        E=2.0e11,
        nu=0.0,
        fixed_faces=[(2, 0)],
        tractions=[TractionLoad(face=(2, 1), vector=(0.0, -500.0, 0.0))],
        thb_refinement={"level": 1, "balance": "two_to_one",
                        "boxes": [((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4))]},
    )
    res = run_linear(spec)
    assert res.converged is True
    import numpy as np

    box_fine = res.keepalive[1]
    u_tip = box_fine.evaluate(np.asarray(res.U).reshape(-1, 3), np.array([[0.0, 0.0, 1.0]]))
    assert abs(abs(float(u_tip[0, 1])) - 1.0e-6) / 1.0e-6 < 0.03, u_tip


def test_thb_pipeline_rejects_set_regions():
    """Named-set refinements are out of scope -> explicit error, no silent fallback."""
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.validate import ValidationError

    spec = JobSpec(
        **_base_job(thb_refinement={"level": 1, "sets": ["patch_a"]}))
    with pytest.raises(ValidationError):
        run_linear(spec)


def test_thb_pipeline_sphere_region_exact():
    """Sphere refinement path solves exactly on uniaxial tension."""
    import numpy as np

    from wht_iga.solver.pipeline import run_linear

    spec = JobSpec(
        **_base_job(thb_refinement={"level": 1, "balance": "two_to_one",
                                    "spheres": [((5.0, 0.5, 0.5), 1.5)]}))
    res = run_linear(spec)
    assert res.converged is True
    box_fine = res.keepalive[1]
    u_tip = box_fine.evaluate(np.asarray(res.U).reshape(-1, 3), np.array([[10.0, 0.5, 0.5]]))
    assert abs(float(u_tip[0, 0]) - 2.0 * 10.0 / 2.0e5) / (2.0 * 10.0 / 2.0e5) < 1e-9, u_tip
