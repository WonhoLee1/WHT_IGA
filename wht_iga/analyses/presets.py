# -*- coding: utf-8 -*-
"""Problem+option presets (usability layer, Track H Phase 6).

Each preset returns a ready ``(problem, option, recommendation)`` triple with
:func:`recommend_discretization` applied, so users need not hand-tune degree,
method, or budget caps. Pure constructors (no Kratos, no solves).
"""

from __future__ import annotations

from typing import Any

from wht_iga.analyses.methods_options import (
    DiscretizationOption,
    FullIgaOption,
    NativeFeOption,
    ThbAdaptiveIgaOption,
)
from wht_iga.analyses.methods_recommend import (
    DiscretizationRecommendation,
    recommend_discretization,
)
from wht_iga.analyses.methods_solve import BoxMethodProblem


def contact_iga(
    E: float,
    nu: float = 0.0,
    size: float = 1.0,
    traction: float = 2.0,
    n_cells: tuple[int, int, int] = (4, 2, 2),
    degree: int | None = None,
    budget_bytes: float | None = None,
) -> tuple[BoxMethodProblem, DiscretizationOption, DiscretizationRecommendation]:
    """Uniaxial contact-patch block: full IGA body-fit, p=2 default."""
    rec = recommend_discretization("contact", degree, n_cells=n_cells, budget_bytes=budget_bytes)
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0), hi=(10.0 * size, size, size), E=E, nu=nu,
        fixed_face=(0, 0), load_face=(0, 1), traction=(traction, 0.0, 0.0))
    option = FullIgaOption(kind="full_iga_body_fit", degree=rec.degree, n_cells=n_cells)
    return problem, option, rec


def bending_thb(
    E: float,
    nu: float = 0.0,
    degree: int | None = None,
    budget_bytes: float | None = None,
) -> tuple[BoxMethodProblem, DiscretizationOption, DiscretizationRecommendation]:
    """EB cantilever: THB adaptive, refinement at the clamped root."""
    rec = recommend_discretization("bending", degree, n_cells=(1, 1, 4), budget_bytes=budget_bytes)
    lo, hi = (-0.05, -0.05, 0.0), (0.05, 0.05, 1.0)
    problem = BoxMethodProblem(
        lo=lo, hi=hi, E=E, nu=nu,
        fixed_face=(2, 0), load_face=(2, 1), traction=(0.0, -500.0, 0.0))
    option = ThbAdaptiveIgaOption(
        degree=rec.degree, coarse_n_cells=(1, 1, 4), refinement_level=1,
        refine_box=((-0.05, -0.05, 0.0), (0.05, 0.05, 0.4)))
    return problem, option, rec


def bulk_linear(
    E: float,
    nu: float = 0.0,
    span: float = 10.0,
    traction: float = 10.0,
    fe_type: str = "hex20",
    degree: int | None = None,
) -> tuple[BoxMethodProblem, DiscretizationOption, DiscretizationRecommendation]:
    """Bulk uniaxial bar: native FE (cheapest verified path)."""
    rec = recommend_discretization("bulk", degree, fe_element_type=fe_type)  # type: ignore[arg-type]
    problem = BoxMethodProblem(
        lo=(0.0, 0.0, 0.0), hi=(span, 1.0, 1.0), E=E, nu=nu,
        fixed_face=(0, 0), load_face=(0, 1), traction=(traction, 0.0, 0.0))
    option = NativeFeOption(element_type=fe_type, n_cells=(4, 2, 2))  # type: ignore[arg-type]
    return problem, option, rec


def describe(preset: tuple[Any, Any, DiscretizationRecommendation]) -> str:
    """One-line human summary of a preset triple (logs / reports)."""
    problem, option, rec = preset
    return (f"{type(option).__name__}(p={rec.degree}) on box {problem.lo}->{problem.hi} "
            f"E={problem.E:g} [{'; '.join(rec.reasons)}]")
