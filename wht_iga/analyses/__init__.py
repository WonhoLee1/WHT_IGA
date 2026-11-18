"""High-level Phase 4 WHT_IGA analysis workflows and selectable discretization methods."""

from __future__ import annotations

from wht_iga.analyses import presets as presets
from wht_iga.analyses.methods_options import (
    DiscretizationOption,
    FullIgaOption,
    IgaFeHybridOption,
    MethodDescriptor,
    MethodId,
    NativeFeOption,
    ThbAdaptiveIgaOption,
    list_discretization_methods,
    parse_discretization_option,
)
from wht_iga.analyses.methods_recommend import (
    DiscretizationRecommendation,
    PhysicsKind,
    recommend_discretization,
)
from wht_iga.analyses.methods_solve import (
    BoxMethodProblem,
    MethodSolution,
    compare_method_solutions,
    solve_with_discretization_method,
)
from wht_iga.analyses.presets import (
    bending_thb,
    bulk_linear,
    contact_iga,
    describe,
)

__all__: list[str] = [
    "BoxMethodProblem",
    "DiscretizationOption",
    "DiscretizationRecommendation",
    "FullIgaOption",
    "IgaFeHybridOption",
    "MethodDescriptor",
    "MethodId",
    "MethodSolution",
    "NativeFeOption",
    "PhysicsKind",
    "ThbAdaptiveIgaOption",
    "bending_thb",
    "bulk_linear",
    "compare_method_solutions",
    "contact_iga",
    "describe",
    "list_discretization_methods",
    "parse_discretization_option",
    "presets",
    "recommend_discretization",
    "solve_with_discretization_method",
]
