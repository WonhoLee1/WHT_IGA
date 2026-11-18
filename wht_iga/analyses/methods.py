"""Selectable multi-method discretization options and unified solver facade.

Supports five first-class structural discretization strategies over a unified
specification and solver interface, without importing Kratos at module load time:

* ``"full_iga_body_fit"`` — Untrimmed body-fit B-spline IGA via ``wht_iga.solver.pipeline``.
* ``"full_iga_immersed"`` — Immersed/trimmed B-spline IGA (FCM / cut-cell quadrature).
* ``"thb_adaptive_iga"`` — Locally refined Truncated Hierarchical B-splines (THB).
* ``"iga_fe_hybrid"`` — Coupled Native FEM + IGA domains across a split interface.
* ``"native_fe"`` — Standard Kratos solid finite elements (hex8/hex20/tet4/tet10/wedge6).
"""

from wht_iga.analyses.methods_options import (
    Degree,
    DiscretizationOption,
    FeElementType,
    FullIgaOption,
    HybridCouplingKind,
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


__all__ = [
    "BoxMethodProblem",
    "Degree",
    "DiscretizationOption",
    "DiscretizationRecommendation",
    "FeElementType",
    "FullIgaOption",
    "HybridCouplingKind",
    "IgaFeHybridOption",
    "MethodDescriptor",
    "MethodId",
    "MethodSolution",
    "NativeFeOption",
    "PhysicsKind",
    "ThbAdaptiveIgaOption",
    "compare_method_solutions",
    "list_discretization_methods",
    "parse_discretization_option",
    "recommend_discretization",
    "solve_with_discretization_method",
]
