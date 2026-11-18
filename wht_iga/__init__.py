"""WHT_IGA: immersed isogeometric analysis on Kratos Multiphysics.

Kratos is never imported at package import time; it is loaded lazily through
``wht_iga.kratos_backend``. The user API (plan v2.1 §8, D18) is exported lazily
(PEP 562), so ``import wht_iga`` stays light and does not import pydantic until
an API name is used::

    import wht_iga as wi
    study = wi.StudySpec(...)
    plan = wi.compile_plan(study)
    with wi.AnalysisSession(plan, backend=wi.KratosPythonBackend()) as s:
        results = s.run()
"""

from __future__ import annotations

import importlib
from typing import Any

__version__ = "0.0.1"

_LAZY: dict[str, str] = {
    # Layer 1: specs
    **{n: "wht_iga.spec" for n in (
        "StudySpec", "PartSpec", "MaterialSpec", "LinearElasticSpec", "ImmersedGridSpec", "BodyFitSpec",
        "RefinementSpec", "CadEntitiesRegionSpec", "SphereRegionSpec", "BoxRegionSpec", "GridFaceSetSpec",
        "BoxGeometrySpec", "StepGeometrySpec", "StaticStepSpec", "DisplacementBCSpec", "PressureLoadSpec",
        "TractionLoadSpec", "GravityLoadSpec", "PointProbeSpec", "LineProbeSpec", "MultiPointProbeSpec")},
    # Layer 2: plan
    **{n: "wht_iga.plan" for n in ("compile_plan", "Plan", "ExecutionPlan", "ResolveReport", "PlanError",
                                   "ArtifactPaths")},
    # Layer 3: ports / backends
    **{n: "wht_iga.ports" for n in ("SolverBackendPort", "MeshEnginePort", "CancelToken", "AnalysisCancelled",
                                    "ProgressEvent")},
    **{n: "wht_iga.backends" for n in ("KratosPythonBackend", "ReferenceNumbaBackend", "GridMeshEngine")},
    # Layer 4: runner
    **{n: "wht_iga.runner" for n in ("AnalysisSession", "StudyResults", "StepResult")},
    # Layer 5: query
    "ResultView": "wht_iga.query",
    # Selectable multi-method discretization options
    **{n: "wht_iga.analyses.methods" for n in (
        "BoxMethodProblem", "DiscretizationOption", "DiscretizationRecommendation",
        "FullIgaOption", "IgaFeHybridOption",
        "MethodDescriptor", "MethodId", "MethodSolution", "NativeFeOption",
        "ObbFrame", "PhysicsKind", "ThbAdaptiveIgaOption",
        "compare_method_solutions", "compute_pca_obb",
        "list_discretization_methods", "parse_discretization_option",
        "recommend_discretization", "solve_immersed_obb",
        "solve_with_discretization_method")},
    **{n: "wht_iga.analyses.presets" for n in (
        "bending_thb", "bulk_linear", "contact_iga", "describe")},
}

__all__ = ["__version__", *sorted(_LAZY)]


def __getattr__(name: str) -> Any:
    """Import user-API names on first access (PEP 562)."""
    mod = _LAZY.get(name)
    if mod is None:
        raise AttributeError(f"module 'wht_iga' has no attribute {name!r}")
    value = getattr(importlib.import_module(mod), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Module attributes including lazy user-API names."""
    return sorted(set(globals()) | set(_LAZY))
