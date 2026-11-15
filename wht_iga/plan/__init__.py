"""Layer 2: pre-flight execution plan and diagnostics (plan v2.1 §8.1, W3)."""

from __future__ import annotations

from wht_iga.plan.artifacts import ArtifactPaths, StepArtifacts
from wht_iga.plan.compiler import (
    COMPUTABLE_PROBE_VARIABLES,
    Diagnostic,
    ExecutionPlan,
    GridResolution,
    Plan,
    PlanError,
    ResolveReport,
    StepPrediction,
    compile_plan,
)

__all__ = [
    "ArtifactPaths", "COMPUTABLE_PROBE_VARIABLES", "Diagnostic", "ExecutionPlan", "GridResolution",
    "Plan", "PlanError", "ResolveReport", "StepArtifacts", "StepPrediction", "compile_plan",
]
