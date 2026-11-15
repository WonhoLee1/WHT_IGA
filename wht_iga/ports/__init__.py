"""Layer 3: Hexagonal ports (``typing.Protocol``) and the data contracts they exchange."""

from __future__ import annotations

from wht_iga.ports.contracts import (
    AnalysisCancelled,
    BackendCapabilities,
    BackendSolution,
    CancelToken,
    MeshModel,
    ProgressCallback,
    ProgressEvent,
    StaticJob,
)
from wht_iga.ports.mesh_engine import MeshEnginePort
from wht_iga.ports.solver_backend import SolverBackendPort

__all__ = [
    "AnalysisCancelled", "BackendCapabilities", "BackendSolution", "CancelToken", "MeshEnginePort",
    "MeshModel", "ProgressCallback", "ProgressEvent", "SolverBackendPort", "StaticJob",
]
