"""``SolverBackendPort``: the swappable solver adapter (plan v2.1 §8.1 item 3).

The same compiled plan runs on ``KratosPythonBackend`` (production) or
``ReferenceNumbaBackend`` (oracle) without code changes; the session only talks
to this protocol.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from wht_iga.ports.contracts import BackendCapabilities, BackendSolution, StaticJob


@runtime_checkable
class SolverBackendPort(Protocol):
    """Solver adapter contract."""

    name: str

    def capabilities(self) -> BackendCapabilities:
        """Static description of supported features."""
        ...

    def solve_static(self, job: StaticJob) -> BackendSolution:
        """Solve one linear static step."""
        ...

    def close(self) -> None:
        """Release native resources (idempotent)."""
        ...
