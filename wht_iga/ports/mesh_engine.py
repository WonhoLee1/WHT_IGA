"""``MeshEnginePort``: builds the background grid and cut-cell quadrature (W2)."""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

from wht_iga.ports.contracts import CancelToken, MeshModel, ProgressCallback


@runtime_checkable
class MeshEnginePort(Protocol):
    """Mesh engine contract."""

    def build(self, plan: Any, *, on_progress: ProgressCallback | None = None,
              cancel: CancelToken | None = None) -> MeshModel:
        """Return the :class:`MeshModel` for a compiled plan."""
        ...
