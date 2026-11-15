"""Layer 3 adapters: solver backends and the default mesh engine.

Importing this package never imports Kratos; ``KratosPythonBackend`` loads it
lazily in ``solve_static``.
"""

from __future__ import annotations

from wht_iga.backends.kratos_python import KratosPythonBackend
from wht_iga.backends.mesh_engine import GridMeshEngine, box_triangle_soup
from wht_iga.backends.reference_numba import BackendUnsupported, ReferenceNumbaBackend, body_force_vector

__all__ = ["BackendUnsupported", "GridMeshEngine", "KratosPythonBackend", "ReferenceNumbaBackend",
           "body_force_vector", "box_triangle_soup"]
