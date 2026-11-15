"""Kratos-independent Phase 3 query package for WHT_IGA.

Provides :class:`ResultView` for unified read-only querying of ``VTKHDF``
and ``SQLite`` simulation outputs.
"""

from __future__ import annotations

from wht_iga.query.results import (
    FIELD_ALIASES,
    EquilibriumReport,
    ResultView,
    StressLinearizationResult,
    resolve_field_name,
    stress_to_matrix,
    von_mises_stress,
)

__all__ = [
    "FIELD_ALIASES",
    "EquilibriumReport",
    "ResultView",
    "StressLinearizationResult",
    "resolve_field_name",
    "stress_to_matrix",
    "von_mises_stress",
]
