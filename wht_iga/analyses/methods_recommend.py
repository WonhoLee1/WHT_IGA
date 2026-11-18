"""Physics-based degree recommendation and RAM-budget guard (split from methods.py)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from wht_iga.analyses.methods_options import Degree, FeElementType


PhysicsKind = Literal[
    "contact",
    "bending",
    "thermal",
    "bulk",
    "explicit_dynamics",
    "general",
]

_FE_IMPLICIT_ORDER: dict[str, int] = {
    "hex8": 1,
    "tet4": 1,
    "wedge6": 1,
    "hex20": 2,
    "tet10": 2,
}

_PHYSICS_BASE_DEGREE: dict[str, int] = {
    "contact": 2,
    "bending": 2,
    "thermal": 1,
    "bulk": 1,
    "explicit_dynamics": 1,
    "general": 2,
}


@dataclass(frozen=True)
class DiscretizationRecommendation:
    """Result of :func:`recommend_discretization`.

    Attributes
    ----------
    degree:
        Recommended (and capped) spline degree in ``[1, 3]``.
    reasons:
        Human-readable reasons ordered by application sequence.
    capped:
        Whether a user-requested degree was lowered.
    warnings:
        Non-fatal advisories (hybrid mismatch, budget overrun at p=1, ...).
    """

    degree: int
    reasons: tuple[str, ...] = ()
    capped: bool = False
    warnings: tuple[str, ...] = ()


def _estimate_peak_bytes(n_cells: tuple[int, int, int], degree: int) -> int:
    """Conservative direct-solver peak estimate from the ``(2p+1)^3`` stencil."""
    from wht_iga.solver.memory import est_pardiso_bytes

    n_dof = 3
    for nc in n_cells:
        n_dof *= int(nc) + int(degree)
    nnz = n_dof * (2 * int(degree) + 1) ** 3
    return int(est_pardiso_bytes(n_dof, nnz, degree))


def recommend_discretization(
    physics: PhysicsKind,
    requested_degree: int | None = None,
    *,
    fe_element_type: FeElementType | None = None,
    n_cells: tuple[int, int, int] | None = None,
    budget_bytes: float | None = None,
) -> DiscretizationRecommendation:
    """Recommend (and cap) the spline degree for a physics + resource context.

    Rules (all pure, no Kratos):
    1. Base degree by physics: contact/bending/general → 2, thermal/bulk/explicit → 1.
    2. Hard cap ``[1, 3]``: requests above 3 are capped with ``capped=True``;
       requests below 1 raise ``ValueError``.
    3. Hybrid pairing guard: when ``fe_element_type`` is given, the implicit FE
       order gap ``degree - fe_order`` may not exceed 1 (linear FE caps IGA at
       p=2); larger gaps are capped with a warning.
    4. RAM budget guard: when ``n_cells`` and ``budget_bytes`` are both given,
       the ``(2p+1)^3``-stencil Pardiso peak estimate must fit; otherwise the
       degree steps down 3 → 2 → 1 (``capped=True``). Overrun even at p=1 yields
       a warning suggesting THB or a coarser mesh.
    """
    if physics not in _PHYSICS_BASE_DEGREE:
        raise ValueError(
            f"unknown physics={physics!r}; expected one of {sorted(_PHYSICS_BASE_DEGREE)}"
        )
    reasons: list[str] = [f"base degree for physics={physics}"]
    warnings: list[str] = []
    degree = _PHYSICS_BASE_DEGREE[physics]
    capped = False

    if requested_degree is not None:
        requested_degree = int(requested_degree)
        if requested_degree < 1:
            raise ValueError(f"requested_degree must be >= 1, got {requested_degree}")
        if requested_degree > 3:
            reasons.append(f"requested p={requested_degree} capped to hard maximum 3")
            requested_degree = 3
            capped = True
        degree = requested_degree

    if fe_element_type is not None:
        if fe_element_type not in _FE_IMPLICIT_ORDER:
            raise ValueError(
                f"unknown fe_element_type={fe_element_type!r}; "
                f"expected one of {sorted(_FE_IMPLICIT_ORDER)}"
            )
        gap = int(degree) - _FE_IMPLICIT_ORDER[fe_element_type]
        if gap >= 2:
            new_degree = _FE_IMPLICIT_ORDER[fe_element_type] + 1
            reasons.append(
                f"hybrid pairing: {fe_element_type} (order "
                f"{_FE_IMPLICIT_ORDER[fe_element_type]}) caps IGA degree "
                f"{degree} -> {new_degree} (gap {gap} >= 2)"
            )
            warnings.append(
                f"IGA p={degree} with {fe_element_type} wastes DOFs: interface "
                "accuracy is limited by the linear FE side; THB local refinement "
                "is cheaper than global p elevation"
            )
            degree = new_degree
            capped = True
        elif gap < 0:
            warnings.append(
                f"IGA p={degree} below {fe_element_type} order "
                f"{_FE_IMPLICIT_ORDER[fe_element_type]}: consider raising IGA "
                "degree to match the FE side"
            )

    if n_cells is not None and budget_bytes is not None:
        budget_bytes = float(budget_bytes)
        if not np.isfinite(budget_bytes) or budget_bytes <= 0.0:
            raise ValueError(f"budget_bytes must be finite and > 0, got {budget_bytes}")
        while int(degree) > 1 and _estimate_peak_bytes(n_cells, int(degree)) > budget_bytes:
            reasons.append(
                f"p={degree} peak estimate exceeds budget {budget_bytes:.3g} B; stepping down"
            )
            degree = int(degree) - 1
            capped = True
        if _estimate_peak_bytes(n_cells, int(degree)) > budget_bytes:
            warnings.append(
                "peak estimate exceeds budget even at p=1; use THB local "
                "refinement or a coarser mesh instead of global p elevation"
            )

    return DiscretizationRecommendation(
        degree=int(degree),
        reasons=tuple(reasons),
        capped=capped,
        warnings=tuple(warnings),
    )
