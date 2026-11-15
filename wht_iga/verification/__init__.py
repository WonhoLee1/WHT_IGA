"""Kratos-independent verification schema, QoI, comparator, and reporting tools."""

from __future__ import annotations

from .comparator import (
    CheckItem,
    ComparisonReport,
    compare_against_tolerances,
    compare_csv_rows,
)
from .qoi import (
    LameSphereExact,
    analyze_dynamic_history,
    filter_participating_modes,
    lame_sphere_exact,
    line_profile_rel_l2,
    log_log_convergence_rate,
    principal_stresses,
    spherical_projections,
    time_window_rms_error,
)
from .reporter import (
    BIMETALLIC_TIE_HEADERS,
    E03_CONV_HEADERS,
    MULTIPLANAR_RESPONSES_HEADERS,
    ROARK_QOI_HEADERS,
    ROARK_REFINEMENT_HEADERS,
    SPHERE_QOI_HEADERS,
    V09_EIGEN_HEADERS,
    render_verification_markdown,
    write_coreform_csv,
)
from .schema import (
    VerificationCaseSpec,
    list_available_cases,
    load_case_spec,
    load_reference_csv,
)

__all__ = [
    "BIMETALLIC_TIE_HEADERS",
    "CheckItem",
    "ComparisonReport",
    "E03_CONV_HEADERS",
    "LameSphereExact",
    "MULTIPLANAR_RESPONSES_HEADERS",
    "ROARK_QOI_HEADERS",
    "ROARK_REFINEMENT_HEADERS",
    "SPHERE_QOI_HEADERS",
    "V09_EIGEN_HEADERS",
    "VerificationCaseSpec",
    "analyze_dynamic_history",
    "compare_against_tolerances",
    "compare_csv_rows",
    "filter_participating_modes",
    "lame_sphere_exact",
    "line_profile_rel_l2",
    "list_available_cases",
    "load_case_spec",
    "load_reference_csv",
    "log_log_convergence_rate",
    "principal_stresses",
    "render_verification_markdown",
    "spherical_projections",
    "time_window_rms_error",
    "write_coreform_csv",
]
