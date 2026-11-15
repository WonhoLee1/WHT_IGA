"""Coreform-compatible CSV writers and Markdown verification reporters."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Mapping, Sequence

from .comparator import ComparisonReport

ROARK_QOI_HEADERS: tuple[str, ...] = (
    "test_case",
    "h",
    "p",
    "peak_stress",
    "peak_rel_err_vs_roark",
    "peak_rel_err_vs_reference",
    "stress_profile_rel_l2",
    "ux_profile_rel_l2",
    "uy_profile_rel_l2",
)

ROARK_REFINEMENT_HEADERS: tuple[str, ...] = (
    "level",
    "finest_h",
    "local_test_case",
    "uniform_test_case",
    "local_dofs",
    "uniform_dofs",
    "local_peak_rel_err_vs_reference",
    "uniform_peak_rel_err_vs_reference",
    "local_peak_rel_err_vs_roark",
    "uniform_peak_rel_err_vs_roark",
    "local_stress_profile_rel_l2",
    "uniform_stress_profile_rel_l2",
    "num_cpus",
    "local_trim_wall_s",
    "uniform_trim_wall_s",
    "local_solve_wall_s",
    "uniform_solve_wall_s",
)

BIMETALLIC_TIE_HEADERS: tuple[str, ...] = (
    "p",
    "interface_workflow",
    "tip_uy",
    "rel_err_vs_transformed_section",
)

MULTIPLANAR_RESPONSES_HEADERS: tuple[str, ...] = (
    "cross_section",
    "p",
    "level",
    "h",
    "quantity",
    "measured_value",
    "reference_value",
    "relative_difference",
)

SPHERE_QOI_HEADERS: tuple[str, ...] = (
    "p",
    "level",
    "h",
    "quantity",
    "rel_l2_error_vs_lame",
    "tolerance",
)

V09_EIGEN_HEADERS: tuple[str, ...] = (
    "mode",
    "body_fitted_hz",
    "immersed_hz",
    "rel_diff",
    "participation_y",
    "participating",
)

E03_CONV_HEADERS: tuple[str, ...] = (
    "step",
    "step_name",
    "increment",
    "attempt",
    "accepted",
    "total_time",
    "trial_total_time",
    "step_time",
    "trial_step_time",
    "increment_size",
    "equilibrium_iterations",
    "severe_discontinuity_iterations",
    "total_iterations",
)


def write_coreform_csv(
    path: str | Path,
    headers: Sequence[str],
    rows: Sequence[Mapping[str, Any] | Sequence[Any]],
) -> Path:
    """Write a UTF-8 (no BOM) CSV with exact Coreform-compatible headers.

    Parameters
    ----------
    path:
        Output CSV file path. Parent directories are created automatically.
    headers:
        Ordered column names for the header row.
    rows:
        Sequence of row mappings (keyed by ``headers``) or positional sequences
        of the same length as ``headers``.

    Returns
    -------
    Path
        Resolved output path written to disk.
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header_list = list(headers)

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header_list)
        for row in rows:
            if isinstance(row, Mapping):
                writer.writerow([row.get(col, "") for col in header_list])
            else:
                row_list = list(row)
                if len(row_list) != len(header_list):
                    raise ValueError(
                        f"Row length {len(row_list)} does not match "
                        f"headers length {len(header_list)}."
                    )
                writer.writerow(row_list)

    return out_path


def render_verification_markdown(
    report: ComparisonReport,
    metadata: Mapping[str, Any] | None = None,
) -> str:
    """Render a standard verification report in GitHub-flavored Markdown.

    Parameters
    ----------
    report:
        Comparison report produced by :mod:`wht_iga.verification.comparator`.
    metadata:
        Optional run metadata (e.g. ``title``, ``dofs``, ``wall_s``,
        ``phase_gate``, ``mesh_size``).

    Returns
    -------
    str
        Formatted Markdown report suitable for ``docs/verification/<ID>.md``.
    """
    status_label = "PASS" if report.passed else "FAIL"
    lines: list[str] = [
        f"# Verification Report: {report.case_id}",
        "",
        f"- **Status:** {status_label}",
    ]
    if report.summary:
        lines.append(f"- **Summary:** {report.summary}")

    if metadata:
        lines.extend(["", "## Run Metadata", "", "| Key | Value |", "|---|---|"])
        for key, value in metadata.items():
            lines.append(f"| {key} | {value} |")

    lines.extend(
        [
            "",
            "## Comparison Checks",
            "",
            "| Check | Measured | Reference | Abs Diff | Rel Diff | Tolerance | Unit | Status |",
            "|---|---|---|---|---|---|---|---|",
        ]
    )
    for item in report.checks:
        item_status = "PASS" if item.passed else "FAIL"
        unit_str = item.unit if item.unit else "-"
        lines.append(
            f"| {item.name} "
            f"| {item.measured:.6g} "
            f"| {item.reference:.6g} "
            f"| {item.abs_diff:.6g} "
            f"| {item.rel_diff:.6g} "
            f"| {item.tol:.6g} "
            f"| {unit_str} "
            f"| {item_status} |"
        )

    lines.append("")
    return "\n".join(lines)
