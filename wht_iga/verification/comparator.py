"""Verification comparison engine for analytical, FE, and Coreform CSV checks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class CheckItem:
    """Single scalar verification check result.

    Parameters
    ----------
    name:
        Identifier of the checked quantity or row-column pair.
    measured:
        Measured numerical value from WHT_IGA.
    reference:
        Accepted analytical, FE, or Coreform reference value.
    rel_diff:
        Relative difference ``|measured - reference| / |reference|``.
    abs_diff:
        Absolute difference ``|measured - reference|``.
    tol:
        Acceptance tolerance threshold.
    passed:
        Whether the check satisfied ``tol``.
    unit:
        Optional physical unit label.
    """

    name: str
    measured: float
    reference: float
    rel_diff: float
    abs_diff: float
    tol: float
    passed: bool
    unit: str = ""


@dataclass(frozen=True)
class ComparisonReport:
    """Aggregated verification comparison report for a case.

    Parameters
    ----------
    case_id:
        Verification or example case identifier (e.g. ``"V06"``).
    passed:
        ``True`` if every item in ``checks`` passed.
    checks:
        List of individual :class:`CheckItem` results.
    summary:
        Human-readable summary string.
    """

    case_id: str
    passed: bool
    checks: list[CheckItem]
    summary: str = ""


def _compute_diffs(measured: float, reference: float) -> tuple[float, float]:
    abs_diff = abs(measured - reference)
    if reference == 0.0:
        rel_diff = 0.0 if abs_diff == 0.0 else float("inf")
    else:
        rel_diff = abs_diff / abs(reference)
    return rel_diff, abs_diff


def compare_against_tolerances(
    case_id: str,
    measured: Mapping[str, float],
    reference: Mapping[str, float],
    tolerances: Mapping[str, float],
    *,
    mode: str = "rel",
    units: Mapping[str, str] | None = None,
) -> ComparisonReport:
    """Compare measured scalar metrics against reference values and tolerances.

    Parameters
    ----------
    case_id:
        Case identifier for the resulting report.
    measured:
        Mapping of metric name to measured float value.
    reference:
        Mapping of metric name to reference float value.
    tolerances:
        Mapping of metric name to tolerance threshold.
    mode:
        Comparison mode: ``"rel"`` (relative difference ``<= tol``) or
        ``"abs"`` (absolute difference ``<= tol``).
    units:
        Optional mapping of metric name to unit string.

    Returns
    -------
    ComparisonReport
        Report containing per-metric :class:`CheckItem` objects and overall
        pass/fail status.
    """
    if mode not in ("rel", "abs"):
        raise ValueError(f"Unsupported tolerance mode '{mode}'; expected 'rel' or 'abs'.")

    checks: list[CheckItem] = []
    unit_map = units or {}

    for name, tol in tolerances.items():
        if name not in measured:
            raise KeyError(f"Missing measured value for metric '{name}'.")
        if name not in reference:
            raise KeyError(f"Missing reference value for metric '{name}'.")

        m_val = float(measured[name])
        r_val = float(reference[name])
        tol_val = float(tol)
        rel_diff, abs_diff = _compute_diffs(m_val, r_val)
        err_metric = rel_diff if mode == "rel" else abs_diff
        passed = err_metric <= tol_val

        checks.append(
            CheckItem(
                name=name,
                measured=m_val,
                reference=r_val,
                rel_diff=rel_diff,
                abs_diff=abs_diff,
                tol=tol_val,
                passed=passed,
                unit=str(unit_map.get(name, "")),
            )
        )

    n_pass = sum(1 for c in checks if c.passed)
    all_passed = n_pass == len(checks)
    summary = f"{case_id}: {n_pass}/{len(checks)} checks passed ({mode} tolerance)"
    return ComparisonReport(
        case_id=case_id,
        passed=all_passed,
        checks=checks,
        summary=summary,
    )


def compare_csv_rows(
    measured_rows: Sequence[Mapping[str, Any]],
    reference_rows: Sequence[Mapping[str, Any]],
    key_cols: Sequence[str],
    value_cols: Sequence[str],
    rel_tol: float | Mapping[str, float],
    *,
    case_id: str = "csv_comparison",
) -> ComparisonReport:
    """Compare measured CSV rows 1:1 against reference CSV rows by composite key.

    Parameters
    ----------
    measured_rows:
        Measured row dictionaries (e.g. from WHT_IGA runner output).
    reference_rows:
        Reference row dictionaries (e.g. from Coreform ``generated-downloads/``).
    key_cols:
        Column names uniquely identifying each row.
    value_cols:
        Numeric column names to compare.
    rel_tol:
        Scalar relative tolerance or per-column mapping of relative tolerances.
    case_id:
        Case identifier for the returned :class:`ComparisonReport`.

    Returns
    -------
    ComparisonReport
        Detailed per-key, per-column comparison report.
    """
    ref_index: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    for row in reference_rows:
        key = tuple(row[k] for k in key_cols)
        ref_index[key] = row

    checks: list[CheckItem] = []
    for m_row in measured_rows:
        key = tuple(m_row[k] for k in key_cols)
        key_label = ",".join(f"{k}={v}" for k, v in zip(key_cols, key))
        r_row = ref_index.get(key)
        if r_row is None:
            for col in value_cols:
                tol_val = (
                    float(rel_tol[col])
                    if isinstance(rel_tol, Mapping)
                    else float(rel_tol)
                )
                checks.append(
                    CheckItem(
                        name=f"{key_label}:{col}",
                        measured=float(m_row[col]),
                        reference=float("nan"),
                        rel_diff=float("inf"),
                        abs_diff=float("inf"),
                        tol=tol_val,
                        passed=False,
                    )
                )
            continue

        for col in value_cols:
            m_val = float(m_row[col])
            r_val = float(r_row[col])
            tol_val = (
                float(rel_tol[col])
                if isinstance(rel_tol, Mapping)
                else float(rel_tol)
            )
            rel_diff, abs_diff = _compute_diffs(m_val, r_val)
            checks.append(
                CheckItem(
                    name=f"{key_label}:{col}",
                    measured=m_val,
                    reference=r_val,
                    rel_diff=rel_diff,
                    abs_diff=abs_diff,
                    tol=tol_val,
                    passed=rel_diff <= tol_val,
                )
            )

    n_pass = sum(1 for c in checks if c.passed)
    all_passed = n_pass == len(checks) and len(checks) > 0
    summary = f"{case_id}: {n_pass}/{len(checks)} CSV cell checks passed"
    return ComparisonReport(
        case_id=case_id,
        passed=all_passed,
        checks=checks,
        summary=summary,
    )
