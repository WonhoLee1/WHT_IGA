"""Issue log shared by S1/S3/S4 geometry stages.

Severity ``error`` blocks :func:`validate`; ``warning`` is recorded only.
Mirrors the ``issue`` table of ``<job>.wiga.sqlite`` (v2 §4.2).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Issue:
    """Single validation issue."""

    severity: str  # "error" | "warning" | "info"
    code: str  # e.g. "NON_WATERTIGHT", "QUAD_FALLBACK", "THIN_FEATURE"
    message: str
    entity: str = ""


class IssueLog:
    """Ordered container of :class:`Issue` with a validate gate."""

    def __init__(self) -> None:
        self._issues: list[Issue] = []

    def add(self, severity: str, code: str, message: str, entity: str = "") -> None:
        """Append an issue entry."""
        self._issues.append(Issue(severity, code, message, entity))

    def __len__(self) -> int:
        return len(self._issues)

    def __iter__(self):  # type: ignore[no-untyped-def]
        return iter(self._issues)

    def errors(self) -> list[Issue]:
        """Return entries with severity ``error``."""
        return [i for i in self._issues if i.severity == "error"]

    def validate(self) -> None:
        """Raise ``ValueError`` if any error-severity issue exists."""
        errs = self.errors()
        if errs:
            detail = "; ".join(f"[{e.code}] {e.message}" for e in errs)
            raise ValueError(f"geometry validation failed: {detail}")

    def to_rows(self) -> list[tuple[str, str, str, str]]:
        """Return rows matching the ``issue`` DB table layout."""
        return [(i.severity, i.code, i.message, i.entity) for i in self._issues]
