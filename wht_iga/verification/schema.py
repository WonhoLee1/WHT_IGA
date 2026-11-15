"""Case specification and reference CSV loaders for WHT_IGA verification.

Provides :class:`VerificationCaseSpec`, :func:`load_case_spec`,
:func:`list_available_cases`, and :func:`load_reference_csv` with full
backward compatibility for the existing ``v04``–``v12`` and ``e01``–``e05``
JSON specifications as well as ``v01``–``v03``.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

_DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "tests" / "verification" / "data"
_INT_PATTERN = re.compile(r"^[+-]?\d+$")


@dataclass(frozen=True)
class VerificationCaseSpec:
    """Parsed verification or example problem specification.

    Parameters
    ----------
    case_id:
        Normalized case identifier (e.g. ``"V04"`` or ``"E01"``).
    title:
        Human-readable problem title.
    category:
        Problem category (``"verification"`` or ``"example"``).
    status:
        Specification completeness status (e.g. ``"complete"``).
    raw:
        Full raw JSON dictionary for lossless access to problem-specific keys.
    geometry:
        Geometry specification dictionary.
    materials:
        Material properties dictionary.
    boundary_conditions:
        Boundary conditions specification dictionary.
    loads:
        Load definitions dictionary.
    mesh_sweep:
        Mesh and polynomial degree sweep specification dictionary.
    acceptance:
        Acceptance criteria and gate thresholds dictionary.
    reference_values:
        Analytical or Coreform reference values dictionary.
    provenance:
        Source documentation and reference CSV paths dictionary.
    meta:
        Optional standardized ``meta`` block.
    probes:
        Optional standardized ``probes`` block.
    sweeps:
        Optional standardized ``sweeps`` block.
    """

    case_id: str
    title: str
    category: str
    status: str
    raw: dict[str, Any] = field(repr=False)
    geometry: dict[str, Any] = field(default_factory=dict)
    materials: dict[str, Any] = field(default_factory=dict)
    boundary_conditions: dict[str, Any] = field(default_factory=dict)
    loads: dict[str, Any] = field(default_factory=dict)
    mesh_sweep: dict[str, Any] = field(default_factory=dict)
    acceptance: dict[str, Any] = field(default_factory=dict)
    reference_values: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)
    probes: dict[str, Any] = field(default_factory=dict)
    sweeps: dict[str, Any] = field(default_factory=dict)

    def get(self, key: str, default: Any = None) -> Any:
        """Return a top-level field from the raw JSON specification."""
        return self.raw.get(key, default)


def _resolve_data_dir(data_dir: str | Path | None) -> Path:
    if data_dir is None:
        return _DEFAULT_DATA_DIR
    return Path(data_dir)


def _normalize_case_stem(case_id: str) -> str:
    stem = case_id.strip()
    if stem.lower().endswith(".json"):
        stem = stem[:-5]
    return stem.lower()


def list_available_cases(data_dir: str | Path | None = None) -> list[str]:
    """List available verification and example case IDs in ``data_dir``.

    Parameters
    ----------
    data_dir:
        Directory containing ``*.json`` case specifications. Defaults to
        ``WHT_IGA/tests/verification/data``.

    Returns
    -------
    list[str]
        Sorted list of lowercase case stems (e.g. ``["e01", ..., "v04", ...]``).
    """
    directory = _resolve_data_dir(data_dir)
    if not directory.is_dir():
        return []
    return sorted(p.stem.lower() for p in directory.glob("*.json"))


def load_case_spec(
    case_id: str,
    data_dir: str | Path | None = None,
) -> VerificationCaseSpec:
    """Load a verification or example case specification from JSON.

    Supports both legacy flat-section JSONs (``v04``–``v12``, ``e01``–``e05``)
    and standardized schemas containing ``meta``, ``model``, ``sweeps``, and
    ``probes`` blocks.

    Parameters
    ----------
    case_id:
        Case identifier such as ``"v04"``, ``"V04"``, or ``"e03"``.
    data_dir:
        Optional directory override containing case JSON files.

    Returns
    -------
    VerificationCaseSpec
        Parsed specification dataclass.
    """
    directory = _resolve_data_dir(data_dir)
    stem = _normalize_case_stem(case_id)
    json_path = directory / f"{stem}.json"
    if not json_path.is_file():
        available = list_available_cases(directory)
        raise FileNotFoundError(
            f"Verification case '{case_id}' not found at '{json_path}'. "
            f"Available cases: {available}"
        )

    raw: dict[str, Any] = json.loads(json_path.read_text(encoding="utf-8"))
    meta = dict(raw.get("meta", {}))
    model = dict(raw.get("model", {}))

    resolved_id = str(raw.get("id") or meta.get("id") or stem.upper())
    default_cat = "example" if stem.startswith("e") else "verification"
    category = str(meta.get("category") or raw.get("category") or default_cat)

    return VerificationCaseSpec(
        case_id=resolved_id,
        title=str(raw.get("title", "")),
        category=category,
        status=str(raw.get("status", "unknown")),
        raw=raw,
        geometry=dict(raw.get("geometry") or model.get("geometry") or {}),
        materials=dict(raw.get("materials") or model.get("materials") or {}),
        boundary_conditions=dict(
            raw.get("boundary_conditions") or model.get("boundary_conditions") or {}
        ),
        loads=dict(raw.get("loads") or model.get("loads") or {}),
        mesh_sweep=dict(raw.get("mesh_sweep") or raw.get("sweeps") or {}),
        acceptance=dict(raw.get("acceptance") or {}),
        reference_values=dict(raw.get("reference_values") or {}),
        provenance=dict(raw.get("provenance") or {}),
        meta=meta,
        probes=dict(raw.get("probes") or {}),
        sweeps=dict(raw.get("sweeps") or raw.get("mesh_sweep") or {}),
    )


def _parse_csv_scalar(value: str | None) -> Any:
    """Parse a CSV string cell into ``bool``, ``int``, ``float``, or ``str``."""
    if value is None:
        return ""
    stripped = value.strip()
    if stripped == "":
        return ""
    lower = stripped.lower()
    if lower == "true":
        return True
    if lower == "false":
        return False
    if _INT_PATTERN.match(stripped):
        try:
            return int(stripped)
        except ValueError:
            return stripped
    try:
        return float(stripped)
    except ValueError:
        return stripped


def load_reference_csv(path: str | Path) -> list[dict[str, Any]]:
    """Load a Coreform or WHT_IGA reference CSV with automatic scalar typing.

    Parameters
    ----------
    path:
        Path to the CSV file.

    Returns
    -------
    list[dict[str, Any]]
        List of row dictionaries with values parsed into ``bool``, ``int``,
        ``float``, or ``str``.
    """
    csv_path = Path(path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Reference CSV not found: '{csv_path}'")

    rows: list[dict[str, Any]] = []
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return rows
        for raw_row in reader:
            parsed = {
                str(k).strip(): _parse_csv_scalar(v)
                for k, v in raw_row.items()
                if k is not None
            }
            rows.append(parsed)
    return rows
