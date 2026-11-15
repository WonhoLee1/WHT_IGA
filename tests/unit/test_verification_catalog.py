"""Unit tests for the 17-case WHT_IGA Verification & Example catalog."""

from __future__ import annotations

import pytest

from wht_iga.verification.catalog import (
    CASE_CATALOG,
    CaseCatalogEntry,
    ProbeSpec,
    SweepSpec,
    get_catalog_entry,
    list_cases_by_phase,
    resolve_reference_csvs,
)
from wht_iga.verification.schema import (
    VerificationCaseSpec,
    load_case_spec,
    load_reference_csv,
)

EXPECTED_CASE_IDS: tuple[str, ...] = (
    *(f"V{i:02d}" for i in range(1, 13)),
    *(f"E{i:02d}" for i in range(1, 6)),
)

EXPECTED_PHASE_GATES: dict[str, list[str]] = {
    "G1": ["V01", "V02", "V03"],
    "G2": ["V06"],
    "G3": ["V04", "V05", "E01", "E02"],
    "G4": ["V09", "V10", "V11", "V12"],
    "G5": ["V07", "V08", "E03", "E04", "E05"],
}

VALID_EXTRACTORS: set[str] = {
    "FaceAndGridProbeExtractor",
    "LineProfileExtractor",
    "VolumeQuadratureL2Extractor",
    "ContactHistoryExtractor",
    "ModalAndSpectrumExtractor",
    "DynamicHistoryExtractor",
    "NonlinearThermoExtractor",
}


def test_catalog_contains_all_17_cases() -> None:
    """Verify V01-V12 and E01-E05 are all registered in CASE_CATALOG."""
    assert tuple(CASE_CATALOG.keys()) == EXPECTED_CASE_IDS
    assert len(CASE_CATALOG) == 17


@pytest.mark.parametrize("case_id", EXPECTED_CASE_IDS)
def test_catalog_entries_load_case_spec_cleanly(case_id: str) -> None:
    """Verify every catalog entry matches and loads its JSON spec."""
    entry = get_catalog_entry(case_id)
    assert isinstance(entry, CaseCatalogEntry)
    assert entry.case_id == case_id

    spec = load_case_spec(case_id)
    spec_via_entry = entry.load_spec()
    assert isinstance(spec, VerificationCaseSpec)
    assert spec.case_id == case_id
    assert spec_via_entry.case_id == case_id
    assert spec.category == entry.category
    assert entry.unit_system.strip() != ""
    assert entry.qoi_extractor in VALID_EXTRACTORS
    assert len(entry.probes) >= 1
    assert isinstance(entry.sweeps, SweepSpec)
    assert len(entry.coreform_output_csvs) >= 1


def test_get_catalog_entry_case_insensitive_and_error() -> None:
    """Verify get_catalog_entry normalizes case and rejects unknown IDs."""
    assert get_catalog_entry("v01").case_id == "V01"
    assert get_catalog_entry(" e03 ").case_id == "E03"
    with pytest.raises(KeyError):
        get_catalog_entry("V99")


@pytest.mark.parametrize("gate,expected_ids", EXPECTED_PHASE_GATES.items())
def test_list_cases_by_phase(gate: str, expected_ids: list[str]) -> None:
    """Verify list_cases_by_phase returns exact expected subsets for G1-G5."""
    entries = list_cases_by_phase(gate)
    assert [e.case_id for e in entries] == expected_ids
    assert [e.case_id for e in list_cases_by_phase(gate.lower())] == expected_ids


def test_list_cases_by_phase_invalid_gate() -> None:
    """Verify list_cases_by_phase rejects unknown phase gates."""
    with pytest.raises(ValueError):
        list_cases_by_phase("G9")


def test_probe_specs_have_required_geometry_fields() -> None:
    """Verify ProbeSpec coordinates and line endpoints are populated."""
    for entry in CASE_CATALOG.values():
        for probe in entry.probes:
            assert isinstance(probe, ProbeSpec)
            assert len(probe.fields) >= 1
            if probe.kind == "point":
                assert probe.coordinates is not None
                assert len(probe.coordinates) in (2, 3)
            elif probe.kind == "line" and probe.start is not None:
                assert probe.end is not None
                assert len(probe.start) == 3
                assert len(probe.end) == 3
                assert probe.n_points is not None and probe.n_points >= 2


@pytest.mark.parametrize("case_id", EXPECTED_CASE_IDS)
def test_declared_reference_csvs_exist_and_load(case_id: str) -> None:
    """Verify all declared reference CSVs exist on disk and parse cleanly."""
    resolved = resolve_reference_csvs(case_id)
    entry = get_catalog_entry(case_id)
    assert len(resolved) == len(entry.reference_csv_relpaths)
    for csv_path in resolved:
        assert csv_path.is_file(), f"Missing reference CSV for {case_id}: {csv_path}"
        rows = load_reference_csv(csv_path)
        assert isinstance(rows, list)
        assert len(rows) > 0, f"Empty reference CSV for {case_id}: {csv_path}"


def test_reference_csv_cases_are_populated() -> None:
    """Verify cases with Coreform CSV downloads have non-empty reference_csv_relpaths."""
    for cid in ("V02", "V03", "V04", "V05", "V06", "E03"):
        assert len(resolve_reference_csvs(cid)) > 0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
