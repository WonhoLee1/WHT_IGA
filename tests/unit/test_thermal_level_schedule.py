"""Thermal level schedules (Kratos-free): uniform ``thermal_levels`` and the
explicit-schedule check ``check_levels`` used to replay a reference solver's
accepted increments (E03, Coreform/Abaqus convergence-history CSV)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from wht_iga.analyses import thermal_plasticity as TP

_CSV = (Path(__file__).resolve().parents[2] / "research-cf-IGA" / "user-guide" / "generated-downloads"
        / "example-problems" / "glass-to-metal-seal" / "glass-to-metal-seal-convergence-history.csv")


def test_thermal_levels_uniform() -> None:
    assert TP.thermal_levels(2, 2) == [0.5, 1.0, 1.5, 2.0]
    assert TP.thermal_levels(2, [1, 4]) == [1.0, 1.25, 1.5, 1.75, 2.0]


def test_check_levels_accepts_valid_schedule() -> None:
    assert TP.check_levels([0.01, 0.5, 1.0, 1.2, 2.0], 2) == [0.01, 0.5, 1.0, 1.2, 2.0]


@pytest.mark.parametrize("levels, msg", [
    ([], "empty"),
    ([0.0, 1.0], "strictly increasing"),
    ([0.5, 0.5, 1.0], "strictly increasing"),
    ([0.6, 0.4, 1.0], "strictly increasing"),
    ([0.5, 1.0, 1.5], "missing"),          # step end 2 not reached
    ([0.5, 1.5, 2.0], "missing"),          # step end 1 skipped
    ([1.0, 2.0, 2.5], "beyond"),
])
def test_check_levels_rejects(levels, msg) -> None:
    with pytest.raises(ValueError, match=msg):
        TP.check_levels(levels, 2)


def test_coreform_e03_history_is_a_valid_schedule() -> None:
    """The Coreform E03 CSV (accepted rows) has 15 increments per step and
    ends each step exactly at 1, 2, 3 (reference for the replay test)."""
    rows = [r for r in csv.DictReader(_CSV.open(encoding="utf-8")) if r["accepted"] == "True"]
    per_step: dict[str, int] = {}
    for r in rows:
        per_step[r["step_name"]] = per_step.get(r["step_name"], 0) + 1
    assert per_step == {"cool_to_room": 15, "cool_to_minus_55": 15, "heat_to_200": 15}
    lv = TP.check_levels([float(r["total_time"]) for r in rows], 3)
    assert len(lv) == 45 and lv[-1] == pytest.approx(3.0)
    assert min(b - a for a, b in zip([0.0] + lv, lv)) == pytest.approx(0.01)
