# -*- coding: utf-8 -*-
"""Unit tests for ASCII ACIS ``.sat`` CAD import via :mod:`wht_iga.geometry.sat_io`."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from wht_iga.geometry.cad import brep_volume, shape_bbox, tessellate, validate_shape
from wht_iga.geometry.occ_io import SAT_EXTENSIONS, load_cad_file

_SAT_PATH = Path(__file__).resolve().parents[2] / "wht_iga" / "resources" / "Door_Lock_Assembly.sat"


def test_sat_extension_registered() -> None:
    """Ensure ``.sat`` is registered in ``SAT_EXTENSIONS``."""
    assert ".sat" in SAT_EXTENSIONS


def test_load_door_lock_assembly_sat() -> None:
    """Verify ``Door_Lock_Assembly.sat`` imports as 3 valid watertight B-Rep solids."""
    from OCP.BRepCheck import BRepCheck_Analyzer

    assert _SAT_PATH.is_file(), f"Missing bundled SAT model: {_SAT_PATH}"
    cad = load_cad_file(_SAT_PATH, heal=False)
    assert cad.file_format == "sat"
    assert cad.healing_report.free_edges_after == 0
    assert cad.healing_report.watertight
    assert len(cad.solids) == 3

    vols = [brep_volume(s) for s in cad.solids]
    bboxes = [shape_bbox(s) for s in cad.solids]

    # Link (~3658 mm^3), Hinge (~14662 mm^3), Pin (~1529 mm^3)
    assert vols[0] == pytest.approx(3658.069, rel=1e-3)
    assert vols[1] == pytest.approx(14662.065, rel=1e-3)
    assert vols[2] == pytest.approx(1529.453, rel=1e-3)

    for solid, vol in zip(cad.solids, vols):
        assert BRepCheck_Analyzer(solid).IsValid()
        mesh = tessellate(solid, deflection=0.05)
        rep = validate_shape(solid, mesh, n_samples=5)
        assert rep.watertight
        assert rep.brep_volume == pytest.approx(vol, rel=1e-6)

    # Pin centered along Y at (x=-1.5, z=10.5) spanning y in [-21.15, +21.15]
    p_lo, p_hi = bboxes[2]
    assert np.isclose(0.5 * (p_lo[1] + p_hi[1]), 0.0, atol=0.1)
    assert p_hi[1] - p_lo[1] == pytest.approx(42.29, abs=0.2)
