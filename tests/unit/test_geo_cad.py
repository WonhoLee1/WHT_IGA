"""S1 CAD tests: Roark plate (analytic refs), STEP roundtrip, fuse (v2 §7.1)."""

from __future__ import annotations

import importlib.util
import os

import numpy as np

from wht_iga.geometry import cad

DATA = os.path.join(os.path.dirname(__file__), "..", "..", "data")
STEP = os.path.join(DATA, "roark_plate_hole.step")

L, D, R, T = 200.0, 40.0, 10.0, 1.0
ANA_VOL = L * D * T - np.pi * R * R * T
ANA_AREA = 2 * (L * D - np.pi * R * R) + 2 * T * (L + D) + 2 * np.pi * R * T


def _load_gen():
    spec = importlib.util.spec_from_file_location(
        "gen_roark_plate", os.path.join(DATA, "gen_roark_plate.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_roark_brep_matches_analytic() -> None:
    """BRep volume/area vs closed-form box-minus-cylinder values."""
    shape = cad.read_step(STEP)
    np.testing.assert_allclose(cad.brep_volume(shape), ANA_VOL, rtol=1e-12)
    np.testing.assert_allclose(cad.brep_area(shape), ANA_AREA, rtol=1e-9)
    lo, hi = cad.shape_bbox(shape)
    np.testing.assert_allclose(lo, [-L / 2, -D / 2, -T / 2], atol=1e-6)
    np.testing.assert_allclose(hi, [L / 2, D / 2, T / 2], atol=1e-6)


def test_roark_tessellation_validates() -> None:
    """Tessellation: 7 faces, watertight, validate() green, bounded error."""
    shape = cad.read_step(STEP)
    tess = cad.tessellate(shape, deflection=0.25)
    assert len(set(tess.face_ids.tolist())) == 7
    assert cad.boundary_edge_count(tess) == 0
    rep = cad.validate_shape(shape, tess)
    assert rep.watertight
    assert rep.tess_volume > 0
    np.testing.assert_allclose(rep.tess_volume, rep.brep_volume,
                               rtol=10 * tess.deflection / cad.bbox_diag(shape))
    assert rep.area_err < 5e-4
    assert rep.max_dev < tess.deflection
    rep.issues.validate()  # no error-severity issues


def test_roark_generator_reproducible(tmp_path) -> None:
    """Rerunning the generator reproduces the committed STEP geometry."""
    gen = _load_gen()
    out = str(tmp_path / "roark_repro.step")
    gen.build_plate()
    shape0 = cad.read_step(STEP)
    cad.write_step(gen.build_plate(), out)
    shape1 = cad.read_step(out)
    t0 = cad.tessellate(shape0, deflection=0.25)
    t1 = cad.tessellate(shape1, deflection=0.25)
    assert len(t0.tris) == len(t1.tris) and len(t0.nodes) == len(t1.nodes)
    np.testing.assert_allclose(cad.tess_area(t0), cad.tess_area(t1), rtol=0)


def test_step_roundtrip_preserves_geometry(tmp_path) -> None:
    """Write→read keeps bbox and area."""
    shape = cad.read_step(STEP)
    out = str(tmp_path / "rt.step")
    cad.write_step(shape, out)
    back = cad.read_step(out)
    lo0, hi0 = cad.shape_bbox(shape)
    lo1, hi1 = cad.shape_bbox(back)
    np.testing.assert_allclose(lo0, lo1, rtol=0)
    np.testing.assert_allclose(hi0, hi1, rtol=0)
    np.testing.assert_allclose(cad.brep_area(shape), cad.brep_area(back),
                               rtol=1e-12)


def test_fuse_cells_general_fuse() -> None:
    """BOPAlgo_Builder general fuse of two touching boxes: V = 2000."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    a = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 10.0, 10.0, 10.0).Shape()
    b = BRepPrimAPI_MakeBox(gp_Pnt(10, 0, 0), 10.0, 10.0, 10.0).Shape()
    fused = cad.fuse_cells([cad.heal(a), cad.heal(b)])
    np.testing.assert_allclose(cad.brep_volume(fused), 2000.0, rtol=1e-9)
    tess = cad.tessellate(fused, deflection=0.5)
    rep = cad.validate_shape(fused, tess)
    assert rep.watertight
    rep.issues.validate()
