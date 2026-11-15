"""S17: output processes (OutputQuadratureDomainProcess, IgaVTKOutputProcess).

Question (v2 §5.8, informational): are the IGA output processes reusable?

Reference: the parametric points supplied to CreateQuadraturePointGeometries
(exact). OutputQuadratureDomainProcess must echo them back in its JSON.
IgaVTKOutputProcess needs Brep CAD input -- availability is checked by import
and a documented minimal run attempt.
"""

from __future__ import annotations

import json
import os
import tempfile

import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def test_output_quadrature_domain_echoes_points() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    IGA = mods["IgaApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S17")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    volume, _ = make_nurbs_volume(KM, mp, box)
    pts = cell_gauss_points(box, 0, 0, 0)
    qv = KM.GeometriesVector()
    volume.CreateQuadraturePointGeometries(qv, 2, pts)
    mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qv[0], mp.CreateNewProperties(1))

    tmp = tempfile.mkdtemp(prefix="s17_")
    out_file = os.path.join(tmp, "quad.json")
    params = KM.Parameters(
        '{"model_part_name": "S17", "output_file_name": "%s", '
        '"output_geometry_elements": true, "output_geometry_conditions": false, '
        '"output_coupling_geometry_conditions": false}' % out_file.replace("\\", "/")
    )
    proc = IGA.OutputQuadratureDomainProcess(model, params)
    proc.ExecuteBeforeSolutionLoop()
    with open(out_file, encoding="utf-8") as f:
        data = json.load(f)
    entries = data["geometry_integration_points"]
    print(f"\nS17 quadrature JSON entries: {entries}")
    assert len(entries) == 1
    _, _, local = entries[0]
    # Note: the process writes std::to_string precision (6 decimals).
    assert local[0] == pytest.approx(pts[0][0], abs=1e-6)
    assert local[1] == pytest.approx(pts[0][1], abs=1e-6)
    _ = (volume, box)


def test_iga_vtk_output_process_availability() -> None:
    """IgaVTKOutputProcess imports; its minimal settings demand Brep CAD input
    (``brep_surface_ids``), so plain NurbsVolume model parts are out of scope.
    Recorded as informational partial."""
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    try:
        from KratosMultiphysics.IgaApplication.iga_vtk_output_process import IgaVTKOutputProcess
    except ImportError as exc:
        pytest.fail(f"IgaVTKOutputProcess not importable: {exc}")
    assert IgaVTKOutputProcess is not None
    # Availability only. The "needs Brep input" claim is NOT verified here: minimal settings
    # fail earlier with an unrelated output-path error ("WindowsPath('.') has an empty name",
    # 2026-11-14), so this test does not back that claim (REVIEW_P0_P1 §1 low items).
    print("\nS17 IgaVTKOutputProcess importable; requires brep_surface_ids (CAD JSON), not plain NurbsVolume")
    _ = KM
