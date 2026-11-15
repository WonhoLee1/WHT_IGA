"""S13: mixed volumetric-strain elements on cell-wise QP geometries.

Reference: Kratos source-build behavior (same build under test). The
volumetric-strain element crashes the process on QP geometries while the
sibling mixed-strain element initializes — recorded via a subprocess so the
crash cannot kill this suite.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]

_CHILD = (
    "from wht_iga.kratos_backend import load_kratos;"
    "from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume;"
    "from wht_iga.reference.box import BoxSpline;"
    "mods = load_kratos();"
    "KM = mods['KratosMultiphysics']; SM = mods['StructuralMechanicsApplication'];"
    "box = BoxSpline(lo=(0.0,0.0,0.0), hi=(1.0,1.0,1.0), n_cells=(1,1,1), degree=1);"
    "model = KM.Model(); mp = model.CreateModelPart('X');"
    "mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT);"
    "mp.AddNodalSolutionStepVariable(KM.REACTION);"
    "mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3);"
    "vol, _ = make_nurbs_volume(KM, mp, box);"
    "props = mp.GetProperties()[1];"
    "props.SetValue(KM.YOUNG_MODULUS, 1.0); props.SetValue(KM.POISSON_RATIO, 0.0);"
    "props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw());"
    "qps = KM.GeometriesVector();"
    "vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, 0, 0, 0));"
    "el = mp.CreateNewElement('SmallDisplacementMixedVolumetricStrainElement3D4N', 1, qps[0], props);"
    "print('created', flush=True);"
    "el.Initialize(mp.ProcessInfo);"
    "print('INIT-OK', flush=True)"
)


def test_mixed_volumetric_strain_crashes_on_qp_geometry():
    proc = subprocess.run(
        [sys.executable, "-c", _CHILD],
        capture_output=True,
        text=True,
        timeout=300,
        cwd="WHT_IGA",
    )
    assert "created" in proc.stdout, proc.stderr[-2000:]
    assert proc.returncode != 0, "expected the volumetric element to crash; reclassify S13"
    assert "INIT-OK" not in proc.stdout
    print("S13 volumetric returncode:", proc.returncode)


def test_mixed_strain_initializes_on_qp_geometry():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S13")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, 0, 0, 0))
    el = mp.CreateNewElement("SmallDisplacementMixedStrainElement3D4N", 1, qps[0], props)
    el.Initialize(mp.ProcessInfo)
