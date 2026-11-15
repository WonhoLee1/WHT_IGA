"""S16: model-build scaling — ~1e5 cell-wise QP geometries + SM elements.

Reference: exact cell count (46^3) and node count (48^3); time/memory are
recorded observations against the v2 targets (< 60 s, < 8 GB).
"""

from __future__ import annotations

import time

import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos, pytest.mark.slow]


def test_1e5_qp_elements_creation_time_and_memory():
    import psutil

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    n = 46
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, n), degree=2)
    model = KM.Model()
    mp = model.CreateModelPart("S16")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box)
    assert mp.NumberOfNodes() == (n + 2) ** 3
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    proc = psutil.Process()
    rss0 = proc.memory_info().rss / 1e9
    t0 = time.perf_counter()
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, j, k))
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                eid += 1
    dt = time.perf_counter() - t0
    rss1 = proc.memory_info().rss / 1e9
    n_el = mp.NumberOfElements()
    print(f"elements={n_el} nodes={mp.NumberOfNodes()} time={dt:.1f}s rss_delta={rss1 - rss0:.2f}GB")
    assert n_el == n**3 == 97336
    assert dt < 60.0, dt
    assert rss1 - rss0 < 8.0, rss1 - rss0


def test_p3_qp_elements_small_scale_rate():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    n = 10
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n, n, n), degree=3)
    model = KM.Model()
    mp = model.CreateModelPart("S16p3")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    t0 = time.perf_counter()
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, j, k))
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                eid += 1
    dt = time.perf_counter() - t0
    print(f"p=3: elements={mp.NumberOfElements()} time={dt:.2f}s rate={mp.NumberOfElements() / dt:.0f}/s")
    assert mp.NumberOfElements() == n**3
