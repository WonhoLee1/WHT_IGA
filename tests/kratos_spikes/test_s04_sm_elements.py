"""S04: StructuralMechanics elements on cell-wise QP geometries.

References: oracle ``assemble_stiffness`` (K equality), analytic linear field
u = (x, 2y, 3z) with E=1/nu=0 (patch displacement, Cauchy stress, von Mises).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    solve_linear,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference import assemble_elastic as ref
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _model(KM, name: str = "S04"):
    model = KM.Model()
    mp = model.CreateModelPart(name)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    return model, mp


def _props(KM, SM, mp, E: float = 1.0, nu: float = 0.0):
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    return props


def _fix_patch(KM, mp, box: BoxSpline, faces: str = "all"):
    ids = set()
    for a in range(3):
        for s in (0, 1):
            ids.update(int(v) for v in box.boundary_control_points(a, s))
    for node in mp.Nodes:
        if node.Id - 1 in ids:
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, node.X0)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 2.0 * node.Y0)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 3.0 * node.Z0)


def _nodal_u(KM, mp, n: int) -> np.ndarray:
    u = np.zeros(3 * n)
    for node in mp.Nodes:
        I = node.Id - 1
        u[3 * I] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
        u[3 * I + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
        u[3 * I + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)
    return u


@pytest.mark.parametrize("degree,n_cells", [(1, (1, 1, 1)), (1, (2, 2, 2)), (2, (2, 2, 2)), (3, (1, 1, 1)), (3, (2, 2, 2))])
def test_global_stiffness_matches_reference(degree, n_cells):
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=n_cells, degree=degree)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp))
    strat = solve_linear(KM, mp)
    A = system_matrix_csr(strat).toarray()
    perm = equation_permutation(KM, mp, box)
    Kref = ref.assemble_stiffness(box, 1.0, 0.0).toarray()
    Kk = Kref[np.ix_(np.argsort(perm), np.argsort(perm))]
    rel = np.linalg.norm(A - Kk) / np.linalg.norm(Kk)
    print(f"S04 K rel err p={degree} cells={n_cells}:", rel)
    assert rel < 1e-12, rel


def test_patch_displacement():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2), degree=2)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp))
    _fix_patch(KM, mp, box)
    solve_linear(KM, mp)
    u = _nodal_u(KM, mp, box.n_ctrl)
    X = box.control_points()
    exact = np.column_stack([X[:, 0], 2 * X[:, 1], 3 * X[:, 2]]).ravel()
    print("S04 patch max err:", np.max(np.abs(u - exact)))
    assert np.max(np.abs(u - exact)) < 1e-10, np.max(np.abs(u - exact))


def test_cauchy_stress_and_von_mises_on_patch():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp))
    _fix_patch(KM, mp, box)
    solve_linear(KM, mp)
    sig = np.array(els[0].CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo))
    assert sig.shape == (8, 6), sig.shape
    np.testing.assert_allclose(sig, np.tile([1.0, 2.0, 3.0, 0.0, 0.0, 0.0], (8, 1)), atol=1e-12)
    vm = np.array(els[0].CalculateOnIntegrationPoints(SM.VON_MISES_STRESS, mp.ProcessInfo)).ravel()
    np.testing.assert_allclose(vm, np.full(8, np.sqrt(3.0)), rtol=1e-12)
    print("S04 von Mises:", vm[0], "expected:", np.sqrt(3.0))


@pytest.mark.parametrize("element_name", ["TotalLagrangianElement3D4N", "UpdatedLagrangianElement3D4N"])
def test_nonlinear_elements_initialize(element_name):
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=2)
    _, mp = _model(KM)
    vol, _ = make_nurbs_volume(KM, mp, box)
    els = attach_cell_elements(KM, mp, vol, box, _props(KM, SM, mp), element_name=element_name)
    for el in els:
        el.Initialize(mp.ProcessInfo)
