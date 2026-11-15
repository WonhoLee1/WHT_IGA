"""S06: eigensolver under MPC — Kratos path vs T-projection math.

References: Euler-Bernoulli f1 (analytic, loose: coarse linear mesh),
Kratos hex-mesh eigen (tight cross-check), local consistent mass
(textbook Gauss formula, verified sum/symmetry).

Finding: Kratos EigensolverStrategy fails on cell-wise QP geometries
because SM element mass on a QuadraturePointGeometry is identically zero.
The T-projection eigen math itself is verified with scipy.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference import assemble_elastic as ref
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

L, H = 5.0, 1.0
E, NU, RHO = 1.0, 0.0, 1.0
EB_OMEGA1_SQ = ((1.87510407**2) * np.sqrt(E * (H**4 / 12) / (RHO * H * H)) / L**2) ** 2


def _cantilever_box() -> BoxSpline:
    return BoxSpline(lo=(0.0, 0.0, 0.0), hi=(L, H, H), n_cells=(4, 1, 1), degree=1)


def test_hex_cantilever_eigen_order_of_magnitude():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("Structure")
    cl = SM.LinearElastic3DLaw()

    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
        StructuralMechanicsAnalysis,
    )

    params = KM.Parameters(
        """{
        "problem_data": {"problem_name": "s06", "parallel_type": "OpenMP", "echo_level": 0, "start_time": 0.0, "end_time": 1.0},
        "solver_settings": {
            "solver_type": "eigen_value", "model_part_name": "Structure", "domain_size": 3,
            "echo_level": 0, "model_import_settings": {"input_type": "use_input_model_part"},
            "block_builder": true, "material_import_settings": {"materials_filename": ""},
            "time_stepping": {"time_step": 1.1},
            "eigensolver_settings": {"solver_type": "eigen_eigensystem", "max_iteration": 1000,
                "tolerance": 1e-9, "number_of_eigenvalues": 3, "echo_level": 0},
            "rotation_dofs": false
        },
        "processes": {"constraints_process_list": [{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "Parameters": {"model_part_name": "Structure.Dirichlet", "variable_name": "DISPLACEMENT",
                "interval": [0.0, "End"], "constrained": [true, true, true], "value": [0.0, 0.0, 0.0]
        }}], "loads_process_list": []}
    }"""
    )
    analysis = StructuralMechanicsAnalysis(model, params)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    nx = 5
    for k in range(2):
        for j in range(2):
            for i in range(nx):
                mp.CreateNewNode(1 + i + nx * (j + 2 * k), float(i) * L / (nx - 1), float(j) * H, float(k) * H)

    def nid_at(i, j, k):
        return 1 + i + nx * (j + 2 * k)

    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, NU)
    props.SetValue(KM.DENSITY, RHO)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    eid = 1
    for k in range(1):
        for j in range(1):
            for i in range(nx - 1):
                mp.CreateNewElement(
                    "SmallDisplacementElement3D8N", eid,
                    [nid_at(i, j, k), nid_at(i + 1, j, k), nid_at(i + 1, j + 1, k), nid_at(i, j + 1, k),
                     nid_at(i, j, k + 1), nid_at(i + 1, j, k + 1), nid_at(i + 1, j + 1, k + 1), nid_at(i, j + 1, k + 1)],
                    props,
                )
                eid += 1
    fixed = mp.CreateSubModelPart("Dirichlet")
    for node in mp.Nodes:
        if node.X0 == 0.0:
            fixed.AddNodes([node.Id])
    analysis.Run()
    ev = mp.ProcessInfo[SM.EIGENVALUE_VECTOR]
    vals = sorted(float(ev[i]) for i in range(ev.Size()))
    assert vals[0] > 0.0
    assert abs(vals[0] - vals[1]) / vals[0] < 1e-6, "two bending planes must coincide"
    assert vals[0] / EB_OMEGA1_SQ < 3.0 and EB_OMEGA1_SQ / vals[0] < 3.0, (vals[0], EB_OMEGA1_SQ)
    print("hex eigvals:", vals, "EB:", EB_OMEGA1_SQ)


def test_qp_mass_matrix_is_zero():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("S06M")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(1, 1, 1), degree=1)
    vol, _ = make_nurbs_volume(KM, mp, box)
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, 1.0)
    props.SetValue(KM.POISSON_RATIO, 0.0)
    props.SetValue(KM.DENSITY, 2.0)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    from wht_iga.kratos_backend.geometry import cell_gauss_points

    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, 0, 0, 0))
    el = mp.CreateNewElement("SmallDisplacementElement3D4N", 1, qps[0], props)
    el.Initialize(mp.ProcessInfo)
    M = KM.Matrix(0, 0)
    el.CalculateMassMatrix(M, mp.ProcessInfo)
    assert float(np.array(M).sum()) == 0.0, "QP consistent mass must be exactly zero (root cause)"
    lump = KM.Vector(0)
    el.CalculateLumpedMassVector(lump, mp.ProcessInfo)
    assert abs(float(np.array(lump).sum()) - 6.0) < 1e-12, "lumped mass works: 3 x rho*V"
    # control: standard hex element gives rho*V per component block
    mp2 = model.CreateModelPart("S06H")
    mp2.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp2.AddNodalSolutionStepVariable(KM.REACTION)
    mp2.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    for n, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], start=1):
        mp2.CreateNewNode(n, float(x), float(y), float(z))
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp2)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp2)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp2)
    props2 = mp2.GetProperties()[1]
    props2.SetValue(KM.YOUNG_MODULUS, 1.0)
    props2.SetValue(KM.POISSON_RATIO, 0.0)
    props2.SetValue(KM.DENSITY, 2.0)
    props2.SetValue(KM.CONSTITUTIVE_LAW, cl)
    elh = mp2.CreateNewElement("SmallDisplacementElement3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], props2)
    elh.Initialize(mp2.ProcessInfo)
    Mh = KM.Matrix(0, 0)
    elh.CalculateMassMatrix(Mh, mp2.ProcessInfo)
    assert abs(float(np.array(Mh).sum()) - 6.0) < 1e-12, "3 blocks x rho*V=2"


def test_qp_cantilever_eigen_fails():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    model = KM.Model()
    mp = model.CreateModelPart("Structure")
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
        StructuralMechanicsAnalysis,
    )

    params = KM.Parameters(
        """{
        "problem_data": {"problem_name": "s06", "parallel_type": "OpenMP", "echo_level": 0, "start_time": 0.0, "end_time": 1.0},
        "solver_settings": {
            "solver_type": "eigen_value", "model_part_name": "Structure", "domain_size": 3,
            "echo_level": 0, "model_import_settings": {"input_type": "use_input_model_part"},
            "block_builder": true, "material_import_settings": {"materials_filename": ""},
            "time_stepping": {"time_step": 1.1},
            "eigensolver_settings": {"solver_type": "eigen_eigensystem", "max_iteration": 1000,
                "tolerance": 1e-9, "number_of_eigenvalues": 3, "echo_level": 0},
            "rotation_dofs": false
        },
        "processes": {"constraints_process_list": [{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "Parameters": {"model_part_name": "Structure.Dirichlet", "variable_name": "DISPLACEMENT",
                "interval": [0.0, "End"], "constrained": [true, true, true], "value": [0.0, 0.0, 0.0]
        }}], "loads_process_list": []}
    }"""
    )
    analysis = StructuralMechanicsAnalysis(model, params)
    box = _cantilever_box()
    vol, _ = make_nurbs_volume(KM, mp, box)
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, E)
    props.SetValue(KM.POISSON_RATIO, NU)
    props.SetValue(KM.DENSITY, RHO)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    from wht_iga.kratos_backend.geometry import cell_gauss_points

    eid = 1
    for i in range(4):
        qps = KM.GeometriesVector()
        vol.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, 0, 0))
        mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
        eid += 1
    fixed = mp.CreateSubModelPart("Dirichlet")
    for node in mp.Nodes:
        if node.X0 == 0.0:
            fixed.AddNodes([node.Id])
    analysis.Run()
    ev = mp.ProcessInfo[SM.EIGENVALUE_VECTOR]
    vals = sorted(float(ev[i]) for i in range(ev.Size()))
    print("qp eigvals:", vals)
    assert not (vals[0] > 0.0 and vals[0] / EB_OMEGA1_SQ < 3.0 and EB_OMEGA1_SQ / vals[0] < 3.0), (
        "QP eigen unexpectedly sane — reclassify S06",
        vals,
    )


def _consistent_mass(box: BoxSpline, rho: float) -> sp.csr_matrix:
    """Textbook consistent mass for the tensor B-spline solid (Gauss p+1)."""
    from wht_iga.grid.bspline import basis_1d

    p = box.degree
    xg, wg = np.polynomial.legendre.leggauss(p + 1)
    tabs = [box.axis_tables(d, xg) for d in range(3)]
    (fu, Nu, _, _, _), (fv, Nv, _, _, _), (fw, Nw, _, _, _) = tabs
    ncu, ncv, ncw = Nu.shape[0], Nv.shape[0], Nw.shape[0]
    rows, cols, vals = [], [], []
    nu, nv, _ = box.n_cp
    for k in range(ncw):
        for j in range(ncv):
            for i in range(ncu):
                loc_ids = [
                    (fu[i] + a) + nu * ((fv[j] + b) + nv * (fw[k] + c))
                    for c in range(p + 1)
                    for b in range(p + 1)
                    for a in range(p + 1)
                ]
                for gz in range(p + 1):
                    for gy in range(p + 1):
                        for gx in range(p + 1):
                            wdet = tabs[0][4][i] * tabs[1][4][j] * tabs[2][4][k] * wg[gx] * wg[gy] * wg[gz]
                            nvec = np.array([Nu[i, gx, a] * Nv[j, gy, b] * Nw[k, gz, c]
                                             for c in range(p + 1) for b in range(p + 1) for a in range(p + 1)])
                            for r, I in enumerate(loc_ids):
                                for s, J in enumerate(loc_ids):
                                    for d in range(3):
                                        rows.append(3 * I + d)
                                        cols.append(3 * J + d)
                                        vals.append(rho * wdet * nvec[r] * nvec[s])
    M = sp.coo_matrix((vals, (rows, cols)), shape=(box.n_dofs, box.n_dofs)).tocsr()
    assert abs(M.sum() - 3.0 * rho * L * H * H) < 1e-9 * rho
    assert abs((M - M.T).max()) < 1e-15
    return M


def test_reduced_eigen_matches_beam_theory_and_tie():
    box = _cantilever_box()
    K = ref.assemble_stiffness(box, E, NU).tocsr()
    M = _consistent_mass(box, RHO)
    nu_ = box.n_cp[0]
    fixed = sorted([3 * I + d for I in range(box.n_ctrl) if I % nu_ == 0 for d in range(3)])
    free = [d for d in range(box.n_dofs) if d not in set(fixed)]
    w, _ = eigsh(K[free][:, free], M=M[free][:, free], k=3, which="SM")
    w = np.sort(w)
    assert w[0] > 0.0
    assert abs(w[0] - w[1]) / w[0] < 1e-6
    assert w[0] / EB_OMEGA1_SQ < 3.0 and EB_OMEGA1_SQ / w[0] < 3.0, (w[0], EB_OMEGA1_SQ)
    print("scipy eigvals:", w.tolist(), "EB:", EB_OMEGA1_SQ)
    # T-projection: duplicate tip-center node, tie it 1:1 (extension-style MPC)
    tip = int(nu_ - 1 + nu_ * (1 + 2 * 1))
    n = box.n_dofs
    T = sp.lil_matrix((n + 3, n))
    T[:n, :] = sp.eye(n)
    for d in range(3):
        T[n + d, 3 * tip + d] = 1.0
    T = T.tocsr()
    Ke = sp.block_diag((K, sp.csr_matrix((3, 3)))).tocsr()
    Me = sp.block_diag((M, sp.csr_matrix((3, 3)))).tocsr()
    Kred = (T.T @ Ke @ T).tocsc()
    Mred = (T.T @ Me @ T).tocsc()
    free_e = [d for d in range(n) if d not in set(fixed)]
    w2, _ = eigsh(Kred[free_e][:, free_e], M=Mred[free_e][:, free_e], k=3, which="SM")
    w2 = np.sort(w2)
    np.testing.assert_allclose(w2, w, rtol=1e-8)
