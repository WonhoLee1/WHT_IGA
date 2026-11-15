"""S20: Temperature-dependent E(T) and alpha(T) via TableAccessor vs Properties cloning.

Tests:
1. ``TableAccessor`` is not directly exposed in Python ``KratosMultiphysics``,
   but can be attached via ``ReadMaterialsUtility`` JSON ``"accessors"``.
2. ``SM.LinearElastic3DLaw`` bypasses ``Accessor`` (uses ``rProps[YOUNG_MODULUS]``).
3. ``CL.ThermalElasticIsotropic3D`` honors ``TableAccessor`` for ``YOUNG_MODULUS``
   (via ``rProps.GetValue(YOUNG_MODULUS, geom, N, info)``) and matches the
   analytical clamped thermal-stress solution to < 1e-8 under p=4 refinement.
4. ``CL.ThermalElasticIsotropic3D`` free thermal expansion under linear T(x)
   matches u_x(L) = 0.5 * alpha * dT * L to < 1e-15.
5. ``CL.ThermalElasticIsotropic3D`` bypasses ``Accessor`` for
   ``THERMAL_EXPANSION_COEFFICIENT`` (``SubstractThermalStrain`` reads
   ``rProps[THERMAL_EXPANSION_COEFFICIENT]`` directly).
6. Cell-wise ``Properties`` cloning works for both ``LinearElastic3DLaw`` and
   ``ThermalElasticIsotropic3D`` (both E(T) and alpha(T)) with O(h^2) convergence.
"""

from __future__ import annotations

import json
import math
import os
import tempfile

import numpy as np
import pytest

from wht_iga.analyses.thermal import (
    assemble_stiffness_tdep,
    recover_thermoelastic_stress_tdep,
    solve_thermoelastic_tdep,
    temperature_dependent_eigenstrain_load,
)
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import attach_cell_elements, solve_linear
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


def _load_apps():
    return load_kratos(
        (
            "IgaApplication",
            "StructuralMechanicsApplication",
            "LinearSolversApplication",
            "ConstitutiveLawsApplication",
        )
    )


def _setup_bar_model(KM, n_x: int, degree: int):
    model = KM.Model()
    mp = model.CreateModelPart("Bar")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(n_x, 1, 1), degree=degree)
    vol, _ = make_nurbs_volume(KM, mp, box)
    return model, mp, box, vol


def _read_material_json(KM, model, mat_dict: dict):
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump({"properties": [{"model_part_name": "Bar", "properties_id": 1, "Material": mat_dict}]}, f)
        fname = f.name
    try:
        KM.ReadMaterialsUtility(KM.Parameters(json.dumps({"Parameters": {"materials_filename": fname}})), model)
    finally:
        os.remove(fname)


def test_table_accessor_attached_via_read_materials():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    assert not hasattr(KM, "TableAccessor")

    model, mp, _, _ = _setup_bar_model(KM, n_x=2, degree=1)
    _read_material_json(
        KM,
        model,
        {
            "name": "ThermoElastic",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": 100.0,
                "POISSON_RATIO": 0.0,
                "THERMAL_EXPANSION_COEFFICIENT": 1e-5,
                "REFERENCE_TEMPERATURE": 0.0,
            },
            "Tables": {
                "E_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "YOUNG_MODULUS",
                    "data": [[0.0, 100.0], [100.0, 200.0]],
                }
            },
            "accessors": {
                "E_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "YOUNG_MODULUS",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    props = mp.GetProperties()[1]
    assert props.HasAccessor(KM.YOUNG_MODULUS)
    assert props.HasTable(KM.TEMPERATURE, KM.YOUNG_MODULUS)


def test_linear_elastic_3d_ignores_table_accessor():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    model, mp, box, vol = _setup_bar_model(KM, n_x=4, degree=2)

    _read_material_json(
        KM,
        model,
        {
            "name": "LinearWithTable",
            "constitutive_law": {"name": "LinearElastic3DLaw"},
            "Variables": {"YOUNG_MODULUS": 100.0, "POISSON_RATIO": 0.0},
            "Tables": {
                "E_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "YOUNG_MODULUS",
                    "data": [[0.0, 100.0], [100.0, 200.0]],
                }
            },
            "accessors": {
                "E_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "YOUNG_MODULUS",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    els = attach_cell_elements(KM, mp, vol, box, mp.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    u_L = 0.01
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, 100.0 * node.X0)
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        elif abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, u_L)
    solve_linear(KM, mp)

    sigs = np.vstack([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)) for el in els])
    sig_const = 100.0 * u_L
    err_vs_const = float(np.max(np.abs(sigs[:, 0] - sig_const)))
    print(f"S20 LinearElastic3DLaw with TableAccessor: sig_xx={sigs[0, 0]:.12f} (const E0={sig_const}), err={err_vs_const:.3e}")
    assert err_vs_const < 1e-12


def test_thermal_elastic_free_expansion_exact():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    model, mp, box, vol = _setup_bar_model(KM, n_x=4, degree=2)

    alpha0, dT = 1e-5, 100.0
    _read_material_json(
        KM,
        model,
        {
            "name": "ThermoElastic",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": 100.0,
                "POISSON_RATIO": 0.0,
                "THERMAL_EXPANSION_COEFFICIENT": alpha0,
                "REFERENCE_TEMPERATURE": 0.0,
            },
        },
    )
    els = attach_cell_elements(KM, mp, vol, box, mp.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, dT * node.X0)
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
    solve_linear(KM, mp)

    end_u = np.array([n.GetSolutionStepValue(KM.DISPLACEMENT_X) for n in mp.Nodes if abs(n.X0 - 1.0) < 1e-12])
    exact_u = 0.5 * alpha0 * dT
    err = float(np.max(np.abs(end_u - exact_u)))
    sigs = np.vstack([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)) for el in els])
    max_sig = float(np.max(np.abs(sigs[:, 0])))
    print(f"S20 ThermalElastic free expansion: end_u={end_u[0]:.16e}, exact={exact_u:.16e}, err={err:.3e}, max_sig_xx={max_sig:.3e}")
    assert err < 1e-15
    assert max_sig < 1e-13


def test_thermal_elastic_table_accessor_E_clamped_bar_analytic():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    E0, E1, dT, alpha0 = 100.0, 1.0, 100.0, 1e-5
    sig_exact = -(alpha0 * E1 * dT * dT) / (2.0 * math.log(1.0 + E1 * dT / E0))

    model, mp, box, vol = _setup_bar_model(KM, n_x=24, degree=4)
    _read_material_json(
        KM,
        model,
        {
            "name": "ThermoElastic",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": E0,
                "POISSON_RATIO": 0.0,
                "THERMAL_EXPANSION_COEFFICIENT": alpha0,
                "REFERENCE_TEMPERATURE": 0.0,
            },
            "Tables": {
                "E_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "YOUNG_MODULUS",
                    "data": [[0.0, E0], [dT, E0 + E1 * dT]],
                }
            },
            "accessors": {
                "E_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "YOUNG_MODULUS",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    els = attach_cell_elements(KM, mp, vol, box, mp.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, dT * node.X0)
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12 or abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
    solve_linear(KM, mp)

    sigs = np.array([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo))[:, 0] for el in els])
    sig_mean = float(np.mean(sigs))
    rel_err = abs(sig_mean - sig_exact) / abs(sig_exact)
    print(f"S20 ThermalElastic TableAccessor E(T) p=4 n_x=24: sig={sig_mean:.12f}, exact={sig_exact:.12f}, rel_err={rel_err:.3e}")
    assert rel_err < 1e-8


def test_thermal_elastic_ignores_table_accessor_for_alpha():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    alpha0, dT = 1e-5, 100.0
    model, mp, box, vol = _setup_bar_model(KM, n_x=4, degree=2)
    _read_material_json(
        KM,
        model,
        {
            "name": "ThermoElasticAlphaTable",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": 100.0,
                "POISSON_RATIO": 0.0,
                "THERMAL_EXPANSION_COEFFICIENT": alpha0,
                "REFERENCE_TEMPERATURE": 0.0,
            },
            "Tables": {
                "alpha_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "THERMAL_EXPANSION_COEFFICIENT",
                    "data": [[0.0, alpha0], [dT, 3.0 * alpha0]],
                }
            },
            "accessors": {
                "alpha_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "THERMAL_EXPANSION_COEFFICIENT",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    attach_cell_elements(KM, mp, vol, box, mp.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, dT * node.X0)
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
    solve_linear(KM, mp)

    end_u = float(np.mean([n.GetSolutionStepValue(KM.DISPLACEMENT_X) for n in mp.Nodes if abs(n.X0 - 1.0) < 1e-12]))
    u_const_alpha = 0.5 * alpha0 * dT
    err_vs_const = abs(end_u - u_const_alpha)
    print(f"S20 ThermalElastic TableAccessor alpha(T): end_u={end_u:.12e}, const_alpha_u={u_const_alpha:.12e}, diff={err_vs_const:.3e}")
    assert err_vs_const < 1e-15


def test_cell_properties_cloning_handles_both_E_and_alpha():
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]

    E0, E1, alpha0, alpha1, dT = 100.0, 1.0, 1e-5, 2e-7, 100.0
    # Clamped bar with E(x) = E0 + E1*dT*x and alpha(x) = alpha0 + alpha1*dT*x:
    # int_0^1 alpha(x)*T(x) dx = dT*(alpha0/2 + alpha1*dT/3)
    # int_0^1 1/E(x) dx = ln(1 + E1*dT/E0) / (E1*dT)
    int_eps_th = dT * (0.5 * alpha0 + (1.0 / 3.0) * alpha1 * dT)
    int_inv_E = math.log(1.0 + E1 * dT / E0) / (E1 * dT)
    sig_exact = -int_eps_th / int_inv_E

    errs = []
    for n_x in (8, 16, 32):
        model, mp, box, vol = _setup_bar_model(KM, n_x=n_x, degree=2)
        els = []
        for i in range(n_x):
            xc = (i + 0.5) / n_x
            Tc = dT * xc
            props_i = mp.GetProperties()[i + 1]
            props_i.SetValue(KM.YOUNG_MODULUS, E0 + E1 * Tc)
            props_i.SetValue(KM.POISSON_RATIO, 0.0)
            props_i.SetValue(KM.THERMAL_EXPANSION_COEFFICIENT, alpha0 + alpha1 * Tc)
            props_i.SetValue(KM.REFERENCE_TEMPERATURE, 0.0)
            props_i.SetValue(KM.CONSTITUTIVE_LAW, CL.ThermalElasticIsotropic3D())
            els.extend(
                attach_cell_elements(
                    KM,
                    mp,
                    vol,
                    box,
                    props_i,
                    element_name="SmallDisplacementElement3D4N",
                    start_id=i + 1,
                    cells=[(i, 0, 0)],
                )
            )
        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, dT * node.X0)
            node.Fix(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
            node.Fix(KM.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
            if abs(node.X0) < 1e-12 or abs(node.X0 - 1.0) < 1e-12:
                node.Fix(KM.DISPLACEMENT_X)
                node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        solve_linear(KM, mp)
        sigs = np.array([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo))[:, 0] for el in els])
        sig_mean = float(np.mean(sigs))
        rel_err = abs(sig_mean - sig_exact) / abs(sig_exact)
        errs.append(rel_err)
        print(f"S20 cell-wise Properties clone n_x={n_x}: sig={sig_mean:.12f}, exact={sig_exact:.12f}, rel_err={rel_err:.3e}")

    # Midpoint-per-cell approximation converges at O(h^2): ratio ~ 4 per halving
    assert errs[0] / errs[1] > 3.5
    assert errs[1] / errs[2] > 3.5
    assert errs[2] < 5e-4


def _quad_x_bspline_coeffs(box: BoxSpline, c0: float, c1: float, c2: float) -> np.ndarray:
    """Exact B-spline control-point values for T(x) = c0 + c1*x + c2*x^2 (degree >= 2)."""
    p = box.degree
    kx = box.knots[0]
    cx = np.empty(box.n_cp[0], dtype=float)
    for i in range(box.n_cp[0]):
        t = kx[i + 1 : i + p + 1]
        xi1 = float(t.mean())
        xi2 = (2.0 / (p * (p - 1))) * sum(float(t[r] * t[s]) for r in range(p) for s in range(r + 1, p))
        cx[i] = c0 + c1 * xi1 + c2 * xi2
    return np.tile(cx, box.n_cp[1] * box.n_cp[2])


def test_s20_quadratic_T_table_accessor_E_kratos_python_parity():
    """Quadratic T(x) + affine E(T) via TableAccessor matches Python Gauss tdep K, f, u, sigma and analytic."""
    import scipy.sparse as sp

    mods = _load_apps()
    KM = mods["KratosMultiphysics"]

    # Part 1: 3D (nu=0.25) element-by-element K, f, solved u and Gauss-point stress parity
    E0, E1, nu, alpha0, T_ref, dT = 100.0, 2.0, 0.25, 1.0e-5, 20.0, 100.0
    model, mp, box, vol = _setup_bar_model(KM, n_x=6, degree=2)
    _read_material_json(
        KM,
        model,
        {
            "name": "ThermoElasticQuadT",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": E0,
                "POISSON_RATIO": nu,
                "THERMAL_EXPANSION_COEFFICIENT": alpha0,
                "REFERENCE_TEMPERATURE": T_ref,
            },
            "Tables": {
                "E_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "YOUNG_MODULUS",
                    "data": [[0.0, E0], [200.0, E0 + E1 * 200.0]],
                }
            },
            "accessors": {
                "E_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "YOUNG_MODULUS",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    els = attach_cell_elements(KM, mp, vol, box, mp.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    # T(x) = T_ref + 20*x + 80*x^2 in [20, 120]
    t_cp = _quad_x_bspline_coeffs(box, T_ref, 20.0, 80.0)
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, float(t_cp[node.Id - 1]))
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12 or abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)

    # Assemble Kratos local LHS and thermal RHS at u = 0 before solving
    rows, cols, k_vals = [], [], []
    f_kr = np.zeros(box.n_dofs, dtype=float)
    for el in els:
        el.Initialize(mp.ProcessInfo)
        lhs, rhs = KM.Matrix(0, 0), KM.Vector(0)
        el.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        ke = np.asarray(lhs, dtype=float)
        fe = np.asarray(rhs, dtype=float)
        ids = np.array([g.Id - 1 for g in el.GetGeometry()], dtype=np.int64)
        dofs = (3 * ids[:, None] + np.arange(3)[None, :]).reshape(-1)
        rows.append(np.repeat(dofs, dofs.size))
        cols.append(np.tile(dofs, dofs.size))
        k_vals.append(ke.reshape(-1))
        f_kr[dofs] += fe
    K_kr = sp.coo_matrix((np.concatenate(k_vals), (np.concatenate(rows), np.concatenate(cols))), shape=(box.n_dofs, box.n_dofs)).tocsr()

    E_fn = lambda T: E0 + E1 * T  # noqa: E731
    K_py = assemble_stiffness_tdep(box, t_cp, E_fn, nu)
    f_py = temperature_dependent_eigenstrain_load(box, t_cp, E_fn, nu, alpha0, T_ref=T_ref)
    assert abs(K_kr - K_py).max() < 1e-12 * abs(K_py).max()
    assert np.max(np.abs(f_kr - f_py)) < 1e-12 * np.max(np.abs(f_py))

    solve_linear(KM, mp)
    u_kr = np.zeros(box.n_dofs, dtype=float)
    for node in sorted(mp.Nodes, key=lambda n: n.Id):
        idx = node.Id - 1
        u_kr[3 * idx + 0] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
        u_kr[3 * idx + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
        u_kr[3 * idx + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)

    all_cps = np.arange(box.n_ctrl, dtype=np.int64)
    ends = np.concatenate([box.boundary_control_points(0, 0), box.boundary_control_points(0, 1)])
    fixed = np.unique(np.concatenate([3 * ends + 0, 3 * all_cps + 1, 3 * all_cps + 2]))
    res_py = solve_thermoelastic_tdep(box, t_cp, E_fn, nu, alpha0, fixed, T_ref=T_ref)
    assert np.max(np.abs(u_kr - res_py.u)) < 1e-12 * np.max(np.abs(res_py.u))

    gp_pts = np.vstack([np.asarray(cell_gauss_points(box, i, 0, 0))[:, :3] for i in range(6)])
    sig_kr = np.vstack([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo)) for el in els])
    sig_py = recover_thermoelastic_stress_tdep(box, res_py.u, t_cp, E_fn, nu, alpha0, gp_pts, T_ref=T_ref)
    assert np.max(np.abs(sig_kr - sig_py)) < 1e-11 * np.max(np.abs(sig_py))

    # Part 2: p=4 clamped bar (nu=0) vs closed-form arctan analytic solution for T(x) = dT * x^2
    E1_bar = 0.5
    model2, mp2, box2, vol2 = _setup_bar_model(KM, n_x=24, degree=4)
    _read_material_json(
        KM,
        model2,
        {
            "name": "ThermoElasticQuadTAnalytic",
            "constitutive_law": {"name": "ThermalElasticIsotropic3D"},
            "Variables": {
                "YOUNG_MODULUS": E0,
                "POISSON_RATIO": 0.0,
                "THERMAL_EXPANSION_COEFFICIENT": alpha0,
                "REFERENCE_TEMPERATURE": 0.0,
            },
            "Tables": {
                "E_vs_T": {
                    "input_variable": "TEMPERATURE",
                    "output_variable": "YOUNG_MODULUS",
                    "data": [[0.0, E0], [dT, E0 + E1_bar * dT]],
                }
            },
            "accessors": {
                "E_acc": {
                    "accessor_type": "table_accessor",
                    "properties": {
                        "table_input_variable": "TEMPERATURE",
                        "table_output_variable": "YOUNG_MODULUS",
                        "table_input_variable_type": "node_historical",
                    },
                }
            },
        },
    )
    els2 = attach_cell_elements(KM, mp2, vol2, box2, mp2.GetProperties()[1], element_name="SmallDisplacementElement3D4N")
    t_cp2 = _quad_x_bspline_coeffs(box2, 0.0, 0.0, dT)
    for node in mp2.Nodes:
        node.SetSolutionStepValue(KM.TEMPERATURE, float(t_cp2[node.Id - 1]))
        node.Fix(KM.DISPLACEMENT_Y)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        node.Fix(KM.DISPLACEMENT_Z)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
        if abs(node.X0) < 1e-12 or abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.DISPLACEMENT_X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
    solve_linear(KM, mp2)

    beta = E1_bar * dT / E0
    sig_exact = -(E0 * math.sqrt(beta) / math.atan(math.sqrt(beta))) * (alpha0 * dT / 3.0)
    gp2 = np.vstack([np.asarray(cell_gauss_points(box2, i, 0, 0)) for i in range(24)])
    sigs2 = np.concatenate([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp2.ProcessInfo))[:, 0] for el in els2])
    sig_mean = float(np.sum(gp2[:, 3] * sigs2) / np.sum(gp2[:, 3]))
    rel_err = abs(sig_mean - sig_exact) / abs(sig_exact)
    print(f"S20 Quad T(x) + TableAccessor E(T) p=4 n_x=24: sig={sig_mean:.12f}, exact={sig_exact:.12f}, rel_err={rel_err:.3e}")
    assert rel_err < 1e-10


def test_s20_dual_tdep_E_and_alpha_kratos_properties_vs_python_gauss():
    """Quadratic T(x) + affine E(T), alpha(T): Kratos cell-wise alpha converges O(h^2) to Python Gauss & exact."""
    mods = _load_apps()
    KM = mods["KratosMultiphysics"]
    CL = mods["ConstitutiveLawsApplication"]

    E0, E1, alpha0, alpha1, dT = 100.0, 0.5, 1.0e-5, 2.5e-7, 100.0
    E_fn = lambda T: E0 + E1 * T  # noqa: E731
    alpha_fn = lambda T: alpha0 + alpha1 * T  # noqa: E731

    beta = E1 * dT / E0
    sig_exact = -(E0 * math.sqrt(beta) / math.atan(math.sqrt(beta))) * (alpha0 * dT / 3.0 + alpha1 * dT**2 / 5.0)

    # Continuous Python Gauss-quadrature reference on p=4, n_x=24 matches exact closed form to < 1e-12
    box_ref = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(24, 1, 1), degree=4)
    t_cp_ref = _quad_x_bspline_coeffs(box_ref, 0.0, 0.0, dT)
    all_cps = np.arange(box_ref.n_ctrl, dtype=np.int64)
    ends_ref = np.concatenate([box_ref.boundary_control_points(0, 0), box_ref.boundary_control_points(0, 1)])
    fixed_ref = np.unique(np.concatenate([3 * ends_ref + 0, 3 * all_cps + 1, 3 * all_cps + 2]))
    res_ref = solve_thermoelastic_tdep(box_ref, t_cp_ref, E_fn, 0.0, alpha_fn, fixed_ref, T_ref=0.0, n_gauss=6)
    gp_ref = np.vstack([np.asarray(cell_gauss_points(box_ref, i, 0, 0)) for i in range(24)])
    sig_py_vals = recover_thermoelastic_stress_tdep(box_ref, res_ref.u, t_cp_ref, E_fn, 0.0, alpha_fn, gp_ref[:, :3])[:, 0]
    sig_py_mean = float(np.sum(gp_ref[:, 3] * sig_py_vals) / np.sum(gp_ref[:, 3]))
    assert abs(sig_py_mean - sig_exact) / abs(sig_exact) < 1e-12

    # Kratos: cell-wise Properties cloning for E(T_c), alpha(T_c) converges O(h^2) to Python Gauss & exact
    errs_vs_py = []
    for n_x in (8, 16, 32):
        model, mp, box, vol = _setup_bar_model(KM, n_x=n_x, degree=2)
        t_cp = _quad_x_bspline_coeffs(box, 0.0, 0.0, dT)
        els = []
        for i in range(n_x):
            xc = (i + 0.5) / n_x
            Tc = dT * xc * xc
            props_i = mp.GetProperties()[i + 1]
            props_i.SetValue(KM.YOUNG_MODULUS, E_fn(Tc))
            props_i.SetValue(KM.POISSON_RATIO, 0.0)
            props_i.SetValue(KM.THERMAL_EXPANSION_COEFFICIENT, alpha_fn(Tc))
            props_i.SetValue(KM.REFERENCE_TEMPERATURE, 0.0)
            props_i.SetValue(KM.CONSTITUTIVE_LAW, CL.ThermalElasticIsotropic3D())
            els.extend(
                attach_cell_elements(
                    KM,
                    mp,
                    vol,
                    box,
                    props_i,
                    element_name="SmallDisplacementElement3D4N",
                    start_id=i + 1,
                    cells=[(i, 0, 0)],
                )
            )
        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, float(t_cp[node.Id - 1]))
            node.Fix(KM.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
            node.Fix(KM.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
            if abs(node.X0) < 1e-12 or abs(node.X0 - 1.0) < 1e-12:
                node.Fix(KM.DISPLACEMENT_X)
                node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        solve_linear(KM, mp)
        sigs = np.array([np.array(el.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo))[:, 0] for el in els])
        sig_mean = float(np.mean(sigs))
        rel_err = abs(sig_mean - sig_py_mean) / abs(sig_py_mean)
        errs_vs_py.append(rel_err)
        print(f"S20 dual quad T(x) Kratos clone n_x={n_x}: sig={sig_mean:.12f}, py_gauss={sig_py_mean:.12f}, rel_err={rel_err:.3e}")

    assert errs_vs_py[0] / errs_vs_py[1] > 3.5
    assert errs_vs_py[1] / errs_vs_py[2] > 3.5
    assert errs_vs_py[2] < 5e-4
