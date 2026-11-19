# -*- coding: utf-8 -*-
"""Phase 6 SBM Comparison Track (Master Plan §10.5).

Compares Kratos IgaApplication Shifted Boundary Method (``NurbsGeometryModelerSbm`` +
``IgaModelerSbm`` + ``SolidElement`` + ``SbmSolidCondition`` / ``SbmLoadSolidCondition``)
against WHT_IGA's ``cut + extension`` and ``FCM alpha`` paths at matched DOF counts.

Exact Kratos 10.4.3 runtime findings verified by this suite (no mocking):
1. ``SolidElement``, ``SbmSolidCondition``, and ``SbmLoadSolidCondition`` are
   2D-only in Kratos 10.4.3 (raising ``RuntimeError`` for 3D ``NurbsVolume``).
   Consequently, V06 (3D pressurized sphere) cannot be solved with Kratos SBM,
   so cross-method head-to-head comparison is performed on:
   - V02 Roark plate with hole (2D plane-stress SBM vs 3D nu=0 nz=1 cut+ext / FCM alpha)
   - 2D Lamé pressurized ring (a=0.5, b=1.0, E=1, nu=0, p_i=1e-6) across all 3 methods
2. ``NurbsCurveOnSurfaceGeometry::ComputeAxisIntersectionSBM`` uses a hardcoded
   ``1e-15`` intersection tolerance, causing non-dyadic knot spans (e.g. ``200/48``)
   to fail with ``Wrong number of intersection found (<2)``. Dyadic/integer-step
   grids (e.g. ``(40,8)``, ``(50,10)``, ``(80,16)``, ``(100,20)``) avoid this.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
import time
from typing import Any

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import system_matrix_csr
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import bc_roark as _bc
from wht_iga.solver import probes as PR
from wht_iga.solver import trimmed as T
from wht_iga.solver.active import fix_inactive
from wht_iga.solver.penalty import rule_for_degree
from wht_iga.solver.pipeline import run_linear
from wht_iga.solver.trimmed import _VOIGT, _d_matrix, _strain_at_points
from wht_iga.solver.validate import JobSpec
from wht_iga.stabilize import cond as Q

pytestmark = pytest.mark.kratos

ROARK_STEP = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "roark_plate_hole.step"
)
ROARK_LO = (-100.0, -20.0, -0.5)
ROARK_HI = (100.0, 20.0, 0.5)
ROARK_PEAK = 4.3175e-6
ROARK_ABAQUS = 4.34796947956784e-6


def _extract_free_system_matrix(strategy: Any) -> Any:
    """Extract free-DOF sparse system matrix from a Kratos strategy."""
    K = system_matrix_csr(strategy).tocsr()
    K.eliminate_zeros()
    nnz = np.diff(K.indptr)
    rows = np.repeat(np.arange(K.shape[0]), nnz)
    only_diag = np.bincount(rows[K.indices == rows], minlength=K.shape[0]) == nnz
    fixed = (nnz == 1) & only_diag
    free = np.flatnonzero(~fixed)
    return K[free][:, free]


def _compute_condition_number(K_sparse: Any, symmetric: bool = True) -> float:
    """Compute 2-norm condition number using exact SVD (<=1500 symmetric / <=2200 unsym) or eigsh."""
    n = int(K_sparse.shape[0])
    if symmetric and n > 1500:
        return Q.cond_eigsh(K_sparse)[0]
    if n <= 2200:
        return Q.cond_dense(K_sparse.toarray())
    if symmetric:
        return Q.cond_eigsh(K_sparse)[0]
    return float("nan")


def _solve_kratos_sbm_2d(
    *,
    lo_xy: tuple[float, float],
    hi_xy: tuple[float, float],
    n_cells_xy: tuple[int, int],
    degree: int,
    outer_lines: list[dict[str, Any]],
    outer_layer_map: list[dict[str, str]],
    inner_lines: list[dict[str, Any]] | None,
    inner_layer_map: list[dict[str, str]] | None,
    bc_setter: Any,
    E: float = 1.0,
    nu: float = 0.0,
    penalty_factor: float = -1.0,
    lambda_inner: float = 0.5,
    lambda_outer: float = 0.5,
    n_import_pts: int = 400,
) -> tuple[BoxSpline, np.ndarray, Any, float]:
    """Build and solve a 2D Kratos SBM solid elasticity model.

    Parameters
    ----------
    lo_xy, hi_xy : tuple[float, float]
        Background NURBS bounding box in XY.
    n_cells_xy : tuple[int, int]
        Number of knot spans in (u, v).
    degree : int
        Polynomial order along both axes.
    outer_lines, inner_lines : list[dict]
        NURBS curve JSON entries for ``ImportNurbsSbmModeler``.
    outer_layer_map, inner_layer_map : list[dict]
        Layer-to-condition mappings for ``ImportNurbsSbmModeler``.
    bc_setter : callable
        Callback ``bc_setter(KM, IGA, skin_model_part)`` setting ``DISPLACEMENT``/``DIRECTION``/``FORCE``.
    E, nu, penalty_factor : float
        Linear elastic plane-stress parameters and SBM ``PENALTY_FACTOR`` (``-1.0`` = penalty-free).
    lambda_inner, lambda_outer : float
        SBM surrogate offset factors passed to ``NurbsGeometryModelerSbm``.
    n_import_pts : int
        Sampling points per imported NURBS curve in ``SnakeSbmProcess``.

    Returns
    -------
    tuple[BoxSpline, np.ndarray, csr_matrix, float]
        Lifted 3D ``BoxSpline`` (nz=1), full DOF coefficient vector, system CSR matrix, solve wall time.
    """
    import KratosMultiphysics.LinearSolversApplication as LA
    import KratosMultiphysics.scipy_conversion_tools as sct
    from KratosMultiphysics.modeler_factory import KratosModelerFactory

    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    IGA = mods["IgaApplication"]

    t0 = time.perf_counter()
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f_out:
        json.dump({"Lines": outer_lines}, f_out)
        p_out = f_out.name.replace("\\", "/")

    p_in: str | None = None
    if inner_lines:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f_in:
            json.dump({"Lines": inner_lines}, f_in)
            p_in = f_in.name.replace("\\", "/")

    try:
        model = KM.Model()
        iga_mp = model.CreateModelPart("IgaModelPart")
        iga_mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        iga_mp.AddNodalSolutionStepVariable(KM.REACTION)
        iga_mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        iga_mp.SetBufferSize(2)

        modelers_cfg: list[dict[str, Any]] = [
            {
                "modeler_name": "ImportNurbsSbmModeler",
                "Parameters": {
                    "input_filename": p_out,
                    "model_part_name": "skinModelPart_outer_initial",
                    "link_layer_to_condition_name": outer_layer_map,
                },
            }
        ]
        if p_in is not None and inner_layer_map is not None:
            modelers_cfg.append(
                {
                    "modeler_name": "ImportNurbsSbmModeler",
                    "Parameters": {
                        "input_filename": p_in,
                        "model_part_name": "skinModelPart_inner_initial",
                        "link_layer_to_condition_name": inner_layer_map,
                    },
                }
            )

        nurbs_params: dict[str, Any] = {
            "model_part_name": "IgaModelPart",
            "lower_point_xyz": [float(lo_xy[0]), float(lo_xy[1]), 0.0],
            "upper_point_xyz": [float(hi_xy[0]), float(hi_xy[1]), 0.0],
            "lower_point_uvw": [float(lo_xy[0]), float(lo_xy[1]), 0.0],
            "upper_point_uvw": [float(hi_xy[0]), float(hi_xy[1]), 0.0],
            "polynomial_order": [int(degree), int(degree)],
            "number_of_knot_spans": [int(n_cells_xy[0]), int(n_cells_xy[1])],
            "lambda_outer": float(lambda_outer),
            "number_of_inner_loops": 1 if inner_lines else 0,
            "number_initial_points_if_importing_nurbs": int(n_import_pts),
            "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",
            "skin_model_part_name": "skin_model_part",
            "echo_level": 0,
        }
        if inner_lines:
            nurbs_params["lambda_inner"] = float(lambda_inner)
            nurbs_params["skin_model_part_inner_initial_name"] = "skinModelPart_inner_initial"
        modelers_cfg.append({"modeler_name": "NurbsGeometryModelerSbm", "Parameters": nurbs_params})

        elem_cond_list: list[dict[str, Any]] = [
            {
                "geometry_type": "GeometrySurface",
                "iga_model_part": "ComputationalDomain",
                "type": "element",
                "name": "SolidElement",
                "shape_function_derivatives_order": 2,
                "variables": [{"variable_name": "BODY_FORCE", "value": ["0.0", "0.0", "0.0"]}],
            },
            {
                "geometry_type": "SurfaceEdge",
                "iga_model_part": "SBM_Support_outer",
                "type": "condition",
                "name": "SbmCondition",
                "shape_function_derivatives_order": 2 * int(degree) + 1,
                "sbm_parameters": {"is_inner": False},
            },
        ]
        if inner_lines:
            elem_cond_list.append(
                {
                    "geometry_type": "SurfaceEdge",
                    "iga_model_part": "SBM_Support_inner",
                    "type": "condition",
                    "name": "SbmCondition",
                    "shape_function_derivatives_order": 2 * int(degree) + 1,
                    "sbm_parameters": {"is_inner": True},
                }
            )
        modelers_cfg.append(
            {
                "modeler_name": "IgaModelerSbm",
                "Parameters": {
                    "echo_level": 0,
                    "skin_model_part_name": "skin_model_part",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": elem_cond_list,
                },
            }
        )

        factory = KratosModelerFactory()
        modelers = factory.ConstructListOfModelers(model, KM.Parameters(json.dumps(modelers_cfg)))
        for mod in modelers:
            mod.SetupGeometryModel()
        for mod in modelers:
            mod.PrepareGeometryModel()
        for mod in modelers:
            mod.SetupModelPart()

        skin = model["skin_model_part"]
        bc_setter(KM, IGA, skin)

        props_json = KM.Parameters(
            json.dumps(
                {
                    "properties": [
                        {
                            "model_part_name": "IgaModelPart",
                            "properties_id": 0,
                            "Material": {
                                "name": "lin_el",
                                "constitutive_law": {"name": "LinearElasticPlaneStress2DLaw"},
                                "Variables": {
                                    "PENALTY_FACTOR": float(penalty_factor),
                                    "THICKNESS": 1.0,
                                    "YOUNG_MODULUS": float(E),
                                    "POISSON_RATIO": float(nu),
                                    "DENSITY": 0.0,
                                },
                                "Tables": {},
                            },
                        }
                    ]
                }
            )
        )
        KM.ReadMaterialsUtility(model).ReadMaterials(props_json)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, iga_mp)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, iga_mp)

        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        solver = LA.PardisoLUSolver()
        builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
        strategy = KM.ResidualBasedLinearStrategy(iga_mp, scheme, builder, False, False, False, False)
        strategy.SetEchoLevel(0)
        strategy.Initialize()
        iga_mp.CloneTimeStep(1.0)
        strategy.Solve()

        box = BoxSpline(
            lo=(float(lo_xy[0]), float(lo_xy[1]), -0.5),
            hi=(float(hi_xy[0]), float(hi_xy[1]), 0.5),
            n_cells=(int(n_cells_xy[0]), int(n_cells_xy[1]), 1),
            degree=int(degree),
        )
        ncp_x = int(n_cells_xy[0]) + int(degree)
        ncp_y = int(n_cells_xy[1]) + int(degree)
        coeffs = np.zeros(box.n_dofs, dtype=float)
        coeffs3 = coeffs.reshape(-1, 3)
        for node in iga_mp.Nodes:
            nid = int(node.Id) - 1
            ix, iy = nid % ncp_x, nid // ncp_x
            ux = float(node.GetSolutionStepValue(KM.DISPLACEMENT_X))
            uy = float(node.GetSolutionStepValue(KM.DISPLACEMENT_Y))
            for iz in range(box.n_cp[2]):
                idx3 = (iz * ncp_y + iy) * ncp_x + ix
                coeffs3[idx3, 0] = ux
                coeffs3[idx3, 1] = uy

        K_csr = sct.to_csr(strategy.GetSystemMatrix()).tocsr()
        dt = time.perf_counter() - t0
        return box, coeffs, K_csr, dt
    finally:
        if os.path.exists(p_out):
            os.remove(p_out)
        if p_in is not None and os.path.exists(p_in):
            os.remove(p_in)


def _ring_quarter_tris(
    a: float, b: float, z0: float, z1: float, n_th: int = 64
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build watertight outward-oriented surface triangles for a quarter annulus."""
    ths = np.linspace(0.0, 0.5 * math.pi, n_th + 1)
    c, s = np.cos(ths), np.sin(ths)
    nodes: list[tuple[float, float, float]] = []

    def add_ring(r: float, z: float) -> int:
        idx0 = len(nodes)
        for ci, si in zip(c, s):
            nodes.append((float(r * ci), float(r * si), float(z)))
        return idx0

    i_az0 = add_ring(a, z0)
    i_bz0 = add_ring(b, z0)
    i_bz1 = add_ring(b, z1)
    i_az1 = add_ring(a, z1)
    tris: list[tuple[int, int, int]] = []
    inner_tris: list[tuple[int, int, int]] = []
    for k in range(n_th):
        tris.extend([(i_az0 + k, i_bz0 + k + 1, i_bz0 + k), (i_az0 + k, i_az0 + k + 1, i_bz0 + k + 1)])
        tris.extend([(i_az1 + k, i_bz1 + k, i_bz1 + k + 1), (i_az1 + k, i_bz1 + k + 1, i_az1 + k + 1)])
        tris.extend([(i_bz0 + k, i_bz0 + k + 1, i_bz1 + k + 1), (i_bz0 + k, i_bz1 + k + 1, i_bz1 + k)])
        t1 = (i_az0 + k, i_az1 + k + 1, i_az0 + k + 1)
        t2 = (i_az0 + k, i_az1 + k, i_az1 + k + 1)
        tris.extend([t1, t2])
        inner_tris.extend([t1, t2])
    tris.extend([(i_az0, i_bz0, i_bz1), (i_az0, i_bz1, i_az1)])
    tris.extend([(i_az0 + n_th, i_bz1 + n_th, i_bz0 + n_th), (i_az0 + n_th, i_az1 + n_th, i_bz1 + n_th)])
    return np.asarray(nodes, dtype=float), np.asarray(tris, dtype=np.int64), np.asarray(inner_tris, dtype=np.int64)


def test_kratos_sbm_3d_solid_limitation(kratos):
    """Verify exact Kratos 10.4.3 2D-only limitation for Solid SBM & 1e-15 intersection edge case."""
    KM = kratos["KratosMultiphysics"]
    IGA = kratos["IgaApplication"]
    SM = kratos["StructuralMechanicsApplication"]
    from KratosMultiphysics.modeler_factory import KratosModelerFactory

    # 1. 3D NurbsVolumeGeometry + SolidElement raises explicit 2D-only RuntimeError.
    m3 = KM.Model()
    mp3 = m3.CreateModelPart("M3")
    mp3.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp3.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    nodes = KM.NodesVector()
    nid = 1
    for z in (0.0, 1.0):
        for y in (0.0, 1.0):
            for x in (0.0, 1.0):
                n = mp3.CreateNewNode(nid, x, y, z)
                n.AddDof(KM.DISPLACEMENT_X)
                n.AddDof(KM.DISPLACEMENT_Y)
                n.AddDof(KM.DISPLACEMENT_Z)
                nodes.append(n)
                nid += 1
    kv = KM.Vector([0.0, 0.0, 1.0, 1.0])
    vol = KM.NurbsVolumeGeometry(nodes, 1, 1, 1, kv, kv, kv)
    qps = KM.GeometriesVector()
    vol.CreateQuadraturePointGeometries(qps, 2)
    prop3 = mp3.CreateNewProperties(1)
    prop3.SetValue(KM.YOUNG_MODULUS, 1.0)
    prop3.SetValue(KM.POISSON_RATIO, 0.0)
    prop3.SetValue(KM.THICKNESS, 1.0)
    prop3.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElasticPlaneStrain2DLaw())
    el3 = mp3.CreateNewElement("SolidElement", 1, qps[0], prop3)
    el3.Initialize(mp3.ProcessInfo)
    lhs, rhs = KM.Matrix(), KM.Vector()
    with pytest.raises(RuntimeError, match=r"SolidElement momentarily only supports 2D elements"):
        el3.CalculateLocalSystem(lhs, rhs, mp3.ProcessInfo)

    # 2. 3D NurbsGeometryModelerSbm + IgaModelerSbm produces a 2D-only failure on SolidElement
    #    and an unpopulated (0x0) local system on SbmSolidCondition; direct 3D NurbsVolume QP
    #    on SbmSolidCondition raises the explicit 2D-only RuntimeError.
    m_sbm3 = KM.Model()
    iga3 = m_sbm3.CreateModelPart("IgaModelPart")
    iga3.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    iga3.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    skin3 = m_sbm3.CreateModelPart("skin_outer_3d")
    skin3.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    skin3.CreateNewProperties(1)
    p = skin3.GetProperties()[1]
    coords = [
        (0.25, 0.25, 0.25), (1.75, 0.25, 0.25), (1.75, 1.75, 0.25), (0.25, 1.75, 0.25),
        (0.25, 0.25, 1.75), (1.75, 0.25, 1.75), (1.75, 1.75, 1.75), (0.25, 1.75, 1.75),
    ]
    for i, (cx, cy, cz) in enumerate(coords, start=1):
        skin3.CreateNewNode(i, cx, cy, cz)
    tris = [
        [1, 3, 2], [1, 4, 3], [5, 6, 7], [5, 7, 8],
        [1, 2, 6], [1, 6, 5], [4, 8, 7], [4, 7, 3],
        [1, 5, 8], [1, 8, 4], [2, 7, 6], [2, 3, 7],
    ]
    for cid, conn in enumerate(tris, start=1):
        skin3.CreateNewCondition("SurfaceCondition3D3N", cid, conn, p)
    cfg3 = KM.Parameters(
        json.dumps(
            [
                {
                    "modeler_name": "NurbsGeometryModelerSbm",
                    "Parameters": {
                        "model_part_name": "IgaModelPart",
                        "lower_point_xyz": [0.0, 0.0, 0.0],
                        "upper_point_xyz": [2.0, 2.0, 2.0],
                        "lower_point_uvw": [0.0, 0.0, 0.0],
                        "upper_point_uvw": [2.0, 2.0, 2.0],
                        "polynomial_order": [1, 1, 1],
                        "number_of_knot_spans": [4, 4, 4],
                        "lambda_outer": 0.5,
                        "number_of_inner_loops": 0,
                        "skin_model_part_outer_initial_name": "skin_outer_3d",
                        "skin_model_part_name": "skin_3d",
                        "echo_level": 0,
                    },
                },
                {
                    "modeler_name": "IgaModelerSbm",
                    "Parameters": {
                        "echo_level": 0,
                        "skin_model_part_name": "skin_3d",
                        "analysis_model_part_name": "IgaModelPart",
                        "element_condition_list": [
                            {
                                "geometry_type": "GeometrySurface",
                                "iga_model_part": "ComputationalDomain",
                                "type": "element",
                                "name": "SolidElement",
                                "shape_function_derivatives_order": 2,
                                "variables": [{"variable_name": "BODY_FORCE", "value": ["0.0", "0.0", "0.0"]}],
                            },
                            {
                                "geometry_type": "SurfaceEdge",
                                "iga_model_part": "SBM_Support_outer",
                                "type": "condition",
                                "name": "SbmSolidCondition",
                                "shape_function_derivatives_order": 2,
                                "sbm_parameters": {"is_inner": False},
                            },
                        ],
                    },
                },
            ]
        )
    )
    modelers = KratosModelerFactory().ConstructListOfModelers(m_sbm3, cfg3)
    for mod in modelers:
        mod.SetupGeometryModel()
    for mod in modelers:
        mod.PrepareGeometryModel()
    for mod in modelers:
        mod.SetupModelPart()
    prop_sbm3 = iga3.CreateNewProperties(1)
    prop_sbm3.SetValue(KM.YOUNG_MODULUS, 1.0)
    prop_sbm3.SetValue(KM.POISSON_RATIO, 0.0)
    prop_sbm3.SetValue(KM.THICKNESS, 1.0)
    prop_sbm3.SetValue(IGA.PENALTY_FACTOR, -1.0)
    prop_sbm3.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElasticPlaneStrain2DLaw())
    el_sbm3 = m_sbm3.GetModelPart("IgaModelPart.ComputationalDomain").Elements[1]
    el_sbm3.Properties = prop_sbm3
    el_sbm3.Initialize(iga3.ProcessInfo)
    with pytest.raises(RuntimeError, match=r"SolidElement momentarily only supports 2D elements"):
        el_sbm3.CalculateLocalSystem(lhs, rhs, iga3.ProcessInfo)
    cond3 = m_sbm3.GetModelPart("IgaModelPart.SBM_Support_outer").Conditions[1]
    cond3.Properties = prop_sbm3
    cond3.Initialize(iga3.ProcessInfo)
    cond3.CalculateLocalSystem(lhs, rhs, iga3.ProcessInfo)
    assert lhs.Size1() == 0 and rhs.Size() == 0
    c_vol3 = mp3.CreateNewCondition("SbmSolidCondition", 1, qps[0], prop_sbm3)
    with pytest.raises(RuntimeError, match=r"SbmSolidCondition momentarily only supports 2D conditions"):
        c_vol3.Initialize(mp3.ProcessInfo)

    # 3. Non-dyadic knot spans trigger ComputeAxisIntersectionSBM (<2) 1e-15 tolerance bug.
    outer_lines = [
        {"Layer": "bottom", "Weights": [1, 1], "CPCoordinates": [[-100.0, -20.0, 0.0], [100.0, -20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 1, "pDegree": 1},
        {"Layer": "right", "Weights": [1, 1], "CPCoordinates": [[100.0, -20.0, 0.0], [100.0, 20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 2, "pDegree": 1},
        {"Layer": "top", "Weights": [1, 1], "CPCoordinates": [[100.0, 20.0, 0.0], [-100.0, 20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 3, "pDegree": 1},
        {"Layer": "left", "Weights": [1, 1], "CPCoordinates": [[-100.0, 20.0, 0.0], [-100.0, -20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 4, "pDegree": 1},
    ]
    w = math.sqrt(0.5)
    R = 10.0
    inner_lines = [
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[R, 0.0, 0.0], [R, R, 0.0], [0.0, R, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 1, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[0.0, R, 0.0], [-R, R, 0.0], [-R, 0.0, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 2, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[-R, 0.0, 0.0], [-R, -R, 0.0], [0.0, -R, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 3, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[0.0, -R, 0.0], [R, -R, 0.0], [R, 0.0, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 4, "pDegree": 2},
    ]
    with pytest.raises(RuntimeError, match=r"ComputeAxisIntersectionSBM :: Wrong number of intersection found"):
        _solve_kratos_sbm_2d(
            lo_xy=(-100.0, -20.0),
            hi_xy=(100.0, 20.0),
            n_cells_xy=(48, 10),
            degree=2,
            outer_lines=outer_lines,
            outer_layer_map=[
                {"layer_name": "left", "condition_name": "SbmSolidCondition"},
                {"layer_name": "right", "condition_name": "SbmLoadSolidCondition"},
                {"layer_name": "top", "condition_name": "SbmLoadSolidCondition"},
                {"layer_name": "bottom", "condition_name": "SbmLoadSolidCondition"},
            ],
            inner_lines=inner_lines,
            inner_layer_map=[{"layer_name": "hole", "condition_name": "SbmLoadSolidCondition"}],
            bc_setter=lambda *_: None,
        )


def test_kratos_sbm_2d_uniaxial_patch_exact():
    """Immersed 2D [0,2]^2 uniaxial tension patch inside [-0.5,2.5]^2 reproduces linear field."""
    outer_lines = [
        {"Layer": "bottom", "Weights": [1, 1], "CPCoordinates": [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 1, "pDegree": 1},
        {"Layer": "right", "Weights": [1, 1], "CPCoordinates": [[2.0, 0.0, 0.0], [2.0, 2.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 2, "pDegree": 1},
        {"Layer": "top", "Weights": [1, 1], "CPCoordinates": [[2.0, 2.0, 0.0], [0.0, 2.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 3, "pDegree": 1},
        {"Layer": "left", "Weights": [1, 1], "CPCoordinates": [[0.0, 2.0, 0.0], [0.0, 0.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 4, "pDegree": 1},
    ]
    layer_map = [
        {"layer_name": "left", "condition_name": "SbmSolidCondition"},
        {"layer_name": "right", "condition_name": "SbmLoadSolidCondition"},
        {"layer_name": "top", "condition_name": "SbmLoadSolidCondition"},
        {"layer_name": "bottom", "condition_name": "SbmLoadSolidCondition"},
    ]

    def set_bcs(KM, _IGA, skin):
        for n in skin.GetSubModelPart("outer.left").Nodes:
            n.SetValue(KM.DISPLACEMENT, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.right").Nodes:
            n.SetValue(KM.FORCE, [1.0e-3, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.top").Nodes:
            n.SetValue(KM.FORCE, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.bottom").Nodes:
            n.SetValue(KM.FORCE, [0.0, 0.0, 0.0])

    box, coeffs, K_csr, dt = _solve_kratos_sbm_2d(
        lo_xy=(-0.5, -0.5),
        hi_xy=(2.5, 2.5),
        n_cells_xy=(10, 10),
        degree=2,
        outer_lines=outer_lines,
        outer_layer_map=layer_map,
        inner_lines=None,
        inner_layer_map=None,
        bc_setter=set_bcs,
        E=1.0,
        nu=0.0,
        penalty_factor=-1.0,
    )
    gx = np.linspace(0.1, 1.9, 5)
    gy = np.linspace(0.1, 1.9, 5)
    XX, YY = np.meshgrid(gx, gy)
    pts = np.column_stack([XX.ravel(), YY.ravel(), np.zeros(XX.size)])
    u = box.evaluate(coeffs.reshape(-1, 3), pts)
    err_x = float(np.max(np.abs(u[:, 0] - 1.0e-3 * pts[:, 0])))
    err_y = float(np.max(np.abs(u[:, 1])))
    kappa = _compute_condition_number(K_csr, symmetric=False)
    print(f"[SBM Patch] dofs={K_csr.shape[0]} cond={kappa:.2f} max|ux-1e-3x|={err_x:.3e} max|uy|={err_y:.3e} t={dt:.2f}s", flush=True)
    assert err_x < 1e-13
    assert err_y < 1e-13
    assert np.isfinite(kappa) and kappa < 1e4


def test_lame_ring_sbm_vs_cut_ext_vs_fcm(kratos):
    """2D Lame pressurized ring (a=0.5, b=1.0, E=1, nu=0, pi=1e-6): cut+ext vs FCM alpha vs SBM."""
    import KratosMultiphysics.LinearSolversApplication as LA

    KM = kratos["KratosMultiphysics"]
    IGA = kratos["IgaApplication"]
    SM = kratos["StructuralMechanicsApplication"]

    a, b, E, nu, pi = 0.5, 1.0, 1.0, 0.0, 1.0e-6
    c_lame = pi * a**2 / (b**2 - a**2)
    k_lame = pi * a**2 * b**2 / (b**2 - a**2)
    D_mat = _d_matrix(E, nu)

    rs = np.linspace(a + 0.05 * (b - a), b - 0.05 * (b - a), 12)
    ths = np.linspace(0.05, 0.5 * math.pi - 0.05, 25)
    RR, TT = np.meshgrid(rs, ths)
    pts_eval = np.column_stack([RR.ravel() * np.cos(TT.ravel()), RR.ravel() * np.sin(TT.ravel()), np.zeros(RR.size)])
    r_eval = np.hypot(pts_eval[:, 0], pts_eval[:, 1])
    er_eval = pts_eval[:, :2] / r_eval[:, None]
    ur_ref = (c_lame * r_eval + k_lame / r_eval) / E
    u_ref = ur_ref[:, None] * er_eval
    srr_ref = c_lame - k_lame / r_eval**2
    stt_ref = c_lame + k_lame / r_eval**2

    def eval_ring_errors(box: BoxSpline, coeffs: np.ndarray) -> tuple[float, float, float]:
        u_num = box.evaluate(coeffs.reshape(-1, 3), pts_eval)[:, :2]
        u_l2 = float(np.linalg.norm(u_num - u_ref) / np.linalg.norm(u_ref))
        S = _strain_at_points(box, coeffs, pts_eval) @ D_mat.T
        srr_num = S[:, _VOIGT["xx"]] * er_eval[:, 0] ** 2 + S[:, _VOIGT["yy"]] * er_eval[:, 1] ** 2 + 2.0 * S[:, _VOIGT["xy"]] * er_eval[:, 0] * er_eval[:, 1]
        stt_num = S[:, _VOIGT["xx"]] * er_eval[:, 1] ** 2 + S[:, _VOIGT["yy"]] * er_eval[:, 0] ** 2 - 2.0 * S[:, _VOIGT["xy"]] * er_eval[:, 0] * er_eval[:, 1]
        srr_l2 = float(np.linalg.norm(srr_num - srr_ref) / np.linalg.norm(srr_ref))
        stt_l2 = float(np.linalg.norm(stt_num - stt_ref) / np.linalg.norm(stt_ref))
        return u_l2, srr_l2, stt_l2

    nodes_m, tris_m, inner_tris = _ring_quarter_tris(a, b, -0.05, 0.05, n_th=64)
    in_tris_xyz = [nodes_m[t] for t in inner_tris]

    rows: list[dict[str, Any]] = []
    for nc3 in [(6, 6, 1), (10, 10, 1), (14, 14, 1)]:
        box = BoxSpline(lo=(0.0, 0.0, -0.05), hi=(1.0, 1.0, 0.05), n_cells=nc3, degree=2)
        qp_table, out_cells, _issues, cell_records, _diag = T.build_trimmed_table(nodes_m, tris_m, box, q=4)
        out_set = {tuple(c) for c in out_cells}
        for method, use_ext, out_mode in [("cut+ext", True, "drop"), ("FCM_alpha", False, "ghost")]:
            t0 = time.perf_counter()
            model = KM.Model()
            mp = model.CreateModelPart("Ring")
            mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            mp.AddNodalSolutionStepVariable(KM.REACTION)
            mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            volume, _ = make_nurbs_volume(KM, mp, box)
            vol_props = mp.CreateNewProperties(1)
            vol_props.SetValue(KM.YOUNG_MODULUS, E)
            vol_props.SetValue(KM.POISSON_RATIO, nu)
            vol_props.SetValue(KM.DENSITY, 0.0)
            vol_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
            ghost_props = None
            if out_cells and out_mode == "ghost":
                ghost_props = mp.CreateNewProperties(2)
                ghost_props.SetValue(KM.YOUNG_MODULUS, E * 1e-8)
                ghost_props.SetValue(KM.POISSON_RATIO, nu)
                ghost_props.SetValue(KM.DENSITY, 0.0)
                ghost_props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
            nx, ny, nz = box.n_cells
            eid = 1
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        if out_mode == "drop" and (i, j, k) in out_set:
                            eid += 1
                            continue
                        pts = (
                            np.asarray(qp_table[(i, j, k)], dtype=float).tolist()
                            if (i, j, k) in qp_table
                            else cell_gauss_points(box, i, j, k)
                        )
                        qps = KM.GeometriesVector()
                        volume.CreateQuadraturePointGeometries(qps, 2, pts)
                        mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], ghost_props if (i, j, k) in out_set else vol_props)
                        eid += 1
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, mp)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, mp)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, mp)
            comps = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
            for axis in range(3):
                for cp_id in box.boundary_control_points(axis, 0):
                    mp.GetNode(int(cp_id) + 1).Fix(comps[axis])
            if out_mode == "drop" and out_cells:
                fix_inactive(KM, mp, box, out_set)
            order = int(rule_for_degree(2))
            qp_geoms, wrappers, _, _ = make_surface_quadrature(
                KM, mp, volume, in_tris_xyz, rule_order=order, first_node_id=box.n_ctrl + 1
            )
            pid = 10
            for qp in qp_geoms:
                props = mp.CreateNewProperties(pid)
                pid += 1
                cond = mp.CreateNewCondition("LoadCondition", 10_000 + pid, qp, props)
                cond.SetValue(KM.PRESSURE, float(-pi))
            if use_ext:
                start_cid = max((int(c.Id) for c in mp.MasterSlaveConstraints), default=0) + 1
                _tk = T.apply_extension_ties(KM, mp, box, cell_records, theta=0.3, start_cid=start_cid)
            mp.SetBufferSize(2)
            scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
            builder = KM.ResidualBasedBlockBuilderAndSolver(LA.PardisoLDLTSolver())
            strategy = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
            strategy.SetEchoLevel(0)
            strategy.Initialize()
            mp.CloneTimeStep(1.0)
            strategy.Solve()
            U = np.zeros(box.n_dofs, dtype=float)
            for node in mp.Nodes:
                idx = int(node.Id) - 1
                if 0 <= idx < box.n_ctrl:
                    U[3 * idx] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
                    U[3 * idx + 1] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
                    U[3 * idx + 2] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)
            u_l2, srr_l2, stt_l2 = eval_ring_errors(box, U)
            Kf = _extract_free_system_matrix(strategy)
            kappa = _compute_condition_number(Kf, symmetric=True)
            dt = time.perf_counter() - t0
            rows.append(
                {
                    "method": method,
                    "nc": str(nc3),
                    "dofs": int(Kf.shape[0]),
                    "u_L2": u_l2,
                    "srr_L2": srr_l2,
                    "stt_L2": stt_l2,
                    "cond": kappa,
                    "t_s": dt,
                }
            )
            del wrappers, qp_geoms, volume

    w = math.sqrt(0.5)
    ring_lines = [
        {"Layer": "bottom", "Weights": [1.0, 1.0], "CPCoordinates": [[a, 0.0, 0.0], [b, 0.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 1, "pDegree": 1},
        {"Layer": "outer", "Weights": [1.0, w, 1.0], "CPCoordinates": [[b, 0.0, 0.0], [b, b, 0.0], [0.0, b, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 2, "pDegree": 2},
        {"Layer": "left", "Weights": [1.0, 1.0], "CPCoordinates": [[0.0, b, 0.0], [0.0, a, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 3, "pDegree": 1},
        {"Layer": "inner", "Weights": [1.0, w, 1.0], "CPCoordinates": [[0.0, a, 0.0], [a, a, 0.0], [a, 0.0, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 4, "pDegree": 2},
    ]
    ring_map = [
        {"layer_name": "bottom", "condition_name": "SbmSolidCondition"},
        {"layer_name": "left", "condition_name": "SbmSolidCondition"},
        {"layer_name": "outer", "condition_name": "SbmLoadSolidCondition"},
        {"layer_name": "inner", "condition_name": "SbmLoadSolidCondition"},
    ]

    def set_ring_bcs(KM_m, _IGA_m, skin):
        for n in skin.GetSubModelPart("outer.bottom").Nodes:
            n.SetValue(KM_m.DIRECTION, [0.0, 1.0, 0.0])
            n.SetValue(KM_m.DISPLACEMENT, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.left").Nodes:
            n.SetValue(KM_m.DIRECTION, [1.0, 0.0, 0.0])
            n.SetValue(KM_m.DISPLACEMENT, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.outer").Nodes:
            n.SetValue(KM_m.FORCE, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.inner").Nodes:
            r = math.hypot(n.X, n.Y)
            n.SetValue(KM_m.FORCE, [pi * n.X / r, pi * n.Y / r, 0.0])

    for nc2 in [(15, 15), (25, 25), (31, 31)]:
        box_s, coeffs_s, K_s, dt_s = _solve_kratos_sbm_2d(
            lo_xy=(0.0, 0.0),
            hi_xy=(1.0, 1.0),
            n_cells_xy=nc2,
            degree=2,
            outer_lines=ring_lines,
            outer_layer_map=ring_map,
            inner_lines=None,
            inner_layer_map=None,
            bc_setter=set_ring_bcs,
            E=E,
            nu=nu,
            penalty_factor=-1.0,
        )
        u_l2, srr_l2, stt_l2 = eval_ring_errors(box_s, coeffs_s)
        kappa_s = _compute_condition_number(K_s, symmetric=False)
        rows.append(
            {
                "method": "Kratos_SBM",
                "nc": str(nc2),
                "dofs": int(K_s.shape[0]),
                "u_L2": u_l2,
                "srr_L2": srr_l2,
                "stt_L2": stt_l2,
                "cond": kappa_s,
                "t_s": dt_s,
            }
        )

    print("\n[P6 Lame Ring 2D] method | nc | dofs | u_L2 | srr_L2 | stt_L2 | cond | t_s", flush=True)
    for r in rows:
        print(
            f"[P6 Lame Ring 2D] {r['method']:10s} | {r['nc']:11s} | {r['dofs']:4d} | "
            f"{r['u_L2']:.3e} | {r['srr_L2']:.3e} | {r['stt_L2']:.3e} | {r['cond']:.3e} | {r['t_s']:.2f}",
            flush=True,
        )

    cut_rows = [r for r in rows if r["method"] == "cut+ext"]
    fcm_rows = [r for r in rows if r["method"] == "FCM_alpha"]
    sbm_rows = [r for r in rows if r["method"] == "Kratos_SBM"]
    assert all(r["u_L2"] < 1e-3 for r in cut_rows)
    assert all(r["cond"] < 1e6 for r in cut_rows)
    assert all(r["cond"] > 1e9 for r in fcm_rows)
    assert all(r["u_L2"] < 1e-2 for r in sbm_rows)
    # At matched DOF tiers (~400 and ~900), cut+ext is more accurate than SBM.
    assert cut_rows[0]["u_L2"] < sbm_rows[0]["u_L2"]
    assert cut_rows[1]["u_L2"] < sbm_rows[1]["u_L2"]


def test_v02_roark_plate_hole_sbm_vs_cut_ext_vs_fcm():
    """V02 Roark plate with hole: cut+ext vs FCM alpha vs 2D Kratos SBM at matched DOFs."""
    from wht_iga.geometry import cad as G

    shape = G.read_step(ROARK_STEP)
    rows: list[dict[str, Any]] = []

    for nc3 in [(20, 4, 1), (33, 7, 1)]:
        h_ip = 200.0 / nc3[0]
        box = BoxSpline(lo=ROARK_LO, hi=ROARK_HI, n_cells=nc3, degree=2)
        tess = G.tessellate(shape, deflection=float(h_ip / 40.0))
        qp_table, out_cells, issues, cell_records, _diag = T.build_trimmed_table(tess.nodes, tess.tris, box, q=4)
        for method, ext, out_mode in [("cut+ext", {"theta": 0.3}, "drop"), ("FCM_alpha", None, "ghost")]:
            t0 = time.perf_counter()
            spec = JobSpec(
                name=f"p6_v02_{method}_{nc3[0]}",
                lo=ROARK_LO,
                hi=ROARK_HI,
                n_cells=nc3,
                degree=2,
                E=1.0,
                nu=0.0,
                out_cells=out_cells,
                out_mode=out_mode,
                ghost_alpha=1e-8,
                **_bc.roark_left_pinned(),
            )
            res = run_linear(spec, geo_issues=issues, qp_table=qp_table, cell_records=cell_records, extension=ext)
            Kf = _extract_free_system_matrix(res.strategy)
            kappa = _compute_condition_number(Kf, symmetric=True)
            peak_A = float(
                PR.max_principal_profile(box, res.U, 1.0, 0.0, np.array([[0.0, 10.0, 0.0], [0.0, -10.0, 0.0]])).max()
            )
            u_right = float(box.evaluate(res.U.reshape(-1, 3), np.array([[100.0, 0.0, 0.0]]))[0, 0])
            dt = time.perf_counter() - t0
            rows.append(
                {
                    "method": method,
                    "nc": str(nc3),
                    "dofs": int(Kf.shape[0]),
                    "peak_A": peak_A,
                    "err_roark": abs(peak_A - ROARK_PEAK) / ROARK_PEAK,
                    "err_abq": abs(peak_A - ROARK_ABAQUS) / ROARK_ABAQUS,
                    "u_right": u_right,
                    "cond": kappa,
                    "t_s": dt,
                }
            )

    w = math.sqrt(0.5)
    R = 10.0
    outer_lines = [
        {"Layer": "bottom", "Weights": [1, 1], "CPCoordinates": [[-100.0, -20.0, 0.0], [100.0, -20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 1, "pDegree": 1},
        {"Layer": "right", "Weights": [1, 1], "CPCoordinates": [[100.0, -20.0, 0.0], [100.0, 20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 2, "pDegree": 1},
        {"Layer": "top", "Weights": [1, 1], "CPCoordinates": [[100.0, 20.0, 0.0], [-100.0, 20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 3, "pDegree": 1},
        {"Layer": "left", "Weights": [1, 1], "CPCoordinates": [[-100.0, 20.0, 0.0], [-100.0, -20.0, 0.0]], "knotVector": [0.0, 0.0, 1.0, 1.0], "id": 4, "pDegree": 1},
    ]
    inner_lines = [
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[R, 0.0, 0.0], [R, R, 0.0], [0.0, R, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 1, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[0.0, R, 0.0], [-R, R, 0.0], [-R, 0.0, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 2, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[-R, 0.0, 0.0], [-R, -R, 0.0], [0.0, -R, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 3, "pDegree": 2},
        {"Layer": "hole", "Weights": [1.0, w, 1.0], "CPCoordinates": [[0.0, -R, 0.0], [R, -R, 0.0], [R, 0.0, 0.0]], "knotVector": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], "id": 4, "pDegree": 2},
    ]

    def set_v02_bcs(KM_m, _IGA_m, skin):
        for n in skin.GetSubModelPart("outer.left").Nodes:
            n.SetValue(KM_m.DISPLACEMENT, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.right").Nodes:
            n.SetValue(KM_m.FORCE, [1.0e-6, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.top").Nodes:
            n.SetValue(KM_m.FORCE, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("outer.bottom").Nodes:
            n.SetValue(KM_m.FORCE, [0.0, 0.0, 0.0])
        for n in skin.GetSubModelPart("inner.hole").Nodes:
            n.SetValue(KM_m.FORCE, [0.0, 0.0, 0.0])

    for nc2 in [(40, 8), (50, 10), (66, 14), (100, 20)]:
        box_s, coeffs_s, K_s, dt_s = _solve_kratos_sbm_2d(
            lo_xy=(-100.0, -20.0),
            hi_xy=(100.0, 20.0),
            n_cells_xy=nc2,
            degree=2,
            outer_lines=outer_lines,
            outer_layer_map=[
                {"layer_name": "left", "condition_name": "SbmSolidCondition"},
                {"layer_name": "right", "condition_name": "SbmLoadSolidCondition"},
                {"layer_name": "top", "condition_name": "SbmLoadSolidCondition"},
                {"layer_name": "bottom", "condition_name": "SbmLoadSolidCondition"},
            ],
            inner_lines=inner_lines,
            inner_layer_map=[{"layer_name": "hole", "condition_name": "SbmLoadSolidCondition"}],
            bc_setter=set_v02_bcs,
            E=1.0,
            nu=0.0,
            penalty_factor=-1.0,
        )
        peak_A = float(
            PR.max_principal_profile(box_s, coeffs_s, 1.0, 0.0, np.array([[0.0, 10.0, 0.0], [0.0, -10.0, 0.0]])).max()
        )
        u_right = float(box_s.evaluate(coeffs_s.reshape(-1, 3), np.array([[100.0, 0.0, 0.0]]))[0, 0])
        kappa_s = _compute_condition_number(K_s, symmetric=False)
        rows.append(
            {
                "method": "Kratos_SBM",
                "nc": str(nc2),
                "dofs": int(K_s.shape[0]),
                "peak_A": peak_A,
                "err_roark": abs(peak_A - ROARK_PEAK) / ROARK_PEAK,
                "err_abq": abs(peak_A - ROARK_ABAQUS) / ROARK_ABAQUS,
                "u_right": u_right,
                "cond": kappa_s,
                "t_s": dt_s,
            }
        )

    print("\n[P6 V02 Roark] method | nc | dofs | peak_A | err_roark | err_abq | u_right | cond | t_s", flush=True)
    for r in rows:
        print(
            f"[P6 V02 Roark] {r['method']:10s} | {r['nc']:11s} | {r['dofs']:4d} | "
            f"{r['peak_A']:.6e} | {r['err_roark']:.2%} | {r['err_abq']:.2%} | "
            f"{r['u_right']:.6e} | {r['cond']:.3e} | {r['t_s']:.2f}",
            flush=True,
        )

    for r in rows:
        assert np.isfinite(r["peak_A"]) and r["peak_A"] > 0.0
        assert np.isfinite(r["u_right"]) and r["u_right"] > 0.0
