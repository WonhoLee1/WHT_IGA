"""S11: steady heat conduction (IgaApp LaplacianElement + SupportLaplacianCondition).

Question (v2 §5.8): do ``LaplacianElement`` + ``SupportLaplacianCondition``
(3D), BDF transient, and slave-surface ``ThermalFace3D3N`` film/radiation work?

References: 1D analytic profiles (steady Dirichlet, steady film, steady
radiation solved by scalar Newton in-test) and a geometry-based consistent
mass check exposing the transient limitation (mass matrix carries
conductivity, no rho*c_p: ``laplacian_element.cpp`` l.158-166).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.kratos_backend.thermal import (
    add_temperature_dofs,
    convection_diffusion_settings,
    ensure_thermal_variables,
    make_conduction_properties,
    make_thermal_face_properties,
)
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos

_THERMAL_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
    "ConvectionDiffusionApplication",
)


def _load():  # type: ignore[no-untyped-def]
    from wht_iga.kratos_backend import load_kratos

    return load_kratos(_THERMAL_APPS)


def _bar_model(KM, SM, name: str):  # type: ignore[no-untyped-def]
    model = KM.Model()
    mp = model.CreateModelPart(name)
    ensure_thermal_variables(mp, KM)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings(KM))
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 0.1, 0.1), n_cells=(2, 1, 1), degree=1)
    volume, _ = make_nurbs_volume(KM, mp, box, add_dofs=False)
    add_temperature_dofs(KM, mp)
    props = make_conduction_properties(mp, KM, 1.0)
    eid = 1
    for i in range(2):
        qv = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(qv, 2, cell_gauss_points(box, i, 0, 0))
        mp.CreateNewElement("LaplacianElement", eid, qv[0], props)
        eid += 1
    return mp, box, volume


def _solve_thermal(KM, mp) -> None:  # type: ignore[no-untyped-def]
    mp.SetBufferSize(2)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    crit = KM.ResidualCriteria(1e-9, 1e-12)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(mp, scheme, crit, builder, 15, True, False, False)
    strategy.SetEchoLevel(0)
    mp.CloneTimeStep(1.0)
    strategy.Solve()


def _face_triangles(x: float) -> list[np.ndarray]:
    y0, y1, z0, z1 = 0.0, 0.1, 0.0, 0.1
    if x == 0.0:
        return [np.array([[x, y0, z0], [x, y0, z1], [x, y1, z1]]),
                np.array([[x, y0, z0], [x, y1, z1], [x, y1, z0]])]
    return [np.array([[x, y0, z0], [x, y1, z0], [x, y1, z1]]),
            np.array([[x, y0, z0], [x, y1, z1], [x, y0, z1]])]


def _carrier_grid(x: float, n: int) -> list[np.ndarray]:
    """n x n subgrid of the x-face, each quad split into 2 triangles (2n^2)."""
    ys = np.linspace(0.0, 0.1, n + 1)
    zs = np.linspace(0.0, 0.1, n + 1)
    tris = []
    for j in range(n):
        for k in range(n):
            p00 = np.array([x, ys[j], zs[k]])
            p10 = np.array([x, ys[j + 1], zs[k]])
            p11 = np.array([x, ys[j + 1], zs[k + 1]])
            p01 = np.array([x, ys[j], zs[k + 1]])
            tris.append(np.array([p00, p10, p11]))
            tris.append(np.array([p00, p11, p01]))
    return tris


def _tie_carriers(KM, mp, box, carriers, first_dof_id: int = 1) -> None:  # type: ignore[no-untyped-def]
    """Tie carrier TEMPERATURE dofs to the spline field (general weights)."""
    from wht_iga.grid.bspline import basis_1d

    nu, nv, _ = box.n_cp
    for node in carriers:
        node.AddDof(KM.TEMPERATURE, KM.REACTION_FLUX)
    for node in carriers:
        weights: dict[int, float] = {}
        spans = []
        for d, xd in enumerate((node.X0, node.Y0, node.Z0)):
            span, ders = basis_1d(box.knots[d], box.degree, xd, 0)
            spans.append((span - box.degree, ders[0]))
        (ox, wx), (oy, wy), (oz, wz) = spans
        for c in range(box.degree + 1):
            for b in range(box.degree + 1):
                for a in range(box.degree + 1):
                    I = (ox + a) + nu * ((oy + b) + nv * (oz + c))
                    weights[I] = wx[a] * wy[b] * wz[c]
        assert abs(sum(weights.values()) - 1.0) < 1e-12
        masters = sorted(weights)
        master_dofs = [next(n for n in mp.Nodes if n.Id == idx + 1).GetDof(KM.TEMPERATURE) for idx in masters]
        rel = KM.Matrix(1, len(masters))
        for j, idx in enumerate(masters):
            rel[0, j] = weights[idx]
        const = KM.Vector(1)
        const[0] = 0.0
        cid = mp.NumberOfMasterSlaveConstraints() + 1
        mp.AddMasterSlaveConstraint(
            KM.LinearMasterSlaveConstraint(cid, master_dofs, [node.GetDof(KM.TEMPERATURE)], rel, const)
        )


def test_support_laplacian_condition_raises() -> None:
    """Negative result (deterministic): SupportLaplacianCondition requires
    ``Normal(0, method)``, which per-point QP surface geometry does not define
    (``[QUADRATURE_POINT_SURFACE_IN_VOLUME_GEOMETRY] Normal not defined``).
    On the node-less wrapper geometry it segfaults instead (evidence in S11.md,
    not executed here). Either way the condition cannot back thermal Dirichlet
    BCs on immersed Python-built geometry in this build."""
    mods = _load()
    KM = mods["KratosMultiphysics"]
    IGA = mods["IgaApplication"]
    mp, box, volume = _bar_model(KM, None, "S11n")
    qp_geoms, wrappers, _, _ = make_surface_quadrature(
        KM, mp, volume, _face_triangles(0.0), rule_order=1, first_node_id=10_000
    )
    props = mp.CreateNewProperties(2)
    props.SetValue(IGA.PENALTY_FACTOR, 1.0e10)
    cond = mp.CreateNewCondition("SupportLaplacianCondition", 1, qp_geoms[0], props)
    cond.SetValue(KM.TEMPERATURE, 0.0)
    with pytest.raises(RuntimeError, match="Normal not defined"):
        cond.CalculateLocalSystem(KM.Matrix(), KM.Vector(), mp.ProcessInfo)
    _ = (box, wrappers)


def _fix_face_temperature(KM, mp, x: float, value: float) -> None:  # type: ignore[no-untyped-def]
    for node in mp.Nodes:
        if abs(node.X0 - x) < 1e-12:
            node.Fix(KM.TEMPERATURE)
            node.SetSolutionStepValue(KM.TEMPERATURE, value)


def test_steady_dirichlet_linear_profile() -> None:
    mods = _load()
    KM = mods["KratosMultiphysics"]
    IGA = mods["IgaApplication"]
    mp, box, volume = _bar_model(KM, None, "S11d")
    # SupportLaplacianCondition is unusable on immersed QP geometry (see
    # test_support_laplacian_condition_raises); impose T strongly on faces.
    for node in mp.Nodes:
        if abs(node.X0) < 1e-12:
            node.Fix(KM.TEMPERATURE)
            node.SetSolutionStepValue(KM.TEMPERATURE, 0.0)
        elif abs(node.X0 - 1.0) < 1e-12:
            node.Fix(KM.TEMPERATURE)
            node.SetSolutionStepValue(KM.TEMPERATURE, 100.0)
    _solve_thermal(KM, mp)
    for node in mp.Nodes:
        if abs(node.X0 - 0.5) < 1e-12:
            t = node.GetSolutionStepValue(KM.TEMPERATURE)
            print(f"\nS11 Dirichlet: T(x=0.5) = {t:.10f} (expect 50.0)")
            np.testing.assert_allclose(t, 50.0, rtol=0.0, atol=1e-6)
            break
    _ = (volume, box, IGA)


def _film_run(KM, n: int):  # type: ignore[no-untyped-def]
    """BC-C film model with n x n carrier subgrid; returns (avg_T1, q_cond, q_film)."""
    from wht_iga.kratos_backend.surface import triangle_area

    mp, box, volume = _bar_model(KM, None, f"S11f{n}")
    _fix_face_temperature(KM, mp, 0.0, 0.0)
    tris = _carrier_grid(1.0, n)
    _, wrappers, _, _ = make_surface_quadrature(KM, mp, volume, tris, rule_order=1, first_node_id=20_000)
    carriers = sorted([nd for nd in mp.Nodes if nd.Id >= 20_000], key=lambda nd: nd.Id)
    assert len(carriers) == 6 * n * n  # 3 coincident carriers per triangle
    _tie_carriers(KM, mp, box, carriers)
    props = make_thermal_face_properties(mp, KM, convection_coeff=2.0, ambient_temperature=100.0, prop_id=5)
    # rebuild triangles in carrier order: _carrier_grid emits 2 quads-split tris
    # per quad, vertices shared via coincident carriers; use index triples:
    ys = np.linspace(0.0, 0.1, n + 1)
    idx = {(1.0, round(float(nd.Y0), 12), round(float(nd.Z0), 12)): nd for nd in carriers}
    areas = []
    for tri in tris:
        tnodes = [idx[(1.0, round(float(v[1]), 12), round(float(v[2]), 12))] for v in tri]
        mp.CreateNewCondition(
            "ThermalFace3D3N", mp.NumberOfConditions() + 1,
            KM.Triangle3D3(tnodes[0], tnodes[1], tnodes[2]), props,
        )
        areas.append(triangle_area(np.array([[nd.X0, nd.Y0, nd.Z0] for nd in tnodes])))
    _solve_thermal(KM, mp)
    t1 = [nd.GetSolutionStepValue(KM.TEMPERATURE) for nd in mp.Nodes if abs(nd.X0 - 1.0) < 1e-12 and nd.Id < 20_000]
    avg = float(sum(t1) / len(t1))
    q_cond = float(sum(nd.GetSolutionStepValue(KM.REACTION_FLUX) for nd in mp.Nodes if abs(nd.X0) < 1e-12))
    temps = {nd.Id: nd.GetSolutionStepValue(KM.TEMPERATURE) for nd in carriers}
    q_film = 0.0
    for tri, a in zip(tris, areas):
        tbar = float(sum(temps[idx[(1.0, round(float(v[1]), 12), round(float(v[2]), 12))].Id] for v in tri)) / 3.0
        q_film += 2.0 * (tbar - 100.0) * a  # linear-in-triangle T: centroid rule exact
    _ = (wrappers, box, volume)
    return avg, q_cond, q_film


def test_steady_film_robin() -> None:
    """BC-C film: energy balance (exact) + analytic accuracy + refinement trend.

    Triangulated film loads diagonal corners twice (valence 2 vs 1), so the
    discrete average carries an O(h) non-conforming-coupling error: n=1 gives
    ~0.2%, n=2 must improve on n=1 while both conserve energy to 1e-6.
    """
    mods = _load()
    KM = mods["KratosMultiphysics"]
    avg1, qc1, qf1 = _film_run(KM, 1)
    avg2, qc2, qf2 = _film_run(KM, 2)
    for tag, avg, qc, qf in (("n=1", avg1, qc1, qf1), ("n=2", avg2, qc2, qf2)):
        print(f"\nS11 film {tag}: avg T(x=1) = {avg:.8f} (expect 66.666...), "
              f"Q_cond = {qc:.6e}, Q_film = {qf:.6e}")
    np.testing.assert_allclose(abs(qc1), abs(qf1), rtol=1e-6, atol=0.0)
    np.testing.assert_allclose(abs(qc2), abs(qf2), rtol=1e-6, atol=0.0)
    np.testing.assert_allclose(avg1, 200.0 / 3.0, rtol=0.01, atol=0.0)
    np.testing.assert_allclose(avg2, 200.0 / 3.0, rtol=0.01, atol=0.0)
    assert abs(avg2 - 200.0 / 3.0) < abs(avg1 - 200.0 / 3.0)


def test_radiation_face_converges_to_scalar_root() -> None:
    """k=1, T(0)=0, radiation eps=0.5/Tinf=100 covering x=1: kC = eps*sig*(Tinf^4-C^4).

    Same BC-C carrier setup as the film test; Newton loop handles the quartic.
    """
    from wht_iga.grid.bspline import basis_1d

    mods = _load()
    KM = mods["KratosMultiphysics"]
    mp, box, volume = _bar_model(KM, None, "S11r")
    _fix_face_temperature(KM, mp, 0.0, 0.0)
    _, wrappers, _, _ = make_surface_quadrature(
        KM, mp, volume, _face_triangles(1.0), rule_order=1, first_node_id=20_000
    )
    carriers = sorted([n for n in mp.Nodes if n.Id >= 20_000], key=lambda n: n.Id)
    assert len(carriers) == 6
    _tie_carriers(KM, mp, box, carriers)
    props = make_thermal_face_properties(mp, KM, ambient_temperature=100.0, emissivity=0.5, prop_id=5)
    for t in range(2):
        tri = KM.Triangle3D3(carriers[3 * t], carriers[3 * t + 1], carriers[3 * t + 2])
        mp.CreateNewCondition("ThermalFace3D3N", mp.NumberOfConditions() + 1, tri, props)
    _solve_thermal(KM, mp)
    # Stefan-Boltzmann constant used by ThermalFace (thermal_face.h:92).
    sig = 5.67e-8
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - 0.5 * sig * (100.0**4 - mid**4) > 0.0:
            hi = mid
        else:
            lo = mid
    cref = 0.5 * (lo + hi)
    t1 = [node.GetSolutionStepValue(KM.TEMPERATURE) for node in mp.Nodes if abs(node.X0 - 1.0) < 1e-12 and node.Id < 20_000]
    avg = float(sum(t1) / len(t1))
    print(f"\nS11 radiation: T(x=1) = {avg:.10f} (scalar root {cref:.10f}, sig={sig:.4e})")
    np.testing.assert_allclose(avg, cref, rtol=0.01, atol=0.0)
    _ = (wrappers, box)


def test_mass_matrix_has_no_heat_capacity() -> None:
    """Transient limitation: M == k * N^T N (no rho*c_p), so BDF heat-capacity
    transients are unphysical through LaplacianElement. Evidence, not a pass."""
    mods = _load()
    KM = mods["KratosMultiphysics"]
    mp, box, volume = _bar_model(KM, None, "S11m")
    el = next(iter(mp.Elements))
    el.Initialize(mp.ProcessInfo)
    M = KM.Matrix()
    el.CalculateMassMatrix(M, mp.ProcessInfo)
    geom = el.GetGeometry()
    N = np.array(geom.ShapeFunctionsValues())
    npts = geom.IntegrationPointsNumber()
    # LaplacianElement::CalculateLocalSystem/MassMatrix multiplies the stored
    # Weight() only (no extra detJ): for QP geometries the stored weight is the
    # physical measure (S08 totals are exact). Supplied weights here are uniform
    # physical dV = 0.5 * 0.1 * 0.1 / 8 per point.
    wval = 0.5 * 0.1 * 0.1 / 8.0
    ref = (N.T @ N) * wval  # conductivity k = 1
    got = np.array([[M[i, j] for j in range(M.Size2())] for i in range(M.Size1())])
    print(f"\nS11 mass: M[0,0] = {got[0, 0]:.6e}, k*NtN[0,0] = {ref[0, 0]:.6e}")
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=0.0)
    _ = (box, volume)
