"""Thermal helpers for Track B spikes S11 (conduction) and S12 (expansion).

Ownership: Track B file. May be promoted to shared use by the integrator.
"""

from __future__ import annotations


def ensure_thermal_variables(model_part, KM) -> None:
    """Register nodal solution-step variables needed by thermal spikes."""
    for var in (KM.TEMPERATURE, KM.HEAT_FLUX, KM.CONDUCTIVITY, KM.SPECIFIC_HEAT, KM.DENSITY, KM.REACTION_FLUX):
        model_part.AddNodalSolutionStepVariable(var)
    if hasattr(KM, "FACE_HEAT_FLUX"):
        model_part.AddNodalSolutionStepVariable(KM.FACE_HEAT_FLUX)


def add_temperature_dofs(KM, model_part) -> None:
    """Add the (TEMPERATURE, REACTION_FLUX) DOF pair to every node."""
    KM.VariableUtils().AddDof(KM.TEMPERATURE, KM.REACTION_FLUX, model_part)


def convection_diffusion_settings(KM):
    """ProcessInfo settings selecting T as unknown, k as diffusivity, q as source."""
    settings = KM.ConvectionDiffusionSettings()
    settings.SetUnknownVariable(KM.TEMPERATURE)
    settings.SetDiffusionVariable(KM.CONDUCTIVITY)
    settings.SetVolumeSourceVariable(KM.HEAT_FLUX)
    if hasattr(KM, "TEMPERATURE_GRADIENT"):
        settings.SetGradientVariable(KM.TEMPERATURE_GRADIENT)
    if hasattr(KM, "FACE_HEAT_FLUX"):
        settings.SetSurfaceSourceVariable(KM.FACE_HEAT_FLUX)
    return settings


def make_conduction_properties(model_part, KM, conductivity: float, prop_id: int = 1):
    """Properties carrying CONDUCTIVITY for ``LaplacianElement``."""
    prop = model_part.CreateNewProperties(prop_id)
    prop.SetValue(KM.CONDUCTIVITY, float(conductivity))
    return prop


def make_thermal_face_properties(
    model_part,
    KM,
    convection_coeff: float = 0.0,
    ambient_temperature: float = 0.0,
    emissivity: float = 0.0,
    prop_id: int = 2,
    conductivity: float = 1.0,
):
    """Properties for ConvDiff ``ThermalFace3D3N`` (film and/or radiation).

    ``conductivity`` is stored because the face condition reads it from its
    Properties; pass the solid's value (the 1.0 default is the S11 test value).
    """
    prop = model_part.CreateNewProperties(prop_id)
    prop.SetValue(KM.CONDUCTIVITY, float(conductivity))
    prop.SetValue(KM.CONVECTION_COEFFICIENT, float(convection_coeff))
    prop.SetValue(KM.AMBIENT_TEMPERATURE, float(ambient_temperature))
    prop.SetValue(KM.EMISSIVITY, float(emissivity))
    return prop


def make_thermal_elastic_properties(
    model_part, KM, SM, CLA, E: float, nu: float, alpha: float, T_ref: float, prop_id: int = 1,
    density: float = 1.0,
):
    """Properties with ``ThermalElasticIsotropic3D`` law for S12.

    ``density`` only matters for dynamics/body force; the 1.0 default is the S12 static value.
    """
    prop = model_part.CreateNewProperties(prop_id)
    prop.SetValue(KM.YOUNG_MODULUS, float(E))
    prop.SetValue(KM.POISSON_RATIO, float(nu))
    prop.SetValue(KM.DENSITY, float(density))
    prop.SetValue(KM.THERMAL_EXPANSION_COEFFICIENT, float(alpha))
    prop.SetValue(KM.REFERENCE_TEMPERATURE, float(T_ref))
    prop.SetValue(KM.CONSTITUTIVE_LAW, CLA.ThermalElasticIsotropic3D())
    return prop


def kratos_box_conductivity(box, conductivity: float):
    """Conduction matrix ``int k grad(N)^T grad(N) dV`` of a box B-spline assembled from Kratos elements.

    One ``LaplacianElement`` per cell on a quadrature-point geometry (route R1).  The local
    matrices (``CalculateLocalSystem``) are scattered with the control-point numbering
    ``node id - 1``, so no Dirichlet treatment touches the matrix: the full ``K_th`` (including
    the coupling to later-prescribed temperatures) is returned.  Control point ``I`` carries
    one TEMPERATURE DOF.

    Args:
        box: Box B-spline (untrimmed).
        conductivity: Isotropic thermal conductivity ``k``.

    Returns:
        ``scipy.sparse.csr_matrix`` of shape ``(box.n_ctrl, box.n_ctrl)``.
    """
    import numpy as np
    import scipy.sparse as sp

    from wht_iga.kratos_backend import load_kratos
    from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume

    mods = load_kratos(
        ("IgaApplication", "StructuralMechanicsApplication", "LinearSolversApplication",
         "ConstitutiveLawsApplication", "ConvectionDiffusionApplication")
    )
    KM = mods["KratosMultiphysics"]
    model = KM.Model()
    mp = model.CreateModelPart("P4Thermal")
    ensure_thermal_variables(mp, KM)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    mp.ProcessInfo.SetValue(KM.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings(KM))
    volume, _nodes = make_nurbs_volume(KM, mp, box, add_dofs=False)
    add_temperature_dofs(KM, mp)
    props = make_conduction_properties(mp, KM, conductivity)
    nx, ny, nz = box.n_cells
    rows, cols, vals = [], [], []
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                qv = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(qv, 2, cell_gauss_points(box, i, j, k))
                el = mp.CreateNewElement("LaplacianElement", eid, qv[0], props)
                eid += 1
                el.Initialize(mp.ProcessInfo)
                lhs, rhs = KM.Matrix(0, 0), KM.Vector(0)
                el.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
                local = np.array(lhs)
                ids = np.array([g.Id for g in el.GetGeometry()], dtype=np.int64) - 1
                if local.shape != (ids.size, ids.size):
                    raise RuntimeError(f"unexpected local matrix shape {local.shape} for {ids.size} nodes")
                rows.append(np.repeat(ids, ids.size))
                cols.append(np.tile(ids, ids.size))
                vals.append(local.reshape(-1))
    n = box.n_ctrl
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
