"""Layer 1: immutable domain specs (Pydantic v2, ``frozen=True``; plan v2.1 §8.1, D18).

Importing this package needs only ``pydantic`` and ``numpy``; Kratos, OCP and
the geometry kernels are not loaded.
"""

from __future__ import annotations

from wht_iga.spec.assembly import (
    DEFAULT_GRID_FACE_SETS,
    BoxGeometrySpec,
    GeometrySpec,
    GridFaceSetSpec,
    PartSpec,
    SetSpec,
    StepGeometrySpec,
)
from wht_iga.spec.base import SpecBase
from wht_iga.spec.mesh import (
    BodyFitSpec,
    BoxRegionSpec,
    CadEntitiesRegionSpec,
    ImmersedGridSpec,
    MeshSpec,
    RefinementSpec,
    RegionSpec,
    SphereRegionSpec,
)
from wht_iga.spec.physics import (
    BCSpec,
    DisplacementBCSpec,
    GravityLoadSpec,
    LinearElasticSpec,
    LoadSpec,
    MaterialSpec,
    PressureLoadSpec,
    StaticStepSpec,
    StepSpec,
    TractionLoadSpec,
)
from wht_iga.spec.study import (
    PROBE_VARIABLE_CANONICAL,
    LineProbeSpec,
    MultiPointProbeSpec,
    PointProbeSpec,
    ProbeSpec,
    StudySpec,
)
from wht_iga.spec.units import UNIT_SYSTEMS, UnitSystem, check_material_units

__all__ = [
    "BCSpec", "BodyFitSpec", "BoxGeometrySpec", "BoxRegionSpec", "CadEntitiesRegionSpec",
    "DEFAULT_GRID_FACE_SETS", "DisplacementBCSpec", "GeometrySpec", "GravityLoadSpec",
    "GridFaceSetSpec", "ImmersedGridSpec", "LineProbeSpec", "LinearElasticSpec", "LoadSpec",
    "MaterialSpec", "MeshSpec", "MultiPointProbeSpec", "PROBE_VARIABLE_CANONICAL", "PartSpec",
    "PointProbeSpec", "PressureLoadSpec", "ProbeSpec", "RefinementSpec", "RegionSpec", "SetSpec",
    "SpecBase", "SphereRegionSpec", "StaticStepSpec", "StepGeometrySpec", "StepSpec", "StudySpec",
    "TractionLoadSpec", "UNIT_SYSTEMS", "UnitSystem", "check_material_units",
]
