"""Mesh specs: immersed rectilinear grid, body-fit grid, refinement and regions.

Field semantics follow the Coreform portable mesh spec (``ImmersedRectilinearMeshSpec``:
``continuity < degree``, ``small_cell_volume_ratio`` in ``[0, 0.5)``, nonnegative
``padding``), renamed to WHT_IGA vocabulary. Sizing is either a uniform
``element_size`` or per-axis ``n_cells`` (exactly one).
"""

from __future__ import annotations

from typing import Annotated, Literal, Self, Union

from pydantic import Field, model_validator

from wht_iga.spec.base import Name, PositiveFinite, SpecBase, Vec3

Degree = Annotated[int, Field(ge=1, le=3)]


class CadEntitiesRegionSpec(SpecBase):
    """Region given by CAD entity pointers (e.g. ``"@face:f003"``)."""

    kind: Literal["cad_entities"] = "cad_entities"
    entities: tuple[Annotated[str, Field(min_length=1)], ...] = Field(min_length=1)


class SphereRegionSpec(SpecBase):
    """Ball region (center, radius)."""

    kind: Literal["sphere"] = "sphere"
    center: Vec3
    radius: PositiveFinite


class BoxRegionSpec(SpecBase):
    """Axis-aligned box region (``hi`` strictly above ``lo`` on every axis)."""

    kind: Literal["box"] = "box"
    lo: Vec3
    hi: Vec3

    @model_validator(mode="after")
    def _check_box(self) -> Self:
        if any(h <= l for h, l in zip(self.hi, self.lo)):
            raise ValueError(f"box region hi must exceed lo on every axis, got lo={self.lo} hi={self.hi}")
        return self


RegionSpec = Annotated[Union[CadEntitiesRegionSpec, SphereRegionSpec, BoxRegionSpec],
                       Field(discriminator="kind")]


class RefinementSpec(SpecBase):
    """Local (THB) refinement request on named sets or on a region (exactly one)."""

    name: Name
    level: Annotated[int, Field(ge=1, le=8)]
    sets: tuple[Name, ...] = ()
    region: RegionSpec | None = None
    balance: Literal["two_to_one", "none"] = "two_to_one"

    @model_validator(mode="after")
    def _one_target(self) -> Self:
        if bool(self.sets) == (self.region is not None):
            raise ValueError(f"refinement {self.name!r}: give exactly one of sets=... or region=...")
        return self


class _GridCommon(SpecBase):
    """Fields shared by immersed and body-fit grids."""

    degree: Degree = 2
    continuity: Annotated[int, Field(ge=0)] | None = None  # None -> degree - 1 (max smoothness)
    element_size: PositiveFinite | None = None
    n_cells: tuple[Annotated[int, Field(ge=1)], Annotated[int, Field(ge=1)],
                   Annotated[int, Field(ge=1)]] | None = None
    bc_penalty: PositiveFinite | None = None  # absolute penalty beta; None -> c_pen*E*p^2/h
    refinements: tuple[RefinementSpec, ...] = ()

    @model_validator(mode="after")
    def _check_grid(self) -> Self:
        if self.continuity is not None and self.continuity >= self.degree:
            raise ValueError(f"continuity must be < degree, got continuity={self.continuity}, "
                             f"degree={self.degree}")
        if (self.element_size is None) == (self.n_cells is None):
            raise ValueError("give exactly one of element_size=... or n_cells=(nx, ny, nz)")
        names = [r.name for r in self.refinements]
        if len(names) != len(set(names)):
            raise ValueError(f"refinement names must be unique, got {names}")
        return self

    @property
    def resolved_continuity(self) -> int:
        """Continuity actually used (``degree - 1`` when not given)."""
        return self.degree - 1 if self.continuity is None else self.continuity


class ImmersedGridSpec(_GridCommon):
    """Immersed (trimmed) rectilinear B-spline grid around the part geometry."""

    kind: Literal["immersed"] = "immersed"
    tessellation_size: PositiveFinite | None = None  # None -> min(h)/4 (recorded as assumed)
    small_cell_volume_ratio: Annotated[float, Field(ge=0.0, lt=0.5)] = 0.2  # 0 disables extension
    padding: tuple[Annotated[int, Field(ge=0)], Annotated[int, Field(ge=0)],
                   Annotated[int, Field(ge=0)]] = (0, 0, 0)
    bounds: tuple[Vec3, Vec3] | None = None  # explicit grid box (before padding)
    quadrature_degree: Annotated[int, Field(ge=1, le=6)] | None = None  # None -> min(2p, 4)
    out_mode: Literal["drop", "ghost"] = "drop"
    ghost_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = 1e-8

    @model_validator(mode="after")
    def _check_bounds(self) -> Self:
        if self.bounds is not None and any(h <= l for h, l in zip(self.bounds[1], self.bounds[0])):
            raise ValueError(f"bounds hi must exceed lo on every axis, got {self.bounds}")
        return self

    @property
    def resolved_quadrature_degree(self) -> int:
        """Cut-cell moment degree: ``min(2p, 4)`` by default (Strang lemma, V06 evidence)."""
        return min(2 * self.degree, 4) if self.quadrature_degree is None else self.quadrature_degree


class BodyFitSpec(_GridCommon):
    """Body-fitted B-spline grid: the grid box IS the (box-shaped) part, no trimming."""

    kind: Literal["body_fit"] = "body_fit"


MeshSpec = Annotated[Union[ImmersedGridSpec, BodyFitSpec], Field(discriminator="kind")]
