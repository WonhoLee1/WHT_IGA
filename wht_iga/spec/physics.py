"""Physics specs: materials, boundary conditions, loads and analysis steps.

Sign conventions are those of the solver pipeline (``wht_iga.solver.loads``):
pressure gives traction ``+magnitude * n_out`` on the selected face (S08), a
traction vector is applied as given, gravity is an acceleration multiplied by
the material density.
"""

from __future__ import annotations

from typing import Annotated, Literal, Self, Union

from pydantic import Field, model_validator

from wht_iga.spec.base import Finite, Name, PositiveFinite, SpecBase, Vec3


class MaterialSpec(SpecBase):
    """Base class of material specs (discriminated by ``kind``)."""


class LinearElasticSpec(MaterialSpec):
    """Isotropic linear elastic material."""

    kind: Literal["linear_elastic"] = "linear_elastic"
    E: PositiveFinite
    nu: Annotated[float, Field(gt=-1.0, lt=0.5)]
    density: Annotated[float, Field(ge=0.0, allow_inf_nan=False)] = 0.0


MaterialUnion = LinearElasticSpec


class DisplacementBCSpec(SpecBase):
    """Prescribed displacement components on a named set (``None`` = free)."""

    kind: Literal["displacement"] = "displacement"
    set: Name
    ux: Finite | None = None
    uy: Finite | None = None
    uz: Finite | None = None

    @model_validator(mode="after")
    def _at_least_one(self) -> Self:
        if self.ux is None and self.uy is None and self.uz is None:
            raise ValueError(f"displacement BC on {self.set!r} constrains no component (give ux/uy/uz)")
        return self

    @property
    def components(self) -> tuple[bool, bool, bool]:
        """Which components are constrained."""
        return (self.ux is not None, self.uy is not None, self.uz is not None)

    @property
    def value(self) -> tuple[float, float, float]:
        """Prescribed values (free components reported as 0.0)."""
        return (self.ux or 0.0, self.uy or 0.0, self.uz or 0.0)


BCSpec = DisplacementBCSpec


class PressureLoadSpec(SpecBase):
    """Uniform pressure on a named set: traction ``+magnitude * n_out`` (repo convention)."""

    kind: Literal["pressure"] = "pressure"
    set: Name
    magnitude: Finite
    follower: bool = False


class TractionLoadSpec(SpecBase):
    """Uniform traction vector (force / area) on a named set."""

    kind: Literal["traction"] = "traction"
    set: Name
    vector: Vec3


class GravityLoadSpec(SpecBase):
    """Body acceleration; body force = density * acceleration."""

    kind: Literal["gravity"] = "gravity"
    acceleration: Vec3


LoadSpec = Annotated[Union[PressureLoadSpec, TractionLoadSpec, GravityLoadSpec],
                     Field(discriminator="kind")]

LinearSolverName = Literal["auto", "pardiso", "pardiso_ldlt", "pardiso_lu", "pardiso_llt",
                           "amgcl", "amgcl_cg", "cg", "splu", "skyline_lu"]


class StaticStepSpec(SpecBase):
    """Linear static step (independent solve; steps do not accumulate)."""

    kind: Literal["static"] = "static"
    name: Name
    bcs: tuple[BCSpec, ...] = Field(min_length=1)
    loads: tuple[LoadSpec, ...] = ()
    linear_solver: LinearSolverName = "auto"
    memory_budget_ratio: Annotated[float, Field(gt=0.0, le=1.0)] = 0.70
    num_threads_mkl: Annotated[int, Field(ge=1)] | None = None

    @model_validator(mode="after")
    def _check_loads(self) -> Self:
        n_grav = sum(1 for ld in self.loads if isinstance(ld, GravityLoadSpec))
        if n_grav > 1:
            raise ValueError(f"step {self.name!r}: at most one gravity load, got {n_grav}")
        return self


#: step union (single member today; Modal/Implicit/Explicit/Thermal steps extend it)
StepSpec = StaticStepSpec
