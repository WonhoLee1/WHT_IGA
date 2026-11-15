"""Part / geometry / set specs (plan v2.1 §5.9 vocabulary, single-part subset).

Geometry sources:

* :class:`BoxGeometrySpec` -- analytic axis-aligned box (body-fit or immersed).
* :class:`StepGeometrySpec` -- STEP file read through ``wht_iga.geometry.cad``.

Named sets map to background-grid faces (:class:`GridFaceSetSpec`). Six default
names always exist: ``xmin, xmax, ymin, ymax, zmin, zmax``. Boundary conditions
and loads in the current solver pipeline act on grid faces only
(``wht_iga.solver.penalty`` / ``loads``); CAD-face sets are rejected by
:func:`wht_iga.plan.compile_plan` with a diagnostic instead of being ignored.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, Self, Union

from pydantic import Field, model_validator

from wht_iga.spec.base import Name, SpecBase, Vec3


class BoxGeometrySpec(SpecBase):
    """Analytic axis-aligned box ``[lo, hi]``."""

    kind: Literal["box"] = "box"
    lo: Vec3
    hi: Vec3

    @model_validator(mode="after")
    def _check_box(self) -> Self:
        if any(h <= l for h, l in zip(self.hi, self.lo)):
            raise ValueError(f"box geometry hi must exceed lo on every axis, got lo={self.lo} hi={self.hi}")
        return self


class StepGeometrySpec(SpecBase):
    """STEP file geometry (read and tessellated at ``build_mesh`` time)."""

    kind: Literal["step"] = "step"
    path: Annotated[str, Field(min_length=1)]
    heal: bool = True


GeometrySpec = Annotated[Union[BoxGeometrySpec, StepGeometrySpec], Field(discriminator="kind")]


class PartSpec(SpecBase):
    """One solid part: geometry source + material name."""

    name: Name
    geometry: GeometrySpec
    material: Name

    @classmethod
    def from_step(cls, path: str | Path, material: str, *, heal: bool = True,
                  name: str | None = None) -> PartSpec:
        """Part from a STEP file; ``name`` defaults to the file stem."""
        p = Path(path)
        return cls(name=name or p.stem, geometry=StepGeometrySpec(path=str(p), heal=heal),
                   material=material)

    @classmethod
    def from_box(cls, lo: Vec3, hi: Vec3, material: str, *, name: str = "box") -> PartSpec:
        """Part from an analytic axis-aligned box."""
        return cls(name=name, geometry=BoxGeometrySpec(lo=lo, hi=hi), material=material)


class GridFaceSetSpec(SpecBase):
    """Named set = one background-grid face (``axis`` 0/1/2, ``side`` 0 = lo, 1 = hi)."""

    kind: Literal["grid_face"] = "grid_face"
    axis: Annotated[int, Field(ge=0, le=2)]
    side: Literal[0, 1]


SetSpec = GridFaceSetSpec

#: always-available set names -> (axis, side)
DEFAULT_GRID_FACE_SETS: dict[str, tuple[int, int]] = {
    "xmin": (0, 0), "xmax": (0, 1),
    "ymin": (1, 0), "ymax": (1, 1),
    "zmin": (2, 0), "zmax": (2, 1),
}
