"""Top-level study spec and probe specs.

:class:`StudySpec` is the single immutable description of an analysis. Its
validator checks every cross reference at construction (part -> material,
BC/load/refinement -> set, unique names) and the material numbers against the
declared unit system (:mod:`wht_iga.spec.units`).
"""

from __future__ import annotations

from typing import Annotated, Literal, Self, Union

import numpy as np
from pydantic import Field, field_validator, model_validator

from wht_iga.spec.assembly import DEFAULT_GRID_FACE_SETS, PartSpec, SetSpec
from wht_iga.spec.base import Name, SpecBase, Vec3
from wht_iga.spec.mesh import MeshSpec
from wht_iga.spec.physics import MaterialUnion, StepSpec
from wht_iga.spec.units import UnitSystemName, check_material_units

#: accepted probe variable names (Coreform aliases and Kratos names) -> canonical name
PROBE_VARIABLE_CANONICAL: dict[str, str] = {
    "U": "U", "DISPLACEMENT": "U",
    "S": "S", "CAUCHY_STRESS_TENSOR": "S",
    "E": "E", "STRAIN": "E", "GREEN_LAGRANGE_STRAIN_TENSOR": "E",
    "MISES": "MISES", "VON_MISES_STRESS": "MISES",
    "NT11": "NT11", "TEMPERATURE": "NT11",
    "RF": "RF", "REACTION": "RF",
}


def _check_variables(values: tuple[str, ...]) -> tuple[str, ...]:
    """Reject unknown probe variable names (case-sensitive, like Coreform/Kratos)."""
    bad = [v for v in values if v not in PROBE_VARIABLE_CANONICAL]
    if bad:
        raise ValueError(f"unknown probe variable(s) {bad}; use one of {sorted(PROBE_VARIABLE_CANONICAL)}")
    if len(set(values)) != len(values):
        raise ValueError(f"duplicate probe variables {list(values)}")
    return values


class _ProbeBase(SpecBase):
    """Fields shared by all probes."""

    name: Name
    variables: tuple[str, ...] = ("U",)

    @field_validator("variables")
    @classmethod
    def _vars(cls, v: tuple[str, ...]) -> tuple[str, ...]:
        return _check_variables(v)

    def points(self) -> np.ndarray:
        """Evaluation points as an ``(n, 3)`` array."""
        raise NotImplementedError  # every concrete probe overrides this


class PointProbeSpec(_ProbeBase):
    """Single-point probe."""

    kind: Literal["point"] = "point"
    location: Vec3

    def points(self) -> np.ndarray:
        """One evaluation point."""
        return np.asarray([self.location], dtype=float)


class LineProbeSpec(_ProbeBase):
    """Uniformly sampled segment probe (endpoints included)."""

    kind: Literal["line"] = "line"
    start: Vec3
    end: Vec3
    n_points: Annotated[int, Field(ge=2, le=100_000)] = 11

    @model_validator(mode="after")
    def _distinct(self) -> Self:
        if tuple(self.start) == tuple(self.end):
            raise ValueError(f"line probe {self.name!r}: start and end must differ")
        return self

    def points(self) -> np.ndarray:
        """``n_points`` evenly spaced points from start to end."""
        t = np.linspace(0.0, 1.0, self.n_points)[:, None]
        a, b = np.asarray(self.start, float), np.asarray(self.end, float)
        return a + t * (b - a)


class MultiPointProbeSpec(_ProbeBase):
    """Explicit list of evaluation points."""

    kind: Literal["multi_point"] = "multi_point"
    locations: tuple[Vec3, ...] = Field(min_length=1)

    def points(self) -> np.ndarray:
        """All listed points."""
        return np.asarray(self.locations, dtype=float)


ProbeSpec = Annotated[Union[PointProbeSpec, LineProbeSpec, MultiPointProbeSpec],
                      Field(discriminator="kind")]


class StudySpec(SpecBase):
    """Complete immutable analysis description (plan v2.1 §8.2)."""

    name: Name
    units: UnitSystemName = "mm-N-tonne-s"
    description: str = ""
    parts: tuple[PartSpec, ...] = Field(min_length=1)
    materials: dict[Name, MaterialUnion] = Field(min_length=1)
    mesh: MeshSpec
    steps: tuple[StepSpec, ...] = Field(min_length=1)
    probes: tuple[ProbeSpec, ...] = ()
    sets: dict[Name, SetSpec] = Field(default_factory=dict)
    output_dir: str | None = None

    @model_validator(mode="after")
    def _cross_references(self) -> Self:
        errors: list[str] = []

        def dupes(names: list[str]) -> list[str]:
            return sorted({n for n in names if names.count(n) > 1})

        for label, names in (("part", [p.name for p in self.parts]),
                             ("step", [s.name for s in self.steps]),
                             ("probe", [p.name for p in self.probes])):
            if dupes(names):
                errors.append(f"duplicate {label} names {dupes(names)}")
        for p in self.parts:
            if p.material not in self.materials:
                errors.append(f"part {p.name!r} references unknown material {p.material!r} "
                              f"(defined: {sorted(self.materials)})")
        shadow = sorted(set(self.sets) & set(DEFAULT_GRID_FACE_SETS))
        if shadow:
            errors.append(f"set names {shadow} are reserved default grid-face sets")
        known = set(self.sets) | set(DEFAULT_GRID_FACE_SETS)
        for s in self.steps:
            for bc in s.bcs:
                if bc.set not in known:
                    errors.append(f"step {s.name!r}: BC references unknown set {bc.set!r}")
            for ld in s.loads:
                tgt = getattr(ld, "set", None)
                if tgt is not None and tgt not in known:
                    errors.append(f"step {s.name!r}: {ld.kind} load references unknown set {tgt!r}")
        for r in self.mesh.refinements:
            for name in r.sets:
                if name not in known:
                    errors.append(f"refinement {r.name!r} references unknown set {name!r}")
        for mname, mat in self.materials.items():
            errors.extend(check_material_units(self.units, mat.E, mat.density, label=f"material {mname!r}"))
        if errors:
            raise ValueError("; ".join(errors))
        return self

    def material_of(self, part: PartSpec) -> MaterialUnion:
        """Material spec of a part."""
        return self.materials[part.material]

    def set_face(self, name: str) -> tuple[int, int]:
        """Grid face ``(axis, side)`` of a named set (user sets, then defaults)."""
        if name in self.sets:
            s = self.sets[name]
            return (int(s.axis), int(s.side))
        return DEFAULT_GRID_FACE_SETS[name]

    def step(self, name: str) -> StepSpec:
        """Step by name (``KeyError`` if absent)."""
        for s in self.steps:
            if s.name == name:
                return s
        raise KeyError(f"no step named {name!r}; steps: {[s.name for s in self.steps]}")
