"""Layer 1 base: immutable, validated spec models (Pydantic v2, D18 / plan v2.1 §8.1).

Every user-facing spec derives from :class:`SpecBase`:

* ``frozen=True``: attribute assignment raises; sequences are stored as tuples.
* ``extra="forbid"``: misspelled keywords fail at construction instead of being
  silently ignored.
* :meth:`SpecBase.with_updates` returns a re-validated copy (the original is
  never touched), so parametric sweeps cannot corrupt a base study.
* JSON round trip: ``model_dump_json`` / ``model_validate_json`` plus the
  short aliases :meth:`to_json` / :meth:`from_json` and :meth:`json_schema`.
"""

from __future__ import annotations

from typing import Annotated, Any, Self

from pydantic import BaseModel, ConfigDict, Field

#: finite float (rejects nan / inf)
Finite = Annotated[float, Field(allow_inf_nan=False)]
#: strictly positive finite float
PositiveFinite = Annotated[float, Field(gt=0.0, allow_inf_nan=False)]
#: 3-vector of finite floats (lists are coerced to tuples)
Vec3 = tuple[Finite, Finite, Finite]
#: identifier used for names that end up in file names and set lookups
Name = Annotated[str, Field(min_length=1, max_length=64, pattern=r"^[A-Za-z0-9_][A-Za-z0-9_.\-]*$")]


class SpecBase(BaseModel):
    """Frozen Pydantic v2 base model shared by all WHT_IGA specs."""

    model_config = ConfigDict(frozen=True, extra="forbid", validate_default=True,
                              str_strip_whitespace=True)

    def with_updates(self, **changes: Any) -> Self:
        """Return a validated copy with ``changes`` applied (original unchanged).

        Unknown field names raise ``ValueError``. Cross-field validators run
        again on the merged data, so an update cannot produce an invalid spec.
        """
        fields = type(self).model_fields
        unknown = sorted(set(changes) - set(fields))
        if unknown:
            raise ValueError(f"{type(self).__name__}.with_updates: unknown field(s) {unknown}")
        data = {name: getattr(self, name) for name in fields}
        data.update(changes)
        return type(self).model_validate(data)

    def to_json(self, *, indent: int | None = None) -> str:
        """Serialize to a JSON document (alias of ``model_dump_json``)."""
        return self.model_dump_json(indent=indent)

    @classmethod
    def from_json(cls, text: str | bytes) -> Self:
        """Validate and build an instance from a JSON document."""
        return cls.model_validate_json(text)

    @classmethod
    def json_schema(cls) -> dict[str, Any]:
        """JSON Schema of this spec (alias of ``model_json_schema``) for AI/MCP/UI use."""
        return cls.model_json_schema()
