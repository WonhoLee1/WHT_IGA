"""Selection sets, boolean set algebra, and named set registry (v2.1 §8.5).

Provides :class:`SelectionSet` for node, element, surface, and CAD topology
entity ID sets with full boolean algebra (``|``, ``&``, ``-``, ``^``, ``~``)
and :class:`SetRegistry` for named registration and derived set management.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

SetKind = Literal[
    "node",
    "element",
    "surface_tri",
    "cad_face",
    "cad_edge",
    "cad_vertex",
    "cad_solid",
]

VALID_SET_KINDS: tuple[SetKind, ...] = (
    "node",
    "element",
    "surface_tri",
    "cad_face",
    "cad_edge",
    "cad_vertex",
    "cad_solid",
)


@dataclass(frozen=True)
class SelectionSet:
    """Immutable set of 0-based entity IDs for a single CAD or mesh domain kind.

    Parameters
    ----------
    kind:
        Entity domain kind (``"node"``, ``"element"``, ``"surface_tri"``,
        ``"cad_face"``, ``"cad_edge"``, ``"cad_vertex"``, ``"cad_solid"``).
    ids:
        Iterable of integer entity IDs. Stored deduplicated and sorted.
    name:
        Optional human-readable set name.
    universe:
        Optional full set of valid entity IDs for complement (``~``) operations.
    metadata:
        Optional auxiliary metadata dictionary.
    """

    kind: SetKind
    ids: tuple[int, ...] = ()
    name: str = ""
    universe: frozenset[int] | None = None
    metadata: dict[str, Any] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if self.kind not in VALID_SET_KINDS:
            raise ValueError(
                f"Invalid SelectionSet kind {self.kind!r}; expected one of {VALID_SET_KINDS}"
            )
        normalized = tuple(sorted({int(i) for i in self.ids}))
        object.__setattr__(self, "ids", normalized)
        if self.universe is not None and not isinstance(self.universe, frozenset):
            object.__setattr__(self, "universe", frozenset(int(i) for i in self.universe))

    @classmethod
    def from_ids(
        cls,
        kind: SetKind,
        ids: Iterable[int],
        *,
        name: str = "",
        universe: Iterable[int] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> SelectionSet:
        """Construct a :class:`SelectionSet` from any iterable of integer IDs.

        Parameters
        ----------
        kind:
            Entity domain kind.
        ids:
            Iterable of integer IDs.
        name:
            Optional set name.
        universe:
            Optional universe of valid IDs for complement operations.
        metadata:
            Optional metadata dictionary.
        """
        uni = frozenset(int(i) for i in universe) if universe is not None else None
        return cls(
            kind=kind,
            ids=tuple(int(i) for i in ids),
            name=name,
            universe=uni,
            metadata=dict(metadata) if metadata else {},
        )

    @property
    def id_set(self) -> frozenset[int]:
        """Return entity IDs as a ``frozenset[int]``."""
        return frozenset(self.ids)

    @property
    def is_empty(self) -> bool:
        """Return ``True`` if the set contains no entity IDs."""
        return len(self.ids) == 0

    def __len__(self) -> int:
        return len(self.ids)

    def __iter__(self):
        return iter(self.ids)

    def __contains__(self, item: object) -> bool:
        return isinstance(item, (int, np.integer)) and int(item) in self.id_set

    def __bool__(self) -> bool:
        return len(self.ids) > 0

    def to_list(self) -> list[int]:
        """Return sorted list of entity IDs."""
        return list(self.ids)

    def to_array(self) -> np.ndarray:
        """Return sorted ``int64`` NumPy array of entity IDs."""
        return np.asarray(self.ids, dtype=np.int64)

    def with_name(self, name: str) -> SelectionSet:
        """Return a copy of this set with an updated name."""
        return SelectionSet(
            kind=self.kind,
            ids=self.ids,
            name=name,
            universe=self.universe,
            metadata=dict(self.metadata),
        )

    def with_universe(self, universe: Iterable[int]) -> SelectionSet:
        """Return a copy of this set with an attached universe."""
        return SelectionSet(
            kind=self.kind,
            ids=self.ids,
            name=self.name,
            universe=frozenset(int(i) for i in universe),
            metadata=dict(self.metadata),
        )

    def _require_compatible(self, other: SelectionSet) -> None:
        if not isinstance(other, SelectionSet):
            raise TypeError(f"Expected SelectionSet, got {type(other).__name__}")
        if self.kind != other.kind:
            raise ValueError(
                f"Cannot combine SelectionSets of different kinds: {self.kind!r} vs {other.kind!r}"
            )

    def _merged_universe(self, other: SelectionSet) -> frozenset[int] | None:
        if self.universe is not None and other.universe is not None:
            return self.universe | other.universe
        if self.universe is not None:
            return self.universe
        return other.universe

    def union(self, other: SelectionSet) -> SelectionSet:
        """Return the set union ``self | other``.

        Parameters
        ----------
        other:
            Another :class:`SelectionSet` of the same ``kind``.
        """
        self._require_compatible(other)
        return SelectionSet(
            kind=self.kind,
            ids=tuple(sorted(self.id_set | other.id_set)),
            universe=self._merged_universe(other),
        )

    def intersection(self, other: SelectionSet) -> SelectionSet:
        """Return the set intersection ``self & other``.

        Parameters
        ----------
        other:
            Another :class:`SelectionSet` of the same ``kind``.
        """
        self._require_compatible(other)
        return SelectionSet(
            kind=self.kind,
            ids=tuple(sorted(self.id_set & other.id_set)),
            universe=self._merged_universe(other),
        )

    def difference(self, other: SelectionSet) -> SelectionSet:
        """Return the set difference ``self - other``.

        Parameters
        ----------
        other:
            Another :class:`SelectionSet` of the same ``kind``.
        """
        self._require_compatible(other)
        return SelectionSet(
            kind=self.kind,
            ids=tuple(sorted(self.id_set - other.id_set)),
            universe=self._merged_universe(other),
        )

    def symmetric_difference(self, other: SelectionSet) -> SelectionSet:
        """Return the symmetric difference ``self ^ other``.

        Parameters
        ----------
        other:
            Another :class:`SelectionSet` of the same ``kind``.
        """
        self._require_compatible(other)
        return SelectionSet(
            kind=self.kind,
            ids=tuple(sorted(self.id_set ^ other.id_set)),
            universe=self._merged_universe(other),
        )

    def complement(self, universe: Iterable[int] | None = None) -> SelectionSet:
        """Return the complement ``~self`` within ``universe`` (or ``self.universe``).

        Parameters
        ----------
        universe:
            Optional iterable of all valid entity IDs. Required if ``self.universe``
            is ``None``.
        """
        uni = frozenset(int(i) for i in universe) if universe is not None else self.universe
        if uni is None:
            raise ValueError(
                "complement (~) requires a universe of all valid IDs; "
                "pass universe=... or construct the SelectionSet with universe"
            )
        return SelectionSet(
            kind=self.kind,
            ids=tuple(sorted(uni - self.id_set)),
            universe=uni,
        )

    def __or__(self, other: SelectionSet) -> SelectionSet:
        return self.union(other)

    def __and__(self, other: SelectionSet) -> SelectionSet:
        return self.intersection(other)

    def __sub__(self, other: SelectionSet) -> SelectionSet:
        return self.difference(other)

    def __xor__(self, other: SelectionSet) -> SelectionSet:
        return self.symmetric_difference(other)

    def __invert__(self) -> SelectionSet:
        return self.complement()

    def to_dict(self) -> dict[str, Any]:
        """Serialize this set to a JSON-compatible dictionary."""
        return {
            "kind": self.kind,
            "name": self.name,
            "ids": list(self.ids),
            "universe": sorted(self.universe) if self.universe is not None else None,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SelectionSet:
        """Deserialize a :class:`SelectionSet` from a dictionary."""
        uni = data.get("universe")
        return cls(
            kind=data["kind"],
            ids=tuple(int(i) for i in data.get("ids", [])),
            name=str(data.get("name", "")),
            universe=frozenset(int(i) for i in uni) if uni is not None else None,
            metadata=dict(data.get("metadata", {})),
        )


class SetRegistry:
    """Registry managing named :class:`SelectionSet` instances and boolean derivations."""

    def __init__(self) -> None:
        self._sets: dict[str, SelectionSet] = {}

    def define_set(
        self,
        name: str,
        selection_set: SelectionSet,
        *,
        overwrite: bool = True,
    ) -> SelectionSet:
        """Register a :class:`SelectionSet` under ``name``.

        Parameters
        ----------
        name:
            Non-empty registry key for the set.
        selection_set:
            The :class:`SelectionSet` (or object with ``to_set()``) to register.
        overwrite:
            If ``False`` and ``name`` already exists, raise :class:`KeyError`.
        """
        if not name or not str(name).strip():
            raise ValueError("Set name must be a non-empty string")
        if not overwrite and name in self._sets:
            raise KeyError(f"SelectionSet {name!r} is already defined")
        if not isinstance(selection_set, SelectionSet):
            if hasattr(selection_set, "to_set"):
                selection_set = selection_set.to_set(name=name)
            else:
                raise TypeError(
                    f"Expected SelectionSet or query with to_set(), got {type(selection_set).__name__}"
                )
        stored = selection_set.with_name(name)
        self._sets[name] = stored
        return stored

    def get(self, name: str) -> SelectionSet:
        """Retrieve a registered :class:`SelectionSet` by name.

        Parameters
        ----------
        name:
            Registered set name.
        """
        if name not in self._sets:
            raise KeyError(f"Unknown SelectionSet {name!r}")
        return self._sets[name]

    def remove(self, name: str) -> SelectionSet:
        """Remove and return a registered :class:`SelectionSet` by name."""
        if name not in self._sets:
            raise KeyError(f"Unknown SelectionSet {name!r}")
        return self._sets.pop(name)

    def names(self) -> list[str]:
        """Return list of all registered set names."""
        return list(self._sets.keys())

    def by_kind(self, kind: SetKind) -> dict[str, SelectionSet]:
        """Return all registered sets matching ``kind``."""
        return {k: v for k, v in self._sets.items() if v.kind == kind}

    def derive_union(self, name: str, *set_names: str) -> SelectionSet:
        """Create and register a union of two or more registered sets."""
        if not set_names:
            raise ValueError("derive_union requires at least one set name")
        acc = self.get(set_names[0])
        for sname in set_names[1:]:
            acc = acc | self.get(sname)
        return self.define_set(name, acc)

    def derive_intersection(self, name: str, *set_names: str) -> SelectionSet:
        """Create and register an intersection of two or more registered sets."""
        if not set_names:
            raise ValueError("derive_intersection requires at least one set name")
        acc = self.get(set_names[0])
        for sname in set_names[1:]:
            acc = acc & self.get(sname)
        return self.define_set(name, acc)

    def derive_difference(self, name: str, base_name: str, *subtract_names: str) -> SelectionSet:
        """Create and register ``base_name - (s1 | s2 | ...)``."""
        acc = self.get(base_name)
        for sname in subtract_names:
            acc = acc - self.get(sname)
        return self.define_set(name, acc)

    def derive_symmetric_difference(self, name: str, a_name: str, b_name: str) -> SelectionSet:
        """Create and register ``a_name ^ b_name``."""
        return self.define_set(name, self.get(a_name) ^ self.get(b_name))

    def derive_complement(
        self,
        name: str,
        source_name: str,
        universe: Iterable[int] | None = None,
    ) -> SelectionSet:
        """Create and register the complement of ``source_name``."""
        return self.define_set(name, self.get(source_name).complement(universe=universe))

    def __getitem__(self, name: str) -> SelectionSet:
        return self.get(name)

    def __contains__(self, name: object) -> bool:
        return name in self._sets

    def __len__(self) -> int:
        return len(self._sets)

    def __iter__(self):
        return iter(self._sets)

    def items(self):
        """Return ``(name, SelectionSet)`` pairs."""
        return self._sets.items()

    def to_dict(self) -> dict[str, Any]:
        """Serialize all registered sets to a JSON-compatible dictionary."""
        return {k: v.to_dict() for k, v in self._sets.items()}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SetRegistry:
        """Deserialize a :class:`SetRegistry` from a dictionary."""
        reg = cls()
        for k, v in data.items():
            reg.define_set(k, SelectionSet.from_dict(v))
        return reg
