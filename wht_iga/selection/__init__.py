"""Unified CAD + Mesh selection, geometric measurement, and contact check API (v2.1 §8.5)."""

from __future__ import annotations

from .contact_check import ContactGapReport, inspect_contact
from .measure import BBoxResult, MeasureResult, measure_cad, measure_mesh
from .selector import CadAdapter, EntityQuery, MeshAdapter
from .sets import VALID_SET_KINDS, SelectionSet, SetKind, SetRegistry

__all__ = [
    "BBoxResult",
    "CadAdapter",
    "ContactGapReport",
    "EntityQuery",
    "MeasureResult",
    "MeshAdapter",
    "SelectionSet",
    "SetKind",
    "SetRegistry",
    "VALID_SET_KINDS",
    "inspect_contact",
    "measure_cad",
    "measure_mesh",
]
