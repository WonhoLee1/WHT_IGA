"""Mesh readers, writers, body-fit NURBS converter, and internal model (v2.1 §5.10). Kratos-free.

Scope is mesh + sets ONLY (K25): loads, materials, steps and other
analysis-deck keywords are reported as issues and never interpreted.
``unit_system`` is mandatory on every read (K26). Shell/beam/spring
elements are read, issued and excluded (K27).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from wht_iga.meshio_ext.bodyfit import (
    BodyFitConversionResult,
    BodyFitSplinePatch,
    BodyFitValidationError,
    PatchInterface,
    ValidationIssue,
    build_kratos_nurbs_volumes,
    mesh_to_bodyfit_patches,
)
from wht_iga.meshio_ext.model import ElemBlock, MeshModel
from wht_iga.meshio_ext.read_inp import read_inp
from wht_iga.meshio_ext.read_k import read_k
from wht_iga.meshio_ext.read_msh import read_msh
from wht_iga.meshio_ext.read_nastran import read_nastran
from wht_iga.meshio_ext.read_rad import read_rad
from wht_iga.meshio_ext.write_inp import write_inp
from wht_iga.meshio_ext.write_k import write_k
from wht_iga.meshio_ext.write_msh import write_msh
from wht_iga.meshio_ext.write_nastran import write_nastran
from wht_iga.meshio_ext.write_rad import write_rad

_READERS = {
    ".inp": read_inp,
    "inp": read_inp,
    ".bdf": read_nastran,
    ".fem": read_nastran,
    ".nas": read_nastran,
    "bdf": read_nastran,
    "fem": read_nastran,
    "nastran": read_nastran,
    "optistruct": read_nastran,
    ".rad": read_rad,
    "rad": read_rad,
    ".k": read_k,
    ".key": read_k,
    ".dyn": read_k,
    "k": read_k,
    "key": read_k,
    ".msh": read_msh,
    "msh": read_msh,
}

_WRITERS = {
    ".inp": write_inp,
    "inp": write_inp,
    ".bdf": write_nastran,
    ".fem": write_nastran,
    ".nas": write_nastran,
    "bdf": write_nastran,
    "fem": write_nastran,
    "nastran": write_nastran,
    "optistruct": write_nastran,
    ".rad": write_rad,
    "rad": write_rad,
    ".k": write_k,
    ".key": write_k,
    ".dyn": write_k,
    "k": write_k,
    "key": write_k,
    ".msh": write_msh,
    "msh": write_msh,
}


def read_mesh(
    path: str | Path,
    unit_system: str,
    file_format: str | None = None,
) -> tuple[MeshModel, list[dict]]:
    """Read a solid mesh from any supported format dispatched by extension or ``file_format``.

    Supported formats: ``.inp`` (Abaqus), ``.bdf``/``.fem``/``.nas`` (Nastran/OptiStruct),
    ``.rad`` (Radioss), ``.k``/``.key``/``.dyn`` (LS-DYNA), ``.msh`` (Gmsh).
    """
    p = Path(path)
    key = (file_format.lower().strip() if file_format else p.suffix.lower())
    reader = _READERS.get(key)
    if reader is None:
        raise ValueError(
            f"Unsupported mesh format {key!r} for {str(p)!r}; "
            "expected .inp, .bdf, .fem, .rad, .k, or .msh"
        )
    return reader(str(p), unit_system)


def write_mesh(
    model: MeshModel,
    path: str | Path,
    file_format: str | None = None,
    extra_sets: Any = None,
    **kwargs: Any,
) -> Path:
    """Write a :class:`MeshModel` to any supported format dispatched by extension or ``file_format``."""
    p = Path(path)
    key = (file_format.lower().strip() if file_format else p.suffix.lower())
    writer = _WRITERS.get(key)
    if writer is None:
        raise ValueError(
            f"Unsupported mesh format {key!r} for {str(p)!r}; "
            "expected .inp, .bdf, .fem, .rad, .k, or .msh"
        )
    if writer is write_nastran and "format_dialect" not in kwargs and file_format:
        if key in ("optistruct", "fem", ".fem"):
            kwargs["format_dialect"] = "optistruct"
        elif key in ("nastran", "bdf", ".bdf"):
            kwargs["format_dialect"] = "nastran"
    return writer(model, p, extra_sets=extra_sets, **kwargs)


__all__ = [
    "BodyFitConversionResult",
    "BodyFitSplinePatch",
    "BodyFitValidationError",
    "ElemBlock",
    "MeshModel",
    "PatchInterface",
    "ValidationIssue",
    "build_kratos_nurbs_volumes",
    "mesh_to_bodyfit_patches",
    "read_inp",
    "read_k",
    "read_mesh",
    "read_msh",
    "read_nastran",
    "read_rad",
    "write_inp",
    "write_k",
    "write_mesh",
    "write_msh",
    "write_nastran",
    "write_rad",
]
