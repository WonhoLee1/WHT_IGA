"""Unified OpenCASCADE STEP / IGES / BREP loader and B-Rep healing pipeline (T-P3.5).

Loads ``.step``/``.stp``, ``.iges``/``.igs``, and ``.brep`` CAD files into a
:class:`CADGeometry` container with topological entity lists (solids, faces,
edges, vertices) and a structured :class:`HealingReport`. When ``heal=True``,
free faces/shells (such as trimmed NURBS surfaces loaded from IGES mode 0) are
sewn into watertight shells via ``BRepBuilderAPI_Sewing``, upgraded to solids
via ``BRepBuilderAPI_MakeSolid``, and oriented with ``ShapeFix_Solid``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from wht_iga.geometry.cad import (
    bbox_diag,
    brep_area,
    brep_volume,
    shape_bbox,
)

STEP_EXTENSIONS: frozenset[str] = frozenset({".step", ".stp"})
IGES_EXTENSIONS: frozenset[str] = frozenset({".iges", ".igs"})
BREP_EXTENSIONS: frozenset[str] = frozenset({".brep"})
SAT_EXTENSIONS: frozenset[str] = frozenset({".sat"})
SUPPORTED_CAD_EXTENSIONS: frozenset[str] = (
    STEP_EXTENSIONS | IGES_EXTENSIONS | BREP_EXTENSIONS | SAT_EXTENSIONS
)


@dataclass
class HealingReport:
    """Summary of B-Rep healing and solid recovery performed by :func:`load_cad_file`.

    Parameters
    ----------
    healed:
        Whether healing (``ShapeFix_Shape``, sewing, solid upgrade) was executed.
    sewing_tol:
        Sewing tolerance passed to ``BRepBuilderAPI_Sewing`` (0.0 if ``heal=False``).
    solids_before, solids_after:
        Unique solid counts before and after healing.
    shells_before, shells_after:
        Unique shell counts before and after healing.
    faces_before, faces_after:
        Unique face counts before and after healing.
    edges_before, edges_after:
        Unique non-degenerated edge counts before and after healing.
    free_edges_after:
        Free (unshared) non-degenerated edges remaining after sewing.
    solids_constructed:
        Number of solids newly promoted from free shells via ``BRepBuilderAPI_MakeSolid``.
    volume_before, volume_after:
        Exact B-Rep volumes computed via ``BRepGProp.VolumeProperties_s``.
    area_before, area_after:
        Exact B-Rep surface areas computed via ``BRepGProp.SurfaceProperties_s``.
    watertight:
        ``True`` when at least one solid exists and zero free edges remain.
    issues:
        Non-fatal diagnostic records produced during healing/validation.
    """

    healed: bool = False
    sewing_tol: float = 0.0
    solids_before: int = 0
    solids_after: int = 0
    shells_before: int = 0
    shells_after: int = 0
    faces_before: int = 0
    faces_after: int = 0
    edges_before: int = 0
    edges_after: int = 0
    free_edges_after: int = 0
    solids_constructed: int = 0
    volume_before: float = 0.0
    volume_after: float = 0.0
    area_before: float = 0.0
    area_after: float = 0.0
    watertight: bool = False
    issues: list[dict[str, str]] = field(default_factory=list)


@dataclass
class CADGeometry:
    """Unified CAD geometry wrapper compatible with :mod:`wht_iga.geometry.cad` and ``CadAdapter``.

    Parameters
    ----------
    shape:
        OpenCASCADE ``TopoDS_Shape`` instance (typically a ``TopoDS_Solid`` or
        ``TopoDS_Compound`` after healing).
    solids:
        Unique ``TopoDS_Solid`` subshapes extracted from ``shape``.
    faces:
        Unique ``TopoDS_Face`` subshapes extracted from ``shape``.
    edges:
        Unique ``TopoDS_Edge`` subshapes extracted from ``shape``.
    vertices:
        Unique ``TopoDS_Vertex`` subshapes extracted from ``shape``.
    healing_report:
        Detailed :class:`HealingReport` recording pre/post healing topology and volume.
    file_format:
        Normalized format tag: ``"step"``, ``"iges"``, or ``"brep"``.
    source_path:
        String representation of the loaded source CAD file path.
    unit_system:
        Optional unit label (e.g. ``"mm"``, ``"m"``) associated with the geometry.
    """

    shape: Any
    solids: list[Any]
    faces: list[Any]
    edges: list[Any]
    vertices: list[Any]
    healing_report: HealingReport
    file_format: str
    source_path: str
    unit_system: str | None = None

    @property
    def volume(self) -> float:
        """Exact B-Rep volume via ``BRepGProp``."""
        return brep_volume(self.shape)

    @property
    def area(self) -> float:
        """Exact B-Rep surface area via ``BRepGProp``."""
        return brep_area(self.shape)

    @property
    def bbox(self) -> tuple[np.ndarray, np.ndarray]:
        """Axis-aligned bounding box ``(lo, hi)`` of ``self.shape``."""
        return shape_bbox(self.shape)

    @property
    def has_solid(self) -> bool:
        """Whether at least one ``TopAbs_SOLID`` exists in the loaded shape."""
        return len(self.solids) > 0


def _collect_unique(shape: Any, shape_type: Any, avoid_type: Any = None) -> list[Any]:
    """Collect topologically unique subshapes using ``IsSame`` deduplication."""
    from OCP.TopExp import TopExp_Explorer

    exp = (
        TopExp_Explorer(shape, shape_type)
        if avoid_type is None
        else TopExp_Explorer(shape, shape_type, avoid_type)
    )
    items: list[Any] = []
    while exp.More():
        cur = exp.Current()
        if not any(prev.IsSame(cur) for prev in items):
            items.append(cur)
        exp.Next()
    return items


def _cast_topods(items: list[Any], kind: str) -> list[Any]:
    """Downcast collected ``TopoDS_Shape`` items to typed OCC wrappers."""
    from OCP.TopoDS import TopoDS

    caster = {
        "solid": TopoDS.Solid,
        "shell": TopoDS.Shell,
        "face": TopoDS.Face,
        "edge": TopoDS.Edge,
        "vertex": TopoDS.Vertex,
    }[kind]
    return [caster(item) for item in items]


def _collect_non_degenerated_edges(shape: Any) -> list[Any]:
    """Collect unique non-degenerated edges in ``shape``."""
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_EDGE

    raw_edges = _cast_topods(_collect_unique(shape, TopAbs_EDGE), "edge")
    return [e for e in raw_edges if not BRep_Tool.Degenerated_s(e)]


def _count_free_edges(shape: Any) -> int:
    """Count non-degenerated edges shared by fewer than 2 unique faces."""
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    faces = _cast_topods(_collect_unique(shape, TopAbs_FACE), "face")
    unique_edges: list[Any] = []
    edge_counts: list[int] = []

    for face in faces:
        seen_in_face: list[Any] = []
        exp = TopExp_Explorer(face, TopAbs_EDGE)
        while exp.More():
            edge = TopoDS.Edge(exp.Current())
            if (
                not BRep_Tool.Degenerated_s(edge)
                and not BRep_Tool.IsClosed_s(edge, face)
                and not any(prev.IsSame(edge) for prev in seen_in_face)
            ):
                seen_in_face.append(edge)
                matched = False
                for idx, ref in enumerate(unique_edges):
                    if ref.IsSame(edge):
                        edge_counts[idx] += 1
                        matched = True
                        break
                if not matched:
                    unique_edges.append(edge)
                    edge_counts.append(1)
            exp.Next()

    return sum(1 for count in edge_counts if count < 2)


def _read_raw_cad_shape(path: Path) -> tuple[Any, str]:
    """Read raw ``TopoDS_Shape`` from ``.step``/``.stp``, ``.iges``/``.igs``, or ``.brep``."""
    ext = path.suffix.lower()
    if ext in STEP_EXTENSIONS:
        from OCP.STEPControl import STEPControl_Reader

        reader = STEPControl_Reader()
        if int(reader.ReadFile(str(path))) != 1:
            raise IOError(f"STEP read failed: {path}")
        reader.TransferRoots()
        shape = reader.OneShape()
        if shape.IsNull():
            raise IOError(f"STEP produced empty shape: {path}")
        return shape, "step"

    if ext in IGES_EXTENSIONS:
        from OCP.IGESControl import IGESControl_Reader

        reader = IGESControl_Reader()
        if int(reader.ReadFile(str(path))) != 1:
            raise IOError(f"IGES read failed: {path}")
        reader.TransferRoots()
        shape = reader.OneShape()
        if shape.IsNull():
            raise IOError(f"IGES produced empty shape: {path}")
        return shape, "iges"

    if ext in BREP_EXTENSIONS:
        from OCP.BRep import BRep_Builder
        from OCP.BRepTools import BRepTools
        from OCP.TopoDS import TopoDS_Shape

        shape = TopoDS_Shape()
        if not BRepTools.Read_s(shape, str(path), BRep_Builder()) or shape.IsNull():
            raise IOError(f"BREP read failed: {path}")
        return shape, "brep"

    if ext in SAT_EXTENSIONS:
        from wht_iga.geometry.sat_io import read_sat_shape

        return read_sat_shape(path), "sat"

    raise ValueError(
        f"Unsupported CAD extension {ext!r}; expected one of {sorted(SUPPORTED_CAD_EXTENSIONS)}"
    )


def _heal_and_make_solids(
    raw_shape: Any,
    sewing_tol: float,
) -> tuple[Any, int, int, list[dict[str, str]]]:
    """Run ``ShapeFix_Shape``, ``BRepBuilderAPI_Sewing``, ``MakeSolid``, and ``ShapeFix_Solid``."""
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Sewing
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Shell, ShapeFix_Solid
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID
    from OCP.TopoDS import TopoDS, TopoDS_Compound

    issues: list[dict[str, str]] = []

    fixer = ShapeFix_Shape(raw_shape)
    fixer.SetPrecision(sewing_tol)
    fixer.SetMinTolerance(min(1e-9, sewing_tol))
    fixer.SetMaxTolerance(max(sewing_tol * 10.0, 1e-4))
    fixer.Perform()
    fixed_raw = fixer.Shape()

    existing_solids = _collect_unique(fixed_raw, TopAbs_SOLID)
    free_shells_raw = _collect_unique(fixed_raw, TopAbs_SHELL, TopAbs_SOLID)
    free_faces_raw = _collect_unique(fixed_raw, TopAbs_FACE, TopAbs_SHELL)

    # Fast path when the input already consists purely of watertight solids.
    if existing_solids and not free_shells_raw and not free_faces_raw and _count_free_edges(fixed_raw) == 0:
        healed_solids: list[Any] = []
        for solid_raw in existing_solids:
            sfs = ShapeFix_Solid(TopoDS.Solid(solid_raw))
            sfs.SetPrecision(sewing_tol)
            sfs.Perform()
            healed_solid = TopoDS.Solid(sfs.Solid())
            if brep_volume(healed_solid) < 0.0:
                healed_solid.Reverse()
            healed_solids.append(healed_solid)
        if len(healed_solids) == 1:
            return healed_solids[0], 0, 0, issues
        comp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(comp)
        for s in healed_solids:
            builder.Add(comp, s)
        return comp, 0, 0, issues

    sewing = BRepBuilderAPI_Sewing(sewing_tol)
    sewing.Add(fixed_raw)
    sewing.Perform()
    sewed = sewing.SewedShape()
    free_edges = int(sewing.NbFreeEdges())

    result_parts: list[Any] = []
    solids_constructed = 0

    for solid_raw in _collect_unique(sewed, TopAbs_SOLID):
        solid = TopoDS.Solid(solid_raw)
        sfs = ShapeFix_Solid(solid)
        sfs.SetPrecision(sewing_tol)
        sfs.Perform()
        healed_solid = TopoDS.Solid(sfs.Solid())
        if brep_volume(healed_solid) < 0.0:
            healed_solid.Reverse()
        result_parts.append(healed_solid)

    for shell_raw in _collect_unique(sewed, TopAbs_SHELL, TopAbs_SOLID):
        shell = TopoDS.Shell(shell_raw)
        sf_shell = ShapeFix_Shell(shell)
        sf_shell.SetPrecision(sewing_tol)
        sf_shell.Perform()
        clean_shell = TopoDS.Shell(sf_shell.Shell())

        shell_free = _count_free_edges(clean_shell)
        if shell_free == 0:
            maker = BRepBuilderAPI_MakeSolid(clean_shell)
            if maker.IsDone():
                new_solid = maker.Solid()
                sfs = ShapeFix_Solid(new_solid)
                sfs.SetPrecision(sewing_tol)
                sfs.Perform()
                healed_solid = TopoDS.Solid(sfs.Solid())
                if brep_volume(healed_solid) < 0.0:
                    healed_solid.Reverse()
                if BRepCheck_Analyzer(healed_solid).IsValid() or brep_volume(healed_solid) > 0.0:
                    result_parts.append(healed_solid)
                    solids_constructed += 1
                    continue
        issues.append(
            {
                "severity": "warning",
                "code": "OPEN_SHELL_REMAINS",
                "message": f"Shell with {shell_free} free edges could not be upgraded to a solid.",
                "entity": "CADGeometry",
            }
        )
        result_parts.append(clean_shell)

    # Preserve any unsewn loose faces (not belonging to any shell or solid).
    for loose_face in _collect_unique(sewed, TopAbs_FACE, TopAbs_SHELL):
        issues.append(
            {
                "severity": "warning",
                "code": "UNSEWN_FACE_REMAINS",
                "message": "Loose face remained outside shells after sewing.",
                "entity": "CADGeometry",
            }
        )
        result_parts.append(TopoDS.Face(loose_face))

    if not result_parts:
        return sewed, free_edges, solids_constructed, issues

    if len(result_parts) == 1:
        final_shape = result_parts[0]
    else:
        comp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(comp)
        for part in result_parts:
            builder.Add(comp, part)
        final_shape = comp

    return final_shape, free_edges, solids_constructed, issues


def load_cad_file(
    path: str | Path,
    heal: bool = True,
    sewing_tol: float | None = 1e-6,
    unit_system: str | None = None,
) -> CADGeometry:
    """Load a STEP, IGES, or BREP CAD file and optionally heal shells into solids.

    Parameters
    ----------
    path:
        Input CAD file path (``.step``, ``.stp``, ``.iges``, ``.igs``, or ``.brep``).
    heal:
        When ``True``, run ``ShapeFix_Shape``, ``BRepBuilderAPI_Sewing``,
        ``BRepBuilderAPI_MakeSolid`` on closed shells, and ``ShapeFix_Solid``.
    sewing_tol:
        Sewing tolerance for ``BRepBuilderAPI_Sewing``. Pass ``None`` to use
        ``1e-6 * bbox_diag(raw_shape)``.
    unit_system:
        Optional unit system tag stored on the returned :class:`CADGeometry`.

    Returns
    -------
    CADGeometry
        Healed (or raw if ``heal=False``) geometry container with extracted
        topology lists and :class:`HealingReport`.
    """
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX

    src = Path(path)
    if not src.is_file():
        raise FileNotFoundError(f"CAD file not found: {src}")

    raw_shape, fmt = _read_raw_cad_shape(src)

    solids_before = len(_collect_unique(raw_shape, TopAbs_SOLID))
    shells_before = len(_collect_unique(raw_shape, TopAbs_SHELL))
    faces_before = len(_collect_unique(raw_shape, TopAbs_FACE))
    edges_before = len(_collect_non_degenerated_edges(raw_shape))
    vol_before = brep_volume(raw_shape)
    area_before = brep_area(raw_shape)

    issues: list[dict[str, str]] = []
    if heal:
        eff_tol = 1e-6 * bbox_diag(raw_shape) if sewing_tol is None else float(sewing_tol)
        if eff_tol <= 0.0:
            eff_tol = 1e-6
        shape, sewing_free_edges, solids_constructed, heal_issues = _heal_and_make_solids(
            raw_shape, eff_tol
        )
        issues.extend(heal_issues)
        free_edges_after = max(sewing_free_edges, _count_free_edges(shape))
    else:
        eff_tol = 0.0
        shape = raw_shape
        solids_constructed = 0
        free_edges_after = _count_free_edges(shape)

    solids = _cast_topods(_collect_unique(shape, TopAbs_SOLID), "solid")
    shells = _cast_topods(_collect_unique(shape, TopAbs_SHELL), "shell")
    faces = _cast_topods(_collect_unique(shape, TopAbs_FACE), "face")
    edges = _collect_non_degenerated_edges(shape)
    vertices = _cast_topods(_collect_unique(shape, TopAbs_VERTEX), "vertex")

    vol_after = brep_volume(shape)
    area_after = brep_area(shape)
    watertight = len(solids) > 0 and free_edges_after == 0 and vol_after > 0.0

    report = HealingReport(
        healed=heal,
        sewing_tol=eff_tol,
        solids_before=solids_before,
        solids_after=len(solids),
        shells_before=shells_before,
        shells_after=len(shells),
        faces_before=faces_before,
        faces_after=len(faces),
        edges_before=edges_before,
        edges_after=len(edges),
        free_edges_after=free_edges_after,
        solids_constructed=solids_constructed,
        volume_before=vol_before,
        volume_after=vol_after,
        area_before=area_before,
        area_after=area_after,
        watertight=watertight,
        issues=issues,
    )

    return CADGeometry(
        shape=shape,
        solids=solids,
        faces=faces,
        edges=edges,
        vertices=vertices,
        healing_report=report,
        file_format=fmt,
        source_path=str(src),
        unit_system=unit_system,
    )


__all__ = [
    "BREP_EXTENSIONS",
    "CADGeometry",
    "HealingReport",
    "IGES_EXTENSIONS",
    "SAT_EXTENSIONS",
    "STEP_EXTENSIONS",
    "SUPPORTED_CAD_EXTENSIONS",
    "load_cad_file",
]
