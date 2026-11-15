"""Multi-part, dependent-instance, multi-cell, and THB assembly builder for Kratos.

Implements Track P3-C (`T-P3.1` & `T-P3.2`) over a single root ``Structure``
``ModelPart`` (S15 compliance), managing:

* Monotonic non-overlapping Node, Element, Condition, Properties, and Constraint
  ID ranges via :class:`IdOffsetTable` serialized to ``<job>.resolved.json``.
* Dependent vs. independent part instances with rigid coordinate
  transformations ``x' = R x + t`` (and consistent vector/tensor transforms).
* Multi-cell shared B-spline fields (`T-P3.1` Path B) where adjacent CAD cells
  or material domains share common Greville control points (``C^{p-1}``
  continuity) while each cell's quadrature points receive distinct ``Properties``.
* Multi-level THB assembly (`T-P3.1` Path A) integrating :class:`THBHierarchy`
  with ``C_THB`` and optional trimming extension ``E`` (``T = C_THB @ E``),
  assembling Kratos ``SmallDisplacementElement3D4N`` elements on fine
  sub-spans of leaf cells and supporting both algebraic reduction
  (``K_red = T^T K_hier T``) and Kratos ``LinearMasterSlaveConstraint`` solve.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import scipy.sparse as sp

from wht_iga.grid.hierarchy import LeafCell, THBHierarchy
from wht_iga.grid.knots import breakpoints, uniform_clamped
from wht_iga.meshio_ext.model import MeshModel
from wht_iga.operators.compose import (
    apply_TKT,
    apply_Tt_vec,
    compose_thb,
    expand_dof_operator,
    solve_reduced,
)
from wht_iga.reference.box import BoxSpline


_ProvenanceKind = Literal[
    "control_point",
    "thb_master",
    "fe_node",
    "surface_layer",
    "virtual",
    "reference_point",
]


def validate_rotation_matrix(
    rotation_matrix: np.ndarray | None,
    tol: float = 1e-10,
) -> np.ndarray:
    """Validate and return a 3x3 proper orthogonal rotation matrix.

    Parameters
    ----------
    rotation_matrix:
        Optional 3x3 rotation matrix. If ``None``, returns ``np.eye(3)``.
    tol:
        Tolerance for orthogonality (``||R^T R - I||``) and determinant check.

    Returns
    -------
    np.ndarray
        Validated ``(3, 3)`` ``float64`` rotation matrix with ``det(R) == +1``.
    """
    if rotation_matrix is None:
        return np.eye(3, dtype=np.float64)
    R = np.asarray(rotation_matrix, dtype=np.float64)
    if R.shape != (3, 3):
        raise ValueError(f"rotation_matrix must have shape (3, 3), got {R.shape}")
    ortho_err = float(np.max(np.abs(R.T @ R - np.eye(3))))
    if ortho_err > tol:
        raise ValueError(f"rotation_matrix is not orthogonal (max err={ortho_err:.3e})")
    det = float(np.linalg.det(R))
    if abs(det - 1.0) > tol:
        raise ValueError(f"rotation_matrix must be proper (det=+1), got det={det:.6f}")
    return R.copy()


def transform_stress_voigt(
    stress_voigt: np.ndarray,
    rotation_matrix: np.ndarray,
) -> np.ndarray:
    """Rotate 3D Voigt stress vector ``[xx, yy, zz, xy, yz, xz]`` via ``R @ S @ R.T``.

    Parameters
    ----------
    stress_voigt:
        Array of shape ``(..., 6)`` in Kratos Voigt order ``[xx, yy, zz, xy, yz, xz]``.
    rotation_matrix:
        Proper ``(3, 3)`` rotation matrix mapping local coordinates to global.

    Returns
    -------
    np.ndarray
        Rotated Voigt stress array of the same shape ``(..., 6)``.
    """
    s = np.asarray(stress_voigt, dtype=np.float64)
    R = validate_rotation_matrix(rotation_matrix)
    single = s.ndim == 1
    if single:
        s = s[None, :]
    if s.shape[-1] != 6:
        raise ValueError(f"Expected Voigt stress last dimension 6, got {s.shape}")
    out = np.empty_like(s)
    for idx in range(s.shape[0]):
        xx, yy, zz, xy, yz, xz = s[idx]
        mat = np.array(
            [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
            dtype=np.float64,
        )
        rot = R @ mat @ R.T
        out[idx] = [
            rot[0, 0],
            rot[1, 1],
            rot[2, 2],
            rot[0, 1],
            rot[1, 2],
            rot[0, 2],
        ]
    return out[0] if single else out


def compose_thb_operator(
    extension_op: sp.spmatrix | None,
    c_thb: sp.spmatrix | None,
) -> sp.csr_matrix | None:
    """Compose ``T = C_THB @ E`` supporting both square and rectangular ``C_THB``.

    Calls :func:`wht_iga.operators.compose.compose_thb` when ``C_THB`` is absent
    or square, and handles rectangular ``(n_fine, n_thb)`` ``C_THB`` matrices
    returned by :meth:`THBHierarchy.build_c_thb_matrix`.

    Parameters
    ----------
    extension_op:
        Optional extension operator ``E`` (or ``None`` if no extension is used).
    c_thb:
        Optional THB prolongation matrix ``C_THB`` mapping active THB DOFs to
        fine control-point DOFs.

    Returns
    -------
    sp.csr_matrix | None
        Composed scalar basis prolongation operator ``T``, or ``None`` if both
        inputs are ``None``.
    """
    if c_thb is None and extension_op is None:
        return None
    if c_thb is None:
        assert extension_op is not None
        return compose_thb(extension_op, None)
    C = c_thb.tocsr()
    if extension_op is None:
        return C
    E = extension_op.tocsr()
    if C.shape[0] == C.shape[1] and C.shape[1] == E.shape[0]:
        return compose_thb(E, C)
    if C.shape[1] == E.shape[0]:
        return (C @ E).tocsr()
    raise ValueError(
        f"Incompatible shapes for C_thb {C.shape} and extension_op {E.shape}"
    )


@dataclass(frozen=True)
class CellDomainSpec:
    """Material and domain descriptor for one CAD cell inside a shared field.

    Parameters
    ----------
    cell_id:
        Unique positive integer cell/material domain ID used as Kratos ``Properties.Id``.
    name:
        Human-readable name for the cell or material phase.
    young_modulus:
        Young's modulus ``E`` assigned to quadrature points in this cell.
    poisson_ratio:
        Poisson's ratio ``nu`` assigned to quadrature points in this cell.
    thickness:
        Element thickness property passed to Kratos ``Properties``.
    density:
        Mass density assigned to Kratos ``Properties``.
    bounds:
        Optional physical/local axis-aligned box
        ``((xmin, xmax), (ymin, ymax), (zmin, zmax))`` selecting this cell.
    cells:
        Optional explicit list of spline grid cell indices ``(i, j, k)``.
    predicate:
        Optional point predicate ``fn(xyz_local) -> bool`` for arbitrary shapes.
    """

    cell_id: int
    name: str
    young_modulus: float
    poisson_ratio: float
    thickness: float = 1.0
    density: float = 0.0
    bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None = None
    cells: list[tuple[int, int, int]] | None = None
    predicate: Callable[[np.ndarray], bool] | None = None

    def contains_point(self, xyz_local: np.ndarray, tol: float = 1e-12) -> bool:
        """Return whether a local coordinate lies in this cell domain."""
        if self.bounds is not None:
            x, y, z = float(xyz_local[0]), float(xyz_local[1]), float(xyz_local[2])
            (x0, x1), (y0, y1), (z0, z1) = self.bounds
            in_bounds = (
                (x0 - tol <= x <= x1 + tol)
                and (y0 - tol <= y <= y1 + tol)
                and (z0 - tol <= z <= z1 + tol)
            )
            if not in_bounds:
                return False
        if self.predicate is not None:
            return bool(self.predicate(np.asarray(xyz_local, dtype=np.float64)))
        return self.bounds is not None or self.predicate is not None


@dataclass
class PartSpec:
    """Template part specification for IGA, THB, or Native FE geometry.

    Parameters
    ----------
    name:
        Part identifier referenced by :class:`InstanceSpec`.
    kind:
        Part discretization type: ``"iga"``, ``"thb"``, or ``"native_fe"``.
    box:
        Uniform :class:`BoxSpline` template (required when ``kind == "iga"``).
    thb:
        :class:`THBHierarchy` template (required when ``kind == "thb"``).
    mesh:
        :class:`MeshModel` template (required when ``kind == "native_fe"``).
    young_modulus:
        Default Young's modulus when ``cell_domains`` is not specified.
    poisson_ratio:
        Default Poisson's ratio when ``cell_domains`` is not specified.
    thickness:
        Default Kratos ``THICKNESS`` property.
    density:
        Default Kratos ``DENSITY`` property.
    cell_domains:
        Optional list of :class:`CellDomainSpec` for multi-cell shared fields.
    extension_T:
        Optional trimming extension matrix ``E`` composed with ``C_THB``.
    n_gauss:
        Gauss quadrature order per axis for IGA/THB cells (defaults to ``degree + 1``).
    fe_block_names:
        Optional subset of :attr:`MeshModel.blocks` names to populate.
    """

    name: str
    kind: Literal["iga", "thb", "native_fe"]
    box: BoxSpline | None = None
    thb: THBHierarchy | None = None
    mesh: MeshModel | None = None
    young_modulus: float = 1e5
    poisson_ratio: float = 0.3
    thickness: float = 1.0
    density: float = 0.0
    cell_domains: list[CellDomainSpec] = field(default_factory=list)
    extension_T: sp.spmatrix | None = None
    n_gauss: int | None = None
    fe_block_names: list[str] | None = None


@dataclass
class InstanceSpec:
    """Placement and mesh-ownership specification for one part instance.

    Parameters
    ----------
    part_name:
        Key matching :attr:`AssemblySpec.parts`.
    instance_name:
        Unique instance suffix within the part (full submodelpart ``<part>_<instance>``).
    dependent:
        ``True`` to share the parent part's template grid/operators and only
        rigidly transform coordinates; ``False`` for an independent per-instance
        discretization.
    translation:
        Global translation vector ``t = (tx, ty, tz)``.
    rotation_matrix:
        Optional ``(3, 3)`` proper rotation matrix ``R`` applied before ``t``
        (``x' = R x + t``).
    young_modulus_override:
        Optional Young's modulus override for this instance.
    poisson_ratio_override:
        Optional Poisson's ratio override for this instance.
    box_override:
        Optional per-instance :class:`BoxSpline` (allowed only when ``dependent=False``).
    thb_override:
        Optional per-instance :class:`THBHierarchy` (allowed only when ``dependent=False``).
    mesh_override:
        Optional per-instance :class:`MeshModel` (allowed only when ``dependent=False``).
    """

    part_name: str
    instance_name: str
    dependent: bool = True
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation_matrix: np.ndarray | None = None
    young_modulus_override: float | None = None
    poisson_ratio_override: float | None = None
    box_override: BoxSpline | None = None
    thb_override: THBHierarchy | None = None
    mesh_override: MeshModel | None = None

    def full_name(self) -> str:
        """Return canonical ``<part_name>_<instance_name>`` identifier."""
        return f"{self.part_name}_{self.instance_name}"

    def rotation(self) -> np.ndarray:
        """Return validated ``(3, 3)`` rotation matrix."""
        return validate_rotation_matrix(self.rotation_matrix)

    def translation_vector(self) -> np.ndarray:
        """Return ``(3,)`` translation vector."""
        t = np.asarray(self.translation, dtype=np.float64)
        if t.shape != (3,):
            raise ValueError(f"translation must have shape (3,), got {t.shape}")
        return t

    def transform_points(self, points: np.ndarray) -> np.ndarray:
        """Transform local coordinates ``(..., 3)`` to global ``x' = R x + t``."""
        pts = np.asarray(points, dtype=np.float64)
        R = self.rotation()
        t = self.translation_vector()
        return pts @ R.T + t

    def transform_vector(self, vec: np.ndarray) -> np.ndarray:
        """Rotate local vector quantity ``(..., 3)`` to global ``v' = R v``."""
        v = np.asarray(vec, dtype=np.float64)
        R = self.rotation()
        return v @ R.T

    def inverse_transform_vector(self, vec: np.ndarray) -> np.ndarray:
        """Rotate global vector quantity ``(..., 3)`` back to local ``v = R^T v'``."""
        v = np.asarray(vec, dtype=np.float64)
        R = self.rotation()
        return v @ R


@dataclass
class AssemblySpec:
    """Top-level multi-part assembly specification.

    Parameters
    ----------
    job_name:
        Name used for metadata and ``<job_name>.resolved.json`` output.
    parts:
        Mapping of part name to :class:`PartSpec`.
    instances:
        Ordered list of :class:`InstanceSpec` to instantiate into ``Structure``.
    surfaces:
        Optional metadata dictionary for named surfaces under ``Surfaces.<name>``.
    couplings:
        Optional list of coupling specifications (e.g., Tie, RBE2, RBE3).
    use_pardiso:
        Prefer ``PardisoLUSolver`` when solving the assembled system.
    """

    job_name: str
    parts: dict[str, PartSpec]
    instances: list[InstanceSpec]
    surfaces: dict[str, dict[str, Any]] = field(default_factory=dict)
    couplings: list[Any] = field(default_factory=list)
    use_pardiso: bool = False


@dataclass(frozen=True)
class IdRangeEntry:
    """Single contiguous inclusive ID allocation record.

    Parameters
    ----------
    entity_type:
        One of ``"nodes"``, ``"elements"``, ``"conditions"``, ``"properties"``,
        or ``"constraints"``.
    scope:
        Owning scope name (e.g., ``"Parts.beam_left"``, ``"Couplings.tie1"``).
    start_id:
        First inclusive ID in range (1-based).
    end_id:
        Last inclusive ID in range (``start_id + count - 1``).
    count:
        Number of IDs allocated.
    provenance:
        Semantic role or element/constraint classification string.
    """

    entity_type: str
    scope: str
    start_id: int
    end_id: int
    count: int
    provenance: str

    def to_dict(self) -> dict[str, Any]:
        """Serialize entry to JSON-safe dictionary."""
        return {
            "entity_type": self.entity_type,
            "scope": self.scope,
            "start_id": self.start_id,
            "end_id": self.end_id,
            "count": self.count,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "IdRangeEntry":
        """Reconstruct :class:`IdRangeEntry` from dictionary."""
        return cls(
            entity_type=str(data["entity_type"]),
            scope=str(data["scope"]),
            start_id=int(data["start_id"]),
            end_id=int(data["end_id"]),
            count=int(data["count"]),
            provenance=str(data["provenance"]),
        )


class IdOffsetTable:
    """Monotonic non-overlapping ID allocator for Kratos assembly entities.

    Tracks separate 1-based monotonic counters for ``nodes``, ``elements``,
    ``conditions``, ``properties``, and ``constraints``, validates zero range
    overlaps, and serializes cleanly to UTF-8 JSON (``<job>.resolved.json``).
    """

    ENTITY_TYPES = ("nodes", "elements", "conditions", "properties", "constraints")

    def __init__(self, job_name: str = "assembly") -> None:
        """Initialize an empty table with all ID counters starting at ``1``.

        Parameters
        ----------
        job_name:
            Job identifier stored alongside serialized ID ranges.
        """
        self.job_name = job_name
        self._next_id: dict[str, int] = {k: 1 for k in self.ENTITY_TYPES}
        self.entries: list[IdRangeEntry] = []
        self.metadata: dict[str, Any] = {}

    def _allocate(
        self,
        entity_type: str,
        scope: str,
        count: int,
        provenance: str,
    ) -> IdRangeEntry:
        """Allocate a contiguous block of ``count`` IDs for ``entity_type``."""
        if entity_type not in self._next_id:
            raise KeyError(f"Unknown entity_type '{entity_type}'")
        if count <= 0:
            raise ValueError(f"Allocation count must be positive, got {count}")
        start = self._next_id[entity_type]
        end = start + count - 1
        self._next_id[entity_type] = end + 1
        entry = IdRangeEntry(
            entity_type=entity_type,
            scope=scope,
            start_id=start,
            end_id=end,
            count=count,
            provenance=provenance,
        )
        self.entries.append(entry)
        return entry

    def allocate_nodes(
        self,
        scope: str,
        count: int,
        provenance: _ProvenanceKind = "control_point",
    ) -> IdRangeEntry:
        """Allocate ``count`` monotonic node IDs with provenance metadata."""
        return self._allocate("nodes", scope, count, provenance)

    def allocate_elements(
        self,
        scope: str,
        count: int,
        kind: str = "SmallDisplacementElement3D4N",
    ) -> IdRangeEntry:
        """Allocate ``count`` monotonic element IDs."""
        return self._allocate("elements", scope, count, kind)

    def allocate_conditions(
        self,
        scope: str,
        count: int,
        kind: str = "SurfaceCondition3D3N",
    ) -> IdRangeEntry:
        """Allocate ``count`` monotonic condition IDs."""
        return self._allocate("conditions", scope, count, kind)

    def allocate_properties(
        self,
        scope: str,
        count: int = 1,
        provenance: str = "material",
    ) -> IdRangeEntry:
        """Allocate ``count`` monotonic Kratos ``Properties`` IDs."""
        return self._allocate("properties", scope, count, provenance)

    def allocate_constraints(
        self,
        scope: str,
        count: int,
        provenance: str = "LinearMasterSlaveConstraint",
    ) -> IdRangeEntry:
        """Allocate ``count`` monotonic Kratos constraint IDs."""
        return self._allocate("constraints", scope, count, provenance)

    def next_id(self, entity_type: str) -> int:
        """Return next free ID for ``entity_type`` without advancing it."""
        return self._next_id[entity_type]

    def validate_non_overlapping(self) -> None:
        """Assert that all allocated ID ranges per entity type are strictly disjoint."""
        for etype in self.ENTITY_TYPES:
            ranges = sorted(
                (e for e in self.entries if e.entity_type == etype),
                key=lambda e: e.start_id,
            )
            prev_end = 0
            for r in ranges:
                if r.start_id <= prev_end or r.end_id < r.start_id:
                    raise ValueError(
                        f"Overlapping or non-monotonic {etype} ID range: {r}"
                    )
                prev_end = r.end_id

    def to_dict(self) -> dict[str, Any]:
        """Serialize ID offset table and counters to a JSON-compatible dict."""
        self.validate_non_overlapping()
        return {
            "schema": "wht_iga.resolved_assembly.v1",
            "job_name": self.job_name,
            "next_ids": dict(self._next_id),
            "ranges": [e.to_dict() for e in self.entries],
            "metadata": self.metadata,
        }

    def to_json(self, path: str | Path) -> Path:
        """Write resolved ID table to ``<job>.resolved.json`` with UTF-8 encoding."""
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(self.to_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        return out_path

    @classmethod
    def from_json(cls, path: str | Path) -> "IdOffsetTable":
        """Read resolved ID table from UTF-8 JSON file and validate disjointness."""
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        table = cls(job_name=str(data.get("job_name", "assembly")))
        table.entries = [IdRangeEntry.from_dict(item) for item in data.get("ranges", [])]
        table.metadata = dict(data.get("metadata", {}))
        for etype in cls.ENTITY_TYPES:
            stored_next = data.get("next_ids", {}).get(etype, 1)
            max_end = max(
                (e.end_id for e in table.entries if e.entity_type == etype),
                default=0,
            )
            table._next_id[etype] = max(int(stored_next), max_end + 1)
        table.validate_non_overlapping()
        return table


def validate_assembly_spec(spec: AssemblySpec) -> None:
    """Validate part references, instance uniqueness, rotation matrices, and ownership.

    Enforces:
    * Every instance references a defined :class:`PartSpec`.
    * No duplicate instance names (``<part>_<instance>`` or ``instance_name``).
    * All rotation matrices are proper orthogonal (``det=+1``).
    * Part mesh ownership rule: within any part, instances must either ALL be
      ``dependent=True`` or ALL be ``dependent=False``. Dependent instances may
      not specify per-instance mesh/grid overrides.
    """
    seen_full: set[str] = set()
    ownership: dict[str, set[bool]] = {}
    for inst in spec.instances:
        if inst.part_name not in spec.parts:
            raise KeyError(
                f"Instance '{inst.instance_name}' references unknown part '{inst.part_name}'"
            )
        full = inst.full_name()
        if full in seen_full:
            raise ValueError(f"Duplicate instance identifier '{full}'")
        seen_full.add(full)
        validate_rotation_matrix(inst.rotation_matrix)
        inst.translation_vector()
        ownership.setdefault(inst.part_name, set()).add(bool(inst.dependent))
        if inst.dependent and (
            inst.box_override is not None
            or inst.thb_override is not None
            or inst.mesh_override is not None
        ):
            raise ValueError(
                f"Dependent instance '{full}' cannot override template grid/mesh"
            )
        part = spec.parts[inst.part_name]
        if part.kind == "iga" and part.box is None and inst.box_override is None:
            raise ValueError(f"IGA part '{part.name}' requires a BoxSpline")
        if part.kind == "thb" and part.thb is None and inst.thb_override is None:
            raise ValueError(f"THB part '{part.name}' requires a THBHierarchy")
        if part.kind == "native_fe" and part.mesh is None and inst.mesh_override is None:
            raise ValueError(f"Native FE part '{part.name}' requires a MeshModel")

    for pname, flags in ownership.items():
        if len(flags) > 1:
            raise ValueError(
                f"Part '{pname}' mixes dependent and independent instances, "
                "violating single mesh-ownership policy"
            )


@dataclass
class InstanceRuntime:
    """Runtime handles created for a single instantiated part instance."""

    spec: InstanceSpec
    part_spec: PartSpec
    sub_model_part: Any
    node_ids: list[int]
    element_ids: list[int]
    properties_map: dict[int, Any]
    box: BoxSpline | None = None
    thb: THBHierarchy | None = None
    volume: Any | None = None
    T_basis: sp.csr_matrix | None = None
    thb_master_node_ids: list[int] | None = None
    fe_handle: Any | None = None
    keepalive: list[Any] = field(default_factory=list)


@dataclass
class BuiltAssembly:
    """Complete assembled Kratos ``ModelPart("Structure")`` and runtime metadata."""

    model: Any
    structure: Any
    parts_smp: Any
    surfaces_smp: Any
    couplings_smp: Any
    id_table: IdOffsetTable
    instances: dict[str, InstanceRuntime]
    keepalive: list[Any] = field(default_factory=list)

    def get_instance(self, name: str) -> InstanceRuntime:
        """Retrieve :class:`InstanceRuntime` by ``<part>_<instance>`` or ``instance_name``."""
        if name in self.instances:
            return self.instances[name]
        matches = [
            inst for inst in self.instances.values() if inst.spec.instance_name == name
        ]
        if len(matches) == 1:
            return matches[0]
        raise KeyError(f"Unknown instance '{name}'")

    def save_resolved_json(self, path: str | Path) -> Path:
        """Serialize the assembly's ID offset table and instance metadata to UTF-8 JSON."""
        return self.id_table.to_json(path)


def _create_properties(
    KM: Any,
    SM: Any,
    model_part: Any,
    prop_id: int,
    young_modulus: float,
    poisson_ratio: float,
    thickness: float = 1.0,
    density: float = 0.0,
) -> Any:
    """Create a linear-elastic 3D Kratos ``Properties`` object."""
    props = model_part.CreateNewProperties(int(prop_id))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())
    props.SetValue(KM.YOUNG_MODULUS, float(young_modulus))
    props.SetValue(KM.POISSON_RATIO, float(poisson_ratio))
    props.SetValue(KM.THICKNESS, float(thickness))
    if density > 0.0:
        props.SetValue(KM.DENSITY, float(density))
    return props


def _build_nurbs_volume_transformed(
    KM: Any,
    sub_model_part: Any,
    box: BoxSpline,
    inst: InstanceSpec,
    first_node_id: int,
    add_dofs: bool = True,
) -> tuple[Any, Any, list[int]]:
    """Build Kratos ``NurbsVolumeGeometry`` with rigidly transformed control points.

    Parametric knots stay untouched in the part-local space while Greville
    control-point physical coordinates undergo ``x' = R x + t``.
    """
    from wht_iga.grid.knots import to_kratos

    ku_full, kv_full, kw_full = box.knots
    p = int(box.degree)
    nu, nv, nw = (int(v) for v in box.n_cp)

    cp_local = box.control_points()
    cp_global = inst.transform_points(cp_local)

    nodes_vec = KM.NodesVector()
    node_ids: list[int] = []
    nid = int(first_node_id)
    for k in range(nw):
        for j in range(nv):
            for i in range(nu):
                flat = i + nu * (j + nv * k)
                x, y, z = cp_global[flat]
                node = sub_model_part.CreateNewNode(nid, float(x), float(y), float(z))
                if add_dofs:
                    node.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
                    node.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
                    node.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)
                nodes_vec.append(node)
                node_ids.append(nid)
                nid += 1

    vol = KM.NurbsVolumeGeometry(
        nodes_vec,
        p,
        p,
        p,
        KM.Vector(to_kratos(ku_full, p).tolist()),
        KM.Vector(to_kratos(kv_full, p).tolist()),
        KM.Vector(to_kratos(kw_full, p).tolist()),
    )
    sub_model_part.AddGeometry(vol)
    return vol, nodes_vec, node_ids


def _cell_gauss_subbox(
    bounds_u: tuple[float, float],
    bounds_v: tuple[float, float],
    bounds_w: tuple[float, float],
    n_gauss_u: int,
    n_gauss_v: int,
    n_gauss_w: int,
) -> list[list[float]]:
    """Build tensor-product Legendre Gauss quadrature points ``[u, v, w, weight]``."""
    u0, u1 = bounds_u
    v0, v1 = bounds_v
    w0, w1 = bounds_w
    gxu, gwu = np.polynomial.legendre.leggauss(n_gauss_u)
    gxv, gwv = np.polynomial.legendre.leggauss(n_gauss_v)
    gxw, gww = np.polynomial.legendre.leggauss(n_gauss_w)

    pts_u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gxu
    wts_u = 0.5 * (u1 - u0) * gwu
    pts_v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gxv
    wts_v = 0.5 * (v1 - v0) * gwv
    pts_w = 0.5 * (w0 + w1) + 0.5 * (w1 - w0) * gxw
    wts_w = 0.5 * (w1 - w0) * gww

    out: list[list[float]] = []
    for wk, zw in zip(wts_w, pts_w):
        for wj, yv in zip(wts_v, pts_v):
            for wi, xu in zip(wts_u, pts_u):
                out.append([float(xu), float(yv), float(zw), float(wi * wj * wk)])
    return out


def _subdivide_box_by_domains(
    cell_bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    domains: list[CellDomainSpec],
    tol: float = 1e-11,
) -> list[tuple[CellDomainSpec, tuple[tuple[float, float], tuple[float, float], tuple[float, float]]]]:
    """Partition an axis-aligned spline cell against ``CellDomainSpec.bounds``."""
    (cx0, cx1), (cy0, cy1), (cz0, cz1) = cell_bounds
    center = np.array([0.5 * (cx0 + cx1), 0.5 * (cy0 + cy1), 0.5 * (cz0 + cz1)])

    # Collect explicit non-straddling axis-aligned overlaps first.
    axis_overlaps: list[
        tuple[CellDomainSpec, tuple[tuple[float, float], tuple[float, float], tuple[float, float]]]
    ] = []
    for dom in domains:
        if dom.bounds is None:
            continue
        (dx0, dx1), (dy0, dy1), (dz0, dz1) = dom.bounds
        ox0, ox1 = max(cx0, dx0), min(cx1, dx1)
        oy0, oy1 = max(cy0, dy0), min(cy1, dy1)
        oz0, oz1 = max(cz0, dz0), min(cz1, dz1)
        if ox1 - ox0 > tol and oy1 - oy0 > tol and oz1 - oz0 > tol:
            axis_overlaps.append((dom, ((ox0, ox1), (oy0, oy1), (oz0, oz1))))
    if axis_overlaps:
        return axis_overlaps

    # Fallback: predicate check at cell center.
    for dom in domains:
        if dom.contains_point(center, tol=tol):
            return [(dom, cell_bounds)]
    return []


class AssemblyBuilder:
    """Compiles an :class:`AssemblySpec` into a unified Kratos ``ModelPart("Structure")``.

    Supports IGA, multi-cell shared B-spline fields (`T-P3.1` Path B),
    multi-level THB hierarchies (`T-P3.1` Path A), dependent vs. independent
    part instances (`T-P3.2`), and Native FE parts (`T-P3.3`).
    """

    def __init__(self, KM: Any | None = None, SM: Any | None = None) -> None:
        """Initialize builder with optional pre-loaded Kratos modules.

        Parameters
        ----------
        KM:
            Optional ``KratosMultiphysics`` module handle. Loaded lazily if omitted.
        SM:
            Optional ``StructuralMechanicsApplication`` handle. Loaded lazily if omitted.
        """
        if KM is None or SM is None:
            from wht_iga.kratos_backend import load_kratos

            apps = load_kratos()
            KM = apps["KratosMultiphysics"]
            SM = apps["StructuralMechanicsApplication"]
        self.KM = KM
        self.SM = SM

    @staticmethod
    def _ensure_nodal_variables(KM: Any, SM: Any, structure: Any) -> None:
        """Register standard structural nodal solution-step variables on ``Structure``."""
        for var in (
            KM.DISPLACEMENT,
            KM.REACTION,
            KM.ROTATION,
            KM.REACTION_MOMENT,
            SM.POINT_LOAD,
            SM.POINT_MOMENT,
            SM.SURFACE_LOAD,
            KM.VOLUME_ACCELERATION,
        ):
            if not structure.HasNodalSolutionStepVariable(var):
                structure.AddNodalSolutionStepVariable(var)

    def build(self, spec: AssemblySpec) -> BuiltAssembly:
        """Instantiate all parts, submodelparts, properties, and ID ranges.

        Parameters
        ----------
        spec:
            Validated multi-part assembly specification.

        Returns
        -------
        BuiltAssembly
            Root ``Structure`` model part, ID table, and per-instance handles.
        """
        validate_assembly_spec(spec)
        KM, SM = self.KM, self.SM

        model = KM.Model()
        structure = model.CreateModelPart("Structure")
        self._ensure_nodal_variables(KM, SM, structure)
        parts_smp = structure.CreateSubModelPart("Parts")
        surfaces_smp = structure.CreateSubModelPart("Surfaces")
        couplings_smp = structure.CreateSubModelPart("Couplings")

        id_table = IdOffsetTable(job_name=spec.job_name)
        instances: dict[str, InstanceRuntime] = {}
        keepalive: list[Any] = []

        # Cache properties for shared template dependent instances or create per-instance.
        shared_props_cache: dict[str, dict[int, Any]] = {}

        for inst in spec.instances:
            part = spec.parts[inst.part_name]
            scope = f"Parts.{inst.full_name()}"
            inst_smp = parts_smp.CreateSubModelPart(inst.full_name())

            e_default = (
                inst.young_modulus_override
                if inst.young_modulus_override is not None
                else part.young_modulus
            )
            nu_default = (
                inst.poisson_ratio_override
                if inst.poisson_ratio_override is not None
                else part.poisson_ratio
            )

            # Reuse Properties across dependent instances unless instance overrides material.
            use_shared_props = (
                inst.dependent
                and inst.young_modulus_override is None
                and inst.poisson_ratio_override is None
                and part.name in shared_props_cache
            )
            if use_shared_props:
                props_map = shared_props_cache[part.name]
            else:
                props_map = {}
                if part.cell_domains:
                    for dom in part.cell_domains:
                        pid_entry = id_table.allocate_properties(
                            scope, 1, provenance=f"cell_{dom.cell_id}_{dom.name}"
                        )
                        props_map[dom.cell_id] = _create_properties(
                            KM,
                            SM,
                            structure,
                            pid_entry.start_id,
                            dom.young_modulus,
                            dom.poisson_ratio,
                            dom.thickness,
                            dom.density,
                        )
                else:
                    pid_entry = id_table.allocate_properties(scope, 1, provenance="part_material")
                    props_map[0] = _create_properties(
                        KM,
                        SM,
                        structure,
                        pid_entry.start_id,
                        e_default,
                        nu_default,
                        part.thickness,
                        part.density,
                    )
                if inst.dependent and inst.young_modulus_override is None and inst.poisson_ratio_override is None:
                    shared_props_cache[part.name] = props_map

            if part.kind == "iga":
                runtime = self._build_iga_instance(
                    part=part,
                    inst=inst,
                    scope=scope,
                    inst_smp=inst_smp,
                    props_map=props_map,
                    id_table=id_table,
                )
            elif part.kind == "thb":
                runtime = self._build_thb_instance(
                    part=part,
                    inst=inst,
                    scope=scope,
                    inst_smp=inst_smp,
                    props_map=props_map,
                    id_table=id_table,
                )
            elif part.kind == "native_fe":
                runtime = self._build_fe_instance(
                    part=part,
                    inst=inst,
                    scope=scope,
                    inst_smp=inst_smp,
                    props_map=props_map,
                    id_table=id_table,
                )
            else:
                raise ValueError(f"Unsupported part kind '{part.kind}'")

            keepalive.extend(runtime.keepalive)
            instances[inst.full_name()] = runtime

        id_table.metadata["parts"] = {
            name: {"kind": p.kind, "n_domains": len(p.cell_domains)}
            for name, p in spec.parts.items()
        }
        id_table.metadata["instances"] = {
            full: {
                "part_name": rt.spec.part_name,
                "instance_name": rt.spec.instance_name,
                "dependent": rt.spec.dependent,
                "translation": list(rt.spec.translation),
                "rotation_matrix": rt.spec.rotation().tolist(),
                "n_nodes": len(rt.node_ids),
                "n_elements": len(rt.element_ids),
            }
            for full, rt in instances.items()
        }
        id_table.validate_non_overlapping()

        return BuiltAssembly(
            model=model,
            structure=structure,
            parts_smp=parts_smp,
            surfaces_smp=surfaces_smp,
            couplings_smp=couplings_smp,
            id_table=id_table,
            instances=instances,
            keepalive=keepalive,
        )

    def _build_iga_instance(
        self,
        part: PartSpec,
        inst: InstanceSpec,
        scope: str,
        inst_smp: Any,
        props_map: dict[int, Any],
        id_table: IdOffsetTable,
    ) -> InstanceRuntime:
        """Instantiate an IGA BoxSpline instance (single or multi-cell shared field)."""
        KM = self.KM
        box = inst.box_override if (not inst.dependent and inst.box_override is not None) else part.box
        assert box is not None
        node_entry = id_table.allocate_nodes(scope, box.n_ctrl, provenance="control_point")
        vol, nodes_vec, node_ids = _build_nurbs_volume_transformed(
            KM,
            inst_smp,
            box,
            inst,
            first_node_id=node_entry.start_id,
            add_dofs=True,
        )

        bu = box.breakpoints(0)
        bv = box.breakpoints(1)
        bw = box.breakpoints(2)
        nx, ny, nz = box.n_cells
        p = int(box.degree)
        ng_u = part.n_gauss or (p + 1)
        ng_v = part.n_gauss or (p + 1)
        ng_w = part.n_gauss or (p + 1)

        elem_ids: list[int] = []
        created_qps: list[Any] = []
        next_eid = id_table.next_id("elements")
        cell_specs_to_add: list[tuple[Any, list[list[float]]]] = []

        # Map explicit cell assignments if specified.
        explicit_map: dict[tuple[int, int, int], CellDomainSpec] = {}
        for dom in part.cell_domains:
            if dom.cells is not None:
                for cidx in dom.cells:
                    explicit_map[tuple(cidx)] = dom

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    cbounds = (
                        (float(bu[i]), float(bu[i + 1])),
                        (float(bv[j]), float(bv[j + 1])),
                        (float(bw[k]), float(bw[k + 1])),
                    )
                    if not part.cell_domains:
                        pts = _cell_gauss_subbox(*cbounds, ng_u, ng_v, ng_w)
                        cell_specs_to_add.append((props_map[0], pts))
                    elif (i, j, k) in explicit_map:
                        dom = explicit_map[(i, j, k)]
                        pts = _cell_gauss_subbox(*cbounds, ng_u, ng_v, ng_w)
                        cell_specs_to_add.append((props_map[dom.cell_id], pts))
                    else:
                        sub_domains = _subdivide_box_by_domains(cbounds, part.cell_domains)
                        for dom, sub_bounds in sub_domains:
                            pts = _cell_gauss_subbox(*sub_bounds, ng_u, ng_v, ng_w)
                            cell_specs_to_add.append((props_map[dom.cell_id], pts))

        if cell_specs_to_add:
            elem_entry = id_table.allocate_elements(
                scope, len(cell_specs_to_add), kind="SmallDisplacementElement3D4N"
            )
            next_eid = elem_entry.start_id
            for props, pts in cell_specs_to_add:
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                if len(qps) != 1:
                    raise RuntimeError(
                        f"Expected 1 QP geometry inside cell subbox, got {len(qps)}"
                    )
                inst_smp.CreateNewElement("SmallDisplacementElement3D4N", next_eid, qps[0], props)
                elem_ids.append(next_eid)
                created_qps.append(qps[0])
                next_eid += 1

        T_basis = compose_thb_operator(part.extension_T, None)
        return InstanceRuntime(
            spec=inst,
            part_spec=part,
            sub_model_part=inst_smp,
            node_ids=node_ids,
            element_ids=elem_ids,
            properties_map=props_map,
            box=box,
            volume=vol,
            T_basis=T_basis,
            keepalive=[vol, nodes_vec, *created_qps],
        )

    def _build_thb_instance(
        self,
        part: PartSpec,
        inst: InstanceSpec,
        scope: str,
        inst_smp: Any,
        props_map: dict[int, Any],
        id_table: IdOffsetTable,
    ) -> InstanceRuntime:
        """Instantiate a multi-level THB instance integrated over leaf cells."""
        KM = self.KM
        thb = inst.thb_override if (not inst.dependent and inst.thb_override is not None) else part.thb
        assert thb is not None

        L_max = thb.max_level
        nx0, ny0, nz0 = thb.cell_shape(0)
        factor = 2 ** L_max
        pu, pv, pw = thb.degrees
        if pu != pv or pu != pw:
            raise ValueError("THB BoxSpline builder requires isotropic degree pu==pv==pw")

        ku_max = thb._knots_u[L_max]
        kv_max = thb._knots_v[L_max]
        kw_max = thb._knots_w[L_max]
        box_fine = BoxSpline(
            lo=(float(ku_max[0]), float(kv_max[0]), float(kw_max[0])),
            hi=(float(ku_max[-1]), float(kv_max[-1]), float(kw_max[-1])),
            n_cells=(nx0 * factor, ny0 * factor, nz0 * factor),
            degree=pu,
        )

        C_thb = thb.build_c_thb_matrix()
        T_basis = compose_thb_operator(part.extension_T, C_thb)

        node_entry = id_table.allocate_nodes(scope, box_fine.n_ctrl, provenance="control_point")
        vol, nodes_vec, node_ids = _build_nurbs_volume_transformed(
            KM,
            inst_smp,
            box_fine,
            inst,
            first_node_id=node_entry.start_id,
            add_dofs=True,
        )

        # Fine knot breakpoints inside leaf cells ensure exact spline polynomial
        # integration even in transition cells where fine THB functions enter coarse cells.
        bu_fine = box_fine.breakpoints(0)
        bv_fine = box_fine.breakpoints(1)
        bw_fine = box_fine.breakpoints(2)
        ng = part.n_gauss or (pu + 1)

        def _subspans(bps: np.ndarray, lo: float, hi: float) -> list[tuple[float, float]]:
            tol = 1e-12 * max(1.0, abs(hi - lo))
            interior = [float(b) for b in bps if lo + tol < float(b) < hi - tol]
            pts = [float(lo), *interior, float(hi)]
            return [(pts[m], pts[m + 1]) for m in range(len(pts) - 1)]

        subcell_queue: list[tuple[Any, list[list[float]]]] = []
        for lc in thb.leaf_cells():
            su_list = _subspans(bu_fine, *lc.bounds_u)
            sv_list = _subspans(bv_fine, *lc.bounds_v)
            sw_list = _subspans(bw_fine, *lc.bounds_w)
            for sw in sw_list:
                for sv in sv_list:
                    for su in su_list:
                        if not part.cell_domains:
                            subcell_queue.append(
                                (props_map[0], _cell_gauss_subbox(su, sv, sw, ng, ng, ng))
                            )
                        else:
                            for dom, sbox in _subdivide_box_by_domains((su, sv, sw), part.cell_domains):
                                subcell_queue.append(
                                    (
                                        props_map[dom.cell_id],
                                        _cell_gauss_subbox(sbox[0], sbox[1], sbox[2], ng, ng, ng),
                                    )
                                )

        elem_ids: list[int] = []
        created_qps: list[Any] = []
        if subcell_queue:
            elem_entry = id_table.allocate_elements(
                scope, len(subcell_queue), kind="SmallDisplacementElement3D4N"
            )
            eid = elem_entry.start_id
            for props, pts in subcell_queue:
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                if len(qps) != 1:
                    raise RuntimeError(
                        f"Expected 1 QP geometry per THB integration subspan, got {len(qps)}"
                    )
                inst_smp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                elem_ids.append(eid)
                created_qps.append(qps[0])
                eid += 1

        return InstanceRuntime(
            spec=inst,
            part_spec=part,
            sub_model_part=inst_smp,
            node_ids=node_ids,
            element_ids=elem_ids,
            properties_map=props_map,
            box=box_fine,
            thb=thb,
            volume=vol,
            T_basis=T_basis,
            keepalive=[vol, nodes_vec, *created_qps],
        )

    def _build_fe_instance(
        self,
        part: PartSpec,
        inst: InstanceSpec,
        scope: str,
        inst_smp: Any,
        props_map: dict[int, Any],
        id_table: IdOffsetTable,
    ) -> InstanceRuntime:
        """Instantiate a Native FE MeshModel instance via :mod:`native_fe`."""
        from wht_iga.kratos_backend.native_fe import populate_native_fe_part

        mesh = inst.mesh_override if (not inst.dependent and inst.mesh_override is not None) else part.mesh
        assert mesh is not None
        default_props = props_map[0] if 0 in props_map else next(iter(props_map.values()))
        handle = populate_native_fe_part(
            KM=self.KM,
            sub_model_part=inst_smp,
            mesh=mesh,
            props=default_props,
            id_table=id_table,
            scope=scope,
            node_transform=inst.transform_points,
            block_names=part.fe_block_names,
            props_by_part_id=props_map if part.cell_domains else None,
        )
        return InstanceRuntime(
            spec=inst,
            part_spec=part,
            sub_model_part=inst_smp,
            node_ids=handle.node_ids,
            element_ids=handle.element_ids,
            properties_map=props_map,
            fe_handle=handle,
        )

    def attach_thb_master_constraints(
        self,
        built: BuiltAssembly,
        instance_name: str,
        tol: float = 1e-13,
    ) -> list[int]:
        """Create THB master virtual nodes and wire ``u_fine = T @ u_red`` in Kratos.

        Allocates ``n_red`` virtual master nodes under ``Parts.<inst>_thb_masters``
        and creates zero-constant ``LinearMasterSlaveConstraint``s on every fine
        control-point node so Kratos ``ResidualBasedBlockBuilderAndSolver`` solves
        directly in the reduced THB (or extended THB) basis.

        Parameters
        ----------
        built:
            Assembled :class:`BuiltAssembly`.
        instance_name:
            Key of the THB instance inside ``built``.
        tol:
            Minimum absolute weight included in sparse MPC entries.

        Returns
        -------
        list[int]
            Allocated Kratos node IDs for the reduced THB master DOFs.
        """
        from wht_iga.kratos_backend.constraints import add_linear_tie

        KM = self.KM
        rt = built.get_instance(instance_name)
        if rt.T_basis is None:
            raise ValueError(f"Instance '{instance_name}' has no THB/extension T_basis")
        T = rt.T_basis.tocsr()
        n_fine, n_red = T.shape
        if n_fine != len(rt.node_ids):
            raise ValueError(
                f"T_basis rows ({n_fine}) != fine nodes ({len(rt.node_ids)})"
            )

        scope = f"Parts.{rt.spec.full_name()}.thb_masters"
        master_smp = rt.sub_model_part.CreateSubModelPart("thb_masters")
        master_entry = built.id_table.allocate_nodes(
            scope, n_red, provenance="thb_master"
        )
        master_ids: list[int] = []
        master_nodes: list[Any] = []
        # Place virtual master nodes at weighted average coordinates for clean inspection.
        fine_coords = np.array(
            [[built.structure.GetNode(nid).X0,
              built.structure.GetNode(nid).Y0,
              built.structure.GetNode(nid).Z0] for nid in rt.node_ids],
            dtype=np.float64,
        )
        col_sums = np.maximum(np.asarray(T.sum(axis=0)).ravel(), 1e-15)
        approx_coords = (T.T @ fine_coords) / col_sums[:, None]

        for j in range(n_red):
            nid = master_entry.start_id + j
            x, y, z = approx_coords[j]
            nd = master_smp.CreateNewNode(nid, float(x), float(y), float(z))
            nd.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
            nd.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
            nd.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)
            master_ids.append(nid)
            master_nodes.append(nd)

        disp_comps = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
        cid_entry = built.id_table.allocate_constraints(
            scope, n_fine * 3, provenance="thb_prolongation"
        )
        cid = cid_entry.start_id
        for i, fine_nid in enumerate(rt.node_ids):
            row_start, row_end = T.indptr[i], T.indptr[i + 1]
            cols = T.indices[row_start:row_end]
            vals = T.data[row_start:row_end]
            keep = np.abs(vals) > tol
            cols = cols[keep]
            vals = vals[keep]
            if len(cols) == 0:
                raise ValueError(f"Fine control point row {i} has zero THB support")
            fine_node = built.structure.GetNode(fine_nid)
            rmat = vals.reshape(1, -1)
            for comp in disp_comps:
                m_dofs = [master_nodes[int(c)].GetDof(comp) for c in cols]
                s_dofs = [fine_node.GetDof(comp)]
                add_linear_tie(KM, built.structure, cid, s_dofs, m_dofs, rmat, constant=None)
                cid += 1

        rt.thb_master_node_ids = master_ids
        built.id_table.validate_non_overlapping()
        return master_ids

    def assemble_and_solve_thb_reduced(
        self,
        built: BuiltAssembly,
        instance_name: str,
        use_pardiso: bool = False,
    ) -> dict[str, Any]:
        """Assemble hierarchical ``(K_hier, f_hier)`` in Kratos and solve reduced THB system.

        Extracts Kratos global ``(K_hier, f_hier)`` matching ``rt.node_ids`` DOF order,
        computes ``K_red = T_dof^T K_hier T_dof`` and ``f_red = T_dof^T f_hier`` where
        ``T_dof = kron(T_basis, I_3)``, solves ``K_red u_red = f_red`` via
        :func:`wht_iga.operators.compose.solve_reduced` (after applying any penalty
        conditions or fixed fine Dirichlet values), expands ``u_hier = T_dof @ u_red``,
        and writes displacements back onto Kratos control-point nodes.

        Parameters
        ----------
        built:
            Assembled :class:`BuiltAssembly` with loads/Dirichlet/penalty BCs attached.
        instance_name:
            Name of the THB (or extended IGA) instance.
        use_pardiso:
            Ignored here when solving via dense/sparse LU in ``solve_reduced``.

        Returns
        -------
        dict[str, Any]
            Dictionary containing ``K_hier``, ``f_hier``, ``T_dof``, ``K_red``,
            ``f_red``, ``u_red``, and ``u_hier``.
        """
        from wht_iga.kratos_backend.elements import solve_linear, system_matrix_csr

        del use_pardiso
        KM = self.KM
        rt = built.get_instance(instance_name)
        if rt.T_basis is None:
            raise ValueError(f"Instance '{instance_name}' has no T_basis for THB reduction")

        # Assemble K_hier and RHS using Kratos BlockBuilder without constraints
        # if none attached, or assemble element-by-element + condition RHS directly.
        n_fine = len(rt.node_ids)
        T_dof = expand_dof_operator(rt.T_basis, n_comp=3)

        # Check whether Structure only contains this instance's fine nodes & no constraints:
        # assemble element & condition local contributions explicitly into fine DOF ordering.
        node_to_local = {nid: idx for idx, nid in enumerate(rt.node_ids)}
        ndof = 3 * n_fine
        K_coo_rows: list[int] = []
        K_coo_cols: list[int] = []
        K_coo_vals: list[float] = []
        f_hier = np.zeros(ndof, dtype=np.float64)

        disp_names = ("DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z")
        pi = built.structure.ProcessInfo
        for obj in (*built.structure.Elements, *built.structure.Conditions):
            obj.Initialize(pi)
            Kl = KM.Matrix()
            fr = KM.Vector()
            obj.CalculateLocalSystem(Kl, fr, pi)
            geom = obj.GetGeometry()
            local_dofs: list[int] = []
            for idx_n in range(geom.PointsNumber()):
                lidx = node_to_local.get(int(geom[idx_n].Id))
                if lidx is not None:
                    local_dofs.extend([3 * lidx, 3 * lidx + 1, 3 * lidx + 2])
            if len(local_dofs) == len(fr):
                Ka = np.array(Kl, dtype=np.float64)
                fa = np.array(fr, dtype=np.float64)
                if len(local_dofs) > 0:
                    idx_arr = np.asarray(local_dofs, dtype=np.int64)
                    rr, cc = np.meshgrid(idx_arr, idx_arr, indexing="ij")
                    K_coo_rows.extend(rr.ravel().tolist())
                    K_coo_cols.extend(cc.ravel().tolist())
                    K_coo_vals.extend(Ka.ravel().tolist())
                    f_hier[idx_arr] += fa

        # Add nodal POINT_LOAD if no PointLoadCondition covered them (already in Conditions above)
        # Apply Dirichlet fixed DOFs by identifying fixed DOFs on fine nodes:
        # Note: in THB/penalty formulations, SupportPenaltyCondition contributes directly to K_hier/f_hier.
        # If strong nodal Fix() was applied to fine DOFs, fold via symmetric Dirichlet elimination on K_red
        # ONLY for pure nodal values when compatible, or validate exact penalty/Dirichlet consistency.
        K_hier = sp.coo_matrix((K_coo_vals, (K_coo_rows, K_coo_cols)), shape=(ndof, ndof)).tocsr()
        K_red = apply_TKT(K_hier, T_dof)
        f_red = apply_Tt_vec(f_hier, T_dof)

        # Check if strong nodal Fix() was set on any fine node (e.g. SupportPenalty not used)
        fixed_fine: dict[int, float] = {}
        disp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
        for idx, nid in enumerate(rt.node_ids):
            nd = built.structure.GetNode(nid)
            for c, dv in enumerate(disp_vars):
                if nd.IsFixed(dv):
                    fixed_fine[3 * idx + c] = float(nd.GetSolutionStepValue(dv, 0))

        if fixed_fine:
            # Constrain all reduced DOFs having non-negligible support on fixed fine DOFs when u_presc==0
            fixed_indices = sorted(fixed_fine.keys())
            T_fixed = T_dof[fixed_indices, :]
            active_cols = np.where(np.max(np.abs(T_fixed.toarray()), axis=0) > 1e-12)[0]
            # For homogeneous boundary (all 0) or direct projection:
            u_fixed_vals = np.array([fixed_fine[i] for i in fixed_indices], dtype=np.float64)
            u_red_bc = np.zeros(T_dof.shape[1], dtype=np.float64)
            if np.any(u_fixed_vals != 0.0) and len(active_cols) > 0:
                T_sub = T_fixed[:, active_cols].toarray()
                u_sub, *_ = np.linalg.lstsq(T_sub, u_fixed_vals, rcond=None)
                u_red_bc[active_cols] = u_sub
            free_cols = np.setdiff1d(np.arange(T_dof.shape[1]), active_cols)
            rhs_mod = f_red - K_red @ u_red_bc
            u_red = u_red_bc.copy()
            if len(free_cols) > 0:
                K_ff = K_red[free_cols, :][:, free_cols]
                u_red[free_cols] = solve_reduced(K_ff, rhs_mod[free_cols])
        else:
            u_red = solve_reduced(K_red, f_red)

        u_hier = T_dof @ u_red
        for idx, nid in enumerate(rt.node_ids):
            nd = built.structure.GetNode(nid)
            nd.SetSolutionStepValue(KM.DISPLACEMENT_X, 0, float(u_hier[3 * idx]))
            nd.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0, float(u_hier[3 * idx + 1]))
            nd.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0, float(u_hier[3 * idx + 2]))

        return {
            "K_hier": K_hier,
            "f_hier": f_hier,
            "T_dof": T_dof,
            "K_red": K_red,
            "f_red": f_red,
            "u_red": u_red,
            "u_hier": u_hier,
        }


__all__ = [
    "AssemblyBuilder",
    "AssemblySpec",
    "BuiltAssembly",
    "CellDomainSpec",
    "IdOffsetTable",
    "IdRangeEntry",
    "InstanceRuntime",
    "InstanceSpec",
    "PartSpec",
    "compose_thb_operator",
    "transform_stress_voigt",
    "validate_assembly_spec",
    "validate_rotation_matrix",
]
