"""Selectable discretization option models (split from methods.py)."""

from __future__ import annotations

from typing import Annotated, Literal, Self, Union

from pydantic import Field, TypeAdapter, model_validator

from wht_iga.spec.base import PositiveFinite, SpecBase, Vec3


MethodId = Literal[
    "full_iga_body_fit",
    "full_iga_immersed",
    "thb_adaptive_iga",
    "iga_fe_hybrid",
    "native_fe",
]

FeElementType = Literal["hex8", "hex20", "tet4", "tet10", "wedge6"]
HybridCouplingKind = Literal["auto", "penalty_spring", "mortar_l2", "barycentric"]
LockingMitigation = Literal["auto", "none", "fbar"]
Degree = Annotated[int, Field(ge=1, le=3)]
CellCounts = tuple[
    Annotated[int, Field(ge=1)],
    Annotated[int, Field(ge=1)],
    Annotated[int, Field(ge=1)],
]


class MethodDescriptor(SpecBase):
    """Metadata describing one selectable discretization method.

    Attributes
    ----------
    id:
        Canonical method identifier literal.
    label:
        Human-readable English label for UI and reports.
    korean_summary:
        Concise Korean summary of the discretization strategy.
    uses_bspline:
        Whether the method uses B-spline or THB-spline basis functions.
    supports_adaptive_refinement:
        Whether local hierarchical adaptive refinement is supported.
    supports_fe_bulk:
        Whether standard Kratos FE elements are used in all or part of the domain.
    engine_module:
        Primary WHT_IGA engine module implementing the discretization.
    """

    id: MethodId
    label: str
    korean_summary: str
    uses_bspline: bool
    supports_adaptive_refinement: bool
    supports_fe_bulk: bool
    engine_module: str


_METHOD_DESCRIPTORS: tuple[MethodDescriptor, ...] = (
    MethodDescriptor(
        id="full_iga_body_fit",
        label="Full IGA (Body-Fit)",
        korean_summary="전 영역 경계적합(Body-Fit) B-스플라인 등기하해석",
        uses_bspline=True,
        supports_adaptive_refinement=False,
        supports_fe_bulk=False,
        engine_module="wht_iga.solver.pipeline",
    ),
    MethodDescriptor(
        id="full_iga_immersed",
        label="Full IGA (Immersed / Trimmed)",
        korean_summary="비적합 침지(Immersed/FCM) B-스플라인 등기하해석",
        uses_bspline=True,
        supports_adaptive_refinement=False,
        supports_fe_bulk=False,
        engine_module="wht_iga.solver.pipeline",
    ),
    MethodDescriptor(
        id="thb_adaptive_iga",
        label="THB Adaptive IGA",
        korean_summary="절단 계층 B-스플라인(THB) 국소 적응 세분화 등기하해석",
        uses_bspline=True,
        supports_adaptive_refinement=True,
        supports_fe_bulk=False,
        engine_module="wht_iga.grid.hierarchy",
    ),
    MethodDescriptor(
        id="iga_fe_hybrid",
        label="IGA-FE Hybrid",
        korean_summary="FEM 벌크와 정밀 IGA 도메인을 인터페이스 결합한 하이브리드 해석",
        uses_bspline=True,
        supports_adaptive_refinement=False,
        supports_fe_bulk=True,
        engine_module="wht_iga.kratos_backend.native_fe",
    ),
    MethodDescriptor(
        id="native_fe",
        label="Native FEM",
        korean_summary="표준 Kratos 구조 유한요소해석 (Hex8/Hex20/Tet4/Tet10/Wedge6)",
        uses_bspline=False,
        supports_adaptive_refinement=False,
        supports_fe_bulk=True,
        engine_module="wht_iga.kratos_backend.native_fe",
    ),
)


def list_discretization_methods() -> tuple[MethodDescriptor, ...]:
    """Return metadata descriptors for all five selectable discretization methods."""
    return _METHOD_DESCRIPTORS


class FullIgaOption(SpecBase):
    """Option model for body-fit or immersed full-domain B-spline IGA.

    Attributes
    ----------
    kind:
        ``"full_iga_body_fit"`` or ``"full_iga_immersed"``.
    degree:
        Spline polynomial degree (1 to 3).
    continuity:
        Inter-element continuity order; defaults to ``degree - 1``.
    n_cells:
        Grid cell counts ``(nx, ny, nz)``.
    quadrature_degree:
        Cut-cell quadrature degree for immersed mode; defaults to ``min(2*degree, 4)``.
    out_mode:
        Handling of outside cells in immersed mode (``"drop"`` or ``"ghost"``).
    ghost_alpha:
        Stiffness scaling factor for outside cells when ``out_mode == "ghost"``.
    bc_penalty:
        Optional explicit Dirichlet penalty parameter ``beta``.
    trim_deflection:
        Optional CAD tessellation chord deflection for immersed trimming.
    locking_mitigation:
        Volumetric locking mitigation for near-incompressible materials:
        ``"none"`` (default linear elasticity), ``"fbar"`` (F-bar mean-dilatation
        for body-fit boxes), or ``"auto"`` (uses F-bar when ``nu >= 0.49``).
    """

    kind: Literal["full_iga_body_fit", "full_iga_immersed"] = "full_iga_body_fit"
    degree: Degree = 2
    continuity: Annotated[int, Field(ge=0)] | None = None
    n_cells: CellCounts = (4, 2, 2)
    quadrature_degree: Annotated[int, Field(ge=1, le=6)] | None = None
    out_mode: Literal["drop", "ghost"] = "drop"
    ghost_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = 1e-8
    bc_penalty: PositiveFinite | None = None
    trim_deflection: PositiveFinite | None = None
    locking_mitigation: LockingMitigation = "none"

    @model_validator(mode="after")
    def _validate_continuity(self) -> Self:
        if self.continuity is not None and self.continuity >= self.degree:
            raise ValueError(
                f"continuity must be < degree, got continuity={self.continuity}, degree={self.degree}"
            )
        if self.locking_mitigation == "fbar" and self.kind != "full_iga_body_fit":
            raise ValueError(
                "locking_mitigation='fbar' is only supported for kind='full_iga_body_fit'"
            )
        return self

    @property
    def resolved_continuity(self) -> int:
        """Continuity order actually passed to the spline builder."""
        return self.degree - 1 if self.continuity is None else self.continuity

    @property
    def resolved_quadrature_degree(self) -> int:
        """Cut-cell moment quadrature degree used for immersed trimming."""
        return min(2 * self.degree, 4) if self.quadrature_degree is None else self.quadrature_degree


class ThbAdaptiveIgaOption(SpecBase):
    """Option model for locally refined Truncated Hierarchical B-spline (THB) IGA.

    Attributes
    ----------
    kind:
        Discriminator literal ``"thb_adaptive_iga"``.
    degree:
        Isotropic spline polynomial degree (1 to 3).
    coarse_n_cells:
        Level-0 coarse background cell counts ``(nx, ny, nz)``.
    refinement_level:
        Target hierarchical refinement depth (1 to 5).
    refine_box:
        Optional axis-aligned ``(lo, hi)`` box region to refine.
    refine_sphere:
        Optional ``(center, radius)`` spherical region to refine.
    balance:
        2:1 level balance enforcement policy (``"two_to_one"`` or ``"none"``).
    penalty_beta:
        Optional Dirichlet penalty factor; if ``None``, strong nodal clamping is used.
    """

    kind: Literal["thb_adaptive_iga"] = "thb_adaptive_iga"
    degree: Degree = 2
    coarse_n_cells: CellCounts = (2, 2, 4)
    refinement_level: Annotated[int, Field(ge=1, le=5)] = 1
    refine_box: tuple[Vec3, Vec3] | None = None
    refine_sphere: tuple[Vec3, PositiveFinite] | None = None
    balance: Literal["two_to_one", "none"] = "two_to_one"
    penalty_beta: PositiveFinite | None = None

    @model_validator(mode="after")
    def _validate_regions(self) -> Self:
        if self.refine_box is not None:
            r_lo, r_hi = self.refine_box
            if any(h <= l for l, h in zip(r_lo, r_hi)):
                raise ValueError(f"refine_box hi must exceed lo on every axis, got {self.refine_box}")
        return self


class NativeFeOption(SpecBase):
    """Option model for standard Kratos Native FEM solid discretization.

    Attributes
    ----------
    kind:
        Discriminator literal ``"native_fe"``.
    element_type:
        Canonical solid element topology (``hex8``, ``hex20``, ``tet4``, ``tet10``, ``wedge6``).
    n_cells:
        Structured grid subdivision counts ``(nx, ny, nz)``.
    """

    kind: Literal["native_fe"] = "native_fe"
    element_type: FeElementType = "hex20"
    n_cells: CellCounts = (4, 2, 2)


class IgaFeHybridOption(SpecBase):
    """Option model for coupled Native FEM + IGA hybrid domain discretization.

    Attributes
    ----------
    kind:
        Discriminator literal ``"iga_fe_hybrid"``.
    fe_element_type:
        Solid element type used on the Native FEM subdomain.
    fe_n_cells:
        Subdivision counts for the FEM subdomain (defaults to ``(2, 2, 2)``).
    iga_degree:
        B-spline degree used on the IGA subdomain.
    iga_n_cells:
        Cell counts for the IGA subdomain (defaults to ``(2, 2, 2)``).
    split_axis:
        Axis index ``0``, ``1``, or ``2`` along which the box is partitioned.
    split_fraction:
        Normalized interface coordinate fraction in ``(0, 1)``.
    fe_side:
        ``"lo"`` places Native FEM on ``[lo, split]`` and IGA on ``[split, hi]``;
        ``"hi"`` places IGA on ``[lo, split]`` and Native FEM on ``[split, hi]``.
    coupling_kind:
        Interface tie coupling strategy (``"auto"``, ``"penalty_spring"``, ``"mortar_l2"``, ``"barycentric"``).
    penalty_beta:
        Dimensionless penalty factor used when tie resolves to ``"penalty_spring"``.
    """

    kind: Literal["iga_fe_hybrid"] = "iga_fe_hybrid"
    fe_element_type: FeElementType = "hex20"
    fe_n_cells: CellCounts | None = None
    iga_degree: Degree = 2
    iga_n_cells: CellCounts | None = None
    split_axis: Literal[0, 1, 2] = 2
    split_fraction: Annotated[float, Field(gt=0.0, lt=1.0)] = 0.5
    fe_side: Literal["lo", "hi"] = "lo"
    coupling_kind: HybridCouplingKind = "auto"
    penalty_beta: PositiveFinite = 1e3


DiscretizationOption = Annotated[
    Union[FullIgaOption, ThbAdaptiveIgaOption, NativeFeOption, IgaFeHybridOption],
    Field(discriminator="kind"),
]

_OPTION_ADAPTER = TypeAdapter(DiscretizationOption)


def parse_discretization_option(
    data: str | bytes | dict[str, Any] | SpecBase,
) -> DiscretizationOption:
    """Parse and validate a JSON string, bytes, mapping, or model into a ``DiscretizationOption``."""
    if isinstance(data, (str, bytes)):
        return _OPTION_ADAPTER.validate_json(data)
    return _OPTION_ADAPTER.validate_python(data)
