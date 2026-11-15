"""Complete 17-case Verification & Example catalog for WHT_IGA.

Maps Coreform IGA for Abaqus Verification Manual cases (``V01``–``V12``) and
Example Problems (``E01``–``E05``) to standardized probe definitions, parameter
sweeps, QoI extractors, phase gates, Coreform output CSV names, and reference
CSV paths under ``research-coreform-IGA/user-guide/``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from .schema import VerificationCaseSpec, load_case_spec

ProbeKind = Literal[
    "point",
    "face_grid",
    "line",
    "volume_quadrature",
    "history",
]
CaseCategory = Literal["verification", "example"]
PhaseGate = Literal["G1", "G2", "G3", "G4", "G5"]

_DEFAULT_USER_GUIDE_ROOT = (
    Path(__file__).resolve().parents[3] / "research-coreform-IGA" / "user-guide"
)


@dataclass(frozen=True)
class ProbeSpec:
    """Spatial, quadrature, or time-history QoI probe specification.

    Parameters
    ----------
    name:
        Unique identifier for the probe within its case.
    kind:
        Probe category: ``"point"``, ``"face_grid"``, ``"line"``,
        ``"volume_quadrature"``, or ``"history"``.
    fields:
        Tuple of physical field or QoI names extracted at this probe.
    coordinates:
        Point coordinates ``(x, y, z)`` or ``(x, y)`` for ``"point"`` probes.
    start:
        Start coordinate ``(x, y, z)`` for ``"line"`` probes.
    end:
        End coordinate ``(x, y, z)`` for ``"line"`` probes.
    n_points:
        Number of sampling points along a ``"line"`` probe.
    grid_shape:
        Grid dimensions (e.g. ``(21, 21)`` or ``(3, 3, 3)``) for
        ``"face_grid"`` probes.
    domain:
        Optional region, boundary face, time window, or entity set name.
    subcase:
        Optional subcase tag (e.g. ``"T2"``, ``"T3"``, ``"T4"`` for ``V12``).
    notes:
        Additional specification notes.
    """

    name: str
    kind: ProbeKind
    fields: tuple[str, ...]
    coordinates: tuple[float, ...] | None = None
    start: tuple[float, float, float] | None = None
    end: tuple[float, float, float] | None = None
    n_points: int | None = None
    grid_shape: tuple[int, ...] | None = None
    domain: str | None = None
    subcase: str | None = None
    notes: str = ""


@dataclass(frozen=True)
class SweepSpec:
    """Parametric sweep matrix for a verification or example case.

    Parameters
    ----------
    degrees:
        Polynomial spline degrees ``p`` to sweep.
    continuities:
        Inter-element continuity levels (e.g. ``"C^0"``, ``"C^(p-1)"``).
    mesh_sizes:
        Background or body-fitted element sizes ``h`` or ``(hx, hy, hz)`` tuples.
    workflows:
        Discretization or modeling workflows (e.g. ``"immersed"``,
        ``"body_fit"``, ``"native_fe"``, ``"boolean_merge"``, ``"tie"``,
        ``"local_refinement"``).
    load_cases:
        Load case variants for multi-load studies.
    time_steps:
        Time step sizes ``dt`` or end times for transient/nonlinear studies.
    padding:
        Background grid padding specification, if applicable.
    notes:
        Additional sweep configuration notes.
    """

    degrees: tuple[int, ...] = ()
    continuities: tuple[str, ...] = ()
    mesh_sizes: tuple[Any, ...] = ()
    workflows: tuple[str, ...] = ()
    load_cases: tuple[str, ...] = ()
    time_steps: tuple[float, ...] = ()
    padding: str | tuple[int, ...] | None = None
    notes: str = ""


@dataclass(frozen=True)
class CaseCatalogEntry:
    """Catalog metadata entry for a single verification or example case.

    Parameters
    ----------
    case_id:
        Normalized uppercase identifier (``"V01"``–``"V12"``, ``"E01"``–``"E05"``).
    slug:
        Coreform documentation directory slug.
    title:
        Human-readable case title.
    category:
        ``"verification"`` or ``"example"``.
    phase_gate:
        Development gate milestone (``"G1"``–``"G5"``).
    unit_system:
        Consistent unit system label.
    qoi_extractor:
        Name of the standard QoI extractor class in ``wht_iga.verification.qoi``.
    probes:
        Tuple of :class:`ProbeSpec` definitions for this case.
    sweeps:
        :class:`SweepSpec` parameter matrix for this case.
    reference_csv_relpaths:
        Tuple of reference CSV paths relative to
        ``research-coreform-IGA/user-guide/``.
    coreform_output_csvs:
        Tuple of Coreform-compatible output CSV/metric filenames produced by the
        verification pipeline.
    reference_doc_relpath:
        Source Markdown documentation path relative to
        ``research-coreform-IGA/user-guide/``.
    notes:
        Optional remarks or cross-gate requirements.
    """

    case_id: str
    slug: str
    title: str
    category: CaseCategory
    phase_gate: PhaseGate
    unit_system: str
    qoi_extractor: str
    probes: tuple[ProbeSpec, ...]
    sweeps: SweepSpec
    reference_csv_relpaths: tuple[str, ...] = field(default_factory=tuple)
    coreform_output_csvs: tuple[str, ...] = field(default_factory=tuple)
    reference_doc_relpath: str = ""
    notes: str = ""

    def load_spec(self, data_dir: str | Path | None = None) -> VerificationCaseSpec:
        """Load the corresponding JSON :class:`VerificationCaseSpec`.

        Parameters
        ----------
        data_dir:
            Optional directory override containing ``*.json`` case specs.

        Returns
        -------
        VerificationCaseSpec
            Parsed specification loaded via :func:`load_case_spec`.
        """
        return load_case_spec(self.case_id, data_dir=data_dir)


_V03_LOAD_CASES: tuple[str, ...] = (
    "gravity",
    "body_force",
    "pressure_uniform",
    "pressure_total_force",
    "traction",
    "moment",
)

_V03_BASE_REL = (
    "generated-downloads/verification-manual/problems/solid_mechanics/"
    "linear_elastic_stress/beam-bending"
)

_V04_LOAD_CASES: tuple[str, ...] = (
    "gravity",
    "body_force",
    "pressure",
    "traction",
    "moment",
)

_V04_BASE_REL = (
    "generated-downloads/verification-manual/problems/solid_mechanics/"
    "linear_elastic_stress/bimetallic-beam-bending"
)

_V05_SECTIONS: tuple[str, ...] = (
    "bottom-center-third",
    "bottom-left-half",
    "centered-third",
    "nested-l",
    "tetris",
    "zig-zag",
)

_V05_BASE_REL = (
    "generated-downloads/verification-manual/problems/solid_mechanics/"
    "linear_elastic_stress/multiplanar-material-interfaces"
)

_V06_BASE_REL = (
    "generated-downloads/verification-manual/problems/solid_mechanics/"
    "linear_elastic_stress/pressurized-sphere"
)

_V02_BASE_REL = (
    "generated-downloads/verification-manual/problems/solid_mechanics/"
    "linear_elastic_stress/roark-17p1-7a"
)

_E03_GEN_REL = "generated-downloads/example-problems/glass-to-metal-seal"
_E03_DATA_REL = "example-problems/glass-to-metal-seal/data"


def _build_v03_ref_csvs() -> tuple[str, ...]:
    slugs = (
        "gravity",
        "body-force",
        "pressure-uniform",
        "pressure-total-force",
        "traction",
        "moment",
    )
    paths = [
        f"{_V03_BASE_REL}/beam-bending-accepted.csv",
        f"{_V03_BASE_REL}/beam-bending-traction-variants.csv",
    ]
    for slug in slugs:
        paths.append(f"{_V03_BASE_REL}/beam-bending-{slug}-tip-deflection.csv")
        paths.append(f"{_V03_BASE_REL}/beam-bending-{slug}-trim-performance.csv")
        paths.append(f"{_V03_BASE_REL}/beam-bending-{slug}-solve-performance.csv")
    return tuple(paths)


def _build_v04_ref_csvs() -> tuple[str, ...]:
    slugs = ("gravity", "body-force", "pressure", "traction", "moment")
    paths = [
        f"{_V04_BASE_REL}/bimetallic-beam-accepted.csv",
        f"{_V04_BASE_REL}/bimetallic-beam-tie-comparison.csv",
    ]
    for slug in slugs:
        paths.append(f"{_V04_BASE_REL}/bimetallic-beam-{slug}-tip-deflection.csv")
        paths.append(f"{_V04_BASE_REL}/bimetallic-beam-{slug}-trim-performance.csv")
        paths.append(f"{_V04_BASE_REL}/bimetallic-beam-{slug}-solve-performance.csv")
    return tuple(paths)


def _build_v05_ref_csvs() -> tuple[str, ...]:
    paths = [
        f"{_V05_BASE_REL}/multiplanar-material-interface-monolithic-trim-performance.csv",
        f"{_V05_BASE_REL}/multiplanar-material-interface-monolithic-solve-performance.csv",
    ]
    for sec in _V05_SECTIONS:
        paths.append(f"{_V05_BASE_REL}/multiplanar-material-interface-{sec}-responses.csv")
        paths.append(f"{_V05_BASE_REL}/multiplanar-material-interface-{sec}-trim-performance.csv")
        paths.append(f"{_V05_BASE_REL}/multiplanar-material-interface-{sec}-solve-performance.csv")
    return tuple(paths)


CASE_CATALOG: dict[str, CaseCatalogEntry] = {
    "V01": CaseCatalogEntry(
        case_id="V01",
        slug="prescribed-displacement-cube",
        title="Prescribed displacement of a unit cube",
        category="verification",
        phase_gate="G1",
        unit_system="dimensionless (E=1, nu=0, L=1)",
        qoi_extractor="FaceAndGridProbeExtractor",
        probes=(
            ProbeSpec(
                name="six_faces_grid",
                kind="face_grid",
                fields=("U", "S_nn"),
                grid_shape=(21, 21),
                domain="x in {0,1}, y in {0,1}, z in {0,1} on [0,1]^3",
                notes="441 uniform points per face; checks probe coord err = 0",
            ),
            ProbeSpec(
                name="surface_layer_nodes",
                kind="face_grid",
                fields=("U",),
                domain="Triangle3D3 surface layer (Coreform SFM3D3 equivalent)",
            ),
            ProbeSpec(
                name="interior_grid",
                kind="face_grid",
                fields=("S11", "S22", "S33", "S12", "S13", "S23"),
                grid_shape=(3, 3, 3),
                domain="{0.1, 0.5, 0.9}^3",
                notes="27 interior points verifying all 6 stress components",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^0", "C^1", "C^2"),
            mesh_sizes=(0.6,),
            workflows=("immersed",),
            notes="6 (p, C) pairs: (1,0), (2,0), (2,1), (3,0), (3,1), (3,2); center (0.5,0.5,0.5)",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("prescribed-displacement-cube-summary.csv",),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "prescribed-displacement-cube/prescribed-displacement-cube.md"
        ),
        notes="Must also be re-passed with extension enabled at G2.",
    ),
    "V02": CaseCatalogEntry(
        case_id="V02",
        slug="roark-17p1-7a",
        title="Stress concentration in a thin plate with hole (Roark 17.1-7a)",
        category="verification",
        phase_gate="G1",
        unit_system="dimensionless (L=200, D=40, t=1, r=10, E=1, nu=0)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="point_A_peak",
                kind="point",
                fields=("S_max_principal",),
                coordinates=(100.0, 10.0, 0.5),
                domain="hole edge (L/2, r, t/2)",
            ),
            ProbeSpec(
                name="transverse_ligament",
                kind="line",
                fields=("S_max_principal", "Ux", "Uy"),
                start=(100.0, 10.0, 0.5),
                end=(100.0, 20.0, 0.5),
                n_points=100,
                domain="x=L/2, y in [10, 20], z=t/2",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^(p-1)",),
            mesh_sizes=(6.0, 3.0, 1.5, 0.75, 0.375),
            workflows=("immersed", "local_refinement"),
            notes="Global sweep k=0..4 (6/2^k); Phase 3 local THB L=1,2,3 from h=3.0, p=2",
        ),
        reference_csv_relpaths=(
            f"{_V02_BASE_REL}/roark-17p1-7a-accepted-peak.csv",
            f"{_V02_BASE_REL}/roark-17p1-7a-qoi.csv",
            f"{_V02_BASE_REL}/roark-17p1-7a-reference-performance.csv",
            f"{_V02_BASE_REL}/roark-17p1-7a-refinement.csv",
        ),
        coreform_output_csvs=(
            "roark-17p1-7a-qoi.csv",
            "roark-17p1-7a-refinement.csv",
            "roark-17p1-7a-trim-performance.csv",
            "roark-17p1-7a-solve-performance.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "roark-17p1-7a/roark-17p1-7a.md"
        ),
        notes="Global sweep at G1, re-passed at G2 with extension, local THB at G3.",
    ),
    "V03": CaseCatalogEntry(
        case_id="V03",
        slug="beam-bending",
        title="Cantilever beam bending",
        category="verification",
        phase_gate="G1",
        unit_system="SI (m, kg, s, N, Pa)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="tip_deflection",
                kind="point",
                fields=("Uy",),
                coordinates=(0.0, 0.0, 1.0),
                domain="free end (0, 0, L)",
            ),
            ProbeSpec(
                name="beam_axis",
                kind="line",
                fields=("Uy",),
                start=(0.0, 0.0, 0.0),
                end=(0.0, 0.0, 1.0),
                n_points=100,
                domain="x=0, y=0, z in [0, L]",
            ),
            ProbeSpec(
                name="top_fiber",
                kind="line",
                fields=("Szz",),
                start=(0.0, 0.05, 0.0),
                end=(0.0, 0.05, 1.0),
                n_points=100,
                domain="x=0, y=h/2, z in [0, L]",
            ),
            ProbeSpec(
                name="midspan_section",
                kind="line",
                fields=("Szz", "Szy"),
                start=(0.0, -0.05, 0.5),
                end=(0.0, 0.05, 0.5),
                n_points=100,
                domain="x=0, y in [-h/2, h/2], z=L/2",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^(p-1)",),
            mesh_sizes=(
                (0.1, 0.1, 0.25),
                (0.05, 0.05, 0.125),
                (0.025, 0.025, 0.0625),
                (0.0125, 0.0125, 0.03125),
                (0.00625, 0.00625, 0.015625),
            ),
            workflows=("immersed", "native_fe"),
            load_cases=_V03_LOAD_CASES,
            padding="p elements per side",
            notes="moment load gated at G3; follower/resultant traction variants at G5",
        ),
        reference_csv_relpaths=_build_v03_ref_csvs(),
        coreform_output_csvs=(
            "beam-bending-<load>-tip-deflection.csv",
            "beam-bending-<load>-trim-performance.csv",
            "beam-bending-<load>-solve-performance.csv",
            "beam-bending-traction-variants.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "beam-bending/beam-bending.md"
        ),
    ),
    "V04": CaseCatalogEntry(
        case_id="V04",
        slug="bimetallic-beam-bending",
        title="Cantilever bimetallic beam bending",
        category="verification",
        phase_gate="G3",
        unit_system="MMTS (mm, tonne, s, N, MPa)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="tip_deflection",
                kind="point",
                fields=("Uy",),
                coordinates=(0.0, 0.0, 100.0),
                domain="free end (0, 0, L)",
            ),
            ProbeSpec(
                name="interface_deflection",
                kind="line",
                fields=("Uy",),
                start=(0.0, 0.0, 0.0),
                end=(0.0, 0.0, 100.0),
                n_points=100,
                domain="Cu/W interface y=0, z in [0, L]",
            ),
            ProbeSpec(
                name="midspan_section",
                kind="line",
                fields=("Szz", "Szy"),
                start=(0.0, -5.0, 50.0),
                end=(0.0, 5.0, 50.0),
                n_points=100,
                domain="x=0, y in [-5, 5], z=50 (sigma_zz jump at y=0)",
            ),
            ProbeSpec(
                name="fixed_root_section",
                kind="line",
                fields=("Szz", "Szy"),
                start=(0.0, -5.0, 0.0),
                end=(0.0, 5.0, 0.0),
                n_points=100,
                domain="x=0, y in [-5, 5], z=0",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^(p-1) within each layer",),
            mesh_sizes=(
                (5.0, 5.0, 10.0),
                (2.5, 2.5, 5.0),
                (1.25, 1.25, 2.5),
                (0.625, 0.625, 1.25),
            ),
            workflows=("boolean_merge", "tie", "iga_fe_hybrid"),
            load_cases=_V04_LOAD_CASES,
            padding="p elements per side",
        ),
        reference_csv_relpaths=_build_v04_ref_csvs(),
        coreform_output_csvs=(
            "bimetallic-beam-tie-comparison.csv",
            "bimetallic-beam-<load>-tip-deflection.csv",
            "bimetallic-beam-<load>-trim-performance.csv",
            "bimetallic-beam-<load>-solve-performance.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "bimetallic-beam-bending/bimetallic-beam-bending.md"
        ),
    ),
    "V05": CaseCatalogEntry(
        case_id="V05",
        slug="multiplanar-material-interfaces",
        title="Multi-planar material interfaces",
        category="verification",
        phase_gate="G3",
        unit_system="unknown (magnitudes consistent with SI: E=100e9, P=1)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="tip_deflection",
                kind="point",
                fields=("Uy",),
                coordinates=(0.0, 0.0, 100.0),
                domain="free end (0, 0, L)",
            ),
            ProbeSpec(
                name="axial_displacement_line",
                kind="line",
                fields=("Uy",),
                start=(0.0, 0.0, 0.0),
                end=(0.0, 0.0, 100.0),
                n_points=100,
                domain="beam axis z in [0, L]",
            ),
            ProbeSpec(
                name="vertical_section_stress",
                kind="line",
                fields=("Szz", "Szy"),
                start=(0.0, -5.0, 50.0),
                end=(0.0, 5.0, 50.0),
                n_points=100,
                domain="midspan vertical line x=0, y in [-5, 5], z=50",
            ),
            ProbeSpec(
                name="horizontal_section_stress",
                kind="line",
                fields=("Szz", "Szx"),
                start=(-5.0, 0.0, 50.0),
                end=(5.0, 0.0, 50.0),
                n_points=100,
                domain="midspan horizontal line x in [-5, 5], y=0, z=50",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^(p-1)",),
            mesh_sizes=(
                (5.0, 5.0, 10.0),
                (2.5, 2.5, 5.0),
                (1.25, 1.25, 2.5),
                (0.625, 0.625, 1.25),
            ),
            workflows=("monolithic", *_V05_SECTIONS),
            notes="7 cross-sections x 4 grids x 3 degrees = 84 simulations",
        ),
        reference_csv_relpaths=_build_v05_ref_csvs(),
        coreform_output_csvs=(
            "multiplanar-material-interface-<case>-responses.csv",
            "multiplanar-material-interface-<case>-trim-performance.csv",
            "multiplanar-material-interface-<case>-solve-performance.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "multiplanar-material-interfaces/multiplanar-material-interfaces.md"
        ),
    ),
    "V06": CaseCatalogEntry(
        case_id="V06",
        slug="pressurized-sphere",
        title="Pressurized sphere (thick-walled, Lame)",
        category="verification",
        phase_gate="G2",
        unit_system="dimensionless (a=0.9, b=1.0, E=1, nu=0, p_i=1e-6)",
        qoi_extractor="VolumeQuadratureL2Extractor",
        probes=(
            ProbeSpec(
                name="volume_quadrature_l2",
                kind="volume_quadrature",
                fields=("U", "Ur", "S_rr", "S_theta_theta"),
                domain="1/8 hollow sphere 0.9 <= r <= 1.0, x,y,z >= 0",
                notes="Spherical projection at quadrature points vs Lame exact solution",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1, 2, 3),
            continuities=("C^(p-1)",),
            mesh_sizes=(0.1, 0.05, 0.025),
            workflows=("immersed",),
            notes="h = 0.1/2^k for k = 0, 1, 2",
        ),
        reference_csv_relpaths=(
            f"{_V06_BASE_REL}/pressurized-sphere-qoi.csv",
            f"{_V06_BASE_REL}/pressurized-sphere-solve-performance.csv",
            f"{_V06_BASE_REL}/pressurized-sphere-trim-performance.csv",
        ),
        coreform_output_csvs=(
            "pressurized-sphere-qoi.csv",
            "pressurized-sphere-trim-performance.csv",
            "pressurized-sphere-solve-performance.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "pressurized-sphere/pressurized-sphere.md"
        ),
    ),
    "V07": CaseCatalogEntry(
        case_id="V07",
        slug="infinite-plate-with-hole",
        title="Infinite plate with hole (Kirsch, nearly incompressible)",
        category="verification",
        phase_gate="G5",
        unit_system="dimensionless (a=1, E=1e5, nu=0.49999, T=10)",
        qoi_extractor="VolumeQuadratureL2Extractor",
        probes=(
            ProbeSpec(
                name="volume_quadrature_kirsch",
                kind="volume_quadrature",
                fields=("U", "Sxx", "volume"),
                domain="1/4 plate [0,4]^2 x [0,0.25] minus r=1 hole",
                notes="Global L2 errors (e_u, e_sigma_xx), domain volume error, log-log slopes",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=(0.5, 0.25, 0.125),
            workflows=("immersed_mixed_vms",),
            notes="Mixed U(2,C^1)/P(2,C^1) continuous pressure + VMS stabilization",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("global_l2_metrics.json",),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/linear_elastic_stress/"
            "infinite-plate-with-hole/infinite-plate-with-hole.md"
        ),
    ),
    "V08": CaseCatalogEntry(
        case_id="V08",
        slug="contact-indenter",
        title="Elastic spherical indenter sliding on block (Hertz-Coulomb contact)",
        category="verification",
        phase_gate="G5",
        unit_system="IPS (in, lbf, psi, s)",
        qoi_extractor="ContactHistoryExtractor",
        probes=(
            ProbeSpec(
                name="contact_reaction_history",
                kind="history",
                fields=("RF1", "RF2", "RF3", "p0_max_contact_pressure"),
                domain="t in [0.0, 1.0]; indentation at t=0.1, sliding window [0.15, 0.70]",
                notes="Compares F_N, F_T, p0 against Hertz-Coulomb analytical values",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=(0.5, 0.25, 0.125),
            workflows=("immersed", "native_fe"),
            time_steps=(0.1, 1.0),
            notes="NLGEOM=YES, mu=0.3; Step 1 indentation (0..0.1), Step 2 sliding (0.1..1.0)",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("contact-indenter-history.csv",),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/contact/"
            "contact-indenter/contact-indenter.md"
        ),
    ),
    "V09": CaseCatalogEntry(
        case_id="V09",
        slug="cantilever-shock-vibration",
        title="Cantilever beam shock and vibration",
        category="verification",
        phase_gate="G4",
        unit_system="SI (m, kg, s, N, Pa, Hz)",
        qoi_extractor="ModalAndSpectrumExtractor",
        probes=(
            ProbeSpec(
                name="modal_and_tip_response",
                kind="history",
                fields=("eigenfrequency", "gamma_y", "Uy_psd_rms", "Uy_srss_peak"),
                coordinates=(0.0, 0.0, 1.0),
                domain="20 modes; random 10-2000 Hz PSD 0.01 g^2/Hz; SRSS xi=0.02",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1,),
            continuities=("C^0",),
            mesh_sizes=(0.025, 0.05),
            workflows=("immersed", "body_fit", "native_fe"),
            notes="Body-fit C3D8 50 mm vs Immersed IGA 25 mm p=1",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=(
            "eigenfrequencies.csv",
            "implementation.csv",
            "response.csv",
        ),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/shock_and_vibration/"
            "cantilever-shock-vibration/cantilever-shock-vibration.md"
        ),
    ),
    "V10": CaseCatalogEntry(
        case_id="V10",
        slug="vibrating-beam",
        title="Implicit-dynamic vibration of 3D beam (Bossak)",
        category="verification",
        phase_gate="G4",
        unit_system="dimensionless (L=10, b=2, h=2, E=1e6, rho=1e-6, nu=0)",
        qoi_extractor="DynamicHistoryExtractor",
        probes=(
            ProbeSpec(
                name="free_end_axial_history",
                kind="history",
                fields=("U_axial", "V_axial", "A_axial"),
                domain="free end average axial history over t in [0, 2e-3] s",
                notes="Extracts f1, norm disp RMS, phase RMS, amp/energy decay, step ratio, samples/period",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=(1.0,),
            workflows=("native_fe", "body_fit", "immersed"),
            load_cases=("TRANSIENT_FIDELITY", "MODERATE_DISSIPATION"),
            time_steps=(1.0e-6, 2.0e-6),
            notes="3 discretizations x 2 Bossak alpha_B profiles = 6 cases, 26 criteria",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("vibrating-beam-metrics.csv",),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/implicit_dynamics/"
            "vibrating-beam/vibrating-beam.md"
        ),
    ),
    "V11": CaseCatalogEntry(
        case_id="V11",
        slug="nafems-p18-t1",
        title="NAFEMS P18 T1 — membrane with circular hot spot (thermal strain)",
        category="verification",
        phase_gate="G4",
        unit_system="MMTS (mm, N, MPa)",
        qoi_extractor="FaceAndGridProbeExtractor",
        probes=(
            ProbeSpec(
                name="point_D_stress",
                kind="point",
                fields=("Syy",),
                coordinates=(1.0, 0.0, 0.125),
                domain="outer interface point D(1, 0) mm",
                notes="Reference sigma_yy(D) = 50.0 MPa (< 4% rel error)",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1,),
            continuities=("C^0",),
            mesh_sizes=((0.2272727, 0.2272727, 0.25),),
            workflows=("immersed",),
            notes="Uniform 44 x 44 x 1 background grid, 2-cell multicell shared field",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("nafems-p18-t1-qoi.csv",),
        reference_doc_relpath=(
            "verification-manual/problems/solid_mechanics/thermal_strain/"
            "nafems-p18-t1/nafems-p18-t1.md"
        ),
    ),
    "V12": CaseCatalogEntry(
        case_id="V12",
        slug="nafems-p18-t2-t3-t4",
        title="NAFEMS P18 T2/T3/T4 — heat transfer (radiation, transient, convection)",
        category="verification",
        phase_gate="G4",
        unit_system="SI (m, kg, s, W, K, degC)",
        qoi_extractor="FaceAndGridProbeExtractor",
        probes=(
            ProbeSpec(
                name="T2_end_B_temperature",
                kind="point",
                fields=("NT",),
                coordinates=(0.1, 0.0, 0.0),
                domain="end B (x = 0.1 m), steady radiation",
                subcase="T2",
                notes="Reference T_B = 927.0 K (< 1% rel error)",
            ),
            ProbeSpec(
                name="T3_point_C_history",
                kind="history",
                fields=("NT",),
                coordinates=(0.08, 0.0, 0.0),
                domain="x = 0.08 m at t = 32 s, transient conduction",
                subcase="T3",
                notes="Reference T_C(32 s) = 36.60 C (< 2% rel error)",
            ),
            ProbeSpec(
                name="T4_point_E_temperature",
                kind="point",
                fields=("NT",),
                coordinates=(0.6, 0.2, 0.025),
                domain="(0.6, 0.2) m on 0.6 x 1.0 x 0.05 m slab, steady convection",
                subcase="T4",
                notes="Reference T_E = 18.3 C (< 3% rel error)",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(1,),
            continuities=("C^0",),
            mesh_sizes=((0.01, 0.01, 0.01), (0.05, 0.05, 0.05)),
            workflows=("body_fit",),
            load_cases=("T2_steady_radiation", "T3_transient_conduction", "T4_steady_convection"),
            time_steps=(0.1, 0.5),
            notes="T2/T3: 10x1x1 on L=0.1m; T4: 12x20x1 on 0.6x1.0x0.05m",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("nafems-p18-heat-transfer-qoi.csv",),
        reference_doc_relpath=(
            "verification-manual/problems/heat_transfer/steady_state/"
            "nafems-p18-t2/nafems-p18-t2.md"
        ),
    ),
    "E01": CaseCatalogEntry(
        case_id="E01",
        slug="thin-plate-with-hole",
        title="Thin plate with hole (IPS) — stress concentration",
        category="example",
        phase_gate="G3",
        unit_system="IPS (in, slinch, s, lbf, psi)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="Probe-Stress-Profile",
                kind="line",
                fields=("S_max_principal", "S_mises", "U"),
                start=(0.0, 1.0, 0.05),
                end=(0.0, 2.5, 0.05),
                n_points=100,
                domain="hole edge y=1.0 to free edge y=2.5 at x=0, z=0.05",
                notes="Roark Kt=2.23568 -> sigma_max=3.7261e-6 psi at (0, 1.0, 0.05)",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=((0.3, 0.3, 0.1),),
            workflows=("immersed",),
            notes="10 x 5 x 0.1 in plate, r=1 in hole, X_max P=-1e-6 psi",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("thin-plate-with-hole-probe.csv",),
        reference_doc_relpath="example-problems/thin-plate-with-hole/thin-plate-with-hole.md",
    ),
    "E02": CaseCatalogEntry(
        case_id="E02",
        slug="plate-with-hole-adaptivity",
        title="Plate with hole — local refinement toward a feature (THB)",
        category="example",
        phase_gate="G3",
        unit_system="IPS (in, slinch, s, lbf, psi)",
        qoi_extractor="LineProfileExtractor",
        probes=(
            ProbeSpec(
                name="Probe-Stress-Profile",
                kind="line",
                fields=("S_max_principal", "S_mises", "U"),
                start=(0.0, 1.0, 0.05),
                end=(0.0, 2.5, 0.05),
                n_points=100,
                domain="hole edge y=1.0 to free edge y=2.5 at x=0, z=0.05",
                notes="Compares THB level=2 (0.3 -> 0.15 -> 0.075 in) vs E01 uniform",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=((0.3, 0.3, 0.1), (0.15, 0.15, 0.1), (0.075, 0.075, 0.1)),
            workflows=("local_refinement",),
            notes="CadEntitiesRegion(sets=['hole'], offset=0.1), level=2, two_to_one balancing",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("plate-with-hole-adaptivity-probe.csv",),
        reference_doc_relpath=(
            "example-problems/plate-with-hole-adaptivity/plate-with-hole-adaptivity.md"
        ),
    ),
    "E03": CaseCatalogEntry(
        case_id="E03",
        slug="glass-to-metal-seal",
        title="Glass-to-metal seal thermal cycling (13-cell multicell, T-dependent J2)",
        category="example",
        phase_gate="G5",
        unit_system="MMTS (mm, tonne, s, N, MPa, degC)",
        qoi_extractor="NonlinearThermoExtractor",
        probes=(
            ProbeSpec(
                name="glass_max_tensile_principal",
                kind="volume_quadrature",
                fields=("S_max_principal_positive", "U"),
                domain="S-8061 glass cell; sigma_1+ = max(sigma_1, 0)",
            ),
            ProbeSpec(
                name="metal_plastic_strain",
                kind="volume_quadrature",
                fields=("PEEQ", "PE", "LE"),
                domain="304L VAR header + 11 Alloy 52 pins",
            ),
            ProbeSpec(
                name="section_y0_viewcut",
                kind="face_grid",
                fields=("S_max_principal_positive", "PEEQ", "U"),
                domain="longitudinal cross-section y = 0",
            ),
            ProbeSpec(
                name="interface_paths",
                kind="line",
                fields=("S_max_principal_positive", "PEEQ"),
                domain="glass-pin and glass-header interfaces",
            ),
            ProbeSpec(
                name="solver_convergence_history",
                kind="history",
                fields=("accepted_increments", "iterations", "wall_time_s", "dofs"),
                domain="3 thermal steps: 450 -> 20 -> -55 -> 200 degC",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=((0.5, 0.5, 0.5),),
            workflows=("boolean_merge", "immersed"),
            load_cases=("baseline", "P0", "P1", "P2", "P3", "P4"),
            notes="13-solid Boolean-merged multicell; NLGEOM=ON; 15 increments/step target",
        ),
        reference_csv_relpaths=(
            f"{_E03_DATA_REL}/alloy52_2007_powerlaw_digitized.csv",
            f"{_E03_DATA_REL}/s8061_modulus_digitized.csv",
            f"{_E03_DATA_REL}/ss304l_var_mlep_digitized.csv",
            f"{_E03_DATA_REL}/thermal_strain_recommended.csv",
            f"{_E03_GEN_REL}/glass-to-metal-seal-convergence-history.csv",
            f"{_E03_GEN_REL}/glass-to-metal-seal-performance.csv",
            f"{_E03_GEN_REL}/glass-to-metal-seal-problem-size.csv",
        ),
        coreform_output_csvs=(
            "glass-to-metal-seal-convergence-history.csv",
            "glass-to-metal-seal-performance.csv",
            "glass-to-metal-seal-problem-size.csv",
        ),
        reference_doc_relpath="example-problems/glass-to-metal-seal/glass-to-metal-seal.md",
    ),
    "E04": CaseCatalogEntry(
        case_id="E04",
        slug="nearly-incompressible-elasticity",
        title="Cook's membrane (nearly incompressible, locking-free, NLGEOM)",
        category="example",
        phase_gate="G5",
        unit_system="MMTS (mm, N, MPa)",
        qoi_extractor="FaceAndGridProbeExtractor",
        probes=(
            ProbeSpec(
                name="top_right_corner",
                kind="point",
                fields=("Uy",),
                coordinates=(48.0, 60.0, 0.95758),
                domain="top-right corner (L, H_R + Delta_y + H_L?, max-y tip)",
                notes="Target final u_y = 6.98946 mm within 1% rel error",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^0",),
            mesh_sizes=(),
            workflows=("body_fit",),
            notes="U(p=2, C^0) / P0(C^-1) static condensation, Neo-Hookean nu=0.49999, NLGEOM=ON",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("cooks-membrane-tip-uy.csv",),
        reference_doc_relpath=(
            "example-problems/nearly-incompressible-elasticity/"
            "nearly-incompressible-elasticity.md"
        ),
    ),
    "E05": CaseCatalogEntry(
        case_id="E05",
        slug="door-lock-assembly",
        title="Door-lock assembly (3-part frictional contact, elastic-plastic)",
        category="example",
        phase_gate="G5",
        unit_system="MMTS (mm, N, MPa)",
        qoi_extractor="ContactHistoryExtractor",
        probes=(
            ProbeSpec(
                name="hold_link_rf3_history",
                kind="history",
                fields=("RF3", "S", "U", "PEEQ"),
                domain="hold_link_surfaces (summed link reaction RF3)",
                notes="Target final RF3 = -46.946 kN within 5% rel error",
            ),
        ),
        sweeps=SweepSpec(
            degrees=(2,),
            continuities=("C^1",),
            mesh_sizes=((2.0, 2.0, 2.0),),
            workflows=("immersed",),
            padding=(2, 2, 2),
            notes="3 parts (Hinge, Link, Pin), 3 Triangle3D3 contact pairs mu=0.3, NLGEOM=ON",
        ),
        reference_csv_relpaths=(),
        coreform_output_csvs=("door-lock-assembly-rf3-history.csv",),
        reference_doc_relpath="example-problems/door-lock-assembly/door-lock-assembly.md",
    ),
}


def get_catalog_entry(case_id: str) -> CaseCatalogEntry:
    """Return the :class:`CaseCatalogEntry` for ``case_id``.

    Parameters
    ----------
    case_id:
        Case identifier (case-insensitive, e.g. ``"v01"``, ``"V01"``, ``"e03"``).

    Returns
    -------
    CaseCatalogEntry
        Catalog metadata for the requested case.

    Raises
    ------
    KeyError
        If ``case_id`` is not present in :data:`CASE_CATALOG`.
    """
    normalized = case_id.strip().upper()
    if normalized not in CASE_CATALOG:
        available = ", ".join(CASE_CATALOG.keys())
        raise KeyError(
            f"Unknown verification/example case '{case_id}'. Available: {available}"
        )
    return CASE_CATALOG[normalized]


def list_cases_by_phase(phase_gate: str) -> list[CaseCatalogEntry]:
    """List catalog entries belonging to a specific phase gate.

    Parameters
    ----------
    phase_gate:
        Gate identifier (e.g. ``"G1"``..``"G5"``, case-insensitive).

    Returns
    -------
    list[CaseCatalogEntry]
        Ordered list of catalog entries assigned to ``phase_gate``.

    Raises
    ------
    ValueError
        If ``phase_gate`` is not one of ``G1``..``G5``.
    """
    normalized = phase_gate.strip().upper()
    valid_gates = ("G1", "G2", "G3", "G4", "G5")
    if normalized not in valid_gates:
        raise ValueError(
            f"Unknown phase gate '{phase_gate}'. Expected one of {valid_gates}."
        )
    return [entry for entry in CASE_CATALOG.values() if entry.phase_gate == normalized]


def resolve_reference_csvs(
    case_id: str,
    user_guide_root: str | Path | None = None,
) -> list[Path]:
    """Resolve Coreform reference CSV paths for ``case_id``.

    Parameters
    ----------
    case_id:
        Case identifier (e.g. ``"V02"`` or ``"e03"``).
    user_guide_root:
        Optional root directory overriding
        ``research-coreform-IGA/user-guide``.

    Returns
    -------
    list[Path]
        Absolute or root-joined :class:`~pathlib.Path` objects for each
        declared entry in ``reference_csv_relpaths``.
    """
    entry = get_catalog_entry(case_id)
    root = (
        Path(user_guide_root)
        if user_guide_root is not None
        else _DEFAULT_USER_GUIDE_ROOT
    )
    return [root / rel for rel in entry.reference_csv_relpaths]


__all__ = [
    "CASE_CATALOG",
    "CaseCatalogEntry",
    "CaseCategory",
    "PhaseGate",
    "ProbeKind",
    "ProbeSpec",
    "SweepSpec",
    "get_catalog_entry",
    "list_cases_by_phase",
    "resolve_reference_csvs",
]
