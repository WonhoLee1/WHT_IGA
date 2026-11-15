"""External-process Gmsh mesh generator (v2.1 §5.10.5, D13).

GPL isolation rule: this module MUST NOT ``import gmsh``. Instead, it generates
a ``.geo`` control script and executes the external ``gmsh`` CLI binary via
:func:`subprocess.run`.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Mapping

from wht_iga.meshio_ext.model import make_issue

ElementType = Literal["tet4", "tet10", "hex8"]
MshFormat = Literal["msh2", "msh4"]

SUPPORTED_CAD_EXTENSIONS: frozenset[str] = frozenset(
    {".step", ".stp", ".iges", ".igs", ".brep"}
)
_MSH_VERSION_MAP: dict[str, str] = {
    "msh2": "2.2",
    "msh4": "4.1",
}


@dataclass
class GmshOptions:
    """Meshing options passed to the external Gmsh process.

    Parameters
    ----------
    element_type:
        Target 3D element type: ``"tet4"``, ``"tet10"``, or ``"hex8"``.
    mesh_size_min:
        Minimum characteristic element length (``Mesh.CharacteristicLengthMin``).
        If ``None``, Gmsh's default is left unchanged.
    mesh_size_max:
        Maximum characteristic element length (``Mesh.CharacteristicLengthMax``).
        If ``None``, Gmsh's default is left unchanged.
    order:
        Element polynomial order (``1`` or ``2``). When ``None``, defaults to
        ``2`` for ``"tet10"`` and ``1`` for ``"tet4"`` / ``"hex8"``.
    algorithm_3d:
        Gmsh 3D meshing algorithm ID (``Mesh.Algorithm3D``; e.g. ``1`` = Delaunay,
        ``4`` = Frontal, ``10`` = HXT).
    optimize:
        Whether to enable mesh quality optimization (``Mesh.Optimize``).
    threads:
        Number of OpenMP threads (``General.NumThreads``).
    format:
        Output ``.msh`` format version: ``"msh2"`` (2.2) or ``"msh4"`` (4.1).
    """

    element_type: ElementType = "tet4"
    mesh_size_min: float | None = None
    mesh_size_max: float | None = None
    order: int | None = None
    algorithm_3d: int = 1
    optimize: bool = True
    threads: int = 4
    format: MshFormat = "msh2"

    def __post_init__(self) -> None:
        if self.element_type not in ("tet4", "tet10", "hex8"):
            raise ValueError(
                f"unsupported element_type {self.element_type!r}; "
                "expected 'tet4', 'tet10', or 'hex8'"
            )
        if self.order is None:
            self.order = 2 if self.element_type == "tet10" else 1
        elif self.order not in (1, 2):
            raise ValueError(f"order must be 1 or 2, got {self.order!r}")

        expected_order = 2 if self.element_type == "tet10" else 1
        if self.order != expected_order:
            raise ValueError(
                f"element_type={self.element_type!r} requires order={expected_order}, "
                f"got order={self.order}"
            )

        if self.mesh_size_min is not None and self.mesh_size_min <= 0.0:
            raise ValueError(
                f"mesh_size_min must be positive, got {self.mesh_size_min!r}"
            )
        if self.mesh_size_max is not None and self.mesh_size_max <= 0.0:
            raise ValueError(
                f"mesh_size_max must be positive, got {self.mesh_size_max!r}"
            )
        if (
            self.mesh_size_min is not None
            and self.mesh_size_max is not None
            and self.mesh_size_min > self.mesh_size_max
        ):
            raise ValueError(
                f"mesh_size_min ({self.mesh_size_min}) must not exceed "
                f"mesh_size_max ({self.mesh_size_max})"
            )
        if self.algorithm_3d <= 0:
            raise ValueError(
                f"algorithm_3d must be positive, got {self.algorithm_3d!r}"
            )
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads!r}")
        if self.format not in _MSH_VERSION_MAP:
            raise ValueError(
                f"unsupported format {self.format!r}; expected 'msh2' or 'msh4'"
            )


@dataclass
class GmshResult:
    """Outcome of an external Gmsh mesh-generation run.

    Parameters
    ----------
    msh_path:
        Target ``.msh`` output file path.
    command:
        Exact subprocess argument list invoked (empty if not executed).
    returncode:
        Process exit code (``0`` on success; non-zero on failure or missing exe).
    issues:
        List of issue dicts (``{"severity", "code", "message", "entity"}``).
    stdout:
        Captured standard output from Gmsh.
    stderr:
        Captured standard error from Gmsh.
    geo_path:
        Path of the generated ``.geo`` script (if persisted).
    """

    msh_path: str
    command: list[str]
    returncode: int
    issues: list[dict[str, str]] = field(default_factory=list)
    stdout: str = ""
    stderr: str = ""
    geo_path: str = ""

    @property
    def ok(self) -> bool:
        """Return ``True`` when Gmsh succeeded without error-severity issues."""
        return self.returncode == 0 and not any(
            issue.get("severity") == "error" for issue in self.issues
        )


def find_gmsh_executable(gmsh_exe: str | Path | None = None) -> str | None:
    """Locate the external ``gmsh`` executable without importing ``gmsh``.

    Parameters
    ----------
    gmsh_exe:
        Optional explicit executable path or command name. When ``None``,
        searches ``PATH`` for ``"gmsh"``.

    Returns
    -------
    str | None
        Resolved executable path if found, otherwise ``None``.
    """
    if gmsh_exe is not None:
        candidate = Path(gmsh_exe)
        if candidate.is_file():
            return str(candidate)
        return shutil.which(str(gmsh_exe))
    return shutil.which("gmsh")


def _format_physical_groups(physical_groups: Mapping[str, Any] | None) -> list[str]:
    """Format Physical Volume / Physical Surface directives for a ``.geo`` script."""
    if not physical_groups:
        return []

    lines: list[str] = []
    if any(k in physical_groups for k in ("volume", "volumes", "surface", "surfaces")):
        vol_map = physical_groups.get("volume") or physical_groups.get("volumes") or {}
        surf_map = physical_groups.get("surface") or physical_groups.get("surfaces") or {}
    else:
        vol_map = physical_groups
        surf_map = {}

    for name, tags in vol_map.items():
        tag_list = ", ".join(str(int(t)) for t in tags)
        if tag_list:
            lines.append(f'Physical Volume("{name}") = {{{tag_list}}};')

    for name, tags in surf_map.items():
        tag_list = ", ".join(str(int(t)) for t in tags)
        if tag_list:
            lines.append(f'Physical Surface("{name}") = {{{tag_list}}};')

    return lines


def build_geo_script(
    cad_path: str | Path,
    out_msh_path: str | Path,
    options: GmshOptions | None = None,
    physical_groups: Mapping[str, Any] | None = None,
) -> str:
    """Build a Gmsh ``.geo`` script for meshing a STEP/IGES/BREP CAD file.

    Parameters
    ----------
    cad_path:
        Path to the input CAD file (``.step``, ``.stp``, ``.iges``, ``.igs``,
        or ``.brep``).
    out_msh_path:
        Output ``.msh`` file path.
    options:
        Meshing configuration. Defaults to :class:`GmshOptions()`.
    physical_groups:
        Optional physical group mapping. Supports either
        ``{"volume": {"PART": [1]}, "surface": {"LOAD": [2, 3]}}`` or a flat
        volume mapping ``{"PART": [1]}``.

    Returns
    -------
    str
        Complete ``.geo`` script text ready for ``gmsh <script.geo> -``.
    """
    opts = options if options is not None else GmshOptions()
    cad = Path(cad_path)
    if cad.suffix and cad.suffix.lower() not in SUPPORTED_CAD_EXTENSIONS:
        raise ValueError(
            f"unsupported CAD extension {cad.suffix!r}; "
            f"expected one of {sorted(SUPPORTED_CAD_EXTENSIONS)}"
        )

    out_msh = Path(out_msh_path)
    lines: list[str] = [
        'SetFactory("OpenCASCADE");',
        f'Merge "{cad.as_posix()}";',
        '_wht_surfs() = Surface "*";',
        '_wht_vols() = Volume "*";',
        "If (#_wht_vols() == 0 && #_wht_surfs() > 0)",
        "  Coherence;",
        '  _wht_surfs() = Surface "*";',
        "  _wht_sl = newsl;",
        "  Surface Loop(_wht_sl) = _wht_surfs();",
        "  Volume(1) = {_wht_sl};",
        "EndIf",
        f"Mesh.Algorithm3D = {opts.algorithm_3d};",
        f"Mesh.ElementOrder = {opts.order};",
        f"Mesh.Optimize = {1 if opts.optimize else 0};",
        f"General.NumThreads = {opts.threads};",
        f"Mesh.MshFileVersion = {_MSH_VERSION_MAP[opts.format]};",
    ]

    if opts.mesh_size_min is not None:
        lines.append(f"Mesh.CharacteristicLengthMin = {opts.mesh_size_min};")
    if opts.mesh_size_max is not None:
        lines.append(f"Mesh.CharacteristicLengthMax = {opts.mesh_size_max};")

    if opts.element_type == "hex8":
        # SubdivisionAlgorithm = 2 subdivides tetrahedra into all-hexahedra,
        # ensuring pure hex8 output on arbitrary 3D BREP/STEP solids.
        lines.append("Mesh.SubdivisionAlgorithm = 2;")

    lines.extend(_format_physical_groups(physical_groups))
    lines.append("Mesh 3;")
    lines.append(f'Save "{out_msh.as_posix()}";')
    return "\n".join(lines) + "\n"


def generate_mesh_with_gmsh(
    cad_path: str | Path,
    out_msh_path: str | Path,
    options: GmshOptions | None = None,
    gmsh_exe: str | Path | None = None,
    physical_groups: Mapping[str, Any] | None = None,
    geo_path: str | Path | None = None,
    timeout: float | None = 300.0,
) -> GmshResult:
    """Generate a 3D mesh from a CAD file by invoking the external ``gmsh`` CLI.

    Parameters
    ----------
    cad_path:
        Input STEP/IGES/BREP file path.
    out_msh_path:
        Desired output ``.msh`` file path.
    options:
        Meshing options. Defaults to :class:`GmshOptions()`.
    gmsh_exe:
        Optional explicit path/name of the ``gmsh`` executable. Defaults to
        ``shutil.which("gmsh")``.
    physical_groups:
        Optional physical volume/surface groups to include in the ``.geo`` script.
    geo_path:
        Optional path to save the generated ``.geo`` script. When ``None``, a
        temporary directory is used and cleaned up automatically.
    timeout:
        Subprocess timeout in seconds.

    Returns
    -------
    GmshResult
        Execution result including ``msh_path``, ``command``, ``returncode``,
        and any ``issues`` (e.g. ``GMSH_NOT_FOUND`` when ``gmsh`` is missing).
    """
    opts = options if options is not None else GmshOptions()
    issues: list[dict[str, str]] = []
    out_msh = Path(out_msh_path)
    cad = Path(cad_path)

    exe = find_gmsh_executable(gmsh_exe)
    if exe is None:
        target = str(gmsh_exe) if gmsh_exe is not None else "gmsh"
        issues.append(
            make_issue(
                "error",
                "GMSH_NOT_FOUND",
                f"External gmsh executable not found ({target!r}). "
                "Install gmsh CLI externally (GPL isolated; never imported).",
                entity=target,
            )
        )
        return GmshResult(
            msh_path=str(out_msh),
            command=[],
            returncode=127,
            issues=issues,
        )

    if not cad.is_file():
        issues.append(
            make_issue(
                "error",
                "CAD_NOT_FOUND",
                f"Input CAD file does not exist: {cad}",
                entity=str(cad),
            )
        )
        return GmshResult(
            msh_path=str(out_msh),
            command=[],
            returncode=2,
            issues=issues,
        )

    script = build_geo_script(
        cad_path=cad,
        out_msh_path=out_msh,
        options=opts,
        physical_groups=physical_groups,
    )
    out_msh.parent.mkdir(parents=True, exist_ok=True)

    if geo_path is not None:
        geo_file = Path(geo_path)
        geo_file.parent.mkdir(parents=True, exist_ok=True)
        geo_file.write_text(script, encoding="utf-8")
        return _run_gmsh_subprocess(exe, geo_file, out_msh, issues, timeout, keep_geo=True)

    with tempfile.TemporaryDirectory(prefix="wht_iga_gmsh_") as tmp_dir:
        geo_file = Path(tmp_dir) / "mesh_job.geo"
        geo_file.write_text(script, encoding="utf-8")
        return _run_gmsh_subprocess(exe, geo_file, out_msh, issues, timeout, keep_geo=False)


def _run_gmsh_subprocess(
    exe: str,
    geo_file: Path,
    out_msh: Path,
    issues: list[dict[str, str]],
    timeout: float | None,
    keep_geo: bool,
) -> GmshResult:
    """Execute ``gmsh <geo_file> -`` and collect status and issues."""
    cmd = [exe, str(geo_file), "-"]
    geo_record = str(geo_file) if keep_geo else ""
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        issues.append(
            make_issue(
                "error",
                "GMSH_TIMEOUT",
                f"Gmsh process timed out after {timeout}s",
                entity=str(out_msh),
            )
        )
        return GmshResult(
            msh_path=str(out_msh),
            command=cmd,
            returncode=-1,
            issues=issues,
            stdout=exc.stdout or "" if isinstance(exc.stdout, str) else "",
            stderr=exc.stderr or "" if isinstance(exc.stderr, str) else "",
            geo_path=geo_record,
        )
    except OSError as exc:
        issues.append(
            make_issue(
                "error",
                "GMSH_EXEC_ERROR",
                f"Failed to launch gmsh executable: {exc}",
                entity=exe,
            )
        )
        return GmshResult(
            msh_path=str(out_msh),
            command=cmd,
            returncode=-1,
            issues=issues,
            geo_path=geo_record,
        )

    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout).strip().splitlines()
        tail = detail[-1] if detail else f"exit code {proc.returncode}"
        issues.append(
            make_issue(
                "error",
                "GMSH_FAILED",
                f"Gmsh exited with code {proc.returncode}: {tail}",
                entity=str(out_msh),
            )
        )
    elif not out_msh.is_file():
        issues.append(
            make_issue(
                "error",
                "MSH_NOT_CREATED",
                f"Gmsh returned 0 but output file was not created: {out_msh}",
                entity=str(out_msh),
            )
        )

    return GmshResult(
        msh_path=str(out_msh),
        command=cmd,
        returncode=proc.returncode,
        issues=issues,
        stdout=proc.stdout,
        stderr=proc.stderr,
        geo_path=geo_record,
    )
