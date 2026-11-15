"""Unit and end-to-end tests for external-process Gmsh mesh generator (D13)."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import meshio
import pytest

from wht_iga.meshgen import (
    GmshOptions,
    GmshResult,
    build_geo_script,
    generate_mesh_with_gmsh,
)
from wht_iga.meshio_ext.read_msh import read_msh


def _write_ocp_box_step(step_path: Path, lx: float = 10.0, ly: float = 10.0, lz: float = 10.0) -> None:
    """Create a simple STEP box using OCP (no gmsh import)."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

    shape = BRepPrimAPI_MakeBox(lx, ly, lz).Shape()
    writer = STEPControl_Writer()
    writer.Transfer(shape, STEPControl_AsIs)
    status = writer.Write(str(step_path))
    assert status == 1, f"OCP STEP write failed with status {status}"


def test_gmsh_options_defaults_and_validation() -> None:
    """Verify GmshOptions defaults, order inference, and invalid parameter guards."""
    tet4 = GmshOptions(element_type="tet4")
    assert tet4.order == 1
    assert tet4.format == "msh2"

    tet10 = GmshOptions(element_type="tet10")
    assert tet10.order == 2

    hex8 = GmshOptions(element_type="hex8", format="msh4")
    assert hex8.order == 1
    assert hex8.format == "msh4"

    with pytest.raises(ValueError, match="element_type"):
        GmshOptions(element_type="prism6")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="order"):
        GmshOptions(element_type="tet4", order=2)

    with pytest.raises(ValueError, match="order"):
        GmshOptions(element_type="tet10", order=1)

    with pytest.raises(ValueError, match="mesh_size_min"):
        GmshOptions(mesh_size_min=5.0, mesh_size_max=2.0)

    with pytest.raises(ValueError, match="format"):
        GmshOptions(format="msh3")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="threads"):
        GmshOptions(threads=0)


def test_build_geo_script_element_types_and_physical_groups() -> None:
    """Verify .geo script generation for tet4, tet10, hex8, and physical groups."""
    geo_tet4 = build_geo_script(
        cad_path="parts/block.step",
        out_msh_path="out/block_tet4.msh",
        options=GmshOptions(
            element_type="tet4",
            mesh_size_min=1.5,
            mesh_size_max=4.0,
            algorithm_3d=4,
            optimize=True,
            threads=2,
            format="msh2",
        ),
        physical_groups={"volume": {"BODY": [1]}, "surface": {"FIXED": [1, 2]}},
    )
    assert 'Merge "parts/block.step";' in geo_tet4
    assert "Mesh.Algorithm3D = 4;" in geo_tet4
    assert "Mesh.ElementOrder = 1;" in geo_tet4
    assert "Mesh.Optimize = 1;" in geo_tet4
    assert "General.NumThreads = 2;" in geo_tet4
    assert "Mesh.MshFileVersion = 2.2;" in geo_tet4
    assert "Mesh.CharacteristicLengthMin = 1.5;" in geo_tet4
    assert "Mesh.CharacteristicLengthMax = 4.0;" in geo_tet4
    assert "Mesh.SubdivisionAlgorithm" not in geo_tet4
    assert 'Physical Volume("BODY") = {1};' in geo_tet4
    assert 'Physical Surface("FIXED") = {1, 2};' in geo_tet4
    assert "Mesh 3;" in geo_tet4
    assert 'Save "out/block_tet4.msh";' in geo_tet4

    geo_tet10 = build_geo_script(
        cad_path="parts/block.iges",
        out_msh_path="out/block_tet10.msh",
        options=GmshOptions(element_type="tet10", optimize=False, format="msh4"),
    )
    assert "Mesh.ElementOrder = 2;" in geo_tet10
    assert "Mesh.Optimize = 0;" in geo_tet10
    assert "Mesh.MshFileVersion = 4.1;" in geo_tet10

    geo_hex8 = build_geo_script(
        cad_path="parts/block.brep",
        out_msh_path="out/block_hex8.msh",
        options=GmshOptions(element_type="hex8"),
        physical_groups={"SOLID": [1]},
    )
    assert "Mesh.ElementOrder = 1;" in geo_hex8
    assert "Mesh.SubdivisionAlgorithm = 2;" in geo_hex8
    assert 'Physical Volume("SOLID") = {1};' in geo_hex8

    with pytest.raises(ValueError, match="unsupported CAD extension"):
        build_geo_script("parts/block.stl", "out/block.msh")


def test_gpl_isolation_no_gmsh_import() -> None:
    """Ensure importing wht_iga.meshgen never imports the GPL gmsh Python module."""
    assert "gmsh" not in sys.modules


def test_missing_gmsh_executable_returns_issue(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify clear error issue and non-zero status when gmsh binary is absent."""
    cad = tmp_path / "dummy.step"
    cad.write_text("ISO-10303-21;", encoding="utf-8")
    out_msh = tmp_path / "out.msh"

    res_explicit = generate_mesh_with_gmsh(
        cad_path=cad,
        out_msh_path=out_msh,
        gmsh_exe=tmp_path / "nonexistent_gmsh_bin",
    )
    assert isinstance(res_explicit, GmshResult)
    assert not res_explicit.ok
    assert res_explicit.returncode != 0
    assert len(res_explicit.issues) == 1
    assert res_explicit.issues[0]["severity"] == "error"
    assert res_explicit.issues[0]["code"] == "GMSH_NOT_FOUND"

    monkeypatch.setattr(shutil, "which", lambda _name: None)
    res_path = generate_mesh_with_gmsh(cad_path=cad, out_msh_path=out_msh)
    assert not res_path.ok
    assert res_path.issues[0]["code"] == "GMSH_NOT_FOUND"


def test_missing_cad_file_returns_issue(tmp_path: Path) -> None:
    """Verify CAD_NOT_FOUND issue when input CAD file does not exist."""
    res = generate_mesh_with_gmsh(
        cad_path=tmp_path / "missing.step",
        out_msh_path=tmp_path / "out.msh",
        gmsh_exe=sys.executable,  # any existing file to pass exe check
    )
    assert not res.ok
    assert res.issues[0]["code"] == "CAD_NOT_FOUND"


@pytest.mark.skipif(shutil.which("gmsh") is None, reason="external gmsh binary not installed")
@pytest.mark.parametrize(
    ("element_type", "expected_meshio_type", "expected_canon"),
    [
        ("tet4", "tetra", "tet4"),
        ("tet10", "tetra10", "tet10"),
        ("hex8", "hexahedron", "hex8"),
    ],
)
def test_end_to_end_ocp_step_to_gmsh_msh(
    tmp_path: Path,
    element_type: str,
    expected_meshio_type: str,
    expected_canon: str,
) -> None:
    """Generate a STEP box via OCP, mesh with external gmsh, and read via meshio."""
    step_path = tmp_path / "box.step"
    _write_ocp_box_step(step_path, 10.0, 10.0, 10.0)

    out_msh = tmp_path / f"box_{element_type}.msh"
    geo_path = tmp_path / f"box_{element_type}.geo"
    opts = GmshOptions(
        element_type=element_type,  # type: ignore[arg-type]
        mesh_size_min=3.0,
        mesh_size_max=5.0,
        threads=2,
        format="msh2",
    )
    res = generate_mesh_with_gmsh(
        cad_path=step_path,
        out_msh_path=out_msh,
        options=opts,
        physical_groups={"volume": {"SOLID": [1]}},
        geo_path=geo_path,
    )
    assert res.ok, f"gmsh failed: issues={res.issues}, stderr={res.stderr}"
    assert res.returncode == 0
    assert res.issues == []
    assert out_msh.is_file()
    assert geo_path.is_file()

    m = meshio.read(str(out_msh))
    cell_types = {cb.type: len(cb.data) for cb in m.cells}
    assert expected_meshio_type in cell_types
    assert cell_types[expected_meshio_type] > 0
    assert "SOLID" in m.field_data

    model, read_issues = read_msh(str(out_msh), unit_system="mm")
    assert not any(i["severity"] == "error" for i in read_issues)
    assert any(b.elem_type == expected_canon and b.conn.shape[0] > 0 for b in model.blocks)
    assert "SOLID" in model.parts.values()
