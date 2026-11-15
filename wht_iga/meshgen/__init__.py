"""External-process Gmsh mesh generator package (D13, GPL-isolated).

This package never imports ``gmsh`` directly; it invokes the external ``gmsh``
executable via ``subprocess.run``.
"""

from __future__ import annotations

from wht_iga.meshgen.gmsh_runner import (
    GmshOptions,
    GmshResult,
    build_geo_script,
    find_gmsh_executable,
    generate_mesh_with_gmsh,
)

__all__ = [
    "GmshOptions",
    "GmshResult",
    "build_geo_script",
    "find_gmsh_executable",
    "generate_mesh_with_gmsh",
]
