"""Locate and load the Kratos source build used by WHT_IGA.

Only this module imports Kratos. It pins the runtime to the source build at
``E:\\kratos-build\\install`` (Kratos 10.4.3 + IgaApplication); a pip Kratos
wheel found first on ``sys.path`` is rejected rather than silently used.

Environment variables:
    WHT_IGA_KRATOS_PATH      install root (default ``E:\\kratos-build\\install``)
    WHT_IGA_KRATOS_DLL_PATH  ``;``-separated DLL directories
                             (default: ``<root>\\libs`` and the conda deps used by the build)

Adapted from the pattern in ``aieng/src/aieng/topopt/adapters/kratos.py``
(copied, not imported: WHT_IGA must not import aieng).
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path
from types import ModuleType

DEFAULT_KRATOS_PATH = r"E:\kratos-build\install"
DEFAULT_DLL_PATHS = (
    r"E:\kratos-build\install\libs",
    r"D:\PythonCodeStudy\kratos-deps\condadeps\Library\bin",
)

DEFAULT_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
)

_dll_handles: list[object] = []
_configured = False


def kratos_root() -> Path:
    return Path(os.environ.get("WHT_IGA_KRATOS_PATH", DEFAULT_KRATOS_PATH))


def dll_dirs() -> list[Path]:
    raw = os.environ.get("WHT_IGA_KRATOS_DLL_PATH")
    if raw:
        return [Path(p) for p in raw.split(";") if p.strip()]
    root = kratos_root()
    return [root / "libs", Path(DEFAULT_DLL_PATHS[1])]


def configure_paths() -> None:
    """Put the Kratos build first on sys.path and register its DLL directories."""
    global _configured
    if _configured:
        return
    root = kratos_root()
    if not (root / "KratosMultiphysics" / "__init__.py").is_file():
        raise RuntimeError(
            f"Kratos build not found at {root}. Set WHT_IGA_KRATOS_PATH "
            "(see WHT_IGA/docs/ENVIRONMENT.md)."
        )
    for entry in (str(root / "libs"), str(root)):
        while entry in sys.path:
            sys.path.remove(entry)
        sys.path.insert(0, entry)
    if hasattr(os, "add_dll_directory"):
        for d in dll_dirs():
            if d.is_dir():
                _dll_handles.append(os.add_dll_directory(str(d)))
    _configured = True


def _check_origin(module: ModuleType) -> None:
    origin = Path(module.__file__).resolve()
    root = kratos_root().resolve()
    if root not in origin.parents:
        raise RuntimeError(
            f"wrong Kratos loaded: {origin} (expected under {root}). "
            "A pip Kratos wheel must not be installed in the whtiga env."
        )


def load_kratos(apps: tuple[str, ...] = DEFAULT_APPS) -> dict[str, ModuleType]:
    """Import Kratos core and the requested applications from the pinned build.

    Returns a dict mapping ``"KratosMultiphysics"`` and each app name to its module.
    """
    configure_paths()
    km = importlib.import_module("KratosMultiphysics")
    _check_origin(km)
    modules: dict[str, ModuleType] = {"KratosMultiphysics": km}
    for app in apps:
        modules[app] = importlib.import_module(f"KratosMultiphysics.{app}")
    return modules
