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
    dirs = [str(d) for d in dll_dirs() if d.is_dir()]
    if hasattr(os, "add_dll_directory"):
        for d in dirs:
            _dll_handles.append(os.add_dll_directory(d))
    # Kratos delay-loads MKL (Pardiso) with the legacy search order, which
    # reads PATH, not add_dll_directory. Whatever mkl_rt PATH finds first wins:
    # - a shell without the env dirs on PATH: 0xC06D007F on the first Pardiso
    #   call (delay-load procedure not found), silent exit;
    # - condadeps first: Kratos gets MKL 2026.1 (libomp) while pypardiso globs
    #   sys.prefix and loads MKL 2026.2 (libiomp5md) -> "OMP: Error #15", abort.
    # So the env MKL dir leads: one mkl_rt, one OpenMP runtime per process,
    # shared by Kratos and pypardiso (2026-11-15).
    lead = [mkl_dir()] if mkl_dir() else []
    order = lead + [d for d in dirs if d not in lead]
    rest = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p and p not in order]
    os.environ["PATH"] = os.pathsep.join(order + rest)
    if lead and os.name == "nt":
        # PATH alone is not enough: Kratos' import resolves mkl_rt.3.dll through
        # the add_dll_directory dirs (condadeps) first. Preloading by full path
        # makes Windows reuse this module for Kratos' by-name import, and
        # PYPARDISO_MKL_RT points pypardiso at the same file.
        import ctypes

        rt = os.path.join(lead[0], "mkl_rt.3.dll")
        _dll_handles.append(ctypes.CDLL(rt))
        prev = os.environ.get("PYPARDISO_MKL_RT")
        if prev and os.path.normcase(os.path.abspath(prev)) != os.path.normcase(rt):
            import warnings

            warnings.warn(f"PYPARDISO_MKL_RT={prev!r} differs from the MKL Kratos uses ({rt}); "
                          "overriding it — two MKL copies in one process abort with OMP Error #15. "
                          "Set WHT_IGA_MKL_DIR to choose the MKL for both.", RuntimeWarning)
        os.environ["PYPARDISO_MKL_RT"] = rt
    _configured = True


def mkl_dir() -> str | None:
    """Directory whose mkl_rt.3.dll Kratos and pypardiso both load.

    ``WHT_IGA_MKL_DIR`` overrides; default ``<sys.prefix>/Library/bin`` (the
    whtiga env MKL, the one pypardiso finds by globbing ``sys.prefix``).
    """
    raw = os.environ.get("WHT_IGA_MKL_DIR")
    d = Path(raw) if raw else Path(sys.prefix) / "Library" / "bin"
    return str(d) if (d / "mkl_rt.3.dll").is_file() else None


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
