"""S01: environment spike — pytest port of ``scripts/check_env.py`` (T-P0.3).

Question (v2 §5.8): does the ``whtiga`` env import Kratos, Iga, SM,
LinearSolvers, ConstitutiveLaws, ConvDiff, ContactSM, with ``HasMKL()`` True?
"""

from __future__ import annotations

import importlib
import sys

import pytest

pytestmark = pytest.mark.kratos

NUMERIC_STACK = ["numpy", "scipy", "numba", "OCP", "h5py", "meshio", "pyvista", "pypardiso"]

KRATOS_APPS = (
    "IgaApplication",
    "StructuralMechanicsApplication",
    "LinearSolversApplication",
    "ConstitutiveLawsApplication",
    "ConvectionDiffusionApplication",
    "ContactStructuralMechanicsApplication",
)


def test_python_version() -> None:
    assert (sys.version_info.major, sys.version_info.minor) == (3, 12)


@pytest.mark.parametrize("mod", NUMERIC_STACK)
def test_numeric_stack_imports(mod: str) -> None:
    assert importlib.import_module(mod) is not None


def test_numba_njit_runs() -> None:
    import numba
    import numpy as np

    @numba.njit(cache=False)
    def _sum(a):  # type: ignore[no-untyped-def]
        s = 0.0
        for x in a:
            s += x
        return s

    assert _sum(np.arange(10.0)) == 45.0


def test_kratos_origin_pinned() -> None:
    from wht_iga.kratos_backend import kratos_root, load_kratos

    mods = load_kratos(KRATOS_APPS)
    km = mods["KratosMultiphysics"]
    assert "E:\\kratos-build\\install" in str(km.__file__).replace("/", "\\").lower() or "e:\\kratos-build\\install" in str(
        km.__file__
    ).lower()
    assert kratos_root().is_dir()


def test_kratos_apps_import() -> None:
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos(KRATOS_APPS)
    for app in KRATOS_APPS:
        assert app in mods


def test_has_mkl_recorded() -> None:
    """Records HasMKL(); passes only when the MKL build flag is set.

    Run output: HasMKL() == True (2026-11-14, recorded in docs/spikes/S01.md).
    """
    from wht_iga.kratos_backend import load_kratos

    mods = load_kratos(("LinearSolversApplication",))
    has_mkl = bool(mods["LinearSolversApplication"].HasMKL())
    print(f"\nHasMKL() = {has_mkl}")
    assert has_mkl
