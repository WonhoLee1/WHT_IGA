"""Check that the whtiga environment can run WHT_IGA.

Usage:
    micromamba run -n whtiga python WHT_IGA/scripts/check_env.py

Exit code 0 when every check passes, 1 otherwise.
"""

from __future__ import annotations

import importlib
import sys

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))


def main() -> int:
    v = sys.version_info
    check("python 3.12", (v.major, v.minor) == (3, 12), sys.version.split()[0])

    # Import the numeric stack first: if a second MKL/HDF5 copy clashed with the
    # Kratos build's DLLs, the Kratos import below would fail.
    for mod in ("numpy", "scipy", "numba", "OCP", "h5py", "meshio", "pyvista", "pypardiso"):
        try:
            m = importlib.import_module(mod)
            check(f"import {mod}", True, str(getattr(m, "__version__", "")))
        except Exception as exc:  # noqa: BLE001 - report every failure
            check(f"import {mod}", False, f"{type(exc).__name__}: {exc}")

    try:
        import numba
        import numpy as np

        @numba.njit(cache=False)
        def _sum(a):
            s = 0.0
            for x in a:
                s += x
            return s

        val = _sum(np.arange(10.0))
        check("numba njit", val == 45.0, f"sum={val}")
    except Exception as exc:  # noqa: BLE001
        check("numba njit", False, f"{type(exc).__name__}: {exc}")

    try:
        from wht_iga.kratos_backend import kratos_root, load_kratos

        apps = (
            "IgaApplication",
            "StructuralMechanicsApplication",
            "LinearSolversApplication",
            "ConstitutiveLawsApplication",
            "ConvectionDiffusionApplication",
            "ContactStructuralMechanicsApplication",
        )
        mods = load_kratos(apps)
        km = mods["KratosMultiphysics"]
        check("kratos origin", True, f"{km.__file__} (root {kratos_root()})")
        check("kratos apps", True, ", ".join(apps))
        ls = mods["LinearSolversApplication"]
        has_mkl = bool(ls.HasMKL()) if hasattr(ls, "HasMKL") else False
        check("LinearSolvers HasMKL", has_mkl, str(has_mkl))
    except Exception as exc:  # noqa: BLE001
        check("kratos", False, f"{type(exc).__name__}: {exc}")

    width = max(len(n) for n, _, _ in RESULTS)
    failed = 0
    for name, ok, detail in RESULTS:
        failed += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {name:<{width}}  {detail}")
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
