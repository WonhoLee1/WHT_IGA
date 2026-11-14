"""WHT_IGA must not import wht_topopt or aieng (implementation request X3)."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1] / "wht_iga"
FORBIDDEN = re.compile(r"^\s*(import|from)\s+(wht_topopt|aieng)\b", re.MULTILINE)


def test_no_forbidden_imports_in_source():
    offenders = [str(p) for p in PKG.rglob("*.py") if FORBIDDEN.search(p.read_text(encoding="utf-8"))]
    assert offenders == []


def test_import_does_not_load_forbidden_or_kratos():
    code = (
        "import sys, wht_iga, wht_iga.grid.bspline, wht_iga.kratos_backend;"
        "bad=[m for m in sys.modules if m.split('.')[0] in ('wht_topopt','aieng','KratosMultiphysics')];"
        "print(bad); sys.exit(1 if bad else 0)"
    )
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr
