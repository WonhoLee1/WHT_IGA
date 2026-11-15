"""WHT_IGA must not import wht_topopt, aieng, or gmsh (X3, X10)."""

from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1] / "wht_iga"
FORBIDDEN = re.compile(r"^\s*(import|from)\s+(wht_topopt|aieng|gmsh)\b", re.MULTILINE)
FORBIDDEN_TOP_LEVEL = ("wht_topopt", "aieng", "gmsh")


def _find_forbidden_ast_imports(
    source: str, forbidden: tuple[str, ...] = FORBIDDEN_TOP_LEVEL
) -> list[tuple[int, str]]:
    """Return ``(lineno, module_root)`` pairs for forbidden imports in AST."""
    tree = ast.parse(source)
    hits: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".")[0]
                if root in forbidden:
                    hits.append((node.lineno, alias.name))
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            root = node.module.split(".")[0]
            if root in forbidden:
                hits.append((node.lineno, node.module))
    return hits


def test_no_forbidden_imports_in_source():
    offenders = [str(p) for p in PKG.rglob("*.py") if FORBIDDEN.search(p.read_text(encoding="utf-8"))]
    assert offenders == []


def test_no_forbidden_ast_imports_in_source():
    offenders: list[tuple[str, int, str]] = []
    for path in sorted(PKG.rglob("*.py")):
        hits = _find_forbidden_ast_imports(path.read_text(encoding="utf-8"))
        for lineno, mod in hits:
            offenders.append((str(path), lineno, mod))
    assert offenders == []


def test_ast_detector_catches_gmsh():
    sample = "import gmsh\nfrom gmsh import model\nimport wht_iga.meshgen.gmsh_runner\n"
    hits = _find_forbidden_ast_imports(sample)
    assert [(lineno, mod) for lineno, mod in hits] == [(1, "gmsh"), (2, "gmsh")]


def test_import_does_not_load_forbidden_or_kratos():
    code = (
        "import sys, wht_iga, wht_iga.grid.bspline, wht_iga.kratos_backend, "
        "wht_iga.meshgen.gmsh_runner;"
        "bad=[m for m in sys.modules if m.split('.')[0] in "
        "('wht_topopt','aieng','KratosMultiphysics','gmsh')];"
        "print(bad); sys.exit(1 if bad else 0)"
    )
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr
