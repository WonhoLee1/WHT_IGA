"""Kratos backend. Kratos is imported only through ``load_kratos``."""

from ._bootstrap import load_kratos, kratos_root

__all__ = ["load_kratos", "kratos_root"]
