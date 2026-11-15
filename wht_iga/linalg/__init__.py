"""Linear-algebra helpers (tuned MKL Pardiso via pypardiso)."""

from wht_iga.linalg.pardiso import PardisoFactor, pardiso_solve

__all__ = ["PardisoFactor", "pardiso_solve"]
