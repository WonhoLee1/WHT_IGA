"""Phase-1 linear static solver package (Track SOLVE).

Owns ``wht_iga/solver/*`` only. Reuses (read-only): ``wht_iga.grid``,
``wht_iga.reference``, ``wht_iga.kratos_backend`` (``load_kratos``,
``make_nurbs_volume`` + keepalive, ``elements.py``, ``constraints.py``,
``surface.py``). No Kratos import at module level anywhere in this package;
Kratos objects are passed in by the caller (pipeline/tests).
"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.1.0"
