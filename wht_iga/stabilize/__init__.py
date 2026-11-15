"""Stabilization operators for immersed IGA (Track STAB, Phase 2).

Owns ``wht_iga/stabilize/*`` only: small-cell basis-support classification
(T-P2.1), extended B-spline extension operator (T-P2.2), FCM ghost-alpha
baseline (T-P2.4) and the conditioning/V-COND tool (T-P2.5). Read-only reuse
of ``wht_iga.grid`` (knot conventions) and ``wht_iga.reference`` (oracle
matrices in tests). No Kratos import anywhere in this package.
"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.1.0"
