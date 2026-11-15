"""Operator composition + constraint numbering (Track STAB, Phase 2).

Owns ``wht_iga/operators/*`` only: the T = C_THB · E composition path
(T-P2.3) with Phase-2 C_THB = I (no THB yet — the composition path is built,
not the hierarchy), and batch ``LinearMasterSlaveConstraint`` spec generation
with c = 0 enforcement and chained-constraint rejection. Sparse-only
throughout (COO → CSR; dense temporaries are refused).
"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.1.0"
