"""T-P2.4: FCM ghost-alpha baseline + penalty beta sweep helpers.

FCM baseline (v2.1 §7.6): OUT-region elements carry a fictitious stiffness
``alpha * E`` with ``FCM_ALPHA = 1e-8``. This module owns the constant and the
scale law; the Kratos wiring stays Phase-1 owned (``solver.pipeline`` applies
``spec.ghost_alpha`` to OUT cells — read-only reuse, no modification here).

Beta law (v2 §5.4, ADR-002): beta = c_pen · E · p² / h with
``DEFAULT_C_PEN = 3e5`` fixed by the Track-SOLVE patch-cube sweep (exact 1/β,
C ≈ 4–8, slope −1.00 ± 0.01; V01 margins ≥ 3.5×). The canonical single source
for the solver default is ``wht_iga.solver.penalty.default_beta``; the copy
here exists so STAB-side sweep analysis never imports the SOLVE package
(cross-track import firewall — values asserted equal in tests).
"""

from __future__ import annotations

import numpy as np

#: Fictitious OUT-region stiffness scale (v2.1 §7.6 baseline).
FCM_ALPHA = 1e-8

#: Penalty constant fixed by ADR-002 (mirrors solver.penalty.DEFAULT_C_PEN).
DEFAULT_C_PEN = 3.0e5

#: Beta sweep multipliers (phase2 T-P2.4: {1e2 … 1e6}·E/h).
SWEEP_MULTS = (1e2, 1e3, 1e4, 1e5, 1e6)


def fictitious_modulus(E: float, alpha: float = FCM_ALPHA) -> float:
    """Ghost OUT-region modulus alpha·E (v2.1 §7.6)."""
    if E <= 0.0:
        raise ValueError(f"E must be positive, got {E}")
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    return float(alpha * E)


def default_beta(E: float, h: float, p: int,
                 c_pen: float = DEFAULT_C_PEN) -> float:
    """Penalty beta = c_pen · E · p² / h (mirrors solver.penalty.default_beta)."""
    if h <= 0.0:
        raise ValueError(f"h must be positive, got {h}")
    return float(c_pen * E * p * p / h)


def beta_sweep(E: float, h: float, p: int,
               mults: tuple[float, ...] = SWEEP_MULTS) -> dict[float, float]:
    """Sweep betas {mult}·E/h (T-P2.4; plateau/1/β analysis input)."""
    return {float(m): float(m * E / h) for m in mults}


def loglog_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Slope of log10(y) vs log10(x) (least squares; all entries > 0)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.any(x <= 0.0) or np.any(y <= 0.0):
        raise ValueError("log-log slope needs strictly positive x and y")
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])
