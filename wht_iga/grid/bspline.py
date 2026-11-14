"""B-spline basis evaluation (Numba kernels).

Piegl & Tiller, The NURBS Book, Alg. A2.1 (FindSpan) and A2.3 (DersBasisFuns).
All kernels take the standard clamped knot vector (see ``grid.knots``).
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def find_span(n_cp: int, p: int, u: float, knots: np.ndarray) -> int:
    """Index i of the knot span [knots[i], knots[i+1]) containing u (Alg. A2.1)."""
    if u >= knots[n_cp]:
        return n_cp - 1
    if u <= knots[p]:
        return p
    low = p
    high = n_cp
    mid = (low + high) // 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) // 2
    return mid


@njit(cache=True)
def ders_basis_funs(span: int, u: float, p: int, n_ders: int, knots: np.ndarray) -> np.ndarray:
    """Nonzero basis functions and derivatives at u (Alg. A2.3).

    Returns array (n_ders+1, p+1); row k is the k-th derivative of
    N_{span-p..span, p}(u).
    """
    ndu = np.zeros((p + 1, p + 1))
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)
    ndu[0, 0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            ndu[j, r] = right[r + 1] + left[j - r]
            temp = ndu[r, j - 1] / ndu[j, r]
            ndu[r, j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        ndu[j, j] = saved

    nd = min(n_ders, p)
    ders = np.zeros((n_ders + 1, p + 1))
    for j in range(p + 1):
        ders[0, j] = ndu[j, p]

    a = np.zeros((2, p + 1))
    for r in range(p + 1):
        s1 = 0
        s2 = 1
        a[0, 0] = 1.0
        for k in range(1, nd + 1):
            d = 0.0
            rk = r - k
            pk = p - k
            if r >= k:
                a[s2, 0] = a[s1, 0] / ndu[pk + 1, rk]
                d = a[s2, 0] * ndu[rk, pk]
            j1 = 1 if rk >= -1 else -rk
            j2 = k - 1 if r - 1 <= pk else p - r
            for j in range(j1, j2 + 1):
                a[s2, j] = (a[s1, j] - a[s1, j - 1]) / ndu[pk + 1, rk + j]
                d += a[s2, j] * ndu[rk + j, pk]
            if r <= pk:
                a[s2, k] = -a[s1, k - 1] / ndu[pk + 1, r]
                d += a[s2, k] * ndu[r, pk]
            ders[k, r] = d
            s1, s2 = s2, s1

    fac = float(p)
    for k in range(1, nd + 1):
        for j in range(p + 1):
            ders[k, j] *= fac
        fac *= p - k
    return ders


def basis_1d(knots: np.ndarray, p: int, u: float, n_ders: int = 0) -> tuple[int, np.ndarray]:
    """Convenience wrapper: (span, ders) for one parameter value."""
    knots = np.asarray(knots, dtype=float)
    n_cp = len(knots) - p - 1
    span = find_span(n_cp, p, float(u), knots)
    return span, ders_basis_funs(span, float(u), p, n_ders, knots)
