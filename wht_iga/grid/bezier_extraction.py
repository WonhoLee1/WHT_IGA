"""Bézier extraction operators.

Borden, Scott, Evans, Hughes (2011), IJNME 87:15-47, Algorithm 1 (1D), and the
tensor-product rule C^e = C_w ⊗ C_v ⊗ C_u for trivariate splines.

Convention: for element e, ``N_e(xi) = C[e] @ B(xi)`` where ``N_e`` are the p+1
B-spline functions supported on the element and ``B`` the Bernstein basis on the
element's parametric interval mapped to [0, 1].
"""

from __future__ import annotations

import numpy as np


def extraction_1d(knots: np.ndarray, degree: int) -> np.ndarray:
    """1D extraction operators, shape (n_elements, p+1, p+1)."""
    U = np.asarray(knots, dtype=float)
    p = int(degree)
    m = len(U)
    a = p + 1
    b = a + 1
    nb = 0
    C = [np.eye(p + 1)]
    alphas = np.zeros(p + 1)
    while b < m:
        C.append(np.eye(p + 1))
        i = b
        while b < m and U[b] == U[b - 1]:
            b += 1
        mult = b - i + 1
        if mult < p:
            numer = U[b - 1] - U[a - 1]
            for j in range(p, mult, -1):
                alphas[j - mult - 1] = numer / (U[a + j - 1] - U[a - 1])
            r = p - mult
            for j in range(1, r + 1):
                save = r - j + 1
                s = mult + j
                for k in range(p + 1, s, -1):
                    alpha = alphas[k - s - 1]
                    C[nb][:, k - 1] = alpha * C[nb][:, k - 1] + (1.0 - alpha) * C[nb][:, k - 2]
                if b < m:
                    C[nb + 1][save - 1 : save + j, save - 1] = C[nb][p - j : p + 1, p]
        nb += 1
        if b < m:
            a = b
            b += 1
    n_el = len(np.unique(U)) - 1
    return np.array(C[:n_el])


def tensor_extraction(cu: np.ndarray, cv: np.ndarray, cw: np.ndarray) -> np.ndarray:
    """Trivariate element operator for one element: kron(C_w, kron(C_v, C_u)).

    Local function ordering: index = iu + (p_u+1)*(iv + (p_v+1)*iw), same for Bernstein.
    """
    return np.kron(cw, np.kron(cv, cu))


def bernstein(p: int, t: float, n_ders: int = 0) -> np.ndarray:
    """Bernstein polynomials of degree p on [0, 1] and derivatives, shape (n_ders+1, p+1)."""
    from math import comb

    out = np.zeros((n_ders + 1, p + 1))
    for k in range(n_ders + 1):
        if k > p:
            break
        for i in range(p + 1):
            val = 0.0
            # d^k/dt^k of C(p,i) t^i (1-t)^(p-i) via finite-difference formula on B_{p-k}
            fac = 1.0
            for q in range(k):
                fac *= p - q
            for j in range(k + 1):
                ii = i - j
                if 0 <= ii <= p - k:
                    sign = (-1) ** (k - j)
                    val += sign * comb(k, j) * comb(p - k, ii) * t**ii * (1.0 - t) ** (p - k - ii)
            out[k, i] = fac * val
    return out
