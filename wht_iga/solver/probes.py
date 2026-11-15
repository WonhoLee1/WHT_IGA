"""Line-probe stress profiles and self-convergence (Coreform V02 convention).

Coreform roark-17p1-7a §7.5 assesses convergence on the transverse
max-principal-stress profile along the shortest span from the hole edge to
the free edge, as the relative L2 error toward the FINEST immersed solution
(self-convergence), because the body-fitted reference is itself approximate.

All levels are sampled at the same composite Gauss points so the discrete L2
norms are directly comparable.
"""

from __future__ import annotations

import numpy as np


def line_gauss_points(p0, p1, n_seg: int = 20, n_gauss: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Composite Gauss–Legendre points/weights on the segment p0→p1 (weights in length)."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    z, w = np.polynomial.legendre.leggauss(n_gauss)
    s = np.concatenate([(k + 0.5 * (z + 1.0)) / n_seg for k in range(n_seg)])
    ws = np.concatenate([0.5 * w / n_seg for _ in range(n_seg)]) * float(np.linalg.norm(p1 - p0))
    return p0[None, :] + s[:, None] * (p1 - p0)[None, :], ws


def max_principal_profile(box, U: np.ndarray, E: float, nu: float, P: np.ndarray) -> np.ndarray:
    """Max principal Cauchy stress of the discrete solution at points P (N,3)."""
    from wht_iga.solver.trimmed import _d_matrix, _strain_at_points, _tensor_from_voigt

    D = _d_matrix(float(E), float(nu))
    S = _strain_at_points(box, np.asarray(U, float), np.asarray(P, float))
    return np.array([float(np.linalg.eigvalsh(_tensor_from_voigt(D @ S[k]))[-1])
                     for k in range(S.shape[0])])


def rel_l2(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> float:
    """||a - b||_L2 / ||b||_L2 with quadrature weights w."""
    a, b, w = (np.asarray(v, float) for v in (a, b, w))
    return float(np.sqrt(np.sum(w * (a - b) ** 2) / np.sum(w * b ** 2)))


def self_convergence(hs, profiles, w) -> dict:
    """Errors of each level vs the finest profile and the fitted log-log rate.

    ``hs``/``profiles`` ordered coarse → fine; the finest is the reference and
    is excluded from the fit. Returns {errors, rate_all, rate_last3}.
    """
    hs = np.asarray(hs, float)
    ref = profiles[-1]
    errs = np.array([rel_l2(p, ref, w) for p in profiles[:-1]])
    out = {"errors": errs.tolist(), "rate_all": float("nan"), "rate_last3": float("nan")}
    if len(errs) >= 2 and np.all(errs > 0):
        out["rate_all"] = float(np.polyfit(np.log(hs[:-1]), np.log(errs), 1)[0])
    if len(errs) >= 3 and np.all(errs[-3:] > 0):
        out["rate_last3"] = float(np.polyfit(np.log(hs[:-1][-3:]), np.log(errs[-3:]), 1)[0])
    return out


def voigt_to_tensors(V: np.ndarray) -> np.ndarray:
    """(N,6) Kratos-order Voigt -> (N,3,3) symmetric tensors (no 1/2 on shear, F1)."""
    from wht_iga.solver.trimmed import _VOIGT as o

    V = np.asarray(V, float).reshape(-1, 6)
    T = np.empty((V.shape[0], 3, 3))
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = V[:, o["xx"]], V[:, o["yy"]], V[:, o["zz"]]
    T[:, 0, 1] = T[:, 1, 0] = V[:, o["xy"]]
    T[:, 1, 2] = T[:, 2, 1] = V[:, o["yz"]]
    T[:, 0, 2] = T[:, 2, 0] = V[:, o["xz"]]
    return T


def peak_max_principal_fast(cauchy_list) -> float:
    """Batched equivalent of ``trimmed.peak_max_principal`` (one eigvalsh call)."""
    parts = [np.asarray(a, float).reshape(-1, 6) for a in cauchy_list if np.size(a)]
    if not parts:
        return float("-inf")
    return float(np.linalg.eigvalsh(voigt_to_tensors(np.vstack(parts)))[:, -1].max())
