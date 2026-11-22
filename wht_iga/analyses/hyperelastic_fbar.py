# -*- coding: utf-8 -*-
"""Mean-dilatation (F-bar / P0) Total-Lagrangian polynomial hyperelasticity, Numba (S48).

Abaqus convention, compressible, isochoric/volumetric split::

    W(F) = W_iso(I1b, I2b) + U(J)
    W_iso = sum_{1 <= i+j <= 3} C_ij (I1b - 3)^i (I2b - 3)^j
    U(J)  = sum_{k=1..3} (1/D_k) (J - 1)^(2k)

    I1b = J^(-2/3) I1,  I2b = J^(-4/3) I2  (isochoric invariants of C = F^T F)

* Mooney-Rivlin: C10, C01, D1 (``mooney_rivlin``).
* Yeoh (reduced polynomial N=3): C10, C20, C30, D1, D2, D3 (``yeoh``).
* Neo-Hookean: Mooney-Rivlin with C01 = 0 (``neo_hookean``).
* Initial moduli: mu0 = 2 (C10 + C01), K0 = 2 / D1.

Mean dilatation (Simo-Taylor-Pister, same structure as ``fbar.py`` / S40):
the isochoric energy is integrated pointwise (I1b, I2b are invariant under
the F-bar substitution F -> (Jbar/J)^(1/3) F), the volumetric energy sees
only the cell average ``Jbar = (1/V_c) sum_q w_q J_q``::

    Pi_cell = sum_q w_q W_iso(F_q) + V_c U(Jbar)
    R_I     = sum_q w_q P_iso(F_q) G_qI + U'(Jbar) V_I,   V_I = sum_q w_q J_q F_q^-T G_qI
    K       = sum_q w_q G_qI . A_iso(F_q) . G_qJ
              + U''(Jbar) V_I (x) V_J / V_c + U'(Jbar) dV_I/du_J

The pressure P0 is condensed in closed form (p = U'(Jbar)); the tangent is
the exact Hessian of Pi, hence symmetric. ``A_iso = d^2 W_iso / dF dF`` is
built from the closed-form second derivatives of I1b and I2b (verified
against central finite differences and against the complex-step gradient
of an independently coded energy in ``tests/unit/test_hyperelastic_fbar.py``).

Cauchy stress at a QP (the stress the weak form equilibrates)::

    sigma_q = P_iso(F_q) F_q^T / J_q + U'(Jbar) I

For homogeneous deformation Jbar = J and this is the Abaqus hyperelastic
Cauchy stress.

Geometry, quadrature and the CSR pattern come from ``fbar.build_model``
(imported, not modified). ``solve_static`` is the adaptive Newton of
``fbar.solve_static`` (copied and adapted: this assembler, plus a
linearised Dirichlet predictor for prescribed displacements).
Scope: untrimmed single-patch boxes, any clamped continuity C^(p-1)..C^0.
Ogden is not implemented (S48: open).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from wht_iga.analyses import fbar as _fbar
from wht_iga.analyses.fbar import _det3, _inv3

_MAX_ORDER = 3  # i + j <= 3 and k <= 3 (covers MR, Yeoh, Abaqus polynomial N <= 3)


@dataclass(frozen=True)
class HyperelasticMaterial:
    """Polynomial hyperelastic law (Abaqus convention).

    ``C[i, j]`` multiplies ``(I1b - 3)^i (I2b - 3)^j`` (``C[0, 0]`` unused);
    ``inv_D[k-1] = 1/D_k`` (0 means the term is absent).
    """

    name: str
    C: np.ndarray
    inv_D: np.ndarray

    @property
    def mu0(self) -> float:
        """Initial shear modulus 2 (C10 + C01)."""
        return 2.0 * float(self.C[1, 0] + self.C[0, 1])

    @property
    def K0(self) -> float:
        """Initial bulk modulus 2 / D1."""
        return 2.0 * float(self.inv_D[0])


def polynomial(coeffs: dict, D: tuple, name: str = "polynomial") -> HyperelasticMaterial:
    """General polynomial law, ``coeffs = {(i, j): C_ij}`` with 1 <= i+j <= 3.

    ``D = (D1, D2, D3)`` (shorter tuples allowed); D1 must be finite and > 0
    (no hybrid/incompressible formulation: a small D1 approximates it),
    D2/D3 may be ``None`` (term absent).
    """
    C = np.zeros((_MAX_ORDER + 1, _MAX_ORDER + 1), dtype=np.float64)
    for (i, j), v in coeffs.items():
        i, j = int(i), int(j)
        if i < 0 or j < 0 or not 1 <= i + j <= _MAX_ORDER:
            raise ValueError(f"C{i}{j}: need 1 <= i + j <= {_MAX_ORDER}")
        C[i, j] = float(v)
    D = tuple(D)
    if not 1 <= len(D) <= _MAX_ORDER:
        raise ValueError("D must hold 1..3 values")
    inv_D = np.zeros(_MAX_ORDER, dtype=np.float64)
    for k, d in enumerate(D):
        if d is None:
            if k == 0:
                raise ValueError("D1 is required (compressible formulation)")
            continue
        d = float(d)
        if not np.isfinite(d) or d <= 0.0:
            raise ValueError(f"D{k + 1} must be finite and > 0, got {d}")
        inv_D[k] = 1.0 / d
    mat = HyperelasticMaterial(name=name, C=C, inv_D=inv_D)
    if not mat.mu0 > 0.0:
        raise ValueError(f"initial shear modulus 2(C10 + C01) = {mat.mu0} must be > 0")
    return mat


def mooney_rivlin(C10: float, C01: float, D1: float) -> HyperelasticMaterial:
    """W = C10 (I1b-3) + C01 (I2b-3) + (1/D1)(J-1)^2."""
    return polynomial({(1, 0): C10, (0, 1): C01}, (D1,), name="mooney_rivlin")


def neo_hookean(C10: float, D1: float) -> HyperelasticMaterial:
    """W = C10 (I1b-3) + (1/D1)(J-1)^2 (Abaqus NEO HOOKE)."""
    return polynomial({(1, 0): C10}, (D1,), name="neo_hookean")


def yeoh(C10: float, C20: float, C30: float, D1: float,
         D2: float | None = None, D3: float | None = None) -> HyperelasticMaterial:
    """W = sum_i C_i0 (I1b-3)^i + sum_k (1/D_k)(J-1)^(2k), i, k = 1..3."""
    return polynomial({(1, 0): C10, (2, 0): C20, (3, 0): C30}, (D1, D2, D3), name="yeoh")


# --------------------------------------------------------------------------
# Numba kernels
# --------------------------------------------------------------------------

@njit(cache=True)
def _iso_derivs(Cc: np.ndarray, I1b: float, I2b: float):
    """dW/dI1b, dW/dI2b and the three second derivatives of W_iso."""
    a = I1b - 3.0
    b = I2b - 3.0
    W1 = 0.0
    W2 = 0.0
    W11 = 0.0
    W12 = 0.0
    W22 = 0.0
    n = Cc.shape[0]
    for i in range(n):
        for j in range(n):
            c = Cc[i, j]
            if c == 0.0 or (i == 0 and j == 0):
                continue
            if i >= 1:
                W1 += c * i * a ** (i - 1) * b ** j
            if j >= 1:
                W2 += c * j * a ** i * b ** (j - 1)
            if i >= 2:
                W11 += c * i * (i - 1) * a ** (i - 2) * b ** j
            if i >= 1 and j >= 1:
                W12 += c * i * j * a ** (i - 1) * b ** (j - 1)
            if j >= 2:
                W22 += c * j * (j - 1) * a ** i * b ** (j - 2)
    return W1, W2, W11, W12, W22


@njit(cache=True)
def _vol_derivs(inv_D: np.ndarray, J: float):
    """U'(J), U''(J) for U = sum_k inv_D[k-1] (J-1)^(2k)."""
    d = J - 1.0
    U1 = 0.0
    U2 = 0.0
    for k in range(inv_D.shape[0]):
        c = inv_D[k]
        if c == 0.0:
            continue
        e = 2 * (k + 1)
        U1 += c * e * d ** (e - 1)
        U2 += c * e * (e - 1) * d ** (e - 2)
    return U1, U2


@njit(cache=True)
def _iso_PA(F: np.ndarray, Cc: np.ndarray, want_A: bool):
    """First Piola P_iso[a, M] and tangent A[a, M, b, N] = d^2 W_iso / dF_aM dF_bN.

    Closed-form derivatives (H = F^-T, dH_aM/dF_bN = -H_aN H_bM)::

        dI1b = J^-2/3 (2F - 2/3 I1 H)
        dI2  = 2 (I1 F - F C),  dI2b = J^-4/3 (dI2 - 4/3 I2 H)
        d2I1b = J^-2/3 [2 d_ab d_MN - 4/3 (F_aM H_bN + H_aM F_bN)
                        + 4/9 I1 H_aM H_bN + 2/3 I1 H_aN H_bM]
        d2I2  = 2 [2 F_aM F_bN + I1 d_ab d_MN - d_ab C_MN - F_aN F_bM - B_ab d_MN]
        d2I2b = J^-4/3 [d2I2 - 4/3 (dI2_aM H_bN + H_aM dI2_bN)
                        + 16/9 I2 H_aM H_bN + 4/3 I2 H_aN H_bM]
    """
    J = _det3(F)
    Finv = _inv3(F)
    H = np.empty((3, 3))
    Cm = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            H[i, j] = Finv[j, i]
            for k in range(3):
                Cm[i, j] += F[k, i] * F[k, j]
                Bm[i, j] += F[i, k] * F[j, k]
    FC = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                FC[i, j] += F[i, k] * Cm[k, j]
    I1 = Cm[0, 0] + Cm[1, 1] + Cm[2, 2]
    trC2 = 0.0
    for i in range(3):
        for j in range(3):
            trC2 += Cm[i, j] * Cm[i, j]
    I2 = 0.5 * (I1 * I1 - trC2)
    j23 = J ** (-2.0 / 3.0)
    j43 = j23 * j23
    I1b = j23 * I1
    I2b = j43 * I2
    W1, W2, W11, W12, W22 = _iso_derivs(Cc, I1b, I2b)
    dI1b = np.empty((3, 3))
    dI2 = np.empty((3, 3))
    dI2b = np.empty((3, 3))
    P = np.empty((3, 3))
    for a in range(3):
        for M in range(3):
            dI1b[a, M] = j23 * (2.0 * F[a, M] - (2.0 / 3.0) * I1 * H[a, M])
            dI2[a, M] = 2.0 * (I1 * F[a, M] - FC[a, M])
            dI2b[a, M] = j43 * (dI2[a, M] - (4.0 / 3.0) * I2 * H[a, M])
            P[a, M] = W1 * dI1b[a, M] + W2 * dI2b[a, M]
    A = np.zeros((3, 3, 3, 3))
    if not want_A:
        return P, A
    for a in range(3):
        for M in range(3):
            for b in range(3):
                for N in range(3):
                    dab = 1.0 if a == b else 0.0
                    dMN = 1.0 if M == N else 0.0
                    hh = H[a, M] * H[b, N]
                    hx = H[a, N] * H[b, M]
                    d1 = (2.0 * dab * dMN
                          - (4.0 / 3.0) * (F[a, M] * H[b, N] + H[a, M] * F[b, N])
                          + (4.0 / 9.0) * I1 * hh + (2.0 / 3.0) * I1 * hx)
                    d2I2 = 2.0 * (2.0 * F[a, M] * F[b, N] + I1 * dab * dMN - dab * Cm[M, N]
                                  - F[a, N] * F[b, M] - Bm[a, b] * dMN)
                    d2 = (d2I2 - (4.0 / 3.0) * (dI2[a, M] * H[b, N] + H[a, M] * dI2[b, N])
                          + (16.0 / 9.0) * I2 * hh + (4.0 / 3.0) * I2 * hx)
                    A[a, M, b, N] = (W1 * j23 * d1 + W2 * j43 * d2
                                     + W11 * dI1b[a, M] * dI1b[b, N]
                                     + W12 * (dI1b[a, M] * dI2b[b, N] + dI2b[a, M] * dI1b[b, N])
                                     + W22 * dI2b[a, M] * dI2b[b, N])
    return P, A


@njit(cache=True)
def _cell_residual_tangent(G: np.ndarray, W: np.ndarray, u_cell: np.ndarray,
                           Cc: np.ndarray, inv_D: np.ndarray):
    """Residual (3*nen,) and tangent (3*nen, 3*nen) for one cell, and Jbar.

    G: (nq, nen, 3) reference gradients dN/dX. W: (nq,) reference weights.
    Inadmissible state (any J_q <= 0 or Jbar <= 0) returns a NaN residual so
    the Newton guard in ``solve_static`` shortens the step.
    """
    nq, nen, _ = G.shape
    nd = 3 * nen
    R = np.zeros(nd, dtype=np.float64)
    K = np.zeros((nd, nd), dtype=np.float64)
    F = np.zeros((nq, 3, 3), dtype=np.float64)
    detJ = np.zeros(nq, dtype=np.float64)
    for q in range(nq):
        for a in range(3):
            F[q, a, a] = 1.0
            for I in range(nen):
                for m in range(3):
                    F[q, a, m] += u_cell[I, a] * G[q, I, m]
        detJ[q] = _det3(F[q])
    Wcell = 0.0
    Jbar = 0.0
    bad = False
    for q in range(nq):
        Wcell += W[q]
        Jbar += W[q] * detJ[q]
        if not detJ[q] > 0.0:
            bad = True
    Jbar /= Wcell
    if bad or not Jbar > 0.0:
        R[:] = np.nan
        return R, K, Jbar
    U1, U2 = _vol_derivs(inv_D, Jbar)
    # H[q, I] = F^-T G_qI (V, dJbar and the dV/du term, as in fbar)
    H = np.zeros((nq, nen, 3), dtype=np.float64)
    for q in range(nq):
        Finv = _inv3(F[q])
        for I in range(nen):
            for a in range(3):
                h = 0.0
                for m in range(3):
                    h += Finv[m, a] * G[q, I, m]  # (F^-T)[a, m] = Finv[m, a]
                H[q, I, a] = h
    V = np.zeros((nen, 3), dtype=np.float64)
    for I in range(nen):
        for q in range(nq):
            wJ = W[q] * detJ[q]
            for a in range(3):
                V[I, a] += wJ * H[q, I, a]
    # isochoric part, pointwise
    T = np.empty((nen, 3, 3, 3), dtype=np.float64)
    for q in range(nq):
        P, A = _iso_PA(F[q], Cc, True)
        wq = W[q]
        for I in range(nen):
            for a in range(3):
                s = 0.0
                for M in range(3):
                    s += P[a, M] * G[q, I, M]
                R[3 * I + a] += wq * s
        for I in range(nen):
            for a in range(3):
                for b in range(3):
                    for N in range(3):
                        s = 0.0
                        for M in range(3):
                            s += G[q, I, M] * A[a, M, b, N]
                        T[I, a, b, N] = wq * s
        for I in range(nen):
            for J in range(nen):
                for a in range(3):
                    for b in range(3):
                        s = 0.0
                        for N in range(3):
                            s += T[I, a, b, N] * G[q, J, N]
                        K[3 * I + a, 3 * J + b] += s
    # volumetric part with the cell Jbar
    for I in range(nen):
        for a in range(3):
            R[3 * I + a] += U1 * V[I, a]
    for I in range(nen):
        for J in range(nen):
            for a in range(3):
                for b in range(3):
                    K[3 * I + a, 3 * J + b] += U2 * V[I, a] * V[J, b] / Wcell
    # U'(Jbar) dV/du: block (I, J) = sum_q w_q J_q (H_I (x) H_J - H_J (x) H_I),
    # antisymmetric 3x3 with zero diagonal (fbar.py, FD-verified there and here)
    for I in range(nen):
        r0, r1, r2 = 3 * I, 3 * I + 1, 3 * I + 2
        for J in range(I + 1, nen):
            c0, c1, c2 = 3 * J, 3 * J + 1, 3 * J + 2
            for q in range(nq):
                cq = U1 * W[q] * detJ[q]
                hi0, hi1, hi2 = H[q, I, 0], H[q, I, 1], H[q, I, 2]
                hj0, hj1, hj2 = H[q, J, 0], H[q, J, 1], H[q, J, 2]
                v01 = cq * (hi0 * hj1 - hj0 * hi1)
                v02 = cq * (hi0 * hj2 - hj0 * hi2)
                v12 = cq * (hi1 * hj2 - hj1 * hi2)
                K[r0, c1] += v01
                K[r1, c0] -= v01
                K[r0, c2] += v02
                K[r2, c0] -= v02
                K[r1, c2] += v12
                K[r2, c1] -= v12
                K[c1, r0] += v01
                K[c0, r1] -= v01
                K[c2, r0] += v02
                K[c0, r2] -= v02
                K[c2, r1] += v12
                K[c1, r2] -= v12
    return R, K, Jbar


@njit(cache=True)
def _assemble_csr(G_all: np.ndarray, W_all: np.ndarray, conn: np.ndarray, u: np.ndarray,
                  Cc: np.ndarray, inv_D: np.ndarray, indptr: np.ndarray,
                  node_pos: np.ndarray, nnz: int):
    """Global R and CSR data of K (pattern from ``fbar._csr_pattern``; same
    scatter as ``fbar._assemble_csr``, only the cell kernel differs)."""
    nc, nen = conn.shape
    R = np.zeros(u.shape[0], dtype=np.float64)
    data = np.zeros(nnz, dtype=np.float64)
    uc = np.empty((nen, 3), dtype=np.float64)
    for c in range(nc):
        for I in range(nen):
            n = conn[c, I]
            for a in range(3):
                uc[I, a] = u[3 * n + a]
        Rc, Kc, _ = _cell_residual_tangent(G_all[c], W_all[c], uc, Cc, inv_D)
        for I in range(nen):
            n = conn[c, I]
            for a in range(3):
                R[3 * n + a] += Rc[3 * I + a]
                rs = indptr[3 * n + a]
                for J in range(nen):
                    p0 = rs + 3 * node_pos[c, I, J]
                    for b in range(3):
                        data[p0 + b] += Kc[3 * I + a, 3 * J + b]
    return R, data


@njit(cache=True)
def _cauchy_all(G_all: np.ndarray, W_all: np.ndarray, conn: np.ndarray, u: np.ndarray,
                Cc: np.ndarray, inv_D: np.ndarray):
    """(nc, nq, 6) Cauchy stress xx, yy, zz, xy, yz, xz and (nc,) Jbar."""
    nc, nen = conn.shape
    nq = W_all.shape[1]
    S = np.zeros((nc, nq, 6), dtype=np.float64)
    Jb = np.zeros(nc, dtype=np.float64)
    F = np.zeros((nq, 3, 3), dtype=np.float64)
    detJ = np.zeros(nq, dtype=np.float64)
    for c in range(nc):
        F[:] = 0.0
        Wcell = 0.0
        Jbar = 0.0
        for q in range(nq):
            for a in range(3):
                F[q, a, a] = 1.0
                for I in range(nen):
                    n = conn[c, I]
                    for m in range(3):
                        F[q, a, m] += u[3 * n + a] * G_all[c, q, I, m]
            detJ[q] = _det3(F[q])
            Wcell += W_all[c, q]
            Jbar += W_all[c, q] * detJ[q]
        Jbar /= Wcell
        Jb[c] = Jbar
        U1, _ = _vol_derivs(inv_D, Jbar)
        sig = np.zeros((3, 3), dtype=np.float64)
        for q in range(nq):
            P, _A = _iso_PA(F[q], Cc, False)
            for a in range(3):
                for b in range(3):
                    s = 0.0
                    for m in range(3):
                        s += P[a, m] * F[q, b, m]
                    sig[a, b] = s / detJ[q]
                sig[a, a] += U1
            S[c, q, 0] = sig[0, 0]
            S[c, q, 1] = sig[1, 1]
            S[c, q, 2] = sig[2, 2]
            S[c, q, 3] = sig[0, 1]
            S[c, q, 4] = sig[1, 2]
            S[c, q, 5] = sig[0, 2]
    return S, Jb


# --------------------------------------------------------------------------
# Python API
# --------------------------------------------------------------------------

def build_model(box, degree: int, material: HyperelasticMaterial,
                n_gauss: int | None = None) -> dict:
    """Per-cell QP data + CSR pattern from ``fbar.build_model`` plus *material*.

    ``fbar.build_model`` needs (E, nu) only for its own NH law; the dummy
    values are dropped so nothing here can read them.
    """
    if not isinstance(material, HyperelasticMaterial):
        raise TypeError("material must be a HyperelasticMaterial (mooney_rivlin / yeoh / ...)")
    model = _fbar.build_model(box, degree, E=1.0, nu=0.0, n_gauss=n_gauss)
    model.pop("mu")
    model.pop("lam")
    model["material"] = material
    model["C"] = np.ascontiguousarray(material.C, dtype=np.float64)
    model["inv_D"] = np.ascontiguousarray(material.inv_D, dtype=np.float64)
    model["formulation"] = "fbar_mean_dilatation_polynomial_hyperelastic"
    return model


def residual_tangent_sparse(model: dict, u: np.ndarray):
    """Assemble global internal force R (ndarray) and K (scipy CSR, sorted indices)."""
    import scipy.sparse as sp

    ndof = int(model["n_dof"])
    s = model["csr"]
    R, data = _assemble_csr(model["G_all"], model["W_all"], model["conn"],
                            np.ascontiguousarray(u, dtype=np.float64).reshape(ndof),
                            model["C"], model["inv_D"],
                            s["indptr"], s["node_pos"], s["nnz"])
    K = sp.csr_matrix((data, s["indices"], s["indptr"]), shape=(ndof, ndof))
    K.has_sorted_indices = True
    return R, K


def cauchy_stress(model: dict, u: np.ndarray) -> list:
    """Per-QP Cauchy stress, list over cells of (nq, 6) Voigt xx,yy,zz,xy,yz,xz
    (QP order q = gx + ng*(gy + ng*gz), as ``fbar.cauchy_stress``)."""
    S, _ = _cauchy_all(model["G_all"], model["W_all"], model["conn"],
                       np.ascontiguousarray(u, dtype=np.float64).reshape(int(model["n_dof"])),
                       model["C"], model["inv_D"])
    return list(S)


def cell_jbar(model: dict, u: np.ndarray) -> np.ndarray:
    """Cell-averaged Jacobian Jbar per cell (pressure p = -U'(Jbar))."""
    _, Jb = _cauchy_all(model["G_all"], model["W_all"], model["conn"],
                        np.ascontiguousarray(u, dtype=np.float64).reshape(int(model["n_dof"])),
                        model["C"], model["inv_D"])
    return Jb


def solve_static(model: dict, f_ext: np.ndarray, fixed: np.ndarray,
                 n_steps: int = 1, rtol: float = 1e-6, atol: float = 1e-9,
                 max_it: int = 25, prescribed: dict | None = None,
                 dtol: float = 1e-10, min_frac: float = 1e-4,
                 history: list | None = None,
                 linear_solver: str = "auto", predictor: bool = True) -> np.ndarray:
    """Adaptive load-stepped Newton (copied from ``fbar.solve_static``; see its
    docstring for the convergence tests, the guard-only line search and the
    increment halving). ``linear_solver``: ``"auto"`` | ``"pardiso"`` |
    ``"superlu"`` on the CSR tangent. Returns the full displacement vector.

    Difference to fbar: with ``prescribed`` DOFs each increment starts from
    the linearised Dirichlet predictor ``K_ff du_f = -(R_f - f_f) - K_fp du_p``
    (``predictor=False`` restores the fbar start: new boundary values on the
    old interior). The residual gate is then relative to that predictor
    right-hand side. Measured 2026-11-20 (S48, 2x2x1 cells p=2, homogeneous
    uniaxial, Mooney-Rivlin): lam = 0.7 (n_steps 4) 15 -> 4 increments
    (88 -> 12 iterations); K0/mu0 = 1e4 (n_steps 8) 80 -> 8 increments
    (623 -> 24 iterations); converged states equal to <= 1.7e-13 relative.
    """
    ndof = int(model["n_dof"])
    u = np.zeros(ndof)
    fixed = np.asarray(fixed, dtype=np.int64)
    prescribed = {int(k): float(v) for k, v in (prescribed or {}).items()}
    free = np.ones(ndof, dtype=bool)
    free[fixed] = False
    for k in list(prescribed.keys()):
        if k < 0 or k >= ndof:
            raise ValueError(f"prescribed dof {k} out of range")
        free[k] = False
    free_idx = np.flatnonzero(free)
    f = np.asarray(f_ext, dtype=float)
    kind = str(linear_solver).lower()
    if kind not in ("auto", "pardiso", "superlu"):
        raise ValueError(f"unknown linear_solver {linear_solver!r}")
    pardiso = None
    if kind != "superlu":
        try:
            import pypardiso
            pardiso = pypardiso.PyPardisoSolver()
        except ImportError:
            if kind == "pardiso":
                raise

    def lin_solve(K, b):
        A = K[free_idx][:, free_idx]
        if pardiso is not None:
            import pypardiso
            return np.asarray(pypardiso.spsolve(A, b, solver=pardiso), dtype=np.float64)
        import scipy.sparse.linalg as spla
        return np.asarray(spla.spsolve(A.tocsc(), b), dtype=np.float64)

    d_nom = 1.0 / max(1, int(n_steps))
    lam, dlam = 0.0, d_nom

    pres_idx = np.array(sorted(prescribed), dtype=np.int64)
    pres_val = np.array([prescribed[k] for k in pres_idx], dtype=np.float64)

    def _increment(lam_new: float) -> tuple[bool, int]:
        target = f * lam_new
        gate_ref = None
        if predictor and pres_idx.size and free_idx.size:
            # Linearised Dirichlet predictor from the converged state:
            # K_ff du_f = -(R_f - f_f) - K_fp du_p. Without it Newton starts
            # from the kinked state (new boundary values, old interior).
            R, K = residual_tangent_sparse(model, u)
            dp = pres_val * lam_new - u[pres_idx]
            rhs = -(R[free_idx] - target[free_idx]) - K[free_idx][:, pres_idx] @ dp
            gate_ref = float(np.max(np.abs(rhs)))
            u_try = u.copy()
            u_try[pres_idx] += dp
            u_try[free_idx] += lin_solve(K, rhs)
            R, K = residual_tangent_sparse(model, u_try)
            if not np.all(np.isfinite(R)):  # inadmissible prediction: plain start
                u_try[free_idx] = u[free_idx]
                R, K = residual_tangent_sparse(model, u_try)
            u[:] = u_try
        else:
            u[pres_idx] = pres_val * lam_new
            R, K = residual_tangent_sparse(model, u)
        if not np.all(np.isfinite(R)):
            return False, 0
        rr = R[free] - target[free]
        if rr.size == 0:
            return True, 0
        if gate_ref is None:
            gate_ref = float(np.max(np.abs(rr)))
        gate = atol + rtol * max(gate_ref, float(np.max(np.abs(target[free]))), 1e-300)
        for it in range(1, max_it + 1):
            if float(np.max(np.abs(rr))) < gate:
                return True, it - 1
            du = lin_solve(K, -rr)
            step = 1.0
            while True:
                u_try = u.copy()
                u_try[free] += step * du
                R, K = residual_tangent_sparse(model, u_try)
                if np.all(np.isfinite(R)):
                    break
                step *= 0.5
                if step < 1.0 / 1024.0:
                    return False, it
            u[:] = u_try
            rr = R[free] - target[free]
            if float(np.max(np.abs(step * du))) <= dtol * max(float(np.max(np.abs(u[free]))), 1e-300):
                return True, it
        return float(np.max(np.abs(rr))) < gate, max_it

    try:
        while lam < 1.0 - 1e-14:
            lam_new = min(1.0, lam + dlam)
            u_save = u.copy()
            ok, its = _increment(lam_new)
            if ok:
                lam = lam_new
                if history is not None:
                    history.append((lam, its))
                if its <= 6:
                    dlam = min(d_nom, 1.5 * dlam)
            else:
                u[:] = u_save
                dlam *= 0.5
                if dlam < min_frac:
                    raise RuntimeError(
                        f"hyperelastic fbar Newton failed: increment below {min_frac} at load factor {lam:.6g}")
    finally:
        if pardiso is not None:
            pardiso.free_memory(everything=True)
    return u
