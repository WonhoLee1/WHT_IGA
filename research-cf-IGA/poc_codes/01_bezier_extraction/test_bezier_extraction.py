# -*- coding: utf-8 -*-
"""PoC 01: 1D p=2 B-spline Bezier extraction (splinepy 대체 검증).

Knots=[0,0,0,0.5,1,1,1], 요소 e0=[0,0.5] 에서
  N(x) = C^e B(t)  복원 -> K^e = C K_B C^T congruence -> Greville 패치테스트.
numpy만 사용. 전부 assert 통과해야 함.
"""
import numpy as np

KNOTS = np.array([0., 0., 0., 0.5, 1., 1., 1.])
P = 2
N_BASIS = len(KNOTS) - P - 1  # 4


def cox_deboor(x, i, p, knots=KNOTS):
    if p == 0:
        return 1.0 if (knots[i] <= x < knots[i + 1]
                       or (x == knots[-1] and knots[i + 1] == knots[-1])) else 0.0
    left = 0.0
    d1 = knots[i + p] - knots[i]
    if d1 > 0:
        left = (x - knots[i]) / d1 * cox_deboor(x, i, p - 1, knots)
    right = 0.0
    d2 = knots[i + p + 1] - knots[i + 1]
    if d2 > 0:
        right = (knots[i + p + 1] - x) / d2 * cox_deboor(x, i + 1, p - 1, knots)
    return left + right


def basis_vec(x):
    return np.array([cox_deboor(x, i, P) for i in range(N_BASIS)])


def bernstein(t):
    return np.array([(1 - t) ** 2, 2 * t * (1 - t), t ** 2])


def test_partition_of_unity():
    xs = np.linspace(0, 1, 401)
    s = np.array([basis_vec(x).sum() for x in xs])
    err = np.abs(s - 1.0).max()
    print(f"[partition-of-unity] max err = {err:.2e}")
    assert err < 1e-12


def extraction_matrix(e0=0.0, e1=0.5, ns=60):
    xs = np.linspace(e0, e1, ns)
    Npts = np.stack([basis_vec(x) for x in xs])          # ns x 4
    ts = (xs - e0) / (e1 - e0)
    Bpts = np.stack([bernstein(t) for t in ts])          # ns x 3
    C, *_ = np.linalg.lstsq(Bpts, Npts, rcond=None)      # B C^T = N^T
    return C.T                                           # 4 x 3


def test_extraction_reconstruction():
    C = extraction_matrix()
    assert C.shape == (4, 3)
    xs = np.linspace(0.0, 0.5, 121)
    err = max(np.abs(basis_vec(x) - C @ bernstein((x) / 0.5)).max() for x in xs)
    pou = np.abs(C.sum(axis=0) - np.ones(3)).max()  # 열합=1 (분할단위 보존)
    print(f"[extraction] C shape={C.shape}, recon err={err:.2e}, rowsum err={pou:.2e}")
    assert err < 1e-10 and pou < 1e-10
    return C


def test_stiffness_congruence(C):
    he = 0.5
    tg, wg = np.polynomial.legendre.leggauss(4)
    t = 0.5 * (tg + 1)            # [0,1]
    w = 0.5 * wg
    dB = np.stack([2 * t - 2, 2 - 4 * t, 2 * t]) / he   # dB/dx, 3 x nq
    K_B = (dB * w) @ dB.T * he
    K_cong = C @ K_B @ C.T
    # 직접 적분 (수치미분)
    eps = 1e-7
    xs = 0.5 * t
    Np = np.stack([basis_vec(x) for x in xs])
    dN = np.stack([(basis_vec(x + eps) - basis_vec(x - eps)) / (2 * eps) for x in xs])
    K_dir = ((dN * w[:, None]).T @ dN) * he
    err = np.abs(K_cong - K_dir).max()
    print(f"[stiffness] congruence err = {err:.2e}")
    assert err < 1e-6


def test_greville_patch():
    g = np.array([(KNOTS[i + 1] + KNOTS[i + 2]) / 2 for i in range(N_BASIS)])
    xs = np.linspace(0, 1, 401)
    err = max(abs(float(g @ basis_vec(x)) - x) for x in xs)
    print(f"[patch] Greville linear reproduction err = {err:.2e}")
    assert err < 1e-12


if __name__ == "__main__":
    test_partition_of_unity()
    C = test_extraction_reconstruction()
    test_stiffness_congruence(C)
    test_greville_patch()
    print("PoC-01 PASS")
