# -*- coding: utf-8 -*-
"""PoC 04: 1D 대칭 Nitsche 약형식 Dirichlet + 감마 민감도.

-u''=1 on (0,L), u(0)=0(강부과), u(L)=0(Nitsche, L=1+ηh 절단점).
γ 스위프: λmin 부호, L2 오차, κ 추적.
"""
import numpy as np

H = 0.25
XBG = np.arange(0, 1.2501, H)  # 배경노드 0..1.25 (6개)
ETA = 1e-3
L = 1.0 + ETA * H  # 절단점 (배경요소 [1.0,1.25] 내부)


def shape_at_LL():
    # 절단요소 [1.0, 1.25] 선형형상함수 at L (전역노드 4, 5)
    N4 = (XBG[-1] - L) / H
    N5 = (L - XBG[-2]) / H
    dN = np.array([-1 / H, 1 / H])
    return np.array([N4, N5]), dN


def assemble_nitsche(gamma):
    """비적합 조립: 배경요소 [1.0,1.25] 형상함수를 물리부분 [1.0,L]에서만 적분."""
    n = len(XBG)
    K = np.zeros((n, n))
    F = np.zeros(n)
    for i in range(n - 2):  # 온전 요소
        K[i:i + 2, i:i + 2] += np.array([[1, -1], [-1, 1]]) / H
        F[i:i + 2] += np.array([H / 2, H / 2])
    eta = (L - XBG[-2]) / H  # 활성 체적비 (형상함수 평가에 사용)
    K[-2:, -2:] += eta * np.array([[1, -1], [-1, 1]]) / H
    F[-2:] += np.array([eta * H / 2, eta * H / 2])
    N, dN = shape_at_LL()
    # 대칭 Nitsche: -u'(L)v(L) - u(L)v'(L) + γ/h u(L)v(L), 전역노드 4,5
    idx = [4, 5]
    for a in range(2):
        for b_ in range(2):
            K[idx[a], idx[b_]] += -dN[a] * N[b_] - N[a] * dN[b_] \
                + gamma / H * N[a] * N[b_]
    return K[1:, 1:], F[1:]


def l2_error(u_full, nq=400):
    xs = np.linspace(0, L, nq)
    uh = np.interp(xs, XBG, np.append([0.0], u_full))
    ue = xs * (L - xs) / 2
    return float(np.sqrt(np.trapezoid((uh - ue) ** 2, xs)))


def test_gamma_sweep():
    print(f"{'gamma':>8} {'lambda_min':>12} {'L2err':>10} {'kappa':>10}")
    rows = []
    for g in [1, 10, 100, 2000, 20000]:
        K, F = assemble_nitsche(g)
        w = np.linalg.eigvalsh(K)
        coer = w[0] > 0
        if coer:
            u = np.linalg.solve(K, F)
            e = l2_error(u)
            k = w[-1] / w[0]
        else:
            e, k = float("nan"), float("inf")
        rows.append((g, w[0], e, k))
        print(f"{g:>8} {w[0]:>12.3e} {e:>10.2e} {k:>10.2e}")
    assert rows[0][1] < 0, "소γ에서 강압성 상실(음 고유값)이 재현되어야 함"
    good = [r for r in rows if r[1] > 0]
    assert good, "충분한 γ에서 강압성 회복되어야 함"
    assert min(r[2] for r in good) < 1e-2, "적정 γ에서 오차가 작아야 함"
    assert rows[-1][3] > rows[-2][3], "과대 γ에서 조건수 단조 악화"


if __name__ == "__main__":
    test_gamma_sweep()
    print("PoC-04 PASS")
