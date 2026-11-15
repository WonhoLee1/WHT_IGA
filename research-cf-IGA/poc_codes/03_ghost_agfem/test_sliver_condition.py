# -*- coding: utf-8 -*-
"""PoC 03: 소형 절단셀 조건수 — 미안정화 폭발 vs Ghost vs AgFEM.

-u''=1 on (0,L), u(0)=0, 우단 자유(Neumann). 배경격자 h=0.25, L=1+θ,
마지막 요소 활성비율 η=θ/h 스위프. 4x4(좌단 소거) SPD 행렬의 κ 추적.
"""
import numpy as np

H = 0.25
XBG = np.arange(0, 1.2501, H)  # 배경노드 0..1.25 (6개), 마지막 요소가 절단됨


def assemble(L):
    """비적합 조립: 배경요소 형상함수(노드 1.0/1.25)를 물리부분 [1.0,L]에서만 적분.
    ke_cut = η·[[1,-1],[-1,1]]/h (η→0에서 소멸 → 고전적 small-cut 병목)."""
    n = len(XBG)
    K = np.zeros((n, n))
    F = np.zeros(n)
    for i in range(n - 2):  # 온전 요소
        K[i:i + 2, i:i + 2] += np.array([[1, -1], [-1, 1]]) / H
        F[i:i + 2] += np.array([H / 2, H / 2])
    eta = (L - XBG[-2]) / H  # 활성 체적비
    K[-2:, -2:] += eta * np.array([[1, -1], [-1, 1]]) / H
    F[-2:] += np.array([eta * H / 2, eta * H / 2])
    return K[1:, 1:], F[1:]  # u(0)=0 소거 → 5x5


def cond(K):
    w = np.linalg.eigvalsh(K)
    return w[-1] / w[0]


def ghost_matrix(gamma):
    """면 x=1.0 도약 γh⟦u'⟧². 국소 순서: 노드1..5 (전역 1..5)."""
    # jump = (u5-u4)/h - (u4-u3)/h ; 국소 idx: u3=2,u4=3,u5=4
    j = np.array([0.0, 0.0, 1 / H, -2 / H, 1 / H])
    return gamma * H * np.outer(j, j)


def solve_agfem(L):
    """외부노드 u5를 내부 선형외삽 u5=2u4-u3으로 소거 (T 행렬)."""
    K5, F5 = assemble(L)
    T = np.zeros((5, 4))
    T[0, 0] = T[1, 1] = T[2, 2] = T[3, 3] = 1.0
    T[4, 2] = -1.0
    T[4, 3] = 2.0
    Kr = T.T @ K5 @ T
    Fr = T.T @ F5
    return Kr, Fr, T


def l2_error(u_full, L, nq=400):
    xs = np.linspace(0, L, nq)
    uh = np.interp(xs, XBG, np.append([0.0], u_full))
    ue = xs * (2 * L - xs) / 2  # 해석해: -u''=1, u(0)=0, u'(L)=0
    return float(np.sqrt(np.trapezoid((uh - ue) ** 2, xs)))


def test_sliver_sweep():
    print(f"{'eta':>10} {'kappa_plain':>12} {'kappa_ghost':>12} {'kappa_agfem':>12}")
    kappas = []
    for eta in [0.5, 1e-1, 1e-3, 1e-6]:
        L = 1.0 + eta * H
        K, F = assemble(L)
        kp = cond(K)
        kg = cond(K + ghost_matrix(0.05))
        Kr, _, _ = solve_agfem(L)
        ka = cond(Kr)
        kappas.append((kp, kg, ka))
        print(f"{eta:>10.1e} {kp:>12.2e} {kg:>12.2e} {ka:>12.2e}")
    assert kappas[-1][0] > 1e4 * kappas[0][0], "미안정화 폭발이 재현되어야 함"
    assert kappas[-1][1] < 10 * kappas[0][1], "Ghost 바운딩 실패"
    assert kappas[-1][2] < 10 * kappas[0][2], "AgFEM 바운딩 실패"


def test_agfem_accuracy():
    L = 1.0 + 1e-3 * H
    K, F = assemble(L)
    u_plain = np.linalg.solve(K, F)
    Kr, Fr, T = solve_agfem(L)
    u_ag = T @ np.linalg.solve(Kr, Fr)
    e_plain = l2_error(u_plain, L)
    e_ag = l2_error(u_ag, L)
    print(f"[accuracy] L2 plain={e_plain:.2e} agfem={e_ag:.2e}")
    assert e_ag < 1e-2, "안정화 후에도 coarse-mesh 오차 수준을 유지해야 함"


if __name__ == "__main__":
    test_sliver_sweep()
    test_agfem_accuracy()
    print("PoC-03 PASS")
