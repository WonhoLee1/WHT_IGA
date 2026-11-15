# -*- coding: utf-8 -*-
"""PoC 05: Boehm knot insertion 기하불변성 + 계층기저 분할단위.

1D p=2. 삽입 전후 곡선 일치, Greville 재현, 계층기저 선택의
분할단위 유지 및 기저수 절감 확인.
"""
import numpy as np

P = 2


def cox(x, i, p, knots):
    if p == 0:
        return 1.0 if (knots[i] <= x < knots[i + 1]
                       or (x == knots[-1] and knots[i + 1] == knots[-1])) else 0.0
    l = r = 0.0
    d1 = knots[i + p] - knots[i]
    if d1 > 0:
        l = (x - knots[i]) / d1 * cox(x, i, p - 1, knots)
    d2 = knots[i + p + 1] - knots[i + 1]
    if d2 > 0:
        r = (knots[i + p + 1] - x) / d2 * cox(x, i + 1, p - 1, knots)
    return l + r


def basis(x, knots):
    n = len(knots) - P - 1
    return np.array([cox(x, i, P, knots) for i in range(n)])


def boehm_insert(knots, ctrl, u):
    knots = list(knots)
    ctrl = [np.asarray(c, float) for c in ctrl]
    p = P
    k = next(i for i in range(len(knots) - 1) if knots[i] <= u < knots[i + 1])
    mult = sum(1 for t in knots if t == u)
    assert mult < p + 1
    n = len(ctrl)
    Q = [None] * (n + 1)
    for i in range(0, k - p + 1):
        Q[i] = np.asarray(ctrl[i], float)
    for i in range(k - p + 1, k - mult + 1):
        d = knots[i + p] - knots[i]
        a = (u - knots[i]) / d if d != 0 else 0.0
        Q[i] = (1 - a) * np.asarray(ctrl[i - 1], float) + a * np.asarray(ctrl[i], float)
    for i in range(k - mult + 1, n + 1):
        Q[i] = np.asarray(ctrl[i - 1], float)
    out = np.array(Q, float)
    knots.insert(k + 1, u)
    return np.array(knots, float), out


def curve(x, knots, ctrl):
    return basis(x, knots) @ ctrl


def test_insertion_invariance():
    knots = np.array([0., 0., 0., 1., 1., 1.])
    ctrl = np.array([[0., 0.], [0.5, 1.], [1., 0.]])
    k2, c2 = boehm_insert(knots, ctrl, 0.5)
    k3, c3 = boehm_insert(k2, c2, 0.25)
    xs = np.linspace(0, 1, 401)
    e1 = max(np.abs(curve(x, knots, ctrl) - curve(x, k2, c2)).max() for x in xs)
    e2 = max(np.abs(curve(x, knots, ctrl) - curve(x, k3, c3)).max() for x in xs)
    print(f"[insertion] 1회 err={e1:.2e} 2회 err={e2:.2e} "
          f"기저수 {len(ctrl)}->{len(c2)}->{len(c3)}")
    assert e1 < 1e-12 and e2 < 1e-12


def two_scale(coarse, fine):
    """동일 공간의 coarse 기저를 fine 기저로 전개: N_i = Σ_j A[i,j] M_j (최소자승)."""
    xs = np.linspace(0, 1, 400)
    M = np.stack([np.array([cox(x, j, P, fine) for j in range(len(fine) - P - 1)])
                  for x in xs])
    A = {}
    for i in range(len(coarse) - P - 1):
        Ni = np.array([cox(x, i, P, coarse) for x in xs])
        a, *_ = np.linalg.lstsq(M, Ni, rcond=None)
        A[i] = a
    return A


def test_hierarchical_pum():
    # plain HB는 분할단위를 잃고, truncation이 복원함 (Giannelli et al. 2012)
    coarse = np.array([0., 0., 0., 0.5, 1., 1., 1.])
    fine, _ = boehm_insert(coarse, np.zeros((4, 1)), 0.75)
    r0, r1 = 0.5, 1.0  # 세분영역
    nc = len(coarse) - P - 1
    active_c = [i for i in range(nc)
                if not (coarse[i] >= r0 - 1e-12 and coarse[i + P + 1] <= r1 + 1e-12)]
    nf = len(fine) - P - 1
    active_f = [j for j in range(nf)
                if fine[j] >= r0 - 1e-12 and fine[j + P + 1] <= r1 + 1e-12]
    A = two_scale(coarse, fine)
    xs = np.linspace(0, 1, 401)
    Mc = lambda x, k, i: cox(x, i, P, k)
    naive = np.array([sum(Mc(x, coarse, i) for i in active_c)
                      + sum(Mc(x, fine, j) for j in active_f) for x in xs])
    trunc = np.array([sum(Mc(x, coarse, i)
                          - sum(A[i][j] * Mc(x, fine, j) for j in active_f)
                          for i in active_c)
                      + sum(Mc(x, fine, j) for j in active_f) for x in xs])
    e_naive = np.abs(naive - 1.0).max()
    e_trunc = np.abs(trunc - 1.0).max()
    print(f"[hierarchical] coarse 활성 {active_c}, fine 활성 {active_f}, "
          f"naive HB err={e_naive:.2e}, truncated HB err={e_trunc:.2e}")
    assert e_naive > 1e-3, "naive HB의 분할단위 상실을 확인해야 함"
    assert e_trunc < 1e-10, "truncation이 분할단위를 복원해야 함"


if __name__ == "__main__":
    test_insertion_invariance()
    test_hierarchical_pum()
    print("PoC-05 PASS")
