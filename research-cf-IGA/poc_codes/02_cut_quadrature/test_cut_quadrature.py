# -*- coding: utf-8 -*-
"""PoC 02: 절단셀 적분 — 모멘트피팅 vs 높이함수 축소(Saye식 2D).

원형 레벨셋 phi<=0 영역(단위셀 내부)과 단항식 적분을
 (a) nnls 모멘트피팅, (b) x-스트립 GL + y 해석적분으로 계산해
조밀격자 기준값과 비교. numpy/scipy만 사용.
"""
import numpy as np
from scipy.optimize import nnls

CX, CY, R = 0.35, 0.35, 0.45


def inside(x, y):
    return (x - CX) ** 2 + (y - CY) ** 2 <= R ** 2


def reference(n=1600):
    xs = (np.arange(n) + 0.5) / n
    ys = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(xs, ys)
    M = inside(X, Y)
    dA = 1.0 / n ** 2
    out = {"area": M.sum() * dA}
    need = [(a, b) for a in range(6) for b in range(6) if a + b <= 5]
    for (a, b) in need:
        out[(a, b)] = ((X ** a * Y ** b)[M]).sum() * dA
    return out


def monomials():
    return [(a, b) for a in range(5) for b in range(5 - a) if a + b <= 4]


def test_moment_fitting(ref):
    n = 80
    xs = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(xs, ys := xs)
    M = inside(X, Y)
    px, py = X[M], Y[M]
    monos = monomials()
    A = np.stack([px ** a * py ** b for (a, b) in monos], axis=0)  # 10 x npts
    b = np.array([ref[m] for m in monos])
    w, _ = nnls(A, b)
    nnz = int((w > 1e-12 * w.max()).sum())
    errs = {}
    for (a, b_) in [(4, 0), (0, 4), (3, 2), (2, 2)]:
        errs[(a, b_)] = abs(float((w * px ** a * py ** b_).sum()) - ref[(a, b_)])
    print(f"[moment-fitting] pts={len(px)} nnz={nnz} " +
          " ".join(f"M{k}={v:.2e}" for k, v in errs.items()))
    assert nnz <= 40, "too many points for degree-3 fit"
    assert max(errs.values()) < 5e-4
    return nnz


def y_intervals(x):
    d2 = R ** 2 - (x - CX) ** 2
    if d2 < 0:
        return []
    s = np.sqrt(d2)
    lo, hi = max(CY - s, 0.0), min(CY + s, 1.0)
    return [(lo, hi)] if hi > lo else []


def int_y(a, b, lo, hi):
    from math import factorial
    # ∫∫ x^a y^b dy dx: y 원시함수 계수 b+1
    return (hi ** (b + 1) - lo ** (b + 1)) / (b + 1)


def test_height_function(ref, nx=16, nGL=12):
    tg, wg = np.polynomial.legendre.leggauss(nGL)
    edges = np.linspace(0, 1, nx + 1)
    res = {(4, 0): 0.0, (0, 4): 0.0, (3, 2): 0.0, (2, 2): 0.0, "area": 0.0}
    for i in range(nx):
        x0, x1 = edges[i], edges[i + 1]
        for t, w in zip(tg, wg):
            x = 0.5 * ((x1 - x0) * t + x0 + x1)
            jac = 0.5 * (x1 - x0) * w
            for (lo, hi) in y_intervals(x):
                res["area"] += jac * (hi - lo)
                for (a, b) in [(4, 0), (0, 4), (3, 2), (2, 2)]:
                    res[(a, b)] += jac * x ** a * int_y(a, b, lo, hi)
    errs = {k: abs(res[k] - ref[k]) for k in res}
    print(f"[height-function] nx={nx} nGL={nGL} " +
          " ".join(f"{k}={v:.2e}" for k, v in errs.items()))
    # 접선 스트립의 단점 특이성 → 순수 GL은 대수 수렴 (Algoim은 Tanh-Sinh로 해결)
    assert errs["area"] < 1e-4 and max(errs.values()) < 1e-4
    return errs


if __name__ == "__main__":
    ref = reference()
    print(f"[reference] area={ref['area']:.10f}")
    test_moment_fitting(ref)
    test_height_function(ref)
    print("PoC-02 PASS")
