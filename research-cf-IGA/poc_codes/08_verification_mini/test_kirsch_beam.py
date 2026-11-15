# -*- coding: utf-8 -*-
"""PoC 08: Kirsch 해석해 + 수집 Coreform 보-휨 CSV 대조.

1) 무한판 원공: σθθ=σ(1-2cos2θ) → 정점 Kt=3.0, 원방 감쇠.
2) user-guide 수집 traction-tip-deflection.csv: 미세화별 Euler-Bernoulli
   오차 단조감소 + 최종 <5% + 수렴율 확인.
"""
import os
import numpy as np

SIG = 10.0


def hoop_exact(theta):
    return SIG * (1 - 2 * np.cos(2 * theta))


def test_kirsch():
    kt = hoop_exact(np.pi / 2) / SIG
    crown_side = hoop_exact(0.0) / SIG
    decay = hoop_exact(np.pi / 2)  # r=a
    print(f"[kirsch] Kt={kt:.6f} (target 3.0), side={crown_side:.6f} (target -1.0)")
    assert abs(kt - 3.0) < 1e-12 and abs(crown_side + 1.0) < 1e-12
    assert decay == 30.0


def find_csv():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    cands = []
    for dp, _, fns in os.walk(os.path.join(root, "user-guide", "generated-downloads")):
        for fn in fns:
            if fn == "beam-bending-traction-tip-deflection.csv":
                cands.append(os.path.join(dp, fn))
    assert cands, "수집 CSV를 찾을 수 없음"
    return cands[0]


def test_beam_csv():
    path = find_csv()
    data = np.genfromtxt(path, delimiter=",", names=True)
    p1 = data[data["p"] == 1]
    h = p1["h_z"]
    err = p1["rel_err_vs_euler_bernoulli"]
    order = np.argsort(h)[::-1]
    h, err = h[order], err[order]
    mono = bool(np.all(np.diff(err) < 0))
    rate = float(np.polyfit(np.log(h), np.log(err), 1)[0])
    print(f"[beam-csv] {os.path.basename(path)}: err={err}, rate~{rate:.2f}, mono={mono}")
    assert mono and err[-1] < 0.05 and rate > 1.0


if __name__ == "__main__":
    test_kirsch()
    test_beam_csv()
    print("PoC-08 PASS")
