# -*- coding: utf-8 -*-
"""PoC 07: 2차 텐서곱 패치 Newton 역매핑 + hdf5/hvascii 후처리 스케치.

S(u,v) 쌍2차 패치, 제어망에 비아핀 왜곡 포함. 물리점->모수점 Newton 복원,
control_points.hdf5 라운드트립, 최소 hvascii 블록 작성·파싱.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KNOTS = np.array([0., 0., 0., 1., 1., 1.])
P = 2
# 3x3 제어망 (x, y는 비아핀, z는 범프)
CP = np.zeros((3, 3, 3))
for i in range(3):
    for j in range(3):
        CP[i, j] = [i / 2 + 0.08 * np.sin(j * 1.7),
                    j / 2 + 0.08 * np.cos(i * 1.3),
                    0.15 * np.sin(np.pi * i / 2) * np.sin(np.pi * j / 2)]


def b1(t):
    B = np.array([(1 - t) ** 2, 2 * t * (1 - t), t ** 2])
    dB = np.array([2 * t - 2, 2 - 4 * t, 2 * t])
    return B, dB


def surf(u, v):
    Bu, dBu = b1(u)
    Bv, dBv = b1(v)
    S = np.einsum("i,j,ijk->k", Bu, Bv, CP)
    dSdu = np.einsum("i,j,ijk->k", dBu, Bv, CP)
    dSdv = np.einsum("i,j,ijk->k", Bu, dBv, CP)
    return S, dSdu, dSdv


def invert(pxy, u0=0.4, v0=0.6, tol=1e-13):
    uv = np.array([u0, v0])
    for _ in range(20):
        S, dSu, dSv = surf(*uv)
        r = S[:2] - pxy
        if np.linalg.norm(r) < tol:
            break
        J = np.column_stack([dSu[:2], dSv[:2]])
        uv -= np.linalg.solve(J, r)
    S, _, _ = surf(*uv)
    return uv, float(np.linalg.norm(S[:2] - pxy))


def test_inversion():
    rng = np.random.default_rng(3)
    errs = []
    for (u, v) in rng.uniform(0.05, 0.95, size=(12, 2)):
        S, _, _ = surf(u, v)
        uv, e = invert(S[:2])
        errs.append(max(e, abs(uv[0] - u), abs(uv[1] - v)))
    print(f"[inversion] max roundtrip err = {max(errs):.2e}")
    assert max(errs) < 1e-10


def test_hdf5():
    import h5py
    p = os.path.join(HERE, "control_points.hdf5")
    U = np.arange(27, dtype=float).reshape(9, 3) * 1e-4
    with h5py.File(p, "w") as f:
        f.create_dataset("displacement", data=U,
                         compression="gzip", chunks=(9, 3))
    with h5py.File(p, "r") as f:
        V = f["displacement"][:]
    print(f"[hdf5] shape={V.shape} roundtrip={np.abs(V - U).max():.1e}")
    assert V.shape == (9, 3) and np.abs(V - U).max() == 0.0


def test_hvascii():
    p = os.path.join(HERE, "probe.hvascii")
    with open(p, "w", encoding="utf-8") as f:
        f.write("# PoC-07 probe export (hvascii-style)\nNODES 4\n")
        for i, s in enumerate([30.0e6, 28.5e6, 10.0e6, 9.8e6], start=1):
            f.write(f"{i} VON_MISES {s!r}\n")
    rows = [l.split() for l in open(p, encoding="utf-8") if l[0].isdigit()]
    vals = [float(r[2]) for r in rows]
    print(f"[hvascii] rows={len(rows)} max={max(vals):.3e}")
    assert len(rows) == 4 and abs(max(vals) - 30.0e6) < 1.0


if __name__ == "__main__":
    test_inversion()
    test_hdf5()
    test_hvascii()
    print("PoC-07 PASS")
