"""V06 sphere-octant tessellation + Lame reference (Kratos-free helpers).

Analytic octant tessellation (inner r=a, outer r=b, three symmetry
planes) and the Lame closed form from
``WHT_IGA/tests/verification/data/v06.json``. Imported by
``test_p2mesh_v06.py`` (driver); the tessellation self-test lives here
so it runs without Kratos.
"""

from __future__ import annotations

import numpy as np
import pytest

A = 0.9
B = 1.0
PI = 1e-6
E_MOD = 1.0
NU = 0.0

THICK = B**3 - A**3
C_LAME = PI * A**3 / THICK
KAPPA = PI * A**3 * B**3 / THICK


def lame_ur(r: np.ndarray) -> np.ndarray:
    """Radial displacement, E=1 (v06.json Lame closed form)."""
    r = np.asarray(r, float)
    return (1.0 - 2.0 * NU) * C_LAME * r + (1.0 + NU) * KAPPA / (2.0 * r**2)


def lame_srr(r: np.ndarray) -> np.ndarray:
    """Radial Cauchy stress."""
    return C_LAME - KAPPA / np.asarray(r, float) ** 3


def lame_stt(r: np.ndarray) -> np.ndarray:
    """Tangential Cauchy stress."""
    return C_LAME + KAPPA / (2.0 * np.asarray(r, float) ** 3)


def _fix_winding(tris: np.ndarray, nodes: np.ndarray, want: np.ndarray) -> np.ndarray:
    """Flip triangles whose normal opposes the desired outward direction."""
    out = tris.copy()
    for k, t in enumerate(out):
        a, b, c = nodes[t]
        n = np.cross(b - a, c - a)
        if float(n @ want[k]) < 0.0:
            out[k] = (t[0], t[2], t[1])
    return out


def sphere_octant_tris(a: float, b: float, n: int,
                       compensate: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Analytic octant tessellation. Returns (nodes, tris, face_ids).

    face_ids: 7 = inner sphere, 8 = outer sphere, 1/2/3 = x/y/z planes.
    All windings point outward from the solid.

    When ``compensate`` is True (default), radial coordinates on the UV sphere
    grid are scaled by ``1 + (delta**2 / 12) * (1 + sin(theta)**2)`` with
    ``delta = pi / (2 * n)`` to cancel the leading-order ``O(n**-2)`` inward
    chord-sagitta bias on triangle-averaged radius (important for ``p=3``).
    """
    th = np.linspace(0.0, np.pi / 2.0, n + 1)  # from +z
    ph = np.linspace(0.0, np.pi / 2.0, n + 1)  # in xy from +x
    delta = np.pi / (2.0 * n)

    def patch(radius: float) -> tuple[np.ndarray, np.ndarray]:
        idx: dict[tuple[int, int], int] = {}
        xyz: list[list[float]] = []
        tri: list[tuple[int, int, int]] = []

        def vid(i: int, j: int) -> int:
            key = (0, 0) if i == 0 else (i, j)  # pole collapse: one node
            if key not in idx:
                st, ct = np.sin(th[i]), np.cos(th[i])
                rr = radius * (1.0 + (delta ** 2 / 12.0) * (1.0 + st ** 2)) if compensate else radius
                idx[key] = len(xyz)
                xyz.append([rr * st * np.cos(ph[j]),
                            rr * st * np.sin(ph[j]),
                            rr * ct])
            return idx[key]

        for i in range(n):
            for j in range(n):
                v00, v10 = vid(i, j), vid(i + 1, j)
                v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
                if i == 0 and v00 == v01:
                    tri.append((v00, v10, v11))  # pole fan: single tri
                else:
                    tri.append((v00, v10, v11))
                    tri.append((v00, v11, v01))
        return np.asarray(xyz, dtype=float), np.asarray(tri, dtype=np.int64)

    in_xyz, in_tri = patch(a)
    out_xyz, out_tri = patch(b)

    def plane(axis: int) -> tuple[np.ndarray, np.ndarray]:
        rr = np.linspace(a, b, n + 1)
        aa = np.linspace(0.0, np.pi / 2.0, n + 1)
        idx: dict[tuple[int, int], int] = {}
        xyz: list[list[float]] = []
        tri: list[tuple[int, int, int]] = []

        def vid(i: int, j: int) -> int:
            if (i, j) not in idx:
                p = [0.0, 0.0, 0.0]
                u = (axis + 1) % 3
                v = (axis + 2) % 3
                if axis == 2:
                    s2t = 1.0
                elif axis == 0:
                    s2t = float(np.cos(aa[j]) ** 2)
                else:
                    s2t = float(np.sin(aa[j]) ** 2)
                scale = (1.0 + (delta ** 2 / 12.0) * (1.0 + s2t)) if compensate else 1.0
                r_val = rr[i] * scale
                p[u] = r_val * np.cos(aa[j])
                p[v] = r_val * np.sin(aa[j])
                idx[(i, j)] = len(xyz)
                xyz.append(p)
            return idx[(i, j)]

        for i in range(n):
            for j in range(n):
                v00, v10 = vid(i, j), vid(i + 1, j)
                v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
                if i == 0 and v00 == v01:
                    tri.append((v00, v10, v11))  # pole fan: single tri
                else:
                    tri.append((v00, v10, v11))
                    tri.append((v00, v11, v01))
        return np.asarray(xyz, dtype=float), np.asarray(tri, dtype=np.int64)

    parts = [(in_xyz, in_tri, 7, -1.0), (out_xyz, out_tri, 8, 1.0)]
    for axis in range(3):
        p_xyz, p_tri = plane(axis)
        parts.append((p_xyz, p_tri, axis + 1, 0.0))
    nodes_list: list[np.ndarray] = []
    tris_list: list[np.ndarray] = []
    fids_list: list[np.ndarray] = []
    for xyz, tri, fid, sign in parts:
        off = sum(len(v) for v in nodes_list)
        ref = np.vstack(nodes_list + [xyz])
        tri = np.asarray(tri)
        want = np.zeros((len(tri), 3))
        if sign != 0.0:
            for k, t in enumerate(tri):
                cen = xyz[list(t)].mean(axis=0)
                want[k] = sign * cen / float(np.linalg.norm(cen))
        else:
            want[:, fid - 1] = -1.0  # symmetry plane: outward is -axis
        tris_list.append(_fix_winding(tri + off, ref, want))
        nodes_list.append(xyz)
        fids_list.append(np.full(len(tri), fid, dtype=np.int64))
    return np.vstack(nodes_list), np.vstack(tris_list), np.concatenate(fids_list)


def tri_area(nodes: np.ndarray, tris: np.ndarray) -> float:
    """Total area of a triangle soup."""
    area = 0.0
    for t in tris:
        a, b, c = nodes[t]
        area += float(np.linalg.norm(np.cross(b - a, c - a))) / 2.0
    return area


def lame_errors(box, U: np.ndarray, stresses: list[np.ndarray],
                qp_table: dict, out_set: set) -> dict[str, float]:
    """Volume-weighted rel-L2 errors vs Lame (CUT QP tables + IN Gauss, OUT skipped).

    ``box`` needs ``.evaluate(coeffs, X) -> (M, 3)`` and ``cell_gauss_points``-
    compatible IN-cell points; kept duck-typed (Kratos-free).
    """
    from wht_iga.kratos_backend.geometry import cell_gauss_points

    coeffs = np.asarray(U, float).reshape(-1, 3)
    num_u = den_u = num_ur = den_ur = 0.0
    num_sr = den_sr = num_st = den_st = 0.0
    ci = 0
    nx, ny, nz = tuple(int(v) for v in box.n_cells)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                svals = np.asarray(stresses[ci], float)
                if (i, j, k) in out_set:
                    ci += 1
                    continue
                if (i, j, k) in qp_table:
                    P = np.asarray(qp_table[(i, j, k)], dtype=float)
                    X, W = P[:, :3], P[:, 3]
                else:
                    full = np.asarray(cell_gauss_points(box, i, j, k), dtype=float)
                    X, W = full[:, :3], full[:, 3]
                uh = box.evaluate(coeffs, X)
                r = np.linalg.norm(X, axis=1)
                rhat = X / r[:, None]
                u_ex = (lame_ur(r)[:, None] * rhat)
                num_u += float(np.sum(W * np.sum((uh - u_ex) ** 2, axis=1)))
                den_u += float(np.sum(W * np.sum(u_ex**2, axis=1)))
                ur_n = np.sum(uh * rhat, axis=1)
                ur_e = lame_ur(r)
                num_ur += float(np.sum(W * (ur_n - ur_e) ** 2))
                den_ur += float(np.sum(W * ur_e**2))
                # Kratos CAUCHY_STRESS_VECTOR order is [xx, yy, zz, xy, yz, xz]
                # (verified against Kratos output: test_s04_voigt_order_kratos).
                # The former unpacking read [xx, yy, zz, yz, xz, xy] -> wrong srr/stt.
                sxx, syy, szz, sxy, syz, sxz = (svals[:, 0], svals[:, 1], svals[:, 2],
                                                svals[:, 3], svals[:, 4], svals[:, 5])
                tr = sxx + syy + szz
                srr_n = (sxx * rhat[:, 0]**2 + syy * rhat[:, 1]**2 + szz * rhat[:, 2]**2
                         + 2.0 * (sxy * rhat[:, 0] * rhat[:, 1]
                                  + sxz * rhat[:, 0] * rhat[:, 2]
                                  + syz * rhat[:, 1] * rhat[:, 2]))
                stt_n = (tr - srr_n) / 2.0
                srr_e, stt_e = lame_srr(r), lame_stt(r)
                num_sr += float(np.sum(W * (srr_n - srr_e) ** 2))
                den_sr += float(np.sum(W * srr_e**2))
                num_st += float(np.sum(W * (stt_n - stt_e) ** 2))
                den_st += float(np.sum(W * stt_e**2))
                ci += 1
    return {"u": float(np.sqrt(num_u / den_u)), "ur": float(np.sqrt(num_ur / den_ur)),
            "srr": float(np.sqrt(num_sr / den_sr)), "stt": float(np.sqrt(num_st / den_st))}


def test_v06_sphere_tessellation() -> None:
    """Analytic octant tessellation: patch areas + outward winding."""
    from wht_iga.kratos_backend.surface import triangle_normal

    nodes, tris, fids = sphere_octant_tris(A, B, 8)
    inner, outer = tris[fids == 7], tris[fids == 8]
    assert abs(tri_area(nodes, inner) - np.pi * A**2 / 2.0) / (np.pi * A**2 / 2.0) < 0.02
    assert abs(tri_area(nodes, outer) - np.pi * B**2 / 2.0) / (np.pi * B**2 / 2.0) < 0.02
    for ax in range(3):
        exact = np.pi * (B**2 - A**2) / 4.0
        got = tri_area(nodes, tris[fids == ax + 1])
        assert abs(got - exact) / exact < 0.02, (ax, got, exact)
    for t, fid in zip(tris, fids):
        n = triangle_normal(nodes[t])
        cen = nodes[t].mean(axis=0)
        if fid == 7:
            assert float(n @ cen) < 0.0  # toward the cavity
        elif fid == 8:
            assert float(n @ cen) > 0.0
        else:
            assert float(n[fid - 1]) < -0.999


def test_v06_sphere_watertight_and_sagitta_compensated() -> None:
    """Watertightness of compensated octant + O(n^-2) bias cancellation check."""
    nodes, tris, fids = sphere_octant_tris(A, B, 12, compensate=True)
    keys = [tuple(np.round(p, 12)) for p in nodes]
    edges: dict[tuple[tuple[float, float, float], tuple[float, float, float]], int] = {}
    for a, b, c in tris:
        for u, v in ((keys[a], keys[b]), (keys[b], keys[c]), (keys[c], keys[a])):
            edges[(u, v)] = edges.get((u, v), 0) + 1
    assert all(c == 1 and edges.get((v, u), 0) == 1 for (u, v), c in edges.items())

    def signed_vol(nd: np.ndarray, tr: np.ndarray) -> float:
        p0, p1, p2 = nd[tr[:, 0]], nd[tr[:, 1]], nd[tr[:, 2]]
        return float(np.sum(np.einsum("ij,ij->i", p0, np.cross(p1, p2))) / 6.0)

    n_c, t_c, f_c = sphere_octant_tris(A, B, 24, compensate=True)
    n_r, t_r, f_r = sphere_octant_tris(A, B, 24, compensate=False)
    v_ex = np.pi * (B**3 - A**3) / 6.0
    a_ex = np.pi * A**2 / 2.0
    assert abs(signed_vol(n_c, t_c) - v_ex) < 0.005 * abs(signed_vol(n_r, t_r) - v_ex)
    assert abs(tri_area(n_c, t_c[f_c == 7]) - a_ex) < 0.5 * abs(tri_area(n_r, t_r[f_r == 7]) - a_ex)


def _exact_voigt(X: np.ndarray) -> np.ndarray:
    """Exact Lame Cartesian stress in Kratos Voigt order [xx,yy,zz,xy,yz,xz] at X."""
    r = np.linalg.norm(X, axis=1)
    rhat = X / r[:, None]
    srr, stt = lame_srr(r), lame_stt(r)
    S = np.zeros((X.shape[0], 6))
    for d in range(3):
        for e in range(d, 3):
            comp = srr * rhat[:, d] * rhat[:, e]
            if d == e:
                comp = comp + stt * (1.0 - rhat[:, d] ** 2)
            else:
                comp = comp + stt * (-rhat[:, d] * rhat[:, e])
            voigt = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3,
                     (1, 2): 4, (0, 2): 5}[(d, e)]  # Kratos order (was xx,yy,zz,yz,xz,xy)
            S[:, voigt] += comp
    return S


class _FakeBox:
    """Duck-typed box stub: evaluate returns factor * exact Lame displacement."""

    def __init__(self, n_cells: tuple[int, int, int], factor: float = 1.0) -> None:
        self.n_cells = n_cells
        self._factor = factor

    def evaluate(self, coeffs: np.ndarray, X: np.ndarray) -> np.ndarray:
        r = np.linalg.norm(X, axis=1)
        return self._factor * (lame_ur(r)[:, None] * (X / r[:, None]))


def test_lame_errors_exact_and_scaled() -> None:
    """Error evaluator: ~0 on the exact field, 1% on a 1%-scaled field."""
    rng = np.random.default_rng(7)
    X = 0.9 + 0.1 * rng.random((40, 3))
    X = X / np.linalg.norm(X, axis=1, keepdims=True) * (0.9 + 0.1 * rng.random((40, 1)))
    W = np.full(40, 0.001)
    box = _FakeBox((1, 1, 1))
    qp_table = {(0, 0, 0): np.column_stack([X, W])}
    S = _exact_voigt(X)
    err = lame_errors(box, np.zeros(3), [S], qp_table, set())
    for qty, val in err.items():
        assert val < 1e-12, (qty, val)
    box2 = _FakeBox((1, 1, 1), factor=1.01)
    err2 = lame_errors(box2, np.zeros(3), [1.01 * S], qp_table, set())
    for qty, val in err2.items():
        assert abs(val - 0.01) < 1e-12, (qty, val)


@pytest.mark.parametrize("cell", [(3, 7, 0), (7, 2, 2), (6, 4, 2), (6, 6, 1)])
def test_cell_moment_volume_matches_sampling_and_rotation(cell) -> None:
    """Regression (2026-11-15): cross-section rect-edge splits used projected
    triangle edges with the wrong parameter -> negative / axis-dependent M000.
    Exact moment volume must match dense winding sampling and be invariant
    under cyclic axis rotations (Kratos-free)."""
    from wht_iga.geometry import quadrature as Q
    from wht_iga.geometry.classify import winding_inside

    nodes, tris, _ = sphere_octant_tris(A, B, 24)
    soup = nodes[tris]
    h = np.full(3, 0.1)
    blo = h * np.array(cell)
    g = (np.arange(60) + 0.5) / 60
    G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    sampled = winding_inside(blo + G * h, nodes, tris).mean()
    vols = []
    for perm in ([0, 1, 2], [1, 2, 0], [2, 0, 1]):
        sp, bp = soup[:, :, perm], blo[perm]
        hit = np.array([Q._tri_box_overlap(t[0], t[1], t[2], bp, bp + h) for t in sp])
        vols.append(Q._cell_moments(sp[hit], bp, h, 2, whole=sp)[0, 0, 0] / float(np.prod(h)))
    np.testing.assert_allclose(vols, vols[0], rtol=0, atol=1e-8)  # coplanar tol 1e-9*h level
    assert abs(vols[0] - sampled) < 3e-3  # 60^3 midpoint sampling error
