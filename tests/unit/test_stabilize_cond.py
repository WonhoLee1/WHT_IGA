"""T-P2.5 module tests: kappa tool vs analytic oracles + cut-strip study (Kratos-free).

Oracles: 1D Dirichlet Laplacian eigenvalues (exact), dense SVD reference,
reference-assembler stiffness (T-P0.5) for the h-slope leg (untrimmed analog
of the V-COND slope acceptance) and the eta-sweep conditioning study
(unstabilized blow-up vs extension-flat vs FCM floor).
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline
from wht_iga.stabilize import cond as Q


def _laplacian_1d(n: int) -> sp.csr_matrix:
    d = np.full(n, 2.0)
    o = np.full(n - 1, -1.0)
    return sp.diags([o, d, o], [-1, 0, 1], format="csr")


def test_cond_matches_laplacian_spectrum():
    n = 20
    K = _laplacian_1d(n)
    lam = np.array([2 - 2 * np.cos(k * np.pi / (n + 1))
                    for k in range(1, n + 1)])
    ref = lam.max() / lam.min()
    kappa, lmin, lmax = Q.cond_eigsh(K)
    assert kappa == pytest.approx(ref, rel=0.01)
    assert lmin == pytest.approx(lam.min(), rel=0.01)
    assert lmax == pytest.approx(lam.max(), rel=0.01)


def test_cond_eigsh_vs_dense_reference():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(60, 60))
    K = sp.csr_matrix(A @ A.T + 60 * np.eye(60))
    ref = Q.cond_dense(K.toarray())
    kappa, _lmin, _lmax = Q.cond_eigsh(K)
    assert abs(kappa - ref) / ref < 0.05


def test_reduced_matrix_identity_and_selection():
    K = sp.diags([1.0, 2.0, 4.0], format="csr")
    out = Q.reduced_matrix(K, None)
    assert (out != K).nnz == 0
    T = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]]))
    Kr = Q.reduced_matrix(K, T)
    ref = T.toarray().T @ np.diag([1.0, 2.0, 4.0]) @ T.toarray()
    np.testing.assert_allclose(Kr.toarray(), ref, rtol=0, atol=1e-14)


def test_h_slope_math_and_full_dirichlet_minus_two():
    # Tool math on exact h^-2 data.
    s = Q.h_slope(np.array([1.0, 0.5, 0.25]), np.array([1.0, 4.0, 16.0]))
    assert abs(s + 2.0) < 1e-12
    with pytest.raises(ValueError):
        Q.h_slope(np.array([1.0, 0.0]), np.array([1.0, 2.0]))
    # End-to-end on a real system: full-Dirichlet p=1 cube attains the
    # textbook -2 (validates the tool; the V-COND -2±0.3 bar is attainable
    # under full constraint — measured -2.02 here).
    hs, kaps = [], []
    for nc in (4, 6, 8, 10, 12):
        box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (nc, nc, nc), 1)
        K = assemble_stiffness(box, 1.0, 0.0)
        ncp = box.n_cp[0]
        bd = {a + ncp * (b + ncp * c)
              for c in range(ncp) for b in range(ncp) for a in range(ncp)
              if a in (0, ncp - 1) or b in (0, ncp - 1)
              or c in (0, ncp - 1)}
        fixed = np.sort(np.concatenate([np.array(sorted(bd)) * 3 + d
                                        for d in range(3)]))
        free = np.setdiff1d(np.arange(box.n_dofs), fixed)
        kaps.append(Q.cond_eigsh(K[free][:, free])[0])
        hs.append(1.0 / nc)
    slope = Q.h_slope(np.array(hs), np.array(kaps))
    print("full-dirichlet p=1 h:", hs,
          "kappa:", [f"{k:.3e}" for k in kaps], "slope:", slope)
    assert abs(slope + 2.0) <= 0.3


def test_h_slope_single_face_pre_asymptotic_pin():
    # Regression pin (NOT the V-COND bar): single-face-fixed p=2 boxes are
    # pre-asymptotic at nc <= 12 (measured -1.35 raw, -1.67 Jacobi): mixed
    # Dirichlet/Neumann coarse grids show flatter slopes; textbook -2 needs
    # full constraint (above) or finer h. Recorded for the V-COND design
    # note (P25): the -2±0.3 G2 bar must be read at V-COND's own levels.
    hs, kaps = [], []
    for nc in (4, 6, 8, 12):
        box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (nc, nc, nc), 2)
        K = assemble_stiffness(box, 1.0, 0.0)
        ncp = box.n_cp[0]
        fixed_cps = [a + ncp * (b + ncp * 0) for b in range(ncp)
                     for a in range(ncp)]
        fixed = np.sort(np.concatenate([np.array(fixed_cps) * 3 + d
                                        for d in range(3)]))
        free = np.setdiff1d(np.arange(box.n_dofs), fixed)
        kaps.append(Q.cond_eigsh(K[free][:, free])[0])
        hs.append(1.0 / nc)
    slope = Q.h_slope(np.array(hs), np.array(kaps))
    print("single-face p=2 h:", hs, "kappa:", [f"{k:.3e}" for k in kaps],
          "slope:", slope)
    assert abs(slope + 1.35) <= 0.25


def test_eta_sweep_unstabilized_blowup_vs_extension_flat():
    # Synthetic small-cut analog on a 4^3 p=2 box (one interior basis scaled
    # by eta): unstabilized kappa grows ~1/eta; neighbor-average extension
    # rows (unit row sums, exact on constants) stay flat; FCM alpha floor
    # sits between. Documents the three-way pattern V-COND repeats on V02.
    box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (4, 4, 4), 2)
    K = assemble_stiffness(box, 1.0, 0.0).tocsr()
    ncp = box.n_cp[0]
    fixed_cps = [a + ncp * (b + ncp * 0) for b in range(ncp)
                 for a in range(ncp)]
    fixed = np.sort(np.concatenate([np.array(fixed_cps) * 3 + d
                                    for d in range(3)]))
    free = np.setdiff1d(np.arange(box.n_dofs), fixed)
    bad_g = 2 + ncp * (2 + ncp * 2)
    bad = np.array([3 * bad_g + d for d in range(3)])
    assert all(b in set(free.tolist()) for b in bad)
    pos = {d: k for k, d in enumerate(free.tolist())}
    n_f = len(free)
    # extension T: drop bad DOFs, replace by 6-neighbor averages (same comp).
    nbd = []
    for da, db, dc in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                       (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        g = (2 + da) + ncp * ((2 + db) + ncp * (2 + dc))
        nbd.append(g)
    keep = [d for d in free if d not in set(bad.tolist())]
    kpos = {d: k for k, d in enumerate(keep)}
    rows, cols, dat = [], [], []
    for d in keep:
        rows.append(pos[d])
        cols.append(kpos[d])
        dat.append(1.0)
    for dd, b in enumerate(bad):
        for g in nbd:
            rows.append(pos[b])
            cols.append(kpos[3 * g + (dd % 3)])
            dat.append(1.0 / 6.0)
    Text = sp.csr_matrix((dat, (rows, cols)), shape=(n_f, len(keep)))
    etas = [1.0, 1e-2, 1e-4, 1e-6, 1e-8]
    ku, ke, kf = [], [], []
    for eta in etas:
        s = np.ones(box.n_dofs)
        s[bad] = np.sqrt(eta)
        Smat = sp.diags(s, format="csr")
        Ke = (Smat @ K @ Smat).tocsr()
        Kf = Ke[free][:, free]
        ku.append(Q.cond_eigsh(Kf)[0])
        ke.append(Q.cond_eigsh(Q.reduced_matrix(Kf, Text))[0])
        Kfcm = (1 - 1e-8) * Ke + 1e-8 * K
        kf.append(Q.cond_eigsh(Kfcm[free][:, free])[0])
    ku, ke, kf = (np.array(v) for v in (ku, ke, kf))
    print("unstab:", [f"{v:.3e}" for v in ku])
    print("extend:", [f"{v:.3e}" for v in ke])
    print("fcm   :", [f"{v:.3e}" for v in kf])
    assert ku[-1] / ku[0] >= 1e4  # unstabilized blow-up
    assert Q.flatness(ke) <= 10.0  # extension flat (V-COND <= 10x analog)
    assert kf[-1] < ku[-1]  # FCM floor beats unstabilized
