# -*- coding: utf-8 -*-
"""D1 verification: energy-norm & L2 projection bounds for V06 L0 p=1 (10x10x10).

Demonstrates mathematically why the unpadded 10x10x10 p=1 pure-Q1 Galerkin
discretization (1 cell across wall thickness b-a=0.1, 1,527 active master DOFs)
hits a variational locking floor ``u_E ~ 1.009e-2`` > Coreform ``0.008`` tolerance,
whereas L2 interpolation reaches ``8.64e-4``.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from tests.kratos_spikes.test_p2mesh_v06geo import (
    A,
    B,
    lame_errors,
    lame_srr,
    lame_stt,
    lame_ur,
    sphere_octant_tris,
)
from wht_iga.analyses.bimaterial import _cell_shape_and_strain, isotropic_D
from wht_iga.analyses.thermal import solve_with_dirichlet
from wht_iga.kratos_backend.geometry import cell_gauss_points
from wht_iga.kratos_backend.mass import reduce_by_operator
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.stabilize.extension import build_extension
from wht_iga.stabilize.support import classify_basis


def compute_v06_l0_p1_projections() -> dict:
    """Compute energy-norm and L2-norm best approximations on V06 L0 p=1."""
    t0 = time.perf_counter()
    nodes, tris, _ = sphere_octant_tris(A, B, 24)
    box = BoxSpline((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (10, 10, 10), degree=1)
    qp_table, out_cells, _issues, _recs, diag = T.build_trimmed_table(nodes, tris, box, q=2)
    out_set = {tuple(c) for c in out_cells}

    nx, ny, nz = box.n_cells
    bx, by, bz = (box.breakpoints(d) for d in range(3))
    hx, hy, hz = np.diff(bx), np.diff(by), np.diff(bz)
    cvol = np.einsum("i,j,k->kji", hx, hy, hz).reshape(-1)
    v_act = np.zeros_like(cvol)

    D = isotropic_D(1.0, 0.0)
    rows, cols, kvals, mvals = [], [], [], []
    b_e = np.zeros(box.n_dofs)
    b_l2 = np.zeros(box.n_dofs)

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = (i, j, k)
                cid = i + nx * (j + ny * k)
                if c in out_set:
                    continue
                pts_w = np.asarray(qp_table[c] if c in qp_table else cell_gauss_points(box, i, j, k), dtype=float)
                v_act[cid] = float(pts_w[:, 3].sum())
                n_u, b_mat, edofs = _cell_shape_and_strain(box, c, pts_w[:, :3])
                x_pts, w_pts = pts_w[:, :3], pts_w[:, 3]
                r = np.linalg.norm(x_pts, axis=1)
                rhat = x_pts / r[:, None]
                u_ex = lame_ur(r)[:, None] * rhat
                srr, stt = lame_srr(r), lame_stt(r)
                d_s = srr - stt
                sig_ex = np.column_stack([
                    stt + d_s * rhat[:, 0] ** 2,
                    stt + d_s * rhat[:, 1] ** 2,
                    stt + d_s * rhat[:, 2] ** 2,
                    d_s * rhat[:, 0] * rhat[:, 1],
                    d_s * rhat[:, 1] * rhat[:, 2],
                    d_s * rhat[:, 0] * rhat[:, 2],
                ])
                k_el = np.einsum("q,qai,ab,qbj->ij", w_pts, b_mat, D, b_mat)
                m_el = np.einsum("q,qca,qcb->ab", w_pts, n_u, n_u)
                b_e[edofs] += np.einsum("q,qsi,qs->i", w_pts, b_mat, sig_ex)
                b_l2[edofs] += np.einsum("q,qci,qc->i", w_pts, n_u, u_ex)
                m_dofs = edofs.size
                rows.append(np.repeat(edofs, m_dofs))
                cols.append(np.tile(edofs, m_dofs))
                kvals.append(k_el.reshape(-1))
                mvals.append(m_el.reshape(-1))

    k_mat = sp.coo_matrix(
        (np.concatenate(kvals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(box.n_dofs, box.n_dofs),
    ).tocsr()
    m_mat = sp.coo_matrix(
        (np.concatenate(mvals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(box.n_dofs, box.n_dofs),
    ).tocsr()

    cls = classify_basis(box.n_cp, list(box.knots), v_act, cvol, theta=0.2)
    ext = build_extension(box.n_cp, cls.labels, box.degree, knots=list(box.knots))
    t_op = sp.kron(ext.E, sp.eye(3, format="csr"), format="csr")
    k_red = reduce_by_operator(k_mat, t_op)
    m_red = reduce_by_operator(m_mat, t_op)
    b_e_red = np.asarray(t_op.T @ b_e, dtype=float)
    b_l2_red = np.asarray(t_op.T @ b_l2, dtype=float)

    mset = {int(b): m_idx for m_idx, b in enumerate(ext.master_ids)}
    fix_red = sorted({
        3 * mset[b] + ax
        for ax in range(3)
        for b in box.boundary_control_points(ax, 0)
        if int(b) in mset
    })
    u_e = np.asarray(
        t_op @ solve_with_dirichlet(k_red, b_e_red, np.array(fix_red, dtype=np.int64), np.zeros(len(fix_red))),
        dtype=float,
    )
    u_l2 = np.asarray(
        t_op @ solve_with_dirichlet(m_red, b_l2_red, np.array(fix_red, dtype=np.int64), np.zeros(len(fix_red))),
        dtype=float,
    )
    u_l2_free = np.asarray(t_op @ spla.spsolve(sp.csc_matrix(m_red), b_l2_red), dtype=float)

    def _eval(u_vec: np.ndarray) -> dict[str, float]:
        stresses = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    c = (i, j, k)
                    if c in out_set:
                        stresses.append(np.zeros((0, 6)))
                        continue
                    pts_w = np.asarray(qp_table[c] if c in qp_table else cell_gauss_points(box, i, j, k), dtype=float)
                    _, b_mat, edofs = _cell_shape_and_strain(box, c, pts_w[:, :3])
                    stresses.append(np.einsum("ab,qbi,i->qa", D, b_mat, u_vec[edofs]))
        return lame_errors(box, u_vec, stresses, qp_table, out_set)

    return {
        "n_cells": [10, 10, 10],
        "degree": 1,
        "q": 2,
        "n_geo": 24,
        "n_master": int(cls.n_master),
        "n_extended": int(cls.n_extended),
        "n_reduced_dofs": int(k_red.shape[0]),
        "n_fallback": int(diag.get("n_fallback", 0)),
        "energy_projection_sym": _eval(u_e),
        "l2_projection_sym": _eval(u_l2),
        "l2_projection_unconstrained": _eval(u_l2_free),
        "runtime_s": float(time.perf_counter() - t0),
    }


def main() -> None:
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "docs", "phase2", "evidence", "v06_nnlsfix", "L0_p1_energy_projection.json",
    )
    res = compute_v06_l0_p1_projections()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
