"""V06 diagnostic (Kratos-free): on our trimmed B-spline space, compare
  (1) L2 best approximation of the exact Lame displacement (weighted lstsq),
  (2) energy projection u_E = K^-1 b_E, b_E = int sigma_ex : eps(v)  (Galerkin
      solution if load / BC / geometry are consistent),
  (3) independent Galerkin solve with the pressure integrated on the EXACT
      inner sphere (Gauss in theta/phi), strong symmetry (u_d = 0 on plane d=0
      control points, = Coreform XSYMM/YSYMM/ZSYMM),
against Kratos's u error (docs/phase2/evidence/v06_drop). nu = 0, E = 1:
sigma = eps, a(u,v) = int eps(u):eps(v).
"""
import sys, json, time
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.reference.box import BoxSpline
from wht_iga.grid.bspline import basis_1d
from wht_iga.solver import trimmed as T
from wht_iga.kratos_backend.geometry import cell_gauss_points
from tests.kratos_spikes.test_p2mesh_v06geo import (sphere_octant_tris, lame_ur, lame_srr, lame_stt,
                                                    A, B, PI)
from tests.kratos_spikes.test_p2mesh_v06 import LO, HI

OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_energy_proj.json"


def basis_rows(box, X):
    """Values and gradients of all nonzero basis functions at points X."""
    p = box.degree
    nu, nv, nw = box.n_cp
    nloc = (p + 1) ** 3
    idx = np.zeros((len(X), nloc), dtype=np.int64)
    val = np.zeros((len(X), nloc))
    grd = np.zeros((len(X), nloc, 3))
    for a, x in enumerate(X):
        s, D = [], []
        for d in range(3):
            span, ders = basis_1d(box.knots[d], p, float(x[d]), 1)
            s.append(span - p); D.append(np.asarray(ders, float))
        m = 0
        for c in range(p + 1):
            for b in range(p + 1):
                for aa in range(p + 1):
                    idx[a, m] = (s[0] + aa) + nu * ((s[1] + b) + nv * (s[2] + c))
                    val[a, m] = D[0][0, aa] * D[1][0, b] * D[2][0, c]
                    grd[a, m] = (D[0][1, aa] * D[1][0, b] * D[2][0, c],
                                 D[0][0, aa] * D[1][1, b] * D[2][0, c],
                                 D[0][0, aa] * D[1][0, b] * D[2][1, c])
                    m += 1
    return idx, val, grd


def errors(uh, X, W):
    r = np.linalg.norm(X, axis=1); rhat = X / r[:, None]
    u_ex = lame_ur(r)[:, None] * rhat
    e_u = np.sqrt(np.sum(W * np.sum((uh - u_ex) ** 2, axis=1)) / np.sum(W * np.sum(u_ex ** 2, axis=1)))
    urh = np.sum(uh * rhat, axis=1); ure = lame_ur(r)
    e_ur = np.sqrt(np.sum(W * (urh - ure) ** 2) / np.sum(W * ure ** 2))
    return float(e_u), float(e_ur)


def run(n, p, n_geo=24, q=None, n_sph=48):
    t0 = time.perf_counter()
    box = BoxSpline(lo=LO, hi=HI, n_cells=(n, n, n), degree=p)
    q = min(2 * p, 4) if q is None else q
    nodes, tris, _ = sphere_octant_tris(A, B, n_geo)
    qp_table, out_cells, *_ = T.build_trimmed_table(nodes, tris, box, q=q)
    out_set = {tuple(c) for c in out_cells}
    Xs, Ws = [], []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if (i, j, k) in out_set:
                    continue
                P = np.asarray(qp_table[(i, j, k)] if (i, j, k) in qp_table
                               else cell_gauss_points(box, i, j, k), dtype=float)
                Xs.append(P[:, :3]); Ws.append(P[:, 3])
    X = np.vstack(Xs); W = np.concatenate(Ws)
    idx, val, grd = basis_rows(box, X)
    nq, nl = val.shape
    nc = box.n_ctrl
    # scalar mass / Laplacian-type blocks via sparse products
    rows = np.repeat(np.arange(nq), nl)
    Nm = sp.csr_matrix((val.ravel(), (rows, idx.ravel())), shape=(nq, nc))
    G = [sp.csr_matrix((grd[:, :, d].ravel(), (rows, idx.ravel())), shape=(nq, nc)) for d in range(3)]
    Wd = sp.diags(W)
    GG = [[(G[a].T @ Wd @ G[b]).tocsr() for b in range(3)] for a in range(3)]
    lap = GG[0][0] + GG[1][1] + GG[2][2]
    # K[(i,d),(j,e)] = 0.5 (delta_de grad_i.grad_j + g_i,e g_j,d); DOF ordering d-major blocks
    blocks = [[0.5 * ((lap if d == e else 0) + GG[e][d]) for e in range(3)] for d in range(3)]
    K = sp.bmat(blocks).tocsr()
    # exact stress at QPs
    r = np.linalg.norm(X, axis=1); rh = X / r[:, None]
    srr, stt = lame_srr(r), lame_stt(r)
    sig = stt[:, None, None] * np.eye(3)[None] + (srr - stt)[:, None, None] * np.einsum("ni,nj->nij", rh, rh)
    bE = np.concatenate([sum(G[k].T @ (W * sig[:, d, k]) for k in range(3)) for d in range(3)])
    # pressure on the exact inner sphere r = A, traction t = PI * rhat (sigma_rr(A) = -PI)
    gt, gw = np.polynomial.legendre.leggauss(n_sph)
    th = 0.25 * np.pi * (gt + 1.0); wt = 0.25 * np.pi * gw
    TH, PH = np.meshgrid(th, th, indexing="ij"); WW = np.outer(wt, wt)
    rs = np.column_stack([np.sin(TH).ravel() * np.cos(PH).ravel(), np.sin(TH).ravel() * np.sin(PH).ravel(),
                          np.cos(TH).ravel()])
    dA = (A ** 2 * np.sin(TH) * WW).ravel()
    Xs_ = A * rs
    sidx, sval, _ = basis_rows(box, Xs_)
    Ns = sp.csr_matrix((sval.ravel(), (np.repeat(np.arange(len(Xs_)), sval.shape[1]), sidx.ravel())),
                       shape=(len(Xs_), nc))
    bT = np.concatenate([Ns.T @ (dA * PI * rs[:, d]) for d in range(3)])
    area_rel = (dA.sum() - np.pi * A ** 2 / 2) / (np.pi * A ** 2 / 2)
    # active dofs and strong symmetry
    act = np.flatnonzero(lap.diagonal() > 0)
    free = np.zeros(3 * nc, bool)
    for d in range(3):
        free[d * nc + act] = True
        free[d * nc + box.boundary_control_points(d, 0)] = False
    fi = np.flatnonzero(free)
    Kf = K[fi][:, fi].tocsc()
    lu = spla.splu(Kf)
    out = {"n": n, "p": p, "q": q, "n_geo": n_geo, "n_free": int(fi.size), "n_qp": int(nq),
           "sphere_area_rel": float(area_rel)}
    for tag, b in (("energy_proj", bE), ("galerkin_exact_load", bT)):
        c = np.zeros(3 * nc); c[fi] = lu.solve(b[fi])
        uh = np.column_stack([Nm @ c[d * nc:(d + 1) * nc] for d in range(3)])
        out[tag] = errors(uh, X, W)
    # weighted least-squares best approximation (all active functions, no BC)
    u_ex = lame_ur(r)[:, None] * rh
    sw = np.sqrt(W)
    Na = (sp.diags(sw) @ Nm[:, act]).toarray()
    coef, *_ = np.linalg.lstsq(Na, sw[:, None] * u_ex, rcond=1e-12)
    c = np.zeros((nc, 3)); c[act] = coef
    out["l2_best"] = errors(Nm @ c, X, W)
    out["t"] = time.perf_counter() - t0
    return out


if __name__ == "__main__":
    res = []
    for n, p in [(10, 1), (20, 1), (10, 2)]:
        r = run(n, p)
        print(r, flush=True)
        res.append(r)
        json.dump(res, open(OUT, "w"), indent=1)
    print("DONE", flush=True)
