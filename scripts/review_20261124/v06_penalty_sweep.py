"""V06 L0: symmetry by penalty (Coreform mesh API: essential BCs via bc_penalty_value)
instead of strong control-point fixing. Kratos-free, same space/quadrature as
v06_energy_proj. Penalty term beta * int_{plane d} u_d v_d dA over the exact
quarter annulus a <= rho <= b on each symmetry plane. Reports u, ur, srr, stt
relative L2 errors (sigma = eps for E = 1, nu = 0)."""
import sys, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import v06_energy_proj as EP
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.kratos_backend.geometry import cell_gauss_points
from tests.kratos_spikes.test_p2mesh_v06geo import sphere_octant_tris, lame_ur, lame_srr, lame_stt, A, B, PI
from tests.kratos_spikes.test_p2mesh_v06 import LO, HI

OUT = r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_penalty_sweep.json"


def setup(n, p, n_geo=24):
    box = BoxSpline(lo=LO, hi=HI, n_cells=(n, n, n), degree=p)
    nodes, tris, _ = sphere_octant_tris(A, B, n_geo)
    qp, out, *_ = T.build_trimmed_table(nodes, tris, box, q=min(2 * p, 4))
    out_set = {tuple(c) for c in out}
    Xs, Ws = [], []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if (i, j, k) in out_set:
                    continue
                P = np.asarray(qp[(i, j, k)] if (i, j, k) in qp else cell_gauss_points(box, i, j, k), float)
                Xs.append(P[:, :3]); Ws.append(P[:, 3])
    X = np.vstack(Xs); W = np.concatenate(Ws)
    idx, val, grd = EP.basis_rows(box, X)
    nq, nl = val.shape; nc = box.n_ctrl; rows = np.repeat(np.arange(nq), nl)
    Nm = sp.csr_matrix((val.ravel(), (rows, idx.ravel())), shape=(nq, nc))
    G = [sp.csr_matrix((grd[:, :, d].ravel(), (rows, idx.ravel())), shape=(nq, nc)) for d in range(3)]
    Wd = sp.diags(W); GG = [[(G[a].T @ Wd @ G[b]).tocsr() for b in range(3)] for a in range(3)]
    lap = GG[0][0] + GG[1][1] + GG[2][2]
    K = sp.bmat([[0.5 * ((lap if d == e else 0) + GG[e][d]) for e in range(3)] for d in range(3)]).tocsr()
    # exact-sphere pressure load
    gt, gw = np.polynomial.legendre.leggauss(48)
    th = 0.25 * np.pi * (gt + 1.0); wt = 0.25 * np.pi * gw
    TH, PH = np.meshgrid(th, th, indexing="ij"); WW = np.outer(wt, wt)
    rs = np.column_stack([np.sin(TH).ravel() * np.cos(PH).ravel(), np.sin(TH).ravel() * np.sin(PH).ravel(), np.cos(TH).ravel()])
    dA = (A ** 2 * np.sin(TH) * WW).ravel()
    si, sv, _ = EP.basis_rows(box, A * rs)
    Ns = sp.csr_matrix((sv.ravel(), (np.repeat(np.arange(len(rs)), sv.shape[1]), si.ravel())), shape=(len(rs), nc))
    f = np.concatenate([Ns.T @ (dA * PI * rs[:, d]) for d in range(3)])
    # penalty mass on each symmetry plane (quarter annulus, polar Gauss)
    gr, grw = np.polynomial.legendre.leggauss(12)
    rho = A + 0.5 * (B - A) * (gr + 1); wr = 0.5 * (B - A) * grw
    gp, gpw = np.polynomial.legendre.leggauss(48)
    phi = 0.25 * np.pi * (gp + 1); wp = 0.25 * np.pi * gpw
    RR, PP = np.meshgrid(rho, phi, indexing="ij"); wA = (np.outer(wr, wp) * RR).ravel()
    s1, s2 = (RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel()
    Mp = []
    for d in range(3):
        Xp = np.zeros((len(s1), 3)); o = [e for e in range(3) if e != d]
        Xp[:, o[0]], Xp[:, o[1]] = s1, s2
        pi_, pv, _ = EP.basis_rows(box, Xp)
        Np = sp.csr_matrix((pv.ravel(), (np.repeat(np.arange(len(Xp)), pv.shape[1]), pi_.ravel())), shape=(len(Xp), nc))
        Mp.append((Np.T @ sp.diags(wA) @ Np).tocsr())
    act = np.flatnonzero(lap.diagonal() > 0)
    return dict(box=box, X=X, W=W, Nm=Nm, G=G, K=K, f=f, Mp=Mp, act=act, nc=nc)


def errors(S, c):
    X, W, Nm, G, nc = S["X"], S["W"], S["Nm"], S["G"], S["nc"]
    r = np.linalg.norm(X, axis=1); rh = X / r[:, None]
    uh = np.column_stack([Nm @ c[d * nc:(d + 1) * nc] for d in range(3)])
    gu = np.stack([np.column_stack([G[k] @ c[d * nc:(d + 1) * nc] for k in range(3)]) for d in range(3)], axis=1)  # (nq, d, k)
    eps = 0.5 * (gu + np.transpose(gu, (0, 2, 1)))
    srr_h = np.einsum("ni,nij,nj->n", rh, eps, rh); stt_h = (np.trace(eps, axis1=1, axis2=2) - srr_h) / 2
    uex = lame_ur(r)[:, None] * rh
    rel = lambda a, b: float(np.sqrt(np.sum(W * a) / np.sum(W * b)))
    return {"u": rel(np.sum((uh - uex) ** 2, 1), np.sum(uex ** 2, 1)),
            "ur": rel((np.sum(uh * rh, 1) - lame_ur(r)) ** 2, lame_ur(r) ** 2),
            "srr": rel((srr_h - lame_srr(r)) ** 2, lame_srr(r) ** 2),
            "stt": rel((stt_h - lame_stt(r)) ** 2, lame_stt(r) ** 2)}


def solve(S, beta=None):
    nc, act, K, f = S["nc"], S["act"], S["K"], S["f"]
    free = np.zeros(3 * nc, bool)
    for d in range(3):
        free[d * nc + act] = True
    Kt = K.copy()
    if beta is None:  # strong
        for d in range(3):
            free[d * nc + S["box"].boundary_control_points(d, 0)] = False
    else:
        Kt = K + sp.block_diag([beta * S["Mp"][d] for d in range(3)]).tocsr()
    fi = np.flatnonzero(free)
    c = np.zeros(3 * nc); c[fi] = spla.spsolve(Kt[fi][:, fi].tocsc(), f[fi])
    return errors(S, c)


if __name__ == "__main__":
    coreform = {"u": 0.00584, "ur": 0.00568, "srr": 0.285, "stt": 0.0189}
    S = setup(10, 1)
    h, E = 0.1, 1.0
    res = [{"beta": "strong", **solve(S)}]
    for bh in (1, 3, 10, 30, 100, 300, 1e3, 1e4, 1e6):
        res.append({"beta": bh * E / h, "beta_h_over_E": bh, **solve(S, bh * E / h)})
    for row in res:
        print({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print("Coreform h0p1:", coreform)
    json.dump({"rows": res, "coreform_h0p1": coreform}, open(OUT, "w"), indent=1)
