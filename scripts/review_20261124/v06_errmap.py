import sys
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import v06_energy_proj as EP
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.kratos_backend.geometry import cell_gauss_points
from tests.kratos_spikes.test_p2mesh_v06geo import sphere_octant_tris, lame_ur, lame_srr, lame_stt, A, B
from tests.kratos_spikes.test_p2mesh_v06 import LO, HI
n = 10; h = 0.1
box = BoxSpline(lo=LO, hi=HI, n_cells=(n, n, n), degree=1)
nodes, tris, _ = sphere_octant_tris(A, B, 24)
qp, out, *_ = T.build_trimmed_table(nodes, tris, box, q=2)
out_set = {tuple(c) for c in out}
Xs, Ws, cid = [], [], []
for k in range(n):
    for j in range(n):
        for i in range(n):
            if (i, j, k) in out_set: continue
            P = np.asarray(qp[(i, j, k)] if (i, j, k) in qp else cell_gauss_points(box, i, j, k), float)
            Xs.append(P[:, :3]); Ws.append(P[:, 3]); cid += [(i, j, k)] * len(P)
X = np.vstack(Xs); W = np.concatenate(Ws); cid = np.array(cid)
idx, val, grd = EP.basis_rows(box, X); nq, nl = val.shape; nc = box.n_ctrl; rows = np.repeat(np.arange(nq), nl)
Nm = sp.csr_matrix((val.ravel(), (rows, idx.ravel())), shape=(nq, nc))
G = [sp.csr_matrix((grd[:, :, d].ravel(), (rows, idx.ravel())), shape=(nq, nc)) for d in range(3)]
Wd = sp.diags(W); GG = [[(G[a].T @ Wd @ G[b]).tocsr() for b in range(3)] for a in range(3)]
lap = GG[0][0] + GG[1][1] + GG[2][2]
K = sp.bmat([[0.5 * ((lap if d == e else 0) + GG[e][d]) for e in range(3)] for d in range(3)]).tocsr()
r = np.linalg.norm(X, axis=1); rh = X / r[:, None]
srr, stt = lame_srr(r), lame_stt(r)
sig = stt[:, None, None] * np.eye(3)[None] + (srr - stt)[:, None, None] * np.einsum("ni,nj->nij", rh, rh)
bE = np.concatenate([sum(G[k].T @ (W * sig[:, d, k]) for k in range(3)) for d in range(3)])
act = np.flatnonzero(lap.diagonal() > 0); free = np.zeros(3 * nc, bool)
for d in range(3):
    free[d * nc + act] = True; free[d * nc + box.boundary_control_points(d, 0)] = False
fi = np.flatnonzero(free); c = np.zeros(3 * nc); c[fi] = spla.spsolve(K[fi][:, fi].tocsc(), bE[fi])
uh = np.column_stack([Nm @ c[d * nc:(d + 1) * nc] for d in range(3)])
uex = lame_ur(r)[:, None] * rh
e2 = W * np.sum((uh - uex) ** 2, axis=1); den = np.sum(W * np.sum(uex ** 2, axis=1))
print("total rel L2", np.sqrt(e2.sum() / den))
# per cell
cells = {}
for a in range(nq):
    t = tuple(cid[a]); cells.setdefault(t, [0.0, 0.0]); cells[t][0] += e2[a]; cells[t][1] += W[a]
fr = np.array([v[1] / h**3 for v in cells.values()]); ce = np.array([v[0] for v in cells.values()])
order = np.argsort(-ce)
keys = list(cells.keys())
print("share of error in top 10 / 30 cells:", ce[order[:10]].sum() / ce.sum(), ce[order[:30]].sum() / ce.sum())
for o in order[:12]:
    i, j, k = keys[o]; cen = (np.array(keys[o]) + 0.5) * h
    print(f"  cell {keys[o]} fill={fr[o]:.3f} r_c={np.linalg.norm(cen):.3f} share={ce[o]/ce.sum():.3f}")
for lo_, hi_ in [(0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 1.01)]:
    m = (fr >= lo_) & (fr < hi_)
    print(f"fill [{lo_},{hi_}): cells {m.sum()}, error share {ce[m].sum()/ce.sum():.3f}")
# error by component: radial vs tangential
du = uh - uex; dr = np.sum(du * rh, axis=1); dt = du - dr[:, None] * rh
print("radial share", np.sum(W * dr**2) / e2.sum(), "tangential share", np.sum(W * np.sum(dt**2, axis=1)) / e2.sum())
mind = X.min(axis=1)
for a_, b_ in [(0, 0.1), (0.1, 0.3), (0.3, 1)]:
    m = (mind >= a_) & (mind < b_)
    print(f"min-coordinate [{a_},{b_}) error share {e2[m].sum()/e2.sum():.3f}, vol share {W[m].sum()/W.sum():.3f}")
print("--- pole cells: ur_h/ur_ex at QPs")
for cell in [(9, 0, 0), (8, 0, 0), (9, 1, 0), (5, 5, 5), (6, 6, 3)]:
    m = np.all(cid == np.array(cell), axis=1)
    if not m.any():
        print(cell, "inactive"); continue
    ratio = np.sum(uh[m] * rh[m], axis=1) / lame_ur(r[m])
    print(cell, f"n_qp={m.sum()} ur ratio min {ratio.min():.4f} max {ratio.max():.4f} mean {ratio.mean():.4f} r range {r[m].min():.3f}-{r[m].max():.3f}")
# control point values near the x-pole: u_x at ctrl points (i,0,0)
nu = box.n_cp[0]
for i in range(7, 11):
    I = i
    print(f" ctrl ({i},0,0) x={i*h:.1f}: u_x = {c[I]:.4e}, exact ur(x) = {lame_ur(np.array([max(i*h,1e-9)]))[0]:.4e}")
