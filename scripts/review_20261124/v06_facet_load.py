"""V06: effect of the facetted inner-sphere pressure (n_geo) on the p=2 L0 Galerkin error (Kratos-free)."""
import sys, json
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import v06_energy_proj as EP
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.kratos_backend.geometry import cell_gauss_points
from tests.kratos_spikes.test_p2mesh_v06geo import sphere_octant_tris, A, B, PI
from tests.kratos_spikes.test_p2mesh_v06 import LO, HI
# Dunavant degree-5, 7 points (barycentric, weights sum 1)
a1, b1 = 0.059715871789770, 0.470142064105115
a2, b2 = 0.797426985353087, 0.101286507323456
BARY = np.array([[1/3, 1/3, 1/3], [a1, b1, b1], [b1, a1, b1], [b1, b1, a1], [a2, b2, b2], [b2, a2, b2], [b2, b2, a2]])
WB = np.array([0.225] + [0.132394152788506] * 3 + [0.125939180544827] * 3)
def facet_load(box, n_geo_load):
    nodes, tris, fids = sphere_octant_tris(A, B, n_geo_load)
    inner = tris[fids == 7]
    P = nodes[inner]                       # (nt,3,3)
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]); area = 0.5 * np.linalg.norm(nrm, axis=1)
    m = nrm / (2 * area)[:, None]
    cen = P.mean(axis=1); m *= np.sign(np.sum(m * cen, axis=1))[:, None]   # radially outward = into solid at r=a
    X = np.einsum("qk,tkd->tqd", BARY, P).reshape(-1, 3)
    w = (area[:, None] * WB[None]).ravel()
    t = np.repeat(PI * m, len(WB), axis=0)
    idx, val, _ = EP.basis_rows(box, X)
    Ns = sp.csr_matrix((val.ravel(), (np.repeat(np.arange(len(X)), val.shape[1]), idx.ravel())), shape=(len(X), box.n_ctrl))
    res = np.array([np.sum(w * t[:, d]) for d in range(3)])
    return np.concatenate([Ns.T @ (w * t[:, d]) for d in range(3)]), res
n, p = 10, 2
box = BoxSpline(lo=LO, hi=HI, n_cells=(n, n, n), degree=p)
nodes, tris, _ = sphere_octant_tris(A, B, 24)
qp_table, out_cells, *_ = T.build_trimmed_table(nodes, tris, box, q=4)
out_set = {tuple(c) for c in out_cells}
Xs, Ws = [], []
for k in range(n):
    for j in range(n):
        for i in range(n):
            if (i, j, k) in out_set: continue
            Pq = np.asarray(qp_table[(i, j, k)] if (i, j, k) in qp_table else cell_gauss_points(box, i, j, k), float)
            Xs.append(Pq[:, :3]); Ws.append(Pq[:, 3])
X = np.vstack(Xs); W = np.concatenate(Ws)
idx, val, grd = EP.basis_rows(box, X)
nq, nl = val.shape; nc = box.n_ctrl; rows = np.repeat(np.arange(nq), nl)
Nm = sp.csr_matrix((val.ravel(), (rows, idx.ravel())), shape=(nq, nc))
G = [sp.csr_matrix((grd[:, :, d].ravel(), (rows, idx.ravel())), shape=(nq, nc)) for d in range(3)]
Wd = sp.diags(W); GG = [[(G[a].T @ Wd @ G[b]).tocsr() for b in range(3)] for a in range(3)]
lap = GG[0][0] + GG[1][1] + GG[2][2]
K = sp.bmat([[0.5 * ((lap if d == e else 0) + GG[e][d]) for e in range(3)] for d in range(3)]).tocsr()
act = np.flatnonzero(lap.diagonal() > 0); free = np.zeros(3 * nc, bool)
for d in range(3):
    free[d * nc + act] = True; free[d * nc + box.boundary_control_points(d, 0)] = False
fi = np.flatnonzero(free); lu = spla.splu(K[fi][:, fi].tocsc())
exact_res = PI * np.pi * A ** 2 / 4
out = []
for ng in (24, 48, 96, 192):
    b, res = facet_load(box, ng)
    c = np.zeros(3 * nc); c[fi] = lu.solve(b[fi])
    uh = np.column_stack([Nm @ c[d * nc:(d + 1) * nc] for d in range(3)])
    e = EP.errors(uh, X, W)
    row = {"n_geo_load": ng, "resultant_rel": float((res[0] - exact_res) / exact_res), "u": e[0], "ur": e[1]}
    print(row, flush=True); out.append(row)
json.dump(out, open(r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight/v06_facet_load.json", "w"), indent=1)
