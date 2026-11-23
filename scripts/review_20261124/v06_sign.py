import sys
sys.path.insert(0, r"D:/PythonCodeStudy/CAD-CAE-COPILOT/WHT_IGA")
sys.path.insert(0, r"C:/Users/GOODMAN/AppData/Local/Temp/claude/694f846d5550b343/57a213e2-2165-46f7-af77-ce5cc827bb7c/scratchpad/overnight")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import v06_energy_proj as EP
from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T
from wht_iga.kratos_backend.geometry import cell_gauss_points
from tests.kratos_spikes.test_p2mesh_v06geo import sphere_octant_tris, lame_ur, A, B, PI
from tests.kratos_spikes.test_p2mesh_v06 import LO, HI
for n in (10, 20):
    box = BoxSpline(lo=LO, hi=HI, n_cells=(n, n, n), degree=1)
    nodes, tris, _ = sphere_octant_tris(A, B, 24)
    qp, out, *_ = T.build_trimmed_table(nodes, tris, box, q=2)
    print("n", n, "active cells", n**3 - len(out), "cut", len(qp))
    out_set = {tuple(c) for c in out}
    Xs, Ws = [], []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if (i, j, k) in out_set: continue
                P = np.asarray(qp[(i, j, k)] if (i, j, k) in qp else cell_gauss_points(box, i, j, k), float)
                Xs.append(P[:, :3]); Ws.append(P[:, 3])
    X = np.vstack(Xs); W = np.concatenate(Ws)
    idx, val, grd = EP.basis_rows(box, X); nq, nl = val.shape; nc = box.n_ctrl; rows = np.repeat(np.arange(nq), nl)
    Nm = sp.csr_matrix((val.ravel(), (rows, idx.ravel())), shape=(nq, nc))
    G = [sp.csr_matrix((grd[:, :, d].ravel(), (rows, idx.ravel())), shape=(nq, nc)) for d in range(3)]
    Wd = sp.diags(W); GG = [[(G[a].T @ Wd @ G[b]).tocsr() for b in range(3)] for a in range(3)]
    lap = GG[0][0] + GG[1][1] + GG[2][2]
    K = sp.bmat([[0.5 * ((lap if d == e else 0) + GG[e][d]) for e in range(3)] for d in range(3)]).tocsr()
    r = np.linalg.norm(X, axis=1); rh = X / r[:, None]
    from tests.kratos_spikes.test_p2mesh_v06geo import lame_srr, lame_stt
    srr, stt = lame_srr(r), lame_stt(r)
    sig = stt[:, None, None] * np.eye(3)[None] + (srr - stt)[:, None, None] * np.einsum("ni,nj->nij", rh, rh)
    bE = np.concatenate([sum(G[k].T @ (W * sig[:, d, k]) for k in range(3)) for d in range(3)])
    act = np.flatnonzero(lap.diagonal() > 0); free = np.zeros(3 * nc, bool)
    for d in range(3):
        free[d * nc + act] = True; free[d * nc + box.boundary_control_points(d, 0)] = False
    fi = np.flatnonzero(free); c = np.zeros(3 * nc); c[fi] = spla.spsolve(K[fi][:, fi].tocsc(), bE[fi])
    uh = np.column_stack([Nm @ c[d * nc:(d + 1) * nc] for d in range(3)])
    urh = np.sum(uh * rh, axis=1); ure = lame_ur(r)
    print("  mean ur_h/ur_ex (vol-weighted):", float(np.sum(W * urh) / np.sum(W * ure)),
          " strain energy ratio U_h/U_ex:", float(c[fi] @ bE[fi]) / float(np.sum(W * np.einsum('nij,nij->n', sig, sig))))
