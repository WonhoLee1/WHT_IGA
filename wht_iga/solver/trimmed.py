# -*- coding: utf-8 -*-
"""Trimmed-cell bridge (integration-owned): GEO S3/S4 → Kratos per-cell QP tables.

Uses ``wht_iga.geometry`` (Kratos-free) to classify the background grid and fit
cut-cell rules, then feeds them to ``solver.pipeline.run_linear`` via its
``qp_table`` hook (CUT cells only; IN cells keep full Gauss, OUT cells use
``spec.out_cells`` + FCM ghost-alpha). Created during Phase-1 integration;
ownership: integrator (not Track GEO/SOLVE).
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from wht_iga.grid.bspline import ders_basis_funs as _bs_ders, find_span as _bs_span


_W: dict = {}  # per-worker geometry (set by _init_worker)

#: Below this many non-OUT cells the serial loop wins (process spawn + numba
#: cache load cost ~ seconds) -- unless the per-cell work is large:
_PARALLEL_MIN_CELLS = 400
#: ... cells x surface triangles above this go parallel anyway (V06 level 0 with
#: a 40k-triangle sphere: 355 cells ran serially for 348 s).
_PARALLEL_MIN_WORK = 2_000_000


def _init_worker(nodes: np.ndarray, tris: np.ndarray, h: np.ndarray, q: int, p: int) -> None:
    import os

    # One thread per worker process: the pool supplies the parallelism; nested
    # numba/MKL/OpenMP threads would oversubscribe the cores.
    for var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS"):
        os.environ[var] = "1"
    from wht_iga.geometry import accel as A
    from wht_iga.geometry import quadrature as Q  # noqa: F401  (import after env)

    _W.update(nodes=nodes, tris=tris, h=h, q=q, p=p,
              accel=A.TriSoupAccel(nodes, tris))


def _rule_chunk(blos: np.ndarray) -> list[tuple]:
    from wht_iga.geometry import quadrature as Q

    ac = _W.get("accel")
    out = []
    for blo in blos:
        r = Q.cell_quadrature(_W["nodes"], _W["tris"], blo, _W["h"],
                              q=_W["q"], p=_W["p"], accel=ac)
        out.append((r.x, r.w, bool(r.fallback), float(r.residual), float(r.volume),
                    bool(r.relaxed)))
    return out


def _default_workers() -> int:
    import os

    env = os.environ.get("WHT_GEO_WORKERS", "")
    if env:
        return max(1, int(env))
    return max(1, (os.cpu_count() or 2) // 2)  # physical cores (SMT machines)


def _cell_rules(nodes, tris, blos: np.ndarray, h, q: int, p: int, workers: int,
                accel=None) -> list[tuple]:
    """cell_quadrature for every row of ``blos`` (same order), optionally in a process pool."""
    from wht_iga.geometry import accel as A
    from wht_iga.progress import Ticker, log

    if len(blos) == 0:
        return []
    small = len(blos) < _PARALLEL_MIN_CELLS and len(blos) * len(tris) < _PARALLEL_MIN_WORK
    if workers <= 1 or small:
        # serial: same worker function in-process (no thread-env changes here)
        log(f"cell rules: {len(blos)} active cells, serial")
        ac = A.get_accel(nodes, tris) if accel is None else accel
        _W.update(nodes=nodes, tris=tris, h=h, q=q, p=p, accel=ac)
        return _rule_chunk(blos)
    from concurrent.futures import ProcessPoolExecutor

    n_chunks = workers * 8
    chunks = [c for c in np.array_split(blos, n_chunks) if len(c)]
    log(f"cell rules: {len(blos)} active cells, {workers} processes, {len(chunks)} chunks")
    tick = Ticker("cell rules", len(blos))
    parts = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(nodes, tris, h, q, p)) as ex:
        for part in ex.map(_rule_chunk, chunks):  # in order; progress per chunk
            parts.append(part)
            tick.advance(len(part))
    return [r for part in parts for r in part]


def _cache_key(nodes: np.ndarray, tris: np.ndarray, box, q: int) -> str:
    """Content hash of every input + the source of the code that computes it."""
    import hashlib
    import os

    hsh = hashlib.sha256()
    hsh.update(np.ascontiguousarray(nodes, dtype=np.float64).tobytes())
    hsh.update(np.ascontiguousarray(tris, dtype=np.int64).tobytes())
    hsh.update(repr((tuple(float(v) for v in box.lo), tuple(float(v) for v in box.hi),
                     tuple(int(v) for v in box.n_cells), int(box.degree), int(q))).encode())
    import inspect

    here = os.path.dirname(os.path.abspath(__file__))
    geo = os.path.join(here, "..", "geometry")
    for f in (os.path.join(geo, "accel.py"),
              os.path.join(geo, "classify.py"),
              os.path.join(geo, "quadrature.py")):
        with open(f, "rb") as fh:
            hsh.update(fh.read())
    # only the result-defining functions of this module (scheduling, logging and
    # post-processing edits elsewhere in trimmed.py must not invalidate caches)
    for fn in (_build_trimmed_table, _rule_chunk):
        hsh.update(inspect.getsource(fn).encode())
    return hsh.hexdigest()


def _cache_dir(cache_dir: str | None) -> str | None:
    """None/"" -> env WHT_GEO_CACHE; "0"/"off" disables; default WHT_IGA/.cache/trimmed."""
    import os

    val = cache_dir if cache_dir is not None else os.environ.get("WHT_GEO_CACHE", "")
    if val in ("0", "off", "false"):
        return None
    if val:
        return val
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "..", ".cache", "trimmed"))


def build_trimmed_table(tess_nodes: np.ndarray, tess_tris: np.ndarray,
                        box, q: int, workers: int | None = None,
                        cache_dir: str | None = None):
    """Cached front of :func:`_build_trimmed_table` (same return value).

    The key hashes the tessellation, grid (lo/hi/n_cells/degree), q AND the
    source of classify.py / quadrature.py / trimmed.py, so any code change
    invalidates old entries. Sweeps that only change BC / extension / alpha /
    solver reuse the table. ``cache_dir``: path, or "0"/"off" to disable;
    default env ``WHT_GEO_CACHE`` or ``WHT_IGA/.cache/trimmed``.
    """
    import os
    import pickle

    nodes = np.asarray(tess_nodes, float)
    tris = np.asarray(tess_tris, int)
    cdir = _cache_dir(cache_dir)
    path = None
    if cdir is not None:
        path = os.path.join(cdir, _cache_key(nodes, tris, box, q) + ".pkl")
        if os.path.exists(path):
            from wht_iga.progress import log

            log(f"trimmed table: cache hit {os.path.basename(path)[:12]}")
            with open(path, "rb") as fh:  # local, self-written cache only
                return pickle.load(fh)
    out = _build_trimmed_table(nodes, tris, box, q, workers=workers)
    if path is not None:
        os.makedirs(cdir, exist_ok=True)
        tmp = path + f".{os.getpid()}.tmp"
        with open(tmp, "wb") as fh:
            pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)  # atomic: concurrent runs never read a partial file
    return out


def _build_trimmed_table(tess_nodes: np.ndarray, tess_tris: np.ndarray,
                         box, q: int, workers: int | None = None):
    """Return (qp_table, out_cells, geo_issues, cell_records, diag).

    qp_table: {(i,j,k): (N,4) [[x,y,z,w],...]} for CUT cells only.
    out_cells: [(i,j,k),...] with no physical volume.
    cell_records: mesh-DB rows with true eta/n_qp (v2 §4.2).
    workers: process count for the per-cell rule fits (independent cells);
        None = env ``WHT_GEO_WORKERS`` or physical cores. Results and their
        order are identical to the serial loop.
    """
    from wht_iga.geometry import accel as A
    from wht_iga.geometry import classify as C
    from wht_iga.geometry import quadrature as Q

    lo = np.asarray(box.lo, float)
    hi = np.asarray(box.hi, float)
    n = tuple(int(v) for v in box.n_cells)
    h = (hi - lo) / np.asarray(n)
    nodes = np.asarray(tess_nodes, float)
    tris = np.asarray(tess_tris, int)

    ac = A.get_accel(nodes, tris)
    state = C.classify_cells(nodes, tris, lo, h, n, accel=ac)["state"].reshape(n)
    qp_table: dict[tuple[int, int, int], np.ndarray] = {}
    out_cells: list[tuple[int, int, int]] = []
    issues: list[dict] = []
    cell_records: list[dict] = []
    nx, ny, nz = n
    n_cut = n_fallback = n_relaxed = 0

    # Strict interior IN cells (no triangle AABB overlap) have exact full Gauss rule;
    # only CUT cells and near-surface IN cells need cell_quadrature.
    in_ijk = np.argwhere(state == C.IN)
    strict_in = np.zeros(n, dtype=bool)
    if len(in_ijk):
        # argwhere returns (i, j, k) sorted
        in_blos = lo + h * in_ijk
        has_tri = ac.has_boxes_aabb_overlap(in_blos, h)
        free = in_ijk[~has_tri]
        if len(free):
            strict_in[free[:, 0], free[:, 1], free[:, 2]] = True

    active = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)
              if int(state[i, j, k]) != C.OUT and not strict_in[i, j, k]]
    blos = (lo + h * np.asarray(active, dtype=float)) if active else np.zeros((0, 3))
    rules = iter(_cell_rules(nodes, tris, blos, h, q, int(box.degree),
                             _default_workers() if workers is None else int(workers),
                             accel=ac))
    full_qp_count = int((q + 1) ** 3)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                s = int(state[i, j, k])
                if s == C.OUT:
                    out_cells.append((i, j, k))
                    cell_records.append({"i": i, "j": j, "k": k, "status": 0,
                                         "eta": 0.0, "n_qp": 0})
                    continue
                if strict_in[i, j, k]:
                    cell_records.append({"i": i, "j": j, "k": k, "status": 1,
                                         "eta": 1.0, "n_qp": full_qp_count})
                    continue
                rx, rw, rfb, rres, rvol, rrel = next(rules)
                rule = Q.CellQuadrature(rx, rw, np.zeros(0), fallback=rfb,
                                        residual=rres, volume=rvol, relaxed=rrel)
                if rule.relaxed:
                    n_relaxed += 1
                    issues.append({"severity": "warning", "code": "QUAD_RELAXED",
                                   "message": f"cell {(i, j, k)} moment-fit residual "
                                              f"{rule.residual / float(np.prod(h)):.1e}·|cell| "
                                              f"(> {Q._NNLS_TOL:g}, <= {Q._NNLS_RELAXED:g})",
                                   "entity": f"cell:{i},{j},{k}"})
                if rule.fallback:
                    n_fallback += 1
                    issues.append({"severity": "warning", "code": "QUAD_FALLBACK",
                                   "message": f"cell {(i, j, k)} residual {rule.residual:.2e}",
                                   "entity": f"cell:{i},{j},{k}"})
                if rule.x.shape[0] == 0:
                    # Touch-only degenerate cell: no volume, ghost-alpha is exact
                    # for all practical purposes; record and move on.
                    out_cells.append((i, j, k))
                    issues.append({"severity": "info", "code": "EMPTY_RULE",
                                   "message": f"cell {(i, j, k)} has no quadrature points",
                                   "entity": f"cell:{i},{j},{k}"})
                    cell_records.append({"i": i, "j": j, "k": k, "status": 0,
                                         "eta": 0.0, "n_qp": 0})
                    continue
                if s == C.CUT:
                    n_cut += 1
                    qp_table[(i, j, k)] = np.column_stack(
                        [rule.x, rule.w]).astype(float)
                eta = float(rule.volume / float(np.prod(h))) if float(np.prod(h)) > 0 else 1.0
                cell_records.append({"i": i, "j": j, "k": k, "status": 1,
                                     "eta": min(max(eta, 0.0), 1.0),
                                     "n_qp": int(rule.x.shape[0])})
    diag = {"n_cut": n_cut, "n_fallback": n_fallback, "n_relaxed": n_relaxed,
            "n_qp_cut": int(sum(v.shape[0] for v in qp_table.values()))}
    return qp_table, out_cells, issues, cell_records, diag


def apply_extension_ties(KM, mp, box, cell_records, theta: float = 0.2,
                         start_cid: int = 1) -> dict:
    """AgFEM-style extension as Kratos MPC ties (R1 path, v2 §7.5/§5.2.3).

    cell_records: [{i,j,k,eta,...}] (eta = inside fraction; OUT cells eta 0).
    Builds support classification → extension operator → one
    ``LinearMasterSlaveConstraint`` per slave DOF component (c = 0 enforced).
    Returns keep dict with counts (record in G2 evidence).
    """
    import numpy as np

    from wht_iga.kratos_backend.constraints import add_linear_tie
    from wht_iga.operators import compose as C
    from wht_iga.operators import numbering as N
    from wht_iga.stabilize import extension as E
    from wht_iga.stabilize import support as S

    nx, ny, nz = (int(v) for v in box.n_cells)
    ncell = nx * ny * nz
    h = (np.asarray(box.hi, float) - np.asarray(box.lo, float)) / np.asarray([nx, ny, nz])
    cv = float(np.prod(h))
    eta = np.zeros(ncell)
    if cell_records:
        for c in cell_records:
            eta[int(c["i"]) + nx * (int(c["j"]) + ny * int(c["k"]))] = float(c.get("eta", 1.0))
    else:
        eta[:] = 1.0
    knots = [np.asarray(k, dtype=float) for k in box.knots]
    cls = S.classify_basis(tuple(int(v) for v in box.n_cp), knots,
                           eta * cv, np.full(ncell, cv), theta=float(theta))
    op = E.build_extension(tuple(int(v) for v in box.n_cp), cls.labels,
                           int(box.degree), knots=knots)
    Td = C.expand_dof_operator(C.compose_thb(op.E), n_comp=3).tocsr()
    master_dofs = (np.asarray(op.master_ids).reshape(-1, 1) * 3
                   + np.arange(3)).ravel()
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    dof_of: dict[int, object] = {}
    n_ctrl = int(box.n_ctrl)
    for node in mp.Nodes:
        I = int(node.Id) - 1
        if not 0 <= I < n_ctrl:
            continue  # wrapper/virtual nodes are not spline DOFs
        for d, var in enumerate(comp_vars):
            dof_of[3 * I + d] = node.GetDof(var)
    cid, n_ties = int(start_cid), 0
    for g in [int(v) for v in op.slave_ids]:
        for d in range(3):
            sd = 3 * g + d
            r = Td.getrow(sd)
            mdofs = [int(master_dofs[c]) for c in r.indices]
            N.require_zero_constant(np.zeros(1))
            add_linear_tie(KM, mp, cid, [dof_of[sd]],
                           [dof_of[m] for m in mdofs],
                           np.asarray(r.data).reshape(1, -1), np.zeros(1))
            cid += 1
            n_ties += 1
    return {"n_ties": n_ties, "n_slave": int(len(op.slave_ids)),
            "n_master": int(len(op.master_ids)), "theta": float(theta),
            "keep": (Td, master_dofs)}


# Kratos Voigt order: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz] (F2).
# Single source of truth for Voigt <-> tensor mapping in this module.
_VOIGT = {"xx": 0, "yy": 1, "zz": 2, "xy": 3, "yz": 4, "xz": 5}


def _tensor_from_voigt(v: np.ndarray) -> np.ndarray:
    """Symmetric 3x3 Cauchy tensor from a Kratos-order Voigt vector.

    No 1/2 scaling: Voigt shear entries ARE tensor shear components
    (F1: the old code halved them a second time).
    """
    v = np.asarray(v, dtype=float).ravel()
    if v.shape != (6,):
        raise ValueError(f"Voigt vector must have 6 entries, got {v.shape}")
    o = _VOIGT
    return np.array([[v[o["xx"]], v[o["xy"]], v[o["xz"]]],
                     [v[o["xy"]], v[o["yy"]], v[o["yz"]]],
                     [v[o["xz"]], v[o["yz"]], v[o["zz"]]]])


@njit(cache=True, inline="always")
def _sym3_max_eig_njit(xx: float, yy: float, zz: float,
                       xy: float, yz: float, xz: float) -> float:
    """Largest eigenvalue of a real symmetric 3x3 tensor (exact Smith/Kopp trig)."""
    q = (xx + yy + zz) / 3.0
    b00 = xx - q
    b11 = yy - q
    b22 = zz - q
    p2 = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * (xy * xy + yz * yz + xz * xz)
    if p2 <= 0.0:
        return q
    p = np.sqrt(p2 / 6.0)
    inv_p = 1.0 / p
    c00, c11, c22 = b00 * inv_p, b11 * inv_p, b22 * inv_p
    c01, c12, c02 = xy * inv_p, yz * inv_p, xz * inv_p
    det_b = (c00 * (c11 * c22 - c12 * c12)
             - c01 * (c01 * c22 - c12 * c02)
             + c02 * (c01 * c12 - c11 * c02))
    r = 0.5 * det_b
    if r <= -1.0:
        return q + p
    if r >= 1.0:
        return q + 2.0 * p
    return q + 2.0 * p * np.cos(np.arccos(r) / 3.0)


@njit(cache=True)
def _max_principal_voigt_njit(sig: np.ndarray) -> float:
    """Max over rows of (N, 6) Voigt [xx, yy, zz, xy, yz, xz] of max principal stress."""
    n = sig.shape[0]
    best = -np.inf
    for i in range(n):
        e = _sym3_max_eig_njit(sig[i, 0], sig[i, 1], sig[i, 2],
                               sig[i, 3], sig[i, 4], sig[i, 5])
        if e > best:
            best = e
    return best


def peak_max_principal(cauchy_list: list[np.ndarray]) -> float:
    """Max over QPs of the max principal stress. Kratos Voigt order."""
    best = -np.inf
    for arr in cauchy_list:
        a = np.ascontiguousarray(arr, dtype=np.float64).reshape(-1, 6)
        if a.shape[0]:
            e = float(_max_principal_voigt_njit(a))
            if e > best:
                best = e
    return float(best)


def _d_matrix(E: float, nu: float) -> np.ndarray:
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)) if abs(1.0 - 2.0 * nu) > 1e-12 else 0.0
    mu = E / (2.0 * (1.0 + nu))
    D = np.zeros((6, 6))
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2.0 * mu
    D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = lam
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def hole_face_by_cylinder_wall(tess_nodes: np.ndarray, tess_tris: np.ndarray,
                               face_ids: np.ndarray, radius: float,
                               axis: int = 2, rtol: float = 0.25) -> int:
    """Find the cylindrical-wall face id by geometry, not by heuristic rank.

    Keeps faces whose triangles all have normals perpendicular to ``axis``
    and centroid radius within ``rtol`` of ``radius``. Raises if the match
    is not unique. (F6: explicit face selection.)
    """
    nodes = np.asarray(tess_nodes, float)
    tris = np.asarray(tess_tris, int)
    fids = np.asarray(face_ids)
    e1 = nodes[tris][:, 1] - nodes[tris][:, 0]
    e2 = nodes[tris][:, 2] - nodes[tris][:, 0]
    n = np.cross(e1, e2)
    area = np.linalg.norm(n, axis=1)
    ok = area > 0
    n = n[ok] / area[ok, None]
    cents = nodes[tris][ok].mean(axis=1)
    rr = np.hypot(cents[:, 0], cents[:, 1]) if axis == 2 else np.hypot(
        cents[:, [d for d in range(3) if d != axis][0]],
        cents[:, [d for d in range(3) if d != axis][1]])
    wall = np.abs(n[:, axis]) < 0.1
    tri_face = fids[ok]
    cands = []
    for f in sorted(set(tri_face.tolist())):
        m = tri_face == f
        if m.sum() == 0:
            continue
        if bool(np.all(wall[m])) and abs(float(np.median(rr[m])) - radius) <= rtol * radius:
            cands.append(f)
    if len(cands) != 1:
        raise ValueError(f"cylindrical wall not unique (radius~{radius}): {cands}")
    return int(cands[0])


@njit(cache=True, parallel=True)
def _strain_at_points_njit(pts: np.ndarray, k0: np.ndarray, k1: np.ndarray, k2: np.ndarray,
                           p: int, nu_: int, nv_: int, nw_: int,
                           uh: np.ndarray) -> np.ndarray:
    """Compiled parallel evaluation of engineering strain at physical points."""
    n_pts = pts.shape[0]
    out = np.zeros((n_pts, 6), dtype=np.float64)
    for r in prange(n_pts):
        x0, x1, x2 = pts[r, 0], pts[r, 1], pts[r, 2]
        sp = _bs_span(nu_, p, x0, k0)
        sq = _bs_span(nv_, p, x1, k1)
        sr = _bs_span(nw_, p, x2, k2)
        d0 = _bs_ders(sp, x0, p, 1, k0)
        d1 = _bs_ders(sq, x1, p, 1, k1)
        d2 = _bs_ders(sr, x2, p, 1, k2)
        exx = eyy = ezz = exy = eyz = exz = 0.0
        for c in range(p + 1):
            Nc, Dc = d2[0, c], d2[1, c]
            kc = sr - p + c
            for b in range(p + 1):
                Nb, Db = d1[0, b], d1[1, b]
                base_bc = nu_ * (sq - p + b + nv_ * kc)
                NbNc = Nb * Nc
                DbNc = Db * Nc
                NbDc = Nb * Dc
                for a in range(p + 1):
                    I = (sp - p + a) + base_bc
                    Na, Da = d0[0, a], d0[1, a]
                    ux, uy, uz = uh[3 * I], uh[3 * I + 1], uh[3 * I + 2]
                    dNdx = Da * NbNc
                    dNdy = Na * DbNc
                    dNdz = Na * NbDc
                    exx += dNdx * ux
                    eyy += dNdy * uy
                    ezz += dNdz * uz
                    exy += dNdx * uy + dNdy * ux
                    eyz += dNdy * uz + dNdz * uy
                    exz += dNdx * uz + dNdz * ux
        out[r, 0] = exx
        out[r, 1] = eyy
        out[r, 2] = ezz
        out[r, 3] = exy
        out[r, 4] = eyz
        out[r, 5] = exz
    return out


def _strain_at_points(box, U: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Engineering strain (N,6) in Kratos Voigt order at physical points.

    Kratos-free; sums the tensor-product basis derivatives at each point.
    """
    p = int(box.degree)
    k0 = np.ascontiguousarray(box.knots[0], dtype=np.float64)
    k1 = np.ascontiguousarray(box.knots[1], dtype=np.float64)
    k2 = np.ascontiguousarray(box.knots[2], dtype=np.float64)
    nu_, nv_, nw_ = (int(v) for v in box.n_cp)
    uh = np.ascontiguousarray(U, dtype=np.float64).ravel()
    if uh.shape[0] != int(box.n_ctrl) * 3:
        raise ValueError(f"U length {uh.shape[0]} != 3*n_ctrl")
    pts = np.ascontiguousarray(P, dtype=np.float64).reshape(-1, 3)
    if pts.shape[0] == 0:
        return np.zeros((0, 6), dtype=np.float64)
    return _strain_at_points_njit(pts, k0, k1, k2, p, nu_, nv_, nw_, uh)


def boundary_peak_stress(box, U: np.ndarray, E: float, nu: float,
                         tess_nodes: np.ndarray, tess_tris: np.ndarray,
                         face_ids: np.ndarray, order: int = 6,
                         face_id: int | None = None) -> dict:
    """Peak max-principal stress sampled ON a CAD face (explicit face id).

    Duffy-mapped points on the face triangles + discrete-solution strain
    from spline basis derivatives (grid bspline, arbitrary-point
    evaluation). ``face_id`` is required; use
    :func:`hole_face_by_cylinder_wall` to determine it (F6: no
    in-function heuristic). Returns {peak, n_points, face_id}.
    """
    from wht_iga.geometry.quadrature import _duffy

    if face_id is None:
        raise ValueError("boundary_peak_stress requires an explicit face_id "
                         "(see hole_face_by_cylinder_wall)")
    nodes = np.asarray(tess_nodes, float)
    tris = np.asarray(tess_tris, int)
    fids = np.asarray(face_ids)
    keep = tris[fids == int(face_id)]
    if len(keep) == 0:
        raise ValueError(f"face_id {face_id} matches no triangles")
    DP, _ = _duffy(order)
    pts = []
    for t in nodes[keep]:
        v0, e1, e2 = t[0], t[1] - t[0], t[2] - t[0]
        if float(np.linalg.norm(np.cross(e1, e2))) == 0.0:
            continue
        pts.append(v0[None, :] + e1[None, :] * DP[:, :1] + e2[None, :] * DP[:, 1:])
    P = np.vstack(pts)
    D = _d_matrix(float(E), float(nu))
    S = _strain_at_points(box, U, P)
    sig = S @ D.T
    peak = float(_max_principal_voigt_njit(sig)) if sig.shape[0] else -np.inf
    return {"peak": peak, "n_points": int(P.shape[0]), "face_id": int(face_id)}
