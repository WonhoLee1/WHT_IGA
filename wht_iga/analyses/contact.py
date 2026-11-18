# -*- coding: utf-8 -*-
"""Surface-layer contact (Track 3, Phase 5 Python precursor).

Pattern (S22-verified, S25-verified for explicit): IGA bodies keep their
QuadraturePointGeometry elements; each contact face gets a Triangle3D3
surface layer whose nodes are tied 1:1 to the spline control points via
``LinearMasterSlaveConstraint``. Contact traction is applied as nodal
``POINT_LOAD`` pairs updated from the measured gap (penalty), so the
high-level Kratos contact processes (which require ``GenerateFaces`` on the
solid elements) are never needed.

Friction (Coulomb, mu > 0) is supported through the same per-step/per-iteration
update with a tangential stick predictor + slip return mapping (prototyped in
the Track 3 spike tests; production hardening is Phase 5 scope).
"""

from __future__ import annotations

import numpy as np
from numba import njit

from wht_iga.grid.bspline import ders_basis_funs, find_span


def tied_surface_layer(KM, SM, mp, box, *, axis: int, side: int,
                       node_id_start: int, cond_id_start: int,
                       tie_id_start: int, props) -> dict:
    """Build a 1:1 tied Triangle3D3 layer over a box face.

    Layer nodes coincide with the face control points (coordinates read from
    the existing spline nodes, id ``cp_index + 1``), connected as two
    triangles per face quad with outward orientation, and tied component-wise
    to the matching spline node. Returns dict with ``nodes``, ``conditions``,
    ``next_node_id``, ``next_cond_id``, ``next_tie_id``.
    """
    from wht_iga.kratos_backend.constraints import add_linear_tie

    n = list(box.n_cells)
    nu, nv, nw = (int(v) for v in box.n_cp)
    if not (0 <= axis <= 2 and side in (0, 1)):
        raise ValueError(f"axis must be 0..2 and side 0/1, got {axis}, {side}")
    fixed = {0: nu, 1: nv, 2: nw}[axis]
    idx = fixed - 1 if side == 1 else 0
    others = [d for d in range(3) if d != axis]

    def cp_index(q):
        ijk = [0, 0, 0]
        ijk[axis] = idx
        ijk[others[0]] = q[0]
        ijk[others[1]] = q[1]
        return ijk[0] + nu * (ijk[1] + nv * ijk[2])

    counts = [{0: nu, 1: nv, 2: nw}[d] for d in others]
    smp = mp.CreateSubModelPart(f"ContactSurf_{axis}{side}_{node_id_start}")
    grid = {}
    nid = int(node_id_start)
    for r in range(counts[1]):
        for c in range(counts[0]):
            cp = cp_index((c, r))
            ref = mp.GetNode(cp + 1)
            grid[(c, r)] = smp.CreateNewNode(nid, float(ref.X0), float(ref.Y0), float(ref.Z0))
            nid += 1
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, smp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, smp)
    KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, smp)
    # outward unit normal from the box corner ordering
    outward = [0.0, 0.0, 0.0]
    outward[axis] = 1.0 if side == 1 else -1.0
    cid = int(cond_id_start)
    for r in range(counts[1] - 1):
        for c in range(counts[0] - 1):
            a = grid[(c, r)]
            b = grid[(c + 1, r)]
            d = grid[(c, r + 1)]
            e = grid[(c + 1, r + 1)]
            ab = np.array([b.X0 - a.X0, b.Y0 - a.Y0, b.Z0 - a.Z0])
            ad = np.array([d.X0 - a.X0, d.Y0 - a.Y0, d.Z0 - a.Z0])
            flip = bool(np.dot(np.cross(ab, ad), np.asarray(outward, dtype=float)) < 0.0)
            quads = ([a.Id, b.Id, e.Id], [a.Id, e.Id, d.Id])
            if flip:
                quads = ([a.Id, e.Id, b.Id], [a.Id, d.Id, e.Id])
            for tri in quads:
                smp.CreateNewCondition("SurfaceCondition3D3N", cid, list(tri), props)
                cid += 1
    tid = int(tie_id_start)
    for (c, r), nd in grid.items():
        cp = cp_index((c, r))
        master = mp.GetNode(cp + 1)
        for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
            add_linear_tie(KM, mp, tid, [nd.GetDof(var)],
                           [master.GetDof(var)], np.array([[1.0]]))
            tid += 1
    return {"nodes": [grid[(c, r)] for r in range(counts[1]) for c in range(counts[0])],
            "sub_model_part": smp,
            "next_node_id": nid,
            "next_cond_id": cid,
            "next_tie_id": tid}


def paired_gap_z(layer_a, layer_b, KM) -> np.ndarray:
    """Signed z-gaps (b - a) from current nodal positions."""
    ga = np.array([float(n.GetSolutionStepValue(KM.DISPLACEMENT_Z)) + float(n.Z0)
                   for n in layer_a])
    gb = np.array([float(n.GetSolutionStepValue(KM.DISPLACEMENT_Z)) + float(n.Z0)
                   for n in layer_b])
    return gb - ga


# --------------------------------------------------------------------------------------
# Gauss-weighted B-spline surface Semi-Smooth Newton (SSN) contact solver
# --------------------------------------------------------------------------------------


def bspline_basis_matrix(box, points: np.ndarray) -> "scipy.sparse.csr_matrix":
    """Evaluate exact scalar B-spline basis ``N_I(x)`` at physical ``points (M, 3)``.

    Uses :func:`wht_iga.grid.bspline.find_span` and
    :func:`wht_iga.grid.bspline.ders_basis_funs` (exact for any degree ``p >= 1``).
    Control-point column order matches ``BoxSpline``: ``I = a + nu*(b + nv*c)``.
    Returns a canonical CSR matrix of shape ``(M, box.n_ctrl)``.
    """
    import scipy.sparse as sp

    from wht_iga.grid.bspline import ders_basis_funs, find_span

    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim == 1:
        pts = pts.reshape(1, -1)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"points must have shape (M, 3), got {pts.shape}")

    M = pts.shape[0]
    if M == 0:
        return sp.csr_matrix((0, box.n_ctrl), dtype=np.float64)

    p = int(box.degree)
    nu, nv, nw = (int(v) for v in box.n_cp)
    kx, ky, kz = (np.asarray(box.knots[d], dtype=np.float64) for d in range(3))
    n_loc = (p + 1) ** 3

    indices = np.empty(M * n_loc, dtype=np.int64)
    data = np.empty(M * n_loc, dtype=np.float64)
    for q in range(M):
        px, py, pz = float(pts[q, 0]), float(pts[q, 1]), float(pts[q, 2])
        sx = find_span(nu, p, px, kx)
        dx = ders_basis_funs(sx, px, p, 0, kx)[0]
        sy = find_span(nv, p, py, ky)
        dy = ders_basis_funs(sy, py, p, 0, ky)[0]
        sz = find_span(nw, p, pz, kz)
        dz = ders_basis_funs(sz, pz, p, 0, kz)[0]

        ix = np.arange(sx - p, sx + 1, dtype=np.int64)
        iy = np.arange(sy - p, sy + 1, dtype=np.int64)
        iz = np.arange(sz - p, sz + 1, dtype=np.int64)
        I_blk = (ix[None, None, :] + nu * (iy[None, :, None] + nv * iz[:, None, None])).ravel()
        W_blk = (dz[:, None, None] * dy[None, :, None] * dx[None, None, :]).ravel()

        s0 = q * n_loc
        indices[s0 : s0 + n_loc] = I_blk
        data[s0 : s0 + n_loc] = W_blk

    indptr = np.arange(M + 1, dtype=np.int64) * n_loc
    mat = sp.csr_matrix((data, indices, indptr), shape=(M, box.n_ctrl), dtype=np.float64)
    mat.has_sorted_indices = True
    return mat


def expand_vector_basis(B_scalar, ndim: int = 3):
    """Expand scalar basis ``(M, n_ctrl)`` to vector basis ``(ndim*M, ndim*n_ctrl)``."""
    import scipy.sparse as sp

    return sp.kron(B_scalar, sp.eye(int(ndim), format="csr", dtype=np.float64), format="csr")


def sample_box_face_gauss(
    box,
    axis: int,
    side: int,
    n_gp: int = 3,
    bounds: tuple[tuple[float, float], tuple[float, float]] | None = None,
    n_sub: int = 1,
):
    """Sample Gauss-Legendre quadrature points on a ``BoxSpline`` boundary face.

    Iterates over face cells of ``axis`` (``side in (0, 1)``), optionally clipped
    to an active tangential window ``bounds = ((u_lo, u_hi), (v_lo, v_hi))`` and
    subdivided into ``n_sub x n_sub`` uniform sub-cells per cell for accurate
    patch integration, preserving exact physical weights ``w_q = jac_a * jac_b * wa * wb``.

    Returns ``(points (M, 3), weights (M,), B_scalar (M, box.n_ctrl))``.
    """
    if not (0 <= axis <= 2 and side in (0, 1)):
        raise ValueError(f"axis must be 0..2 and side 0/1, got {axis}, {side}")
    if int(n_gp) < 1 or int(n_sub) < 1:
        raise ValueError(f"n_gp and n_sub must be >= 1, got {n_gp}, {n_sub}")

    xg, wg = np.polynomial.legendre.leggauss(int(n_gp))
    others = [d for d in range(3) if d != axis]
    x_fixed = float(box.lo[axis] if side == 0 else box.hi[axis])

    pts_1d: list[np.ndarray] = []
    wts_1d: list[np.ndarray] = []
    for idx, d in enumerate(others):
        br = box.breakpoints(d)
        p_list, w_list = [], []
        for a, b in zip(br[:-1], br[1:]):
            aa, bb = float(a), float(b)
            if bounds is not None:
                lo_b, hi_b = float(bounds[idx][0]), float(bounds[idx][1])
                aa = max(aa, lo_b)
                bb = min(bb, hi_b)
            if bb <= aa + 1e-15:
                continue
            edges = np.linspace(aa, bb, int(n_sub) + 1)
            for s0, s1 in zip(edges[:-1], edges[1:]):
                jac = 0.5 * (s1 - s0)
                p_list.append(0.5 * (s0 + s1) + jac * xg)
                w_list.append(jac * wg)
        if not p_list:
            empty_pts = np.empty((0, 3), dtype=np.float64)
            empty_wts = np.empty(0, dtype=np.float64)
            return empty_pts, empty_wts, bspline_basis_matrix(box, empty_pts)
        pts_1d.append(np.concatenate(p_list))
        wts_1d.append(np.concatenate(w_list))

    U0, U1 = np.meshgrid(pts_1d[0], pts_1d[1], indexing="ij")
    W0, W1 = np.meshgrid(wts_1d[0], wts_1d[1], indexing="ij")
    weights = (W0 * W1).ravel()
    points = np.empty((weights.size, 3), dtype=np.float64)
    points[:, axis] = x_fixed
    points[:, others[0]] = U0.ravel()
    points[:, others[1]] = U1.ravel()

    B_scalar = bspline_basis_matrix(box, points)
    return points, weights, B_scalar


class RigidSphereObstacle:
    """Rigid sphere obstacle ``||x - center|| >= radius``."""

    def __init__(self, center, radius: float) -> None:
        c = np.asarray(center, dtype=np.float64).reshape(3)
        r = float(radius)
        if r <= 0.0:
            raise ValueError(f"radius must be > 0, got {radius}")
        self.center = c
        self.radius = r

    def gap_and_normal(self, x_q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(g_n, n, curv)`` with ``g_n = ||x - c|| - R`` and ``n = (x - c) / r``."""
        d_vec = x_q - self.center[None, :]
        r = np.linalg.norm(d_vec, axis=1)
        if np.any(r <= 1e-15):
            raise ValueError("contact point coincides with rigid sphere center")
        n = d_vec / r[:, None]
        g_n = r - self.radius
        curv = g_n / r
        return g_n, n, curv


class RigidPlaneObstacle:
    """Rigid half-space obstacle ``(x - point) . normal >= 0`` (``normal`` into body)."""

    def __init__(self, point, normal) -> None:
        p0 = np.asarray(point, dtype=np.float64).reshape(3)
        n = np.asarray(normal, dtype=np.float64).reshape(3)
        norm = float(np.linalg.norm(n))
        if norm <= 1e-15:
            raise ValueError("plane normal must be nonzero")
        self.point = p0
        self.normal = n / norm

    def gap_and_normal(self, x_q: np.ndarray) -> tuple[np.ndarray, np.ndarray, None]:
        """Return ``(g_n, n, None)`` with ``g_n = (x - point) . normal``."""
        g_n = (x_q - self.point[None, :]) @ self.normal
        n = np.broadcast_to(self.normal, x_q.shape).copy()
        return g_n, n, None


def _solve_linear_step(A, b: np.ndarray, use_pardiso: bool = True) -> np.ndarray:
    """Solve ``A x = b`` via pypardiso (if available) or SciPy sparse LU."""
    if use_pardiso:
        try:
            import pypardiso

            return np.asarray(pypardiso.spsolve(A.tocsr(), b), dtype=np.float64)
        except ImportError:
            pass
    import scipy.sparse.linalg as spla

    return np.asarray(spla.spsolve(A.tocsc(), b), dtype=np.float64)


def solve_ssn_contact(
    Kf,
    X_q: np.ndarray,
    w_q: np.ndarray,
    B_scalar,
    obstacle,
    eps_n: float,
    *,
    free_mask: np.ndarray | None = None,
    u_fixed: np.ndarray | None = None,
    f_ext: np.ndarray | None = None,
    u_init: np.ndarray | None = None,
    mu: float = 0.0,
    eps_t: float | None = None,
    slip_state: dict[str, np.ndarray] | None = None,
    obstacle_disp=None,
    use_curvature: bool = True,
    max_iter: int = 50,
    rtol: float = 1e-10,
    use_pardiso: bool = True,
) -> dict:
    """Solve Gauss-weighted penalty contact equilibrium via Semi-Smooth Newton.

    Solves ``Kf @ u - f_ext + f_c(u) = 0`` on free DOFs where ``f_c(u)`` is the
    Gauss-weighted B-spline surface integral of penalty normal traction and
    optional Coulomb friction (stick predictor + radial return mapping).

    Curvature ``(g_n / r) * (I - n x n)`` is enabled once the active set settles
    and ``rnorm < 0.1`` (with automatic fallback to the SPD ``n x n`` tangent if
    a step increases the residual or cycles), preventing 3-cycle overshoot from
    negative curvature at large initial penetrations.
    """
    import scipy.sparse as sp

    if float(eps_n) <= 0.0:
        raise ValueError(f"eps_n must be > 0, got {eps_n}")
    mu = float(mu)
    if mu < 0.0:
        raise ValueError(f"mu must be >= 0, got {mu}")
    if mu > 0.0:
        eps_t = float(eps_n if eps_t is None else eps_t)
        if eps_t <= 0.0:
            raise ValueError(f"eps_t must be > 0 when mu > 0, got {eps_t}")

    X_q = np.asarray(X_q, dtype=np.float64)
    w_q = np.asarray(w_q, dtype=np.float64)
    M, n_ctrl = B_scalar.shape
    n_dofs = 3 * n_ctrl
    if free_mask is None:
        free = np.ones(n_dofs, dtype=bool)
    else:
        free = np.asarray(free_mask, dtype=bool)
        if free.shape != (n_dofs,):
            raise ValueError(f"free_mask shape {free.shape} != ({n_dofs},)")
    n_free = int(free.sum())
    Kf_csr = Kf.tocsr()
    if Kf_csr.shape != (n_free, n_free):
        raise ValueError(f"Kf shape {Kf_csr.shape} != ({n_free}, {n_free})")

    B_vec = expand_vector_basis(B_scalar, ndim=3)
    B_free = B_vec[:, free].tocsr()

    if u_fixed is not None:
        u_fix_full = np.asarray(u_fixed, dtype=np.float64).reshape(n_dofs)
        x_base = X_q + (B_vec[:, ~free] @ u_fix_full[~free]).reshape(M, 3)
    else:
        u_fix_full = None
        x_base = X_q

    f_ext_vec = np.zeros(n_free, dtype=np.float64) if f_ext is None else np.asarray(f_ext, dtype=np.float64).reshape(n_free)
    u = np.zeros(n_free, dtype=np.float64) if u_init is None else np.asarray(u_init, dtype=np.float64).reshape(n_free).copy()

    if slip_state is None:
        x_prev = X_q.copy()
        t_t_prev = np.zeros((M, 3), dtype=np.float64)
    else:
        x_prev = np.asarray(slip_state["x_prev"], dtype=np.float64).reshape(M, 3)
        t_t_prev = np.asarray(slip_state["t_t_prev"], dtype=np.float64).reshape(M, 3)
    dx_obs = np.zeros(3, dtype=np.float64) if obstacle_disp is None else np.asarray(obstacle_disp, dtype=np.float64).reshape(3)

    gap_fn = obstacle.gap_and_normal if hasattr(obstacle, "gap_and_normal") else obstacle
    I3 = np.eye(3, dtype=np.float64)

    def _evaluate(u_vec: np.ndarray, allow_curv: bool):
        x_q = x_base + (B_free @ u_vec).reshape(M, 3)
        g_n, n, curv = gap_fn(x_q)
        active = g_n < 0.0
        act_idx = np.flatnonzero(active)
        t_res = np.zeros((M, 3), dtype=np.float64)
        t_t_res = np.zeros((M, 3), dtype=np.float64)
        slip_mask = np.zeros(M, dtype=bool)
        if act_idx.size == 0:
            f_c = np.zeros(n_free, dtype=np.float64)
            res = Kf_csr @ u_vec - f_ext_vec
            denom = max(1.0, float(np.linalg.norm(f_ext_vec)))
            rnorm = float(np.linalg.norm(res) / denom)
            return x_q, g_n, n, curv, active, act_idx, slip_mask, t_res, t_t_res, f_c, res, rnorm, None

        w_act = w_q[act_idx]
        g_act = g_n[act_idx]
        n_act = n[act_idx]
        nn = np.einsum("qi,qj->qij", n_act, n_act)
        P_t = I3[None, :, :] - nn
        t_n = (eps_n * g_act)[:, None] * n_act

        if allow_curv and curv is not None:
            k_loc = (w_act * eps_n)[:, None, None] * (nn + curv[act_idx, None, None] * P_t)
        else:
            k_loc = (w_act * eps_n)[:, None, None] * nn

        if mu > 0.0:
            du_rel = (x_q[act_idx] - x_prev[act_idx]) - dx_obs[None, :]
            du_t = np.einsum("qij,qj->qi", P_t, du_rel)
            t_prev_proj = np.einsum("qij,qj->qi", P_t, t_t_prev[act_idx])
            t_trial = t_prev_proj + float(eps_t) * du_t
            tn = np.linalg.norm(t_trial, axis=1)
            cap = -mu * eps_n * g_act
            slip = tn > cap
            slip_mask[act_idx] = slip

            t_t_act = t_trial.copy()
            k_t = (w_act * float(eps_t))[:, None, None] * P_t
            if np.any(slip):
                s_idx = np.flatnonzero(slip)
                d_hat = t_trial[s_idx] / tn[s_idx, None]
                t_t_act[s_idx] = cap[s_idx, None] * d_hat
                dd = np.einsum("qi,qj->qij", d_hat, d_hat)
                dn = np.einsum("qi,qj->qij", d_hat, n_act[s_idx])
                k_t[s_idx] = w_act[s_idx, None, None] * (
                    (float(eps_t) * (cap[s_idx] / tn[s_idx]))[:, None, None] * (P_t[s_idx] - dd)
                    - (mu * eps_n) * dn
                )
            t_t_res[act_idx] = t_t_act
            t_res[act_idx] = w_act[:, None] * (t_n + t_t_act)
            k_loc = k_loc + k_t
        else:
            t_res[act_idx] = w_act[:, None] * t_n

        rows = (3 * act_idx[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
        B_act = B_free[rows, :]
        f_c = B_act.T @ t_res[act_idx].ravel()
        res = Kf_csr @ u_vec - f_ext_vec + f_c
        denom = max(1.0, float(np.linalg.norm(f_ext_vec)), float(np.linalg.norm(f_c)))
        rnorm = float(np.linalg.norm(res) / denom)
        return x_q, g_n, n, curv, active, act_idx, slip_mask, t_res, t_t_res, f_c, res, rnorm, (B_act, k_loc)

    prev_active = None
    prev_slip = None
    allow_curv = False
    converged = False
    n_solves = 0

    for _ in range(1, int(max_iter) + 1):
        (
            x_q, g_n, n, curv, active, act_idx, slip_mask,
            t_res, t_t_res, _f_c, res, rnorm, tan_data
        ) = _evaluate(u, allow_curv)

        same_active = prev_active is not None and np.array_equal(active, prev_active)
        same_slip = prev_slip is not None and np.array_equal(slip_mask, prev_slip)
        if same_active and same_slip and rnorm <= rtol:
            converged = True
            break

        # Enable exact curvature term only inside the active-set basin; fall back to SPD if needed.
        if not (use_curvature and curv is not None and same_active and rnorm < 0.1):
            allow_curv = False
        elif not allow_curv:
            allow_curv = True
            _, _, _, _, _, _, _, _, _, _, _, _, tan_data = _evaluate(u, True)

        prev_active = active.copy()
        prev_slip = slip_mask.copy()

        if tan_data is None:
            K_tan = Kf_csr
        else:
            B_act, k_loc = tan_data
            n_a = act_idx.size
            K_blk = sp.bsr_matrix(
                (k_loc, np.arange(n_a, dtype=np.int64), np.arange(n_a + 1, dtype=np.int64)),
                shape=(3 * n_a, 3 * n_a),
            ).tocsr()
            K_tan = Kf_csr + B_act.T @ K_blk @ B_act

        du = _solve_linear_step(K_tan, -res, use_pardiso=use_pardiso)
        n_solves += 1
        u_next = u + du

        # Guard against curvature-induced overshoot when active set was assumed settled.
        if allow_curv:
            _, _, _, _, act_next, _, _, _, _, _, _, rn_next, _ = _evaluate(u_next, False)
            if not np.array_equal(act_next, active) or rn_next > max(rnorm, rtol):
                allow_curv = False
                _, _, _, _, _, _, _, _, _, _, _, _, tan_spd = _evaluate(u, False)
                if tan_spd is not None:
                    B_act, k_loc = tan_spd
                    n_a = act_idx.size
                    K_blk = sp.bsr_matrix(
                        (k_loc, np.arange(n_a, dtype=np.int64), np.arange(n_a + 1, dtype=np.int64)),
                        shape=(3 * n_a, 3 * n_a),
                    ).tocsr()
                    du = _solve_linear_step(Kf_csr + B_act.T @ K_blk @ B_act, -res, use_pardiso=use_pardiso)
                    n_solves += 1
                    u_next = u + du

        # Residual backtracking when friction is active to step smoothly across stick/slip kinks.
        if mu > 0.0 and float(np.linalg.norm(du)) > 0.0:
            res_norm0 = float(np.linalg.norm(res))
            alpha = 1.0
            while alpha >= 1.0 / 64.0:
                _, _, _, _, _, _, _, _, _, _, res_try, _, _ = _evaluate(u + alpha * du, False)
                if float(np.linalg.norm(res_try)) < (1.0 - 1e-4 * alpha) * res_norm0:
                    break
                alpha *= 0.5
            u_next = u + alpha * du

        u = u_next

    if not converged:
        raise RuntimeError(
            f"SSN contact did not converge in {max_iter} iterations "
            f"(rnorm={rnorm:.3e}, n_active={int(act_idx.size)})"
        )

    u_full = np.zeros(n_dofs, dtype=np.float64)
    u_full[free] = u
    if u_fix_full is not None:
        u_full[~free] = u_fix_full[~free]

    return {
        "u_free": u,
        "u_ctrl": u_full.reshape(n_ctrl, 3),
        "active": active,
        "n_active": int(act_idx.size),
        "n_slip": int(slip_mask.sum()),
        "gap": g_n,
        "normal": n,
        "traction_res": t_res,
        "force_total": t_res.sum(axis=0),
        "normal_pressure_integral": float(np.sum(-w_q[act_idx] * eps_n * g_n[act_idx])) if act_idx.size else 0.0,
        "max_pressure": float(np.max(-eps_n * g_n[act_idx])) if act_idx.size else 0.0,
        "iterations": int(n_solves),
        "rnorm": float(rnorm),
        "slip_state": {"x_prev": x_q.copy(), "t_t_prev": t_t_res.copy()},
    }


# --------------------------------------------------------------------------------------
# Deformable-deformable surface Gauss SSN contact (V08 Cases R1 & R2)
# --------------------------------------------------------------------------------------


def sample_surface_triangles_gauss(
    nodes: np.ndarray,
    tris: np.ndarray,
    box_slave,
    order: int = 2,
):
    """Sample Dunavant triangle quadrature on physical slave triangles.

    Parameters
    ----------
    nodes:
        Vertex coordinates ``(N, 3)``.
    tris:
        Triangle vertex connectivity ``(K, 3)`` indexing ``nodes``.
    box_slave:
        ``BoxSpline`` discretizing the slave body.
    order:
        Dunavant triangle quadrature order (1, 2, 4, or 6).

    Returns
    -------
    tuple[np.ndarray, np.ndarray, scipy.sparse.csr_matrix]
        ``(X_s (M, 3), w_q (M,), B_s (M, box_slave.n_ctrl))`` with exact
        physical area weights ``w_q = 2 * tri_area * ref_w`` satisfying
        ``sum(w_q) == sum(tri_area)``.
    """
    from wht_iga.kratos_backend.surface import triangle_quadrature

    pts_arr = np.asarray(nodes, dtype=np.float64)
    tri_arr = np.asarray(tris, dtype=np.int64)
    if tri_arr.ndim != 2 or tri_arr.shape[1] != 3:
        raise ValueError(f"tris must have shape (K, 3), got {tri_arr.shape}")
    if tri_arr.shape[0] == 0:
        empty_pts = np.empty((0, 3), dtype=np.float64)
        empty_wts = np.empty(0, dtype=np.float64)
        return empty_pts, empty_wts, bspline_basis_matrix(box_slave, empty_pts)

    ref_pts, ref_w = triangle_quadrature(int(order))
    tri_xyz = pts_arr[tri_arr]  # (K, 3, 3)
    v0 = tri_xyz[:, 0, :]
    e1 = tri_xyz[:, 1, :] - v0
    e2 = tri_xyz[:, 2, :] - v0
    tri_area = 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1)

    X_s = (
        v0[:, None, :]
        + ref_pts[None, :, 0:1] * e1[:, None, :]
        + ref_pts[None, :, 1:2] * e2[:, None, :]
    ).reshape(-1, 3)
    w_q = (2.0 * tri_area[:, None] * ref_w[None, :]).ravel()
    B_s = bspline_basis_matrix(box_slave, X_s)
    return X_s, w_q, B_s


@njit(cache=True)
def _eval_face_1st_njit(
    a: float,
    b: float,
    X_face: np.ndarray,
    p: int,
    n_a: int,
    n_b: int,
    ka: np.ndarray,
    kb: np.ndarray,
):
    sa = find_span(n_a, p, a, ka)
    sb = find_span(n_b, p, b, kb)
    da = ders_basis_funs(sa, a, p, 1, ka)
    db = ders_basis_funs(sb, b, p, 1, kb)
    ia0 = sa - p
    jb0 = sb - p
    xm0 = xm1 = xm2 = 0.0
    ta0 = ta1 = ta2 = 0.0
    tb0 = tb1 = tb2 = 0.0
    for i in range(p + 1):
        na0 = da[0, i]
        na1 = da[1, i]
        ia = ia0 + i
        for j in range(p + 1):
            nb0 = db[0, j]
            nb1 = db[1, j]
            jb = jb0 + j
            w = na0 * nb0
            wa = na1 * nb0
            wb = na0 * nb1
            cx = X_face[ia, jb, 0]
            cy = X_face[ia, jb, 1]
            cz = X_face[ia, jb, 2]
            xm0 += w * cx
            xm1 += w * cy
            xm2 += w * cz
            ta0 += wa * cx
            ta1 += wa * cy
            ta2 += wa * cz
            tb0 += wb * cx
            tb1 += wb * cy
            tb2 += wb * cz
    return sa, sb, da, db, xm0, xm1, xm2, ta0, ta1, ta2, tb0, tb1, tb2


@njit(cache=True)
def _eval_face_0th_njit(
    a: float,
    b: float,
    X_face: np.ndarray,
    p: int,
    n_a: int,
    n_b: int,
    ka: np.ndarray,
    kb: np.ndarray,
):
    sa = find_span(n_a, p, a, ka)
    sb = find_span(n_b, p, b, kb)
    na = ders_basis_funs(sa, a, p, 0, ka)[0]
    nb = ders_basis_funs(sb, b, p, 0, kb)[0]
    ia0 = sa - p
    jb0 = sb - p
    xm0 = xm1 = xm2 = 0.0
    for i in range(p + 1):
        nai = na[i]
        ia = ia0 + i
        for j in range(p + 1):
            w = nai * nb[j]
            jb = jb0 + j
            xm0 += w * X_face[ia, jb, 0]
            xm1 += w * X_face[ia, jb, 1]
            xm2 += w * X_face[ia, jb, 2]
    return xm0, xm1, xm2


@njit(cache=True)
def _project_box_face_njit(
    x_s: np.ndarray,
    X_face: np.ndarray,
    p: int,
    n_a: int,
    n_b: int,
    ka: np.ndarray,
    kb: np.ndarray,
    lo_a: float,
    hi_a: float,
    lo_b: float,
    hi_b: float,
    d_a: int,
    d_b: int,
    orient: float,
):
    M = x_s.shape[0]
    ab_out = np.empty((M, 2), dtype=np.float64)
    x_m_out = np.empty((M, 3), dtype=np.float64)
    n_m_out = np.empty((M, 3), dtype=np.float64)
    g_n_out = np.empty(M, dtype=np.float64)
    sa_out = np.empty(M, dtype=np.int64)
    sb_out = np.empty(M, dtype=np.int64)
    Na_out = np.empty((M, p + 1), dtype=np.float64)
    Nb_out = np.empty((M, p + 1), dtype=np.float64)

    span_len = max(1.0, max(hi_a - lo_a, hi_b - lo_b))
    tol = 1e-12 * span_len

    for q in range(M):
        xs0 = x_s[q, 0]
        xs1 = x_s[q, 1]
        xs2 = x_s[q, 2]
        a = x_s[q, d_a]
        b = x_s[q, d_b]
        if a < lo_a:
            a = lo_a
        elif a > hi_a:
            a = hi_a
        if b < lo_b:
            b = lo_b
        elif b > hi_b:
            b = hi_b

        for _ in range(20):
            sa, sb, da, db, xm0, xm1, xm2, ta0, ta1, ta2, tb0, tb1, tb2 = _eval_face_1st_njit(
                a, b, X_face, p, n_a, n_b, ka, kb
            )
            r0 = xm0 - xs0
            r1 = xm1 - xs1
            r2 = xm2 - xs2
            r2_sq = r0 * r0 + r1 * r1 + r2 * r2
            g0 = ta0 * r0 + ta1 * r1 + ta2 * r2
            g1 = tb0 * r0 + tb1 * r1 + tb2 * r2
            H00 = ta0 * ta0 + ta1 * ta1 + ta2 * ta2
            H01 = ta0 * tb0 + ta1 * tb1 + ta2 * tb2
            H11 = tb0 * tb0 + tb1 * tb1 + tb2 * tb2
            det = H00 * H11 - H01 * H01
            if det <= 1e-24 * max(1.0, H00 * H11):
                break
            da_step = (-g0 * H11 + g1 * H01) / det
            db_step = (-H00 * g1 + H01 * g0) / det

            alpha = 1.0
            a_next = a
            b_next = b
            for _ls in range(5):
                a_try = a + alpha * da_step
                b_try = b + alpha * db_step
                if a_try < lo_a:
                    a_try = lo_a
                elif a_try > hi_a:
                    a_try = hi_a
                if b_try < lo_b:
                    b_try = lo_b
                elif b_try > hi_b:
                    b_try = hi_b
                if alpha == 1.0 and abs(a_try - a) <= tol and abs(b_try - b) <= tol:
                    a_next = a_try
                    b_next = b_try
                    break
                tx, ty, tz = _eval_face_0th_njit(a_try, b_try, X_face, p, n_a, n_b, ka, kb)
                d0 = tx - xs0
                d1 = ty - xs1
                d2 = tz - xs2
                if d0 * d0 + d1 * d1 + d2 * d2 <= r2_sq + 1e-18 or alpha <= 0.125:
                    a_next = a_try
                    b_next = b_try
                    break
                alpha *= 0.5

            step_diff = max(abs(a_next - a), abs(b_next - b))
            a = a_next
            b = b_next
            if step_diff <= tol:
                break

        sa, sb, da, db, xm0, xm1, xm2, ta0, ta1, ta2, tb0, tb1, tb2 = _eval_face_1st_njit(
            a, b, X_face, p, n_a, n_b, ka, kb
        )
        nx = orient * (ta1 * tb2 - ta2 * tb1)
        ny = orient * (ta2 * tb0 - ta0 * tb2)
        nz = orient * (ta0 * tb1 - ta1 * tb0)
        nlen = np.sqrt(nx * nx + ny * ny + nz * nz)
        if nlen > 1e-15:
            nx /= nlen
            ny /= nlen
            nz /= nlen
        g_n = (xs0 - xm0) * nx + (xs1 - xm1) * ny + (xs2 - xm2) * nz

        ab_out[q, 0] = a
        ab_out[q, 1] = b
        x_m_out[q, 0] = xm0
        x_m_out[q, 1] = xm1
        x_m_out[q, 2] = xm2
        n_m_out[q, 0] = nx
        n_m_out[q, 1] = ny
        n_m_out[q, 2] = nz
        g_n_out[q] = g_n
        sa_out[q] = sa
        sb_out[q] = sb
        for i in range(p + 1):
            Na_out[q, i] = da[0, i]
            Nb_out[q, i] = db[0, i]

    return ab_out, x_m_out, n_m_out, g_n_out, sa_out, sb_out, Na_out, Nb_out


class DeformableBoxFaceMaster:
    """Closest-point projector and basis evaluator on a deformed ``BoxSpline`` face.

    Parameters
    ----------
    box_master:
        ``BoxSpline`` discretizing the master body.
    axis:
        Face normal coordinate index (``0``, ``1``, or ``2``).
    side:
        Face side (``0`` for ``lo[axis]``, ``1`` for ``hi[axis]``).
    """

    def __init__(self, box_master, axis: int = 1, side: int = 1) -> None:
        if not (0 <= axis <= 2 and side in (0, 1)):
            raise ValueError(f"axis must be 0..2 and side 0/1, got {axis}, {side}")
        self.box = box_master
        self.axis = int(axis)
        self.side = int(side)
        self.p = int(box_master.degree)
        self.n_ctrl = int(box_master.n_ctrl)

        others = [d for d in range(3) if d != self.axis]
        self.d_a, self.d_b = others[0], others[1]
        nu, nv, nw = (int(v) for v in box_master.n_cp)
        n_cp = (nu, nv, nw)
        self.n_a, self.n_b = n_cp[self.d_a], n_cp[self.d_b]
        k_fixed = 0 if self.side == 0 else n_cp[self.axis] - 1

        self.ka = np.ascontiguousarray(box_master.knots[self.d_a], dtype=np.float64)
        self.kb = np.ascontiguousarray(box_master.knots[self.d_b], dtype=np.float64)
        self.lo_a = float(box_master.lo[self.d_a])
        self.hi_a = float(box_master.hi[self.d_a])
        self.lo_b = float(box_master.lo[self.d_b])
        self.hi_b = float(box_master.hi[self.d_b])
        # Outward orientation: e_{d_a} x e_{d_b} = (-1)^axis e_{axis}.
        self.orient = (1.0 if self.side == 1 else -1.0) * (-1.0 if self.axis == 1 else 1.0)

        ia = np.arange(self.n_a, dtype=np.int64)[:, None]
        ib = np.arange(self.n_b, dtype=np.int64)[None, :]
        ijk = [None, None, None]
        ijk[self.axis] = np.full((self.n_a, self.n_b), k_fixed, dtype=np.int64)
        ijk[self.d_a] = np.broadcast_to(ia, (self.n_a, self.n_b))
        ijk[self.d_b] = np.broadcast_to(ib, (self.n_a, self.n_b))
        self.face_cp_idx = (ijk[0] + nu * (ijk[1] + nv * ijk[2])).astype(np.int64)
        self.X_face_0 = np.ascontiguousarray(
            box_master.control_points()[self.face_cp_idx], dtype=np.float64
        )

    def evaluate(
        self,
        x_s: np.ndarray,
        U_m: np.ndarray,
    ):
        """Project deformed slave points ``x_s (M, 3)`` onto deformed master face.

        Returns ``(g_n (M,), n_m (M, 3), x_m (M, 3), ab_m (M, 2), B_m (M, n_ctrl_m))``.
        ``n_m`` is the outward unit normal of ``box_master``; ``g_n = (x_s - x_m) . n_m``
        is negative on penetration.
        """
        import scipy.sparse as sp

        xs = np.ascontiguousarray(x_s, dtype=np.float64)
        if xs.ndim == 1:
            xs = xs.reshape(1, 3)
        M = xs.shape[0]
        if M == 0:
            return (
                np.empty(0, dtype=np.float64),
                np.empty((0, 3), dtype=np.float64),
                np.empty((0, 3), dtype=np.float64),
                np.empty((0, 2), dtype=np.float64),
                sp.csr_matrix((0, self.n_ctrl), dtype=np.float64),
            )

        um = np.asarray(U_m, dtype=np.float64).reshape(self.n_ctrl, 3)
        X_face = np.ascontiguousarray(self.X_face_0 + um[self.face_cp_idx], dtype=np.float64)

        ab_m, x_m, n_m, g_n, sa, sb, Na, Nb = _project_box_face_njit(
            xs,
            X_face,
            self.p,
            self.n_a,
            self.n_b,
            self.ka,
            self.kb,
            self.lo_a,
            self.hi_a,
            self.lo_b,
            self.hi_b,
            self.d_a,
            self.d_b,
            self.orient,
        )

        p1 = self.p + 1
        n_loc = p1 * p1
        off = np.arange(-self.p, 1, dtype=np.int64)
        ia_mat = sa[:, None] + off[None, :]  # (M, p+1)
        ib_mat = sb[:, None] + off[None, :]  # (M, p+1)
        col_blk = self.face_cp_idx[ia_mat[:, :, None], ib_mat[:, None, :]]  # (M, p+1, p+1)
        val_blk = Na[:, :, None] * Nb[:, None, :]
        indptr = np.arange(M + 1, dtype=np.int64) * n_loc
        B_m = sp.csr_matrix(
            (val_blk.ravel(), col_blk.ravel(), indptr),
            shape=(M, self.n_ctrl),
            dtype=np.float64,
        )
        B_m.sort_indices()
        return g_n, n_m, x_m, ab_m, B_m


def solve_deformable_ssn_contact(
    Kf_master,
    Kf_slave,
    X_s: np.ndarray,
    w_q: np.ndarray,
    B_s,
    master_face: DeformableBoxFaceMaster,
    eps_n: float,
    *,
    free_m: np.ndarray | None = None,
    u_fixed_m: np.ndarray | None = None,
    free_s: np.ndarray | None = None,
    u_fixed_s: np.ndarray | None = None,
    f_ext_m: np.ndarray | None = None,
    f_ext_s: np.ndarray | None = None,
    u_init_m: np.ndarray | None = None,
    u_init_s: np.ndarray | None = None,
    mu: float = 0.0,
    eps_t: float | None = None,
    slip_state: dict[str, np.ndarray] | None = None,
    max_iter: int = 50,
    rtol: float = 1e-10,
    use_pardiso: bool = True,
) -> dict:
    """Solve deformable-deformable surface Gauss SSN contact equilibrium.

    Assembles the coupled equilibrium system on free DOFs ``(u_m, u_s)``:

    ``[Kf_m  0  ] [u_m]   [-B_m_free^T t_res]   [f_ext_m]``
    ``[  0  Kf_s] [u_s] + [+B_s_free^T t_res] = [f_ext_s]``

    with Gauss-weighted penalty normal traction ``t_n = eps_n * g_n * n_m`` on
    active slave Gauss points (``g_n = (x_s - x_m) . n_m < 0``) and optional
    Coulomb friction (``mu > 0``). Here ``t_res = w_q * (t_n + t_t)`` is the
    force exerted **by the slave on the master** (so ``-t_res`` is exerted on the
    slave by the master).

    Convected tangential slip increment (``mu > 0``):
    At each slave Gauss point ``q``, the slave material displacement increment is
    ``x_s - x_s_prev``, while the master surface particle currently under ``x_s``
    at closest-point parameter ``xi_curr = (a*, b*)`` has Lagrangian material
    displacement increment ``B_m(xi_curr) @ (U_m - U_m_prev)``. Thus the relative
    tangential slip increment between the contacting material particles is
    ``du_rel = (x_s - x_s_prev) - B_m(xi_curr) @ (U_m - U_m_prev)``, which
    vanishes under simultaneous rigid translation and stays exact under finite
    sliding across the master face (using ``x_m(xi_curr) - x_m(xi_prev)`` would
    spuriously cancel rigid sliding of the slave).
    """
    import scipy.sparse as sp

    if float(eps_n) <= 0.0:
        raise ValueError(f"eps_n must be > 0, got {eps_n}")
    mu = float(mu)
    if mu < 0.0:
        raise ValueError(f"mu must be >= 0, got {mu}")
    if mu > 0.0:
        eps_t = float(eps_n if eps_t is None else eps_t)
        if eps_t <= 0.0:
            raise ValueError(f"eps_t must be > 0 when mu > 0, got {eps_t}")

    X_s = np.asarray(X_s, dtype=np.float64)
    w_q = np.asarray(w_q, dtype=np.float64)
    M, n_ctrl_s = B_s.shape
    n_ctrl_m = int(master_face.n_ctrl)
    n_dofs_m = 3 * n_ctrl_m
    n_dofs_s = 3 * n_ctrl_s

    free_m_arr = np.ones(n_dofs_m, dtype=bool) if free_m is None else np.asarray(free_m, dtype=bool)
    free_s_arr = np.ones(n_dofs_s, dtype=bool) if free_s is None else np.asarray(free_s, dtype=bool)
    if free_m_arr.shape != (n_dofs_m,):
        raise ValueError(f"free_m shape {free_m_arr.shape} != ({n_dofs_m},)")
    if free_s_arr.shape != (n_dofs_s,):
        raise ValueError(f"free_s shape {free_s_arr.shape} != ({n_dofs_s},)")

    n_free_m = int(free_m_arr.sum())
    n_free_s = int(free_s_arr.sum())
    Kf_m_csr = Kf_master.tocsr()
    Kf_s_csr = Kf_slave.tocsr()
    if Kf_m_csr.shape != (n_free_m, n_free_m):
        raise ValueError(f"Kf_master shape {Kf_m_csr.shape} != ({n_free_m}, {n_free_m})")
    if Kf_s_csr.shape != (n_free_s, n_free_s):
        raise ValueError(f"Kf_slave shape {Kf_s_csr.shape} != ({n_free_s}, {n_free_s})")

    K_diag = sp.block_diag([Kf_m_csr, Kf_s_csr], format="csr")

    u_m_base = np.zeros(n_dofs_m, dtype=np.float64)
    if u_fixed_m is not None:
        u_fix_m_in = np.asarray(u_fixed_m, dtype=np.float64).reshape(n_dofs_m)
        u_m_base[~free_m_arr] = u_fix_m_in[~free_m_arr]

    u_s_base = np.zeros(n_dofs_s, dtype=np.float64)
    if u_fixed_s is not None:
        u_fix_s_in = np.asarray(u_fixed_s, dtype=np.float64).reshape(n_dofs_s)
        u_s_base[~free_s_arr] = u_fix_s_in[~free_s_arr]

    B_s_vec = expand_vector_basis(B_s, ndim=3)
    B_s_free = B_s_vec[:, free_s_arr].tocsr()
    x_s_base = X_s + (B_s_vec[:, ~free_s_arr] @ u_s_base[~free_s_arr]).reshape(M, 3)

    f_ext_m_vec = np.zeros(n_free_m, dtype=np.float64) if f_ext_m is None else np.asarray(f_ext_m, dtype=np.float64).reshape(n_free_m)
    f_ext_s_vec = np.zeros(n_free_s, dtype=np.float64) if f_ext_s is None else np.asarray(f_ext_s, dtype=np.float64).reshape(n_free_s)

    u_m = np.zeros(n_free_m, dtype=np.float64) if u_init_m is None else np.asarray(u_init_m, dtype=np.float64).reshape(n_free_m).copy()
    u_s = np.zeros(n_free_s, dtype=np.float64) if u_init_s is None else np.asarray(u_init_s, dtype=np.float64).reshape(n_free_s).copy()
    u = np.concatenate([u_m, u_s])

    if slip_state is None:
        x_s_prev = X_s.copy()
        U_m_prev = np.zeros((n_ctrl_m, 3), dtype=np.float64)
        t_t_prev = np.zeros((M, 3), dtype=np.float64)
    else:
        x_s_prev = np.asarray(slip_state["x_s_prev"], dtype=np.float64).reshape(M, 3)
        U_m_prev = np.asarray(slip_state["U_m_prev"], dtype=np.float64).reshape(n_ctrl_m, 3)
        t_t_prev = np.asarray(slip_state["t_t_prev"], dtype=np.float64).reshape(M, 3)

    I3 = np.eye(3, dtype=np.float64)

    def _unpack_full_m(um_free: np.ndarray) -> np.ndarray:
        full = u_m_base.copy()
        full[free_m_arr] = um_free
        return full.reshape(n_ctrl_m, 3)

    def _unpack_full_s(us_free: np.ndarray) -> np.ndarray:
        full = u_s_base.copy()
        full[free_s_arr] = us_free
        return full.reshape(n_ctrl_s, 3)

    def _evaluate(u_vec: np.ndarray, need_tangent: bool = True):
        um_f = u_vec[:n_free_m]
        us_f = u_vec[n_free_m:]
        U_m_curr = _unpack_full_m(um_f)
        x_s = x_s_base + (B_s_free @ us_f).reshape(M, 3)
        g_n, n_m, x_m, ab_m, B_m = master_face.evaluate(x_s, U_m_curr)

        active = g_n < 0.0
        act_idx = np.flatnonzero(active)
        t_res = np.zeros((M, 3), dtype=np.float64)
        t_t_res = np.zeros((M, 3), dtype=np.float64)
        slip_mask = np.zeros(M, dtype=bool)

        if act_idx.size == 0:
            res_m = Kf_m_csr @ um_f - f_ext_m_vec
            res_s = Kf_s_csr @ us_f - f_ext_s_vec
            res = np.concatenate([res_m, res_s])
            denom = max(1.0, float(np.linalg.norm(f_ext_m_vec)), float(np.linalg.norm(f_ext_s_vec)))
            rnorm = float(np.linalg.norm(res) / denom)
            return x_s, x_m, ab_m, B_m, g_n, n_m, active, act_idx, slip_mask, t_res, t_t_res, res, rnorm, None

        w_act = w_q[act_idx]
        g_act = g_n[act_idx]
        n_act = n_m[act_idx]
        nn = np.einsum("qi,qj->qij", n_act, n_act)
        P_t = I3[None, :, :] - nn
        t_n = (eps_n * g_act)[:, None] * n_act
        k_loc = (w_act * eps_n)[:, None, None] * nn if need_tangent else None

        B_m_act_scalar = B_m[act_idx, :]
        if mu > 0.0:
            du_m_act = B_m_act_scalar @ (U_m_curr - U_m_prev)
            du_rel = (x_s[act_idx] - x_s_prev[act_idx]) - du_m_act
            du_t = np.einsum("qij,qj->qi", P_t, du_rel)
            t_prev_proj = np.einsum("qij,qj->qi", P_t, t_t_prev[act_idx])
            t_trial = t_prev_proj + float(eps_t) * du_t
            tn = np.linalg.norm(t_trial, axis=1)
            cap = -mu * eps_n * g_act
            slip = tn > cap
            slip_mask[act_idx] = slip

            t_t_act = t_trial.copy()
            if need_tangent:
                k_t = (w_act * float(eps_t))[:, None, None] * P_t
            if np.any(slip):
                s_idx = np.flatnonzero(slip)
                d_hat = t_trial[s_idx] / tn[s_idx, None]
                t_t_act[s_idx] = cap[s_idx, None] * d_hat
                if need_tangent:
                    dd = np.einsum("qi,qj->qij", d_hat, d_hat)
                    dn = np.einsum("qi,qj->qij", d_hat, n_act[s_idx])
                    k_t[s_idx] = w_act[s_idx, None, None] * (
                        (float(eps_t) * (cap[s_idx] / tn[s_idx]))[:, None, None] * (P_t[s_idx] - dd)
                        - (mu * eps_n) * dn
                    )
            t_t_res[act_idx] = t_t_act
            t_res[act_idx] = w_act[:, None] * (t_n + t_t_act)
            if need_tangent:
                k_loc = k_loc + k_t
        else:
            t_res[act_idx] = w_act[:, None] * t_n

        rows = (3 * act_idx[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
        B_s_act = B_s_free[rows, :]
        B_m_act = expand_vector_basis(B_m_act_scalar, ndim=3)[:, free_m_arr].tocsr()

        t_act_flat = t_res[act_idx].ravel()
        f_c_m = -(B_m_act.T @ t_act_flat)
        f_c_s = +(B_s_act.T @ t_act_flat)
        res_m = Kf_m_csr @ um_f - f_ext_m_vec + f_c_m
        res_s = Kf_s_csr @ us_f - f_ext_s_vec + f_c_s
        res = np.concatenate([res_m, res_s])
        denom = max(
            1.0,
            float(np.linalg.norm(f_ext_m_vec)),
            float(np.linalg.norm(f_ext_s_vec)),
            float(np.linalg.norm(f_c_m)),
            float(np.linalg.norm(f_c_s)),
        )
        rnorm = float(np.linalg.norm(res) / denom)
        if not need_tangent:
            return x_s, x_m, ab_m, B_m, g_n, n_m, active, act_idx, slip_mask, t_res, t_t_res, res, rnorm, None
        B_pair = sp.hstack([-B_m_act, B_s_act], format="csr")
        return x_s, x_m, ab_m, B_m, g_n, n_m, active, act_idx, slip_mask, t_res, t_t_res, res, rnorm, (B_pair, k_loc)

    prev_active = None
    prev_slip = None
    converged = False
    n_solves = 0

    for _ in range(1, int(max_iter) + 1):
        (
            x_s, x_m, ab_m, B_m, g_n, n_m, active, act_idx, slip_mask,
            t_res, t_t_res, res, rnorm, tan_data
        ) = _evaluate(u)

        same_active = prev_active is not None and np.array_equal(active, prev_active)
        same_slip = prev_slip is not None and np.array_equal(slip_mask, prev_slip)
        if same_active and same_slip and rnorm <= rtol:
            converged = True
            break

        prev_active = active.copy()
        prev_slip = slip_mask.copy()

        if tan_data is None:
            K_tan = K_diag
        else:
            B_pair, k_loc = tan_data
            n_a = act_idx.size
            K_blk = sp.bsr_matrix(
                (k_loc, np.arange(n_a, dtype=np.int64), np.arange(n_a + 1, dtype=np.int64)),
                shape=(3 * n_a, 3 * n_a),
            ).tocsr()
            K_tan = K_diag + B_pair.T @ K_blk @ B_pair

        du = _solve_linear_step(K_tan, -res, use_pardiso=use_pardiso)
        n_solves += 1
        u_next = u + du

        if mu > 0.0 and float(np.linalg.norm(du)) > 0.0:
            res_norm0 = float(np.linalg.norm(res))
            alpha = 1.0
            while alpha >= 1.0 / 64.0:
                _, _, _, _, _, _, _, _, _, _, _, res_try, _, _ = _evaluate(u + alpha * du, need_tangent=False)
                if float(np.linalg.norm(res_try)) < (1.0 - 1e-4 * alpha) * res_norm0:
                    break
                alpha *= 0.5
            u_next = u + alpha * du

        u = u_next

    if not converged:
        raise RuntimeError(
            f"Deformable SSN contact did not converge in {max_iter} iterations "
            f"(rnorm={rnorm:.3e}, n_active={int(act_idx.size)})"
        )

    u_free_m = u[:n_free_m]
    u_free_s = u[n_free_m:]
    U_m_ctrl = _unpack_full_m(u_free_m)
    U_s_ctrl = _unpack_full_s(u_free_s)
    if act_idx.size > 0:
        B_m_vec_act = expand_vector_basis(B_m[act_idx, :], ndim=3)
        f_c_m_full = -(B_m_vec_act.T @ t_res[act_idx].ravel())
        rows_full = (3 * act_idx[:, None] + np.arange(3, dtype=np.int64)[None, :]).ravel()
        f_c_s_full = +(B_s_vec[rows_full, :].T @ t_res[act_idx].ravel())
    else:
        f_c_m_full = np.zeros(n_dofs_m, dtype=np.float64)
        f_c_s_full = np.zeros(n_dofs_s, dtype=np.float64)

    return {
        "u_free_m": u_free_m,
        "u_ctrl_m": U_m_ctrl,
        "u_free_s": u_free_s,
        "u_ctrl_s": U_s_ctrl,
        "active": active,
        "n_active": int(act_idx.size),
        "n_slip": int(slip_mask.sum()),
        "gap": g_n,
        "normal": n_m,
        "x_s": x_s,
        "x_m": x_m,
        "ab_m": ab_m,
        "B_m": B_m,
        "traction_res": t_res,
        "force_total": t_res.sum(axis=0),
        "f_c_m_full": f_c_m_full,
        "f_c_s_full": f_c_s_full,
        "normal_pressure_integral": float(np.sum(-w_q[act_idx] * eps_n * g_n[act_idx])) if act_idx.size else 0.0,
        "max_pressure": float(np.max(-eps_n * g_n[act_idx])) if act_idx.size else 0.0,
        "iterations": int(n_solves),
        "rnorm": float(rnorm),
        "slip_state": {
            "x_s_prev": x_s.copy(),
            "U_m_prev": U_m_ctrl.copy(),
            "t_t_prev": t_t_res.copy(),
        },
    }
