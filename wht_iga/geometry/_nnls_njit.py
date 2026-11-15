"""Lawson-Hanson active-set NNLS using thin incremental QR in Numba ``@njit``.

Implements genuine Lawson-Hanson NNLS with:
- Thin QR factorization ``A_P = Q^T R`` (``Q`` stored row-wise for contiguous access)
- Two-pass Classical Gram-Schmidt with re-orthogonalization (CGS2)
- Givens rotation QR downdate on column drop
- Two-pass ``Q``-projected residual inside the active-set loop plus exact
  ``b - A_P x_P`` residual on exit
- Optional warm-start passive-set reuse (``prefix_len > 0``)
- Explicit ``ok=False`` signaling on stagnation, numerical rank loss, or iteration
  exhaustion so callers can fall back deterministically to ``scipy.optimize.nnls``
"""

from __future__ import annotations

import math
import threading

import numpy as np
from numba import njit  # type: ignore[import-untyped]

_EMPTY_I64 = np.empty(0, dtype=np.int64)
_EPS = 2.220446049250313e-16
_RCOND_TOL = 1e-12
_FASTMATH = {"reassoc", "contract"}

_WS_LOCAL = threading.local()


def alloc_nnls_workspace(max_m: int = 125, max_n: int = 2048) -> tuple[np.ndarray, ...]:
    """Allocate reusable workspace arrays for :func:`nnls_qr_ws`."""
    Q = np.empty((max_m, max_m), dtype=np.float64)
    R = np.empty((max_m, max_m), dtype=np.float64)
    qtb = np.empty(max_m, dtype=np.float64)
    z = np.empty(max_m, dtype=np.float64)
    v = np.empty(max_m, dtype=np.float64)
    r = np.empty(max_m, dtype=np.float64)
    P = np.empty(max_m, dtype=np.int64)
    in_P = np.empty(max_n, dtype=np.uint8)
    x = np.empty(max_n, dtype=np.float64)
    g = np.empty(max_n, dtype=np.float64)
    c_buf = np.empty(max_m, dtype=np.float64)
    return (Q, R, qtb, z, v, r, P, in_P, x, g, c_buf)


def get_nnls_workspace(min_m: int = 125, min_n: int = 2048) -> tuple[np.ndarray, ...]:
    """Return a per-thread reusable NNLS workspace large enough for ``(min_m, min_n)``."""
    ws = getattr(_WS_LOCAL, "ws", None)
    cur_m = ws[0].shape[0] if ws is not None else 0
    cur_n = ws[7].shape[0] if ws is not None else 0
    if ws is None or cur_m < min_m or cur_n < min_n:
        ws = alloc_nnls_workspace(max(125, int(min_m), cur_m), max(2048, int(min_n), cur_n))
        _WS_LOCAL.ws = ws
    return ws


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _add_col_cgs2(
    A: np.ndarray,
    b: np.ndarray,
    m: int,
    k: int,
    t: int,
    Q: np.ndarray,
    R: np.ndarray,
    qtb: np.ndarray,
    v: np.ndarray,
    c_buf: np.ndarray,
    P: np.ndarray,
    in_P: np.ndarray,
    max_diag_R: float,
) -> tuple[int, float, bool]:
    """Append column ``t`` of ``A`` to thin QR ``(Q, R)`` via 2-pass Gram-Schmidt."""
    if k >= m:
        return k, max_diag_R, False

    v_norm2 = 0.0
    for row in range(m):
        val = A[row, t]
        v[row] = val
        v_norm2 += val * val
    if v_norm2 <= 0.0:
        return k, max_diag_R, False

    if k > 0:
        for i in range(k):
            qi = Q[i]
            dot = 0.0
            for row in range(m):
                dot += qi[row] * v[row]
            c_buf[i] = dot
        for i in range(k):
            qi = Q[i]
            ci = c_buf[i]
            R[i, k] = ci
            for row in range(m):
                v[row] -= ci * qi[row]

        for i in range(k):
            qi = Q[i]
            dot = 0.0
            for row in range(m):
                dot += qi[row] * v[row]
            c_buf[i] = dot
        for i in range(k):
            qi = Q[i]
            ci = c_buf[i]
            R[i, k] += ci
            for row in range(m):
                v[row] -= ci * qi[row]

    r_kk2 = 0.0
    for row in range(m):
        val = v[row]
        r_kk2 += val * val
    r_kk = math.sqrt(r_kk2)

    new_max = r_kk if r_kk > max_diag_R else max_diag_R
    if new_max <= 0.0 or r_kk < _RCOND_TOL * new_max:
        return k, new_max, False

    if new_max > max_diag_R:
        thresh = _RCOND_TOL * new_max
        for i in range(k):
            if abs(R[i, i]) < thresh:
                return k, new_max, False

    inv_r = 1.0 / r_kk
    qk = Q[k]
    qtb_k = 0.0
    for row in range(m):
        qval = v[row] * inv_r
        qk[row] = qval
        qtb_k += qval * b[row]

    R[k, k] = r_kk
    qtb[k] = qtb_k
    P[k] = t
    in_P[t] = 1
    return k + 1, new_max, True


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _solve_upper(
    R: np.ndarray,
    qtb: np.ndarray,
    z: np.ndarray,
    k: int,
    max_diag_R: float,
) -> bool:
    """Solve upper-triangular system ``R[:k, :k] z[:k] = qtb[:k]`` by back-substitution."""
    thresh = _RCOND_TOL * max_diag_R
    for i in range(k - 1, -1, -1):
        rii = R[i, i]
        if rii == 0.0 or abs(rii) < thresh:
            return False
        s = qtb[i]
        ri = R[i]
        for j in range(i + 1, k):
            s -= ri[j] * z[j]
        z[i] = s / rii
    return True


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _drop_col_givens(
    pos: int,
    m: int,
    k: int,
    Q: np.ndarray,
    R: np.ndarray,
    qtb: np.ndarray,
    P: np.ndarray,
    in_P: np.ndarray,
) -> tuple[int, float]:
    """Remove column at position ``pos`` from passive set and restore QR via Givens rotations."""
    in_P[P[pos]] = 0
    for idx in range(pos, k - 1):
        P[idx] = P[idx + 1]

    for j in range(pos, k - 1):
        for i in range(j + 2):
            R[i, j] = R[i, j + 1]

    for i in range(pos, k - 1):
        a = R[i, i]
        bv = R[i + 1, i]
        if bv != 0.0:
            if abs(bv) > abs(a):
                tau = a / bv
                s = 1.0 / math.sqrt(1.0 + tau * tau)
                if bv < 0.0:
                    s = -s
                c = s * tau
            else:
                tau = bv / a
                c = 1.0 / math.sqrt(1.0 + tau * tau)
                if a < 0.0:
                    c = -c
                s = c * tau

            R[i, i] = c * a + s * bv
            R[i + 1, i] = 0.0

            ri = R[i]
            ri1 = R[i + 1]
            for j in range(i + 1, k - 1):
                t1 = ri[j]
                t2 = ri1[j]
                ri[j] = c * t1 + s * t2
                ri1[j] = -s * t1 + c * t2

            qi = Q[i]
            qi1 = Q[i + 1]
            for col in range(m):
                q1 = qi[col]
                q2 = qi1[col]
                qi[col] = c * q1 + s * q2
                qi1[col] = -s * q1 + c * q2

            qb1 = qtb[i]
            qb2 = qtb[i + 1]
            qtb[i] = c * qb1 + s * qb2
            qtb[i + 1] = -s * qb1 + c * qb2

    new_k = k - 1
    max_diag = 0.0
    for i in range(new_k):
        d = abs(R[i, i])
        if d > max_diag:
            max_diag = d
    return new_k, max_diag


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _compute_qr_residual_2pass(
    b: np.ndarray,
    m: int,
    k: int,
    Q: np.ndarray,
    qtb: np.ndarray,
    c_buf: np.ndarray,
    r: np.ndarray,
) -> float:
    """Compute ``r = (I - Q^T Q)^2 b`` against ``Q[:k]`` and return ``||r||_2``."""
    for row in range(m):
        r[row] = b[row]
    for i in range(k):
        qi = Q[i]
        qbi = qtb[i]
        for row in range(m):
            r[row] -= qbi * qi[row]
    if k > 0:
        for i in range(k):
            qi = Q[i]
            dot = 0.0
            for row in range(m):
                dot += qi[row] * r[row]
            c_buf[i] = dot
        for i in range(k):
            qi = Q[i]
            ci = c_buf[i]
            for row in range(m):
                r[row] -= ci * qi[row]
    norm2 = 0.0
    for row in range(m):
        val = r[row]
        norm2 += val * val
    return math.sqrt(norm2)


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _compute_exact_residual(
    A: np.ndarray,
    b: np.ndarray,
    m: int,
    P: np.ndarray,
    x: np.ndarray,
    k: int,
    r: np.ndarray,
) -> float:
    """Compute exact ``r = b - A_P x_P`` and return ``||r||_2``."""
    for i in range(m):
        r[i] = b[i]
    for idx in range(k):
        col = P[idx]
        xc = x[col]
        if xc != 0.0:
            for i in range(m):
                r[i] -= A[i, col] * xc
    norm2 = 0.0
    for i in range(m):
        ri = r[i]
        norm2 += ri * ri
    return math.sqrt(norm2)


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _compute_grad_argmax(
    A: np.ndarray,
    r: np.ndarray,
    g: np.ndarray,
    in_P: np.ndarray,
    m: int,
    n: int,
) -> tuple[int, float]:
    """Compute ``g = A^T r`` and return ``(argmax_{j not in P} g[j], g_max)``."""
    for j in range(n):
        g[j] = 0.0

    i = 0
    while i + 4 <= m:
        r0 = r[i]
        r1 = r[i + 1]
        r2 = r[i + 2]
        r3 = r[i + 3]
        a0 = A[i]
        a1 = A[i + 1]
        a2 = A[i + 2]
        a3 = A[i + 3]
        for j in range(n):
            g[j] += r0 * a0[j] + r1 * a1[j] + r2 * a2[j] + r3 * a3[j]
        i += 4

    while i < m:
        ri = r[i]
        ai = A[i]
        for j in range(n):
            g[j] += ri * ai[j]
        i += 1

    t = -1
    g_max = 0.0
    for j in range(n):
        if in_P[j] == 0:
            gj = g[j]
            if t == -1 or gj > g_max:
                g_max = gj
                t = j
    return t, g_max


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def _nnls_qr_core(
    A: np.ndarray,
    b: np.ndarray,
    maxiter: int,
    warm_active: np.ndarray,
    prefix_len: int,
    Q: np.ndarray,
    R: np.ndarray,
    qtb: np.ndarray,
    z: np.ndarray,
    v: np.ndarray,
    r: np.ndarray,
    P: np.ndarray,
    in_P: np.ndarray,
    x: np.ndarray,
    g: np.ndarray,
    c_buf: np.ndarray,
) -> tuple[np.ndarray, float, bool]:
    """Core Lawson-Hanson QR NNLS loop using caller-provided workspace buffers."""
    m = A.shape[0]
    n = A.shape[1]

    if maxiter <= 0:
        maxiter = 10 * n

    b_norm2 = 0.0
    for i in range(m):
        bi = b[i]
        b_norm2 += bi * bi
    b_norm = math.sqrt(b_norm2)

    w_out = np.zeros(n, dtype=np.float64)
    if m == 0 or n == 0 or b_norm == 0.0:
        return w_out, b_norm, True

    for j in range(n):
        in_P[j] = 0
        x[j] = 0.0

    k = 0
    max_diag_R = 0.0
    iters = 0
    obj = b_norm
    for i in range(m):
        r[i] = b[i]

    n_warm = warm_active.shape[0]
    if prefix_len > 0 and n_warm > 0:
        eff_pref = prefix_len if prefix_len < n else n
        warm_ok = True
        for wi in range(n_warm):
            j = warm_active[wi]
            if 0 <= j < eff_pref and in_P[j] == 0:
                k, max_diag_R, add_ok = _add_col_cgs2(
                    A, b, m, k, j, Q, R, qtb, v, c_buf, P, in_P, max_diag_R
                )
                if not add_ok:
                    warm_ok = False
                    break

        if warm_ok and k > 0:
            if not _solve_upper(R, qtb, z, k, max_diag_R):
                warm_ok = False

        if warm_ok and k > 0:
            while k > 0:
                all_pos = True
                for i in range(k):
                    if z[i] <= 0.0:
                        all_pos = False
                        break
                if all_pos:
                    break

                iters += 1
                if iters > maxiter:
                    warm_ok = False
                    break

                alpha = 2.0
                min_pos = -1
                for i in range(k):
                    if z[i] <= 0.0:
                        xp = x[P[i]]
                        denom = xp - z[i]
                        ratio = xp / denom if denom > 0.0 else 0.0
                        if ratio < alpha:
                            alpha = ratio
                            min_pos = i
                if alpha > 1.0:
                    alpha = 1.0

                for i in range(k):
                    col = P[i]
                    val = x[col] + alpha * (z[i] - x[col])
                    x[col] = val if val > 0.0 else 0.0
                if min_pos >= 0:
                    x[P[min_pos]] = 0.0

                for i in range(k - 1, -1, -1):
                    col = P[i]
                    if x[col] <= 0.0:
                        x[col] = 0.0
                        k, max_diag_R = _drop_col_givens(i, m, k, Q, R, qtb, P, in_P)

                if k > 0:
                    if not _solve_upper(R, qtb, z, k, max_diag_R):
                        warm_ok = False
                        break

        if warm_ok and k > 0:
            for i in range(k):
                x[P[i]] = z[i]
            obj = _compute_qr_residual_2pass(b, m, k, Q, qtb, c_buf, r)
        elif not warm_ok:
            for j in range(n):
                in_P[j] = 0
                x[j] = 0.0
            k = 0
            max_diag_R = 0.0
            iters = 0
            obj = b_norm
            for i in range(m):
                r[i] = b[i]

    obj_tol = 1e-14 * b_norm
    kkt_scale = 8.0 * _EPS
    while True:
        if k == m or obj <= obj_tol:
            break

        t, g_max = _compute_grad_argmax(A, r, g, in_P, m, n)
        if t == -1 or g_max <= 0.0:
            break

        col_norm2 = 0.0
        for i in range(m):
            val = A[i, t]
            col_norm2 += val * val
        if g_max <= kkt_scale * math.sqrt(col_norm2) * obj:
            break

        iters += 1
        if iters > maxiter:
            obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
            for i in range(k):
                w_out[P[i]] = x[P[i]]
            return w_out, obj_ex, False

        k, max_diag_R, add_ok = _add_col_cgs2(
            A, b, m, k, t, Q, R, qtb, v, c_buf, P, in_P, max_diag_R
        )
        if not add_ok:
            obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
            for i in range(k):
                w_out[P[i]] = x[P[i]]
            return w_out, obj_ex, False

        if not _solve_upper(R, qtb, z, k, max_diag_R):
            obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
            for i in range(k):
                w_out[P[i]] = x[P[i]]
            return w_out, obj_ex, False

        while True:
            all_pos = True
            for i in range(k):
                if z[i] <= 0.0:
                    all_pos = False
                    break
            if all_pos:
                break

            iters += 1
            if iters > maxiter:
                obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
                for i in range(k):
                    w_out[P[i]] = x[P[i]]
                return w_out, obj_ex, False

            alpha = 2.0
            min_pos = -1
            for i in range(k):
                if z[i] <= 0.0:
                    xp = x[P[i]]
                    denom = xp - z[i]
                    ratio = xp / denom if denom > 0.0 else 0.0
                    if ratio < alpha:
                        alpha = ratio
                        min_pos = i
            if alpha > 1.0:
                alpha = 1.0
            if alpha <= 0.0:
                obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
                for i in range(k):
                    w_out[P[i]] = x[P[i]]
                return w_out, obj_ex, False

            for i in range(k):
                col = P[i]
                val = x[col] + alpha * (z[i] - x[col])
                x[col] = val if val > 0.0 else 0.0
            if min_pos >= 0:
                x[P[min_pos]] = 0.0

            dropped_any = False
            for i in range(k - 1, -1, -1):
                col = P[i]
                if x[col] <= 0.0:
                    x[col] = 0.0
                    k, max_diag_R = _drop_col_givens(i, m, k, Q, R, qtb, P, in_P)
                    dropped_any = True

            if not dropped_any:
                obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
                for i in range(k):
                    w_out[P[i]] = x[P[i]]
                return w_out, obj_ex, False

            if k == 0:
                break

            if not _solve_upper(R, qtb, z, k, max_diag_R):
                obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
                for i in range(k):
                    w_out[P[i]] = x[P[i]]
                return w_out, obj_ex, False

        for i in range(k):
            x[P[i]] = z[i]

        obj_new = _compute_qr_residual_2pass(b, m, k, Q, qtb, c_buf, r)
        if obj_new >= obj:
            obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
            for i in range(k):
                w_out[P[i]] = x[P[i]]
            return w_out, obj_ex, False
        obj = obj_new

    if k > m:
        obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
        return w_out, obj_ex, False

    for i in range(k):
        col = P[i]
        val = x[col]
        if val > 0.0:
            w_out[col] = val

    obj_ex = _compute_exact_residual(A, b, m, P, x, k, r)
    return w_out, obj_ex, True


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def nnls_qr_ws(
    A: np.ndarray,
    b: np.ndarray,
    maxiter: int,
    warm_active: np.ndarray,
    prefix_len: int,
    Q: np.ndarray,
    R: np.ndarray,
    qtb: np.ndarray,
    z: np.ndarray,
    v: np.ndarray,
    r: np.ndarray,
    P: np.ndarray,
    in_P: np.ndarray,
    x: np.ndarray,
    g: np.ndarray,
    c_buf: np.ndarray,
) -> tuple[np.ndarray, float, bool]:
    """Solve ``argmin_{x >= 0} ||A x - b||_2`` via Lawson-Hanson QR with preallocated buffers."""
    m = A.shape[0]
    n = A.shape[1]
    if m <= Q.shape[0] and n <= in_P.shape[0]:
        return _nnls_qr_core(
            A, b, maxiter, warm_active, prefix_len,
            Q, R, qtb, z, v, r, P, in_P, x, g, c_buf,
        )
    return _nnls_qr_core(
        A, b, maxiter, warm_active, prefix_len,
        np.empty((m, m), dtype=np.float64),
        np.empty((m, m), dtype=np.float64),
        np.empty(m, dtype=np.float64),
        np.empty(m, dtype=np.float64),
        np.empty(m, dtype=np.float64),
        np.empty(m, dtype=np.float64),
        np.empty(m, dtype=np.int64),
        np.empty(n, dtype=np.uint8),
        np.empty(n, dtype=np.float64),
        np.empty(n, dtype=np.float64),
        np.empty(m, dtype=np.float64),
    )


@njit(cache=True, nogil=True, fastmath=_FASTMATH)
def nnls_qr(
    A: np.ndarray,
    b: np.ndarray,
    maxiter: int = -1,
    warm_active: np.ndarray = _EMPTY_I64,
    prefix_len: int = 0,
) -> tuple[np.ndarray, float, bool]:
    """Standalone ``@njit(cache=True, nogil=True)`` Lawson-Hanson QR NNLS solver."""
    m = A.shape[0]
    n = A.shape[1]
    Q = np.empty((m, m), dtype=np.float64)
    R = np.empty((m, m), dtype=np.float64)
    qtb = np.empty(m, dtype=np.float64)
    z = np.empty(m, dtype=np.float64)
    v = np.empty(m, dtype=np.float64)
    r = np.empty(m, dtype=np.float64)
    P = np.empty(m, dtype=np.int64)
    in_P = np.empty(n, dtype=np.uint8)
    x = np.empty(n, dtype=np.float64)
    g = np.empty(n, dtype=np.float64)
    c_buf = np.empty(m, dtype=np.float64)
    return _nnls_qr_core(
        A, b, maxiter, warm_active, prefix_len,
        Q, R, qtb, z, v, r, P, in_P, x, g, c_buf,
    )
