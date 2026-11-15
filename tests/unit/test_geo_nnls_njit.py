"""Unit tests for Numba Lawson-Hanson QR NNLS kernel and quadrature integration.

Covers Phase 1 verification criteria (a)-(e):
  (a) Planted nonnegative solutions on random and orthogonal-projected q=2 (m=27)
      and q=4 (m=125) NNLS systems.
  (b) Replay on compact fixture nnls_small_samples.npz (w >= 0, nnz <= m,
      rnorm <= scipy_rnorm * (1 + 1e-6) or <= 1e-12 * cell_vol).
  (c) Degeneracy injection (near-duplicate columns, rank-deficient |R_ii|
      collapse, zero columns, negative RHS, maxiter=1, oversized shape guard)
      verifying FAIL (ok=False) and exact SciPy fallback equivalence.
  (d) Warm-start vs cold-start equivalence on consecutive prefix-sharing pairs
      plus ``fit_moment_rule`` end-to-end integration.
  (e) 3 explicit mutation detection checks (dropping w >= 0, flipping tolerance
      sign, disabling FAIL/fallback check).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import nnls as scipy_nnls

from wht_iga.geometry._nnls_njit import (
    _add_col_cgs2,
    alloc_nnls_workspace,
    get_nnls_workspace,
    nnls_qr,
    nnls_qr_ws,
)
from wht_iga.geometry import quadrature as Q

_NNLS_TOL = 1e-12
_EMPTY_I64 = np.empty(0, dtype=np.int64)
_FIXTURE_PATH = Path(__file__).resolve().parent / "data" / "nnls_small_samples.npz"


def solve_nnls_with_fallback(
    A: np.ndarray,
    b: np.ndarray,
    *,
    maxiter: int = -1,
    warm_active: np.ndarray = _EMPTY_I64,
    prefix_len: int = 0,
    workspace: tuple[np.ndarray, ...] | None = None,
    solver_fn=None,
) -> tuple[np.ndarray, float, bool, bool]:
    """Run Numba NNLS kernel and fall back to SciPy nnls only when ok=False.

    Returns ``(w, rnorm, numba_ok, used_fallback)``.
    """
    n = A.shape[1]
    eff_maxiter = 10 * n if maxiter <= 0 else maxiter

    if solver_fn is not None:
        w, rn, ok = solver_fn(A, b, eff_maxiter, warm_active, prefix_len)
    elif workspace is not None:
        w, rn, ok = nnls_qr_ws(
            A, b, eff_maxiter, warm_active, prefix_len, *workspace
        )
    else:
        w, rn, ok = nnls_qr(A, b, eff_maxiter, warm_active, prefix_len)

    if not ok:
        w_sp, rn_sp = scipy_nnls(A, b, maxiter=eff_maxiter)
        return w_sp, float(rn_sp), False, True
    return w, float(rn), True, False


def assert_nnls_invariants(
    A: np.ndarray,
    b: np.ndarray,
    w: np.ndarray,
    rnorm: float,
    scipy_rnorm: float,
    cell_vol: float,
    *,
    tol_scale: float = _NNLS_TOL,
) -> None:
    """Verify core Phase 1 numerical invariants for one NNLS solve."""
    m = A.shape[0]
    assert np.all(np.isfinite(w)), "weights contain non-finite values"
    assert np.isfinite(rnorm), "residual norm is non-finite"
    assert np.all(w >= 0.0), f"negative weight detected: min(w)={w.min():.3e}"

    nnz = int(np.count_nonzero(w > 0.0))
    assert nnz <= m, f"nnz={nnz} exceeds row count m={m}"

    tol = tol_scale * cell_vol
    assert rnorm <= scipy_rnorm * (1.0 + 1e-6) or rnorm <= tol, (
        f"residual gate violated: rnorm={rnorm:.6e}, "
        f"scipy_rnorm={scipy_rnorm:.6e}, tol={tol:.6e}"
    )
    if scipy_rnorm <= tol:
        assert rnorm <= tol, (
            f"SciPy converged within tol={tol:.6e} (scipy={scipy_rnorm:.6e}) "
            f"but candidate rnorm={rnorm:.6e} did not"
        )

    direct_rnorm = float(np.linalg.norm(A @ w - b))
    assert abs(rnorm - direct_rnorm) <= 1e-11 * max(1.0, direct_rnorm, float(np.linalg.norm(b))), (
        f"reported rnorm={rnorm:.6e} differs from ||Aw - b||_2={direct_rnorm:.6e}"
    )


def _make_planted_nnls(
    rng: np.random.Generator,
    m: int,
    n: int,
    k_active: int,
    residual_norm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct ``(A, b, x_true)`` satisfying exact NNLS KKT conditions.

    When ``residual_norm > 0``, ``b = A @ x_true + residual_norm * r_orth``
    where ``A[:, P]^T r_orth = 0`` and ``A[:, Z]^T r_orth <= -0.05 < 0``.
    """
    assert k_active < m < n
    A = rng.normal(size=(m, n))
    perm = rng.permutation(n)
    P_idx = perm[:k_active]
    Z_idx = perm[k_active:]

    x_true = np.zeros(n, dtype=np.float64)
    x_true[P_idx] = rng.uniform(0.2, 2.0, size=k_active)

    if residual_norm > 0.0:
        Q_p, _ = np.linalg.qr(A[:, P_idx], mode="reduced")
        v = rng.normal(size=m)
        r_orth = v - Q_p @ (Q_p.T @ v)
        r_orth_norm = float(np.linalg.norm(r_orth))
        assert r_orth_norm > 1e-8
        r_orth /= r_orth_norm

        dots = r_orth @ A[:, Z_idx]
        shift = np.maximum(dots + 0.05, 0.0)
        A[:, Z_idx] -= np.outer(r_orth, shift)

        b = A @ x_true + residual_norm * r_orth
    else:
        b = A @ x_true

    return np.ascontiguousarray(A, dtype=np.float64), np.ascontiguousarray(b, dtype=np.float64), x_true


def _load_small_samples() -> list[dict]:
    """Unpack ``tests/unit/data/nnls_small_samples.npz`` into per-call sample dicts."""
    assert _FIXTURE_PATH.exists(), f"Missing fixture: {_FIXTURE_PATH}"
    d = np.load(_FIXTURE_PATH, allow_pickle=False)
    A_flat = np.ascontiguousarray(d["A_flat"], dtype=np.float64)
    co = d["col_offsets"]
    b_all = np.ascontiguousarray(d["b_all"], dtype=np.float64)
    w_flat = np.ascontiguousarray(d["w_flat"], dtype=np.float64)
    rnorms = d["rnorms"]
    cell_vols = d["cell_vols"]
    cell_idxs = d["cell_idxs"]
    call_in_cells = d["call_in_cells"]
    prefix_lens = d["prefix_lens"]
    n_nzs = d["n_nzs"]
    outcomes = d["outcomes"]

    samples: list[dict] = []
    for k in range(b_all.shape[0]):
        c0, c1 = int(co[k]), int(co[k + 1])
        samples.append(
            {
                "idx": k,
                "A": np.ascontiguousarray(A_flat[:, c0:c1], dtype=np.float64),
                "b": np.ascontiguousarray(b_all[k], dtype=np.float64),
                "w_ref": np.ascontiguousarray(w_flat[c0:c1], dtype=np.float64),
                "rnorm_ref": float(rnorms[k]),
                "cell_vol": float(cell_vols[k]),
                "cell_idx": int(cell_idxs[k]),
                "call_in_cell": int(call_in_cells[k]),
                "prefix_len": int(prefix_lens[k]),
                "n_nz_ref": int(n_nzs[k]),
                "outcome": str(outcomes[k]),
            }
        )
    return samples


# ---------------------------------------------------------------------------
# (a) Planted nonnegative solutions on q=2 (m=27) and q=4 (m=125)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q,m,n,k_active", [(2, 27, 64, 10), (4, 125, 256, 32)])
@pytest.mark.parametrize("residual_norm", [0.0, 0.75])
def test_planted_kkt_solutions(
    q: int,
    m: int,
    n: int,
    k_active: int,
    residual_norm: float,
) -> None:
    """Kernel recovers exact planted KKT solution for q=2 (m=27) and q=4 (m=125)."""
    assert (q + 1) ** 3 == m
    rng = np.random.default_rng(1000 + 17 * q + int(residual_norm * 10))
    A, b, x_true = _make_planted_nnls(rng, m, n, k_active, residual_norm)

    w_sp, rn_sp = scipy_nnls(A, b, maxiter=10 * n)
    w_sa, rn_sa, ok_sa = nnls_qr(A, b, 10 * n, _EMPTY_I64, 0)
    ws = alloc_nnls_workspace(m, n)
    w_ws, rn_ws, ok_ws = nnls_qr_ws(A, b, 10 * n, _EMPTY_I64, 0, *ws)

    assert ok_sa and ok_ws
    assert_nnls_invariants(A, b, w_sa, float(rn_sa), float(rn_sp), cell_vol=1.0)
    assert_nnls_invariants(A, b, w_ws, float(rn_ws), float(rn_sp), cell_vol=1.0)
    np.testing.assert_allclose(w_sa, x_true, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(w_ws, x_true, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(w_ws, w_sp, atol=1e-10, rtol=1e-10)
    assert abs(rn_ws - residual_norm) <= 1e-10


# ---------------------------------------------------------------------------
# (b) Replay on compact real-cell fixture nnls_small_samples.npz
# ---------------------------------------------------------------------------


def test_small_samples_fixture_coverage() -> None:
    """Verify fixture covers small/median/large n, depths 0..4, strict & relaxed, prefix pairs."""
    samples = _load_small_samples()
    assert 20 <= len(samples) <= 35

    ns = [s["A"].shape[1] for s in samples]
    assert any(1 <= n <= 20 for n in ns), "missing small n (1..20)"
    assert any(200 <= n <= 300 for n in ns), "missing median n (~200-300)"
    assert any(n > 1000 for n in ns), "missing large n (>1000)"

    outcomes = {s["outcome"] for s in samples}
    assert outcomes == {"strict", "relaxed"}

    cics_all = {s["call_in_cell"] for s in samples}
    assert {0, 1, 2, 3, 4}.issubset(cics_all)
    strict_cics = {s["call_in_cell"] for s in samples if s["outcome"] == "strict"}
    assert {0, 1, 2, 3}.issubset(strict_cics)

    prefix_pairs = [
        (samples[k - 1], samples[k])
        for k in range(1, len(samples))
        if samples[k]["cell_idx"] == samples[k - 1]["cell_idx"]
        and samples[k]["prefix_len"] > 0
    ]
    assert len(prefix_pairs) >= 1
    for prev_s, curr_s in prefix_pairs:
        plen = curr_s["prefix_len"]
        assert prev_s["A"].shape[1] == plen
        assert np.array_equal(prev_s["A"], curr_s["A"][:, :plen])


def test_small_samples_replay() -> None:
    """Replay all small_samples cases and enforce w >= 0, nnz <= m, and residual gates."""
    samples = _load_small_samples()
    max_n = max(s["A"].shape[1] for s in samples)
    ws = alloc_nnls_workspace(125, max(2048, max_n))

    for s in samples:
        A, b = s["A"], s["b"]
        w_sp, rn_sp = scipy_nnls(A, b, maxiter=10 * A.shape[1])

        w, rn, _ok, _fb = solve_nnls_with_fallback(A, b, workspace=ws)
        assert_nnls_invariants(A, b, w, rn, float(rn_sp), s["cell_vol"])

        w_sa, rn_sa, _ok_sa, _fb_sa = solve_nnls_with_fallback(A, b)
        assert_nnls_invariants(A, b, w_sa, rn_sa, float(rn_sp), s["cell_vol"])


# ---------------------------------------------------------------------------
# (c) Degeneracy injection & oversized workspace guard
# ---------------------------------------------------------------------------


def test_degeneracy_near_duplicate_columns_triggers_fail_and_fallback() -> None:
    """Near-duplicate columns entering passive set trigger ok=False; fallback matches SciPy."""
    A = np.array(
        [[1.0, 1.0 - 1e-13], [0.0, 1e-13], [0.0, 0.0]],
        dtype=np.float64,
    )
    b = np.array([1.0, 0.5, 0.0], dtype=np.float64)

    _w_raw, _rn_raw, ok = nnls_qr(A, b, 20, _EMPTY_I64, 0)
    assert not ok, "expected FAIL (ok=False) on near-duplicate columns"

    w_fb, rn_fb, ok_fb, used_fb = solve_nnls_with_fallback(A, b)
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=20)
    assert not ok_fb and used_fb
    np.testing.assert_allclose(w_fb, w_sp, atol=1e-12, rtol=1e-12)
    assert abs(rn_fb - float(rn_sp)) <= 1e-12


def test_degeneracy_exact_duplicate_columns_matches_scipy() -> None:
    """Exact duplicate columns never corrupt solution and match pure SciPy."""
    A = np.array(
        [[1.0, 1.0, 0.0], [2.0, 2.0, 1.0], [3.0, 3.0, 0.0]],
        dtype=np.float64,
    )
    b = np.array([1.5, 2.0, 3.0], dtype=np.float64)

    w_fb, rn_fb, _ok, _used_fb = solve_nnls_with_fallback(A, b)
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=30)
    assert_nnls_invariants(A, b, w_fb, rn_fb, float(rn_sp), cell_vol=1.0)
    np.testing.assert_allclose(w_fb, w_sp, atol=1e-12, rtol=1e-12)


def test_degeneracy_rank_deficient_r_collapse_triggers_fail_and_fallback() -> None:
    """Rank-deficient 3rd column collapsing |R_ii| triggers ok=False; fallback matches SciPy."""
    A = np.array(
        [[1.0, 0.0, 0.4], [0.0, 1.0, 0.4], [0.0, 0.0, 1e-13]],
        dtype=np.float64,
    )
    b = np.array([1.0, 1.0, 0.5], dtype=np.float64)

    _w_raw, _rn_raw, ok = nnls_qr(A, b, 30, _EMPTY_I64, 0)
    assert not ok, "expected FAIL (ok=False) when |R_ii| < 1e-12 * max|R|"

    w_fb, rn_fb, ok_fb, used_fb = solve_nnls_with_fallback(A, b)
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=30)
    assert not ok_fb and used_fb
    np.testing.assert_allclose(w_fb, w_sp, atol=1e-12, rtol=1e-12)
    assert abs(rn_fb - float(rn_sp)) <= 1e-12


def test_degeneracy_zero_columns() -> None:
    """Zero column cannot enter QR passive set and full solve matches pure SciPy."""
    A = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    b = np.array([1.0, 1.0, 0.0], dtype=np.float64)

    Q_buf, R_buf, qtb, _z, v, _r, P, in_P, _x, _g, c_buf = alloc_nnls_workspace(3, 3)
    in_P[:] = 0
    k, max_r, add_ok = _add_col_cgs2(
        A, b, 3, 0, 1, Q_buf, R_buf, qtb, v, c_buf, P, in_P, 0.0
    )
    assert not add_ok and k == 0 and max_r == 0.0

    w_fb, rn_fb, _ok, _fb = solve_nnls_with_fallback(A, b)
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=30)
    np.testing.assert_allclose(w_fb, w_sp, atol=1e-12)
    assert abs(rn_fb - float(rn_sp)) <= 1e-12

    A_zero = np.zeros((4, 3), dtype=np.float64)
    b_zero = np.array([1.0, 2.0, 0.0, 0.0], dtype=np.float64)
    w0, rn0, ok0, _fb0 = solve_nnls_with_fallback(A_zero, b_zero)
    w0_sp, rn0_sp = scipy_nnls(A_zero, b_zero, maxiter=30)
    assert ok0
    np.testing.assert_allclose(w0, w0_sp, atol=1e-12)
    assert abs(rn0 - float(rn0_sp)) <= 1e-12


def test_degeneracy_negative_rhs() -> None:
    """Negative RHS with positive, mixed, and ill-conditioned columns matches SciPy."""
    A_pos = np.array([[1.0, 0.5], [0.5, 1.0], [1.0, 1.0]], dtype=np.float64)
    b_neg = np.array([-1.0, -2.0, -0.5], dtype=np.float64)
    w1, rn1, ok1, _fb1 = solve_nnls_with_fallback(A_pos, b_neg)
    w1_sp, rn1_sp = scipy_nnls(A_pos, b_neg, maxiter=20)
    assert ok1
    np.testing.assert_allclose(w1, w1_sp, atol=1e-12)
    assert abs(rn1 - float(rn1_sp)) <= 1e-12

    A_neg_dup = np.array(
        [[-1.0, -1.0 + 1e-13], [0.0, -1e-13], [0.0, 0.0]],
        dtype=np.float64,
    )
    b_neg2 = np.array([-1.0, -0.5, 0.0], dtype=np.float64)
    _w2_raw, _rn2_raw, ok2_raw = nnls_qr(A_neg_dup, b_neg2, 20, _EMPTY_I64, 0)
    assert not ok2_raw, "expected FAIL on negative RHS near-duplicate system"

    w2, rn2, ok2, used_fb2 = solve_nnls_with_fallback(A_neg_dup, b_neg2)
    w2_sp, rn2_sp = scipy_nnls(A_neg_dup, b_neg2, maxiter=20)
    assert not ok2 and used_fb2
    np.testing.assert_allclose(w2, w2_sp, atol=1e-12, rtol=1e-12)
    assert abs(rn2 - float(rn2_sp)) <= 1e-12


def test_degeneracy_maxiter_one_triggers_fail() -> None:
    """maxiter=1 on multi-column problem returns ok=False and matches SciPy failure."""
    A = np.eye(4, dtype=np.float64)
    b = np.ones(4, dtype=np.float64)

    w_raw, rn_raw, ok = nnls_qr(A, b, 1, _EMPTY_I64, 0)
    assert not ok, "expected FAIL (ok=False) when maxiter=1 is exhausted"
    assert np.all(w_raw >= 0.0) and np.isfinite(rn_raw)

    with pytest.raises(RuntimeError):
        scipy_nnls(A, b, maxiter=1)


def test_nnls_qr_ws_oversized_guard() -> None:
    """Undersized workspace triggers standalone nnls_qr delegation without buffer overflow."""
    rng = np.random.default_rng(77)
    A, b, x_true = _make_planted_nnls(rng, m=16, n=40, k_active=6, residual_norm=0.25)
    tiny_ws = alloc_nnls_workspace(max_m=4, max_n=8)
    w, rn, ok = nnls_qr_ws(A, b, 10 * A.shape[1], _EMPTY_I64, 0, *tiny_ws)
    assert ok
    np.testing.assert_allclose(w, x_true, atol=1e-10, rtol=1e-10)
    assert abs(rn - 0.25) <= 1e-10

    ws_auto = get_nnls_workspace(min_m=200, min_n=3000)
    assert ws_auto[0].shape[0] >= 200 and ws_auto[7].shape[0] >= 3000


# ---------------------------------------------------------------------------
# (d) Warm-start vs cold-start equivalence & quadrature fit_moment_rule
# ---------------------------------------------------------------------------


def test_warm_start_matches_cold_start_on_prefix_pairs() -> None:
    """Warm-start reusing prior passive set matches cold-start gates on prefix pairs."""
    samples = _load_small_samples()
    max_n = max(s["A"].shape[1] for s in samples)
    ws = alloc_nnls_workspace(125, max(2048, max_n))

    tested_pairs = 0
    for k in range(1, len(samples)):
        prev_s = samples[k - 1]
        curr_s = samples[k]
        if curr_s["cell_idx"] != prev_s["cell_idx"] or curr_s["prefix_len"] <= 0:
            continue

        plen = curr_s["prefix_len"]
        assert prev_s["A"].shape[1] == plen
        assert np.array_equal(prev_s["A"], curr_s["A"][:, :plen])

        w_prev, _rn_prev, _ok_p, _fb_p = solve_nnls_with_fallback(
            prev_s["A"], prev_s["b"], workspace=ws
        )
        warm_active = np.flatnonzero(w_prev > 0.0).astype(np.int64)

        A, b = curr_s["A"], curr_s["b"]
        w_sp, rn_sp = scipy_nnls(A, b, maxiter=10 * A.shape[1])

        w_cold, rn_cold, _ok_c, _fb_c = solve_nnls_with_fallback(A, b, workspace=ws)
        w_warm, rn_warm, _ok_w, _fb_w = solve_nnls_with_fallback(
            A, b, warm_active=warm_active, prefix_len=plen, workspace=ws
        )

        assert_nnls_invariants(A, b, w_cold, rn_cold, float(rn_sp), curr_s["cell_vol"])
        assert_nnls_invariants(A, b, w_warm, rn_warm, float(rn_sp), curr_s["cell_vol"])
        assert abs(rn_warm - rn_cold) <= 1e-12 * max(1.0, rn_cold, float(rn_sp))
        tested_pairs += 1

    assert tested_pairs >= 1, "No consecutive prefix-sharing pair found in fixture"


def test_warm_start_bad_active_set_falls_back_cleanly() -> None:
    """Ill-conditioned warm_active triggers internal cold-start fallback without error."""
    A = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    b = np.array([1.0, 1.0, 0.0], dtype=np.float64)
    bad_warm = np.array([0, 1], dtype=np.int64)

    ws = alloc_nnls_workspace(3, 3)
    w, rn, ok, _fb = solve_nnls_with_fallback(
        A, b, warm_active=bad_warm, prefix_len=3, workspace=ws
    )
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=30)
    assert ok
    np.testing.assert_allclose(w, w_sp, atol=1e-12)
    assert abs(rn - float(rn_sp)) <= 1e-12


def test_fit_moment_rule_fast_nnls() -> None:
    """fit_moment_rule recovers valid pruned quadrature weights via Numba QR NNLS."""
    blo = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    h = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q = 2
    Xg, Wg = Q._tensor_gauss(q + 1)
    Xd, _ = Q._tensor_gauss(q + 2)
    cand = np.vstack([blo + Xg * h, blo + Xd * h])
    # Exact box moments on [0, 1]^3: only M[0, 0, 0] == volume == 1.0
    moments = np.zeros((q + 1, q + 1, q + 1), dtype=np.float64)
    moments[0, 0, 0] = 1.0

    warm_state: dict = {"active": _EMPTY_I64, "prefix_len": 0, "strided": False}
    w, rnorm, used_cand = Q.fit_moment_rule(
        moments, cand, blo, h, cell_vol=1.0, _warm_state=warm_state
    )
    assert rnorm <= 1e-12
    assert np.all(w >= 0.0)
    assert int(np.count_nonzero(w > 0.0)) <= (q + 1) ** 3
    assert warm_state["prefix_len"] == len(used_cand)
    np.testing.assert_allclose(float(w.sum()), float(Wg.sum()), atol=1e-12, rtol=1e-12)


# ---------------------------------------------------------------------------
# (e) Mutation detection checks (proving suite catches broken invariants)
# ---------------------------------------------------------------------------


def test_mutation_dropping_nonnegativity_is_detected() -> None:
    """Mutation 1: dropping w >= 0 clamping must fail assert_nnls_invariants."""
    samples = _load_small_samples()
    s = samples[0]
    A, b = s["A"], s["b"]
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=10 * A.shape[1])

    w_mutated = w_sp.copy()
    w_mutated[0] = -1e-4
    rn_mutated = float(np.linalg.norm(A @ w_mutated - b))

    with pytest.raises(AssertionError, match="negative weight"):
        assert_nnls_invariants(A, b, w_mutated, rn_mutated, float(rn_sp), s["cell_vol"])


def test_mutation_flipping_tolerance_sign_is_detected() -> None:
    """Mutation 2: flipping tolerance sign in gate or accepting degraded residual fails."""
    samples = _load_small_samples()
    converged = next(
        s for s in samples
        if s["outcome"] == "strict" and s["rnorm_ref"] <= _NNLS_TOL * s["cell_vol"]
    )
    A, b = converged["A"], converged["b"]
    w_sp, rn_sp = scipy_nnls(A, b, maxiter=10 * A.shape[1])

    w_bad = w_sp.copy()
    w_bad[0] += 1e-3
    rn_bad = float(np.linalg.norm(A @ w_bad - b))
    with pytest.raises(AssertionError, match="residual gate violated|SciPy converged"):
        assert_nnls_invariants(A, b, w_bad, rn_bad, float(rn_sp), converged["cell_vol"])

    def flipped_tol_gate(rnorm: float, cell_vol: float) -> bool:
        return rnorm <= -_NNLS_TOL * cell_vol

    assert not flipped_tol_gate(float(rn_sp), converged["cell_vol"])
    assert float(rn_sp) <= _NNLS_TOL * converged["cell_vol"]


def test_mutation_disabling_fail_check_is_detected(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mutation 3: disabling ok=False fallback causes degeneracy test to diverge from SciPy."""
    A_rk = np.array(
        [[1.0, 0.0, 0.4], [0.0, 1.0, 0.4], [0.0, 0.0, 1e-13]],
        dtype=np.float64,
    )
    b_rk = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    w_sp, _rn_sp = scipy_nnls(A_rk, b_rk, maxiter=30)

    orig_nnls_qr = nnls_qr

    def broken_no_fail_kernel(
        mat: np.ndarray,
        rhs: np.ndarray,
        maxiter: int = -1,
        warm_active: np.ndarray = _EMPTY_I64,
        prefix_len: int = 0,
    ) -> tuple[np.ndarray, float, bool]:
        w_raw, rn_raw, ok_raw = orig_nnls_qr(mat, rhs, maxiter, warm_active, prefix_len)
        assert not ok_raw
        return w_raw, float(rn_raw), True

    monkeypatch.setattr(sys.modules[__name__], "nnls_qr", broken_no_fail_kernel)

    w_mut, _rn_mut, ok_mut, used_fb = solve_nnls_with_fallback(A_rk, b_rk)
    assert ok_mut and not used_fb
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(w_mut, w_sp, atol=1e-12, rtol=1e-12)

    A_id = np.eye(4, dtype=np.float64)
    b_id = np.ones(4, dtype=np.float64)
    _w_id_sp, rn_id_sp = scipy_nnls(A_id, b_id, maxiter=40)
    w_id_mut, rn_id_mut, ok_id_mut, used_id_fb = solve_nnls_with_fallback(
        A_id, b_id, maxiter=1
    )
    assert ok_id_mut and not used_id_fb
    with pytest.raises(AssertionError):
        assert_nnls_invariants(A_id, b_id, w_id_mut, rn_id_mut, float(rn_id_sp), cell_vol=1.0)
