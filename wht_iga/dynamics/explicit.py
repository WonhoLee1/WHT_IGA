# -*- coding: utf-8 -*-
"""Explicit dynamics helpers (Phase 6 Alpha).

Kratos ``MechanicalExplicitStrategy`` drives the time loop; this module owns
the immersed-IGA specifics that Kratos cannot do alone on trimmed grids:

- HRZ (Hinton-Rock-Zienkiewicz) diagonal lumping, including of the reduced
  ``T^T M_c T`` system (variational lumped mass for extension-stabilized cells).
- Polynomial x-extension operator ``T`` eliminating the last x-layer control
  points (cut-cell stabilization restoring ``dt_crit ~ h/c``).
- ``NODAL_MASS`` overwrite helpers for the hybrid path: after Kratos
  ``Initialize()``, master-node masses are replaced by HRZ values so that
  negative extension weights do not produce zero/negative lumped masses.
- Small-model ``dt_crit`` estimator (generalized eigenproblem, SciPy).

Ownership: Track 2 (Phase 6 Explicit). New file; no existing module modified.
"""

from __future__ import annotations

import numpy as np


def hrz_lump_diagonal(M_consistent: np.ndarray, total_mass: float) -> np.ndarray:
    """Hinton-Rock-Zienkiewicz diagonal scaling, per spatial component.

    Row-sum the consistent mass diagonal per component (x, y, z) and rescale
    each component so its sum equals ``total_mass``. Raises ``ValueError`` on
    non-positive component sums (rather than silently producing bad masses).
    """
    diag = np.diag(np.asarray(M_consistent, dtype=float)).copy()
    if diag.shape[0] % 3 != 0:
        raise ValueError(f"mass DOFs {diag.shape[0]} not a multiple of 3")
    out = diag.copy()
    for d in range(3):
        s = float(np.sum(diag[d::3]))
        if not s > 0.0:
            raise ValueError(f"HRZ component-{d} diagonal sum must be positive, got {s}")
        out[d::3] *= total_mass / s
    return out


def _lagrange_weights_1d(x_eval: float, x_nodes: np.ndarray) -> np.ndarray:
    """1D Lagrange interpolation/extrapolation weights at ``x_eval``."""
    x_nodes = np.asarray(x_nodes, dtype=float)
    n = len(x_nodes)
    w = np.ones(n, dtype=float)
    for j in range(n):
        for k in range(n):
            if j != k:
                w[j] *= (x_eval - x_nodes[k]) / (x_nodes[j] - x_nodes[k])
    return w


def build_x_extension(box) -> dict:
    """Polynomial extension eliminating the last x-layer control points.

    Returns dict with ``T`` (u_full = T @ u_master), ``n_good_x``,
    ``slave_cp`` (eliminated control-point indices), and per-slave
    ``(master_idx, weights)`` in x (identical for every (y, z) row).
    Masters are the ``p + 1`` nearest interior Greville abscissae.
    """
    p = int(box.degree)
    nu_cp, nv_cp, nw_cp = (int(v) for v in box.n_cp)
    kv = np.asarray(box.knots[0], dtype=float)
    gx = np.array([float(np.mean(kv[i + 1 : i + 1 + p])) for i in range(nu_cp)])
    n_good_x = nu_cp - 1
    Tx = np.zeros((nu_cp, n_good_x), dtype=float)
    Tx[:n_good_x, :n_good_x] = np.eye(n_good_x)
    src = np.arange(n_good_x - (p + 1), n_good_x)
    weights = _lagrange_weights_1d(float(gx[nu_cp - 1]), gx[src])
    Tx[nu_cp - 1, src] = weights
    T_scalar = np.kron(np.eye(nw_cp * nv_cp), Tx)
    T = np.kron(T_scalar, np.eye(3))
    slave_cp = [nu_cp - 1 + nu_cp * (j + nv_cp * k)
                for k in range(nw_cp) for j in range(nv_cp)]
    masters = [(int(s), np.asarray(weights, dtype=float)) for s in src]
    return {"T": T, "n_good_x": n_good_x, "slave_cp": slave_cp,
            "masters": masters, "src": np.asarray(src, dtype=int)}


def read_nodal_mass(KM, mp, node_ids: list[int]) -> np.ndarray:
    """Read non-historical ``NODAL_MASS`` for the given node ids (in order)."""
    nodes = {nd.Id: nd for nd in mp.Nodes}
    return np.array([float(nodes[i].GetValue(KM.NODAL_MASS)) for i in node_ids],
                    dtype=float)


def overwrite_nodal_mass(KM, mp, masses_by_id: dict[int, float]) -> dict[int, float]:
    """Overwrite non-historical ``NODAL_MASS``; returns previous values.

    Intended for the hybrid path: call after ``strategy.Initialize()`` so the
    Python-computed HRZ values survive (Kratos computes its own masses once in
    ``Initialize`` and does not refresh them every step — verified by probe).
    """
    nodes = {nd.Id: nd for nd in mp.Nodes}
    prev = {i: float(nodes[i].GetValue(KM.NODAL_MASS)) for i in masses_by_id}
    for i, m in masses_by_id.items():
        nodes[i].SetValue(KM.NODAL_MASS, float(m))
    return prev


def estimate_dt_crit(K: np.ndarray, m_diag: np.ndarray,
                     free_mask: np.ndarray) -> tuple[float, float]:
    """``dt_crit = 2 / wmax`` from the generalized eigenproblem (small models).

    Solves ``K v = w^2 M v`` on free DOFs with diagonal ``M`` via SciPy
    ``eigh``. Returns ``(dt_crit, wmax)``. Dense solver: small models only.
    """
    from scipy.linalg import eigh

    Kf = np.asarray(K, dtype=float)[np.ix_(free_mask, free_mask)]
    mf = np.asarray(m_diag, dtype=float)[free_mask]
    if np.any(mf <= 0.0):
        raise ValueError("non-positive lumped mass in free set; explicit unstable")
    vals = eigh(Kf, np.diag(mf), eigvals_only=True)
    wmax = float(np.sqrt(max(float(np.max(vals)), 0.0)))
    if wmax <= 0.0:
        raise ValueError("non-positive max eigenfrequency")
    return 2.0 / wmax, wmax


def floor_lumped_mass(m_diag: np.ndarray, m_min: float) -> np.ndarray:
    """Selective floor for lumped masses (mass scaling of the stiff tail only).

    Entries below ``m_min`` are raised to ``m_min``; the rest untouched.
    Returns a new array.
    """
    m = np.asarray(m_diag, dtype=float).copy()
    m[m < m_min] = float(m_min)
    return m


def run_explicit_history(KM, mp, strategy, dt: float, n_steps: int,
                         probe_fn=None, step_hook=None) -> dict:
    """Drive a Kratos explicit strategy with absolute-time stepping.

    ``CloneTimeStep`` takes absolute TIME (not an increment); this helper owns
    the ``t += dt`` bookkeeping (a past spike froze by passing ``dt`` every step).
    ``probe_fn(mp) -> float`` records one scalar per step; ``step_hook(mp, i)``
    runs before each solve (e.g. penalty updates). Returns history dict.
    Raises on non-finite probe values.
    """
    hist = []
    t_abs = 0.0
    for i in range(n_steps):
        if step_hook is not None:
            step_hook(mp, i)
        t_abs += dt
        mp.CloneTimeStep(t_abs)
        mp.ProcessInfo[KM.STEP] = i + 1
        strategy.Solve()
        if probe_fn is not None:
            v = float(probe_fn(mp))
            if not np.isfinite(v):
                raise RuntimeError(f"non-finite probe at step {i}")
            hist.append(v)
    return {"dt": float(dt), "n_steps": int(n_steps), "history": hist}


def paired_penalty_forces(gaps: np.ndarray, k_pen: float) -> np.ndarray:
    """Linear penalty force magnitudes for paired surface nodes.

    ``gaps > 0`` means separated (no force); penetration ``-gaps`` gives
    ``k_pen * penetration``. The caller applies ``-f`` on the lower node and
    ``+f`` on the upper node (internal pair: momentum-conserving by
    construction). Pure function: unit-testable without Kratos.
    """
    gaps = np.asarray(gaps, dtype=float)
    if float(k_pen) <= 0.0:
        raise ValueError(f"k_pen must be positive, got {k_pen}")
    return float(k_pen) * np.maximum(0.0, -gaps)


def coulomb_friction_update(g_t, k_t: float, mu: float, f_n: float,
                            ref=None):
    """Tangential stick predictor + radial return mapping (per contact pair).

    - ``g_t``: current tangential relative displacement (2-vector).
    - ``ref``: stick reference = accumulated slip (``None`` = virgin contact).
    - ``f_n``: current normal force magnitude (>= 0).
    Returns ``(f_t, ref_new, slipped)`` with ``|f_t| <= mu * f_n``.
    Pure function: unit-testable without Kratos. Sign-convention free
    (works with either normal-gap convention; only magnitudes enter).
    """
    g_t = np.asarray(g_t, dtype=float).reshape(2)
    ref = np.zeros(2, dtype=float) if ref is None else np.asarray(ref, dtype=float).reshape(2)
    trial = float(k_t) * (g_t - ref)
    tn = float(np.linalg.norm(trial))
    cap = float(mu) * float(f_n)
    if tn <= cap or tn == 0.0:
        return trial, ref.copy(), False
    direction = trial / tn
    f_t = cap * direction
    ref_new = g_t - f_t / float(k_t)
    return f_t, ref_new, True
