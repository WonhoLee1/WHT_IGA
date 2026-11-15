"""Pure NumPy/SciPy Quantity-of-Interest (QoI) extractors and calculators.

Implements Kratos-independent analytical and post-processing metrics for
WHT_IGA verification problems V01–V12 and example problems E01–E05.
"""

from __future__ import annotations

from typing import NamedTuple, Sequence

import numpy as np


class LameSphereExact(NamedTuple):
    """Exact Lamé solution for a thick-walled pressurized hollow sphere.

    Attributes
    ----------
    u_vec:
        Cartesian displacement vectors, shape ``(N, 3)``.
    u_r:
        Radial displacement, shape ``(N,)``.
    sigma_rr:
        Radial stress, shape ``(N,)``.
    sigma_tt:
        Tangential (hoop) stress, shape ``(N,)``.
    stress_voigt:
        Cartesian stress in 6-component Voigt order
        ``[sxx, syy, szz, sxy, syz, sxz]``, shape ``(N, 6)``.
    """

    u_vec: np.ndarray
    u_r: np.ndarray
    sigma_rr: np.ndarray
    sigma_tt: np.ndarray
    stress_voigt: np.ndarray


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Integrate ``y`` along ``x`` via trapezoidal rule (NumPy 1.x/2.x safe)."""
    trap_fn = getattr(np, "trapezoid", None) or np.trapz
    return float(trap_fn(y, x))


def _to_stress_matrices(stress_voigt_or_mat: np.ndarray) -> tuple[np.ndarray, bool]:
    """Normalize ``(6,)``, ``(N, 6)``, ``(3, 3)``, or ``(N, 3, 3)`` to ``(N, 3, 3)``."""
    arr = np.asarray(stress_voigt_or_mat, dtype=np.float64)
    if arr.shape == (3, 3):
        return arr[np.newaxis, :, :], True
    if arr.ndim == 3 and arr.shape[1:] == (3, 3):
        return arr, False
    if arr.shape == (6,):
        arr = arr[np.newaxis, :]
        was_single = True
    elif arr.ndim == 2 and arr.shape[1] == 6:
        was_single = False
    else:
        raise ValueError(
            "Expected stress shape (6,), (N, 6), (3, 3), or (N, 3, 3), "
            f"got {arr.shape}"
        )

    n = arr.shape[0]
    mats = np.empty((n, 3, 3), dtype=np.float64)
    mats[:, 0, 0] = arr[:, 0]
    mats[:, 1, 1] = arr[:, 1]
    mats[:, 2, 2] = arr[:, 2]
    mats[:, 0, 1] = mats[:, 1, 0] = arr[:, 3]
    mats[:, 1, 2] = mats[:, 2, 1] = arr[:, 4]
    mats[:, 0, 2] = mats[:, 2, 0] = arr[:, 5]
    return mats, was_single


def principal_stresses(
    stress_voigt_or_mat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute sorted principal stresses and positive maximum principal stress.

    Parameters
    ----------
    stress_voigt_or_mat:
        Symmetric stress tensor(s) as ``(6,)`` or ``(N, 6)`` Voigt vectors
        ``[sxx, syy, szz, sxy, syz, sxz]``, or ``(3, 3)`` / ``(N, 3, 3)``
        matrices.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        ``(principals, sigma1_plus)`` where ``principals`` has descending
        eigenvalues ``(sigma_1 >= sigma_2 >= sigma_3)`` of shape ``(3,)`` or
        ``(N, 3)``, and ``sigma1_plus = max(sigma_1, 0)`` has shape ``()`` or
        ``(N,)``.
    """
    mats, was_single = _to_stress_matrices(stress_voigt_or_mat)
    # Enforce symmetry in case of minor floating-point asymmetry in 3x3 input
    sym = 0.5 * (mats + np.swapaxes(mats, -1, -2))
    eigvals = np.linalg.eigvalsh(sym)
    principals = eigvals[:, ::-1].copy()
    sigma1_plus = np.maximum(principals[:, 0], 0.0)
    if was_single:
        return principals[0], sigma1_plus[0]
    return principals, sigma1_plus


def line_profile_rel_l2(
    s_coords: np.ndarray | Sequence[float],
    values_num: np.ndarray | Sequence[float],
    values_ref: np.ndarray | Sequence[float],
) -> float:
    """Compute relative L2 error along a 1D arc-length profile using trapezoids.

    Parameters
    ----------
    s_coords:
        1D arc-length coordinates, shape ``(N,)`` with ``N >= 2``.
    values_num:
        Numerical values along the line, shape ``(N,)`` or ``(N, D)``.
    values_ref:
        Reference values along the line, matching ``values_num`` shape.

    Returns
    -------
    float
        Relative L2 error ``sqrt(int ||v_num - v_ref||^2 ds / int ||v_ref||^2 ds)``.
    """
    s = np.asarray(s_coords, dtype=np.float64).ravel()
    if s.size < 2:
        raise ValueError(f"Expected at least 2 points along line, got {s.size}")

    v_num = np.asarray(values_num, dtype=np.float64)
    v_ref = np.asarray(values_ref, dtype=np.float64)
    if v_num.shape[0] != s.size or v_ref.shape != v_num.shape:
        raise ValueError(
            f"Shape mismatch: s={s.shape}, num={v_num.shape}, ref={v_ref.shape}"
        )

    order = np.argsort(s)
    s_sorted = s[order]
    diff = v_num[order] - v_ref[order]
    ref_sorted = v_ref[order]

    if diff.ndim == 1:
        diff_sq = diff * diff
        ref_sq = ref_sorted * ref_sorted
    else:
        axes = tuple(range(1, diff.ndim))
        diff_sq = np.sum(diff * diff, axis=axes)
        ref_sq = np.sum(ref_sorted * ref_sorted, axis=axes)

    err_int = max(_trapezoid(diff_sq, s_sorted), 0.0)
    ref_int = max(_trapezoid(ref_sq, s_sorted), 0.0)
    if ref_int == 0.0:
        return 0.0 if err_int == 0.0 else float("inf")
    return float(np.sqrt(err_int / ref_int))


def spherical_projections(
    points: np.ndarray,
    u: np.ndarray,
    stress_voigt: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project Cartesian displacement and stress onto spherical coordinates.

    Parameters
    ----------
    points:
        Evaluation points ``(x, y, z)``, shape ``(N, 3)``.
    u:
        Cartesian displacement vectors ``(ux, uy, uz)``, shape ``(N, 3)``.
    stress_voigt:
        Cartesian stress in Voigt ``(N, 6)`` or tensor ``(N, 3, 3)`` format.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        ``(u_r, sigma_rr, sigma_tt)``, each of shape ``(N,)``.
    """
    pts = np.asarray(points, dtype=np.float64)
    disp = np.asarray(u, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"Expected points shape (N, 3), got {pts.shape}")
    if disp.shape != pts.shape:
        raise ValueError(f"Expected u shape {pts.shape}, got {disp.shape}")

    radii = np.linalg.norm(pts, axis=1)
    if np.any(radii == 0.0):
        raise ValueError("Spherical projection is undefined at the origin r = 0.")
    e_r = pts / radii[:, np.newaxis]

    u_r = np.sum(disp * e_r, axis=1)

    mats, was_single = _to_stress_matrices(stress_voigt)
    if was_single and pts.shape[0] != 1:
        raise ValueError("Single stress tensor provided for multiple points.")
    if mats.shape[0] != pts.shape[0]:
        raise ValueError(
            f"Stress count {mats.shape[0]} does not match points {pts.shape[0]}"
        )

    sigma_rr = np.einsum("ni,nij,nj->n", e_r, mats, e_r)
    trace = np.trace(mats, axis1=1, axis2=2)
    sigma_tt = 0.5 * (trace - sigma_rr)
    return u_r, sigma_rr, sigma_tt


def lame_sphere_exact(
    points: np.ndarray,
    a: float,
    b: float,
    p_i: float,
    p_o: float,
    E: float,
    nu: float,
) -> LameSphereExact:
    """Evaluate exact Lamé solution for a pressurized thick-walled sphere.

    Parameters
    ----------
    points:
        Cartesian coordinates ``(x, y, z)``, shape ``(N, 3)``.
    a:
        Inner radius (``0 < a < b``).
    b:
        Outer radius.
    p_i:
        Internal pressure.
    p_o:
        External pressure.
    E:
        Young's modulus.
    nu:
        Poisson's ratio.

    Returns
    -------
    LameSphereExact
        NamedTuple containing ``u_vec``, ``u_r``, ``sigma_rr``, ``sigma_tt``,
        and ``stress_voigt``.
    """
    if not (0.0 < a < b):
        raise ValueError(f"Require 0 < a < b, got a={a}, b={b}")
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"Expected points shape (N, 3), got {pts.shape}")

    radii = np.linalg.norm(pts, axis=1)
    if np.any(radii == 0.0):
        raise ValueError("Lamé solution is singular at r = 0.")
    e_r = pts / radii[:, np.newaxis]

    denom = b**3 - a**3
    c = (p_i * a**3 - p_o * b**3) / denom
    kappa = (p_i - p_o) * (a**3) * (b**3) / denom

    r2 = radii * radii
    r3 = r2 * radii
    sigma_rr = c - kappa / r3
    sigma_tt = c + kappa / (2.0 * r3)
    u_r = (1.0 / E) * ((1.0 - 2.0 * nu) * c * radii + (1.0 + nu) * kappa / (2.0 * r2))
    u_vec = u_r[:, np.newaxis] * e_r

    # Cartesian stress: sigma = sigma_rr * (e_r x e_r) + sigma_tt * (I - e_r x e_r)
    diff_s = sigma_rr - sigma_tt
    sxx = sigma_tt + diff_s * e_r[:, 0] * e_r[:, 0]
    syy = sigma_tt + diff_s * e_r[:, 1] * e_r[:, 1]
    szz = sigma_tt + diff_s * e_r[:, 2] * e_r[:, 2]
    sxy = diff_s * e_r[:, 0] * e_r[:, 1]
    syz = diff_s * e_r[:, 1] * e_r[:, 2]
    sxz = diff_s * e_r[:, 0] * e_r[:, 2]
    stress_voigt = np.column_stack([sxx, syy, szz, sxy, syz, sxz])

    return LameSphereExact(
        u_vec=u_vec,
        u_r=u_r,
        sigma_rr=sigma_rr,
        sigma_tt=sigma_tt,
        stress_voigt=stress_voigt,
    )


def log_log_convergence_rate(
    h_values: np.ndarray | Sequence[float],
    errors: np.ndarray | Sequence[float],
) -> float:
    """Compute least-squares convergence rate ``alpha`` in ``log(e) = alpha * log(h) + c``.

    Parameters
    ----------
    h_values:
        Positive mesh sizes, length ``>= 2``.
    errors:
        Positive error norms corresponding to ``h_values``.

    Returns
    -------
    float
        Estimated algebraic convergence rate ``alpha``.
    """
    h = np.asarray(h_values, dtype=np.float64).ravel()
    err = np.asarray(errors, dtype=np.float64).ravel()
    if h.size < 2 or h.size != err.size:
        raise ValueError(
            f"Require equal-length arrays with >= 2 points, got {h.size} and {err.size}"
        )
    if np.any(h <= 0.0) or np.any(err <= 0.0):
        raise ValueError("h_values and errors must all be strictly positive.")

    x = np.log(h)
    y = np.log(err)
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    dx = x - x_mean
    denom = float(np.dot(dx, dx))
    if denom == 0.0:
        raise ValueError("h_values must contain at least two distinct values.")
    return float(np.dot(dx, y - y_mean) / denom)


def time_window_rms_error(
    t: np.ndarray | Sequence[float],
    q_num: np.ndarray | Sequence[float],
    q_ref: np.ndarray | Sequence[float],
    t_min: float,
    t_max: float,
) -> float:
    """Compute relative RMS error within a closed time window ``[t_min, t_max]``.

    Parameters
    ----------
    t:
        Time stamps, shape ``(N,)``.
    q_num:
        Numerical time-history values, shape ``(N,)`` or ``(N, D)``.
    q_ref:
        Reference time-history values, matching ``q_num`` shape.
    t_min:
        Start of comparison window (inclusive).
    t_max:
        End of comparison window (inclusive).

    Returns
    -------
    float
        Relative RMS error ``sqrt(sum(err^2)) / sqrt(sum(ref^2))`` in window.
    """
    t_arr = np.asarray(t, dtype=np.float64).ravel()
    num_arr = np.asarray(q_num, dtype=np.float64)
    ref_arr = np.asarray(q_ref, dtype=np.float64)
    if num_arr.shape[0] != t_arr.size or ref_arr.shape != num_arr.shape:
        raise ValueError(
            f"Shape mismatch: t={t_arr.shape}, num={num_arr.shape}, ref={ref_arr.shape}"
        )

    mask = (t_arr >= t_min) & (t_arr <= t_max)
    if not np.any(mask):
        raise ValueError(f"No time samples found in window [{t_min}, {t_max}].")

    diff_sq = float(np.sum((num_arr[mask] - ref_arr[mask]) ** 2))
    ref_sq = float(np.sum(ref_arr[mask] ** 2))
    if ref_sq == 0.0:
        return 0.0 if diff_sq == 0.0 else float("inf")
    return float(np.sqrt(diff_sq / ref_sq))


def filter_participating_modes(
    freqs: np.ndarray | Sequence[float],
    participation_factors: np.ndarray | Sequence[float],
    threshold: float = 0.01,
) -> np.ndarray:
    """Return indices of modes whose participation exceeds ``threshold * max(|Gamma|)``.

    Parameters
    ----------
    freqs:
        Modal frequencies, shape ``(M,)``.
    participation_factors:
        Modal participation factors ``Gamma_i``, shape ``(M,)``.
    threshold:
        Relative cutoff fraction (default ``0.01`` per V09).

    Returns
    -------
    np.ndarray
        1D integer array of participating mode indices.
    """
    f_arr = np.asarray(freqs, dtype=np.float64).ravel()
    gamma = np.asarray(participation_factors, dtype=np.float64).ravel()
    if f_arr.size != gamma.size:
        raise ValueError(
            f"Length mismatch: freqs={f_arr.size}, participation={gamma.size}"
        )
    if gamma.size == 0:
        return np.array([], dtype=np.intp)

    abs_gamma = np.abs(gamma)
    max_gamma = float(np.max(abs_gamma))
    if max_gamma == 0.0:
        return np.array([], dtype=np.intp)
    return np.flatnonzero(abs_gamma >= threshold * max_gamma)


def analyze_dynamic_history(
    t: np.ndarray | Sequence[float],
    u: np.ndarray | Sequence[float],
    t_free_start: float,
    ref_freq: float | None = None,
) -> tuple[float, float, float, float]:
    """Analyze free-vibration dynamic response history (V10).

    Parameters
    ----------
    t:
        Time stamps, shape ``(N,)``.
    u:
        Displacement history, shape ``(N,)``.
    t_free_start:
        Start time of the free-vibration regime.
    ref_freq:
        Optional reference frequency in Hz used to compute ``samples_per_period``
        and normalized reference-frequency residual; if ``None``, uses the
        extracted fundamental frequency ``f1_hz``.

    Returns
    -------
    tuple[float, float, float, float]
        ``(f1_hz, norm_rms, amp_decay_ratio, samples_per_period)`` where:
        - ``f1_hz``: fundamental frequency in Hz via zero-crossing / FFT.
        - ``norm_rms``: RMS normalized by peak free-vibration amplitude.
        - ``amp_decay_ratio``: ratio of last envelope peak to first envelope peak.
        - ``samples_per_period``: average number of samples per fundamental period.
    """
    t_arr = np.asarray(t, dtype=np.float64).ravel()
    u_arr = np.asarray(u, dtype=np.float64).ravel()
    if t_arr.size != u_arr.size or t_arr.size < 4:
        raise ValueError("Require matching time/displacement arrays with >= 4 samples.")

    mask = t_arr >= t_free_start
    tf = t_arr[mask]
    uf = u_arr[mask]
    if tf.size < 4:
        raise ValueError(f"Insufficient samples for t >= {t_free_start}.")

    uf_centered = uf - float(np.mean(uf))

    # Upward zero-crossings with linear interpolation
    cross_times: list[float] = []
    for i in range(uf_centered.size - 1):
        y0, y1 = uf_centered[i], uf_centered[i + 1]
        if y0 < 0.0 <= y1 or (y0 == 0.0 and y1 > 0.0):
            denom = y1 - y0
            frac = 0.0 if denom == 0.0 else -y0 / denom
            cross_times.append(float(tf[i] + frac * (tf[i + 1] - tf[i])))

    if len(cross_times) >= 2:
        periods = np.diff(cross_times)
        f1_hz = float(1.0 / np.mean(periods))
    else:
        # FFT fallback on uniform grid
        n_fft = max(tf.size, 256)
        t_uni = np.linspace(tf[0], tf[-1], n_fft)
        u_uni = np.interp(t_uni, tf, uf_centered)
        dt_uni = float(t_uni[1] - t_uni[0])
        spec = np.abs(np.fft.rfft(u_uni))
        freqs = np.fft.rfftfreq(n_fft, d=dt_uni)
        if spec.size <= 1:
            raise ValueError("Unable to determine fundamental frequency from signal.")
        peak_idx = int(np.argmax(spec[1:]) + 1)
        f1_hz = float(freqs[peak_idx])

    peak_amp = float(np.max(np.abs(uf_centered)))
    if peak_amp == 0.0:
        norm_rms = 0.0
    else:
        norm_rms = float(np.sqrt(np.mean(uf_centered**2)) / peak_amp)

    # Local positive peaks of |uf_centered| for amplitude decay ratio
    abs_u = np.abs(uf_centered)
    peak_indices = [
        i
        for i in range(1, abs_u.size - 1)
        if abs_u[i] >= abs_u[i - 1] and abs_u[i] > abs_u[i + 1]
    ]
    if len(peak_indices) >= 2:
        amp_decay_ratio = float(abs_u[peak_indices[-1]] / abs_u[peak_indices[0]])
    else:
        amp_decay_ratio = 1.0

    dt_mean = float(np.mean(np.diff(tf)))
    target_freq = ref_freq if ref_freq is not None else f1_hz
    samples_per_period = float((1.0 / target_freq) / dt_mean)

    return f1_hz, norm_rms, amp_decay_ratio, samples_per_period
