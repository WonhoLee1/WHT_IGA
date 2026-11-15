"""Coreform-conforming response spectrum and random response (plan add_1 G-3, V09).

Follows the Coreform IGA for Abaqus ``response-spectrum`` user-guide pages:
log-log spectra with damping interpolation (see :mod:`wht_iga.dynamics.case_file`),
the seven modal summation rules, the four directional rules, the Abaqus
frequency sweep and Segalman's RMS von Mises stress.

``modal_superposition.modal_response_spectrum`` / ``modal_random_psd`` keep their
original single-direction behaviour (linear interpolation, ABS/SRSS/CQC); this
module is the case-file driven path used by V09.

Documented assumptions (the user guide does not define them): the exact Abaqus
bias formula of the frequency sweep, and the textbook (NRC RG 1.92) forms of
``TENP``, ``GRP`` and ``DSC``.  Both are isolated in single functions so the
sensitivity of a result to them can be measured.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from wht_iga.dynamics.case_file import (
    CaseFileError,
    PSDType,
    RandomResponseCase,
    ResponseSpectrumCase,
)

#: Frequency ratio below which two modes are "closely spaced" (ten percent).
CLOSE_MODE_RATIO: float = 0.10

#: ASCE 4-98 secondary-direction fractions for the R40 / R30 rules.
_DIRECTIONAL_FRACTION = {"R40": 0.4, "R30": 0.3}

#: Von Mises quadratic form ``sigma_vm**2 = s^T A s`` for Voigt [xx,yy,zz,xy,yz,xz].
VON_MISES_A: np.ndarray = np.array(
    [
        [1.0, -0.5, -0.5, 0.0, 0.0, 0.0],
        [-0.5, 1.0, -0.5, 0.0, 0.0, 0.0],
        [-0.5, -0.5, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 3.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 3.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 3.0],
    ]
)

_QUANTITY_KIND = {"U": "D", "V": "V", "A": "A", "S": "D"}
_LAMBDA = {"ACCELERATION": 0, "VELOCITY": 1, "DISPLACEMENT": 2}


@dataclass(frozen=True)
class ResponseSpectrumResult:
    """Peak response of a case-file response-spectrum analysis.

    Attributes:
        peaks: Combined peak per output, shape ``(n_out,)`` for a 2D shape array
            or ``(n_out, n_comp)`` for a 3D one.  Always non-negative.
        per_direction: Modal-combined peak of every excitation before the
            directional rule, shape ``(n_exc, *peaks.shape)`` (one entry for
            ``ALGEBRAIC``, whose directions merge before the modal rule).
        frequencies_hz: Modal frequencies, shape ``(n_modes,)``.
        damping: Modal damping ratio used for the spectrum lookup and CQC/DSC.
        spectrum_values: ``scale * S(f_k, xi_k)`` of every excitation in the unit
            demanded by the output quantity, shape ``(n_exc, n_modes)``.
    """

    peaks: np.ndarray
    per_direction: np.ndarray
    frequencies_hz: np.ndarray
    damping: np.ndarray
    spectrum_values: np.ndarray


@dataclass(frozen=True)
class RandomResponseResult:
    """RMS results of a case-file random-response analysis.

    Attributes:
        frequencies_hz: Sweep frequencies used by the trapezoidal integration.
        ru: RMS displacement per DOF, shape ``(n_dofs,)``.
        rv: RMS velocity per DOF, shape ``(n_dofs,)``.
        ra: RMS acceleration per DOF, shape ``(n_dofs,)``.
        modal_variance: Real part of the generalized-coordinate variance matrix
            ``V_ab = int Re S^q_ab df``, shape ``(n_modes, n_modes)``.
        rmises: Segalman RMS von Mises stress per output point, or ``None``
            when no stress mode shapes were supplied.
    """

    frequencies_hz: np.ndarray
    ru: np.ndarray
    rv: np.ndarray
    ra: np.ndarray
    modal_variance: np.ndarray
    rmises: np.ndarray | None


def cqc_correlation(omega: np.ndarray, zeta: np.ndarray) -> np.ndarray:
    """Der Kiureghian (1981) CQC correlation matrix ``rho_ij``.

    ``rho_ij = 8 sqrt(xi_i xi_j) (xi_i + r xi_j) r**1.5 /
    [(1 - r**2)**2 + 4 xi_i xi_j r (1 + r**2) + 4 (xi_i**2 + xi_j**2) r**2]``
    with ``r = w_j / w_i``.  The diagonal is exactly one; an undamped pair with
    a vanishing denominator is treated as uncorrelated.

    Args:
        omega: Modal circular frequencies, shape ``(n,)``.
        zeta: Modal damping ratios, shape ``(n,)``.

    Returns:
        Symmetric matrix of shape ``(n, n)``.
    """
    w = np.asarray(omega, dtype=float).reshape(-1)
    z = np.asarray(zeta, dtype=float).reshape(-1)
    r = w[None, :] / w[:, None]
    zi, zj = z[:, None], z[None, :]
    num = 8.0 * np.sqrt(zi * zj) * (zi + r * zj) * r**1.5
    den = (1.0 - r * r) ** 2 + 4.0 * zi * zj * r * (1.0 + r * r) + 4.0 * (zi * zi + zj * zj) * r * r
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(den > 0.0, num / den, 0.0)
    np.fill_diagonal(rho, 1.0)
    return rho


def dsc_correlation(omega: np.ndarray, zeta: np.ndarray, duration: float) -> np.ndarray:
    """Double-sum (Rosenblueth-Elorduy / RG 1.92) coefficients ``eps_ij``.

    ``eps_ij = 1 / (1 + ((w'_i - w'_j) / (xi'_i w_i + xi'_j w_j))**2)`` with the
    damped frequency ``w' = w sqrt(1 - xi**2)`` and the duration-corrected damping
    ``xi'_i = xi_i + 2 / (w_i t_d)``.

    Args:
        omega: Modal circular frequencies, shape ``(n,)``.
        zeta: Modal damping ratios, shape ``(n,)``.
        duration: Strong-motion duration ``t_d`` (positive).

    Returns:
        Symmetric matrix of shape ``(n, n)`` with unit diagonal.
    """
    w = np.asarray(omega, dtype=float).reshape(-1)
    z = np.asarray(zeta, dtype=float).reshape(-1)
    wd = w * np.sqrt(1.0 - z * z)
    zp = z + 2.0 / (w * float(duration))
    eps = 1.0 / (1.0 + ((wd[:, None] - wd[None, :]) / (zp[:, None] * w[:, None] + zp[None, :] * w[None, :])) ** 2)
    np.fill_diagonal(eps, 1.0)
    return eps


def combine_modal_peaks(
    r: np.ndarray,
    omega: np.ndarray,
    zeta: np.ndarray,
    rule: str,
    duration: float | None = None,
) -> np.ndarray:
    """Combine signed modal peak amplitudes with one of the seven Abaqus rules.

    Rules (``R_k`` = modal peak of one output component, modes sorted by frequency):

    * ``ABS``: ``sum |R_k|``
    * ``SRSS``: ``sqrt(sum R_k**2)``
    * ``NRL``: ``max|R_k| + sqrt(sum of the squares of the others)``
    * ``TENP``: ``sqrt(sum R_k**2 + 2 sum_{i<j, (w_j-w_i)/w_i <= 0.1} |R_i R_j|)``
    * ``CQC``: ``sqrt(R^T rho R)`` with signed ``R`` and :func:`cqc_correlation`
    * ``GRP``: groups start at the lowest unassigned mode and take every mode with
      ``w <= 1.1 w_start``; ``sqrt(sum_groups (sum_{k in group} |R_k|)**2)``
    * ``DSC``: ``sqrt(sum_ij |R_i R_j| eps_ij)`` with :func:`dsc_correlation`

    Args:
        r: Signed modal amplitudes, shape ``(n_out, n_modes)``.
        omega: Modal circular frequencies, shape ``(n_modes,)``.
        zeta: Modal damping ratios, shape ``(n_modes,)``.
        rule: One of the rule names above (case-insensitive).
        duration: Strong-motion duration, required for ``DSC``.

    Returns:
        Combined peak per output, shape ``(n_out,)``.
    """
    rule = str(rule).upper()
    r = np.atleast_2d(np.asarray(r, dtype=float))
    w = np.asarray(omega, dtype=float).reshape(-1)
    z = np.asarray(zeta, dtype=float).reshape(-1)
    if r.shape[1] != w.size or z.size != w.size:
        raise ValueError(f"r has {r.shape[1]} modes, omega {w.size}, zeta {z.size}")
    if np.any(w <= 0.0):
        raise ValueError("response spectrum combination needs positive modal frequencies")
    order = np.argsort(w, kind="stable")
    r, w, z = r[:, order], w[order], z[order]
    a = np.abs(r)
    if rule == "ABS":
        return a.sum(axis=1)
    if rule == "SRSS":
        return np.sqrt((r * r).sum(axis=1))
    if rule == "NRL":
        sq = (r * r).sum(axis=1)
        top = a.max(axis=1)
        return top + np.sqrt(np.maximum(sq - top * top, 0.0))
    if rule == "TENP":
        close = (w[None, :] - w[:, None]) / w[:, None] <= CLOSE_MODE_RATIO
        upper = np.triu(close, k=1).astype(float)
        return np.sqrt((r * r).sum(axis=1) + 2.0 * np.einsum("ni,ij,nj->n", a, upper, a))
    if rule == "CQC":
        return np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", r, cqc_correlation(w, z), r), 0.0))
    if rule == "GRP":
        total = np.zeros(r.shape[0])
        start = 0
        while start < w.size:
            end = start + 1
            while end < w.size and w[end] <= (1.0 + CLOSE_MODE_RATIO) * w[start]:
                end += 1
            total += a[:, start:end].sum(axis=1) ** 2
            start = end
        return np.sqrt(total)
    if rule == "DSC":
        if duration is None or duration <= 0.0:
            raise ValueError("DSC needs a positive strong-motion duration")
        return np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", a, dsc_correlation(w, z, duration), a), 0.0))
    raise ValueError(f"unknown modal summation rule '{rule}'")


def combine_directions(per_direction: np.ndarray, rule: str) -> np.ndarray:
    """Combine per-direction peaks ``R_d >= 0`` (ALGEBRAIC is applied before the modal rule).

    ``SRSS``: ``sqrt(sum R_d**2)``.  ``R40`` / ``R30`` (ASCE 4-98): the largest of
    ``R_a + f (R_b + R_c)`` over the choice of the dominant direction ``a``,
    ``f = 0.4`` or ``0.3``, which equals ``f sum R + (1 - f) max R``.

    Args:
        per_direction: Shape ``(n_dir, ...)``.
        rule: ``SRSS``, ``R40`` or ``R30``.

    Returns:
        Combined peak, shape ``per_direction.shape[1:]``.
    """
    rule = str(rule).upper()
    p = np.abs(np.asarray(per_direction, dtype=float))
    if rule == "SRSS":
        return np.sqrt((p * p).sum(axis=0))
    if rule in _DIRECTIONAL_FRACTION:
        f = _DIRECTIONAL_FRACTION[rule]
        return f * p.sum(axis=0) + (1.0 - f) * p.max(axis=0)
    raise ValueError(f"unknown directional summation rule '{rule}'")


def response_spectrum_from_case(
    omega_modes: np.ndarray,
    shapes: np.ndarray,
    participation_factors: np.ndarray,
    case: ResponseSpectrumCase,
    quantity: str = "U",
) -> ResponseSpectrumResult:
    """Run a Coreform-style response-spectrum case on a mass-normalised modal basis.

    For excitation ``e`` (unit direction ``d_e``, scale ``s_e``) the signed modal
    amplitude of output component ``N`` is
    ``R_ek,N = Gamma_k(d_e) * s_e * S_q(f_k, xi_k) * shape_kN`` with
    ``Gamma_k(d) = sum_i d_i Gamma_k,i`` and ``S_q`` the spectrum in the unit of
    the requested quantity (displacement, velocity, acceleration; stress mode
    shapes use the displacement spectrum).  ``ALGEBRAIC`` sums the excitations
    first and applies the modal rule once; the other directional rules apply
    the modal rule per excitation and combine afterwards.

    Args:
        omega_modes: Modal circular frequencies, shape ``(n_modes,)``.
        shapes: Mode shapes of the output quantity, ``(n_out, n_modes)`` or
            ``(n_out, n_modes, n_comp)`` (e.g. ``Phi`` or stress shapes).
        participation_factors: Participation factors ``(n_modes, 3)``.
        case: A :class:`~wht_iga.dynamics.case_file.ResponseSpectrumCase`.
        quantity: ``"U"``, ``"V"``, ``"A"`` or ``"S"``.

    Returns:
        :class:`ResponseSpectrumResult`.
    """
    q = quantity.upper()
    if q not in _QUANTITY_KIND:
        raise ValueError(f"quantity must be U, V, A or S, got {quantity!r}")
    kind = _QUANTITY_KIND[q]
    w = np.asarray(omega_modes, dtype=float).reshape(-1)
    n_modes = w.size
    gam = np.asarray(participation_factors, dtype=float)
    if gam.shape != (n_modes, 3):
        raise ValueError(f"participation_factors must have shape ({n_modes}, 3), got {gam.shape}")
    sh = np.asarray(shapes, dtype=float)
    if sh.ndim not in (2, 3) or sh.shape[1] != n_modes:
        raise ValueError(f"shapes must be (n_out, {n_modes}) or (n_out, {n_modes}, n_comp), got {sh.shape}")
    out_shape = (sh.shape[0],) if sh.ndim == 2 else (sh.shape[0], sh.shape[2])
    flat = (sh[:, :, None] if sh.ndim == 2 else sh).transpose(0, 2, 1).reshape(-1, n_modes)

    f_hz = w / (2.0 * np.pi)
    xi = case.damping.ratios(f_hz)
    spec_vals, amps = [], []
    for e in case.excitations:
        s = e.scale * case.spectra[e.spectrum].at(kind, f_hz, xi)
        spec_vals.append(s)
        amps.append((gam @ np.asarray(e.direction)) * s * flat)

    rule = case.modal_summation.value
    if case.directional_summation.value == "ALGEBRAIC":
        if rule == "DSC" and len({e.duration for e in case.excitations}) > 1:
            raise ValueError("ALGEBRAIC + DSC needs one common duration on all excitations")
        merged = combine_modal_peaks(np.sum(amps, axis=0), w, xi, rule, case.excitations[0].duration)
        per_dir = merged.reshape(1, *out_shape)
        peaks = per_dir[0]
    else:
        per_dir = np.asarray(
            [
                combine_modal_peaks(a, w, xi, rule, e.duration).reshape(out_shape)
                for a, e in zip(amps, case.excitations, strict=True)
            ]
        )
        peaks = combine_directions(per_dir, case.directional_summation.value)
    return ResponseSpectrumResult(
        peaks=peaks, per_direction=per_dir, frequencies_hz=f_hz, damping=xi, spectrum_values=np.asarray(spec_vals)
    )


def _bias_map(xi: np.ndarray, bias: float, model: str) -> np.ndarray:
    """Map uniform ``xi in [-1, 1]`` to biased ``t in [-1, 1]`` clustering toward the ends.

    ``abaqus``: ``t = sign(xi) |xi|**(1/bias)`` -- the Abaqus bias formula
    ``f_k = (f1 + f2)/2 + (f2 - f1)/2 |y|^(1/p) sign(y)``, ``y = -1 + 2 (k-1)/(n-1)``
    (Abaqus Analysis User's Guide, mode-based steady-state dynamics / random response;
    ``f`` is the frequency or its logarithm depending on the scale).
    ``sinh_ratio``: ``t = asinh(xi sinh c) / c`` with ``cosh(c) = bias`` (mid/end spacing
    ratio = ``bias``); kept only to measure sensitivity to the sweep.
    ``bias == 1`` is uniform in both.
    """
    if bias == 1.0:
        return xi
    if model in ("abaqus", "power"):
        return np.sign(xi) * np.abs(xi) ** (1.0 / bias)
    if model == "sinh_ratio":
        c = np.arccosh(bias)
        return np.arcsinh(xi * np.sinh(c)) / c
    raise ValueError(f"unknown bias model '{model}'")


def abaqus_frequency_sweep(
    f_modes: np.ndarray,
    f_lower: float,
    f_upper: float,
    n_points: int = 20,
    bias: float = 3.0,
    scale: str = "LOGARITHMIC",
    model: str = "abaqus",
) -> np.ndarray:
    """Frequency points of a random-response sweep (Abaqus convention by default).

    The range ``[f_lower, f_upper]`` is divided at every eigenfrequency inside it.
    Each interval gets ``n_points`` points, both bounding frequencies included (the
    Abaqus manual: "the user-defined number of points (which includes the bounding
    frequencies for the interval)"), clustered toward the interval ends by the Abaqus
    bias formula (see :func:`_bias_map`).  The formula was read from the Abaqus
    documentation, not from the Coreform guide (which only names it); it is the
    convention the Coreform ``respspec`` sweep claims to match point for point.

    Args:
        f_modes: Eigenfrequencies in Hz.
        f_lower: Lower sweep bound in Hz (> 0).
        f_upper: Upper sweep bound in Hz.
        n_points: Points per interval (>= 2).
        bias: Clustering parameter (>= 1; 1 = uniform).
        scale: ``"LINEAR"`` or ``"LOGARITHMIC"`` spacing inside an interval.
        model: Bias mapping, ``"abaqus"`` (default, alias ``"power"``) or
            ``"sinh_ratio"`` (sensitivity only).

    Returns:
        Sorted, unique frequencies in Hz.
    """
    if not 0.0 < f_lower < f_upper:
        raise ValueError("need 0 < f_lower < f_upper")
    if n_points < 2 or bias < 1.0:
        raise ValueError("need n_points >= 2 and bias >= 1")
    modes = np.asarray(f_modes, dtype=float).reshape(-1)
    inner = np.sort(modes[(modes > f_lower) & (modes < f_upper)])
    edges = np.concatenate([[f_lower], inner, [f_upper]])
    edges = edges[np.concatenate([[True], np.diff(edges) > 1e-9 * edges[1:]])]
    t = 0.5 * (1.0 + _bias_map(np.linspace(-1.0, 1.0, int(n_points)), float(bias), model))
    pts = []
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        seg = a * (b / a) ** t if scale.upper() == "LOGARITHMIC" else a + (b - a) * t
        seg[0], seg[-1] = a, b  # shared ends must be bit-identical so np.unique merges them
        pts.append(seg)
    return np.unique(np.concatenate(pts))


def trapezoid_weights(x: np.ndarray) -> np.ndarray:
    """Weights ``w_i`` such that ``sum_i w_i g_i`` is the trapezoidal rule on ``x``."""
    d = np.diff(x)
    w = np.zeros_like(x)
    w[:-1] += 0.5 * d
    w[1:] += 0.5 * d
    return w


def random_response_from_case(
    omega_modes: np.ndarray,
    modes: np.ndarray,
    participation_factors: np.ndarray,
    case: RandomResponseCase,
    stress_modes: np.ndarray | None = None,
    frequencies_hz: np.ndarray | None = None,
    bias_model: str = "abaqus",
) -> RandomResponseResult:
    """RMS displacement/velocity/acceleration and Segalman RMISES for a random case.

    Implements the Coreform guide equations: the base-motion cross-spectral
    density ``S^f_ab = sum_J g^2 P^J(f) sum_I [sum_ij Gamma^i_a Gamma^j_b Psi_ij] (2 pi f)^(2 lambda_I)``
    (``CORRELATED``: participation factors of one load case add before the outer
    product; ``UNCORRELATED``: diagonal only), the modal transfer function
    ``H_a = 1 / (w_a**2 - w**2 + 2 i xi_a w_a w)``, ``S^q_ab = H_a S^f_ab H_b^*``
    and trapezoidal integration over the sweep.  Velocity and acceleration carry
    ``(2 pi f)**2`` and ``(2 pi f)**4``.  RMISES is ``sqrt(sum_ab V_ab T_ab)`` with
    ``T_ab = psi_a^T A psi_b`` (Segalman et al.).

    Args:
        omega_modes: Modal circular frequencies, shape ``(n_modes,)``.
        modes: Mass-orthonormal mode shapes, shape ``(n_dofs, n_modes)``.
        participation_factors: Shape ``(n_modes, 3)``.
        case: A :class:`~wht_iga.dynamics.case_file.RandomResponseCase`.
        stress_modes: Optional modal stress vectors ``(n_pts, n_modes, 6)`` in
            Voigt order ``[xx, yy, zz, xy, yz, xz]`` for RMISES.
        frequencies_hz: Override of the sweep (default: :func:`abaqus_frequency_sweep`).
        bias_model: Bias mapping forwarded to :func:`abaqus_frequency_sweep`.

    Returns:
        :class:`RandomResponseResult`.
    """
    w = np.asarray(omega_modes, dtype=float).reshape(-1)
    phi = np.asarray(modes, dtype=float)
    gam = np.asarray(participation_factors, dtype=float)
    if phi.shape[1] != w.size or gam.shape != (w.size, 3):
        raise ValueError("modes / participation_factors do not match omega_modes")
    if np.any(w <= 0.0):
        raise ValueError("random response needs positive modal frequencies")
    f_modes = w / (2.0 * np.pi)
    xi = case.damping.ratios(f_modes)
    freqs = (
        np.asarray(frequencies_hz, dtype=float)
        if frequencies_hz is not None
        else abaqus_frequency_sweep(
            f_modes, case.lower_frequency_hz, case.upper_frequency_hz, case.n_points, case.bias,
            case.frequency_scale, bias_model,
        )
    )
    wf = 2.0 * np.pi * freqs

    sources: list[tuple[np.ndarray, np.ndarray]] = []  # (gamma (n_modes,), power(f) (n_f,))
    for corr in case.correlations:
        psd = case.psds[corr.psd]
        if psd.type is PSDType.FORCE:
            raise CaseFileError("FORCE PSDs need applied-force loads; only base-motion excitation is implemented")
        pf = psd.value(freqs)
        if np.any(pf.imag != 0.0):
            raise CaseFileError(f"psd '{psd.name}': complex (cross-spectral) ordinates are not supported")
        pf = pf.real * (psd.gravity**2 if psd.type is PSDType.BASE else 1.0)
        for lc, (wre, wim) in corr.weights.items():
            moves = [m for m in case.base_motions if m.load_case == lc]
            power = pf * (wre * wre + wim * wim) * wf ** (2 * _LAMBDA[moves[0].kind])
            if corr.type == "CORRELATED":
                sources.append((sum(gam[:, m.dof - 1] for m in moves), power))
            else:
                sources.extend((gam[:, m.dof - 1], power) for m in moves)

    h = 1.0 / (w[None, :] ** 2 - wf[:, None] ** 2 + 2j * xi[None, :] * w[None, :] * wf[:, None])
    tw = trapezoid_weights(freqs)
    var = []
    for freq_factor in (np.ones_like(wf), wf**2, wf**4):
        v = np.zeros((w.size, w.size))
        for gvec, power in sources:
            z = h * gvec[None, :]
            v += np.real(np.einsum("f,fa,fb->ab", tw * freq_factor * power, z, np.conj(z)))
        var.append(v)
    ru, rv, ra = (np.sqrt(np.maximum(np.einsum("na,ab,nb->n", phi, v, phi), 0.0)) for v in var)

    rmises = None
    if stress_modes is not None:
        psi = np.asarray(stress_modes, dtype=float)
        if psi.ndim != 3 or psi.shape[1] != w.size or psi.shape[2] != 6:
            raise ValueError(f"stress_modes must be (n_pts, {w.size}, 6), got {psi.shape}")
        t_ab = np.einsum("pak,kl,pbl->pab", psi, VON_MISES_A, psi)
        rmises = np.sqrt(np.maximum(np.einsum("pab,ab->p", t_ab, var[0]), 0.0))
    return RandomResponseResult(frequencies_hz=freqs, ru=ru, rv=rv, ra=ra, modal_variance=var[0], rmises=rmises)
