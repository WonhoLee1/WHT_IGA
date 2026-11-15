"""General-purpose PSD random vibration analysis engine (Abaqus / ANSYS / Nastran parity).

Provides a pure-Python / NumPy / SciPy object-oriented framework for stationary
random vibration analysis without requiring JSON case files or Kratos runtime
dependencies, while remaining directly compatible with :class:`EigenmodeResult`,
:class:`RandomResponseCase`, and precomputed direct harmonic FRF tensors.

Key capabilities:
    * **PSD curves**: :class:`PSDTable` (``log-log``, ``log-lin``, ``lin-log``,
      ``lin-lin`` interpolation; ``"zero"``, ``"clamp"``, ``"extrapolate"``
      out-of-range policies), :class:`WhiteNoisePSD`, and :class:`OctaveBandPSD`
      (IEC 61260 / ANSI S1.11 1/1 and 1/3 octave bands in dB).
    * **Damping models**: :class:`ConstantModalDamping`,
      :class:`ModalDampingList`, :class:`RayleighDamping`, and
      :class:`StructuralDamping`.
    * **Excitations & Multi-input Cross-PSD**: :class:`BaseExcitation`,
      :class:`NodalForceExcitation`, and :class:`CrossPSDMatrix` supporting
      uncorrelated, fully correlated, coherence + phase, co/quad (``C - i Q``),
      magnitude + phase, and explicit complex Hermitian PSD validation.
    * **Frequency grid builders**: uniform, logarithmic, Abaqus biased sweeps,
      and half-power adaptive modal clustering around ``f_n (1 +/- c * zeta_n)``.
    * **Solvers**: :class:`ModalRandomVibration` (with optional static residual
      flexibility) and :class:`DirectFRFRandomVibration`.
    * **Post-processing & Fatigue**: spectral moments (``lambda_0, 1, 2, 4``),
      zero-crossing/peak rates, irregularity factor, Segalman-Fulcher RMS von
      Mises stress, and Basquin S-N fatigue via :class:`NarrowBandRayleigh`,
      :class:`Steinberg3Band`, and :class:`Dirlik`.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import scipy.integrate as integrate
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import gamma as gamma_func

from wht_iga.dynamics.case_file import BaseKind, PSDType, RandomResponseCase
from wht_iga.dynamics.consistent_mass import EigenmodeResult
from wht_iga.dynamics.response_spectrum import VON_MISES_A, abaqus_frequency_sweep

InterpolationScale = Literal["log-log", "log-lin", "lin-log", "lin-lin"]
ExtrapolationPolicy = Literal["zero", "clamp", "extrapolate"]
ExcitationKind = Literal["acceleration", "velocity", "displacement"]
ResponseFrame = Literal["relative", "absolute"]
StressDefinition = Literal["amplitude", "range"]

#: Base-10 IEC 61260 exact octave ratio ``G = 10**(3/10)``.
IEC61260_BASE10_G: float = 10.0 ** (3.0 / 10.0)

_KIND_LAMBDA: dict[str, int] = {
    "acceleration": 0,
    "velocity": 1,
    "displacement": 2,
}

_STEINBERG_WEIGHTS: tuple[tuple[float, float], ...] = (
    (1.0, 0.6831),
    (2.0, 0.2710),
    (3.0, 0.0433),
)


# ---------------------------------------------------------------------------
# Internal utilities
# ---------------------------------------------------------------------------


def _as_1d_frequencies(freqs_hz: Sequence[float] | np.ndarray, *, allow_zero: bool = True) -> np.ndarray:
    """Validate and return a 1D float array of strictly increasing frequencies."""
    f = np.atleast_1d(np.asarray(freqs_hz, dtype=float))
    if f.ndim != 1 or f.size == 0:
        raise ValueError("frequencies_hz must be a non-empty 1D array.")
    if not np.all(np.isfinite(f)):
        raise ValueError("frequencies_hz must contain only finite values.")
    if allow_zero:
        if np.any(f < 0.0):
            raise ValueError("frequencies_hz must be non-negative.")
    elif np.any(f <= 0.0):
        raise ValueError("frequencies_hz must be strictly positive.")
    if f.size > 1 and np.any(np.diff(f) <= 0.0):
        raise ValueError("frequencies_hz must be strictly increasing.")
    return f


def _validate_hermitian_psd_slice(mat: np.ndarray, freq_index: int, freq_hz: float, tol: float) -> None:
    """Verify that a single ``(n_in, n_in)`` CPSD slice is Hermitian and PSD."""
    herm_err = float(np.max(np.abs(mat - mat.conj().T)))
    scale = max(1.0, float(np.linalg.norm(mat, ord="fro")))
    if herm_err > tol * scale:
        raise ValueError(
            f"Cross-PSD matrix at frequency index {freq_index} ({freq_hz:.6g} Hz) is not Hermitian: "
            f"max|S - S^H| = {herm_err:.3e} > tol*scale = {tol * scale:.3e}."
        )
    sym = 0.5 * (mat + mat.conj().T)
    eigvals = np.linalg.eigvalsh(sym)
    min_eig = float(np.min(eigvals))
    max_eig = max(1.0, float(np.max(np.abs(eigvals))))
    if min_eig < -tol * max_eig:
        raise ValueError(
            f"Cross-PSD matrix at frequency index {freq_index} ({freq_hz:.6g} Hz) is not positive "
            f"semi-definite: min eigenvalue = {min_eig:.3e} < -{tol * max_eig:.3e}."
        )


def _direction_to_vector(dof_or_direction: int | str | Sequence[float] | np.ndarray, *, one_based: bool = False) -> np.ndarray:
    """Resolve a translational DOF index, axis name, or 3-vector to a unit 3D vector."""
    if isinstance(dof_or_direction, (int, np.integer)):
        idx = int(dof_or_direction)
        if one_based:
            if idx not in (1, 2, 3):
                raise ValueError(f"1-based translational DOF must be 1, 2, or 3, got {idx}.")
            idx -= 1
        elif idx not in (0, 1, 2):
            raise ValueError(
                f"0-based translational axis index must be 0, 1, or 2 (got {idx}); "
                "pass 'x'/'y'/'z', a 3-vector, or one_based=True for 1-based DOFs."
            )
        vec = np.zeros(3, dtype=float)
        vec[idx] = 1.0
        return vec

    if isinstance(dof_or_direction, str):
        token = dof_or_direction.strip().lower()
        mapping = {"x": 0, "ux": 0, "1": 0, "y": 1, "uy": 1, "2": 1, "z": 2, "uz": 2, "3": 2}
        if token not in mapping:
            raise ValueError(f"Unknown direction string {dof_or_direction!r}; expected 'x', 'y', or 'z'.")
        vec = np.zeros(3, dtype=float)
        vec[mapping[token]] = 1.0
        return vec

    vec = np.asarray(dof_or_direction, dtype=float).reshape(-1)
    if vec.size != 3:
        raise ValueError(f"Direction vector must have length 3, got {vec.size}.")
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0:
        raise ValueError("Direction vector must be non-zero.")
    return vec / norm


# ---------------------------------------------------------------------------
# 1. PSD Curves
# ---------------------------------------------------------------------------


class PSDCurve(ABC):
    """Abstract base class for frequency-dependent Power Spectral Density functions."""

    @abstractmethod
    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        """Evaluate the PSD ordinate at ``freqs_hz`` (in Hz).

        Args:
            freqs_hz: Scalar or array of non-negative frequencies in Hz.

        Returns:
            Complex or real array matching the shape of ``freqs_hz`` (1D if scalar).
        """

    def __call__(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        """Alias for :meth:`evaluate`."""
        return self.evaluate(freqs_hz)


class PSDTable(PSDCurve):
    """Tabulated PSD curve supporting log/linear interpolation and boundary policies.

    Args:
        freqs_hz: Strictly increasing table frequencies in Hz.
        values: Tabulated PSD ordinates matching ``freqs_hz``. May be real or
            complex (or combine with ``imag_values``).
        scale: Interpolation scale in ``{"log-log", "log-lin", "lin-log", "lin-lin"}``.
        extrapolation: Out-of-range policy in ``{"zero", "clamp", "extrapolate"}``.
            Alias ``"constant"`` is accepted for ``"clamp"``.
        imag_values: Optional imaginary-part table added to real ``values``.
        name: Optional human-readable label.
    """

    def __init__(
        self,
        freqs_hz: Sequence[float] | np.ndarray,
        values: Sequence[float | complex] | np.ndarray,
        scale: InterpolationScale = "log-log",
        extrapolation: ExtrapolationPolicy = "zero",
        imag_values: Sequence[float] | np.ndarray | None = None,
        name: str | None = None,
    ) -> None:
        self.name = name
        scale_norm = scale.strip().lower()
        if scale_norm not in ("log-log", "log-lin", "lin-log", "lin-lin"):
            raise ValueError(f"scale must be 'log-log', 'log-lin', 'lin-log', or 'lin-lin', got {scale!r}.")
        self.scale: InterpolationScale = scale_norm  # type: ignore[assignment]

        extrap_norm = extrapolation.strip().lower()
        if extrap_norm == "constant":
            extrap_norm = "clamp"
        if extrap_norm not in ("zero", "clamp", "extrapolate"):
            raise ValueError(f"extrapolation must be 'zero', 'clamp', or 'extrapolate', got {extrapolation!r}.")
        self.extrapolation: ExtrapolationPolicy = extrap_norm  # type: ignore[assignment]

        x_log, y_log = self.scale.split("-")
        self._x_is_log = x_log == "log"
        self._y_is_log = y_log == "log"

        f_arr = _as_1d_frequencies(freqs_hz, allow_zero=not self._x_is_log)
        v_raw = np.asarray(values)
        if v_raw.ndim != 1 or v_raw.size != f_arr.size:
            raise ValueError("values must be a 1D array with the same length as freqs_hz.")
        if imag_values is not None:
            im_arr = np.asarray(imag_values, dtype=float)
            if im_arr.shape != f_arr.shape:
                raise ValueError("imag_values must match freqs_hz shape.")
            if np.iscomplexobj(v_raw):
                raise ValueError("Pass either complex values or separate imag_values, not both.")
            v_arr = v_raw.astype(float) + 1j * im_arr
        else:
            v_arr = np.asarray(v_raw, dtype=complex if np.iscomplexobj(v_raw) else float)

        if not np.all(np.isfinite(v_arr)):
            raise ValueError("PSDTable values must be finite.")

        if self._y_is_log:
            if np.iscomplexobj(v_arr):
                re = v_arr.real
                im = v_arr.imag
                if np.any(re <= 0.0):
                    raise ValueError("log-y interpolation requires strictly positive real parts.")
                if not (np.all(im == 0.0) or np.all(im > 0.0) or np.all(im < 0.0)):
                    raise ValueError("log-y interpolation requires imaginary parts to be all zero or strictly one-signed.")
            elif np.any(v_arr <= 0.0):
                raise ValueError("log-y interpolation requires strictly positive PSD values.")

        self.freqs_hz: np.ndarray = f_arr
        self.values: np.ndarray = v_arr

    @property
    def is_complex(self) -> bool:
        """Whether the table carries non-zero imaginary ordinates."""
        return bool(np.iscomplexobj(self.values) and np.any(self.values.imag != 0.0))

    def _interp_real_component(self, f_eval: np.ndarray, y_tab: np.ndarray, log_y: bool) -> np.ndarray:
        """Interpolate and extrapolate a single real-valued ordinate vector."""
        x_tab = np.log(self.freqs_hz) if self._x_is_log else self.freqs_hz
        y_trans = np.log(y_tab) if log_y else y_tab

        out = np.zeros(f_eval.shape, dtype=float)
        f_lo, f_hi = float(self.freqs_hz[0]), float(self.freqs_hz[-1])
        in_mask = (f_eval >= f_lo) & (f_eval <= f_hi)
        lo_mask = f_eval < f_lo
        hi_mask = f_eval > f_hi

        if np.any(in_mask):
            x_in = np.log(f_eval[in_mask]) if self._x_is_log else f_eval[in_mask]
            y_in = np.interp(x_in, x_tab, y_trans)
            out[in_mask] = np.exp(y_in) if log_y else y_in

        if self.extrapolation == "clamp":
            out[lo_mask] = float(y_tab[0])
            out[hi_mask] = float(y_tab[-1])
        elif self.extrapolation == "extrapolate":
            if np.any(lo_mask):
                if self._x_is_log and np.any(f_eval[lo_mask] <= 0.0):
                    raise ValueError("Log-frequency extrapolation requires strictly positive evaluation frequencies.")
                if x_tab.size == 1:
                    out[lo_mask] = float(y_tab[0])
                else:
                    slope_lo = (y_trans[1] - y_trans[0]) / (x_tab[1] - x_tab[0])
                    x_lo = np.log(f_eval[lo_mask]) if self._x_is_log else f_eval[lo_mask]
                    y_lo = y_trans[0] + slope_lo * (x_lo - x_tab[0])
                    out[lo_mask] = np.exp(y_lo) if log_y else y_lo
            if np.any(hi_mask):
                if x_tab.size == 1:
                    out[hi_mask] = float(y_tab[-1])
                else:
                    slope_hi = (y_trans[-1] - y_trans[-2]) / (x_tab[-1] - x_tab[-2])
                    x_hi = np.log(f_eval[hi_mask]) if self._x_is_log else f_eval[hi_mask]
                    y_hi = y_trans[-1] + slope_hi * (x_hi - x_tab[-1])
                    out[hi_mask] = np.exp(y_hi) if log_y else y_hi

        return out

    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        """Evaluate the tabulated PSD at ``freqs_hz``."""
        f_eval = np.asarray(freqs_hz, dtype=float)
        scalar_input = f_eval.ndim == 0
        f_flat = np.atleast_1d(f_eval).reshape(-1)
        if np.any(f_flat < 0.0) or not np.all(np.isfinite(f_flat)):
            raise ValueError("Evaluation frequencies must be finite and non-negative.")

        if np.iscomplexobj(self.values):
            re = self._interp_real_component(f_flat, self.values.real, self._y_is_log)
            im_vals = self.values.imag
            if np.all(im_vals == 0.0):
                out_flat = re.astype(float)
            elif self._y_is_log:
                sign = 1.0 if im_vals[0] > 0.0 else -1.0
                im = sign * self._interp_real_component(f_flat, np.abs(im_vals), True)
                out_flat = re + 1j * im
            else:
                im = self._interp_real_component(f_flat, im_vals, False)
                out_flat = re + 1j * im
        else:
            out_flat = self._interp_real_component(f_flat, self.values, self._y_is_log)

        if scalar_input:
            return out_flat.reshape(())
        return out_flat.reshape(f_eval.shape)


class WhiteNoisePSD(PSDCurve):
    """Band-limited constant (white-noise) one-sided PSD.

    Args:
        value: Constant PSD level ``W0 >= 0`` (in ``unit^2 / Hz``).
        f_min: Lower active frequency bound in Hz (default ``0.0``).
        f_max: Upper active frequency bound in Hz (default ``np.inf``).
    """

    def __init__(self, value: float, f_min: float = 0.0, f_max: float = float("inf")) -> None:
        val = float(value)
        flo = float(f_min)
        fhi = float(f_max)
        if not math.isfinite(val) or val < 0.0:
            raise ValueError(f"WhiteNoisePSD value must be finite and >= 0, got {value}.")
        if not math.isfinite(flo) or flo < 0.0:
            raise ValueError(f"WhiteNoisePSD f_min must be finite and >= 0, got {f_min}.")
        if fhi <= flo:
            raise ValueError(f"WhiteNoisePSD requires f_max > f_min, got [{flo}, {fhi}].")
        self.value = val
        self.f_min = flo
        self.f_max = fhi

    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        """Return ``value`` inside ``[f_min, f_max]`` and ``0.0`` outside."""
        f_arr = np.asarray(freqs_hz, dtype=float)
        if np.any(f_arr < 0.0) or not np.all(np.isfinite(f_arr)):
            raise ValueError("Evaluation frequencies must be finite and non-negative.")
        active = (f_arr >= self.f_min) & (f_arr <= self.f_max)
        return np.where(active, self.value, 0.0)


def iec61260_center_frequencies(
    f_min: float,
    f_max: float,
    fraction: int = 3,
    f_ref: float = 1000.0,
    base: Literal[10, 2] = 10,
) -> np.ndarray:
    """Compute exact IEC 61260 / ANSI S1.11 1/1 or 1/3 octave center frequencies.

    Args:
        f_min: Minimum center frequency bound in Hz (inclusive within tolerance).
        f_max: Maximum center frequency bound in Hz (inclusive within tolerance).
        fraction: Band division ``b`` (``1`` for 1/1 octave, ``3`` for 1/3 octave).
        f_ref: Reference center frequency in Hz (default ``1000.0`` Hz).
        base: Octave ratio base: ``10`` for IEC 61260 ``G = 10**(3/10)`` (default),
            or ``2`` for exact binary octave ``G = 2``.

    Returns:
        Sorted 1D array of exact center frequencies in Hz.
    """
    if f_min <= 0.0 or f_max < f_min:
        raise ValueError("Require 0 < f_min <= f_max.")
    if fraction < 1:
        raise ValueError(f"fraction must be >= 1, got {fraction}.")
    g = IEC61260_BASE10_G if base == 10 else 2.0
    if base not in (10, 2):
        raise ValueError(f"base must be 10 or 2, got {base!r}.")
    step = g ** (1.0 / float(fraction))
    k_lo = math.ceil(math.log(f_min / f_ref, step) - 1e-12)
    k_hi = math.floor(math.log(f_max / f_ref, step) + 1e-12)
    ks = np.arange(k_lo, k_hi + 1, dtype=float)
    return f_ref * (step**ks)


def iec61260_band_edges(
    center_freqs_hz: Sequence[float] | np.ndarray,
    fraction: int = 3,
    base: Literal[10, 2] = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute lower and upper IEC 61260 band-edge frequencies for center frequencies.

    Args:
        center_freqs_hz: Center frequencies in Hz.
        fraction: Octave fraction ``b`` (``1`` for 1/1 octave, ``3`` for 1/3 octave).
        base: Octave ratio base (``10`` or ``2``).

    Returns:
        Tuple ``(f_lower, f_upper)`` of arrays matching ``center_freqs_hz``.
    """
    fc = _as_1d_frequencies(center_freqs_hz, allow_zero=False)
    if fraction < 1:
        raise ValueError(f"fraction must be >= 1, got {fraction}.")
    g = IEC61260_BASE10_G if base == 10 else 2.0
    if base not in (10, 2):
        raise ValueError(f"base must be 10 or 2, got {base!r}.")
    half = g ** (1.0 / (2.0 * float(fraction)))
    return fc / half, fc * half


class OctaveBandPSD(PSDCurve):
    """IEC 61260 / ANSI S1.11 1/1 or 1/3 octave-band dB PSD definition.

    Converts decibel levels into linear PSD values (``unit^2 / Hz``) using either:

    * ``level_kind="psd_db"`` (default): dB specifies PSD level directly relative
      to ``ref_psd`` (so ``0 dB = ref_psd``, default ``1.0 unit^2/Hz``):
      ``S(f_c) = ref_psd * 10**(L_dB / 10)``.
    * ``level_kind="band_db"``: dB specifies band-integrated mean-square / SPL
      power relative to ``ref_psd`` (in ``unit^2``), dividing by the exact IEC
      61260 bandwidth ``df = f_upper - f_lower`` so that integrating ``S(f)``
      over the band equals ``ref_psd * 10**(L_dB / 10)``.

    Args:
        center_freqs_hz: Strictly increasing positive center frequencies in Hz.
        levels_db: Decibel levels corresponding to ``center_freqs_hz``.
        fraction: Octave fraction (``1`` for 1/1 octave, ``3`` for 1/3 octave).
        ref_psd: Reference level (``1.0`` so ``0 dB = 1 unit^2/Hz`` for
            ``"psd_db"``, or reference mean-square ``p_ref**2`` for ``"band_db"``).
        level_kind: ``"psd_db"`` or ``"band_db"``.
        base: IEC 61260 octave base ``10`` (default, ``G = 10**0.3``) or ``2``.
        mode: ``"continuous"`` (interpolates center-frequency table as
            :class:`PSDTable`) or ``"band_constant"`` (piecewise-constant PSD
            across each octave band ``[f_lo, f_hi]``).
        scale: Interpolation scale used in ``"continuous"`` mode and
            :meth:`to_psd_table`.
        extrapolation: Boundary policy used in ``"continuous"`` mode.
    """

    def __init__(
        self,
        center_freqs_hz: Sequence[float] | np.ndarray,
        levels_db: Sequence[float] | np.ndarray,
        fraction: int = 3,
        ref_psd: float = 1.0,
        level_kind: Literal["psd_db", "band_db"] = "psd_db",
        base: Literal[10, 2] = 10,
        mode: Literal["continuous", "band_constant"] = "continuous",
        scale: InterpolationScale = "log-log",
        extrapolation: ExtrapolationPolicy = "zero",
    ) -> None:
        fc = _as_1d_frequencies(center_freqs_hz, allow_zero=False)
        db = np.asarray(levels_db, dtype=float).reshape(-1)
        if db.size != fc.size or not np.all(np.isfinite(db)):
            raise ValueError("levels_db must be a finite 1D array matching center_freqs_hz.")
        if fraction not in (1, 3) and fraction < 1:
            raise ValueError(f"fraction must be >= 1 (typically 1 or 3), got {fraction}.")
        ref = float(ref_psd)
        if not math.isfinite(ref) or ref <= 0.0:
            raise ValueError(f"ref_psd must be positive and finite, got {ref_psd}.")
        lkind = level_kind.strip().lower()
        if lkind not in ("psd_db", "band_db"):
            raise ValueError(f"level_kind must be 'psd_db' or 'band_db', got {level_kind!r}.")
        m_norm = mode.strip().lower()
        if m_norm not in ("continuous", "band_constant"):
            raise ValueError(f"mode must be 'continuous' or 'band_constant', got {mode!r}.")

        f_lo, f_hi = iec61260_band_edges(fc, fraction=int(fraction), base=base)
        bandwidths = f_hi - f_lo
        power = ref * (10.0 ** (db / 10.0))
        psd_vals = power if lkind == "psd_db" else power / bandwidths

        self.center_freqs_hz = fc
        self.levels_db = db
        self.fraction = int(fraction)
        self.ref_psd = ref
        self.level_kind = lkind
        self.base = base
        self.mode = m_norm
        self.f_lower_hz = f_lo
        self.f_upper_hz = f_hi
        self.bandwidths_hz = bandwidths
        self.psd_values = psd_vals
        self._table = PSDTable(fc, psd_vals, scale=scale, extrapolation=extrapolation)

    @classmethod
    def from_ansi_bands(
        cls,
        band_numbers: Sequence[int] | np.ndarray,
        levels_db: Sequence[float] | np.ndarray,
        fraction: int = 3,
        ref_psd: float = 1.0,
        level_kind: Literal["psd_db", "band_db"] = "psd_db",
        **kwargs: Any,
    ) -> OctaveBandPSD:
        """Build an :class:`OctaveBandPSD` from ANSI S1.11 band indices.

        For 1/3-octave (``fraction=3``), band ``30`` is ``1000 Hz`` (``10**(B/10)``).
        For 1/1-octave (``fraction=1``), band ``10`` is ``1000 Hz``
        (``1000 * 10**(0.3 * (B - 10))``).
        """
        bands = np.asarray(band_numbers, dtype=int).reshape(-1)
        if bands.size == 0:
            raise ValueError("band_numbers must be non-empty.")
        if fraction == 3:
            fc = 1000.0 * (IEC61260_BASE10_G ** ((bands - 30.0) / 3.0))
        elif fraction == 1:
            fc = 1000.0 * (IEC61260_BASE10_G ** (bands - 10.0))
        else:
            raise ValueError("from_ansi_bands supports fraction=1 or fraction=3.")
        return cls(
            center_freqs_hz=fc,
            levels_db=levels_db,
            fraction=fraction,
            ref_psd=ref_psd,
            level_kind=level_kind,
            **kwargs,
        )

    def to_psd_table(self) -> PSDTable:
        """Convert the octave-band PSD at center frequencies to a :class:`PSDTable`."""
        return self._table

    def band_mean_squares(self) -> np.ndarray:
        """Return the integrated mean-square value ``S(f_c) * df`` per band."""
        return self.psd_values * self.bandwidths_hz

    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        """Evaluate the octave-band PSD at ``freqs_hz``."""
        if self.mode == "continuous":
            return self._table.evaluate(freqs_hz)

        f_eval = np.asarray(freqs_hz, dtype=float)
        scalar_input = f_eval.ndim == 0
        f_flat = np.atleast_1d(f_eval).reshape(-1)
        if np.any(f_flat < 0.0) or not np.all(np.isfinite(f_flat)):
            raise ValueError("Evaluation frequencies must be finite and non-negative.")

        out = np.zeros(f_flat.shape, dtype=float)
        # Assign band value where f lies in [f_lower, f_upper]; shared edges resolve consistently
        for k in range(self.center_freqs_hz.size):
            if k == self.center_freqs_hz.size - 1:
                mask = (f_flat >= self.f_lower_hz[k]) & (f_flat <= self.f_upper_hz[k])
            else:
                mask = (f_flat >= self.f_lower_hz[k]) & (f_flat < self.f_upper_hz[k])
            out[mask] = self.psd_values[k]

        if scalar_input:
            return out.reshape(())
        return out.reshape(f_eval.shape)


# ---------------------------------------------------------------------------
# 2. Damping Models
# ---------------------------------------------------------------------------


class DampingModel(ABC):
    """Abstract base class for modal/structural damping models."""

    @abstractmethod
    def equivalent_viscous_ratios(self, omega_modes: np.ndarray) -> np.ndarray:
        """Return equivalent viscous damping ratios ``zeta_k`` per mode."""

    @abstractmethod
    def modal_transfer_function(
        self,
        omega_modes: np.ndarray,
        omega_grid: np.ndarray,
        modal_masses: np.ndarray | None = None,
    ) -> np.ndarray:
        """Evaluate modal FRF ``H_k(omega)`` of shape ``(n_modes, n_freq)``."""


def _as_damping_model(damping: DampingModel | float | Sequence[float] | np.ndarray) -> DampingModel:
    """Normalize a damping specification into a :class:`DampingModel` instance."""
    if isinstance(damping, DampingModel):
        return damping
    if isinstance(damping, (int, float, np.floating, np.integer)):
        return ConstantModalDamping(float(damping))
    return ModalDampingList(damping)


class ConstantModalDamping(DampingModel):
    """Uniform viscous modal damping ratio ``zeta`` across all modes."""

    def __init__(self, zeta: float) -> None:
        z = float(zeta)
        if not math.isfinite(z) or z < 0.0:
            raise ValueError(f"Modal damping ratio must be finite and >= 0, got {zeta}.")
        self.zeta = z

    def equivalent_viscous_ratios(self, omega_modes: np.ndarray) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        return np.full(w.shape, self.zeta, dtype=float)

    def modal_transfer_function(
        self,
        omega_modes: np.ndarray,
        omega_grid: np.ndarray,
        modal_masses: np.ndarray | None = None,
    ) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        wg = np.asarray(omega_grid, dtype=float).reshape(-1)
        m = np.ones_like(w) if modal_masses is None else np.asarray(modal_masses, dtype=float).reshape(-1)
        denom = m[:, None] * (w[:, None] ** 2 - wg[None, :] ** 2 + 2j * self.zeta * w[:, None] * wg[None, :])
        return 1.0 / denom


class ModalDampingList(DampingModel):
    """Explicit per-mode viscous damping ratios ``[zeta_1, ..., zeta_m]``."""

    def __init__(self, zetas: Sequence[float] | np.ndarray) -> None:
        z = np.asarray(zetas, dtype=float).reshape(-1)
        if z.size == 0 or not np.all(np.isfinite(z)) or np.any(z < 0.0):
            raise ValueError("ModalDampingList requires a non-empty array of finite zetas >= 0.")
        self.zetas = z

    def equivalent_viscous_ratios(self, omega_modes: np.ndarray) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        if self.zetas.size != w.size:
            raise ValueError(f"ModalDampingList length ({self.zetas.size}) must match n_modes ({w.size}).")
        return self.zetas.copy()

    def modal_transfer_function(
        self,
        omega_modes: np.ndarray,
        omega_grid: np.ndarray,
        modal_masses: np.ndarray | None = None,
    ) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        z = self.equivalent_viscous_ratios(w)
        wg = np.asarray(omega_grid, dtype=float).reshape(-1)
        m = np.ones_like(w) if modal_masses is None else np.asarray(modal_masses, dtype=float).reshape(-1)
        denom = m[:, None] * (w[:, None] ** 2 - wg[None, :] ** 2 + 2j * z[:, None] * w[:, None] * wg[None, :])
        return 1.0 / denom


class RayleighDamping(DampingModel):
    """Proportional Rayleigh damping ``C = alpha * M + beta * K``.

    Gives modal damping ratio ``zeta_k = alpha / (2 * omega_k) + beta * omega_k / 2``.
    """

    def __init__(self, alpha: float = 0.0, beta: float = 0.0) -> None:
        a = float(alpha)
        b = float(beta)
        if not (math.isfinite(a) and math.isfinite(b)) or a < 0.0 or b < 0.0:
            raise ValueError(f"RayleighDamping requires finite alpha >= 0 and beta >= 0, got ({alpha}, {beta}).")
        self.alpha = a
        self.beta = b

    def equivalent_viscous_ratios(self, omega_modes: np.ndarray) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        if np.any(w <= 0.0):
            raise ValueError("RayleighDamping requires strictly positive modal circular frequencies.")
        return 0.5 * self.alpha / w + 0.5 * self.beta * w

    def modal_transfer_function(
        self,
        omega_modes: np.ndarray,
        omega_grid: np.ndarray,
        modal_masses: np.ndarray | None = None,
    ) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        z = self.equivalent_viscous_ratios(w)
        wg = np.asarray(omega_grid, dtype=float).reshape(-1)
        m = np.ones_like(w) if modal_masses is None else np.asarray(modal_masses, dtype=float).reshape(-1)
        denom = m[:, None] * (w[:, None] ** 2 - wg[None, :] ** 2 + 2j * z[:, None] * w[:, None] * wg[None, :])
        return 1.0 / denom


class StructuralDamping(DampingModel):
    """Hysteretic / structural loss-factor damping ``K* = K * (1 + i * g_loss)``.

    Uses the exact frequency-independent hysteretic modal dynamic stiffness
    ``m_k * (omega_k**2 - omega**2 + i * g_k * omega_k**2)``.
    """

    def __init__(self, g_loss: float | Sequence[float] | np.ndarray) -> None:
        g = np.asarray(g_loss, dtype=float)
        if not np.all(np.isfinite(g)) or np.any(g < 0.0):
            raise ValueError(f"StructuralDamping loss factor g_loss must be finite and >= 0, got {g_loss!r}.")
        self.g_loss = float(g) if g.ndim == 0 else g.reshape(-1)

    def _loss_array(self, n_modes: int) -> np.ndarray:
        if isinstance(self.g_loss, float):
            return np.full(n_modes, self.g_loss, dtype=float)
        if self.g_loss.size != n_modes:
            raise ValueError(f"StructuralDamping g_loss length ({self.g_loss.size}) must match n_modes ({n_modes}).")
        return self.g_loss

    def equivalent_viscous_ratios(self, omega_modes: np.ndarray) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        return 0.5 * self._loss_array(w.size)

    def modal_transfer_function(
        self,
        omega_modes: np.ndarray,
        omega_grid: np.ndarray,
        modal_masses: np.ndarray | None = None,
    ) -> np.ndarray:
        w = np.asarray(omega_modes, dtype=float).reshape(-1)
        g = self._loss_array(w.size)
        wg = np.asarray(omega_grid, dtype=float).reshape(-1)
        m = np.ones_like(w) if modal_masses is None else np.asarray(modal_masses, dtype=float).reshape(-1)
        denom = m[:, None] * (w[:, None] ** 2 - wg[None, :] ** 2 + 1j * g[:, None] * (w[:, None] ** 2))
        return 1.0 / denom


# ---------------------------------------------------------------------------
# 3. Excitations & Multi-input Cross-PSD Matrix
# ---------------------------------------------------------------------------


def _as_psd_curve(psd: PSDCurve | float | Callable[[np.ndarray], np.ndarray]) -> PSDCurve:
    """Wrap scalar or callable PSD inputs into a :class:`PSDCurve`."""
    if isinstance(psd, PSDCurve):
        return psd
    if isinstance(psd, (int, float, np.floating, np.integer)):
        return WhiteNoisePSD(float(psd))
    if callable(psd):
        return _CallablePSD(psd)
    raise TypeError(f"psd must be a PSDCurve, float, or callable, got {type(psd).__name__}.")


class _CallablePSD(PSDCurve):
    """Lightweight wrapper for user-supplied PSD callables ``f(freqs) -> values``."""

    def __init__(self, fn: Callable[[np.ndarray], np.ndarray | float]) -> None:
        self._fn = fn

    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray | float) -> np.ndarray:
        f_arr = np.asarray(freqs_hz, dtype=float)
        vals = np.asarray(self._fn(f_arr))
        if vals.ndim == 0:
            return np.full(f_arr.shape, vals.item())
        return vals.reshape(f_arr.shape)


class ExcitationBase(ABC):
    """Abstract base class for random vibration input excitations."""

    psd: PSDCurve

    @abstractmethod
    def amplitude_scale(self) -> float:
        """Linear scale factor multiplying the excitation amplitude (PSD scales by ``scale**2``)."""


@dataclass(frozen=True)
class BaseExcitation(ExcitationBase):
    """Translational base-motion PSD excitation.

    Args:
        dof_or_direction: Global translational axis index ``0, 1, 2`` (or ``1, 2, 3``
            when ``one_based=True``), string ``'x' | 'y' | 'z'``, or 3D vector
            ``(dx, dy, dz)``.
        psd: Auto-PSD curve (or constant float) in units corresponding to ``kind``
            divided by ``gravity_scale``.
        gravity_scale: Linear multiplier on base motion amplitude (e.g., ``9.81``
            to convert a ``g^2/Hz`` acceleration PSD into ``(m/s^2)^2/Hz``).
        kind: Kinematic type of the input base motion: ``"acceleration"``,
            ``"velocity"``, or ``"displacement"``.
        one_based: Interpret integer ``dof_or_direction`` as 1-based (``1, 2, 3``).
    """

    dof_or_direction: int | str | Sequence[float] | np.ndarray
    psd: PSDCurve
    gravity_scale: float = 1.0
    kind: ExcitationKind = "acceleration"
    one_based: bool = False

    def __init__(
        self,
        dof_or_direction: int | str | Sequence[float] | np.ndarray,
        psd: PSDCurve | float | Callable[[np.ndarray], np.ndarray],
        gravity_scale: float = 1.0,
        kind: ExcitationKind = "acceleration",
        one_based: bool = False,
    ) -> None:
        vec = _direction_to_vector(dof_or_direction, one_based=one_based)
        psd_obj = _as_psd_curve(psd)
        g = float(gravity_scale)
        if not math.isfinite(g) or g <= 0.0:
            raise ValueError(f"gravity_scale must be finite and positive, got {gravity_scale}.")
        kind_norm = kind.strip().lower()
        if kind_norm not in _KIND_LAMBDA:
            raise ValueError(f"kind must be 'acceleration', 'velocity', or 'displacement', got {kind!r}.")
        object.__setattr__(self, "dof_or_direction", dof_or_direction)
        object.__setattr__(self, "direction", vec)
        object.__setattr__(self, "psd", psd_obj)
        object.__setattr__(self, "gravity_scale", g)
        object.__setattr__(self, "kind", kind_norm)
        object.__setattr__(self, "one_based", bool(one_based))
        object.__setattr__(self, "lambda_exp", _KIND_LAMBDA[kind_norm])

    def amplitude_scale(self) -> float:
        return self.gravity_scale


@dataclass(frozen=True)
class NodalForceExcitation(ExcitationBase):
    """Spatially distributed nodal force PSD excitation ``f(t) = force_vector * p(t)``.

    Args:
        force_vector: Real or complex spatial load pattern of shape ``(n_dofs,)``.
        psd: Auto-PSD curve (or constant float) of the scalar multiplier ``p(t)``.
        scale: Optional linear multiplier on ``force_vector``.
    """

    force_vector: np.ndarray
    psd: PSDCurve
    scale: float = 1.0

    def __init__(
        self,
        force_vector: Sequence[float | complex] | np.ndarray,
        psd: PSDCurve | float | Callable[[np.ndarray], np.ndarray],
        scale: float = 1.0,
    ) -> None:
        fvec = np.asarray(force_vector)
        if fvec.ndim != 1 or fvec.size == 0:
            raise ValueError("force_vector must be a non-empty 1D array.")
        fvec = fvec.astype(complex if np.iscomplexobj(fvec) else float, copy=True)
        sc = float(scale)
        if not math.isfinite(sc) or sc <= 0.0:
            raise ValueError(f"scale must be finite and positive, got {scale}.")
        object.__setattr__(self, "force_vector", fvec)
        object.__setattr__(self, "psd", _as_psd_curve(psd))
        object.__setattr__(self, "scale", sc)

    def amplitude_scale(self) -> float:
        return self.scale


class CrossPSDMatrix:
    """Multi-input Hermitian positive semi-definite Cross-PSD matrix ``S_in(f)``.

    Evaluates to a complex array of shape ``(n_in, n_in, n_f)`` satisfying:
    ``S_in[i, j, k] = conj(S_in[j, i, k])`` and ``lambda_min(S_in[:, :, k]) >= -tol``.

    Supports:
        1. Uncorrelated (diagonal) inputs (:meth:`uncorrelated`).
        2. Fully correlated inputs with optional phases (:meth:`fully_correlated`).
        3. Coherence ``gamma_ij^2`` + phase (:meth:`from_coherence`).
        4. Co-spectrum / quad-spectrum ``C - i Q`` (:meth:`from_co_quad`).
        5. Magnitude + phase (:meth:`from_magnitude_phase`).
        6. Explicit cross-term curves (:meth:`from_cross_terms`) or prebuilt
           arrays (:meth:`from_array`).
    """

    def __init__(
        self,
        excitations: Sequence[ExcitationBase] | None,
        evaluator: Callable[[np.ndarray], np.ndarray],
        n_inputs: int,
        tol: float = 1e-10,
    ) -> None:
        if n_inputs < 1:
            raise ValueError("CrossPSDMatrix requires at least one input.")
        if tol < 0.0 or not math.isfinite(tol):
            raise ValueError(f"tol must be finite and >= 0, got {tol}.")
        self.excitations = tuple(excitations) if excitations is not None else ()
        self._evaluator = evaluator
        self.n_inputs = int(n_inputs)
        self.tol = float(tol)

    @classmethod
    def uncorrelated(
        cls,
        excitations: Sequence[ExcitationBase],
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Create a diagonal (mutually uncorrelated) input Cross-PSD matrix."""
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_in, n_f = len(excs), freqs.size
            s = np.zeros((n_in, n_in, n_f), dtype=complex)
            for i, exc in enumerate(excs):
                sc = exc.amplitude_scale()
                s[i, i, :] = (sc * sc) * np.real(exc.psd.evaluate(freqs))
            return s

        return cls(excs, _eval, len(excs), tol=tol)

    @classmethod
    def fully_correlated(
        cls,
        excitations: Sequence[ExcitationBase],
        phases: Sequence[float] | np.ndarray | None = None,
        phase_unit: Literal["rad", "deg"] = "rad",
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Create a rank-1 fully correlated input Cross-PSD matrix.

        ``S_ij(f) = sqrt(S_ii(f) * S_jj(f)) * exp(i * (phi_i - phi_j))``.
        """
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")
        n_in = len(excs)
        if phases is None:
            phi = np.zeros(n_in, dtype=float)
        else:
            phi = np.asarray(phases, dtype=float).reshape(-1)
            if phi.size != n_in:
                raise ValueError(f"phases length ({phi.size}) must equal n_inputs ({n_in}).")
            if phase_unit == "deg":
                phi = np.deg2rad(phi)
            elif phase_unit != "rad":
                raise ValueError(f"phase_unit must be 'rad' or 'deg', got {phase_unit!r}.")
        phasors = np.exp(1j * phi)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_f = freqs.size
            amps = np.empty((n_in, n_f), dtype=complex)
            for i, exc in enumerate(excs):
                sc = exc.amplitude_scale()
                s_ii = (sc * sc) * np.real(exc.psd.evaluate(freqs))
                if np.any(s_ii < -tol * max(1.0, float(np.max(np.abs(s_ii))))):
                    raise ValueError(f"Auto-PSD of excitation {i} contains negative values.")
                amps[i, :] = phasors[i] * np.sqrt(np.maximum(s_ii, 0.0))
            return np.einsum("if,jf->ijf", amps, amps.conj())

        return cls(excs, _eval, n_in, tol=tol)

    @classmethod
    def from_coherence(
        cls,
        excitations: Sequence[ExcitationBase],
        coherence: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable],
        phase: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable] | None = None,
        phase_unit: Literal["rad", "deg"] = "rad",
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Build a Cross-PSD matrix from ordinary squared coherence ``gamma_ij^2`` and phase.

        ``|S_ij(f)| = sqrt(gamma_ij^2(f) * S_ii(f) * S_jj(f))``,
        ``S_ij(f) = |S_ij(f)| * exp(i * theta_ij(f))``.
        If ``coherence`` or ``phase`` is given as a dict keyed by ``(i, j)``,
        symmetric/Hermitian ``(j, i)`` conjugates are filled automatically.
        """
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")
        n_in = len(excs)
        if phase_unit not in ("rad", "deg"):
            raise ValueError(f"phase_unit must be 'rad' or 'deg', got {phase_unit!r}.")
        deg_factor = math.pi / 180.0 if phase_unit == "deg" else 1.0

        def _eval_entry(
            spec: float | PSDCurve | Callable | np.ndarray, freqs: np.ndarray
        ) -> np.ndarray:
            if isinstance(spec, PSDCurve):
                return np.asarray(spec.evaluate(freqs), dtype=float)
            if callable(spec):
                val = np.asarray(spec(freqs), dtype=float)
                return np.full(freqs.shape, float(val)) if val.ndim == 0 else val.reshape(freqs.shape)
            arr = np.asarray(spec, dtype=float)
            if arr.ndim == 0:
                return np.full(freqs.shape, float(arr))
            return arr.reshape(freqs.shape)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_f = freqs.size
            auto = np.empty((n_in, n_f), dtype=float)
            for i, exc in enumerate(excs):
                sc = exc.amplitude_scale()
                auto[i, :] = (sc * sc) * np.real(exc.psd.evaluate(freqs))
                if np.any(auto[i] < -tol * max(1.0, float(np.max(np.abs(auto[i]))))):
                    raise ValueError(f"Auto-PSD of excitation {i} has negative values.")
                auto[i] = np.maximum(auto[i], 0.0)

            s = np.zeros((n_in, n_in, n_f), dtype=complex)
            for i in range(n_in):
                s[i, i, :] = auto[i, :]

            if isinstance(coherence, Mapping):
                for (i, j), coh_spec in coherence.items():
                    if i == j:
                        continue
                    g2 = _eval_entry(coh_spec, freqs)
                    if np.any(g2 < -tol) or np.any(g2 > 1.0 + tol):
                        raise ValueError(f"Coherence gamma^2[({i},{j})] must lie in [0, 1].")
                    g2 = np.clip(g2, 0.0, 1.0)
                    ph_spec = 0.0
                    if phase is not None and isinstance(phase, Mapping):
                        if (i, j) in phase:
                            ph_spec = phase[(i, j)]
                        elif (j, i) in phase:
                            ph_spec = -np.asarray(_eval_entry(phase[(j, i)], freqs))
                    th = _eval_entry(ph_spec, freqs) * deg_factor
                    val = np.sqrt(g2 * auto[i] * auto[j]) * np.exp(1j * th)
                    s[i, j, :] = val
                    s[j, i, :] = val.conj()
            else:
                coh_arr = np.asarray(coherence, dtype=float)
                if coh_arr.shape == (n_in, n_in):
                    coh_grid = np.broadcast_to(coh_arr[:, :, None], (n_in, n_in, n_f))
                elif coh_arr.shape == (n_in, n_in, n_f):
                    coh_grid = coh_arr
                else:
                    raise ValueError(f"coherence array shape {coh_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")

                if phase is None:
                    ph_grid = np.zeros((n_in, n_in, n_f), dtype=float)
                else:
                    ph_arr = np.asarray(phase, dtype=float)
                    if ph_arr.shape == (n_in, n_in):
                        ph_grid = np.broadcast_to(ph_arr[:, :, None], (n_in, n_in, n_f))
                    elif ph_arr.shape == (n_in, n_in, n_f):
                        ph_grid = ph_arr
                    else:
                        raise ValueError(f"phase array shape {ph_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")
                ph_grid = ph_grid * deg_factor

                if np.any(coh_grid < -tol) or np.any(coh_grid > 1.0 + tol):
                    raise ValueError("Coherence gamma^2 entries must lie in [0, 1].")
                for i in range(n_in):
                    for j in range(n_in):
                        if i == j:
                            continue
                        g2 = np.clip(coh_grid[i, j], 0.0, 1.0)
                        s[i, j, :] = np.sqrt(g2 * auto[i] * auto[j]) * np.exp(1j * ph_grid[i, j])
            return s

        return cls(excs, _eval, n_in, tol=tol)

    @classmethod
    def from_co_quad(
        cls,
        excitations: Sequence[ExcitationBase],
        co_matrix: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable],
        quad_matrix: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable] | None = None,
        sign: Literal[-1, 1] = -1,
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Build a Cross-PSD matrix from co-spectrum ``C`` and quad-spectrum ``Q``.

        Uses ``S_ij(f) = C_ij(f) + sign * i * Q_ij(f)`` (default ``sign=-1``, i.e.
        ``C - i Q`` standard in Bendat & Piersol / random vibration).
        Note: Entries are multiplied by ``scale_i * scale_j`` from ``excitations``.
        """
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")
        if sign not in (-1, 1):
            raise ValueError("sign must be -1 (C - iQ) or +1 (C + iQ).")
        n_in = len(excs)
        scales = np.array([e.amplitude_scale() for e in excs], dtype=float)

        def _eval_entry(spec: float | PSDCurve | Callable | np.ndarray, freqs: np.ndarray) -> np.ndarray:
            if isinstance(spec, PSDCurve):
                return np.asarray(spec.evaluate(freqs), dtype=float)
            if callable(spec):
                v = np.asarray(spec(freqs), dtype=float)
                return np.full(freqs.shape, float(v)) if v.ndim == 0 else v.reshape(freqs.shape)
            a = np.asarray(spec, dtype=float)
            return np.full(freqs.shape, float(a)) if a.ndim == 0 else a.reshape(freqs.shape)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_f = freqs.size
            s = np.zeros((n_in, n_in, n_f), dtype=complex)
            for i, exc in enumerate(excs):
                s[i, i, :] = (scales[i] ** 2) * np.real(exc.psd.evaluate(freqs))

            if isinstance(co_matrix, Mapping):
                pairs = set(co_matrix.keys())
                if isinstance(quad_matrix, Mapping):
                    pairs |= set(quad_matrix.keys())
                for i, j in pairs:
                    c_val = _eval_entry(co_matrix.get((i, j), 0.0), freqs)
                    q_val = _eval_entry(quad_matrix.get((i, j), 0.0), freqs) if isinstance(quad_matrix, Mapping) else 0.0
                    val = (scales[i] * scales[j]) * (c_val + sign * 1j * q_val)
                    s[i, j, :] = val
                    if i != j:
                        s[j, i, :] = val.conj()
            else:
                co_arr = np.asarray(co_matrix, dtype=float)
                co_grid = (
                    np.broadcast_to(co_arr[:, :, None], (n_in, n_in, n_f))
                    if co_arr.shape == (n_in, n_in)
                    else co_arr
                )
                if co_grid.shape != (n_in, n_in, n_f):
                    raise ValueError(f"co_matrix shape {co_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")
                if quad_matrix is None:
                    q_grid = np.zeros((n_in, n_in, n_f), dtype=float)
                else:
                    q_arr = np.asarray(quad_matrix, dtype=float)  # type: ignore[arg-type]
                    q_grid = (
                        np.broadcast_to(q_arr[:, :, None], (n_in, n_in, n_f))
                        if q_arr.shape == (n_in, n_in)
                        else q_arr
                    )
                    if q_grid.shape != (n_in, n_in, n_f):
                        raise ValueError(f"quad_matrix shape {q_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")
                scale_mat = np.outer(scales, scales)[:, :, None]
                s = scale_mat * (co_grid + sign * 1j * q_grid)
            return s

        return cls(excs, _eval, n_in, tol=tol)

    @classmethod
    def from_magnitude_phase(
        cls,
        excitations: Sequence[ExcitationBase],
        magnitude: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable],
        phase: np.ndarray | Mapping[tuple[int, int], float | PSDCurve | Callable] = 0.0,
        phase_unit: Literal["rad", "deg"] = "rad",
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Build a Cross-PSD matrix from off-diagonal magnitude ``|S_ij|`` and phase."""
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")
        if phase_unit not in ("rad", "deg"):
            raise ValueError(f"phase_unit must be 'rad' or 'deg', got {phase_unit!r}.")
        deg_factor = math.pi / 180.0 if phase_unit == "deg" else 1.0
        n_in = len(excs)
        scales = np.array([e.amplitude_scale() for e in excs], dtype=float)

        def _eval_entry(spec: float | PSDCurve | Callable | np.ndarray, freqs: np.ndarray) -> np.ndarray:
            if isinstance(spec, PSDCurve):
                return np.asarray(spec.evaluate(freqs), dtype=float)
            if callable(spec):
                v = np.asarray(spec(freqs), dtype=float)
                return np.full(freqs.shape, float(v)) if v.ndim == 0 else v.reshape(freqs.shape)
            a = np.asarray(spec, dtype=float)
            return np.full(freqs.shape, float(a)) if a.ndim == 0 else a.reshape(freqs.shape)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_f = freqs.size
            s = np.zeros((n_in, n_in, n_f), dtype=complex)
            for i, exc in enumerate(excs):
                s[i, i, :] = (scales[i] ** 2) * np.real(exc.psd.evaluate(freqs))

            if isinstance(magnitude, Mapping):
                for (i, j), mag_spec in magnitude.items():
                    mag_v = _eval_entry(mag_spec, freqs)
                    if np.any(mag_v < 0.0):
                        raise ValueError(f"Magnitude entry ({i}, {j}) must be non-negative.")
                    ph_spec: Any = 0.0
                    if isinstance(phase, Mapping):
                        if (i, j) in phase:
                            ph_spec = phase[(i, j)]
                        elif (j, i) in phase:
                            ph_spec = -_eval_entry(phase[(j, i)], freqs)
                    th = _eval_entry(ph_spec, freqs) * deg_factor
                    val = (scales[i] * scales[j]) * mag_v * np.exp(1j * th)
                    s[i, j, :] = val
                    if i != j:
                        s[j, i, :] = val.conj()
            else:
                m_arr = np.asarray(magnitude, dtype=float)
                m_grid = np.broadcast_to(m_arr[:, :, None], (n_in, n_in, n_f)) if m_arr.shape == (n_in, n_in) else m_arr
                if m_grid.shape != (n_in, n_in, n_f):
                    raise ValueError(f"magnitude shape {m_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")
                if np.any(m_grid < 0.0):
                    raise ValueError("Cross-PSD magnitudes must be non-negative.")
                p_arr = np.asarray(phase, dtype=float)
                if p_arr.ndim == 0:
                    p_grid = np.full((n_in, n_in, n_f), float(p_arr))
                elif p_arr.shape == (n_in, n_in):
                    p_grid = np.broadcast_to(p_arr[:, :, None], (n_in, n_in, n_f))
                else:
                    p_grid = p_arr
                if p_grid.shape != (n_in, n_in, n_f):
                    raise ValueError(f"phase shape {p_arr.shape} incompatible with ({n_in}, {n_in}, {n_f}).")
                scale_mat = np.outer(scales, scales)[:, :, None]
                s = scale_mat * m_grid * np.exp(1j * (p_grid * deg_factor))
            return s

        return cls(excs, _eval, n_in, tol=tol)

    @classmethod
    def from_cross_terms(
        cls,
        excitations: Sequence[ExcitationBase],
        cross_terms: Mapping[tuple[int, int], PSDCurve | float | complex | Callable],
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Build a Cross-PSD matrix from explicit ``(i, j)`` complex cross-PSD curves."""
        excs = tuple(excitations)
        if not excs:
            raise ValueError("At least one excitation is required.")
        n_in = len(excs)
        scales = np.array([e.amplitude_scale() for e in excs], dtype=float)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            n_f = freqs.size
            s = np.zeros((n_in, n_in, n_f), dtype=complex)
            for i, exc in enumerate(excs):
                s[i, i, :] = (scales[i] ** 2) * np.real(exc.psd.evaluate(freqs))
            for (i, j), term in cross_terms.items():
                if not (0 <= i < n_in and 0 <= j < n_in):
                    raise IndexError(f"Cross-term index ({i}, {j}) out of range for {n_in} inputs.")
                curve = _as_psd_curve(term) if not isinstance(term, complex) else _CallablePSD(lambda f, v=term: v)
                val = (scales[i] * scales[j]) * np.asarray(curve.evaluate(freqs), dtype=complex)
                s[i, j, :] = val
                if i != j:
                    s[j, i, :] = val.conj()
            return s

        return cls(excs, _eval, n_in, tol=tol)

    @classmethod
    def from_array(
        cls,
        freqs_hz: Sequence[float] | np.ndarray,
        s_in: np.ndarray,
        excitations: Sequence[ExcitationBase] | None = None,
        tol: float = 1e-10,
    ) -> CrossPSDMatrix:
        """Wrap and validate a precomputed ``(n_in, n_in, n_f)`` Cross-PSD array."""
        f_ref = _as_1d_frequencies(freqs_hz, allow_zero=True)
        arr = np.asarray(s_in, dtype=complex)
        if arr.ndim != 3 or arr.shape[0] != arr.shape[1] or arr.shape[2] != f_ref.size:
            raise ValueError(f"s_in must have shape (n_in, n_in, {f_ref.size}), got {arr.shape}.")
        n_in = arr.shape[0]
        for k in range(f_ref.size):
            _validate_hermitian_psd_slice(arr[:, :, k], k, float(f_ref[k]), tol)

        def _eval(freqs: np.ndarray) -> np.ndarray:
            f_req = _as_1d_frequencies(freqs, allow_zero=True)
            if f_req.shape == f_ref.shape and np.allclose(f_req, f_ref, rtol=1e-12, atol=1e-15):
                return arr.copy()
            # Linearly interpolate each entry onto requested frequencies
            out = np.empty((n_in, n_in, f_req.size), dtype=complex)
            for i in range(n_in):
                for j in range(n_in):
                    out[i, j] = np.interp(f_req, f_ref, arr[i, j].real) + 1j * np.interp(f_req, f_ref, arr[i, j].imag)
            return out

        return cls(excitations, _eval, n_in, tol=tol)

    def evaluate(self, freqs_hz: Sequence[float] | np.ndarray, *, validate: bool = True) -> np.ndarray:
        """Evaluate and validate the ``(n_in, n_in, n_f)`` Hermitian Cross-PSD matrix."""
        f = _as_1d_frequencies(freqs_hz, allow_zero=True)
        s = np.asarray(self._evaluator(f), dtype=complex)
        expected = (self.n_inputs, self.n_inputs, f.size)
        if s.shape != expected:
            raise ValueError(f"CrossPSDMatrix evaluator returned shape {s.shape}, expected {expected}.")
        if validate:
            for k in range(f.size):
                _validate_hermitian_psd_slice(s[:, :, k], k, float(f[k]), self.tol)
        return s


# ---------------------------------------------------------------------------
# 4. Frequency Grid Builders
# ---------------------------------------------------------------------------


def uniform_frequency_grid(f_min: float, f_max: float, n_points: int = 500) -> np.ndarray:
    """Build a linearly spaced frequency grid on ``[f_min, f_max]`` in Hz."""
    if f_min < 0.0 or f_max <= f_min:
        raise ValueError(f"Require 0 <= f_min < f_max, got [{f_min}, {f_max}].")
    if n_points < 2:
        raise ValueError("n_points must be >= 2.")
    return np.linspace(float(f_min), float(f_max), int(n_points))


def log_frequency_grid(f_min: float, f_max: float, n_points: int = 500) -> np.ndarray:
    """Build a logarithmically spaced frequency grid on ``[f_min, f_max]`` in Hz."""
    if f_min <= 0.0 or f_max <= f_min:
        raise ValueError(f"Require 0 < f_min < f_max, got [{f_min}, {f_max}].")
    if n_points < 2:
        raise ValueError("n_points must be >= 2.")
    return np.geomspace(float(f_min), float(f_max), int(n_points))


def abaqus_frequency_grid(
    f_modes_hz: Sequence[float] | np.ndarray,
    f_min: float,
    f_max: float,
    n_points: int = 20,
    bias: float = 3.0,
    scale: Literal["LOGARITHMIC", "LINEAR"] = "LOGARITHMIC",
    model: str = "abaqus",
) -> np.ndarray:
    """Build an Abaqus-style biased frequency sweep split at eigenfrequencies."""
    return abaqus_frequency_sweep(
        np.asarray(f_modes_hz, dtype=float),
        float(f_min),
        float(f_max),
        n_points=int(n_points),
        bias=float(bias),
        scale=scale,
        model=model,
    )


def adaptive_modal_frequency_grid(
    f_modes_hz: Sequence[float] | np.ndarray,
    damping: DampingModel | float | Sequence[float] | np.ndarray,
    f_min: float,
    f_max: float,
    n_base: int = 81,
    scale: Literal["log", "linear"] = "log",
    half_power_multiplier: float = 35.0,
    points_per_mode: int = 121,
    min_zeta: float = 1e-4,
) -> np.ndarray:
    """Build a frequency grid with hyperbolic-sine half-power clustering at resonances.

    Around each mode ``f_n`` with equivalent viscous damping ``zeta_n``, maps
    ``u = sinh(s)`` uniformly in ``s in [-asinh(c), +asinh(c)]`` so
    ``f = f_n * (1 + zeta_eff * u)`` covers ``f_n * (1 +/- c * zeta_eff)``.
    The spacing is dense and uniform across the half-power core ``f_n (1 +/- zeta_n)``
    and transitions logarithmically into the wings, resolving sharp resonance
    peaks (``zeta ~ 0.005``) accurately with far fewer points than uniform grids.

    Args:
        f_modes_hz: Modal natural frequencies in Hz.
        damping: Damping model, scalar ``zeta``, or per-mode ``zeta`` sequence.
        f_min: Lower frequency bound in Hz (``>= 0``; ``> 0`` for ``scale="log"``).
        f_max: Upper frequency bound in Hz (``> f_min``).
        n_base: Number of background grid points spanning ``[f_min, f_max]``.
        scale: Background grid spacing: ``"log"`` (default) or ``"linear"``.
        half_power_multiplier: Half-width factor ``c`` in ``f_n (1 +/- c * zeta_n)``.
        points_per_mode: Number of clustered points per mode (odd preferred so ``f_n``
            is hit exactly).
        min_zeta: Floor on ``zeta_n`` for clustering width when damping is zero.

    Returns:
        Strictly increasing 1D frequency grid in Hz.
    """
    if f_min < 0.0 or f_max <= f_min:
        raise ValueError(f"Require 0 <= f_min < f_max, got [{f_min}, {f_max}].")
    if n_base < 2 or points_per_mode < 3:
        raise ValueError("Require n_base >= 2 and points_per_mode >= 3.")
    if half_power_multiplier <= 0.0:
        raise ValueError("half_power_multiplier must be positive.")

    f_modes = np.asarray(f_modes_hz, dtype=float).reshape(-1)
    omega_modes = 2.0 * np.pi * f_modes
    d_model = _as_damping_model(damping)
    zetas = d_model.equivalent_viscous_ratios(omega_modes) if f_modes.size > 0 else np.empty(0)

    scale_l = scale.strip().lower()
    if scale_l == "log":
        base = log_frequency_grid(max(float(f_min), 1e-12) if f_min == 0.0 else f_min, f_max, n_base)
        if f_min == 0.0:
            base = np.concatenate([[0.0], base])
    elif scale_l == "linear":
        base = uniform_frequency_grid(f_min, f_max, n_base)
    else:
        raise ValueError(f"scale must be 'log' or 'linear', got {scale!r}.")

    clusters = [base]
    s_max = math.asinh(float(half_power_multiplier))
    s_vals = np.linspace(-s_max, s_max, int(points_per_mode))
    u_sinh = np.sinh(s_vals)

    for fn, zn in zip(f_modes, zetas, strict=True):
        if fn <= 0.0:
            continue
        z_eff = max(float(zn), float(min_zeta))
        f_lo = fn * (1.0 - half_power_multiplier * z_eff)
        f_hi = fn * (1.0 + half_power_multiplier * z_eff)
        if f_hi < f_min or f_lo > f_max:
            continue
        f_c = fn * (1.0 + z_eff * u_sinh)
        f_c = f_c[(f_c >= f_min) & (f_c <= f_max)]
        if f_c.size > 0:
            clusters.append(f_c)
        if f_min <= fn <= f_max:
            clusters.append(np.array([fn]))

    merged = np.sort(np.concatenate(clusters))
    # Deduplicate near-coincident points to avoid degenerate trapezoid intervals
    keep = np.ones(merged.size, dtype=bool)
    tol_step = 1e-11 * max(float(f_max), 1.0)
    last_kept = merged[0]
    for i in range(1, merged.size):
        if merged[i] - last_kept <= max(1e-12 * max(last_kept, 1.0), tol_step * 1e-3):
            keep[i] = False
        else:
            last_kept = merged[i]
    out = merged[keep]
    out[0] = float(f_min)
    if out[-1] < float(f_max):
        out = np.append(out, float(f_max))
    else:
        out[-1] = float(f_max)
    return out


# ---------------------------------------------------------------------------
# 5. Post-processing & Fatigue
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpectralMoments:
    """Spectral moments ``lambda_0, lambda_1, lambda_2, lambda_4`` and statistics.

    Moments are defined over frequency ``f`` in Hz:
    ``lambda_n = integral(f**n * S(f) df)``.
    """

    m0: np.ndarray | float
    m1: np.ndarray | float
    m2: np.ndarray | float
    m4: np.ndarray | float
    rms: np.ndarray | float
    peak_2sigma: np.ndarray | float
    peak_3sigma: np.ndarray | float
    f0: np.ndarray | float
    fp: np.ndarray | float
    gamma2: np.ndarray | float

    @property
    def lambda_0(self) -> np.ndarray | float:
        return self.m0

    @property
    def lambda_1(self) -> np.ndarray | float:
        return self.m1

    @property
    def lambda_2(self) -> np.ndarray | float:
        return self.m2

    @property
    def lambda_4(self) -> np.ndarray | float:
        return self.m4

    @property
    def zero_crossing_rate(self) -> np.ndarray | float:
        return self.f0

    @property
    def peak_rate(self) -> np.ndarray | float:
        return self.fp

    @property
    def irregularity_factor(self) -> np.ndarray | float:
        return self.gamma2


def spectral_moments(
    freqs_hz: Sequence[float] | np.ndarray,
    psd: np.ndarray,
    axis: int = -1,
) -> SpectralMoments:
    """Compute spectral moments ``lambda_0, lambda_1, lambda_2, lambda_4`` and metrics.

    Args:
        freqs_hz: Frequency grid in Hz, shape ``(n_f,)``.
        psd: Non-negative real one-sided PSD array with length ``n_f`` along ``axis``.
        axis: Frequency axis in ``psd`` (default ``-1``).

    Returns:
        :class:`SpectralMoments` with scalars (for 1D ``psd``) or arrays.
    """
    f = _as_1d_frequencies(freqs_hz, allow_zero=True)
    if f.size < 2:
        raise ValueError("spectral_moments requires at least 2 frequency points.")
    s_raw = np.asarray(psd)
    if np.iscomplexobj(s_raw):
        max_im = float(np.max(np.abs(s_raw.imag)))
        max_re = max(1.0, float(np.max(np.abs(s_raw.real))))
        if max_im > 1e-10 * max_re:
            raise ValueError(f"Auto-PSD must be real; found max|imag| = {max_im:.3e}.")
        s_real = s_raw.real.astype(float)
    else:
        s_real = s_raw.astype(float)

    s_moved = np.moveaxis(s_real, axis, -1)
    if s_moved.shape[-1] != f.size:
        raise ValueError(f"psd length along axis {axis} ({s_moved.shape[-1]}) != len(freqs_hz) ({f.size}).")
    s_clamped = np.maximum(s_moved, 0.0)

    m0 = integrate.trapezoid(s_clamped, f, axis=-1)
    m1 = integrate.trapezoid(s_clamped * f, f, axis=-1)
    f2 = f * f
    m2 = integrate.trapezoid(s_clamped * f2, f, axis=-1)
    m4 = integrate.trapezoid(s_clamped * (f2 * f2), f, axis=-1)

    m0 = np.maximum(m0, 0.0)
    m1 = np.maximum(m1, 0.0)
    m2 = np.maximum(m2, 0.0)
    m4 = np.maximum(m4, 0.0)

    rms = np.sqrt(m0)
    with np.errstate(divide="ignore", invalid="ignore"):
        f0 = np.where(m0 > 0.0, np.sqrt(m2 / m0), 0.0)
        fp = np.where(m2 > 0.0, np.sqrt(m4 / m2), 0.0)
        denom_g = np.sqrt(m0 * m4)
        gamma2 = np.where(denom_g > 0.0, np.clip(m2 / denom_g, 0.0, 1.0), 0.0)

    is_scalar = s_real.ndim == 1
    if is_scalar:
        return SpectralMoments(
            m0=float(m0),
            m1=float(m1),
            m2=float(m2),
            m4=float(m4),
            rms=float(rms),
            peak_2sigma=float(2.0 * rms),
            peak_3sigma=float(3.0 * rms),
            f0=float(f0),
            fp=float(fp),
            gamma2=float(gamma2),
        )
    return SpectralMoments(
        m0=m0,
        m1=m1,
        m2=m2,
        m4=m4,
        rms=rms,
        peak_2sigma=2.0 * rms,
        peak_3sigma=3.0 * rms,
        f0=f0,
        fp=fp,
        gamma2=gamma2,
    )


def von_mises_psd(stress_cpsd: np.ndarray) -> np.ndarray:
    """Compute Segalman equivalent von Mises stress PSD ``S_vM(f) = Re tr(A_vM S_ss(f))``.

    Args:
        stress_cpsd: Complex Hermitian Voigt stress CPSD of shape ``(6, 6, n_f)``
            or ``(..., 6, 6, n_f)``.

    Returns:
        Non-negative real von Mises stress PSD of shape ``(n_f,)`` or ``(..., n_f)``.
    """
    s = np.asarray(stress_cpsd)
    if s.ndim < 3 or s.shape[-3:-1] != (6, 6):
        raise ValueError(f"stress_cpsd must have shape (..., 6, 6, n_f), got {s.shape}.")
    eq = np.einsum("ij,...jif->...f", VON_MISES_A, s)
    return np.maximum(np.real(eq), 0.0)


def segalman_von_mises_rms(
    stress_cpsd_or_cov: np.ndarray,
    freqs_hz: Sequence[float] | np.ndarray | None = None,
) -> np.ndarray | float:
    """Compute Segalman-Fulcher RMS von Mises stress ``sqrt(tr(A_vM * Lambda_0^sigma))``.

    Accepts either:
        * Stress covariance matrix ``Lambda_0^sigma`` of shape ``(6, 6)`` or
          ``(..., 6, 6)`` (when ``freqs_hz is None``), or
        * Frequency-dependent complex Hermitian stress CPSD ``S_ss(f)`` of shape
          ``(6, 6, n_f)`` or ``(..., 6, 6, n_f)`` (when ``freqs_hz`` is supplied).
    """
    arr = np.asarray(stress_cpsd_or_cov)
    if freqs_hz is not None:
        s_vm = von_mises_psd(arr)
        f = _as_1d_frequencies(freqs_hz, allow_zero=True)
        if s_vm.shape[-1] != f.size:
            raise ValueError(f"stress_cpsd frequency length ({s_vm.shape[-1]}) != len(freqs_hz) ({f.size}).")
        var = integrate.trapezoid(s_vm, f, axis=-1)
        rms = np.sqrt(np.maximum(var, 0.0))
        return float(rms) if arr.ndim == 3 else rms

    if arr.ndim < 2 or arr.shape[-2:] != (6, 6):
        raise ValueError(f"stress covariance must have shape (..., 6, 6), got {arr.shape}.")
    cov_real = np.real(arr)
    var = np.einsum("ij,...ji->...", VON_MISES_A, cov_real)
    rms = np.sqrt(np.maximum(var, 0.0))
    return float(rms) if arr.ndim == 2 else rms


@dataclass(frozen=True)
class SNCurve:
    """Basquin S-N fatigue curve ``N * S**b = C``.

    Args:
        C: Fatigue strength constant ``C > 0``.
        b: Basquin slope exponent ``b > 0``.
        stress_definition: Whether ``S`` in ``N * S**b = C`` denotes zero-to-peak
            stress ``"amplitude"`` (default) or peak-to-peak stress ``"range"``
            (``2 * amplitude``).
    """

    C: float
    b: float
    stress_definition: StressDefinition = "amplitude"

    def __post_init__(self) -> None:
        if not math.isfinite(self.C) or self.C <= 0.0:
            raise ValueError(f"S-N constant C must be positive and finite, got {self.C}.")
        if not math.isfinite(self.b) or self.b <= 0.0:
            raise ValueError(f"S-N exponent b must be positive and finite, got {self.b}.")
        sdef = self.stress_definition.strip().lower()
        if sdef not in ("amplitude", "range"):
            raise ValueError(f"stress_definition must be 'amplitude' or 'range', got {self.stress_definition!r}.")
        object.__setattr__(self, "stress_definition", sdef)

    @property
    def amplitude_scale(self) -> float:
        """Multiplier converting stress amplitude ``S_a`` into the curve's ``S``."""
        return 2.0 if self.stress_definition == "range" else 1.0


@dataclass(frozen=True)
class FatigueResult:
    """Random vibration fatigue damage and life summary."""

    method: str
    damage_rate: np.ndarray | float
    life_seconds: np.ndarray | float
    life_cycles: np.ndarray | float
    cycle_rate_hz: np.ndarray | float
    damage_in_duration: np.ndarray | float | None = None


class FatigueModel(ABC):
    """Abstract base class for frequency-domain random vibration fatigue models."""

    @abstractmethod
    def evaluate(
        self,
        moments: SpectralMoments,
        sn_curve: SNCurve,
        duration_seconds: float | None = None,
    ) -> FatigueResult:
        """Compute expected fatigue damage rate and life from spectral moments."""

    def damage_rate(self, moments: SpectralMoments, sn_curve: SNCurve) -> np.ndarray | float:
        """Return expected Miner's damage per second (``1 / s``)."""
        return self.evaluate(moments, sn_curve).damage_rate


class NarrowBandRayleigh(FatigueModel):
    """Narrow-band Rayleigh peak fatigue model (Bendat / Miles).

    Assumes Rayleigh distribution of stress amplitudes ``S_a`` with rate ``f_0``:
    ``D_rate = (f_0 / C) * (sqrt(2 * lambda_0) * scale)**b * Gamma(1 + b / 2)``.
    """

    @staticmethod
    def pdf(stress_amplitudes: Sequence[float] | np.ndarray, moments: SpectralMoments) -> np.ndarray:
        """Evaluate Rayleigh amplitude PDF ``p(S_a) = (S_a / m0) * exp(-S_a^2 / (2 m0))``."""
        s = np.asarray(stress_amplitudes, dtype=float)
        m0 = float(moments.m0)
        if m0 <= 0.0:
            return np.zeros_like(s)
        return np.where(s >= 0.0, (s / m0) * np.exp(-0.5 * (s * s) / m0), 0.0)

    def evaluate(
        self,
        moments: SpectralMoments,
        sn_curve: SNCurve,
        duration_seconds: float | None = None,
    ) -> FatigueResult:
        rms = np.asarray(moments.rms, dtype=float)
        f0 = np.asarray(moments.f0, dtype=float)
        b = sn_curve.b
        scale = sn_curve.amplitude_scale
        expected_sb = ((math.sqrt(2.0) * scale * rms) ** b) * float(gamma_func(1.0 + 0.5 * b))
        d_rate = (f0 / sn_curve.C) * expected_sb
        return _package_fatigue_result("NarrowBandRayleigh", d_rate, f0, duration_seconds, np.ndim(moments.m0) == 0)


class Steinberg3Band(FatigueModel):
    """Steinberg three-band Gaussian fatigue approximation.

    Uses discrete 1-sigma (68.31%), 2-sigma (27.10%), and 3-sigma (4.33%)
    amplitude bands at cycle rate ``f_0`` (or ``f_p`` when configured).
    """

    def __init__(self, cycle_rate: Literal["f0", "fp"] = "f0") -> None:
        cr = cycle_rate.strip().lower()
        if cr not in ("f0", "fp"):
            raise ValueError(f"cycle_rate must be 'f0' or 'fp', got {cycle_rate!r}.")
        self.cycle_rate = cr

    def evaluate(
        self,
        moments: SpectralMoments,
        sn_curve: SNCurve,
        duration_seconds: float | None = None,
    ) -> FatigueResult:
        rms = np.asarray(moments.rms, dtype=float)
        rate = np.asarray(moments.f0 if self.cycle_rate == "f0" else moments.fp, dtype=float)
        b = sn_curve.b
        scale = sn_curve.amplitude_scale
        band_factor = sum(w * ((k * scale) ** b) for k, w in _STEINBERG_WEIGHTS)
        d_rate = (rate / sn_curve.C) * (rms**b) * band_factor
        return _package_fatigue_result("Steinberg3Band", d_rate, rate, duration_seconds, np.ndim(moments.m0) == 0)


class Dirlik(FatigueModel):
    """Dirlik (1985) empirical rainflow amplitude PDF and closed-form damage model.

    Combines an exponential term and two Rayleigh terms parameterized by
    ``m0, m1, m2, m4``, accurately capturing wide-band and multi-modal spectra
    while reducing smoothly to :class:`NarrowBandRayleigh` in the narrow-band
    limit ``gamma_2 -> 1``.
    """

    @staticmethod
    def parameters(moments: SpectralMoments) -> tuple[np.ndarray, ...]:
        """Compute Dirlik parameters ``(D1, D2, D3, Q, R)`` from spectral moments."""
        m0 = np.asarray(moments.m0, dtype=float)
        m1 = np.asarray(moments.m1, dtype=float)
        m2 = np.asarray(moments.m2, dtype=float)
        m4 = np.asarray(moments.m4, dtype=float)

        with np.errstate(divide="ignore", invalid="ignore"):
            gamma2 = np.where((m0 > 0.0) & (m4 > 0.0), np.clip(m2 / np.sqrt(m0 * m4), 0.0, 1.0), 0.0)
            xm = np.where((m0 > 0.0) & (m4 > 0.0), (m1 / m0) * np.sqrt(m2 / m4), 0.0)

            d1 = 2.0 * (xm - gamma2 * gamma2) / np.maximum(1.0 + gamma2 * gamma2, 1e-30)
            d1 = np.clip(d1, 0.0, 1.0)
            den_r = 1.0 - gamma2 - d1 + d1 * d1
            narrow = (gamma2 >= 1.0 - 1e-10) | (np.abs(den_r) <= 1e-12) | (d1 <= 1e-12)

            r_val = np.where(narrow, 1.0, (gamma2 - xm - d1 * d1) / np.where(narrow, 1.0, den_r))
            r_val = np.clip(r_val, 1e-12, 1.0 - 1e-12)
            r_val = np.where(narrow, 1.0, r_val)

            d2 = np.where(narrow, 0.0, den_r / np.maximum(1.0 - r_val, 1e-30))
            d2 = np.clip(d2, 0.0, 1.0 - d1)
            d3 = np.maximum(1.0 - d1 - d2, 0.0)
            # Algebraically gamma2 - D3 - D2*R == D1**2, avoiding 0/0 cancellation near narrowband
            q_val = np.where(narrow, 1.0, 1.25 * d1)
            d1 = np.where(narrow, 0.0, d1)
            d3 = np.where(narrow, 1.0, d3)

        return d1, d2, d3, q_val, r_val

    @classmethod
    def pdf(cls, stress_amplitudes: Sequence[float] | np.ndarray, moments: SpectralMoments) -> np.ndarray:
        """Evaluate the Dirlik rainflow stress-amplitude PDF ``p(S_a)``."""
        s = np.asarray(stress_amplitudes, dtype=float)
        rms = float(moments.rms)
        if rms <= 0.0:
            return np.zeros_like(s)
        d1, d2, d3, q_val, r_val = cls.parameters(moments)
        d1_f, d2_f, d3_f, q_f, r_f = float(d1), float(d2), float(d3), float(q_val), float(r_val)
        z = np.maximum(s, 0.0) / rms
        pz = np.zeros_like(z)
        if d1_f > 0.0 and q_f > 0.0:
            pz += (d1_f / q_f) * np.exp(-z / q_f)
        if d2_f > 0.0:
            pz += (d2_f * z / (r_f * r_f)) * np.exp(-0.5 * (z * z) / (r_f * r_f))
        if d3_f > 0.0:
            pz += d3_f * z * np.exp(-0.5 * z * z)
        return np.where(s >= 0.0, pz / rms, 0.0)

    def evaluate(
        self,
        moments: SpectralMoments,
        sn_curve: SNCurve,
        duration_seconds: float | None = None,
    ) -> FatigueResult:
        rms = np.asarray(moments.rms, dtype=float)
        fp = np.asarray(moments.fp, dtype=float)
        b = sn_curve.b
        scale = sn_curve.amplitude_scale

        d1, d2, d3, q_val, r_val = self.parameters(moments)
        g_exp = float(gamma_func(1.0 + b))
        g_ray = float(gamma_func(1.0 + 0.5 * b))
        two_half_b = 2.0 ** (0.5 * b)

        term1 = np.where(d1 > 0.0, d1 * (q_val**b) * g_exp, 0.0)
        term2 = two_half_b * g_ray * (d2 * (r_val**b) + d3)
        expected_zb = term1 + term2
        d_rate = (fp / sn_curve.C) * ((scale * rms) ** b) * expected_zb
        return _package_fatigue_result("Dirlik", d_rate, fp, duration_seconds, np.ndim(moments.m0) == 0)


def _package_fatigue_result(
    method: str,
    d_rate: np.ndarray,
    cycle_rate: np.ndarray,
    duration_seconds: float | None,
    is_scalar: bool,
) -> FatigueResult:
    """Wrap fatigue arrays into :class:`FatigueResult` with safe zero-damage handling."""
    d_arr = np.maximum(np.asarray(d_rate, dtype=float), 0.0)
    r_arr = np.maximum(np.asarray(cycle_rate, dtype=float), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        life_s = np.where(d_arr > 0.0, 1.0 / d_arr, np.inf)
        life_cyc = np.where(d_arr > 0.0, r_arr * life_s, np.inf)
    dmg_dur: np.ndarray | float | None = None
    if duration_seconds is not None:
        dur = float(duration_seconds)
        if dur < 0.0 or not math.isfinite(dur):
            raise ValueError(f"duration_seconds must be finite and >= 0, got {duration_seconds}.")
        dmg_dur = d_arr * dur
        if is_scalar:
            dmg_dur = float(dmg_dur)
    if is_scalar:
        return FatigueResult(
            method=method,
            damage_rate=float(d_arr),
            life_seconds=float(life_s),
            life_cycles=float(life_cyc),
            cycle_rate_hz=float(r_arr),
            damage_in_duration=dmg_dur,
        )
    return FatigueResult(
        method=method,
        damage_rate=d_arr,
        life_seconds=life_s,
        life_cycles=life_cyc,
        cycle_rate_hz=r_arr,
        damage_in_duration=dmg_dur,
    )


# ---------------------------------------------------------------------------
# 6. Solvers & Result Container
# ---------------------------------------------------------------------------


def compute_residual_flexibility(
    K_full: sp.spmatrix | np.ndarray,
    phi: np.ndarray,
    omega_modes: np.ndarray,
    M_modal: np.ndarray | None = None,
    free_dofs: Sequence[int] | np.ndarray | None = None,
) -> np.ndarray:
    """Compute static residual flexibility matrix ``G_res = K^{-1} - Phi diag(1 / (m_k w_k^2)) Phi^T``.

    Constrained DOFs (outside ``free_dofs``) have zero rows and columns in ``G_res``.
    """
    phi_arr = np.asarray(phi, dtype=float)
    w = np.asarray(omega_modes, dtype=float).reshape(-1)
    n_dofs, n_modes = phi_arr.shape
    if w.size != n_modes or np.any(w <= 0.0):
        raise ValueError("omega_modes must match phi.shape[1] and be strictly positive.")
    m = np.ones(n_modes, dtype=float) if M_modal is None else np.asarray(M_modal, dtype=float).reshape(n_modes)
    if np.any(m <= 0.0):
        raise ValueError("M_modal entries must be strictly positive.")

    free = np.arange(n_dofs, dtype=np.int64) if free_dofs is None else np.unique(np.asarray(free_dofs, dtype=np.int64))
    K_csr = sp.csr_matrix(K_full, dtype=float)
    if K_csr.shape != (n_dofs, n_dofs):
        raise ValueError(f"K_full shape {K_csr.shape} does not match n_dofs ({n_dofs}).")

    K_ff = K_csr[free, :][:, free]
    K_ff = 0.5 * (K_ff + K_ff.T)
    eye_ff = np.eye(free.size, dtype=float)
    if free.size <= 512:
        k_inv_ff = la.solve(K_ff.toarray(), eye_ff, assume_a="sym")
    else:
        lu = spla.splu(K_ff.tocsc())
        k_inv_ff = lu.solve(eye_ff)

    G_res = np.zeros((n_dofs, n_dofs), dtype=float)
    G_res[np.ix_(free, free)] = k_inv_ff
    modal_flex = (phi_arr / (m * (w * w))[None, :]) @ phi_arr.T
    G_res -= modal_flex
    return 0.5 * (G_res + G_res.T)


@dataclass(frozen=True)
class RandomVibrationResult:
    """Complete output of a PSD random vibration solve.

    Attributes:
        frequencies_hz: Frequency grid used in Hz, shape ``(n_f,)``.
        disp_psd: Output auto-PSD of primary response (displacement), shape ``(n_out, n_f)``.
        vel_psd: Output velocity auto-PSD ``(2 pi f)^2 S_u``, shape ``(n_out, n_f)``.
        acc_psd: Output acceleration auto-PSD ``(2 pi f)^4 S_u`` (or absolute
            acceleration PSD when ``response_frame="absolute"``), shape ``(n_out, n_f)``.
        modal_psd: Complex Hermitian modal generalized coordinate PSD ``S_qq(f)``
            of shape ``(n_modes, n_modes, n_f)``, or ``None`` for direct FRF solves.
        stress_psd: Optional complex Hermitian Voigt stress CPSD ``S_ss(f)``
            of shape ``(n_pts, 6, 6, n_f)``.
        rmises: Optional Segalman-Fulcher RMS von Mises stress, shape ``(n_pts,)``.
        von_mises_psd: Optional equivalent von Mises stress PSD, shape ``(n_pts, n_f)``.
        input_cross_psd: Evaluated input Cross-PSD ``S_in(f)``, shape ``(n_in, n_in, n_f)``.
        response_frame: ``"relative"`` or ``"absolute"``.
    """

    frequencies_hz: np.ndarray
    disp_psd: np.ndarray
    vel_psd: np.ndarray
    acc_psd: np.ndarray
    modal_psd: np.ndarray | None
    stress_psd: np.ndarray | None
    rmises: np.ndarray | None
    von_mises_psd: np.ndarray | None
    input_cross_psd: np.ndarray
    response_frame: ResponseFrame = "relative"
    _cross_evaluator: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None

    @property
    def ru(self) -> np.ndarray:
        """RMS displacement per output DOF, shape ``(n_out,)``."""
        return self.displacement_moments().rms  # type: ignore[return-value]

    @property
    def rv(self) -> np.ndarray:
        """RMS velocity per output DOF, shape ``(n_out,)``."""
        return self.velocity_moments().rms  # type: ignore[return-value]

    @property
    def ra(self) -> np.ndarray:
        """RMS acceleration per output DOF, shape ``(n_out,)``."""
        return self.acceleration_moments().rms  # type: ignore[return-value]

    def output_psd(self, quantity: Literal["displacement", "velocity", "acceleration"] = "displacement") -> np.ndarray:
        """Return auto-PSD array of shape ``(n_out, n_f)`` for the requested quantity."""
        q = quantity.strip().lower()
        if q in ("displacement", "u", "d"):
            return self.disp_psd
        if q in ("velocity", "v"):
            return self.vel_psd
        if q in ("acceleration", "a"):
            return self.acc_psd
        raise ValueError(f"Unknown quantity {quantity!r}; expected 'displacement', 'velocity', or 'acceleration'.")

    def cross_psd(
        self,
        dofs_i: Sequence[int] | np.ndarray,
        dofs_j: Sequence[int] | np.ndarray | None = None,
    ) -> np.ndarray:
        """Compute displacement cross-PSD ``S_uu[dofs_i, dofs_j, f]`` on demand.

        Args:
            dofs_i: Row output indices of length ``ni``.
            dofs_j: Column output indices of length ``nj`` (defaults to ``dofs_i``).

        Returns:
            Complex Hermitian Cross-PSD array of shape ``(ni, nj, n_f)``.
        """
        idx_i = np.asarray(dofs_i, dtype=np.int64).reshape(-1)
        idx_j = idx_i if dofs_j is None else np.asarray(dofs_j, dtype=np.int64).reshape(-1)
        if self._cross_evaluator is None:
            raise RuntimeError("Cross-PSD evaluator is not available for this result.")
        return self._cross_evaluator(idx_i, idx_j)

    def displacement_moments(self, dof: int | None = None) -> SpectralMoments:
        """Compute spectral moments of displacement PSD (all DOFs or single ``dof``)."""
        data = self.disp_psd[dof] if dof is not None else self.disp_psd
        return spectral_moments(self.frequencies_hz, data, axis=-1)

    def velocity_moments(self, dof: int | None = None) -> SpectralMoments:
        """Compute spectral moments of velocity PSD (all DOFs or single ``dof``)."""
        data = self.vel_psd[dof] if dof is not None else self.vel_psd
        return spectral_moments(self.frequencies_hz, data, axis=-1)

    def acceleration_moments(self, dof: int | None = None) -> SpectralMoments:
        """Compute spectral moments of acceleration PSD (all DOFs or single ``dof``)."""
        data = self.acc_psd[dof] if dof is not None else self.acc_psd
        return spectral_moments(self.frequencies_hz, data, axis=-1)

    def stress_moments(self, point_index: int | None = None, component: int | None = None) -> SpectralMoments:
        """Compute spectral moments of stress PSD.

        If ``component is None``, uses Segalman equivalent von Mises stress PSD
        ``S_vM(f)``. If ``component in 0..5`` is specified, uses the Voigt
        component auto-PSD ``S_ss[point, component, component, f]``.
        """
        if self.stress_psd is None or self.von_mises_psd is None:
            raise ValueError("No stress_modes were provided during the random vibration solve.")
        if component is None:
            src = self.von_mises_psd[point_index] if point_index is not None else self.von_mises_psd
        else:
            if not 0 <= int(component) < 6:
                raise ValueError(f"Voigt stress component must be in 0..5, got {component}.")
            c = int(component)
            src = (
                np.real(self.stress_psd[point_index, c, c, :])
                if point_index is not None
                else np.real(self.stress_psd[:, c, c, :])
            )
        return spectral_moments(self.frequencies_hz, src, axis=-1)

    def fatigue(
        self,
        sn_curve: SNCurve,
        method: FatigueModel | Literal["dirlik", "rayleigh", "steinberg"] = "dirlik",
        point_index: int | None = None,
        component: int | None = None,
        duration_seconds: float | None = None,
    ) -> FatigueResult:
        """Evaluate random vibration fatigue damage and life on stress results."""
        moments = self.stress_moments(point_index=point_index, component=component)
        if isinstance(method, FatigueModel):
            model = method
        else:
            m = method.strip().lower()
            if m == "dirlik":
                model = Dirlik()
            elif m in ("rayleigh", "narrowband", "narrow_band"):
                model = NarrowBandRayleigh()
            elif m in ("steinberg", "steinberg3band"):
                model = Steinberg3Band()
            else:
                raise ValueError(f"Unknown fatigue method {method!r}.")
        return model.evaluate(moments, sn_curve, duration_seconds=duration_seconds)


class ModalRandomVibration:
    """Modal superposition PSD random vibration solver.

    Accepts modal frequencies, mode shapes, optional modal or full mass matrix,
    damping model, optional Voigt stress mode shapes, optional modal
    participation factors, and optional static residual flexibility.

    Args:
        frequencies_hz: Modal natural frequencies in Hz, shape ``(n_modes,)``.
        phi: Mode shapes ``Phi``, shape ``(n_dofs, n_modes)``.
        M_modal: Optional generalized modal masses ``m_k = phi_k^T M phi_k``
            (1D ``(n_modes,)`` or 2D ``(n_modes, n_modes)``). If ``None``,
            computed from ``M_full`` if given, otherwise assumed ``1.0``.
        M_full: Optional full mass matrix ``(n_dofs, n_dofs)`` (dense or sparse).
        damping: :class:`DampingModel`, scalar ``zeta``, or per-mode ``zetas``.
        stress_modes: Optional Voigt stress mode shapes ``(n_pts, n_modes, 6)``.
        participation_factors: Optional modal participation factors ``(n_modes, 3)``
            or ``(n_modes,)`` for :class:`BaseExcitation`.
        residual_flexibility: Optional static residual flexibility matrix
            ``(n_dofs, n_dofs)`` or callable ``(F_mat) -> u_res_mat``.
        influence_vectors: Optional rigid-body influence vectors ``(n_dofs, 3)``
            used for absolute base-excitation acceleration and residual inertial loads.
    """

    def __init__(
        self,
        frequencies_hz: Sequence[float] | np.ndarray,
        phi: np.ndarray,
        M_modal: np.ndarray | None = None,
        M_full: sp.spmatrix | np.ndarray | None = None,
        damping: DampingModel | float | Sequence[float] | np.ndarray = 0.02,
        stress_modes: np.ndarray | None = None,
        participation_factors: np.ndarray | None = None,
        residual_flexibility: np.ndarray | sp.spmatrix | Callable[[np.ndarray], np.ndarray] | None = None,
        influence_vectors: np.ndarray | None = None,
    ) -> None:
        f_modes = _as_1d_frequencies(frequencies_hz, allow_zero=False)
        phi_arr = np.asarray(phi, dtype=float)
        if phi_arr.ndim != 2 or phi_arr.shape[1] != f_modes.size:
            raise ValueError(f"phi must have shape (n_dofs, {f_modes.size}), got {phi_arr.shape}.")
        n_dofs, n_modes = phi_arr.shape

        self.frequencies_hz = f_modes
        self.omega_modes = 2.0 * np.pi * f_modes
        self.phi = phi_arr
        self.n_dofs = n_dofs
        self.n_modes = n_modes
        self.damping: DampingModel = _as_damping_model(damping)
        self.M_full = M_full

        if M_modal is not None:
            m_arr = np.asarray(M_modal, dtype=float)
            m_vec = np.diag(m_arr) if m_arr.ndim == 2 else m_arr.reshape(-1)
            if m_vec.size != n_modes or np.any(m_vec <= 0.0):
                raise ValueError("M_modal must have length n_modes with strictly positive entries.")
            self.M_modal = m_vec
        elif M_full is not None:
            m_phi = M_full @ phi_arr
            m_vec = np.sum(phi_arr * m_phi, axis=0)
            if np.any(m_vec <= 0.0):
                raise ValueError("Computed modal masses Phi^T M Phi must be strictly positive.")
            self.M_modal = m_vec
        else:
            self.M_modal = np.ones(n_modes, dtype=float)

        if stress_modes is not None:
            psi = np.asarray(stress_modes, dtype=float)
            if psi.ndim != 3 or psi.shape[1] != n_modes or psi.shape[2] != 6:
                raise ValueError(f"stress_modes must have shape (n_pts, {n_modes}, 6), got {psi.shape}.")
            self.stress_modes: np.ndarray | None = psi
        else:
            self.stress_modes = None

        self.influence_vectors = None if influence_vectors is None else np.asarray(influence_vectors, dtype=float)
        if self.influence_vectors is not None and self.influence_vectors.shape[0] != n_dofs:
            raise ValueError("influence_vectors row count must equal n_dofs.")

        if participation_factors is not None:
            gam = np.asarray(participation_factors, dtype=float)
            if gam.ndim == 1:
                if gam.size != n_modes:
                    raise ValueError(f"1D participation_factors length ({gam.size}) != n_modes ({n_modes}).")
                gam_2d = np.zeros((n_modes, 3), dtype=float)
                gam_2d[:, 0] = gam
                self.participation_factors: np.ndarray | None = gam_2d
            elif gam.ndim == 2 and gam.shape[0] == n_modes:
                self.participation_factors = gam
            else:
                raise ValueError(f"participation_factors must have shape ({n_modes}, k), got {gam.shape}.")
        elif M_full is not None:
            r_mat = self._default_influence_vectors()
            if r_mat is not None:
                self.participation_factors = (phi_arr.T @ (M_full @ r_mat)) / self.M_modal[:, None]
            else:
                self.participation_factors = None
        else:
            self.participation_factors = None

        self.residual_flexibility = residual_flexibility

    @classmethod
    def from_eigenmodes(
        cls,
        eigenmodes: EigenmodeResult,
        damping: DampingModel | float | Sequence[float] | np.ndarray = 0.02,
        stress_modes: np.ndarray | None = None,
        residual_flexibility: np.ndarray | sp.spmatrix | Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> ModalRandomVibration:
        """Construct solver directly from an :class:`EigenmodeResult`."""
        return cls(
            frequencies_hz=eigenmodes.frequencies_hz,
            phi=eigenmodes.modes,
            M_modal=np.ones(eigenmodes.frequencies_hz.size, dtype=float),
            damping=damping,
            stress_modes=stress_modes,
            participation_factors=eigenmodes.participation_factors,
            residual_flexibility=residual_flexibility,
        )

    def _default_influence_vectors(self) -> np.ndarray | None:
        """Build standard 3D translational rigid-body vectors ``R`` when inferable."""
        if self.influence_vectors is not None:
            return self.influence_vectors
        if self.n_dofs == 1:
            return np.ones((1, 1), dtype=float)
        if self.n_dofs % 3 == 0:
            r = np.zeros((self.n_dofs, 3), dtype=float)
            r[0::3, 0] = 1.0
            r[1::3, 1] = 1.0
            r[2::3, 2] = 1.0
            return r
        return None

    def _build_modal_load_tensor(
        self,
        excitations: Sequence[ExcitationBase],
        omega_grid: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
        """Build generalized modal load tensor ``L[mode, in, freq]`` and optional residual/base terms."""
        n_in = len(excitations)
        n_f = omega_grid.size
        L = np.zeros((self.n_modes, n_in, n_f), dtype=complex)
        has_force = False
        force_cols = np.zeros((self.n_dofs, n_in), dtype=complex)
        base_r_cols = np.zeros((self.n_dofs, n_in, n_f), dtype=complex)
        has_base = False

        for i, exc in enumerate(excitations):
            if isinstance(exc, NodalForceExcitation):
                if exc.force_vector.size != self.n_dofs:
                    raise ValueError(
                        f"NodalForceExcitation[{i}] force_vector size ({exc.force_vector.size}) != n_dofs ({self.n_dofs})."
                    )
                has_force = True
                force_cols[:, i] = exc.force_vector
                modal_f = self.phi.T.conj() @ exc.force_vector
                L[:, i, :] = modal_f[:, None]
            elif isinstance(exc, BaseExcitation):
                has_base = True
                if self.participation_factors is None:
                    raise ValueError(
                        "BaseExcitation requires either participation_factors or M_full to be provided."
                    )
                n_dir = self.participation_factors.shape[1]
                if n_dir == 1:
                    gamma_i = self.participation_factors[:, 0]
                    r_dir = np.ones(self.n_dofs, dtype=float)
                else:
                    d_vec = exc.direction[:n_dir]
                    gamma_i = self.participation_factors @ d_vec
                    r_mat = self._default_influence_vectors()
                    if r_mat is not None:
                        r_dir = r_mat @ exc.direction[: r_mat.shape[1]]
                    else:
                        r_dir = self.phi @ gamma_i

                lam = exc.lambda_exp
                kin_factor = (1j * omega_grid) ** lam if lam > 0 else np.ones_like(omega_grid, dtype=complex)
                # Modal RHS = -M_modal * gamma_i * (i*w)^lam * a_base
                L[:, i, :] = -(self.M_modal * gamma_i)[:, None] * kin_factor[None, :]
                base_r_cols[:, i, :] = r_dir[:, None] * kin_factor[None, :]
            else:
                raise TypeError(f"Unsupported excitation type {type(exc).__name__}.")

        res_disp: np.ndarray | None = None
        if self.residual_flexibility is not None:
            res_disp = np.zeros((self.n_dofs, n_in, n_f), dtype=complex)
            if has_force:
                if callable(self.residual_flexibility):
                    r_f = np.asarray(self.residual_flexibility(force_cols), dtype=complex)
                else:
                    r_f = np.asarray(self.residual_flexibility @ force_cols, dtype=complex)
                res_disp += r_f[:, :, None]
            if has_base and self.M_full is not None:
                # Quasi-static inertial load residual for base acceleration
                for i, exc in enumerate(excitations):
                    if isinstance(exc, BaseExcitation):
                        r_mat = self._default_influence_vectors()
                        if r_mat is not None:
                            r_dir = r_mat @ exc.direction[: r_mat.shape[1]]
                            f_inert = -(self.M_full @ r_dir)
                            if callable(self.residual_flexibility):
                                u_res_i = np.asarray(self.residual_flexibility(f_inert[:, None]), dtype=complex).reshape(-1)
                            else:
                                u_res_i = np.asarray(self.residual_flexibility @ f_inert, dtype=complex).reshape(-1)
                            lam = exc.lambda_exp
                            kin = (1j * omega_grid) ** lam if lam > 0 else np.ones_like(omega_grid, dtype=complex)
                            res_disp[:, i, :] += u_res_i[:, None] * kin[None, :]

        return L, res_disp, (base_r_cols if has_base else None)

    def solve(
        self,
        excitation: ExcitationBase | Sequence[ExcitationBase] | CrossPSDMatrix,
        freqs_hz: Sequence[float] | np.ndarray | None = None,
        *,
        f_min: float | None = None,
        f_max: float | None = None,
        grid_type: Literal["adaptive", "abaqus", "log", "uniform"] = "adaptive",
        n_points: int = 200,
        response_frame: ResponseFrame = "relative",
        validate_input_psd: bool = True,
    ) -> RandomVibrationResult:
        """Solve the modal random vibration problem.

        Args:
            excitation: A single excitation, a sequence of uncorrelated excitations,
                or a :class:`CrossPSDMatrix`.
            freqs_hz: Explicit frequency grid in Hz. If ``None``, built automatically
                from ``[f_min, f_max]`` using ``grid_type``.
            f_min: Lower frequency bound when ``freqs_hz is None``.
            f_max: Upper frequency bound when ``freqs_hz is None``.
            grid_type: Grid builder when ``freqs_hz is None`` (default ``"adaptive"``).
            n_points: Base point count for automatic grid generation.
            response_frame: ``"relative"`` (default, standard relative DOF motion)
                or ``"absolute"`` (adds base motion to acceleration/velocity/disp
                for :class:`BaseExcitation`).
            validate_input_psd: Whether to verify Hermitian + PSD eigenvalues on ``S_in``.

        Returns:
            :class:`RandomVibrationResult`.
        """
        frame = response_frame.strip().lower()
        if frame not in ("relative", "absolute"):
            raise ValueError(f"response_frame must be 'relative' or 'absolute', got {response_frame!r}.")

        if isinstance(excitation, CrossPSDMatrix):
            cpsd = excitation
            excs = cpsd.excitations
        elif isinstance(excitation, ExcitationBase):
            excs = (excitation,)
            cpsd = CrossPSDMatrix.uncorrelated(excs)
        else:
            excs = tuple(excitation)
            cpsd = CrossPSDMatrix.uncorrelated(excs)

        if not excs:
            raise ValueError("ModalRandomVibration requires explicit excitation objects.")

        if freqs_hz is not None:
            f_grid = _as_1d_frequencies(freqs_hz, allow_zero=True)
        else:
            if f_min is None or f_max is None:
                raise ValueError("Provide either freqs_hz or both f_min and f_max.")
            if grid_type == "adaptive":
                f_grid = adaptive_modal_frequency_grid(self.frequencies_hz, self.damping, f_min, f_max, n_base=n_points)
            elif grid_type == "abaqus":
                f_grid = abaqus_frequency_grid(self.frequencies_hz, f_min, f_max, n_points=n_points)
            elif grid_type == "log":
                f_grid = log_frequency_grid(f_min, f_max, n_points)
            elif grid_type == "uniform":
                f_grid = uniform_frequency_grid(f_min, f_max, n_points)
            else:
                raise ValueError(f"Unknown grid_type {grid_type!r}.")

        omega_grid = 2.0 * np.pi * f_grid
        s_in = cpsd.evaluate(f_grid, validate=validate_input_psd)  # (n_in, n_in, n_f)
        h_modal = self.damping.modal_transfer_function(self.omega_modes, omega_grid, self.M_modal)  # (n_modes, n_f)
        L, res_disp, base_r_cols = self._build_modal_load_tensor(excs, omega_grid)  # (n_modes, n_in, n_f)

        # Modal transfer from inputs: T_q[m, i, f] = H_m(f) * L[m, i, f]
        t_q = h_modal[:, None, :] * L
        # Modal generalized PSD: S_qq[a, b, f] = sum_ij T_q[a, i, f] * S_in[i, j, f] * conj(T_q[b, j, f])
        s_qq = np.einsum("aif,ijf,bjf->abf", t_q, s_in, t_q.conj())
        # Enforce exact Hermitian symmetry against floating-point roundoff
        s_qq = 0.5 * (s_qq + np.swapaxes(s_qq.conj(), 0, 1))

        if res_disp is None and frame == "relative":
            # Fast O(n_dofs * n_modes^2 * n_f) auto-PSD without materializing (n_dofs, n_dofs, n_f)
            phi_s = np.einsum("dm,mnf->dnf", self.phi, s_qq)
            disp_psd = np.maximum(np.real(np.sum(phi_s * self.phi[:, :, None], axis=1)), 0.0)
            vel_psd = (omega_grid[None, :] ** 2) * disp_psd
            acc_psd = (omega_grid[None, :] ** 4) * disp_psd

            def _cross_eval(idx_i: np.ndarray, idx_j: np.ndarray) -> np.ndarray:
                return np.einsum("im,mnf,jn->ijf", self.phi[idx_i], s_qq, self.phi[idx_j].conj())
        else:
            # Build input-to-displacement FRF H_u[d, i, f] = sum_m Phi[d, m] * t_q[m, i, f] (+ residual)
            h_u = np.einsum("dm,mif->dif", self.phi, t_q)
            if res_disp is not None:
                h_u = h_u + res_disp

            if frame == "absolute" and base_r_cols is not None:
                # For base acceleration inputs (lam=0), total acceleration is a_rel + r_dir * a_base:
                # -w^2 u_rel + r_dir * a_base = (-w^2 H_u + base_r_cols) * a_base.
                # For velocity/displacement base inputs, base_r_cols = r_dir * (i*w)^lam, while u_base = r_dir * (i*w)^(lam-2)
                # so total displacement = H_u - base_r_cols / w^2 (except at w=0 where u_abs = u_base).
                h_acc = -(omega_grid[None, None, :] ** 2) * h_u + base_r_cols
                with np.errstate(divide="ignore", invalid="ignore"):
                    iw = 1j * omega_grid[None, None, :]
                    h_vel = np.where(omega_grid[None, None, :] > 0.0, h_acc / iw, 0.0)
                    h_disp = np.where(omega_grid[None, None, :] > 0.0, -h_acc / (omega_grid[None, None, :] ** 2), h_u)
                disp_psd = np.maximum(np.real(np.einsum("dif,ijf,djf->df", h_disp, s_in, h_disp.conj())), 0.0)
                vel_psd = np.maximum(np.real(np.einsum("dif,ijf,djf->df", h_vel, s_in, h_vel.conj())), 0.0)
                acc_psd = np.maximum(np.real(np.einsum("dif,ijf,djf->df", h_acc, s_in, h_acc.conj())), 0.0)
                h_cross = h_disp
            else:
                disp_psd = np.maximum(np.real(np.einsum("dif,ijf,djf->df", h_u, s_in, h_u.conj())), 0.0)
                vel_psd = (omega_grid[None, :] ** 2) * disp_psd
                acc_psd = (omega_grid[None, :] ** 4) * disp_psd
                h_cross = h_u

            def _cross_eval(idx_i: np.ndarray, idx_j: np.ndarray) -> np.ndarray:
                return np.einsum("iaf,abf,jbf->ijf", h_cross[idx_i], s_in, h_cross[idx_j].conj())

        stress_psd: np.ndarray | None = None
        s_vm_psd: np.ndarray | None = None
        rmises: np.ndarray | None = None
        if self.stress_modes is not None:
            stress_psd = np.einsum("pak,abf,pbl->pklf", self.stress_modes, s_qq, self.stress_modes.conj())
            stress_psd = 0.5 * (stress_psd + np.swapaxes(stress_psd.conj(), 1, 2))
            s_vm_psd = von_mises_psd(stress_psd)
            rmises = np.asarray(segalman_von_mises_rms(stress_psd, f_grid), dtype=float)

        return RandomVibrationResult(
            frequencies_hz=f_grid,
            disp_psd=disp_psd,
            vel_psd=vel_psd,
            acc_psd=acc_psd,
            modal_psd=s_qq,
            stress_psd=stress_psd,
            rmises=rmises,
            von_mises_psd=s_vm_psd,
            input_cross_psd=s_in,
            response_frame=frame,  # type: ignore[arg-type]
            _cross_evaluator=_cross_eval,
        )


class DirectFRFRandomVibration:
    """Direct Frequency Response Function (FRF) PSD random vibration solver.

    Computes exact response CPSD ``S_yy(f) = H(f) S_in(f) H(f)^H`` from a
    precomputed complex FRF tensor ``H(f)`` of shape ``(n_out, n_in, n_f)``
    (e.g., from direct harmonic analysis or analytical closed-form FRFs).

    Args:
        frequencies_hz: Frequency grid in Hz, shape ``(n_f,)``.
        frf: Complex FRF tensor of shape ``(n_out, n_in, n_f)`` (or ``(n_out, n_f)``
            for single-input systems).
        stress_frf: Optional complex Voigt stress FRF tensor of shape
            ``(n_pts, 6, n_in, n_f)``.
    """

    def __init__(
        self,
        frequencies_hz: Sequence[float] | np.ndarray,
        frf: np.ndarray,
        stress_frf: np.ndarray | None = None,
    ) -> None:
        self.frequencies_hz = _as_1d_frequencies(frequencies_hz, allow_zero=True)
        h = np.asarray(frf, dtype=complex)
        if h.ndim == 2:
            if h.shape[1] != self.frequencies_hz.size:
                raise ValueError(f"2D frf shape {h.shape} must be (n_out, {self.frequencies_hz.size}).")
            h = h[:, None, :]
        elif h.ndim != 3 or h.shape[2] != self.frequencies_hz.size:
            raise ValueError(f"frf must have shape (n_out, n_in, {self.frequencies_hz.size}), got {h.shape}.")
        self.frf = h
        self.n_out, self.n_in, self.n_freq = h.shape

        if stress_frf is not None:
            sf = np.asarray(stress_frf, dtype=complex)
            if sf.ndim != 4 or sf.shape[1:] != (6, self.n_in, self.n_freq):
                raise ValueError(
                    f"stress_frf must have shape (n_pts, 6, {self.n_in}, {self.n_freq}), got {sf.shape}."
                )
            self.stress_frf: np.ndarray | None = sf
        else:
            self.stress_frf = None

    def solve(
        self,
        cross_psd: CrossPSDMatrix | np.ndarray,
        *,
        tol: float = 1e-10,
        validate_input_psd: bool = True,
    ) -> RandomVibrationResult:
        """Evaluate output displacement/velocity/acceleration/stress PSDs from ``S_in(f)``."""
        if isinstance(cross_psd, CrossPSDMatrix):
            if cross_psd.n_inputs != self.n_in:
                raise ValueError(f"CrossPSDMatrix has {cross_psd.n_inputs} inputs, FRF expects {self.n_in}.")
            s_in = cross_psd.evaluate(self.frequencies_hz, validate=validate_input_psd)
        else:
            s_in = np.asarray(cross_psd, dtype=complex)
            if s_in.shape != (self.n_in, self.n_in, self.n_freq):
                raise ValueError(
                    f"cross_psd array must have shape ({self.n_in}, {self.n_in}, {self.n_freq}), got {s_in.shape}."
                )
            if validate_input_psd:
                for k in range(self.n_freq):
                    _validate_hermitian_psd_slice(s_in[:, :, k], k, float(self.frequencies_hz[k]), tol)

        omega = 2.0 * np.pi * self.frequencies_hz
        disp_psd = np.maximum(np.real(np.einsum("iaf,abf,ibf->if", self.frf, s_in, self.frf.conj())), 0.0)
        vel_psd = (omega[None, :] ** 2) * disp_psd
        acc_psd = (omega[None, :] ** 4) * disp_psd

        def _cross_eval(idx_i: np.ndarray, idx_j: np.ndarray) -> np.ndarray:
            return np.einsum("iaf,abf,jbf->ijf", self.frf[idx_i], s_in, self.frf[idx_j].conj())

        stress_psd: np.ndarray | None = None
        s_vm_psd: np.ndarray | None = None
        rmises: np.ndarray | None = None
        if self.stress_frf is not None:
            stress_psd = np.einsum("pkaf,abf,plbf->pklf", self.stress_frf, s_in, self.stress_frf.conj())
            stress_psd = 0.5 * (stress_psd + np.swapaxes(stress_psd.conj(), 1, 2))
            s_vm_psd = von_mises_psd(stress_psd)
            rmises = np.asarray(segalman_von_mises_rms(stress_psd, self.frequencies_hz), dtype=float)

        return RandomVibrationResult(
            frequencies_hz=self.frequencies_hz,
            disp_psd=disp_psd,
            vel_psd=vel_psd,
            acc_psd=acc_psd,
            modal_psd=None,
            stress_psd=stress_psd,
            rmises=rmises,
            von_mises_psd=s_vm_psd,
            input_cross_psd=s_in,
            response_frame="relative",
            _cross_evaluator=_cross_eval,
        )


# ---------------------------------------------------------------------------
# 7. Lightweight Case-File Bridge Helper
# ---------------------------------------------------------------------------


def random_vibration_from_case(
    omega_modes: np.ndarray,
    modes: np.ndarray,
    participation_factors: np.ndarray,
    case: RandomResponseCase,
    stress_modes: np.ndarray | None = None,
    frequencies_hz: np.ndarray | None = None,
    bias_model: str = "abaqus",
) -> RandomVibrationResult:
    """Bridge a validated :class:`RandomResponseCase` to :class:`ModalRandomVibration`.

    Preserves exact :func:`~wht_iga.dynamics.response_spectrum.random_response_from_case`
    semantics (log-log clamped PSD tables, ``BASE`` gravity scaling, and
    ``CORRELATED`` / ``UNCORRELATED`` load-case weights) while returning the full
    :class:`RandomVibrationResult` (PSDs, spectral moments, cross-PSD, and fatigue).
    """
    w = np.asarray(omega_modes, dtype=float).reshape(-1)
    f_modes = w / (2.0 * np.pi)
    damping = ModalDampingList(case.damping.ratios(f_modes))
    freqs = (
        _as_1d_frequencies(frequencies_hz, allow_zero=False)
        if frequencies_hz is not None
        else abaqus_frequency_grid(
            f_modes,
            case.lower_frequency_hz,
            case.upper_frequency_hz,
            n_points=case.n_points,
            bias=case.bias,
            scale=case.frequency_scale,  # type: ignore[arg-type]
            model=bias_model,
        )
    )

    solver = ModalRandomVibration(
        frequencies_hz=f_modes,
        phi=modes,
        damping=damping,
        stress_modes=stress_modes,
        participation_factors=participation_factors,
    )

    # Expand RandomResponseCase correlations into physical uncorrelated modal-weighted inputs
    excitations: list[BaseExcitation] = []
    for corr in case.correlations:
        psd_def = case.psds[corr.psd]
        if psd_def.type is PSDType.FORCE:
            raise ValueError("FORCE PSDs need applied-force vectors; only base-motion cases are bridged here.")
        # Validate via PSD.value(freqs) to preserve CaseFileError checks on DB/complex tables
        vals = psd_def.value(freqs)
        if np.any(vals.imag != 0.0):
            raise ValueError(f"psd '{psd_def.name}': complex ordinates are not supported in case bridge.")
        tab = np.asarray(psd_def.rows, dtype=float)
        curve = PSDTable(tab[:, 2], tab[:, 0], scale="log-log", extrapolation="clamp", name=psd_def.name)
        g_scale = float(psd_def.gravity) if psd_def.type is PSDType.BASE else 1.0

        for lc, (wre, wim) in corr.weights.items():
            w_amp = math.hypot(wre, wim)
            if w_amp == 0.0:
                continue
            moves = [m for m in case.base_motions if m.load_case == lc]
            kind_str = BaseKind(moves[0].kind).value.lower()
            if corr.type == "CORRELATED":
                d_vec = np.zeros(3, dtype=float)
                for m in moves:
                    d_vec[m.dof - 1] += 1.0
                d_norm = float(np.linalg.norm(d_vec))
                if d_norm > 0.0:
                    excitations.append(
                        BaseExcitation(
                            dof_or_direction=d_vec,
                            psd=curve,
                            gravity_scale=g_scale * w_amp * d_norm,
                            kind=kind_str,  # type: ignore[arg-type]
                        )
                    )
            else:
                for m in moves:
                    excitations.append(
                        BaseExcitation(
                            dof_or_direction=m.dof - 1,
                            psd=curve,
                            gravity_scale=g_scale * w_amp,
                            kind=kind_str,  # type: ignore[arg-type]
                        )
                    )

    if not excitations:
        raise ValueError("RandomResponseCase produced no active excitations.")

    return solver.solve(CrossPSDMatrix.uncorrelated(excitations), freqs_hz=freqs)
