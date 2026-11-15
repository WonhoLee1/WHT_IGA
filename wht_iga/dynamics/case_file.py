"""Coreform-style response-spectrum / random-response case files (plan add_1 G-3).

The JSON schemas follow the Coreform IGA for Abaqus user guide pages
``response-spectrum/case-files`` and ``response-spectrum/random-response``:

* a response-spectrum case has exactly the keys ``spectra``, ``excitations``,
  ``modal_summation`` and ``directional_summation`` plus ``damping``;
* a random-response case has ``psds``, ``base_motions``, ``correlations``,
  ``lower_frequency_hz``, ``upper_frequency_hz``, ``n_points``, ``bias``,
  ``frequency_scale`` and ``damping``.

Everything here is Kratos-free.  Anything the user guide does not define is
rejected with :class:`CaseFileError` instead of being guessed (for example the
band-number table of ``DB`` power spectral densities).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np


class CaseFileError(ValueError):
    """Raised when a case file violates the documented schema or its rules."""


class SpectrumType(str, Enum):
    """Physical type of the magnitudes tabulated in a spectrum."""

    ACCELERATION = "ACCELERATION"
    VELOCITY = "VELOCITY"
    DISPLACEMENT = "DISPLACEMENT"
    G = "G"


class ModalSummation(str, Enum):
    """Rule combining the peak responses of the individual modes."""

    ABS = "ABS"
    SRSS = "SRSS"
    NRL = "NRL"
    TENP = "TENP"
    CQC = "CQC"
    GRP = "GRP"
    DSC = "DSC"


class DirectionalSummation(str, Enum):
    """Rule combining the excitation directions."""

    ALGEBRAIC = "ALGEBRAIC"
    SRSS = "SRSS"
    R40 = "R40"
    R30 = "R30"


class PSDType(str, Enum):
    """Kind of a power-spectral-density definition (``*PSD-DEFINITION``)."""

    BASE = "BASE"
    FORCE = "FORCE"
    DB = "DB"


class BaseKind(str, Enum):
    """Quantity a base motion prescribes (sets the ``(2 pi f)**(2 lambda)`` PSD factor)."""

    ACCELERATION = "ACCELERATION"
    VELOCITY = "VELOCITY"
    DISPLACEMENT = "DISPLACEMENT"


def _enum(cls: type[Enum], value: Any, what: str) -> Any:
    """Parse ``value`` (case-insensitive) into ``cls`` or raise ``CaseFileError``."""
    if isinstance(value, cls):
        return value
    try:
        return cls(str(value).strip().upper())
    except ValueError:
        allowed = ", ".join(m.value for m in cls)
        raise CaseFileError(f"{what}: '{value}' is not one of {allowed}") from None


def _number(value: Any, what: str) -> float:
    """Return ``value`` as a finite float or raise ``CaseFileError``."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise CaseFileError(f"{what}: expected a number, got {value!r}") from None
    if not math.isfinite(out):
        raise CaseFileError(f"{what}: expected a finite number, got {value!r}")
    return out


# ---------------------------------------------------------------------------
# Response spectrum
# ---------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class SpectrumCurve:
    """Magnitude-versus-frequency curve at one damping ratio.

    Attributes:
        damping: Fraction of critical damping of this curve.
        freqs_hz: Strictly increasing positive frequencies, shape ``(n,)``.
        mags: Positive magnitudes in the units of the spectrum type, shape ``(n,)``.
    """

    damping: float
    freqs_hz: np.ndarray
    mags: np.ndarray


@dataclass(frozen=True, eq=False)
class Spectrum:
    """A named response spectrum: curves sorted by damping ratio.

    Attributes:
        name: Name used by excitations to refer to the spectrum.
        type: Physical type of the tabulated magnitudes.
        gravity: Acceleration of gravity; only used (and required) for ``G``.
        curves: One :class:`SpectrumCurve` per damping ratio, ascending damping.
    """

    name: str
    type: SpectrumType
    gravity: float | None
    curves: tuple[SpectrumCurve, ...]

    @classmethod
    def from_rows(
        cls,
        name: str,
        rows: Any,
        type: SpectrumType | str = SpectrumType.ACCELERATION,  # noqa: A002 - mirrors the JSON key
        gravity: float | None = None,
    ) -> Spectrum:
        """Build a spectrum from ``[magnitude, frequency, damping]`` rows.

        Args:
            name: Spectrum name.
            rows: Iterable of ``(magnitude, frequency_hz, damping_ratio)`` triples
                in the order of a ``*SPECTRUM`` data line.  Rows of one damping
                ratio must be contiguous and form one curve.
            type: Spectrum type (``ACCELERATION`` default).
            gravity: Gravity, required when ``type`` is ``G``.

        Returns:
            The validated :class:`Spectrum`.
        """
        stype = _enum(SpectrumType, type, f"spectrum '{name}' type")
        if stype is SpectrumType.G:
            if gravity is None:
                raise CaseFileError(f"spectrum '{name}': type G needs a 'gravity' value")
            gravity = _number(gravity, f"spectrum '{name}' gravity")
            if gravity <= 0.0:
                raise CaseFileError(f"spectrum '{name}': gravity must be positive")
        table = list(rows)
        if not table:
            raise CaseFileError(f"spectrum '{name}': 'rows' is empty")
        order: list[float] = []
        groups: dict[float, list[tuple[float, float]]] = {}
        for k, row in enumerate(table):
            if len(row) != 3:
                raise CaseFileError(f"spectrum '{name}' row {k}: expected [magnitude, frequency, damping]")
            mag = _number(row[0], f"spectrum '{name}' row {k} magnitude")
            freq = _number(row[1], f"spectrum '{name}' row {k} frequency")
            xi = _number(row[2], f"spectrum '{name}' row {k} damping")
            if xi < 0.0:
                raise CaseFileError(f"spectrum '{name}' row {k}: damping ratio must be >= 0")
            if not order or order[-1] != xi:
                if xi in groups:
                    raise CaseFileError(f"spectrum '{name}': damping ratio {xi} is repeated in separate row blocks")
                order.append(xi)
                groups[xi] = []
            groups[xi].append((freq, mag))
        curves = []
        for xi in sorted(groups):
            pts = groups[xi]
            f = np.array([p[0] for p in pts], dtype=float)
            m = np.array([p[1] for p in pts], dtype=float)
            if np.any(f <= 0.0) or np.any(np.diff(f) <= 0.0):
                raise CaseFileError(f"spectrum '{name}' damping {xi}: frequencies must be positive and strictly increasing")
            if np.any(m <= 0.0):
                raise CaseFileError(f"spectrum '{name}' damping {xi}: magnitudes must be positive")
            curves.append(SpectrumCurve(damping=xi, freqs_hz=f, mags=m))
        return cls(name=name, type=stype, gravity=gravity, curves=tuple(curves))

    def tabulated(self, freq_hz: np.ndarray, damping: np.ndarray | float) -> np.ndarray:
        """Interpolate the magnitude in the *supplied* spectrum units.

        Within a curve the magnitude is interpolated linearly on a log-log scale
        and held constant beyond the table ends.  Between two damping curves the
        weight is linear in the damping ratio and the magnitude interpolates
        logarithmically; damping outside the tabulated range is clamped.  A
        ``G`` spectrum is multiplied by ``gravity`` (it becomes an acceleration).

        Args:
            freq_hz: Frequencies in Hz (positive), shape ``(n,)``.
            damping: Damping ratio per frequency, shape ``(n,)``, or a scalar.

        Returns:
            Magnitudes, shape ``(n,)``.
        """
        f = np.atleast_1d(np.asarray(freq_hz, dtype=float))
        if np.any(f <= 0.0):
            raise CaseFileError("spectrum evaluation needs positive frequencies")
        xi = np.broadcast_to(np.asarray(damping, dtype=float), f.shape)
        logf = np.log(f)
        per_curve = np.stack(
            [np.interp(logf, np.log(c.freqs_hz), np.log(c.mags)) for c in self.curves]
        )  # (n_curves, n): log-magnitude of every curve at every frequency
        xis = np.array([c.damping for c in self.curves])
        if len(self.curves) == 1:
            logs = per_curve[0]
        else:
            xi_c = np.clip(xi, xis[0], xis[-1])
            hi = np.clip(np.searchsorted(xis, xi_c, side="right"), 1, len(xis) - 1)
            lo = hi - 1
            w = (xi_c - xis[lo]) / (xis[hi] - xis[lo])
            idx = np.arange(f.size)
            logs = (1.0 - w) * per_curve[lo, idx] + w * per_curve[hi, idx]
        out = np.exp(logs)
        if self.type is SpectrumType.G:
            out = out * float(self.gravity)  # type: ignore[arg-type]
        return out

    def at(self, kind: str, freq_hz: np.ndarray, damping: np.ndarray | float) -> np.ndarray:
        """Spectral value of ``kind`` (``'D'``, ``'V'`` or ``'A'``) at ``freq_hz``.

        Uses the undamped relations ``S_D = S_V / w = S_A / w**2`` (Abaqus), with
        ``w = 2 pi f``, to convert from whichever type was supplied.
        """
        base = self.tabulated(freq_hz, damping)
        w = 2.0 * np.pi * np.atleast_1d(np.asarray(freq_hz, dtype=float))
        supplied = "A" if self.type in (SpectrumType.ACCELERATION, SpectrumType.G) else self.type.value[0]
        s_a = {"A": base, "V": base * w, "D": base * w * w}[supplied]
        kind = kind.upper()
        if kind == "A":
            return s_a
        if kind == "V":
            return s_a / w
        if kind == "D":
            return s_a / (w * w)
        raise CaseFileError(f"spectrum kind must be 'D', 'V' or 'A', got {kind!r}")


@dataclass(frozen=True)
class Excitation:
    """One base-excitation line of a response-spectrum case.

    Attributes:
        spectrum: Name of an entry of the case's ``spectra``.
        direction: Unit direction cosines ``(dx, dy, dz)``.
        scale: Factor multiplying the spectrum magnitudes.
        duration: Strong-motion duration, required by the ``DSC`` rule.
    """

    spectrum: str
    direction: tuple[float, float, float]
    scale: float = 1.0
    duration: float | None = None


@dataclass(frozen=True)
class ModalDamping:
    """Modal damping given as a uniform ratio, mode ranges or a frequency table.

    Exactly one of the three fields is set.  Modes not covered by ``mode_ranges``
    get zero damping, as in Abaqus; ``frequency_table`` is linear in frequency
    and constant beyond its ends.
    """

    uniform: float | None = None
    mode_ranges: tuple[tuple[int, int, float], ...] | None = None
    frequency_table: tuple[tuple[float, float], ...] | None = None

    @classmethod
    def parse(cls, raw: Any) -> ModalDamping:
        """Parse the ``damping`` JSON value (number or one-key dict)."""
        if isinstance(raw, bool) or raw is None:
            raise CaseFileError("damping: a number or {'mode_ranges'|'frequency_table': ...} is required")
        if isinstance(raw, (int, float)):
            xi = _number(raw, "damping")
            if xi < 0.0:
                raise CaseFileError("damping: ratio must be >= 0")
            return cls(uniform=xi)
        if isinstance(raw, dict) and len(raw) == 1:
            key, val = next(iter(raw.items()))
            if key == "mode_ranges":
                rngs = []
                for r in val:
                    if len(r) != 3:
                        raise CaseFileError("damping.mode_ranges: entries are [first, last, ratio]")
                    first, last = int(r[0]), int(r[1])
                    xi = _number(r[2], "damping.mode_ranges ratio")
                    if first < 1 or last < first or xi < 0.0:
                        raise CaseFileError(f"damping.mode_ranges: invalid entry {list(r)}")
                    rngs.append((first, last, xi))
                return cls(mode_ranges=tuple(rngs))
            if key == "frequency_table":
                tab = [(_number(r[0], "damping.frequency_table f"), _number(r[1], "damping.frequency_table ratio")) for r in val]
                if not tab or any(b[0] <= a[0] for a, b in zip(tab, tab[1:], strict=False)):
                    raise CaseFileError("damping.frequency_table: frequencies must be strictly increasing")
                return cls(frequency_table=tuple(tab))
        raise CaseFileError("damping: expected a number, {'mode_ranges': ...} or {'frequency_table': ...}")

    def ratios(self, freqs_hz: np.ndarray) -> np.ndarray:
        """Damping ratio of every mode (mode numbers are ``1..len(freqs_hz)``)."""
        f = np.asarray(freqs_hz, dtype=float).reshape(-1)
        if self.uniform is not None:
            return np.full(f.shape, self.uniform)
        if self.mode_ranges is not None:
            out = np.zeros(f.shape)
            for first, last, xi in self.mode_ranges:
                out[first - 1 : min(last, f.size)] = xi
            return out
        tab = np.asarray(self.frequency_table, dtype=float)
        return np.interp(f, tab[:, 0], tab[:, 1])

    def to_json(self) -> Any:
        """Inverse of :meth:`parse`."""
        if self.uniform is not None:
            return self.uniform
        if self.mode_ranges is not None:
            return {"mode_ranges": [list(r) for r in self.mode_ranges]}
        return {"frequency_table": [list(r) for r in self.frequency_table]}  # type: ignore[union-attr]


@dataclass(frozen=True, eq=False)
class ResponseSpectrumCase:
    """A complete response-spectrum analysis definition."""

    spectra: dict[str, Spectrum]
    excitations: tuple[Excitation, ...]
    modal_summation: ModalSummation
    directional_summation: DirectionalSummation
    damping: ModalDamping


# ---------------------------------------------------------------------------
# Random response
# ---------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class PSD:
    """A frequency function ``P(f)`` (``*PSD-DEFINITION``).

    Attributes:
        name: PSD name.
        type: ``BASE`` (ordinates in ``g**2`` per Hz), ``FORCE`` (power units) or ``DB``.
        gravity: Acceleration of gravity, required for ``BASE``.
        db_reference: Reference level for ``DB`` tables (not evaluable here).
        rows: ``(real, imag, frequency_hz)`` triples, ascending frequency.
    """

    name: str
    type: PSDType
    gravity: float | None
    db_reference: float | None
    rows: tuple[tuple[float, float, float], ...]

    def value(self, freq_hz: np.ndarray) -> np.ndarray:
        """Complex ordinate at ``freq_hz``: log-log interpolation, constant beyond the table.

        The imaginary part must be identically zero or strictly one-signed
        (otherwise a log-log interpolation is undefined and the PSD is rejected).
        """
        if self.type is PSDType.DB:
            raise CaseFileError(
                f"psd '{self.name}': type DB needs the band-number to frequency table, "
                "which the Coreform user guide does not tabulate; supply a BASE or FORCE PSD"
            )
        f = np.atleast_1d(np.asarray(freq_hz, dtype=float))
        tab = np.asarray(self.rows, dtype=float)
        logf = np.log(tab[:, 2])
        re = np.exp(np.interp(np.log(f), logf, np.log(tab[:, 0])))
        if np.all(tab[:, 1] == 0.0):
            return re.astype(complex)
        if np.all(tab[:, 1] > 0.0) or np.all(tab[:, 1] < 0.0):
            sign = 1.0 if tab[0, 1] > 0.0 else -1.0
            im = sign * np.exp(np.interp(np.log(f), logf, np.log(np.abs(tab[:, 1]))))
            return re + 1j * im
        raise CaseFileError(f"psd '{self.name}': imaginary parts must be all zero or one-signed")


@dataclass(frozen=True)
class BaseMotion:
    """One base-motion degree of freedom of a load case.

    Attributes:
        dof: Global degree of freedom 1-3 (translations).
        load_case: Load case number.
        kind: ``ACCELERATION``, ``VELOCITY`` or ``DISPLACEMENT``.
    """

    dof: int
    load_case: int
    kind: str


@dataclass(frozen=True)
class Correlation:
    """Ties one PSD to load cases with complex weights.

    Attributes:
        psd: Name of an entry of ``psds``.
        type: ``CORRELATED`` or ``UNCORRELATED``.
        weights: ``load_case -> (real, imag)`` weight.
    """

    psd: str
    type: str
    weights: dict[int, tuple[float, float]]


@dataclass(frozen=True, eq=False)
class RandomResponseCase:
    """A complete random-response analysis definition."""

    psds: dict[str, PSD]
    base_motions: tuple[BaseMotion, ...]
    correlations: tuple[Correlation, ...]
    lower_frequency_hz: float
    upper_frequency_hz: float
    n_points: int
    bias: float
    frequency_scale: str
    damping: ModalDamping


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_RS_KEYS = {"spectra", "excitations", "modal_summation", "directional_summation", "damping"}
_RANDOM_KEYS = {
    "psds", "base_motions", "correlations", "lower_frequency_hz", "upper_frequency_hz",
    "n_points", "bias", "frequency_scale", "damping",
}


def _check_keys(raw: dict, required: set[str], what: str) -> None:
    """Reject missing and unknown top-level keys."""
    missing = sorted(required - set(raw))
    extra = sorted(set(raw) - required)
    if missing:
        raise CaseFileError(f"{what}: missing required key(s) {missing}")
    if extra:
        raise CaseFileError(f"{what}: unknown key(s) {extra}")


def _parse_response_spectrum(raw: dict) -> ResponseSpectrumCase:
    """Parse and validate a response-spectrum case dictionary."""
    _check_keys(raw, _RS_KEYS, "response-spectrum case")
    spectra = {
        name: Spectrum.from_rows(name, body.get("rows", []), body.get("type", "ACCELERATION"), body.get("gravity"))
        for name, body in raw["spectra"].items()
    }
    if not spectra:
        raise CaseFileError("spectra: at least one spectrum is required")
    mod = _enum(ModalSummation, raw["modal_summation"], "modal_summation")
    dirsum = _enum(DirectionalSummation, raw["directional_summation"], "directional_summation")
    exc_raw = raw["excitations"]
    if not 1 <= len(exc_raw) <= 3:
        raise CaseFileError("excitations: between one and three entries are required")
    excs = []
    for k, e in enumerate(exc_raw):
        if e.get("spectrum") not in spectra:
            raise CaseFileError(f"excitation {k}: spectrum '{e.get('spectrum')}' is not defined in spectra")
        d = tuple(_number(x, f"excitation {k} direction") for x in e.get("direction", ()))
        if len(d) != 3 or abs(math.sqrt(sum(x * x for x in d)) - 1.0) > 1e-8:
            raise CaseFileError(f"excitation {k}: direction must be a unit vector of length 3")
        scale = _number(e.get("scale", 1.0), f"excitation {k} scale")
        dur = e.get("duration")
        dur = None if dur is None else _number(dur, f"excitation {k} duration")
        if mod is ModalSummation.DSC and (dur is None or dur <= 0.0):
            raise CaseFileError(f"excitation {k}: DSC summation needs a positive duration on every excitation")
        excs.append(Excitation(spectrum=e["spectrum"], direction=d, scale=scale, duration=dur))
    for a in range(len(excs)):
        for b in range(a + 1, len(excs)):
            if abs(sum(x * y for x, y in zip(excs[a].direction, excs[b].direction, strict=True))) > 1e-8:
                raise CaseFileError(f"excitations {a} and {b}: directions must be mutually orthogonal")
    return ResponseSpectrumCase(
        spectra=spectra, excitations=tuple(excs), modal_summation=mod,
        directional_summation=dirsum, damping=ModalDamping.parse(raw["damping"]),
    )


def _parse_random(raw: dict) -> RandomResponseCase:
    """Parse and validate a random-response case dictionary."""
    _check_keys(raw, _RANDOM_KEYS, "random-response case")
    psds: dict[str, PSD] = {}
    for name, body in raw["psds"].items():
        ptype = _enum(PSDType, body.get("type", "BASE"), f"psd '{name}' type")
        grav = body.get("gravity")
        if ptype is PSDType.BASE and (grav is None or _number(grav, f"psd '{name}' gravity") <= 0.0):
            raise CaseFileError(f"psd '{name}': type BASE needs a positive 'gravity'")
        rows = tuple((_number(r[0], f"psd '{name}' re"), _number(r[1], f"psd '{name}' im"), _number(r[2], f"psd '{name}' f"))
                     for r in body.get("rows", []))
        if not rows:
            raise CaseFileError(f"psd '{name}': 'rows' is empty")
        if any(b[2] <= a[2] for a, b in zip(rows, rows[1:], strict=False)) or rows[0][2] <= 0.0:
            raise CaseFileError(f"psd '{name}': frequencies must be positive and strictly increasing")
        if any(r[0] <= 0.0 for r in rows) and ptype is not PSDType.DB:
            raise CaseFileError(f"psd '{name}': real parts must be positive (log-log interpolation)")
        psds[name] = PSD(name=name, type=ptype, gravity=None if grav is None else float(grav),
                         db_reference=body.get("db_reference"), rows=rows)
    if not psds:
        raise CaseFileError("psds: at least one PSD is required")
    motions = []
    for k, m in enumerate(raw["base_motions"]):
        dof, lc = int(m.get("dof", 0)), int(m.get("load_case", 0))
        kind = _enum(BaseKind, m.get("kind", "ACCELERATION"), f"base_motion {k} kind").value
        if not 1 <= dof <= 6:
            raise CaseFileError(f"base_motion {k}: dof must be 1..6")
        if dof > 3:
            raise CaseFileError(f"base_motion {k}: rotational base motion (dof {dof}) is not supported")
        motions.append(BaseMotion(dof=dof, load_case=lc, kind=kind))
    if not motions:
        raise CaseFileError("base_motions: at least one entry is required")
    for lc in {m.load_case for m in motions}:
        if len({m.kind for m in motions if m.load_case == lc}) > 1:
            raise CaseFileError(f"load case {lc}: all base motions of one load case must share a kind")
    corrs = []
    for k, c in enumerate(raw["correlations"]):
        if c.get("psd") not in psds:
            raise CaseFileError(f"correlation {k}: psd '{c.get('psd')}' is not defined in psds")
        ctype = str(c.get("type", "")).upper()
        if ctype not in ("CORRELATED", "UNCORRELATED"):
            raise CaseFileError(f"correlation {k}: type must be CORRELATED or UNCORRELATED")
        w = {int(lc): (_number(v[0], f"correlation {k} weight"), _number(v[1], f"correlation {k} weight"))
             for lc, v in c.get("weights", {}).items()}
        if not w:
            raise CaseFileError(f"correlation {k}: 'weights' is empty")
        for lc in w:
            if lc not in {m.load_case for m in motions}:
                raise CaseFileError(f"correlation {k}: load case {lc} has no base motion")
        corrs.append(Correlation(psd=c["psd"], type=ctype, weights=w))
    if not corrs:
        raise CaseFileError("correlations: at least one entry is required")
    lo = _number(raw["lower_frequency_hz"], "lower_frequency_hz")
    hi = _number(raw["upper_frequency_hz"], "upper_frequency_hz")
    if not 0.0 < lo < hi:
        raise CaseFileError("lower_frequency_hz / upper_frequency_hz: need 0 < lower < upper")
    n_pts = int(raw["n_points"])
    if n_pts < 2:
        raise CaseFileError("n_points must be >= 2")
    bias = _number(raw["bias"], "bias")
    if bias < 1.0:
        raise CaseFileError("bias must be >= 1 (1 = uniform spacing)")
    scale = str(raw["frequency_scale"]).upper()
    if scale not in ("LINEAR", "LOGARITHMIC"):
        raise CaseFileError("frequency_scale must be LINEAR or LOGARITHMIC")
    return RandomResponseCase(
        psds=psds, base_motions=tuple(motions), correlations=tuple(corrs), lower_frequency_hz=lo,
        upper_frequency_hz=hi, n_points=n_pts, bias=bias, frequency_scale=scale,
        damping=ModalDamping.parse(raw["damping"]),
    )


def parse_case(raw: dict) -> ResponseSpectrumCase | RandomResponseCase:
    """Parse a case dictionary; the key set decides spectrum versus random.

    Args:
        raw: Decoded JSON object.

    Returns:
        A :class:`ResponseSpectrumCase` or :class:`RandomResponseCase`.

    Raises:
        CaseFileError: On any schema or rule violation.
    """
    if not isinstance(raw, dict):
        raise CaseFileError("a case file must be a JSON object")
    if "psds" in raw:
        return _parse_random(raw)
    if "spectra" in raw:
        return _parse_response_spectrum(raw)
    raise CaseFileError("case file has neither 'spectra' (response spectrum) nor 'psds' (random response)")


def load_case_file(path: str | Path) -> ResponseSpectrumCase | RandomResponseCase:
    """Read and validate a case file from disk (UTF-8 JSON).

    Args:
        path: File path.

    Returns:
        The parsed case.
    """
    with open(path, encoding="utf-8") as fh:
        return parse_case(json.load(fh))


def case_to_dict(case: ResponseSpectrumCase | RandomResponseCase) -> dict[str, Any]:
    """Serialise a case to the documented JSON structure (inverse of :func:`parse_case`)."""
    if isinstance(case, ResponseSpectrumCase):
        spectra = {}
        for name, sp in case.spectra.items():
            rows = [[float(m), float(f), c.damping] for c in sp.curves for f, m in zip(c.freqs_hz, c.mags, strict=True)]
            spectra[name] = {"type": sp.type.value, "gravity": sp.gravity, "rows": rows}
        return {
            "spectra": spectra,
            "excitations": [
                {"spectrum": e.spectrum, "direction": list(e.direction), "scale": e.scale, "duration": e.duration}
                for e in case.excitations
            ],
            "modal_summation": case.modal_summation.value,
            "directional_summation": case.directional_summation.value,
            "damping": case.damping.to_json(),
        }
    return {
        "psds": {
            n: {"type": p.type.value, "gravity": p.gravity, "db_reference": p.db_reference, "rows": [list(r) for r in p.rows]}
            for n, p in case.psds.items()
        },
        "base_motions": [{"dof": m.dof, "load_case": m.load_case, "kind": m.kind} for m in case.base_motions],
        "correlations": [
            {"psd": c.psd, "type": c.type, "weights": {str(k): list(v) for k, v in c.weights.items()}}
            for c in case.correlations
        ],
        "lower_frequency_hz": case.lower_frequency_hz,
        "upper_frequency_hz": case.upper_frequency_hz,
        "n_points": case.n_points,
        "bias": case.bias,
        "frequency_scale": case.frequency_scale,
        "damping": case.damping.to_json(),
    }


def dump_case_file(case: ResponseSpectrumCase | RandomResponseCase, path: str | Path) -> None:
    """Write ``case`` as UTF-8 JSON (no BOM) that :func:`load_case_file` reads back."""
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(case_to_dict(case), fh, indent=1, ensure_ascii=False)
