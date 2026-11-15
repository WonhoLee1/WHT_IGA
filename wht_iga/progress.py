"""Timestamped stage/progress messages for long runs (stdout, flushed).

``WHT_PROGRESS=0`` silences them. Kept dependency-free so geometry, solver and
scripts can all use it without importing Kratos.
"""

from __future__ import annotations

import os
import sys
import time

_T0 = time.perf_counter()


def enabled() -> bool:
    return os.environ.get("WHT_PROGRESS", "1") not in ("0", "off", "false")


def log(msg: str) -> None:
    """Print ``[HH:MM:SS +elapsed] msg`` immediately."""
    if not enabled():
        return
    el = time.perf_counter() - _T0
    sys.stdout.write(f"[{time.strftime('%H:%M:%S')} +{el:7.1f}s] {msg}\n")
    sys.stdout.flush()


class Ticker:
    """Progress for a loop of ``total`` items: prints at ~``every_frac`` steps."""

    def __init__(self, label: str, total: int, every_frac: float = 0.1) -> None:
        self.label, self.total = label, max(int(total), 1)
        self.step = max(1, int(self.total * every_frac))
        self.done = 0
        self.t = time.perf_counter()
        self._next = self.step

    def advance(self, n: int = 1) -> None:
        self.done += n
        if self.done >= self._next or self.done >= self.total:
            dt = time.perf_counter() - self.t
            rate = self.done / dt if dt > 0 else 0.0
            eta = (self.total - self.done) / rate if rate > 0 else float("nan")
            log(f"{self.label}: {self.done}/{self.total} ({100 * self.done / self.total:.0f}%), "
                f"{dt:.0f}s elapsed, ETA {eta:.0f}s")
            while self._next <= self.done:
                self._next += self.step
