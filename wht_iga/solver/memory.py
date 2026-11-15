"""RAM probing, solver memory estimation, and automatic solver selection."""

from __future__ import annotations

import ctypes
import importlib.util
import math
import os
import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class SystemMemoryInfo:
    """Snapshot of host physical and available RAM."""

    total_bytes: int
    available_bytes: int
    used_percent: float
    source: str

    @property
    def total_gb(self) -> float:
        return self.total_bytes / (1024.0 ** 3)

    @property
    def available_gb(self) -> float:
        return self.available_bytes / (1024.0 ** 3)


@dataclass(frozen=True)
class SolverMemoryEstimate:
    """Estimated peak memory (bytes) for candidate linear solvers."""

    n_dof: int
    nnz: int
    degree: int
    pardiso_bytes: int
    amgcl_bytes: int

    @property
    def pardiso_gb(self) -> float:
        return self.pardiso_bytes / (1024.0 ** 3)

    @property
    def amgcl_gb(self) -> float:
        return self.amgcl_bytes / (1024.0 ** 3)


@dataclass(frozen=True)
class SolverSelectionDecision:
    """Result of memory-aware linear solver selection."""

    requested_solver: str
    selected_solver: str
    n_dof: int
    nnz: int
    degree: int
    est_pardiso_bytes: int
    est_amgcl_bytes: int
    avail_bytes: int
    total_bytes: int
    memory_budget_ratio: float
    budget_bytes: int
    dof_threshold: int
    switched: bool
    switched_reason: str | None = None

    @property
    def est_pardiso_gb(self) -> float:
        return self.est_pardiso_bytes / (1024.0 ** 3)

    @property
    def est_amgcl_gb(self) -> float:
        return self.est_amgcl_bytes / (1024.0 ** 3)

    @property
    def avail_gb(self) -> float:
        return self.avail_bytes / (1024.0 ** 3)


class _WindowsMemoryStatusEx(ctypes.Structure):
    """Win32 MEMORYSTATUSEX layout for ctypes fallback."""

    _fields_ = [
        ("dwLength", ctypes.c_uint32),
        ("dwMemoryLoad", ctypes.c_uint32),
        ("ullTotalPhys", ctypes.c_uint64),
        ("ullAvailPhys", ctypes.c_uint64),
        ("ullTotalPageFile", ctypes.c_uint64),
        ("ullAvailPageFile", ctypes.c_uint64),
        ("ullTotalVirtual", ctypes.c_uint64),
        ("ullAvailVirtual", ctypes.c_uint64),
        ("sullAvailExtendedVirtual", ctypes.c_uint64),
    ]


def _query_os_memory() -> tuple[int, int, str]:
    """Query total and available physical RAM from psutil or OS fallback."""
    try:
        import psutil  # type: ignore[import-untyped]

        vm = psutil.virtual_memory()
        return int(vm.total), int(vm.available), "psutil"
    except Exception:
        pass

    if sys.platform == "win32":
        status = _WindowsMemoryStatusEx()
        status.dwLength = ctypes.sizeof(_WindowsMemoryStatusEx)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):  # type: ignore[attr-defined]
            return int(status.ullTotalPhys), int(status.ullAvailPhys), "win32_global_memory_status_ex"

    meminfo_path = "/proc/meminfo"
    if os.path.exists(meminfo_path):
        total_kb = 0
        avail_kb = 0
        free_kb = 0
        buffers_kb = 0
        cached_kb = 0
        with open(meminfo_path, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2:
                    continue
                key = parts[0].rstrip(":")
                try:
                    val = int(parts[1])
                except ValueError:
                    continue
                if key == "MemTotal":
                    total_kb = val
                elif key == "MemAvailable":
                    avail_kb = val
                elif key == "MemFree":
                    free_kb = val
                elif key == "Buffers":
                    buffers_kb = val
                elif key == "Cached":
                    cached_kb = val
        if total_kb > 0:
            eff_avail_kb = avail_kb if avail_kb > 0 else (free_kb + buffers_kb + cached_kb)
            return int(total_kb * 1024), int(eff_avail_kb * 1024), "proc_meminfo"

    raise RuntimeError("Unable to query host physical memory via psutil, Win32 API, or /proc/meminfo")


def _env_avail_ram_bytes() -> int | None:
    """Return override from ``WHT_SOLVER_AVAIL_RAM_GB`` if set and positive."""
    raw = os.environ.get("WHT_SOLVER_AVAIL_RAM_GB")
    if raw is None or not raw.strip():
        return None
    val_gb = float(raw)
    if val_gb <= 0.0:
        raise ValueError(f"WHT_SOLVER_AVAIL_RAM_GB must be > 0, got {raw!r}")
    return int(val_gb * (1024.0 ** 3))


def get_total_ram_bytes() -> int:
    """Return total host physical RAM in bytes."""
    total, _, _ = _query_os_memory()
    return total


def get_available_ram_bytes() -> int:
    """Return available physical RAM in bytes (honoring ``WHT_SOLVER_AVAIL_RAM_GB``)."""
    override = _env_avail_ram_bytes()
    if override is not None:
        return override
    _, avail, _ = _query_os_memory()
    return avail


def probe_system_memory() -> SystemMemoryInfo:
    """Collect a full snapshot of host total and available memory."""
    total, avail, source = _query_os_memory()
    override = _env_avail_ram_bytes()
    if override is not None:
        avail = override
        source = f"{source}+WHT_SOLVER_AVAIL_RAM_GB"
    used_pct = 0.0 if total <= 0 else max(0.0, min(100.0, 100.0 * (1.0 - avail / float(total))))
    return SystemMemoryInfo(
        total_bytes=int(total),
        available_bytes=int(avail),
        used_percent=float(used_pct),
        source=source,
    )


def est_pardiso_bytes(n_dof: int, nnz: int, degree: int = 2) -> int:
    """Estimate peak sparse direct LU/LDLT factor and workspace memory in bytes.

    Uses the empirical 3D immersed-IGA fill-in scaling:
    ``nnz_LU = nnz * max(8.0, (degree + 1) ** 2 * max(n_dof, 1) ** 0.25)``
    together with 8-byte float64 values, 4-byte int32 column indices, CSR/DOF
    metadata, and a 1.35x METIS reordering + OOC workspace factor.
    """
    if n_dof <= 0 or nnz <= 0:
        return 0
    deg = max(int(degree), 1)
    fill_factor = max(8.0, float((deg + 1) ** 2) * (float(max(n_dof, 1)) ** 0.25))
    nnz_lu = float(nnz) * fill_factor
    # 8 B float64 factor entries + 4 B int32 indices = 12 B per LU non-zero.
    factor_bytes = nnz_lu * 12.0
    csr_and_dof_bytes = float(nnz) * 12.0 + float(n_dof) * 64.0
    return int(math.ceil((factor_bytes + csr_and_dof_bytes) * 1.35))


def est_amgcl_bytes(n_dof: int, nnz: int) -> int:
    """Estimate peak AMGCL + Krylov hierarchy memory in bytes."""
    if n_dof <= 0 or nnz <= 0:
        return 0
    csr_bytes = float(nnz) * 12.0 + float(max(n_dof, 1)) * 8.0
    # Matrix + aggregation coarsened levels (~1.8x CSR) + ILU0 smoother + Krylov vectors.
    return int(math.ceil(3.2 * csr_bytes + float(n_dof) * 256.0))


def estimate_solver_memory(n_dof: int, nnz: int, degree: int = 2) -> SolverMemoryEstimate:
    """Estimate memory requirements for both Pardiso and AMGCL solvers."""
    return SolverMemoryEstimate(
        n_dof=int(n_dof),
        nnz=int(nnz),
        degree=int(degree),
        pardiso_bytes=est_pardiso_bytes(n_dof, nnz, degree),
        amgcl_bytes=est_amgcl_bytes(n_dof, nnz),
    )


def _is_pypardiso_importable() -> bool:
    """Check whether ``pypardiso`` can be discovered without forcing heavy imports."""
    return importlib.util.find_spec("pypardiso") is not None


def select_linear_solver(
    solver_spec: str,
    n_dof: int,
    nnz: int,
    degree: int = 2,
    memory_budget_ratio: float = 0.70,
    avail_bytes: int | None = None,
    dof_threshold: int = 150_000,
    pypardiso_available: bool | None = None,
) -> SolverSelectionDecision:
    """Choose a linear solver respecting user preference and RAM budget.

    Parameters
    ----------
    solver_spec:
        Requested solver (e.g., ``"auto"``, ``"pardiso"``, ``"amgcl"``,
        ``"cg"``, ``"splu"``, or Kratos direct names ``"pardiso_ldlt"``,
        ``"pardiso_lu"``, ``"pardiso_llt"``, ``"amgcl_cg"``, ``"skyline_lu"``).
    n_dof:
        Number of active degrees of freedom in the linear system.
    nnz:
        Number of non-zeros in the assembled sparse system matrix.
    degree:
        B-spline polynomial degree used for 3D fill-in estimation.
    memory_budget_ratio:
        Fraction of available RAM allowed for sparse direct factorization.
    avail_bytes:
        Optional override for available host RAM in bytes. When ``None``,
        host memory (or ``WHT_SOLVER_AVAIL_RAM_GB``) is queried.
    dof_threshold:
        Maximum DOF count above which ``"auto"`` prefers iterative ``"amgcl"``.
    pypardiso_available:
        Optional test override indicating whether ``pypardiso`` is importable.
    """
    if not 0.0 < float(memory_budget_ratio) <= 1.0:
        raise ValueError(f"memory_budget_ratio must be in (0, 1], got {memory_budget_ratio}")

    total_ram = get_total_ram_bytes()
    eff_avail = int(avail_bytes) if avail_bytes is not None else get_available_ram_bytes()
    budget_bytes = int(math.floor(float(memory_budget_ratio) * float(max(eff_avail, 0))))
    est_p = est_pardiso_bytes(n_dof, nnz, degree)
    est_a = est_amgcl_bytes(n_dof, nnz)

    req = str(solver_spec).strip()
    if req != "auto":
        return SolverSelectionDecision(
            requested_solver=req,
            selected_solver=req,
            n_dof=int(n_dof),
            nnz=int(nnz),
            degree=int(degree),
            est_pardiso_bytes=est_p,
            est_amgcl_bytes=est_a,
            avail_bytes=eff_avail,
            total_bytes=total_ram,
            memory_budget_ratio=float(memory_budget_ratio),
            budget_bytes=budget_bytes,
            dof_threshold=int(dof_threshold),
            switched=False,
            switched_reason=None,
        )

    has_pypardiso = _is_pypardiso_importable() if pypardiso_available is None else bool(pypardiso_available)
    reasons: list[str] = []
    if est_p > budget_bytes:
        reasons.append(
            f"est_pardiso={est_p / (1024.0 ** 3):.3f} GB > "
            f"{memory_budget_ratio * 100.0:.0f}% of avail RAM "
            f"({eff_avail / (1024.0 ** 3):.3f} GB, budget={budget_bytes / (1024.0 ** 3):.3f} GB)"
        )
    if int(n_dof) > int(dof_threshold):
        reasons.append(f"n_dof={int(n_dof)} > dof_threshold={int(dof_threshold)}")

    if reasons:
        reason_msg = (
            "[HWT-IGA SOLVER][AUTO-SOLVER] Switched direct->iterative: "
            + "; ".join(reasons)
            + f". Selected 'amgcl' (est_amgcl={est_a / (1024.0 ** 3):.3f} GB)."
        )
        return SolverSelectionDecision(
            requested_solver="auto",
            selected_solver="amgcl",
            n_dof=int(n_dof),
            nnz=int(nnz),
            degree=int(degree),
            est_pardiso_bytes=est_p,
            est_amgcl_bytes=est_a,
            avail_bytes=eff_avail,
            total_bytes=total_ram,
            memory_budget_ratio=float(memory_budget_ratio),
            budget_bytes=budget_bytes,
            dof_threshold=int(dof_threshold),
            switched=True,
            switched_reason=reason_msg,
        )

    selected = "pardiso" if has_pypardiso else "splu"
    return SolverSelectionDecision(
        requested_solver="auto",
        selected_solver=selected,
        n_dof=int(n_dof),
        nnz=int(nnz),
        degree=int(degree),
        est_pardiso_bytes=est_p,
        est_amgcl_bytes=est_a,
        avail_bytes=eff_avail,
        total_bytes=total_ram,
        memory_budget_ratio=float(memory_budget_ratio),
        budget_bytes=budget_bytes,
        dof_threshold=int(dof_threshold),
        switched=False,
        switched_reason=None,
    )
