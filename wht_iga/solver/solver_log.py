"""Commercial-style 6-stage ASCII/UTF-8 solver logger for WHT-IGA."""

from __future__ import annotations

import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, TextIO

if TYPE_CHECKING:
    from wht_iga.solver.memory import SolverSelectionDecision


SolverLogLevel = Literal[
    "summary",
    "detailed",
    "debug",
    "minimal",
    "standard",
    "verbose",
]

_LEVEL_RANK = {
    "summary": 1,
    "minimal": 1,
    "detailed": 2,
    "standard": 2,
    "debug": 3,
    "verbose": 3,
}

_CANONICAL_LEVEL = {
    "summary": "summary",
    "minimal": "summary",
    "detailed": "detailed",
    "standard": "detailed",
    "debug": "debug",
    "verbose": "debug",
}


def resolve_solver_log_enabled(enable_solver_log: bool = True) -> bool:
    """Resolve effective logging switch considering ``WHT_SOLVER_LOG`` env override."""
    if not enable_solver_log:
        return False
    env = os.environ.get("WHT_SOLVER_LOG")
    if env is not None:
        val = env.strip().lower()
        if val in ("0", "false", "no", "off"):
            return False
        if val in ("1", "true", "yes", "on"):
            return True
    return True


def _fmt_seq(val: Any) -> str:
    """Format vector/tuple quantities compactly for log tables."""
    if val is None:
        return "N/A"
    if hasattr(val, "__iter__") and not isinstance(val, (str, bytes)):
        return "[" + ", ".join(f"{float(x):.6e}" for x in val) + "]"
    return str(val)


class SolverLogger:
    """6-stage Abaqus/ANSYS-style solver diagnostic logger.

    Parameters
    ----------
    enable_solver_log:
        Master enable/disable flag. When ``False`` (or ``WHT_SOLVER_LOG=0``),
        no file is created and no stdout lines are emitted.
    solver_log_file:
        Optional output path for the solver log (written with ``encoding="utf-8"``).
    solver_log_stdout:
        If ``True``, mirror log lines to stdout (or ``stdout_stream``).
    solver_log_level:
        Verbosity tier: ``"summary"`` (or ``"minimal"``), ``"detailed"``
        (or ``"standard"``), or ``"debug"`` (or ``"verbose"``).
    """

    BANNER_TITLE = "HWT-IGA SOLVER — Immersed Isogeometric Analysis Engine"

    def __init__(
        self,
        enable_solver_log: bool = True,
        solver_log_file: str | Path | None = None,
        solver_log_stdout: bool = False,
        solver_log_level: SolverLogLevel = "detailed",
        *,
        version: str = "2.1",
        stdout_stream: TextIO | None = None,
    ) -> None:
        lvl_key = str(solver_log_level).strip().lower()
        if lvl_key not in _LEVEL_RANK:
            raise ValueError(
                f"solver_log_level must be one of {sorted(_LEVEL_RANK)}, got {solver_log_level!r}"
            )
        self.enabled = resolve_solver_log_enabled(bool(enable_solver_log))
        self.solver_log_file = Path(solver_log_file) if solver_log_file is not None else None
        self.solver_log_stdout = bool(solver_log_stdout)
        self.solver_log_level = _CANONICAL_LEVEL[lvl_key]
        self._rank = _LEVEL_RANK[lvl_key]
        self.version = version
        self._stdout = stdout_stream if stdout_stream is not None else sys.stdout
        self._fh: TextIO | None = None
        self.lines: list[str] = []
        self._stage5_header_written = False
        self._banner_written = False

    def __enter__(self) -> SolverLogger:
        self.open()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()

    def open(self) -> None:
        """Open the backing log file in UTF-8 mode if logging is enabled."""
        if not self.enabled or self.solver_log_file is None or self._fh is not None:
            return
        self.solver_log_file.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.solver_log_file, "w", encoding="utf-8", newline="\n")

    def close(self) -> None:
        """Flush and close any open backing file handle."""
        if self._fh is not None:
            try:
                self._fh.flush()
            finally:
                self._fh.close()
                self._fh = None

    def _emit(self, line: str, min_rank: int = 1) -> None:
        """Emit one formatted line to memory buffer, file, and/or stdout."""
        if not self.enabled or self._rank < min_rank:
            return
        if self.solver_log_file is not None and self._fh is None:
            self.open()
        self.lines.append(line)
        if self._fh is not None:
            self._fh.write(line + "\n")
            self._fh.flush()
        if self.solver_log_stdout:
            self._stdout.write(line + "\n")
            self._stdout.flush()

    def log_line(self, line: str, *, level: str = "detailed") -> None:
        """Write an arbitrary diagnostic line at the given verbosity level."""
        rank = _LEVEL_RANK.get(str(level).lower(), 2)
        self._emit(line, min_rank=rank)

    def log_warning(self, code: str, message: str) -> None:
        """Emit an explicit solver warning tag."""
        self._emit(f"  [WARNING:{code}] {message}", min_rank=1)

    def write_banner(self, job_name: str = "") -> None:
        """Emit the top solver banner once."""
        if not self.enabled or self._banner_written:
            return
        self._banner_written = True
        rule = "=" * 84
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        self._emit(rule, min_rank=1)
        self._emit(f"  {self.BANNER_TITLE} (v{self.version})", min_rank=1)
        if job_name:
            self._emit(f"  Job: {job_name}    Started: {ts}    Level: {self.solver_log_level}", min_rank=1)
        else:
            self._emit(f"  Started: {ts}    Level: {self.solver_log_level}", min_rank=1)
        self._emit(rule, min_rank=1)

    def log_stage_1_environment(
        self,
        job_name: str = "",
        cpu_count: int | None = None,
        total_ram_gb: float = 0.0,
        avail_ram_gb: float = 0.0,
        mkl_threads: int | str | None = None,
        python_version: str | None = None,
        kratos_version: str = "",
        extra: dict[str, Any] | None = None,
    ) -> None:
        """Emit Stage [1/6] JOB & SYSTEM ENVIRONMENT."""
        if not self.enabled:
            return
        self.write_banner(job_name)
        cpus = os.cpu_count() if cpu_count is None else cpu_count
        pyver = python_version or platform.python_version()
        host = platform.node() or "unknown-host"
        os_str = f"{platform.system()} {platform.release()} ({platform.machine()})"
        self._emit("", min_rank=1)
        self._emit("[1/6] JOB & SYSTEM ENVIRONMENT", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(f"  Job Name         : {job_name or '(unnamed)'}", min_rank=1)
        self._emit(f"  Host / OS        : {host} / {os_str}", min_rank=2)
        self._emit(f"  Python           : {pyver}", min_rank=2)
        if kratos_version:
            self._emit(f"  Kratos Kernel    : {kratos_version}", min_rank=2)
        self._emit(
            f"  CPU / MKL        : {cpus} logical cores | mkl_threads={mkl_threads if mkl_threads is not None else 'default'}",
            min_rank=2,
        )
        self._emit(
            f"  Host RAM         : total={total_ram_gb:.2f} GB | available={avail_ram_gb:.2f} GB",
            min_rank=1,
        )
        if extra:
            for k, v in extra.items():
                self._emit(f"  {k:<17}: {v}", min_rank=3)

    def log_stage_2_model(
        self,
        n_cells: tuple[int, int, int],
        degree: int,
        continuity: int | None,
        n_ctrl_pts: int,
        n_total_dofs: int,
        n_active_dofs: int,
        n_in_cells: int,
        n_cut_cells: int,
        n_out_cells: int,
        out_mode: str,
        n_quad_pts: int = 0,
        nnz: int = 0,
    ) -> None:
        """Emit Stage [2/6] MODEL & IMMERSED GRID SUMMARY."""
        if not self.enabled:
            return
        self._emit("", min_rank=1)
        self._emit("[2/6] MODEL & IMMERSED GRID SUMMARY", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(
            f"  Background Grid  : n_cells={tuple(n_cells)} | degree=p{degree} | continuity={continuity}",
            min_rank=1,
        )
        self._emit(
            f"  Immersed Cells   : IN={n_in_cells} | CUT={n_cut_cells} | OUT={n_out_cells} | out_mode={out_mode}",
            min_rank=1,
        )
        self._emit(
            f"  Control Points   : {n_ctrl_pts} | Total DOFs={n_total_dofs} | Active DOFs={n_active_dofs}",
            min_rank=1,
        )
        self._emit(
            f"  Quadrature / NNZ : n_qp={n_quad_pts} | est_nnz={nnz}",
            min_rank=2,
        )

    def log_stage_3_constraints(
        self,
        n_dirichlet_bcs: int = 0,
        n_fixed_faces: int = 0,
        n_inactive_fixed: int = 0,
        n_ties: int = 0,
        n_slave_dofs: int = 0,
        n_master_dofs: int = 0,
        beta: float | None = None,
        c_pen: float | None = None,
        theta_thb: float | None = None,
        extra_notes: list[str] | None = None,
    ) -> None:
        """Emit Stage [3/6] CONSTRAINT & INTERACTION SUMMARY."""
        if not self.enabled:
            return
        self._emit("", min_rank=1)
        self._emit("[3/6] CONSTRAINT & INTERACTION SUMMARY", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(
            f"  Dirichlet / Fix  : bcs={n_dirichlet_bcs} | fixed_faces={n_fixed_faces} | inactive_fixed={n_inactive_fixed}",
            min_rank=1,
        )
        beta_str = f"{beta:.6e}" if beta is not None else "N/A"
        cpen_str = f"{c_pen:.6e}" if c_pen is not None else "N/A"
        theta_str = f"{theta_thb:.6e}" if theta_thb is not None else "N/A"
        self._emit(
            f"  Penalty Params   : beta={beta_str} | C_pen={cpen_str}",
            min_rank=2,
        )
        self._emit(
            f"  Extension Ties   : n_ties={n_ties} | slave_dofs={n_slave_dofs} | master_dofs={n_master_dofs} | theta={theta_str}",
            min_rank=1 if n_ties > 0 else 2,
        )
        if extra_notes:
            for note in extra_notes:
                self._emit(f"  Note             : {note}", min_rank=2)

    def log_stage_4_solver_selection(
        self,
        decision: SolverSelectionDecision | None = None,
        *,
        requested_solver: str = "auto",
        selected_solver: str = "",
        est_pardiso_gb: float = 0.0,
        est_amgcl_gb: float = 0.0,
        avail_ram_gb: float = 0.0,
        budget_ratio: float = 0.70,
        switched: bool = False,
        switched_reason: str | None = None,
        solver_params: dict[str, Any] | None = None,
    ) -> None:
        """Emit Stage [4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION."""
        if not self.enabled:
            return
        if decision is not None:
            requested_solver = decision.requested_solver
            selected_solver = decision.selected_solver
            est_pardiso_gb = decision.est_pardiso_gb
            est_amgcl_gb = decision.est_amgcl_gb
            avail_ram_gb = decision.avail_gb
            budget_ratio = decision.memory_budget_ratio
            switched = decision.switched
            switched_reason = decision.switched_reason

        self._emit("", min_rank=1)
        self._emit("[4/6] MEMORY ESTIMATION & LINEAR SOLVER SELECTION", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(
            f"  Solver Request   : requested={requested_solver} -> selected={selected_solver}",
            min_rank=1,
        )
        self._emit(
            f"  Memory Estimate  : est_pardiso={est_pardiso_gb:.3f} GB | est_amgcl={est_amgcl_gb:.3f} GB | "
            f"avail={avail_ram_gb:.3f} GB (budget={budget_ratio * 100.0:.0f}%)",
            min_rank=1,
        )
        if switched and switched_reason:
            self._emit(f"  Auto-Switch      : {switched_reason}", min_rank=1)
        if solver_params:
            self._emit(f"  Solver Settings  : {solver_params}", min_rank=3)

    def log_stage_5_header(self) -> None:
        """Emit Stage [5/6] STEP & ITERATION PROGRESS TABLE header."""
        if not self.enabled or self._stage5_header_written:
            return
        self._stage5_header_written = True
        self._emit("", min_rank=1)
        self._emit("[5/6] STEP & ITERATION PROGRESS TABLE", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(
            "Step | Inc | Att |      Time |        dt | LinIt |      RelRes |  EnergyNorm | Wall(s)",
            min_rank=1,
        )
        self._emit(
            "-----+-----+-----+-----------+-----------+-------+-------------+-------------+--------",
            min_rank=1,
        )

    def log_stage_5_row(
        self,
        step: int = 1,
        inc: int = 1,
        att: int = 1,
        time: float = 1.0,
        dt: float = 1.0,
        lin_it: int = 1,
        rel_res: float = 0.0,
        energy_norm: float = 0.0,
        wall_s: float = 0.0,
    ) -> None:
        """Emit one increment/iteration row in Stage [5/6]."""
        if not self.enabled:
            return
        if not self._stage5_header_written:
            self.log_stage_5_header()
        row = (
            f"{step:4d} | {inc:3d} | {att:3d} | {time:9.4f} | {dt:9.4e} | "
            f"{lin_it:5d} | {rel_res:11.3e} | {energy_norm:11.3e} | {wall_s:7.3f}"
        )
        self._emit(row, min_rank=1)

    def log_stage_6_summary(
        self,
        converged: bool = True,
        max_disp: float = 0.0,
        total_force: Any = None,
        reaction_residual: float | None = None,
        setup_s: float = 0.0,
        solve_s: float = 0.0,
        post_s: float = 0.0,
        total_s: float | None = None,
        peak_ram_gb: float | None = None,
        warnings: list[str] | None = None,
        error_message: str | None = None,
    ) -> None:
        """Emit Stage [6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY."""
        if not self.enabled:
            return
        eff_total = (setup_s + solve_s + post_s) if total_s is None else total_s
        self._emit("", min_rank=1)
        self._emit("[6/6] SOLUTION QUALITY, REACTION BALANCE & TIMING SUMMARY", min_rank=1)
        self._emit("-" * 84, min_rank=1)
        self._emit(
            f"  Solution Norms   : max|U|={max_disp:.6e} | total_force={_fmt_seq(total_force)}",
            min_rank=1,
        )
        if reaction_residual is not None:
            self._emit(f"  Reaction Balance : ||R_sum + F_ext|| = {reaction_residual:.6e}", min_rank=1)
        self._emit(
            f"  Timings (s)      : setup={setup_s:.4f} | solve={solve_s:.4f} | post={post_s:.4f} | total={eff_total:.4f}",
            min_rank=1,
        )
        if peak_ram_gb is not None:
            self._emit(f"  Peak RSS RAM     : {peak_ram_gb:.3f} GB", min_rank=2)
        if warnings:
            self._emit(f"  Warnings ({len(warnings)})     :", min_rank=1)
            for w in warnings:
                self._emit(f"    - {w}", min_rank=1)
        if error_message:
            self._emit(f"  Error Detail     : {error_message}", min_rank=1)

        self._emit("-" * 84, min_rank=1)
        if converged:
            self._emit("  THE ANALYSIS HAS COMPLETED SUCCESSFULLY", min_rank=1)
        else:
            self._emit("  THE ANALYSIS DID NOT CONVERGE", min_rank=1)
        self._emit("=" * 84, min_rank=1)

    def get_text(self) -> str:
        """Return all logged lines with trailing newline matching the log file (empty if disabled)."""
        if not self.enabled or not self.lines:
            return ""
        return "\n".join(self.lines) + "\n"
