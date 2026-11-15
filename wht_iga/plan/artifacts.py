"""Artifact file layout of a study (plan v2.1 §4: ``.wiga.sqlite``, ``_iga_results.sqlite``,
``.vtkhdf``, ``.spec.json``, ``_probe_data.json``, ``_solver.log``).

A single-step study uses the job name as file prefix (``roark.wiga.sqlite``);
a multi-step study prefixes per step (``roark_load.wiga.sqlite``). The file
names of the solver pipeline (``wht_iga.solver.pipeline._write_outputs``) are
reproduced exactly so the session can find what ``run_linear`` wrote.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class StepArtifacts:
    """Files belonging to one step."""

    output_dir: Path
    prefix: str

    @property
    def mesh_db(self) -> Path:
        """Mesh / quadrature database written by the Kratos pipeline."""
        return self.output_dir / f"{self.prefix}.wiga.sqlite"

    @property
    def results_db(self) -> Path:
        """Results database (frames, control-point and QP fields)."""
        return self.output_dir / f"{self.prefix}_iga_results.sqlite"

    @property
    def vtu(self) -> Path:
        """Single-shot VTU proxy written by the Kratos pipeline."""
        return self.output_dir / f"{self.prefix}.vtu"

    @property
    def vtkhdf(self) -> Path:
        """VTKHDF proxy fields written by the session (all backends)."""
        return self.output_dir / f"{self.prefix}.vtkhdf"

    @property
    def probe_json(self) -> Path:
        """Probe values (JSON)."""
        return self.output_dir / f"{self.prefix}_probe_data.json"

    @property
    def solver_log(self) -> Path:
        """HWT-IGA SOLVER 6-stage log."""
        return self.output_dir / f"{self.prefix}_solver.log"

    def all(self) -> tuple[Path, ...]:
        """Every file this step may produce."""
        return (self.mesh_db, self.results_db, self.vtu, self.vtkhdf, self.probe_json, self.solver_log)


@dataclass(frozen=True)
class ArtifactPaths:
    """All artifact paths of a study (resolved before any solver runs)."""

    output_dir: Path
    job: str
    step_names: tuple[str, ...]

    @property
    def spec_json(self) -> Path:
        """Input-reproduction copy of the ``StudySpec``."""
        return self.output_dir / f"{self.job}.spec.json"

    def step(self, name: str) -> StepArtifacts:
        """Artifacts of the named step."""
        if name not in self.step_names:
            raise KeyError(f"no step named {name!r}; steps: {list(self.step_names)}")
        prefix = self.job if len(self.step_names) == 1 else f"{self.job}_{name}"
        return StepArtifacts(self.output_dir, prefix)

    def all(self) -> tuple[Path, ...]:
        """Every file the study may produce."""
        out: list[Path] = [self.spec_json]
        for s in self.step_names:
            out.extend(self.step(s).all())
        return tuple(out)

    def existing(self) -> tuple[Path, ...]:
        """Artifacts that exist on disk now."""
        return tuple(p for p in self.all() if p.exists())
