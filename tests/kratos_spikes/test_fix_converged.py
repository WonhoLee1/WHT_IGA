# -*- coding: utf-8 -*-
"""F3: solver convergence is recorded, not fabricated.

run_linear must report the real Kratos IsConverged / NL_ITERATION_NUMBER
and write them into the results DB frames table (default frames mean
"not recorded", never "converged").
"""

from __future__ import annotations

import pytest

from wht_iga.solver.validate import DirichletBC, JobSpec, PressureLoad

pytestmark = pytest.mark.kratos


def _untrimmed() -> JobSpec:
    return JobSpec(
        name="f3_conv", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
        n_cells=(2, 2, 2), degree=2, E=1.0, nu=0.0,
        bcs=[DirichletBC(face=(0, 0), components=(True, False, False))],
        pressures=[PressureLoad(face=(0, 1), magnitude=1.0)],
        output_dir=None,
    )


def test_converged_and_niter_recorded() -> None:
    from wht_iga.solver.pipeline import run_linear

    res = run_linear(_untrimmed())
    assert res.converged is True
    assert isinstance(res.n_iter, int) and res.n_iter >= 1


def test_converged_frame_written_to_db(tmp_path) -> None:  # type: ignore[no-untyped-def]
    import os

    from wht_iga.solver import db as meshdb
    from wht_iga.solver.pipeline import run_linear

    outdir = str(tmp_path)
    import dataclasses

    spec = dataclasses.replace(_untrimmed(), output_dir=outdir)
    res = run_linear(spec)
    back = meshdb.read_results_db(os.path.join(outdir, "f3_conv_iga_results.sqlite"))
    assert back["frames"][0]["converged"] == 1
    assert back["frames"][0]["n_iter"] == int(res.n_iter) >= 1


def test_nonconverged_raises_validation_error(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    """Non-convergence must raise ValidationError(SOLVER_NOT_CONVERGED), not NameError.

    Regression (2026-11-15): the branch referenced ValidationIssue without
    importing it, so a real non-convergence surfaced as NameError.
    A 1-iteration CG solver forces the Newton loop to exhaust its iterations.
    """
    from wht_iga.solver import pipeline as PL
    from wht_iga.solver.validate import ValidationError

    monkeypatch.setattr(
        PL, "_make_linear_solver",
        lambda KM, name, thr=None, amgcl_params=None: KM.CGSolver(1e-30, 1),
    )
    spec = _untrimmed()
    spec.bcs = [DirichletBC(face=(0, 0))]  # well-posed: CG alone is the failure cause
    with pytest.warns(RuntimeWarning), pytest.raises(ValidationError) as ei:
        PL.run_linear(spec)
    codes = [i.code for i in ei.value.issues]
    assert "SOLVER_NOT_CONVERGED" in codes
    # matrix-only callers (V-COND) can still read the assembled K from the error
    from wht_iga.kratos_backend.elements import system_matrix_csr

    K = system_matrix_csr(ei.value.strategy).tocsr()
    assert K.shape[0] > 0 and K.nnz > 0
    assert abs(K - K.T).max() <= 1e-10 * abs(K).max()
