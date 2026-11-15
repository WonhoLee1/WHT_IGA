"""Kratos CAUCHY_STRESS_VECTOR Voigt order, checked against Kratos itself.

The unit test test_trimmed_probe.py::test_voigt_order_matches_kratos only checks
our own convention. Here a pure-shear affine field is imposed on a 1-cell cube
(E=2, nu=0 -> mu=1, sigma_ij = gamma) and the slot that carries the stress is
read from Kratos: expected order [xx, yy, zz, xy, yz, xz] = trimmed._VOIGT.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.solver.validate import AffineField, DirichletBC, JobSpec

pytestmark = pytest.mark.kratos

G = 1e-3
CASES = {"xy": ((0, G, 0), (0, 0, 0), (0, 0, 0)),   # u = g*y e_x
         "yz": ((0, 0, 0), (0, 0, G), (0, 0, 0)),   # u = g*z e_y
         "xz": ((0, 0, G), (0, 0, 0), (0, 0, 0))}   # u = g*z e_x


@pytest.mark.parametrize("comp", list(CASES))
def test_shear_slot_matches_voigt_map(comp: str) -> None:
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.trimmed import _VOIGT

    A = CASES[comp]
    bcs = [DirichletBC(face=(a, s), field=AffineField(A=A, b=(0.0, 0.0, 0.0)))
           for a in range(3) for s in (0, 1)]
    spec = JobSpec(name=f"voigt_{comp}", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                   n_cells=(1, 1, 1), degree=1, E=2.0, nu=0.0, bcs=bcs)
    s = np.asarray(run_linear(spec).cauchy_stress[0], float)[0] / G
    expected = np.zeros(6)
    expected[_VOIGT[comp]] = 1.0
    np.testing.assert_allclose(s, expected, atol=1e-4)
