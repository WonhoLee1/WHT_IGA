"""P0 c!=0 guard on real Kratos DOFs (REVIEW_P0_P1 §1 #1, v2.1 §5.2.3).

The pure-Python guard is covered in ``tests/unit/test_constraints_guard.py``.
This Kratos-marked test proves the same ``ValueError`` fires when real
``Dof`` objects are passed (guard runs before any Kratos call), and that
zero/``None`` constants still create the constraint.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie

pytestmark = pytest.mark.kratos


def _slave_master_dofs(KM):  # type: ignore[no-untyped-def]
    """One slave + one master X-DOF on two nodes of a fresh ModelPart."""
    model = KM.Model()
    mp = model.CreateModelPart("P0Guard")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
    n2 = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
    n1.AddDof(KM.DISPLACEMENT_X)
    n2.AddDof(KM.DISPLACEMENT_X)
    return mp, [n1.GetDof(KM.DISPLACEMENT_X)], [n2.GetDof(KM.DISPLACEMENT_X)]


def test_nonzero_constant_rejected_on_kratos_dofs() -> None:
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mp, slave, master = _slave_master_dofs(KM)
    with pytest.raises(ValueError):
        add_linear_tie(KM, mp, 1, slave, master, np.array([[2.0]]), np.array([1e-6]))
    print("\nP0 guard: nonzero constant raised ValueError on Kratos DOFs")
    assert mp.NumberOfMasterSlaveConstraints() == 0


@pytest.mark.parametrize("constant", [None, np.zeros(1)])
def test_zero_and_none_constant_accepted_on_kratos_dofs(constant) -> None:  # type: ignore[no-untyped-def]
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    mp, slave, master = _slave_master_dofs(KM)
    add_linear_tie(KM, mp, 1, slave, master, np.array([[2.0]]), constant)
    print(f"\nP0 guard: constant={constant} accepted, constraint created")
    assert mp.NumberOfMasterSlaveConstraints() == 1
