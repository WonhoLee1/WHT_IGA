"""Master-slave constraint helpers for the R1 ``u_geom = T u_master`` route.

Kratos pitfall K11: a master DOF may never be the slave of another
constraint, so every composition (THB o extension o interpolation) is
collapsed into a single relation row here in Python (v2 §5.2.3 rule 5).

Binding order note: ``CreateNewMasterSlaveConstraint`` takes the *master*
DOF list before the *slave* DOF list
(``add_model_part_to_python.cpp`` ``CreateNewMasterSlaveConstraint1``),
while the ``LinearMasterSlaveConstraint`` constructor takes them in the
same master-first order. Both enforce ``slave = R * master + c``.

Dirichlet rule (spike S05 finding): a master DOF that carries a prescribed
(Dirichlet) value must stay *in* the relation with ``c = 0`` — Kratos folds
the prescribed value itself. Pre-folding it into ``c`` double-counts the
constant (observed: slave = R u + 2c).

Ownership: Track A file (see implementation-request_v1.md parallel tracks).
"""

from __future__ import annotations

from typing import Any

import numpy as np


def _as_dof_lists(KM: Any, slave_dofs: list[Any], master_dofs: list[Any]) -> tuple[list[Any], list[Any]]:
    return list(slave_dofs), list(master_dofs)


def add_linear_tie(
    KM: Any,
    model_part: Any,
    constraint_id: int,
    slave_dofs: list[Any],
    master_dofs: list[Any],
    relation: np.ndarray,
    constant: np.ndarray | None = None,
) -> Any:
    """Add one ``LinearMasterSlaveConstraint``: ``slave = R @ master + c``.

    ``relation`` is (n_slave, n_master); ``constant`` defaults to zeros.

    The constant must be exactly zero (c = 0 rule, v2.1 §5.2.3, K-pitfall:
    the block builder re-adds ``c`` on every build). Express offsets through a
    fixed master DOF instead. A nonzero entry raises ``ValueError``.
    """
    n_s, n_m = len(slave_dofs), len(master_dofs)
    relation = np.asarray(relation, dtype=float).reshape(n_s, n_m)
    if constant is None:
        constant = np.zeros(n_s)
    constant = np.asarray(constant, dtype=float).reshape(-1)
    if constant.shape != (n_s,):
        raise ValueError(f"add_linear_tie: constant has shape {constant.shape}, expected ({n_s},)")
    if np.any(constant != 0.0):
        raise ValueError(
            "add_linear_tie: nonzero constraint constant rejected (c = 0 rule); "
            "express the offset through a fixed master DOF"
        )
    rmat = KM.Matrix(n_s, n_m)
    for i in range(n_s):
        for j in range(n_m):
            rmat[i, j] = float(relation[i, j])
    cvec = KM.Vector(n_s)
    for i in range(n_s):
        cvec[i] = float(constant[i])
    s_dofs, m_dofs = _as_dof_lists(KM, slave_dofs, master_dofs)
    return model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", constraint_id, m_dofs, s_dofs, rmat, cvec
    )


def add_identity_tie(
    KM: Any,
    model_part: Any,
    constraint_id: int,
    master_node: Any,
    slave_node: Any,
    variables: tuple[Any, ...] | None = None,
) -> Any:
    """Tie one slave node to one master node component-wise (weight 1, c 0)."""
    if variables is None:
        variables = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    slave_dofs = [slave_node.GetDof(v) for v in variables]
    master_dofs = [master_node.GetDof(v) for v in variables]
    return add_linear_tie(
        KM, model_part, constraint_id, slave_dofs, master_dofs, np.eye(len(variables))
    )


def nodal_displacements(KM: Any, node: Any) -> np.ndarray:
    """Current solution-step DISPLACEMENT of a node as a 3-vector."""
    return np.array(
        [
            node.GetSolutionStepValue(KM.DISPLACEMENT_X),
            node.GetSolutionStepValue(KM.DISPLACEMENT_Y),
            node.GetSolutionStepValue(KM.DISPLACEMENT_Z),
        ]
    )
