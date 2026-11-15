"""T-P2.3: reduced numbering + batch LinearMasterSlaveConstraint specs.

Constraint model (v2 §5.2.3): one relation per slave control point, expanded
to 3 scalar relations (one per displacement component) with a 3×3m block-
diagonal relation matrix (m = nonzeros of that T row). Thermal problems share
the same T with 1 component (caller selects ``n_comp``).

Hard rules (own implementation — ``kratos_backend/constraints.py`` is
Phase-1 owned and NOT modified here):
* ``require_zero_constant``: any nonzero constant c is rejected with an error
  (v2 §5.2.3: c is always 0; Kratos BlockBuilder re-assembles c every Build
  and double-counts it in incremental schemes — Track A F2).
* ``validate_no_chains``: a master of one relation may never be the slave of
  another (pitfall K22); the batch builder refuses chained sets before any
  Kratos object is created.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp


def require_zero_constant(constant: np.ndarray | None, tol: float = 0.0,
                          context: str = "") -> np.ndarray:
    """Return zeros; raise if any |c_i| > tol (c≠0 refusal, STAB-owned check)."""
    if constant is None:
        return np.zeros(0)
    c = np.asarray(constant, dtype=float)
    bad = np.flatnonzero(np.abs(c) > tol)
    if len(bad):
        raise ValueError(
            f"nonzero master-slave constant rejected{': ' + context if context else ''} "
            f"(max |c| = {np.abs(c).max():.3e} at {bad[:8].tolist()}; "
            "v2 §5.2.3: c is always 0)")
    return np.zeros_like(c)


def validate_no_chains(slave_ids: np.ndarray,
                       master_lists: list[np.ndarray]) -> None:
    """Raise if any master id also appears as a slave (chained constraint)."""
    slaves = set(int(s) for s in np.asarray(slave_ids))
    offenders: list[int] = []
    for ml in master_lists:
        for m in np.asarray(ml, dtype=np.int64).ravel():
            if int(m) in slaves and int(m) not in offenders:
                offenders.append(int(m))
    if offenders:
        raise ValueError(
            f"chained constraints rejected: {len(offenders)} master(s) are "
            f"also slaves (e.g. {offenders[:8]}); pre-compose in Python "
            "(v2 §5.2.3 rule 5, pitfall K22)")


@dataclass
class ReducedNumbering:
    """Master control-point ids ↔ reduced columns; slaves listed separately."""

    master_ids: np.ndarray  # (n_m,) flat basis ids, column order
    slave_ids: np.ndarray  # (n_s,) flat basis ids
    n_full: int

    @property
    def n_master(self) -> int:
        return len(self.master_ids)

    def column_of(self, master_id: int) -> int:
        pos = int(np.searchsorted(self.master_ids, int(master_id)))
        if pos >= len(self.master_ids) or self.master_ids[pos] != master_id:
            raise KeyError(f"basis {master_id} is not a master")
        return pos


@dataclass
class ConstraintSpec:
    """One scalar LinearMasterSlaveConstraint: slave = R @ master (c = 0)."""

    slave_dof: int
    master_dofs: np.ndarray  # (m,) scalar DOFs, reference ordering 3*I+d
    row: np.ndarray  # (m,) relation weights
    constant: np.ndarray  # zeros (enforced)


def build_numbering(n_full: int, master_ids: np.ndarray,
                    slave_ids: np.ndarray) -> ReducedNumbering:
    """Validate disjointness (chain check) and build the numbering."""
    m = np.asarray(master_ids, dtype=np.int64)
    s = np.asarray(slave_ids, dtype=np.int64)
    validate_no_chains(s, [m])
    order = np.argsort(m)
    return ReducedNumbering(master_ids=m[order], slave_ids=np.sort(s),
                            n_full=int(n_full))


def make_scalar_specs(numbering: ReducedNumbering, T: sp.spmatrix,
                      slave_cp_ids: np.ndarray, n_comp: int = 3,
                      dof_stride: int = 3,
                      constant: np.ndarray | None = None) -> list[ConstraintSpec]:
    """Expand control-point relations to scalar-DOF specs (block-diagonal).

    :param T: full (n_full × n_m) CSR operator (columns follow
        ``numbering.master_ids`` order).
    :param slave_cp_ids: flat basis ids of slaves, row order in T.
    :param n_comp: components per node (3 structure, 1 thermal).
    :param constant: optional constant vector; any nonzero entry is rejected
        (c = 0 rule, v2 §5.2.3). ``None`` (default) stores exact zeros.
    """
    if not sp.issparse(T):
        raise TypeError("T must be scipy.sparse (sparse-only rule)")
    T = T.tocsr()
    c0 = require_zero_constant(
        np.zeros(1) if constant is None else constant, tol=0.0,
        context="make_scalar_specs")
    specs: list[ConstraintSpec] = []
    for g in np.asarray(slave_cp_ids, dtype=np.int64).ravel():
        r = T.getrow(int(g))
        masters = numbering.master_ids[r.indices]
        for d in range(n_comp):
            mdofs = (masters * dof_stride + d)
            specs.append(ConstraintSpec(
                slave_dof=int(g) * dof_stride + int(d),
                master_dofs=np.asarray(mdofs, dtype=np.int64),
                row=np.asarray(r.data, dtype=float),
                constant=np.array([c0[0]])))
    return specs
