# -*- coding: utf-8 -*-
"""Kratos penalty mortar contact (unilateral, frictionless) between tied surface layers.

Recipe from S29/S29b (``docs/spikes/S29.md``). Kratos' ``PenaltyContactProcess``
cannot build its interface on IGA models (``GenerateFaces`` fails on
``QuadraturePointGeometry``). Instead each contact surface gets a *surface layer*:
plain FE nodes + ``SurfaceCondition3D4N`` quads whose nodes are tied 1:1 to
control points by ``LinearMasterSlaveConstraint``. The layers are handed to the
penalty process directly, and the solver is set up the way Kratos'
``ContactStaticMechanicalSolver`` does it, with the changes below.

Everything a run needs, and why (each item was isolated by a differential test
against Kratos' own ``penalty_frictionless_contact_test_3D``):

1. ``MoveMeshFlag = True`` (mandatory). The weighted gap and the active set
   (``MortarExplicitContributionUtilities`` -> ``ComputePenaltyFrictionlessActiveSet``)
   are computed from the *current node coordinates*, not from X0 + u. With a
   fixed mesh the weighted gap stays at its initial value (0 for touching
   layers), so the active set is decided by round-off: matching triangle
   layers kept 1 of 4 slave nodes active (force 0.87 instead of 5, the same
   tension when pulled: the S29 failure), quad layers dropped all nodes (no
   force at all). The condition's own residual uses X0 + u, so the force is
   right whenever the active set is right.
   ``SmallDisplacementElement`` on QP geometry uses the reference
   configuration, so moving the mesh does not change its stiffness (FE hexa8
   and IGA give identical reactions, 4.990020 at penalty 1e6).
2. ``MortarAndConvergenceCriteria(<displacement criterion>,
   PenaltyFrictionlessMortarConvergenceCriteria)``: the mortar criterion
   recomputes the weighted gap and updates ``ACTIVE`` after each iteration and
   reports active-set convergence. Plain Newton never updates the active set.
3. The base criterion must be displacement-based. Kratos' contact factory sets
   ``SetActualizeRHSFlag(True)`` and the block builder's ``BuildRHS`` does not
   apply the master-slave transformation, so with tied layers the residual
   still contains the layer (slave) DOF forces and a residual criterion never
   converges (quad layers: 30 iterations, although the reaction is already
   the converged 4.990020). A residual criterion with ``SetActualizeRHSFlag(False)``
   also converges (2 iterations); the displacement criterion is used here.
4. Quad layers (``SurfaceCondition3D4N``). Matching *triangle* layers with the
   same diagonal on both sides make the mortar clipping of identical triangles
   degenerate: one pair intermittently gets no weighted gap, its node drops out
   of the active set and the Newton loop chatters (no convergence in 30 its).
5. ``ProcessInfo[IS_RESTARTED] = 1`` only during ``ExecuteInitialize`` (skips
   ``NormalCheckProcess``, which needs element faces), and
   ``_interface_preprocess`` overridden to a no-op (the layers already are the
   interface conditions). The layer normals must point outward; this module
   orients the quads from a given outward vector.
6. ``model_part.Set(INTERACTION, True)`` (official default: semi-smooth active
   set update in every iteration) and ``ResidualBasedNewtonRaphsonContactStrategy``
   with ``ContactResidualBasedBlockBuilderAndSolver``. These are what the
   official solver uses; for penalty contact on this model they were not
   decisive (plain Newton + block builder gave the same numbers).

Penalty scaling: each contact condition multiplies the penalty by its own
*weighted* gap (gap x integral of N_i over that condition) and integrates the
resulting pressure again, so for a uniform gap the spring of one quad pair of
area A_e is ``penalty * A_e**2 / 4`` and the interface spring is the sum over
conditions (measured to 1e-14: 1x1 layer -> penalty/4, 2x2 layer -> penalty/16
on a unit face). Per unit area that is penalty*A_e/4, so a target c*E/h per
unit area needs penalty ~ 4*c*E/(h*A_e). The active set, in contrast, uses the
assembled nodal ``WEIGHTED_GAP``.

Not supported: friction. Kratos 10.4's penalty frictional mortar condition
weights the tangential residual by the change of the mortar operators over the
step (M_old - M, D_old - D) instead of by M, D; the friction force is then of
order mu*p_n*(slip increment/element size). Measured: 1.4e-4 tangential
force at 1e-3 tangential travel for mu*N = 1.5 (S29 probe, 3D IGA layers), and
|Rx/Ry| <= 0.03 with mu = 1 in Kratos' own 2D "stick" patch test while the
upper interface nodes move 6e-3..5e-2 more than the lower ones.

Kratos is imported only through ``wht_iga.kratos_backend.load_kratos``.
"""

from __future__ import annotations

import json
from typing import Any, Callable, Sequence

import numpy as np

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.constraints import add_linear_tie

APPS = ("StructuralMechanicsApplication", "ContactStructuralMechanicsApplication")

_LAYER_PROCESS_CLASS = None


def _kratos() -> tuple[Any, Any]:
    mods = load_kratos(APPS)
    return mods["KratosMultiphysics"], mods["ContactStructuralMechanicsApplication"]


def add_contact_variables(model_part: Any) -> None:
    """Register the nodal variables the penalty mortar contact needs.

    Call before any node is created (Kratos pitfall K13). Adds DISPLACEMENT,
    REACTION, NORMAL, NODAL_H and WEIGHTED_GAP (``AuxiliaryAddVariables`` for
    ``PenaltyContactFrictionless``).
    """
    KM, CSMA = _kratos()
    for var in (KM.DISPLACEMENT, KM.REACTION, KM.NORMAL, KM.NODAL_H, CSMA.WEIGHTED_GAP):
        model_part.AddNodalSolutionStepVariable(var)


def add_tied_quad_layer(
    model_part: Any,
    name: str,
    master_nodes: Sequence[Sequence[Any]],
    outward: Sequence[float],
    properties: Any,
    first_node_id: int,
    first_condition_id: int,
    first_constraint_id: int,
) -> tuple[Any, int]:
    """Create a contact surface layer tied 1:1 to existing nodes.

    ``master_nodes[i][j]`` is a structured (n_i x n_j) grid of Kratos nodes that
    lie on the contact surface (for p = 1 IGA faces: the face control points,
    which interpolate the surface). For each master node a layer node is created
    at the same position, every grid cell becomes one ``SurfaceCondition3D4N``
    whose geometric normal points along ``outward``, and each layer DOF is tied
    to its master DOF (``slave = master``). A 1:1 tie is exact only where the
    face interpolates its control points; p >= 2 faces need ties through the
    spline basis (not provided here).

    Returns ``(sub_model_part, next_constraint_id)``.
    """
    KM, _ = _kratos()
    grid = [list(row) for row in master_nodes]
    ni, nj = len(grid), len(grid[0])
    if ni < 2 or nj < 2 or any(len(row) != nj for row in grid):
        raise ValueError(f"add_tied_quad_layer: need a rectangular grid of at least 2x2 nodes, got {ni} rows")
    out = np.asarray(outward, dtype=float)
    smp = model_part.CreateSubModelPart(name)
    layer = [[None] * nj for _ in range(ni)]
    nid = first_node_id
    for i in range(ni):
        for j in range(nj):
            m = grid[i][j]
            node = smp.CreateNewNode(nid, m.X0, m.Y0, m.Z0)
            for var, reac in ((KM.DISPLACEMENT_X, KM.REACTION_X), (KM.DISPLACEMENT_Y, KM.REACTION_Y),
                              (KM.DISPLACEMENT_Z, KM.REACTION_Z)):
                node.AddDof(var, reac)
            layer[i][j] = node
            nid += 1
    cid = first_condition_id
    for i in range(ni - 1):
        for j in range(nj - 1):
            quad = [layer[i][j], layer[i + 1][j], layer[i + 1][j + 1], layer[i][j + 1]]
            p = np.array([[n.X0, n.Y0, n.Z0] for n in quad])
            normal = np.cross(p[1] - p[0], p[3] - p[0])
            if float(normal @ out) < 0.0:
                quad = [quad[0], quad[3], quad[2], quad[1]]
            smp.CreateNewCondition("SurfaceCondition3D4N", cid, [n.Id for n in quad], properties)
            cid += 1
    kid = first_constraint_id
    for i in range(ni):
        for j in range(nj):
            for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z):
                add_linear_tie(KM, model_part, kid, [layer[i][j].GetDof(var)], [grid[i][j].GetDof(var)],
                               np.array([[1.0]]))
                kid += 1
    return smp, kid


def _layer_process_class() -> type:
    """``PenaltyContactProcess`` that uses ready-made layer conditions (S29 bypass)."""
    global _LAYER_PROCESS_CLASS
    if _LAYER_PROCESS_CLASS is None:
        _kratos()
        import KratosMultiphysics.ContactStructuralMechanicsApplication.penalty_contact_process as pcp

        class SurfaceLayerPenaltyContactProcess(pcp.PenaltyContactProcess):
            def _interface_preprocess(self, partial_model_part, key="0"):
                pass  # layers already are the interface; GenerateFaces fails on QP geometry

        _LAYER_PROCESS_CLASS = SurfaceLayerPenaltyContactProcess
    return _LAYER_PROCESS_CLASS


class SurfaceLayerContact:
    """Unilateral penalty mortar contact between two surface layers + its solver.

    Construct after the model (elements, layers, ties, fixities) is complete.
    ``slave`` and ``master`` are names of sub model parts of ``model_part``
    created by :func:`add_tied_quad_layer`. ``penalty`` is Kratos'
    ``INITIAL_PENALTY`` (see the module docstring for its scaling).

    Then call :meth:`solve_step` once per load step. The strategy moves the
    mesh: after a step every node (control points included) sits at X0 + u, so
    read reference positions from ``X0``/``Y0``/``Z0``. Verified with
    ``SmallDisplacementElement`` on QP geometry (reference-configuration B).
    """

    def __init__(
        self,
        model: Any,
        model_part: Any,
        slave: str,
        master: str,
        penalty: float,
        contact_property_id: int,
        linear_solver: Any | None = None,
        max_iterations: int = 30,
        rel_tol: float = 1e-8,
        abs_tol: float = 1e-12,
        echo: int = 0,
    ) -> None:
        KM, CSMA = _kratos()
        if model_part.GetBufferSize() < 2:
            raise ValueError("SurfaceLayerContact: model_part buffer size must be >= 2 (previous-step DISPLACEMENT)")
        for name in (slave, master):
            if not model_part.HasSubModelPart(name):
                raise ValueError(f"SurfaceLayerContact: no sub model part {name!r}")
        self.KM, self.CSMA = KM, CSMA
        self.model_part = model_part
        self.slave, self.master = slave, master
        if not model_part.HasProperties(contact_property_id):
            model_part.CreateNewProperties(contact_property_id)
        cfg = {
            "model_part_name": model_part.FullName(),
            "contact_model_part": {"0": [slave, master]},
            "assume_master_slave": {"0": [slave]},
            "contact_property_ids": {"0": int(contact_property_id)},
            "contact_type": "Frictionless",
            "advance_ALM_parameters": {"manual_ALM": True, "penalty": float(penalty)},
        }
        self.process = _layer_process_class()(model, KM.Parameters(json.dumps(cfg)))
        info = model_part.ProcessInfo
        info[KM.IS_RESTARTED] = 1  # skip NormalCheckProcess (needs element faces)
        try:
            self.process.ExecuteInitialize()
        finally:
            info[KM.IS_RESTARTED] = 0  # elements must initialize their constitutive laws
        if info[KM.INITIAL_PENALTY] != float(penalty):
            raise RuntimeError(f"INITIAL_PENALTY is {info[KM.INITIAL_PENALTY]}, expected {penalty}")
        self.process.ExecuteBeforeSolutionLoop()
        model_part.Set(KM.INTERACTION, True)

        base = CSMA.DisplacementContactCriteria(rel_tol, abs_tol, rel_tol, abs_tol, False)
        mortar = CSMA.PenaltyFrictionlessMortarConvergenceCriteria(False, False, False)
        crit = CSMA.MortarAndConvergenceCriteria(base, mortar, False, None)
        crit.SetActualizeRHSFlag(True)
        for c in (base, mortar, crit):
            c.SetEchoLevel(int(echo))
        solver = linear_solver if linear_solver is not None else KM.SkylineLUFactorizationSolver()
        builder = CSMA.ContactResidualBasedBlockBuilderAndSolver(solver)
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        self.strategy = CSMA.ResidualBasedNewtonRaphsonContactStrategy(
            model_part, scheme, crit, builder, int(max_iterations),
            True,   # compute reactions
            True,   # reform DOFs each step (contact solver default)
            True,   # MoveMeshFlag: weighted gap / active set use current coordinates
            KM.Parameters("{}"), None, None)
        self.strategy.SetEchoLevel(int(echo))
        self.strategy.Initialize()
        self._criteria, self._echo = crit, int(echo)
        self.step = 0
        self.iterations = 0

    def solve_step(self, time: float, apply_loads: Callable[[], None]) -> bool:
        """Advance one load step and solve it. Returns whether the step converged.

        Order (as in ``AnalysisStage``): ``CloneTimeStep`` -> ``apply_loads()``
        (set prescribed DISPLACEMENT values / loads; never clone after this) ->
        contact search -> strategy initialize/predict/solve/finalize ->
        process finalize. Convergence is the return value of
        ``SolveSolutionStep``: ``strategy.IsConverged()`` after
        ``FinalizeSolutionStep`` re-runs the criteria on the system that
        ``reform_dofs`` has already cleared (access violation).
        """
        KM = self.KM
        mp = self.model_part
        self.step += 1
        mp.CloneTimeStep(float(time))
        mp.ProcessInfo[KM.STEP] = self.step
        apply_loads()
        self.process.ExecuteInitializeSolutionStep()
        self.strategy.InitializeSolutionStep()
        # The Newton base constructor sets the criteria echo level to 1 and the contact
        # strategy restores that value in InitializeSolutionStep; re-apply ours.
        self._criteria.SetEchoLevel(self._echo)
        self.strategy.Predict()
        converged = bool(self.strategy.SolveSolutionStep())
        self.iterations = int(mp.ProcessInfo[KM.NL_ITERATION_NUMBER])
        self.strategy.FinalizeSolutionStep()
        self.process.ExecuteFinalizeSolutionStep()
        return converged

    def slave_state(self) -> dict[int, tuple[bool, float]]:
        """``{slave layer node id: (ACTIVE, WEIGHTED_GAP)}`` after the last step."""
        KM, CSMA = self.KM, self.CSMA
        return {n.Id: (bool(n.Is(KM.ACTIVE)), float(n.GetSolutionStepValue(CSMA.WEIGHTED_GAP)))
                for n in self.model_part.GetSubModelPart(self.slave).Nodes}
