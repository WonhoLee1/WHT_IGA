"""Cell-wise quadrature-point elements on B-spline background grids (Track A).

Builds one ``QuadraturePointGeometry`` per grid cell (Kratos pitfall K2: a
single call of ``CreateQuadraturePointGeometries`` with ``nderiv=2`` packs all
given points of that cell into exactly one geometry) and attaches a
StructuralMechanics element to it (pitfall K5: the ``3D4N`` suffix is only a
prototype name; the real geometry is supplied at creation).

Ownership: Track A file (see implementation-request_v1.md parallel tracks).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from wht_iga.kratos_backend.geometry import cell_gauss_points


def attach_cell_elements(
    KM: Any,
    model_part: Any,
    volume: Any,
    box: Any,
    props: Any,
    element_name: str = "SmallDisplacementElement3D4N",
    start_id: int = 1,
    cells: list[tuple[int, int, int]] | None = None,
) -> list[Any]:
    """Create one QP geometry + one element per cell.

    Returns the created elements in cell order. ``volume`` is the NurbsVolume
    from ``make_nurbs_volume`` (shared helper, read-only use).
    """
    nx, ny, nz = box.n_cells
    if cells is None:
        cells = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    elements = []
    eid = start_id
    for i, j, k in cells:
        qps = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(qps, 2, cell_gauss_points(box, i, j, k))
        assert len(qps) == 1, f"cell {(i, j, k)}: expected 1 QP geometry, got {len(qps)}"
        elements.append(model_part.CreateNewElement(element_name, eid, qps[0], props))
        eid += 1
    return elements


def solve_linear(
    KM: Any,
    model_part: Any,
    use_pardiso: bool = False,
    echo: int = 0,
) -> Any:
    """Run a linear static strategy to completion and return it.

    The caller fixes/prescribes DOFs before calling. The returned strategy
    exposes ``GetSystemMatrix()`` for K extraction.
    """
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    if use_pardiso:
        import KratosMultiphysics.LinearSolversApplication as LA

        solver = LA.PardisoLUSolver()
    else:
        solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(solver)
    strategy = KM.ResidualBasedLinearStrategy(model_part, scheme, builder, False, False, False, False)
    strategy.SetEchoLevel(echo)
    strategy.Initialize()
    model_part.CloneTimeStep(1.0)
    strategy.Solve()
    return strategy


def system_matrix_csr(strategy: Any) -> Any:
    """Return the assembled global matrix as a scipy CSR matrix."""
    import KratosMultiphysics.scipy_conversion_tools as sct

    return sct.to_csr(strategy.GetSystemMatrix())


def equation_permutation(KM: Any, model_part: Any, box: Any) -> np.ndarray:
    """Reference-DOF index ``3*I+d`` -> Kratos equation id.

    Node ids are ``1 + control-point index`` (see shared ``make_nurbs_volume``).
    Requires the strategy to have run (equation ids are assigned at setup).
    """
    vars3 = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    perm = np.empty(box.n_dofs, dtype=np.int64)
    for node in model_part.Nodes:
        I = node.Id - 1
        if 0 <= I < box.n_ctrl:
            for d, var in enumerate(vars3):
                perm[3 * I + d] = node.GetDof(var).EquationId
    return perm
