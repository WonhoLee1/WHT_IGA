"""Native Kratos FE element population, Jacobian checks, and IGA-FE Hybrid solver.

Implements Track P3-C (`T-P3.3`):

* Populates Kratos ``SmallDisplacementElement3D8N`` (``hex8``), ``3D4N``
  (``tet4``), ``3D10N`` (``tet10``), ``3D20N`` (``hex20``), and ``3D6N``
  (``wedge6``) from :class:`wht_iga.meshio_ext.model.MeshModel` into
  ``Parts.<name>`` submodelparts of the root ``Structure`` ``ModelPart``.
* Handles canonical-to-Kratos node reordering (crucially ``hex20`` midside
  edges ``[0..11, 16..19, 12..15]``) and validates positive Jacobian
  determinants across elements with deterministic auto-flip for inverted
  hex/tet/wedge cells.
* Provides :func:`solve_hybrid_iga_fe` to solve coupled Native FE + IGA (or
  multi-part FE/IGA/THB) systems in a unified Kratos ``Structure`` model part
  via :mod:`wht_iga.kratos_backend.interactions`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from wht_iga.meshio_ext.model import MeshModel


KRATOS_FE_ELEMENT_MAP: dict[str, str] = {
    "hex8": "SmallDisplacementElement3D8N",
    "tet4": "SmallDisplacementElement3D4N",
    "tet10": "SmallDisplacementElement3D10N",
    "hex20": "SmallDisplacementElement3D20N",
    "wedge6": "SmallDisplacementElement3D6N",
}

_EXPECTED_NNODES: dict[str, int] = {
    "hex8": 8,
    "tet4": 4,
    "tet10": 10,
    "hex20": 20,
    "wedge6": 6,
}

# MeshModel hex20: corners 0..7, bottom mids 8..11, top mids 12..15, vertical mids 16..19.
# Kratos Hexahedra3D20: corners 0..7, bottom mids 8..11, vertical mids 12..15, top mids 16..19.
_HEX20_TO_KRATOS = np.array(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15],
    dtype=np.int64,
)

# Reference signs and flip permutations in MeshModel canonical order.
_HEX_CORNER_SIGNS = np.array(
    [
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1],
    ],
    dtype=np.float64,
)
_HEX_GAUSS_PTS = np.array(
    [
        [sx / np.sqrt(3.0), sy / np.sqrt(3.0), sz / np.sqrt(3.0)]
        for sz in (-1.0, 1.0)
        for sy in (-1.0, 1.0)
        for sx in (-1.0, 1.0)
    ],
    dtype=np.float64,
)
_FLIP_HEX8 = np.array([4, 5, 6, 7, 0, 1, 2, 3], dtype=np.int64)
_FLIP_HEX20 = np.array(
    [4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 16, 17, 18, 19],
    dtype=np.int64,
)
_FLIP_TET4 = np.array([0, 2, 1, 3], dtype=np.int64)
_FLIP_TET10 = np.array([0, 2, 1, 3, 6, 5, 4, 7, 9, 8], dtype=np.int64)
_FLIP_WEDGE6 = np.array([3, 4, 5, 0, 1, 2], dtype=np.int64)


def _hex8_detj_samples(pts: np.ndarray) -> np.ndarray:
    """Evaluate trilinear hex8 Jacobian determinant at center + 8 Gauss points."""
    corners = pts[:8]
    sample_pts = np.vstack([np.zeros((1, 3), dtype=np.float64), _HEX_GAUSS_PTS])
    dets = np.empty(len(sample_pts), dtype=np.float64)
    for k, (xi, eta, zeta) in enumerate(sample_pts):
        dN_dxi = 0.125 * _HEX_CORNER_SIGNS[:, 0] * (1.0 + _HEX_CORNER_SIGNS[:, 1] * eta) * (1.0 + _HEX_CORNER_SIGNS[:, 2] * zeta)
        dN_deta = 0.125 * (1.0 + _HEX_CORNER_SIGNS[:, 0] * xi) * _HEX_CORNER_SIGNS[:, 1] * (1.0 + _HEX_CORNER_SIGNS[:, 2] * zeta)
        dN_dzeta = 0.125 * (1.0 + _HEX_CORNER_SIGNS[:, 0] * xi) * (1.0 + _HEX_CORNER_SIGNS[:, 1] * eta) * _HEX_CORNER_SIGNS[:, 2]
        J = np.column_stack([corners.T @ dN_dxi, corners.T @ dN_deta, corners.T @ dN_dzeta])
        dets[k] = float(np.linalg.det(J))
    return dets


def _wedge6_detj_samples(pts: np.ndarray) -> np.ndarray:
    """Evaluate wedge6 prism Jacobian determinant at bottom and top Gauss points."""
    dets = np.empty(2, dtype=np.float64)
    for idx, zeta in enumerate((-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))):
        r, s = 1.0 / 3.0, 1.0 / 3.0
        # N0..2 = (1-r-s, r, s)*(1-zeta)/2; N3..5 = (1-r-s, r, s)*(1+zeta)/2
        f_lo = 0.5 * (1.0 - zeta)
        f_hi = 0.5 * (1.0 + zeta)
        dN_dr = np.array([-f_lo, f_lo, 0.0, -f_hi, f_hi, 0.0], dtype=np.float64)
        dN_ds = np.array([-f_lo, 0.0, f_lo, -f_hi, 0.0, f_hi], dtype=np.float64)
        dN_dz = 0.5 * np.array([-(1.0 - r - s), -r, -s, (1.0 - r - s), r, s], dtype=np.float64)
        J = np.column_stack([pts.T @ dN_dr, pts.T @ dN_ds, pts.T @ dN_dz])
        dets[idx] = float(np.linalg.det(J))
    return dets


def check_element_jacobian(
    elem_type: str,
    elem_coords: np.ndarray,
) -> tuple[float, float]:
    """Return ``(min_detJ, center_detJ)`` for one element in canonical ``MeshModel`` order."""
    pts = np.asarray(elem_coords, dtype=np.float64)
    if elem_type in ("tet4", "tet10"):
        J = np.column_stack([pts[1] - pts[0], pts[2] - pts[0], pts[3] - pts[0]])
        det = float(np.linalg.det(J))
        return det, det
    if elem_type in ("hex8", "hex20"):
        dets = _hex8_detj_samples(pts)
        return float(np.min(dets)), float(dets[0])
    if elem_type == "wedge6":
        dets = _wedge6_detj_samples(pts[:6])
        return float(np.min(dets)), float(np.mean(dets))
    raise ValueError(f"Unsupported FE element type '{elem_type}'")


def ensure_positive_jacobian(
    elem_type: str,
    conn_row: np.ndarray,
    coords: np.ndarray,
    auto_flip: bool = True,
) -> tuple[np.ndarray, bool]:
    """Verify positive Jacobian determinant and optionally flip inverted node winding.

    Parameters
    ----------
    elem_type:
        Canonical element type: ``"hex8"``, ``"tet4"``, ``"tet10"``, ``"hex20"``,
        or ``"wedge6"``.
    conn_row:
        1D integer node indices in canonical :class:`MeshModel` node order.
    coords:
        ``(N, 3)`` mesh nodal coordinate array.
    auto_flip:
        If ``True``, automatically permute inverted elements so ``det(J) > 0``.
        If ``False``, raise ``ValueError`` when an inverted element is encountered.

    Returns
    -------
    tuple[np.ndarray, bool]
        Validated (possibly flipped) connectivity row in canonical order, plus
        boolean indicating whether a flip was performed.
    """
    row = np.asarray(conn_row, dtype=np.int64)
    expected = _EXPECTED_NNODES[elem_type]
    if row.shape != (expected,):
        raise ValueError(
            f"Element '{elem_type}' expects {expected} nodes, got {row.shape}"
        )

    min_det, center_det = check_element_jacobian(elem_type, coords[row])
    if min_det > 0.0:
        return row.copy(), False

    if center_det < 0.0:
        if not auto_flip:
            raise ValueError(
                f"Inverted '{elem_type}' element detected (detJ={center_det:.3e} <= 0)"
            )
        flip_map = {
            "tet4": _FLIP_TET4,
            "tet10": _FLIP_TET10,
            "hex8": _FLIP_HEX8,
            "hex20": _FLIP_HEX20,
            "wedge6": _FLIP_WEDGE6,
        }[elem_type]
        flipped = row[flip_map]
        min_det_after, center_after = check_element_jacobian(elem_type, coords[flipped])
        if min_det_after <= 0.0:
            raise ValueError(
                f"Degenerate '{elem_type}' element remains non-positive after flip "
                f"(min_detJ={min_det_after:.3e}, center_detJ={center_after:.3e})"
            )
        return flipped, True

    raise ValueError(
        f"Tangled/degenerate '{elem_type}' element (center_detJ={center_det:.3e}, min_detJ={min_det:.3e})"
    )


def to_kratos_connectivity(elem_type: str, canonical_ids: np.ndarray | list[int]) -> list[int]:
    """Map canonical :class:`MeshModel` element node IDs to Kratos element node order."""
    arr = np.asarray(canonical_ids, dtype=np.int64)
    if elem_type == "hex20":
        return arr[_HEX20_TO_KRATOS].tolist()
    return arr.tolist()


@dataclass
class NativeFEPartHandle:
    """Runtime metadata produced by :func:`populate_native_fe_part`."""

    sub_model_part: Any
    mesh: MeshModel
    coords_global: np.ndarray
    node_ids: list[int]
    mesh_to_kratos_node: dict[int, int]
    element_ids: list[int]
    flipped_count: int
    element_types_used: list[str] = field(default_factory=list)


def populate_native_fe_part(
    KM: Any,
    sub_model_part: Any,
    mesh: MeshModel,
    props: Any,
    id_table: Any | None = None,
    scope: str = "Parts.fe",
    node_id_start: int = 1,
    elem_id_start: int = 1,
    auto_flip: bool = True,
    node_transform: Callable[[np.ndarray], np.ndarray] | None = None,
    block_names: list[str] | None = None,
    props_by_part_id: dict[int, Any] | None = None,
) -> NativeFEPartHandle:
    """Populate a Kratos ``ModelPart`` from a :class:`MeshModel`.

    Creates nodes once across all selected blocks, registers translational DOFs
    (``DISPLACEMENT`` / ``REACTION``), checks/flips element Jacobians, reorders
    canonical ``hex20`` midside nodes into Kratos ``Hexahedra3D20`` order, and
    creates ``SmallDisplacementElement3D*`` elements.
    """
    coords_global = (
        np.asarray(node_transform(mesh.nodes), dtype=np.float64)
        if node_transform is not None
        else np.asarray(mesh.nodes, dtype=np.float64)
    )

    blocks = [
        b
        for b in mesh.blocks
        if block_names is None
        or b.elem_type in block_names
        or mesh.parts.get(b.part_id) in block_names
    ]
    if not blocks:
        raise ValueError(f"No matching element blocks found in MeshModel (filter={block_names})")

    used_vertices = sorted({int(v) for b in blocks for v in b.conn.ravel()})
    n_nodes = len(used_vertices)
    total_elems = sum(len(b.conn) for b in blocks)

    if id_table is not None:
        nid0 = id_table.allocate_nodes(scope, n_nodes, provenance="fe_node").start_id
        eid0 = id_table.allocate_elements(scope, total_elems, kind="NativeFE").start_id
    else:
        nid0 = int(node_id_start)
        eid0 = int(elem_id_start)

    mesh_to_kratos: dict[int, int] = {}
    node_ids: list[int] = []
    for offset, m_idx in enumerate(used_vertices):
        kid = nid0 + offset
        x, y, z = coords_global[m_idx]
        node = sub_model_part.CreateNewNode(kid, float(x), float(y), float(z))
        node.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
        node.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
        node.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)
        mesh_to_kratos[m_idx] = kid
        node_ids.append(kid)

    element_ids: list[int] = []
    flipped_count = 0
    elem_types_used: list[str] = []
    eid = eid0

    for block in blocks:
        etype = block.elem_type
        if etype not in KRATOS_FE_ELEMENT_MAP:
            raise ValueError(f"Unsupported FE element type '{etype}'")
        kratos_elem_name = KRATOS_FE_ELEMENT_MAP[etype]
        if kratos_elem_name not in elem_types_used:
            elem_types_used.append(kratos_elem_name)

        elem_props = (
            props_by_part_id[block.part_id]
            if (props_by_part_id is not None and block.part_id in props_by_part_id)
            else props
        )

        for conn_row in block.conn:
            valid_row, was_flipped = ensure_positive_jacobian(
                etype, conn_row, coords_global, auto_flip=auto_flip
            )
            if was_flipped:
                flipped_count += 1
            k_canon_ids = [mesh_to_kratos[int(v)] for v in valid_row]
            k_node_ids = to_kratos_connectivity(etype, k_canon_ids)
            sub_model_part.CreateNewElement(kratos_elem_name, eid, k_node_ids, elem_props)
            element_ids.append(eid)
            eid += 1

    return NativeFEPartHandle(
        sub_model_part=sub_model_part,
        mesh=mesh,
        coords_global=coords_global,
        node_ids=node_ids,
        mesh_to_kratos_node=mesh_to_kratos,
        element_ids=element_ids,
        flipped_count=flipped_count,
        element_types_used=elem_types_used,
    )


@dataclass
class HybridSolveResult:
    """Outputs returned by :func:`solve_hybrid_iga_fe`."""

    built: Any
    tie_results: list[Any]
    u_vector: np.ndarray
    node_ids_ordered: list[int]
    max_displacement: float


def _kratos_matrix_to_csr(KM: Any, kmat: Any) -> sp.csr_matrix:
    """Convert Kratos compressed CSR matrix cleanly to SciPy CSR."""
    del KM
    n1 = int(kmat.Size1())
    n2 = int(kmat.Size2())
    if hasattr(kmat, "index1_data") and hasattr(kmat, "index2_data") and hasattr(kmat, "value_data"):
        indptr_raw = kmat.index1_data() if callable(kmat.index1_data) else kmat.index1_data
        indices_raw = kmat.index2_data() if callable(kmat.index2_data) else kmat.index2_data
        data_raw = kmat.value_data() if callable(kmat.value_data) else kmat.value_data
        indptr = np.asarray(indptr_raw, dtype=np.int64)
        indices = np.asarray(indices_raw, dtype=np.int64)
        data = np.asarray(data_raw, dtype=np.float64)
        return sp.csr_matrix((data, indices, indptr), shape=(n1, n2))
    rows, cols, vals = [], [], []
    for i in range(n1):
        for j in range(n2):
            v = float(kmat[i, j])
            if v != 0.0:
                rows.append(i)
                cols.append(j)
                vals.append(v)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n1, n2))


def solve_hybrid_iga_fe(
    KM: Any,
    SM: Any,
    built: Any,
    ties: list[Any] | None = None,
    rbe2s: list[Any] | None = None,
    rbe3_specs_and_loads: list[tuple[Any, tuple[float, float, float], tuple[float, float, float]]] | None = None,
    use_pardiso: bool = False,
) -> HybridSolveResult:
    """Assemble and solve a hybrid Native FE + Immersed IGA (or multi-part) assembly.

    Applies all configured Tie, RBE2, and RBE3 couplings to ``built.structure``.
    If all ties compile to pure ``LinearMasterSlaveConstraint``s (``mpc_matching``
    or ``mortar_l2``), runs Kratos ``ResidualBasedLinearStrategy`` natively via
    :func:`wht_iga.kratos_backend.elements.solve_linear`. If any tie resolves to
    surface quadrature ``penalty_spring`` (e.g. IGA-FE or IGA-IGA interface),
    builds the unconstrained Kratos ``(K_base, b_base)``, injects ``K_tie``,
    applies any Kratos ``MasterSlaveConstraint``s and Dirichlet conditions via
    ``T_mpc`` reduction, solves with sparse direct solver, and writes solutions
    back onto every Kratos node.
    """
    from wht_iga.kratos_backend.elements import solve_linear
    from wht_iga.kratos_backend.interactions import (
        apply_rbe2,
        apply_rbe3,
        apply_tie,
    )

    model_part = built.structure
    tie_results: list[Any] = []
    for ts in (ties or []):
        res = apply_tie(
            KM,
            model_part,
            ts,
            id_table=built.id_table,
            couplings_smp=built.couplings_smp,
        )
        tie_results.append(res)

    for r2 in (rbe2s or []):
        apply_rbe2(
            KM,
            model_part,
            r2,
            id_table=built.id_table,
            couplings_smp=built.couplings_smp,
        )

    for r3_spec, f_vec, m_vec in (rbe3_specs_and_loads or []):
        apply_rbe3(
            KM,
            SM,
            model_part,
            r3_spec,
            force=f_vec,
            moment=m_vec,
            mode="mpc",
            id_table=built.id_table,
            couplings_smp=built.couplings_smp,
        )

    penalty_ties = [tr for tr in tie_results if tr.penalty_matrix is not None]
    node_ids_ordered = sorted(int(nd.Id) for nd in model_part.Nodes)

    if not penalty_ties:
        solve_linear(KM, model_part, use_pardiso=use_pardiso, echo=0)
        disp_arr = np.array(
            [list(model_part.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
             for nid in node_ids_ordered],
            dtype=np.float64,
        )
        return HybridSolveResult(
            built=built,
            tie_results=tie_results,
            u_vector=disp_arr.ravel(),
            node_ids_ordered=node_ids_ordered,
            max_displacement=float(np.max(np.abs(disp_arr))) if len(disp_arr) else 0.0,
        )

    # Zero free DOFs cleanly before assembling base unconstrained (K_base, b_base).
    disp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    rot_vars = (KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z)
    for nd in model_part.Nodes:
        for var in (*disp_vars, *rot_vars):
            if nd.HasDofFor(var) and not nd.IsFixed(var):
                nd.SetSolutionStepValue(var, 0, 0.0)

    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    linear_solver = KM.SkylineLUFactorizationSolver()
    builder = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
    strategy = KM.ResidualBasedLinearStrategy(
        model_part, scheme, builder, False, False, False, False
    )
    strategy.SetEchoLevel(0)
    strategy.Initialize()
    model_part.CloneTimeStep(1.0)
    builder.SetUpDofSet(scheme, model_part)
    builder.SetUpSystem(model_part)

    A_mat = KM.CompressedMatrix()
    b_vec = KM.Vector()
    dx_vec = KM.Vector()
    builder.ResizeAndInitializeVectors(scheme, A_mat, dx_vec, b_vec, model_part)
    builder.Build(scheme, model_part, A_mat, b_vec)

    K_base = _kratos_matrix_to_csr(KM, A_mat)
    rhs_base = np.asarray(b_vec, dtype=np.float64).copy()
    neq = K_base.shape[0]

    # Inject all interface penalty matrices K_tie into global equation numbering.
    K_pen_total = sp.csr_matrix((neq, neq), dtype=np.float64)
    for tr in penalty_ties:
        pen_nids = tr.penalty_node_ids
        K_loc = tr.penalty_matrix.tocoo()
        eq_lookup = np.full(3 * len(pen_nids), -1, dtype=np.int64)
        for pos, nid in enumerate(pen_nids):
            nd = model_part.GetNode(int(nid))
            for c, dv in enumerate(disp_vars):
                eq_lookup[3 * pos + c] = int(nd.GetDof(dv).EquationId)
        gr = eq_lookup[K_loc.row]
        gc = eq_lookup[K_loc.col]
        valid = (gr >= 0) & (gc >= 0)
        K_glob = sp.coo_matrix(
            (K_loc.data[valid], (gr[valid], gc[valid])), shape=(neq, neq)
        ).tocsr()
        K_pen_total = K_pen_total + K_glob

    K_full = (K_base + K_pen_total).tocsr()
    f_full = rhs_base

    # Apply any Kratos MasterSlaveConstraints via exact T_mpc prolongation u_full = T_mpc u_act.
    if model_part.NumberOfMasterSlaveConstraints() > 0:
        pi = model_part.ProcessInfo
        slave_eqs: set[int] = set()
        t_rows: list[int] = []
        t_cols_eq: list[int] = []
        t_vals: list[float] = []
        for msc in model_part.MasterSlaveConstraints:
            rmat = KM.Matrix()
            cvec = KM.Vector()
            msc.CalculateLocalSystem(rmat, cvec, pi)
            s_dofs = msc.GetSlaveDofsVector()
            m_dofs = msc.GetMasterDofsVector()
            R_np = np.array(rmat, dtype=np.float64)
            for ri, sdof in enumerate(s_dofs):
                seq = int(sdof.EquationId)
                slave_eqs.add(seq)
                for ci, mdof in enumerate(m_dofs):
                    w = float(R_np[ri, ci])
                    if w != 0.0:
                        t_rows.append(seq)
                        t_cols_eq.append(int(mdof.EquationId))
                        t_vals.append(w)

        active_eqs = [eq for eq in range(neq) if eq not in slave_eqs]
        eq_to_act = {eq: col for col, eq in enumerate(active_eqs)}
        for col, eq in enumerate(active_eqs):
            t_rows.append(eq)
            t_cols_eq.append(eq)
            t_vals.append(1.0)
        t_cols = [eq_to_act[meq] for meq in t_cols_eq]
        T_mpc = sp.csr_matrix(
            (t_vals, (t_rows, t_cols)), shape=(neq, len(active_eqs))
        )
        K_sys = (T_mpc.T @ K_full @ T_mpc).tocsr()
        f_sys = np.asarray(T_mpc.T @ f_full, dtype=np.float64)
    else:
        active_eqs = list(range(neq))
        eq_to_act = {eq: eq for eq in active_eqs}
        T_mpc = None
        K_sys = K_full
        f_sys = f_full

    # Apply fixed Dirichlet BCs on active DOFs via exact partitioning.
    n_act = len(active_eqs)
    fixed_act: dict[int, float] = {}
    for nd in model_part.Nodes:
        for var in (*disp_vars, *rot_vars):
            if nd.HasDofFor(var) and nd.IsFixed(var):
                eq_id = int(nd.GetDof(var).EquationId)
                if eq_id in eq_to_act:
                    fixed_act[eq_to_act[eq_id]] = float(nd.GetSolutionStepValue(var, 0))

    u_act = np.zeros(n_act, dtype=np.float64)
    if fixed_act:
        fixed_idx = sorted(fixed_act.keys())
        for idx in fixed_idx:
            u_act[idx] = fixed_act[idx]
        free_idx = np.setdiff1d(np.arange(n_act, dtype=np.int64), np.asarray(fixed_idx, dtype=np.int64))
        rhs_mod = f_sys[free_idx] - K_sys[free_idx, :][:, fixed_idx] @ u_act[fixed_idx]
        K_ff = K_sys[free_idx, :][:, free_idx].tocsc()
        u_act[free_idx] = spla.spsolve(K_ff, rhs_mod)
    else:
        u_act = spla.spsolve(K_sys.tocsc(), f_sys)

    u_sol = np.asarray(T_mpc @ u_act, dtype=np.float64) if T_mpc is not None else u_act

    for nd in model_part.Nodes:
        for var in (*disp_vars, *rot_vars):
            if nd.HasDofFor(var) and not nd.IsFixed(var):
                eq_id = int(nd.GetDof(var).EquationId)
                if 0 <= eq_id < neq:
                    nd.SetSolutionStepValue(var, 0, float(u_sol[eq_id]))

    disp_arr = np.array(
        [list(model_part.GetNode(nid).GetSolutionStepValue(KM.DISPLACEMENT, 0))
         for nid in node_ids_ordered],
        dtype=np.float64,
    )
    return HybridSolveResult(
        built=built,
        tie_results=tie_results,
        u_vector=disp_arr.ravel(),
        node_ids_ordered=node_ids_ordered,
        max_displacement=float(np.max(np.abs(disp_arr))) if len(disp_arr) else 0.0,
    )


__all__ = [
    "HybridSolveResult",
    "KRATOS_FE_ELEMENT_MAP",
    "NativeFEPartHandle",
    "check_element_jacobian",
    "ensure_positive_jacobian",
    "populate_native_fe_part",
    "solve_hybrid_iga_fe",
    "to_kratos_connectivity",
]
