"""Kratos-independent Phase 3 post-processing package for WHT_IGA.

Provides visualization proxy mesh generation (hex/tet + B-spline evaluation),
VTKHDF writer/reader with ``/WHT_IGA/control_points`` storage, and
quadrature-based ``L2``, ``H1``, and energy error norm calculators.
"""

from __future__ import annotations

from wht_iga.post.norms import (
    H1NormResult,
    NormResult,
    compute_energy_norm,
    compute_h1_norm,
    compute_l2_norm,
    gauss_legendre_hex_rule,
    integrate_error_against_analytic,
    integrate_error_on_tets,
)
from wht_iga.post.proxy_mesh import (
    CELL_CUT,
    CELL_IN,
    CELL_OUT,
    VTK_HEXAHEDRON,
    VTK_TETRA,
    BSplineFieldEvaluator,
    ProxyMesh,
    build_proxy_mesh,
    decompose_hex_to_6_tets,
    generate_cut_tet_proxy,
    generate_hex_proxy_for_cells,
)
from wht_iga.post.probes import (
    LineProbe,
    MultiPointProbe,
    PointProbe,
    ProbeSample,
    create_zero_weight_probe_element,
    evaluate_zero_weight_kratos_probes,
    export_probes_json,
    export_probes_sqlite,
)
from wht_iga.post.vtk_writer import read_vtkhdf, write_vtkhdf

__all__ = [
    "CELL_CUT",
    "CELL_IN",
    "CELL_OUT",
    "VTK_HEXAHEDRON",
    "VTK_TETRA",
    "BSplineFieldEvaluator",
    "H1NormResult",
    "LineProbe",
    "MultiPointProbe",
    "NormResult",
    "PointProbe",
    "ProbeSample",
    "ProxyMesh",
    "build_proxy_mesh",
    "compute_energy_norm",
    "compute_h1_norm",
    "compute_l2_norm",
    "create_zero_weight_probe_element",
    "decompose_hex_to_6_tets",
    "evaluate_zero_weight_kratos_probes",
    "export_probes_json",
    "export_probes_sqlite",
    "gauss_legendre_hex_rule",
    "generate_cut_tet_proxy",
    "generate_hex_proxy_for_cells",
    "integrate_error_against_analytic",
    "integrate_error_on_tets",
    "read_vtkhdf",
    "write_vtkhdf",
]
