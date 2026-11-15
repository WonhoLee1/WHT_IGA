"""Build Kratos geometries from WHT_IGA grid objects (shared by all spikes).

Ownership: shared module. Spike tracks must not modify it; request changes via the
integrator (see implementation-request_v1.md, parallel tracks).
"""

from __future__ import annotations

import numpy as np

from wht_iga.grid.knots import to_kratos
from wht_iga.reference.box import BoxSpline


def make_nurbs_volume(KM, model_part, box: BoxSpline, first_node_id: int = 1, add_dofs: bool = True):
    """Create one Kratos node per control point (Greville positions) and a NurbsVolumeGeometry.

    Node order matches ``BoxSpline``: control point I = a + n_u*(b + n_v*c) gets node id
    ``first_node_id + I``. The model part must already have DISPLACEMENT (and REACTION when
    add_dofs=True) registered as nodal solution-step variables (Kratos pitfall K13).

    Returns (volume_geometry, nodes_vector).

    Lifetime: the volume is registered on ``model_part`` (``AddGeometry``) so the model
    part owns it. Quadrature-point geometries keep a raw parent pointer; when the last
    Python reference to the volume was dropped, ``TotalLagrangianElement3D4N`` (whose
    reference Jacobian re-evaluates the parent at the point) read freed memory and gave
    intermittent ``DETJ0: -0.5`` / NaN stresses (S19 note, REVIEW_P0_P1 §1 #8; 2026-11-14
    integration run: 6/12 NaN without the reference, 0/12 with it).
    """
    X = box.control_points()
    nodes = KM.NodesVector()
    for I, (x, y, z) in enumerate(X):
        node = model_part.CreateNewNode(first_node_id + I, float(x), float(y), float(z))
        nodes.append(node)
    if add_dofs:
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
    p = box.degree
    kk = [KM.Vector(to_kratos(k, p).tolist()) for k in box.knots]
    volume = KM.NurbsVolumeGeometry(nodes, p, p, p, kk[0], kk[1], kk[2])
    model_part.AddGeometry(volume)
    return volume, nodes


def cell_gauss_points(box: BoxSpline, i: int, j: int, k: int, n_gauss: int | None = None) -> list[list[float]]:
    """Tensor Gauss points of cell (i, j, k) as [[x, y, z, w], ...] in physical = parameter space."""
    ng = box.degree + 1 if n_gauss is None else n_gauss
    xg, wg = np.polynomial.legendre.leggauss(ng)
    pts = []
    br = [box.breakpoints(d) for d in range(3)]
    lo = [br[0][i], br[1][j], br[2][k]]
    hi = [br[0][i + 1], br[1][j + 1], br[2][k + 1]]
    half = [(hi[d] - lo[d]) / 2 for d in range(3)]
    for gz, wz in zip(xg, wg):
        for gy, wy in zip(xg, wg):
            for gx, wx in zip(xg, wg):
                x = lo[0] + (gx + 1) * half[0]
                y = lo[1] + (gy + 1) * half[1]
                z = lo[2] + (gz + 1) * half[2]
                pts.append([x, y, z, wx * wy * wz * half[0] * half[1] * half[2]])
    return pts
