"""Boundary-surface quadrature on immersed volumes (Track B).

Builds ``SurfaceInNurbsVolumeGeometry`` wrappers over physical triangles and
creates boundary quadrature-point geometries for IgaApplication
``LoadCondition`` / ``SupportPenaltyCondition``. ``SupportNitscheCondition`` is NOT
usable in 3D on these geometries (S18: segfault, membrane-only) -- C++ track.

Lifetime: each ``SurfaceInNurbsVolumeGeometry`` wrapper is registered on the model
part (``AddGeometry``) so its QP geometries never point at a freed parent (same
cause as the S19 TL NaN, see ``geometry.make_nurbs_volume``). The returned
``wrappers`` list is kept for backward compatibility.

Ownership: Track B file. May be promoted to shared use by the integrator.
"""

from __future__ import annotations

import numpy as np


def triangle_quadrature(order: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Dunavant rules on the reference triangle {(xi, eta): xi>=0, eta>=0, xi+eta<=1}.

    Returns (points (n, 2), weights (n,)) with weights summing to 1/2
    (the reference-triangle area). Only non-negative-weight rules.
    """
    if order <= 1:
        pts = np.array([[1.0 / 3.0, 1.0 / 3.0]])
        wts = np.array([0.5])
    elif order == 2:
        pts = np.array([[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]])
        wts = np.full(3, 1.0 / 6.0)
    elif order == 4:
        # Dunavant degree-4 6-point rule (barycentric -> (xi, eta)).
        a, wa = 0.445948490915965, 0.11169079483905
        b, wb = 0.091576213509771, 0.0549758718227661
        bary = [(a, a, 1.0 - 2.0 * a), (a, 1.0 - 2.0 * a, a), (1.0 - 2.0 * a, a, a),
                (b, b, 1.0 - 2.0 * b), (b, 1.0 - 2.0 * b, b), (1.0 - 2.0 * b, b, b)]
        pts = np.array([[l1, l2] for (l1, l2, _l3) in bary])
        wts = np.array([wa, wa, wa, wb, wb, wb])
    elif order == 6:
        # Dunavant degree-6 12-point rule (all weights positive; Dunavant 1985 table).
        a1, w1 = 0.249286745170910, 0.116786275726379
        a2, w2 = 0.063089014491502, 0.050844906370207
        b1, b2, w3 = 0.053145049844817, 0.310352451033784, 0.082851075618374
        bary = [(a1, a1, 1.0 - 2.0 * a1), (a1, 1.0 - 2.0 * a1, a1), (1.0 - 2.0 * a1, a1, a1),
                (a2, a2, 1.0 - 2.0 * a2), (a2, 1.0 - 2.0 * a2, a2), (1.0 - 2.0 * a2, a2, a2)]
        b3 = 1.0 - b1 - b2
        bary += [(b1, b2, b3), (b2, b1, b3), (b1, b3, b2), (b3, b1, b2), (b2, b3, b1), (b3, b2, b1)]
        pts = np.array([[l1, l2] for (l1, l2, _l3) in bary])
        wts = 0.5 * np.array([w1] * 3 + [w2] * 3 + [w3] * 6)
    else:
        raise ValueError(f"triangle_quadrature: order {order} not available (use 1, 2, 4 or 6)")
    return pts, wts


def triangle_area(tri: np.ndarray) -> float:
    """Physical area of a (3, 3) triangle."""
    e1 = tri[1] - tri[0]
    e2 = tri[2] - tri[0]
    return 0.5 * float(np.linalg.norm(np.cross(e1, e2)))


def triangle_normal(tri: np.ndarray) -> np.ndarray:
    """Unit normal from winding (right-hand rule)."""
    e1 = tri[1] - tri[0]
    e2 = tri[2] - tri[0]
    n = np.cross(e1, e2)
    return n / float(np.linalg.norm(n))


def make_surface_quadrature(
    KM,
    model_part,
    volume,
    triangles: list[np.ndarray],
    rule_order: int = 2,
    first_node_id: int = 10_000,
    shift_inward: float = 0.0,
) -> tuple[list, list, list[np.ndarray], list[float]]:
    """Create surface quadrature-point geometries for physical triangles.

    Each triangle gets carrier ``Triangle3D3`` nodes (created *after* the caller
    added DOFs, so carriers stay DOF-free and never enter the DofSet), one
    ``SurfaceInNurbsVolumeGeometry`` wrapper, and quadrature-point geometries
    from ``CreateQuadraturePointGeometries(g, 2, [[xi, eta, 0, w], ...])``.

    LIFETIME REQUIREMENT (spike S08 finding): the caller must keep the
    ``volume`` (``NurbsVolumeGeometry``) Python object alive for as long as the
    returned QP geometries are assembled/solved. Dropping the last Python
    reference silently zeroes the assembled condition RHS (zero-residual false
    convergence, reactions NaN) with no exception. The wrapper geometries are
    therefore also returned so the caller can retain them.

    ``shift_inward`` optionally pulls triangle vertices towards the volume
    interior along the inverted face normal (K3 safeguard for trimmed faces);
    0.0 keeps exact positions (used for untrimmed box faces).

    (``SurfaceInNurbsVolumeGeometry::CreateQuadraturePointGeometries`` override:
    ``surface_in_nurbs_volume_geometry.h`` l.278 — one QP geometry **per**
    integration point, each with its own nonzero control-point set.)

    Returns (qp_geometries, surface_wrappers, unit_normals, triangle_areas)
    with one QP geometry per quadrature point
    (``len == n_triangles * n_rule_points``).
    """
    ref_pts, ref_w = triangle_quadrature(rule_order)
    qp_geoms: list = []
    wrappers: list = []
    normals: list[np.ndarray] = []
    areas: list[float] = []
    node_id = first_node_id
    for tri in triangles:
        tri = np.asarray(tri, dtype=float)
        if shift_inward > 0.0:
            n = triangle_normal(tri)
            tri = tri - shift_inward * n
        tri_nodes = []
        for v in tri:
            tri_nodes.append(model_part.CreateNewNode(node_id, float(v[0]), float(v[1]), float(v[2])))
            node_id += 1
        flat = KM.Triangle3D3(tri_nodes[0], tri_nodes[1], tri_nodes[2])
        surf = KM.SurfaceInNurbsVolumeGeometry(volume, flat)
        qp_list = [[float(xi), float(eta), 0.0, float(w)] for (xi, eta), w in zip(ref_pts, ref_w)]
        out = KM.GeometriesVector()
        surf.CreateQuadraturePointGeometries(out, 2, qp_list)
        if len(out) != len(qp_list):
            raise RuntimeError(f"expected {len(qp_list)} surface QP geometries, got {len(out)}")
        model_part.AddGeometry(surf)
        qp_geoms.extend(out[i] for i in range(len(out)))
        wrappers.append(surf)
        normals.append(triangle_normal(np.asarray(tri, dtype=float)))
        areas.append(triangle_area(np.asarray(tri, dtype=float)))
    return qp_geoms, wrappers, normals, areas
