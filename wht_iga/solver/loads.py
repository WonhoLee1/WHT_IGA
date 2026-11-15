"""Load conditions: pressure / traction / follower + gravity (Track SOLVE).

Builders wrap the S08-proven ``LoadCondition`` patterns (``surface.py``
reused read-only, never modified):
- ``PRESSURE`` (core): traction = +p·n_out — S08: +p pushes along the
  outward face normal (documented sign convention).
- ``SURFACE_LOAD`` (SM): traction vector integrated as given.
- ``PRESSURE_FOLLOWER_LOAD`` (IgaApp): needs move_mesh Newton solve.
- Total-force pressure: p = F / Σw with Σw = tessellation face area.
- Gravity/body force: ``Properties`` ``VOLUME_ACCELERATION`` × ``DENSITY``
  (``structural_mechanics_element_utilities.cpp`` l.88-96).
"""

from __future__ import annotations

import numpy as np

from wht_iga.kratos_backend.surface import make_surface_quadrature
from wht_iga.solver.penalty import _qp_physical_xy, box_face_triangles, face_subdiv, rule_for_degree
from wht_iga.solver.validate import PressureLoad, TractionLoad


def resultant_pressure(force: float, triangles: list[np.ndarray]) -> float:
    """Uniform pressure giving total normal force ``force`` (p = F / Σw)."""
    from wht_iga.kratos_backend.surface import triangle_area

    total = sum(triangle_area(np.asarray(t)) for t in triangles)
    return float(force) / float(total)


def attach_loads(
    KM: object,
    SM: object,
    IGA: object,
    model_part: object,
    volume: object,
    box: object,
    pressures: list[PressureLoad],
    tractions: list[TractionLoad],
    first_cond_id: int = 10_000,
    first_node_id: int = 50_000,
    first_prop_id: int = 500,
) -> dict:
    """Attach ``LoadCondition``s on box faces. Returns keepalive refs + totals."""
    keep: dict = {"wrappers": [], "qp_geoms": [], "n_cond": 0, "total_force": np.zeros(3)}
    cid = first_cond_id
    nid = first_node_id
    pid = first_prop_id
    order = rule_for_degree(box.degree)
    jobs: list[tuple[tuple[int, int], str, np.ndarray | float, bool]] = [
        (p.face, "pressure", float(p.magnitude), bool(p.follower)) for p in pressures
    ] + [(t.face, "traction", np.asarray(t.vector, dtype=float), False) for t in tractions]
    for face, kind, val, follower in jobs:
        tris = box_face_triangles(np.asarray(box.lo), np.asarray(box.hi), *face,
                                  subdiv=face_subdiv(box, face[0]))
        qp_geoms, wrappers, normals, _ = make_surface_quadrature(
            KM, model_part, volume, tris, rule_order=order, first_node_id=nid)
        nid += 3 * len(tris)
        xs, dAs = _qp_physical_xy(tris, order)
        assert len(qp_geoms) == len(xs), "one surface QP geometry per integration point"
        n_out = np.asarray(normals[0], dtype=float)  # planar face: single outward normal
        for qp, x, dA in zip(qp_geoms, xs, dAs):
            props = model_part.CreateNewProperties(pid)
            pid += 1
            cond = model_part.CreateNewCondition("LoadCondition", cid, qp, props)
            cid += 1
            if kind == "pressure":
                var = IGA.PRESSURE_FOLLOWER_LOAD if follower else KM.PRESSURE
                cond.SetValue(var, float(val))
                keep["total_force"] = keep["total_force"] + float(val) * float(dA) * n_out
            else:
                v = KM.Vector(3)
                v[0], v[1], v[2] = float(val[0]), float(val[1]), float(val[2])
                cond.SetValue(SM.SURFACE_LOAD, v)
                keep["total_force"] = keep["total_force"] + np.asarray(val, dtype=float) * float(dA)
        keep["wrappers"].extend(wrappers)
        keep["qp_geoms"].extend(qp_geoms)
        keep["n_cond"] += len(qp_geoms)
    keep["total_force"] = np.asarray(keep["total_force"], dtype=float)
    return keep


def apply_gravity(props: object, KM: object, density: float, gravity: tuple[float, float, float]) -> None:
    """Set DENSITY + VOLUME_ACCELERATION on volume properties (body force)."""
    props.SetValue(KM.DENSITY, float(density))
    vec = KM.Vector(3)
    vec[0], vec[1], vec[2] = float(gravity[0]), float(gravity[1]), float(gravity[2])
    props.SetValue(KM.VOLUME_ACCELERATION, vec)
