# -*- coding: utf-8 -*-
"""ASCII ACIS ``.sat`` reader producing OpenCASCADE B-Rep solids.

Supports ACIS version-700 (and 400+ unindexed/indexed) ASCII ``.sat`` files
composed of analytical surfaces (``plane-surface``, ``cone-surface`` cylinder
and conical sheets, ``torus-surface`` including negative major/minor radii via
finite-UV ``Geom_SurfaceOfRevolution``) and B-spline/NURBS surfaces & curves
(``spline-surface``, ``intcurve-curve`` including ``exactcur``, ``parcur``,
``surfintcur``, ``straight-curve``, ``ellipse-curve``, and shared ``{ ref N }``
subtype blocks).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np


def _ptr(token: str) -> int:
    """Parse an ACIS entity pointer token ``$N`` (returning ``-1`` for ``$-1``)."""
    if token.startswith("$"):
        return int(token[1:])
    return -1


def _split_acis_records(text: str) -> tuple[int, list[str]]:
    """Split an ASCII ``.sat`` file into its ACIS version and ``#``-terminated records."""
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError("Invalid ASCII ACIS .sat file: fewer than 4 header lines")
    ver_tokens = lines[0].split()
    version = int(ver_tokens[0])

    recs: list[str] = []
    cur: list[str] = []
    for ln in lines[3:]:
        cur.append(ln)
        if "#" in ln:
            rec = "\n".join(cur).strip()
            cur = []
            if not rec:
                continue
            first = rec.split(None, 1)[0]
            # Strip optional '-<id>' record prefix used in older ACIS exports.
            if first.startswith("-") and first[1:].isdigit():
                rec = rec.split(None, 1)[1].strip()
            if rec.startswith("End-of-ACIS"):
                break
            recs.append(rec)
    return version, recs


def _parse_spline_block(
    toks: list[str],
    start_idx: int,
    subtype_tbl: list[Any],
) -> tuple[Any, int]:
    """Recursively parse a braced ACIS subtype block ``{ ... }`` in pre-order."""
    kind = toks[start_idx + 1]
    if kind == "ref":
        ref_id = int(toks[start_idx + 2])
        return subtype_tbl[ref_id], start_idx + 4

    ref_id = len(subtype_tbl)
    subtype_tbl.append(None)  # reserve pre-order slot before descending nested blocks

    depth = 1
    j = start_idx + 1
    inner: list[Any] = []
    while j < len(toks) and depth > 0:
        if toks[j] == "{":
            sub_val, j = _parse_spline_block(toks, j, subtype_tbl)
            inner.append(("SUB", sub_val))
            continue
        if toks[j] == "}":
            depth -= 1
            j += 1
            break
        inner.append(toks[j])
        j += 1

    res: Any = None
    if kind == "exactsur":
        mode = str(inner[2])
        du, dv = int(inner[3]), int(inner[4])
        idx = 5
        rational = mode == "nurbs"
        if rational and str(inner[idx]) in ("both", "u", "v", "none"):
            idx += 1
        idx += 4  # u_closure, v_closure, singleton flags
        nku, nkv = int(inner[idx]), int(inner[idx + 1])
        idx += 2
        uk: list[float] = []
        um: list[int] = []
        for _ in range(nku):
            uk.append(float(inner[idx]))
            um.append(int(inner[idx + 1]))
            idx += 2
        vk: list[float] = []
        vm: list[int] = []
        for _ in range(nkv):
            vk.append(float(inner[idx]))
            vm.append(int(inner[idx + 1]))
            idx += 2
        nu = sum(um) - du + 1
        nv = sum(vm) - dv + 1
        poles = np.zeros((nu, nv, 3), dtype=np.float64)
        weights = np.ones((nu, nv), dtype=np.float64)
        for iv in range(nv):
            for iu in range(nu):
                poles[iu, iv, 0] = float(inner[idx])
                poles[iu, iv, 1] = float(inner[idx + 1])
                poles[iu, iv, 2] = float(inner[idx + 2])
                idx += 3
                if rational:
                    weights[iu, iv] = float(inner[idx])
                    idx += 1
        res = ("exactsur", du, dv, uk, um, vk, vm, poles, weights, rational)
    elif kind in ("exactcur", "parcur", "surfintcur"):
        mode = str(inner[2])
        deg = int(inner[3])
        idx = 4
        rational = mode == "nurbs"
        if rational and str(inner[idx]) in ("both", "u", "v", "none"):
            idx += 1
        idx += 1  # open | periodic
        nk = int(inner[idx])
        idx += 1
        ck: list[float] = []
        cm: list[int] = []
        for _ in range(nk):
            ck.append(float(inner[idx]))
            cm.append(int(inner[idx + 1]))
            idx += 2
        nc = sum(cm) - deg + 1
        poles = np.zeros((nc, 3), dtype=np.float64)
        weights = np.ones(nc, dtype=np.float64)
        for ic in range(nc):
            poles[ic, 0] = float(inner[idx])
            poles[ic, 1] = float(inner[idx + 1])
            poles[ic, 2] = float(inner[idx + 2])
            idx += 3
            if rational:
                weights[ic] = float(inner[idx])
                idx += 1
        res = ("cur", deg, ck, cm, poles, weights, rational)

    subtype_tbl[ref_id] = res
    return res, j


def _make_occt_curve(
    recs: list[str],
    parsed_blocks: dict[int, Any],
    ci: int,
) -> tuple[Any, float]:
    """Convert an ACIS curve entity record to ``(Geom_Curve, param_scale)``."""
    from OCP.collections import Array1_double, Array1_gp_Pnt, Array1_int
    from OCP.Geom import Geom_BSplineCurve, Geom_Ellipse, Geom_Line
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    toks = recs[ci].split()
    typ = toks[0]
    if typ == "straight-curve":
        origin = np.array([float(toks[4]), float(toks[5]), float(toks[6])], dtype=np.float64)
        direction = np.array([float(toks[7]), float(toks[8]), float(toks[9])], dtype=np.float64)
        scale = float(np.linalg.norm(direction))
        unit_dir = direction / max(scale, 1e-30)
        return Geom_Line(gp_Pnt(*origin), gp_Dir(*unit_dir)), scale

    if typ == "ellipse-curve":
        center = np.array([float(toks[4]), float(toks[5]), float(toks[6])], dtype=np.float64)
        normal = np.array([float(toks[7]), float(toks[8]), float(toks[9])], dtype=np.float64)
        normal /= float(np.linalg.norm(normal))
        major = np.array([float(toks[10]), float(toks[11]), float(toks[12])], dtype=np.float64)
        r_maj = float(np.linalg.norm(major))
        r_min = r_maj * float(toks[13])
        return Geom_Ellipse(
            gp_Ax2(gp_Pnt(*center), gp_Dir(*normal), gp_Dir(*(major / r_maj))),
            r_maj,
            r_min,
        ), 1.0

    if typ == "intcurve-curve":
        sense = toks[4]
        _, deg, ck, cm, poles, weights, rational = parsed_blocks[ci]
        if sense == "reversed":
            ck = [-float(k) for k in reversed(ck)]
            cm = list(reversed(cm))
            poles = poles[::-1].copy()
            weights = weights[::-1].copy()
        nc = len(poles)
        cm2 = list(cm)
        cm2[0] = deg + 1
        cm2[-1] = deg + 1
        pa = Array1_gp_Pnt(1, nc)
        for i, p in enumerate(poles, 1):
            pa.SetValue(i, gp_Pnt(float(p[0]), float(p[1]), float(p[2])))
        cka = Array1_double(1, len(ck))
        cma = Array1_int(1, len(ck))
        for i, (k, m) in enumerate(zip(ck, cm2), 1):
            cka.SetValue(i, float(k))
            cma.SetValue(i, int(m))
        if rational:
            wa = Array1_double(1, nc)
            for i, w in enumerate(weights, 1):
                wa.SetValue(i, float(w))
            return Geom_BSplineCurve(pa, wa, cka, cma, deg), 1.0
        return Geom_BSplineCurve(pa, cka, cma, deg), 1.0

    raise ValueError(f"Unsupported ACIS curve record type: {typ!r}")


def _face_loop_points(recs: list[str], vpts: dict[int, Any], fi: int) -> np.ndarray:
    """Collect 3D vertex coordinates around all loops of ACIS face ``fi``."""
    ft = recs[fi].split()
    lp = _ptr(ft[5])
    pts: list[list[float]] = []
    while lp >= 0:
        lt = recs[lp].split()
        ce = _ptr(lt[5])
        seen: set[int] = set()
        while ce >= 0 and ce not in seen:
            seen.add(ce)
            cet = recs[ce].split()
            et = recs[_ptr(cet[7])].split()
            for vi in (_ptr(et[4]), _ptr(et[6])):
                if vi >= 0:
                    p = vpts[vi]
                    pts.append([float(p.X()), float(p.Y()), float(p.Z())])
            ce = _ptr(cet[4])
        lp = _ptr(lt[4])
    return np.asarray(pts, dtype=np.float64)


def _make_occt_surface(
    recs: list[str],
    parsed_blocks: dict[int, Any],
    vpts: dict[int, Any],
    fi: int,
) -> tuple[Any, str, tuple[float, float] | None]:
    """Build ``(Geom_Surface, surf_type, optional_torus_v_range)`` for face ``fi``."""
    from OCP.collections import Array1_double, Array1_int, Array2_double, Array2_gp_Pnt
    from OCP.Geom import (
        Geom_BSplineSurface,
        Geom_Circle,
        Geom_ConicalSurface,
        Geom_CylindricalSurface,
        Geom_Plane,
        Geom_SurfaceOfRevolution,
    )
    from OCP.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pnt

    si = _ptr(recs[fi].split()[8])
    toks = recs[si].split()
    typ = toks[0]

    if typ == "plane-surface":
        origin = np.array([float(toks[4]), float(toks[5]), float(toks[6])], dtype=np.float64)
        normal = np.array([float(toks[7]), float(toks[8]), float(toks[9])], dtype=np.float64)
        normal /= float(np.linalg.norm(normal))
        u_dir = np.array([float(toks[10]), float(toks[11]), float(toks[12])], dtype=np.float64)
        u_dir /= float(np.linalg.norm(u_dir))
        return Geom_Plane(gp_Ax3(gp_Pnt(*origin), gp_Dir(*normal), gp_Dir(*u_dir))), typ, None

    if typ == "cone-surface":
        origin = np.array([float(toks[4]), float(toks[5]), float(toks[6])], dtype=np.float64)
        axis = np.array([float(toks[7]), float(toks[8]), float(toks[9])], dtype=np.float64)
        axis /= float(np.linalg.norm(axis))
        u_dir = np.array([float(toks[10]), float(toks[11]), float(toks[12])], dtype=np.float64)
        u_dir /= float(np.linalg.norm(u_dir))
        sin_a = float(toks[16])
        cos_a = float(toks[17])
        radius = float(toks[18])
        ax3 = gp_Ax3(gp_Pnt(*origin), gp_Dir(*axis), gp_Dir(*u_dir))
        if abs(sin_a) < 1e-12:
            return Geom_CylindricalSurface(ax3, radius), typ, None
        semi_angle = math.atan2(sin_a, -cos_a)
        return Geom_ConicalSurface(ax3, semi_angle, radius), typ, None

    if typ == "torus-surface":
        center = np.array([float(toks[4]), float(toks[5]), float(toks[6])], dtype=np.float64)
        axis = np.array([float(toks[7]), float(toks[8]), float(toks[9])], dtype=np.float64)
        axis /= float(np.linalg.norm(axis))
        r_maj, r_min = float(toks[10]), float(toks[11])
        u_dir = np.array([float(toks[12]), float(toks[13]), float(toks[14])], dtype=np.float64)
        u_dir /= float(np.linalg.norm(u_dir))
        v_dir = np.cross(axis, u_dir)
        circle_center = center + r_maj * u_dir
        circ = Geom_Circle(
            gp_Ax2(gp_Pnt(*circle_center), gp_Dir(*v_dir), gp_Dir(*u_dir)),
            abs(r_min),
        )
        surf = Geom_SurfaceOfRevolution(circ, gp_Ax1(gp_Pnt(*center), gp_Dir(*axis)))
        pts = _face_loop_points(recs, vpts, fi)
        d = pts - center
        z = d @ axis
        rad = np.linalg.norm(d - np.outer(z, axis), axis=1)
        if rad.size and rad.min() < 1e-4:
            th = np.unwrap(np.arctan2(-z, rad - r_maj))
            return surf, typ, (float(th.min()), float(th.max()))
        return surf, typ, None

    if typ == "spline-surface":
        _, du, dv, uk, um, vk, vm, poles, weights, rational = parsed_blocks[si]
        nu, nv = poles.shape[:2]
        um2 = list(um)
        um2[0] = du + 1
        um2[-1] = du + 1
        vm2 = list(vm)
        vm2[0] = dv + 1
        vm2[-1] = dv + 1
        pa = Array2_gp_Pnt(1, nu, 1, nv)
        for iu in range(nu):
            for iv in range(nv):
                p = poles[iu, iv]
                pa.SetValue(iu + 1, iv + 1, gp_Pnt(float(p[0]), float(p[1]), float(p[2])))
        uka = Array1_double(1, len(uk))
        uma = Array1_int(1, len(uk))
        for i, (k, m) in enumerate(zip(uk, um2), 1):
            uka.SetValue(i, float(k))
            uma.SetValue(i, int(m))
        vka = Array1_double(1, len(vk))
        vma = Array1_int(1, len(vk))
        for i, (k, m) in enumerate(zip(vk, vm2), 1):
            vka.SetValue(i, float(k))
            vma.SetValue(i, int(m))
        if rational:
            wa = Array2_double(1, nu, 1, nv)
            for iu in range(nu):
                for iv in range(nv):
                    wa.SetValue(iu + 1, iv + 1, float(weights[iu, iv]))
            return Geom_BSplineSurface(pa, wa, uka, vka, uma, vma, du, dv), typ, None
        return Geom_BSplineSurface(pa, uka, vka, uma, vma, du, dv), typ, None

    raise ValueError(f"Unsupported ACIS surface record type: {typ!r}")


def _apply_acis_transform(shape: Any, rec: str) -> Any:
    """Apply an ACIS ``transform`` record if it is non-identity."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Trsf

    toks = rec.split()
    vals = [float(x) for x in toks[3:16]]  # 3x3 linear (9) + translation (3) + scale (1)
    mat = np.array(vals[:9], dtype=np.float64).reshape(3, 3)
    trans = np.array(vals[9:12], dtype=np.float64)
    scale = float(vals[12])
    # ACIS stores row-vector M^T -> column vector transform uses mat.T * scale
    linear = mat.T * scale
    if np.allclose(linear, np.eye(3), atol=1e-12) and np.allclose(trans, 0.0, atol=1e-12):
        return shape
    tr = gp_Trsf()
    tr.SetValues(
        linear[0, 0], linear[0, 1], linear[0, 2], trans[0],
        linear[1, 0], linear[1, 1], linear[1, 2], trans[1],
        linear[2, 0], linear[2, 1], linear[2, 2], trans[2],
    )
    return BRepBuilderAPI_Transform(shape, tr, True).Shape()


def read_sat_shape(path: str | Path, sewing_tol: float = 5e-3) -> Any:
    """Read an ASCII ACIS ``.sat`` file and return an OpenCASCADE ``TopoDS_Shape``.

    Parameters
    ----------
    path:
        Path to an ASCII ACIS ``.sat`` file.
    sewing_tol:
        Sewing tolerance passed to ``BRepBuilderAPI_Sewing`` when stitching
        the reconstructed trimmed B-Rep faces of each ACIS body into a solid.

    Returns
    -------
    TopoDS_Shape
        A ``TopoDS_Solid`` when the file contains 1 body or a ``TopoDS_Compound``
        of reconstructed solids otherwise.
    """
    from OCP.BRep import BRep_Builder, BRep_Tool
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeSolid,
        BRepBuilderAPI_MakeVertex,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.gp import gp_Pnt
    from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Solid
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS, TopoDS_Compound

    from wht_iga.geometry.cad import bbox_diag

    src = Path(path)
    text = src.read_text(encoding="utf-8", errors="replace")
    _, recs = _split_acis_records(text)

    subtype_tbl: list[Any] = []
    parsed_blocks: dict[int, Any] = {}
    for ri, rec in enumerate(recs):
        if "{" in rec:
            toks = rec.replace("{", " { ").replace("}", " } ").split()
            val, _ = _parse_spline_block(toks, toks.index("{"), subtype_tbl)
            parsed_blocks[ri] = val

    builder = BRep_Builder()
    verts: dict[int, Any] = {}
    vpts: dict[int, Any] = {}
    for ri, rec in enumerate(recs):
        toks = rec.split()
        if toks[0] in ("vertex", "tvertex-vertex"):
            pt_toks = recs[_ptr(toks[5])].split()
            p = gp_Pnt(float(pt_toks[4]), float(pt_toks[5]), float(pt_toks[6]))
            v = BRepBuilderAPI_MakeVertex(p).Vertex()
            tol = float(toks[6]) if toks[0] == "tvertex-vertex" else 1e-6
            builder.UpdateVertex(v, max(tol, 1e-4))
            verts[ri] = v
            vpts[ri] = p

    edge_specs: list[tuple[int, Any]] = []
    for ri, rec in enumerate(recs):
        et = rec.split()
        if et[0] not in ("edge", "tedge-edge"):
            continue
        v0i, t0, v1i, t1, cp, sense = (
            _ptr(et[4]),
            float(et[5]),
            _ptr(et[6]),
            float(et[7]),
            _ptr(et[9]),
            et[10],
        )
        if cp < 0:
            edge_specs.append((ri, None))
            continue
        crv, scale = _make_occt_curve(recs, parsed_blocks, cp)
        u0, u1 = t0 * scale, t1 * scale
        d0 = crv.Value(u0).Distance(vpts[v0i])
        d1 = crv.Value(u1).Distance(vpts[v1i])
        builder.UpdateVertex(verts[v0i], max(BRep_Tool.Tolerance_s(verts[v0i]), d0 + 1e-4))
        builder.UpdateVertex(verts[v1i], max(BRep_Tool.Tolerance_s(verts[v1i]), d1 + 1e-4))
        etol = float(et[13]) if et[0] == "tedge-edge" else 1e-5
        edge_specs.append((ri, (crv, verts[v0i], verts[v1i], u0, u1, sense, etol)))

    edges: dict[int, Any] = {}
    for ri, spec in edge_specs:
        if spec is None:
            edges[ri] = None
            continue
        crv, va, vb, u0, u1, sense, etol = spec
        me = BRepBuilderAPI_MakeEdge(crv, va, vb, u0, u1)
        if not me.IsDone():
            edges[ri] = None
            continue
        e = me.Edge()
        builder.UpdateEdge(e, max(etol, 1e-4))
        edges[ri] = TopoDS.Edge(e.Reversed()) if sense == "reversed" else e

    solids: list[Any] = []
    for bi, rec in enumerate(recs):
        bt = rec.split()
        if bt[0] != "body":
            continue
        lump_id = _ptr(bt[4])
        while lump_id >= 0:
            lump_toks = recs[lump_id].split()
            sh_id = _ptr(lump_toks[5])
            sew = BRepBuilderAPI_Sewing(float(sewing_tol))
            while sh_id >= 0:
                sh_toks = recs[sh_id].split()
                cur_f = _ptr(sh_toks[6])
                while cur_f >= 0:
                    ft = recs[cur_f].split()
                    surf, _, v_cap = _make_occt_surface(recs, parsed_blocks, vpts, cur_f)
                    fsense = ft[9]
                    if v_cap is not None:
                        face = BRepBuilderAPI_MakeFace(
                            surf, 0.0, 2.0 * math.pi, v_cap[0], v_cap[1], 1e-6
                        ).Face()
                        if fsense == "reversed":
                            face = TopoDS.Face(face.Reversed())
                        sf = ShapeFix_Face(face)
                        sf.Perform()
                        sew.Add(sf.Face())
                    else:
                        wires: list[Any] = []
                        lp = _ptr(ft[5])
                        while lp >= 0:
                            lt = recs[lp].split()
                            ce = _ptr(lt[5])
                            seen: set[int] = set()
                            mw = BRepBuilderAPI_MakeWire()
                            cnt = 0
                            while ce >= 0 and ce not in seen:
                                seen.add(ce)
                                cet = recs[ce].split()
                                e = edges.get(_ptr(cet[7]))
                                if e is not None:
                                    e_use = TopoDS.Edge(e.Reversed()) if cet[8] == "reversed" else e
                                    mw.Add(e_use)
                                    cnt += 1
                                ce = _ptr(cet[4])
                            if cnt > 0 and mw.IsDone():
                                wires.append(mw.Wire())
                            lp = _ptr(lt[4])
                        if wires:
                            if len(wires) > 1:
                                wires.sort(key=lambda w: bbox_diag(w), reverse=True)
                            mf = BRepBuilderAPI_MakeFace(surf, wires[0], True)
                            for hw in wires[1:]:
                                mf.Add(hw)
                            if mf.IsDone():
                                face = mf.Face()
                                sf = ShapeFix_Face(face)
                                sf.FixOrientationMode = 1
                                sf.Perform()
                                face = sf.Face()
                                if fsense == "reversed":
                                    face = TopoDS.Face(face.Reversed())
                                sew.Add(face)
                    cur_f = _ptr(ft[4])
                sh_id = _ptr(sh_toks[4])

            sew.Perform()
            sewn = sew.SewedShape()
            ex = TopExp_Explorer(sewn, TopAbs_SHELL)
            body_solids: list[Any] = []
            while ex.More():
                sh = TopoDS.Shell(ex.Current())
                ms = BRepBuilderAPI_MakeSolid(sh)
                if ms.IsDone():
                    sfs = ShapeFix_Solid(ms.Solid())
                    sfs.Perform()
                    solid = TopoDS.Solid(sfs.Solid())
                    props = GProp_GProps()
                    BRepGProp.VolumeProperties_s(solid, props)
                    if float(props.Mass()) < 0.0:
                        solid = TopoDS.Solid(solid.Reversed())
                    body_solids.append(solid)
                ex.Next()

            tr_id = _ptr(bt[6])
            for s in body_solids:
                if tr_id >= 0:
                    s = _apply_acis_transform(s, recs[tr_id])
                solids.append(s)
            lump_id = _ptr(lump_toks[4])

    if not solids:
        raise IOError(f"SAT read produced no valid solids: {src}")
    if len(solids) == 1:
        return solids[0]
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for s in solids:
        builder.Add(comp, s)
    return comp


__all__ = ["read_sat_shape"]
