# -*- coding: utf-8 -*-
"""P5 / S43: E05 (door-lock-assembly) B-Rep reconstruction & reduced 3-part contact precursor.

References:
- ``WHT_IGA/research-cf-IGA/user-guide/example-problems/door-lock-assembly/door-lock-assembly.md``
- ``WHT_IGA/tests/verification/data/e05.json``

Tests:
1. Documented contract check (MMTS units, weak/strong steel, 3 contact pairs, RF3 = -46.946 kN).
2. Parametric 3-part B-Rep reconstruction in OCP (hinge, link, and 3-zone partitioned pin:
   strong-steel cylindrical shaft + 2 weak-steel end caps), verifying ``watertight=True``,
   ``free_edges=0``, and positive B-Rep/tessellated volumes.
3. Reduced 3-part deformable IGA contact load-transfer precursor (coarse ``p=2`` BoxSplines in
   Kratos ``SmallDisplacementElement3D4N``, coupled via ``wht_iga.analyses.contact`` SSN penalty
   for ``mu=0.0`` and ``mu=0.3``): verifies load transfer ``hinge -> pin -> link``, opposing
   ``RF3 < 0`` on ``link``, exact global force balance ``sum_hinge_RF + sum_link_RF == 0``,
   and measurable Coulomb friction response.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pytest

from wht_iga.analyses.contact import (
    DeformableBoxFaceMaster,
    _solve_linear_step,
    expand_vector_basis,
    sample_box_face_gauss,
)
from wht_iga.geometry.cad import brep_area, brep_volume, tessellate, validate_shape
from wht_iga.geometry.occ_io import _count_free_edges
from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import equation_permutation
from wht_iga.kratos_backend.geometry import cell_gauss_points, make_nurbs_volume
from wht_iga.reference.box import BoxSpline

pytestmark = [pytest.mark.kratos]

_WHT_ROOT = Path(__file__).resolve().parents[2]
_E05_JSON = _WHT_ROOT / "tests" / "verification" / "data" / "e05.json"
_E05_DOC = (
    _WHT_ROOT
    / "research-cf-IGA"
    / "user-guide"
    / "example-problems"
    / "door-lock-assembly"
    / "door-lock-assembly.md"
)


# --------------------------------------------------------------------------------------
# Test 1: Documented model & verification contract
# --------------------------------------------------------------------------------------


def test_e05_doc_and_json_contract() -> None:
    """Verify E05 JSON and source doc contain all documented parameters and RF3 target."""
    data = json.loads(_E05_JSON.read_text(encoding="utf-8"))
    doc_text = _E05_DOC.read_text(encoding="utf-8")

    assert data["id"] == "E05"
    assert data["geometry"]["units"] == "MMTS (mm, N, MPa)"
    assert data["materials"]["weak_steel"]["E_MPa"] == 200000.0
    assert data["materials"]["weak_steel"]["nu"] == 0.29
    assert data["materials"]["strong_steel"]["E_MPa"] == 200000.0
    assert data["materials"]["strong_steel"]["nu"] == 0.29
    assert data["boundary_conditions"]["hinge"] == "prescribed displacement (U1, U2, U3) = (0, 0, 5) mm"
    assert data["boundary_conditions"]["link"] == "pinned"

    pairs = {p["pair"]: p["main"] for p in data["boundary_conditions"]["contact_pairs"]}
    assert pairs == {
        "pin_to_hinge_contact": "Pin",
        "pin_to_link_contact": "Pin",
        "hinge_to_link_contact": "Hinge",
    }
    assert data["reference_values"]["summed_link_RF3_kN"] == pytest.approx(-46.946)
    assert data["acceptance"]["rel_tol"] == pytest.approx(0.05)

    for token in ("pin_to_hinge_contact", "pin_to_link_contact", "hinge_to_link_contact", "-46.946"):
        assert token in doc_text


# --------------------------------------------------------------------------------------
# Test 2: Parametric 3-part B-Rep reconstruction in OCP
# --------------------------------------------------------------------------------------


def _make_cylinder(x0: float, y0: float, z0: float, dx: float, r: float):
    """Create an X-aligned cylinder starting at ``(x0, y0, z0)`` with length ``dx`` and radius ``r``."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    ax = gp_Ax2(gp_Pnt(x0, y0, z0), gp_Dir(1.0, 0.0, 0.0))
    return BRepPrimAPI_MakeCylinder(ax, float(r), float(dx)).Shape()


def _make_z_cylinder(x0: float, y0: float, z0: float, dz: float, r: float):
    """Create a Z-aligned cylinder starting at ``(x0, y0, z0)`` with height ``dz`` and radius ``r``."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    ax = gp_Ax2(gp_Pnt(x0, y0, z0), gp_Dir(0.0, 0.0, 1.0))
    return BRepPrimAPI_MakeCylinder(ax, float(r), float(dz)).Shape()


def _make_y_cylinder(x0: float, y0: float, z0: float, dy: float, r: float):
    """Create a Y-aligned cylinder starting at ``(x0, y0, z0)`` with length ``dy`` and radius ``r``."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    ax = gp_Ax2(gp_Pnt(x0, y0, z0), gp_Dir(0.0, 1.0, 0.0))
    return BRepPrimAPI_MakeCylinder(ax, float(r), float(dy)).Shape()


def _make_box(xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float):
    """Create an axis-aligned B-Rep box."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(xmin, ymin, zmin), gp_Pnt(xmax, ymax, zmax)).Shape()


def _fuse(a, b):
    """Boolean union of two B-Rep shapes."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

    op = BRepAlgoAPI_Fuse(a, b)
    op.Build()
    assert op.IsDone(), "BRepAlgoAPI_Fuse failed"
    return op.Shape()


def _cut(a, b):
    """Boolean difference ``a - b`` of two B-Rep shapes."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    op = BRepAlgoAPI_Cut(a, b)
    op.Build()
    assert op.IsDone(), "BRepAlgoAPI_Cut failed"
    return op.Shape()


def _split_closed_faces(shape):
    """Split periodic closed faces (cylinders) into two half-faces so seam edges are 2-manifold."""
    from OCP.ShapeUpgrade import ShapeUpgrade_ShapeDivideClosed

    up = ShapeUpgrade_ShapeDivideClosed(shape)
    up.SetNbSplitPoints(1)
    up.Perform()
    return up.Result()


def build_door_lock_brep_parts() -> dict[str, object]:
    """Build parametric watertight B-Rep solids for the E05 door-lock assembly (mm).

    Pin axis is along X, centered at ``(y=0, z=0)``:
    - ``pin_cyl``: strong-steel shaft, ``r=4.0 mm``, ``x in [-16, 16]``
    - ``pin_cap_l`` / ``pin_cap_r``: weak-steel end caps, ``r=5.5 mm``, ``x in [-19, -16]`` & ``[16, 19]``
    - ``hinge``: vertical backplate (2 mounting holes) + 2 clevis ears with ``r=4.0 mm`` pin bores
    - ``link``: U-shaped clevis saddle with ``r=4.0 mm`` pin bore + 2 pinned mounting feet
    """
    r_pin = 4.0
    r_cap = 5.5

    # 1. Connecting pin (3 partitioned solids per doc §1.4)
    pin_cyl = _make_cylinder(-16.0, 0.0, 0.0, 32.0, r_pin)
    pin_cap_l = _make_cylinder(-19.0, 0.0, 0.0, 3.0, r_cap)
    pin_cap_r = _make_cylinder(16.0, 0.0, 0.0, 3.0, r_cap)

    # 2. Hinge: backplate + two clevis ears + rounded noses - pin bore - mounting holes
    backplate = _make_box(-24.0, 24.0, -24.0, -18.0, -8.0, 36.0)
    ear_l_arm = _make_box(-16.0, -11.0, -20.0, 0.0, -7.8, 7.8)
    ear_l_nose = _make_cylinder(-16.0, 0.0, 0.0, 5.0, 8.0)
    ear_r_arm = _make_box(11.0, 16.0, -20.0, 0.0, -7.8, 7.8)
    ear_r_nose = _make_cylinder(11.0, 0.0, 0.0, 5.0, 8.0)

    hinge = _fuse(backplate, ear_l_arm)
    hinge = _fuse(hinge, ear_l_nose)
    hinge = _fuse(hinge, ear_r_arm)
    hinge = _fuse(hinge, ear_r_nose)
    # Coaxial pin bore through both ears + 2 mounting holes in backplate
    hinge = _cut(hinge, _make_cylinder(-17.0, 0.0, 0.0, 34.0, r_pin))
    hinge = _cut(hinge, _make_y_cylinder(-12.0, -25.0, 22.0, 8.0, 4.5))
    hinge = _cut(hinge, _make_y_cylinder(12.0, -25.0, 22.0, 8.0, 4.5))

    # 3. Link: U-shaped clevis saddle around pin + two mounting feet with holes
    saddle_block = _make_box(-10.0, 10.0, -7.8, 7.8, -8.0, 0.0)
    saddle_top = _make_cylinder(-10.0, 0.0, 0.0, 20.0, 8.0)
    foot_l = _make_box(-10.0, -4.0, 6.0, 26.0, -8.0, -4.0)
    foot_r = _make_box(4.0, 10.0, 6.0, 26.0, -8.0, -4.0)
    link = _fuse(saddle_block, saddle_top)
    link = _fuse(link, foot_l)
    link = _fuse(link, foot_r)
    link = _cut(link, _make_cylinder(-11.0, 0.0, 0.0, 22.0, r_pin))
    link = _cut(link, _make_z_cylinder(-7.0, 18.0, -9.0, 6.0, 2.5))
    link = _cut(link, _make_z_cylinder(7.0, 18.0, -9.0, 6.0, 2.5))

    return {
        "hinge": _split_closed_faces(hinge),
        "link": _split_closed_faces(link),
        "pin_cyl": _split_closed_faces(pin_cyl),
        "pin_cap_l": _split_closed_faces(pin_cap_l),
        "pin_cap_r": _split_closed_faces(pin_cap_r),
    }


def test_e05_brep_reconstruction_watertight() -> None:
    """Verify all 5 reconstructed B-Rep solids are watertight with 0 free edges and positive volume."""
    parts = build_door_lock_brep_parts()
    print("\n[E05 B-Rep Validation]")
    for name, shape in parts.items():
        n_free = _count_free_edges(shape)
        tess = tessellate(shape, deflection=0.2)
        rep = validate_shape(shape, tess)
        v_brep = brep_volume(shape)
        a_brep = brep_area(shape)
        print(
            f"  {name:10s}: V_brep={v_brep:10.3f} mm^3, A_brep={a_brep:10.3f} mm^2, "
            f"tris={len(tess.tris):5d}, free_edges={n_free}, watertight={rep.watertight}"
        )
        assert n_free == 0, f"{name} has {n_free} free edges"
        assert rep.watertight, f"{name} is not watertight"
        assert v_brep > 0.0 and rep.tess_volume > 0.0
        assert rep.area_err < 0.05


# --------------------------------------------------------------------------------------
# Test 3: Reduced 3-part deformable IGA contact load-transfer precursor
# --------------------------------------------------------------------------------------


def _assemble_kratos_body(
    KM,
    SM,
    box: BoxSpline,
    E: float,
    nu: float,
    *,
    zero_fn=None,
    presc_fn=None,
):
    """Assemble unconstrained Kratos stiffness ``K_full`` and partition DOFs for ``box``."""
    import KratosMultiphysics.scipy_conversion_tools as sct

    model = KM.Model()
    mp = model.CreateModelPart("DoorLockBody")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    vol, _ = make_nurbs_volume(KM, mp, box)
    props = mp.CreateNewProperties(1)
    props.SetValue(KM.YOUNG_MODULUS, float(E))
    props.SetValue(KM.POISSON_RATIO, float(nu))
    props.SetValue(KM.CONSTITUTIVE_LAW, SM.LinearElastic3DLaw())

    nx, ny, nz = (int(v) for v in box.n_cells)
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                pts = cell_gauss_points(box, i, j, k)
                qps = KM.GeometriesVector()
                vol.CreateQuadraturePointGeometries(qps, 2, pts)
                mp.CreateNewElement("SmallDisplacementElement3D4N", eid, qps[0], props)
                eid += 1

    for v, r in (
        (KM.DISPLACEMENT_X, KM.REACTION_X),
        (KM.DISPLACEMENT_Y, KM.REACTION_Y),
        (KM.DISPLACEMENT_Z, KM.REACTION_Z),
    ):
        KM.VariableUtils().AddDof(v, r, mp)

    builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    strat = KM.ResidualBasedLinearStrategy(mp, scheme, builder, False, False, False, False)
    strat.SetEchoLevel(0)
    strat.Initialize()
    mp.CloneTimeStep(1.0)
    strat.InitializeSolutionStep()

    A = strat.GetSystemMatrix()
    b = strat.GetSystemVector()
    space = KM.UblasSparseSpace()
    space.SetToZeroMatrix(A)
    space.SetToZeroVector(b)
    builder.BuildLHS(scheme, mp, A)

    perm = equation_permutation(KM, mp, box)
    inv = np.argsort(perm)
    K_full = sct.to_csr(A)[inv, :][:, inv].tocsr()

    Xcp = box.control_points()
    zero = np.zeros(box.n_dofs, dtype=bool)
    presc = np.zeros(box.n_dofs, dtype=bool)
    if zero_fn is not None:
        for I in range(box.n_ctrl):
            for d in zero_fn(Xcp[I]):
                zero[3 * I + d] = True
    if presc_fn is not None:
        for I in range(box.n_ctrl):
            for d in presc_fn(Xcp[I]):
                if not zero[3 * I + d]:
                    presc[3 * I + d] = True
    free = ~(zero | presc)

    return {
        "K_full": K_full,
        "Kf": K_full[free, :][:, free].tocsr(),
        "free": free,
        "zero": zero,
        "presc": presc,
        "n_dofs": box.n_dofs,
    }


def _curved_crown_box(
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    n_cells: tuple[int, int, int],
    z_bot_fn,
    z_top_fn,
) -> BoxSpline:
    """Create a ``p=2`` BoxSpline whose Z control points are lofted between ``z_bot_fn(y)`` and ``z_top_fn(y)``."""
    box = BoxSpline(lo=lo, hi=hi, n_cells=n_cells, degree=2)
    cp = box.control_points().copy()
    z0, z1 = float(lo[2]), float(hi[2])
    for I in range(box.n_ctrl):
        y = float(cp[I, 1])
        t = (float(cp[I, 2]) - z0) / (z1 - z0)
        cp[I, 2] = (1.0 - t) * float(z_bot_fn(y)) + t * float(z_top_fn(y))
    box.cp = cp
    return box


def _solve_3body_ssn_contact(
    bodies: list[dict],
    u_fixed_list: list[np.ndarray],
    pairs: list[dict],
    eps_n: float,
    *,
    mu: float = 0.0,
    eps_t: float | None = None,
    max_iter: int = 60,
    rtol: float = 1e-7,
) -> dict:
    """Coupled 3-body Gauss-weighted SSN contact solver with Coulomb friction.

    Each pair in ``pairs`` specifies:
    - ``name``: contact pair identifier
    - ``m_idx``, ``s_idx``: body indices in ``bodies``
    - ``master_face``: ``DeformableBoxFaceMaster``
    - ``X_s``, ``w_q``, ``B_s``: slave surface Gauss points, weights, and scalar basis
    """
    import scipy.sparse as sp

    n_bodies = len(bodies)
    n_free = [int(b["free"].sum()) for b in bodies]
    offsets = np.cumsum([0] + n_free)
    n_tot = int(offsets[-1])

    K_diag = sp.block_diag([b["Kf"] for b in bodies], format="csr")
    f_ext = []
    u_base = []
    for b, u_fix in zip(bodies, u_fixed_list):
        u_b = np.zeros(b["n_dofs"], dtype=np.float64)
        u_b[~b["free"]] = u_fix[~b["free"]]
        u_base.append(u_b)
        # -(K_fp @ u_presc)
        f_ext.append(-(b["K_full"][b["free"], :][:, ~b["free"]] @ u_b[~b["free"]]))
    f_ext_vec = np.concatenate(f_ext)

    pair_prep = []
    for p in pairs:
        b_s = bodies[p["s_idx"]]
        B_s_vec = expand_vector_basis(p["B_s"], ndim=3)
        B_s_free = B_s_vec[:, b_s["free"]].tocsr()
        u_sb = u_base[p["s_idx"]]
        M = len(p["X_s"])
        x_s_base = p["X_s"] + (B_s_vec[:, ~b_s["free"]] @ u_sb[~b_s["free"]]).reshape(M, 3)
        pair_prep.append({
            **p,
            "M": M,
            "B_s_free": B_s_free,
            "x_s_base": x_s_base,
            "x_s_prev": p["X_s"].copy(),
            "U_m_prev": np.zeros((p["master_face"].n_ctrl, 3), dtype=np.float64),
        })

    if mu > 0.0 and eps_t is None:
        eps_t = float(eps_n)

    I3 = np.eye(3, dtype=np.float64)
    u = np.zeros(n_tot, dtype=np.float64)

    def _unpack_full(b_idx: int, u_vec: np.ndarray) -> np.ndarray:
        b = bodies[b_idx]
        full = u_base[b_idx].copy()
        full[b["free"]] = u_vec[offsets[b_idx] : offsets[b_idx + 1]]
        return full

    def _eval(u_vec: np.ndarray, initial_touch: bool = False):
        f_c_blocks = [np.zeros(nf, dtype=np.float64) for nf in n_free]
        K_blocks = [[None] * n_bodies for _ in range(n_bodies)]
        pair_stats = {}

        for pp in pair_prep:
            mi, si = pp["m_idx"], pp["s_idx"]
            M = pp["M"]
            um_f = u_vec[offsets[mi] : offsets[mi + 1]]
            us_f = u_vec[offsets[si] : offsets[si + 1]]
            U_m_curr = _unpack_full(mi, u_vec).reshape(-1, 3)
            x_s = pp["x_s_base"] + (pp["B_s_free"] @ us_f).reshape(M, 3)

            g_n, n_m, _x_m, _ab_m, B_m = pp["master_face"].evaluate(x_s, U_m_curr)
            B_m_vec = expand_vector_basis(B_m, ndim=3)
            B_m_free = B_m_vec[:, bodies[mi]["free"]].tocsr()

            active = (g_n <= 1e-9) if initial_touch else (g_n < 0.0)
            t_n = np.zeros((M, 3), dtype=np.float64)
            t_t = np.zeros((M, 3), dtype=np.float64)
            D_loc = np.zeros((M, 3, 3), dtype=np.float64)

            if np.any(active):
                t_n[active] = (eps_n * g_n[active])[:, None] * n_m[active]
                nn = n_m[active, :, None] * n_m[active, None, :]
                D_loc[active] = eps_n * nn

                if mu > 0.0:
                    du_m = B_m @ (U_m_curr - pp["U_m_prev"])
                    du_rel = (x_s - pp["x_s_prev"]) - du_m
                    for q in np.flatnonzero(active):
                        n_q = n_m[q]
                        P_T = I3 - np.outer(n_q, n_q)
                        du_t = P_T @ du_rel[q]
                        t_tr = float(eps_t) * du_t
                        norm_tr = float(np.linalg.norm(t_tr))
                        cap = mu * abs(eps_n * float(g_n[q]))
                        if norm_tr <= cap or norm_tr < 1e-30:
                            t_t[q] = t_tr
                            D_loc[q] += float(eps_t) * P_T
                        else:
                            t_t[q] = (cap / norm_tr) * t_tr
                            d_dir = t_tr / norm_tr
                            D_loc[q] += (
                                (cap / norm_tr) * (P_T - np.outer(d_dir, d_dir))
                                - mu * eps_n * np.outer(d_dir, n_q)
                            )

            w_q = pp["w_q"]
            t_res = (w_q[:, None] * (t_n + t_t)).ravel()
            f_c_blocks[mi] -= pp["B_m_free" if "B_m_free" in pp else None] if False else (B_m_free.T @ t_res)
            f_c_blocks[si] += pp["B_s_free"].T @ t_res

            w_D = w_q[:, None, None] * D_loc
            rows = np.repeat(np.arange(3 * M, dtype=np.int64), 3)
            cols = (np.arange(M, dtype=np.int64)[:, None, None] * 3 + np.arange(3)[None, None, :])
            cols = np.broadcast_to(cols, (M, 3, 3)).ravel()
            D_sp = sp.csr_matrix((w_D.ravel(), (rows, cols)), shape=(3 * M, 3 * M))

            K_mm = (B_m_free.T @ D_sp @ B_m_free).tocsr()
            K_ms = -(B_m_free.T @ D_sp @ pp["B_s_free"]).tocsr()
            K_sm = -(pp["B_s_free"].T @ D_sp @ B_m_free).tocsr()
            K_ss = (pp["B_s_free"].T @ D_sp @ pp["B_s_free"]).tocsr()

            for r_i, c_i, blk in ((mi, mi, K_mm), (mi, si, K_ms), (si, mi, K_sm), (si, si, K_ss)):
                K_blocks[r_i][c_i] = blk if K_blocks[r_i][c_i] is None else (K_blocks[r_i][c_i] + blk)

            pair_stats[pp["name"]] = {
                "n_active": int(active.sum()),
                "F_total": (w_q[:, None] * (t_n + t_t)).sum(axis=0),
                "F_t": (w_q[:, None] * t_t).sum(axis=0),
            }

        for r_i in range(n_bodies):
            for c_i in range(n_bodies):
                if K_blocks[r_i][c_i] is None:
                    K_blocks[r_i][c_i] = sp.csr_matrix((n_free[r_i], n_free[c_i]), dtype=np.float64)

        K_c = sp.bmat(K_blocks, format="csr")
        f_c = np.concatenate(f_c_blocks)
        res = K_diag @ u_vec + f_c - f_ext_vec
        return res, (K_diag + K_c).tocsr(), pair_stats

    res, K_T, pair_stats = _eval(u, initial_touch=True)
    r0 = float(np.linalg.norm(res))
    scale = max(r0, float(np.linalg.norm(f_ext_vec)), 1.0)
    converged = False
    for it in range(max_iter):
        if float(np.linalg.norm(res)) / scale < rtol and it > 0:
            converged = True
            break
        du = _solve_linear_step(K_T, -res, use_pardiso=False)
        alpha = 1.0
        cur_norm = float(np.linalg.norm(res))
        while alpha > 1e-4:
            u_try = u + alpha * du
            res_try, _, _ = _eval(u_try)
            if float(np.linalg.norm(res_try)) < cur_norm or float(np.linalg.norm(res_try)) / scale < rtol:
                u = u_try
                break
            alpha *= 0.5
        else:
            u = u + alpha * du
        res, K_T, pair_stats = _eval(u)
        if float(np.linalg.norm(res)) / scale < rtol:
            converged = True
            break

    # Recover full displacements and total nodal reactions R = K_full @ u_full
    u_full_list = [_unpack_full(i, u) for i in range(n_bodies)]
    reactions = [b["K_full"] @ uf for b, uf in zip(bodies, u_full_list)]
    return {
        "converged": converged,
        "iters": it + 1,
        "rel_res": float(np.linalg.norm(res)) / scale,
        "u_full": u_full_list,
        "reactions": reactions,
        "pair_stats": pair_stats,
    }


def _sum_rf(body: dict, R_full: np.ndarray, mask_key: str = "presc") -> np.ndarray:
    """Sum ``(RF1, RF2, RF3)`` over ``zero`` or ``presc`` DOFs of a body."""
    m = body[mask_key]
    R_m = R_full[m]
    comp = np.flatnonzero(m) % 3
    return np.array(
        [float(R_m[comp == 0].sum()), float(R_m[comp == 1].sum()), float(R_m[comp == 2].sum())],
        dtype=np.float64,
    )


def test_e05_reduced_3part_contact_load_transfer() -> None:
    """Reduced 3-part contact precursor: hinge (+Z pull) -> pin -> link (pinned), mu=0.0 & 0.3."""
    t0 = time.perf_counter()
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]

    E = 200000.0
    nu = 0.29
    u_pull_z = 0.2  # mm (+Z prescribed on hinge, same direction as doc (0, 0, 5) mm)
    eps_n = 5e5

    # Symmetric parabolic crown z(y) = 0.02 * y^2 so +Z pull generates both normal & tangential slip
    crown = lambda y: 0.02 * (y ** 2)
    box_hinge = _curved_crown_box(
        (-6.0, -5.0, -10.0), (6.0, 5.0, -4.0), (2, 2, 1),
        lambda _y: -10.0,
        lambda y: -4.0 + crown(y),
    )
    box_pin = _curved_crown_box(
        (-4.0, -4.0, -4.0), (4.0, 4.0, 0.0), (2, 2, 1),
        lambda y: -4.0 + crown(y),
        lambda y: 0.0 + crown(y),
    )
    box_link = _curved_crown_box(
        (-6.0, -5.0, 0.0), (6.0, 5.0, 6.0), (2, 2, 1),
        lambda y: 0.0 + crown(y),
        lambda _y: 6.0,
    )

    # BCs: hinge pulled at z=-10; link pinned at z=6; pin has NO Z-BC (only tiny X/Y rigid-mode reg)
    b_hinge = _assemble_kratos_body(
        KM, SM, box_hinge, E, nu,
        presc_fn=lambda x: (0, 1, 2) if abs(x[2] + 10.0) < 1e-10 else (),
    )
    b_pin = _assemble_kratos_body(KM, SM, box_pin, E, nu)
    # Add tiny X/Y-only diagonal stabilization (1e-2 N/mm) to pin so X/Y/Rz rigid modes are
    # non-singular at g_n=0 while Z remains 100% unconstrained (exact global Z force balance).
    import scipy.sparse as sp

    xy_diag = np.zeros(b_pin["n_dofs"], dtype=np.float64)
    xy_diag[0::3] = 1e-2
    xy_diag[1::3] = 1e-2
    K_xy_reg = sp.diags(xy_diag, format="csr")
    b_pin["K_full"] = (b_pin["K_full"] + K_xy_reg).tocsr()
    b_pin["Kf"] = (b_pin["Kf"] + K_xy_reg[b_pin["free"], :][:, b_pin["free"]]).tocsr()

    b_link = _assemble_kratos_body(
        KM, SM, box_link, E, nu,
        zero_fn=lambda x: (0, 1, 2) if abs(x[2] - 6.0) < 1e-10 else (),
    )

    u_fix_hinge = np.zeros(b_hinge["n_dofs"], dtype=np.float64)
    pm_h = np.flatnonzero(b_hinge["presc"])
    u_fix_hinge[pm_h[pm_h % 3 == 2]] = u_pull_z
    u_fix_pin = np.zeros(b_pin["n_dofs"], dtype=np.float64)
    u_fix_link = np.zeros(b_link["n_dofs"], dtype=np.float64)

    # 3 documented contact pairs: pin_to_hinge, pin_to_link, hinge_to_link
    # Map parametric Gauss points to physical curved crown coordinates X_s = B_s @ cp
    def _curved_face_gauss(box: BoxSpline, axis: int, side: int, n_gp: int, bounds):
        _, w_par, B_s = sample_box_face_gauss(box, axis=axis, side=side, n_gp=n_gp, bounds=bounds)
        X_phys = B_s @ box.control_points()
        jac = np.sqrt(1.0 + (0.04 * X_phys[:, 1]) ** 2)
        return X_phys, w_par * jac, B_s

    X_ph, w_ph, B_ph = _curved_face_gauss(box_hinge, axis=2, side=1, n_gp=3, bounds=((-3.5, 3.5), (-3.5, 3.5)))
    X_pl, w_pl, B_pl = _curved_face_gauss(box_link, axis=2, side=0, n_gp=3, bounds=((-3.5, 3.5), (-3.5, 3.5)))
    X_hl, w_hl, B_hl = _curved_face_gauss(box_link, axis=2, side=0, n_gp=2, bounds=((4.2, 5.8), (-4.0, 4.0)))

    pairs = [
        {
            "name": "pin_to_hinge_contact",
            "m_idx": 1,
            "s_idx": 0,
            "master_face": DeformableBoxFaceMaster(box_pin, axis=2, side=0),
            "X_s": X_ph,
            "w_q": w_ph,
            "B_s": B_ph,
        },
        {
            "name": "pin_to_link_contact",
            "m_idx": 1,
            "s_idx": 2,
            "master_face": DeformableBoxFaceMaster(box_pin, axis=2, side=1),
            "X_s": X_pl,
            "w_q": w_pl,
            "B_s": B_pl,
        },
        {
            "name": "hinge_to_link_contact",
            "m_idx": 0,
            "s_idx": 2,
            "master_face": DeformableBoxFaceMaster(box_hinge, axis=2, side=1),
            "X_s": X_hl,
            "w_q": w_hl,
            "B_s": B_hl,
        },
    ]

    rf3_by_mu = {}
    print("\n[E05 Reduced 3-Part Contact Precursor]")
    for mu in (0.0, 0.3):
        res = _solve_3body_ssn_contact(
            [b_hinge, b_pin, b_link],
            [u_fix_hinge, u_fix_pin, u_fix_link],
            pairs,
            eps_n=eps_n,
            mu=mu,
            eps_t=0.2 * eps_n,
        )
        assert res["converged"], f"mu={mu} did not converge: rel_res={res['rel_res']}"
        rf_hinge = _sum_rf(b_hinge, res["reactions"][0], "presc")
        rf_link = _sum_rf(b_link, res["reactions"][2], "zero")
        bal_z = abs(rf_hinge[2] + rf_link[2]) / max(abs(rf_link[2]), 1.0)
        rf3_by_mu[mu] = rf_link[2]

        ft_ph = float(np.linalg.norm(res["pair_stats"]["pin_to_hinge_contact"]["F_t"]))
        ft_pl = float(np.linalg.norm(res["pair_stats"]["pin_to_link_contact"]["F_t"]))
        print(
            f"  mu={mu:.1f}: iters={res['iters']:2d}, rel_res={res['rel_res']:.2e}, "
            f"hinge_RF3={rf_hinge[2]/1e3:+8.3f} kN, link_RF3={rf_link[2]/1e3:+8.3f} kN, "
            f"Z_bal={bal_z:.2e}, |Ft_ph|={ft_ph:.1f} N, |Ft_pl|={ft_pl:.1f} N"
        )

        # Link RF3 opposes +Z pull on hinge
        assert rf_hinge[2] > 0.0
        assert rf_link[2] < 0.0
        assert bal_z < 1e-7
        assert res["pair_stats"]["pin_to_hinge_contact"]["n_active"] > 0
        assert res["pair_stats"]["pin_to_link_contact"]["n_active"] > 0

    # Friction increases load transfer resistance compared to frictionless
    assert abs(rf3_by_mu[0.3]) > abs(rf3_by_mu[0.0])
    elapsed = time.perf_counter() - t0
    print(f"  elapsed={elapsed:.2f}s")
    assert elapsed < 30.0
