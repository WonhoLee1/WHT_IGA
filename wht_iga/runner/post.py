"""Backend-independent post-processing of a control-point solution.

Probe values are evaluated from the control-point displacements with the exact
B-spline basis (value and first derivatives), so both backends report probes
through the same code path:

* ``U``: displacement ``Σ N_I u_I``;
* ``E``: small strain, Voigt ``[xx, yy, zz, xy, yz, xz]`` with engineering shear;
* ``S``: ``D(E, nu) · E`` (isotropic linear elastic, same Voigt order as Kratos
  ``CAUCHY_STRESS_VECTOR``);
* ``MISES``: von Mises of ``S``.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from wht_iga.grid.bspline import ders_basis_funs, find_span
from wht_iga.spec.study import PROBE_VARIABLE_CANONICAL


@njit(cache=True)
def _grad_box_field_njit(kx, ky, kz, p, nu, nv, nw, coeffs, points):
    """Gradient ``d u_c / d x_d`` (M, 3, 3) of a B-spline vector field (I = i + nu*(j + nv*k))."""
    m = points.shape[0]
    out = np.zeros((m, 3, 3))
    for n in range(m):
        px = points[n, 0]
        py = points[n, 1]
        pz = points[n, 2]
        sx = find_span(nu, p, px, kx)
        sy = find_span(nv, p, py, ky)
        sz = find_span(nw, p, pz, kz)
        Dx = ders_basis_funs(sx, px, p, 1, kx)
        Dy = ders_basis_funs(sy, py, p, 1, ky)
        Dz = ders_basis_funs(sz, pz, p, 1, kz)
        for c in range(p + 1):
            for b in range(p + 1):
                for a in range(p + 1):
                    I = (sx - p + a) + nu * ((sy - p + b) + nv * (sz - p + c))
                    gx = Dx[1, a] * Dy[0, b] * Dz[0, c]
                    gy = Dx[0, a] * Dy[1, b] * Dz[0, c]
                    gz = Dx[0, a] * Dy[0, b] * Dz[1, c]
                    for comp in range(3):
                        u = coeffs[I, comp]
                        out[n, comp, 0] += gx * u
                        out[n, comp, 1] += gy * u
                        out[n, comp, 2] += gz * u
    return out


def displacement_gradient(box, U: np.ndarray, points: np.ndarray) -> np.ndarray:
    """``(n, 3, 3)`` displacement gradient at physical points."""
    pts = np.ascontiguousarray(np.atleast_2d(points), dtype=np.float64)
    cf = np.ascontiguousarray(np.asarray(U, float).reshape(-1, 3))
    nu, nv, nw = (int(v) for v in box.n_cp)
    return _grad_box_field_njit(np.ascontiguousarray(box.knots[0], dtype=np.float64),
                                np.ascontiguousarray(box.knots[1], dtype=np.float64),
                                np.ascontiguousarray(box.knots[2], dtype=np.float64),
                                int(box.degree), nu, nv, nw, cf, pts)


def strain_voigt(grad: np.ndarray) -> np.ndarray:
    """Small strain ``[xx, yy, zz, xy, yz, xz]`` (engineering shear) from gradients."""
    g = np.asarray(grad, float)
    return np.stack([g[:, 0, 0], g[:, 1, 1], g[:, 2, 2],
                     g[:, 0, 1] + g[:, 1, 0], g[:, 1, 2] + g[:, 2, 1], g[:, 0, 2] + g[:, 2, 0]], axis=1)


def von_mises_voigt(s: np.ndarray) -> np.ndarray:
    """Von Mises stress of Voigt stresses ``[xx, yy, zz, xy, yz, xz]``."""
    s = np.asarray(s, float)
    return np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2)
                   + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))


def evaluate_point_fields(box, U: np.ndarray, points: np.ndarray, E: float, nu: float,
                          variables) -> dict[str, np.ndarray]:
    """Requested variables (keys exactly as requested, e.g. ``"S"`` or ``"CAUCHY_STRESS_TENSOR"``)."""
    from wht_iga.reference.assemble_elastic import isotropic_D

    pts = np.ascontiguousarray(np.atleast_2d(points), dtype=np.float64)
    canon = {v: PROBE_VARIABLE_CANONICAL[v] for v in variables}
    out: dict[str, np.ndarray] = {}
    eps = sig = None
    if any(c in ("E", "S", "MISES") for c in canon.values()):
        eps = strain_voigt(displacement_gradient(box, U, pts))
        sig = eps @ isotropic_D(E, nu).T
    for name, c in canon.items():
        if c == "U":
            out[name] = box.evaluate(np.asarray(U, float).reshape(-1, 3), pts)
        elif c == "E":
            out[name] = eps
        elif c == "S":
            out[name] = sig
        elif c == "MISES":
            out[name] = von_mises_voigt(sig)
        else:
            raise ValueError(f"variable {name!r} ({c}) is not a point field of a static structural solve")
    return out


def write_step_vtkhdf(path, box, U: np.ndarray, out_cells, E: float, nu: float) -> str:
    """Hex proxy over all non-OUT cells with point fields ``U`` and ``S`` + control-point ``U``.

    CUT cells are drawn as full hexahedra (fields there are the B-spline
    extension); the exact cut geometry is not reproduced in this proxy.
    """
    from wht_iga.post.proxy_mesh import generate_hex_proxy_for_cells
    from wht_iga.post.vtk_writer import write_vtkhdf

    out = {tuple(int(v) for v in c) for c in out_cells}
    nx, ny, nz = box.n_cells
    br = [box.breakpoints(d) for d in range(3)]
    bounds = [[[br[0][i], br[1][j], br[2][k]], [br[0][i + 1], br[1][j + 1], br[2][k + 1]]]
              for k in range(nz) for j in range(ny) for i in range(nx) if (i, j, k) not in out]
    pm = generate_hex_proxy_for_cells(np.asarray(bounds, float), n_sub=1, block_name="active_cells")
    fields = evaluate_point_fields(box, U, pm.points, E, nu, ("U", "S"))
    pm.point_data["U"] = fields["U"]
    pm.point_data["S"] = fields["S"]
    return write_vtkhdf(path, pm, control_points={0: {"U": np.asarray(U, float).reshape(-1, 3)}},
                        times=[1.0], block_names=["active_cells"])
