"""Unit and Kratos-backed verification tests for bimaterial Nitsche + AgFEM (V11 / P4).

Covers:
1. Quarter-disc and planar cut-region volume conservation & interface arc quadrature accuracy.
2. Homogeneous and dissimilar-modulus free thermal expansion patch tests (stress == 0 to 1e-11).
3. Cut-fraction sliver stability (``eta in {1e-2, 1e-4, 1e-6}``): bounded condition number and
   eta-independent Nitsche coercivity with AgFEM extension + volume-fraction weighting.
4. Negative / mutation control: insufficient Nitsche penalty (``C_gamma = 2``) fails coercivity.
5. V11 quarter-disc coercivity sweep across ``theta`` and ``C_gamma``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from wht_iga.analyses.bimaterial import (
    interface_coercivity_diagnostic,
    isotropic_D,
    kratos_region_stiffness,
    planar_bimaterial_regions,
    quarter_disc_bimaterial_regions,
    solve_bimaterial_thermal_stress,
)
from wht_iga.analyses.thermal import eigenstrain_load, plane_stress_as_plane_strain
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.trimmed import _strain_at_points


def _symmetry_fixed_dofs(box: BoxSpline) -> list[int]:
    """Return symmetry DOFs (u_x(x=0)=0, u_y(y=0)=0, u_z=0 everywhere)."""
    fix = (
        [3 * i for i in box.boundary_control_points(0, 0)]
        + [3 * i + 1 for i in box.boundary_control_points(1, 0)]
        + [3 * i + 2 for i in range(box.n_ctrl)]
    )
    return sorted(set(fix))


def test_quarter_disc_and_planar_region_conservation() -> None:
    """Volume conservation vol_0 + vol_1 == cell_vol and exact quarter-arc geometry."""
    box = BoxSpline((0.0, 0.0, 0.0), (10.0, 10.0, 0.25), (16, 16, 1), degree=1)
    reg = quarter_disc_bimaterial_regions(box, radius=1.0, q=4)

    bx, by, bz = (box.breakpoints(d) for d in range(3))
    hx, hy, hz = np.diff(bx), np.diff(by), np.diff(bz)
    cvol = np.einsum("i,j,k->kji", hx, hy, hz).reshape(-1)
    np.testing.assert_allclose(reg.cell_volumes[0] + reg.cell_volumes[1], cvol, rtol=1e-13, atol=1e-14)
    assert reg.cell_volumes[0].sum() == pytest.approx(0.25 * math.pi * 1.0**2 * 0.25, rel=1e-12)

    arc_area = sum(float(pts_w[:, 3].sum()) for pts_w, _nrm in reg.interface.values())
    assert arc_area == pytest.approx(0.5 * math.pi * 1.0 * 0.25, rel=1e-10)
    for pts_w, nrms in reg.interface.values():
        r = np.hypot(pts_w[:, 0], pts_w[:, 1])
        np.testing.assert_allclose(r, 1.0, atol=1e-13)
        np.testing.assert_allclose(np.linalg.norm(nrms, axis=1), 1.0, atol=1e-13)

    preg = planar_bimaterial_regions(box, axis=0, interface_pos=2.1, q=2)
    np.testing.assert_allclose(preg.cell_volumes[0] + preg.cell_volumes[1], cvol, rtol=1e-13, atol=1e-14)
    assert preg.cell_volumes[0].sum() == pytest.approx(2.1 * 10.0 * 0.25, rel=1e-12)


@pytest.mark.parametrize("e2_factor, nu1, nu2", [(1.0, 0.3, 0.3), (2.5, 0.3, 0.22)])
def test_bimaterial_free_thermal_expansion_patch_test(
    e2_factor: float, nu1: float, nu2: float
) -> None:
    """Constant isotropic eigenstrain across a bimaterial interface yields zero stress."""
    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 1.0, 0.25), (4, 2, 1), degree=1)
    e1 = 1e5
    e2 = e2_factor * e1
    eps0 = np.array([1e-3, 1e-3, 1e-3, 0.0, 0.0, 0.0])

    reg = planar_bimaterial_regions(box, axis=0, interface_pos=1.05, q=2)
    fix = sorted(
        set(
            [3 * i for i in box.boundary_control_points(0, 0)]
            + [3 * i + 1 for i in box.boundary_control_points(1, 0)]
            + [3 * i + 2 for i in box.boundary_control_points(2, 0)]
        )
    )
    res = solve_bimaterial_thermal_stress(
        reg, e1, nu1, e2, nu2, fix, fix, eigenstrain1=eps0, eigenstrain2=eps0, theta=0.2, C_gamma=10.0
    )

    d1, d2 = isotropic_D(e1, nu1), isotropic_D(e2, nu2)
    pts1 = np.array([[0.25, 0.5, 0.125], [1.04, 0.5, 0.125]])
    pts2 = np.array([[1.06, 0.5, 0.125], [1.75, 0.5, 0.125]])
    sig1 = (_strain_at_points(box, res.u1, pts1) - eps0) @ d1.T
    sig2 = (_strain_at_points(box, res.u2, pts2) - eps0) @ d2.T
    assert float(np.max(np.abs(sig1))) < 1e-10
    assert float(np.max(np.abs(sig2))) < 1e-10

    ax, ay, az = box.greville_axes()
    zz, yy, xx = np.meshgrid(az, ay, ax, indexing="ij")
    u_lin = (1e-3 * np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])).reshape(-1)
    k1 = kratos_region_stiffness(box, e1, nu1, reg.cut_tables[0], reg.full_cells[0])
    f1 = eigenstrain_load(box, d1, eps0, reg.cut_tables[0], reg.full_cells[0])
    assert float(np.linalg.norm(k1 @ u_lin - f1)) < 1e-11


def test_sliver_eta_sweep_coercivity_and_mutation_control() -> None:
    """Sliver cut fractions eta down to 1e-6 keep uniform coercivity; low C_gamma breaks it."""
    box = BoxSpline((0.0, 0.0, 0.0), (2.0, 1.0, 0.25), (4, 2, 1), degree=1)
    e1, e2, nu = 1e5, 2e5, 0.3
    eps0 = np.array([1e-3, 1e-3, 1e-3, 0.0, 0.0, 0.0])
    fix = sorted(
        set(
            [3 * i for i in box.boundary_control_points(0, 0)]
            + [3 * i + 1 for i in box.boundary_control_points(1, 0)]
            + [3 * i + 2 for i in box.boundary_control_points(2, 0)]
        )
    )

    ratios = []
    for eta in (1e-2, 1e-4, 1e-6):
        reg = planar_bimaterial_regions(box, axis=0, interface_pos=0.5 + 0.5 * eta, q=2)

        res_ok = solve_bimaterial_thermal_stress(
            reg, e1, nu, e2, nu, fix, fix, eigenstrain1=eps0, eigenstrain2=eps0, theta=0.25, C_gamma=10.0
        )
        d_ok = interface_coercivity_diagnostic(res_ok)
        assert d_ok["coercive"] == 1.0
        assert d_ok["lambda_min_coercivity_ratio"] > 0.5
        assert d_ok["kappa_A"] < 1e3
        ratios.append(d_ok["lambda_min_coercivity_ratio"])

        res_bad = solve_bimaterial_thermal_stress(
            reg, e1, nu, e2, nu, fix, fix, eigenstrain1=eps0, eigenstrain2=eps0, theta=0.25, C_gamma=2.0
        )
        d_bad = interface_coercivity_diagnostic(res_bad)
        assert d_bad["coercive"] == 0.0
        assert d_bad["lambda_min_coercivity_ratio"] < 0.0

    # Nearly eta-independent coercivity across 4 orders of magnitude of sliver cut fraction.
    assert max(ratios) - min(ratios) < 0.05


def test_quarter_disc_coercivity_sweep() -> None:
    """Coercivity holds on the V11 quarter-disc geometry across (theta, C_gamma)."""
    e_p, nu_p, eps_p = plane_stress_as_plane_strain(1e5, 0.3, 1e-3)
    eps1 = np.array([eps_p, eps_p, eps_p, 0.0, 0.0, 0.0])
    box = BoxSpline((0.0, 0.0, 0.0), (10.0, 10.0, 0.25), (16, 16, 1), degree=1)
    reg = quarter_disc_bimaterial_regions(box, radius=1.0, q=4)
    fix = _symmetry_fixed_dofs(box)

    for theta in (0.1, 0.2, 0.35):
        for c_gam in (5.0, 10.0, 20.0):
            res = solve_bimaterial_thermal_stress(
                reg, e_p, nu_p, e_p, nu_p, fix, fix, eigenstrain1=eps1, eigenstrain2=None,
                theta=theta, C_gamma=c_gam,
            )
            diag = interface_coercivity_diagnostic(res)
            assert diag["coercive"] == 1.0
            assert diag["lambda_min_coercivity_ratio"] > 0.7
            assert diag["kappa_A"] < 5e3
