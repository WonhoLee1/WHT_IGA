"""P4-A: consistent mass, Kratos stiffness bridge, eigenmodes, modal stress recovery.

Oracles
-------
* Total mass ``rho * V`` and the linear-field energy ``c^T M c = rho int x^2 dV``:
  B-splines with Greville control points reproduce linear fields exactly, so both
  have closed forms on any axis-aligned sub-box (no use of the code under test).
* ``assemble_stiffness`` of the reference stack for the Kratos stiffness (independent
  assembler, same discretisation).
* Euler-Bernoulli cantilever frequency and the uniform-strain stress field.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from wht_iga.dynamics import assemble_consistent_mass_from_bspline
from wht_iga.kratos_backend.mass import (
    assemble_box_mass,
    face_dofs,
    kratos_box_stiffness,
    kratos_recover_modal_stresses,
    kratos_total_mass,
    recover_modal_stresses,
    reduce_by_operator,
    solve_kratos_eigenmodes,
)
from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline

RHO = 7.5


def _box(p: int = 2, n=(3, 2, 2), hi=(3.0, 1.0, 2.0)) -> BoxSpline:
    return BoxSpline(lo=(0.0, 0.0, 0.0), hi=hi, n_cells=n, degree=p)


def _linear_x_coeffs(box: BoxSpline, comp: int) -> np.ndarray:
    """DOF vector of the field u_comp = x (Greville control points reproduce it exactly)."""
    c = np.zeros(box.n_dofs)
    c[comp::3] = box.control_points()[:, 0]
    return c


# ---------------------------------------------------------------------------
# Kratos-free mass math
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [1, 2, 3])
def test_general_path_equals_fast_path_when_ghost_alpha_is_one(p: int) -> None:
    """A ghost cell with alpha = 1 is a full cell: the per-cell path must equal the separable path."""
    box = _box(p, n=(4, 2, 2))
    fast = assemble_consistent_mass_from_bspline(box.knots, degrees=p, rho=RHO)
    general = assemble_box_mass(box, RHO, out_cells=[(1, 1, 0)], ghost_alpha=1.0)
    assert abs(fast - general).max() < 1e-12 * abs(fast).max()


def test_drop_and_ghost_remove_or_scale_out_cell_mass() -> None:
    """Total mass: full = rho V; drop removes the OUT cell volume; ghost keeps alpha of it."""
    box = _box(2, n=(4, 2, 2), hi=(4.0, 2.0, 2.0))  # cells 1 x 1 x 1
    out = [(0, 0, 0), (3, 1, 1)]
    r = np.zeros(box.n_dofs)
    r[0::3] = 1.0
    total = lambda M: float(r @ (M @ r))  # noqa: E731 - one-liner oracle helper
    assert total(assemble_box_mass(box, RHO)) == pytest.approx(RHO * 16.0, rel=1e-13)
    assert total(assemble_box_mass(box, RHO, out_cells=out, ghost_alpha=None)) == pytest.approx(RHO * 14.0, rel=1e-13)
    assert total(assemble_box_mass(box, RHO, out_cells=out, ghost_alpha=1e-3)) == pytest.approx(RHO * (14.0 + 2e-3), rel=1e-13)


@pytest.mark.parametrize("p", [1, 2, 3])
def test_cut_cell_table_mass_and_linear_energy_closed_forms(p: int) -> None:
    """A cut cell given as its sub-box Gauss rule [0.4, c]: mass and c^T M c match analytic integrals."""
    box = _box(p, n=(3, 2, 2), hi=(3.0, 2.0, 2.0))  # unit cells
    cut = (1, 0, 1)
    c_cut = 1.35  # cut plane inside cell x in [1, 2]
    xg, wg = np.polynomial.legendre.leggauss(p + 2)
    px = 1.0 + 0.5 * (c_cut - 1.0) * (xg + 1.0)
    wx = 0.5 * (c_cut - 1.0) * wg
    py, wy = 0.5 * (xg + 1.0), 0.5 * wg  # y in [0, 1]
    pz, wz = 1.0 + 0.5 * (xg + 1.0), 0.5 * wg  # z in [1, 2]
    X, Y, Z = np.meshgrid(px, py, pz, indexing="ij")
    W = np.einsum("i,j,k->ijk", wx, wy, wz)
    table = np.column_stack([X.ravel(), Y.ravel(), Z.ravel(), W.ravel()])
    M = assemble_box_mass(box, RHO, qp_table={cut: table})

    # exact: all cells except `cut` are full unit cells, the cut cell keeps x in [1, c_cut]
    vol = 12.0 - (2.0 - c_cut)
    r = np.zeros(box.n_dofs)
    r[2::3] = 1.0
    assert float(r @ (M @ r)) == pytest.approx(RHO * vol, rel=1e-12)
    c = _linear_x_coeffs(box, 0)
    int_x2 = 0.0
    for k in range(2):
        for j in range(2):
            for i in range(3):
                if (i, j, k) == cut:
                    int_x2 += (c_cut**3 - 1.0) / 3.0 * 1.0 * 1.0
                else:
                    int_x2 += ((i + 1) ** 3 - i**3) / 3.0
    assert float(c @ (M @ c)) == pytest.approx(RHO * int_x2, rel=1e-12)
    assert abs(M - M.T).max() < 1e-13 * abs(M).max()


def test_mass_is_positive_definite_on_active_support_only() -> None:
    """With dropped OUT cells, DOFs whose support lies entirely in OUT cells carry zero mass."""
    box = _box(1, n=(3, 1, 1), hi=(3.0, 1.0, 1.0))
    M = assemble_box_mass(box, RHO, out_cells=[(2, 0, 0)], ghost_alpha=None)
    diag = M.diagonal().reshape(-1, 3)[:, 0]
    nu, nv, nw = box.n_cp
    last_layer = [a + nu * (b + nv * c) for c in range(nw) for b in range(nv) for a in (nu - 1,)]
    assert np.all(diag[last_layer] == 0.0)  # x-layer 3 touches only the dropped cell 2
    assert np.all(np.delete(diag, last_layer) > 0.0)


def test_reduce_by_operator_matches_dense_triple_product() -> None:
    """T^T A T equals the dense product and stays symmetric."""
    rng = np.random.default_rng(3)
    A = rng.normal(size=(7, 7))
    A = A + A.T
    T = sp.csr_matrix(rng.normal(size=(7, 4)))
    got = reduce_by_operator(sp.csr_matrix(A), T).toarray()
    np.testing.assert_allclose(got, T.toarray().T @ A @ T.toarray(), rtol=1e-13, atol=1e-13)
    assert np.allclose(got, got.T, atol=0.0)


def test_modal_stress_recovery_uniform_and_linear_strain_fields() -> None:
    """Field u = (a x, 0, 0): sigma_xx = (lam + 2 mu) a, sigma_yy = sigma_zz = lam a, shear 0 (Hooke)."""
    E, nu = 3.0, 0.25
    lam, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    for p in (1, 2, 3):
        box = _box(p, n=(3, 2, 2), hi=(3.0, 1.0, 2.0))
        a = 0.7
        modes = (a * _linear_x_coeffs(box, 0))[:, None]
        pts = np.array([[0.2, 0.3, 0.4], [1.0, 0.5, 1.0], [2.9, 0.99, 1.7], [3.0, 1.0, 2.0], [0.0, 0.0, 0.0]])
        psi = recover_modal_stresses(box, E, nu, modes, pts)
        want = np.array([(lam + 2 * mu) * a, lam * a, lam * a, 0.0, 0.0, 0.0])
        np.testing.assert_allclose(psi[:, 0, :], np.tile(want, (len(pts), 1)), rtol=0, atol=1e-12)
    # pure shear u_x = g y: gamma_xy = g -> sigma_xy = mu g in Kratos Voigt slot 3
    box = _box(2)
    c = np.zeros(box.n_dofs)
    c[0::3] = 0.4 * box.control_points()[:, 1]
    psi = recover_modal_stresses(box, E, nu, c[:, None], np.array([[1.3, 0.5, 0.9]]))
    np.testing.assert_allclose(psi[0, 0], [0.0, 0.0, 0.0, mu * 0.4, 0.0, 0.0], atol=1e-12)
    # rigid rotation u = theta (-y, x, 0) is strain-free: needs the symmetric part du_x/dy + du_y/dx
    rot = np.zeros(box.n_dofs)
    cp = box.control_points()
    rot[0::3], rot[1::3] = -0.3 * cp[:, 1], 0.3 * cp[:, 0]
    psi_rot = recover_modal_stresses(box, E, nu, rot[:, None], np.array([[1.3, 0.5, 0.9], [0.1, 0.9, 1.8]]))
    np.testing.assert_allclose(psi_rot, 0.0, atol=1e-12)
    # mixed shear u = (g y, h x, 0): gamma_xy = g + h
    mix = np.zeros(box.n_dofs)
    mix[0::3], mix[1::3] = 0.4 * cp[:, 1], 0.25 * cp[:, 0]
    psi_mix = recover_modal_stresses(box, E, nu, mix[:, None], np.array([[1.3, 0.5, 0.9]]))
    np.testing.assert_allclose(psi_mix[0, 0], [0.0, 0.0, 0.0, mu * 0.65, 0.0, 0.0], atol=1e-12)


def test_face_dofs_components_and_sides() -> None:
    """face_dofs returns 3 I + d for the control points of the requested face only."""
    box = _box(2, n=(3, 2, 2))
    nu, nv, nw = box.n_cp
    d = face_dofs(box, 0, 0)
    assert d.size == 3 * nv * nw
    assert set(d % 3) == {0, 1, 2}
    assert np.all((d // 3) % nu == 0)
    dy = face_dofs(box, 0, 1, components=(1,))
    assert dy.size == nv * nw and np.all(dy % 3 == 1) and np.all((dy // 3) % nu == nu - 1)


# ---------------------------------------------------------------------------
# Kratos contract
# ---------------------------------------------------------------------------
@pytest.mark.kratos
@pytest.mark.parametrize("side", [0, 1])
@pytest.mark.parametrize("p", [1, 2])
def test_kratos_stiffness_equals_reference_assembler(p: int, side: int) -> None:
    """Kratos route-R1 K (permuted to reference order) equals the independent assembler's K.

    The clamped face is both the low (first control points) and the high (last control
    points) x face, so a wrong equation-id ordering assumption would show up.
    """
    box = _box(p, n=(4, 2, 2), hi=(2.0, 0.5, 0.5))
    E, nu = 2.1, 0.3
    fixed = face_dofs(box, 0, side)
    K = kratos_box_stiffness(box, E, nu, fixed)
    Kref = assemble_stiffness(box, E, nu)
    free = np.setdiff1d(np.arange(box.n_dofs), fixed)
    a, b = K[free][:, free], Kref[free][:, free]
    assert abs(a - b).max() < 1e-12 * abs(b).max()
    assert a.nnz == b.nnz


@pytest.mark.kratos
def test_kratos_lumped_total_mass_equals_rho_volume() -> None:
    """Kratos' (valid) lumped mass sums to rho V; the python consistent mass agrees."""
    box = _box(2, n=(3, 2, 2), hi=(3.0, 1.0, 2.0))
    assert kratos_total_mass(box, RHO) == pytest.approx(RHO * 6.0, rel=1e-13)
    r = np.zeros(box.n_dofs)
    r[1::3] = 1.0
    M = assemble_box_mass(box, RHO)
    assert float(r @ (M @ r)) == pytest.approx(kratos_total_mass(box, RHO), rel=1e-13)


@pytest.mark.kratos
def test_kratos_cantilever_modes_euler_bernoulli_and_orthonormality() -> None:
    """Square cantilever: degenerate bending pair near the Euler-Bernoulli f1; M-orthonormal; mass sums."""
    L, b, E, rho = 1.0, 0.1, 2.0e11, 8000.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(b, b, L), n_cells=(2, 2, 16), degree=2)
    model = solve_kratos_eigenmodes(box, E, 0.0, rho, face_dofs(box, 2, 0), n_modes=6)
    eig = model.eig
    f_eb = (1.8751040687119611**2 / (2.0 * np.pi * L**2)) * np.sqrt(E * (b**4 / 12.0) / (rho * b * b))
    assert eig.frequencies_hz[0] == pytest.approx(f_eb, rel=0.01)  # shear flexibility keeps it within ~1 %
    assert eig.frequencies_hz[1] == pytest.approx(eig.frequencies_hz[0], rel=1e-6)  # square section
    gram = eig.modes.T @ (model.M @ eig.modes)
    assert np.abs(gram - np.eye(6)).max() < 1e-10
    assert eig.total_mass[0] == pytest.approx(rho * b * b * L, rel=1e-12)
    assert np.all(eig.effective_modal_masses.sum(axis=0) <= eig.total_mass * (1.0 + 1e-12))


def _independent_DB_modal_stress(
    box: BoxSpline, E: float, nu: float, modes: np.ndarray, points: np.ndarray
) -> np.ndarray:
    """Independent sigma_n(x) = D B(x) phi_n evaluator using basis_1d + Lame Hooke's law."""
    from wht_iga.grid.bspline import basis_1d

    p = box.degree
    nu_cp, nv_cp, _nw_cp = box.n_cp
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    D = np.array(
        [
            [lam + 2 * mu, lam, lam, 0.0, 0.0, 0.0],
            [lam, lam + 2 * mu, lam, 0.0, 0.0, 0.0],
            [lam, lam, lam + 2 * mu, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, mu, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, mu, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, mu],
        ],
        dtype=float,
    )
    n_modes = modes.shape[1]
    out = np.zeros((points.shape[0], n_modes, 6), dtype=float)
    for q_idx, pt in enumerate(points):
        s_x, dx = basis_1d(box.knots[0], p, float(pt[0]), 1)
        s_y, dy = basis_1d(box.knots[1], p, float(pt[1]), 1)
        s_z, dz = basis_1d(box.knots[2], p, float(pt[2]), 1)
        n_loc = (p + 1) ** 3
        B = np.zeros((6, 3 * n_loc), dtype=float)
        dofs = np.empty(3 * n_loc, dtype=np.int64)
        loc = 0
        for c in range(p + 1):
            for b in range(p + 1):
                for a in range(p + 1):
                    cp = (int(s_x) - p + a) + nu_cp * ((int(s_y) - p + b) + nv_cp * (int(s_z) - p + c))
                    gx = dx[1, a] * dy[0, b] * dz[0, c]
                    gy = dx[0, a] * dy[1, b] * dz[0, c]
                    gz = dx[0, a] * dy[0, b] * dz[1, c]
                    c0 = 3 * loc
                    B[0, c0] = gx
                    B[1, c0 + 1] = gy
                    B[2, c0 + 2] = gz
                    B[3, c0] = gy
                    B[3, c0 + 1] = gx
                    B[4, c0 + 1] = gz
                    B[4, c0 + 2] = gy
                    B[5, c0] = gz
                    B[5, c0 + 2] = gx
                    dofs[c0 : c0 + 3] = (3 * cp, 3 * cp + 1, 3 * cp + 2)
                    loc += 1
        out[q_idx] = (D @ (B @ modes[dofs, :])).T
    return out


@pytest.mark.kratos
def test_kratos_zero_weight_modal_stress_parity_on_beam_modes() -> None:
    """Zero-weight Kratos probe modal stresses match independent D*B*phi_n on bending and axial modes."""
    from wht_iga.kratos_backend.elements import solve_linear, system_matrix_csr
    from wht_iga.kratos_backend.mass import _build_kratos_box

    L, b, E, nu, rho = 1.0, 0.1, 2.0e11, 0.3, 8000.0
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(b, b, L), n_cells=(2, 2, 8), degree=2)
    fixed = face_dofs(box, 2, 0)
    modal = solve_kratos_eigenmodes(box, E, nu, rho, fixed, n_modes=6)
    modes = modal.eig.modes
    pf = modal.eig.participation_factors
    # Verify both bending modes (x/y participation) and axial mode (z participation) are tested
    assert np.max(np.abs(pf[:, 0])) > 1.0 and np.max(np.abs(pf[:, 1])) > 1.0
    assert np.max(np.abs(pf[:, 2])) > 1.0

    pts = np.array(
        [
            [0.02, 0.03, 0.15],
            [0.08, 0.07, 0.55],
            [0.03, 0.08, 0.80],
            [0.05, 0.05, 0.50],  # knot-line point
            [0.00, 0.04, 0.25],  # boundary face point
        ],
        dtype=float,
    )

    psi_kr = kratos_recover_modal_stresses(box, E, nu, modes, pts, verify_zero_contribution=True)
    psi_kw = recover_modal_stresses(box, E, nu, modes, pts, backend="kratos")
    psi_py = recover_modal_stresses(box, E, nu, modes, pts)
    psi_ref = _independent_DB_modal_stress(box, E, nu, modes, pts)

    np.testing.assert_allclose(psi_kr, psi_kw, rtol=0.0, atol=0.0)
    for k in range(modes.shape[1]):
        scale = float(np.max(np.abs(psi_ref[:, k, :])))
        assert scale > 0.0
        rel_err_ref = float(np.max(np.abs(psi_kr[:, k, :] - psi_ref[:, k, :])) / scale)
        rel_err_py = float(np.max(np.abs(psi_kr[:, k, :] - psi_py[:, k, :])) / scale)
        assert rel_err_ref < 1e-10
        assert rel_err_py < 1e-10

    # Verify adding zero-weight probe elements into a solved Kratos ModelPart leaves global K untouched
    KM, mp_base, _vol_b, _els_b = _build_kratos_box(box, E, nu, rho)
    KM2, mp_probe, vol_p, _els_p = _build_kratos_box(box, E, nu, rho)
    comp_vars = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)
    for dof in fixed:
        mp_base.GetNode(int(dof) // 3 + 1).Fix(comp_vars[int(dof) % 3])
        mp_probe.GetNode(int(dof) // 3 + 1).Fix(comp_vars[int(dof) % 3])
    props_p = mp_probe.GetProperties()[1]
    for idx, pt in enumerate(pts):
        qps = KM2.GeometriesVector()
        vol_p.CreateQuadraturePointGeometries(qps, 2, [[float(pt[0]), float(pt[1]), float(pt[2]), 0.0]])
        mp_probe.CreateNewElement("SmallDisplacementElement3D4N", 900001 + idx, qps[0], props_p)
    K_base = system_matrix_csr(solve_linear(KM, mp_base))
    K_probe = system_matrix_csr(solve_linear(KM2, mp_probe))
    assert abs(K_probe - K_base).max() == 0.0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
