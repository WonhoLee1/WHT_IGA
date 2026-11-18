# -*- coding: utf-8 -*-
"""Unit tests for wht_iga.dynamics.explicit pure-python helpers (Kratos-free)."""

import numpy as np

from wht_iga.dynamics import explicit as EX
from wht_iga.reference.box import BoxSpline


def test_hrz_lump_diagonal_rescales_per_component():
    M = np.diag([1.0, 2.0, 3.0, 2.0, 1.0, 1.0])
    out = EX.hrz_lump_diagonal(M, total_mass=12.0)
    assert out.shape == (6,)
    for d in range(3):
        assert abs(out[d::3].sum() - 12.0) < 1e-12
    assert np.all(out > 0.0)


def test_build_x_extension_partition_of_unity():
    for degree in (1, 2, 3):
        box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(4.0, 1.0, 1.0),
                        n_cells=(4, 1, 1), degree=degree)
        info = EX.build_x_extension(box)
        T = info["T"]
        assert T.shape[0] == box.n_dofs
        assert T.shape[1] == (box.n_ctrl - box.n_cp[1] * box.n_cp[2]) * 3
        assert np.allclose(T.sum(axis=1), 1.0, atol=1e-12)
        assert len(info["slave_cp"]) == box.n_cp[1] * box.n_cp[2]


def test_floor_lumped_mass():
    m = EX.floor_lumped_mass(np.array([0.0, 0.05, 2.0]), 0.1)
    assert list(m) == [0.1, 0.1, 2.0]


def test_paired_penalty_forces():
    f = EX.paired_penalty_forces(np.array([0.1, 0.0, -0.02]), 1000.0)
    assert list(f) == [0.0, 0.0, 20.0]


def test_coulomb_stick_and_slip():
    import numpy as np
    f, ref, slipped = EX.coulomb_friction_update(np.array([0.01, 0.0]), 1000.0, 0.3, 100.0)
    assert slipped is False
    assert list(f) == [10.0, 0.0]
    assert list(ref) == [0.0, 0.0]
    f2, ref2, slipped2 = EX.coulomb_friction_update(np.array([0.1, 0.0]), 1000.0, 0.3, 100.0)
    assert slipped2 is True
    assert abs(float(np.linalg.norm(f2)) - 30.0) < 1e-12
    assert abs(float(np.linalg.norm(np.array([0.1, 0.0]) - ref2)) - 0.03) < 1e-12
    f3, _, slipped3 = EX.coulomb_friction_update(np.array([0.1, 0.0]), 1000.0, 0.3, 100.0, ref2)
    assert slipped3 is False
    assert abs(float(np.linalg.norm(f3)) - 30.0) < 1e-12
