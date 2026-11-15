"""inactive_ctrl_mask vs an independent oracle (basis evaluated at cell centres).

A control point is inactive iff its basis function vanishes on every non-OUT
cell. The oracle evaluates the unit-coefficient field with BoxSpline.evaluate
at each active cell's centre (B-splines are strictly positive inside their
support), so it does not share the index arithmetic under test.
"""

import numpy as np
import pytest

from wht_iga.solver.active import inactive_ctrl_mask
from wht_iga.reference.box import BoxSpline


@pytest.mark.parametrize("p", [1, 2, 3])
def test_inactive_mask_matches_basis_support(p):
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 0.75, 0.5), n_cells=(4, 3, 2), degree=p)
    rng = np.random.default_rng(p)
    cells = [(i, j, k) for k in range(2) for j in range(3) for i in range(4)]
    out_set = {c for c in cells if rng.random() < 0.6}
    active = [c for c in cells if c not in out_set]
    h = (np.asarray(box.hi) - np.asarray(box.lo)) / np.asarray(box.n_cells)
    centres = np.array([np.asarray(box.lo) + (np.asarray(c) + 0.5) * h for c in active])

    mask = inactive_ctrl_mask(box, out_set)
    assert mask.shape == (box.n_ctrl,)
    for I in range(box.n_ctrl):
        coeffs = np.zeros((box.n_ctrl, 3))
        coeffs[I, 0] = 1.0
        vals = box.evaluate(coeffs, centres)[:, 0] if len(centres) else np.zeros(0)
        oracle_inactive = not np.any(np.abs(vals) > 1e-14)
        assert mask[I] == oracle_inactive, (p, I)
    assert 0 < mask.sum() < box.n_ctrl
