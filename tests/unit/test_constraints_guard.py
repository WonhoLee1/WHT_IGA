"""c = 0 rule guard on ``add_linear_tie`` (REVIEW_P0_P1 §1 #1, v2.1 §5.2.3).

The guard runs before any Kratos object is touched, so no Kratos import is
needed: a nonzero or mis-shaped constant must raise ``ValueError``.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend.constraints import add_linear_tie


@pytest.mark.parametrize("constant", [np.array([1e-6]), np.array([-2.0]), [0.0, 0.0]])
def test_add_linear_tie_rejects_nonzero_or_misshaped_constant(constant) -> None:
    with pytest.raises(ValueError):
        add_linear_tie(None, None, 1, ["s"], ["m"], np.array([[1.0]]), constant)
