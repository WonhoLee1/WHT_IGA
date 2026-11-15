"""P4-C: Kratos conduction matrix versus the independent python assembly (contract).

``LaplacianElement`` local matrices scattered by control-point id must equal
``assemble_conduction`` (a separate tensor-product Gauss assembly) to round-off, and
must expose the constant null space of pure conduction.  The temperature-dependent
material path of spike S20 is covered by ``test_s20_accessor_temp.py`` (not repeated here).
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.analyses.thermal import assemble_conduction
from wht_iga.kratos_backend.thermal import kratos_box_conductivity
from wht_iga.reference.box import BoxSpline

pytestmark = pytest.mark.kratos


@pytest.mark.parametrize("p", [1, 2, 3])
def test_kratos_conductivity_equals_python_assembly(p: int) -> None:
    """Full (BC-free) K_th from Kratos equals the python assembly entry by entry."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(0.6, 1.0, 0.05), n_cells=(3, 4, 2), degree=p)
    k = 52.0
    K_kr = kratos_box_conductivity(box, k)
    K_py = assemble_conduction(box, k)
    assert K_kr.shape == (box.n_ctrl, box.n_ctrl)
    assert abs(K_kr - K_py).max() < 1e-12 * abs(K_py).max()
    assert np.abs(K_kr @ np.ones(box.n_ctrl)).max() < 1e-11 * abs(K_kr).max()  # constants are in the null space


def test_kratos_conductivity_scales_linearly_with_k() -> None:
    """Doubling the conductivity doubles the matrix (the property reaches the element)."""
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(0.1, 0.01, 0.01), n_cells=(4, 1, 1), degree=2)
    a, b = kratos_box_conductivity(box, 10.0), kratos_box_conductivity(box, 20.0)
    assert abs(b - 2.0 * a).max() < 1e-13 * abs(b).max()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
