# -*- coding: utf-8 -*-
"""QP element-name guard: unverified (mixed / B-bar / Q1P0) elements are rejected.

S13: Mixed*VolumetricStrain crashes on QP geometry; S39: TotalLagrangianQ1P0Mixed*
gives NaN / detF < 0 on the Cook model. The guard turns these silent or fatal
failures into an immediate ValueError before any element is created.
"""

from __future__ import annotations

import pytest

from wht_iga.kratos_backend.elements import check_qp_element_name


@pytest.mark.parametrize("name", [
    "SmallDisplacementElement3D4N",
    "SmallDisplacementElement3D8N",
    "TotalLagrangianElement3D4N",
    "UpdatedLagrangianElement3D4N",
])
def test_verified_families_pass(name):
    assert check_qp_element_name(name) == name


@pytest.mark.parametrize("name", [
    "TotalLagrangianQ1P0MixedElement3D8N",
    "SmallDisplacementMixedVolumetricStrainElement3D4N",
    "SmallDisplacementMixedStrainElement3D4N",
    "SmallDisplacementBbarElement3D8N",
    "TotalLagrangianMixedVolumetricStrainElement3D4N",
])
def test_unverified_elements_rejected(name):
    with pytest.raises(ValueError, match="not verified on QuadraturePointGeometry"):
        check_qp_element_name(name)


def test_attach_cell_elements_rejects_before_creating_anything():
    """No Kratos objects are touched: the guard runs first."""
    from wht_iga.kratos_backend.elements import attach_cell_elements

    class _Boom:
        def __getattr__(self, item):  # any Kratos call would hit this
            raise AssertionError(f"Kratos touched: {item}")

    with pytest.raises(ValueError):
        attach_cell_elements(_Boom(), _Boom(), _Boom(), _Boom(), _Boom(),
                             element_name="TotalLagrangianQ1P0MixedElement3D8N")
