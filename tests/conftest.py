"""Shared fixtures.

Tests marked ``kratos`` need the Kratos source build. If Kratos cannot be loaded
they FAIL (no skip): deselect them explicitly with ``-m "not kratos"``.
"""

from __future__ import annotations

import pytest


@pytest.fixture(scope="session")
def kratos():
    from wht_iga.kratos_backend import load_kratos

    return load_kratos()
