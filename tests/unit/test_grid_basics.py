"""Unit tests for grid basics (no Kratos). Implementation request T-P0.4."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from wht_iga.grid.bezier_extraction import bernstein, extraction_1d
from wht_iga.grid.bspline import basis_1d
from wht_iga.grid.greville import greville
from wht_iga.grid.knots import from_kratos, n_control_points, to_kratos, uniform_clamped

CASES = [(p, k, n) for p in (1, 2, 3) for k in range(p) for n in (1, 3, 5)]


def all_basis(knots, p, u, n_ders=0):
    n_cp = n_control_points(knots, p)
    span, ders = basis_1d(knots, p, u, n_ders)
    full = np.zeros((n_ders + 1, n_cp))
    full[:, span - p : span + 1] = ders
    return full


@pytest.mark.parametrize("p,k,n", CASES)
def test_partition_of_unity(p, k, n):
    knots = uniform_clamped(-1.0, 2.0, n, p, k)
    for u in np.linspace(-1.0, 2.0, 37):
        N = all_basis(knots, p, u, 1)
        np.testing.assert_allclose(N[0].sum(), 1.0, atol=1e-14)
        np.testing.assert_allclose(N[1].sum(), 0.0, atol=1e-12)


@pytest.mark.parametrize("p,k,n", CASES)
def test_derivatives_match_finite_difference(p, k, n):
    knots = uniform_clamped(0.0, 1.0, n, p, k)
    breaks = np.unique(knots)
    h = 1e-6
    rng = np.random.default_rng(0)
    for _ in range(20):
        u = rng.uniform(0.0, 1.0)
        if np.min(np.abs(breaks - u)) < 10 * h:
            continue  # derivative may be discontinuous at C0 knots
        N = all_basis(knots, p, u, 1)
        fd = (all_basis(knots, p, u + h)[0] - all_basis(knots, p, u - h)[0]) / (2 * h)
        np.testing.assert_allclose(N[1], fd, atol=1e-6)


@pytest.mark.parametrize("p,k,n", CASES)
def test_kratos_knot_roundtrip(p, k, n):
    knots = uniform_clamped(0.0, 4.0, n, p, k)
    kr = to_kratos(knots, p)
    assert len(kr) == n_control_points(knots, p) + p - 1
    np.testing.assert_array_equal(from_kratos(kr, p), knots)


def test_kratos_single_span_convention():
    # IgaApplication/tests/test_nurbs_volume_element.py: one span of degree p -> [0]*p + [1]*p
    for p in (1, 2, 3):
        np.testing.assert_array_equal(to_kratos(uniform_clamped(0.0, 1.0, 1, p), p), [0.0] * p + [1.0] * p)


@pytest.mark.parametrize("p,k,n", CASES)
def test_identity_map_with_greville_points(p, k, n):
    knots = uniform_clamped(-2.0, 3.0, n, p, k)
    g = greville(knots, p)
    for u in np.linspace(-2.0, 3.0, 41):
        N = all_basis(knots, p, u, 1)
        np.testing.assert_allclose(N[0] @ g, u, atol=1e-13)
        np.testing.assert_allclose(N[1] @ g, 1.0, atol=1e-12)  # det J = 1


@pytest.mark.parametrize("p,k,n", CASES)
def test_bezier_extraction_reproduces_bsplines(p, k, n):
    knots = uniform_clamped(0.0, 2.0, n, p, k)
    C = extraction_1d(knots, p)
    breaks = np.unique(knots)
    assert C.shape == (len(breaks) - 1, p + 1, p + 1)
    for e, (a, b) in enumerate(zip(breaks[:-1], breaks[1:])):
        for t in np.linspace(0.0, 1.0, 9)[1:-1]:
            u = a + t * (b - a)
            span, ders = basis_1d(knots, p, u, 0)
            np.testing.assert_allclose(C[e] @ bernstein(p, t)[0], ders[0], atol=1e-14)


def test_bernstein_derivative_matches_fd():
    for p, t in itertools.product((1, 2, 3, 4), (0.13, 0.5, 0.87)):
        d = bernstein(p, t, 2)
        h = 1e-6
        fd1 = (bernstein(p, t + h)[0] - bernstein(p, t - h)[0]) / (2 * h)
        np.testing.assert_allclose(d[1], fd1, atol=1e-7)
        np.testing.assert_allclose(d[0].sum(), 1.0, atol=1e-14)
