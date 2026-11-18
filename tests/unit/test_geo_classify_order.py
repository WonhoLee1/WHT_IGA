"""Regression tests for C-order flat cell indexing in ``classify._find_cut_cells_njit``."""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.backends.mesh_engine import box_triangle_soup
from wht_iga.geometry import classify
from wht_iga.reference import hertz as H
from wht_iga.reference.box import BoxSpline
from wht_iga.solver.trimmed import build_trimmed_table


def _trimmed_volume(
    qp: dict[tuple[int, int, int], np.ndarray],
    out_cells: list[tuple[int, int, int]],
    lo: tuple[float, float, float],
    hi: tuple[float, float, float],
    n: tuple[int, int, int],
) -> float:
    """Integrate total volume across all active (IN + CUT) grid cells."""
    out = set(out_cells)
    cell_vol = float(np.prod((np.asarray(hi, float) - np.asarray(lo, float)) / np.asarray(n, float)))
    vol = 0.0
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                if (i, j, k) in out:
                    continue
                if (i, j, k) in qp:
                    vol += float(np.asarray(qp[(i, j, k)])[:, 3].sum())
                else:
                    vol += cell_vol
    return vol


@pytest.mark.parametrize(
    "n_cells, box_lo, box_hi, expected_ijk",
    [
        (
            (4, 2, 3),
            (0.01, 0.01, 0.01),
            (0.24, 0.49, 0.74),
            {(0, j, k) for j in range(2) for k in range(3)},
        ),
        (
            (2, 3, 5),
            (0.26, 0.01, 0.26),
            (0.49, 0.74, 1.24),
            {(1, j, k) for j in range(3) for k in range(1, 5)},
        ),
    ],
)
def test_sat_only_box_c_order_noncubic_grids(
    n_cells: tuple[int, int, int],
    box_lo: tuple[float, float, float],
    box_hi: tuple[float, float, float],
    expected_ijk: set[tuple[int, int, int]],
) -> None:
    """Interior box hitting zero grid vertices (SAT-only) marks exact C-order cells."""
    h = np.array([0.25, 0.25, 0.25], dtype=float)
    lo = np.zeros(3, dtype=float)
    hi = tuple(0.25 * v for v in n_cells)
    nodes, tris = box_triangle_soup(box_lo, box_hi)

    cut = classify.find_cut_cells(nodes, tris, lo, h, n_cells)
    nx, ny, nz = n_cells
    expected_flat = sorted((i * ny + j) * nz + k for i, j, k in expected_ijk)
    assert cut.shape == (nx * ny * nz,)
    assert np.flatnonzero(cut).tolist() == expected_flat

    res = classify.classify_cells(nodes, tris, lo, h, n_cells)
    assert not res["vertex_inside"].any()
    state_grid = res["state"].reshape(n_cells)
    cut_ijk = {tuple(int(v) for v in row) for row in np.argwhere(state_grid == classify.CUT)}
    assert cut_ijk == expected_ijk
    assert not np.any(state_grid == classify.IN)

    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=hi, n_cells=n_cells, degree=2)
    qp, out_cells, _issues, _rec, diag = build_trimmed_table(
        nodes, tris, box, q=3, cache_dir="off"
    )
    assert set(qp) == expected_ijk
    assert diag["n_fallback"] == 0
    v_quad = _trimmed_volume(qp, out_cells, (0.0, 0.0, 0.0), hi, n_cells)
    v_exact = float(np.prod(np.subtract(box_hi, box_lo)))
    assert v_quad == pytest.approx(v_exact, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("h", [0.5, 0.25])
def test_half_hollow_hemisphere_trimmed_volume(h: float) -> None:
    """Non-symmetric half V08 indenter grid (nx != nz) recovers exact soup volume."""
    center = (-1.5, 1.525, 0.0)
    r_out, r_in = 0.5, 0.25
    nodes, tris, _ = H.hollow_hemisphere_soup(r_out, r_in, center, n=16, sym="half")
    lo, hi, n_cells = H.v08_indenter_grid(h, center, r_out, sym="half")
    assert n_cells[0] != n_cells[2]

    box = BoxSpline(lo=lo, hi=hi, n_cells=n_cells, degree=2)
    qp, out_cells, _issues, _rec, diag = build_trimmed_table(
        nodes, tris, box, q=4, cache_dir="off"
    )
    assert diag["n_fallback"] == 0
    v_quad = _trimmed_volume(qp, out_cells, lo, hi, n_cells)
    v_soup = H.soup_signed_volume(nodes, tris)
    assert v_quad == pytest.approx(v_soup, rel=1e-8)


def test_off_center_box_nonuniform_noncubic_grid_volume() -> None:
    """Off-center box on a (3, 5, 4) grid with anisotropic cell sizes recovers exact volume."""
    n_cells = (3, 5, 4)
    lo = (0.0, 0.0, 0.0)
    hi = (0.9, 1.0, 1.0)  # h = (0.3, 0.2, 0.25)
    box_lo = (0.11, 0.23, 0.08)
    box_hi = (0.79, 0.87, 0.92)
    nodes, tris = box_triangle_soup(box_lo, box_hi)

    box = BoxSpline(lo=lo, hi=hi, n_cells=n_cells, degree=2)
    qp, out_cells, _issues, _rec, diag = build_trimmed_table(
        nodes, tris, box, q=3, cache_dir="off"
    )
    assert diag["n_fallback"] == 0
    v_quad = _trimmed_volume(qp, out_cells, lo, hi, n_cells)
    v_exact = float(np.prod(np.subtract(box_hi, box_lo)))
    assert v_quad == pytest.approx(v_exact, rel=1e-12, abs=1e-12)
