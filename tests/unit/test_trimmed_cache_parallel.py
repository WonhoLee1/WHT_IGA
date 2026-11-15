"""build_trimmed_table: process-pool path and content-hash cache give identical results."""

from __future__ import annotations

import os

import numpy as np

from wht_iga.reference.box import BoxSpline
from wht_iga.solver import trimmed as T


def _box_surface(lo=0.23, hi=0.81):
    """Closed, outward-wound triangulated cube [lo, hi]^3 (12 triangles)."""
    v = np.array([[x, y, z] for z in (lo, hi) for y in (lo, hi) for x in (lo, hi)], float)
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    tris = []
    for a, b, c, d in quads:
        tris += [(a, b, c), (a, c, d)]
    return v, np.array(tris, int)


def _same(a, b):
    assert a[0].keys() == b[0].keys()
    assert all(np.array_equal(a[0][k], b[0][k]) for k in a[0])
    assert a[1] == b[1] and a[2] == b[2] and a[3] == b[3] and a[4] == b[4]


def _box():
    return BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(3, 3, 3), degree=2)


def test_pool_matches_serial(monkeypatch):
    nodes, tris = _box_surface()
    ser = T.build_trimmed_table(nodes, tris, _box(), q=2, workers=1, cache_dir="off")
    assert ser[4]["n_cut"] > 0
    monkeypatch.setattr(T, "_PARALLEL_MIN_CELLS", 1)
    par = T.build_trimmed_table(nodes, tris, _box(), q=2, workers=2, cache_dir="off")
    _same(ser, par)


def test_cache_hit_identical_and_keyed(tmp_path):
    nodes, tris = _box_surface()
    ref = T.build_trimmed_table(nodes, tris, _box(), q=2, workers=1, cache_dir="off")
    d = str(tmp_path)
    first = T.build_trimmed_table(nodes, tris, _box(), q=2, workers=1, cache_dir=d)
    assert len(os.listdir(d)) == 1
    second = T.build_trimmed_table(nodes, tris, _box(), q=2, workers=1, cache_dir=d)
    _same(ref, first)
    _same(ref, second)
    T.build_trimmed_table(nodes, tris, _box(), q=3, workers=1, cache_dir=d)  # new key
    nodes2 = nodes.copy()
    nodes2[0, 0] += 1e-3
    T.build_trimmed_table(nodes2, tris, _box(), q=2, workers=1, cache_dir=d)  # new key
    assert len([f for f in os.listdir(d) if f.endswith(".pkl")]) == 3


def test_cache_off_writes_nothing(tmp_path, monkeypatch):
    nodes, tris = _box_surface()
    monkeypatch.setenv("WHT_GEO_CACHE", "off")
    T.build_trimmed_table(nodes, tris, _box(), q=2, workers=1)
    assert T._cache_dir(None) is None
