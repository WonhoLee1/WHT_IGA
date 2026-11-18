# -*- coding: utf-8 -*-
"""Regression: ``solver.trimmed._cell_rules`` serial path is re-entrant.

Review 2026-11-17 (dev_log/track2_review_20261116_opus.md, N1): the serial
branch passed geometry through a module-global dict, so two threads (e.g.
``compat.coreform.Job.submit``) could compute cut-cell rules with each
other's surface. Reference: the same calls made sequentially.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np

from wht_iga.backends.mesh_engine import box_triangle_soup
from wht_iga.solver import trimmed as TR


def _case(lo, hi):
    nodes, tris = box_triangle_soup(np.array(lo, float), np.array(hi, float))
    h = np.array([0.5, 0.5, 0.5])
    blos = np.array([[i * 0.5, j * 0.5, k * 0.5]
                     for k in range(2) for j in range(2) for i in range(2)])
    return nodes, tris, blos, h


def _rules(case):
    nodes, tris, blos, h = case
    return TR._cell_rules(nodes, tris, blos, h, q=2, p=1, workers=1)


def _same(a, b):
    return len(a) == len(b) and all(
        np.array_equal(ra[0], rb[0]) and np.array_equal(ra[1], rb[1]) for ra, rb in zip(a, b))


def test_serial_cell_rules_have_no_shared_state():
    assert not hasattr(TR, "_W")
    case_a = _case((0.05, 0.05, 0.05), (0.95, 0.95, 0.95))
    case_b = _case((0.10, 0.20, 0.30), (0.70, 0.80, 0.90))
    ref_a, ref_b = _rules(case_a), _rules(case_b)
    assert not _same(ref_a, ref_b)  # the two geometries give different rules
    jobs = [case_a, case_b] * 12
    with ThreadPoolExecutor(max_workers=4) as ex:
        outs = list(ex.map(_rules, jobs))
    for case, out in zip(jobs, outs):
        assert _same(out, ref_a if case is case_a else ref_b)
