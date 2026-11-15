# -*- coding: utf-8 -*-
"""F5: monotonic ID allocator never collides, even at 40k+ constraints."""

from __future__ import annotations

from wht_iga.solver.pipeline import _next_id


class _O:
    def __init__(self, i: int):
        self.Id = i


def test_empty_container_returns_start() -> None:
    assert _next_id([], 100) == 100
    assert _next_id(iter([]), 7) == 7


def test_monotonic_past_max() -> None:
    objs = [_O(i) for i in [1, 5, 9]]
    assert _next_id(objs, 100) == 10


def test_no_collision_at_scale() -> None:
    """Simulate BC + load + extension allocations on one model part."""
    used: set[int] = set()
    allocs: dict[str, int] = {}
    # control nodes 1..N already exist; each block starts past the max
    nodes = [_O(i) for i in range(1, 101)]
    nxt = _next_id(nodes, 1)
    for n_new in (30, 40000, 7):  # BC-B virtuals, ties, loads
        ids = list(range(nxt, nxt + n_new))
        assert not (set(ids) & used), "collision"
        used.update(ids)
        nodes.extend(_O(i) for i in ids)
        nxt = _next_id(nodes, nxt)
    assert len(used) == 40037
