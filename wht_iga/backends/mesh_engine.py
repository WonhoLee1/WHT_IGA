"""Default ``MeshEnginePort`` adapter: geometry -> classification -> cut-cell quadrature.

Calls only public geometry / bridge functions (``wht_iga.geometry.occ_io.load_cad_file``,
``wht_iga.geometry.cad.tessellate``, ``wht_iga.solver.trimmed.build_trimmed_table``);
it never modifies them. Untrimmed grids skip classification entirely.
"""

from __future__ import annotations

import numpy as np

from wht_iga.ports.contracts import CancelToken, MeshModel, ProgressCallback, ProgressEvent
from wht_iga.spec.assembly import BoxGeometrySpec, StepGeometrySpec

# outward-oriented triangles of a box with vertex v(i, j, k) = i + 2 j + 4 k
_BOX_TRIS = np.array([[0, 2, 1], [1, 2, 3],   # z = lo (-z)
                      [4, 5, 6], [5, 7, 6],   # z = hi (+z)
                      [0, 1, 5], [0, 5, 4],   # y = lo (-y)
                      [2, 6, 7], [2, 7, 3],   # y = hi (+y)
                      [0, 4, 6], [0, 6, 2],   # x = lo (-x)
                      [1, 3, 7], [1, 7, 5]],  # x = hi (+x)
                     dtype=np.int64)


def box_triangle_soup(lo, hi) -> tuple[np.ndarray, np.ndarray]:
    """Closed, outward-oriented 12-triangle surface of the box ``[lo, hi]``."""
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    nodes = np.array([[(lo, hi)[i][0], (lo, hi)[j][1], (lo, hi)[k][2]]
                      for k in (0, 1) for j in (0, 1) for i in (0, 1)], dtype=float)
    return nodes, _BOX_TRIS.copy()


class GridMeshEngine:
    """Builds the :class:`MeshModel` of a compiled plan.

    ``workers`` / ``cache_dir`` are passed to ``build_trimmed_table`` (``cache_dir="off"``
    disables the on-disk trimmed-table cache).
    """

    def __init__(self, *, workers: int | None = None, cache_dir: str | None = None) -> None:
        self.workers = workers
        self.cache_dir = cache_dir

    def _surface(self, plan) -> tuple[np.ndarray, np.ndarray]:
        """Triangle soup of the (single) part."""
        geom = plan.spec.parts[0].geometry
        if isinstance(geom, BoxGeometrySpec):
            return box_triangle_soup(geom.lo, geom.hi)
        assert isinstance(geom, StepGeometrySpec)
        from wht_iga.geometry import cad, occ_io

        shape = occ_io.load_cad_file(geom.path, heal=geom.heal).shape
        g = plan.report.grid
        tess = cad.tessellate(shape, deflection=float(g.tessellation_size))
        return np.asarray(tess.nodes, float), np.asarray(tess.tris, np.int64)

    def build(self, plan, *, on_progress: ProgressCallback | None = None,
              cancel: CancelToken | None = None) -> MeshModel:
        """Classify the grid against the part and fit cut-cell rules (if trimmed)."""
        plan.report.raise_if_invalid()
        g = plan.report.grid
        emit = on_progress or (lambda ev: None)
        if not g.trimmed:
            emit(ProgressEvent("mesh", 1.0, "untrimmed grid: no classification needed"))
            return MeshModel(lo=g.lo, hi=g.hi, n_cells=g.n_cells, degree=g.degree,
                             continuity=g.continuity, trimmed=False)
        from wht_iga.reference.box import BoxSpline
        from wht_iga.solver.trimmed import build_trimmed_table

        emit(ProgressEvent("mesh", 0.0, "tessellating part"))
        nodes, tris = self._surface(plan)
        if cancel is not None:
            cancel.raise_if_cancelled("mesh: after tessellation")
        emit(ProgressEvent("mesh", 0.3, f"classifying {np.prod(g.n_cells)} cells against {len(tris)} triangles"))
        box = BoxSpline(lo=g.lo, hi=g.hi, n_cells=g.n_cells, degree=g.degree, continuity=g.continuity)
        qp_table, out_cells, issues, cell_records, diag = build_trimmed_table(
            nodes, tris, box, q=g.quadrature_degree, workers=self.workers, cache_dir=self.cache_dir)
        emit(ProgressEvent("mesh", 1.0, f"cut={diag.get('n_cut')} out={len(out_cells)} "
                                        f"fallback={diag.get('n_fallback')}"))
        return MeshModel(lo=g.lo, hi=g.hi, n_cells=g.n_cells, degree=g.degree, continuity=g.continuity,
                         trimmed=True, qp_table=qp_table,
                         out_cells=tuple(tuple(int(v) for v in c) for c in out_cells),
                         cell_records=list(cell_records), geo_issues=list(issues), diag=dict(diag),
                         tess_nodes=nodes, tess_tris=tris)
