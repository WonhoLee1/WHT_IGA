"""Quadrature DB writer/reader per v2 §4.4 (`<job>_iga_quadrature.sqlite`).

Tables: ``meta``, ``qp`` (qp_id, cell_id, material, x, y, z, dV),
``qp_val`` (frame_id, qp_id, U, S blobs). Arrays are little-endian
float64/int64 BLOBs (v2 §4.1).
"""

from __future__ import annotations

import sqlite3

import numpy as np

SCHEMA_VERSION = "1"

_META = {"schema_version": SCHEMA_VERSION, "units": "input-consistent"}


def _blob(a: np.ndarray) -> bytes:
    return np.ascontiguousarray(a, dtype=np.float64).tobytes()


def _from_blob(b: bytes, shape: tuple[int, ...]) -> np.ndarray:
    return np.frombuffer(b, dtype=np.float64).reshape(shape)


def write_quad_db(path: str, cells: list[dict], meta: dict | None = None,
                  frame_id: int = 0) -> int:
    """Write per-cell rules to a quadrature DB. Returns total point count.

    Each cell dict: ``cell_id``, ``material``, ``x`` (N,3), ``w`` (N,),
    optional ``U`` (N,3), ``S`` (N,6) field blobs for ``qp_val``.
    """
    con = sqlite3.connect(path)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    cur.execute("""CREATE TABLE qp (qp_id INTEGER PRIMARY KEY, cell_id INTEGER,
                   material INTEGER, x REAL, y REAL, z REAL, dV REAL)""")
    cur.execute("""CREATE TABLE qp_val (frame_id INTEGER, qp_id INTEGER,
                   U BLOB, S BLOB)""")
    merged = dict(_META)
    if meta:
        merged.update({k: str(v) for k, v in meta.items()})
    cur.executemany("INSERT INTO meta VALUES (?, ?)", list(merged.items()))
    qid, n = 0, 0
    for c in cells:
        x = np.asarray(c["x"], float)
        w = np.asarray(c["w"], float)
        u = c.get("U")
        s = c.get("S")
        for k in range(len(x)):
            cur.execute("INSERT INTO qp VALUES (?, ?, ?, ?, ?, ?, ?)",
                        (qid, int(c["cell_id"]), int(c.get("material", 0)),
                         float(x[k, 0]), float(x[k, 1]), float(x[k, 2]),
                         float(w[k])))
            cur.execute("INSERT INTO qp_val VALUES (?, ?, ?, ?)",
                        (frame_id, qid,
                         _blob(np.zeros(3) if u is None else np.asarray(u)[k]),
                         _blob(np.zeros(6) if s is None else np.asarray(s)[k])))
            qid += 1
        n += len(x)
    con.commit()
    con.close()
    return n


def read_quad_db(path: str) -> dict:
    """Read back meta, qp table and qp_val blobs (roundtrip check)."""
    con = sqlite3.connect(path)
    cur = con.cursor()
    meta = dict(cur.execute("SELECT key, value FROM meta").fetchall())
    qp = cur.execute(
        "SELECT qp_id, cell_id, material, x, y, z, dV FROM qp").fetchall()
    vals = cur.execute("SELECT frame_id, qp_id, U, S FROM qp_val").fetchall()
    con.close()
    return {"meta": meta, "qp": qp, "qp_val": vals}


def qp_totals(path: str) -> tuple[float, int]:
    """(ΣdV, point count) straight from the DB."""
    con = sqlite3.connect(path)
    cur = con.cursor()
    tot, cnt = cur.execute("SELECT SUM(dV), COUNT(*) FROM qp").fetchone()
    con.close()
    return float(tot), int(cnt)
