"""Mesh + results SQLite DBs per v2 §4.2 / §4.5 (mesh stage subset) (Track SOLVE).

BLOB convention (v2 §4.1): little-endian float64/int64 ``ndarray.tobytes()``,
array shape stored in a sibling column. Read-back is bit-exact (roundtrip
test asserts exact equality). ``meta`` carries schema/provenance; every
``validate()`` issue is stored in ``issue``.
"""

from __future__ import annotations

import json
import sqlite3

import numpy as np

MESH_SCHEMA_VERSION = "wiga.mesh/1"
RESULTS_SCHEMA_VERSION = "wiga.results/1"

_MESH_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE grid_level (level INTEGER PRIMARY KEY, degree INTEGER, continuity INTEGER,
                         lo BLOB, hi BLOB, n_cells BLOB);
CREATE TABLE cell (cell_id INTEGER PRIMARY KEY, level INTEGER, i INTEGER, j INTEGER, k INTEGER,
                   status INTEGER, eta REAL, n_qp INTEGER);
CREATE TABLE qp_volume (cell_id INTEGER, material INTEGER, block INTEGER, shape TEXT, xw BLOB);
CREATE TABLE qp_surface (face_id INTEGER, cell_id INTEGER, block INTEGER, shape TEXT, xwn BLOB);
CREATE TABLE probe (probe_id INTEGER PRIMARY KEY, name TEXT, kind TEXT, xyz BLOB, variables TEXT);
CREATE TABLE issue (severity TEXT, code TEXT, message TEXT, entity TEXT);
"""

_RESULTS_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE step (step_id INTEGER PRIMARY KEY, name TEXT, kind TEXT, t_end REAL);
CREATE TABLE frame (frame_id INTEGER PRIMARY KEY, step_id INTEGER, time REAL,
                    increment INTEGER, converged INTEGER, n_iter INTEGER);
CREATE TABLE cp_field (frame_id INTEGER, var TEXT, shape TEXT, data BLOB);
CREATE TABLE qp_field (frame_id INTEGER, var TEXT, shape TEXT, data BLOB);
"""


def encode(a: np.ndarray) -> bytes:
    return np.ascontiguousarray(a, dtype=np.float64).tobytes()


def decode(blob: bytes, shape: tuple[int, ...]) -> np.ndarray:
    return np.frombuffer(blob, dtype=np.float64).reshape(shape).copy()


def _write_meta(con: sqlite3.Connection, meta: dict) -> None:
    con.executemany("INSERT INTO meta (key, value) VALUES (?, ?)", list(meta.items()))


def read_meta(path: str) -> dict:
    con = sqlite3.connect(path)
    try:
        return dict(con.execute("SELECT key, value FROM meta").fetchall())
    finally:
        con.close()


def write_mesh_db(
    path: str,
    spec_json: dict,
    meta_extra: dict | None = None,
    degree: int = 1,
    continuity: int | None = None,
    lo: tuple = (0.0, 0.0, 0.0),
    hi: tuple = (1.0, 1.0, 1.0),
    n_cells: tuple = (1, 1, 1),
    cells: list[dict] | None = None,
    qp_volumes: dict[int, np.ndarray] | None = None,
    qp_surfaces: dict[tuple[int, int], np.ndarray] | None = None,
    probes: list[dict] | None = None,
    issues: list[dict] | None = None,
) -> str:
    """Write ``<job>.wiga.sqlite`` mesh stage (v2 §4.2 subset for Phase 1)."""
    import os

    if os.path.exists(path):
        os.remove(path)
    con = sqlite3.connect(path)
    try:
        con.executescript(_MESH_SCHEMA)
        meta = {"schema_version": MESH_SCHEMA_VERSION, "spec_json": json.dumps(spec_json)}
        meta.update(meta_extra or {})
        _write_meta(con, meta)
        con.execute(
            "INSERT INTO grid_level (level, degree, continuity, lo, hi, n_cells) VALUES (0, ?, ?, ?, ?, ?)",
            (int(degree), -1 if continuity is None else int(continuity),
             encode(np.asarray(lo, dtype=float)), encode(np.asarray(hi, dtype=float)),
             encode(np.asarray(n_cells, dtype=float))),
        )
        nx, ny, nz = (int(v) for v in n_cells)
        for c in cells or []:
            ci, cj, ck = int(c["i"]), int(c["j"]), int(c["k"])
            con.execute(
                "INSERT INTO cell (cell_id, level, i, j, k, status, eta, n_qp) VALUES (?, 0, ?, ?, ?, ?, ?, ?)",
                ((ck * ny + cj) * nx + ci, ci, cj, ck, int(c.get("status", 1)),
                 float(c.get("eta", 1.0)), int(c.get("n_qp", 0))),
            )
        for cell_id, xw in (qp_volumes or {}).items():
            xw = np.asarray(xw, dtype=float)
            con.execute("INSERT INTO qp_volume (cell_id, material, block, shape, xw) VALUES (?, 0, 0, ?, ?)",
                        (int(cell_id), json.dumps(list(xw.shape)), encode(xw)))
        for (axis, side), xwn in (qp_surfaces or {}).items():
            xwn = np.asarray(xwn, dtype=float)
            con.execute("INSERT INTO qp_surface (face_id, cell_id, block, shape, xwn) VALUES (?, -1, 0, ?, ?)",
                        (int(axis * 2 + side), json.dumps(list(xwn.shape)), encode(xwn)))
        for k, p in enumerate(probes or []):
            con.execute("INSERT INTO probe (probe_id, name, kind, xyz, variables) VALUES (?, ?, ?, ?, ?)",
                        (k, p.get("name", f"p{k}"), p.get("kind", "point"),
                         encode(np.asarray(p.get("xyz", (0.0, 0.0, 0.0)), dtype=float)),
                         json.dumps(p.get("variables", []))))
        for i in issues or []:
            con.execute("INSERT INTO issue (severity, code, message, entity) VALUES (?, ?, ?, ?)",
                        (i.get("severity", "info"), i.get("code", ""), i.get("message", ""), i.get("entity", "")))
        con.commit()
    finally:
        con.close()
    return path


def read_mesh_db(path: str) -> dict:
    """Read back all mesh-stage tables (exact arrays)."""
    con = sqlite3.connect(path)
    try:
        con.row_factory = sqlite3.Row
        grid = con.execute("SELECT * FROM grid_level").fetchone()
        cells = [dict(r) for r in con.execute("SELECT * FROM cell ORDER BY rowid")]
        qp_vol = [(r["cell_id"], decode(r["xw"], tuple(json.loads(r["shape"]))))
                  for r in con.execute("SELECT * FROM qp_volume")]
        qp_surf = [(r["face_id"], decode(r["xwn"], tuple(json.loads(r["shape"]))))
                   for r in con.execute("SELECT * FROM qp_surface")]
        probes = [dict(r) for r in con.execute("SELECT * FROM probe ORDER BY probe_id")]
        issues = [dict(r) for r in con.execute("SELECT * FROM issue")]
        return {"meta": read_meta(path), "grid": dict(grid), "cells": cells,
                "qp_volumes": qp_vol, "qp_surfaces": qp_surf, "probes": probes, "issues": issues}
    finally:
        con.close()


def write_results_db(
    path: str,
    spec_json: dict,
    meta_extra: dict | None = None,
    steps: list[dict] | None = None,
    frames: list[dict] | None = None,
    cp_fields: dict[tuple[int, str], np.ndarray] | None = None,
    qp_fields: dict[tuple[int, str], np.ndarray] | None = None,
) -> str:
    """Write ``<job>_iga_results.sqlite`` (v2 §4.3)."""
    import os

    if os.path.exists(path):
        os.remove(path)
    con = sqlite3.connect(path)
    try:
        con.executescript(_RESULTS_SCHEMA)
        meta = {"schema_version": RESULTS_SCHEMA_VERSION, "spec_json": json.dumps(spec_json)}
        meta.update(meta_extra or {})
        _write_meta(con, meta)
        for s in steps or [{"step_id": 0, "name": "load", "kind": "static", "t_end": 1.0}]:
            con.execute("INSERT INTO step (step_id, name, kind, t_end) VALUES (?, ?, ?, ?)",
                        (s["step_id"], s["name"], s["kind"], s["t_end"]))
        # NOTE: default converged=0 means "not recorded" (unknown), never
        # "converged". Callers with a real solver verdict must pass frames
        # explicitly (pipeline.run_linear does).
        for f in frames or [{"frame_id": 0, "step_id": 0, "time": 1.0, "increment": 0, "converged": 0, "n_iter": 0}]:
            con.execute("INSERT INTO frame (frame_id, step_id, time, increment, converged, n_iter) VALUES (?, ?, ?, ?, ?, ?)",
                        (f["frame_id"], f["step_id"], f["time"], f["increment"], f["converged"], f["n_iter"]))
        for (frame_id, var), data in (cp_fields or {}).items():
            data = np.asarray(data, dtype=float)
            con.execute("INSERT INTO cp_field (frame_id, var, shape, data) VALUES (?, ?, ?, ?)",
                        (int(frame_id), var, json.dumps(list(data.shape)), encode(data)))
        for (frame_id, var), data in (qp_fields or {}).items():
            data = np.asarray(data, dtype=float)
            con.execute("INSERT INTO qp_field (frame_id, var, shape, data) VALUES (?, ?, ?, ?)",
                        (int(frame_id), var, json.dumps(list(data.shape)), encode(data)))
        con.commit()
    finally:
        con.close()
    return path


def read_results_db(path: str) -> dict:
    con = sqlite3.connect(path)
    try:
        con.row_factory = sqlite3.Row
        steps = [dict(r) for r in con.execute("SELECT * FROM step")]
        frames = [dict(r) for r in con.execute("SELECT * FROM frame")]
        cp = {(r["frame_id"], r["var"]): decode(r["data"], tuple(json.loads(r["shape"])))
              for r in con.execute("SELECT * FROM cp_field")}
        qp = {(r["frame_id"], r["var"]): decode(r["data"], tuple(json.loads(r["shape"])))
              for r in con.execute("SELECT * FROM qp_field")}
        return {"meta": read_meta(path), "steps": steps, "frames": frames, "cp_fields": cp, "qp_fields": qp}
    finally:
        con.close()
