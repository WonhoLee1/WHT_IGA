# -*- coding: utf-8 -*-
"""PoC 06: mesh_data.sqlite BLOB 버스 + bezier_grid.mdpa 라운드트립.

master-plan §2.1 스키마 그대로 생성, 요소별 가변 C_tilde BLOB·적분점 기록 후
역조회 일치 확인. 표준 sqlite3+numpy만 사용.
"""
import os
import sqlite3
import numpy as np

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mesh_data.sqlite")
MDPA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bezier_grid.mdpa")


def build_db():
    if os.path.exists(DB):
        os.remove(DB)
    con = sqlite3.connect(DB)
    cur = con.cursor()
    cur.execute("""CREATE TABLE elements (
        elem_id INTEGER PRIMARY KEY, status INTEGER NOT NULL,
        num_qp INTEGER NOT NULL, c_tilde BLOB NOT NULL)""")
    cur.execute("""CREATE TABLE quadrature_points (
        elem_id INTEGER NOT NULL, qp_id INTEGER NOT NULL,
        x REAL NOT NULL, y REAL NOT NULL, z REAL NOT NULL,
        weight REAL NOT NULL, nx REAL NOT NULL, ny REAL NOT NULL,
        nz REAL NOT NULL, PRIMARY KEY (elem_id, qp_id))""")
    cur.execute("CREATE INDEX idx_qp_elem ON quadrature_points(elem_id)")
    cur.execute("CREATE TABLE model_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    rng = np.random.default_rng(7)
    # elem 1: Interior (4x4 C_tilde), elem 2: Cut (4x6 + 5qp), elem 3: Exterior
    C1 = np.eye(4)
    C2 = rng.normal(size=(4, 6)) * 0.5
    C3 = np.zeros((4, 4))
    for eid, status, C, nqp in [(1, 1, C1, 8), (2, 2, C2, 5), (3, 0, C3, 0)]:
        cur.execute("INSERT INTO elements VALUES (?,?,?,?)",
                    (eid, status, nqp, C.astype(np.float64).tobytes()))
        for q in range(nqp):
            cur.execute("INSERT INTO quadrature_points VALUES (?,?,?,?,?,?,?,?,?)",
                        (eid, q, 0.1 * q, 0.2 * q, 0.0, 0.01 * (q + 1), 0, 0, 1))
    cur.execute("INSERT INTO model_metadata VALUES (?,?)", ("degree", "[2,2,2]"))
    cur.execute("INSERT INTO model_metadata VALUES (?,?)", ("eta_crit", "1e-2"))
    con.commit()
    con.close()


def test_roundtrip():
    build_db()
    con = sqlite3.connect(DB)
    cur = con.cursor()
    n_elem = cur.execute("SELECT COUNT(*) FROM elements").fetchone()[0]
    n_qp = cur.execute("SELECT COUNT(*) FROM quadrature_points").fetchone()[0]
    blob = cur.execute("SELECT c_tilde FROM elements WHERE elem_id=2").fetchone()[0]
    C2 = np.frombuffer(blob, dtype=np.float64).reshape(4, 6)
    meta = dict(cur.execute("SELECT key, value FROM model_metadata").fetchall())
    con.close()
    print(f"[sqlite] elements={n_elem} qp_rows={n_qp} C2.shape={C2.shape} meta={meta}")
    assert n_elem == 3 and n_qp == 13
    assert C2.shape == (4, 6) and np.abs(C2).max() > 0
    assert meta["eta_crit"] == "1e-2"


def write_mdpa():
    with open(MDPA, "w", encoding="utf-8") as f:
        f.write("Begin ModelPartData\nDENSITY 7850.0\nYOUNG_MODULUS 2.1e11\n"
                "POISSON_RATIO 0.3\nEnd ModelPartData\n\nBegin Nodes\n")
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
            f.write(f" {i} {x} {y} 0.0\n")
        f.write("End Nodes\n\nBegin Elements Element3D8N\n 1 1 1 2 3 4 4 3 2 1\n"
                "End Elements\n")


def test_mdpa():
    write_mdpa()
    nodes = elems = 0
    mode = None
    for line in open(MDPA, encoding="utf-8"):
        s = line.split()
        if s[:2] == ["Begin", "Nodes"]:
            mode = "n"
        elif s[:2] == ["Begin", "Elements"]:
            mode = "e"
        elif s[:2] == ["End", "Nodes"] or s[:2] == ["End", "Elements"]:
            mode = None
        elif mode == "n" and s:
            nodes += 1
        elif mode == "e" and s:
            elems += 1
    print(f"[mdpa] nodes={nodes} elements={elems}")
    assert nodes == 4 and elems == 1


if __name__ == "__main__":
    test_roundtrip()
    test_mdpa()
    print("PoC-06 PASS")
