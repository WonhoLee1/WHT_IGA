# -*- coding: utf-8 -*-
"""verify_report determinism + artifact completeness (Kratos-free)."""

from __future__ import annotations

import json
import os


def test_verify_report_artifacts_complete(tmp_path):
    from wht_iga.post.verify_report import build_report

    out = str(tmp_path / "rep")
    manifest = build_report(out)
    for name in ("report.html", "fig_compare_accuracy.png",
                 "fig_compare_cost.png", "fig_hertz_convergence.png",
                 "fig_coreform_duel.png", "data.json"):
        p = os.path.join(out, name)
        assert os.path.isfile(p), name
        assert os.path.getsize(p) > 0, name
    assert manifest["gates_total"] == 20
    with open(os.path.join(out, "data.json"), encoding="utf-8") as fh:
        data = json.load(fh)
    assert len(data["gates"]) == 20
    assert any(g["id"] == "V08" for g in data["gates"])
    html = open(os.path.join(out, "report.html"), encoding="utf-8").read()
    assert "WHT_IGA Verification Report" in html
    assert "Anti-mocking rule" in html


def test_verify_report_deterministic(tmp_path):
    from wht_iga.post.verify_report import build_report

    a = str(tmp_path / "a")
    b = str(tmp_path / "b")
    build_report(a)
    build_report(b)
    da = json.load(open(os.path.join(a, "data.json"), encoding="utf-8"))
    db = json.load(open(os.path.join(b, "data.json"), encoding="utf-8"))
    da.pop("generated_utc", None)
    db.pop("generated_utc", None)
    assert da == db
