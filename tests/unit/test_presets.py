# -*- coding: utf-8 -*-
"""Presets wiring tests (Kratos-free)."""

from __future__ import annotations


def test_contact_iga_preset_wires_recommendation():
    from wht_iga.analyses import presets as P

    problem, opt, rec = P.contact_iga(2.0e5, 0.0)
    assert rec.degree == 2 and not rec.capped
    assert opt.kind == "full_iga_body_fit" and opt.degree == 2
    assert problem.E == 2.0e5 and problem.traction == (2.0, 0.0, 0.0)


def test_bending_thb_preset_selects_thb_p2():
    from wht_iga.analyses import presets as P

    problem, opt, rec = P.bending_thb(2.0e11, 0.0)
    assert rec.degree == 2
    assert opt.kind == "thb_adaptive_iga" and opt.refinement_level == 1


def test_bulk_linear_preset_caps_with_linear_fe():
    from wht_iga.analyses import presets as P

    _, opt, rec = P.bulk_linear(1.0e5, 0.0, fe_type="hex8", degree=3)
    assert rec.degree == 2 and rec.capped is True
    assert opt.kind == "native_fe"


def test_describe_mentions_method_and_degree():
    from wht_iga.analyses import presets as P

    s = P.describe(P.contact_iga(2.0e5, 0.0))
    assert "FullIgaOption" in s and "p=2" in s
