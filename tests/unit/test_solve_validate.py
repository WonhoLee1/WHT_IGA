"""Unit tests for solver.validate (Kratos-free)."""

from __future__ import annotations

from wht_iga.solver.validate import DirichletBC, JobSpec, ValidationError, validate


def _clean() -> JobSpec:
    return JobSpec(name="t", lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0), n_cells=(2, 2, 2),
                   degree=2, E=1.0, nu=0.0,
                   bcs=[DirichletBC(face=(2, 0))])


def test_clean_spec_passes() -> None:
    assert validate(_clean(), []) == []


def test_no_dirichlet_blocks() -> None:
    spec = _clean()
    spec.bcs = []
    issues = validate(spec, [])
    assert any(i.severity == "error" and i.code == "NO_DIRICHLET" for i in issues)


def test_watertight_and_volume_sign_block() -> None:
    geo = [{"severity": "error", "code": "WATERTIGHT", "message": "open shell", "entity": "cell:3"},
           {"severity": "error", "code": "VOLUME_SIGN", "message": "negative vol", "entity": "cell:3"}]
    codes = [i.code for i in validate(_clean(), geo) if i.severity == "error"]
    assert "WATERTIGHT" in codes and "VOLUME_SIGN" in codes


def test_quad_fallback_face_entity_blocks_on_loaded_face() -> None:
    from wht_iga.solver.validate import PressureLoad

    spec = _clean()
    spec.pressures = [PressureLoad(face=(2, 1), magnitude=1.0)]
    loaded_geo = [{"severity": "warning", "code": "QUAD_FALLBACK", "message": "octree", "entity": "face:2:1"}]
    other_geo = [{"severity": "warning", "code": "QUAD_FALLBACK", "message": "octree", "entity": "face:0:0"}]
    assert any(i.severity == "error" and i.code == "QUAD_FALLBACK_LOADED" for i in validate(spec, loaded_geo))
    got = validate(spec, other_geo)
    assert not any(i.severity == "error" for i in got)
    assert any(i.code == "QUAD_FALLBACK" and i.severity == "warning" for i in got)


def test_quad_fallback_cell_entity_adjacency() -> None:
    """Cell-tagged QUAD_FALLBACK blocks when the cell touches a loaded face."""
    from wht_iga.solver.validate import PressureLoad

    spec = _clean()
    spec.pressures = [PressureLoad(face=(0, 1), magnitude=1.0)]
    # cell (1,0,0) of a (2,2,2) grid touches the +x loaded face
    touching = [{"severity": "warning", "code": "QUAD_FALLBACK", "message": "octree",
                 "entity": "cell:1,0,0"}]
    assert any(i.severity == "error" and i.code == "QUAD_FALLBACK_LOADED"
               for i in validate(spec, touching))
    # interior cell: warning only
    interior = [{"severity": "warning", "code": "QUAD_FALLBACK", "message": "octree",
                 "entity": "cell:0,0,0"}]
    got = validate(spec, interior)
    assert any(i.code == "QUAD_FALLBACK" and i.severity == "warning" for i in got)
    assert not any(i.code == "QUAD_FALLBACK_LOADED" for i in got)


def test_all_false_bcs_count_as_no_dirichlet() -> None:
    from wht_iga.solver.validate import DirichletBC

    spec = _clean()
    spec.bcs = [DirichletBC(face=(0, 0), components=(False, False, False))]
    assert any(i.code == "NO_DIRICHLET" and i.severity == "error"
               for i in validate(spec, []))


def test_bad_material_and_beta_block() -> None:
    spec = _clean()
    spec.E = -1.0
    spec.nu = 0.5
    spec.beta = 0.0
    codes = {i.code for i in validate(spec, []) if i.severity == "error"}
    assert {"MAT_E", "MAT_NU", "BETA_VALUE"} <= codes


def test_gravity_without_density_warns() -> None:
    spec = _clean()
    spec.gravity = (0.0, 0.0, -9.81)
    assert any(i.code == "GRAVITY_NO_DENSITY" for i in validate(spec, []))


def test_json_roundtrip() -> None:
    spec = _clean()
    spec2 = JobSpec.from_json(spec.to_json())
    assert spec2.to_json() == spec.to_json()


def test_validation_error_carries_blockers() -> None:
    spec = _clean()
    spec.bcs = []
    try:
        raise ValidationError(validate(spec, []))
    except ValidationError as e:
        assert any(i.code == "NO_DIRICHLET" for i in e.issues)


def test_out_mode_default_roundtrip_and_reject() -> None:
    spec = _clean()
    assert spec.out_mode == "ghost"  # V02 evidence ran with the FCM ghost fill
    spec.out_mode = "drop"
    back = JobSpec.from_json(spec.to_json())
    assert back.out_mode == "drop"
    assert validate(back, []) == []
    spec.out_mode = "zero"
    issues = validate(spec, [])
    assert any(i.severity == "error" and i.code == "OUT_MODE" for i in issues)
