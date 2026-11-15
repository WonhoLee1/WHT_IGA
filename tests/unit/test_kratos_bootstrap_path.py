"""configure_paths must put the Kratos DLL directories first on PATH.

Regression: a session whose PATH lacked the build deps made the first Pardiso
call exit with 0xC06D007F (delay-load procedure not found), silently.
"""

import os

from wht_iga.kratos_backend import _bootstrap as B


def test_dll_dirs_lead_path(monkeypatch, tmp_path):
    root = tmp_path / "install"
    (root / "KratosMultiphysics").mkdir(parents=True)
    (root / "KratosMultiphysics" / "__init__.py").write_text("")
    (root / "libs").mkdir()
    deps = tmp_path / "deps"
    deps.mkdir()
    other = str(tmp_path / "other")
    monkeypatch.setenv("WHT_IGA_KRATOS_PATH", str(root))
    monkeypatch.setenv("WHT_IGA_KRATOS_DLL_PATH", f"{root / 'libs'};{deps}")
    mkl = tmp_path / "mkl"
    mkl.mkdir()
    (mkl / "mkl_rt.3.dll").write_bytes(b"")
    monkeypatch.setenv("WHT_IGA_MKL_DIR", str(mkl))
    monkeypatch.setenv("PATH", os.pathsep.join([other, str(deps)]))
    monkeypatch.setattr(B, "_configured", False)
    monkeypatch.setattr(B, "_dll_handles", [])
    monkeypatch.setattr(B.sys, "path", list(B.sys.path))
    # setenv first so teardown restores the original state (delenv of an
    # absent key records nothing, and configure_paths sets it)
    monkeypatch.setenv("PYPARDISO_MKL_RT", "x")
    monkeypatch.delenv("PYPARDISO_MKL_RT")
    loaded = []
    import ctypes

    monkeypatch.setattr(ctypes, "CDLL", lambda path: loaded.append(path) or object())
    B.configure_paths()
    parts = os.environ["PATH"].split(os.pathsep)
    # one MKL per process: the env MKL dir leads, then the Kratos DLL dirs
    assert parts[:4] == [str(mkl), str(root / "libs"), str(deps), other]
    assert parts.count(str(deps)) == 1
    if os.name == "nt":  # preloaded by full path, and pypardiso pointed at the same file
        rt = os.path.join(str(mkl), "mkl_rt.3.dll")
        assert loaded == [rt]
        assert os.environ["PYPARDISO_MKL_RT"] == rt


def test_mkl_dir_default_is_env_and_missing_is_none(monkeypatch, tmp_path):
    monkeypatch.delenv("WHT_IGA_MKL_DIR", raising=False)
    d = B.mkl_dir()
    if d is not None:  # whtiga env ships mkl_rt.3.dll under Library/bin
        assert os.path.normcase(d) == os.path.normcase(os.path.join(B.sys.prefix, "Library", "bin"))
    monkeypatch.setenv("WHT_IGA_MKL_DIR", str(tmp_path))  # no mkl_rt.3.dll there
    assert B.mkl_dir() is None


def test_conflicting_pypardiso_mkl_is_overridden_with_warning(monkeypatch, tmp_path):
    import ctypes
    import warnings

    root = tmp_path / "install"
    (root / "KratosMultiphysics").mkdir(parents=True)
    (root / "KratosMultiphysics" / "__init__.py").write_text("")
    mkl = tmp_path / "mkl"
    mkl.mkdir()
    (mkl / "mkl_rt.3.dll").write_bytes(b"")
    monkeypatch.setenv("WHT_IGA_KRATOS_PATH", str(root))
    monkeypatch.setenv("WHT_IGA_KRATOS_DLL_PATH", str(root))
    monkeypatch.setenv("WHT_IGA_MKL_DIR", str(mkl))
    monkeypatch.setenv("PYPARDISO_MKL_RT", str(tmp_path / "other" / "mkl_rt.3.dll"))
    monkeypatch.setenv("PATH", os.environ.get("PATH", ""))
    monkeypatch.setattr(B, "_configured", False)
    monkeypatch.setattr(B, "_dll_handles", [])
    monkeypatch.setattr(B.sys, "path", list(B.sys.path))
    monkeypatch.setattr(ctypes, "CDLL", lambda path: object())
    if os.name != "nt":
        return
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        B.configure_paths()
    assert os.environ["PYPARDISO_MKL_RT"] == os.path.join(str(mkl), "mkl_rt.3.dll")
    assert any("PYPARDISO_MKL_RT" in str(w.message) for w in rec)
