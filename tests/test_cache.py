"""Tests for the pointwise-table disk cache."""

import pytest

from fission.data import EndfLibrary
from fission.endf import Evaluation, cache
from fission.endf.pointwise import pointwise


@pytest.fixture(scope="module")
def lib():
    lib = EndfLibrary()
    if not (lib.root / "archives" / "neutrons-version.VIII.1.tar.gz").exists():
        pytest.skip("ENDF/B-VIII.1 neutron data not downloaded")
    return lib


@pytest.fixture
def build_counter(monkeypatch):
    """Count how often a table is actually built (rather than loaded)."""
    calls = []
    real = cache.pointwise

    def counting(*args, **kwargs):
        calls.append(args)
        return real(*args, **kwargs)

    monkeypatch.setattr(cache, "pointwise", counting)
    return calls


def test_second_call_loads_identical_table(lib, tmp_path, build_counter):
    first = cache.cached_pointwise("B10", lib=lib, cache_dir=tmp_path)
    second = cache.cached_pointwise("B10", lib=lib, cache_dir=tmp_path)
    assert len(build_counter) == 1
    assert second.energies == first.energies and second.xs == first.xs
    fresh = pointwise(Evaluation.from_file(lib.neutron("B10")))
    assert second.energies == fresh.energies and second.xs == fresh.xs


def test_parameters_are_part_of_the_key(lib, tmp_path, build_counter):
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path, tol=1e-2)
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path, hi=1e6)
    assert len(build_counter) == 3
    assert len(list(tmp_path.glob("H1-*.npz"))) == 3


def test_refresh_rebuilds(lib, tmp_path, build_counter):
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path, refresh=True)
    assert len(build_counter) == 2


def test_corrupt_file_is_rebuilt(lib, tmp_path, build_counter):
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    (path,) = tmp_path.glob("H1-*.npz")
    path.write_bytes(b"not an npz file")
    table = cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    assert len(build_counter) == 2 and len(table) > 100


def test_code_change_invalidates_and_prunes(lib, tmp_path, build_counter, monkeypatch):
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    (old,) = tmp_path.glob("H1-*.npz")
    monkeypatch.setattr(cache, "code_fingerprint", lambda: "edited-code")
    cache.cached_pointwise("H1", lib=lib, cache_dir=tmp_path)
    assert len(build_counter) == 2
    assert not old.exists()  # entry from the old code was deleted
    assert len(list(tmp_path.glob("H1-*.npz"))) == 1


def test_fingerprint_covers_sources_and_extension():
    cache.code_fingerprint.cache_clear()
    assert len(cache.code_fingerprint()) == 64
