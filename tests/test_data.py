"""Offline tests for nodal.data: a tiny fake library stands in for the NNDC."""

import hashlib
import io
import json
import tarfile
from pathlib import Path
from typing import ClassVar

import pytest

from nodal import data
from nodal.data import EndfLibrary, _neutron_file_stem

FILES = {
    "neutrons": {
        "n-001_H_001.endf": b"hydrogen\n",
        "n-040_Zr_090.endf": b"zr90\n",
        "n-040_Zr_091.endf": b"zr91\n",
        "n-092_U_238.endf": b"u238\n",
        "n-095_Am_242m1.endf": b"am242m\n",
        "README.md": b"readme\n",
    },
    "thermal_scatt": {"tsl-HinH2O.endf": b"h in water\n", "tsl-HinH2O.leapr": b"input\n"},
}


@pytest.fixture
def lib(tmp_path: Path) -> EndfLibrary:
    """A data directory that already holds metadata and archives (no downloads)."""
    subs = {}
    for sub, files in FILES.items():
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            for name, content in files.items():
                info = tarfile.TarInfo(f"{sub}-version.VIII.1/{name}")
                info.size = len(content)
                tar.addfile(info, io.BytesIO(content))
        archive = f"{sub}-version.VIII.1.tar.gz"
        (tmp_path / "archives").mkdir(exist_ok=True)
        (tmp_path / "archives" / archive).write_bytes(buf.getvalue())
        subs[sub] = {"files": {"ENDF-6": {
            "name": archive, "size": "1 KB",
            "checksums": {"SHA1": hashlib.sha1(buf.getvalue()).hexdigest()}}}}
    (tmp_path / "metadata.json").write_text(json.dumps({"sublibraries": subs}))
    return EndfLibrary(tmp_path)


@pytest.mark.parametrize("name, key", [
    ("U238", ("U", 238, "")), ("u-238", ("U", 238, "")), ("H1", ("H", 1, "")),
    ("Am242m", ("Am", 242, "m1")), ("am242m1", ("Am", 242, "m1")),
])
def test_name_parsing(name, key):
    assert _neutron_file_stem(name) == key


def test_bad_name():
    with pytest.raises(ValueError):
        _neutron_file_stem("uranium")


def test_neutron_lookup_extracts_file(lib):
    path = lib.neutron("U-238")
    assert path == lib.root / "neutrons" / "n-092_U_238.endf"
    assert path.read_bytes() == b"u238\n"


def test_batch_and_metastable(lib):
    paths = lib.neutrons(["H1", "Am242m"])
    assert [p.read_bytes() for p in paths] == [b"hydrogen\n", b"am242m\n"]


def test_isotopes_and_missing_natural_element(lib):
    assert lib.isotopes("zr") == ["Zr90", "Zr91"]
    with pytest.raises(KeyError, match="Zr90"):
        lib.neutron("Zr92")


def test_thermal(lib):
    assert lib.thermal_materials() == ["HinH2O"]
    assert lib.thermal("hinh2o").read_bytes() == b"h in water\n"
    with pytest.raises(KeyError):
        lib.thermal("HinD2O")


def test_download_verifies_sha1(tmp_path, monkeypatch) -> None:
    payload = b"x" * 100

    class FakeResponse(io.BytesIO):
        headers: ClassVar[dict[str, str]] = {"Content-Length": str(len(payload))}

    monkeypatch.setattr(data.urllib.request, "urlopen", lambda req: FakeResponse(payload))
    good = tmp_path / "good.bin"
    data._download("https://example.invalid/f", good, sha1=hashlib.sha1(payload).hexdigest())
    assert good.read_bytes() == payload

    bad = tmp_path / "bad.bin"
    with pytest.raises(OSError, match="SHA-1 mismatch"):
        data._download("https://example.invalid/f", bad, sha1="0" * 40)
    assert not bad.exists() and not bad.with_name("bad.bin.part").exists()
