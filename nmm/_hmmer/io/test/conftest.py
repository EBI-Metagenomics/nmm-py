import pytest


@pytest.fixture
def tblout(tmp_path):
    return _write_file(tmp_path, "tblout")


def _write_file(path, filename):
    import importlib_resources as pkg_resources
    import nmm

    text = pkg_resources.read_text(nmm._hmmer.io.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
