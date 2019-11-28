import pytest


@pytest.fixture
def PF03373(tmp_path):
    return _write_file(tmp_path, "PF03373.hmm")


def _write_file(path, filename):
    import importlib_resources as pkg_resources
    import nmm

    text = pkg_resources.read_text(nmm._hmmer.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
