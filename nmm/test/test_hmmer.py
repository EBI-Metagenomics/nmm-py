import pytest
import importlib_resources as pkg_resources
from numpy.testing import assert_equal, assert_allclose

from nmm import read_hmmer


def test_read_hmmer(PF03373_path):
    hmmer = read_hmmer(PF03373_path)


@pytest.fixture
def PF03373_path(tmp_path):
    import nmm

    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    return tmp_path / filename
