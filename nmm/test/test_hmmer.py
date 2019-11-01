import pytest
import importlib_resources as pkg_resources
from numpy.testing import assert_allclose

from nmm import read_hmmer


def test_read_hmmer(PF03373_path):
    hmmer = read_hmmer(PF03373_path)
    most_likely_seq = "PGKEDNNK"
    score = hmmer.lr(most_likely_seq)
    assert_allclose(score, 14.300392905370323)


@pytest.fixture
def PF03373_path(tmp_path):
    import nmm

    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    return tmp_path / filename
