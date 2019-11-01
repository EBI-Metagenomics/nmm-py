import importlib_resources as pkg_resources
from numpy.testing import assert_allclose

from nmm import read_hmmer


def test_read_hmmer(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    hmmer = read_hmmer(filepath)
    most_likely_seq = "PGKEDNNK"
    lr = hmmer.lr(most_likely_seq)
    assert_allclose(lr.score, 11.867796719423442)


def write_file(path, filename):
    import nmm

    text = pkg_resources.read_text(nmm.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
