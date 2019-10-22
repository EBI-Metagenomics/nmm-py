import pytest
import importlib_resources as pkg_resources
from numpy.testing import assert_equal
from numpy.random import RandomState

import nmm


def test_read_hmmer_1(tmp_path):
    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    hmmfile = nmm.hmmer.read_file(tmp_path / filename)
    hmm = nmm.hmmer.create_profile(hmmfile)

    assert_equal(set("ACDEFGHIKLMNPQRSTVWY"), set(hmm.alphabet.symbols))
