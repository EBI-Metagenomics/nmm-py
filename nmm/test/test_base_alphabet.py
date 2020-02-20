import pytest
from numpy.testing import assert_equal

from nmm.alphabet import BaseAlphabet


def test_base_alphabet():
    base = BaseAlphabet.create(b"ACGT", b"X")

    assert_equal(base.symbols, b"ACGT")
    assert_equal(str(base), "{ACGT}")
    assert_equal(repr(base), "<BaseAlphabet:{ACGT}>")

    with pytest.raises(RuntimeError):
        BaseAlphabet.create(b"ACGTK", b"X")
