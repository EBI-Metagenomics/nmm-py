import pytest
from numpy.testing import assert_equal

from nmm.alphabet import Alphabet, BaseAlphabet


def test_base():
    base = BaseAlphabet(Alphabet(b"ACGT", b"X"))

    assert_equal(base.symbols, b"ACGT")
    assert_equal(str(base), "{ACGT}")
    assert_equal(repr(base), "<BaseAlphabet:{ACGT}>")

    with pytest.raises(RuntimeError):
        BaseAlphabet(Alphabet(b"ACGTK", b"X"))
