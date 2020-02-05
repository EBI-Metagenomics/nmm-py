import pytest
from numpy.testing import assert_equal

from nmm import Alphabet, Base


def test_base():
    alphabet = Alphabet(b"ACGT", b"X")
    base = Base(alphabet)

    assert_equal(base.symbols, b"ACGT")
    assert_equal(str(base), "{ACGT}")
    assert_equal(repr(base), "<Base:{ACGT}>")

    with pytest.raises(RuntimeError):
        alphabet = Alphabet(b"ACGTK", b"X")
        Base(alphabet)
