import pytest
from numpy.testing import assert_equal

from nmm import Alphabet, Base


def test_base():
    base = Base(Alphabet(b"ACGT", b"X"))

    assert_equal(base.symbols, b"ACGT")
    assert_equal(str(base), "{ACGT}")
    assert_equal(repr(base), "<Base:{ACGT}>")

    with pytest.raises(RuntimeError):
        Base(Alphabet(b"ACGTK", b"X"))
