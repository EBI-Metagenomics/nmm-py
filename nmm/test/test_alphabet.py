import pytest
from numpy.testing import assert_equal

from nmm import Alphabet


def test_alphabet():
    abc = Alphabet(b"ACGT")
    assert_equal(abc.length, 4)

    assert_equal(abc.has_symbol(b"A"), True)
    assert_equal(abc.has_symbol(b"C"), True)
    assert_equal(abc.has_symbol(b"G"), True)
    assert_equal(abc.has_symbol(b"T"), True)

    assert_equal(abc.symbol_idx(b"A"), 0)
    assert_equal(abc.symbol_idx(b"C"), 1)
    assert_equal(abc.symbol_idx(b"G"), 2)
    assert_equal(abc.symbol_idx(b"T"), 3)

    assert_equal(abc.symbol_id(0), b"A")
    assert_equal(abc.symbol_id(1), b"C")
    assert_equal(abc.symbol_id(2), b"G")
    assert_equal(abc.symbol_id(3), b"T")

    assert_equal(abc.symbols, b"ACGT")

    assert_equal(str(abc), "{ACGT}")
    assert_equal(repr(abc), "{Alphabet:ACGT}")

    with pytest.raises(TypeError):
        Alphabet("ACGTç")

    with pytest.raises(RuntimeError):
        Alphabet("ACGTç".encode())
