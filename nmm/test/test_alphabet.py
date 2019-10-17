from numpy.testing import assert_equal
from nmm import Alphabet


def test_alphabet():
    abc = Alphabet("ACGT")
    assert_equal(abc.length, 4)

    assert_equal(abc.has_symbol("A"), True)
    assert_equal(abc.has_symbol("C"), True)
    assert_equal(abc.has_symbol("G"), True)
    assert_equal(abc.has_symbol("T"), True)

    assert_equal(abc.symbol_idx("A"), 0)
    assert_equal(abc.symbol_idx("C"), 1)
    assert_equal(abc.symbol_idx("G"), 2)
    assert_equal(abc.symbol_idx("T"), 3)

    assert_equal(abc.symbol_id(0), "A")
    assert_equal(abc.symbol_id(1), "C")
    assert_equal(abc.symbol_id(2), "G")
    assert_equal(abc.symbol_id(3), "T")

