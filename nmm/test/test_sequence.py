import pytest
from numpy.testing import assert_equal

from nmm.sequence import Sequence
from nmm.alphabet import Alphabet


def test_sequence():
    alphabet = Alphabet(b"ACGT", b"X")
    seq = Sequence(b"ACAAAGATX", alphabet)

    assert_equal(seq.length, 9)
    assert_equal(seq.symbols, b"ACAAAGATX")

    assert_equal(str(seq), "[ACAAAGATX]")
    assert_equal(repr(seq), "<Sequence:[ACAAAGATX]>")

    Sequence(b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence(b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence("ACGTç".encode(), alphabet)
