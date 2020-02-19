import pytest
from numpy.testing import assert_equal

from nmm import Interval
from nmm.sequence import Sequence
from nmm.alphabet import Alphabet, BaseAlphabet


def test_sequence():
    alphabet = Alphabet(b"ACGT", b"X")
    seq = Sequence(b"ACAAAGATX", alphabet)

    assert_equal(seq.length, 9)
    assert_equal(bytes(seq), b"ACAAAGATX")

    assert_equal(str(seq), "ACAAAGATX")
    assert_equal(repr(seq), "<Sequence:ACAAAGATX>")

    Sequence(b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence(b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence("ACGTç".encode(), alphabet)


def test_sequence_base():
    alphabet = BaseAlphabet(b"ACGT", b"X")
    seq = Sequence[BaseAlphabet](b"ACAAAGATX", alphabet)

    assert_equal(seq.length, 9)
    assert_equal(bytes(seq), b"ACAAAGATX")

    assert_equal(str(seq), "ACAAAGATX")
    assert_equal(repr(seq), "<Sequence:ACAAAGATX>")

    subseq = seq[1:7]
    assert_equal(str(subseq), "CAAAGA")
    subseq = subseq[Interval(0, 5)]
    assert_equal(str(subseq), "CAAAG")
    assert_equal(subseq.alphabet.symbols, b"ACGT")

    del subseq
    assert_equal(seq.alphabet.symbols, b"ACGT")

    Sequence[BaseAlphabet](b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence[BaseAlphabet](b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence[BaseAlphabet]("ACGTç".encode(), alphabet)
