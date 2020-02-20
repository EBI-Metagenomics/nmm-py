import pytest
from numpy.testing import assert_equal

from nmm import Interval
from nmm.alphabet import Alphabet, BaseAlphabet
from nmm.sequence import Sequence


def test_sequence():
    alphabet = Alphabet.create(b"ACGT", b"X")
    seq = Sequence.create(b"ACAAAGATX", alphabet)

    assert_equal(len(seq), 9)
    assert_equal(bytes(seq), b"ACAAAGATX")

    assert_equal(str(seq), "ACAAAGATX")
    assert_equal(repr(seq), "<Sequence:ACAAAGATX>")

    Sequence.create(b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create(b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create("ACGTç".encode(), alphabet)


def test_sequence_base():
    alphabet = BaseAlphabet.create(b"ACGT", b"X")
    seq = Sequence.create(b"ACAAAGATX", alphabet)

    assert_equal(len(seq), 9)
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

    Sequence.create(b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create(b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create("ACGTç".encode(), alphabet)
