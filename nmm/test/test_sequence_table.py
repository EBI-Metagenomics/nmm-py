from math import log
import pytest
from numpy.testing import assert_equal, assert_allclose

from nmm import Sequence, Alphabet, SequenceTable, lprob_is_zero


def test_sequence_table():
    alphabet = Alphabet(b"ACGT", b"X")
    seqt = SequenceTable(alphabet)

    with pytest.raises(RuntimeError):
        seqt.normalize()

    seqt.add(Sequence(b"AGTG", alphabet), log(0.2))
    seqt.add(Sequence(b"T", alphabet), log(1.2))

    assert_allclose(seqt.lprob(Sequence(b"AGTG", alphabet)), log(0.2))
    assert_allclose(seqt.lprob(Sequence(b"T", alphabet)), log(1.2))
    assert_equal(lprob_is_zero(seqt.lprob(Sequence(b"", alphabet))), True)

    with pytest.raises(RuntimeError):
        seqt.lprob(Sequence(b"AT", Alphabet(b"AT", b"X")))

    with pytest.raises(RuntimeError):
        seqt.add(Sequence(b"AT", Alphabet(b"AT", b"X")), log(0.2))

    seqt.normalize()

    assert_allclose(seqt.lprob(Sequence(b"AGTG", alphabet)), log(0.2 / 1.4))
    assert_allclose(seqt.lprob(Sequence(b"T", alphabet)), log(1.2 / 1.4))
