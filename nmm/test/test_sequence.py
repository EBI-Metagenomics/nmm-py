import pytest
from imm import Interval, Sequence

from nmm import BaseAlphabet


def test_sequence():
    alphabet = BaseAlphabet.create(b"ACGT", b"X")
    seq = Sequence.create(b"ACAAAGATX", alphabet)

    assert len(seq) == 9
    assert bytes(seq) == b"ACAAAGATX"

    assert str(seq) == "ACAAAGATX"
    assert repr(seq) == "<Sequence:ACAAAGATX>"

    subseq = seq[1:7]
    assert str(subseq) == "CAAAGA"
    subseq = subseq[Interval(0, 5)]
    assert str(subseq) == "CAAAG"
    assert subseq.alphabet.symbols == b"ACGT"

    del subseq
    assert seq.alphabet.symbols == b"ACGT"

    Sequence.create(b"ACGXXT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create(b"ACGWT", alphabet)

    with pytest.raises(RuntimeError):
        Sequence.create("ACGTç".encode(), alphabet)
