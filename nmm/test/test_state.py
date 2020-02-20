from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm.alphabet import Alphabet, BaseAlphabet
from nmm.codon import Codon
from nmm.prob import BaseTable, CodonProb, CodonTable, SequenceTable, lprob_is_zero
from nmm.sequence import Sequence
from nmm.state import FrameState, MuteState, NormalState, TableState


def test_normal_state():
    alphabet = Alphabet.create(b"ACGT", b"X")

    state = NormalState(b"M0", alphabet, [log(0.1), log(0.2), log(0.3), log(0.3)],)
    assert_equal(state.name, b"M0")
    assert_equal(state.lprob(Sequence.create(b"A", alphabet)), log(0.1))
    assert_equal(state.lprob(Sequence.create(b"C", alphabet)), log(0.2))
    assert_equal(state.lprob(Sequence.create(b"G", alphabet)), log(0.3))
    assert_equal(state.lprob(Sequence.create(b"T", alphabet)), log(0.3))
    assert_equal(state.min_seq, 1)
    assert_equal(state.max_seq, 1)

    with pytest.raises(RuntimeError):
        state.lprob(Sequence.create(b"T", Alphabet.create(b"ACGT", b"X")))

    assert_equal(lprob_is_zero(state.lprob(Sequence.create(b"AC", alphabet))), True)

    assert_equal(str(state), "M0")
    assert_equal(repr(state), "<NormalState:M0>")


def test_mute_state():
    alphabet = Alphabet.create(b"ACGU", b"X")
    state = MuteState(b"S", alphabet)

    assert_equal(state.name, b"S")
    assert_equal(state.lprob(Sequence.create(b"", alphabet)), log(1.0))
    assert_equal(lprob_is_zero(state.lprob(Sequence.create(b"AC", alphabet))), True)
    assert_equal(state.min_seq, 0)
    assert_equal(state.max_seq, 0)
    assert_equal(str(state), "S")
    assert_equal(repr(state), "<MuteState:S>")


def test_table_state():
    alphabet = Alphabet.create(b"ACGU", b"X")
    seqt = SequenceTable(alphabet)
    seqt.add(Sequence.create(b"AUG", alphabet), log(0.8))
    seqt.add(Sequence.create(b"AUU", alphabet), log(0.4))

    state = TableState(b"M2", seqt)
    assert_equal(state.name, b"M2")
    assert_allclose(state.lprob(Sequence.create(b"AUG", alphabet)), log(0.8))
    assert_allclose(state.lprob(Sequence.create(b"AUU", alphabet)), log(0.4))
    assert_equal(str(state), "M2")
    assert_equal(repr(state), "<TableState:M2>")


def test_frame_state():
    base = BaseAlphabet.create(b"ACGU", b"X")
    baset = BaseTable(base, (log(0.25), log(0.25), log(0.25), log(0.25)))

    codonp = CodonProb(base)
    codonp.set_lprob(Codon.create(b"AUG", base), log(0.8))
    codonp.set_lprob(Codon.create(b"AUU", base), log(0.1))

    frame_state = FrameState(b"M1", baset, CodonTable(codonp), 0.0)

    assert_equal(lprob_is_zero(frame_state.lprob(Sequence.create(b"AUA", base))), True)
    assert_allclose(frame_state.lprob(Sequence.create(b"AUG", base)), log(0.8))
    assert_allclose(frame_state.lprob(Sequence.create(b"AUU", base)), log(0.1))
    assert_equal(lprob_is_zero(frame_state.lprob(Sequence.create(b"AU", base))), True)
    assert_equal(lprob_is_zero(frame_state.lprob(Sequence.create(b"A", base))), True)
    assert_equal(lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUA", base))), True)
    assert_equal(
        lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUAA", base))), True
    )

    codonp.normalize()
    frame_state = FrameState(b"M1", baset, CodonTable(codonp), 0.1)

    assert_allclose(
        frame_state.lprob(Sequence.create(b"AUA", base)), -6.905597115665666
    )
    assert_allclose(
        frame_state.lprob(Sequence.create(b"AUG", base)), -0.5347732882047062
    )
    assert_allclose(
        frame_state.lprob(Sequence.create(b"AUU", base)), -2.5902373304999466
    )
    assert_allclose(
        frame_state.lprob(Sequence.create(b"AU", base)), -2.9158434238698336
    )
    assert_allclose(frame_state.lprob(Sequence.create(b"A", base)), -5.914503505971854)
    assert_allclose(
        frame_state.lprob(Sequence.create(b"AUUA", base)), -6.881032208841384
    )
    assert_allclose(
        frame_state.lprob(Sequence.create(b"AUUAA", base)), -12.08828960987379
    )
    assert_equal(
        lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUAAA", base))), True
    )

    codon = Codon.create(b"XXX", base)
    lprob = frame_state.decode(Sequence.create(b"AUA", base), codon)
    assert_allclose(lprob, -7.128586690537968)
    assert_equal(codon.symbols, b"AUG")

    lprob = frame_state.decode(Sequence.create(b"AUAG", base), codon)
    assert_allclose(lprob, -4.813151489562624)
    assert_equal(codon.symbols, b"AUG")

    lprob = frame_state.decode(Sequence.create(b"A", base), codon)
    assert_allclose(lprob, -6.032286541628237)
    assert_equal(codon.symbols, b"AUG")

    lprob = frame_state.decode(Sequence.create(b"UUU", base), codon)
    assert_allclose(lprob, -8.110186062956258)
    assert_equal(codon.symbols, b"AUU")
