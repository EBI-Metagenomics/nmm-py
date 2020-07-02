from math import inf, log

from imm import Sequence, lprob_is_zero
from imm.testing import assert_allclose

from nmm import (
    BaseAlphabet,
    BaseTable,
    Codon,
    CodonProb,
    CodonState,
    CodonTable,
    FrameState,
)


def test_codon_state():
    base = BaseAlphabet.create(b"ACGU", b"X")
    codonp = CodonProb.create(base)
    codonp.set_lprob(Codon.create(b"AUG", base), log(0.8))
    codonp.set_lprob(Codon.create(b"AUU", base), log(0.1))
    state = CodonState.create(b"M1", codonp)
    assert state.name == b"M1"
    assert_allclose(state.lprob(Sequence.create(b"AUG", base)), log(0.8))
    assert_allclose(state.lprob(Sequence.create(b"AUU", base)), log(0.1))
    assert_allclose(state.lprob(Sequence.create(b"ACU", base)), -inf)


def test_frame_state():
    base = BaseAlphabet.create(b"ACGU", b"X")
    baset = BaseTable.create(base, (log(0.25), log(0.25), log(0.25), log(0.25)))

    codonp = CodonProb.create(base)
    codonp.set_lprob(Codon.create(b"AUG", base), log(0.8))
    codonp.set_lprob(Codon.create(b"AUU", base), log(0.1))

    frame_state = FrameState.create(b"M1", baset, CodonTable.create(codonp), 0.0)

    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"AUA", base)))
    assert_allclose(frame_state.lprob(Sequence.create(b"AUG", base)), log(0.8))
    assert_allclose(frame_state.lprob(Sequence.create(b"AUU", base)), log(0.1))
    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"AU", base)))
    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"A", base)))
    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUA", base)))
    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUAA", base)))

    codonp.normalize()
    frame_state = FrameState.create(b"M1", baset, CodonTable.create(codonp), 0.1)

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
    assert lprob_is_zero(frame_state.lprob(Sequence.create(b"AUUAAA", base)))

    lprob, codon = frame_state.decode(Sequence.create(b"AUA", base))
    assert_allclose(lprob, -7.128586690537968)
    assert codon.symbols == b"AUG"

    lprob, codon = frame_state.decode(Sequence.create(b"AUAG", base))
    assert_allclose(lprob, -4.813151489562624)
    assert codon.symbols == b"AUG"

    lprob, codon = frame_state.decode(Sequence.create(b"A", base))
    assert_allclose(lprob, -6.032286541628237)
    assert codon.symbols == b"AUG"

    lprob, codon = frame_state.decode(Sequence.create(b"UUU", base))
    assert_allclose(lprob, -8.110186062956258)
    assert codon.symbols == b"AUU"
