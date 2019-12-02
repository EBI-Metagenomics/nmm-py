from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import (
    LOG0,
    Alphabet,
    Base,
    BaseTable,
    Codon,
    CodonTable,
    FrameState,
    MuteState,
    NormalState,
    TableState,
)


def test_normal_state():
    alphabet = Alphabet(b"ACGT")
    with pytest.raises(ValueError):
        NormalState(b"M0", alphabet, {b"A": log(0.1), b"C": log(0.2), b"X": log(0.3)})

    state = NormalState(
        b"M0",
        alphabet,
        {b"A": log(0.1), b"C": log(0.2), b"G": log(0.3), b"T": log(0.3)},
    )
    assert_equal(state.name, b"M0")
    assert_equal(state.lprob(b"A"), log(0.1))
    assert_equal(state.lprob(b"C"), log(0.2))
    assert_equal(state.lprob(b"G"), log(0.3))
    assert_equal(state.lprob(b"T"), log(0.3))
    assert_equal(state.min_seq, 1)
    assert_equal(state.max_seq, 1)

    state.normalize()

    assert_equal(state.lprob(b"A"), log(0.1) - log(0.9))
    assert_equal(state.lprob(b"C"), log(0.2) - log(0.9))
    assert_equal(state.lprob(b"G"), log(0.3) - log(0.9))
    assert_equal(state.lprob(b"T"), log(0.3) - log(0.9))

    assert_equal(state.alphabet.length, alphabet.length)
    assert_equal(state.alphabet.imm_abc, alphabet.imm_abc)

    assert_equal(str(state), "<M0>")
    assert_equal(repr(state), "<NormalState:M0>")

    state = NormalState(b"M", alphabet, {b"A": LOG0})
    with pytest.raises(RuntimeError):
        state.normalize()


def test_mute_state():
    alphabet = Alphabet(b"ACGU")
    state = MuteState(b"S", alphabet)
    assert_equal(state.name, b"S")
    assert_equal(state.alphabet.symbols, b"ACGU")
    assert_equal(state.lprob(b""), log(1.0))
    assert_equal(state.lprob(b"A"), LOG0)
    assert_equal(str(state), "<S>")
    assert_equal(repr(state), "<MuteState:S>")


def test_table_state():
    alphabet = Alphabet(b"ACGU")
    state = TableState(b"M2", alphabet, {b"AUG": log(0.8), b"AUU": log(0.4)})
    assert_equal(state.name, b"M2")
    assert_equal(set(state.alphabet.symbols), set(b"ACGU"))
    assert_allclose(state.lprob(b"AUG"), log(0.8))
    assert_allclose(state.lprob(b"AUU"), log(0.4))
    assert_equal(state.lprob(b"AGU"), LOG0)
    assert_equal(str(state), "<M2>")
    assert_equal(repr(state), "<TableState:M2>")
    state.normalize()
    assert_allclose(state.lprob(b"AUG"), log(0.8) - log(1.2))
    assert_allclose(state.lprob(b"AUU"), log(0.4) - log(1.2))
    assert_equal(state.lprob(b"AGU"), LOG0)

    state = TableState(b"M", alphabet, {b"A": LOG0})
    with pytest.raises(RuntimeError):
        state.normalize()


def test_frame_state():
    alphabet = Alphabet(b"ACGU")
    base = BaseTable.create(
        alphabet,
        {
            Base("A"): log(0.25),
            Base("C"): log(0.25),
            Base("G"): log(0.25),
            Base("U"): log(0.25),
        },
    )
    codon = CodonTable.create(
        alphabet, {Codon("AUG"): log(0.8), Codon("AUU"): log(0.1)}
    )

    frame_state = FrameState(b"M1", base, codon, epsilon=0.0)
    assert_allclose(frame_state.lprob(b"AUA"), LOG0)
    assert_allclose(frame_state.lprob(b"AUG"), log(0.8))
    assert_allclose(frame_state.lprob(b"AUU"), log(0.1))
    assert_allclose(frame_state.lprob(b"AU"), LOG0)
    assert_allclose(frame_state.lprob(b"A"), LOG0)
    assert_allclose(frame_state.lprob(b"AUUA"), LOG0)
    assert_allclose(frame_state.lprob(b"AUUAA"), LOG0)

    codon.normalize()
    frame_state = FrameState(b"M1", base, codon, 0.1)
    assert_allclose(frame_state.lprob(b"AUA"), -6.905597115665666)
    assert_allclose(frame_state.lprob(b"AUG"), -0.5347732882047062)
    assert_allclose(frame_state.lprob(b"AUU"), -2.5902373304999466)
    assert_allclose(frame_state.lprob(b"AU"), -2.9158434238698336)
    assert_allclose(frame_state.lprob(b"A"), -5.914503505971854)
    assert_allclose(frame_state.lprob(b"AUUA"), -6.881032208841384)
    assert_allclose(frame_state.lprob(b"AUUAA"), -12.08828960987379)
    assert_allclose(frame_state.lprob(b"AUUAAA"), LOG0)

    alphabet = Alphabet(b"ACGT")
    base = BaseTable.create(
        alphabet,
        {
            Base("A"): log(0.1),
            Base("C"): log(0.2),
            Base("G"): log(0.3),
            Base("T"): log(0.4),
        },
    )
    codon = CodonTable.create(
        alphabet,
        {Codon("ATG"): log(0.8), Codon("ATT"): log(0.1), Codon("GTC"): log(0.4)},
    )
    codon.normalize()
    frame_state = FrameState(b"M2", base, codon, 0.1)
    assert_allclose(frame_state.lprob(b"A"), -6.282228286097171)
    assert_allclose(frame_state.lprob(b"C"), -7.0931585023135)
    assert_allclose(frame_state.lprob(b"G"), -5.99454621364539)
    assert_allclose(frame_state.lprob(b"T"), -5.840395533818132)
    assert_allclose(frame_state.lprob(b"AT"), -3.283414346005771)
    assert_allclose(frame_state.lprob(b"CG"), -9.395743595307545)
    assert_allclose(frame_state.lprob(b"ATA"), -8.18911998648269)
    assert_allclose(frame_state.lprob(b"ATG"), -0.9021560981322401)
    assert_allclose(frame_state.lprob(b"ATT"), -2.9428648000333952)
    assert_allclose(frame_state.lprob(b"ATC"), -7.314811395663229)
    assert_allclose(frame_state.lprob(b"GTC"), -1.5951613351178675)
    assert_allclose(frame_state.lprob(b"ATTA"), -8.157369364264277)
    assert_allclose(frame_state.lprob(b"GTTC"), -4.711642430498609)
    assert_allclose(frame_state.lprob(b"ACTG"), -5.404361876760574)
    assert_allclose(frame_state.lprob(b"ATTAA"), -14.288595853747417)
    assert_allclose(frame_state.lprob(b"GTCAA"), -12.902301492627526)
    assert_allclose(frame_state.lprob(b"ATTAAA"), LOG0)

    assert_equal(str(frame_state), "<M2>")
    assert_equal(repr(frame_state), "<FrameState:M2>")

    codon, lprob = frame_state.decode(b"ATG")
    assert_allclose(lprob, -0.9025667061364364)
    assert_equal(bytes(codon), b"ATG")

    codon, lprob = frame_state.decode(b"ATGGG")
    assert_allclose(lprob, -8.913317446063251)
    assert_equal(bytes(codon), b"ATG")

    codon, lprob = frame_state.decode(b"T")
    assert_allclose(lprob, -6.400011321753555)
    assert_equal(bytes(codon), b"ATG")

    codon, lprob = frame_state.decode(b"TT")
    assert_allclose(lprob, -5.579253016600963)
    assert_equal(bytes(codon), b"ATT")

    codon, lprob = frame_state.decode(b"GC")
    assert_allclose(lprob, -4.199705077879926)
    assert_equal(bytes(codon), b"GTC")
