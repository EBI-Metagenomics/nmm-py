import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import (
    LOG,
    Alphabet,
    FrameState,
    MuteState,
    NormalState,
    TableState,
    Base,
    Codon,
)


def test_normal_state():
    alphabet = Alphabet("ACGT")
    with pytest.raises(ValueError):
        NormalState("M0", alphabet, [LOG(v) for v in [0.1, 0.2, 0.3]])

    state = NormalState("M0", alphabet, [LOG(v) for v in [0.1, 0.2, 0.3, 0.3]])
    assert_equal(state.name, "M0")
    assert_equal(state.lprob("A"), LOG(0.1))
    assert_equal(state.lprob("C"), LOG(0.2))
    assert_equal(state.lprob("G"), LOG(0.3))
    assert_equal(state.lprob("T"), LOG(0.3))
    assert_equal(state.min_seq, 1)
    assert_equal(state.max_seq, 1)

    state.normalize()

    assert_equal(state.lprob("A"), LOG(0.1) - LOG(0.9))
    assert_equal(state.lprob("C"), LOG(0.2) - LOG(0.9))
    assert_equal(state.lprob("G"), LOG(0.3) - LOG(0.9))
    assert_equal(state.lprob("T"), LOG(0.3) - LOG(0.9))

    assert_equal(state.alphabet.length, alphabet.length)
    assert_equal(state.alphabet.cdata, alphabet.cdata)

    table = state.emission()
    assert_equal(table[0][0], "G")
    assert_allclose(table[0][1], LOG(0.3) - LOG(0.9))
    assert_equal(table[1][0], "T")
    assert_allclose(table[1][1], LOG(0.3) - LOG(0.9))
    assert_equal(table[2][0], "C")
    assert_allclose(table[2][1], LOG(0.2) - LOG(0.9))
    assert_equal(table[3][0], "A")
    assert_allclose(table[3][1], LOG(0.1) - LOG(0.9))

    assert_equal(str(state), "<M0>")
    assert_equal(repr(state), "<NormalState:M0>")


def test_mute_state():
    alphabet = Alphabet("ACGU")
    state = MuteState("S", alphabet)
    assert_equal(state.name, "S")
    assert_equal(state.alphabet.symbols, "ACGU")
    assert_equal(state.lprob(""), LOG(1.0))
    assert_equal(state.lprob("A"), LOG(0.0))
    assert_equal(str(state), "<S>")
    assert_equal(repr(state), "<MuteState:S>")

    table = state.emission()
    assert_equal(len(table), 1)
    assert_equal(table[0][0], "")
    assert_equal(table[0][1], LOG(1.0))


def test_table_state():
    alphabet = Alphabet("ACGU")
    state = TableState("M2", alphabet, {"AUG": LOG(0.8), "AUU": LOG(0.4)})
    assert_equal(state.name, "M2")
    assert_equal(set(state.alphabet.symbols), set("ACGU"))
    assert_allclose(state.lprob("AUG"), LOG(0.8))
    assert_allclose(state.lprob("AUU"), LOG(0.4))
    assert_equal(state.lprob("AGU"), LOG(0.0))
    assert_equal(str(state), "<M2>")
    assert_equal(repr(state), "<TableState:M2>")
    state.normalize()
    assert_allclose(state.lprob("AUG"), LOG(0.8) - LOG(1.2))
    assert_allclose(state.lprob("AUU"), LOG(0.4) - LOG(1.2))
    assert_equal(state.lprob("AGU"), LOG(0.0))


def test_frame_state():
    alphabet = Alphabet("ACGU")
    base = Base(
        alphabet, {"A": LOG(0.25), "C": LOG(0.25), "G": LOG(0.25), "U": LOG(0.25)}
    )
    codon = Codon(alphabet, {"AUG": LOG(0.8), "AUU": LOG(0.1)})

    frame_state = FrameState("M1", base, codon, epsilon=0.0)
    assert_allclose(frame_state.lprob("AUA"), LOG(0.0))
    assert_allclose(frame_state.lprob("AUG"), LOG(0.8))
    assert_allclose(frame_state.lprob("AUU"), LOG(0.1))
    assert_allclose(frame_state.lprob("AU"), LOG(0.0))
    assert_allclose(frame_state.lprob("A"), LOG(0.0))
    assert_allclose(frame_state.lprob("AUUA"), LOG(0.0))
    assert_allclose(frame_state.lprob("AUUAA"), LOG(0.0))

    codon.normalize()
    frame_state = FrameState("M1", base, codon, 0.1)
    assert_allclose(frame_state.lprob("AUA"), -6.905597115665666)
    assert_allclose(frame_state.lprob("AUG"), -0.5347732882047062)
    assert_allclose(frame_state.lprob("AUU"), -2.5902373304999466)
    assert_allclose(frame_state.lprob("AU"), -2.9158434238698336)
    assert_allclose(frame_state.lprob("A"), -5.914503505971854)
    assert_allclose(frame_state.lprob("AUUA"), -6.881032208841384)
    assert_allclose(frame_state.lprob("AUUAA"), -12.08828960987379)
    assert_allclose(frame_state.lprob("AUUAAA"), LOG(0.0))

    alphabet = Alphabet("ACGT")
    base = Base(alphabet, {"A": LOG(0.1), "C": LOG(0.2), "G": LOG(0.3), "T": LOG(0.4)})
    codon = Codon(alphabet, {"ATG": LOG(0.8), "ATT": LOG(0.1), "GTC": LOG(0.4)})
    codon.normalize()
    frame_state = FrameState("M2", base, codon, 0.1)
    assert_allclose(frame_state.lprob("A"), -6.282228286097171)
    assert_allclose(frame_state.lprob("C"), -7.0931585023135)
    assert_allclose(frame_state.lprob("G"), -5.99454621364539)
    assert_allclose(frame_state.lprob("T"), -5.840395533818132)
    assert_allclose(frame_state.lprob("AT"), -3.283414346005771)
    assert_allclose(frame_state.lprob("CG"), -9.395743595307545)
    assert_allclose(frame_state.lprob("ATA"), -8.18911998648269)
    assert_allclose(frame_state.lprob("ATG"), -0.9021560981322401)
    assert_allclose(frame_state.lprob("ATT"), -2.9428648000333952)
    assert_allclose(frame_state.lprob("ATC"), -7.314811395663229)
    assert_allclose(frame_state.lprob("GTC"), -1.5951613351178675)
    assert_allclose(frame_state.lprob("ATTA"), -8.157369364264277)
    assert_allclose(frame_state.lprob("GTTC"), -4.711642430498609)
    assert_allclose(frame_state.lprob("ACTG"), -5.404361876760574)
    assert_allclose(frame_state.lprob("ATTAA"), -14.288595853747417)
    assert_allclose(frame_state.lprob("GTCAA"), -12.902301492627526)
    assert_allclose(frame_state.lprob("ATTAAA"), LOG(0.0))
