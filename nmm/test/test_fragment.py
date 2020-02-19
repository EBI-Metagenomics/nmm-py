from typing import Union
from numpy.testing import assert_equal
from math import log

from nmm.sequence import Sequence
from nmm.alphabet import BaseAlphabet
from nmm.fragment import Fragment
from nmm.path import Path, Step
from nmm.state import MuteState, NormalState


def test_fragment():
    alphabet = BaseAlphabet(b"ACGT", b"X")
    seq = Sequence[BaseAlphabet](b"ACAAAGATX", alphabet)

    S = MuteState(b"S", alphabet)
    E = MuteState(b"E", alphabet)
    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2), log(0.01), log(0.01)],)
    M2 = NormalState(b"M2", alphabet, [log(0.4), log(0.6), log(0.1), log(0.6)])

    path = Path([Step(S, 0), Step(M1, 1), Step(M2, 1), Step(E, 0)])

    fragment = Fragment[BaseAlphabet, Union[MuteState, NormalState]](seq, path)
    i = iter(fragment)

    frag_step = next(i)
    assert_equal(frag_step.sequence.symbols, b"")
    assert_equal(frag_step.step.seq_len, 0)
    assert_equal(frag_step.step.state.name, S.name)

    frag_step = next(i)
    assert_equal(frag_step.sequence.symbols, b"A")
    assert_equal(frag_step.step.seq_len, 1)
    assert_equal(frag_step.step.state.name, M1.name)
