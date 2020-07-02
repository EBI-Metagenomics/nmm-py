from math import log
from typing import Union

from imm import Fragment, MuteState, NormalState, Path, Sequence, Step

from nmm import BaseAlphabet


def test_fragment():
    alphabet = BaseAlphabet.create(b"ACGT", b"X")
    seq = Sequence.create(b"ACAAAGATX", alphabet)

    S = MuteState.create(b"S", alphabet)
    E = MuteState.create(b"E", alphabet)
    M1 = NormalState.create(
        b"M1", alphabet, [log(0.8), log(0.2), log(0.01), log(0.01)],
    )
    M2 = NormalState.create(b"M2", alphabet, [log(0.4), log(0.6), log(0.1), log(0.6)])

    path = Path.create(
        [Step.create(S, 0), Step.create(M1, 1), Step.create(M2, 1), Step.create(E, 0)]
    )

    fragment = Fragment[BaseAlphabet, Union[MuteState, NormalState]](seq, path)
    i = iter(fragment)

    frag_step = next(i)
    assert bytes(frag_step.sequence) == b""
    assert frag_step.step.seq_len == 0
    assert frag_step.step.state.name == S.name

    frag_step = next(i)
    assert bytes(frag_step.sequence) == b"A"
    assert frag_step.step.seq_len == 1
    assert frag_step.step.state.name == M1.name
