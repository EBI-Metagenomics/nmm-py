from math import log

import pytest

from nmm import Alphabet, MuteState, NormalState, Path


def test_path():
    alphabet = Alphabet(b"AC")
    S = MuteState(b"S", alphabet)
    E = MuteState(b"E", alphabet)
    M = NormalState(b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2)})

    path = Path()

    step = path.append(S, 0)
    assert step.state.name == b"S"
    assert step.seq_len == 0

    step = path.append(E, 0)
    assert step.state.name == b"E"
    assert step.seq_len == 0

    step = path.append(M, 1)
    assert step.state.name == b"M1"
    assert step.seq_len == 1

    assert len(list(path.steps())) == 3

    path = Path()
    with pytest.raises(RuntimeError):
        path.append(M, 0)

    with pytest.raises(RuntimeError):
        path.append(M, -1)

    path.append(M, 1)
    assert len(list(path.steps())) == 1
