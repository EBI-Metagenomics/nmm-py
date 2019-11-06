from math import log

import pytest

from nmm import Alphabet, MuteState, NormalState, Path, Step


def test_path():
    alphabet = Alphabet(b"AC")
    S = MuteState(b"S", alphabet)
    E = MuteState(b"E", alphabet)
    M = NormalState(b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2)})

    with pytest.raises(RuntimeError):
        Path([Step(S, 0), Step(E, 0), Step(M, -1)])
