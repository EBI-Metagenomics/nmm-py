from math import log

import pytest

from nmm import Alphabet, MuteState, NormalState, Path, Step


def test_path():
    alphabet = Alphabet("AC")
    S = MuteState("S", alphabet)
    E = MuteState("E", alphabet)
    M = NormalState("M1", alphabet, {"A": log(0.8), "C": log(0.2)})

    with pytest.raises(ValueError):
        Path([Step(S, 0), Step(E, 0), Step(M, -1)])
