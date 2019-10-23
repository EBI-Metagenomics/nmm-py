from math import log

import pytest

from nmm import Alphabet, MuteState, NormalState, Path


def test_path():
    alphabet = Alphabet("AC")
    S = MuteState("S", alphabet)
    E = MuteState("E", alphabet)
    M = NormalState("M1", alphabet, {"A": log(0.8), "C": log(0.2)})

    with pytest.raises(ValueError):
        Path([(S, 0), (E, 0), (M, -1)])
