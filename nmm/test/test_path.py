from math import log

import pytest
from numpy.testing import assert_equal

from nmm import Alphabet, CPath, MuteState, NormalState


def test_path():
    alphabet = Alphabet(b"AC")
    S = MuteState(b"S", alphabet)
    E = MuteState(b"E", alphabet)
    M = NormalState(b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2)})

    path = CPath()

    step = path.append_cstep(S, 0)
    assert_equal(step.state.name, b"S")
    assert_equal(step.seq_len, 0)

    step = path.append_cstep(E, 0)
    assert_equal(step.state.name, b"E")
    assert_equal(step.seq_len, 0)

    step = path.append_cstep(M, 1)
    assert_equal(step.state.name, b"M1")
    assert_equal(step.seq_len, 1)

    assert_equal(len(list(path.steps())), 3)

    path = CPath()
    with pytest.raises(RuntimeError):
        path.append_cstep(M, 0)

    with pytest.raises(RuntimeError):
        path.append_cstep(M, -1)

    path.append_cstep(M, 1)
    assert_equal(len(list(path.steps())), 1)

    assert_equal(str(path), "<CPath:<M1:1>>")
    assert_equal(str(next(path.steps())), "<M1,1>")
