import pytest
from math import log
from pathlib import Path
from numpy.testing import assert_allclose, assert_equal

from nmm import HMM
from nmm.alphabet import Alphabet
from nmm.io import Model, Output, Input
from nmm.prob import lprob_zero
from nmm.state import MuteState, NormalState
from nmm.sequence import Sequence


@pytest.fixture
def example():
    alphabet = Alphabet.create(b"AC", b"X")
    hmm = HMM.create(alphabet)

    S = MuteState.create(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState.create(b"E", alphabet)
    hmm.add_state(E, lprob_zero())

    M1 = NormalState.create(b"M1", alphabet, [log(0.8), log(0.2)])
    hmm.add_state(M1, lprob_zero())

    M2 = NormalState.create(b"M2", alphabet, [log(0.4), log(0.6)])
    hmm.add_state(M2, lprob_zero())

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()
    hmm.set_transition(E, E, lprob_zero())

    dp = hmm.create_dp(E)

    return {"hmm": hmm, "dp": dp, "alphabet": alphabet}


def test_imm_model(tmpdir, example):
    alphabet = example["alphabet"]
    hmm = example["hmm"]
    dp = example["dp"]

    score = dp.viterbi(Sequence.create(b"AC", alphabet))[0].loglikelihood
    assert_allclose(score, log(0.48))

    filepath = Path(tmpdir / "model.nmm")
    output = Output.create(bytes(filepath))
    output.write(Model.create(hmm, dp))
    output.close()

    input = Input.create(bytes(filepath))
    model = input.read()
    input.close()
    alphabet = model.alphabet
    score = model.dp.viterbi(Sequence.create(b"AC", alphabet))[0].loglikelihood
    assert_allclose(score, log(0.48))


def test_imm_model_iter(tmpdir, example):
    alphabet = example["alphabet"]
    hmm = example["hmm"]
    dp = example["dp"]

    score = dp.viterbi(Sequence.create(b"AC", alphabet))[0].loglikelihood
    assert_allclose(score, log(0.48))

    filepath = Path(tmpdir / "model.nmm")
    with Output.create(bytes(filepath)) as output:
        output.write(Model.create(hmm, dp))
        output.write(Model.create(hmm, dp))
        output.write(Model.create(hmm, dp))

    input = Input.create(bytes(filepath))
    nmodels = 0
    for model in input:
        alphabet = model.alphabet
        seq = Sequence.create(b"AC", alphabet)
        score = model.dp.viterbi(seq)[0].loglikelihood
        assert_allclose(score, log(0.48))
        nmodels += 1
    input.close()
    assert_equal(nmodels, 3)

    with Input.create(bytes(filepath)) as input:
        nmodels = 0
        for model in input:
            alphabet = model.alphabet
            seq = Sequence.create(b"AC", alphabet)
            score = model.dp.viterbi(seq)[0].loglikelihood
            assert_allclose(score, log(0.48))
            nmodels += 1
        assert_equal(nmodels, 3)


# def test_nmm_model(tmpdir):
#     alphabet = Alphabet.create(b"AC", b"X")
#     hmm = HMM.create(alphabet)

#     S = MuteState.create(b"S", alphabet)
#     hmm.add_state(S, log(1.0))

#     E = MuteState.create(b"E", alphabet)
#     hmm.add_state(E, lprob_zero())

#     M1 = NormalState.create(b"M1", alphabet, [log(0.8), log(0.2)])
#     hmm.add_state(M1, lprob_zero())

#     M2 = NormalState.create(b"M2", alphabet, [log(0.4), log(0.6)])
#     hmm.add_state(M2, lprob_zero())

#     hmm.set_transition(S, M1, log(1.0))
#     hmm.set_transition(M1, M2, log(1.0))
#     hmm.set_transition(M2, E, log(1.0))
#     hmm.set_transition(E, E, log(1.0))
#     hmm.normalize()
#     hmm.set_transition(E, E, lprob_zero())

#     dp = hmm.create_dp(E)

#     score = dp.viterbi(Sequence.create(b"AC", alphabet))[0].loglikelihood
#     assert_allclose(score, log(0.48))

#     filepath = Path(tmpdir / "model.nmm")
#     output = Output.create(bytes(filepath))
#     output.write(Model.create(hmm, dp))
#     output.close()

#     input = Input.create(bytes(filepath))
# model = input.read()
# alphabet = model.alphabet
# score = model.dp.viterbi(Sequence.create(b"AC", alphabet))[0].loglikelihood
# assert_allclose(score, log(0.48))
