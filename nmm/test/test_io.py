from math import log
from pathlib import Path

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import HMM
from nmm.alphabet import Alphabet, BaseAlphabet
from nmm.codon import Codon
from nmm.io import Input, Model, Output
from nmm.prob import BaseTable, CodonProb, CodonTable, lprob_zero
from nmm.sequence import Sequence
from nmm.state import FrameState, MuteState, NormalState


@pytest.fixture
def imm_example():
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


@pytest.fixture
def nmm_example():
    abc = BaseAlphabet.create(b"ACGU", b"X")
    baset = BaseTable.create(abc, (log(0.25), log(0.25), log(0.25), log(0.25)))

    codonp = CodonProb.create(abc)
    codonp.set_lprob(Codon.create(b"AUG", abc), log(0.8))
    codonp.set_lprob(Codon.create(b"AUU", abc), log(0.1))

    B = MuteState.create(b"B", abc)
    M1 = FrameState.create(b"M1", baset, CodonTable.create(codonp), 0.02)
    M2 = FrameState.create(b"M2", baset, CodonTable.create(codonp), 0.01)
    E = MuteState.create(b"E", abc)

    hmm = HMM.create(abc)
    hmm.add_state(B, log(0.5))
    hmm.add_state(M1)
    hmm.add_state(M2)
    hmm.add_state(E)

    hmm.set_transition(B, M1, log(0.8))
    hmm.set_transition(B, M2, log(0.2))
    hmm.set_transition(M1, M2, log(0.1))
    hmm.set_transition(M1, E, log(0.4))
    hmm.set_transition(M2, E, log(0.3))

    dp = hmm.create_dp(E)

    return {"hmm": hmm, "dp": dp, "alphabet": abc}


def test_imm_model(tmpdir, imm_example):
    alphabet = imm_example["alphabet"]
    hmm = imm_example["hmm"]
    dp = imm_example["dp"]

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


def test_imm_model_iter(tmpdir, imm_example):
    alphabet = imm_example["alphabet"]
    hmm = imm_example["hmm"]
    dp = imm_example["dp"]

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


def test_nmm_model(tmpdir, nmm_example):
    alphabet = nmm_example["alphabet"]
    hmm = nmm_example["hmm"]
    dp = nmm_example["dp"]

    seq = Sequence.create(b"AUGAUU", alphabet)
    results = dp.viterbi(seq)
    assert_equal(len(results), 1)
    assert_allclose(results[0].loglikelihood, -7.069201008427531)

    filepath = Path(tmpdir / "model.nmm")
    with Output.create(bytes(filepath)) as output:
        output.write(Model.create(hmm, dp))
        output.write(Model.create(hmm, dp))
        output.write(Model.create(hmm, dp))

    with Input.create(bytes(filepath)) as input:
        nmodels = 0
        for model in input:
            alphabet = model.alphabet
            seq = Sequence.create(b"AUGAUU", alphabet)
            score = model.dp.viterbi(seq)[0].loglikelihood
            assert_allclose(score, -7.069201008427531)
            nmodels += 1
        assert_equal(nmodels, 3)
