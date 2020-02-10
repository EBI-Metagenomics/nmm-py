from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import (
    HMM,
    Alphabet,
    MuteState,
    NormalState,
    TableState,
    SequenceTable,
    Sequence,
    LPROB_ZERO,
    LPROB_INVALID,
    Path,
)


def test_hmm_states():
    alphabet = Alphabet(b"ACGU", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S)

    seqt = SequenceTable(alphabet)
    seqt.add(Sequence(b"AGU", alphabet), log(0.8))
    seqt.add(Sequence(b"AGG", alphabet), log(0.2))

    M = TableState(b"M", seqt)
    hmm.add_state(M)

    with pytest.raises(ValueError):
        hmm.add_state(S)

    with pytest.raises(ValueError):
        hmm.add_state(M)

    assert_equal(len(hmm.states()), 2)


def test_hmm_trans_prob():
    alphabet = Alphabet(b"ACGU", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    with pytest.raises(RuntimeError):
        hmm.set_start_lprob(S, log(0.4))
    hmm.add_state(S)

    E = MuteState(b"E", alphabet)
    with pytest.raises(RuntimeError):
        hmm.transition(S, E)

    with pytest.raises(ValueError):
        hmm.set_transition(S, E, LPROB_ZERO)

    with pytest.raises(ValueError):
        hmm.set_transition(E, S, LPROB_ZERO)

    with pytest.raises(ValueError):
        hmm.del_state(E)

    hmm.add_state(E)

    with pytest.raises(RuntimeError):
        hmm.set_transition(E, S, LPROB_INVALID)

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_transition(S, E, log(0.5))

    assert_allclose(hmm.transition(S, S), LPROB_ZERO)
    assert_allclose(hmm.transition(S, E), log(0.5))
    assert_allclose(hmm.transition(E, S), LPROB_ZERO)
    assert_allclose(hmm.transition(E, E), LPROB_ZERO)

    with pytest.raises(ValueError):
        hmm.normalize()

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_start_lprob(S, log(0.4))
    hmm.set_transition(E, E, log(0.1))

    hmm.normalize()

    assert_allclose(hmm.transition(S, E), log(1.0))
    assert_allclose(hmm.transition(E, S), LPROB_ZERO)
    assert_allclose(hmm.transition(S, S), LPROB_ZERO)
    assert_allclose(hmm.transition(E, E), log(1.0))


def test_hmm_likelihood():
    alphabet = Alphabet(b"ACGU", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LPROB_ZERO)

    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2), LPROB_ZERO, LPROB_ZERO],)
    hmm.add_state(M1, LPROB_ZERO)

    M2 = NormalState(
        b"M2", alphabet, [log(0.4 / 1.6), log(0.6 / 1.6), LPROB_ZERO, log(0.6 / 1.6)]
    )
    hmm.add_state(M2, LPROB_ZERO)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()

    p = hmm.likelihood(
        Sequence(b"AC", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(
        Sequence(b"AA", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.2))

    p = hmm.likelihood(
        Sequence(b"AG", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"AU", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(
        Sequence(b"CC", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(
        Sequence(b"CA", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.05))

    p = hmm.likelihood(
        Sequence(b"CG", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"CG", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"CU", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(
        Sequence(b"GC", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"GA", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"GG", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"GU", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"UC", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"UA", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"UG", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    p = hmm.likelihood(
        Sequence(b"UU", alphabet), Path([(S, 0), (M1, 1), (M2, 1), (E, 0)])
    )
    assert_allclose(p, LPROB_ZERO)

    M3 = NormalState(b"M2", alphabet, [log(0.4), log(0.6), LPROB_ZERO, log(0.6)],)

    with pytest.raises(ValueError):
        hmm.likelihood(
            Sequence(b"UU", alphabet), Path([(S, 0), (M1, 1), (M3, 1), (E, 0)])
        )


def test_hmm_viterbi_1():
    alphabet = Alphabet(b"ACGU", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LPROB_ZERO)

    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2), LPROB_ZERO, LPROB_ZERO],)
    hmm.add_state(M1, LPROB_ZERO)

    M2 = NormalState(
        b"M2", alphabet, [log(0.4 / 1.6), log(0.6 / 1.6), LPROB_ZERO, log(0.6 / 1.6)],
    )
    hmm.add_state(M2, LPROB_ZERO)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()

    hmm.set_transition(E, E, LPROB_ZERO)
    assert_allclose(hmm.transition(E, E), LPROB_ZERO)
    assert_allclose(hmm.transition(S, S), LPROB_ZERO)
    assert_allclose(hmm.transition(S, E), LPROB_ZERO)
    assert_allclose(hmm.transition(E, S), LPROB_ZERO)

    results = hmm.viterbi(Sequence(b"AC", alphabet), E)
    assert_equal(len(results), 1)
    assert_allclose(results[0].loglikelihood, log(0.3))


def test_hmm_viterbi_2():
    alphabet = Alphabet(b"AC", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LPROB_ZERO)

    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2)])
    hmm.add_state(M1, LPROB_ZERO)

    M2 = NormalState(b"M2", alphabet, [log(0.4), log(0.6)])
    hmm.add_state(M2, LPROB_ZERO)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()
    hmm.set_transition(E, E, LPROB_ZERO)

    score = hmm.viterbi(Sequence(b"AC", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.48))

    score = hmm.viterbi(Sequence(b"AA", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.32))

    score = hmm.viterbi(Sequence(b"CA", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.08))

    score = hmm.viterbi(Sequence(b"CC", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.12))

    hmm.set_transition(M1, E, log(1.0))

    score = hmm.viterbi(Sequence(b"AC", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.48))

    score = hmm.viterbi(Sequence(b"AA", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.32))


def test_hmm_viterbi_3():
    alphabet = Alphabet(b"AC", b"X")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LPROB_ZERO)

    M1 = NormalState(b"M1", alphabet, [log(0.8), log(0.2)])
    hmm.add_state(M1, LPROB_ZERO)

    D1 = MuteState(b"D1", alphabet)
    hmm.add_state(D1, LPROB_ZERO)

    M2 = NormalState(b"M2", alphabet, [log(0.4), log(0.6)])
    hmm.add_state(M2, LPROB_ZERO)

    D2 = MuteState(b"D2", alphabet)
    hmm.add_state(D2, LPROB_ZERO)

    hmm.set_transition(S, M1, log(0.8))
    hmm.set_transition(S, D1, log(0.2))

    hmm.set_transition(M1, M2, log(0.8))
    hmm.set_transition(M1, D2, log(0.2))

    hmm.set_transition(D1, D2, log(0.2))
    hmm.set_transition(D1, M2, log(0.8))

    hmm.set_transition(D2, E, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()
    hmm.set_transition(E, E, LPROB_ZERO)

    results = hmm.viterbi(Sequence(b"AC", alphabet), E)
    score = results[0].loglikelihood
    assert_equal(results[0].sequence.symbols, b"AC")
    path = results[0].path
    steps = list(path.steps())
    assert_equal(steps[0].seq_len, 0)
    assert_equal(steps[1].seq_len, 1)
    assert_equal(steps[2].seq_len, 1)
    assert_equal(steps[3].seq_len, 0)

    assert_allclose(score, log(0.3072))

    score = hmm.viterbi(Sequence(b"AA", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.2048))

    score = hmm.viterbi(Sequence(b"A", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.128))

    score = hmm.viterbi(Sequence(b"AC", alphabet), E)[0].loglikelihood
    assert_allclose(score, log(0.3072))

    score = hmm.viterbi(Sequence(b"AC", alphabet), M2)[0].loglikelihood
    assert_allclose(score, log(0.3072))

    hmm.del_state(E)

    score = hmm.viterbi(Sequence(b"AC", alphabet), M2)[0].loglikelihood
    assert_allclose(score, log(0.3072))
