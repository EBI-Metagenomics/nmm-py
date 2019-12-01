from math import log, nan

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import HMM, LOG0, Alphabet, CPath, MuteState, NormalState, TableState


def test_hmm_states():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S)
    M = TableState(b"M", alphabet, {b"AGU": log(0.8), b"AGG": log(0.2)})
    hmm.add_state(M, LOG0)

    with pytest.raises(ValueError):
        hmm.add_state(S)

    with pytest.raises(ValueError):
        hmm.add_state(M)

    assert_equal(len(hmm.states()), 2)


def test_hmm_trans_prob():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    with pytest.raises(RuntimeError):
        hmm.set_start_lprob(S, log(0.4))
    hmm.add_state(S)

    E = MuteState(b"E", alphabet)
    with pytest.raises(RuntimeError):
        hmm.transition(S, E)

    with pytest.raises(ValueError):
        hmm.set_transition(S, E, LOG0)

    with pytest.raises(ValueError):
        hmm.set_transition(E, S, LOG0)

    with pytest.raises(ValueError):
        hmm.del_state(E)

    hmm.add_state(E)

    with pytest.raises(RuntimeError):
        hmm.set_transition(E, S, nan)

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_transition(S, E, log(0.5))

    assert_allclose(hmm.transition(S, S), LOG0)
    assert_allclose(hmm.transition(S, E), log(0.5))
    assert_allclose(hmm.transition(E, S), LOG0)
    assert_allclose(hmm.transition(E, E), LOG0)

    with pytest.raises(ValueError):
        hmm.normalize()

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_start_lprob(S, log(0.4))
    hmm.set_transition(E, E, log(0.1))

    hmm.normalize()

    assert_allclose(hmm.transition(S, E), log(1.0))
    assert_allclose(hmm.transition(E, S), LOG0)
    assert_allclose(hmm.transition(S, S), LOG0)
    assert_allclose(hmm.transition(E, E), log(1.0))


def test_hmm_lik_1():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    M1 = NormalState(
        b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2), b"G": LOG0, b"U": LOG0}
    )
    hmm.add_state(M1, LOG0)

    M2 = NormalState(
        b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6), b"G": LOG0, b"U": log(0.6)}
    )
    M2.normalize()
    hmm.add_state(M2, LOG0)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()

    p = hmm.likelihood(b"AC", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(b"AA", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    p = hmm.likelihood(b"AG", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"AU", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(b"CC", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(b"CA", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.05))

    p = hmm.likelihood(b"CG", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"CG", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"CU", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(b"GC", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GA", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GG", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GU", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UC", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UA", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UG", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UU", CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    M3 = NormalState(
        b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6), b"G": LOG0, b"U": log(0.6)}
    )
    with pytest.raises(ValueError):
        hmm.likelihood(b"UU", CPath.create_cpath([(S, 0), (M1, 1), (M3, 1), (E, 0)]))


def test_hmm_lik_2():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = NormalState(
        b"S", alphabet, {b"A": log(0.8), b"C": log(0.2), b"G": LOG0, b"U": LOG0}
    )
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    hmm.set_transition(S, E, log(1.0))

    p = hmm.likelihood(b"A", CPath.create_cpath([(S, 1), (E, 0)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", CPath.create_cpath([(S, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    p = hmm.likelihood(b"A", CPath.create_cpath([(S, 1)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", CPath.create_cpath([(S, 1)]))
    assert_allclose(p, log(0.2))

    with pytest.raises(RuntimeError):
        CPath.create_cpath([(S, 2)])

    with pytest.raises(RuntimeError):
        CPath.create_cpath([(S, 0)])

    with pytest.raises(RuntimeError):
        CPath.create_cpath([(E, 1)])

    with pytest.raises(ValueError):
        hmm.likelihood(b"", CPath.create_cpath([]))

    with pytest.raises(ValueError):
        hmm.likelihood(b"A", CPath.create_cpath([]))


def test_hmm_lik_3():
    alphabet = Alphabet(b"AC")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    M1 = MuteState(b"M1", alphabet)
    hmm.add_state(M1, LOG0)

    M2 = NormalState(b"M2", alphabet, {b"A": log(0.8), b"C": log(0.2)})
    hmm.add_state(M2, LOG0)

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()

    with pytest.raises(RuntimeError):
        CPath.create_cpath([(S, 1), (E, 0)])

    p = hmm.likelihood(b"A", CPath.create_cpath([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", CPath.create_cpath([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    with pytest.raises(RuntimeError):
        CPath.create_cpath([(S, 0), (M1, 1), (M2, 1), (E, 0)])

    hmm.set_transition(M1, E, log(1.0))
    hmm.normalize()

    p = hmm.likelihood(b"A", CPath.create_cpath([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.4))

    p = hmm.likelihood(b"C", CPath.create_cpath([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.1))

    p = hmm.likelihood(b"", CPath.create_cpath([(S, 0), (M1, 0), (E, 0)]))
    assert_allclose(p, log(0.5))


def test_hmm_viterbi_1():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    M1 = NormalState(
        b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2), b"G": LOG0, b"U": LOG0}
    )
    hmm.add_state(M1, LOG0)

    M2 = NormalState(
        b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6), b"G": LOG0, b"U": log(0.6)}
    )
    M2.normalize()
    hmm.add_state(M2, LOG0)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()

    hmm.set_transition(E, E, LOG0)
    assert_allclose(hmm.transition(E, E), LOG0)
    assert_allclose(hmm.transition(S, S), LOG0)
    assert_allclose(hmm.transition(S, E), LOG0)
    assert_allclose(hmm.transition(E, S), LOG0)

    score = hmm.viterbi(b"AC", E)[0]
    assert_allclose(score, log(0.3))


def test_hmm_viterbi_2():
    alphabet = Alphabet(b"AC")

    hmm = HMM(alphabet)
    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    M1 = NormalState(b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2)})
    hmm.add_state(M1, LOG0)

    M2 = NormalState(b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6)})
    hmm.add_state(M2, LOG0)

    hmm.set_transition(S, M1, log(1.0))
    hmm.set_transition(M1, M2, log(1.0))
    hmm.set_transition(M2, E, log(1.0))
    hmm.set_transition(E, E, log(1.0))
    hmm.normalize()
    hmm.set_transition(E, E, LOG0)

    score = hmm.viterbi(b"AC", E)[0]
    assert_allclose(score, log(0.48))

    score = hmm.viterbi(b"AA", E)[0]
    assert_allclose(score, log(0.32))

    score = hmm.viterbi(b"CA", E)[0]
    assert_allclose(score, log(0.08))

    score = hmm.viterbi(b"CC", E)[0]
    assert_allclose(score, log(0.12))

    hmm.set_transition(M1, E, log(1.0))

    score = hmm.viterbi(b"AC", E)[0]
    assert_allclose(score, log(0.48))

    score = hmm.viterbi(b"AA", E)[0]
    assert_allclose(score, log(0.32))


def test_hmm_viterbi_3():
    alphabet = Alphabet(b"AC")

    hmm = HMM(alphabet)
    assert_equal(alphabet, hmm.alphabet)

    S = MuteState(b"S", alphabet)
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    M1 = NormalState(b"M1", alphabet, {b"A": log(0.8), b"C": log(0.2)})
    hmm.add_state(M1, LOG0)

    D1 = MuteState(b"D1", alphabet)
    hmm.add_state(D1, LOG0)

    M2 = NormalState(b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6)})
    hmm.add_state(M2, LOG0)

    D2 = MuteState(b"D2", alphabet)
    hmm.add_state(D2, LOG0)

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
    hmm.set_transition(E, E, LOG0)

    score = hmm.viterbi(b"AC", E)[0]
    assert_allclose(score, log(0.3072))

    score = hmm.viterbi(b"AA", E)[0]
    assert_allclose(score, log(0.2048))

    score = hmm.viterbi(b"A", E)[0]
    assert_allclose(score, log(0.128))

    score = hmm.viterbi(b"AC", E)[0]
    assert_allclose(score, log(0.3072))

    score = hmm.viterbi(b"AC", M2)[0]
    assert_allclose(score, log(0.3072))

    hmm.del_state(E)

    score = hmm.viterbi(b"AC", M2)[0]
    assert_allclose(score, log(0.3072))
