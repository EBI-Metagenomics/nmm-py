from math import log, nan

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import HMM, LOG0, Alphabet, MuteState, NormalState, Path, TableState


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

    assert_equal(len(hmm.states), 2)


def test_hmm_trans_prob():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = MuteState(b"S", alphabet)
    with pytest.raises(RuntimeError):
        hmm.set_start_lprob(S, log(0.4))
    hmm.add_state(S)

    E = MuteState(b"E", alphabet)
    with pytest.raises(RuntimeError):
        hmm.trans(S, E)

    with pytest.raises(ValueError):
        hmm.set_trans(S, E, LOG0)

    with pytest.raises(ValueError):
        hmm.set_trans(E, S, LOG0)

    with pytest.raises(ValueError):
        hmm.del_state(E)

    hmm.add_state(E)

    with pytest.raises(RuntimeError):
        hmm.set_trans(E, S, nan)

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_trans(S, E, log(0.5))

    assert_allclose(hmm.trans(S, S), LOG0)
    assert_allclose(hmm.trans(S, E), log(0.5))
    assert_allclose(hmm.trans(E, S), LOG0)
    assert_allclose(hmm.trans(E, E), LOG0)

    with pytest.raises(ValueError):
        hmm.normalize()

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_start_lprob(S, log(0.4))
    hmm.set_trans(E, E, log(0.1))

    hmm.normalize()

    assert_allclose(hmm.trans(S, E), log(1.0))
    assert_allclose(hmm.trans(E, S), LOG0)
    assert_allclose(hmm.trans(S, S), LOG0)
    assert_allclose(hmm.trans(E, E), log(1.0))


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

    hmm.set_trans(S, M1, log(1.0))
    hmm.set_trans(M1, M2, log(1.0))
    hmm.set_trans(M2, E, log(1.0))
    hmm.set_trans(E, E, log(1.0))
    hmm.normalize()

    p = hmm.likelihood(b"AC", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(b"AA", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    p = hmm.likelihood(b"AG", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"AU", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.3))

    p = hmm.likelihood(b"CC", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(b"CA", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.05))

    p = hmm.likelihood(b"CG", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"CG", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"CU", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.075))

    p = hmm.likelihood(b"GC", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GA", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GG", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"GU", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UC", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UA", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UG", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    p = hmm.likelihood(b"UU", Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG0)

    M3 = NormalState(
        b"M2", alphabet, {b"A": log(0.4), b"C": log(0.6), b"G": LOG0, b"U": log(0.6)}
    )
    with pytest.raises(ValueError):
        hmm.likelihood(b"UU", Path.create([(S, 0), (M1, 1), (M3, 1), (E, 0)]))


def test_hmm_lik_2():
    alphabet = Alphabet(b"ACGU")
    hmm = HMM(alphabet)

    S = NormalState(
        b"S", alphabet, {b"A": log(0.8), b"C": log(0.2), b"G": LOG0, b"U": LOG0}
    )
    hmm.add_state(S, log(1.0))

    E = MuteState(b"E", alphabet)
    hmm.add_state(E, LOG0)

    hmm.set_trans(S, E, log(1.0))

    p = hmm.likelihood(b"A", Path.create([(S, 1), (E, 0)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", Path.create([(S, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    p = hmm.likelihood(b"A", Path.create([(S, 1)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", Path.create([(S, 1)]))
    assert_allclose(p, log(0.2))

    with pytest.raises(RuntimeError):
        Path.create([(S, 2)])

    with pytest.raises(RuntimeError):
        Path.create([(S, 0)])

    with pytest.raises(RuntimeError):
        Path.create([(E, 1)])

    with pytest.raises(ValueError):
        hmm.likelihood(b"", Path.create([]))

    with pytest.raises(ValueError):
        hmm.likelihood(b"A", Path.create([]))


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

    hmm.set_trans(S, M1, log(1.0))
    hmm.set_trans(M1, M2, log(1.0))
    hmm.set_trans(M2, E, log(1.0))
    hmm.set_trans(E, E, log(1.0))
    hmm.normalize()

    with pytest.raises(RuntimeError):
        Path.create([(S, 1), (E, 0)])

    p = hmm.likelihood(b"A", Path.create([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.8))

    p = hmm.likelihood(b"C", Path.create([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.2))

    with pytest.raises(RuntimeError):
        Path.create([(S, 0), (M1, 1), (M2, 1), (E, 0)])

    hmm.set_trans(M1, E, log(1.0))
    hmm.normalize()

    p = hmm.likelihood(b"A", Path.create([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.4))

    p = hmm.likelihood(b"C", Path.create([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, log(0.1))

    p = hmm.likelihood(b"", Path.create([(S, 0), (M1, 0), (E, 0)]))
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

    hmm.set_trans(S, M1, log(1.0))
    hmm.set_trans(M1, M2, log(1.0))
    hmm.set_trans(M2, E, log(1.0))
    hmm.set_trans(E, E, log(1.0))
    hmm.normalize()

    hmm.set_trans(E, E, LOG0)
    assert_allclose(hmm.trans(E, E), LOG0)
    assert_allclose(hmm.trans(S, S), LOG0)
    assert_allclose(hmm.trans(S, E), LOG0)
    assert_allclose(hmm.trans(E, S), LOG0)

    lik = hmm.viterbi(b"AC", E)
    assert_allclose(lik.score, log(0.3))


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

    hmm.set_trans(S, M1, log(1.0))
    hmm.set_trans(M1, M2, log(1.0))
    hmm.set_trans(M2, E, log(1.0))
    hmm.set_trans(E, E, log(1.0))
    hmm.normalize()
    hmm.set_trans(E, E, LOG0)

    lik = hmm.viterbi(b"AC", E)
    assert_allclose(lik.score, log(0.48))

    lik = hmm.viterbi(b"AA", E)
    assert_allclose(lik.score, log(0.32))

    lik = hmm.viterbi(b"CA", E)
    assert_allclose(lik.score, log(0.08))

    lik = hmm.viterbi(b"CC", E)
    assert_allclose(lik.score, log(0.12))

    hmm.set_trans(M1, E, log(1.0))

    lik = hmm.viterbi(b"AC", E)
    assert_allclose(lik.score, log(0.48))

    lik = hmm.viterbi(b"AA", E)
    assert_allclose(lik.score, log(0.32))


def test_hmm_viterbi_3():
    alphabet = Alphabet(b"AC")

    hmm = HMM(alphabet)
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

    hmm.set_trans(S, M1, log(0.8))
    hmm.set_trans(S, D1, log(0.2))

    hmm.set_trans(M1, M2, log(0.8))
    hmm.set_trans(M1, D2, log(0.2))

    hmm.set_trans(D1, D2, log(0.2))
    hmm.set_trans(D1, M2, log(0.8))

    hmm.set_trans(D2, E, log(1.0))
    hmm.set_trans(M2, E, log(1.0))
    hmm.set_trans(E, E, log(1.0))
    hmm.normalize()
    hmm.set_trans(E, E, LOG0)

    lik = hmm.viterbi(b"AC", E)
    assert_allclose(lik.score, log(0.3072))

    lik = hmm.viterbi(b"AA", E)
    assert_allclose(lik.score, log(0.2048))

    lik = hmm.viterbi(b"A", E)
    assert_allclose(lik.score, log(0.128))

    lik = hmm.viterbi(b"AC", E)
    assert_allclose(lik.score, log(0.3072))

    lik = hmm.viterbi(b"AC", M2)
    assert_allclose(lik.score, log(0.3072))

    hmm.del_state(E)

    lik = hmm.viterbi(b"AC", M2)
    assert_allclose(lik.score, log(0.3072))
