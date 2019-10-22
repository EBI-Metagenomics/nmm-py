import pytest
from numpy.testing import assert_allclose
from nmm import MuteState, NormalState, TableState, Alphabet
from nmm import HMM, LOG, Path


def test_hmm_states():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = MuteState("S", alphabet)
    hmm.add_state(S)
    M = TableState("M", alphabet, {"AGU": LOG(0.8), "AGG": LOG(0.2)})
    hmm.add_state(M, LOG(0.0))

    with pytest.raises(ValueError):
        hmm.add_state(S)

    with pytest.raises(ValueError):
        hmm.add_state(M)


def test_hmm_trans_prob():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = MuteState("S", alphabet)
    hmm.add_state(S)

    E = MuteState("E", alphabet)
    hmm.add_state(E)

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_trans(S, E, LOG(0.5))

    assert_allclose(hmm.trans(S, S), LOG(0.0))
    assert_allclose(hmm.trans(S, E), LOG(0.5))
    assert_allclose(hmm.trans(E, S), LOG(0.0))
    assert_allclose(hmm.trans(E, E), LOG(0.0))

    with pytest.raises(ValueError):
        hmm.normalize()

    with pytest.raises(ValueError):
        hmm.normalize()

    hmm.set_start_lprob(S, LOG(0.4))
    hmm.set_trans(E, E, LOG(0.1))

    hmm.normalize()

    assert_allclose(hmm.trans(S, E), LOG(1.0))
    assert_allclose(hmm.trans(E, S), LOG(0.0))
    assert_allclose(hmm.trans(S, S), LOG(0.0))
    assert_allclose(hmm.trans(E, E), LOG(1.0))


def test_hmm_lik_1():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = MuteState("S", alphabet)
    hmm.add_state(S, LOG(1.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    M1 = NormalState(
        "M1", alphabet, {"A": LOG(0.8), "C": LOG(0.2), "G": LOG(0.0), "U": LOG(0.0)}
    )
    hmm.add_state(M1, LOG(0.0))

    M2 = NormalState(
        "M2", alphabet, {"A": LOG(0.4), "C": LOG(0.6), "G": LOG(0.0), "U": LOG(0.6)}
    )
    M2.normalize()
    hmm.add_state(M2, LOG(0.0))

    hmm.set_trans(S, M1, LOG(1.0))
    hmm.set_trans(M1, M2, LOG(1.0))
    hmm.set_trans(M2, E, LOG(1.0))
    hmm.set_trans(E, E, LOG(1.0))
    hmm.normalize()

    p = hmm.likelihood("AC", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.3))

    p = hmm.likelihood("AA", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.2))

    p = hmm.likelihood("AG", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0))

    p = hmm.likelihood("AU", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.3))

    p = hmm.likelihood("CC", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.075))

    p = hmm.likelihood("CA", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.05))

    p = hmm.likelihood("CG", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("CG", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("CU", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.075))

    p = hmm.likelihood("GC", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("GA", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("GG", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("GU", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("UC", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("UA", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("UG", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("UU", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    M3 = NormalState(
        "M2", alphabet, {"A": LOG(0.4), "C": LOG(0.6), "G": LOG(0.0), "U": LOG(0.6)}
    )
    with pytest.raises(ValueError):
        hmm.likelihood("UU", Path([(S, 0), (M1, 1), (M3, 1), (E, 0)]))


def test_hmm_lik_2():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = NormalState(
        "S", alphabet, {"A": LOG(0.8), "C": LOG(0.2), "G": LOG(0.0), "U": LOG(0.0)}
    )
    hmm.add_state(S, LOG(1.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    hmm.set_trans(S, E, LOG(1.0))

    p = hmm.likelihood("A", Path([(S, 1), (E, 0)]))
    assert_allclose(p, LOG(0.8))

    p = hmm.likelihood("C", Path([(S, 1), (E, 0)]))
    assert_allclose(p, LOG(0.2))

    p = hmm.likelihood("A", Path([(S, 1)]))
    assert_allclose(p, LOG(0.8))

    p = hmm.likelihood("C", Path([(S, 1)]))
    assert_allclose(p, LOG(0.2))

    with pytest.raises(ValueError):
        hmm.likelihood("C", Path([(S, 2)]))

    with pytest.raises(ValueError):
        hmm.likelihood("C", Path([(S, 0)]))

    p = hmm.likelihood("C", Path([(E, 1)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("", Path([]))
    assert_allclose(p, LOG(1.0))

    p = hmm.likelihood("A", Path([]))
    assert_allclose(p, LOG(0.0))


def test_hmm_lik_3():
    alphabet = Alphabet("AC")
    hmm = HMM(alphabet)

    S = MuteState("S", alphabet)
    hmm.add_state(S, LOG(1.0))

    M1 = MuteState("M1", alphabet)
    hmm.add_state(M1, LOG(0.0))

    M2 = NormalState("M2", alphabet, {"A": LOG(0.8), "C": LOG(0.2)})
    hmm.add_state(M2, LOG(0.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    hmm.set_trans(S, M1, LOG(1.0))
    hmm.set_trans(M1, M2, LOG(1.0))
    hmm.set_trans(M2, E, LOG(1.0))
    hmm.set_trans(E, E, LOG(1.0))
    hmm.normalize()

    p = hmm.likelihood("A", Path([(S, 1), (E, 0)]))
    assert_allclose(p, LOG(0.0))

    p = hmm.likelihood("A", Path([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.8))

    p = hmm.likelihood("C", Path([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.2))

    with pytest.raises(ValueError):
        p = hmm.likelihood("C", Path([(S, 0), (M1, 1), (M2, 1), (E, 0)]))

    hmm.set_trans(M1, E, LOG(1.0))
    hmm.normalize()

    p = hmm.likelihood("A", Path([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.4))

    p = hmm.likelihood("C", Path([(S, 0), (M1, 0), (M2, 1), (E, 0)]))
    assert_allclose(p, LOG(0.1))

    p = hmm.likelihood("", Path([(S, 0), (M1, 0), (E, 0)]))
    assert_allclose(p, LOG(0.5))


def test_hmm_viterbi_1():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = MuteState("S", alphabet)
    hmm.add_state(S, LOG(1.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    M1 = NormalState(
        "M1", alphabet, {"A": LOG(0.8), "C": LOG(0.2), "G": LOG(0.0), "U": LOG(0.0)}
    )
    hmm.add_state(M1, LOG(0.0))

    M2 = NormalState(
        "M2", alphabet, {"A": LOG(0.4), "C": LOG(0.6), "G": LOG(0.0), "U": LOG(0.6)}
    )
    M2.normalize()
    hmm.add_state(M2, LOG(0.0))

    hmm.set_trans(S, M1, LOG(1.0))
    hmm.set_trans(M1, M2, LOG(1.0))
    hmm.set_trans(M2, E, LOG(1.0))
    hmm.set_trans(E, E, LOG(1.0))
    hmm.normalize()

    with pytest.raises(ValueError):
        hmm.viterbi("AC", E)

    hmm.set_trans(E, E, LOG(0.0))
    assert_allclose(hmm.trans(E, E), LOG(0.0))
    assert_allclose(hmm.trans(S, S), LOG(0.0))
    assert_allclose(hmm.trans(S, E), LOG(0.0))
    assert_allclose(hmm.trans(E, S), LOG(0.0))

    lik = hmm.viterbi("AC", E)
    assert_allclose(lik, LOG(0.3))


def test_hmm_viterbi_2():
    alphabet = Alphabet("AC")

    hmm = HMM(alphabet)
    S = MuteState("S", alphabet)
    hmm.add_state(S, LOG(1.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    M1 = NormalState("M1", alphabet, {"A": LOG(0.8), "C": LOG(0.2)})
    hmm.add_state(M1, LOG(0.0))

    M2 = NormalState("M2", alphabet, {"A": LOG(0.4), "C": LOG(0.6)})
    hmm.add_state(M2, LOG(0.0))

    hmm.set_trans(S, M1, LOG(1.0))
    hmm.set_trans(M1, M2, LOG(1.0))
    hmm.set_trans(M2, E, LOG(1.0))
    hmm.set_trans(E, E, LOG(1.0))
    hmm.normalize()
    hmm.set_trans(E, E, LOG(0.0))

    lik = hmm.viterbi("AC", E)
    assert_allclose(lik, LOG(0.48))

    lik = hmm.viterbi("AA", E)
    assert_allclose(lik, LOG(0.32))

    lik = hmm.viterbi("CA", E)
    assert_allclose(lik, LOG(0.08))

    lik = hmm.viterbi("CC", E)
    assert_allclose(lik, LOG(0.12))

    hmm.set_trans(M1, E, LOG(1.0))

    lik = hmm.viterbi("AC", E)
    assert_allclose(lik, LOG(0.48))

    lik = hmm.viterbi("AA", E)
    assert_allclose(lik, LOG(0.32))


def test_hmm_viterbi_3():
    alphabet = Alphabet("AC")

    hmm = HMM(alphabet)
    S = MuteState("S", alphabet)
    hmm.add_state(S, LOG(1.0))

    E = MuteState("E", alphabet)
    hmm.add_state(E, LOG(0.0))

    M1 = NormalState("M1", alphabet, {"A": LOG(0.8), "C": LOG(0.2)})
    hmm.add_state(M1, LOG(0.0))

    D1 = MuteState("D1", alphabet)
    hmm.add_state(D1, LOG(0.0))

    M2 = NormalState("M2", alphabet, {"A": LOG(0.4), "C": LOG(0.6)})
    hmm.add_state(M2, LOG(0.0))

    D2 = MuteState("D2", alphabet)
    hmm.add_state(D2, LOG(0.0))

    hmm.set_trans(S, M1, LOG(0.8))
    hmm.set_trans(S, D1, LOG(0.2))

    hmm.set_trans(M1, M2, LOG(0.8))
    hmm.set_trans(M1, D2, LOG(0.2))

    hmm.set_trans(D1, D2, LOG(0.2))
    hmm.set_trans(D1, M2, LOG(0.8))

    hmm.set_trans(D2, E, LOG(1.0))
    hmm.set_trans(M2, E, LOG(1.0))
    hmm.set_trans(E, E, LOG(1.0))
    hmm.normalize()
    hmm.set_trans(E, E, LOG(0.0))

    lik = hmm.viterbi("AC", E)
    assert_allclose(lik, LOG(0.3072))

    lik = hmm.viterbi("AA", E)
    assert_allclose(lik, LOG(0.2048))

    lik = hmm.viterbi("A", E)
    assert_allclose(lik, LOG(0.128))

    lik = hmm.viterbi("AC", E)
    assert_allclose(lik, LOG(0.3072))

    lik = hmm.viterbi("AC", M2)
    assert_allclose(lik, LOG(0.3072))

    hmm.del_state(E)

    lik = hmm.viterbi("AC", M2)
    assert_allclose(lik, LOG(0.3072))
