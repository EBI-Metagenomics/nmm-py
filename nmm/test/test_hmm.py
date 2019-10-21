from math import isinf, nan
import pytest
from numpy.testing import assert_allclose, assert_equal
from nmm import MuteState, NormalState, TableState, FrameState, Alphabet
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

    M1 = NormalState("M1", alphabet, [LOG(0.8), LOG(0.2), LOG(0.0), LOG(0.0)])
    hmm.add_state(M1, LOG(0.0))

    M2 = NormalState("M2", alphabet, [LOG(0.4), LOG(0.6), LOG(0.0), LOG(0.6)])
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

    M3 = NormalState("M2", alphabet, [LOG(0.4), LOG(0.6), LOG(0.0), LOG(0.6)])
    with pytest.raises(ValueError):
        hmm.likelihood("UU", Path([(S, 0), (M1, 1), (M3, 1), (E, 0)]))


def test_hmm_lik_2():
    alphabet = Alphabet("ACGU")
    hmm = HMM(alphabet)

    S = NormalState("S", alphabet, [LOG(0.8), LOG(0.2), LOG(0.0), LOG(0.0)])
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

    M2 = NormalState("M2", alphabet, [LOG(0.8), LOG(0.2)])
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


# def test_hmm_viterbi_1():
#     alphabet = "ACGU"

#     hmm = HMM(alphabet)
#     start_state = MuteState("S", alphabet)
#     hmm.add_state(start_state, LOG(1.0))

#     end_state = MuteState("E", alphabet)
#     hmm.add_state(end_state, LOG(0.0))

#     M1 = NormalState("M1", {"A": LOG(0.8), "C": LOG(0.2), "G": LOG(0.0), "U": LOG(0.0)})
#     hmm.add_state(M1, LOG(0.0))

#     M2 = NormalState("M2", {"A": LOG(0.4), "C": LOG(0.6), "G": LOG(0.0), "U": LOG(0.6)})
#     hmm.add_state(M2, LOG(0.0))

#     hmm.set_trans("S", "M1", LOG(1.0))
#     hmm.set_trans("M1", "M2", LOG(1.0))
#     hmm.set_trans("M2", "E", LOG(1.0))
#     hmm.normalize()

#     lik, path = hmm.viterbi("AC", "E")
#     assert_allclose(lik, 0.3)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert path[0][1] == 0
#     assert path[1][1] == 1
#     assert path[2][1] == 1
#     assert path[3][1] == 0

#     lik = hmm.viterbi("AC", "E", True)[0]
#     assert_allclose(lik, LOG(0.3))

#     p = hmm.likelihood("AC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.3)
#     logp = hmm.likelihood("AC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)], True)
#     assert_allclose(logp, LOG(0.3))


# def test_hmm_viterbi_2():
#     alphabet = "AC"

#     hmm = HMM(alphabet)
#     start_state = MuteState("S", alphabet)
#     hmm.add_state(start_state, LOG(1.0))

#     end_state = MuteState("E", alphabet)
#     hmm.add_state(end_state, LOG(0.0))

#     M1 = NormalState("M1", {"A": LOG(0.8), "C": LOG(0.2)})
#     hmm.add_state(M1, LOG(0.0))

#     M2 = NormalState("M2", {"A": LOG(0.4), "C": LOG(0.6)})
#     hmm.add_state(M2, LOG(0.0))

#     hmm.set_trans("S", "M1", LOG(1.0))
#     hmm.set_trans("M1", "M2", LOG(1.0))
#     hmm.set_trans("M2", "E", LOG(1.0))
#     hmm.normalize()

#     lik, path = hmm.viterbi("AC", "E")
#     assert_allclose(lik, 0.48)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.48)

#     lik, path = hmm.viterbi("AA", "E")
#     assert_allclose(lik, 0.32)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AA", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.32)

#     lik, path = hmm.viterbi("CA", "E")
#     assert_allclose(lik, 0.08)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("CA", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.08)

#     lik, path = hmm.viterbi("CC", "E")
#     assert_allclose(lik, 0.12)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("CC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.12)

#     hmm.set_trans("M1", "E", LOG(1.0))
#     hmm.normalize()

#     lik, path = hmm.viterbi("AC", "E")
#     assert_allclose(lik, 0.48)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.24)

#     lik, path = hmm.viterbi("AA", "E")
#     assert_allclose(lik, 0.32)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AA", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.16)

#     lik, path = hmm.viterbi("CA", "E")
#     assert_allclose(lik, 0.08)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("CA", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.04)

#     lik, path = hmm.viterbi("CC", "E")
#     assert_allclose(lik, 0.12)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("CC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.06)


# def test_hmm_viterbi_3():
#     alphabet = "AC"

#     hmm = HMM(alphabet)
#     start_state = MuteState("S", alphabet)
#     hmm.add_state(start_state, LOG(1.0))

#     end_state = MuteState("E", alphabet)
#     hmm.add_state(end_state, LOG(0.0))

#     M1 = NormalState("M1", {"A": LOG(0.8), "C": LOG(0.2)})
#     hmm.add_state(M1, LOG(0.0))

#     D1 = MuteState("D1", alphabet)
#     hmm.add_state(D1, LOG(0.0))

#     M2 = NormalState("M2", {"A": LOG(0.4), "C": LOG(0.6)})
#     hmm.add_state(M2, LOG(0.0))

#     D2 = MuteState("D2", alphabet)
#     hmm.add_state(D2, LOG(0.0))

#     hmm.set_trans("S", "M1", LOG(0.8))
#     hmm.set_trans("S", "D1", LOG(0.2))

#     hmm.set_trans("M1", "M2", LOG(0.8))
#     hmm.set_trans("M1", "D2", LOG(0.2))

#     hmm.set_trans("D1", "D2", LOG(0.2))
#     hmm.set_trans("D1", "M2", LOG(0.8))

#     hmm.set_trans("D2", "E", LOG(1.0))
#     hmm.set_trans("M2", "E", LOG(1.0))
#     hmm.normalize()

#     lik, path = hmm.viterbi("AC", "E")
#     assert_allclose(lik, 0.3072)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AC", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.3072)

#     lik, path = hmm.viterbi("AA", "E")
#     assert_allclose(lik, 0.2048)
#     assert [p[0] for p in path] == ["S", "M1", "M2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 1, 0]
#     p = hmm.likelihood("AA", [("S", 0), ("M1", 1), ("M2", 1), ("E", 0)])
#     assert_allclose(p, 0.2048)

#     lik, path = hmm.viterbi("A", "E")
#     assert_allclose(lik, 0.128)
#     assert [p[0] for p in path] == ["S", "M1", "D2", "E"]
#     assert list(p[1] for p in path) == [0, 1, 0, 0]
#     p = hmm.likelihood("A", [("S", 0), ("M1", 1), ("D2", 0), ("E", 0)])
#     assert_allclose(p, 0.128)


# def test_hmm_draw(tmp_path):
#     hmm = _create_hmm()

#     base_emission = {"A": LOG(0.5), "C": LOG(0.5)}
#     codon_emission = {"ACC": LOG(0.8), "AAA": LOG(0.2)}
#     epsilon = 0.1
#     M0 = FrameState("M0", base_emission, codon_emission, epsilon)
#     hmm.add_state(M0, LOG(1.0))
#     hmm.set_trans("M0", "E", LOG(1.0))
#     hmm.normalize()

#     hmm.draw(tmp_path / "test.pdf", emissions=5, init_prob=True)
#     hmm.draw(tmp_path / "test.pdf", emissions=3, init_prob=False)
#     hmm.draw(tmp_path / "test.pdf", emissions=0, init_prob=True)
#     hmm.draw(tmp_path / "test.pdf", emissions=0, init_prob=False)
#     hmm.draw(tmp_path / "test.pdf", emissions=50, init_prob=True)


# def test_hmm_rename_state():
#     hmm = _create_hmm()

#     with pytest.raises(ValueError):
#         hmm.rename_state("DD", "D1")

#     with pytest.raises(ValueError):
#         hmm.rename_state("E", "M1")

#     p = hmm.trans("D1", "D2")
#     hmm.rename_state("S", "B")
#     hmm.rename_state("D2", "DD")
#     assert "B" in hmm.states
#     assert "S" not in hmm.states
#     assert "DD" in hmm.states
#     assert "D2" not in hmm.states
#     assert hmm.trans("D1", "DD") == p


# def test_hmm_delete_state():
#     hmm = _create_hmm()

#     with pytest.raises(ValueError):
#         hmm.delete_state("D22")

#     hmm.delete_state("D2")
#     assert "D2" not in hmm.states


# def test_hmm_single_state():
#     alphabet = "ACGU"
#     hmm = HMM(alphabet)
#     state = NormalState(
#         "I", {"A": LOG(0.8), "C": LOG(0.2), "G": LOG(0.0), "U": LOG(0.0)}
#     )
#     hmm.add_state(state)
#     hmm.set_trans("I", "I", LOG(1.0))
#     hmm.normalize()
#     lik, path = hmm.viterbi("ACC", "I")
#     assert abs(lik - 0.032) < 1e-7
#     assert path == [("I", 1), ("I", 1), ("I", 1)]


# def _create_hmm():
#     alphabet = "AC"

#     hmm = HMM(alphabet)
#     start_state = MuteState("S", alphabet)
#     hmm.add_state(start_state, LOG(1.0))

#     end_state = MuteState("E", alphabet)
#     hmm.add_state(end_state, LOG(0.0))

#     M1 = NormalState("M1", {"A": LOG(0.8), "C": LOG(0.2)})
#     hmm.add_state(M1, LOG(0.0))

#     D1 = MuteState("D1", alphabet)
#     hmm.add_state(D1, LOG(0.0))

#     M2 = NormalState("M2", {"A": LOG(0.4), "C": LOG(0.6)})
#     hmm.add_state(M2, LOG(0.0))

#     D2 = MuteState("D2", alphabet)
#     hmm.add_state(D2, LOG(0.0))

#     hmm.set_trans("S", "M1", LOG(0.8))
#     hmm.set_trans("S", "D1", LOG(0.2))

#     hmm.set_trans("M1", "M2", LOG(0.8))
#     hmm.set_trans("M1", "D2", LOG(0.2))

#     hmm.set_trans("D1", "D2", LOG(0.2))
#     hmm.set_trans("D1", "M2", LOG(0.8182787382))

#     hmm.set_trans("D2", "E", LOG(1.0))
#     hmm.set_trans("M2", "E", LOG(1.0))

#     return hmm
