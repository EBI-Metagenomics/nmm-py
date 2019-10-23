import pytest
import importlib_resources as pkg_resources
from numpy.testing import assert_equal, assert_allclose

import nmm
from nmm import LOG0


def test_hmmer_create_profile(PF03373):
    hmm = nmm.hmmer.create_profile(PF03373)

    assert_equal(set("ACDEFGHIKLMNPQRSTVWY"), set(hmm.alphabet.symbols))

    path = [(hmm.find_state("B"), 0)]
    path += [(hmm.find_state(f"M{i}"), 1) for i in range(1, 9)]
    path += [(hmm.find_state("E"), 0)]
    most_likely_seq = "PGKEDNNK"
    lik = hmm.likelihood(most_likely_seq, nmm.Path(path))
    assert_allclose(lik, -3.910270475034747)
    score = hmm.viterbi(most_likely_seq, hmm.find_state("E"))
    assert_allclose(score, -3.9102704750347463)

    seq = "PGKEDNNSQ"
    path = [
        ("B", 0),
        ("M1", 1),
        ("M2", 1),
        ("M3", 1),
        ("M4", 1),
        ("M5", 1),
        ("M6", 1),
        ("M7", 1),
        ("I7", 1),
        ("M8", 1),
        ("E", 0),
    ]
    path = [(hmm.find_state(step[0]), step[1]) for step in path]

    assert_allclose(hmm.likelihood(seq, nmm.Path(path)), -14.730593544566743)


def test_hmmer_create_core_profile(PF03373):
    hmm = nmm.hmmer.create_core_profile(PF03373)

    assert_equal(set("ACDEFGHIKLMNPQRSTVWY"), set(hmm.alphabet.symbols))

    path = [(hmm.find_state("B"), 0)]
    path += [(hmm.find_state(f"M{i}"), 1) for i in range(1, 2)]
    lik = hmm.likelihood("P", nmm.Path(path))
    assert_allclose(lik, -0.09290742187949368)

    path = [(hmm.find_state("B"), 0)]
    path += [(hmm.find_state(f"M{i}"), 1) for i in range(1, 9)]
    path += [(hmm.find_state("E"), 0)]
    most_likely_seq = "PGKEDNNK"
    lik = hmm.likelihood(most_likely_seq, nmm.Path(path))
    assert_allclose(lik, -3.910270475034747)
    score = hmm.viterbi(most_likely_seq, hmm.find_state("E"))
    assert_allclose(score, -3.9102704750347463)

    seq = "PGKEDNNSQ"
    path = [
        ("B", 0),
        ("M1", 1),
        ("M2", 1),
        ("M3", 1),
        ("M4", 1),
        ("M5", 1),
        ("M6", 1),
        ("M7", 1),
        ("I7", 1),
        ("M8", 1),
        ("E", 0),
    ]
    path = [(hmm.find_state(step[0]), step[1]) for step in path]

    assert_allclose(hmm.likelihood(seq, nmm.Path(path)), -14.730593544566743)


def test_hmmer_core_profile_bg(PF03373):
    hmm = nmm.hmmer.create_core_profile(PF03373)

    bg_hmm = nmm.HMM(hmm.alphabet)
    N = hmm.find_state("I1")
    bg_hmm.add_state(N, 0.0, "N")
    bg_hmm.set_trans(N, N, 0.0)
    most_likely_seq = "PGKEDNNK"
    score = bg_hmm.viterbi(most_likely_seq, N)
    assert_allclose(score, -21.810690069999648)


def test_hmmer_global_profile(tmp_path):
    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    hmmer = nmm.hmmer.read_file2(tmp_path / filename)
    most_likely_seq = "PGKEDNNK"

    path = [(hmmer._special_states["S"], 0)]
    lik = hmmer.hmm.likelihood("", nmm.Path(path))
    assert_allclose(lik, 0.0, atol=1e-7)

    path = [(hmmer._special_states["S"], 0)]
    path += [(hmmer._special_states["B"], 0)]
    lik = hmmer.hmm.likelihood("", nmm.Path(path))
    assert_allclose(lik, 0.0, atol=1e-7)

    path = [(hmmer._special_states["S"], 0)]
    path += [(hmmer._special_states["N"], 0)]
    lik = hmmer.hmm.likelihood("", nmm.Path(path))
    assert_equal(lik, LOG0)

    path = [(hmmer._special_states["S"], 0)]
    path += [(hmmer._special_states["N"], 1)]
    with pytest.raises(ValueError):
        hmmer.hmm.likelihood("", nmm.Path(path))

    path = [(hmmer._special_states["S"], 0)]
    path += [(hmmer._special_states["B"], 0)]
    path += [(hmmer._core_nodes[0].M, 1)]
    lik = hmmer.hmm.likelihood(most_likely_seq[:1], nmm.Path(path))
    assert_allclose(lik, -0.08208032322059217)

    path = [(hmmer._special_states["S"], 0)]
    path += [(hmmer._special_states["B"], 0)]
    path += [(hmmer._core_nodes[i].M, 1) for i in range(0, 8)]
    path += [(hmmer._special_states["E"], 0)]
    path += [(hmmer._special_states["T"], 0)]

    lik = hmmer.hmm.likelihood(most_likely_seq, nmm.Path(path))
    assert_allclose(lik, -3.892142512894315)

    score = hmmer.viterbi(most_likely_seq)
    assert_allclose(score, -3.8921425128943143)

    lr = hmmer.lr(most_likely_seq)
    print(lr)

    lr = hmmer.lr("ACDEFGHI")
    print(lr)

    # loc, scale = (-11.570078552677991, 2.579038545186777)
    loc, scale = (-5.4106, 0.77110)
    pv = hmmer.gumbel(most_likely_seq, loc, scale)
    print(pv)
    pv = hmmer.gumbel("ACDEFGHI", loc, scale)
    print(pv)
    pv = hmmer.gumbel("ACDEFGHI"[::-1], loc, scale)
    print(pv)


@pytest.fixture
def PF03373(tmp_path):
    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    return nmm.hmmer.read_file(tmp_path / filename)
