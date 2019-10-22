import importlib_resources as pkg_resources
from numpy.testing import assert_equal, assert_allclose

import nmm


def test_read_hmmer_1(tmp_path):
    filename = "PF03373.hmm"
    text = pkg_resources.read_text(nmm.test, filename)

    with open(tmp_path / filename, "w") as f:
        f.write(text)

    hmmfile = nmm.hmmer.read_file(tmp_path / filename)
    hmm = nmm.hmmer.create_profile(hmmfile)

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
