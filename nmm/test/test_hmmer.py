import importlib_resources as pkg_resources
from numpy.testing import assert_allclose, assert_equal

from nmm import read_hmmer, create_hmmer_profile


def test_read_hmmer_unihit_homologous_1(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_hmmer_profile(reader)
    most_likely_seq = b"PGKEDNNK"
    r = hmmer.lr(most_likely_seq)
    assert_allclose(r.score, 11.867796719423442)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence, most_likely_seq)


def test_read_hmmer_unihit_homologous_2(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_hmmer_profile(reader)
    seq = b"PGKENNK"
    r = hmmer.lr(seq)
    assert_allclose(r.score, 3.299501501364073)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence, seq)


def test_read_hmmer_unihit_homologous_3(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_hmmer_profile(reader)
    seq = b"PGKEPNNK"
    r = hmmer.lr(seq)
    assert_allclose(r.score, 6.883636719423446)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence, seq)


def test_read_hmmer_nonhomo_and_homologous(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_hmmer_profile(reader)
    seq = b"KKKPGKEDNNK"
    r = hmmer.lr(seq)
    assert_allclose(r.score, 10.707618955640605)
    frags = r.fragments
    assert_equal(len(frags), 2)
    assert_equal(frags[0].homologous, False)
    assert_equal(frags[0].sequence, b"KKK")
    assert_equal(frags[1].homologous, True)
    assert_equal(frags[1].sequence, b"PGKEDNNK")


def test_read_hmmer_multihit_homologous(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_hmmer_profile(reader)
    seq = b"PPPPGKEDNNKDDDPGKEDNNKEEEE"
    r = hmmer.lr(seq)
    assert_allclose(r.score, 20.329227532144742)
    frags = r.fragments
    assert_equal(len(frags), 4)
    assert_equal(frags[0].homologous, False)
    assert_equal(frags[0].sequence, b"PPP")
    assert_equal(frags[1].homologous, True)
    assert_equal(frags[1].sequence, b"PGKEDNNK")
    assert_equal(frags[2].homologous, False)
    assert_equal(frags[2].sequence, b"DDD")
    assert_equal(frags[3].homologous, True)
    assert_equal(frags[3].sequence, b"PGKEDNNK")


def write_file(path, filename):
    import nmm

    text = pkg_resources.read_text(nmm.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
