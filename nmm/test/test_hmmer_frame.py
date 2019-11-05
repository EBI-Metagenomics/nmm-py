import importlib_resources as pkg_resources
from numpy.testing import assert_allclose, assert_equal

from nmm import create_frame_profile, read_hmmer


def test_read_hmmer_frame_1(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader)

    # most_likely_seq = b"PGKEDNNK"
    most_likely_rna_seq = b"CCU GGU AAA GAA GAU AAU AAC AAA"
    most_likely_rna_seq = most_likely_rna_seq.replace(b" ", b"")

    r = hmmer.lr(most_likely_rna_seq)
    assert_allclose(r.score, 125.83363182422178)
    frags = r.fragments
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence, most_likely_rna_seq)


def test_read_hmmer_frame_2(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader, epsilon=0.1)

    # seq = b"KKKPGKEDNNK"
    rna_seq = b"AAA AAA AAA CCU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")

    r = hmmer.lr(rna_seq)
    assert_allclose(r.score, 168.23071232889802)
    frags = r.fragments
    assert_equal(frags[0].homologous, False)
    assert_equal(frags[0].sequence, b"AAAAAAAAA")
    assert_equal(frags[1].homologous, True)
    assert_equal(frags[1].sequence, b"CCUGGUAAAGAAGAUAAUAACAAA")


def test_read_hmmer_frame_3(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader, epsilon=0.0)

    # seq = b"PGKEDNNK"
    rna_seq = b"CCU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")

    r = hmmer.lr(rna_seq)
    frags = r.fragments
    assert_equal(len(frags), 1)
    assert_equal(frags[0].homologous, True)
    assert_equal(frags[0].sequence, b"CCUGGUAAAGAAGAUAAUAACAAA")


def test_read_hmmer_frame_4(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader, epsilon=0.0)

    # seq = b"PGKEDNNK"
    rna_seq = b"CCUU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")

    r = hmmer.lr(rna_seq)
    frags = r.fragments
    assert_equal(len(frags), 0)


def test_read_hmmer_frame_5(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader, epsilon=0.00001)

    # seq = b"PGKEDNNK"
    rna_seq = b"CCUU GGU AAA GAA GAU AAU AAC AAA"
    rna_seq = rna_seq.replace(b" ", b"")

    r = hmmer.lr(rna_seq)
    frags = r.fragments
    assert_equal(len(frags), 1)
    assert_equal(frags[0].homologous, True)
    assert_equal(frags[0].sequence, b"CCUUGGUAAAGAAGAUAAUAACAAA")


def write_file(path, filename):
    import nmm

    text = pkg_resources.read_text(nmm.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
