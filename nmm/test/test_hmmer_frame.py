import importlib_resources as pkg_resources
from numpy.testing import assert_allclose, assert_equal

from nmm import create_frame_profile, read_hmmer


def test_read_hmmer_frame_unihit_homologous_1(tmp_path):
    filepath = write_file(tmp_path, "PF03373.hmm")
    reader = read_hmmer(filepath)
    hmmer = create_frame_profile(reader)
    # most_likely_seq = "PGKEDNNK"
    # r = hmmer.lr(most_likely_seq)
    # assert_allclose(r.score, 11.867796719423442)
    # frags = r.fragments
    # assert_equal(len(frags), 1)
    # frag = frags[0]
    # assert_equal(frag.homologous, True)
    # assert_equal(frag.sequence, most_likely_seq)


def write_file(path, filename):
    import nmm

    text = pkg_resources.read_text(nmm.test, filename)

    with open(path / filename, "w") as f:
        f.write(text)

    return path / filename
