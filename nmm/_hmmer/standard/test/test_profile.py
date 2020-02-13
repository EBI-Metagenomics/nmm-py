from numpy.testing import assert_allclose, assert_equal

from hmmer_reader import open_hmmer

from nmm import create_standard_profile, Sequence


def test_standard_profile_unihit_homologous_1(PF03373):
    with open_hmmer(PF03373) as reader:
        hmmer = create_standard_profile(reader.read_profile())

    alphabet = hmmer.alphabet
    most_likely_seq = Sequence(b"PGKEDNNK", alphabet)
    r = hmmer.search(most_likely_seq)

    assert_allclose(r.loglikelihood, 11.867796719423442)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence.symbols, most_likely_seq.symbols)

    hmmer.multiple_hits = False
    r = hmmer.search(most_likely_seq)
    assert_allclose(r.loglikelihood, 11.94063404337571)
    frags = r.fragments
    assert_equal(len(frags), 1)
    frag = frags[0]
    assert_equal(frag.homologous, True)
    assert_equal(frag.sequence.symbols, most_likely_seq.symbols)


# def test_standard_profile_unihit_homologous_2(PF03373):
#     with open_hmmer(PF03373) as reader:
#         hmmer = create_standard_profile(reader.read_profile())
#     seq = b"PGKENNK"
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 3.299501501364073)
#     frags = r.fragments
#     assert_equal(len(frags), 1)
#     frag = frags[0]
#     assert_equal(frag.homologous, True)
#     assert_equal(frag.sequence, seq)
#     assert_equal(str(frag), "<StandardFragment:PGKENNK>")


# def test_standard_profile_unihit_homologous_3(PF03373):
#     with open_hmmer(PF03373) as reader:
#         hmmer = create_standard_profile(reader.read_profile())
#     seq = b"PGKEPNNK"
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 6.883636719423446)
#     frags = r.fragments
#     assert_equal(len(frags), 1)
#     frag = frags[0]
#     assert_equal(frag.homologous, True)
#     assert_equal(frag.sequence, seq)


# def test_standard_profile_nonhomo_and_homologous(PF03373):
#     with open_hmmer(PF03373) as reader:
#         hmmer = create_standard_profile(reader.read_profile())
#     seq = b"KKKPGKEDNNK"
#     assert_equal(hmmer.multiple_hits, True)
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 10.707618955640605)
#     frags = r.fragments
#     assert_equal(len(frags), 2)
#     assert_equal(frags[0].homologous, False)
#     assert_equal(frags[0].sequence, b"KKK")
#     assert_equal(frags[1].homologous, True)
#     assert_equal(frags[1].sequence, b"PGKEDNNK")

#     hmmer.multiple_hits = False
#     assert_equal(hmmer.multiple_hits, False)
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 10.96037578075283)
#     frags = r.fragments
#     assert_equal(len(frags), 2)
#     assert_equal(frags[0].homologous, False)
#     assert_equal(frags[0].sequence, b"KKK")
#     assert_equal(frags[1].homologous, True)
#     assert_equal(frags[1].sequence, b"PGKEDNNK")


# def test_standard_profile_multihit_homologous1(PF03373):
#     with open_hmmer(PF03373) as reader:
#         hmmer = create_standard_profile(reader.read_profile())
#     seq = b"PPPPGKEDNNKDDDPGKEDNNKEEEE"
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 20.329227532144742)
#     frags = r.fragments
#     assert_equal(len(frags), 5)
#     assert_equal(frags[0].homologous, False)
#     assert_equal(frags[0].sequence, b"PPP")
#     assert_equal(frags[1].homologous, True)
#     assert_equal(frags[1].sequence, b"PGKEDNNK")
#     assert_equal(frags[2].homologous, False)
#     assert_equal(frags[2].sequence, b"DDD")
#     assert_equal(frags[3].homologous, True)
#     assert_equal(frags[3].sequence, b"PGKEDNNK")
#     assert_equal(frags[4].homologous, False)
#     assert_equal(frags[4].sequence, b"EEEE")

#     items = list(frags[0].items())

#     assert_equal(items[0][0], b"")
#     assert_equal(str(items[0][1]), "<S,0>")
#     assert_equal(items[1][0], b"P")
#     assert_equal(str(items[1][1]), "<N,1>")
#     assert_equal(items[2][0], b"P")
#     assert_equal(str(items[2][1]), "<N,1>")
#     assert_equal(items[3][0], b"P")
#     assert_equal(str(items[3][1]), "<N,1>")
#     assert_equal(items[4][0], b"")
#     assert_equal(str(items[4][1]), "<B,0>")

#     items = list(frags[1].items())

#     assert_equal(items[0][0], b"P")
#     assert_equal(str(items[0][1]), "<M1,1>")
#     assert_equal(items[1][0], b"G")
#     assert_equal(str(items[1][1]), "<M2,1>")
#     assert_equal(items[2][0], b"K")
#     assert_equal(str(items[2][1]), "<M3,1>")
#     assert_equal(items[3][0], b"E")
#     assert_equal(str(items[3][1]), "<M4,1>")
#     assert_equal(items[4][0], b"D")
#     assert_equal(str(items[4][1]), "<M5,1>")
#     assert_equal(items[5][0], b"N")
#     assert_equal(str(items[5][1]), "<M6,1>")
#     assert_equal(items[6][0], b"N")
#     assert_equal(str(items[6][1]), "<M7,1>")
#     assert_equal(items[7][0], b"K")
#     assert_equal(str(items[7][1]), "<M8,1>")

#     hmmer.multiple_hits = False
#     r = hmmer.search(seq)
#     assert_allclose(r.loglikelihood, 8.666478660222928)
#     frags = r.fragments
#     assert_equal(len(frags), 3)
#     assert_equal(frags[0].homologous, False)
#     assert_equal(frags[1].homologous, True)
#     assert_equal(frags[1].sequence, b"PGKEDNNK")
#     assert_equal(frags[2].homologous, False)
