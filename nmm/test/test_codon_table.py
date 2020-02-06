from math import log

from numpy.testing import assert_allclose

from nmm import Alphabet, Base, Codon, CodonProb, CodonTable


def test_codon_table():
    base = Base(Alphabet(b"ACGT", b"X"))
    codonp = CodonProb(base)

    codonp.set_lprob(Codon(b"AAA", base), log(0.01))
    codonp.set_lprob(Codon(b"AGA", base), log(0.31))
    codonp.set_lprob(Codon(b"CAA", base), log(0.40))
    codonp.set_lprob(Codon(b"CAT", base), log(0.40))

    codont = CodonTable(codonp)
    assert_allclose(codont.lprob(Codon(b"CAT", base)), log(0.40))
    assert_allclose(codont.lprob(Codon(b"CAX", base)), log(0.80))
    assert_allclose(codont.lprob(Codon(b"XXX", base)), log(1.12))
