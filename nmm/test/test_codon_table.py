from math import log

from numpy.testing import assert_allclose

from nmm.alphabet import BaseAlphabet
from nmm.codon import Codon
from nmm.prob import CodonProb, CodonTable


def test_codon_table():
    base = BaseAlphabet.create(b"ACGT", b"X")
    codonp = CodonProb(base)

    codonp.set_lprob(Codon.create(b"AAA", base), log(0.01))
    codonp.set_lprob(Codon.create(b"AGA", base), log(0.31))
    codonp.set_lprob(Codon.create(b"CAA", base), log(0.40))
    codonp.set_lprob(Codon.create(b"CAT", base), log(0.40))

    codont = CodonTable(codonp)
    assert_allclose(codont.lprob(Codon.create(b"CAT", base)), log(0.40))
    assert_allclose(codont.lprob(Codon.create(b"CAX", base)), log(0.80))
    assert_allclose(codont.lprob(Codon.create(b"XXX", base)), log(1.12))
