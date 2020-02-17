from numpy.testing import assert_equal

from nmm import Alphabet, BaseAlphabet, codon_iter


def test_codon_iter():
    base = BaseAlphabet(Alphabet(b"ACGT", b"X"))

    codons = list(codon_iter(base))
    assert_equal(len(codons), 64)
    assert_equal(codons[0].symbols, b"AAA")
    assert_equal(codons[1].symbols, b"AAC")
