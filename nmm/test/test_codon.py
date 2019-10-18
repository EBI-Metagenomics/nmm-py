from numpy.testing import assert_allclose, assert_equal
from nmm import Codon, LOG, Alphabet


def test_codon():
    bases = Alphabet("ACGT")
    codon = Codon(bases)
    codon.set_lprob("ACT", LOG(0.3))
    codon.set_lprob("CCC", LOG(0.5))
    assert_allclose(codon.get_lprob("ACT"), LOG(0.3))
    assert_allclose(codon.get_lprob("CCC"), LOG(0.5))
    assert_equal(codon.get_lprob("CTC"), LOG(0.0))
    codon.normalize()
    assert_allclose(codon.get_lprob("ACT"), LOG(0.3) - LOG(0.8))
    assert_allclose(codon.get_lprob("CCC"), LOG(0.5) - LOG(0.8))
