from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import Alphabet, Base, Codon, CodonProb, lprob_is_zero


def test_codon_prob():
    base = Base(Alphabet(b"ACGT", b"X"))
    codonp = CodonProb(base)

    with pytest.raises(RuntimeError):
        codonp.normalize()

    codonp.set_lprob(Codon(b"AAA", base), log(0.01))
    assert_allclose(codonp.get_lprob(Codon(b"AAA", base)), log(0.01))

    codonp.normalize()
    assert_allclose(codonp.get_lprob(Codon(b"AAA", base)), log(1.0))

    codonp.set_lprob(Codon(b"AAA", base), log(0.01))
    assert_allclose(codonp.get_lprob(Codon(b"AAA", base)), log(0.01))

    assert_equal(lprob_is_zero(codonp.get_lprob(Codon(b"ACA", base))), True)
    with pytest.raises(RuntimeError):
        codonp.get_lprob(Codon(b"AXA", base))