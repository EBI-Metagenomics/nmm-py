from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from imm import lprob_is_zero
from nmm import BaseAlphabet, Codon, CodonProb


def test_codon_prob():
    base = BaseAlphabet.create(b"ACGT", b"X")
    codonp = CodonProb.create(base)

    with pytest.raises(RuntimeError):
        codonp.normalize()

    codonp.set_lprob(Codon.create(b"AAA", base), log(0.01))
    assert_allclose(codonp.get_lprob(Codon.create(b"AAA", base)), log(0.01))

    codonp.normalize()
    assert_allclose(codonp.get_lprob(Codon.create(b"AAA", base)), log(1.0))

    codonp.set_lprob(Codon.create(b"AAA", base), log(0.01))
    assert_allclose(codonp.get_lprob(Codon.create(b"AAA", base)), log(0.01))

    assert_equal(lprob_is_zero(codonp.get_lprob(Codon.create(b"ACA", base))), True)
    with pytest.raises(RuntimeError):
        codonp.get_lprob(Codon.create(b"AXA", base))
