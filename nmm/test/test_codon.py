from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, Codon


def test_codon():
    bases = Alphabet("ACGT")
    codon = Codon(bases)
    codon.set_lprob("ACT", log(0.3))
    codon.set_lprob("CCC", log(0.5))
    assert_allclose(codon.get_lprob("ACT"), log(0.3))
    assert_allclose(codon.get_lprob("CCC"), log(0.5))
    assert_equal(codon.get_lprob("CTC"), LOG0)
    codon.normalize()
    assert_allclose(codon.get_lprob("ACT"), log(0.3) - log(0.8))
    assert_allclose(codon.get_lprob("CCC"), log(0.5) - log(0.8))

    with pytest.raises(ValueError):
        codon.set_lprob("X", 0.0)

    with pytest.raises(ValueError):
        codon.set_lprob("XX", 0.0)

    with pytest.raises(ValueError):
        codon.get_lprob("XX")

    with pytest.raises(ValueError):
        codon.set_lprob("XXX", 0.0)

    codon.set_lprob("ACT", LOG0)
    codon.set_lprob("CCC", LOG0)

    with pytest.raises(RuntimeError):
        codon.normalize()
