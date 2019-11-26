from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, CodonTable


def test_codon():
    bases = Alphabet(b"ACGT")
    codon = CodonTable(bases)
    codon.set_lprob(b"ACT", log(0.3))
    codon.set_lprob(b"CCC", log(0.5))
    assert_allclose(codon.get_lprob(b"ACT"), log(0.3))
    assert_allclose(codon.get_lprob(b"CCC"), log(0.5))
    assert_equal(codon.get_lprob(b"CTC"), LOG0)
    codon.normalize()
    assert_allclose(codon.get_lprob(b"ACT"), log(0.3) - log(0.8))
    assert_allclose(codon.get_lprob(b"CCC"), log(0.5) - log(0.8))

    with pytest.raises(ValueError):
        codon.set_lprob(b"X", 0.0)

    with pytest.raises(ValueError):
        codon.set_lprob(b"XX", 0.0)

    with pytest.raises(ValueError):
        codon.get_lprob(b"XX")

    with pytest.raises(ValueError):
        codon.set_lprob(b"XXX", 0.0)

    codon.set_lprob(b"ACT", LOG0)
    codon.set_lprob(b"CCC", LOG0)

    with pytest.raises(RuntimeError):
        codon.normalize()
