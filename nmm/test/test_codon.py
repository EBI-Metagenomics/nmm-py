from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, CodonTable, Codon


def test_codon():
    bases = Alphabet(b"ACGT")
    codon = CodonTable.create(bases)
    codon.set_lprob(Codon("ACT"), log(0.3))
    codon.set_lprob(Codon("CCC"), log(0.5))
    assert_allclose(codon.get_lprob(Codon("ACT")), log(0.3))
    assert_allclose(codon.get_lprob(Codon("CCC")), log(0.5))
    assert_equal(codon.get_lprob(Codon("CTC")), LOG0)
    codon.normalize()
    assert_allclose(codon.get_lprob(Codon("ACT")), log(0.3) - log(0.8))
    assert_allclose(codon.get_lprob(Codon("CCC")), log(0.5) - log(0.8))

    with pytest.raises(ValueError):
        codon.set_lprob(Codon("X"), 0.0)

    with pytest.raises(ValueError):
        codon.set_lprob(Codon("XX"), 0.0)

    with pytest.raises(ValueError):
        codon.get_lprob(Codon("XX"))

    with pytest.raises(ValueError):
        codon.set_lprob(Codon("XXX"), 0.0)

    codon.set_lprob(Codon("ACT"), LOG0)
    codon.set_lprob(Codon("CCC"), LOG0)

    with pytest.raises(RuntimeError):
        codon.normalize()
