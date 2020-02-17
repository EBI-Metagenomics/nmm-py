import pytest
from numpy.testing import assert_equal

from nmm import Alphabet, Codon, BaseAlphabet


def test_codon():
    base = BaseAlphabet(Alphabet(b"ACGT", b"X"))

    codon = Codon(b"AAA", base)
    assert_equal(codon.symbols, b"AAA")

    codon.symbols = b"GTX"
    assert_equal(codon.symbols, b"GTX")

    with pytest.raises(ValueError):
        codon.symbols = b"GTGG"

    with pytest.raises(ValueError):
        codon.symbols = b"GT"

    with pytest.raises(ValueError):
        codon.symbols = b"ADA"
