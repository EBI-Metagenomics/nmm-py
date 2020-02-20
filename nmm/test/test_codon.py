import pytest
from numpy.testing import assert_equal

from nmm.alphabet import BaseAlphabet
from nmm.codon import Codon


def test_codon():
    base = BaseAlphabet.create(b"ACGT", b"X")

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
