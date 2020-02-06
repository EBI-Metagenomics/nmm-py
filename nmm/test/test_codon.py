import pytest
from numpy.testing import assert_equal

from nmm import Alphabet, Codon, Base


def test_codon():
    base = Base(Alphabet(b"ACGT", b"X"))

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
