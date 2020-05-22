import pytest

from nmm import BaseAlphabet, Codon, codon_iter


def test_codon():
    base = BaseAlphabet.create(b"ACGT", b"X")

    codon = Codon.create(b"AAA", base)
    assert codon.symbols == b"AAA"

    codon.symbols = b"GTX"
    assert codon.symbols == b"GTX"

    with pytest.raises(ValueError):
        codon.symbols = b"GTGG"

    with pytest.raises(ValueError):
        codon.symbols = b"GT"

    with pytest.raises(ValueError):
        codon.symbols = b"ADA"


def test_codon_iter():
    base = BaseAlphabet.create(b"ACGT", b"X")

    codons = list(codon_iter(base))
    assert len(codons) == 64
    assert codons[0].symbols == b"AAA"
    assert codons[1].symbols == b"AAC"
