import pytest
from numpy.testing import assert_equal

from nmm.alphabet import AminoAlphabet


def test_amino_alphabet():
    amino = AminoAlphabet.create(b"ACDEFGHIKLMNPQRSTVWY", b"X")

    assert_equal(amino.symbols, b"ACDEFGHIKLMNPQRSTVWY")
    assert_equal(str(amino), "{ACDEFGHIKLMNPQRSTVWY}")
    assert_equal(repr(amino), "<AminoAlphabet:{ACDEFGHIKLMNPQRSTVWY}>")

    with pytest.raises(RuntimeError):
        AminoAlphabet.create(b"ACGTK", b"X")
