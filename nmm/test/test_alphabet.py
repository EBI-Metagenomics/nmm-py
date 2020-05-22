import pytest

from nmm import AminoAlphabet, BaseAlphabet


def test_amino_alphabet():
    amino = AminoAlphabet.create(b"ACDEFGHIKLMNPQRSTVWY", b"X")

    assert amino.symbols == b"ACDEFGHIKLMNPQRSTVWY"
    assert str(amino) == "{ACDEFGHIKLMNPQRSTVWY}"
    assert repr(amino) == "<AminoAlphabet:{ACDEFGHIKLMNPQRSTVWY}>"

    with pytest.raises(RuntimeError):
        AminoAlphabet.create(b"ACGTK", b"X")


def test_base_alphabet():
    base = BaseAlphabet.create(b"ACGT", b"X")

    assert base.symbols == b"ACGT"
    assert str(base) == "{ACGT}"
    assert repr(base) == "<BaseAlphabet:{ACGT}>"

    with pytest.raises(RuntimeError):
        BaseAlphabet.create(b"ACGTK", b"X")
