from math import log

import pytest
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, CodonTable
from nmm._codon import CCodonTable
from nmm._ffi import lib, ffi


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


def test_codon_errors():

    with pytest.raises(ValueError):
        CodonTable(None)

    with pytest.raises(ValueError):
        CCodonTable(nmm_codont=None, alphabet=None)

    with pytest.raises(RuntimeError):
        CCodonTable(nmm_codont=ffi.NULL, alphabet=None)

    with pytest.raises(ValueError):
        CCodonTable(nmm_codont=ffi.NULL, alphabet=Alphabet(b"abcd"))

    alphabet = Alphabet(b"abcd")
    nmm_codont = lib.nmm_codont_create(alphabet.imm_abc)
    codont = CCodonTable(nmm_codont=nmm_codont)
    assert_equal(codont.imm_abc, alphabet.imm_abc)
