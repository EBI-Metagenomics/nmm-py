from math import log

import pytest
from imm.testing import assert_allclose

from nmm import BaseAlphabet, BaseTable, Codon, CodonProb, CodonTable


def test_base_table():
    base = BaseAlphabet.create(b"ACGT", b"X")
    baset = BaseTable.create(base, (log(0.1), log(0.2), log(0.3), log(0.4)))
    assert_allclose(baset.lprob(b"A"), log(0.1))
    assert_allclose(baset.lprob(b"C"), log(0.2))
    assert_allclose(baset.lprob(b"G"), log(0.3))
    assert_allclose(baset.lprob(b"T"), log(0.4))

    with pytest.raises(Exception):
        baset = BaseTable.create(base, (log(0.1), log(0.2), log(0.3)))


def test_codon_table():
    base = BaseAlphabet.create(b"ACGT", b"X")
    codonp = CodonProb.create(base)

    codonp.set_lprob(Codon.create(b"AAA", base), log(0.01))
    codonp.set_lprob(Codon.create(b"AGA", base), log(0.31))
    codonp.set_lprob(Codon.create(b"CAA", base), log(0.40))
    codonp.set_lprob(Codon.create(b"CAT", base), log(0.40))

    codont = CodonTable.create(codonp)
    assert_allclose(codont.lprob(Codon.create(b"CAT", base)), log(0.40))
    assert_allclose(codont.lprob(Codon.create(b"CAX", base)), log(0.80))
    assert_allclose(codont.lprob(Codon.create(b"XXX", base)), log(1.12))
