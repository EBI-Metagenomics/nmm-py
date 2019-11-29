from math import log

import pytest
from numpy import nan
from numpy.testing import assert_allclose, assert_equal

from nmm import LOG0, Alphabet, BaseTable
from nmm._base import CBaseTable
from nmm._ffi import lib, ffi


def test_base():
    alphabet = Alphabet(b"ACGT")
    base = BaseTable(alphabet)
    base.set_lprob(b"A", log(0.3))
    base.set_lprob(b"T", log(0.3))

    assert_allclose(base.get_lprob(b"A"), log(0.3))
    assert_allclose(base.get_lprob(b"T"), log(0.3))

    assert_equal(base.get_lprob(b"C"), LOG0)
    assert_equal(base.get_lprob(b"G"), LOG0)

    assert_equal(base.get_lprob(b"X"), nan)

    with pytest.raises(ValueError):
        base.set_lprob(b"X", 0.0)

    with pytest.raises(ValueError):
        base.set_lprob(b"XX", 0.0)

    with pytest.raises(ValueError):
        base.get_lprob(b"XX")

    base.normalize()

    assert_allclose(base.get_lprob(b"A"), log(0.3) - log(0.6))
    assert_allclose(base.get_lprob(b"T"), log(0.3) - log(0.6))

    assert_equal(base.get_lprob(b"C"), LOG0)
    assert_equal(base.get_lprob(b"G"), LOG0)

    assert_equal(base.get_lprob(b"X"), nan)

    assert_equal(set(base.alphabet.symbols), set(b"ACGT"))

    base.set_lprob(b"A", LOG0)
    base.set_lprob(b"C", LOG0)
    base.set_lprob(b"G", LOG0)
    base.set_lprob(b"T", LOG0)

    with pytest.raises(RuntimeError):
        base.normalize()


def test_base_errors():

    with pytest.raises(ValueError):
        BaseTable(None)

    with pytest.raises(ValueError):
        CBaseTable(nmm_baset=None, alphabet=None)

    with pytest.raises(RuntimeError):
        CBaseTable(nmm_baset=ffi.NULL, alphabet=None)

    with pytest.raises(ValueError):
        CBaseTable(nmm_baset=ffi.NULL, alphabet=Alphabet(b"abcd"))

    alphabet = Alphabet(b"abcd")
    nmm_baset = lib.nmm_baset_create(alphabet.imm_abc)
    baset = CBaseTable(nmm_baset=nmm_baset)
    assert_equal(baset.imm_abc, alphabet.imm_abc)
